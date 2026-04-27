#!/usr/bin/env python3
"""End-to-end MEG source-reconstruction pipeline using individual anatomy.

Uses FastSurfer-reconstructed individual surfaces, BEM model, and the
``aparc.DKTatlas`` parcellation (62 cortical regions).

Coregistration relies on subject-specific ``forMEGCoregistration_T1w``
fiducials aligned to the three MEG cardinal points (nasion, LPA, RPA).

Pipeline per subject
--------------------
1. Read CTF .ds data, crop to ``max_minutes``.
2. Pick MEG + ref channels, filter (bandpass + notch), drop ref, resample.
3. Segment into fixed-length epochs; artefact rejection by peak-to-peak
   threshold.
4. Compute epoch covariance.
5. Compute BEM model + solution from individual surfaces.
6. Set up ico-4 source space on individual anatomy.
7. Compute forward solution using saved coregistration transform.
8. Build & apply epoch-based LCMV beamformer.
9. Extract ``aparc.DKTatlas`` label time courses (62 regions, PCA-flip).
10. Save ``label_ts.npy``, ``label_ts_epochs.npy``, ``labels.txt``,
    ``meta.json``.

Prerequisites
-------------
- FastSurfer recon completed (``derivatives/fastsurfer/<sub>/`` with surfaces).
- Coregistration ``<sub>-trans.fif`` saved (``run_coreg_auto.py``).

Usage
-----
    python scripts/source_recon/run_pipeline_individual.py
    python scripts/source_recon/run_pipeline_individual.py --subjects sub-ON02747
    python scripts/source_recon/run_pipeline_individual.py --dry-run
    python scripts/source_recon/run_pipeline_individual.py --force
"""

from __future__ import annotations

import argparse
import gc as gc_mod
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (
    add_common_args,
    apply_cli_overrides,
    dry_run,
    fwd_path,
    load_config,
    load_subjects,
    meg_ds_path,
    save_meta,
    setup_logging,
    trans_path,
)


# ── Preprocessing ──────────────────────────────────────────────────────────

def preprocess_raw(raw, cfg, logger):
    """Crop, pick channels, filter, notch, drop ref_meg, resample.  In-place."""
    import mne

    pp = cfg["preproc"]

    # Crop to max_minutes when configured; null means use the full recording.
    max_minutes = pp.get("max_minutes")
    max_sec = max_minutes * 60.0 if max_minutes is not None else None
    if max_sec is not None and raw.times[-1] > max_sec:
        logger.info("  Cropping to first %.0f s (of %.0f s)", max_sec, raw.times[-1])
        raw.crop(tmax=max_sec)

    # Pick channel types
    picks_kw = pp.get("pick_types", {})
    raw.pick_types(**picks_kw)
    logger.info("  Channels after pick: %d", len(raw.ch_names))

    # Bandpass
    l_freq, h_freq = pp.get("l_freq"), pp.get("h_freq")
    if l_freq or h_freq:
        logger.info("  Bandpass: %.1f – %.1f Hz", l_freq or 0, h_freq or np.inf)
        raw.filter(l_freq=l_freq, h_freq=h_freq, n_jobs=1)

    # Notch
    notch = pp.get("notch_freqs")
    if notch:
        logger.info("  Notch: %s Hz", notch)
        raw.notch_filter(freqs=notch, n_jobs=1)

    # Drop reference channels — they were needed for CTF 3rd-order
    # gradient compensation during filtering but should not enter the
    # beamformer or artefact rejection.
    ref_picks = mne.pick_types(raw.info, meg=False, ref_meg=True)
    if len(ref_picks) > 0:
        raw.drop_channels([raw.ch_names[p] for p in ref_picks])
        logger.info("  Dropped %d ref_meg channels → %d channels remain",
                     len(ref_picks), len(raw.ch_names))

    # Resample
    resample = pp.get("resample_freq")
    if resample:
        logger.info("  Resampling to %.0f Hz", resample)
        raw.resample(resample, n_jobs=1)

    return raw


def make_clean_epochs(raw, cfg, logger):
    """Segment continuous raw into fixed-length epochs and reject artefacts.

    Artefact rejection uses a coarse peak-to-peak amplitude threshold
    following the release configuration.

    Returns
    -------
    epochs : mne.Epochs
        Clean, preloaded epochs.
    qc : dict
        Epoch-level QC summary (counts and fractions kept/dropped).
    """
    import mne

    pp = cfg["preproc"]
    epoch_dur = pp.get("epoch_duration", 2.0)

    # ── Fixed-length epochs ───────────────────────────────────────────────
    events = mne.make_fixed_length_events(raw, duration=epoch_dur)
    epochs = mne.Epochs(
        raw, events,
        tmin=0.0,
        tmax=epoch_dur - 1.0 / raw.info["sfreq"],
        baseline=None,
        preload=True,
    )
    n_total = len(epochs)
    logger.info("  Created %d epochs of %.1f s", n_total, epoch_dur)

    # ── Stage 1: coarse peak-to-peak rejection ───────────────────────────
    reject_ptp = pp.get("reject_ptp")
    if reject_ptp:
        reject_dict = {k: float(v) for k, v in reject_ptp.items()}
        logger.info("  Peak-to-peak rejection thresholds: %s", reject_dict)
        epochs.drop_bad(reject=reject_dict)
        n_after_ptp = len(epochs)
        logger.info("  After PTP: %d / %d kept (%.0f %% dropped)",
                     n_after_ptp, n_total,
                     100 * (1 - n_after_ptp / n_total))
    else:
        n_after_ptp = n_total

    n_final = len(epochs)
    n_dropped_ptp = n_total - n_after_ptp
    n_dropped_total = n_total - n_final
    frac_kept = (n_final / n_total) if n_total else 0.0
    frac_dropped = 1.0 - frac_kept

    logger.info("  Final: %d clean epochs (%.0f s of usable data)",
                n_final, n_final * epoch_dur)

    qc = {
        "n_epochs_total": int(n_total),
        "n_epochs_after_ptp": int(n_after_ptp),
        "n_epochs_clean": int(n_final),
        "n_epochs_dropped_ptp": int(n_dropped_ptp),
        "n_epochs_dropped_total": int(n_dropped_total),
        "frac_epochs_kept": float(frac_kept),
        "frac_epochs_dropped_total": float(frac_dropped),
        "usable_duration_s": float(n_final * epoch_dur),
    }
    return epochs, qc


# ── Per-subject processing ─────────────────────────────────────────────────

def process_subject(sub: str, cfg: dict, logger) -> bool:
    """Full individual-anatomy pipeline for one subject (epoch-based).

    Returns True on success, False on failure.
    """
    import mne
    from mne.beamformer import apply_lcmv_epochs, make_lcmv

    src_cfg = cfg["source"]
    pp = cfg["preproc"]
    subjects_dir = cfg["freesurfer"]["subjects_dir"]
    parc = "aparc.DKTatlas"   # FastSurfer native parcellation (62 regions)
    t0 = time.time()

    # ── Output paths ──────────────────────────────────────────────────────
    results_root = Path(cfg["results"]["source_reconstructed_dir"])
    out_dir = results_root / sub
    ts_file = out_dir / "label_ts.npy"

    if ts_file.exists() and cfg["pipeline"]["skip_existing"]:
        logger.info("SKIP %s — %s exists.", sub, ts_file)
        return True

    # ── Check prerequisites ───────────────────────────────────────────────
    meg = meg_ds_path(cfg, sub)
    tp = trans_path(cfg, sub)
    fs_sub = Path(subjects_dir) / sub

    missing = []
    if not meg.exists():
        missing.append(f"MEG data: {meg}")
    if not tp.exists():
        missing.append(f"Transform: {tp}  (run run_coreg_auto.py first)")
    if not (fs_sub / "surf" / "lh.white").exists():
        missing.append(f"FreeSurfer surfaces: {fs_sub}/surf/")
    annot = fs_sub / "label" / f"lh.{parc}.annot"
    if not annot.exists():
        missing.append(f"Parcellation: {annot}")
    if missing:
        for m in missing:
            logger.error("MISSING %s", m)
        logger.error("Skipping %s.", sub)
        return False

    # ── 1. Read MEG ──────────────────────────────────────────────────────
    logger.info("[%s] Reading CTF: %s", sub, meg)
    raw = mne.io.read_raw_ctf(str(meg), preload=True, system_clock="ignore")
    logger.info("[%s] Raw: %d ch, %.0f s, sfreq=%.0f Hz",
                sub, len(raw.ch_names), raw.times[-1], raw.info["sfreq"])

    # ── 2. Preprocess (filter, resample) ─────────────────────────────────
    raw = preprocess_raw(raw, cfg, logger)

    # ── 3. Epoch + artefact rejection ────────────────────────────────────
    epochs, epoch_qc = make_clean_epochs(raw, cfg, logger)

    # Save sfreq before deleting raw
    sfreq = epochs.info["sfreq"]
    info = epochs.info.copy()
    del raw
    gc_mod.collect()

    # ── 4. Covariance from clean epochs ──────────────────────────────────
    logger.info("[%s] Computing covariance (method=%s) from %d clean epochs …",
                sub, src_cfg["cov_method"], len(epochs))
    data_cov = mne.compute_covariance(
        epochs, method=src_cfg["cov_method"], rank="info",
    )
    logger.info("[%s] Covariance: %s", sub, data_cov.data.shape)

    # ── 5. BEM surfaces + model + solution ───────────────────────────────
    #
    # FastSurfer does not produce BEM surfaces (inner_skull etc.).
    # We use FreeSurfer's mri_watershed via MNE's make_watershed_bem()
    # to create them from the T1 volume.
    #
    # For MEG-only data, a single-shell BEM (inner skull only) is
    # standard and well-validated (Hämäläinen & Sarvas 1989).
    # Conductivity (0.3,) → single compartment.
    #
    inner_skull = Path(subjects_dir) / sub / "bem" / "inner_skull.surf"
    if not inner_skull.exists():
        logger.info("[%s] BEM surfaces missing — running make_watershed_bem …", sub)
        mne.bem.make_watershed_bem(
            subject=sub,
            subjects_dir=subjects_dir,
            overwrite=True,
        )
        if not inner_skull.exists():
            logger.error("[%s] make_watershed_bem did not produce %s", sub, inner_skull)
            return False
        logger.info("[%s] Watershed BEM surfaces created.", sub)
    else:
        logger.info("[%s] BEM surfaces already exist.", sub)

    fwd_dir = Path(cfg["results"]["fwd_dir"]) / sub
    fwd_dir.mkdir(parents=True, exist_ok=True)
    bem_sol_file = fwd_dir / "bem-sol.fif"

    if bem_sol_file.exists() and cfg["pipeline"]["skip_existing"]:
        logger.info("[%s] Reading cached BEM: %s", sub, bem_sol_file)
        bem_sol = mne.read_bem_solution(str(bem_sol_file))
    else:
        # Single-shell BEM for MEG-only data
        conductivity = (0.3,)
        ico = src_cfg.get("bem_ico", 4)
        logger.info("[%s] Making single-shell BEM model (conductivity=%s, ico=%d) …",
                     sub, conductivity, ico)
        model = mne.make_bem_model(
            subject=sub,
            subjects_dir=subjects_dir,
            conductivity=conductivity,
            ico=ico,
        )
        logger.info("[%s] Making BEM solution …", sub)
        bem_sol = mne.make_bem_solution(model)
        mne.write_bem_solution(str(bem_sol_file), bem_sol, overwrite=True)
        logger.info("[%s] Wrote %s", sub, bem_sol_file)
        del model

    # ── 6. Source space ──────────────────────────────────────────────────
    src_file = fwd_dir / "src.fif"
    spacing = src_cfg.get("spacing", "ico4")

    if src_file.exists() and cfg["pipeline"]["skip_existing"]:
        logger.info("[%s] Reading cached source space: %s", sub, src_file)
        src = mne.read_source_spaces(str(src_file))
    else:
        logger.info("[%s] Setting up source space (spacing=%s) …", sub, spacing)
        src = mne.setup_source_space(
            subject=sub,
            subjects_dir=subjects_dir,
            spacing=spacing,
            add_dist="patch",
        )
        mne.write_source_spaces(str(src_file), src, overwrite=True)
        logger.info("[%s] Wrote %s", sub, src_file)

    n_src = sum(s["nuse"] for s in src)
    logger.info("[%s] Source space: %d active vertices", sub, n_src)

    # ── 7. Forward solution ──────────────────────────────────────────────
    fp = fwd_path(cfg, sub)

    if fp.exists() and cfg["pipeline"]["skip_existing"]:
        logger.info("[%s] Reading cached forward: %s", sub, fp)
        fwd = mne.read_forward_solution(str(fp))
    else:
        trans = mne.read_trans(str(tp))
        logger.info("[%s] Computing forward solution (mindist=%.1f) …",
                     sub, src_cfg["mindist"])
        fwd = mne.make_forward_solution(
            info,
            trans=trans,
            src=src,
            bem=bem_sol,
            mindist=src_cfg["mindist"],
            eeg=False,
        )
        mne.write_forward_solution(str(fp), fwd, overwrite=True)
        logger.info("[%s] Wrote %s", sub, fp)

    n_fwd_src = fwd["nsource"]
    logger.info("[%s] Forward: %d sources × %d sensors",
                sub, n_fwd_src, fwd["nchan"])
    del bem_sol
    gc_mod.collect()

    # ── 8. LCMV beamformer ───────────────────────────────────────────────
    reg = src_cfg.get("reg", 0.05)
    logger.info("[%s] Building LCMV (ori=%s, norm=%s, reg=%.2f) …",
                sub, src_cfg["pick_ori"], src_cfg["weight_norm"], reg)
    filters = make_lcmv(
        info,
        fwd,
        data_cov,
        reg=reg,
        pick_ori=src_cfg["pick_ori"],
        weight_norm=src_cfg["weight_norm"],
    )
    del data_cov
    gc_mod.collect()

    # ── 9. Apply LCMV to epochs ──────────────────────────────────────────
    n_epo = len(epochs)
    logger.info("[%s] Applying LCMV to %d epochs …", sub, n_epo)
    stcs = apply_lcmv_epochs(epochs, filters)
    del epochs, filters, fwd
    gc_mod.collect()

    logger.info("[%s] STC per epoch: %d vertices × %d time points",
                sub, stcs[0].data.shape[0], stcs[0].data.shape[1])

    # ── 10. Label time courses ───────────────────────────────────────────
    mode = src_cfg["label_mode"]

    logger.info("[%s] Reading '%s' labels from %s …", sub, parc, subjects_dir)
    labels = mne.read_labels_from_annot(
        sub, parc=parc, subjects_dir=subjects_dir,
    )
    labels = [l for l in labels if "unknown" not in l.name.lower()]
    n_labels = len(labels)
    logger.info("[%s] %d labels (excl. unknown)", sub, n_labels)

    logger.info("[%s] Extracting label time courses (mode=%s) from %d epoch STCs …",
                sub, mode, len(stcs))
    epoch_ts_list = mne.extract_label_time_course(
        stcs, labels, src, mode=mode, allow_empty=True,
    )
    epoch_ts = np.array(epoch_ts_list)  # (n_epochs, n_labels, n_times_epoch)
    del stcs, epoch_ts_list
    gc_mod.collect()

    n_times_epoch = epoch_ts.shape[2]

    # Concatenated view: (n_labels, n_epochs × n_times_epoch)
    label_ts = epoch_ts.transpose(1, 0, 2).reshape(n_labels, -1)
    logger.info("[%s] label_ts concatenated: %s  (labels × time)", sub, label_ts.shape)

    # ── 11. Save outputs ─────────────────────────────────────────────────
    out_dir.mkdir(parents=True, exist_ok=True)

    # Concatenated time series
    np.save(str(ts_file), label_ts)
    logger.info("[%s] Saved %s", sub, ts_file)

    # Per-epoch time series
    epoch_ts_file = out_dir / "label_ts_epochs.npy"
    np.save(str(epoch_ts_file), epoch_ts)
    logger.info("[%s] Saved %s  shape=%s", sub, epoch_ts_file, epoch_ts.shape)

    # Label names
    labels_file = out_dir / "labels.txt"
    with open(labels_file, "w") as fh:
        for l in labels:
            fh.write(l.name + "\n")

    # Metadata / provenance
    epoch_dur = pp.get("epoch_duration", 2.0)
    save_meta(out_dir, sub, sfreq, cfg, extra={
        "template": "individual",
        "anatomy": "fastsurfer",
        "parc": parc,
        "n_labels": n_labels,
        "n_sources": n_src,
        "n_epochs_clean": int(epoch_qc["n_epochs_clean"]),
        "n_epochs_total": int(epoch_qc["n_epochs_total"]),
        "n_epochs_after_ptp": int(epoch_qc["n_epochs_after_ptp"]),
        "n_epochs_dropped_ptp": int(epoch_qc["n_epochs_dropped_ptp"]),
        "n_epochs_dropped_total": int(epoch_qc["n_epochs_dropped_total"]),
        "frac_epochs_kept": float(epoch_qc["frac_epochs_kept"]),
        "frac_epochs_dropped_total": float(epoch_qc["frac_epochs_dropped_total"]),
        "usable_duration_s": float(epoch_qc["usable_duration_s"]),
        "n_times_per_epoch": int(n_times_epoch),
        "n_times_total": int(label_ts.shape[1]),
        "epoch_duration_s": epoch_dur,
        "label_mode": mode,
        "beamformer": "lcmv",
        "reg": reg,
        "pick_ori": src_cfg["pick_ori"],
        "weight_norm": src_cfg["weight_norm"],
        "preproc": cfg["preproc"],
    })

    elapsed = time.time() - t0
    logger.info("[%s] ✓ Done in %.1f min.", sub, elapsed / 60)
    return True


# ── CLI ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_common_args(parser)
    args = parser.parse_args()

    cfg = load_config(args.config)
    cfg = apply_cli_overrides(cfg, args)
    logger = setup_logging(cfg, "run_pipeline_individual")

    subjects = args.subjects or load_subjects(cfg)
    is_dry = dry_run(cfg)

    ok, fail = 0, 0
    for i, sub in enumerate(subjects, 1):
        logger.info("━━━ [%d/%d] %s ━━━", i, len(subjects), sub)
        if is_dry:
            meg = meg_ds_path(cfg, sub)
            tp = trans_path(cfg, sub)
            fs_sub = Path(cfg["freesurfer"]["subjects_dir"]) / sub
            logger.info("[DRY-RUN] Would run individual-anatomy pipeline for %s", sub)
            logger.info("  MEG: %s (exists=%s)", meg, meg.exists())
            logger.info("  trans: %s (exists=%s)", tp, tp.exists())
            logger.info("  FreeSurfer: %s (exists=%s)", fs_sub,
                        (fs_sub / "surf" / "lh.white").exists())
            continue

        try:
            if process_subject(sub, cfg, logger):
                ok += 1
            else:
                fail += 1
        except Exception:
            logger.exception("ERROR processing %s", sub)
            fail += 1

    if not is_dry:
        logger.info("═══ Summary: OK=%d, FAILED=%d ═══", ok, fail)
        if fail:
            sys.exit(1)


if __name__ == "__main__":
    main()
