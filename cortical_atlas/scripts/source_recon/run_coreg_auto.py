#!/usr/bin/env python3
"""Automated MEG-MRI coregistration.

Computes a rigid-body (head -> MRI) transform for each subject using
MNE's standard ``coregister_fiducials()`` function.

MRI-side fiducials are subject-specific landmarks from a
``forMEGCoregistration_T1w`` scan (BIDS ``acq-forMEGCoregistration``).
These provide per-subject anatomical landmark coordinates (LPA, nasion,
RPA) digitised on the MRI via BrainSight neuro-navigation. The coordinates
are stored in voxel space of the localiser image and are transformed to
FreeSurfer surface-RAS via the NIfTI affine and the FreeSurfer ``orig.mgz``
header.

The MEG-side fiducials come from the CTF digitisation: the three HPI coil
positions, which define the MEG head coordinate frame.

Pipeline
--------
1. Read MEG ``info`` -> extract cardinal points (head coords).
2. Obtain subject-specific MRI-side fiducials.
3. Compute the head->MRI transform via ``mne.coreg.coregister_fiducials()``.
4. Sanity-check fiducial geometry and fit residuals.
5. Save ``<sub>-trans.fif`` and a QC JSON report.

Usage
-----
    python scripts/source_recon/run_coreg_auto.py
    python scripts/source_recon/run_coreg_auto.py --subjects sub-ON02747
    python scripts/source_recon/run_coreg_auto.py --dry-run
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (
    add_common_args,
    apply_cli_overrides,
    dry_run,
    load_config,
    load_subjects,
    meg_ds_path,
    setup_logging,
    trans_path,
)


# == MRI-side fiducials ======================================================

def _make_fiducials_list(lpa, nasion, rpa):
    """Build an MNE fiducials list in the MRI coordinate frame.

    Each entry mirrors what ``mne.coreg.get_mni_fiducials()`` returns:
    a dict with keys *kind*, *ident*, *r* (metres), *coord_frame*.

    Parameters
    ----------
    lpa, nasion, rpa : array-like, shape (3,)
        Positions in FreeSurfer surface-RAS (metres).
    """
    def _pt(ident, r):
        return dict(
            kind=1,       # FIFFV_POINT_CARDINAL
            ident=ident,
            r=np.asarray(r, dtype=np.float64),
            coord_frame=5,  # FIFFV_COORD_MRI
        )

    return [_pt(1, lpa), _pt(2, nasion), _pt(3, rpa)]


def get_subject_specific_fiducials(sub, cfg, logger):
    """Load fiducials from a ``forMEGCoregistration_T1w`` sidecar.

    The JSON contains ``AnatomicalLandmarkCoordinates`` with LPA, NAS,
    RPA in **voxel** coordinates of the forMEGCoregistration image.
    These are transformed to FreeSurfer surface-RAS via::

        forMEGCoreg_voxel  --(NIfTI affine)-->  scanner RAS
        scanner RAS  --(Torig @ inv(Norig))-->  surface RAS

    Returns
    -------
    fiducials : list of dict | None
        MNE-format fiducials in FreeSurfer surface-RAS (metres), or
        ``None`` if the forMEGCoregistration data is not available.
    """
    import nibabel as nib

    ses = cfg["dataset"]["session"]
    clone_dir = Path(cfg["dataset"]["clone_dir"])

    # -- JSON sidecar --
    json_name = f"{sub}_{ses}_acq-forMEGCoregistration_T1w.json"
    json_path = clone_dir / sub / ses / "anat" / json_name
    if not json_path.exists():
        return None

    with open(json_path) as fh:
        meta = json.load(fh)

    coords = meta.get("AnatomicalLandmarkCoordinates")
    if coords is None:
        logger.warning("  forMEGCoreg JSON has no AnatomicalLandmarkCoordinates")
        return None

    for key in ("LPA", "NAS", "RPA"):
        if key not in coords:
            logger.warning("  forMEGCoreg JSON missing '%s'", key)
            return None

    lpa_vox = np.array(coords["LPA"], dtype=np.float64)
    nas_vox = np.array(coords["NAS"], dtype=np.float64)
    rpa_vox = np.array(coords["RPA"], dtype=np.float64)

    logger.info("  forMEGCoreg voxel coords: LPA=%s, NAS=%s, RPA=%s",
                lpa_vox.round(1), nas_vox.round(1), rpa_vox.round(1))

    # -- NIfTI affine (forMEGCoreg image) --
    nii_name = f"{sub}_{ses}_acq-forMEGCoregistration_T1w.nii.gz"
    nii_path = clone_dir / sub / ses / "anat" / nii_name
    if not nii_path.exists():
        logger.warning("  forMEGCoreg NIfTI not available: %s", nii_path)
        return None

    try:
        fmc_img = nib.load(str(nii_path))
    except Exception as exc:
        logger.warning("  Cannot load forMEGCoreg NIfTI: %s", exc)
        return None

    fmc_affine = fmc_img.affine  # voxel -> scanner RAS (mm)

    # -- FreeSurfer surface-RAS transform (from MPRAGE orig.mgz) --
    subjects_dir = Path(cfg["freesurfer"]["subjects_dir"])
    orig_path = subjects_dir / sub / "mri" / "orig.mgz"
    if not orig_path.exists():
        logger.warning("  FreeSurfer orig.mgz not found: %s", orig_path)
        return None

    orig_img = nib.load(str(orig_path))
    Norig = orig_img.header.get_vox2ras()       # MPRAGE voxel -> scanner RAS
    Torig = orig_img.header.get_vox2ras_tkr()   # MPRAGE voxel -> surface RAS
    ras_to_surfras = Torig @ np.linalg.inv(Norig)

    # -- Transform each fiducial --
    def vox_to_surfras_m(vox):
        """forMEGCoreg voxel -> surface RAS (metres)."""
        scanner_ras = fmc_affine @ np.append(vox, 1.0)      # mm
        surf_ras = ras_to_surfras @ scanner_ras               # mm
        return surf_ras[:3] / 1000.0                           # -> metres

    lpa_m = vox_to_surfras_m(lpa_vox)
    nas_m = vox_to_surfras_m(nas_vox)
    rpa_m = vox_to_surfras_m(rpa_vox)

    logger.info("  Surface-RAS (m): LPA=[%.4f, %.4f, %.4f], "
                "NAS=[%.4f, %.4f, %.4f], RPA=[%.4f, %.4f, %.4f]",
                *lpa_m, *nas_m, *rpa_m)

    return _make_fiducials_list(lpa_m, nas_m, rpa_m)


# == Sanity checks ============================================================

def check_fiducial_geometry(fids_dict, label, logger):
    """Check that fiducial positions are anatomically plausible.

    Parameters
    ----------
    fids_dict : dict with 'nasion', 'lpa', 'rpa' (in metres)
    label : str
    logger : logging.Logger

    Returns
    -------
    warnings : list of str
    """
    warnings = []
    lpa, rpa, nas = fids_dict["lpa"], fids_dict["rpa"], fids_dict["nasion"]

    # LPA<->RPA distance (head width) ~ 13-18 cm
    lr_dist = np.linalg.norm(lpa - rpa)
    logger.info("  %s LPA-RPA distance: %.1f cm", label, lr_dist * 100)
    if not (0.10 <= lr_dist <= 0.22):
        w = f"{label} LPA-RPA distance {lr_dist*100:.1f} cm outside [10, 22] cm"
        warnings.append(w)
        logger.warning("  WARNING %s", w)

    # Nasion depth (nasion -> midpoint of LPA/RPA) ~ 7-13 cm
    midpoint = (lpa + rpa) / 2
    nas_depth = np.linalg.norm(nas - midpoint)
    logger.info("  %s nasion depth: %.1f cm", label, nas_depth * 100)
    if not (0.05 <= nas_depth <= 0.16):
        w = f"{label} nasion depth {nas_depth*100:.1f} cm outside [5, 16] cm"
        warnings.append(w)
        logger.warning("  WARNING %s", w)

    return warnings


def extract_meg_fids_dict(info):
    """Extract nasion/lpa/rpa from MEG dig points as a dict of arrays (metres)."""
    IDENT_MAP = {1: "lpa", 2: "nasion", 3: "rpa"}
    fids = {}
    for d in info.get("dig", []):
        if d["kind"] == 1:  # FIFFV_POINT_CARDINAL
            name = IDENT_MAP.get(int(d["ident"]))
            if name:
                fids[name] = np.array(d["r"], dtype=float)

    for name in ("nasion", "lpa", "rpa"):
        if name not in fids:
            raise ValueError(f"MEG fiducial '{name}' not found in digitisation")
        if np.any(np.isnan(fids[name])):
            raise ValueError(f"MEG fiducial '{name}' contains NaN")
        if np.linalg.norm(fids[name]) < 1e-6:
            raise ValueError(f"MEG fiducial '{name}' is at origin")

    return fids


def fids_list_to_dict(fids_list):
    """Convert MNE fiducials list to {nasion, lpa, rpa} dict."""
    IDENT_MAP = {1: "lpa", 2: "nasion", 3: "rpa"}
    return {IDENT_MAP[d["ident"]]: np.array(d["r"]) for d in fids_list}


def compute_residuals(trans_obj, meg_fids_dict, mri_fids_list):
    """Per-fiducial distance (metres) after applying the transform."""
    import mne

    mri_dict = fids_list_to_dict(mri_fids_list)
    residuals = {}
    for name in ("nasion", "lpa", "rpa"):
        meg_in_mri = mne.transforms.apply_trans(trans_obj, meg_fids_dict[name])
        residuals[name] = float(np.linalg.norm(meg_in_mri - mri_dict[name]))

    return residuals


# == QC figure =================================================================

def save_qc_figure(subject, subjects_dir, info, trans_obj, out_path, logger):
    """Save a coregistration QC figure (best-effort, headless-safe)."""
    try:
        import mne
        import matplotlib
        matplotlib.use("Agg")

        # Force off-screen rendering so PyVista/VTK never opens an X window
        try:
            import pyvista
            pyvista.OFF_SCREEN = True
        except ImportError:
            pass

        fig = mne.viz.plot_alignment(
            info,
            trans=trans_obj,
            subject=subject,
            subjects_dir=subjects_dir,
            surfaces=["head-dense"],
            meg="sensors",
            dig=True,
            coord_frame="mri",
            verbose=False,
        )

        try:
            img = fig.plotter.screenshot()
            import matplotlib.pyplot as plt
            fig_mpl, ax = plt.subplots(1, 1, figsize=(8, 6))
            ax.imshow(img)
            ax.axis("off")
            ax.set_title(f"Coregistration QC - {subject}")
            fig_mpl.savefig(str(out_path), dpi=150, bbox_inches="tight")
            plt.close(fig_mpl)
            fig.plotter.close()
            logger.info("  QC figure saved: %s", out_path)
        except Exception as exc:
            logger.warning("  Could not save QC screenshot: %s", exc)
            try:
                fig.plotter.close()
            except Exception:
                pass
    except Exception as exc:
        logger.warning("  QC figure generation failed: %s (non-fatal)", exc)


# == Per-subject processing ====================================================

def process_subject(sub: str, cfg: dict, logger) -> bool:
    """Automated coregistration for one subject."""
    import mne

    t0 = time.time()
    subjects_dir = cfg["freesurfer"]["subjects_dir"]
    tp = trans_path(cfg, sub)

    # Skip if done
    if tp.exists() and cfg["pipeline"]["skip_existing"]:
        logger.info("SKIP %s - %s exists.", sub, tp)
        return True

    # -- Check prerequisites ------------------------------------------------
    fs_sub_dir = Path(subjects_dir) / sub
    if not (fs_sub_dir / "surf" / "lh.white").exists():
        logger.error("FreeSurfer surfaces not found for %s in %s", sub, fs_sub_dir)
        return False

    meg = meg_ds_path(cfg, sub)
    if not meg.exists():
        logger.error("MEG data not found: %s", meg)
        return False

    # -- 1. Read MEG info ---------------------------------------------------
    logger.info("[%s] Reading MEG info: %s", sub, meg)
    raw = mne.io.read_raw_ctf(str(meg), preload=False,
                               system_clock="ignore", verbose=False)
    info = raw.info
    del raw

    try:
        meg_fids_dict = extract_meg_fids_dict(info)
    except ValueError as exc:
        logger.error("[%s] MEG fiducial extraction failed: %s", sub, exc)
        return False

    for name, pos in meg_fids_dict.items():
        logger.info("  MEG %-7s [%.4f, %.4f, %.4f] m", name, *pos)

    # -- 2. Obtain MRI-side fiducials (prefer subject-specific) -------------
    logger.info("[%s] Looking for subject-specific fiducials ...", sub)
    mri_fids = get_subject_specific_fiducials(sub, cfg, logger)

    if mri_fids is None:
        logger.error("[%s] forMEGCoregistration fiducials not found in %s — "
                     "cannot proceed without subject-specific fiducials.",
                     sub, cfg["dataset"]["clone_dir"])
        return False
    fid_method = "subject_specific_forMEGCoregistration"
    logger.info("[%s] Using subject-specific fiducials (forMEGCoregistration)", sub)

    mri_fids_dict = fids_list_to_dict(mri_fids)
    for name, pos in mri_fids_dict.items():
        logger.info("  MRI %-7s [%.4f, %.4f, %.4f] m", name, *pos)

    # -- 3. Sanity-check fiducial geometry ----------------------------------
    logger.info("[%s] Checking fiducial geometry ...", sub)
    warnings = []
    warnings += check_fiducial_geometry(meg_fids_dict, "MEG", logger)
    warnings += check_fiducial_geometry(mri_fids_dict, "MRI", logger)

    # -- 4. Compute head->MRI transform ------------------------------------
    logger.info("[%s] Computing head->MRI transform (coregister_fiducials) ...", sub)
    trans_obj = mne.coreg.coregister_fiducials(info, mri_fids)

    # -- 5. Verify fit residuals -------------------------------------------
    residuals = compute_residuals(trans_obj, meg_fids_dict, mri_fids)
    residuals_mm = {k: round(v * 1000, 4) for k, v in residuals.items()}
    max_resid = max(residuals.values())

    logger.info("  Fit residuals: %s",
                {k: f"{v:.4f} mm" for k, v in residuals_mm.items()})

    if max_resid > 0.005:  # > 5 mm
        w = f"Max fiducial residual {max_resid*1000:.2f} mm > 5 mm"
        warnings.append(w)
        logger.warning("  WARNING %s", w)

    if warnings:
        logger.warning("  Sanity checks: %d warning(s)", len(warnings))
    else:
        logger.info("  Sanity checks: all passed")

    # -- 6. Save transform --------------------------------------------------
    tp.parent.mkdir(parents=True, exist_ok=True)
    mne.write_trans(str(tp), trans_obj, overwrite=True)
    logger.info("[%s] Saved transform: %s", sub, tp)

    # -- 7. Save QC report --------------------------------------------------
    qc_report = {
        "subject": sub,
        "method": fid_method,
        "meg_fiducials_m": {k: v.tolist() for k, v in meg_fids_dict.items()},
        "mri_fiducials_m": {k: v.tolist() for k, v in mri_fids_dict.items()},
        "residuals_mm": residuals_mm,
        "warnings": warnings,
    }
    qc_file = tp.parent / "coreg_qc.json"
    with open(qc_file, "w") as fh:
        json.dump(qc_report, fh, indent=2)
    logger.info("[%s] Saved QC report: %s", sub, qc_file)

    # -- 8. QC figure (best-effort) -----------------------------------------
    try:
        import os
        if os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"):
            qc_fig_path = tp.parent / "coreg_qc.png"
            save_qc_figure(sub, subjects_dir, info, trans_obj, qc_fig_path, logger)
        else:
            logger.info("  Skipping QC figure (no display).")
    except Exception as exc:
        logger.warning("  QC figure skipped: %s", exc)

    elapsed = time.time() - t0
    logger.info("[%s] Coregistration done in %.1f s (method: %s)", sub, elapsed, fid_method)
    return True


# == CLI =======================================================================

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    add_common_args(parser)
    args = parser.parse_args()

    cfg = load_config(args.config)
    cfg = apply_cli_overrides(cfg, args)
    logger = setup_logging(cfg, "run_coreg_auto")

    subjects = args.subjects or load_subjects(cfg)
    is_dry = dry_run(cfg)

    ok, fail = 0, 0
    for i, sub in enumerate(subjects, 1):
        logger.info("=== [%d/%d] %s ===", i, len(subjects), sub)

        if is_dry:
            tp = trans_path(cfg, sub)
            meg = meg_ds_path(cfg, sub)
            fs_sub = Path(cfg["freesurfer"]["subjects_dir"]) / sub
            logger.info("[DRY-RUN] Would compute coreg for %s", sub)
            logger.info("  MEG: %s (exists=%s)", meg, meg.exists())
            logger.info("  FreeSurfer: %s (exists=%s)", fs_sub,
                         (fs_sub / "surf" / "lh.white").exists())
            logger.info("  Output: %s (exists=%s)", tp, tp.exists())
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
        logger.info("=== Summary: OK=%d, FAILED=%d ===", ok, fail)
        if fail:
            sys.exit(1)


if __name__ == "__main__":
    main()
