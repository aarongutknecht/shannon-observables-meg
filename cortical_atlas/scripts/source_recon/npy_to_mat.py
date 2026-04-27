#!/usr/bin/env python3
"""Convert source-reconstructed .npy outputs to .mat for MATLAB/MVGC2.

For every subject directory under derivatives/source_timeseries/sub-*/,
reads label_ts.npy, label_ts_epochs.npy, labels.txt, and meta.json, then
writes a single source_epochs.mat (v7.3 / HDF5) that MATLAB can load directly.

Usage
-----
    python scripts/source_recon/npy_to_mat.py                 # all subjects
    python scripts/source_recon/npy_to_mat.py --subjects sub-ON02747
    python scripts/source_recon/npy_to_mat.py --force          # overwrite existing .mat

The output .mat file contains:
    label_ts          (n_labels × n_times)        float64 — concatenated time series
    label_ts_epochs   (n_epochs × n_labels × T)   float64 — per-epoch time series
    labels            {n_labels × 1} cell          — region names
    sfreq             scalar                        — sampling rate (Hz)
    n_labels          scalar
    n_epochs_clean    scalar
    subject           char                          — subject ID
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

# scipy.io.savemat writes v5 by default; for large arrays use hdf5storage
# which writes -v7.3 (HDF5).  Fall back to scipy if hdf5storage is absent.
try:
    import hdf5storage
    HAS_HDF5 = True
except ImportError:
    HAS_HDF5 = False

try:
    import scipy.io as sio
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


def find_project_root() -> Path:
    """Walk upward from this script to find config.yaml."""
    d = Path(__file__).resolve().parent
    for _ in range(5):
        if (d / "config.yaml").exists():
            return d
        d = d.parent
    raise FileNotFoundError("Cannot locate project root (config.yaml)")


def convert_subject(subj_dir: Path, force: bool = False) -> bool:
    """Convert one subject's .npy files to .mat.  Returns True on success."""
    out_mat = subj_dir / "source_epochs.mat"
    if out_mat.exists() and not force:
        print(f"  SKIP {subj_dir.name} (source_epochs.mat exists, use --force)")
        return False

    # --- Load ----------------------------------------------------------------
    label_ts_path = subj_dir / "label_ts.npy"
    epochs_path   = subj_dir / "label_ts_epochs.npy"
    labels_path   = subj_dir / "labels.txt"
    meta_path     = subj_dir / "meta.json"

    if not label_ts_path.exists():
        print(f"  SKIP {subj_dir.name} (label_ts.npy not found)")
        return False

    label_ts = np.load(label_ts_path)                 # (n_labels, n_times)
    label_ts_epochs = np.load(epochs_path) if epochs_path.exists() else np.array([])
    labels = [l.strip() for l in labels_path.read_text().splitlines() if l.strip()]

    meta = {}
    if meta_path.exists():
        with open(meta_path) as f:
            meta = json.load(f)

    sfreq = float(meta.get("sfreq", 300.0))
    n_labels = int(meta.get("n_labels", label_ts.shape[0]))
    n_epochs = int(meta.get("n_epochs_clean", label_ts_epochs.shape[0] if label_ts_epochs.ndim == 3 else 0))
    subject_id = meta.get("subject", subj_dir.name)

    # --- Build dict for MATLAB -----------------------------------------------
    # MATLAB cell arrays: wrap each string in a list
    mat_labels = np.array(labels, dtype=object).reshape(-1, 1)

    mdict = {
        "label_ts":        label_ts,              # (n_labels × T)
        "label_ts_epochs": label_ts_epochs,       # (n_epochs × n_labels × T)
        "labels":          mat_labels,             # {n_labels×1} cell of strings
        "sfreq":           sfreq,
        "n_labels":        n_labels,
        "n_epochs_clean":  n_epochs,
        "subject":         subject_id,
    }

    # --- Save ----------------------------------------------------------------
    if HAS_HDF5:
        hdf5storage.savemat(str(out_mat), mdict, format='7.3',
                            oned_as='column', store_python_metadata=False)
    elif HAS_SCIPY:
        # scipy v5 .mat — works for arrays < 2 GB
        sio.savemat(str(out_mat), mdict, do_compression=True)
    else:
        print("  ERROR: neither hdf5storage nor scipy is installed", file=sys.stderr)
        return False

    mb = out_mat.stat().st_size / 1e6
    print(f"  {subj_dir.name} → source_epochs.mat  ({label_ts.shape}, {mb:.1f} MB)")
    return True


def main():
    parser = argparse.ArgumentParser(description="Convert .npy source data to .mat")
    parser.add_argument("--subjects", nargs="+", default=None,
                        help="Specific subject IDs (default: all)")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing .mat files")
    args = parser.parse_args()

    root = find_project_root()
    src_dir = root / "derivatives" / "source_timeseries"

    if not src_dir.exists():
        print(f"Source directory not found: {src_dir}", file=sys.stderr)
        sys.exit(1)

    # Discover subjects
    if args.subjects:
        subj_dirs = [src_dir / s for s in args.subjects]
    else:
        subj_dirs = sorted(d for d in src_dir.iterdir() if d.is_dir() and d.name.startswith("sub-"))

    if not subj_dirs:
        print("No subject directories found.", file=sys.stderr)
        sys.exit(1)

    print(f"Converting {len(subj_dirs)} subject(s) in {src_dir}\n")
    ok = 0
    for sd in subj_dirs:
        if convert_subject(sd, force=args.force):
            ok += 1
    print(f"\nDone: {ok}/{len(subj_dirs)} converted.")


if __name__ == "__main__":
    main()
