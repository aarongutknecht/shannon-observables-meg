"""Shared helpers for the NIMH MEG Atlas pipeline."""

from __future__ import annotations

import json
import logging
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml


# ── Config loading ──────────────────────────────────────────────────────────

def find_project_root(marker: str = "config.yaml") -> Path:
    """Walk up from CWD to find the project root (dir containing *marker*)."""
    p = Path.cwd().resolve()
    for parent in [p, *p.parents]:
        if (parent / marker).exists():
            return parent
    raise FileNotFoundError(
        f"Could not find project root (looking for '{marker}' in parents of {p})"
    )


def load_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """Load and validate the pipeline YAML config.

    Parameters
    ----------
    config_path : str | None
        Explicit path.  If *None*, searches upward for ``config.yaml``.

    Returns
    -------
    dict
        Parsed configuration with all relative paths resolved to absolute.
    """
    if config_path is None:
        root = find_project_root()
        config_path = str(root / "config.yaml")
    else:
        root = Path(config_path).resolve().parent

    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)

    cfg["_root"] = str(root)

    # Resolve relative paths
    for section_key, keys in [
        ("dataset", ["clone_dir"]),
        ("freesurfer", ["subjects_dir"]),
        ("results", ["root", "transforms_dir", "fwd_dir", "source_reconstructed_dir"]),
        ("logging", ["dir"]),
    ]:
        section = cfg.get(section_key, {})
        for k in keys:
            val = section.get(k)
            if val and not os.path.isabs(val):
                section[k] = str(root / val)

    sfile = cfg.get("subjects", {}).get("subjects_file", "subjects.txt")
    if not os.path.isabs(sfile):
        cfg["subjects"]["subjects_file"] = str(root / sfile)

    return cfg


# ── Subject list ────────────────────────────────────────────────────────────

def load_subjects(cfg: Dict[str, Any]) -> List[str]:
    """Read subject IDs from the subjects file."""
    sf = cfg["subjects"]["subjects_file"]
    if not os.path.isfile(sf):
        raise FileNotFoundError(f"Subject list not found: {sf}")
    with open(sf) as fh:
        subs = [line.strip() for line in fh if line.strip() and not line.startswith("#")]
    if not subs:
        raise ValueError(f"Subject list {sf} is empty.")
    return subs


# ── Logging ─────────────────────────────────────────────────────────────────

def setup_logging(cfg: Dict[str, Any], script_name: str) -> logging.Logger:
    """Configure a file + console logger for *script_name*.

    Log file: ``<log_dir>/<script_name>_<timestamp>.log``
    """
    log_dir = Path(cfg["logging"]["dir"])
    log_dir.mkdir(parents=True, exist_ok=True)

    ts = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    log_file = log_dir / f"{script_name}_{ts}.log"

    level = getattr(logging, cfg["logging"].get("level", "INFO").upper(), logging.INFO)

    logger = logging.getLogger(script_name)
    logger.setLevel(level)

    fmt = logging.Formatter(
        "%(asctime)s | %(name)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    fh = logging.FileHandler(log_file)
    fh.setLevel(level)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(level)
    ch.setFormatter(fmt)
    logger.addHandler(ch)

    logger.info("Log file: %s", log_file)
    return logger


# ── Path builders ───────────────────────────────────────────────────────────

def meg_ds_path(cfg: Dict[str, Any], sub: str) -> Path:
    """Return path of the CTF .ds folder for a subject."""
    ses = cfg["dataset"]["session"]
    task = cfg["dataset"]["task"]
    run = cfg["dataset"]["run"]
    ds_name = f"{sub}_{ses}_task-{task}_{run}_meg.ds"
    clone_dir = Path(cfg["dataset"]["clone_dir"])
    return clone_dir / sub / ses / "meg" / ds_name


def trans_path(cfg: Dict[str, Any], sub: str) -> Path:
    return Path(cfg["results"]["transforms_dir"]) / sub / f"{sub}-trans.fif"


def fwd_path(cfg: Dict[str, Any], sub: str) -> Path:
    return Path(cfg["results"]["fwd_dir"]) / sub / f"{sub}-fwd.fif"


# ── Idempotency helpers ────────────────────────────────────────────────────

def dry_run(cfg: Dict[str, Any]) -> bool:
    return cfg["pipeline"].get("dry_run", False)


# ── Provenance metadata ────────────────────────────────────────────────────

def save_meta(out_dir: Path, sub: str, sfreq: float, cfg: Dict[str, Any],
              extra: Optional[Dict[str, Any]] = None) -> None:
    """Write a ``meta.json`` provenance sidecar."""
    meta: Dict[str, Any] = {
        "subject": sub,
        "sfreq": sfreq,
        "pipeline_config": os.path.basename(cfg.get("_config_path", "config.yaml")),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "mne_version": _mne_version(),
    }
    if extra:
        meta.update(extra)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "meta.json", "w") as fh:
        json.dump(meta, fh, indent=2, default=_json_default)


def _json_default(obj):
    """Handle numpy types that json.dump can't serialize."""
    import numpy as np
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def _mne_version() -> str:
    try:
        import mne
        return mne.__version__
    except ImportError:
        return "N/A"


# ── CLI helpers ─────────────────────────────────────────────────────────────

def add_common_args(parser):
    """Add --config and --dry-run to an argparse parser."""
    parser.add_argument(
        "--config", "-c", default=None,
        help="Path to config.yaml (default: auto-detect from CWD parents)",
    )
    parser.add_argument(
        "--dry-run", action="store_true", default=False,
        help="Print planned actions without executing.",
    )
    parser.add_argument(
        "--force", action="store_true", default=False,
        help="Recompute even if outputs already exist.",
    )
    parser.add_argument(
        "--subjects", nargs="*", default=None,
        help="Override subject list (space-separated IDs).",
    )
    return parser


def apply_cli_overrides(cfg: Dict[str, Any], args) -> Dict[str, Any]:
    """Merge CLI flags into the config dict."""
    if args.dry_run:
        cfg["pipeline"]["dry_run"] = True
    if args.force:
        cfg["pipeline"]["skip_existing"] = False
    return cfg
