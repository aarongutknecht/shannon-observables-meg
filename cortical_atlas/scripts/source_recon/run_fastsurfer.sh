#!/usr/bin/env bash
# ===========================================================================
# Run FastSurfer (Docker, CPU) for study subjects
# ===========================================================================
#
# Usage:
#   bash scripts/source_recon/run_fastsurfer.sh                  # all subjects
#   bash scripts/source_recon/run_fastsurfer.sh sub-ON08710      # single subject
#
# Prerequisites:
#   - Docker installed, image deepmi/fastsurfer:latest pulled
#   - FreeSurfer license path exported as $FREESURFER_LICENSE
#   - T1w NIfTI in data/ds005752/<sub>/ses-01/anat/
#
# Output:
#   derivatives/fastsurfer/<sub>/  — full FreeSurfer-compatible directory structure
#
# FastSurfer on CPU with 16 threads: ~1.5–2.5 h per subject
# ===========================================================================

set -euo pipefail

# ── Project root (directory containing config.yaml) ────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# ── Configuration ──────────────────────────────────────────────────────────
DOCKER_IMAGE="deepmi/fastsurfer:latest"
LICENSE_FILE="${FREESURFER_LICENSE:-}"
RAW_DIR="$PROJECT_ROOT/data/ds005752"
SUBJECTS_DIR="$PROJECT_ROOT/derivatives/fastsurfer"
LOG_DIR="$PROJECT_ROOT/logs"
SUBJECT_LIST="$PROJECT_ROOT/subjects.txt"
THREADS=16   # match CPU thread count

mkdir -p "$SUBJECTS_DIR" "$LOG_DIR"

# ── Subject list ───────────────────────────────────────────────────────────
if [[ $# -gt 0 ]]; then
    SUBJECTS=("$@")
else
    if [[ ! -f "$SUBJECT_LIST" ]]; then
        echo "ERROR: Subject list not found: $SUBJECT_LIST" >&2
        exit 1
    fi
    mapfile -t SUBJECTS < <(grep -v '^\s*#' "$SUBJECT_LIST" | grep -v '^\s*$')
fi

echo "======================================================================"
echo " FastSurfer (Docker CPU) – ${#SUBJECTS[@]} subject(s)"
echo " Image:       $DOCKER_IMAGE"
echo " License:     $LICENSE_FILE"
echo " Output:      $SUBJECTS_DIR"
echo " Threads:     $THREADS"
echo "======================================================================"

# ── Sanity checks ─────────────────────────────────────────────────────────
if [[ -z "$LICENSE_FILE" ]] || [[ ! -f "$LICENSE_FILE" ]]; then
    echo "ERROR: set FREESURFER_LICENSE to a valid FreeSurfer license file." >&2
    echo "  Get one free at: https://surfer.nmr.mgh.harvard.edu/fswiki/License"
    exit 1
fi

if ! sudo docker image inspect "$DOCKER_IMAGE" &>/dev/null; then
    echo "ERROR: Docker image not found: $DOCKER_IMAGE" >&2
    echo "  Pull it:  sudo docker pull $DOCKER_IMAGE"
    exit 1
fi

# ── Process each subject ──────────────────────────────────────────────────
PASS=0
FAIL=0
SKIP=0

for SUB in "${SUBJECTS[@]}"; do
    echo ""
    echo "━━━ $SUB ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

    # Check if already done
    if [[ -f "$SUBJECTS_DIR/$SUB/scripts/recon-all.done" ]] || \
       [[ -f "$SUBJECTS_DIR/$SUB/scripts/recon-surf.done" ]]; then
        echo "SKIP $SUB — already completed."
        (( SKIP++ )) || true
        continue
    fi

    # Locate T1w (MPRAGE) — deterministic BIDS path
    SES="ses-01"
    T1W_DIR="$RAW_DIR/$SUB/$SES/anat"
    T1W_BASENAME="${SUB}_${SES}_acq-MPRAGE_T1w.nii.gz"
    T1W="$T1W_DIR/$T1W_BASENAME"
    if [[ ! -f "$T1W" ]]; then
        echo "ERROR: MPRAGE T1w not found: $T1W" >&2
        ((FAIL++))
        continue
    fi
    echo "T1w: $T1W"

    LOGFILE="$LOG_DIR/fastsurfer_${SUB}_$(date +%Y%m%dT%H%M%S).log"
    echo "Log: $LOGFILE"
    echo "Starting FastSurfer (CPU, $THREADS threads)..."

    # ── Docker run ─────────────────────────────────────────────────────
    # Mount points:
    #   /data      → directory containing the T1w file
    #   /output    → FreeSurfer subjects_dir
    #   /fs_license → directory containing license.txt
    #
    # Flags:
    #   --device cpu      — no GPU
    #   --threads $THREADS — parallel threads for recon-surf
    #   --seg_only is NOT set → runs full pipeline (seg + recon-surf)
    #   --3T               — data acquired on 3T scanner (NIMH MEG dataset)
    #   Hemispheres auto-parallel when --threads > 1 (FastSurfer ≥2.4)
    sudo docker run --rm \
        --user "$(id -u):$(id -g)" \
        -v "$T1W_DIR":/data:ro \
        -v "$SUBJECTS_DIR":/output \
        -v "$(dirname "$LICENSE_FILE")":/fs_license:ro \
        "$DOCKER_IMAGE" \
        --fs_license /fs_license/$(basename "$LICENSE_FILE") \
        --t1 /data/"$T1W_BASENAME" \
        --sid "$SUB" \
        --sd /output \
        --device cpu \
        --threads "$THREADS" \
        --3T \
        2>&1 | tee "$LOGFILE"

    RC=${PIPESTATUS[0]}
    if [[ $RC -eq 0 ]]; then
        echo "✓ $SUB completed successfully."
        (( PASS++ )) || true
    else
        echo "✗ $SUB FAILED (exit code $RC). See $LOGFILE" >&2
        (( FAIL++ )) || true
    fi
done

echo ""
echo "======================================================================"
echo " Done: $PASS passed, $FAIL failed, $SKIP skipped"
echo "======================================================================"
