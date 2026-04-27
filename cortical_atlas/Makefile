# ============================================================================
# MEG Atlas Pipeline – Makefile
# ============================================================================
#
# Full pipeline:
#   make fastsurfer                  # 1. FastSurfer cortical reconstruction
#   make coreg      SUB=sub-ON08710   # 2. Automated MEG↔MRI coregistration
#   make source     SUB=sub-ON08710   # 3. Source reconstruction (BEM→LCMV→labels)
#   make mat                          # 4. Convert .npy → .mat for MATLAB
#
# Shortcuts:
#   make fastsurfer SUB=sub-ON08710   # run one subject
#   make source                       # run all subjects
#   make source ARGS="--force"        # recompute even if outputs exist
#   make dry-run                      # preview pipeline for all subjects
#
# Requires:
#   - conda env: nimh-meg-atlas
#   - FREESURFER_HOME set (for mri_watershed)
#   - Docker (for FastSurfer)
# ============================================================================

SHELL   := /bin/bash
PYTHON  := python3
SRC     := scripts/source_recon
GC      := scripts/gc_analysis
CONFIG  := config.yaml
ARGS    ?=
SUB     ?=

# FreeSurfer home (needed by mri_watershed in the source step).
# Set before running: export FREESURFER_HOME=/path/to/your/freesurfer
export FREESURFER_HOME ?= /path/to/your/freesurfer
export PATH := $(FREESURFER_HOME)/bin:$(PATH)

# Build --subjects flag if SUB is set
ifdef SUB
  SUB_FLAG := --subjects $(SUB)
else
  SUB_FLAG :=
endif

.PHONY: all fastsurfer coreg source mat \
        gc-analysis dry-run clean-logs help

help:
	@echo "Targets:"
	@echo "  fastsurfer [SUB=sub-ON...] – FastSurfer cortical reconstruction"
	@echo "  coreg      SUB=sub-ON... – Automated MEG↔MRI coregistration"
	@echo "  source     [SUB=sub-ON...] – Source reconstruction pipeline"
	@echo "  mat                      – Convert .npy → .mat for MATLAB"
	@echo "  gc-analysis              – Run GC analysis (MATLAB)"
	@echo "  dry-run                  – Preview pipeline (all subjects)"
	@echo "  clean-logs               – Remove log files"

all: source mat

# ── Step 1: FastSurfer (Docker, CPU) ───────────────────────────────────────
fastsurfer:
	bash $(SRC)/run_fastsurfer.sh $(SUB)

# ── Step 2: Coregistration ────────────────────────────────────────────────
coreg:
	$(PYTHON) $(SRC)/run_coreg_auto.py -c $(CONFIG) $(SUB_FLAG) $(ARGS)

# ── Step 3: Source reconstruction (BEM → fwd → LCMV → labels) ─────────────
source:
	$(PYTHON) $(SRC)/run_pipeline_individual.py -c $(CONFIG) $(SUB_FLAG) $(ARGS)

# ── Step 4: Convert to MATLAB ─────────────────────────────────────────────
mat:
	$(PYTHON) $(SRC)/npy_to_mat.py $(ARGS)

# ── Step 5: GC analysis (MATLAB) ─────────────────────────────────────────
# Usage: make gc-analysis [SUB=sub-ON08710]
gc-analysis:
ifdef SUB
	matlab -batch "gc_subjects = {'$(SUB)'}; run('$(GC)/run_gc_analysis_cortical.m')"
else
	matlab -batch "run('$(GC)/run_gc_analysis_cortical.m')"
endif

# ── Dry-run ───────────────────────────────────────────────────────────────────
dry-run:
	$(PYTHON) $(SRC)/run_pipeline_individual.py -c $(CONFIG) --dry-run

# ── Housekeeping ──────────────────────────────────────────────────────────
clean-logs:
	rm -f logs/*.log
