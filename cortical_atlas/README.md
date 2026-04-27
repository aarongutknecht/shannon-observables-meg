# Cortical Atlas Analysis — Shannon Observables MEG

Source reconstruction and information-theoretic analysis of resting-state MEG from the NIMH Healthy Research Volunteer Dataset ([OpenNeuro ds005752](https://openneuro.org/datasets/ds005752)).

This repository contains the Shannon Observables MEG release: cortical maps of **redundancy** R̄ and **vulnerability** V̄ across **62 DKTatlas cortical regions** in **N = 54 healthy adults**.

## Overview

The pipeline has two stages:

1. **Source reconstruction** (**Python / MNE-Python**)  
   Individual cortical anatomy via [FastSurfer](https://github.com/Deep-MI/FastSurfer), MEG↔MRI coregistration, single-shell BEM, and LCMV beamforming to extract 62-region DKTatlas parcel time series from resting-state CTF MEG.

2. **State-space Granger causality** (**MATLAB / MVGC2**)  
   Innovations-form state-space models fitted to the 62-region time series, yielding joint, unconditional pairwise, and conditional pairwise GC, plus derived redundancy and vulnerability measures per region.

## Dataset

| Property | Value |
|----------|-------|
| Source | [OpenNeuro ds005752](https://openneuro.org/datasets/ds005752) |
| Name | NIMH Healthy Research Volunteer Dataset |
| Subjects analysed | 54 |
| MEG system | CTF 275-channel magnetometers, 1200 Hz native sampling |
| MRI scanner | 3T |
| Paradigm | ~6-minute eyes-open resting state |

## Requirements

**Tested operating system:** Ubuntu 24.04.2 LTS.

| Dependency | Tested version | Scope |
|------------|----------------|-------|
| Python | 3.11 | figure, GC, full pipeline |
| Jupyter Notebook | environment-defined | figure |
| MNE-Python | ≥ 1.6 | figure, GC, full pipeline |
| MATLAB | R2025b | GC, full pipeline |
| MVGC2 | commit `b99eadad18253582f0a7e8e52d9b0c4c9c40a02f` | GC, full pipeline |
| Docker + `deepmi/fastsurfer:latest` | Docker 28.2.2, FastSurfer v2.4.2 at time of analysis | full pipeline only |
| FreeSurfer | 7.4.1 (tarball) | full pipeline only |
| Conda | 26.1.0 | figure, GC, full pipeline |


## Reproduction options

Commands below assume the repository root as the starting directory.

### 1) Generate the paper figure

**Requires:** Python environment from `environment.yml`. The notebook uses MNE `fsaverage` surfaces for cortical rendering and may download them automatically on first run if they are not already cached locally.

GC results are already included under:

`results/gc_analysis/gc_dkt62_results/`

Run:

```bash
cd cortical_atlas
conda env create -f environment.yml
conda activate nimh-meg-atlas
jupyter notebook notebooks/paper_figure_cortical_maps.ipynb
```

**Expected output:** `fig_cortical_atlas_maps_spearman.pdf` and `.png`, containing cortical surface maps of R̄ and V̄, an R̄ vs V̄ scatter across the 62 regions, and leave-one-out subject-average correlations.

**Expected runtime:** ~1–5 minutes.

### 2) Rerun the GC analysis

**Requires:** Python environment from `environment.yml`, MATLAB, and the [MVGC2 toolbox](https://github.com/lcbarnett/MVGC2) cloned into the project root.

Download the source-reconstructed time series from [Zenodo](https://doi.org/10.5281/zenodo.19250047). The archive contains a top-level `source_timeseries/` directory with per-subject `source_epochs.mat` files; extract that directory to:

`derivatives/source_timeseries/`

Then run:

```bash
cd cortical_atlas
conda env create -f environment.yml
conda activate nimh-meg-atlas
git clone https://github.com/lcbarnett/MVGC2.git
git -C MVGC2 checkout b99eadad18253582f0a7e8e52d9b0c4c9c40a02f
make gc-analysis
```

**Expected runtime:** ~10–30 minutes per subject; several hours for all 54 subjects.

### 3) Reproduce everything from raw data

**Requires:** Python environment from `environment.yml`, MATLAB, [MVGC2](https://github.com/lcbarnett/MVGC2), Docker with image `deepmi/fastsurfer:latest`, [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) (tarball, for `mri_watershed`), and raw OpenNeuro data. Ensure `FREESURFER_HOME` and `FREESURFER_LICENSE` are set.

Download the dataset from [OpenNeuro ds005752](https://openneuro.org/datasets/ds005752) under `data/ds005752/`. The pipeline only needs these files per subject:

```
data/ds005752/sub-<ID>/ses-01/anat/sub-<ID>_ses-01_acq-MPRAGE_T1w.nii.gz
data/ds005752/sub-<ID>/ses-01/anat/sub-<ID>_ses-01_acq-forMEGCoregistration_T1w.nii.gz
data/ds005752/sub-<ID>/ses-01/anat/sub-<ID>_ses-01_acq-forMEGCoregistration_T1w.json
data/ds005752/sub-<ID>/ses-01/meg/sub-<ID>_ses-01_task-rest_run-01_meg.ds
```

Run:

```bash
cd cortical_atlas
conda env create -f environment.yml
conda activate nimh-meg-atlas
git clone https://github.com/lcbarnett/MVGC2.git
git -C MVGC2 checkout b99eadad18253582f0a7e8e52d9b0c4c9c40a02f
sudo docker pull deepmi/fastsurfer:latest
export FREESURFER_HOME=/path/to/freesurfer   # tarball install; only mri_watershed is used
export FREESURFER_LICENSE=/path/to/freesurfer/license.txt
export PATH="$FREESURFER_HOME/bin:$PATH"
make fastsurfer
make coreg
make source
make mat
make gc-analysis
```

**Expected runtime:** Several hours per subject. The full 54-subject pipeline takes multiple days if run serially.

## Repository layout

```text
cortical_atlas/
├── config.yaml
├── Makefile
├── environment.yml
├── subjects.txt
├── scripts/
│   ├── source_recon/
│   └── gc_analysis/
├── data/
│   └── ds005752/          # OpenNeuro dataset (user-provided)
├── derivatives/           
│   ├── fastsurfer/        # pipeline-generated 
│   ├── transforms/        # pipeline-generated 
│   ├── forward/           # pipeline-generated 
│   └── source_timeseries/ # pipeline-generated or user-provided (Zenodo)
├── results/
│   └── gc_analysis/gc_dkt62_results/
├── notebooks/
│   └── paper_figure_cortical_maps.ipynb
└── logs/
```

## Configuration

Source-reconstruction parameters live in `config.yaml`.

GC analysis parameters are set in the `PARAMETERS` section of
`scripts/gc_analysis/run_gc_analysis_cortical.m`.

## License

Pipeline code is released under the **MIT License**.

The NIMH dataset is distributed separately via OpenNeuro and remains subject to its own data use terms.
