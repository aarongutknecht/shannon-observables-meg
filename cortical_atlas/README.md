# Cortical Atlas Analysis — Shannon Observables MEG

Source reconstruction and information-theoretic analysis of resting-state MEG from the NIMH Healthy Research Volunteer Dataset ([OpenNeuro ds005752](https://openneuro.org/datasets/ds005752)).

This repository contains the Shannon Observables MEG release: cortical maps of **redundancy** r̄ and **vulnerability** v̄ across **62 DKTatlas cortical regions** in **N = 54 healthy adults**.

## Overview

The pipeline has two stages:

1. **Source reconstruction** (**Python / MNE-Python**)  
   Individual cortical anatomy via [FastSurfer](https://github.com/Deep-MI/FastSurfer), MEG↔MRI coregistration, single-shell BEM, and LCMV beamforming to extract 62-region DKTatlas parcel time series from resting-state CTF MEG.

2. **State-space Granger causality** (**MATLAB / [MVGC2](https://github.com/lcbarnett/MVGC2)**)
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

Final Granger causality results used for the paper figure are included under:

`results/gc_analysis/gc_dkt62_results/`

These are final derived results, not raw or source-reconstructed input data.
The source-reconstructed time series needed to rerun the GC analysis are
distributed separately via [Zenodo](https://doi.org/10.5281/zenodo.19250047),
and the raw MEG/MRI data are distributed via OpenNeuro.

## System requirements

**Tested operating system:** Ubuntu 24.04.2 LTS.

The commands are written for Linux/bash. No non-standard hardware is required.
The full raw-data source reconstruction is CPU and disk intensive but does not
require specialized MEG/MRI hardware; those inputs are ordinary data files.

The "Needed for" column uses the first entry point for which the dependency is
needed. The full raw-data pipeline uses the full Python environment defined in
`environment.yml`.

| Dependency | Tested / recommended version | Needed for |
|------------|------------------------------|------------|
| Conda | 26.1.0 | Paper figure notebook |
| Python | 3.11.14 | Paper figure notebook |
| pip | 26.0.1 | Paper figure notebook |
| NumPy | 2.3.5 | Paper figure notebook |
| SciPy | 1.17.0 | Paper figure notebook |
| Matplotlib | 3.10.8 | Paper figure notebook |
| Jupyter Notebook | 7.5.3 | Paper figure notebook |
| MNE-Python | 1.11.0 (`environment.yml` requires `>=1.6`) | Paper figure notebook |
| h5py | 3.15.1 | Paper figure notebook |
| nibabel | 5.3.2 | Paper figure notebook |
| nilearn | 0.13.1 | Paper figure notebook |
| pandas | 3.0.1 | Paper figure notebook |
| adjustText | 1.3.0 | Paper figure notebook |
| MATLAB | R2025b | GC from source-reconstructed data |
| [MVGC2](https://github.com/lcbarnett/MVGC2) | commit `b99eadad18253582f0a7e8e52d9b0c4c9c40a02f` | GC from source-reconstructed data |
| PyVista | 0.47.0 | Full raw-data pipeline only |
| PyYAML | 6.0.3 | Full raw-data pipeline only |
| joblib | 1.5.3 | Full raw-data pipeline only |
| tqdm | 4.67.3 | Full raw-data pipeline only |
| mne-bids | 0.18.0 (`environment.yml` requires `>=0.14`) | Full raw-data pipeline only |
| hdf5storage | 0.2.2 | Full raw-data pipeline only |
| Docker + `deepmi/fastsurfer:cpu-v2.4.2` | Docker 28.2.2, FastSurfer v2.4.2 | Full raw-data pipeline only |
| FreeSurfer | 7.4.1 (tarball) | Full raw-data pipeline only |

## Installation guide

Commands below assume the repository root as the starting directory.

The full installation consists of the Python environment, MATLAB/MVGC2,
Docker/FastSurfer, and FreeSurfer. The later "Instructions for use" section
gives three entry points; not every entry point requires the full stack.

### Python environment

Required for the paper-figure notebook and for full raw-data reproduction:
Install Conda, Miniconda, or another compatible Conda implementation first if
it is not already available.

```bash
cd cortical_atlas
conda env create -f environment.yml
conda activate nimh-meg-atlas
```

Typical install time on a normal desktop computer is ~5-15 minutes, depending
on network speed and Conda solver/cache state. The first figure run may also
download MNE `fsaverage` surfaces if they are not already cached locally.

### MATLAB and MVGC2

Required for rerunning the Granger causality analysis and for full raw-data
reproduction. Install MATLAB separately and ensure `matlab` is available on
the shell `PATH`. Then clone the tested [MVGC2](https://github.com/lcbarnett/MVGC2)
revision into `cortical_atlas/`:

```bash
cd cortical_atlas
git clone https://github.com/lcbarnett/MVGC2.git
git -C MVGC2 checkout b99eadad18253582f0a7e8e52d9b0c4c9c40a02f
```

Typical MVGC2 install time is <1 minute after MATLAB is available.

### Docker and FastSurfer

Required only for reproducing source reconstruction from raw MRI data. Install
Docker separately, then pull the FastSurfer image:

```bash
sudo docker pull deepmi/fastsurfer:cpu-v2.4.2
```

Typical image pull time is ~5-30 minutes, depending on network speed.

### FreeSurfer

Required only for reproducing source reconstruction from raw MRI data. Install
FreeSurfer 7.4.1 from the tarball distribution and obtain a FreeSurfer license
file from the [FreeSurfer license page](https://surfer.nmr.mgh.harvard.edu/fswiki/License).
The pipeline uses FreeSurfer's `mri_watershed` command for BEM surface
generation. Before running the raw-data workflow, set:

```bash
export FREESURFER_HOME=/path/to/freesurfer
export FREESURFER_LICENSE=/path/to/freesurfer/license.txt
export PATH="$FREESURFER_HOME/bin:$PATH"
```

Typical FreeSurfer setup time after downloading the tarball is ~5-20 minutes.

### Required components by entry point

| Entry point | Required installed components |
|-------------|-------------------------------|
| Generate the paper figure | Python environment |
| Rerun GC from source-reconstructed data | MATLAB, MVGC2 |
| Reproduce everything from raw data | Python environment, MATLAB, MVGC2, Docker/FastSurfer, FreeSurfer |

The full installation time excluding installation of external system software
such as Conda, MATLAB, and Docker itself is typically ~15-60 minutes on a
normal desktop computer, dominated by package downloads and the FastSurfer
Docker image pull.

## Demo

This release supports three demonstration modes. The exact commands are given
in the following "Instructions for use" section.

The lightweight demo is to regenerate the paper figure from the final GC
results included in `results/gc_analysis/gc_dkt62_results/`. Expected output is
`figures/fig_cortical_atlas_maps_spearman.pdf` and
`figures/fig_cortical_atlas_maps_spearman.png`, containing cortical surface
maps of r̄ and v̄, an r̄ vs v̄ scatter across the 62 regions, and leave-one-out
subject-average correlations. Expected runtime is ~1-5 minutes on a normal
desktop computer.

The two computational reproduction workflows can also be run as single-subject
demos: rerun the GC analysis for one subject from source-reconstructed data, or
reproduce the full pipeline for one subject from raw OpenNeuro data. The
GC-only single-subject demo takes ~10-30 minutes. The full raw-data
single-subject demo takes several hours.

## Instructions for use

Commands below assume the repository root as the starting directory.

### Generate the paper figure

**Requires:** Python environment from `environment.yml`. The notebook uses MNE
`fsaverage` surfaces for cortical rendering and may download them automatically
on first run if they are not already cached locally.

The notebook reads the included final GC results from
`results/gc_analysis/gc_dkt62_results/`. If that folder is absent, it falls back
to the most recent timestamped `results/gc_analysis/gc_dkt62_*` folder produced
by a local GC rerun.

Run:

```bash
cd cortical_atlas
conda activate nimh-meg-atlas
jupyter notebook notebooks/paper_figure_cortical_maps.ipynb
```

**Expected output:** `figures/fig_cortical_atlas_maps_spearman.pdf` and
`figures/fig_cortical_atlas_maps_spearman.png`, containing cortical surface
maps of r̄ and v̄, an r̄ vs v̄ scatter across the 62 regions, and leave-one-out
subject-average correlations.

**Expected runtime:** ~1–5 minutes.

### Rerun the GC analysis from source-reconstructed data

**Requires:** MATLAB and the [MVGC2 toolbox](https://github.com/lcbarnett/MVGC2) cloned into the project root.

Download the source-reconstructed time series from [Zenodo](https://doi.org/10.5281/zenodo.19250047). The archive contains a top-level `source_timeseries/` directory with per-subject `source_epochs.mat` files; extract that directory to:

`derivatives/source_timeseries/`

Then run:

```bash
cd cortical_atlas
make gc-analysis
```

For a single-subject demo, make sure the subject has
`source_epochs.mat` in the expected location, for example
`derivatives/source_timeseries/sub-ON08710/source_epochs.mat`. Then replace the
final command with:

```bash
make gc-analysis SUB=sub-ON08710
```

**Expected output:** a new timestamped folder
`results/gc_analysis/gc_dkt62_<timestamp>/` containing `gc_results.mat` and
`gc_summary.csv`.

**Expected runtime:** ~10–30 minutes per subject; several hours for all 54 subjects.

### Reproduce everything from raw data

**Requires:** Python environment from `environment.yml`, MATLAB, [MVGC2](https://github.com/lcbarnett/MVGC2), Docker with image `deepmi/fastsurfer:cpu-v2.4.2`, [FreeSurfer](https://surfer.nmr.mgh.harvard.edu/) (tarball, for `mri_watershed`), and raw OpenNeuro data. Ensure `FREESURFER_HOME` and `FREESURFER_LICENSE` are set.

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
conda activate nimh-meg-atlas
sudo docker pull deepmi/fastsurfer:cpu-v2.4.2
export FREESURFER_HOME=/path/to/freesurfer   # tarball install; only mri_watershed is used
export FREESURFER_LICENSE=/path/to/freesurfer/license.txt
export PATH="$FREESURFER_HOME/bin:$PATH"
make fastsurfer
make coreg
make source
make mat
make gc-analysis
```

For a single-subject demo, download only the required raw files for one subject
listed in `subjects.txt` and pass that subject ID to each subject-specific
target:

```bash
cd cortical_atlas
conda activate nimh-meg-atlas
sudo docker pull deepmi/fastsurfer:cpu-v2.4.2
export FREESURFER_HOME=/path/to/freesurfer
export FREESURFER_LICENSE=/path/to/freesurfer/license.txt
export PATH="$FREESURFER_HOME/bin:$PATH"
make fastsurfer SUB=sub-ON08710
make coreg SUB=sub-ON08710
make source SUB=sub-ON08710
make mat SUB=sub-ON08710
make gc-analysis SUB=sub-ON08710
```

**Expected output:** source-reconstructed parcel time series under
`derivatives/source_timeseries/`, MATLAB-ready inputs from `make mat`, and final
GC outputs under a timestamped `results/gc_analysis/gc_dkt62_<timestamp>/`
folder.

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
