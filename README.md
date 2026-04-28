# Shannon Observables MEG

Code release for the Shannon Observables MEG cortical atlas analysis.

The available arXiv preprint uses the earlier title:

> Gutknecht AJ, Rosas FE, Ehrlich DA, Makkeh A, Mediano PAM, Wibral M (2025).
> **Shannon invariants: A scalable approach to information decomposition.**
> arXiv:2504.15779

This public repository contains the cortical atlas MEG analysis used in the
manuscript. It estimates redundancy and vulnerability from state-space Granger
causality on resting-state MEG source time series from healthy adults in the
[NIMH Healthy Research Volunteer Dataset](https://openneuro.org/datasets/ds005752)
(OpenNeuro `ds005752`).

## Contents

### `cortical_atlas/`

Surface-based cortical analysis across 62 DKT atlas regions in 54 healthy
adults using the OpenNeuro dataset `ds005752`.

See `cortical_atlas/README.md` for full system requirements, installation,
final results, demo options, expected outputs, usage, and reproduction
instructions.

## Data

Raw MEG and MRI data are available from OpenNeuro:

- [OpenNeuro ds005752: NIMH Healthy Research Volunteer Dataset](https://openneuro.org/datasets/ds005752)

Source-reconstructed time series used as input to the Granger causality stage
are available from Zenodo:

- [Zenodo source-reconstructed time series](https://doi.org/10.5281/zenodo.19250047)

The repository includes the derived GC results needed to regenerate the paper
figure under `cortical_atlas/results/gc_analysis/gc_dkt62_results/`.

## Citation

Please cite the arXiv preprint by its preprint title:

```bibtex
@misc{gutknecht2025shannoninvariantsscalableapproach,
  title={Shannon invariants: A scalable approach to information decomposition},
  author={Aaron J. Gutknecht and Fernando E. Rosas and David A. Ehrlich
          and Abdullah Makkeh and Pedro A. M. Mediano and Michael Wibral},
  year={2025},
  eprint={2504.15779},
  archivePrefix={arXiv},
  primaryClass={cs.IT},
  url={https://arxiv.org/abs/2504.15779}
}
```

## License

MIT — see `LICENSE`.
