# Shannon Observables MEG

Code release for the Shannon Observables MEG cortical atlas analysis.

The available arXiv preprint uses the earlier title:

> Gutknecht AJ, Rosas FE, Ehrlich DA, Makkeh A, Mediano PAM, Wibral M (2025).
> **Shannon invariants: A scalable approach to information decomposition.**
> arXiv:2504.15779

This public repository contains the cortical atlas MEG analysis used in the
manuscript. It estimates redundancy and vulnerability from state-space Granger
causality on resting-state MEG source time series from healthy adults.

## Contents

### `cortical_atlas/`

Surface-based cortical analysis across 62 DKT atlas regions in 54 healthy
adults using the OpenNeuro dataset `ds005752`.

See `cortical_atlas/README.md` for setup, data requirements, and usage.

## Data

Source-reconstructed time series (Level 2 input):
https://doi.org/10.5281/zenodo.19250047

## Requirements

The full Python environment for this release is defined in:
- `cortical_atlas/environment.yml`

MATLAB is required for the state-space GC stage together with the MVGC2 toolbox:
https://github.com/lcbarnett/MVGC2

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
