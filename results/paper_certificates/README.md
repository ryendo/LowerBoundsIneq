# Certificates used in the paper

This directory archives the completed interval certificates for the
upper area-perimeter inequality.

- `siudeja_upper_global_rayleigh_manifest.json` records the cover used in
  the paper: the one-dimensional Rayleigh certificate, the conforming
  affine atlas, and the residual calculation.
- `siudeja_upper_global_manifest.json` records the earlier successful
  combination based on the analytic thin-sector calculation.
- `omega_up_residual/omega_up_all_manifest.json` records the complete
  residual calculation on the upper region. The corresponding cell
  records are stored in `omega_up_all_cells.csv.gz`.
- `upper_affine_atlas/upper_manifest_interval_20260728_162102.json`
  records the complete conforming affine-atlas calculation. The
  corresponding cell records are stored in
  `upper_atlas_interval_20260728_162102.csv.gz`.

The thin-region Rayleigh certificate is stored separately in
`results/muhat1_rayleigh_t038_hpc_20260730/`.

The full production calculation for the middle region is still running
from its saved checkpoints. Its final manifest and cell records will be
added only after the validator reports complete coverage.
