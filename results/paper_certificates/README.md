# Certificates used in the paper

This directory archives the completed interval certificates for the
upper area-perimeter inequality.

- `siudeja_upper_global_rayleigh_manifest.json` preserves the cover and
  extrema from the earlier one-dimensional Rayleigh run, but is now marked
  incomplete.  That run rounded the rational exponents before INTLAB and
  did not retain replayable frozen trials.  It must not be cited as a proof
  certificate until a corrected v2 artifact has passed
  `validate_muhat1_rayleigh_certificate` and this global manifest has been
  rebuilt against its SHA-256 digest.
- `siudeja_upper_global_manifest.json` records the earlier successful
  combination based on the analytic thin-sector calculation.
- `omega_up_residual/omega_up_all_manifest.json` records the complete
  residual calculation on the upper region. The corresponding cell
  records are stored in `omega_up_all_cells.csv.gz`.
- `upper_affine_atlas/upper_manifest_interval_20260728_162102.json`
  records the complete conforming affine-atlas calculation. The
  corresponding cell records are stored in
  `upper_atlas_interval_20260728_162102.csv.gz`.
- `omega_mid_geometry/omega_mid_geometric_coverage_manifest.json` records
  the independent exact-decimal/INTLAB proof that the 188,623 middle-region
  cells cover the part of shape space with `y >= 0.04` outside the closed
  radius-0.122 ball about the equilateral point.  Its file SHA-256 is
  `7c9bd54f065d598469e987a6066e25a26b18d93a33c310cd4299c75b32b71395`.

The legacy thin-region output is stored separately in
`results/muhat1_rayleigh_t038_hpc_20260730/`; its README explains why it
is regression evidence rather than a current proof certificate.

The full production calculation for the middle region is still running
from its saved checkpoints. Its final manifest and cell records will be
added only after the validator reports complete coverage.
