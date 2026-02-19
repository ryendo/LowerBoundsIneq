# LowerBoundsIneq

MATLAB code for **computer-assisted verification of sharp Dirichlet Laplacian eigenvalue inequalities on planar triangles**, using verified finite element bounds and (optionally) **interval arithmetic (INTLAB)**.

This repository implements the workflow described in:

- **R. Endo, X. Liu, P. Mariano**, *Sharp Dirichlet Eigenvalue Inequalities on Triangles*.

The main entry point is **`VerificationRunner.m`**, which orchestrates the verification over a decomposition (“cells”) of the triangle parameter space.

---

## What is being verified?

Let

- \(\triangle\) be a planar triangle,
- \(\lambda_1(\triangle)\) be the **first Dirichlet eigenvalue** of the Laplacian on \(\triangle\),
- \(|\triangle|\) be the **area**,
- \(|\partial\triangle|\) be the **perimeter**.

This project certifies nonnegativity of the following shape functionals (see `src/algorithms/compute_J_lower_bound.m`):

### J1 (Laugesen–Siudeja type)
\[
J_1(\triangle)=\lambda_1(\triangle)|\triangle|
-\frac{\pi^2}{16}\frac{|\partial\triangle|^2}{|\triangle|}
-\frac{7\sqrt{3}\pi^2}{12}.
\]

### J2 (Cheeger type)
\[
J_2(\triangle)=\lambda_1(\triangle)|\triangle|
-\frac{4\pi^2}{\left(3+\sqrt{\pi\sqrt{3}}\right)^2}
\cdot
\frac{\left(|\partial \triangle|+\sqrt{4\pi|\triangle|}\right)^2}{4|\triangle|}.
\]

Goal: certify **\(J_1(\triangle)\ge 0\)** and/or **\(J_2(\triangle)\ge 0\)** for all triangles in the parameter domain considered in the paper.

---

## High-level approach

The parameter space is split into regions (paper terminology):

- **Ω_up** (near the equilateral triangle)  
  Verified via **Algorithm 2** (`src/algorithms/Algorithm2_VerifyOmegaUp.m`) using bounds on second derivatives and positivity checks.

- **Ω_mid** (intermediate region)  
  Verified via **Algorithm 3** (`src/algorithms/Algorithm3_VerifyOmegaMid.m`) by splitting Ω_mid into many **cells** and proving \(J\ge 0\) on each cell using:
  - rigorous **lower bounds** for \(\lambda_1\) (CR; optionally LG),
  - rigorous **geometry bounds** (area/perimeter),
  - optional **interval propagation** to keep everything certified.

> Note: `VerificationRunner.runCompleteVerification` currently runs **Ω_up** and **Ω_mid**. Any mention of “Ω_down” may appear in comments but is not a separate step here.

---

## Requirements

- **MATLAB** (recent versions recommended)
- **INTLAB** (optional but recommended for fully rigorous interval propagation)
- **Gmsh** (for mesh generation)  
  `my_intlab_config.m` notes the code was tested with **Gmsh 4.8.4**; other versions may change meshes slightly.

---

## Setup

1. **Clone / unzip** the repository so you have the folder `LowerBoundsIneq/`.

2. **Edit** `my_intlab_config.m` to match your environment:
   - `gmsh_command` : path to your `gmsh` executable
   - `mesh_path` : folder for temporary mesh files (typically `.../LowerBoundsIneq/src/mesh`)
   - `intlab_root` : path to your INTLAB installation (if used)
   - `INTERVAL_MODE` : set to `1` (interval) or `0` (floating point)

3. In MATLAB, from the project root:
   ```matlab
   cd LowerBoundsIneq
   my_intlab_config
````

---

## Quick start

### Run the full workflow (Ω_up + Ω_mid)

```matlab
my_intlab_config;

runner = VerificationRunner( ...
    'verbose', true, ...
    'resume_enabled', true, ...
    'save_intermediate', true);

results = runner.runCompleteVerification('J1', ...
    'cell_def_file', 'inputs/cell_def.csv', ...
    'cell_range', []);  % [] = all cells
```

### Run a subset of Ω_mid cells (useful for debugging / parallel runs)

`cell_range` refers to **row indices** in `cell_def.csv` (not the cell IDs themselves).

```matlab
results_mid = runner.verifyOmegaMid('J1', [1 200], 'inputs/cell_def.csv');
```

---

## Inputs

### `inputs/cell_def.csv`

Ω_mid is defined by a list of “cells”. The CSV contains (11 columns):

* `i` : cell ID
* `x_inf`, `x_sup` : bounds for the x-parameter
* `theta_inf`, `theta_sup` : bounds for the angle parameter
* `mesh_size_upper`, `fem_order_upper` : settings for the (upper-bound) FEM solve
* `mesh_size_lower_cr` : CR mesh size for eigenvalue lower bounds
* `isLG` : `1` to use Lehmann–Goerisch refinement, `0` otherwise
* `mesh_size_lower_LG`, `fem_order_lower_LG` : LG settings (if enabled)

---

## Outputs

All outputs go under `results_dir` (default: `results/`).

### Ω_mid per-cell log (CSV, append + resume)

* `results/J1_OmegaMid.csv` (or `J2_OmegaMid.csv`)
* Columns:

  * `conjecture, cell_id, verified, J_lower, status, note, run_timestamp`

If `resume_enabled = true`, existing rows are read and **already-computed cell IDs are skipped**.

### Ω_up summary (CSV)

* `results/J1_OmegaUp.csv` (or `J2_OmegaUp.csv`)
* Includes the minimum lower bounds found in Ω_up checks, plus the configuration used.

> You may see older example CSVs without headers in some folders. The current `VerificationRunner` writes headers for Ω_mid and Ω_up summaries.

---

## Tests

Tests live in `tests/`. From MATLAB, you can run them individually, e.g.:

```matlab
my_intlab_config
run tests/test_verification_runner.m
```

Or use MATLAB’s test runner:

```matlab
results = runtests('tests');
table(results)
```

---

## Repository layout

```text
LowerBoundsIneq/
  VerificationRunner.m
  my_intlab_config.m
  inputs/
    cell_def.csv
  src/
    algorithms/         # Algorithm 2/3 + shape functional bounds
    fem/                # FEM eigenvalue bounds + utilities
    interval/            # I_*.m wrappers controlled by INTERVAL_MODE
    mesh/                # Gmsh integration + mesh readers/scripts
    lib/
      VFEM2D/            # bundled external FEM/eigenvalue library (has its own LICENSE)
      VFEM2D_revised/    # project-specific revisions
      veigs/             # bundled eigen-solver utilities (has its own LICENSE)
      nouse_mesh/        # unused/legacy mesh utilities
  results/               # example outputs
  tests/                 # test scripts + sample CSVs in tests/results/
```

---

## Troubleshooting

* **“gmsh not found” / mesh generation fails**
  Set `gmsh_command` correctly in `my_intlab_config.m`.

* **INTLAB errors**
  Make sure `intlab_root` is correct and that `startintlab` works in your MATLAB session.

* **Different results across machines**
  Mesh generation can vary across Gmsh versions. For reproducibility, match the version used in your experiments as closely as possible.

---

## Licensing / third-party code

This repository bundles third-party components under `src/lib/` (e.g., `VFEM2D`, `veigs`) which include their own license files. Please check those licenses before redistributing or publishing derived code.

---

## Citation

If you use this code in academic work, please cite the associated paper:

* R. Endo, X. Liu, P. Mariano, *Sharp Dirichlet Eigenvalue Inequalities on Triangles*.
