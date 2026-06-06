# LowerBoundsIneq

Computer-assisted verification of **sharp Dirichlet Laplacian eigenvalue inequalities on planar triangles**, using rigorous finite-element bounds with INTLAB interval arithmetic.

Companion code to

> R. Endo, X. Liu, P. Mariano, *Sharp Dirichlet Eigenvalue Inequalities on Triangles*.

---

## What is verified

Let $\triangle$ be a planar triangle, $\lambda_1(\triangle)$ the first Dirichlet eigenvalue of the Laplacian, $|\triangle|$ the area and $|\partial\triangle|$ the perimeter.

This repository certifies non-negativity of the two shape functionals:

### J1 (Laugesen–Siudeja)
$$
J_1(\triangle)\ =\ \lambda_1(\triangle)|\triangle|\ -\ \frac{\pi^{2}}{16}\frac{|\partial\triangle|^{2}}{|\triangle|}\ -\ \frac{7\sqrt{3}\,\pi^{2}}{12}\ \ge\ 0.
$$

### J2 (Cheeger type)
$$
J_2(\triangle)\ =\ \lambda_1(\triangle)|\triangle|\ -\ \frac{4\pi^{2}}{(3+\sqrt{\pi\sqrt{3}})^{2}}\cdot\frac{(|\partial\triangle|+\sqrt{4\pi|\triangle|})^{2}}{4|\triangle|}\ \ge\ 0.
$$

Both inequalities are sharp with equality at the equilateral triangle.

The parameter space of triangles (up to similarity, restricted by symmetry to $x\in[0.5,1]$, and with the equilateral neighborhood handled separately as $\Omega_{\rm up}$) is tiled by a collection of **cells** $\{C_i\}_{i=1}^{N}$ covering the intermediate region $\Omega_{\rm mid}$.  For each cell the verification uses:

- a **rigorous FEM lower bound** for $\lambda_1$ (Crouzeix–Raviart, optionally sharpened by Lehmann–Goerisch);
- **rigorous interval bounds** on area and perimeter;
- **INTLAB** interval arithmetic to propagate all floating-point rounding errors.

If every cell yields $J_k\ge 0$ then $J_k\ge 0$ on all of $\Omega_{\rm mid}$.

---

## Repository layout

```
LowerBoundsIneq/
├── VerifyTriangleInequalities.m    ← MAIN CLASS (public entry-point)
├── my_intlab_config.m              ← edit paths here before running
├── readme.md
├── inputs/
│   └── cell_def.csv                ← the cell tiling used for verification
├── results/                        ← final verification outputs
│   ├── J1_OmegaMid.csv             ← one row per cell with J1_lower
│   ├── J2_OmegaMid.csv             ← one row per cell with J2_lower
│   ├── verification_summary.md     ← human-readable summary
│   └── verification_summary.json
├── src/
│   ├── VerificationRunner.m        ← orchestrator (called by the main class)
│   ├── algorithms/                 ← Algorithm2/3, J bound, verify_J_positive
│   ├── fem/                        ← eigenvalue lower-bound utilities
│   ├── mesh/                       ← gmsh integration
│   ├── interval/                   ← INTLAB wrappers
│   └── lib/                        ← bundled external FEM/eigensolver libs
├── Intlab_V12/                     ← bundled INTLAB 12 (run once to regenerate its .mat)
├── scripts_run/                    ← parallel driver + helpers
├── tests/                          ← unit-level tests
└── tools/                          ← miscellaneous developer tools
```

---

## Core Libraries & Dependencies

This project relies on specialized libraries for verified numerical computation:

1. **INTLAB**: The fundamental toolbox for rigorous interval arithmetic in MATLAB.
   - **Source:** [http://www.tuhh.de/ti3/intlab/](http://www.tuhh.de/ti3/intlab/) [INTLAB_V12, INTLAB_V14 were used for the computation.]
2. Revised version of **VFEM2D**: Used for rigorous finite element matrix assembly and high-precision eigenvalue bounds (Lehmann–Goerisch method).
   - **Source:** [https://github.com/xfliu/VFEM2D](https://github.com/xfliu/VFEM2D) [2025/12/13]
3. **veigs**: Used for solving generalized matrix eigenvalue problems with rigorous error bounds with the information of indices.
   - **Source:** [https://github.com/yuuka-math/veigs](https://github.com/yuuka-math/veigs) [2025/12/13]

---

## Quick start

1. Clone the repository and `cd` into it.
2. Edit the constants at the top of `my_intlab_config.m` (at minimum set `GMSH_COMMAND`).
3. Start MATLAB in the project root and run **one** of:

### (A) Validate the committed results (no recomputation)

```matlab
v = VerifyTriangleInequalities();
v.run();
```

This reads `inputs/cell_def.csv` and `results/J{1,2}_OmegaMid.csv`, checks every row has `verified==1` and `J_lower>0`, and prints a summary.  Expected output:

```
[load] inputs/cell_def.csv
[load] results/J1_OmegaMid.csv
[load] results/J2_OmegaMid.csv
  cell_def: 95226 cells; J1: 95226 rows; J2: 95026 rows
...
========== VERIFICATION SUMMARY ==========
J1: 95226 cells  verified=95226 (100.000%)  J_lower min/med/max = 6.856e-07 / 8.003e+00 / 2.626e+03
J2: 95026 cells  verified=95026 (100.000%)  J_lower min/med/max = 1.232e-05 / 3.067e+01 / 2.633e+03

*** VERIFIED: J1 >= 0 and J2 >= 0 on all cells of Omega_mid. ***
```

### (B) Reproduce the full computation from scratch

```matlab
v = VerifyTriangleInequalities();
v.compute(20);      % 20 parallel workers; raise if you have more cores
```

`.compute(nworkers)` orchestrates the parallel parfor over all cells, auto-aggregates per-worker outputs into `results/J{1,2}_OmegaMid.csv`, and then calls `.run()` to validate.  On a 96-core machine with 20 workers the full run takes roughly 6–8 hours per conjecture.

---

## Input / Output formats

### `inputs/cell_def.csv`

One row per cell. The tiling used here has 95 226 cells. Columns:

| column              | meaning                                                              |
|---------------------|----------------------------------------------------------------------|
| `i`                 | cell identifier (= row index)                                        |
| `x_inf`, `x_sup`    | interval $[x_{\rm inf},x_{\rm sup}]$ for the apex $x$-parameter      |
| `theta_inf`, `theta_sup` | interval for the apex angle parameter $\theta$                  |
| `mesh_size_lower_cr`| mesh size for the CR lower-bound FEM                                  |
| `isLG`              | `1` to refine with Lehmann–Goerisch, `0` for CR-only                  |
| `mesh_size_lower_LG`| (if `isLG=1`) mesh size for the LG solve                              |
| `fem_order_lower_LG`| (if `isLG=1`) FEM order for the LG solve                              |

Example first rows:
```
i,x_inf,x_sup,theta_inf,theta_sup,mesh_size_lower_cr,isLG,mesh_size_lower_LG,fem_order_lower_LG
1,0.7992187500000001,0.8,0.0487077921309046,0.0487324060873969,0.005,0,0.1249,2
2,0.7992187500000001,0.8,0.0487324060873969,0.0487570200438891,0.005,0,0.1249,2
```

### `results/J1_OmegaMid.csv` and `results/J2_OmegaMid.csv`

One row per cell. Columns:

| column          | meaning                                                            |
|-----------------|--------------------------------------------------------------------|
| `conjecture`    | `J1` or `J2`                                                        |
| `cell_id`       | matches `i` in the cell_def                                         |
| `verified`      | `1` iff the rigorous $J_k$-lower bound is positive on this cell     |
| `J_lower`       | rigorous lower bound for $J_k$ over this cell                       |
| `status`        | `ok` or `error`                                                     |
| `note`          | free-text note (usually empty)                                      |
| `run_timestamp` | when this row was computed                                          |

Example:
```
conjecture,cell_id,verified,J_lower,status,note,run_timestamp
J1,1,1,7.09959390563225767e+00,"ok","",2026-04-19 21:29:22
J1,2,1,7.09708917028804009e+00,"ok","",2026-04-19 21:29:29
```

### `results/verification_summary.md`

Human-readable digest of the verification statistics; see the file itself for the current numbers.

---

## Computation pipeline, if you need to reproduce

If you want to regenerate the per-worker raw outputs yourself:

1. `my_intlab_config;` — initialise INTLAB and add paths.
2. `run_parallel_omegamid('J1', [], 20, 'inputs/cell_def.csv', 'results_raw');`
3. `run_parallel_omegamid('J2', [], 20, 'inputs/cell_def.csv', 'results_raw');`
4. `v = VerifyTriangleInequalities(); v.aggregate(); v.run();`

Internally `VerifyTriangleInequalities.compute()` does exactly this sequence.

The rigorous eigenvalue lower-bound routine is in `src/fem/cell_lower_eig_bound.m`; the J-bound is in `src/algorithms/compute_J_lower_bound.m`; the single-cell driver is `src/algorithms/verify_J_positive.m`.

---

## Citation

If you use this code please cite the associated paper:

- R. Endo, X. Liu, P. Mariano, *Sharp Dirichlet Eigenvalue Inequalities on Triangles*.

---

## Licenses

External libraries bundled under `src/lib/` (`VFEM2D`, `veigs`) and `Intlab_V12/` retain their own licenses (see their directories).
