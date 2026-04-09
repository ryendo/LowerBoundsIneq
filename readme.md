# Computer-Assisted Proof for Sharp Dirichlet Eigenvalue Inequalities on Triangles

This repository provides the MATLAB source code and computational framework for the computer-assisted part of the paper:

> **R. Endo, X. Liu, P. Mariano**  
> **Sharp Dirichlet Eigenvalue Inequalities on Triangles**

The project rigorously verifies sharp lower bounds for the first Dirichlet Laplacian eigenvalue on planar triangles by certifying positivity of two shape functionals:

- **J1**: the Laugesen–Siudeja-type functional,
- **J2**: a Cheeger-type functional.

The code implements the computer-assisted parts of Section 4 of the paper:

- **Step 1: $\Omega_{\mathrm{up}}$** — near the equilateral triangle, using certified lower bounds for second directional derivatives and convexity checks;
- **Step 2: $\Omega_{\mathrm{mid}}$** — the intermediate region, using per-cell certified lower bounds for the functionals.

The paper also treats the degenerate region **$\Omega_{\mathrm{down}}$** analytically. That part is **not** implemented here as a separate MATLAB verification routine.

---

## Background

Let $\triangle$ be a planar triangle, $\lambda_1(\triangle)$ its first Dirichlet eigenvalue, $|\triangle|$ its area, and $|\partial \triangle|$ its perimeter.
The paper proves that the equilateral triangle uniquely minimizes the following functionals:

### J1 (Laugesen–Siudeja type)

```math
J_1(\triangle)=\lambda_1(\triangle)|\triangle|-\frac{\pi^2}{16}\frac{|\partial \triangle|^2}{|\triangle|}-\frac{7\sqrt{3}\pi^2}{12}.
```

### J2 (Cheeger type)

```math
J_2(\triangle)=\lambda_1(\triangle)|\triangle|-\frac{4\pi^2}{\left(3+\sqrt{\pi\sqrt{3}}\right)^2}
\frac{\left(|\partial \triangle|+\sqrt{4\pi|\triangle|}\right)^2}{4|\triangle|}.
```

The computational strategy follows the decomposition of the triangle moduli space into three regions:

- **$\Omega_{\mathrm{up}}$**: a neighborhood of the equilateral triangle;
- **$\Omega_{\mathrm{mid}}$**: an intermediate region covered by finitely many cells in $(x,\theta)$ coordinates;
- **$\Omega_{\mathrm{down}}$**: a nearly degenerate region handled analytically in the paper.

Triangles are normalized as

```math
\triangle^{(x,y)} = \mbox{conv}\{(0,0),(1,0),(x,y)\},
```

with

```math
\Omega = \{(x,y) : x\ge 1/2,\ y>0,\ x^2+y^2\le 1\}.
```

---

## Core Libraries & Dependencies

This project relies on specialized libraries for verified numerical computation:

1. **INTLAB**: The fundamental toolbox for rigorous interval arithmetic in MATLAB.  
   **Source:** http://www.tuhh.de/ti3/intlab/  
   **Note:** `INTLAB_V12` and `INTLAB_V14` were used for the computation.

2. **Revised version of VFEM2D**: Used for rigorous finite element matrix assembly and high-precision eigenvalue bounds (Lehmann–Goerisch method).  
   **Source:** https://github.com/xfliu/VFEM2D  
   **Reference version:** 2025/12/13

3. **veigs**: Used for solving generalized matrix eigenvalue problems with rigorous error bounds together with eigenvalue-index information.  
   **Source:** https://github.com/yuuka-math/veigs  
   **Reference version:** 2025/12/13

In addition, the mesh-generation routines require:

- **gmsh** (tested with **gmsh 4.8.4**).

This archive contains the project-specific routines under `src/lib/VFEM2D_revised/`, but the current MATLAB path configuration in `my_intlab_config.m` also expects external installations of **INTLAB**, **VFEM2D**, and **veigs**.
Before running the code, make sure these dependencies are installed locally and that `my_intlab_config.m` points to the correct locations.

---

## Project Structure

```text
.
├── VerificationRunner.m
├── my_intlab_config.m
├── verify_all.m
├── README.md
├── inputs/
│   └── cell_def.csv
├── results/
│   ├── J1_OmegaMid.csv
│   ├── J1_OmegaUp.csv
│   ├── J1_OmegaUp_step1_2_cells.csv
│   ├── J1_OmegaUp_step1_3_axis.csv
│   ├── J2_OmegaMid.csv
│   ├── J2_OmegaUp.csv
│   ├── J2_OmegaUp_step1_2_cells.csv
│   └── J2_OmegaUp_step1_3_axis.csv
├── src/
│   ├── algorithms/
│   │   ├── Algorithm2_VerifyOmegaUp.m
│   │   ├── Algorithm3_VerifyOmegaMid.m
│   │   ├── compute_J_lower_bound.m
│   │   ├── compute_geometry_bounds.m
│   │   └── verify_J_positive.m
│   ├── fem/
│   │   ├── calc_ddlami_lower_bound.m
│   │   ├── calc_eigen_bounds_any_order_1k_wh.m
│   │   ├── calc_grad_error_bounds.m
│   │   ├── cell_lower_eig_bound.m
│   │   ├── cell_upper_eig_bound.m
│   │   └── auto_cluster_eigenvalues.m
│   ├── interval/
│   │   └── I_*.m
│   ├── mesh/
│   │   ├── make_mesh_by_gmsh.m
│   │   ├── make_mesh_by_gmsh_ref.m
│   │   ├── gmshread.m
│   │   └── create_mesh.sh
│   └── lib/
│       ├── VFEM2D_revised/
│       └── nouse_mesh/
├── tests/
└── tools/
```

### Main entry points

- **`VerificationRunner.m`**  
  Main orchestrator for the verified computations.

- **`src/algorithms/Algorithm2_VerifyOmegaUp.m`**  
  Implements the $\Omega_{\mathrm{up}}$ verification:
  - certification of positivity of $\partial_x^2 J_k$,
  - certification of positivity of $\partial_y^2 J_k$ on the symmetry axis.

- **`src/algorithms/Algorithm3_VerifyOmegaMid.m`**  
  Implements the $\Omega_{\mathrm{mid}}$ batch verification over the cells listed in `inputs/cell_def.csv`.

- **`src/algorithms/verify_J_positive.m`**  
  Single-cell verification routine used in Step 2.

- **`src/fem/calc_ddlami_lower_bound.m`**  
  Core certified routine for lower bounds of second directional shape derivatives of $\lambda_1$.

---

## Installation & Setup

Before running the code, make sure you have:

- MATLAB,
- INTLAB,
- gmsh,
- access to the required VFEM2D / veigs routines.

### 1. Clone or unzip the repository

```bash
git clone https://github.com/ryendo/LowerBoundsIneq.git
cd LowerBoundsIneq
```

or unzip the archive and open the project root in MATLAB.

### 2. Configure `my_intlab_config.m`

Edit the following items:

- `gmsh_command` — full path to your `gmsh` executable;
- `mesh_path` — directory for temporary mesh files;
- `intlab_root` — path to your local INTLAB installation;
- additional `addpath(...)` entries if your local copies of **VFEM2D** and **veigs** are stored elsewhere.

The code was tested with:

- `INTERVAL_MODE = 1`,
- `gmsh 4.8.4`.

### 3. Initialize the MATLAB environment

From the project root:

```matlab
my_intlab_config();
```

This function:

- adds the project directories to the MATLAB path,
- adds the configured dependency paths,
- initializes INTLAB.

---

## Usage: The `VerificationRunner` Class

The `VerificationRunner` class is the main interface for this project.
It manages parameters, runs the verification routines, and writes CSV summaries.

### 1. Initialize

```matlab
my_intlab_config();
runner = VerificationRunner();
```

You can also override parameters at construction time:

```matlab
runner = VerificationRunner( ...
    'verbose', true, ...
    'resume_enabled', true, ...
    'save_intermediate', true, ...
    'N_spectral', 1, ...
    'N_LG', 16, ...
    'N_rho', 64, ...
    'ord_LG', 2, ...
    'Nx', 10, ...
    'Ny', 100, ...
    'Ny_axis', 200);
```

### 2. Step 1: Verify `\Omega_{\mathrm{up}}`

Run the near-equilateral verification for either conjecture:

```matlab
results_up_J1 = runner.verifyOmegaUp('J1');
results_up_J2 = runner.verifyOmegaUp('J2');
```

This routine checks:

- positivity of $\partial_x^2 J_k$ on a rectangular cover of $\Omega_{\mathrm{up}}$,
- positivity of $\partial_y^2 J_k$ on the symmetry axis $x=1/2$.

Intermediate CSV files are written to `results/`:

- `J*_OmegaUp_step1_2_cells.csv`
- `J*_OmegaUp_step1_3_axis.csv`
- `J*_OmegaUp.csv`

### 3. Step 2: Verify $\Omega_{\mathrm{mid}}$

Run the batch verification over the cells in `inputs/cell_def.csv`:

```matlab
results_mid_J1 = runner.verifyOmegaMid('J1', [], 'inputs/cell_def.csv');
results_mid_J2 = runner.verifyOmegaMid('J2', [], 'inputs/cell_def.csv');
```

To verify only a subset of rows in `cell_def.csv`:

```matlab
results_mid_subset = runner.verifyOmegaMid('J1', [1 200], 'inputs/cell_def.csv');
```

Notes:

- `cell_range = [a b]` refers to **row indices** in `cell_def.csv`;
- `VerificationRunner` supports CSV-based resume behavior for $\Omega_{\mathrm{mid}}$;
- results are appended to `results/J1_OmegaMid.csv` or `results/J2_OmegaMid.csv`.

### 4. Run the complete implemented workflow

```matlab
results = runner.runCompleteVerification('J1', ...
    'cell_def_file', 'inputs/cell_def.csv', ...
    'cell_range', []);
```

This runs:

1. $\Omega_{\mathrm{up}}$,
2. $\Omega_{\mathrm{mid}}$.

The current implementation does **not** include a separate `verifyOmegaDown(...)` method.

### 5. Single-cell verification (debugging / experimentation)

If you want to inspect one cell manually, use `verify_J_positive`:

```matlab
cell_data = struct();
cell_data.x_inf = '0.799219';
cell_data.x_sup = '0.8';
cell_data.theta_inf = '0.048708';
cell_data.theta_sup = '0.048732';
cell_data.mesh_size_upper = 0.002633;
cell_data.fem_order_upper = 0.1249;
cell_data.mesh_size_lower_cr = 0.005;
cell_data.isLG = 0;
cell_data.mesh_size_lower_LG = 0.1249;
cell_data.fem_order_lower_LG = 2;

[verified, J_lower, diagnostics] = verify_J_positive('J1', cell_data);
```

### 6. Certified second derivative of $\lambda_1$

For lower bounds of second directional derivatives of the first eigenvalue:

```matlab
base_triangle = [0, 0, 1, 0, 0.5, sqrt(3)/2];
triangle      = [0, 0, 1, 0, 0.5, sqrt(3)/2];
e_direction   = [1, 0];

[lami, dlami, ddlami_lb] = calc_ddlami_lower_bound( ...
    1, base_triangle, triangle, e_direction, 1, 16, 64, 2);
```

This is the core routine used by the $\Omega_{\mathrm{up}}$ verification.

---

## Input Data: `inputs/cell_def.csv`

The file `inputs/cell_def.csv` defines the covering of $\Omega_{\mathrm{mid}}$ and the FEM parameters used per cell.
Its columns are:

- `i`
- `x_inf`, `x_sup`
- `theta_inf`, `theta_sup`
- `mesh_size_upper`, `fem_order_upper`
- `mesh_size_lower_cr`
- `isLG`
- `mesh_size_lower_LG`, `fem_order_lower_LG`

Geometrically, a cell is given in $(x,\theta)$ coordinates, with

```math
y = x\tan(\theta).
```

The Step 2 proof uses domain monotonicity on each cell to combine:

- lower bounds for $\lambda_1$,
- lower / upper geometric evaluations of area and perimeter,
- a certified lower evaluation of `J_1` or `J_2`.

---

## Output Files

The main output files are:

- `results/J1_OmegaUp.csv`, `results/J2_OmegaUp.csv`  
  summary of the $\Omega_{\mathrm{up}}$ verification;

- `results/J1_OmegaMid.csv`, `results/J2_OmegaMid.csv`  
  per-cell verification logs for $\Omega_{\mathrm{mid}}$;

- `results/J1_OmegaUp_step1_2_cells.csv`, `results/J2_OmegaUp_step1_2_cells.csv`  
  detailed rectangle-wise lower bounds for $\partial_x^2 J_k$;

- `results/J1_OmegaUp_step1_3_axis.csv`, `results/J2_OmegaUp_step1_3_axis.csv`  
  detailed interval-wise lower bounds for $\partial_y^2 J_k$ on the axis.

The complete workflow also writes a flattened CSV summary under:

```text
results/complete/
```

---

## Tests and Utilities

The repository includes:

- ad hoc MATLAB test scripts in `tests/`,
- visualization utilities for the `\Omega_{\mathrm{mid}}` cell decomposition,
- notebooks / helper scripts in `tools/` for building and refining `cell_def.csv`.

Useful starting points include:

```matlab
run tests/test_main1_OmegaUp.m
run tests/test_main2_OmegaMid.m
run tests/test_calc_ddlami_lower_bound.m
```

---

## Citation

If you use this code in academic work, please cite the associated paper:

> **R. Endo, X. Liu, P. Mariano**, *Sharp Dirichlet Eigenvalue Inequalities on Triangles*.

If you use the verified FEM components or supporting rigorous eigensolver framework, please also cite the relevant dependency projects and papers.
