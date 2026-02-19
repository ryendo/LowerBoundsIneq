# LowerBoundInequalities

Computer-assisted verification code (MATLAB) for sharp **Dirichlet eigenvalue inequalities on planar triangles**, based on verified finite element methods and (optionally) interval arithmetic.

This repository contains the implementation used for the workflow described in the paper:

- **R. Endo, X. Liu, P. Mariano**, *Sharp Dirichlet Eigenvalue Inequalities on Triangles*.

The main entry point is the `VerificationRunner` class, which orchestrates the verification over a decomposition of the triangle parameter space.

---

## What is being verified?

Let:

- $\triangle$ be a planar triangle,
- $\lambda_1(\triangle)$ be the **first Dirichlet eigenvalue** of the Laplacian on $\triangle$,
- $|\triangle|$ be the **area** of $\triangle$,
- $|\partial \triangle|$ be the **perimeter** of $\triangle$.

This project verifies nonnegativity of the following shape functionals (see `src/algorithms/compute_J_lower_bound.m`):

- **J1 (Laugesen–Siudeja type)**  
  $$
  J_1(\triangle) = \lambda_1(\triangle)\,|\triangle| - \frac{\pi^2}{16}\frac{|\partial \triangle|^2}{|\triangle|} - \frac{7\sqrt{3}\,\pi^2}{12}.
  $$

- **J2 (Cheeger type)**  
  $$
  J_2(\triangle) = \lambda_1(\triangle)\,|\triangle|
  - \frac{4\pi^2}{\left(3+\sqrt{\pi\sqrt{3}}\right)^2}
    \cdot
    \frac{\left(|\partial \triangle|+\sqrt{4\pi|\triangle|}\right)^2}{4|\triangle|}.
  $$

The verification goal is to certify **$J_1(\triangle) \ge 0$** (and/or **$J_2(\triangle) \ge 0$**) for all triangles $\triangle$ in the parameter domain considered in the paper.

---

## High-level approach

The code follows the paper’s “computer-assisted proof” structure and splits the parameter space into regions (terminology from the paper):

- **Ω_up**: triangles near the equilateral one  
  Verified via **Algorithm 2** (`src/algorithms/Algorithm2_VerifyOmegaUp.m`) by bounding second derivatives and verifying positive definiteness of the Hessian.

- **Ω_mid**: an intermediate region  
  Verified via **Algorithm 3** (`src/algorithms/Algorithm3_VerifyOmegaMid.m`) by splitting the region into many “cells” and proving $J \ge 0$ on each cell using:
  - rigorous **lower bounds** for $\lambda_1$ (Crouzeix–Raviart and, optionally, Lehmann–Goerisch),
  - rigorous **geometry bounds** (area/perimeter),
  - interval arithmetic to propagate bounds to $J$.

> Note: `VerificationRunner.runCompleteVerification` currently runs Ω_up and Ω_mid. The mention of “Ω_down” exists in comments but is not implemented as an explicit step in this code.

---

## Repository layout

```text
LowerBoundInequalities/
  VerificationRunner.m        # Main driver class
  my_intlab_config.m          # Path setup + INTLAB initialization (edit before use)

  src/
    algorithms/               # Algorithm 2/3 and shape functional evaluation
    fem/                      # Eigenvalue bounds (CR/LG), FEM utilities
    interval/                 # Small wrappers controlled by INTERVAL_MODE
    mesh/                     # gmsh integration + mesh readers
    lib/mesh/                 # helper mesh utilities

  inputs/
    cell_def.csv              # Ω_mid cell definitions (bounds + FEM parameters)

  results/
    J1_OmegaUp.csv            # Example output
    J1_OmegaMid.csv           # Example output

  tests/
    run_tests.m               # Test runner (may require path tweaks in some environments)

  VFEM2D/                     # included third-party FEM/eigenvalue library (see its README)
  VFEM2D_revised/             # revised scripts used by this project
  veigs/                      # verified generalized eigenvalue solver (see its README)
