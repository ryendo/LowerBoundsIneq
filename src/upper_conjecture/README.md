# Siudeja upper-conjecture verifier

This directory implements the ver10 certificate for

\[
\lambda_1(T)\leq \frac{\pi^2L^2}{12A^2}
                 +\frac{\sqrt3\,\pi^2}{3A}.
\]

For
\(T(x,y)=\operatorname{conv}\{(0,0),(1,0),(x,y)\}\), \(A=y/2\), this is

\[
\lambda_1(T(x,y))
\leq \frac{\pi^2L(x,y)^2}{3y^2}
     +\frac{2\sqrt3\,\pi^2}{3y}.
\]

The proof has three pieces.  No completed production extrema are quoted in
this file until their interval manifests have been generated and committed.

## 1. Thin sector: \(0<y\leq0.06\)

`upper_verify_thin_sector.m` uses the Freitas--Siudeja sector inclusion.
Writing

\[
\nu=\frac{\pi}{\arcsin y},\qquad m=\frac45\nu^{2/3},
\]

the sector trial

\[
u(r,\theta)=r^m(1-r^m)\sin(\nu\theta)
\]

gives

\[
\lambda_1(T)\leq
\frac{2Q(m,\nu)}{1+\sqrt{1-y^2}},\qquad
Q(m,\nu)=
\frac{(m+1)(2m+1)(3m+2)(\nu^2+2m^2)}{6m^3}.
\]

After multiplication by \(y^2\), the substitution
\(z=\nu^{-1/3}\) removes the apparent singularity at \(y=0\).  INTLAB then
checks the entire interval with an alternating Taylor enclosure for
\((1-\cos s)/s^2\), \(s=\arcsin y\).  The verifier does not numerically
evaluate a Bessel zero.

## 2. Compact P4 atlas: \(0.06\leq y\leq0.851\)

`run_upper_conjecture_atlas.m` tiles the diameter-normalized domain by a
rectangular superset cover of its curved boundary.  The production backend
uses one fixed P4/Gmsh topology with `trial_order=4`, `N_trial=8`, and
default widths `dx=dy=0.001`.

For each cell:

1. a conforming trial vector is computed at the midpoint;
2. its coefficients are frozen and transported affinely over the cell;
3. its four reference integrals are rigorously enclosed;
4. the Rayleigh quotient and target difference are evaluated with outward
   rounding.

If \(E_{11},E_{12},E_{22},M\) are the fixed trial integrals, the transported
Rayleigh quotient is

\[
R(x,y)=\frac{E_{11}}{M}
       +\frac{x^2E_{11}-2xE_{12}+E_{22}}{y^2M}.
\]

The min--max principle gives \(\lambda_1(T(x,y))\leq R(x,y)\).
A centered mean-value enclosure of `target - R` controls variation across
the full cell and avoids the wrapping of a repeated natural interval
expression.  The floating-point midpoint Ritz value is diagnostic only;
it is never used as the proof bound.

The P1 backend remains available for fast diagnostics and unit tests.  It
is not the production atlas backend.

## 3. Residual overlap: \(0.85\leq y\leq\sqrt3/2\)

The atlas ends at `y_up=0.851`.  The independent ver10 residual-Hessian
driver

`scripts_run/run_omega_up_all_residual_parallel.m`

certifies the upper functional on the band `y>=0.85`, through the
equilateral point, without explicit eigenfunctions or a spectral tail.
Thus

\[
[0.85,0.851]
\]

is a nonempty proof overlap.  The residual production grid is
`Nx=976`, `Ny=488`, `Ny_axis=1952`, with
`functional_scope='split'`.  Its JSON manifest records the precise band,
coverage, extrema, and Git commit.

An atlas manifest by itself proves only the thin-plus-compact scope through
its recorded `y_up`.  The global upper conjecture additionally requires a
complete interval residual manifest over the overlap and the remaining
equilateral band.

## Runtime requirements

INTLAB is external and is not bundled with the repository.  Set these
variables before starting MATLAB:

```sh
export INTLAB_ROOT=/path/to/Intlab_V12
export LOWERBOUNDS_MESH_PATH=/path/to/private/mesh-workspace
export GMSH_COMMAND=/usr/bin/gmsh
```

For a 20-worker HPC run, also prevent nested BLAS oversubscription:

```sh
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
```

Every worker runs `upper_prepare_runtime` and initializes the same external
`INTLAB_ROOT`.

## Smoke tests

A fast double-precision diagnostic is:

```matlab
addpath scripts_run
run_upper_conjecture_atlas( ...
    'mode','double', ...
    'backend','p1', ...
    'p1_n',12, ...
    'max_cells',8);
```

Double mode is always labelled `exploratory_double` and is never accepted
as a certificate.

To run the interval smoke tests:

```matlab
setenv('RUN_INTLAB_SMOKE','1')
test_upper_conjecture_smoke
test_upper_compact_cover
```

These tests require `INTLAB_ROOT` and Gmsh.  They verify the thin formula, a
small interval cell, the P4 path, and rigorous containment of the circular
boundary.

## Production commands

Run the complete thin-plus-compact atlas on a 20-worker node:

```matlab
addpath scripts_run
run_upper_conjecture_atlas( ...
    'mode','interval', ...
    'workers',20, ...
    'backend','p4', ...
    'trial_order',4, ...
    'N_trial',8, ...
    'y_min',0.06, ...
    'y_up',0.851, ...
    'residual_y_start',0.85, ...
    'dx',0.001, ...
    'dy',0.001);
```

The equivalent checked-in entry point is:

```matlab
run('scripts_run/run_upper_conjecture_atlas_full.m')
```

Run the overlapping residual certificate from the repository root:

```matlab
addpath scripts_run
run_omega_up_all_residual_parallel( ...
    'mode','interval', ...
    'workers',20, ...
    'eps_up','0.122', ...
    'Nx',976, ...
    'Ny',488, ...
    'Ny_axis',1952, ...
    'functional_scope','split', ...
    'jup_y_min','0.85');
```

## Outputs and certificate semantics

The atlas driver writes:

- `results/upper_conjecture/upper_atlas_interval_*.csv`, one row per cell;
- `results/upper_conjecture/upper_manifest_interval_*.json`, containing the
  thin result, atlas configuration, cover-completeness flag, verified-cell
  count, minimum margin, overlap, and output paths.

A run restricted by `max_cells` is intentionally incomplete.  It must not
be reported as an atlas certificate even if every selected cell passes.
Only interval mode with a complete cover and all cell margins positive sets
the atlas certificate flag.

The residual driver separately writes its checkpointed CSV/MAT/JSON outputs
below `results/omega_up_all_residual/`.  The theorem-level proof requires
both complete interval manifests and their recorded overlap.
