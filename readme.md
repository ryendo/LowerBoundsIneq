# LowerBoundsIneq

Computer-assisted proofs of sharp Dirichlet-Laplacian inequalities on
planar triangles.  The ver10 code combines verified finite-element
eigenvalue bounds, INTLAB interval arithmetic, and a residual enclosure for
second shape derivatives.

Companion code to

> R. Endo, X. Liu, P. Mariano, *Sharp Dirichlet Eigenvalue Inequalities on
> Triangles*.

## Statements covered by the code

For a triangle \(T\), let \(\lambda _1(T)\), \(A=|T|\), and \(L=|\partial T|\)
denote its first Dirichlet eigenvalue, area, and perimeter.  The lower-bound
part verifies

\[
J_1(T)=\lambda _1(T)A-\frac{\pi^2}{16}\frac{L^2}{A}
       -\frac{7\sqrt3\,\pi^2}{12}\geq0
\]

and

\[
J_2(T)=\lambda _1(T)A-
\frac{4\pi^2}{(3+\sqrt{\pi\sqrt3})^2}
\frac{(L+\sqrt{4\pi A})^2}{4A}\geq0.
\]

Both are sharp at the equilateral triangle.  The upper-bound part addresses
Siudeja's conjecture

\[
\lambda _1(T)\leq \frac{\pi^2L^2}{12A^2}
                       +\frac{\sqrt3\,\pi^2}{3A}.
\]

The normalized parameterization is
\(T(x,y)=\operatorname{conv}\{(0,0),(1,0),(x,y)\}\), with
\(1/2\leq x\leq1\), \(x^2+y^2\leq1\), and \(A=y/2\).

## Ver10 sharp strong-residual Hessian enclosure

Let \(\lambda_i\) be simple and let \(v\in H^2\cap H^1_0\) approximate its
material derivative.  With

\[
\mathcal J_i(v)=2\ell_i(v)-(a-\lambda_i m)(v,v),\qquad
g_{i,v}=\operatorname{div}(P_e\nabla\phi_i+\nabla v)
        +\alpha_i\phi_i+\lambda_i v,
\]

the exact complementarity identity is

\[
\ell_i(w_i)-\mathcal J_i(v)
=\sum_{k\ne i}\frac{|(g_{i,v},\phi_k)|^2}
{\lambda_k-\lambda_i}.
\]

Consequently

\[
c_i-2\mathcal J_i(v)-\frac{2\|g_{i,v}\|_0^2}
{\lambda_{i+1}-\lambda_i}
\leq\lambda_{i,ee}\leq
c_i-2\mathcal J_i(v)+\frac{2\|g_{i,v}\|_0^2}
{\lambda_i-\lambda_{i-1}},
\]

with the right residual term equal to zero for \(i=1\).  Both constants
are sharp: equality is attained by residuals in the adjacent eigenspaces.

The computable version replaces the exact eigenpair by certified errors
\(\varepsilon_\lambda,\varepsilon_a,\varepsilon_0\).  If
\(\widetilde g\) is the polynomial strong residual, the implementation
returns

\[
\widetilde D_i-\varepsilon_D-2E_i^+
\leq\lambda_{i,ee}\leq
\widetilde D_i+\varepsilon_D+2E_i^-,
\]

where

\[
E_i^\pm=
\left(
\frac{\|\widetilde g\|_0+\delta_0}{\sqrt{d_i^\pm}}
+\sqrt{C_i^\pm}\,\delta_a
\right)^2,\qquad E_1^-=0.
\]

It uses no eigenfunction series, high-mode cutoff, spectral tail, or upper
cluster dimension.  For fixed positive adjacent gaps, its width converges
to zero when the eigenvalue/eigenfunction errors and strong residual
converge to zero.

For \(\lambda_1\), the only verified spectral endpoints are

\[
L_1^{\rm LG}\leq\lambda_1\leq U_1^{\rm Ritz}
<L_2^{\rm CR}\leq\lambda_2.
\]

Here \(L_1^{\rm LG}\) is Lehmann--Goerisch, \(U_1^{\rm Ritz}\) is a
verified conforming Rayleigh--Ritz value, and \(L_2^{\rm CR}\) is the
corrected Crouzeix--Raviart--Liu lower bound.  No enclosure of
\(\lambda_3\) is used.  The same certificate routine accepts an arbitrary
index \(i\), returning LG/Ritz enclosures through \(i\) and one CR--Liu
lower endpoint for \(\lambda_{i+1}\); it never asks for
\(\lambda_{i+2}\).  The \(\Omega_{\rm up}\) Hessian trials are global
degree-11 Bernstein bubbles, so their mass, stiffness, Hessian, and strong
residual norms use exact barycentric moments and outward rounding.  The
Hessian step itself has no mesh, numerical quadrature, or flux
reconstruction; CR/LG is confined to the independent scalar endpoint
certificate.

The main entry points are:

- `src/fem/strong_residual_hessian_enclosure.m`
- `src/fem/calc_ddlambda1_bernstein_strong_bounds.m`
- `src/fem/triangle_spectral_certificate_cr_lg.m`
- `src/algorithms/build_omega_up_spectral_atlas_cr_lg.m`
- `src/algorithms/Algorithm2_VerifyOmegaUpResidual.m`

## Certificate decomposition and current committed results

The proof is split into independently checkable regions.

### \(\Omega_{\rm mid}\) certificate

`inputs/cell_def.csv` contains **188,623** cells.  The authoritative ver10
driver is `scripts_run/run_omega_mid_unified_parallel.m`; the production
wrapper is `scripts_run/run_omega_mid_unified_full.m`.  It computes \(J_1\)
and \(J_2\) together, reusing one strict CR/Liu or
Liu--Lehmann--Goerisch spectral enclosure per cell.  In the LG branch it
records and checks the conforming Ritz separation, interval positivity of
the LG Gram matrix, negativity of the transformed eigenvalues, and
finiteness of the resulting lower bound.

The historical committed CSVs are not a ver10 certificate: each contains
two `J_lower=Inf` rows that the former validator accepted.  The current
validator deliberately rejects them.  Final counts and extrema will be
quoted here only after the unified clean-commit interval run has atomically
replaced both CSVs and regenerated `results/verification_summary.json`.

### Ver10 \(\Omega_{\rm up}\) residual production run

`scripts_run/run_omega_up_all_residual_parallel.m` is the resumable
production driver.  Its production grid is

- `eps_up=0.122`, `Nx=976`, and `Ny=488` for \(x\)-direction rectangles;
- `Ny_axis=1952`, hence axis width \(0.122/1952=0.0000625\), for the
  \(y\)-direction symmetry-axis intervals;
- degree-11 Bernstein eigenfunction/material trials and exact polynomial
  strong-residual norms;
- a reusable coarse atlas of
  \(L_1^{\rm LG},U_1^{\rm Ritz},L_2^{\rm CR}\), transported to each fine
  cell by rigorous \(2\times2\) affine metric factors;
- `functional_scope='split'`: \(J_1,J_2\) on all of
  \(\Omega_{\rm up}\), and the Siudeja upper functional only on
  \(y\geq0.85\).

One Hessian-estimator call per directional cell is reused for every
in-scope functional.  Rectangles outside \(x^2+y^2\leq1\) are recorded as
geometric skips.  Each task has an atomic MAT checkpoint, so an interrupted
HPC run can be resumed.  The merged CSV, MAT file, and JSON manifest record
the configuration, Git commit, functional scopes, extrema, failures, and
whether coverage is complete.

No final \(\Omega_{\rm up}\) extremum is quoted here until a complete
interval-mode production manifest has been generated and committed.  A
`max_tasks` or `task_ids` run is a smoke test and has
`coverage_complete=false`; it is not a proof of the whole region.

### Siudeja upper-bound certificate

The upper conjecture is split into:

1. \(0<y\leq0.06\): an analytic sector-inclusion certificate using an
   explicit radial trial function, evaluated after nonsingular \(y^2\)
   scaling; no Bessel zero is numerically evaluated.
2. \(0.06\leq y\leq0.851\): a conforming P4/\(N=8\) fixed-trial atlas,
   affinely transported over each cell and enclosed by a centered
   mean-value formula in INTLAB.
3. \(0.85\leq y\leq\sqrt3/2\): the independent
   \(\Omega_{\rm up}\) residual-Hessian certificate above.

The overlap \(0.85\leq y\leq0.851\) is intentional and is recorded in the
atlas manifest.  See `src/upper_conjecture/README.md` for the proof formulas
and reproduction commands.

Historical `results/J*_OmegaUp*.csv` and
`results/omega_connection_summary.*` predate the ver10 residual production
driver.  Likewise, files below `results/omega_up_taylor_*` are
double-precision exploratory diagnostics.  They are not the ver10
\(\Omega_{\rm up}\) or upper-conjecture proof certificate.

## Dependencies

- MATLAB R2023b or later (the production run used R2023b).
- Parallel Computing Toolbox for multi-worker production runs.
- Gmsh (tested with 4.8.4).
- An external INTLAB installation for every rigorous run.
- The revised VFEM2D and VEIGS/VEIG routines under `src/lib/`.

INTLAB is **not vendored in this repository**.  Set `INTLAB_ROOT` to the
installation directory before starting MATLAB.  The supported reproducible
configuration also sets separate mesh storage.  `my_intlab_config.m`
defaults to `/usr/bin/gmsh`; set `GMSH_COMMAND` when Gmsh is installed
elsewhere:

```sh
export INTLAB_ROOT=/path/to/Intlab_V12
export LOWERBOUNDS_MESH_PATH=/path/to/private/mesh-workspace
export GMSH_COMMAND=/usr/bin/gmsh
```

Each parallel worker initializes the same `INTLAB_ROOT`.  Temporary meshes
are kept outside the checkout so workers do not share generated files.

## Repository layout

```text
LowerBoundsIneq/
├── VerifyTriangleInequalities.m
├── inputs/
│   └── cell_def.csv
├── results/
│   ├── J1_OmegaMid.csv
│   ├── J2_OmegaMid.csv
│   ├── verification_summary.md
│   └── verification_summary.json
├── scripts_run/
│   ├── run_omega_mid_unified_parallel.m
│   ├── run_omega_mid_unified_full.m
│   ├── run_omega_up_all_residual_parallel.m
│   ├── run_omega_up_spectral_atlas.m
│   ├── run_omega_up_residual_full.m
│   ├── run_omega_up_residual_smoke.m
│   ├── run_upper_conjecture_global_finalize.m
│   └── run_upper_conjecture_atlas_full.m
├── src/
│   ├── algorithms/
│   ├── fem/
│   ├── interval/
│   ├── lib/
│   ├── mesh/
│   └── upper_conjecture/
└── tests/
```

## Validate the published \(\Omega_{\rm mid}\) results

This read-only validation does not require INTLAB:

```matlab
v = VerifyTriangleInequalities();
v.run();
```

It checks the exact cell-ID set and row count, `status=ok`, `verified=1`,
and finite strict positivity of all 188,623 lower bounds for each
functional.  It fails closed on the historical pre-ver10 CSVs.

## Rigorous smoke tests

After setting the external runtime variables:

```matlab
my_intlab_config
test_verified_ritz_enclosures(1)
test_bernstein_strong_residual_estimator(1)
```

The last command checks CR/LG endpoint separation and both a point and a
finite-width Bernstein Hessian cell.  It deliberately does not claim full
regional coverage.  Targeted production-width checks are in
`scripts_run/run_omega_up_rectangle_targeted_smoke.m` and
`scripts_run/run_omega_up_axis_targeted_smoke.m`.

## Reproduce the production computations

The authoritative \(\Omega_{\rm mid}\) atlas is recomputed with:

```matlab
run('scripts_run/run_omega_mid_unified_full.m')
```

The wrapper requires a clean Git checkout, uses 20 workers and resumable
chunk checkpoints, validates exact input coverage, and publishes the J1/J2
CSV pair transactionally only after every strict certificate succeeds.
`VerifyTriangleInequalities().compute(20)` delegates to this same wrapper.

First create the reusable CR/LG endpoint atlas outside the checkout:

```sh
export VER10_SPECTRAL_ATLAS=/absolute/certificates/omega_up_spectral_atlas.mat
export VER10_CERTIFICATE_ROOT=/absolute/certificates
matlab -batch "run('scripts_run/run_omega_up_spectral_atlas.m')"
```

Then run the fine Hessian cover:

```sh
matlab -batch "run('scripts_run/run_omega_up_residual_full.m')"
```

The Siudeja compact atlas is:

```matlab
addpath scripts_run
run('scripts_run/run_upper_conjecture_atlas_full.m')
```

After the compact and residual manifests exist, set
`VER10_COMPACT_UPPER_MANIFEST` and `VER10_OMEGA_UP_MANIFEST`, then run
`scripts_run/run_upper_conjecture_global_finalize.m`.  The finalizer
checks both CSV hashes, strict margins, source commits, complete coverage,
and the nonempty \(0.85\leq y\leq0.851\) overlap before emitting a global
certificate.

On a shared 96-core node, keep threaded BLAS single-threaded so that the 20
MATLAB workers do not oversubscribe the node:

```sh
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
matlab -batch "run('scripts_run/run_omega_up_residual_full.m')"
```

Use a stable `run_name` to resume an interrupted \(\Omega_{\rm up}\) run;
the driver rejects checkpoints whose configuration fingerprint differs.

## Input and output formats

`inputs/cell_def.csv` has one row per \(\Omega_{\rm mid}\) cell:

| column | meaning |
|---|---|
| `i` | cell identifier |
| `x_inf`, `x_sup` | apex \(x\)-interval |
| `theta_inf`, `theta_sup` | apex-angle interval |
| `mesh_size_lower_cr` | CR lower-bound mesh size |
| `isLG` | whether Lehmann--Goerisch refinement is used |
| `mesh_size_lower_LG` | LG mesh size |
| `fem_order_lower_LG` | LG finite-element order |

`results/J1_OmegaMid.csv` and `results/J2_OmegaMid.csv` retain the legacy
columns `conjecture`, `cell_id`, `verified`, `J_lower`, `status`, `note`,
and `run_timestamp`, followed by the shared spectral method,
\(\lambda_1\) lower bound, LG shift/Ritz/separation/transformed-eigenvalue
diagnostics, strictness flags, timing, and error fields.

The residual \(\Omega_{\rm up}\) driver writes one checkpoint per task and,
after merging, `omega_up_all_cells.csv`, `omega_up_all_results.mat`, and
`omega_up_all_manifest.json`.  Only an interval manifest with complete
coverage and every requested functional certified is a regional proof.

The upper-conjecture driver writes `upper_atlas_interval_*.csv` and
`upper_manifest_interval_*.json`.  Its atlas manifest covers only through
its recorded `y_up`; completion of the global upper conjecture also
requires the overlapping residual certificate.

## Citation

If you use this code, please cite:

- R. Endo, X. Liu, P. Mariano, *Sharp Dirichlet Eigenvalue Inequalities on
  Triangles*.

## Licenses

The external libraries under `src/lib/` retain their original licenses.
INTLAB is separately installed and remains subject to its own license.
