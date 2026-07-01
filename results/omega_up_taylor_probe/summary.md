# Omega Up Taylor Probe

This is an exploratory double-precision Taylor spectral probe. It does not modify the production Omega_up or Omega_mid verification code.

The spectral data use analytic equilateral triangle eigenfunctions, but spectral inner products and Rxx cell lower bounds are evaluated numerically in double precision. Therefore the rigorous flag is `exploratory_double`.

Known comparison points:

- Current cellwise eigenspace enclosure succeeds only around distance `1e-5`.
- Idealized fixed-equilateral-coefficient model suggests `1e-2` to `1e-1` scale.
- This Taylor probe tests whether the third-derivative Taylor certificate can bridge that gap.

A row that reaches the configured `--max-radius` search cap should be read as a lower diagnostic, not as the true maximal accepted radius.

Taylor models:

- `norm_bound`: existing first-order Taylor lower bound using `-rx M_xxx - ry M_xxy`.
- `symmetry_x_even`: `lambda_xx(0,0)` lower bound with the `rx M_xxx` penalty removed by equilateral x-symmetry, keeping `-ry M_xxy`.
- `signed_xxy`: `lambda_xx(0,0)` lower bound plus the signed affine `lambda_xxy t` minimum on `t in [-ry,0]`.
- `dirichlet_parseval_tail`: the intentionally crude Dirichlet-only Parseval tail diagnostic; it is not used by the two symmetry-improved models.

Rigor status:

- The Neumann-corrected tail formula for `lambda_xx(0,0)` is the rigorous direction mathematically, but this script still evaluates all integrals in double precision.
- The `symmetry_x_even` and `signed_xxy` cell models are exploratory because the fourth-order `s^2`, `st`, and `t^2` remainders are not yet bounded.
- The `M_xxy` norm bound remains conservative; for `neumann_corrected_tail` only the `P_x` tail entering `lambda_xx` is Neumann-corrected.

Finite signed derivative diagnostics:

| J | cluster complete | complete cutoffs | lambda_x | lambda_xx finite | lambda_xxx finite | lambda_xxy finite | M_xxx/abs | M_xxy/abs |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 30 | True | 30/30 | -4.88786096e-13 | 3.63934075e+01 | 7.20151797e-13 | -7.50212457e+01 | 2.10126127e+15 | 2.57191748e+01 |
| 49 | True | 49/49 | -4.88786096e-13 | 3.63259650e+01 | 7.73517510e-13 | -7.51026047e+01 | 1.95760684e+15 | 2.57059623e+01 |
| 51 | True | 51/51 | -4.88786096e-13 | 3.63242373e+01 | 7.36215552e-13 | -7.51243643e+01 | 2.05682825e+15 | 2.56988879e+01 |
| 99 | True | 99/99 | -4.88786096e-13 | 3.62754001e+01 | 6.88287578e-13 | -7.51588564e+01 | 2.20110209e+15 | 2.56974916e+01 |
| 101 | True | 101/101 | -4.88786096e-13 | 3.62754001e+01 | 6.88287578e-13 | -7.51588564e+01 | 2.20110209e+15 | 2.56974916e+01 |
| 199 | True | 199/199 | -4.88786096e-13 | 3.62550042e+01 | 6.35004396e-13 | -7.51707173e+01 | 2.38626478e+15 | 2.56977148e+01 |
| 203 | True | 203/203 | -4.88786096e-13 | 3.62535315e+01 | 6.29675667e-13 | -7.51776145e+01 | 2.40649281e+15 | 2.56956640e+01 |

The finite `lambda_xxx` values are numerically close to zero at the `1e-6` level.

Basis orthogonality diagnostics:

| bc | J | mass max error | stiffness max error | stiffness relative error |
|---|---:|---:|---:|---:|
| D | 30 | 3.24185123e-14 | 3.52429197e-11 | 3.29279601e-14 |
| D | 49 | 3.44169138e-14 | 5.25233190e-11 | 3.21878519e-14 |
| D | 51 | 3.44169138e-14 | 5.25233190e-11 | 3.08605178e-14 |
| D | 99 | 3.44169138e-14 | 1.33240974e-10 | 4.14963110e-14 |
| D | 101 | 3.44169138e-14 | 1.33240974e-10 | 4.01789678e-14 |
| D | 199 | 3.44169138e-14 | 2.68300937e-10 | 4.21248479e-14 |
| D | 203 | 3.44169138e-14 | 2.68300937e-10 | 4.20091203e-14 |
| N | 29 | 4.90533275e-14 | 3.71755959e-11 | 5.72636460e-14 |
| N | 31 | 4.90533275e-14 | 3.71755959e-11 | 5.43270488e-14 |
| N | 50 | 4.90533275e-14 | 5.83213478e-11 | 4.96107200e-14 |
| N | 98 | 4.90533275e-14 | 1.06410880e-10 | 4.21159231e-14 |
| N | 101 | 4.90533275e-14 | 1.25055521e-10 | 4.84851025e-14 |
| N | 200 | 4.90533275e-14 | 2.27373675e-10 | 4.30523109e-14 |

Shell-wise `lambda_xx` finite sum diagnostic:

| shell | q | multiplicity | complete J | contribution | cumulative | lambda_xx finite | fraction of final cumulative |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 7 | 2 | 3 | 1.44874865e+01 | 1.44874865e+01 | 4.12088806e+01 | 8.53955140e-01 |
| 2 | 12 | 1 | 4 | 2.25210728e-27 | 1.44874865e+01 | 4.12088806e+01 | 8.53955140e-01 |
| 3 | 13 | 2 | 6 | 1.81731030e+00 | 1.63047968e+01 | 3.75742600e+01 | 9.61075272e-01 |
| 4 | 19 | 2 | 8 | 2.34673334e-01 | 1.65394701e+01 | 3.71049133e+01 | 9.74907934e-01 |
| 5 | 21 | 2 | 10 | 1.35543607e-27 | 1.65394701e+01 | 3.71049133e+01 | 9.74907934e-01 |
| 6 | 27 | 1 | 11 | 1.90259250e-27 | 1.65394701e+01 | 3.71049133e+01 | 9.74907934e-01 |
| 7 | 28 | 2 | 13 | 1.90178250e-01 | 1.67296483e+01 | 3.67245568e+01 | 9.86117863e-01 |
| 8 | 31 | 2 | 15 | 4.22375699e-02 | 1.67718859e+01 | 3.66400817e+01 | 9.88607528e-01 |
| 9 | 37 | 2 | 17 | 4.03410692e-02 | 1.68122270e+01 | 3.65593995e+01 | 9.90985405e-01 |
| 10 | 39 | 2 | 19 | 1.11757096e-27 | 1.68122270e+01 | 3.65593995e+01 | 9.90985405e-01 |
| 11 | 43 | 2 | 21 | 2.87474835e-03 | 1.68151017e+01 | 3.65536501e+01 | 9.91154856e-01 |
| 12 | 48 | 1 | 22 | 1.87879601e-27 | 1.68151017e+01 | 3.65536501e+01 | 9.91154856e-01 |
- Full shell table is written to `lambda_xx_shell_contributions.csv`.

Neumann-corrected tail diagnostic:

| J | M | old tail energy | corrected tail energy | ratio | old lambda_xx lower | corrected lambda_xx lower |
|---:|---:|---:|---:|---:|---:|---:|
| 30 | 31 | 5.99866146e+01 | 5.78366856e+00 | 9.64159854e-02 | -8.95784831e+01 | 2.42477035e+01 |
| 30 | 50 | 5.99866146e+01 | 4.56451831e+00 | 7.60922807e-02 | -8.95784831e+01 | 2.68079190e+01 |
| 30 | 101 | 5.99866146e+01 | 3.22820126e+00 | 5.38153601e-02 | -8.95784831e+01 | 2.96141848e+01 |
| 30 | 200 | 5.99866146e+01 | 2.17754094e+00 | 3.63004473e-02 | -8.95784831e+01 | 3.18205715e+01 |
| 51 | 31 | 5.99533506e+01 | 5.75040454e+00 | 9.59146485e-02 | -8.71796648e+01 | 2.44784040e+01 |
| 51 | 50 | 5.99533506e+01 | 4.53125430e+00 | 7.55796675e-02 | -8.71796648e+01 | 2.69898535e+01 |
| 51 | 101 | 5.99533506e+01 | 3.19493725e+00 | 5.32903870e-02 | -8.71796648e+01 | 2.97426666e+01 |
| 51 | 200 | 5.99533506e+01 | 2.14427693e+00 | 3.57657564e-02 | -8.71796648e+01 | 3.19070269e+01 |
| 101 | 31 | 5.99295212e+01 | 5.72657516e+00 | 9.55551628e-02 | -8.54861668e+01 | 2.46404537e+01 |
| 101 | 50 | 5.99295212e+01 | 4.50742491e+00 | 7.52120962e-02 | -8.54861668e+01 | 2.71174574e+01 |
| 101 | 101 | 5.99295212e+01 | 3.17110786e+00 | 5.29139529e-02 | -8.54861668e+01 | 2.98325143e+01 |
| 101 | 200 | 5.99295212e+01 | 2.12044754e+00 | 3.53823541e-02 | -8.54861668e+01 | 3.19671892e+01 |
| 203 | 31 | 5.99187262e+01 | 5.71578013e+00 | 9.53922170e-02 | -8.45715922e+01 | 2.47277550e+01 |
| 203 | 50 | 5.99187262e+01 | 4.49662988e+00 | 7.50454853e-02 | -8.45715922e+01 | 2.71861514e+01 |
| 203 | 101 | 5.99187262e+01 | 3.16031284e+00 | 5.27433248e-02 | -8.45715922e+01 | 2.98808127e+01 |
| 203 | 200 | 5.99187262e+01 | 2.10965251e+00 | 3.52085675e-02 | -8.45715922e+01 | 3.19994519e+01 |
- Smallest corrected tail in this run: `2.10965251e+00` at `J=203`, `M=200`.
- Full comparison is written to `neumann_corrected_tail.csv`.

Cell-size bisection results:

| model | J | M | tail mode | cell shape | max rx | max ry | distance | Jxx1 lower | Jxx2 lower | rigorous flag |
|---|---:|---:|---|---|---:|---:|---:|---:|---:|---|
| norm_bound | 30 |  | tail_ignored | x_only | 4.48081823e-03 | 0.00000000e+00 | 4.48081823e-03 | 1.48524942e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 |  | tail_ignored | y_only | 0.00000000e+00 | 3.40823079e-03 | 3.40823079e-03 | 1.55790964e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 |  | tail_ignored | omega_up | 3.96025396e-03 | 3.96025396e-04 | 3.98000597e-03 | 1.50399992e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 |  | tail_ignored | square | 1.93590103e-03 | 1.93590103e-03 | 2.73777749e-03 | 1.55122141e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 |  | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.16096132e+01 | -5.16111829e+01 | exploratory_double |
| norm_bound | 30 |  | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.16096132e+01 | -5.16111829e+01 | exploratory_double |
| norm_bound | 30 |  | dirichlet_parseval_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.16096132e+01 | -5.16111829e+01 | exploratory_double |
| norm_bound | 30 |  | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.16096132e+01 | -5.16111829e+01 | exploratory_double |
| norm_bound | 30 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| norm_bound | 30 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| norm_bound | 30 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| norm_bound | 30 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| norm_bound | 30 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| norm_bound | 30 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| norm_bound | 30 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| norm_bound | 30 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| norm_bound | 30 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.21282276e+00 | -1.21439243e+00 | exploratory_double |
| norm_bound | 30 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.21282276e+00 | -1.21439243e+00 | exploratory_double |
| norm_bound | 30 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.21282276e+00 | -1.21439243e+00 | exploratory_double |
| norm_bound | 30 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.21282276e+00 | -1.21439243e+00 | exploratory_double |
| norm_bound | 30 | 98 | neumann_corrected_tail | x_only | 1.08194902e-06 | 0.00000000e+00 | 1.08194902e-06 | 1.56967631e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 6.91245027e-07 | 6.91245027e-07 | 1.56967725e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | 98 | neumann_corrected_tail | omega_up | 9.35519802e-07 | 9.35519802e-08 | 9.40185765e-07 | 1.56967644e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | 98 | neumann_corrected_tail | square | 4.21776674e-07 | 4.21776674e-07 | 5.96482292e-07 | 1.56967689e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | 101 | neumann_corrected_tail | x_only | 1.08194902e-06 | 0.00000000e+00 | 1.08194902e-06 | 1.56967631e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 6.91245027e-07 | 6.91245027e-07 | 1.56967725e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | 101 | neumann_corrected_tail | omega_up | 9.35519802e-07 | 9.35519802e-08 | 9.40185765e-07 | 1.56967644e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | 101 | neumann_corrected_tail | square | 4.21776674e-07 | 4.21776674e-07 | 5.96482292e-07 | 1.56967689e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | 200 | neumann_corrected_tail | x_only | 1.39584724e-03 | 0.00000000e+00 | 1.39584724e-03 | 1.56148338e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 8.83103419e-04 | 8.83103419e-04 | 1.56978163e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | 200 | neumann_corrected_tail | omega_up | 1.20533756e-03 | 1.20533756e-04 | 1.21134925e-03 | 1.56370937e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | 200 | neumann_corrected_tail | square | 5.40903114e-04 | 5.40903114e-04 | 7.64952520e-04 | 1.56876962e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 30 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.40091519e+00 | 2.56929486e+00 | exploratory_double |
| symmetry_x_even | 30 |  | tail_ignored | y_only | 0.00000000e+00 | 3.40823079e-03 | 3.40823079e-03 | 1.55790964e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 30 |  | tail_ignored | omega_up | 3.39044348e-02 | 3.39044348e-03 | 3.40735353e-02 | 1.77635684e-15 | 3.30583722e-03 | exploratory_double |
| symmetry_x_even | 30 |  | tail_ignored | square | 3.40808954e-03 | 3.40808954e-03 | 4.81976645e-03 | 1.50878081e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 30 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.40091519e+00 | 2.56929486e+00 | exploratory_double |
| signed_xxy | 30 |  | tail_ignored | y_only | 0.00000000e+00 | 9.88637327e-02 | 9.88637327e-02 | 1.77635684e-15 | 1.52399842e-02 | exploratory_double |
| signed_xxy | 30 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.75747456e+00 | 1.93281022e+00 | exploratory_double |
| signed_xxy | 30 |  | tail_ignored | square | 9.24585593e-02 | 9.24585593e-02 | 1.30756149e-01 | 8.88178420e-15 | 5.52099788e-02 | exploratory_double |
| symmetry_x_even | 30 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| symmetry_x_even | 30 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| symmetry_x_even | 30 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| symmetry_x_even | 30 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| signed_xxy | 30 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| signed_xxy | 30 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| signed_xxy | 30 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| signed_xxy | 30 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| symmetry_x_even | 30 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| symmetry_x_even | 30 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| symmetry_x_even | 30 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| symmetry_x_even | 30 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| signed_xxy | 30 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| signed_xxy | 30 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| signed_xxy | 30 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| signed_xxy | 30 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.32142859e+00 | -2.32299827e+00 | exploratory_double |
| symmetry_x_even | 30 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.21282276e+00 | -1.21439243e+00 | exploratory_double |
| symmetry_x_even | 30 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.21282276e+00 | -1.21439243e+00 | exploratory_double |
| symmetry_x_even | 30 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.21282276e+00 | -1.21439243e+00 | exploratory_double |
| symmetry_x_even | 30 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.21282276e+00 | -1.21439243e+00 | exploratory_double |
| signed_xxy | 30 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.21282276e+00 | -1.21439243e+00 | exploratory_double |
| signed_xxy | 30 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.21282276e+00 | -1.21439243e+00 | exploratory_double |
| signed_xxy | 30 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.21282276e+00 | -1.21439243e+00 | exploratory_double |
| signed_xxy | 30 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.21282276e+00 | -1.21439243e+00 | exploratory_double |
| symmetry_x_even | 30 | 98 | neumann_corrected_tail | x_only | 8.60386500e-03 | 0.00000000e+00 | 8.60386500e-03 | 1.25838955e-03 | 5.32907052e-15 | exploratory_double |
| symmetry_x_even | 30 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 6.91245027e-07 | 6.91245027e-07 | 1.56967725e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 30 | 98 | neumann_corrected_tail | omega_up | 6.91244581e-06 | 6.91244581e-07 | 6.94692206e-06 | 1.56967705e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 30 | 98 | neumann_corrected_tail | square | 6.91245023e-07 | 6.91245023e-07 | 9.77568086e-07 | 1.56967725e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 30 | 98 | neumann_corrected_tail | x_only | 8.60386500e-03 | 0.00000000e+00 | 8.60386500e-03 | 1.25838955e-03 | 5.32907052e-15 | exploratory_double |
| signed_xxy | 30 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 2.91875268e-05 | 2.91875268e-05 | 1.56971480e-03 | 5.32907052e-15 | exploratory_double |
| signed_xxy | 30 | 98 | neumann_corrected_tail | omega_up | 2.91540018e-04 | 2.91540018e-05 | 2.92994092e-04 | 1.56935733e-03 | 5.32907052e-15 | exploratory_double |
| signed_xxy | 30 | 98 | neumann_corrected_tail | square | 2.91871908e-05 | 2.91871908e-05 | 4.12769210e-05 | 1.56971121e-03 | 5.32907052e-15 | exploratory_double |
| symmetry_x_even | 30 | 101 | neumann_corrected_tail | x_only | 8.60386500e-03 | 0.00000000e+00 | 8.60386500e-03 | 1.25838955e-03 | 5.32907052e-15 | exploratory_double |
| symmetry_x_even | 30 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 6.91245027e-07 | 6.91245027e-07 | 1.56967725e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 30 | 101 | neumann_corrected_tail | omega_up | 6.91244581e-06 | 6.91244581e-07 | 6.94692206e-06 | 1.56967705e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 30 | 101 | neumann_corrected_tail | square | 6.91245023e-07 | 6.91245023e-07 | 9.77568086e-07 | 1.56967725e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 30 | 101 | neumann_corrected_tail | x_only | 8.60386500e-03 | 0.00000000e+00 | 8.60386500e-03 | 1.25838955e-03 | 5.32907052e-15 | exploratory_double |
| signed_xxy | 30 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 2.91875268e-05 | 2.91875268e-05 | 1.56971480e-03 | 5.32907052e-15 | exploratory_double |
| signed_xxy | 30 | 101 | neumann_corrected_tail | omega_up | 2.91540018e-04 | 2.91540018e-05 | 2.92994092e-04 | 1.56935733e-03 | 5.32907052e-15 | exploratory_double |
| signed_xxy | 30 | 101 | neumann_corrected_tail | square | 2.91871908e-05 | 2.91871908e-05 | 4.12769210e-05 | 1.56971121e-03 | 5.32907052e-15 | exploratory_double |
| symmetry_x_even | 30 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.20819127e-01 | 5.89198799e-01 | exploratory_double |
| symmetry_x_even | 30 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 8.83103419e-04 | 8.83103419e-04 | 1.56978163e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 30 | 200 | neumann_corrected_tail | omega_up | 8.82363988e-03 | 8.82363988e-04 | 8.86764833e-03 | 1.24188900e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 30 | 200 | neumann_corrected_tail | square | 8.83096011e-04 | 8.83096011e-04 | 1.24888636e-03 | 1.56649733e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 30 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.20819127e-01 | 5.89198799e-01 | exploratory_double |
| signed_xxy | 30 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 3.52403067e-02 | 3.52403067e-02 | 7.10542736e-15 | 2.58832763e-04 | exploratory_double |
| signed_xxy | 30 | 200 | neumann_corrected_tail | omega_up | 1.78513708e-01 | 1.78513708e-02 | 1.79404056e-01 | 1.77635684e-15 | 1.38501340e-01 | exploratory_double |
| signed_xxy | 30 | 200 | neumann_corrected_tail | square | 3.44972939e-02 | 3.44972939e-02 | 4.87865409e-02 | 1.77635684e-15 | 5.49107682e-03 | exploratory_double |
| norm_bound | 49 |  | tail_ignored | x_only | 4.43327938e-03 | 0.00000000e+00 | 4.43327938e-03 | 1.48703137e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 |  | tail_ignored | y_only | 0.00000000e+00 | 3.37247287e-03 | 3.37247287e-03 | 1.55820456e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 |  | tail_ignored | omega_up | 3.91829128e-03 | 3.91829128e-04 | 3.93783400e-03 | 1.50539047e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 |  | tail_ignored | square | 1.91549070e-03 | 1.91549070e-03 | 2.70891293e-03 | 1.55163672e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 |  | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.06703179e+01 | -5.06718875e+01 | exploratory_double |
| norm_bound | 49 |  | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.06703179e+01 | -5.06718875e+01 | exploratory_double |
| norm_bound | 49 |  | dirichlet_parseval_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.06703179e+01 | -5.06718875e+01 | exploratory_double |
| norm_bound | 49 |  | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.06703179e+01 | -5.06718875e+01 | exploratory_double |
| norm_bound | 49 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| norm_bound | 49 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| norm_bound | 49 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| norm_bound | 49 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| norm_bound | 49 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| norm_bound | 49 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| norm_bound | 49 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| norm_bound | 49 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| norm_bound | 49 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.14155722e+00 | -1.14312689e+00 | exploratory_double |
| norm_bound | 49 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.14155722e+00 | -1.14312689e+00 | exploratory_double |
| norm_bound | 49 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.14155722e+00 | -1.14312689e+00 | exploratory_double |
| norm_bound | 49 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.14155722e+00 | -1.14312689e+00 | exploratory_double |
| norm_bound | 49 | 98 | neumann_corrected_tail | x_only | 7.32619692e-05 | 0.00000000e+00 | 7.32619692e-05 | 1.56965375e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 4.70265782e-05 | 4.70265782e-05 | 1.56973714e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 | 98 | neumann_corrected_tail | omega_up | 6.33870254e-05 | 6.33870254e-06 | 6.37031721e-05 | 1.56966798e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 | 98 | neumann_corrected_tail | square | 2.86416443e-05 | 2.86416443e-05 | 4.05054018e-05 | 1.56971065e-03 | 3.55271368e-15 | exploratory_double |
| norm_bound | 49 | 101 | neumann_corrected_tail | x_only | 7.32619692e-05 | 0.00000000e+00 | 7.32619692e-05 | 1.56965375e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 4.70265782e-05 | 4.70265782e-05 | 1.56973714e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 | 101 | neumann_corrected_tail | omega_up | 6.33870254e-05 | 6.33870254e-06 | 6.37031721e-05 | 1.56966798e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 | 101 | neumann_corrected_tail | square | 2.86416443e-05 | 2.86416443e-05 | 4.05054018e-05 | 1.56971065e-03 | 3.55271368e-15 | exploratory_double |
| norm_bound | 49 | 200 | neumann_corrected_tail | x_only | 1.44765192e-03 | 0.00000000e+00 | 1.44765192e-03 | 1.56086396e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 9.20519696e-04 | 9.20519696e-04 | 1.56973763e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 | 200 | neumann_corrected_tail | omega_up | 1.25093292e-03 | 1.25093292e-04 | 1.25717203e-03 | 1.56324291e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 | 200 | neumann_corrected_tail | square | 5.62715967e-04 | 5.62715967e-04 | 7.95800553e-04 | 1.56866410e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 49 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.37171176e+00 | 2.54009143e+00 | exploratory_double |
| symmetry_x_even | 49 |  | tail_ignored | y_only | 0.00000000e+00 | 3.37247287e-03 | 3.37247287e-03 | 1.55820456e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 49 |  | tail_ignored | omega_up | 3.35510519e-02 | 3.35510519e-03 | 3.37183899e-02 | 1.77635684e-15 | 3.20435612e-03 | exploratory_double |
| symmetry_x_even | 49 |  | tail_ignored | square | 3.37233468e-03 | 3.37233468e-03 | 4.76920145e-03 | 1.51010413e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 49 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.37171176e+00 | 2.54009143e+00 | exploratory_double |
| signed_xxy | 49 |  | tail_ignored | y_only | 0.00000000e+00 | 9.80019278e-02 | 9.80019278e-02 | 7.10542736e-15 | 1.49158462e-02 | exploratory_double |
| signed_xxy | 49 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.72894555e+00 | 1.90428121e+00 | exploratory_double |
| signed_xxy | 49 |  | tail_ignored | square | 9.17137311e-02 | 9.17137311e-02 | 1.29702802e-01 | 3.55271368e-15 | 5.42180721e-02 | exploratory_double |
| symmetry_x_even | 49 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| symmetry_x_even | 49 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| symmetry_x_even | 49 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| symmetry_x_even | 49 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| signed_xxy | 49 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| signed_xxy | 49 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| signed_xxy | 49 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| signed_xxy | 49 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| symmetry_x_even | 49 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| symmetry_x_even | 49 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| symmetry_x_even | 49 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| symmetry_x_even | 49 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| signed_xxy | 49 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| signed_xxy | 49 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| signed_xxy | 49 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| signed_xxy | 49 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.23106852e+00 | -2.23263820e+00 | exploratory_double |
| symmetry_x_even | 49 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.14155722e+00 | -1.14312689e+00 | exploratory_double |
| symmetry_x_even | 49 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.14155722e+00 | -1.14312689e+00 | exploratory_double |
| symmetry_x_even | 49 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.14155722e+00 | -1.14312689e+00 | exploratory_double |
| symmetry_x_even | 49 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.14155722e+00 | -1.14312689e+00 | exploratory_double |
| signed_xxy | 49 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.14155722e+00 | -1.14312689e+00 | exploratory_double |
| signed_xxy | 49 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.14155722e+00 | -1.14312689e+00 | exploratory_double |
| signed_xxy | 49 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.14155722e+00 | -1.14312689e+00 | exploratory_double |
| signed_xxy | 49 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.14155722e+00 | -1.14312689e+00 | exploratory_double |
| symmetry_x_even | 49 | 98 | neumann_corrected_tail | x_only | 6.06214704e-02 | 0.00000000e+00 | 6.06214704e-02 | 3.55271368e-15 | 1.39009613e-02 | exploratory_double |
| symmetry_x_even | 49 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 4.70265782e-05 | 4.70265782e-05 | 1.56973714e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 49 | 98 | neumann_corrected_tail | omega_up | 4.70244977e-04 | 4.70244977e-05 | 4.72590353e-04 | 1.56880722e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 49 | 98 | neumann_corrected_tail | square | 4.70265574e-05 | 4.70265574e-05 | 6.65055953e-05 | 1.56972784e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 49 | 98 | neumann_corrected_tail | x_only | 6.06214704e-02 | 0.00000000e+00 | 6.06214704e-02 | 3.55271368e-15 | 1.39009613e-02 | exploratory_double |
| signed_xxy | 49 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 1.96685977e-03 | 1.96685977e-03 | 1.56690501e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 49 | 98 | neumann_corrected_tail | omega_up | 1.83328332e-02 | 1.83328332e-03 | 1.84242694e-02 | 1.49536134e-04 | 3.55271368e-15 | exploratory_double |
| signed_xxy | 49 | 98 | neumann_corrected_tail | square | 1.96532209e-03 | 1.96532209e-03 | 2.77938515e-03 | 1.55061444e-03 | 7.10542736e-15 | exploratory_double |
| symmetry_x_even | 49 | 101 | neumann_corrected_tail | x_only | 6.06214704e-02 | 0.00000000e+00 | 6.06214704e-02 | 3.55271368e-15 | 1.39009613e-02 | exploratory_double |
| symmetry_x_even | 49 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 4.70265782e-05 | 4.70265782e-05 | 1.56973714e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 49 | 101 | neumann_corrected_tail | omega_up | 4.70244977e-04 | 4.70244977e-05 | 4.72590353e-04 | 1.56880722e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 49 | 101 | neumann_corrected_tail | square | 4.70265574e-05 | 4.70265574e-05 | 6.65055953e-05 | 1.56972784e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 49 | 101 | neumann_corrected_tail | x_only | 6.06214704e-02 | 0.00000000e+00 | 6.06214704e-02 | 3.55271368e-15 | 1.39009613e-02 | exploratory_double |
| signed_xxy | 49 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 1.96685977e-03 | 1.96685977e-03 | 1.56690501e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 49 | 101 | neumann_corrected_tail | omega_up | 1.83328332e-02 | 1.83328332e-03 | 1.84242694e-02 | 1.49536134e-04 | 3.55271368e-15 | exploratory_double |
| signed_xxy | 49 | 101 | neumann_corrected_tail | square | 1.96532209e-03 | 1.96532209e-03 | 2.77938515e-03 | 1.55061444e-03 | 7.10542736e-15 | exploratory_double |
| symmetry_x_even | 49 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.54699444e-01 | 6.23079116e-01 | exploratory_double |
| symmetry_x_even | 49 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 9.20519696e-04 | 9.20519696e-04 | 1.56973763e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 49 | 200 | neumann_corrected_tail | omega_up | 9.19710759e-03 | 9.19710759e-04 | 9.24297874e-03 | 1.21347726e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 49 | 200 | neumann_corrected_tail | square | 9.20511591e-04 | 9.20511591e-04 | 1.30179998e-03 | 1.56616891e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 49 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.54699444e-01 | 6.23079116e-01 | exploratory_double |
| signed_xxy | 49 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 3.64274829e-02 | 3.64274829e-02 | 3.55271368e-15 | 3.90880446e-04 | exploratory_double |
| signed_xxy | 49 | 200 | neumann_corrected_tail | omega_up | 1.82654449e-01 | 1.82654449e-02 | 1.83565449e-01 | 1.77635684e-15 | 1.45247079e-01 | exploratory_double |
| signed_xxy | 49 | 200 | neumann_corrected_tail | square | 3.56318568e-02 | 3.56318568e-02 | 5.03910551e-02 | 1.77635684e-15 | 5.98059527e-03 | exploratory_double |
| norm_bound | 51 |  | tail_ignored | x_only | 4.43206314e-03 | 0.00000000e+00 | 4.43206314e-03 | 1.48707671e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 |  | tail_ignored | y_only | 0.00000000e+00 | 3.37155785e-03 | 3.37155785e-03 | 1.55821206e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 |  | tail_ignored | omega_up | 3.91721767e-03 | 3.91721767e-04 | 3.93675504e-03 | 1.50542585e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 |  | tail_ignored | square | 1.91496846e-03 | 1.91496846e-03 | 2.70817437e-03 | 1.55164729e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 |  | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.05708944e+01 | -5.05724641e+01 | exploratory_double |
| norm_bound | 51 |  | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.05708944e+01 | -5.05724641e+01 | exploratory_double |
| norm_bound | 51 |  | dirichlet_parseval_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.05708944e+01 | -5.05724641e+01 | exploratory_double |
| norm_bound | 51 |  | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.05708944e+01 | -5.05724641e+01 | exploratory_double |
| norm_bound | 51 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| norm_bound | 51 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| norm_bound | 51 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| norm_bound | 51 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| norm_bound | 51 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| norm_bound | 51 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| norm_bound | 51 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| norm_bound | 51 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| norm_bound | 51 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.13404283e+00 | -1.13561250e+00 | exploratory_double |
| norm_bound | 51 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.13404283e+00 | -1.13561250e+00 | exploratory_double |
| norm_bound | 51 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.13404283e+00 | -1.13561250e+00 | exploratory_double |
| norm_bound | 51 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.13404283e+00 | -1.13561250e+00 | exploratory_double |
| norm_bound | 51 | 98 | neumann_corrected_tail | x_only | 8.08789822e-05 | 0.00000000e+00 | 8.08789822e-05 | 1.56964881e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 5.19418661e-05 | 5.19418661e-05 | 1.56974314e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 | 98 | neumann_corrected_tail | omega_up | 6.99820610e-05 | 6.99820610e-06 | 7.03311008e-05 | 1.56966517e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 | 98 | neumann_corrected_tail | square | 3.16291313e-05 | 3.16291313e-05 | 4.47303464e-05 | 1.56971370e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 | 101 | neumann_corrected_tail | x_only | 8.08789822e-05 | 0.00000000e+00 | 8.08789822e-05 | 1.56964881e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 5.19418661e-05 | 5.19418661e-05 | 1.56974314e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 | 101 | neumann_corrected_tail | omega_up | 6.99820610e-05 | 6.99820610e-06 | 7.03311008e-05 | 1.56966517e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 | 101 | neumann_corrected_tail | square | 3.16291313e-05 | 3.16291313e-05 | 4.47303464e-05 | 1.56971370e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 | 200 | neumann_corrected_tail | x_only | 1.45309638e-03 | 0.00000000e+00 | 1.45309638e-03 | 1.56079755e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 9.24478503e-04 | 9.24478503e-04 | 1.56973274e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 | 200 | neumann_corrected_tail | omega_up | 1.25572928e-03 | 1.25572928e-04 | 1.26199231e-03 | 1.56319282e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 | 200 | neumann_corrected_tail | square | 5.65017932e-04 | 5.65017932e-04 | 7.99056023e-04 | 1.56865265e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 51 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.37096364e+00 | 2.53934331e+00 | exploratory_double |
| symmetry_x_even | 51 |  | tail_ignored | y_only | 0.00000000e+00 | 3.37155785e-03 | 3.37155785e-03 | 1.55821206e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 51 |  | tail_ignored | omega_up | 3.35420084e-02 | 3.35420084e-03 | 3.37093013e-02 | 1.77635684e-15 | 3.20177319e-03 | exploratory_double |
| symmetry_x_even | 51 |  | tail_ignored | square | 3.37141974e-03 | 3.37141974e-03 | 4.76790752e-03 | 1.51013780e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 51 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.37096364e+00 | 2.53934331e+00 | exploratory_double |
| signed_xxy | 51 |  | tail_ignored | y_only | 0.00000000e+00 | 9.79798235e-02 | 9.79798235e-02 | 5.32907052e-15 | 1.49075816e-02 | exploratory_double |
| signed_xxy | 51 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.72821471e+00 | 1.90355037e+00 | exploratory_double |
| signed_xxy | 51 |  | tail_ignored | square | 9.16946123e-02 | 9.16946123e-02 | 1.29675764e-01 | 1.77635684e-15 | 5.41927448e-02 | exploratory_double |
| symmetry_x_even | 51 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| symmetry_x_even | 51 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| symmetry_x_even | 51 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| symmetry_x_even | 51 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| signed_xxy | 51 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| signed_xxy | 51 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| signed_xxy | 51 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| signed_xxy | 51 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| symmetry_x_even | 51 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| symmetry_x_even | 51 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| symmetry_x_even | 51 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| symmetry_x_even | 51 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| signed_xxy | 51 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| signed_xxy | 51 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| signed_xxy | 51 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| signed_xxy | 51 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.22153236e+00 | -2.22310204e+00 | exploratory_double |
| symmetry_x_even | 51 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.13404283e+00 | -1.13561250e+00 | exploratory_double |
| symmetry_x_even | 51 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.13404283e+00 | -1.13561250e+00 | exploratory_double |
| symmetry_x_even | 51 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.13404283e+00 | -1.13561250e+00 | exploratory_double |
| symmetry_x_even | 51 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.13404283e+00 | -1.13561250e+00 | exploratory_double |
| signed_xxy | 51 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.13404283e+00 | -1.13561250e+00 | exploratory_double |
| signed_xxy | 51 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.13404283e+00 | -1.13561250e+00 | exploratory_double |
| signed_xxy | 51 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.13404283e+00 | -1.13561250e+00 | exploratory_double |
| signed_xxy | 51 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.13404283e+00 | -1.13561250e+00 | exploratory_double |
| symmetry_x_even | 51 | 98 | neumann_corrected_tail | x_only | 6.36190005e-02 | 0.00000000e+00 | 6.36190005e-02 | 1.77635684e-15 | 1.54706492e-02 | exploratory_double |
| symmetry_x_even | 51 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 5.19418661e-05 | 5.19418661e-05 | 1.56974314e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 51 | 98 | neumann_corrected_tail | omega_up | 5.19393259e-04 | 5.19393259e-05 | 5.21983766e-04 | 1.56860866e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 51 | 98 | neumann_corrected_tail | square | 5.19418407e-05 | 5.19418407e-05 | 7.34568556e-05 | 1.56973180e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 51 | 98 | neumann_corrected_tail | x_only | 6.36190005e-02 | 0.00000000e+00 | 6.36190005e-02 | 1.77635684e-15 | 1.54706492e-02 | exploratory_double |
| signed_xxy | 51 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 2.17025651e-03 | 2.17025651e-03 | 1.56599426e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 51 | 98 | neumann_corrected_tail | omega_up | 2.00532224e-02 | 2.00532224e-03 | 2.01532391e-02 | 1.77635684e-15 | 1.30306295e-04 | exploratory_double |
| signed_xxy | 51 | 98 | neumann_corrected_tail | square | 2.16838284e-03 | 2.16838284e-03 | 3.06655643e-03 | 1.54615755e-03 | 5.32907052e-15 | exploratory_double |
| symmetry_x_even | 51 | 101 | neumann_corrected_tail | x_only | 6.36190005e-02 | 0.00000000e+00 | 6.36190005e-02 | 1.77635684e-15 | 1.54706492e-02 | exploratory_double |
| symmetry_x_even | 51 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 5.19418661e-05 | 5.19418661e-05 | 1.56974314e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 51 | 101 | neumann_corrected_tail | omega_up | 5.19393259e-04 | 5.19393259e-05 | 5.21983766e-04 | 1.56860866e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 51 | 101 | neumann_corrected_tail | square | 5.19418407e-05 | 5.19418407e-05 | 7.34568556e-05 | 1.56973180e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 51 | 101 | neumann_corrected_tail | x_only | 6.36190005e-02 | 0.00000000e+00 | 6.36190005e-02 | 1.77635684e-15 | 1.54706492e-02 | exploratory_double |
| signed_xxy | 51 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 2.17025651e-03 | 2.17025651e-03 | 1.56599426e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 51 | 101 | neumann_corrected_tail | omega_up | 2.00532224e-02 | 2.00532224e-03 | 2.01532391e-02 | 1.77635684e-15 | 1.30306295e-04 | exploratory_double |
| signed_xxy | 51 | 101 | neumann_corrected_tail | square | 2.16838284e-03 | 2.16838284e-03 | 3.06655643e-03 | 1.54615755e-03 | 5.32907052e-15 | exploratory_double |
| symmetry_x_even | 51 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.58255398e-01 | 6.26635070e-01 | exploratory_double |
| symmetry_x_even | 51 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 9.24478503e-04 | 9.24478503e-04 | 1.56973274e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 51 | 200 | neumann_corrected_tail | omega_up | 9.23662002e-03 | 9.23662002e-04 | 9.28268823e-03 | 1.21040216e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 51 | 200 | neumann_corrected_tail | square | 9.24470322e-04 | 9.24470322e-04 | 1.30739847e-03 | 1.56613323e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 51 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.58255398e-01 | 6.26635070e-01 | exploratory_double |
| signed_xxy | 51 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 3.65518628e-02 | 3.65518628e-02 | 3.55271368e-15 | 4.05000959e-04 | exploratory_double |
| signed_xxy | 51 | 200 | neumann_corrected_tail | omega_up | 1.83085591e-01 | 1.83085591e-02 | 1.83998742e-01 | 1.77635684e-15 | 1.45959090e-01 | exploratory_double |
| signed_xxy | 51 | 200 | neumann_corrected_tail | square | 3.57506145e-02 | 3.57506145e-02 | 5.05590039e-02 | 5.32907052e-15 | 6.03284268e-03 | exploratory_double |
| norm_bound | 99 |  | tail_ignored | x_only | 4.39771885e-03 | 0.00000000e+00 | 4.39771885e-03 | 1.48835190e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 |  | tail_ignored | y_only | 0.00000000e+00 | 3.34571527e-03 | 3.34571527e-03 | 1.55842286e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 |  | tail_ignored | omega_up | 3.88690061e-03 | 3.88690061e-04 | 3.90628677e-03 | 1.50642096e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 |  | tail_ignored | square | 1.90022006e-03 | 1.90022006e-03 | 2.68731698e-03 | 1.55194437e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 |  | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98508756e+01 | -4.98524453e+01 | exploratory_double |
| norm_bound | 99 |  | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98508756e+01 | -4.98524453e+01 | exploratory_double |
| norm_bound | 99 |  | dirichlet_parseval_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98508756e+01 | -4.98524453e+01 | exploratory_double |
| norm_bound | 99 |  | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98508756e+01 | -4.98524453e+01 | exploratory_double |
| norm_bound | 99 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| norm_bound | 99 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| norm_bound | 99 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| norm_bound | 99 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| norm_bound | 99 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| norm_bound | 99 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| norm_bound | 99 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| norm_bound | 99 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| norm_bound | 99 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07978808e+00 | -1.08135776e+00 | exploratory_double |
| norm_bound | 99 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07978808e+00 | -1.08135776e+00 | exploratory_double |
| norm_bound | 99 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07978808e+00 | -1.08135776e+00 | exploratory_double |
| norm_bound | 99 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07978808e+00 | -1.08135776e+00 | exploratory_double |
| norm_bound | 99 | 98 | neumann_corrected_tail | x_only | 1.35907469e-04 | 0.00000000e+00 | 1.35907469e-04 | 1.56959865e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 8.75999944e-05 | 8.75999944e-05 | 1.56978464e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 | 98 | neumann_corrected_tail | omega_up | 1.17654045e-04 | 1.17654045e-05 | 1.18240852e-04 | 1.56963391e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 | 98 | neumann_corrected_tail | square | 5.32667083e-05 | 5.32667083e-05 | 7.53305014e-05 | 1.56973281e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 | 101 | neumann_corrected_tail | x_only | 1.35907469e-04 | 0.00000000e+00 | 1.35907469e-04 | 1.56959865e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 8.75999944e-05 | 8.75999944e-05 | 1.56978464e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 | 101 | neumann_corrected_tail | omega_up | 1.17654045e-04 | 1.17654045e-05 | 1.18240852e-04 | 1.56963391e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 | 101 | neumann_corrected_tail | square | 5.32667083e-05 | 5.32667083e-05 | 7.53305014e-05 | 1.56973281e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 | 200 | neumann_corrected_tail | x_only | 1.49230235e-03 | 0.00000000e+00 | 1.49230235e-03 | 1.56031197e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 9.53138846e-04 | 9.53138846e-04 | 1.56969605e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 | 200 | neumann_corrected_tail | omega_up | 1.29029380e-03 | 1.29029380e-04 | 1.29672922e-03 | 1.56282612e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 | 200 | neumann_corrected_tail | square | 5.81649209e-04 | 5.81649209e-04 | 8.22576201e-04 | 1.56856817e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 99 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34981648e+00 | 2.51819616e+00 | exploratory_double |
| symmetry_x_even | 99 |  | tail_ignored | y_only | 0.00000000e+00 | 3.34571527e-03 | 3.34571527e-03 | 1.55842286e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 99 |  | tail_ignored | omega_up | 3.32865834e-02 | 3.32865834e-03 | 3.34526023e-02 | 1.77635684e-15 | 3.12911161e-03 | exploratory_double |
| symmetry_x_even | 99 |  | tail_ignored | square | 3.34557935e-03 | 3.34557935e-03 | 4.73136369e-03 | 1.51108481e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 99 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34981648e+00 | 2.51819616e+00 | exploratory_double |
| signed_xxy | 99 |  | tail_ignored | y_only | 0.00000000e+00 | 9.73544387e-02 | 9.73544387e-02 | 1.77635684e-15 | 1.46747705e-02 | exploratory_double |
| signed_xxy | 99 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.70755593e+00 | 1.88289159e+00 | exploratory_double |
| signed_xxy | 99 |  | tail_ignored | square | 9.11533838e-02 | 9.11533838e-02 | 1.28910352e-01 | 5.32907052e-15 | 5.34785288e-02 | exploratory_double |
| symmetry_x_even | 99 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| symmetry_x_even | 99 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| symmetry_x_even | 99 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| symmetry_x_even | 99 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| signed_xxy | 99 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| signed_xxy | 99 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| signed_xxy | 99 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| signed_xxy | 99 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| symmetry_x_even | 99 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| symmetry_x_even | 99 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| symmetry_x_even | 99 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| symmetry_x_even | 99 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| signed_xxy | 99 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| signed_xxy | 99 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| signed_xxy | 99 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| signed_xxy | 99 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15263244e+00 | -2.15420212e+00 | exploratory_double |
| symmetry_x_even | 99 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07978808e+00 | -1.08135776e+00 | exploratory_double |
| symmetry_x_even | 99 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07978808e+00 | -1.08135776e+00 | exploratory_double |
| symmetry_x_even | 99 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07978808e+00 | -1.08135776e+00 | exploratory_double |
| symmetry_x_even | 99 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07978808e+00 | -1.08135776e+00 | exploratory_double |
| signed_xxy | 99 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07978808e+00 | -1.08135776e+00 | exploratory_double |
| signed_xxy | 99 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07978808e+00 | -1.08135776e+00 | exploratory_double |
| signed_xxy | 99 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07978808e+00 | -1.08135776e+00 | exploratory_double |
| signed_xxy | 99 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07978808e+00 | -1.08135776e+00 | exploratory_double |
| symmetry_x_even | 99 | 98 | neumann_corrected_tail | x_only | 8.21423323e-02 | 0.00000000e+00 | 8.21423323e-02 | 3.55271368e-15 | 2.68610482e-02 | exploratory_double |
| symmetry_x_even | 99 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 8.75999944e-05 | 8.75999944e-05 | 1.56978464e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 99 | 98 | neumann_corrected_tail | omega_up | 8.75927290e-04 | 8.75927290e-05 | 8.80296032e-04 | 1.56655788e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 99 | 98 | neumann_corrected_tail | square | 8.75999217e-05 | 8.75999217e-05 | 1.23884997e-04 | 1.56975237e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 99 | 98 | neumann_corrected_tail | x_only | 8.21423323e-02 | 0.00000000e+00 | 8.21423323e-02 | 3.55271368e-15 | 2.68610482e-02 | exploratory_double |
| signed_xxy | 99 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 3.63365189e-03 | 3.63365189e-03 | 1.55596605e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 99 | 98 | neumann_corrected_tail | omega_up | 3.13998063e-02 | 3.13998063e-03 | 3.15564148e-02 | 5.32907052e-15 | 2.60975395e-03 | exploratory_double |
| signed_xxy | 99 | 98 | neumann_corrected_tail | square | 3.62836905e-03 | 3.62836905e-03 | 5.13128872e-03 | 1.50030582e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 99 | 101 | neumann_corrected_tail | x_only | 8.21423323e-02 | 0.00000000e+00 | 8.21423323e-02 | 3.55271368e-15 | 2.68610482e-02 | exploratory_double |
| symmetry_x_even | 99 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 8.75999944e-05 | 8.75999944e-05 | 1.56978464e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 99 | 101 | neumann_corrected_tail | omega_up | 8.75927290e-04 | 8.75927290e-05 | 8.80296032e-04 | 1.56655788e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 99 | 101 | neumann_corrected_tail | square | 8.75999217e-05 | 8.75999217e-05 | 1.23884997e-04 | 1.56975237e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 99 | 101 | neumann_corrected_tail | x_only | 8.21423323e-02 | 0.00000000e+00 | 8.21423323e-02 | 3.55271368e-15 | 2.68610482e-02 | exploratory_double |
| signed_xxy | 99 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 3.63365189e-03 | 3.63365189e-03 | 1.55596605e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 99 | 101 | neumann_corrected_tail | omega_up | 3.13998063e-02 | 3.13998063e-03 | 3.15564148e-02 | 5.32907052e-15 | 2.60975395e-03 | exploratory_double |
| signed_xxy | 99 | 101 | neumann_corrected_tail | square | 3.62836905e-03 | 3.62836905e-03 | 5.13128872e-03 | 1.50030582e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 99 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.83836314e-01 | 6.52215987e-01 | exploratory_double |
| symmetry_x_even | 99 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 9.53138846e-04 | 9.53138846e-04 | 1.56969605e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 99 | 200 | neumann_corrected_tail | omega_up | 9.52266373e-03 | 9.52266373e-04 | 9.57015861e-03 | 1.18774556e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 99 | 200 | neumann_corrected_tail | square | 9.53130104e-04 | 9.53130104e-04 | 1.34792952e-03 | 1.56586971e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 99 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.83836314e-01 | 6.52215987e-01 | exploratory_double |
| signed_xxy | 99 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 3.74453933e-02 | 3.74453933e-02 | 1.77635684e-15 | 5.08042377e-04 | exploratory_double |
| signed_xxy | 99 | 200 | neumann_corrected_tail | omega_up | 1.86168430e-01 | 1.86168430e-02 | 1.87096956e-01 | 8.88178420e-15 | 1.51103315e-01 | exploratory_double |
| signed_xxy | 99 | 200 | neumann_corrected_tail | square | 3.66031414e-02 | 3.66031414e-02 | 5.17646591e-02 | 3.55271368e-15 | 6.41354120e-03 | exploratory_double |
| norm_bound | 101 |  | tail_ignored | x_only | 4.39771885e-03 | 0.00000000e+00 | 4.39771885e-03 | 1.48835190e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 |  | tail_ignored | y_only | 0.00000000e+00 | 3.34571527e-03 | 3.34571527e-03 | 1.55842286e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 |  | tail_ignored | omega_up | 3.88690061e-03 | 3.88690061e-04 | 3.90628677e-03 | 1.50642096e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 |  | tail_ignored | square | 1.90022006e-03 | 1.90022006e-03 | 2.68731698e-03 | 1.55194437e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 |  | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98375883e+01 | -4.98391579e+01 | exploratory_double |
| norm_bound | 101 |  | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98375883e+01 | -4.98391579e+01 | exploratory_double |
| norm_bound | 101 |  | dirichlet_parseval_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98375883e+01 | -4.98391579e+01 | exploratory_double |
| norm_bound | 101 |  | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98375883e+01 | -4.98391579e+01 | exploratory_double |
| norm_bound | 101 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| norm_bound | 101 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| norm_bound | 101 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| norm_bound | 101 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| norm_bound | 101 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| norm_bound | 101 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| norm_bound | 101 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| norm_bound | 101 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| norm_bound | 101 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07878871e+00 | -1.08035839e+00 | exploratory_double |
| norm_bound | 101 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07878871e+00 | -1.08035839e+00 | exploratory_double |
| norm_bound | 101 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07878871e+00 | -1.08035839e+00 | exploratory_double |
| norm_bound | 101 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07878871e+00 | -1.08035839e+00 | exploratory_double |
| norm_bound | 101 | 98 | neumann_corrected_tail | x_only | 1.36922004e-04 | 0.00000000e+00 | 1.36922004e-04 | 1.56959748e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 8.82598613e-05 | 8.82598613e-05 | 1.56978537e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 | 98 | neumann_corrected_tail | omega_up | 1.18533392e-04 | 1.18533392e-05 | 1.19124585e-04 | 1.56963315e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 | 98 | neumann_corrected_tail | square | 5.36665362e-05 | 5.36665362e-05 | 7.58959433e-05 | 1.56973312e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 | 101 | neumann_corrected_tail | x_only | 1.36922004e-04 | 0.00000000e+00 | 1.36922004e-04 | 1.56959748e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 8.82598613e-05 | 8.82598613e-05 | 1.56978537e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 | 101 | neumann_corrected_tail | omega_up | 1.18533392e-04 | 1.18533392e-05 | 1.19124585e-04 | 1.56963315e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 | 101 | neumann_corrected_tail | square | 5.36665362e-05 | 5.36665362e-05 | 7.58959433e-05 | 1.56973312e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 | 200 | neumann_corrected_tail | x_only | 1.49302354e-03 | 0.00000000e+00 | 1.49302354e-03 | 1.56030292e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 9.53668443e-04 | 9.53668443e-04 | 1.56969535e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 | 200 | neumann_corrected_tail | omega_up | 1.29093000e-03 | 1.29093000e-04 | 1.29736860e-03 | 1.56281927e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 | 200 | neumann_corrected_tail | square | 5.81955993e-04 | 5.81955993e-04 | 8.23010057e-04 | 1.56856658e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 101 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34981648e+00 | 2.51819616e+00 | exploratory_double |
| symmetry_x_even | 101 |  | tail_ignored | y_only | 0.00000000e+00 | 3.34571527e-03 | 3.34571527e-03 | 1.55842286e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 101 |  | tail_ignored | omega_up | 3.32865834e-02 | 3.32865834e-03 | 3.34526023e-02 | 1.77635684e-15 | 3.12911161e-03 | exploratory_double |
| symmetry_x_even | 101 |  | tail_ignored | square | 3.34557935e-03 | 3.34557935e-03 | 4.73136369e-03 | 1.51108481e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 101 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34981648e+00 | 2.51819616e+00 | exploratory_double |
| signed_xxy | 101 |  | tail_ignored | y_only | 0.00000000e+00 | 9.73544387e-02 | 9.73544387e-02 | 1.77635684e-15 | 1.46747705e-02 | exploratory_double |
| signed_xxy | 101 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.70755593e+00 | 1.88289159e+00 | exploratory_double |
| signed_xxy | 101 |  | tail_ignored | square | 9.11533838e-02 | 9.11533838e-02 | 1.28910352e-01 | 5.32907052e-15 | 5.34785288e-02 | exploratory_double |
| symmetry_x_even | 101 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| symmetry_x_even | 101 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| symmetry_x_even | 101 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| symmetry_x_even | 101 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| signed_xxy | 101 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| signed_xxy | 101 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| signed_xxy | 101 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| signed_xxy | 101 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| symmetry_x_even | 101 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| symmetry_x_even | 101 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| symmetry_x_even | 101 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| symmetry_x_even | 101 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| signed_xxy | 101 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| signed_xxy | 101 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| signed_xxy | 101 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| signed_xxy | 101 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.15136276e+00 | -2.15293244e+00 | exploratory_double |
| symmetry_x_even | 101 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07878871e+00 | -1.08035839e+00 | exploratory_double |
| symmetry_x_even | 101 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07878871e+00 | -1.08035839e+00 | exploratory_double |
| symmetry_x_even | 101 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07878871e+00 | -1.08035839e+00 | exploratory_double |
| symmetry_x_even | 101 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07878871e+00 | -1.08035839e+00 | exploratory_double |
| signed_xxy | 101 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07878871e+00 | -1.08035839e+00 | exploratory_double |
| signed_xxy | 101 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07878871e+00 | -1.08035839e+00 | exploratory_double |
| signed_xxy | 101 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07878871e+00 | -1.08035839e+00 | exploratory_double |
| signed_xxy | 101 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.07878871e+00 | -1.08035839e+00 | exploratory_double |
| symmetry_x_even | 101 | 98 | neumann_corrected_tail | x_only | 8.24457411e-02 | 0.00000000e+00 | 8.24457411e-02 | 1.77635684e-15 | 2.70718860e-02 | exploratory_double |
| symmetry_x_even | 101 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 8.82598613e-05 | 8.82598613e-05 | 1.56978537e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 101 | 98 | neumann_corrected_tail | omega_up | 8.82524854e-04 | 8.82524854e-05 | 8.86926502e-04 | 1.56650982e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 101 | 98 | neumann_corrected_tail | square | 8.82597876e-05 | 8.82597876e-05 | 1.24818189e-04 | 1.56975261e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 101 | 98 | neumann_corrected_tail | x_only | 8.24457411e-02 | 0.00000000e+00 | 8.24457411e-02 | 1.77635684e-15 | 2.70718860e-02 | exploratory_double |
| signed_xxy | 101 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 3.66053325e-03 | 3.66053325e-03 | 1.55572456e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 101 | 98 | neumann_corrected_tail | omega_up | 3.15981647e-02 | 3.15981647e-03 | 3.17557625e-02 | 3.55271368e-15 | 2.66291491e-03 | exploratory_double |
| signed_xxy | 101 | 98 | neumann_corrected_tail | square | 3.65517139e-03 | 3.65517139e-03 | 5.16919295e-03 | 1.49923677e-03 | 3.55271368e-15 | exploratory_double |
| symmetry_x_even | 101 | 101 | neumann_corrected_tail | x_only | 8.24457411e-02 | 0.00000000e+00 | 8.24457411e-02 | 1.77635684e-15 | 2.70718860e-02 | exploratory_double |
| symmetry_x_even | 101 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 8.82598613e-05 | 8.82598613e-05 | 1.56978537e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 101 | 101 | neumann_corrected_tail | omega_up | 8.82524854e-04 | 8.82524854e-05 | 8.86926502e-04 | 1.56650982e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 101 | 101 | neumann_corrected_tail | square | 8.82597876e-05 | 8.82597876e-05 | 1.24818189e-04 | 1.56975261e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 101 | 101 | neumann_corrected_tail | x_only | 8.24457411e-02 | 0.00000000e+00 | 8.24457411e-02 | 1.77635684e-15 | 2.70718860e-02 | exploratory_double |
| signed_xxy | 101 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 3.66053325e-03 | 3.66053325e-03 | 1.55572456e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 101 | 101 | neumann_corrected_tail | omega_up | 3.15981647e-02 | 3.15981647e-03 | 3.17557625e-02 | 3.55271368e-15 | 2.66291491e-03 | exploratory_double |
| signed_xxy | 101 | 101 | neumann_corrected_tail | square | 3.65517139e-03 | 3.65517139e-03 | 5.16919295e-03 | 1.49923677e-03 | 3.55271368e-15 | exploratory_double |
| symmetry_x_even | 101 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.84306453e-01 | 6.52686125e-01 | exploratory_double |
| symmetry_x_even | 101 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 9.53668443e-04 | 9.53668443e-04 | 1.56969535e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 101 | 200 | neumann_corrected_tail | omega_up | 9.52794916e-03 | 9.52794916e-04 | 9.57547040e-03 | 1.18732039e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 101 | 200 | neumann_corrected_tail | square | 9.53659690e-04 | 9.53659690e-04 | 1.34867847e-03 | 1.56586476e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 101 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.84306453e-01 | 6.52686125e-01 | exploratory_double |
| signed_xxy | 101 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 3.74617947e-02 | 3.74617947e-02 | 1.77635684e-15 | 5.09960113e-04 | exploratory_double |
| signed_xxy | 101 | 200 | neumann_corrected_tail | omega_up | 1.86224784e-01 | 1.86224784e-02 | 1.87153592e-01 | 1.77635684e-15 | 1.51198220e-01 | exploratory_double |
| signed_xxy | 101 | 200 | neumann_corrected_tail | square | 3.66187801e-02 | 3.66187801e-02 | 5.17867755e-02 | 3.55271368e-15 | 6.42061717e-03 | exploratory_double |
| norm_bound | 199 |  | tail_ignored | x_only | 4.38339766e-03 | 0.00000000e+00 | 4.38339766e-03 | 1.48888070e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 |  | tail_ignored | y_only | 0.00000000e+00 | 3.33493667e-03 | 3.33493667e-03 | 1.55851022e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 |  | tail_ignored | omega_up | 3.87425839e-03 | 3.87425839e-04 | 3.89358149e-03 | 1.50683362e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 |  | tail_ignored | square | 1.89406933e-03 | 1.89406933e-03 | 2.67861853e-03 | 1.55206753e-03 | 3.55271368e-15 | exploratory_double |
| norm_bound | 199 |  | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94451199e+01 | -4.94466896e+01 | exploratory_double |
| norm_bound | 199 |  | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94451199e+01 | -4.94466896e+01 | exploratory_double |
| norm_bound | 199 |  | dirichlet_parseval_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94451199e+01 | -4.94466896e+01 | exploratory_double |
| norm_bound | 199 |  | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94451199e+01 | -4.94466896e+01 | exploratory_double |
| norm_bound | 199 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| norm_bound | 199 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| norm_bound | 199 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| norm_bound | 199 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| norm_bound | 199 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| norm_bound | 199 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| norm_bound | 199 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| norm_bound | 199 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| norm_bound | 199 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04931005e+00 | -1.05087972e+00 | exploratory_double |
| norm_bound | 199 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04931005e+00 | -1.05087972e+00 | exploratory_double |
| norm_bound | 199 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04931005e+00 | -1.05087972e+00 | exploratory_double |
| norm_bound | 199 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04931005e+00 | -1.05087972e+00 | exploratory_double |
| norm_bound | 199 | 98 | neumann_corrected_tail | x_only | 1.66858018e-04 | 0.00000000e+00 | 1.66858018e-04 | 1.56955924e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 1.07771165e-04 | 1.07771165e-04 | 1.56980653e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 | 98 | neumann_corrected_tail | omega_up | 1.44487653e-04 | 1.44487653e-05 | 1.45208294e-04 | 1.56960788e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 | 98 | neumann_corrected_tail | square | 6.54792313e-05 | 6.54792313e-05 | 9.26016170e-05 | 1.56974128e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 | 101 | neumann_corrected_tail | x_only | 1.66858018e-04 | 0.00000000e+00 | 1.66858018e-04 | 1.56955924e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 1.07771165e-04 | 1.07771165e-04 | 1.56980653e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 | 101 | neumann_corrected_tail | omega_up | 1.44487653e-04 | 1.44487653e-05 | 1.45208294e-04 | 1.56960788e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 | 101 | neumann_corrected_tail | square | 6.54792313e-05 | 6.54792313e-05 | 9.26016170e-05 | 1.56974128e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 | 200 | neumann_corrected_tail | x_only | 1.51427190e-03 | 0.00000000e+00 | 1.51427190e-03 | 1.56003422e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 9.69312971e-04 | 9.69312971e-04 | 1.56967432e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 | 200 | neumann_corrected_tail | omega_up | 1.30968133e-03 | 1.30968133e-04 | 1.31621344e-03 | 1.56261602e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 | 200 | neumann_corrected_tail | square | 5.91009354e-04 | 5.91009354e-04 | 8.35813443e-04 | 1.56851925e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 199 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34098483e+00 | 2.50936450e+00 | exploratory_double |
| symmetry_x_even | 199 |  | tail_ignored | y_only | 0.00000000e+00 | 3.33493667e-03 | 3.33493667e-03 | 1.55851022e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 199 |  | tail_ignored | omega_up | 3.31800413e-02 | 3.31800413e-03 | 3.33455288e-02 | 1.77635684e-15 | 3.09896911e-03 | exploratory_double |
| symmetry_x_even | 199 |  | tail_ignored | square | 3.33480166e-03 | 3.33480166e-03 | 4.71612173e-03 | 1.51147754e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 199 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34098483e+00 | 2.50936450e+00 | exploratory_double |
| signed_xxy | 199 |  | tail_ignored | y_only | 0.00000000e+00 | 9.70929393e-02 | 9.70929393e-02 | 1.77635684e-15 | 1.45780031e-02 | exploratory_double |
| signed_xxy | 199 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.69892823e+00 | 1.87426389e+00 | exploratory_double |
| signed_xxy | 199 |  | tail_ignored | square | 9.09268969e-02 | 9.09268969e-02 | 1.28590051e-01 | 3.55271368e-15 | 5.31812372e-02 | exploratory_double |
| symmetry_x_even | 199 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| symmetry_x_even | 199 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| symmetry_x_even | 199 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| symmetry_x_even | 199 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| signed_xxy | 199 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| signed_xxy | 199 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| signed_xxy | 199 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| signed_xxy | 199 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| symmetry_x_even | 199 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| symmetry_x_even | 199 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| symmetry_x_even | 199 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| symmetry_x_even | 199 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| signed_xxy | 199 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| signed_xxy | 199 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| signed_xxy | 199 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| signed_xxy | 199 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11389922e+00 | -2.11546889e+00 | exploratory_double |
| symmetry_x_even | 199 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04931005e+00 | -1.05087972e+00 | exploratory_double |
| symmetry_x_even | 199 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04931005e+00 | -1.05087972e+00 | exploratory_double |
| symmetry_x_even | 199 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04931005e+00 | -1.05087972e+00 | exploratory_double |
| symmetry_x_even | 199 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04931005e+00 | -1.05087972e+00 | exploratory_double |
| signed_xxy | 199 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04931005e+00 | -1.05087972e+00 | exploratory_double |
| signed_xxy | 199 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04931005e+00 | -1.05087972e+00 | exploratory_double |
| signed_xxy | 199 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04931005e+00 | -1.05087972e+00 | exploratory_double |
| signed_xxy | 199 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04931005e+00 | -1.05087972e+00 | exploratory_double |
| symmetry_x_even | 199 | 98 | neumann_corrected_tail | x_only | 9.09585363e-02 | 0.00000000e+00 | 9.09585363e-02 | 1.77635684e-15 | 3.33069552e-02 | exploratory_double |
| symmetry_x_even | 199 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 1.07771165e-04 | 1.07771165e-04 | 1.56980653e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 199 | 98 | neumann_corrected_tail | omega_up | 1.07760134e-03 | 1.07760134e-04 | 1.08297594e-03 | 1.56492269e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 199 | 98 | neumann_corrected_tail | square | 1.07771055e-04 | 1.07771055e-04 | 1.52411287e-04 | 1.56975768e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 199 | 98 | neumann_corrected_tail | x_only | 9.09585363e-02 | 0.00000000e+00 | 9.09585363e-02 | 1.77635684e-15 | 3.33069552e-02 | exploratory_double |
| signed_xxy | 199 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 4.45213699e-03 | 4.45213699e-03 | 1.54768200e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 199 | 98 | neumann_corrected_tail | omega_up | 3.72986093e-02 | 3.72986093e-03 | 3.74846384e-02 | 1.77635684e-15 | 4.33531666e-03 | exploratory_double |
| signed_xxy | 199 | 98 | neumann_corrected_tail | square | 4.44418062e-03 | 4.44418062e-03 | 6.28502051e-03 | 1.46407926e-03 | 3.55271368e-15 | exploratory_double |
| symmetry_x_even | 199 | 101 | neumann_corrected_tail | x_only | 9.09585363e-02 | 0.00000000e+00 | 9.09585363e-02 | 1.77635684e-15 | 3.33069552e-02 | exploratory_double |
| symmetry_x_even | 199 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 1.07771165e-04 | 1.07771165e-04 | 1.56980653e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 199 | 101 | neumann_corrected_tail | omega_up | 1.07760134e-03 | 1.07760134e-04 | 1.08297594e-03 | 1.56492269e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 199 | 101 | neumann_corrected_tail | square | 1.07771055e-04 | 1.07771055e-04 | 1.52411287e-04 | 1.56975768e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 199 | 101 | neumann_corrected_tail | x_only | 9.09585363e-02 | 0.00000000e+00 | 9.09585363e-02 | 1.77635684e-15 | 3.33069552e-02 | exploratory_double |
| signed_xxy | 199 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 4.45213699e-03 | 4.45213699e-03 | 1.54768200e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 199 | 101 | neumann_corrected_tail | omega_up | 3.72986093e-02 | 3.72986093e-03 | 3.74846384e-02 | 1.77635684e-15 | 4.33531666e-03 | exploratory_double |
| signed_xxy | 199 | 101 | neumann_corrected_tail | square | 4.44418062e-03 | 4.44418062e-03 | 6.28502051e-03 | 1.46407926e-03 | 3.55271368e-15 | exploratory_double |
| symmetry_x_even | 199 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.98151493e-01 | 6.66531165e-01 | exploratory_double |
| symmetry_x_even | 199 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 9.69312971e-04 | 9.69312971e-04 | 1.56967432e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 199 | 200 | neumann_corrected_tail | omega_up | 9.68407969e-03 | 9.68407969e-04 | 9.73237964e-03 | 1.17465404e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 199 | 200 | neumann_corrected_tail | square | 9.69303902e-04 | 9.69303902e-04 | 1.37080272e-03 | 1.56571691e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 199 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.98151493e-01 | 6.66531165e-01 | exploratory_double |
| signed_xxy | 199 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 3.79444710e-02 | 3.79444710e-02 | 8.88178420e-15 | 5.66823123e-04 | exploratory_double |
| signed_xxy | 199 | 200 | neumann_corrected_tail | omega_up | 1.87879567e-01 | 1.87879567e-02 | 1.88816628e-01 | 1.77635684e-15 | 1.53998875e-01 | exploratory_double |
| signed_xxy | 199 | 200 | neumann_corrected_tail | square | 3.70788463e-02 | 3.70788463e-02 | 5.24374073e-02 | 1.77635684e-15 | 6.63027536e-03 | exploratory_double |
| norm_bound | 203 |  | tail_ignored | x_only | 4.38236413e-03 | 0.00000000e+00 | 4.38236413e-03 | 1.48891880e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 |  | tail_ignored | y_only | 0.00000000e+00 | 3.33415874e-03 | 3.33415874e-03 | 1.55851651e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 |  | tail_ignored | omega_up | 3.87334601e-03 | 3.87334601e-04 | 3.89266457e-03 | 1.50686335e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 |  | tail_ignored | square | 1.89362542e-03 | 1.89362542e-03 | 2.67799075e-03 | 1.55207641e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 |  | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94415658e+01 | -4.94431355e+01 | exploratory_double |
| norm_bound | 203 |  | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94415658e+01 | -4.94431355e+01 | exploratory_double |
| norm_bound | 203 |  | dirichlet_parseval_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94415658e+01 | -4.94431355e+01 | exploratory_double |
| norm_bound | 203 |  | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94415658e+01 | -4.94431355e+01 | exploratory_double |
| norm_bound | 203 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| norm_bound | 203 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| norm_bound | 203 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| norm_bound | 203 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| norm_bound | 203 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| norm_bound | 203 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| norm_bound | 203 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| norm_bound | 203 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| norm_bound | 203 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04904333e+00 | -1.05061300e+00 | exploratory_double |
| norm_bound | 203 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04904333e+00 | -1.05061300e+00 | exploratory_double |
| norm_bound | 203 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04904333e+00 | -1.05061300e+00 | exploratory_double |
| norm_bound | 203 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04904333e+00 | -1.05061300e+00 | exploratory_double |
| norm_bound | 203 | 98 | neumann_corrected_tail | x_only | 1.67129016e-04 | 0.00000000e+00 | 1.67129016e-04 | 1.56955886e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 1.07948151e-04 | 1.07948151e-04 | 1.56980672e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 | 98 | neumann_corrected_tail | omega_up | 1.44722670e-04 | 1.44722670e-05 | 1.45444483e-04 | 1.56960762e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 | 98 | neumann_corrected_tail | square | 6.55862983e-05 | 6.55862983e-05 | 9.27530326e-05 | 1.56974135e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 | 101 | neumann_corrected_tail | x_only | 1.67129016e-04 | 0.00000000e+00 | 1.67129016e-04 | 1.56955886e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 1.07948151e-04 | 1.07948151e-04 | 1.56980672e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 | 101 | neumann_corrected_tail | omega_up | 1.44722670e-04 | 1.44722670e-05 | 1.45444483e-04 | 1.56960762e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 | 101 | neumann_corrected_tail | square | 6.55862983e-05 | 6.55862983e-05 | 9.27530326e-05 | 1.56974135e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 | 200 | neumann_corrected_tail | x_only | 1.51446404e-03 | 0.00000000e+00 | 1.51446404e-03 | 1.56003177e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 9.69454778e-04 | 9.69454778e-04 | 1.56967413e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 | 200 | neumann_corrected_tail | omega_up | 1.30985094e-03 | 1.30985094e-04 | 1.31638390e-03 | 1.56261416e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 | 200 | neumann_corrected_tail | square | 5.91091339e-04 | 5.91091339e-04 | 8.35929389e-04 | 1.56851882e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 203 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34034709e+00 | 2.50872677e+00 | exploratory_double |
| symmetry_x_even | 203 |  | tail_ignored | y_only | 0.00000000e+00 | 3.33415874e-03 | 3.33415874e-03 | 1.55851651e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 203 |  | tail_ignored | omega_up | 3.31723515e-02 | 3.31723515e-03 | 3.33378007e-02 | 1.77635684e-15 | 3.09679732e-03 | exploratory_double |
| symmetry_x_even | 203 |  | tail_ignored | square | 3.33402379e-03 | 3.33402379e-03 | 4.71502166e-03 | 1.51150584e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 203 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34034709e+00 | 2.50872677e+00 | exploratory_double |
| signed_xxy | 203 |  | tail_ignored | y_only | 0.00000000e+00 | 9.70740492e-02 | 9.70740492e-02 | 7.10542736e-15 | 1.45710261e-02 | exploratory_double |
| signed_xxy | 203 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.69830523e+00 | 1.87364089e+00 | exploratory_double |
| signed_xxy | 203 |  | tail_ignored | square | 9.09105320e-02 | 9.09105320e-02 | 1.28566907e-01 | 3.55271368e-15 | 5.31597925e-02 | exploratory_double |
| symmetry_x_even | 203 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| symmetry_x_even | 203 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| symmetry_x_even | 203 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| symmetry_x_even | 203 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| signed_xxy | 203 | 29 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| signed_xxy | 203 | 29 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| signed_xxy | 203 | 29 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| signed_xxy | 203 | 29 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| symmetry_x_even | 203 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| symmetry_x_even | 203 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| symmetry_x_even | 203 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| symmetry_x_even | 203 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| signed_xxy | 203 | 31 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| signed_xxy | 203 | 31 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| signed_xxy | 203 | 31 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| signed_xxy | 203 | 31 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.11356018e+00 | -2.11512986e+00 | exploratory_double |
| symmetry_x_even | 203 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04904333e+00 | -1.05061300e+00 | exploratory_double |
| symmetry_x_even | 203 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04904333e+00 | -1.05061300e+00 | exploratory_double |
| symmetry_x_even | 203 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04904333e+00 | -1.05061300e+00 | exploratory_double |
| symmetry_x_even | 203 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04904333e+00 | -1.05061300e+00 | exploratory_double |
| signed_xxy | 203 | 50 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04904333e+00 | -1.05061300e+00 | exploratory_double |
| signed_xxy | 203 | 50 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04904333e+00 | -1.05061300e+00 | exploratory_double |
| signed_xxy | 203 | 50 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04904333e+00 | -1.05061300e+00 | exploratory_double |
| signed_xxy | 203 | 50 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.04904333e+00 | -1.05061300e+00 | exploratory_double |
| symmetry_x_even | 203 | 98 | neumann_corrected_tail | x_only | 9.10320924e-02 | 0.00000000e+00 | 9.10320924e-02 | 1.77635684e-15 | 3.33635216e-02 | exploratory_double |
| symmetry_x_even | 203 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 1.07948151e-04 | 1.07948151e-04 | 1.56980672e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 203 | 98 | neumann_corrected_tail | omega_up | 1.07937083e-03 | 1.07937083e-04 | 1.08475426e-03 | 1.56490682e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 203 | 98 | neumann_corrected_tail | square | 1.07948040e-04 | 1.07948040e-04 | 1.52661582e-04 | 1.56975771e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 203 | 98 | neumann_corrected_tail | x_only | 9.10320924e-02 | 0.00000000e+00 | 9.10320924e-02 | 1.77635684e-15 | 3.33635216e-02 | exploratory_double |
| signed_xxy | 203 | 98 | neumann_corrected_tail | y_only | 0.00000000e+00 | 4.45928907e-03 | 4.45928907e-03 | 1.54760112e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 203 | 98 | neumann_corrected_tail | omega_up | 3.73489296e-02 | 3.73489296e-03 | 3.75352097e-02 | 1.77635684e-15 | 4.35132619e-03 | exploratory_double |
| signed_xxy | 203 | 98 | neumann_corrected_tail | square | 4.45130689e-03 | 4.45130689e-03 | 6.29509858e-03 | 1.46372918e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 203 | 101 | neumann_corrected_tail | x_only | 9.10320924e-02 | 0.00000000e+00 | 9.10320924e-02 | 1.77635684e-15 | 3.33635216e-02 | exploratory_double |
| symmetry_x_even | 203 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 1.07948151e-04 | 1.07948151e-04 | 1.56980672e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 203 | 101 | neumann_corrected_tail | omega_up | 1.07937083e-03 | 1.07937083e-04 | 1.08475426e-03 | 1.56490682e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 203 | 101 | neumann_corrected_tail | square | 1.07948040e-04 | 1.07948040e-04 | 1.52661582e-04 | 1.56975771e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 203 | 101 | neumann_corrected_tail | x_only | 9.10320924e-02 | 0.00000000e+00 | 9.10320924e-02 | 1.77635684e-15 | 3.33635216e-02 | exploratory_double |
| signed_xxy | 203 | 101 | neumann_corrected_tail | y_only | 0.00000000e+00 | 4.45928907e-03 | 4.45928907e-03 | 1.54760112e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 203 | 101 | neumann_corrected_tail | omega_up | 3.73489296e-02 | 3.73489296e-03 | 3.75352097e-02 | 1.77635684e-15 | 4.35132619e-03 | exploratory_double |
| signed_xxy | 203 | 101 | neumann_corrected_tail | square | 4.45130689e-03 | 4.45130689e-03 | 6.29509858e-03 | 1.46372918e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 203 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.98276627e-01 | 6.66656299e-01 | exploratory_double |
| symmetry_x_even | 203 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 9.69454778e-04 | 9.69454778e-04 | 1.56967413e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 203 | 200 | neumann_corrected_tail | omega_up | 9.68549488e-03 | 9.68549488e-04 | 9.73380188e-03 | 1.17453829e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 203 | 200 | neumann_corrected_tail | square | 9.69445706e-04 | 9.69445706e-04 | 1.37100327e-03 | 1.56571556e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 203 | 200 | neumann_corrected_tail | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 4.98276627e-01 | 6.66656299e-01 | exploratory_double |
| signed_xxy | 203 | 200 | neumann_corrected_tail | y_only | 0.00000000e+00 | 3.79488306e-02 | 3.79488306e-02 | 1.77635684e-15 | 5.67340482e-04 | exploratory_double |
| signed_xxy | 203 | 200 | neumann_corrected_tail | omega_up | 1.87894481e-01 | 1.87894481e-02 | 1.88831617e-01 | 1.77635684e-15 | 1.54024239e-01 | exploratory_double |
| signed_xxy | 203 | 200 | neumann_corrected_tail | square | 3.70830003e-02 | 3.70830003e-02 | 5.24432819e-02 | 3.55271368e-15 | 6.63218155e-03 | exploratory_double |

Most important Omega_up-type result:

- best `rx = 2.00000000e-01`, `ry = 2.00000000e-02` with model `signed_xxy`, `J = 30`, `M = `, tail mode `tail_ignored`.
- target `rx >= 1e-2`, `ry = 1e-3`: reached.

Omega_up-type best row by model:

| model | best rx | best ry | J | M | tail mode | Jxx2 lower |
|---|---:|---:|---:|---:|---|---:|
| norm_bound | 3.96025396e-03 | 3.96025396e-04 | 30 |  | tail_ignored | 1.77635684e-15 |
| signed_xxy | 2.00000000e-01 | 2.00000000e-02 | 30 |  | tail_ignored | 1.93281022e+00 |
| symmetry_x_even | 3.39044348e-02 | 3.39044348e-03 | 30 |  | tail_ignored | 3.30583722e-03 |

Target diagnostic at `rx=1e-2`, `ry=1e-3`:

| model | J | M | tail mode | Jxx1 lower | Jxx2 lower | success |
|---|---:|---:|---|---:|---:|---|
| norm_bound | 30 |  | tail_ignored | -4.47236440e+00 | -4.47351279e+00 | False |
| norm_bound | 30 |  | dirichlet_parseval_tail | -6.61635968e+01 | -6.61647452e+01 | False |
| norm_bound | 30 | 29 | neumann_corrected_tail | -1.07557485e+01 | -1.07568969e+01 | False |
| norm_bound | 30 | 31 | neumann_corrected_tail | -1.07557485e+01 | -1.07568969e+01 | False |
| norm_bound | 30 | 50 | neumann_corrected_tail | -9.47781556e+00 | -9.47896396e+00 | False |
| norm_bound | 30 | 98 | neumann_corrected_tail | -8.07413012e+00 | -8.07527851e+00 | False |
| norm_bound | 30 | 101 | neumann_corrected_tail | -8.07413012e+00 | -8.07527851e+00 | False |
| norm_bound | 30 | 200 | neumann_corrected_tail | -6.96821083e+00 | -6.96935922e+00 | False |
| symmetry_x_even | 30 |  | tail_ignored | 2.07253492e+00 | 2.07138652e+00 | True |
| signed_xxy | 30 |  | tail_ignored | 2.90706148e+00 | 2.90591309e+00 | True |
| symmetry_x_even | 30 | 29 | neumann_corrected_tail | -3.44019933e+00 | -3.44134772e+00 | False |
| signed_xxy | 30 | 29 | neumann_corrected_tail | -2.34610976e+00 | -2.34725815e+00 | False |
| symmetry_x_even | 30 | 31 | neumann_corrected_tail | -3.44019933e+00 | -3.44134772e+00 | False |
| signed_xxy | 30 | 31 | neumann_corrected_tail | -2.34610976e+00 | -2.34725815e+00 | False |
| symmetry_x_even | 30 | 50 | neumann_corrected_tail | -2.32004786e+00 | -2.32119626e+00 | False |
| signed_xxy | 30 | 50 | neumann_corrected_tail | -1.23878403e+00 | -1.23993242e+00 | False |
| symmetry_x_even | 30 | 98 | neumann_corrected_tail | -1.09201115e+00 | -1.09315954e+00 | False |
| signed_xxy | 30 | 98 | neumann_corrected_tail | -2.50384233e-02 | -2.61868153e-02 | False |
| symmetry_x_even | 30 | 101 | neumann_corrected_tail | -1.09201115e+00 | -1.09315954e+00 | False |
| signed_xxy | 30 | 101 | neumann_corrected_tail | -2.50384233e-02 | -2.61868153e-02 | False |
| symmetry_x_even | 30 | 200 | neumann_corrected_tail | -1.26303064e-01 | -1.27451456e-01 | False |
| signed_xxy | 30 | 200 | neumann_corrected_tail | 9.29251839e-01 | 9.28103447e-01 | True |
| norm_bound | 49 |  | tail_ignored | -4.50640470e+00 | -4.50755309e+00 | False |
| norm_bound | 49 |  | dirichlet_parseval_tail | -6.49798024e+01 | -6.49809508e+01 | False |
| norm_bound | 49 | 29 | neumann_corrected_tail | -1.06283437e+01 | -1.06294921e+01 | False |
| norm_bound | 49 | 31 | neumann_corrected_tail | -1.06283437e+01 | -1.06294921e+01 | False |
| norm_bound | 49 | 50 | neumann_corrected_tail | -9.37570751e+00 | -9.37685590e+00 | False |
| norm_bound | 49 | 98 | neumann_corrected_tail | -7.99992681e+00 | -8.00107521e+00 | False |
| norm_bound | 49 | 101 | neumann_corrected_tail | -7.99992681e+00 | -8.00107521e+00 | False |
| norm_bound | 49 | 200 | neumann_corrected_tail | -6.91608918e+00 | -6.91723757e+00 | False |
| symmetry_x_even | 49 |  | tail_ignored | 2.04288936e+00 | 2.04174097e+00 | True |
| signed_xxy | 49 |  | tail_ignored | 2.87789177e+00 | 2.87674338e+00 | True |
| symmetry_x_even | 49 | 29 | neumann_corrected_tail | -3.34121358e+00 | -3.34236198e+00 | False |
| signed_xxy | 49 | 29 | neumann_corrected_tail | -2.25585403e+00 | -2.25700242e+00 | False |
| symmetry_x_even | 49 | 31 | neumann_corrected_tail | -3.34121358e+00 | -3.34236198e+00 | False |
| signed_xxy | 49 | 31 | neumann_corrected_tail | -2.25585403e+00 | -2.25700242e+00 | False |
| symmetry_x_even | 49 | 50 | neumann_corrected_tail | -2.24064026e+00 | -2.24178865e+00 | False |
| signed_xxy | 49 | 50 | neumann_corrected_tail | -1.16760078e+00 | -1.16874917e+00 | False |
| symmetry_x_even | 49 | 98 | neumann_corrected_tail | -1.03407793e+00 | -1.03522633e+00 | False |
| signed_xxy | 49 | 98 | neumann_corrected_tail | 2.52393825e-02 | 2.40909905e-02 | True |
| symmetry_x_even | 49 | 101 | neumann_corrected_tail | -1.03407793e+00 | -1.03522633e+00 | False |
| signed_xxy | 49 | 101 | neumann_corrected_tail | 2.52393825e-02 | 2.40909905e-02 | True |
| symmetry_x_even | 49 | 200 | neumann_corrected_tail | -8.52655437e-02 | -8.64139357e-02 | False |
| signed_xxy | 49 | 200 | neumann_corrected_tail | 9.63093034e-01 | 9.61944642e-01 | True |
| norm_bound | 51 |  | tail_ignored | -4.50727550e+00 | -4.50842389e+00 | False |
| norm_bound | 51 |  | dirichlet_parseval_tail | -6.48546952e+01 | -6.48558436e+01 | False |
| norm_bound | 51 | 29 | neumann_corrected_tail | -1.06149293e+01 | -1.06160777e+01 | False |
| norm_bound | 51 | 31 | neumann_corrected_tail | -1.06149293e+01 | -1.06160777e+01 | False |
| norm_bound | 51 | 50 | neumann_corrected_tail | -9.36496571e+00 | -9.36611410e+00 | False |
| norm_bound | 51 | 98 | neumann_corrected_tail | -7.99213277e+00 | -7.99328116e+00 | False |
| norm_bound | 51 | 101 | neumann_corrected_tail | -7.99213277e+00 | -7.99328116e+00 | False |
| norm_bound | 51 | 200 | neumann_corrected_tail | -6.91062744e+00 | -6.91177583e+00 | False |
| symmetry_x_even | 51 |  | tail_ignored | 2.04213004e+00 | 2.04098165e+00 | True |
| signed_xxy | 51 |  | tail_ignored | 2.87714452e+00 | 2.87599613e+00 | True |
| symmetry_x_even | 51 | 29 | neumann_corrected_tail | -3.33076974e+00 | -3.33191813e+00 | False |
| signed_xxy | 51 | 29 | neumann_corrected_tail | -2.24632887e+00 | -2.24747727e+00 | False |
| symmetry_x_even | 51 | 31 | neumann_corrected_tail | -3.33076974e+00 | -3.33191813e+00 | False |
| signed_xxy | 51 | 31 | neumann_corrected_tail | -2.24632887e+00 | -2.24747727e+00 | False |
| symmetry_x_even | 51 | 50 | neumann_corrected_tail | -2.23226873e+00 | -2.23341713e+00 | False |
| signed_xxy | 51 | 50 | neumann_corrected_tail | -1.16009506e+00 | -1.16124345e+00 | False |
| symmetry_x_even | 51 | 98 | neumann_corrected_tail | -1.02797939e+00 | -1.02912779e+00 | False |
| signed_xxy | 51 | 98 | neumann_corrected_tail | 3.05315776e-02 | 2.93831856e-02 | True |
| symmetry_x_even | 51 | 101 | neumann_corrected_tail | -1.02797939e+00 | -1.02912779e+00 | False |
| signed_xxy | 51 | 101 | neumann_corrected_tail | 3.05315776e-02 | 2.93831856e-02 | True |
| symmetry_x_even | 51 | 200 | neumann_corrected_tail | -8.09553203e-02 | -8.21037123e-02 | False |
| signed_xxy | 51 | 200 | neumann_corrected_tail | 9.66644882e-01 | 9.65496490e-01 | True |
| norm_bound | 99 |  | tail_ignored | -4.53186018e+00 | -4.53300857e+00 | False |
| norm_bound | 99 |  | dirichlet_parseval_tail | -6.39498026e+01 | -6.39509510e+01 | False |
| norm_bound | 99 | 29 | neumann_corrected_tail | -1.05181875e+01 | -1.05193358e+01 | False |
| norm_bound | 99 | 31 | neumann_corrected_tail | -1.05181875e+01 | -1.05193358e+01 | False |
| norm_bound | 99 | 50 | neumann_corrected_tail | -9.28754953e+00 | -9.28869792e+00 | False |
| norm_bound | 99 | 98 | neumann_corrected_tail | -7.93602986e+00 | -7.93717825e+00 | False |
| norm_bound | 99 | 101 | neumann_corrected_tail | -7.93602986e+00 | -7.93717825e+00 | False |
| norm_bound | 99 | 200 | neumann_corrected_tail | -6.87138603e+00 | -6.87253442e+00 | False |
| symmetry_x_even | 99 |  | tail_ignored | 2.02066931e+00 | 2.01952092e+00 | True |
| signed_xxy | 99 |  | tail_ignored | 2.85602178e+00 | 2.85487339e+00 | True |
| symmetry_x_even | 99 | 29 | neumann_corrected_tail | -3.25532668e+00 | -3.25647507e+00 | False |
| signed_xxy | 99 | 29 | neumann_corrected_tail | -2.17750851e+00 | -2.17865690e+00 | False |
| symmetry_x_even | 99 | 31 | neumann_corrected_tail | -3.25532668e+00 | -3.25647507e+00 | False |
| signed_xxy | 99 | 31 | neumann_corrected_tail | -2.17750851e+00 | -2.17865690e+00 | False |
| symmetry_x_even | 99 | 50 | neumann_corrected_tail | -2.17183317e+00 | -2.17298156e+00 | False |
| signed_xxy | 99 | 50 | neumann_corrected_tail | -1.10590297e+00 | -1.10705136e+00 | False |
| symmetry_x_even | 99 | 98 | neumann_corrected_tail | -9.84004350e-01 | -9.85152742e-01 | False |
| signed_xxy | 99 | 98 | neumann_corrected_tail | 6.86895550e-02 | 6.75411630e-02 | True |
| symmetry_x_even | 99 | 101 | neumann_corrected_tail | -9.84004350e-01 | -9.85152742e-01 | False |
| signed_xxy | 99 | 101 | neumann_corrected_tail | 6.86895550e-02 | 6.75411630e-02 | True |
| symmetry_x_even | 99 | 200 | neumann_corrected_tail | -4.99307054e-02 | -5.10790974e-02 | False |
| signed_xxy | 99 | 200 | neumann_corrected_tail | 9.92196260e-01 | 9.91047868e-01 | True |
| norm_bound | 101 |  | tail_ignored | -4.53186018e+00 | -4.53300857e+00 | False |
| norm_bound | 101 |  | dirichlet_parseval_tail | -6.39331222e+01 | -6.39342706e+01 | False |
| norm_bound | 101 | 29 | neumann_corrected_tail | -1.05164077e+01 | -1.05175561e+01 | False |
| norm_bound | 101 | 31 | neumann_corrected_tail | -1.05164077e+01 | -1.05175561e+01 | False |
| norm_bound | 101 | 50 | neumann_corrected_tail | -9.28612593e+00 | -9.28727432e+00 | False |
| norm_bound | 101 | 98 | neumann_corrected_tail | -7.93499898e+00 | -7.93614737e+00 | False |
| norm_bound | 101 | 101 | neumann_corrected_tail | -7.93499898e+00 | -7.93614737e+00 | False |
| norm_bound | 101 | 200 | neumann_corrected_tail | -6.87066581e+00 | -6.87181420e+00 | False |
| symmetry_x_even | 101 |  | tail_ignored | 2.02066931e+00 | 2.01952092e+00 | True |
| signed_xxy | 101 |  | tail_ignored | 2.85602178e+00 | 2.85487339e+00 | True |
| symmetry_x_even | 101 | 29 | neumann_corrected_tail | -3.25393673e+00 | -3.25508512e+00 | False |
| signed_xxy | 101 | 29 | neumann_corrected_tail | -2.17624030e+00 | -2.17738869e+00 | False |
| symmetry_x_even | 101 | 31 | neumann_corrected_tail | -3.25393673e+00 | -3.25508512e+00 | False |
| signed_xxy | 101 | 31 | neumann_corrected_tail | -2.17624030e+00 | -2.17738869e+00 | False |
| symmetry_x_even | 101 | 50 | neumann_corrected_tail | -2.17072015e+00 | -2.17186854e+00 | False |
| signed_xxy | 101 | 50 | neumann_corrected_tail | -1.10490475e+00 | -1.10605314e+00 | False |
| symmetry_x_even | 101 | 98 | neumann_corrected_tail | -9.83195071e-01 | -9.84343463e-01 | False |
| signed_xxy | 101 | 98 | neumann_corrected_tail | 6.93918307e-02 | 6.82434387e-02 | True |
| symmetry_x_even | 101 | 101 | neumann_corrected_tail | -9.83195071e-01 | -9.84343463e-01 | False |
| signed_xxy | 101 | 101 | neumann_corrected_tail | 6.93918307e-02 | 6.82434387e-02 | True |
| symmetry_x_even | 101 | 200 | neumann_corrected_tail | -4.93603917e-02 | -5.05087837e-02 | False |
| signed_xxy | 101 | 200 | neumann_corrected_tail | 9.92665856e-01 | 9.91517464e-01 | True |
| norm_bound | 199 |  | tail_ignored | -4.54210700e+00 | -4.54325539e+00 | False |
| norm_bound | 199 |  | dirichlet_parseval_tail | -6.34407396e+01 | -6.34418880e+01 | False |
| norm_bound | 199 | 29 | neumann_corrected_tail | -1.04639429e+01 | -1.04650913e+01 | False |
| norm_bound | 199 | 31 | neumann_corrected_tail | -1.04639429e+01 | -1.04650913e+01 | False |
| norm_bound | 199 | 50 | neumann_corrected_tail | -9.24417205e+00 | -9.24532045e+00 | False |
| norm_bound | 199 | 98 | neumann_corrected_tail | -7.90463539e+00 | -7.90578378e+00 | False |
| norm_bound | 199 | 101 | neumann_corrected_tail | -7.90463539e+00 | -7.90578378e+00 | False |
| norm_bound | 199 | 200 | neumann_corrected_tail | -6.84947023e+00 | -6.85061862e+00 | False |
| symmetry_x_even | 199 |  | tail_ignored | 2.01170876e+00 | 2.01056037e+00 | True |
| signed_xxy | 199 |  | tail_ignored | 2.84720032e+00 | 2.84605193e+00 | True |
| symmetry_x_even | 199 | 29 | neumann_corrected_tail | -3.21292874e+00 | -3.21407714e+00 | False |
| signed_xxy | 199 | 29 | neumann_corrected_tail | -2.13882001e+00 | -2.13996840e+00 | False |
| symmetry_x_even | 199 | 31 | neumann_corrected_tail | -3.21292874e+00 | -3.21407714e+00 | False |
| signed_xxy | 199 | 31 | neumann_corrected_tail | -2.13882001e+00 | -2.13996840e+00 | False |
| symmetry_x_even | 199 | 50 | neumann_corrected_tail | -2.13789170e+00 | -2.13904009e+00 | False |
| signed_xxy | 199 | 50 | neumann_corrected_tail | -1.07546012e+00 | -1.07660852e+00 | False |
| symmetry_x_even | 199 | 98 | neumann_corrected_tail | -9.59337928e-01 | -9.60486320e-01 | False |
| signed_xxy | 199 | 98 | neumann_corrected_tail | 9.00942879e-02 | 8.89458959e-02 | True |
| symmetry_x_even | 199 | 101 | neumann_corrected_tail | -9.59337928e-01 | -9.60486320e-01 | False |
| signed_xxy | 199 | 101 | neumann_corrected_tail | 9.00942879e-02 | 8.89458959e-02 | True |
| symmetry_x_even | 199 | 200 | neumann_corrected_tail | -3.25613440e-02 | -3.37097360e-02 | False |
| signed_xxy | 199 | 200 | neumann_corrected_tail | 1.00649491e+00 | 1.00534652e+00 | True |
| norm_bound | 203 |  | tail_ignored | -4.54284629e+00 | -4.54399468e+00 | False |
| norm_bound | 203 |  | dirichlet_parseval_tail | -6.34362834e+01 | -6.34374318e+01 | False |
| norm_bound | 203 | 29 | neumann_corrected_tail | -1.04634686e+01 | -1.04646170e+01 | False |
| norm_bound | 203 | 31 | neumann_corrected_tail | -1.04634686e+01 | -1.04646170e+01 | False |
| norm_bound | 203 | 50 | neumann_corrected_tail | -9.24379281e+00 | -9.24494120e+00 | False |
| norm_bound | 203 | 98 | neumann_corrected_tail | -7.90436102e+00 | -7.90550941e+00 | False |
| norm_bound | 203 | 101 | neumann_corrected_tail | -7.90436102e+00 | -7.90550941e+00 | False |
| norm_bound | 203 | 200 | neumann_corrected_tail | -6.84927881e+00 | -6.85042720e+00 | False |
| symmetry_x_even | 203 |  | tail_ignored | 2.01106179e+00 | 2.00991340e+00 | True |
| signed_xxy | 203 |  | tail_ignored | 2.84656333e+00 | 2.84541493e+00 | True |
| symmetry_x_even | 203 | 29 | neumann_corrected_tail | -3.21255768e+00 | -3.21370607e+00 | False |
| signed_xxy | 203 | 29 | neumann_corrected_tail | -2.13848137e+00 | -2.13962976e+00 | False |
| symmetry_x_even | 203 | 31 | neumann_corrected_tail | -3.21255768e+00 | -3.21370607e+00 | False |
| signed_xxy | 203 | 31 | neumann_corrected_tail | -2.13848137e+00 | -2.13962976e+00 | False |
| symmetry_x_even | 203 | 50 | neumann_corrected_tail | -2.13759470e+00 | -2.13874310e+00 | False |
| signed_xxy | 203 | 50 | neumann_corrected_tail | -1.07519371e+00 | -1.07634211e+00 | False |
| symmetry_x_even | 203 | 98 | neumann_corrected_tail | -9.59122171e-01 | -9.60270563e-01 | False |
| signed_xxy | 203 | 98 | neumann_corrected_tail | 9.02815254e-02 | 8.91331334e-02 | True |
| symmetry_x_even | 203 | 101 | neumann_corrected_tail | -9.59122171e-01 | -9.60270563e-01 | False |
| signed_xxy | 203 | 101 | neumann_corrected_tail | 9.02815254e-02 | 8.91331334e-02 | True |
| symmetry_x_even | 203 | 200 | neumann_corrected_tail | -3.24094990e-02 | -3.35578910e-02 | False |
| signed_xxy | 203 | 200 | neumann_corrected_tail | 1.00661990e+00 | 1.00547151e+00 | True |

Parseval decomposition diagnostic:

- smallest absolute combined residual in the grid is `2.10965251e+00` at `J=203`, `M=200`.
- Full values are written to `parseval_decomposition.csv`.
