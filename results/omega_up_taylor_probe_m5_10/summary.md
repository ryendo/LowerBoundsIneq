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
| N | 5 | 3.13813978e-14 | 3.06954462e-12 | 4.37357663e-14 |
| N | 6 | 3.13813978e-14 | 3.06954462e-12 | 2.49918665e-14 |
| N | 7 | 3.72659195e-14 | 4.56168436e-12 | 3.71406905e-14 |
| N | 8 | 3.72659195e-14 | 6.30961949e-12 | 3.99561322e-14 |
| N | 9 | 4.62766492e-14 | 6.30961949e-12 | 3.99561322e-14 |
| N | 10 | 4.62766492e-14 | 8.10018719e-12 | 3.84712759e-14 |

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

| J | M | M complete | old tail energy | corrected tail energy | ratio | old lambda_xx lower | corrected lambda_xx lower |
|---:|---:|---|---:|---:|---:|---:|---:|
| 30 | 5 | True | 5.99866146e+01 | 2.54932806e+01 | 4.24982819e-01 | -8.95784831e+01 | -1.71424817e+01 |
| 30 | 6 | False | 5.99866146e+01 | 2.54932806e+01 | 4.24982819e-01 | -8.95784831e+01 | -1.71424817e+01 |
| 30 | 7 | True | 5.99866146e+01 | 1.26713866e+01 | 2.11236902e-01 | -8.95784831e+01 | 9.78349560e+00 |
| 30 | 8 | False | 5.99866146e+01 | 1.26713866e+01 | 2.11236902e-01 | -8.95784831e+01 | 9.78349560e+00 |
| 30 | 9 | True | 5.99866146e+01 | 1.26713866e+01 | 2.11236902e-01 | -8.95784831e+01 | 9.78349560e+00 |
| 30 | 10 | True | 5.99866146e+01 | 1.26713866e+01 | 2.11236902e-01 | -8.95784831e+01 | 9.78349560e+00 |
| 49 | 5 | True | 5.99541877e+01 | 2.54608537e+01 | 4.24671815e-01 | -8.74092734e+01 | -1.62209032e+01 |
| 49 | 6 | False | 5.99541877e+01 | 2.54608537e+01 | 4.24671815e-01 | -8.74092734e+01 | -1.62209032e+01 |
| 49 | 7 | True | 5.99541877e+01 | 1.26389597e+01 | 2.10810291e-01 | -8.74092734e+01 | 1.02413035e+01 |
| 49 | 8 | False | 5.99541877e+01 | 1.26389597e+01 | 2.10810291e-01 | -8.74092734e+01 | 1.02413035e+01 |
| 49 | 9 | True | 5.99541877e+01 | 1.26389597e+01 | 2.10810291e-01 | -8.74092734e+01 | 1.02413035e+01 |
| 49 | 10 | True | 5.99541877e+01 | 1.26389597e+01 | 2.10810291e-01 | -8.74092734e+01 | 1.02413035e+01 |
| 51 | 5 | True | 5.99533506e+01 | 2.54600166e+01 | 4.24663781e-01 | -8.71796648e+01 | -1.61233968e+01 |
| 51 | 6 | False | 5.99533506e+01 | 2.54600166e+01 | 4.24663781e-01 | -8.71796648e+01 | -1.61233968e+01 |
| 51 | 7 | True | 5.99533506e+01 | 1.26381226e+01 | 2.10799271e-01 | -8.71796648e+01 | 1.02897048e+01 |
| 51 | 8 | False | 5.99533506e+01 | 1.26381226e+01 | 2.10799271e-01 | -8.71796648e+01 | 1.02897048e+01 |
| 51 | 9 | True | 5.99533506e+01 | 1.26381226e+01 | 2.10799271e-01 | -8.71796648e+01 | 1.02897048e+01 |
| 51 | 10 | True | 5.99533506e+01 | 1.26381226e+01 | 2.10799271e-01 | -8.71796648e+01 | 1.02897048e+01 |
| 99 | 5 | True | 5.99295212e+01 | 2.54361872e+01 | 4.24435014e-01 | -8.55168526e+01 | -1.54174964e+01 |
| 99 | 6 | False | 5.99295212e+01 | 2.54361872e+01 | 4.24435014e-01 | -8.55168526e+01 | -1.54174964e+01 |
| 99 | 7 | True | 5.99295212e+01 | 1.26142932e+01 | 2.10485466e-01 | -8.55168526e+01 | 1.06399010e+01 |
| 99 | 8 | False | 5.99295212e+01 | 1.26142932e+01 | 2.10485466e-01 | -8.55168526e+01 | 1.06399010e+01 |
| 99 | 9 | True | 5.99295212e+01 | 1.26142932e+01 | 2.10485466e-01 | -8.55168526e+01 | 1.06399010e+01 |
| 99 | 10 | True | 5.99295212e+01 | 1.26142932e+01 | 2.10485466e-01 | -8.55168526e+01 | 1.06399010e+01 |
| 101 | 5 | True | 5.99295212e+01 | 2.54361872e+01 | 4.24435014e-01 | -8.54861668e+01 | -1.54044723e+01 |
| 101 | 6 | False | 5.99295212e+01 | 2.54361872e+01 | 4.24435014e-01 | -8.54861668e+01 | -1.54044723e+01 |
| 101 | 7 | True | 5.99295212e+01 | 1.26142932e+01 | 2.10485466e-01 | -8.54861668e+01 | 1.06463599e+01 |
| 101 | 8 | False | 5.99295212e+01 | 1.26142932e+01 | 2.10485466e-01 | -8.54861668e+01 | 1.06463599e+01 |
| 101 | 9 | True | 5.99295212e+01 | 1.26142932e+01 | 2.10485466e-01 | -8.54861668e+01 | 1.06463599e+01 |
| 101 | 10 | True | 5.99295212e+01 | 1.26142932e+01 | 2.10485466e-01 | -8.54861668e+01 | 1.06463599e+01 |
| 199 | 5 | True | 5.99194565e+01 | 2.54261225e+01 | 4.24338336e-01 | -8.45798000e+01 | -1.50198355e+01 |
| 199 | 6 | False | 5.99194565e+01 | 2.54261225e+01 | 4.24338336e-01 | -8.45798000e+01 | -1.50198355e+01 |
| 199 | 7 | True | 5.99194565e+01 | 1.26042285e+01 | 2.10352851e-01 | -8.45798000e+01 | 1.08370587e+01 |
| 199 | 8 | False | 5.99194565e+01 | 1.26042285e+01 | 2.10352851e-01 | -8.45798000e+01 | 1.08370587e+01 |
| 199 | 9 | True | 5.99194565e+01 | 1.26042285e+01 | 2.10352851e-01 | -8.45798000e+01 | 1.08370587e+01 |
| 199 | 10 | True | 5.99194565e+01 | 1.26042285e+01 | 2.10352851e-01 | -8.45798000e+01 | 1.08370587e+01 |
| 203 | 5 | True | 5.99187262e+01 | 2.54253921e+01 | 4.24331320e-01 | -8.45715922e+01 | -1.50163527e+01 |
| 203 | 6 | False | 5.99187262e+01 | 2.54253921e+01 | 4.24331320e-01 | -8.45715922e+01 | -1.50163527e+01 |
| 203 | 7 | True | 5.99187262e+01 | 1.26034982e+01 | 2.10343226e-01 | -8.45715922e+01 | 1.08387851e+01 |
| 203 | 8 | False | 5.99187262e+01 | 1.26034982e+01 | 2.10343226e-01 | -8.45715922e+01 | 1.08387851e+01 |
| 203 | 9 | True | 5.99187262e+01 | 1.26034982e+01 | 2.10343226e-01 | -8.45715922e+01 | 1.08387851e+01 |
| 203 | 10 | True | 5.99187262e+01 | 1.26034982e+01 | 2.10343226e-01 | -8.45715922e+01 | 1.08387851e+01 |
- Smallest corrected tail in this run: `1.26034982e+01` at `J=203`, `M=7`.
- Full comparison is written to `neumann_corrected_tail.csv`.

Cell-size bisection results:

| model | J | M | tail mode | cell shape | max rx | max ry | distance | Jxx1 lower | Jxx2 lower | rigorous flag |
|---|---:|---:|---|---|---:|---:|---:|---:|---:|---|
| norm_bound | 30 |  | tail_ignored | x_only | 4.48081823e-03 | 0.00000000e+00 | 4.48081823e-03 | 1.48524942e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 |  | tail_ignored | y_only | 0.00000000e+00 | 3.40823079e-03 | 3.40823079e-03 | 1.55790964e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 |  | tail_ignored | omega_up | 3.96025396e-03 | 3.96025396e-04 | 3.98000597e-03 | 1.50399992e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 |  | tail_ignored | square | 1.93590103e-03 | 1.93590103e-03 | 2.73777749e-03 | 1.55122141e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| norm_bound | 30 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| norm_bound | 30 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| norm_bound | 30 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| norm_bound | 30 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| norm_bound | 30 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| norm_bound | 30 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| norm_bound | 30 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| norm_bound | 30 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 30 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.40091519e+00 | 2.56929486e+00 | exploratory_double |
| symmetry_x_even | 30 |  | tail_ignored | y_only | 0.00000000e+00 | 3.40823079e-03 | 3.40823079e-03 | 1.55790964e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 30 |  | tail_ignored | omega_up | 3.39044348e-02 | 3.39044348e-03 | 3.40735353e-02 | 1.77635684e-15 | 3.30583722e-03 | exploratory_double |
| symmetry_x_even | 30 |  | tail_ignored | square | 3.40808954e-03 | 3.40808954e-03 | 4.81976645e-03 | 1.50878081e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 30 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.40091519e+00 | 2.56929486e+00 | exploratory_double |
| signed_xxy | 30 |  | tail_ignored | y_only | 0.00000000e+00 | 9.88637327e-02 | 9.88637327e-02 | 1.77635684e-15 | 1.52399842e-02 | exploratory_double |
| signed_xxy | 30 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.75747456e+00 | 1.93281022e+00 | exploratory_double |
| signed_xxy | 30 |  | tail_ignored | square | 9.24585593e-02 | 9.24585593e-02 | 1.30756149e-01 | 8.88178420e-15 | 5.52099788e-02 | exploratory_double |
| symmetry_x_even | 30 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| symmetry_x_even | 30 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| symmetry_x_even | 30 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| symmetry_x_even | 30 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| signed_xxy | 30 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| signed_xxy | 30 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| signed_xxy | 30 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| signed_xxy | 30 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| symmetry_x_even | 30 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| symmetry_x_even | 30 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| symmetry_x_even | 30 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| symmetry_x_even | 30 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| signed_xxy | 30 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| signed_xxy | 30 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| signed_xxy | 30 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| signed_xxy | 30 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -2.02439045e+01 | -2.02454742e+01 | exploratory_double |
| symmetry_x_even | 30 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| symmetry_x_even | 30 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| signed_xxy | 30 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.58461434e+00 | -8.58618402e+00 | exploratory_double |
| norm_bound | 49 |  | tail_ignored | x_only | 4.43327938e-03 | 0.00000000e+00 | 4.43327938e-03 | 1.48703137e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 |  | tail_ignored | y_only | 0.00000000e+00 | 3.37247287e-03 | 3.37247287e-03 | 1.55820456e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 |  | tail_ignored | omega_up | 3.91829128e-03 | 3.91829128e-04 | 3.93783400e-03 | 1.50539047e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 |  | tail_ignored | square | 1.91549070e-03 | 1.91549070e-03 | 2.70891293e-03 | 1.55163672e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 49 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| norm_bound | 49 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| norm_bound | 49 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| norm_bound | 49 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| norm_bound | 49 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| norm_bound | 49 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| norm_bound | 49 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| norm_bound | 49 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| norm_bound | 49 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 49 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.37171176e+00 | 2.54009143e+00 | exploratory_double |
| symmetry_x_even | 49 |  | tail_ignored | y_only | 0.00000000e+00 | 3.37247287e-03 | 3.37247287e-03 | 1.55820456e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 49 |  | tail_ignored | omega_up | 3.35510519e-02 | 3.35510519e-03 | 3.37183899e-02 | 1.77635684e-15 | 3.20435612e-03 | exploratory_double |
| symmetry_x_even | 49 |  | tail_ignored | square | 3.37233468e-03 | 3.37233468e-03 | 4.76920145e-03 | 1.51010413e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 49 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.37171176e+00 | 2.54009143e+00 | exploratory_double |
| signed_xxy | 49 |  | tail_ignored | y_only | 0.00000000e+00 | 9.80019278e-02 | 9.80019278e-02 | 7.10542736e-15 | 1.49158462e-02 | exploratory_double |
| signed_xxy | 49 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.72894555e+00 | 1.90428121e+00 | exploratory_double |
| signed_xxy | 49 |  | tail_ignored | square | 9.17137311e-02 | 9.17137311e-02 | 1.29702802e-01 | 3.55271368e-15 | 5.42180721e-02 | exploratory_double |
| symmetry_x_even | 49 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| symmetry_x_even | 49 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| symmetry_x_even | 49 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| symmetry_x_even | 49 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| signed_xxy | 49 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| signed_xxy | 49 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| signed_xxy | 49 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| signed_xxy | 49 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| symmetry_x_even | 49 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| symmetry_x_even | 49 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| symmetry_x_even | 49 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| symmetry_x_even | 49 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| signed_xxy | 49 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| signed_xxy | 49 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| signed_xxy | 49 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| signed_xxy | 49 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98448493e+01 | -1.98464190e+01 | exploratory_double |
| symmetry_x_even | 49 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| symmetry_x_even | 49 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| signed_xxy | 49 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.38637772e+00 | -8.38794739e+00 | exploratory_double |
| norm_bound | 51 |  | tail_ignored | x_only | 4.43206314e-03 | 0.00000000e+00 | 4.43206314e-03 | 1.48707671e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 |  | tail_ignored | y_only | 0.00000000e+00 | 3.37155785e-03 | 3.37155785e-03 | 1.55821206e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 |  | tail_ignored | omega_up | 3.91721767e-03 | 3.91721767e-04 | 3.93675504e-03 | 1.50542585e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 |  | tail_ignored | square | 1.91496846e-03 | 1.91496846e-03 | 2.70817437e-03 | 1.55164729e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 51 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| norm_bound | 51 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| norm_bound | 51 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| norm_bound | 51 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| norm_bound | 51 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| norm_bound | 51 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| norm_bound | 51 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| norm_bound | 51 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| norm_bound | 51 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 51 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.37096364e+00 | 2.53934331e+00 | exploratory_double |
| symmetry_x_even | 51 |  | tail_ignored | y_only | 0.00000000e+00 | 3.37155785e-03 | 3.37155785e-03 | 1.55821206e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 51 |  | tail_ignored | omega_up | 3.35420084e-02 | 3.35420084e-03 | 3.37093013e-02 | 1.77635684e-15 | 3.20177319e-03 | exploratory_double |
| symmetry_x_even | 51 |  | tail_ignored | square | 3.37141974e-03 | 3.37141974e-03 | 4.76790752e-03 | 1.51013780e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 51 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.37096364e+00 | 2.53934331e+00 | exploratory_double |
| signed_xxy | 51 |  | tail_ignored | y_only | 0.00000000e+00 | 9.79798235e-02 | 9.79798235e-02 | 5.32907052e-15 | 1.49075816e-02 | exploratory_double |
| signed_xxy | 51 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.72821471e+00 | 1.90355037e+00 | exploratory_double |
| signed_xxy | 51 |  | tail_ignored | square | 9.16946123e-02 | 9.16946123e-02 | 1.29675764e-01 | 1.77635684e-15 | 5.41927448e-02 | exploratory_double |
| symmetry_x_even | 51 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| symmetry_x_even | 51 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| symmetry_x_even | 51 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| symmetry_x_even | 51 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| signed_xxy | 51 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| signed_xxy | 51 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| signed_xxy | 51 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| signed_xxy | 51 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| symmetry_x_even | 51 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| symmetry_x_even | 51 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| symmetry_x_even | 51 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| symmetry_x_even | 51 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| signed_xxy | 51 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| signed_xxy | 51 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| signed_xxy | 51 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| signed_xxy | 51 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.98026278e+01 | -1.98041975e+01 | exploratory_double |
| symmetry_x_even | 51 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| symmetry_x_even | 51 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| signed_xxy | 51 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.36541933e+00 | -8.36698901e+00 | exploratory_double |
| norm_bound | 99 |  | tail_ignored | x_only | 4.39771885e-03 | 0.00000000e+00 | 4.39771885e-03 | 1.48835190e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 |  | tail_ignored | y_only | 0.00000000e+00 | 3.34571527e-03 | 3.34571527e-03 | 1.55842286e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 |  | tail_ignored | omega_up | 3.88690061e-03 | 3.88690061e-04 | 3.90628677e-03 | 1.50642096e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 |  | tail_ignored | square | 1.90022006e-03 | 1.90022006e-03 | 2.68731698e-03 | 1.55194437e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 99 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| norm_bound | 99 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| norm_bound | 99 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| norm_bound | 99 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| norm_bound | 99 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| norm_bound | 99 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| norm_bound | 99 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| norm_bound | 99 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| norm_bound | 99 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 99 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34981648e+00 | 2.51819616e+00 | exploratory_double |
| symmetry_x_even | 99 |  | tail_ignored | y_only | 0.00000000e+00 | 3.34571527e-03 | 3.34571527e-03 | 1.55842286e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 99 |  | tail_ignored | omega_up | 3.32865834e-02 | 3.32865834e-03 | 3.34526023e-02 | 1.77635684e-15 | 3.12911161e-03 | exploratory_double |
| symmetry_x_even | 99 |  | tail_ignored | square | 3.34557935e-03 | 3.34557935e-03 | 4.73136369e-03 | 1.51108481e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 99 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34981648e+00 | 2.51819616e+00 | exploratory_double |
| signed_xxy | 99 |  | tail_ignored | y_only | 0.00000000e+00 | 9.73544387e-02 | 9.73544387e-02 | 1.77635684e-15 | 1.46747705e-02 | exploratory_double |
| signed_xxy | 99 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.70755593e+00 | 1.88289159e+00 | exploratory_double |
| signed_xxy | 99 |  | tail_ignored | square | 9.11533838e-02 | 9.11533838e-02 | 1.28910352e-01 | 5.32907052e-15 | 5.34785288e-02 | exploratory_double |
| symmetry_x_even | 99 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| symmetry_x_even | 99 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| symmetry_x_even | 99 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| symmetry_x_even | 99 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| signed_xxy | 99 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| signed_xxy | 99 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| signed_xxy | 99 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| signed_xxy | 99 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| symmetry_x_even | 99 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| symmetry_x_even | 99 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| symmetry_x_even | 99 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| symmetry_x_even | 99 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| signed_xxy | 99 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| signed_xxy | 99 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| signed_xxy | 99 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| signed_xxy | 99 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94969640e+01 | -1.94985337e+01 | exploratory_double |
| symmetry_x_even | 99 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| symmetry_x_even | 99 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| signed_xxy | 99 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21377994e+00 | -8.21534961e+00 | exploratory_double |
| norm_bound | 101 |  | tail_ignored | x_only | 4.39771885e-03 | 0.00000000e+00 | 4.39771885e-03 | 1.48835190e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 |  | tail_ignored | y_only | 0.00000000e+00 | 3.34571527e-03 | 3.34571527e-03 | 1.55842286e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 |  | tail_ignored | omega_up | 3.88690061e-03 | 3.88690061e-04 | 3.90628677e-03 | 1.50642096e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 |  | tail_ignored | square | 1.90022006e-03 | 1.90022006e-03 | 2.68731698e-03 | 1.55194437e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 101 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| norm_bound | 101 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| norm_bound | 101 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| norm_bound | 101 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| norm_bound | 101 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| norm_bound | 101 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| norm_bound | 101 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| norm_bound | 101 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| norm_bound | 101 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 101 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34981648e+00 | 2.51819616e+00 | exploratory_double |
| symmetry_x_even | 101 |  | tail_ignored | y_only | 0.00000000e+00 | 3.34571527e-03 | 3.34571527e-03 | 1.55842286e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 101 |  | tail_ignored | omega_up | 3.32865834e-02 | 3.32865834e-03 | 3.34526023e-02 | 1.77635684e-15 | 3.12911161e-03 | exploratory_double |
| symmetry_x_even | 101 |  | tail_ignored | square | 3.34557935e-03 | 3.34557935e-03 | 4.73136369e-03 | 1.51108481e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 101 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34981648e+00 | 2.51819616e+00 | exploratory_double |
| signed_xxy | 101 |  | tail_ignored | y_only | 0.00000000e+00 | 9.73544387e-02 | 9.73544387e-02 | 1.77635684e-15 | 1.46747705e-02 | exploratory_double |
| signed_xxy | 101 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.70755593e+00 | 1.88289159e+00 | exploratory_double |
| signed_xxy | 101 |  | tail_ignored | square | 9.11533838e-02 | 9.11533838e-02 | 1.28910352e-01 | 5.32907052e-15 | 5.34785288e-02 | exploratory_double |
| symmetry_x_even | 101 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| symmetry_x_even | 101 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| symmetry_x_even | 101 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| symmetry_x_even | 101 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| signed_xxy | 101 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| signed_xxy | 101 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| signed_xxy | 101 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| signed_xxy | 101 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| symmetry_x_even | 101 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| symmetry_x_even | 101 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| symmetry_x_even | 101 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| symmetry_x_even | 101 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| signed_xxy | 101 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| signed_xxy | 101 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| signed_xxy | 101 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| signed_xxy | 101 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.94913244e+01 | -1.94928940e+01 | exploratory_double |
| symmetry_x_even | 101 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| symmetry_x_even | 101 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| signed_xxy | 101 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.21098314e+00 | -8.21255281e+00 | exploratory_double |
| norm_bound | 199 |  | tail_ignored | x_only | 4.38339766e-03 | 0.00000000e+00 | 4.38339766e-03 | 1.48888070e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 |  | tail_ignored | y_only | 0.00000000e+00 | 3.33493667e-03 | 3.33493667e-03 | 1.55851022e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 |  | tail_ignored | omega_up | 3.87425839e-03 | 3.87425839e-04 | 3.89358149e-03 | 1.50683362e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 199 |  | tail_ignored | square | 1.89406933e-03 | 1.89406933e-03 | 2.67861853e-03 | 1.55206753e-03 | 3.55271368e-15 | exploratory_double |
| norm_bound | 199 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| norm_bound | 199 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| norm_bound | 199 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| norm_bound | 199 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| norm_bound | 199 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| norm_bound | 199 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| norm_bound | 199 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| norm_bound | 199 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| norm_bound | 199 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 199 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34098483e+00 | 2.50936450e+00 | exploratory_double |
| symmetry_x_even | 199 |  | tail_ignored | y_only | 0.00000000e+00 | 3.33493667e-03 | 3.33493667e-03 | 1.55851022e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 199 |  | tail_ignored | omega_up | 3.31800413e-02 | 3.31800413e-03 | 3.33455288e-02 | 1.77635684e-15 | 3.09896911e-03 | exploratory_double |
| symmetry_x_even | 199 |  | tail_ignored | square | 3.33480166e-03 | 3.33480166e-03 | 4.71612173e-03 | 1.51147754e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 199 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34098483e+00 | 2.50936450e+00 | exploratory_double |
| signed_xxy | 199 |  | tail_ignored | y_only | 0.00000000e+00 | 9.70929393e-02 | 9.70929393e-02 | 1.77635684e-15 | 1.45780031e-02 | exploratory_double |
| signed_xxy | 199 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.69892823e+00 | 1.87426389e+00 | exploratory_double |
| signed_xxy | 199 |  | tail_ignored | square | 9.09268969e-02 | 9.09268969e-02 | 1.28590051e-01 | 3.55271368e-15 | 5.31812372e-02 | exploratory_double |
| symmetry_x_even | 199 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| symmetry_x_even | 199 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| symmetry_x_even | 199 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| symmetry_x_even | 199 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| signed_xxy | 199 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| signed_xxy | 199 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| signed_xxy | 199 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| signed_xxy | 199 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| symmetry_x_even | 199 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| symmetry_x_even | 199 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| symmetry_x_even | 199 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| symmetry_x_even | 199 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| signed_xxy | 199 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| signed_xxy | 199 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| signed_xxy | 199 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| signed_xxy | 199 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93247718e+01 | -1.93263414e+01 | exploratory_double |
| symmetry_x_even | 199 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| symmetry_x_even | 199 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| signed_xxy | 199 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12840815e+00 | -8.12997783e+00 | exploratory_double |
| norm_bound | 203 |  | tail_ignored | x_only | 4.38236413e-03 | 0.00000000e+00 | 4.38236413e-03 | 1.48891880e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 |  | tail_ignored | y_only | 0.00000000e+00 | 3.33415874e-03 | 3.33415874e-03 | 1.55851651e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 |  | tail_ignored | omega_up | 3.87334601e-03 | 3.87334601e-04 | 3.89266457e-03 | 1.50686335e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 |  | tail_ignored | square | 1.89362542e-03 | 1.89362542e-03 | 2.67799075e-03 | 1.55207641e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 203 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| norm_bound | 203 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| norm_bound | 203 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| norm_bound | 203 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| norm_bound | 203 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| norm_bound | 203 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| norm_bound | 203 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| norm_bound | 203 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| norm_bound | 203 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| norm_bound | 203 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34034709e+00 | 2.50872677e+00 | exploratory_double |
| symmetry_x_even | 203 |  | tail_ignored | y_only | 0.00000000e+00 | 3.33415874e-03 | 3.33415874e-03 | 1.55851651e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_x_even | 203 |  | tail_ignored | omega_up | 3.31723515e-02 | 3.31723515e-03 | 3.33378007e-02 | 1.77635684e-15 | 3.09679732e-03 | exploratory_double |
| symmetry_x_even | 203 |  | tail_ignored | square | 3.33402379e-03 | 3.33402379e-03 | 4.71502166e-03 | 1.51150584e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy | 203 |  | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34034709e+00 | 2.50872677e+00 | exploratory_double |
| signed_xxy | 203 |  | tail_ignored | y_only | 0.00000000e+00 | 9.70740492e-02 | 9.70740492e-02 | 7.10542736e-15 | 1.45710261e-02 | exploratory_double |
| signed_xxy | 203 |  | tail_ignored | omega_up | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.69830523e+00 | 1.87364089e+00 | exploratory_double |
| signed_xxy | 203 |  | tail_ignored | square | 9.09105320e-02 | 9.09105320e-02 | 1.28566907e-01 | 3.55271368e-15 | 5.31597925e-02 | exploratory_double |
| symmetry_x_even | 203 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| symmetry_x_even | 203 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| symmetry_x_even | 203 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| symmetry_x_even | 203 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| signed_xxy | 203 | 5 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| signed_xxy | 203 | 5 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| signed_xxy | 203 | 5 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| signed_xxy | 203 | 5 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| symmetry_x_even | 203 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| symmetry_x_even | 203 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| symmetry_x_even | 203 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| symmetry_x_even | 203 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| signed_xxy | 203 | 6 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| signed_xxy | 203 | 6 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| signed_xxy | 203 | 6 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| signed_xxy | 203 | 6 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -1.93232637e+01 | -1.93248333e+01 | exploratory_double |
| symmetry_x_even | 203 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 7 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 7 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 7 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 7 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 8 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 8 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 8 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 8 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 9 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 9 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 9 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 9 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| symmetry_x_even | 203 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 10 | neumann_corrected_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 10 | neumann_corrected_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 10 | neumann_corrected_tail | omega_up | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |
| signed_xxy | 203 | 10 | neumann_corrected_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -8.12766057e+00 | -8.12923025e+00 | exploratory_double |

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
| norm_bound | 30 | 5 | neumann_corrected_tail | -3.11530937e+01 | -3.11542421e+01 | False |
| norm_bound | 30 | 6 | neumann_corrected_tail | -3.11530937e+01 | -3.11542421e+01 | False |
| norm_bound | 30 | 7 | neumann_corrected_tail | -1.79339800e+01 | -1.79351284e+01 | False |
| norm_bound | 30 | 8 | neumann_corrected_tail | -1.79339800e+01 | -1.79351284e+01 | False |
| norm_bound | 30 | 9 | neumann_corrected_tail | -1.79339800e+01 | -1.79351284e+01 | False |
| norm_bound | 30 | 10 | neumann_corrected_tail | -1.79339800e+01 | -1.79351284e+01 | False |
| symmetry_x_even | 30 |  | tail_ignored | 2.07253492e+00 | 2.07138652e+00 | True |
| signed_xxy | 30 |  | tail_ignored | 2.90706148e+00 | 2.90591309e+00 | True |
| symmetry_x_even | 30 | 5 | neumann_corrected_tail | -2.15285099e+01 | -2.15296583e+01 | False |
| signed_xxy | 30 | 5 | neumann_corrected_tail | -2.02478906e+01 | -2.02490390e+01 | False |
| symmetry_x_even | 30 | 6 | neumann_corrected_tail | -2.15285099e+01 | -2.15296583e+01 | False |
| signed_xxy | 30 | 6 | neumann_corrected_tail | -2.02478906e+01 | -2.02490390e+01 | False |
| symmetry_x_even | 30 | 7 | neumann_corrected_tail | -9.76531657e+00 | -9.76646496e+00 | False |
| signed_xxy | 30 | 7 | neumann_corrected_tail | -8.60206340e+00 | -8.60321179e+00 | False |
| symmetry_x_even | 30 | 8 | neumann_corrected_tail | -9.76531657e+00 | -9.76646496e+00 | False |
| signed_xxy | 30 | 8 | neumann_corrected_tail | -8.60206340e+00 | -8.60321179e+00 | False |
| symmetry_x_even | 30 | 9 | neumann_corrected_tail | -9.76531657e+00 | -9.76646496e+00 | False |
| signed_xxy | 30 | 9 | neumann_corrected_tail | -8.60206340e+00 | -8.60321179e+00 | False |
| symmetry_x_even | 30 | 10 | neumann_corrected_tail | -9.76531657e+00 | -9.76646496e+00 | False |
| signed_xxy | 30 | 10 | neumann_corrected_tail | -8.60206340e+00 | -8.60321179e+00 | False |
| norm_bound | 49 |  | tail_ignored | -4.50640470e+00 | -4.50755309e+00 | False |
| norm_bound | 49 | 5 | neumann_corrected_tail | -3.06307243e+01 | -3.06318727e+01 | False |
| norm_bound | 49 | 6 | neumann_corrected_tail | -3.06307243e+01 | -3.06318727e+01 | False |
| norm_bound | 49 | 7 | neumann_corrected_tail | -1.76660117e+01 | -1.76671601e+01 | False |
| norm_bound | 49 | 8 | neumann_corrected_tail | -1.76660117e+01 | -1.76671601e+01 | False |
| norm_bound | 49 | 9 | neumann_corrected_tail | -1.76660117e+01 | -1.76671601e+01 | False |
| norm_bound | 49 | 10 | neumann_corrected_tail | -1.76660117e+01 | -1.76671601e+01 | False |
| symmetry_x_even | 49 |  | tail_ignored | 2.04288936e+00 | 2.04174097e+00 | True |
| signed_xxy | 49 |  | tail_ignored | 2.87789177e+00 | 2.87674338e+00 | True |
| symmetry_x_even | 49 | 5 | neumann_corrected_tail | -2.11141828e+01 | -2.11153311e+01 | False |
| signed_xxy | 49 | 5 | neumann_corrected_tail | -1.98492962e+01 | -1.98504446e+01 | False |
| symmetry_x_even | 49 | 6 | neumann_corrected_tail | -2.11141828e+01 | -2.11153311e+01 | False |
| signed_xxy | 49 | 6 | neumann_corrected_tail | -1.98492962e+01 | -1.98504446e+01 | False |
| symmetry_x_even | 49 | 7 | neumann_corrected_tail | -9.55591814e+00 | -9.55706654e+00 | False |
| signed_xxy | 49 | 7 | neumann_corrected_tail | -8.40405568e+00 | -8.40520407e+00 | False |
| symmetry_x_even | 49 | 8 | neumann_corrected_tail | -9.55591814e+00 | -9.55706654e+00 | False |
| signed_xxy | 49 | 8 | neumann_corrected_tail | -8.40405568e+00 | -8.40520407e+00 | False |
| symmetry_x_even | 49 | 9 | neumann_corrected_tail | -9.55591814e+00 | -9.55706654e+00 | False |
| signed_xxy | 49 | 9 | neumann_corrected_tail | -8.40405568e+00 | -8.40520407e+00 | False |
| symmetry_x_even | 49 | 10 | neumann_corrected_tail | -9.55591814e+00 | -9.55706654e+00 | False |
| signed_xxy | 49 | 10 | neumann_corrected_tail | -8.40405568e+00 | -8.40520407e+00 | False |
| norm_bound | 51 |  | tail_ignored | -4.50727550e+00 | -4.50842389e+00 | False |
| norm_bound | 51 | 5 | neumann_corrected_tail | -3.05755599e+01 | -3.05767083e+01 | False |
| norm_bound | 51 | 6 | neumann_corrected_tail | -3.05755599e+01 | -3.05767083e+01 | False |
| norm_bound | 51 | 7 | neumann_corrected_tail | -1.76377425e+01 | -1.76388909e+01 | False |
| norm_bound | 51 | 8 | neumann_corrected_tail | -1.76377425e+01 | -1.76388909e+01 | False |
| norm_bound | 51 | 9 | neumann_corrected_tail | -1.76377425e+01 | -1.76388909e+01 | False |
| norm_bound | 51 | 10 | neumann_corrected_tail | -1.76377425e+01 | -1.76388909e+01 | False |
| symmetry_x_even | 51 |  | tail_ignored | 2.04213004e+00 | 2.04098165e+00 | True |
| signed_xxy | 51 |  | tail_ignored | 2.87714452e+00 | 2.87599613e+00 | True |
| symmetry_x_even | 51 | 5 | neumann_corrected_tail | -2.10703578e+01 | -2.10715062e+01 | False |
| signed_xxy | 51 | 5 | neumann_corrected_tail | -1.98071234e+01 | -1.98082718e+01 | False |
| symmetry_x_even | 51 | 6 | neumann_corrected_tail | -2.10703578e+01 | -2.10715062e+01 | False |
| signed_xxy | 51 | 6 | neumann_corrected_tail | -1.98071234e+01 | -1.98082718e+01 | False |
| symmetry_x_even | 51 | 7 | neumann_corrected_tail | -9.53378680e+00 | -9.53493519e+00 | False |
| signed_xxy | 51 | 7 | neumann_corrected_tail | -8.38312149e+00 | -8.38426989e+00 | False |
| symmetry_x_even | 51 | 8 | neumann_corrected_tail | -9.53378680e+00 | -9.53493519e+00 | False |
| signed_xxy | 51 | 8 | neumann_corrected_tail | -8.38312149e+00 | -8.38426989e+00 | False |
| symmetry_x_even | 51 | 9 | neumann_corrected_tail | -9.53378680e+00 | -9.53493519e+00 | False |
| signed_xxy | 51 | 9 | neumann_corrected_tail | -8.38312149e+00 | -8.38426989e+00 | False |
| symmetry_x_even | 51 | 10 | neumann_corrected_tail | -9.53378680e+00 | -9.53493519e+00 | False |
| signed_xxy | 51 | 10 | neumann_corrected_tail | -8.38312149e+00 | -8.38426989e+00 | False |
| norm_bound | 99 |  | tail_ignored | -4.53186018e+00 | -4.53300857e+00 | False |
| norm_bound | 99 | 5 | neumann_corrected_tail | -3.01767944e+01 | -3.01779428e+01 | False |
| norm_bound | 99 | 6 | neumann_corrected_tail | -3.01767944e+01 | -3.01779428e+01 | False |
| norm_bound | 99 | 7 | neumann_corrected_tail | -1.74335600e+01 | -1.74347084e+01 | False |
| norm_bound | 99 | 8 | neumann_corrected_tail | -1.74335600e+01 | -1.74347084e+01 | False |
| norm_bound | 99 | 9 | neumann_corrected_tail | -1.74335600e+01 | -1.74347084e+01 | False |
| norm_bound | 99 | 10 | neumann_corrected_tail | -1.74335600e+01 | -1.74347084e+01 | False |
| symmetry_x_even | 99 |  | tail_ignored | 2.02066931e+00 | 2.01952092e+00 | True |
| signed_xxy | 99 |  | tail_ignored | 2.85602178e+00 | 2.85487339e+00 | True |
| symmetry_x_even | 99 | 5 | neumann_corrected_tail | -2.07531569e+01 | -2.07543053e+01 | False |
| signed_xxy | 99 | 5 | neumann_corrected_tail | -1.95018125e+01 | -1.95029609e+01 | False |
| symmetry_x_even | 99 | 6 | neumann_corrected_tail | -2.07531569e+01 | -2.07543053e+01 | False |
| signed_xxy | 99 | 6 | neumann_corrected_tail | -1.95018125e+01 | -1.95029609e+01 | False |
| symmetry_x_even | 99 | 7 | neumann_corrected_tail | -9.37370131e+00 | -9.37484970e+00 | False |
| signed_xxy | 99 | 7 | neumann_corrected_tail | -8.23165720e+00 | -8.23280559e+00 | False |
| symmetry_x_even | 99 | 8 | neumann_corrected_tail | -9.37370131e+00 | -9.37484970e+00 | False |
| signed_xxy | 99 | 8 | neumann_corrected_tail | -8.23165720e+00 | -8.23280559e+00 | False |
| symmetry_x_even | 99 | 9 | neumann_corrected_tail | -9.37370131e+00 | -9.37484970e+00 | False |
| signed_xxy | 99 | 9 | neumann_corrected_tail | -8.23165720e+00 | -8.23280559e+00 | False |
| symmetry_x_even | 99 | 10 | neumann_corrected_tail | -9.37370131e+00 | -9.37484970e+00 | False |
| signed_xxy | 99 | 10 | neumann_corrected_tail | -8.23165720e+00 | -8.23280559e+00 | False |
| norm_bound | 101 |  | tail_ignored | -4.53186018e+00 | -4.53300857e+00 | False |
| norm_bound | 101 | 5 | neumann_corrected_tail | -3.01694469e+01 | -3.01705953e+01 | False |
| norm_bound | 101 | 6 | neumann_corrected_tail | -3.01694469e+01 | -3.01705953e+01 | False |
| norm_bound | 101 | 7 | neumann_corrected_tail | -1.74297999e+01 | -1.74309483e+01 | False |
| norm_bound | 101 | 8 | neumann_corrected_tail | -1.74297999e+01 | -1.74309483e+01 | False |
| norm_bound | 101 | 9 | neumann_corrected_tail | -1.74297999e+01 | -1.74309483e+01 | False |
| norm_bound | 101 | 10 | neumann_corrected_tail | -1.74297999e+01 | -1.74309483e+01 | False |
| symmetry_x_even | 101 |  | tail_ignored | 2.02066931e+00 | 2.01952092e+00 | True |
| signed_xxy | 101 |  | tail_ignored | 2.85602178e+00 | 2.85487339e+00 | True |
| symmetry_x_even | 101 | 5 | neumann_corrected_tail | -2.07473056e+01 | -2.07484540e+01 | False |
| signed_xxy | 101 | 5 | neumann_corrected_tail | -1.94961794e+01 | -1.94973278e+01 | False |
| symmetry_x_even | 101 | 6 | neumann_corrected_tail | -2.07473056e+01 | -2.07484540e+01 | False |
| signed_xxy | 101 | 6 | neumann_corrected_tail | -1.94961794e+01 | -1.94973278e+01 | False |
| symmetry_x_even | 101 | 7 | neumann_corrected_tail | -9.37074942e+00 | -9.37189782e+00 | False |
| signed_xxy | 101 | 7 | neumann_corrected_tail | -8.22886363e+00 | -8.23001202e+00 | False |
| symmetry_x_even | 101 | 8 | neumann_corrected_tail | -9.37074942e+00 | -9.37189782e+00 | False |
| signed_xxy | 101 | 8 | neumann_corrected_tail | -8.22886363e+00 | -8.23001202e+00 | False |
| symmetry_x_even | 101 | 9 | neumann_corrected_tail | -9.37074942e+00 | -9.37189782e+00 | False |
| signed_xxy | 101 | 9 | neumann_corrected_tail | -8.22886363e+00 | -8.23001202e+00 | False |
| symmetry_x_even | 101 | 10 | neumann_corrected_tail | -9.37074942e+00 | -9.37189782e+00 | False |
| signed_xxy | 101 | 10 | neumann_corrected_tail | -8.22886363e+00 | -8.23001202e+00 | False |
| norm_bound | 199 |  | tail_ignored | -4.54210700e+00 | -4.54325539e+00 | False |
| norm_bound | 199 | 5 | neumann_corrected_tail | -2.99526200e+01 | -2.99537684e+01 | False |
| norm_bound | 199 | 6 | neumann_corrected_tail | -2.99526200e+01 | -2.99537684e+01 | False |
| norm_bound | 199 | 7 | neumann_corrected_tail | -1.73188814e+01 | -1.73200298e+01 | False |
| norm_bound | 199 | 8 | neumann_corrected_tail | -1.73188814e+01 | -1.73200298e+01 | False |
| norm_bound | 199 | 9 | neumann_corrected_tail | -1.73188814e+01 | -1.73200298e+01 | False |
| norm_bound | 199 | 10 | neumann_corrected_tail | -1.73188814e+01 | -1.73200298e+01 | False |
| symmetry_x_even | 199 |  | tail_ignored | 2.01170876e+00 | 2.01056037e+00 | True |
| signed_xxy | 199 |  | tail_ignored | 2.84720032e+00 | 2.84605193e+00 | True |
| symmetry_x_even | 199 | 5 | neumann_corrected_tail | -2.05745208e+01 | -2.05756692e+01 | False |
| signed_xxy | 199 | 5 | neumann_corrected_tail | -1.93298192e+01 | -1.93309675e+01 | False |
| symmetry_x_even | 199 | 6 | neumann_corrected_tail | -2.05745208e+01 | -2.05756692e+01 | False |
| signed_xxy | 199 | 6 | neumann_corrected_tail | -1.93298192e+01 | -1.93309675e+01 | False |
| symmetry_x_even | 199 | 7 | neumann_corrected_tail | -9.28360668e+00 | -9.28475507e+00 | False |
| signed_xxy | 199 | 7 | neumann_corrected_tail | -8.14638399e+00 | -8.14753238e+00 | False |
| symmetry_x_even | 199 | 8 | neumann_corrected_tail | -9.28360668e+00 | -9.28475507e+00 | False |
| signed_xxy | 199 | 8 | neumann_corrected_tail | -8.14638399e+00 | -8.14753238e+00 | False |
| symmetry_x_even | 199 | 9 | neumann_corrected_tail | -9.28360668e+00 | -9.28475507e+00 | False |
| signed_xxy | 199 | 9 | neumann_corrected_tail | -8.14638399e+00 | -8.14753238e+00 | False |
| symmetry_x_even | 199 | 10 | neumann_corrected_tail | -9.28360668e+00 | -9.28475507e+00 | False |
| signed_xxy | 199 | 10 | neumann_corrected_tail | -8.14638399e+00 | -8.14753238e+00 | False |
| norm_bound | 203 |  | tail_ignored | -4.54284629e+00 | -4.54399468e+00 | False |
| norm_bound | 203 | 5 | neumann_corrected_tail | -2.99506581e+01 | -2.99518065e+01 | False |
| norm_bound | 203 | 6 | neumann_corrected_tail | -2.99506581e+01 | -2.99518065e+01 | False |
| norm_bound | 203 | 7 | neumann_corrected_tail | -1.73178780e+01 | -1.73190264e+01 | False |
| norm_bound | 203 | 8 | neumann_corrected_tail | -1.73178780e+01 | -1.73190264e+01 | False |
| norm_bound | 203 | 9 | neumann_corrected_tail | -1.73178780e+01 | -1.73190264e+01 | False |
| norm_bound | 203 | 10 | neumann_corrected_tail | -1.73178780e+01 | -1.73190264e+01 | False |
| symmetry_x_even | 203 |  | tail_ignored | 2.01106179e+00 | 2.00991340e+00 | True |
| signed_xxy | 203 |  | tail_ignored | 2.84656333e+00 | 2.84541493e+00 | True |
| symmetry_x_even | 203 | 5 | neumann_corrected_tail | -2.05729564e+01 | -2.05741048e+01 | False |
| signed_xxy | 203 | 5 | neumann_corrected_tail | -1.93283128e+01 | -1.93294612e+01 | False |
| symmetry_x_even | 203 | 6 | neumann_corrected_tail | -2.05729564e+01 | -2.05741048e+01 | False |
| signed_xxy | 203 | 6 | neumann_corrected_tail | -1.93283128e+01 | -1.93294612e+01 | False |
| symmetry_x_even | 203 | 7 | neumann_corrected_tail | -9.28281785e+00 | -9.28396624e+00 | False |
| signed_xxy | 203 | 7 | neumann_corrected_tail | -8.14563728e+00 | -8.14678567e+00 | False |
| symmetry_x_even | 203 | 8 | neumann_corrected_tail | -9.28281785e+00 | -9.28396624e+00 | False |
| signed_xxy | 203 | 8 | neumann_corrected_tail | -8.14563728e+00 | -8.14678567e+00 | False |
| symmetry_x_even | 203 | 9 | neumann_corrected_tail | -9.28281785e+00 | -9.28396624e+00 | False |
| signed_xxy | 203 | 9 | neumann_corrected_tail | -8.14563728e+00 | -8.14678567e+00 | False |
| symmetry_x_even | 203 | 10 | neumann_corrected_tail | -9.28281785e+00 | -9.28396624e+00 | False |
| signed_xxy | 203 | 10 | neumann_corrected_tail | -8.14563728e+00 | -8.14678567e+00 | False |

Parseval decomposition diagnostic:

- smallest absolute combined residual in the grid is `1.26034982e+01` at `J=203`, `M=7`.
- Full values are written to `parseval_decomposition.csv`.
