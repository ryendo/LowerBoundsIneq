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
- `symmetry_first_order`: finite `lambda_xx(0,0)` with the `rx M_xxx` penalty removed by equilateral x-symmetry, keeping `-ry M_xxy`.
- `signed_xxy_first_order`: finite `lambda_xx(0,0)` plus the signed affine `lambda_xxy t` minimum on `t in [-ry,0]`.
- `dirichlet_parseval_tail`: the intentionally crude Dirichlet-only Parseval tail diagnostic; it is not used by the two symmetry-improved models.

Finite signed derivative diagnostics:

| J | cluster complete | complete cutoffs | lambda_x | lambda_xx finite | lambda_xxx finite | lambda_xxy finite | M_xxx/abs | M_xxy/abs |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 30 | True | 30/30 | -4.88786096e-13 | 3.63934075e+01 | 7.20151797e-13 | -7.50212457e+01 | 2.10126127e+15 | 2.57191748e+01 |
| 50 | False | 49/51 | -4.88786096e-13 | 3.63242373e+01 | 7.68189976e-13 | -7.51138036e+01 | 1.97121675e+15 | 2.57023867e+01 |
| 100 | False | 99/101 | -4.88786096e-13 | 3.62754001e+01 | 6.93616649e-13 | -7.51588564e+01 | 2.18419098e+15 | 2.56974916e+01 |
| 200 | False | 199/203 | -4.88786096e-13 | 3.62550041e+01 | 6.35004396e-13 | -7.51707171e+01 | 2.38626478e+15 | 2.56977149e+01 |

The finite `lambda_xxx` values are numerically close to zero at the `1e-6` level.

Basis orthogonality diagnostics:

| bc | J | mass max error | stiffness max error | stiffness relative error |
|---|---:|---:|---:|---:|
| D | 30 | 2.81996648e-14 | 3.76303433e-11 | 3.51585638e-14 |
| D | 50 | 2.81996648e-14 | 7.68523023e-11 | 4.51552166e-14 |
| D | 100 | 3.06421555e-14 | 1.74622983e-10 | 5.26577598e-14 |
| D | 200 | 3.24185123e-14 | 3.36513040e-10 | 5.26894051e-14 |
| N | 30 | 4.90533275e-14 | 3.71755959e-11 | 5.43270488e-14 |
| N | 50 | 4.90533275e-14 | 5.83213478e-11 | 4.96107200e-14 |
| N | 100 | 4.90533275e-14 | 1.06410880e-10 | 4.12564145e-14 |
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

Cell-size bisection results:

| model | J | tail mode | cell shape | max rx | max ry | distance | Jxx1 lower | Jxx2 lower | rigorous flag |
|---|---:|---|---|---:|---:|---:|---:|---:|---|
| norm_bound | 30 | tail_ignored | x_only | 4.48081823e-03 | 0.00000000e+00 | 4.48081823e-03 | 1.48524942e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | tail_ignored | y_only | 0.00000000e+00 | 3.40823079e-03 | 3.40823079e-03 | 1.55790964e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | tail_ignored | omega_up_type | 3.96025396e-03 | 3.96025396e-04 | 3.98000597e-03 | 1.50399992e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | tail_ignored | square | 1.93590103e-03 | 1.93590103e-03 | 2.73777749e-03 | 1.55122141e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 30 | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.16096132e+01 | -5.16111829e+01 | exploratory_double |
| norm_bound | 30 | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.16096132e+01 | -5.16111829e+01 | exploratory_double |
| norm_bound | 30 | dirichlet_parseval_tail | omega_up_type | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.16096132e+01 | -5.16111829e+01 | exploratory_double |
| norm_bound | 30 | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.16096132e+01 | -5.16111829e+01 | exploratory_double |
| symmetry_first_order | 30 | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.40091519e+00 | 2.56929486e+00 | exploratory_double |
| symmetry_first_order | 30 | tail_ignored | y_only | 0.00000000e+00 | 3.40823079e-03 | 3.40823079e-03 | 1.55790964e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_first_order | 30 | tail_ignored | omega_up_type | 3.39044348e-02 | 3.39044348e-03 | 3.40735353e-02 | 1.77635684e-15 | 3.30583722e-03 | exploratory_double |
| symmetry_first_order | 30 | tail_ignored | square | 3.40808954e-03 | 3.40808954e-03 | 4.81976645e-03 | 1.50878081e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy_first_order | 30 | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.40091519e+00 | 2.56929486e+00 | exploratory_double |
| signed_xxy_first_order | 30 | tail_ignored | y_only | 0.00000000e+00 | 9.88637327e-02 | 9.88637327e-02 | 1.77635684e-15 | 1.52399842e-02 | exploratory_double |
| signed_xxy_first_order | 30 | tail_ignored | omega_up_type | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.75747456e+00 | 1.93281022e+00 | exploratory_double |
| signed_xxy_first_order | 30 | tail_ignored | square | 9.24585593e-02 | 9.24585593e-02 | 1.30756149e-01 | 8.88178420e-15 | 5.52099788e-02 | exploratory_double |
| norm_bound | 50 | tail_ignored | x_only | 4.43206314e-03 | 0.00000000e+00 | 4.43206314e-03 | 1.48707671e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 50 | tail_ignored | y_only | 0.00000000e+00 | 3.37157240e-03 | 3.37157240e-03 | 1.55821194e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 50 | tail_ignored | omega_up_type | 3.91721963e-03 | 3.91721963e-04 | 3.93675701e-03 | 1.50542579e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 50 | tail_ignored | square | 1.91497315e-03 | 1.91497315e-03 | 2.70818100e-03 | 1.55164719e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 50 | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.06703179e+01 | -5.06718875e+01 | exploratory_double |
| norm_bound | 50 | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.06703179e+01 | -5.06718875e+01 | exploratory_double |
| norm_bound | 50 | dirichlet_parseval_tail | omega_up_type | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.06703179e+01 | -5.06718875e+01 | exploratory_double |
| norm_bound | 50 | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.06703179e+01 | -5.06718875e+01 | exploratory_double |
| symmetry_first_order | 50 | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.37096364e+00 | 2.53934331e+00 | exploratory_double |
| symmetry_first_order | 50 | tail_ignored | y_only | 0.00000000e+00 | 3.37157240e-03 | 3.37157240e-03 | 1.55821194e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_first_order | 50 | tail_ignored | omega_up_type | 3.35421515e-02 | 3.35421515e-03 | 3.37094451e-02 | 1.77635684e-15 | 3.20181406e-03 | exploratory_double |
| symmetry_first_order | 50 | tail_ignored | square | 3.37143429e-03 | 3.37143429e-03 | 4.76792809e-03 | 1.51013726e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy_first_order | 50 | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.37096364e+00 | 2.53934331e+00 | exploratory_double |
| signed_xxy_first_order | 50 | tail_ignored | y_only | 0.00000000e+00 | 9.79798235e-02 | 9.79798235e-02 | 5.32907052e-15 | 1.49075816e-02 | exploratory_double |
| signed_xxy_first_order | 50 | tail_ignored | omega_up_type | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.72821471e+00 | 1.90355037e+00 | exploratory_double |
| signed_xxy_first_order | 50 | tail_ignored | square | 9.16946123e-02 | 9.16946123e-02 | 1.29675764e-01 | 1.77635684e-15 | 5.41927448e-02 | exploratory_double |
| norm_bound | 100 | tail_ignored | x_only | 4.39771885e-03 | 0.00000000e+00 | 4.39771885e-03 | 1.48835190e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 100 | tail_ignored | y_only | 0.00000000e+00 | 3.34571527e-03 | 3.34571527e-03 | 1.55842286e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 100 | tail_ignored | omega_up_type | 3.88690061e-03 | 3.88690061e-04 | 3.90628677e-03 | 1.50642096e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 100 | tail_ignored | square | 1.90022006e-03 | 1.90022006e-03 | 2.68731698e-03 | 1.55194437e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 100 | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98508756e+01 | -4.98524453e+01 | exploratory_double |
| norm_bound | 100 | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98508756e+01 | -4.98524453e+01 | exploratory_double |
| norm_bound | 100 | dirichlet_parseval_tail | omega_up_type | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98508756e+01 | -4.98524453e+01 | exploratory_double |
| norm_bound | 100 | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98508756e+01 | -4.98524453e+01 | exploratory_double |
| symmetry_first_order | 100 | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34981648e+00 | 2.51819616e+00 | exploratory_double |
| symmetry_first_order | 100 | tail_ignored | y_only | 0.00000000e+00 | 3.34571527e-03 | 3.34571527e-03 | 1.55842286e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_first_order | 100 | tail_ignored | omega_up_type | 3.32865834e-02 | 3.32865834e-03 | 3.34526023e-02 | 1.77635684e-15 | 3.12911161e-03 | exploratory_double |
| symmetry_first_order | 100 | tail_ignored | square | 3.34557935e-03 | 3.34557935e-03 | 4.73136369e-03 | 1.51108481e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy_first_order | 100 | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34981648e+00 | 2.51819616e+00 | exploratory_double |
| signed_xxy_first_order | 100 | tail_ignored | y_only | 0.00000000e+00 | 9.73544387e-02 | 9.73544387e-02 | 1.77635684e-15 | 1.46747705e-02 | exploratory_double |
| signed_xxy_first_order | 100 | tail_ignored | omega_up_type | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.70755593e+00 | 1.88289159e+00 | exploratory_double |
| signed_xxy_first_order | 100 | tail_ignored | square | 9.11533838e-02 | 9.11533838e-02 | 1.28910352e-01 | 5.32907052e-15 | 5.34785288e-02 | exploratory_double |
| norm_bound | 200 | tail_ignored | x_only | 4.38339757e-03 | 0.00000000e+00 | 4.38339757e-03 | 1.48888071e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 200 | tail_ignored | y_only | 0.00000000e+00 | 3.33493661e-03 | 3.33493661e-03 | 1.55851022e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 200 | tail_ignored | omega_up_type | 3.87425831e-03 | 3.87425831e-04 | 3.89358141e-03 | 1.50683363e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 200 | tail_ignored | square | 1.89406929e-03 | 1.89406929e-03 | 2.67861848e-03 | 1.55206753e-03 | 3.55271368e-15 | exploratory_double |
| norm_bound | 200 | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94451199e+01 | -4.94466896e+01 | exploratory_double |
| norm_bound | 200 | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94451199e+01 | -4.94466896e+01 | exploratory_double |
| norm_bound | 200 | dirichlet_parseval_tail | omega_up_type | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94451199e+01 | -4.94466896e+01 | exploratory_double |
| norm_bound | 200 | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94451199e+01 | -4.94466896e+01 | exploratory_double |
| symmetry_first_order | 200 | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34098477e+00 | 2.50936444e+00 | exploratory_double |
| symmetry_first_order | 200 | tail_ignored | y_only | 0.00000000e+00 | 3.33493661e-03 | 3.33493661e-03 | 1.55851022e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_first_order | 200 | tail_ignored | omega_up_type | 3.31800406e-02 | 3.31800406e-03 | 3.33455282e-02 | 1.77635684e-15 | 3.09896892e-03 | exploratory_double |
| symmetry_first_order | 200 | tail_ignored | square | 3.33480159e-03 | 3.33480159e-03 | 4.71612164e-03 | 1.51147755e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy_first_order | 200 | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34098477e+00 | 2.50936444e+00 | exploratory_double |
| signed_xxy_first_order | 200 | tail_ignored | y_only | 0.00000000e+00 | 9.70929376e-02 | 9.70929376e-02 | 5.32907052e-15 | 1.45780025e-02 | exploratory_double |
| signed_xxy_first_order | 200 | tail_ignored | omega_up_type | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.69892817e+00 | 1.87426383e+00 | exploratory_double |
| signed_xxy_first_order | 200 | tail_ignored | square | 9.09268955e-02 | 9.09268955e-02 | 1.28590049e-01 | 1.77635684e-15 | 5.31812353e-02 | exploratory_double |

Most important Omega_up-type result:

- best `rx = 2.00000000e-01`, `ry = 2.00000000e-02` with model `signed_xxy_first_order`, `J = 30`, tail mode `tail_ignored`.
- target `rx >= 1e-2`, `ry = 1e-3`: reached.

Omega_up-type best row by model:

| model | best rx | best ry | J | tail mode | Jxx2 lower |
|---|---:|---:|---:|---|---:|
| norm_bound | 3.96025396e-03 | 3.96025396e-04 | 30 | tail_ignored | 1.77635684e-15 |
| signed_xxy_first_order | 2.00000000e-01 | 2.00000000e-02 | 30 | tail_ignored | 1.93281022e+00 |
| symmetry_first_order | 3.39044348e-02 | 3.39044348e-03 | 30 | tail_ignored | 3.30583722e-03 |

Target diagnostic at `rx=1e-2`, `ry=1e-3`:

| model | J | tail mode | Jxx1 lower | Jxx2 lower | success |
|---|---:|---|---:|---:|---|
| norm_bound | 30 | tail_ignored | -4.47236440e+00 | -4.47351279e+00 | False |
| norm_bound | 30 | dirichlet_parseval_tail | -6.61635968e+01 | -6.61647452e+01 | False |
| symmetry_first_order | 30 | tail_ignored | 2.07253492e+00 | 2.07138652e+00 | True |
| signed_xxy_first_order | 30 | tail_ignored | 2.90706148e+00 | 2.90591309e+00 | True |
| norm_bound | 50 | tail_ignored | -4.50727179e+00 | -4.50842018e+00 | False |
| norm_bound | 50 | dirichlet_parseval_tail | -6.49798024e+01 | -6.49809508e+01 | False |
| symmetry_first_order | 50 | tail_ignored | 2.04213376e+00 | 2.04098536e+00 | True |
| signed_xxy_first_order | 50 | tail_ignored | 2.87714452e+00 | 2.87599613e+00 | True |
| norm_bound | 100 | tail_ignored | -4.53186018e+00 | -4.53300857e+00 | False |
| norm_bound | 100 | dirichlet_parseval_tail | -6.39498026e+01 | -6.39509510e+01 | False |
| symmetry_first_order | 100 | tail_ignored | 2.02066931e+00 | 2.01952092e+00 | True |
| signed_xxy_first_order | 100 | tail_ignored | 2.85602178e+00 | 2.85487339e+00 | True |
| norm_bound | 200 | tail_ignored | -4.54210707e+00 | -4.54325546e+00 | False |
| norm_bound | 200 | dirichlet_parseval_tail | -6.34407396e+01 | -6.34418880e+01 | False |
| symmetry_first_order | 200 | tail_ignored | 2.01170871e+00 | 2.01056032e+00 | True |
| signed_xxy_first_order | 200 | tail_ignored | 2.84720026e+00 | 2.84605187e+00 | True |

Parseval decomposition diagnostic:

- smallest absolute combined residual in the grid is `2.11038277e+00` at `J=200`, `M=200`.
- Full values are written to `parseval_decomposition.csv`.
