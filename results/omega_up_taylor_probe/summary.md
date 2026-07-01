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

Finite third derivative diagnostics:

| J | cluster complete | lambda_xx finite | lambda_xxx finite | lambda_xxy finite | M_xxx/abs | M_xxy/abs |
|---:|---|---:|---:|---:|---:|---:|
| 30 | True | 3.63934075e+01 | 7.30809938e-13 | -7.50212457e+01 | 2.07061645e+15 | 2.57191748e+01 |
| 50 | False | 3.63242373e+01 | 7.68189976e-13 | -7.51138036e+01 | 1.97121675e+15 | 2.57023867e+01 |
| 100 | False | 3.62754001e+01 | 6.88287578e-13 | -7.51588564e+01 | 2.20110209e+15 | 2.56974916e+01 |
| 200 | False | 3.62528033e+01 | 6.19017693e-13 | -7.51710501e+01 | 2.44794379e+15 | 2.56980595e+01 |

The finite `lambda_xxx` values are numerically close to zero at the `1e-6` level.

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
| norm_bound | 200 | tail_ignored | x_only | 4.38185316e-03 | 0.00000000e+00 | 4.38185316e-03 | 1.48893763e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 200 | tail_ignored | y_only | 0.00000000e+00 | 3.33377416e-03 | 3.33377416e-03 | 1.55851962e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 200 | tail_ignored | omega_up_type | 3.87289495e-03 | 3.87289495e-04 | 3.89221126e-03 | 1.50687805e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 200 | tail_ignored | square | 1.89340597e-03 | 1.89340597e-03 | 2.67768040e-03 | 1.55208079e-03 | 1.77635684e-15 | exploratory_double |
| norm_bound | 200 | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94182383e+01 | -4.94198080e+01 | exploratory_double |
| norm_bound | 200 | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94182383e+01 | -4.94198080e+01 | exploratory_double |
| norm_bound | 200 | dirichlet_parseval_tail | omega_up_type | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94182383e+01 | -4.94198080e+01 | exploratory_double |
| norm_bound | 200 | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94182383e+01 | -4.94198080e+01 | exploratory_double |
| symmetry_first_order | 200 | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34003180e+00 | 2.50841147e+00 | exploratory_double |
| symmetry_first_order | 200 | tail_ignored | y_only | 0.00000000e+00 | 3.33377416e-03 | 3.33377416e-03 | 1.55851962e-03 | 1.77635684e-15 | exploratory_double |
| symmetry_first_order | 200 | tail_ignored | omega_up_type | 3.31685500e-02 | 3.31685500e-03 | 3.33339802e-02 | 1.77635684e-15 | 3.09572387e-03 | exploratory_double |
| symmetry_first_order | 200 | tail_ignored | square | 3.33363924e-03 | 3.33363924e-03 | 4.71447783e-03 | 1.51151982e-03 | 1.77635684e-15 | exploratory_double |
| signed_xxy_first_order | 200 | tail_ignored | x_only | 2.00000000e-01 | 0.00000000e+00 | 2.00000000e-01 | 2.34003180e+00 | 2.50841147e+00 | exploratory_double |
| signed_xxy_first_order | 200 | tail_ignored | y_only | 0.00000000e+00 | 9.70647096e-02 | 9.70647096e-02 | 1.77635684e-15 | 1.45675772e-02 | exploratory_double |
| signed_xxy_first_order | 200 | tail_ignored | omega_up_type | 2.00000000e-01 | 2.00000000e-02 | 2.00997512e-01 | 1.69799721e+00 | 1.87333287e+00 | exploratory_double |
| signed_xxy_first_order | 200 | tail_ignored | square | 9.09024407e-02 | 9.09024407e-02 | 1.28555464e-01 | 1.77635684e-15 | 5.31491913e-02 | exploratory_double |

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
| norm_bound | 200 | tail_ignored | -4.54321178e+00 | -4.54436017e+00 | False |
| norm_bound | 200 | dirichlet_parseval_tail | -6.34070362e+01 | -6.34081846e+01 | False |
| symmetry_first_order | 200 | tail_ignored | 2.01074193e+00 | 2.00959354e+00 | True |
| signed_xxy_first_order | 200 | tail_ignored | 2.84624840e+00 | 2.84510000e+00 | True |

Parseval decomposition diagnostic:

- smallest absolute combined residual in the grid is `2.09969695e+00` at `J=200`, `M=200`.
- Full values are written to `parseval_decomposition.csv`.
