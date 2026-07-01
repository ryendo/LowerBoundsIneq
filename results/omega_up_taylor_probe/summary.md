# Omega Up Taylor Probe

This is an exploratory double-precision Taylor spectral probe. It does not modify the production Omega_up or Omega_mid verification code.

The spectral data use analytic equilateral triangle eigenfunctions, but spectral inner products and Rxx cell lower bounds are evaluated numerically in double precision. Therefore the rigorous flag is `exploratory_double`.

Known comparison points:

- Current cellwise eigenspace enclosure succeeds only around distance `1e-5`.
- Idealized fixed-equilateral-coefficient model suggests `1e-2` to `1e-1` scale.
- This Taylor probe tests whether the third-derivative Taylor certificate can bridge that gap.

| model | J | tail mode | cell shape | max rx | max ry | distance | Jxx1 lower | Jxx2 lower | rigorous flag |
|---|---:|---|---|---:|---:|---:|---:|---:|---|
| taylor_spectral_probe | 30 | tail_ignored | x_only | 4.48081823e-03 | 0.00000000e+00 | 4.48081823e-03 | 1.48524942e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 30 | tail_ignored | y_only | 0.00000000e+00 | 3.40823079e-03 | 3.40823079e-03 | 1.55790964e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 30 | tail_ignored | omega_up_type | 3.96025396e-03 | 3.96025396e-04 | 3.98000597e-03 | 1.50399992e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 30 | tail_ignored | square | 1.93590103e-03 | 1.93590103e-03 | 2.73777749e-03 | 1.55122141e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 30 | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.16096132e+01 | -5.16111829e+01 | exploratory_double |
| taylor_spectral_probe | 30 | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.16096132e+01 | -5.16111829e+01 | exploratory_double |
| taylor_spectral_probe | 30 | dirichlet_parseval_tail | omega_up_type | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.16096132e+01 | -5.16111829e+01 | exploratory_double |
| taylor_spectral_probe | 30 | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.16096132e+01 | -5.16111829e+01 | exploratory_double |
| taylor_spectral_probe | 50 | tail_ignored | x_only | 4.43206314e-03 | 0.00000000e+00 | 4.43206314e-03 | 1.48707671e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 50 | tail_ignored | y_only | 0.00000000e+00 | 3.37157240e-03 | 3.37157240e-03 | 1.55821194e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 50 | tail_ignored | omega_up_type | 3.91721963e-03 | 3.91721963e-04 | 3.93675701e-03 | 1.50542579e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 50 | tail_ignored | square | 1.91497315e-03 | 1.91497315e-03 | 2.70818100e-03 | 1.55164719e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 50 | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.06703179e+01 | -5.06718875e+01 | exploratory_double |
| taylor_spectral_probe | 50 | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.06703179e+01 | -5.06718875e+01 | exploratory_double |
| taylor_spectral_probe | 50 | dirichlet_parseval_tail | omega_up_type | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.06703179e+01 | -5.06718875e+01 | exploratory_double |
| taylor_spectral_probe | 50 | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -5.06703179e+01 | -5.06718875e+01 | exploratory_double |
| taylor_spectral_probe | 100 | tail_ignored | x_only | 4.39771885e-03 | 0.00000000e+00 | 4.39771885e-03 | 1.48835190e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 100 | tail_ignored | y_only | 0.00000000e+00 | 3.34571527e-03 | 3.34571527e-03 | 1.55842286e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 100 | tail_ignored | omega_up_type | 3.88690061e-03 | 3.88690061e-04 | 3.90628677e-03 | 1.50642096e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 100 | tail_ignored | square | 1.90022006e-03 | 1.90022006e-03 | 2.68731698e-03 | 1.55194437e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 100 | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98508756e+01 | -4.98524453e+01 | exploratory_double |
| taylor_spectral_probe | 100 | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98508756e+01 | -4.98524453e+01 | exploratory_double |
| taylor_spectral_probe | 100 | dirichlet_parseval_tail | omega_up_type | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98508756e+01 | -4.98524453e+01 | exploratory_double |
| taylor_spectral_probe | 100 | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.98508756e+01 | -4.98524453e+01 | exploratory_double |
| taylor_spectral_probe | 200 | tail_ignored | x_only | 4.38185316e-03 | 0.00000000e+00 | 4.38185316e-03 | 1.48893763e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 200 | tail_ignored | y_only | 0.00000000e+00 | 3.33377416e-03 | 3.33377416e-03 | 1.55851962e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 200 | tail_ignored | omega_up_type | 3.87289495e-03 | 3.87289495e-04 | 3.89221126e-03 | 1.50687805e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 200 | tail_ignored | square | 1.89340597e-03 | 1.89340597e-03 | 2.67768040e-03 | 1.55208079e-03 | 1.77635684e-15 | exploratory_double |
| taylor_spectral_probe | 200 | dirichlet_parseval_tail | x_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94182383e+01 | -4.94198080e+01 | exploratory_double |
| taylor_spectral_probe | 200 | dirichlet_parseval_tail | y_only | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94182383e+01 | -4.94198080e+01 | exploratory_double |
| taylor_spectral_probe | 200 | dirichlet_parseval_tail | omega_up_type | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94182383e+01 | -4.94198080e+01 | exploratory_double |
| taylor_spectral_probe | 200 | dirichlet_parseval_tail | square | 0.00000000e+00 | 0.00000000e+00 | 0.00000000e+00 | -4.94182383e+01 | -4.94198080e+01 | exploratory_double |

Most important Omega_up-type result:

- best `rx = 3.96025396e-03`, `ry = 3.96025396e-04` at `J = 30`, tail mode `tail_ignored`.
- target `rx >= 1e-2`, `ry = 1e-3`: not reached.
