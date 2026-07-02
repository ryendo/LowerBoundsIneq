# Omega_up Taylor Descent Probe

Exploratory double-precision diagnostic; not a proof certificate.

- finite basis modes: `360`
- quadrature order: `220`
- finite `lambda_xx(p0)`: `36.2477175625`
- finite signed y slope `-lambda_xxy(p0)`: `75.1791299048`

## Coupled Signed Taylor Model

| region | rx | ry | y_inf | min J1 | min J2 | min J1 point |
|---|---:|---:|---:|---:|---:|---|
| old_eps_0p122 | 0.244 | 0.122 | 0.744025403784 | 1.24017000e+00 | 1.58612012e+00 | (0.744, 0.7440254) |
| coupled_linear_eps_0p152 | 0.304131295384 | 0.152065647692 | 0.713959756092 | 1.12586768e-03 | 5.83264306e-01 | (0.8041313, 0.71395976) |
| mid_top_narrow | 0.00387334601383 | 0.121382060046 | 0.744643343738 | 2.66510740e+00 | 2.69027717e+00 | (0.50387335, 0.74464334) |
| signed_prev_y_0p8472 | 0.1878944815 | 0.0188254037844 | 0.8472 | 2.39678660e+00 | 2.54506620e+00 | (0.68789448, 0.8660254) |

## Direct Finite Residual Check

| region | min residual | residual point | min direct J1 | min direct J1 point |
|---|---:|---|---:|---|
| old_eps_0p122 | 0.00000000e+00 | (0.5, 0.8660254) | 2.71566483e+00 | (0.744, 0.8660254) |
| coupled_linear_eps_0p152 | 0.00000000e+00 | (0.5, 0.8660254) | 2.62933235e+00 | (0.8041313, 0.8660254) |
| mid_top_narrow | 0.00000000e+00 | (0.5, 0.8660254) | 2.87468960e+00 | (0.50387335, 0.8660254) |
| signed_prev_y_0p8472 | 0.00000000e+00 | (0.5, 0.8660254) | 2.78007641e+00 | (0.68789448, 0.8660254) |
