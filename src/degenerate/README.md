# Verified thin-triangle Rayleigh quotient

This directory verifies the nearly degenerate part of the upper
area--perimeter inequality without evaluating Airy functions.

For the triangle with vertices \((-1,0),(1,0),(s,t)\), set
\(a=t^{2/3}\) and use \(z=ax+s\).  Problem 2 in the thin-triangle model is
then posed on the fixed interval \([-1,1]\).  Its Rayleigh quotient is
evaluated with

\[
 \phi(z)=(1-z^2)\sum_{j=1}^{N}q_jP_{j-1}(z),
\]

where \(P_j\) is the Legendre polynomial of degree \(j\).  The vector
\(q\) is computed at the midpoint of an \(s\)-cell and is frozen over the
whole cell.  For replay, the artifact also stores directed hexadecimal
endpoints for the resulting monomial coefficients of
\(p(z)=\sum_jq_jP_{j-1}(z)\).  Those monomial intervals, rather than a
fresh eigensolve or a rebuilt Legendre table, are authoritative in the
certificate validator.

This basis makes every integral algebraic.  On the two sides of the
apex,

\[
 \frac{\phi(z)^2}{(1+z)^2}
  =(1-z)^2\left(\sum_jq_jP_{j-1}(z)\right)^2,
 \qquad
 \frac{\phi(z)^2}{(1-z)^2}
  =(1+z)^2\left(\sum_jq_jP_{j-1}(z)\right)^2.
\]

Thus the mass, kinetic energy, and potential energy are evaluated by
polynomial antiderivatives.  `muhat1_rayleigh_cover.m` uses INTLAB
outward rounding for the coefficients and evaluates the resulting
polynomial in \(s\) in centered form.

## Direct use in the upper conjecture

Let \(R(s)\) be the certified Rayleigh upper bound at \(t_0=0.38\).
Monotonicity of the one-dimensional model and the min--max principle give

\[
 \lambda_1(\triangle^{(x,y)})
 \le \frac{\pi^2}{y^2}
      +\frac{4R(2x-1)}{(2y)^{4/3}},
 \qquad 0<y\le0.19.
\]

It is unnecessary to replace \(R(s)\) by one uniform constant.  After
multiplication by \(12y\), the desired upper inequality follows if

\[
 6\,2^{2/3}R(s)y^{2/3}+6\pi^2
 -2\pi^2\!\left(
 1+\sqrt{x^2+y^2}+\sqrt{(1-x)^2+y^2}
 \right)^2
 -4\sqrt3\,\pi^2y<0,
 \qquad x=\frac{1+s}{2}.
\]

The production run evaluates this expression simultaneously with the
Rayleigh quotient.  This retains the compensating dependence of the
perimeter on \(s\).  The closed computational hull includes \(s=1\) and
\(y=0\); this is harmless and covers the theorem's half-open domain
\(0\le s<1\), \(0<y\le0.19\).

Every exponent used by the interval path is constructed as an INTLAB
enclosure of the exact rational number (for example
`intval('2')/intval('3')`).  Likewise, the decimal parameters `0.38` and
`0.19` and the endpoints of the uniform \(s\)-partition enter INTLAB as
decimal strings, without first being rounded to binary64.

## Reproduction

Set `INTLAB_ROOT` and an output root outside the clean Git checkout, then
run from the repository root:

```sh
export INTLAB_ROOT=/path/to/Intlab_V12
export VER10_CERTIFICATE_ROOT=/path/to/certificate-output
```

```matlab
run scripts_run/run_muhat1_rayleigh_hpc.m
```

The driver refuses to label a production result complete unless the source
tree is clean and unchanged throughout the run.  It writes a manifest,
lossless hexadecimal frozen coefficients, per-\(s\)-cell Rayleigh bounds,
all per-\((s,r)\)-cell upper-functional bounds, and SHA-256 checksums below
`$VER10_CERTIFICATE_ROOT/muhat1_rayleigh_t038/`.  Before returning, it
replays the proof from the frozen trials and validates every hash,
cell ID, endpoint, mass, and strict upper bound with INTLAB.  A copied
certificate can be checked independently with
`validate_muhat1_rayleigh_certificate`.
