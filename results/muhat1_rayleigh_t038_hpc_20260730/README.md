# HPC/INTLAB certificate at \(t_0=0.38\)

The production calculation used MATLAB R2023b Update 5 and INTLAB on
`liulab-hpc2023`.  It covered \(0\le s\le1\) with 200 cells, used ten
Legendre bubble trial functions, and subdivided
\(0\le y^{1/3}\le0.19^{1/3}\) into 40 intervals for each \(s\)-cell.

The certified extrema are:

- maximum Rayleigh upper endpoint: `14.049354897039173`;
- cell attaining that recorded upper endpoint: `[0.995,1]`;
- maximum scaled upper-conjecture endpoint: `-1.5240735212108216`;
- cell attaining that recorded scaled endpoint: `[0.65,0.655]`.

The negative scaled endpoint certifies the upper area--perimeter
inequality for every \(0<y\le0.19\) and \(0\le s<1\).  The computation
does not certify the proposed uniform estimate \(R(s)\le11.5\); instead,
it retains the \(s\)-dependence of the Rayleigh quotient and evaluates it
together with the perimeter.

`summary.json` records the run configuration and extrema.  `cells.csv`
contains the enclosure on every \(s\)-cell.
