# Legacy HPC/INTLAB output at \(t_0=0.38\) (not a proof certificate)

This directory is retained only as historical output.  Its interval run
formed the exponents `1/3` and `2/3` as rounded binary64 numbers before
passing them to INTLAB, and it did not retain the frozen trial vectors or a
replayable per-\((s,r)\)-cell artifact.  Consequently, neither
`summary.json` nor the old `complete_certificate` claim in a downstream
manifest may be used as a formal certificate for the exact mathematical
statement.

The corrected production driver is
`scripts_run/run_muhat1_rayleigh_hpc.m`.  It uses exact rational interval
exponents and produces a fail-closed, replay-validated artifact.  Only a
new artifact accepted by `validate_muhat1_rayleigh_certificate` supersedes
this directory.

The production calculation used MATLAB R2023b Update 5 and INTLAB on
`liulab-hpc2023`.  It covered \(0\le s\le1\) with 200 cells, used ten
Legendre bubble trial functions, and subdivided
\(0\le y^{1/3}\le0.19^{1/3}\) into 40 intervals for each \(s\)-cell.

The certified extrema are:

- maximum Rayleigh upper endpoint: `14.049354897039173`;
- cell attaining that recorded upper endpoint: `[0.995,1]`;
- maximum scaled upper-conjecture endpoint: `-1.5240735212108216`;
- cell attaining that recorded scaled endpoint: `[0.65,0.655]`.

The negative scaled endpoint is useful regression evidence but, for the
reason above, does not by itself certify the upper area--perimeter
inequality.  The calculation also did not establish the proposed uniform
estimate \(R(s)\le11.5\); it retained the \(s\)-dependence of the quotient
and evaluated it together with the perimeter.

`summary.json` records the run configuration and extrema.  `cells.csv`
contains the enclosure on every \(s\)-cell.
