# Omega_mid adaptive cell builder (correct vertex) — preliminary (float) pass

Rebuilds the Omega_mid covering for the `J1,J2 >= 0` proof **from scratch**, using
the DatabaseTriangle **address** quadtree and evaluating J at the **correct
vertex** (Lemma 4.9):

```
J_k(T^p) >= B_k( p_{i,j} ; lambda_1_lb(cell) )      for all p in the cell,
   geometry (|T|,|dT|)  at the INNER corner  p_{i,j}     = (x_inf, theta_inf),
   eigenvalue lower bnd  valid over the cell (<= lambda_1 at the OUTER corner
                          p_{i+1,j+1} = (x_sup, theta_sup)).
```

The two vertices are different on purpose: the geometry term is minimised at the
inner corner (B is increasing in x and theta), the eigenvalue is smallest at the
outer corner (largest triangle).

## Coverage region (matches the paper's three-region split)

Omega is tiled by `Omega_up u Omega_mid u Omega_down`. This builder covers the
**horizontal slab**

```
  Y_FLOOR = eps_down            <= y <=   Y_TOP = sqrt3/2 - eps_up
       0.04                                       0.74403...
```

inside Omega, with `eps_up = 0.122`, `eps_down = 0.04` (paper eq. 43).

* Below `y = eps_down`  -> Omega_down (handled analytically, Lemma 4.13/4.14).
* Above `y = sqrt3/2 - eps_up` -> Omega_up strip (handled by the convexity argument).
* For `y <= Y_TOP`, every point is automatically outside the Omega_up disk
  (`dist(.,p0) >= sqrt3/2 - y >= eps_up`), so the paper's disk constraint holds.

Using the **strip floor** (a horizontal line) as the ceiling instead of the disk
arc is what keeps the quadtree from blowing up: on the slab `min J ~ J(1/2,Y_TOP)`
is bounded away from 0, so cells converge at a finite depth. Equivalently, the
small-J layer next to the equilateral point is delegated to Omega_up. **B is only
ever evaluated at a mid-cell inner corner with `y >= eps_down`, never at a
degenerate triangle.**

## Pipeline

1. **Iteration 1 — reuse** (`reuse_db_iteration1.py`, no FEM): take each leaf cell
   of `ryendo/DatabaseTriangle` (`database_best_cover.csv`) carrying a *rigorous*
   `lambda_1_lb`, and test `B_k(p_{i,j}; lambda_1_lb) > 0`. Passing cells are
   verified essentially for free.
2. **Iteration 2+ — refine** (`build_cells_integrated.py`): every failing cell is
   split into its 4 address-children (`1`=RU, `2`=LU, `3`=LD, `4`=RD); the
   children's `lambda_1` is estimated by a fast vectorised float FEM
   (`omega_mid_lambda1.py`, aspect-aware: Richardson near the equilateral where J
   is small, one coarse solve for thin triangles where J is large). Recurse until
   both J1,J2 pass (leaf). Boundary (straddle) cells are capped at address length
   `STRADDLE_MAXLEN` and emitted separately.
3. **Rigorous pass** (later, MATLAB/INTLAB): feed the produced `cell_def` to the
   *vertex-corrected* `LowerBoundsIneq` (`verify_J_positive` -> CR first, escalate
   to Lehmann-Goerisch only if CR fails; LG directly near the equilateral). It
   recomputes every leaf rigorously. The float pass only decides the cell
   structure; it is never used as proof.

## Result of the preliminary build (this repo's parameters)

```
seeded 9973 DatabaseTriangle leaves (iteration 1, rigorous lambda reused -> 2678 verified)
verified leaf cells    : 161,347   (DB-reused 2,678 ; FEM-refined 158,669)
boundary slivers       :  24,014   (straddle, address length <= 13)
UNVERIFIED at max depth :      0
leaf level min/med/max  :  6 / 11 / 13
build wall time         : ~5 min on 40 cores
inner-corner y (B point): [0.04000, 0.74376]   (never below eps_down)
outer-corner y          : <= 0.74403 = Y_TOP
```

Output `cell_def_omega_mid_rebuilt.csv` columns:
`address, x_inf, x_sup, theta_inf, theta_sup, mesh_size_lower_cr, isLG, mesh_size_lower_LG, fem_order_lower_LG`.

## Files

| file | role |
|---|---|
| `omega_mid_lambda1.py` | vectorised float P1-FEM lambda1 (Richardson) |
| `reuse_db_iteration1.py` | iteration 1 (DatabaseTriangle lambda_1_lb reuse) |
| `build_cells_integrated.py` | full builder (iter1 reuse + iter2 adaptive refine) |
| `analyze_cells.py` | stats + (x,theta)/(x,y) visualization |

Run: `python3 build_cells_integrated.py <nproc> <max_depth>`
(requires `DatabaseTriangle/results/database_best_cover.csv` at the path in the script).
