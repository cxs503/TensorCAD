# The multi-cell shear flow reduced to the merged solver, three ways

Box 4.0 m x 3.0 m, t = 20 mm, shear force 8000 N. Second moment `I = 0.450000 m^4`, centroid `z_c = 1.500000 m`, flow scale `max|q| = 1466.667 N/m`.

## The three routes at four hand-checkable stations

| station | closed form [N/m] | merged [N/m] | graph [N/m] | oracle [N/m] | worst rel. err. |
|---|---|---|---|---|---|
| web, top corner (u=0.2) | 1166.222222 | 1166.222222 | 1166.222222 | 1166.222222 | 1.17e-15 |
| web, mid height (u=1.5) | 1466.666667 | 1466.666667 | 1466.666667 | 1466.666667 | 1.09e-15 |
| top flange, s=1 from centre | 533.333333 | 533.333333 | 533.333333 | 533.333333 | 2.34e-15 |
| bottom flange, s=1.5 from centre | 800.000000 | 800.000000 | 800.000000 | 800.000000 | 1.14e-15 |

## Agreement over the whole perimeter (41 stations per wall)

* oracle vs merged: `2.046e-12` N/m
* oracle vs graph: `2.501e-12` N/m
* merged vs graph: `1.819e-12` N/m
* the oracle's resultant: `8000.000000000` N (applied 8000 N), circulation `1.268e-12` N/m (zero by symmetry)

All three are round-off apart, i.e. the two solvers agree with a route that shares
no algebra with either of them.

## The vanishing internal web

| web thickness [m] | max difference [N/m] | difference / thickness | ratio to previous |
|---|---|---|---|
| 1e-03 | 4.9083e+01 | 49082.67 | — |
| 1e-04 | 5.1494e+00 | 51494.27 | 9.5 |
| 1e-05 | 5.1749e-01 | 51749.28 | 10.0 |
| 1e-06 | 5.1775e-02 | 51774.93 | 10.0 |
| 1e-07 | 5.1777e-03 | 51777.49 | 10.0 |
| 1e-08 | 5.1778e-04 | 51777.75 | 10.0 |
| 1e-09 | 5.1778e-05 | 51777.78 | 10.0 |

A wall whose thickness goes to zero keeps a finite stress and therefore loses its
flow: the perturbation it exerts on the rest of the section is first order in its
thickness. The ratio column is flat to two digits over seven decades, which is the
quantitative form of that statement — and the reason the limit is *not* the
round-off limit a naive test would assert.

## Figures

* `three_routes.png`
* `vanishing_web.png`
