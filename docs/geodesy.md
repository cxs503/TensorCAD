# Ellipsoidal geodesy: the Earth model the ice channel was not using

> Module: `tensorcad.geodesy` · verification run: `examples/icechannel_geodesy/run_example.py`
> artifacts: `artifacts/icechannel_geodesy/{results.json,route_legs_sphere_vs_wgs84.png,earth_model_vs_projection.png}`

## 1. Why this module exists

`tensorcad.icechannel` builds the de Long -> Pevek corridor on a sphere of
`R = 6371008.8 m` and `docs/icechannel.md` budgets its projection
distortion down to 0.46 % at 1061 km. That budget is complete *given* the sphere: the
sphere itself was never a line item. Over 70-77 deg N it is one, and it is not the
small one. Both `M` (meridional) and `N` (prime-vertical) radii of curvature grow towards
the pole while `R` does not, so every distance the corridor is built from is short by
roughly the same amount -- 0.386 % over this route, i.e.
4.12 km of a 1065.9 km track.

This module supplies WGS84 geodesics so the term can be quoted instead of assumed, and it
verifies itself by two routes that share nothing but the ellipsoid constants.

## 2. Method

| ingredient | what it is | why it is here |
|---|---|---|
| `vincenty_inverse` / `vincenty_direct` | classical series solution (Vincenty 1975), iterated on `Delta lambda` to 1e-12 rad | the production route; the direct problem exists so the inverse can be round-tripped |
| `clairaut_distance` | quadrature of `ds/dphi = M / sqrt(1-(c/r)^2)` with `c = r sin(alpha)` shot on `|Delta lambda|`, `r = N cos(phi)` | an algorithmically unrelated second route: series vs quadrature |
| `meridian_arc_m` | Simpson rule on `M(phi)` | third route for `Delta lambda = 0`, where the Clairaut reduction degenerates |
| `great_circle` | haversine on the mean-radius sphere | the datum the campaign currently uses, kept so the difference is measured, not asserted |

Three facts make the cross-check meaningful rather than circular: the anchors are
external (equator arc `= a dlambda`, published WGS84 meridian quadrant
10 001 965.7293 m), the two routes differ in kind, and the quadrature refuses the cases
it cannot do (vertex inside the latitude interval, or close enough outside that Simpson
stops converging -- see §5).

## 3. Results

### 3.1 Anchors (outside the code under test)

| quantity | measured [m] | anchor [m] | difference [m] |
|---|---|---|---|
| equator arc, 1 deg (`vincenty_inverse`) | 111319.490793 | `a*pi/180` = 111319.490793 | -4.72e-08 |
| meridian quadrant 0-90 deg (`meridian_arc_m`) | 10001965.729313 | published 10 001 965.7293 | +1.27e-05 |
| equator circumference `2 pi a` | 40075016.6856 | published 40 075 016.686 | -4.22e-04 |

### 3.2 Two independent routes on every route leg

Worst disagreement over the 7 legs: **4.76e-12**
relative (series vs Clairaut quadrature; the test bound is 1e-9, the random-sample worst
case at 55-84 deg N is the same order).

| leg | sphere [km] | WGS84 [km] | sphere/WGS84-1 | route agreement |
|---|---|---|---|---|
| de Long (Bennett Island) -> East Siberian Sea north | 187.387 | 188.171 | -0.4169 | 4.4e-12 |
| East Siberian Sea north -> Medvezhyi Islands north | 216.283 | 217.169 | -0.4081 | 2.9e-13 |
| Medvezhyi Islands north -> Kolyma drift | 210.884 | 211.724 | -0.3969 | 7.0e-14 |
| Kolyma drift -> Chaun approach | 185.318 | 186.026 | -0.3805 | 1.8e-13 |
| Chaun approach -> Cape Billings south | 138.062 | 138.546 | -0.3487 | 3.4e-12 |
| Cape Billings south -> Pevek roadstead | 77.925 | 78.187 | -0.3341 | 1.0e-12 |
| Pevek roadstead -> Pevek port | 45.932 | 46.085 | -0.3327 | 4.8e-12 |

### 3.3 Totals

| quantity | value |
|---|---|
| spherical total (what `icechannel.total_km` reports) | 1061.791 km |
| WGS84 geodesic total (true surface distance) | 1065.908 km |
| shortfall of the sphere | -4.117 km (-0.3862 %) |
| `icechannel.total_km` vs the spherical sum recomputed here | +2.22e-14 % |

### 3.4 Cross-track: how wide the 250 m corridor really is

Boundary points placed 125 m either side of the centreline in *sphere plane metres*:

| station | sphere width [m] | WGS84 width [m] | relative |
|---|---|---|---|
| East Siberian Sea north | 250.000 | 251.071 | +0.4283 |
| Medvezhyi Islands north | 250.000 | 251.062 | +0.4248 |
| Kolyma drift | 250.000 | 251.052 | +0.4207 |

### 3.5 The error is in the same league as the distortion it was compared with

| latitude | `M/R - 1` | `N/R - 1` |
|---|---|---|
| 0 | -0.558 | +0.112 |
| 15 | -0.491 | +0.134 |
| 30 | -0.308 | +0.196 |
| 45 | -0.057 | +0.280 |
| 60 | +0.195 | +0.364 |
| 65 | +0.268 | +0.388 |
| 70 | +0.330 | +0.409 |
| 75 | +0.381 | +0.426 |
| 77 | +0.397 | +0.432 |
| 80 | +0.418 | +0.438 |
| 85 | +0.441 | +0.446 |

The route lives at 70-77 deg N: 0.39 % along track,
+0.42 % across track. `docs/icechannel.md`'s 0.46 % projection distortion is the same
size. Neither term may be called negligible relative to the other, and the two nearly
cancel in *area* (+0.042 %) --
which is exactly why the error is easy to miss and worth quoting as a *position* error:
an ice feature placed 1000 km along the spherical track sits 4 m
short of where the WGS84 track puts it.

## 4. What this changes downstream

* A DEM ice field or FEM ice belt built on the corridor should take its along-track
  coordinate from `geodesy.polyline_length_m(points, ellipsoid=True)` and state the
  ellipsoid; the spherical 1061.791 km is a plane coordinate, not a true distance.
* Corridor *width* in true metres must be scaled by `N cos(phi) / (R cos(phi)) = N/R`
  (+0.43 % at 75 deg N) on top of the projection factor `plane_to_true_factor` that
  `tensorcad.icechannel` already applies -- the two are multiplicative and independent.
* Nothing in `icechannel` is wrong as spherical geometry; the claim that needs revising
  is only the completeness of its error budget.

## 5. Limitations and the defect found on the way

* **Vertex case (defect found and fixed).** The first version of `clairaut_distance`
  fitted `c` by bisection and integrated one branch. For `(60, 140) -> (60.1, 175)`
  (35 deg of longitude inside 0.1 deg of latitude) the fitted vertex is at the end of the
  interval and the integrand is nearly singular there: Simpson on 512 panels returned
  **1.07 % long** instead of refusing. The fix is two guards -- an analytic vertex-latitude
  test (two branches) and an a-posteriori panel-halving test (quadrature not converged) --
  both of which raise `GeodesicError`. Regression tests:
  `test_quadrature_refuses_a_vertex_inside_the_interval`,
  `test_quadrature_keeps_its_accuracy_just_outside_the_vertex`.
* **Scope.** No datum shift, no geoid, no height above the ellipsoid; distances are
  surface geodesics. The Clairaut route is limited to monotonic-latitude legs and to
  same-hemisphere legs; `vincenty_inverse` has neither restriction.
* **Vincenty non-convergence** is reported (`converged=False`) rather than raised; it only
  occurs near antipodal points, which this route never visits.
* **Not claimed:** that the ice-breaking analysis is unaffected. Only the geometry
  datum is quantified here.

## 6. Reproduce

```bash
python examples/icechannel_geodesy/run_example.py     # writes artifacts/icechannel_geodesy/
python -m pytest tests/test_geodesy.py -q             # 32 tests
```
