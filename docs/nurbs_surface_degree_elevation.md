# NURBS **surface** degree elevation (tensor-product of ISO 10303-42)

## 1. Case introduction

An assembly is a set of patches, and almost every downstream consumer of TensorCAD
(BEM/FEM panel meshing, stitching, a solver that wants one polynomial degree per
face) needs the patches to agree on degree. `NURBSCurve.elevate_degree` landed the
curve-level operation in 0.23.0; the surface half was explicitly left open ("surface
degree elevation, trimming and stitching remain out of scope" — `CHANGELOG.md`).
A degree-`(2,1)` rational cylinder patch therefore still could not be expressed in
the degree-`(3,3)` space of its bicubic neighbour, even though both describe the
same geometry.

This document lands and measures `NURBSSurface.elevate_degree(times_u, times_v)`.

The reference is the standard algorithm, not a copy of one implementation:

* ISO 10303-42 (STEP) degree elevation for B-spline surfaces;
* Piegl & Tiller, *The NURBS Book*, §5.6 (curves) applied per parametric direction —
  the tensor-product surface is elevated one direction at a time;
* OCCT exposes the same per-direction operation as `Geom_BSplineSurface::IncreaseDegree`.

Inputs (`examples/nurbs_surface_degree_elevation/run_example.py`):

| item | value |
|---|---|
| exact rational cylinder | `cylinder_surface(R = 1.37 m, H = 5.2 m)`, degree `(2,1)`, net `7 x 2`, three 120° arcs, analytic lateral area `2 pi R H = 44.76141212834738 m²` |
| exact rational sphere | `sphere_surface(R = 2.2 m)`, degree `(2,2)`, net `7 x 5`, two singular pole rows, analytic area `4 pi R² = 60.821233773498406 m²` |
| non-rational patch | seeded (`seed = 20261111`) bicubic `4 x 4` net, knots clamped, no symmetry |
| sampled grid | `u: 25 x v: 13` check points inside the domain |
| elevation counts | cylinder `(1,0), (2,0), (3,0), (0,1), (0,2), (1,1), (2,1), (3,2)`; sphere `(1,0), (0,1), (1,1), (2,1), (1,2)` |

## 2. Method

1. **Homogeneous coordinates.** Each control point is raised to `(w x, w y, w z, w)`.
   Elevation is a linear operation on the *homogeneous* net; elevating the 3-D net
   and carrying the old weights over is a different, generally wrong surface (§5).
2. **One direction at a time.** For `u`, every fixed-`v` column is a curve for which
   the knots are the `u`-knots; for `v`, every fixed-`u` row is a curve over the
   `v`-knots. Both directions reuse `_elevate_homogeneous`, the same core that
   `NURBSCurve.elevate_degree` now calls (single source of truth, one code path for
   both halves of the operation).
3. **Knot rule.** A knot of multiplicity `m` becomes `m + times`, ends stay `p+1`-fold
   clamped (ISO 10303-42). Because the curve's continuity is `p - m`, the *continuity
   is unchanged*: a `C1` join stays `C1` in the larger degree.
4. **Minimal growth.** With `S` spans in a direction the net grows by `times * S`
   points — `n_new = n + times * S`, the minimum for the target space.
5. **Solve.** Control points are the unique coefficients of the elevated patch,
   obtained by collocation at the Greville abscissae of the new knot vector. That
   matrix is banded and totally positive, so the elimination below the diagonal never
   needs pivoting; the four homogeneous components share one elimination
   (four right-hand sides).
6. **Uniform weights stay uniform.** If every input weight is the same constant `c`
   (in particular: a non-rational patch, `c = 1`), the solve returns `c ± 4e-16`.
   Piecewise that noise would flip `NURBSSurface.rational` to `True` and silently
   move every later evaluation onto the rational path; the homogeneous net is
   therefore rescaled so the constant is exact, which does not move any point.

## 3. Results

All numbers are produced by

```bash
PYTHONPATH=$PWD python examples/nurbs_surface_degree_elevation/run_example.py
# -> examples/nurbs_surface_degree_elevation/outputs/results.json
PYTHONPATH=$PWD python tools/nurbs_surface_elevation_figure.py
# -> examples/nurbs_surface_degree_elevation/outputs/nurbs_surface_degree_elevation.png
```

Deviation is `max |S_elev(u,v) - S(u,v)|` over the 25 x 13 grid; "radius error" is
`max | |(x,y)| - R |` for the cylinder and `max | |S| - R |` for the sphere; "area
error" is `|A_elev / A_analytic - 1|`.

| surface | times `(u,v)` | degree | net | max deviation [m] | radius error [m] | area error |
|---|---|---|---|---|---|---|
| cylinder | (1, 0) | (2,1) → (3,1) | 7x2 → 10x2 | 2.75e-15 | 8.88e-16 | 4.22e-15 |
| cylinder | (2, 0) | (2,1) → (4,1) | 7x2 → 13x2 | 2.67e-15 | 6.66e-16 | 4.44e-15 |
| cylinder | (3, 0) | (2,1) → (5,1) | 7x2 → 16x2 | 3.56e-15 | 6.66e-16 | 4.66e-15 |
| cylinder | (0, 1) | (2,1) → (2,2) | 7x2 → 7x3 | 2.70e-15 | 6.66e-16 | 4.22e-15 |
| cylinder | (0, 2) | (2,1) → (2,3) | 7x2 → 7x4 | 2.70e-15 | 8.88e-16 | 4.66e-15 |
| cylinder | (1, 1) | (2,1) → (3,2) | 7x2 → 10x3 | 2.70e-15 | 8.88e-16 | 4.66e-15 |
| cylinder | (2, 1) | (2,1) → (4,2) | 7x2 → 13x3 | 2.75e-15 | 8.88e-16 | 4.22e-15 |
| cylinder | (3, 2) | (2,1) → (5,3) | 7x2 → 16x4 | 3.56e-15 | 1.33e-15 | 4.66e-15 |
| sphere | (1, 0) | (2,2) → (3,2) | 7x5 → 10x5 | 1.02e-15 | 8.88e-16 | 1.33e-15 |
| sphere | (0, 1) | (2,2) → (2,3) | 7x5 → 7x7 | 1.40e-15 | 1.33e-15 | 1.33e-15 |
| sphere | (1, 1) | (2,2) → (3,3) | 7x5 → 10x7 | 1.35e-15 | 1.33e-15 | 1.22e-15 |
| sphere | (2, 1) | (2,2) → (4,3) | 7x5 → 13x7 | 2.23e-15 | 1.78e-15 | 1.33e-15 |
| sphere | (1, 2) | (2,2) → (3,4) | 7x5 → 10x9 | 1.43e-15 | 1.33e-15 | 1.33e-15 |
| bicubic patch | (1, 1) | (3,3) → (4,4) | 4x4 → 5x5 | 4.34e-16 | — | — |
| bicubic patch | (2, 3) | (3,3) → (5,6) | 4x4 → 6x7 | 5.24e-16 | — | — |
| bicubic patch | (0, 2) | (3,3) → (3,5) | 4x4 → 4x6 | 5.44e-16 | — | — |

Structural invariants measured on the same artifacts:

| property | measured |
|---|---|
| `(2,1)` vs `(0,1)` then `(2,0)` (directions commute) | `2.22e-16 m` |
| elevate by 1 twice vs by 2 (additivity) | `2.54e-16 m` |
| non-rational patch stays non-rational | `rational == False` for all three patch rows |
| interior `C1` join (multiplicity 2, degree 3) after `(2,1)` | multiplicity 4, `continuity(order=1).passed == True` |
| deviating field max (sphere → `(4,3)`, 61 x 41 samples) | `2.23e-15 m` |
| growth rule `n + times * spans` | cylinder `7 → 10 → 13 → 16` (`S = 3`); sphere rows `5 → 7 → 9` (`S = 2`) |

![Surface degree elevation evidence](../examples/nurbs_surface_degree_elevation/outputs/nurbs_surface_degree_elevation.png)

*(Figure: `examples/nurbs_surface_degree_elevation/outputs/nurbs_surface_degree_elevation.png`,
generated by `tools/nurbs_surface_elevation_figure.py`. Panel (c) is the deviation
field over the parametric domain — the geometry-domain analogue of a field contour —
for the exact sphere elevated to `(4,3)`.)*

## 4. Error judgement against the analytic anchors

| key quantity | analytic value | elevated value | error | < 3 %? |
|---|---|---|---|---|
| cylinder radius (worst case, `(3,2)`) | 1.37 m | 1.37 m | 1.33e-15 m | **yes** (1e-13 % of `R`) |
| cylinder lateral area, `(3,2)` | 44.76141212834738 m² | 44.76141212834738 m² | 4.66e-15 relative | **yes** |
| sphere radius (worst case, `(2,1)`) | 2.2 m | 2.2 m | 1.78e-15 m | **yes** |
| sphere area, `(2,1)` | 60.821233773498406 m² | 60.821233773498406 m² | 1.33e-15 relative | **yes** |
| pointwise geometry, worst of all cases | — | — | 3.56e-15 m | **yes** |

Every key quantity is at the floating-point round-off floor of its own scale, i.e.
seven to fifteen orders of magnitude inside the 3 % acceptance band. This is the
expected outcome for an *exact* algebraic operation, not a fit: there is no
discretisation parameter to converge.

## 5. Notes, traps and limitations

* **The rational trap (measured).** Elevating the 3-D control net and dropping the
  weights gives a *different* surface: on the cylinder this loses `3.425e-01 m` of
  radius (25.0 % of `R = 1.37 m`) and `14.1 %` of the lateral area, while the
  homogeneous route errs by `6.66e-16 m`. Both directions are asserted in
  `tests/test_nurbs_surface_degree_elevation.py::TestRationalTrap`.
* **The round-off flag trap (found while building this).** Before the uniform-weight
  restore, elevating the non-rational bicubic patch returned weights
  `0.9999999999999998` / `1.0000000000000004`, so `NURBSSurface.rational` flipped to
  `True` and every later evaluation switched to the (correct but noisier) rational
  path. The regression test asserts `not elevated.rational`.
* **Singular pole rows are fine.** The sphere's `v = 0` and `v = 1` rows are singular
  in 3-D (all `u` map to one point) but *not* in homogeneous coordinates, so
  elevation along `u` is well posed there; the pole weld in `tessellate` still finds
  one vertex.
* **A discontinuous direction is refused, not fudged.** An interior multiplicity of
  `degree + 1` is a break, not a `C0` join; the collocation matrix would be singular
  by construction, so the call raises `NURBSError("... discontinuous ...")` naming
  the knot. The other direction is still elevated normally (tested).
* **Cost.** Elevation is `O(n * times)` in control points and one banded solve per
  column/row. `area()` on the elevated nets is appreciably slower than on the
  originals (the Gauss quadrature sees `13 x 7` instead of `7 x 5`), which is the
  real price of degree alignment, not the elevation itself.
* **Out of scope (unchanged by this work).** Trimming, stitching, knot removal, and
  B-rep representation. This is the *representation* half of patch-degree alignment;
  the topological half (shared edges, sewn shells) is still missing, so a real
  assembly cannot yet be elevated as a unit.

## 6. Reproduction

```bash
python -m pytest tests/test_nurbs_surface_degree_elevation.py -q     # 17 passed
python -m pytest tests/test_nurbs.py tests/test_nurbs_degree_elevation.py \
                 tests/test_hull_nurbs.py tests/test_public_surface.py -q  # 97 passed (no regression)
PYTHONPATH=$PWD python examples/nurbs_surface_degree_elevation/run_example.py
PYTHONPATH=$PWD python tools/nurbs_surface_elevation_figure.py
```

API:

```python
from tensorcad.nurbs import cylinder_surface

patch = cylinder_surface(1.37, 5.2)          # degree (2, 1)
aligned = patch.elevate_degree(1, 2)         # degree (3, 3), same surface
assert aligned.domain == patch.domain        # same parametrisation, same points
```