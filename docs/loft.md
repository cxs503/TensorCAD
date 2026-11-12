# Exact section skinning (`tensorcad.loft`)

`GAP-CAD-016` in `docs/oss_inventory_cad/GAP-ANCHORS.md` asks for **arbitrary section
sequence lofting**: the surface passes through every section (deviation `< 1e-9 L`), the
joins are `G1`/`G2` *and that is verifiable*, and the round trip through the NURBS
representation drifts `< 1e-6 L`.

`tensorcad.parametric_hull.loft_nurbs_surface` does not do this. It skins **hull
stations**: every section is a `(y, z)` polyline in its own `x = const` plane, all
sections are re-expressed in one half-breadth frame, and the transverse direction is a
fit on a shared uniform parametrisation. A tilted section, a non-planar section, a
section in a different plane, or two sections of different degrees cannot be lofted at
all, and the transverse error between sample rows is an interpolation error that only
disappears as the sample count goes to infinity.

```python
from tensorcad.nurbs import circle_curve
from tensorcad.loft import loft_sections

sections = [circle_in_plane(x) for x in (0.0, 0.5, 1.0, 1.5, 2.0)]
loft = loft_sections(sections)
loft.surface         # one tensor-product NURBSSurface, C2 in the loft direction
loft.measure()       # what it actually does, measured on the surface
```

## Method (exact, not a fit)

1. **Normalise** every section's knot vector to `[0, 1]`. This is a parameter change: no
   point moves (`test_knot_normalisation_does_not_move_the_curve`).
2. **Elevate** every section to the highest degree in the sequence (homogeneous
   elevation, so rational weights survive).
3. **Unify the representation**: the knot vector is the union of the elevated section
   knot vectors, each interior value at its maximum multiplicity, and every missing knot
   is inserted into every section (Boehm's single-knot insertion in homogeneous
   coordinates). Insertion is exact: it adds a control point and moves nothing.
4. **Interpolate the homogeneous control nets** along the loft direction on the mean
   chord-length parametrisation (Piegl & Tiller eqs. 9.1–9.8). At `u = u_i` the basis is
   the Kronecker delta, so the surface restricted to that row *is* the section: same
   degree, same knots, same weights.
5. **One patch.** No stitching, no patch seams. With `degree_u = 3` and simple interior
   knots the skin is `C2` in the loft direction.

A skin of polynomial sections is polynomial, so the loft-direction solve is done in 3-D
and the weights are exactly one; a skin with any rational section solves in 4-D
homogeneous coordinates and carries the weights row by row
(`test_weights_are_carried_row_by_row`).

## Measured evidence (`examples/loft/run_example.py`, `docs/loft/results.json`)

| Study | Reference | Measured |
|---|---|---|
| cylinder: 5 exact rational circles, `R = 0.5 m`, `H = 2 m` | `r = R` everywhere | max radial deviation **3.3e-16 m**; area vs `2 pi R H` **1.4e-14** relative |
| cone: radii `0.6 → 0.2 m` | `r(x) = r0 + k x` | max deviation **≤ 1e-13 m**; area vs `pi (r0+r1) sqrt(H^2+(r0-r1)^2)` **≤ 1e-9** relative |
| parabola family | `z = y^2` at dense `(u, v)` | max `|z - y^2|` **5.6e-16 m** |
| section interpolation (all studies) | surface restricted to `u_i` = section | defect **≤ 3e-16 L** |
| continuity, `degree_u = 3` | `C2` at simple interior knots | second-derivative jump **2.8e-16 L**, order **2** |
| continuity, `degree_u = 2` | `C1` only | second-derivative jump **2.8e-4 L**, order **1** |
| continuity, `degree_u = 1` | `C0` only | first-derivative jump **>1e-3 L**, order **0** |
| mixed degrees (polyline + ellipse) | both sections reproduced | defect **≤ 1e-14 L**, transverse degree elevated to 2 |
| rotating ellipses (non-planar-parallel) | both sections reproduced | defect **≤ 1e-14 L**, order **2** |

Acceptance `GAP-CAD-016`: surface passes through the sections at `< 1e-9 L` — measured
`≤ 3e-16 L`; `G2` in the loft direction — measured order 2 for a cubic skin.

The continuity numbers are **one-sided limits taken exactly**: `find_span` returns the
polynomial piece to the right of an interior knot and the piece to the left is one span
earlier, so the derivative nets are evaluated at the knot itself rather than at `knot +/- delta`
(a finite delta would leave an `O(delta)` floor and make a `C2` patch look discontinuous).
The `degree_u = 2` and `degree_u = 1` runs are the **controls**: the same measurement on a
deliberately reduced skin reports the finite jump, so a "pass" is informative.

## Limitations (honest, and pinned by tests)

* **No section reparametrisation.** Each section keeps its own `v` parametrisation, so the
  knot union is a *representation* union, not a geometric correspondence. Two
  differently-parametrised copies of one circle (three arcs vs four arcs) blend
  point-to-point into a non-circular interior
  (`test_sections_with_different_knot_vectors_are_both_reproduced` asserts both sections
  are exact *and* that the interior differs from the circle — the limitation is pinned,
  not hidden). Lofts of identically-parametrised sections are exact throughout
  (`test_identically_parametrised_sections_blend_to_the_same_circle`).
* **No self-intersection or fold detection.** A section sequence that folds the surface
  through itself is returned, not rejected.
* **No trimming.** Sections are curves, not trimmed boundaries.
* **Breaks are rejected, not silently smoothed.** A section with an interior knot of
  multiplicity `degree + 1` raises with the remedy in the message
  (`NURBSCurve.split_at_discontinuities()`), because degree elevation cannot represent a
  break in a single patch.
* **Degenerate strips raise.** Two consecutive sections that coincide within `1e-14 L`
  raise `LoftError` naming the pair, instead of returning a singular surface.

## Reproduce

```bash
python examples/loft/run_example.py --output docs/loft   # exits 1 if a claim fails
python -m pytest tests/test_loft.py -q                   # 31 tests
```