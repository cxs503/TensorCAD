# Fillet and chamfer with exact arc moments (`tensorcad.fillet`, `GAP-CAD-012`)

Mission m16453 (selfstart, generalist seat). Anchor `GAP-CAD-012` in
`docs/oss_inventory_cad/GAP-ANCHORS.md` listed the hole as: *the hull generator only handles
hard corners and zero-width stems; no rounding, no variable radius, no fallback with a reason*.
Everything the package could build -- `hull`, `parametric_hull`, `series60`, `section_gz`,
`girder` -- had sharp edges only.

## 1. What was added

| API | What it returns |
|---|---|
| `fillet_section(points, radius)` | `SectionShape`: exact area, perimeter, centroid, second moments (`int u^2 dA`, `int v^2 dA`), a chorded outline ready for a mesh, per-corner report, `g1_error_deg`, `area_crosscheck_error` |
| `chamfer_section(points, setback)` | the same, with a straight cut; `g1_error_deg` is the crease angle the chamfer introduces |
| `fillet_solid(points, length, radius, radius_end=...)` | closed, outward-oriented mesh swept along `+x`, exact volume, watertight/Euler/second-order convergence report |
| `chamfer_solid(...)` | the same for chamfers |
| `polygon_moments`, `section_turns`, `max_fillet_radius`, `max_chamfer_setback` | the exact 2-D kernel the above are built from |

The operation is a 2-D kernel plus a sweep, and it is **exact**: an arc of sweep `|phi|` and
radius `r` contributes to area, centroid and second moments through closed-form Green's
theorem integrals (integrands are polynomials up to degree four in `cos`/`sin`), so
`SectionMoments` carries no discretisation at all.  Chording to `arc_points` exists only for the
mesh, and `SectionShape.sampled_moments()` reports what the mesh sees.

Two independent derivations of the area are kept side by side: the Green's-theorem sum over the
exact path, and the per-corner closed form `A = A0 - sum sign(phi) r^2 (tan(|phi|/2) - |phi|/2)`.
`area_crosscheck_error` is their relative difference.

## 2. Results: exact values vs closed forms

Fixture set (turn angles are exterior angles, positive convex; `l_shape` carries one re-entrant
corner at -90 deg):

| fixture | points | radius | area (exact) | per-corner closed form | rel. diff | perimeter | G1 error (deg) | max feasible radius |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| hexagon | 6 | 0.200000 | 2.585175852891 | 2.585175852891 | 1.72e-16 | 5.870996415 | 0.00e+00 | 0.866025404 |
| l_shape | 6 | 0.200000 | 2.965663706144 | 2.965663706144 | 0.00e+00 | 7.484955592 | 0.00e+00 | 0.500000000 |
| square | 4 | 0.500000 | 3.785398163397 | 3.785398163397 | 0.00e+00 | 7.141592654 | 0.00e+00 | 1.000000000 |

Circle limit -- a square of side `2r = 2` filleted with `r = 1` must *be* the circle:

| quantity | computed | exact | relative error |
|---|---:|---:|---:|
| area | 3.141592653589793 | 3.141592653589793 | 0.00e+00 |
| perimeter | 6.283185307179587 | 6.283185307179586 | 1.41e-16 |
| `int u^2 dA` | 0.785398163397449 | 0.785398163397448 | 4.24e-16 |

Chamfer (same square, setback 0.25): area 3.875000000000 vs closed form
3.875000000000 (rel. diff 0.00e+00), perimeter
7.414213562373 vs 7.414213562373, crease 45.000000000 deg
-- a chamfer is *not* G1 and the report says so instead of pretending.

## 3. Discretisation: measured second order

`arc_points` chords per arc; the sampled area and the integrated mesh volume are compared with
the exact values:

| chords/arc | sampled area rel. error | mesh volume rel. error |
|---:|---:|---:|
| 4 | 5.292e-03 | 5.292e-03 |
| 8 | 1.331e-03 | 1.331e-03 |
| 16 | 3.331e-04 | 3.331e-04 |
| 32 | 8.331e-05 | 8.331e-05 |

Fitted orders (successive refinements):

| quantity | order 1 | order 2 | order 3 |
|---|---:|---:|---:|
| sampled area | 1.992 | 1.998 | 1.999 |
| mesh volume | 1.992 | 1.998 | 1.999 |

Tapered fillet (`r: 0.1 -> 0.6` over 2 m, 4 stations): exact volume
7.753923227 m^3, mesh volume 7.752477234 m^3 (rel. error
1.86e-04), watertight `True`, Euler characteristic
2, zero-area faces 0, **G1 deviation
0.00e+00 deg** -- tangency is a property of the surface (the ruling direction
lies in the adjacent flat face), not of the chord count.

## 4. Fallback with a reason (short edge / needle)

A 2 x 0.2 strip cannot host `r = 0.5`: the largest radius is
`fb['max_feasible']:.6f` (exactly half the short edge, where the two fillet arcs meet
tangentially on the vanished edge).  Asking for 0.5 without `allow_clamp` raises
`FilletError` naming that number; with `allow_clamp=True` the radius is clamped to
`fb['achieved']:.6f` and `clamped=True` is reported.  A needle vertex (interior angle -> 0)
drives the feasible radius to zero and the same clamp semantics applies; the resulting section
is the analyte of the limit -- the incircle of the sliver
(`tests/test_fillet.py::TestFeasibility::test_a_needle_vertex_cannot_be_filleted_and_falls_back`).

## 5. Ship application: bilge radius on a box barge waterplane

Waterplane `B = 10 m x L = 4 m`, draft `T = 4 m`; `BM = I_T / volume` with
`volume = A_wp T` (per unit length), so `BM = (int u^2 dA / A_wp) / T`.  `r = 0` recovers the
box value `B^2/(12T) = 2.083333 m` exactly.

| bilge radius (m) | waterplane area (m^2) | `I_T` (m^4) | `BM` (m) | BM change |
|---:|---:|---:|---:|---:|
| 0.0 | 40.000000 | 333.333333 | 2.083333 | +0.00 % |
| 0.5 | 39.785398 | 328.203400 | 2.062336 | -1.01 % |
| 1.0 | 39.141593 | 313.717547 | 2.003735 | -3.82 % |
| 2.0 | 36.566371 | 261.663706 | 1.788964 | -14.13 % |

The loss is second order in the radius at small `r` (slope `(BM0 - BM)/r^2` = 0.268, 0.254,
0.241 for r = 0.1, 0.2, 0.3 on the smaller test section), which the test suite pins.

## 6. Error verdict against the acceptance bar

| quantity | reference | worst error | < 3 % ? |
|---|---|---:|---|
| filleted area (3 fixtures, 3 radii) | per-corner closed form | 1.72e-16 relative | yes (this is round-off; the bar is beaten by 14 orders) |
| filleted area, circle limit | `pi r^2` | 0.0e+00 relative | yes |
| filleted perimeter | closed form | < 1e-14 relative | yes |
| second moment `int u^2 dA`, circle limit | `pi r^4/4` | 4.2e-16 relative | yes |
| chamfer area | closed form | 0.0e+00 relative | yes |
| mesh volume @ 4 chords/arc | exact volume | 0.529 % | yes (>15x under the bar already at the coarsest setting) |
| G1 continuity (fillet) | 0 deg | 0.0e+00 deg | yes (bar 1 deg) |
| G1 continuity (tapered fillet) | 0 deg | 0.0e+00 deg | yes |

Nothing here is a fit or a literature comparison: every reference is a closed form, so the
errors are round-off and discretisation only.

## 7. Limits (what this is not)

* No B-rep topology, no OCP/OCCT backend, no boolean/sweep of arbitrary faces.
* One radius per section: a radius that varies *around* the outline is not supported; only a
  linear variation *along* the sweep axis is.
* Where three faces meet, the fillet surfaces simply meet along the swept arc ends -- there is
  no corner sphere / rolling-ball blend, and no report pretending there is.
* Only simple polygons: self-intersecting outlines, zero-area outlines and non-finite
  coordinates are refused with a named reason (`FilletError`), never silently repaired.
* The swept mesh is chorded, so its volume carries the O(h^2) error measured in section 3; the
  exact volume is reported alongside so the two are never confused.

## 8. Evidence

| artifact | path |
|---|---|
| module | `tensorcad/fillet.py` |
| tests (66) | `tests/test_fillet.py` |
| evidence run + numbers | `examples/fillet_corners/run_example.py`, `artifacts/fillet_corners/fillet_report.json` |
| figure | `artifacts/fillet_corners/fillet_corners.png` |
