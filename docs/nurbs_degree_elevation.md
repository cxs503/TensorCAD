# NURBS curve degree elevation (absorbed from ISO 10303-42 / OCCT)

## 1. Case introduction

A CAD assembly is a set of patches, and a solver or a stitching operation usually
demands that they share one degree. TensorCAD could **insert knots** on a surface
(`tensorcad.highorder.nurbs_geom.insert_knot`, Boehm A5.1) but had no operation that
raised a curve's degree: a degree-2 rational arc could not be represented as the
degree-3 curve its neighbour uses, even though the two describe the same circle.
This document lands the curve-level operation, `NURBSCurve.elevate_degree(times)`,
and measures it against the analytic circle and against a deliberately naive
per-span Bezier construction.

The reference is the standard algorithm, not a copy of a particular implementation:

* ISO 10303-42 (STEP) degree elevation of B-spline curves;
* Piegl & Tiller, *The NURBS Book*, §5.6 / algorithm A5.9 (degree elevation);
* OCCT `Geom_BSplineCurve::IncreaseDegree` / `GeomConvert` expose the same operation.

Inputs used throughout the evidence (`examples/nurbs_degree_elevation/run_example.py`):

| item | value |
|---|---|
| exact rational circle | `circle_curve(R = 2.0 m, arcs = 3)` (three 120° rational quadratic arcs) |
| circle degree / control points / knots | 2 / 7 / 10 |
| analytic arc length | `2 pi R = 12.566370614359172 m` |
| continuity test curve | cubic, one interior knot at 0.4 with multiplicity 1 (declared C2) |
| elevation counts | `times = 0, 1, 2, 3` on the circle; `times = 1` on the cubic |

Everything below is reproduced by

```bash
PYTHONPATH=. python examples/nurbs_degree_elevation/run_example.py --out artifacts/m14481
python -m pytest tests/test_nurbs_degree_elevation.py -q
```

## 2. Method

The elevated curve is *not* obtained by piecing independently elevated Bezier spans
and keeping every shared end point. That route preserves the point set but changes the
data's continuity class: a C2 join of multiplicity 1 in a cubic would come back with
multiplicity 4 in the quartic (C0), and with two more control points than necessary.

The implemented route keeps the original continuity:

1. Target degree `q = p + times`.
2. Grow every knot multiplicity by `times`: an interior multiplicity `m` becomes
   `m + times`, and a clamped end `p + 1` becomes `q + 1`. This is forced by
   `continuity = degree - multiplicity`, so a C2 join stays C2.
3. The new control points are the **unique** coefficients of that same curve in the
   enlarged spline space. They are found by collocation at the Greville abscissae of
   the new knot vector. The collocation matrix is banded and totally positive, so
   Gaussian elimination needs no pivoting; the four homogeneous components share one
   elimination.
4. The solve is done in homogeneous coordinates `(w x, w y, w z, w)` and divided by
   the resulting weights, so rational arcs (`weights = (1, cos 60°, ...)`) stay exact
   instead of being reapproximated.
5. `times = 0` returns the identical object; a non-integer or negative `times` raises
   `NURBSError`.

The operation is deliberately additive: `C.elevate_degree(1).elevate_degree(1)`
produces the same degree-5 representation as `C.elevate_degree(2)` (test
`test_elevating_twice_equals_elevating_by_two`).

## 3. Results

Circle metrics (all numbers from `artifacts/m14481/degree_elevation.json`):

| times | degree | control points | knots | interior multiplicity | max point deviation [m] | max radius deviation [m] | arc-length rel. error |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 2 | 7 | 10 | 2 | 0.0 | 6.66e-16 | 1.41e-16 |
| 1 | 3 | 10 | 14 | 3 | 1.40e-15 | 1.55e-15 | 2.83e-16 |
| 2 | 4 | 13 | 18 | 4 | 9.42e-16 | 8.88e-16 | 4.24e-16 |
| 3 | 5 | 16 | 22 | 5 | 1.11e-15 | 8.88e-16 | 2.83e-16 |

Point deviation is sampled at every distinct knot and 32 points inside every span, so
the maximum includes the C0 joins. The elevated circle remains on `|x| = 2.0 m` to
round-off and its arc length remains `2 pi R` to round-off; the representation grows
by exactly `3` control points per elevation (one per span), never by more.

Cubic continuity case:

| quantity | original | elevated (`times = 1`) |
|---|---:|---:|
| degree | 3 | 4 |
| control points | 5 | 7 |
| interior knot 0.4 multiplicity | 1 | 2 |
| declared continuity at 0.4 | C2 | C2 |
| measured max derivative jump | — | 0.0 |
| max point deviation | — | 1.35e-15 m |

Negative control (the naive route the doc argues against): elevating each span as an
independent quartic Bezier and welding only the shared end point yields **9** control
points and interior multiplicity **4** (`C0`) for the same cubic. It evaluates to the
same points, but it is a different data contract for a stitcher/solver. The committed
test `test_continuity_is_preserved_and_the_representation_grows_minimally` pins the
7-point / multiplicity-2 answer, so a regression to the naive construction is red, not
silently accepted.

## 4. Figure

![degree elevation](../artifacts/m14481/degree_elevation.png)

*Left: the degree-2 rational circle's control polygon and the degree-4 elevated
polygon, both over the analytic circle (dashed). The elevated polygon is not the same
polygon -- degree elevation is the unique change of basis, so its points move onto a
different net that describes the same curve. Right: maximum point deviation between
the original and elevated curves versus the number of elevations; it stays at
`1e-15 m` for a `2 m` radius, i.e. at double-precision round-off, rather than growing
with the degree.*

## 5. Notes, traps and limitations

1. **Preserve continuity, not just points.** The obvious per-span Bezier construction
   is easier and geometrically exact, but it demotes an interior C2 join to C0 and
   inflates the control net. For patch stitching the data contract matters as much as
   the point set. This was the trap in the first draft of this work.
2. **Do not insert knots to degree `p` before elevating unless you remove them again.**
   Knot insertion to full multiplicity is a convenient way to extract Beziers, but the
   inserted knots are not part of the standard elevated representation; leaving them in
   is exactly the naive route above.
3. **Collocation needs no pivoting, but relies on total positivity.** The Greville
   collocation matrix of a B-spline is totally positive (Schoenberg-Whitney), so every
   leading minor and every elimination pivot is positive. The code asserts a non-zero
   pivot rather than silently returning `nan` if that assumption is ever broken.
4. **Rational curves are handled in homogeneous space.** Interpolating `(x, y, z)` and
   `w` separately and dividing afterwards would not preserve the conic. The circle
   tests check `|x| = R` directly, which is the observable consequence.
5. **Scope.** Curve-level only. A `NURBSSurface.elevate_degree` for both parameter
   directions is the natural next step (the same 1-D algorithm applied row by row/column
   by column); multi-patch stitching and trimming remain out of scope, as the module
   documentation already states. A curve with an interior knot of multiplicity
   `degree + 1` is genuinely discontinuous (its two sides share no control point), so
   the single-curve collocation matrix is singular; the method refuses with a named
   remedy ("split it at that knot first") rather than returning a wrong curve.
6. **Cost.** The implementation solves a dense `count x count` system with one
   elimination per call. For the curves in this library (tens of control points) that
   is microseconds; a production kernel would exploit the band explicitly. The
   algorithm's result, not its asymptotic cost, is what this increment claims.

## 6. Error judgement against references

| quantity | measured | reference | error | within 3 %? |
|---|---:|---:|---:|---|
| circle radius after 1 elevation | 2.0 m + 1.55e-15 m | 2.0 m | 7.8e-16 | **yes** |
| circle arc length after 1 elevation | 12.566370614359169 m | 12.566370614359172 m | 2.8e-16 | **yes** |
| circle arc length after 3 elevations | 12.566370614359169 m | 12.566370614359172 m | 2.8e-16 | **yes** |
| point deviation, original vs elevated | 1.11e-15 m | 0 (identity) | absolute 1.1e-15 m | **yes** |
| cubic derivative jump at C2 join | 0.0 | 0 (C2) | 0.0 | **yes** |
| cubic control-point growth | 7 points | 7 (minimal for `q = 4`, `m = 2`) | 0 | **yes** |

The point-deviation and jump references are exact identities, not engineering
approximations, so the relevant statement is "at double-precision round-off": every
error above is at least 12 orders of magnitude below the 3 % threshold. No
hydrodynamic or structural suitability is claimed from a geometric rewrite; the
operation is a representation change, and the only claims are geometric.

Regression evidence: `artifacts/m14481/red_new_test_before_implementation.txt`
(7 failed before the method existed) and
`artifacts/m14481/green_new_test_after_implementation.txt` (7 passed), plus the
targeted suite `artifacts/m14481/targeted_suite.txt` (150 passed across
`test_nurbs.py`, `test_nurbs_pole_surface.py`, `test_highorder_geometry.py`,
`test_public_surface.py`, `test_hull_nurbs.py`, and the new module).