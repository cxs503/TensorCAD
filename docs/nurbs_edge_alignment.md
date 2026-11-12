# Compatible edge representations (degree + knot vector) for patch seams

## 1. Case introduction

A CAD model of a hull is a network of patches. Two patches that meet along an edge
have to describe that edge **twice**, once from each side, and a stitcher, a solver or
an export can only weld the two descriptions if they are the same representation. Two
things can disagree:

* the **degree** -- one patch is quadratic, its neighbour cubic;
* the **knot vector** -- same degree, but the two edges carry different breakpoints or
  different multiplicities, so the same curve has five control points on one side and
  nine on the other.

TensorCAD already lands the first half of that operation
(`NURBSCurve.elevate_degree`, `NURBSSurface.elevate_degree`, see
`docs/nurbs_degree_elevation.md` and `docs/nurbs_surface_degree_elevation.md`), and it
had knot insertion for *surfaces* (`tensorcad.highorder.nurbs_geom.insert_knot`,
Boehm A5.1) used internally for Bezier extraction -- but no curve-level insertion and
nothing that took two edges and made them the same representation. So the sentence
that closed the elevation documents -- *"patches of an assembly that disagree in
degree cannot share one edge representation without it"* -- had no operation behind it.
After elevation the two edges of one curve still had 9 and 5 control points.

This document lands that operation (`tensorcad/edge_alignment.py`) and measures it:

| item | value |
|---|---|
| analytic anchor | exact rational circle `circle_curve(R = 2.0 m, arcs = 3)` |
| surface anchor | `sphere_surface(R = 1.0 m)` (exact rational, degree `(2, 2)`) |
| seam construction | `split_patch(sphere, axis=0, value=0.5)` -- the two halves share `u = 0.5` exactly |
| neighbour re-representation | `insert_knot(v, 0.37, 1)` then `elevate_degree(0, 1)` |
| watertight tolerance | `1e-12 m` |
| knot insertions measured | `(0.25, 2)` and `(0.7, 1)` on the circle |

Everything below is reproduced by

```bash
PYTHONPATH=$PWD python examples/nurbs_edge_alignment/run_example.py
# -> examples/nurbs_edge_alignment/outputs/edge_alignment.json
# -> examples/nurbs_edge_alignment/outputs/edge_alignment.png
python -m pytest tests/test_nurbs_edge_alignment.py -q     # 16 passed
```

## 2. Method

Four operations, each with one job (ISO 10303-42 / Piegl & Tiller A5.1 as the
reference, not a copy of a particular implementation):

1. `insert_knot(curve, value, times)` -- Boehm's single-knot insertion. The blend is
   applied to the **homogeneous** control polygon `(w x, w y, w z, w)`, so a rational
   curve stays exactly the same rational curve; inserting into the 3-D net and
   carrying the weights across is a different, generally wrong curve. An all-equal
   weight net is rescaled back to exactly equal weights, so a unit-weight curve does
   not drift onto the rational evaluation path. Insertion past multiplicity `degree`
   would break the curve: it is refused with the fix ("elevate the degree first")
   named.
2. `boundary_curve(patch, side)` -- the four edge curves `v0`, `v1`, `u0`, `u1` of a
   patch, so the caller does not do control-net index arithmetic and the seam of two
   patches can be taken directly.
3. `reverse_curve(edge)` -- the same edge traced the other way. Reversing the control
   polygon alone would not reverse the parametrisation, so the knot vector is
   reflected about the domain midpoint as well. Two neighbours that run in opposite
   directions must be brought to one direction *explicitly*; the alignment does not
   guess (and section 3 shows why guessing would be wrong).
4. `align_edges(first, second)` -- (a) elevate the lower-degree edge by the difference,
   (b) insert the breakpoints that only one side carries, (c) raise every shared
   breakpoint to the larger multiplicity. The result: equal degree, equal knot vector,
   therefore equal control count -- and, if the two edges describe the same
   parametrised curve, equal control points and equal weights.
5. `edge_report(first, second)` -- the measurement: degrees, control counts, knot
   equality, max control-point difference, max **weight** difference, and the
   geometric gap at the same parameters, plus a `watertight` flag that needs *all* of
   them. Reporting the gap is the point: alignment equalises representation, it never
   makes geometry coincide, and a report that only said "the knots now match" would be
   a false pass.

## 3. Results

### 3.1 Knot insertion does not move the curve

The exact rational circle `R = 2 m`, three 120-degree arcs, knots inserted at `0.25`
(twice) and `0.7` (once), all numbers from
`examples/nurbs_edge_alignment/outputs/edge_alignment.json`:

| quantity | before | after |
|---|---|---|
| degree | 2 | 2 |
| control points | 7 | 10 |
| knots | 10 | 13 |
| max point deviation at the same `u` | -- | **6.68e-16 m** |
| radius error against `R = 2 m` | 6.66e-16 m | 8.88e-16 m |
| arc length | 12.56637061435917 m | 12.566370614359172 m |
| arc-length relative error against `2 pi R` | 1.41e-16 | **0.0** |

### 3.2 A real seam: geometrically closed, representationally open

`sphere_surface(1.0)` split at `u = 0.5`; the neighbour's `v`-representation changed
by one inserted knot and one degree elevation (so the shared edge is the same curve,
represented as degree 2 / 5 control points on the left and degree 3 / 9 on the right):

| quantity | before alignment | after alignment |
|---|---|---|
| degrees | 2 and 3 | 3 and 3 |
| control points on the edge | 5 and 9 | **9 and 9** |
| knot vectors equal | no | **yes** |
| max control-point difference | not comparable (counts differ) | **4.45e-16 m** |
| max weight difference | not comparable | **4.44e-16** |
| max geometric gap at the same `u` | 3.51e-16 m | 4.44e-16 m |
| `watertight` | false | **true** |

The gap being round-off *before* alignment is the honest part of the table: the two
edges were always the same curve (they come from one split), the problem was that
there was no shared data to weld. Elevation and insertion cost, respectively,
`4.00e-16 m` and `0.0 m` of deviation on the two edges, and the neighbour's
`elevate_degree(0, 1)` moved the surface by `4.59e-16 m` -- the representation grows,
the geometry does not move.

The aligned edge carries the union of the breakpoints and nothing else: 9 control
points is exactly `len(knots) - degree - 1` for the union vector (asserted in
`tests/test_nurbs_edge_alignment.py::test_alignment_is_minimal`).

### 3.3 Two controls that must not come out watertight

| case | what it is | after alignment |
|---|---|---|
| `circle_curve(1.0, arcs=3)` vs `circle_curve(1.002, arcs=3)` | two edges that are **different curves** | knots equal, gap **2.00e-3 m**, `watertight` false |
| `circle_curve(1.0, arcs=3)` vs `circle_curve(1.0, arcs=6)` | same point set, **different parametrisation** | knots equal, control points agree to 4.58e-16 m, weights differ by **0.25**, gap **3.85e-2 m**, `watertight` false |
| seam edge vs the same edge reversed | neighbours running in opposite directions | knots equal, gap **2.0 m** (= `2 R`) |

The second row is the trap this module exists to avoid. Three 120-degree arcs and six
60-degree arcs are both the unit circle; after alignment even their dehomogenised
control polygons agree to round-off, so a `watertight` test based on control points
alone would return **true** while the two curves visit the circle at different
parameters (0.0385 m apart at the same `u`, 3.9% of the radius). The weight check is
what catches it: weights are 0.75 against `cos 30 deg = 0.8660254`. The third row is
the same lesson for direction -- a reversed neighbour is 2 m away at the same `u`, so
`reverse_curve` has to be called deliberately, which is why `align_edges` does not
guess a direction.

![Edge alignment evidence](../examples/nurbs_edge_alignment/outputs/edge_alignment.png)

*(Figure: `examples/nurbs_edge_alignment/outputs/edge_alignment.png`. (a) the sphere
seam at `u = 0.5` with the left patch's degree-2 control polygon, the neighbour's
degree-3 polygon and the aligned 9-point polygon both edges share; (b) gap at the same
parameter for the coincident seam, the perturbed-radius pair and the 3-vs-6-arc pair,
against the `1e-12 m` tolerance; (c) control points carried by each side of the seam
before and after alignment.)*

## 4. Error judgement

Every quantity is checked against an independent value, and every one of them is far
inside the 3% threshold -- in fact inside `1e-15`:

| key quantity | measured | benchmark | relative error | < 3%? |
|---|---|---|---|---|
| insertion: max point deviation | 6.68e-16 m | 0 (same curve) | 3.3e-16 of `R` | yes |
| insertion: radius error after | 8.88e-16 m | analytic `R = 2 m` | 4.4e-16 | yes |
| insertion: arc length after | 12.566370614359172 m | `2 pi R` | 0.0 | yes |
| seam: control-point difference after alignment | 4.45e-16 m | 0 (same curve) | 4.5e-16 of `R` | yes |
| seam: weight difference after alignment | 4.44e-16 | 0 | 4.4e-16 | yes |
| seam: geometric gap after alignment | 4.44e-16 m | 0 | 4.4e-16 of `R` | yes |
| elevation cost on the neighbour patch | 4.59e-16 m | 0 | 4.6e-16 of `R` | yes |
| negative control: radii 1.000 / 1.002 | 2.00e-3 m | 2.00e-3 m | 1.2e-13 | yes (it *must* be non-zero) |

## 5. Points to note

* **Alignment is not stitching.** It gives two edges one representation; it does not
  build the topology, does not trim and does not touch the patch interiors. A patch
  network still has to be assembled by the caller.
* **Curves only.** What is aligned is the edge curves. Aligning a whole network means
  taking each shared edge, aligning it, and using the aligned edge as the patch data --
  which is what a stitcher needs -- but the module does not rewrite the patches.
* **Parametrisation, not just degree.** Two edges of the same circle cut into three
  and into six arcs are the same point set and a *different* curve as far as a solver
  is concerned; insertion equalises knot vectors, it cannot equalise the `u -> point`
  map. The report measures that difference (3.85e-2 m at `R = 1 m`) instead of
  hiding it. Re-parametrisation of an edge is out of scope.
* **Same parameter domain required.** A domain mismatch is refused, not rescaled: the
  seam only has a meaning at the same parameters.
* **Insertion cannot exceed multiplicity `degree`** (that would break the curve). The
  alignment never runs into the limit, because it elevates to `max(degree)` first and
  every original multiplicity is at most its own degree.
* **Weights are part of the contract.** A watertight test on control points alone
  passes on the 3-vs-6-arc pair, which is wrong; the measured weight difference of
  0.25 is the reason `edge_report` reports it and `watertight` requires it.
* **Direction is explicit.** `reverse_curve` exists so the caller says what it means;
  a reversed neighbour is 2 m away at the same `u`, so a "closest of the two
  orientations" heuristic would hide a genuine modelling error.

## 6. Files

| file | role |
|---|---|
| `tensorcad/edge_alignment.py` | `insert_knot`, `boundary_curve`, `reverse_curve`, `align_edges`, `edge_report` |
| `tests/test_nurbs_edge_alignment.py` | 16 tests: insertion exactness/guards, boundary curves on the surface, reversal, alignment on a real seam, three negative controls |
| `examples/nurbs_edge_alignment/run_example.py` | produces the JSON and the figure above |
| `examples/nurbs_edge_alignment/outputs/` | `edge_alignment.json`, `edge_alignment.png` |

Related: `docs/nurbs_degree_elevation.md` (curve elevation),
`docs/nurbs_surface_degree_elevation.md` (surface elevation),
`docs/highorder-geometry.md` / `docs/CAD_HIGH_ORDER_GEOMETRY.md` (Bezier extraction,
which is knot insertion used internally), `docs/nurbs-surface.md`.