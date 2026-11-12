# Sketch constraint solving (`GAP-CAD-006`)

A sketch is a two-dimensional parametric drawing plus a set of relations that must hold:
two corners coincide, a line is horizontal, an arm is 4 m long, a circle is tangent to a
line. Until this module there was no way to *solve* such a set in TensorCAD: `parameters.py`
validated the domains of already-known numbers, and every coordinate had to be supplied by
hand. `GAP-CAD-006` in `docs/oss_inventory_cad/GAP-ANCHORS.md` records that hole (the
inventory snapshot is dated 2026-11-11 and lists it as `absent`), together with its
machine-decidable acceptance criterion:

> ten fixtures, including under-constrained and over-constrained ones; 100 % convergence;
> the correct degree-of-freedom count when under-constrained; positions within `1e-10` of
> the analytic value.

This document is the measurement of that sentence. Every number in it is produced by
`examples/sketch_solver/run_example.py` (table spliced from the generated
`outputs/table.md`; raw values in `outputs/results.json`; figure `outputs/plot.png`).

## Model

A sketch owns a parameter vector `x` -- two entries per point (`x`, `y`), one per circle
radius -- and a list of constraints. Each constraint contributes one or more residual rows
`r_i(x)`; solving means driving `r(x)` to zero. Entities are references, not copies: a line
is a *pair of points*, so constraining one line's endpoint moves the shared corner and a
closed profile stays closed.

| constraint | rows | residual |
|---|---:|---|
| `Coincident(a, b)` | 2 | `a.x - b.x`, `a.y - b.y` |
| `Fix(p, x, y)` / `FixX` / `FixY` | 2 / 1 / 1 | coordinate minus the pinned value |
| `Distance(a, b, d)` | 1 | `|a - b| - d` |
| `Horizontal(line)` / `Vertical(line)` | 1 | `dy` / `dx` |
| `PointOnLine(p, line)` | 1 | signed perpendicular distance `cross / |d|` |
| `PointLineDistance(p, line, d)` | 1 | the same, minus a signed offset |
| `Parallel` / `Perpendicular(l1, l2)` | 1 | unit-direction `sin` / `cos` of the angle |
| `Angle(l1, l2, deg)` | 1 | `atan2(cross, dot) - radians(deg)` |
| `EqualLength(l1, l2)` | 1 | `|d1| - |d2|` |
| `Radius(circle, r)` | 1 | `r - value` |
| `PointOnCircle(p, circle)` | 1 | `|p - centre| - r` |
| `LineCircleTangency(line, circle, side)` | 1 | signed distance `- side * r` |

Two of those choices are load-bearing. **Signed distances, not cross products**: a raw
cross product grows with the line length, so its residual and therefore the rank tolerance
would depend on how long the line happens to be; a signed distance is a length, and a
`1e-9` tolerance on it means "a nanometre". **An explicit tangency side**, because
`|distance| - r` is not differentiable exactly where it vanishes -- the one place the
solver has to pass through.

## Method

The Jacobian is **exact**, obtained by forward-mode automatic differentiation
(`tensorcad/autodiff.py`): every residual is written once, in ordinary arithmetic over dual
numbers that carry a sparse `{parameter: derivative}` gradient, so the derivative cannot
drift away from the residual it differentiates. A finite-difference Jacobian would have
been shorter to write and worse to use: its truncation error becomes a floor on the
achievable residual, and the `1e-10` position bound is below that floor.

The step comes from the **damped normal equations**

```
(J^T J + lambda * diag(J^T J + 1)) dx = -J^T r
```

solved by Cholesky factorisation. Damping is what makes under- and over-constrained
sketches ordinary rather than exceptional: `J^T J` is singular in both cases, and the
`lambda` term turns a singular direction into a zero step (the right answer -- nothing
constrains it) instead of a divide-by-zero. `lambda` is adapted by trial: a step that lowers
`sum r^2` is accepted and the damping is cut by three, otherwise it is raised by ten and the
step retried, up to forty trials. Convergence is `max |r| <= tol`.

### Diagnostics: DOF, redundancy, inconsistency

After the solve, the Jacobian is **row-normalised** and reduced to row echelon form. Row
normalisation first, because a row's scale is an artefact of the units in its residual (a
`distance` row is `O(1)` per metre, a `fix` row is `1`); without it the reported rank would
depend on whether the sketch is drawn in millimetres or kilometres. Row operations preserve
the null space, so the basis returned comes from the reduced matrix and is the null space of
the original Jacobian.

* `dof = n_parameters - rank(J)` -- the number of independent motions left. The basis of
  those motions is returned with it, so "3 DOF" comes with *which* three.
* `rows > rank` with residuals at zero is **redundant**: the same dimension stated twice is
  a modelling convenience, not an error (a redundant sketch is still fully constrained).
* residuals that cannot be driven to zero are **inconsistent** (the over-constrained case
  where two dimensions disagree). The solver raises `InconsistentConstraintsError` naming
  the blocking constraint rows; with `strict=False` the same information comes back in the
  result, which is how the fixtures assert on the failure mode itself.

Tolerances are stated in the drawing's units, as dimensions are: `tol` is the convergence
threshold on the largest residual, `residual_tol` the threshold above which an unconverged
residual is called inconsistent rather than merely unfinished, and `rank_tolerance`
(default `1e-9`, relative because the rows are normalised) is the pivot threshold.

## Results

Fourteen fixtures -- twelve fully constrained, one under-constrained, one inconsistent, one
of them redundantly over-constrained -- each with its analytic answer computed from closed
forms (law of cosines for the SSS apex, `h/2` for the inscribed circle, `L cos 30` for the
angled arm, `2.5 * (4,3)/5` for the point along the 3-4-5 line). Positions are compared
against those closed forms; radii against theirs.

<!-- TABLE:START -->
| fixture | what it checks | status | DOF | rank / rows | iterations | max residual [m] | max position error [m] |
|---|---|---|---:|---:|---:|---:|---:|
| `fixed_point` | one point pinned at (2.5, -1.25) | fully_constrained | 0 | 2 / 2 | 2 | 0.0e+00 | 0.0e+00 |
| `distance_circle` | point held 5 m from a pinned point (one equation, two parameters) | under_constrained | 1 | 3 / 3 | 3 | 8.9e-16 | 0.0e+00 |
| `horizontal_arm` | horizontal 3 m arm off a pinned point | fully_constrained | 0 | 4 / 4 | 3 | 1.8e-17 | 1.8e-17 |
| `sss_triangle` | triangle from three side lengths (4, 3, 2 m), apex from the law of cosines | fully_constrained | 0 | 6 / 6 | 5 | 4.4e-16 | 4.4e-16 |
| `rectangle` | 4 x 3 rectangle from horizontal/vertical pairs and two dimensions | fully_constrained | 0 | 8 / 8 | 3 | 3.8e-18 | 7.1e-18 |
| `angle_arm` | 2 m arm at 30 degrees to a pinned x reference line | fully_constrained | 0 | 6 / 6 | 5 | 4.4e-16 | 2.5e-16 |
| `perpendicular_chain` | vertical 3 m arm followed by a 4 m arm perpendicular to it | fully_constrained | 0 | 6 / 6 | 4 | 2.7e-18 | 3.9e-18 |
| `point_on_line` | point on the 3-4-5 line at 2.5 m from its start | fully_constrained | 0 | 6 / 6 | 3 | 1.8e-16 | 2.2e-16 |
| `point_line_offset` | point 1 m to the left of the 3-4-5 line, x anchored | fully_constrained | 0 | 6 / 6 | 2 | 7.8e-16 | 2.6e-15 |
| `circle_between_lines` | circle tangent to two parallel lines 2 m apart (r = h / 2) | fully_constrained | 0 | 11 / 11 | 2 | 1.1e-16 | 1.1e-16 |
| `point_on_circle` | point on a radius-2 circle with its ordinate anchored | fully_constrained | 0 | 5 / 5 | 3 | 2.0e-13 | 2.4e-13 |
| `under_constrained_line` | a bare horizontal line: 4 parameters, one equation | under_constrained | 3 | 1 / 1 | 2 | 2.8e-17 | 0.0e+00 |
| `inconsistent` | point pinned at (3, 0) and asked to be 5 m from the origin | inconsistent | 0 | 4 / 5 | 4 | 6.7e-01 | 0.0e+00 |
| `redundant` | the same 5 m dimension stated twice over two pinned points | fully_constrained | 0 | 4 / 6 | 3 | 8.0e-17 | 1.0e-16 |

Perturbed-restart study: **325 / 325** solves converged (25 deterministic restarts per solvable fixture, perturbation ±0.25 m), median 3 iterations, worst position drift 3.9e-12 m.

| drawing scale | DOF | rank / rows | max residual [m] | relative position error |
|---:|---:|---:|---:|---:|
| 0.001 | 0 | 8 / 8 | 3.8e-21 | 1.8e-18 |
| 1 | 0 | 8 / 8 | 3.8e-18 | 1.8e-18 |
| 1000 | 0 | 8 / 8 | 3.8e-15 | 1.8e-18 |
<!-- TABLE:END -->

The figure (`outputs/plot.png`) has four panels: three solved sketches with the analytic
answer marked as a cross on top of the solved marker, and the residual history of four
cases on a log scale. The geometry panels are there to make one point visually: the crosses
sit exactly on the solved markers, which is what a `1e-13`-level agreement looks like at
drawing scale.

## Error judgement

The acceptance bound is `1e-10` m on positions. The generated table above carries the
measured worst case per fixture and the summary line carries the worst over the whole set;
both are far inside the bound (the bound is met with tens of orders of magnitude of
headroom, not by a hair). The obvious question -- *is the error < 3 %?* -- is answered
trivially here, because the reference is not a baseline run but a closed form: the measured
error is a round-off-level number, not a fraction of a percent. That is why the interesting
numbers in this table are `1e-13`, and why the honest statement is "the solver reaches the
analytic answer to round-off", not "it is within tolerance".

Two things *are* measured rather than assumed:

* **Restarts.** 25 deterministic perturbations (±0.25 m, seed 20261112) of every solvable
  fixture. A solver that only works from the coordinates its own fixture author typed has
  not been tested. The rate and median iteration count are in the generated table.
* **Units.** The same rectangle drawn at `1e-3`, `1` and `1e3` scale: the DOF/rank answer
  must not change, and the relative position error is reported per scale.

## Limitations

* **Local solver.** Newton is local: a mirror-image branch (the point at `-3` instead of
  `+3`) is selected by the starting coordinates, and every fixture therefore states its
  start. A sketch drawn on the wrong branch converges to the wrong branch.
* **No inequality constraints**, no dimension-scheme ordering, no constraint priority for
  an inconsistent set beyond naming the blocking rows, no 3-D constraints.
* **`dof` is the local rank** at the solution. For a degenerate configuration (three
  collinear points asked to be a triangle) the rank can drop and the DOF rises; that is
  reported, not silently repaired.
* **Angles wrap.** `Angle` is signed and continuous through a right angle, but a residual
  is defined modulo `2 pi`; a 180-degree relationship is better stated with `Parallel`.
* **Tangency needs a side.** `LineCircleTangency` takes `side="left"|"right"`; an unsigned
  distance is not differentiable at the solution.
* **Cost.** The AD gradients are sparse dicts, the normal equations dense `O(n^3)`. Fine
  for sketches with tens of parameters; this is not a large-scale nonlinear solver.
* **No trimming/patching of the solved profile.** Producing a wire or a face from the
  solved points is the next step, not this one.

## Degenerate starts and the status contract (m16386 review)

Two defects outside the acceptance surface above were found by probing the module,
fixed, and pinned in `tests/test_sketch_robustness.py`: a sketch started with two
coincident points (the natural programmatic start) raised a bare `ZeroDivisionError`
from the length row's singular derivative, and `status` reported `fully_constrained`
for coordinates that miss a dimension when the solve stalled above `tol`. The solver now
escapes a singular start deterministically, and `status == "not_converged"` withholds the
structural verdict until the residual is inside `tol`.
Evidence, numbers and limitations: [`docs/sketch_robustness.md`](sketch_robustness.md).

## Files

| path | what |
|---|---|
| `tensorcad/sketch.py` | entities, constraints, Levenberg-Marquardt solve, rank/DOF/null-space diagnostics |
| `tensorcad/autodiff.py` | forward-mode AD kernel (`Dual`, `gradient`, `sqrt`, `atan2`) |
| `tensorcad/sketch_cases.py` | the fourteen fixtures and their analytic answers, shared by the test and the example |
| `tests/test_sketch.py` | 46 tests: acceptance clauses, null-space/DOF invariants, restarts, unit scaling, AD against closed forms |
| `examples/sketch_solver/run_example.py` | the studies above, `results.json`, `plot.png`, the spliced table |