# Sketch solver: a degenerate start and a status that claimed too much

**Date**: 2026-11-12 · **Mission**: m16386 (selfstart generalist seat, TensorCAD) ·
**Code under test**: `tensorcad/sketch.py` (`solve`, `_evaluate`, `SolveResult`) ·
**Evidence**: [`results.json`](../artifacts/sketch_robustness/results.json),
[`plot.png`](../artifacts/sketch_robustness/plot.png) ·
**Generator**: `python examples/sketch_robustness/run_example.py` ·
**Tests**: `tests/test_sketch_robustness.py` (13 cases), `tests/test_sketch.py` (46, unchanged)

This page is an adversarial review of the sketch constraint solver that landed at m16105
(`docs/sketch_constraints.md`, GAP-CAD-006) — the newest large capability in this
repository and one whose own acceptance run reports round-off-level accuracy on fourteen
fixtures. Both findings below are outside the surface that acceptance run covers: it
starts every fixture from coordinates its author typed and it never asks what the status
means when the residual does not reach `tol`.

## 1. Case: two ways a sketch is reported better than it is

### 1.1 A start on the singular set of a length row is fatal

A sketch written programmatically starts with coincident points -- that is the natural
way to say "these two points are 5 m apart" before there is any geometry to speak of:

```python
sketch.add(Fix(origin, 0.0, 0.0), Distance(free, origin, 5.0))   # free starts at the origin
```

`Distance` is the row `|free - origin| - 5`, whose derivative is `d / |d|`. At `d = 0` the
automatic-differentiation kernel correctly refuses (`sqrt derivative is singular at zero`),
and the constant escaped `Sketch.solve` as a bare `ZeroDivisionError`:

```
$ python -c "...coincident start, Distance 5..."
ZeroDivisionError: sqrt derivative is singular at zero
```

The *geometry* is well posed. Only the start is not. A solver that can only be launched
from a non-degenerate configuration is a solver whose caller must know which
configurations those are -- and the published restart study could not reveal this, because
it perturbs a well-posed start by `±0.25 m`, which never lands on the singular set.

### 1.2 `fully_constrained` did not imply "the dimensions are met"

The status was the structural verdict alone:

```python
if max_residual > residual_tol:  status = "inconsistent"
elif dof > 0:                    status = "under_constrained"
else:                            status = "fully_constrained"
```

`dof` is the rank of the Jacobian *at the coordinates the iteration stopped on*, whether
or not those coordinates satisfy the sketch. When the residual stalls between `tol`
(default `1e-12`) and `residual_tol` (default `1e-9`), or whenever a caller widens
`residual_tol` -- documented as the "unfinished versus inconsistent" threshold -- the
answer was `status == "fully_constrained"` with `converged is False`:

```
status fully_constrained  conv False  dof 0  rank 4  rows 5  maxres 6.667e-01
pos {'P0': (-0.667, 0.0), 'P1': (3.667, 0.0)}      # asked to be 5 m apart
```

A consumer that branches on `result.fully_constrained` -- the public property the API
offers -- accepts geometry that misses a dimension by 0.67 m. `converged` was right there
in the same record and told the truth; the status contradicted it.

## 2. What was changed

Minimal, in `tensorcad/sketch.py`; no residual, no tolerance and no fixture was touched.

| change | where | effect |
|---|---|---|
| `SINGULAR_START_NUDGE = 1e-6`, `_nudge_from_singular` | new | displacements `+1, -2, +3, -4, ...` × nudge: every parameter gets a distinct magnitude, so *any* two coincident parameters are separated. An alternating sign alone is not enough — the two coordinates of one point and of the next point receive the same pattern and move together (measured: the first version did exactly that and still raised) |
| restart once from the nudge | `solve` entry | a singular *start* is escaped; `SolveResult.nudged_from_singular_start` records that it happened, so the escape is never silent |
| a trial step that lands on the singular set is rejected | LM loop | the step cannot be scored, so it is damped rather than crashing the solve |
| `SketchError` naming the constraint | `solve`, `residuals`, `jacobian` | an actually ill-posed sketch (a direction constraint on a line whose endpoints the constraints hold together) fails as an API error, not an arithmetic one |
| `status = "not_converged"` when `max_residual > tol` | status ladder | the structural verdict is withheld until the residual is inside `tol`; `blocking` then names the rows above `tol` |
| `SolveResult.not_converged`, `SolveResult.solved` | new properties | `fully_constrained` can no longer be true for an unsatisfied sketch; `solved` is `converged` |

## 3. Results

Everything below is read from `artifacts/sketch_robustness/results.json` (the script
prints the same table and exits non-zero if a claim stops holding).

### 3.1 Degenerate starts now solve, to round-off

| case | start | status | dof | nudged | iterations | target | solved | error | relative |
|---|---|---|---:|---|---:|---:|---:|---:|---:|
| `single_distance` | two coincident points | under_constrained | 1 | yes | 2 | 5 m | 4.999999999999999 | 8.88e-16 m | 1.8e-16 |
| `sss_triangle` | three coincident points, sides 4/3/2 m | under_constrained | 1 | yes | 10 | 4 m | 4.000000000000001 | 8.88e-16 m | 2.2e-16 |
| `sss_triangle` | " | " | " | " | " | 3 m | 3.000000000000001 | 8.88e-16 m | 3.0e-16 |
| `sss_triangle` | " | " | " | " | " | 2 m | 2.0 | 0 | 0 |
| `zero_radius` | circle of radius 0 | fully_constrained | 0 | no | 2 | 2.5 m | 2.5 | 4.44e-16 m | 1.8e-16 |

The `zero_radius` row is the control: a radius row is *linear* (`r - value`), never
singular, so it must not be nudged -- and it is not. `dof = 1` on both point cases is
correct, not a leftover: three distances fix the shape, one vertex is pinned, and the
triangle can still spin about that vertex (the `single_distance` point slides on a
circle).

### 3.2 Restart study from the singular set

25 deterministic perturbations of the *coincident* start, jitter `±1e-3 m`,
seed `20261112`: **25 / 25 converged**, worst side-length error `9.90e-13 m`. This is the
study the m16105 acceptance run could not make: its `±0.25 m` perturbations of a
well-posed start never approach `|d| = 0`.

### 3.3 The status contract, swept over the band

System: the origin pinned, a point pinned at `(3, 0)`, and `Distance(point, origin, 5)`.
No solution exists; `dof = 0` at every stop, so the old ladder answered
`fully_constrained` in every row.

| `tol` / `residual_tol` / iterations | max residual | status now | `fully_constrained` | blocking | old ladder |
|---|---:|---|---|---|---|
| `1e-12` / `1e-9` / 100 (defaults) | 6.67e-01 | inconsistent | False | `distance`, `fix` | `fully_constrained` |
| `1e-3` / `1.0` / 100 | 6.67e-01 | not_converged | False | `distance` | `fully_constrained` |
| `1e-3` / `10.0` / 100 | 6.67e-01 | not_converged | False | `distance` | `fully_constrained` |
| `1e-30` / `1e3` / 1 | 6.67e-01 | not_converged | False | `fix` | `fully_constrained` |

The "old ladder" column is the previous rule re-derived (`dof == 0 → fully_constrained`
with no convergence test), not a run of the old code; the defect is that the rule had no
convergence test, which is exactly what the re-derivation shows. Even at the defaults the
row differs only because this particular system is inconsistent (`6.67e-1 > 1e-9`); the
dangerous window is the one the caller controls, and rows 2-4 are inside it.

## 4. Figure

![degenerate start and status contract](../artifacts/sketch_robustness/plot.png)

(a) worst position/radius error of the three degenerate-start cases, log scale, against
the closed-form answer; (b) the residual each status was returned for in the
unsatisfiable sweep -- green is `not_converged`, red is `inconsistent`; the bar heights
are identical, which is the point: only the convergence test tells them apart.

## 5. Error judgement: is the error < 3 %?

There is no baseline run to be a percentage of: the reference in every case is a closed
form (a stated dimension, or the law of cosines for the triangle's third side), so the
error is an absolute length, and the honest summary is relative error of `2.2e-16` at
worst -- round-off, sixteen orders inside a 3 % band and six inside the module's own
`1e-10 m` acceptance bound. The status findings are not an accuracy question at all: no
change to the residual moves a *wrong* status, because the residual already said
`6.67e-01`.

## 6. Pitfalls and limits of this work

* **The escape changes the start, so it can change the branch.** Newton is local
  (documented at m16105); a `±1e-6 m` displacement decides nothing a `±0.25 m` restart
  study would not also decide, but a sketch drawn *exactly* on a branch point is still a
  sketch drawn on a branch point. The nudge is recorded in the result for this reason.
* **An ill-posed direction constraint is still an error, not a solution.** The nudge
  separates coincident points; it does not invent a direction for a line whose endpoints
  the constraints hold together. Those sketches now fail with a named constraint, which
  is the recommended outcome, not a capability.
* **`not_converged` is a new status string.** Consumers that compare `status` against the
  three old values see a fourth; consumers that used the `fully_constrained` property
  see a strictly safer answer (`True` only when converged).
* **Not measured here**: the linear solver's conditioning on ill-scaled sketches (m16105
  covers unit scaling), and inequality constraints -- neither exists in this module.
* The probe that found both defects is a five-line script; it is worth noting that the
  acceptance run of m16105 passed 46 tests and 14 fixtures while both defects were
  present, because neither defect is reachable from a fixture that starts at a good point
  and converges.