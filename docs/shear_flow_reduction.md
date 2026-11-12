# The multi-cell shear flow reduced to the merged solver — three routes, one answer

*Verification debt closed: `tests/test_shear_flow_multicell.py` carried three `xfail` items
("the pointwise reduction to the merged solver is not yet wired up correctly in this
test"). The physics was never in doubt; the **comparison** was wrong. This document records
the diagnosis, the third route that makes the claim testable, and the numbers.*

## 1. What is being claimed

`tensorcad.shear_flow` solves one open tree or one closed loop by sweeping a first moment
`Q(s)`; `tensorcad.shear_flow_multicell` solves a graph with one unknown flow constant per
wall plus one no-twist equation per cell. They must agree wherever the bodies coincide:

* a two-cell box whose internal web vanishes must reproduce the merged single-cell box;
* a one-cell box must reproduce the merged box point by point;
* an open tree (an I-section) must reproduce the merged open-section solver.

Agreement of two implementations is necessary, not sufficient — a shared misreading would
survive it. So the claim is checked against a **third route** that shares no algebra with
either: `tests/shear_flow_reduction.py::PerimeterOracle` integrates the differential relation

    dq/ds = -(V / I) t(s) (v(s) - v_c)

around the perimeter from a cut and closes it with Bredt's no-twist condition
`CLOSED_INT q/t ds = 0`. It forms no first moment anywhere, its global sign is fixed by the
resultant (as the modules' is), and because every wall is straight the integrand is
piecewise quadratic, so every integral is closed form: the oracle is round-off accurate, not
discretisation limited. It is invariant to where the cut is placed and to which way the loop
is walked (both asserted).

## 2. Why the comparison failed while the physics held

Three defects in the *test*, each a trap a caller will meet again
(`tests/shear_flow_reduction.py` documents them at the top of the module):

| # | Defect | Symptom |
|---|---|---|
| 1 | The arc sense is not shared between the modules | `q_at(f)` vs `q_at(f)` compares two points mirrored about the wall centre; error looks O(1) |
| 2 | A relative tolerance at a zero crossing | `|q|` is zero at a free edge / a symmetric flange midpoint; `1e-13` against `5e-13` scores as 100 % error |
| 3 | A spliced wall read as an unspliced one | the multi-cell body must split a flange where a web meets it; wall-by-wall matching then finds no partner |

The comparators now map both reports onto a shared convention (`canonical_q`: flow measured
from the wall's lexicographically smaller endpoint), compare as a **field of the physical
point**, and divide every tolerance by the section's own flow scale (`flow_scale`), never by
the value being compared.

## 3. Results — box 4 m × 3 m, t = 20 mm, V = 8 kN

`I = 0.45 m^4`, `z_c = 1.5 m`, `max|q| = 1466.667 N/m`. Four stations a reader can check on
paper (the closed form comes from the symmetry cut, where the circulation is zero):

| station | closed form [N/m] | merged | graph | oracle | worst rel. error |
|---|---|---|---|---|---|
| web at top corner (u = 0.2 m) | 1166.222222 | 1166.222222 | 1166.222222 | 1166.222222 | 1.17e-15 |
| web at mid height (u = 1.5 m) | 1466.666667 | 1466.666667 | 1466.666667 | 1466.666667 | 1.09e-15 |
| top flange, 1.0 m from centre | 533.333333 | 533.333333 | 533.333333 | 533.333333 | 2.34e-15 |
| bottom flange, 1.5 m from centre | 800.000000 | 800.000000 | 800.000000 | 800.000000 | 1.14e-15 |

Over the whole perimeter (41 interior stations per wall, which excludes the corners where two
walls with opposite canonical directions meet):

| pair | max difference [N/m] |
|---|---|
| oracle vs merged | 2.046e-12 |
| oracle vs graph | 2.501e-12 |
| merged vs graph | 1.819e-12 |

Oracle self-checks: resultant `8000.000000000 N` against an applied `8000 N`; circulation
`1.27e-12 N/m` at the symmetry cut (zero, as symmetry requires) and `-533.33 N/m` at an
off-centre cut (the no-twist condition is doing work, not being satisfied by accident).

## 4. The vanishing web is first order, not round-off

The reduction test used to demand `1e-5` agreement at a `1e-3` web. That was the wrong
statement: a wall whose thickness goes to zero keeps a **finite stress**, so it loses its
**flow** (`q_web = tau · t_web`), and the perturbation it exerts on the rest of the section is
first order in its thickness.

| web thickness [m] | max difference on the side webs [N/m] | difference / thickness | ratio per decade |
|---|---|---|---|
| 1e-3 | 4.908e+01 | 49082.67 | — |
| 1e-4 | 5.149e+00 | 51494.27 | 9.5 |
| 1e-5 | 5.175e-01 | 51749.28 | 10.0 |
| 1e-6 | 5.178e-02 | 51774.93 | 10.0 |
| 1e-7 | 5.178e-03 | 51777.49 | 10.0 |
| 1e-8 | 5.178e-04 | 51777.75 | 10.0 |
| 1e-9 | 5.178e-05 | 51777.78 | 10.0 |

`difference / thickness` is flat at `5.178e4` to five significant digits over seven decades:
a measured convergence **rate** (and a prediction that the physics makes before the numbers
do). The test now asserts the rate, not an invented ceiling.

## 5. What changed

* `tests/shear_flow_reduction.py` (new) — the three traps, the geometry-aware comparators
  (`canonical_q`, `flow_scale`, `difference_on_spans`, `field_at`, `relative_error`) and
  `PerimeterOracle`, the independent differential route.
* `tests/test_shear_flow_multicell.py` — the three `xfail` items are green assertions on the
  scale-relative difference; the vanishing-web test asserts first-order convergence.
* `tests/test_shear_flow_reduction.py` (new) — the closed-form anchor (both solvers), the
  oracle against both solvers over the whole perimeter, the oracle's cut/direction/shear-force
  invariance, and one regression test per comparison trap.
* `examples/shear_flow_reduction/` — the runnable study behind this document
  (`outputs/results.json`, `outputs/report.md`, `outputs/three_routes.png`,
  `outputs/vanishing_web.png`).

## 6. Scope and limits

* One closed cell and open trees; the multi-cell case is covered by reduction, not by a
  closed form of its own (an unequal two-cell box has no single-circle solution to quote).
* Constant thickness per wall, straight walls, no product of inertia — the modules'
  own idealisation, unchanged here.
* The oracle's "round-off accurate" claim rests on straight walls; a curved wall would need
  quadrature and would be discretisation limited.

## 7. References

* Bredt, R. (1896) *Kritische Bemerkungen zur Drehungselastizität* — the no-twist circulation.
* Megson, T.H.G. *Aircraft Structures for Engineering Students*, §"shear flow in closed
  sections" — the first-moment-from-a-cut construction used for the closed form.
* Timoshenko, S. *Theory of Elasticity* — thin-walled shear flow for the open tree.
