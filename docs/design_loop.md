# The design loop

Everything else in TensorCAD answers "what is this surface?". This module answers the
question an agent actually receives: *"give me a hull that displaces 12.06 m³ with a
B/T of 2.67 and a wetted area under 33.3 m², and prove it."*

```
targets ─► guess ─► build + measure ─► residuals ─► Newton step ─► bounds ─► repeat
                                 ▲                                          │
                                 └──────────────────────────────────────────┘
```

The loop is a damped Newton iteration on `log` of the free dimensions, with a
finite-difference Jacobian, a backtracking line search, hard bounds from the hull model's
own parameter specs, and **every intermediate design kept** so the record shows the path
rather than only the answer. It is standard-library only, like the hull model it drives.

Everything below is measured by `python examples/design_loop/run_example.py`, which exits
non-zero if any claim stops holding. `tests/test_design.py` pins the same behaviour.

## What a design problem is

| Piece | Meaning |
|---|---|
| equalities | `displacement_m3`, `beam_to_draft`, `length_to_beam`, `draft_m` — the loop drives these to zero residual |
| budget | `wetted_area_m2_max` — measured every iteration, enforced on the verdict, not on the step |
| free variables | `length`, `beam`, `draft` — at least as many as there are equalities, checked up front |
| verdict | `converged` **and** every check and budget green, or the design is not delivered |

Two geometric integrals are deliberately *not* offered as a target pair. Displacement and
wetted area both scale roughly with beam times draft on this family, so asking for both
makes the Jacobian nearly rank deficient: measured, the loop stalls at a residual of
4.0 × 10⁻³ and cannot remove it, because no hull in the family has that pair. A shape
ratio and a displacement are independent, and that is what a naval architect specifies
anyway. The residual table in the example is the evidence.

## Measured: recovering a known hull

The brief is built from a reference hull so the answer is checkable: L = 12 m,
B = 2.4 m, T = 0.9 m → displacement 12.063418 m³, wetted area 31.743484 m². The loop
starts from L = 10 m, B = 2.0 m, T = 0.8 m and is given the reference's displacement,
B/T and L/B.

| it | L (m) | B (m) | T (m) | ‖r‖ | V (m³) | A (m²) | digest |
|---|---|---|---|---|---|---|---|
| 0 | 10.000000 | 2.000000 | 0.800000 | 3.88e-01 | 7.446554 | 22.838010 | — |
| 1 | 12.572023 | 2.514405 | 0.940897 | 1.48e-01 | 13.842608 | 34.802342 | fd3d8c9f |
| 2 | 12.036239 | 2.407248 | 0.902716 | 9.08e-03 | 12.173012 | 31.935460 | 7257c38a |
| 3 | 12.000163 | 2.400033 | 0.900012 | 4.08e-05 | 12.063910 | 31.744347 | 914cfba1 |
| 4 | 12.000000 | 2.400000 | 0.900000 | 8.31e-10 | 12.063418 | 31.743484 | 21dae993 |
| 5 | 12.000000 | 2.400000 | 0.900000 | 2.48e-16 | 12.063418 | 31.743484 | 814f4770 |

Six evaluations, 1.1 s, and the reference hull comes back to 1e-7 relative on all three
dimensions. The last row exists because the loop evaluates before checking convergence;
its digest differs from row 4 by the last bit of a vertex, which is worth seeing — the
mesh is the artifact, and the design digest tracks it, not the parameters.

## Measured: is it really Newton?

Two independent checks, both of which would fail if the Jacobian were wrong:

- **Observed order 1.940** from the residual triple (9.1e-3, 4.1e-5, 8.3e-10). The
  estimate deliberately ignores the last iterations: the central difference has a
  truncation error of order `step² ≈ 1e-10`, so the tail stops contracting quadratically
  at the finite-difference floor rather than at the solution.
- **Sensitivities** `d ln(measured) / d ln(variable)` at the first design:

  | variable | displacement | B/T | L/B |
  |---|---|---|---|
  | length | +1.000000 | 0.000000 | +1.000000 |
  | beam | +1.000000 | +1.000000 | −1.000000 |
  | draft | +1.000000 | −1.000000 | 0.000000 |

  A ratio has sensitivity exactly ±1 to its own numerator and denominator, and a hull
  whose sections are proportional to the beam displaces linearly in beam and draft. Both
  hold to 1e-9, which says the Jacobian is measuring geometry, not grid artefacts.

## Measured: state, replay and re-open

- **Replay is bit-identical**: the same problem solved twice returns the same digest
  `814f4770b0189c58` and the same parameters. The generator is deterministic, so a design
  run is reproducible from its targets alone.
- **Re-opening carries state**: asking for 15 % more displacement and starting from the
  previous answer converges in 5 evaluations to L = 12.572275, B = 2.514455,
  T = 0.942921, V = 13.872930 m³ — bigger in every dimension, and the iteration count is
  one less than a cold start.

## Measured: the verdict the loop must refuse

The case that matters is the one that *converges* to a design the validation layer then
rejects -- a design loop that only ever reports success is not a design loop:

```
target 200 m³ at L = 10 m  ->  status 'converged', |r| = 1.13e-12, solved = False
failed gates: mesh_min_angle, mesh_aspect_ratio, mesh_sliver_share
```

The residual reaches 1e-12 — the Newton step did exactly what it was asked — and the
design is still refused, because at L/B = 0.965 the triangles have a 1.38° minimum angle
at a 41.5 aspect ratio and the validation layer rejects them. The manifest says
`solver.status: "not_ready"` and names the three failed gates. A design loop that reports
this as a success is worse than no design loop.

## Measured: how much of the answer is the grid

The loop converges on the grid it measures, so the grid is part of the answer:

| resolution | displacement (m³) | wetted area (m²) |
|---|---|---|
| (41, 13) | 12.063418 | 31.743484 |
| (81, 25) | 12.133078 | 31.799270 |

+0.577 % on displacement, +0.18 % on area. That is the number to quote next to any design
result, and it is in the record.

## The solver manifest

`solution.manifest()` is what leaves the repository:

```json
{
  "schema": "tensorcad.design/1",
  "status": "converged",
  "solved": true,
  "targets": {"displacement_m3": 12.063418, "beam_to_draft": 2.666667, ...},
  "design_variables": ["length", "beam", "draft"],
  "path": [ ... one full record per evaluation, with hydrostatics, handoff and checks ... ],
  "reference_quantities": {"length_m": 12.0, "wetted_area_m2": 31.743484,
                            "water_density_kg_per_m3": 1025.0, "block_coefficient": 0.465410, ...},
  "handoff": { ... tensorcad.geometry.handoff/1, digest included ... },
  "verification": {"passed": true, "failed": [], "checks": 23},
  "solver": {"status": "ready", "executed_here": false, "reason": "..."},
  "limits": [ ... four sentences about what this is not ... ]
}
```

`solver.executed_here` is `false` in every record this repository writes. TensorCAD
measures and hands off geometry; running the flow or the structure is TensorFEM's and
TensorLBM's job, and pretending otherwise would turn a measured design record into a
marketing claim.

The feedback seam for a coupled run is `DesignProblem.solve(feedback=...)`: the callable
receives the current record and returns adjustments to the targets, which are recorded in
the manifest. That is how a measured sinkage or a measured drag correction re-enters the
next design step instead of the loop guessing it.

## API

```python
from tensorcad import DesignProblem, DesignTargets

problem = DesignProblem(
    base=HullParameters(length=12.0, beam=2.0, draft=0.8),
    targets=DesignTargets(
        displacement_m3=12.063418,
        beam_to_draft=2.4 / 0.9,
        wetted_area_m2_max=33.3,
    ),
    variables=("length", "beam", "draft"),
)
solution = problem.solve()
solution.status          # "converged" | "max_iterations" | "no_descent" | "infeasible" | "singular_jacobian"
solution.solved          # converged AND every gate green
solution.failed_checks   # names, so a caller never has to parse a message
solution.convergence_order
solution.sensitivity()
solution.final.handoff   # tensorcad.geometry.handoff/1
write_solver_manifest(solution, "results/solver_manifest.json")
```

## What this does not do

- It does not optimise a hydrodynamic objective. The residual is displacement and shape
  ratios; there is no resistance, no stability, no seakeeping, no class rule.
- It does not search a topology or a shape family. One hull family, three dimensions.
- It does not run a solver, and it does not claim a design is *good* — only that it meets
  the stated targets and passes the geometric gates, with the measurements attached.
