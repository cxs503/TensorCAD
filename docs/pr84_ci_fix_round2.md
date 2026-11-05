# PR #84 CI round 2: the merge-graph conflict that kept the matrix red

Round 1 (`cf10921`) cleared the 26 `ruff` findings. The checks still could not go green
because, after `main` merged **PR #45** (`tensorcad/gz.py`, a righting-arm/GZ module),
`agent/m3858` no longer merged cleanly into `main`:

```
$ git merge-tree --write-tree HEAD origin/main
CONFLICT (content): Merge conflict in tensorcad/__init__.py
```

GitHub cannot build the pull-request merge commit for a conflicted PR, so the `tests`
matrix never reaches a green verdict: whichever leg fails first cancels its siblings,
which is exactly the `test (3.13)=failure, 3.10/3.11/3.12=cancelled` pattern.

## What collided

Both sides added a righting-arm implementation and exported it from the package root:

| side | module | public names |
|---|---|---|
| PR #84 (`agent/m3858`) | `tensorcad/righting_arm.py` | `DEFAULT_HEELS_DEG`, `HeeledDisplacement`, `RightingArm`, `RightingArmCurve`, `RightingArmError`, `heeled_displacement`, `righting_arm`, `righting_arm_curve`, `wall_sided_righting_arm`, `wall_sided_validity_deg` |
| `main` (PR #45) | `tensorcad/gz.py` | `DEFAULT_ANGLES_DEG`, `FreeSurfaceCorrection`, `GzError`, `HeeledHydrostatics`, `RightingArmCurve`, `box_righting_arms_analytic`, `box_wall_sided_limit_deg`, `deck_edge_immersion_deg`, `free_surface_correction`, `heeled_hydrostatics`, `rectangular_tank_inertia`, `righting_arm_curve`, `wall_sided_righting_arm`, `with_freeboard` |

Three names are genuinely ambiguous at the package root (`RightingArmCurve`,
`righting_arm_curve`, `wall_sided_righting_arm`). Both modules are self-contained (neither
imports the other), so the duplicate implementations can coexist; only the *root exports*
need a single owner.

## Resolution (union, no implementation deleted)

* `main` keeps ownership of the three ambiguous root names. `tensorcad/gz.py` is already
  merged and its `tests/test_gz.py` imports it directly; keeping the established public
  names stable is the conservative choice. The PR's versions remain importable from
  `tensorcad.righting_arm`, which is where `tests/test_righting_arm.py` takes them from.
* All non-colliding names from both sides are exported: the `gz` family *and*
  `DEFAULT_HEELS_DEG`, `HeeledDisplacement`, `RightingArm`, `RightingArmError`,
  `heeled_displacement`, `righting_arm`, `wall_sided_validity_deg`.
* `__all__` is the sorted union (141 unique names); `tests/test_public_surface.py` already
  stores `EXPECTED_API` as a set and auto-merged to the same union.

No public symbol was renamed, no implementation was removed, no numeric behaviour changed.
The collision itself (two GZ modules) is a genuine duplication that should be folded into
one module in a follow-up; this fix only makes the merge graph and the public surface
consistent so the PR can land.

## Verification (local, Python 3.11.6; the runner image's 3.13 is not available offline)

| check | result |
|---|---|
| `ruff check .` | All checks passed |
| `mypy` (34 source files, strict) | clean |
| `pytest tests/test_public_surface.py tests/test_gz.py tests/test_righting_arm.py` | **73 passed** |
| `python -c "import tensorcad"` | 141 unique `__all__` names, `RightingArmCurve is tensorcad.gz.RightingArmCurve`, `RightingArm is tensorcad.righting_arm.RightingArm` |

The merge commit is `Merge origin/main into agent/m3858` (PR #84 cure).
