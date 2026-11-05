# PR #48 cure — merging `agent/m3479` into `main` (mission 3871)

## What was dirty

`agent/m3479` (PR #48, large-angle `GZ`) branched from `e604f02`; `main` has since taken
PR #39 (`2ff80b2`: waterplane second moments, split fore/aft taper exponents, transom-stern
hull family). `mergeable=false`.

## Three-way merge and per-file resolution

`git merge origin/main` from the PR tip produced one textual conflict and four auto-merged
files that were audited rather than trusted.

| file | how it merged | resolution |
| --- | --- | --- |
| `CHANGELOG.md` | conflict (both sides wrote a `## 0.21.0` section) | one `## 0.21.0 — 2026-11-04` heading; PR prose kept verbatim under `### Large-angle stability: the GZ curve`, main's list kept verbatim under `### Added`, plus a new `### Fixed` bullet for the integration defect below. Nothing dropped. |
| `README.md` | auto | PR semantics kept (`large-angle GZ` removed from the "does not do" list, pointed at `tensorcad.gz`); main's transom/asymmetric sections untouched. |
| `tensorcad/__init__.py` | auto (both sides appended) | union checked programmatically: PR's eleven `gz` names on top of main's exports, every `__all__` name resolvable (`missing exports: []`). |
| `tests/test_public_surface.py` | auto (both sides appended) | PR's eleven expected names added, main's entries kept; the set-equality test passes. |
| `pyproject.toml`, `__version__` | auto | both sides bumped to `0.21.0`: one release, no renumbering. |
| 28 added files | clean | main's transom/asymmetric/Wigley artefacts and the PR's `gz` module, docs, example, tests all present. |

## Integration defect the union exposed, and the fix

`gz.solve_waterline` started from the bracket midpoint. For a hull generated over
`[0, draft]` that midpoint is exactly the half-draft ring, and the transom family emits a
ring there — so `submerged_solid` refused the trial waterline ("the waterline z=-0.45 m
meets the surface along mesh edges") and `gz_curve` raised on the family's default grid
`nz=13`. Each feature was correct alone; the pair was not. `tensorcad/gz.py` now starts
with `_start_waterline`, which steps off mesh rings by a fixed fraction of the span
(displacement is continuous in the waterline height, so the converged height and its
conservation are unchanged, still ~9e-13 relative). `submerged_solid` still refuses an
on-edge waterline that is *requested*: the documented guard is unchanged. Two regression
tests in `tests/test_gz.py` pin both halves.

Cross-feature evidence: 28 checks in `tools/verify_union_m3871.py` →
`artifacts/m3871/union_verification.{json,md}`.
GZ and `initial_stability` GM agree to 2e-16…4e-15 m on the transom mesh (tau = 0, 0.25,
0.4); displacement held to 9e-13 relative at every heel; `GZ(0)` equals the mesh's own
lateral centre of buoyancy; on the split-exponent asymmetric hull `GZ(0) = -4.6e-17 m` and
`LCB` mesh-vs-reference 9.3 mm (0.15 %).

## Test evidence

| run | result |
| --- | --- |
| merged tree, full suite (minus two files that fail collection for want of `safetensors`) | **2 failed, 849 passed, 11 skipped** in 1257 s |
| `main` at `2ff80b2`, same suite | the same two `tests/test_geometry_field.py` failures (reproduced in isolation: `2 failed, 47 passed`) — pre-existing, not from the merge |
| merged tree after the start fix: `test_gz`+`test_transom`+`test_public_surface`+`test_cli`+`test_hull`+`test_hydrostatics_inertia`+`test_asymmetric_hull` | **239 passed** |
| `tools/verify_union_m3871.py` | 28/28 checks pass |

## Branch

- old `agent/m3479`: `55f6e7d`
- merged tip: the commit carrying this file (fast-forward from `55f6e7d`; first parent the PR
  tip, second parent `2ff80b2`).

## Limitations carried, not introduced

- The transom mesh is mirror-symmetric in its vertices but not in its triangulation:
  `|y_B| = 6.2e-5 m` at the design waterline versus `6.6e-17 m` for the plain hull, so
  `GZ(0)` there is 6.2e-5 m instead of exactly zero (0.0026 % of beam, far inside the
  family's 0.045–0.88 % discretisation error).
- `tensorcad.analytic_hydrostatics` is the breadth-exponent-only reference; the
  split-exponent one is `tensorcad.hull.analytic_hydrostatics`. One name, two functions — a
  pre-existing collision, left as is.
