# PR #52 re-opened --- the third `GZ` route merged under the current `main` (m4534)

Repo `TensorCAD`. PR #52 (`agent/m3581`, "matured mission deliverables") was closed by the
owner on 2026-11-05 instead of merged; the task was re-opened as m4534 with the instruction
to redevelop the feature on the current `main`. This page records why the branch could not
be merged as it stood, what the redevelopment changed, and the measured numbers that show
all three righting-arm routes working on the merged tree. `agent/m3581` stays in place as the
archive.

## 1. State before the redevelopment

| ref | note |
|---|---|
| `origin/main` | `faa0da6`, "Merge #156 (m4092) via merge-master"; version 0.22.0 |
| `origin/agent/m3581` (PR #52 head) | `b691ce4`, one commit, cut from 0.20.0 |
| merge base | before `gz.py` (0.21.0) and `righting_arm.py` (0.22.0) landed |

PR #52 was a **second semantic duplicate**, the same failure mode that closed PR #40 (see
`docs/m4038_pr40_merge_cure.md`). Its `tensorcad/righting.py` publishes `righting_arm_curve`,
`RightingArmCurve` and `wall_sided_righting_arm`, and by review time `main` bound all three of
those names to other modules:

| name | `main` at 0.22.0 |
|---|---|
| `righting_arm_curve` / `RightingArmCurve` | `tensorcad/gz.py` --- divergence theorem on the rotated closed surface, waterline by Newton on the projected wetted area |
| `bisection_righting_arm_curve` / `BisectionRightingArmCurve` | `tensorcad/righting_arm.py` --- per-face clipped sections, waterline by bracket-and-bisect |
| `wall_sided_righting_arm` | `tensorcad/gz.py` |

Merging as delivered would have bound the third route's names over the first route's, leaving
`gz` unreachable from the package while `tests/test_righting.py` kept passing on a direct
`from tensorcad.righting import ...`. That is a silent loss, not an import error, which is why
the cure is a naming decision plus a test rather than a merge-conflict resolution.

## 2. What the redevelopment did

* **Cut from current `main`** (`faa0da6`) and carried the branch's eleven files forward:
  `tensorcad/righting.py`, `tests/test_righting.py`, `examples/righting/run_example.py`,
  `docs/righting.md` and `docs/righting/`, plus the README/CHANGELOG entries, rewritten
  against 0.23.0 and re-measured on this tree.
* **No route dropped; the overlap became the check.** The third route is published as
  `section_righting_arm_curve` / `SectionRightingArmCurve`; `submerged_body`, `ImmersedBody`,
  `RightingArmPoint` and `RightingError` do not collide and keep their names. The table is
  asserted in both directions by
  `tests/test_righting.py::test_package_aliases_point_here_and_do_not_shadow_the_other_routes`.
* **Two new cross-route tests** (43 in total, the branch's 40 plus these): the barge
  three-route agreement, and the hull agreement against `gz` with a `KG`-shift **negative
  control** so the comparison is shown to be able to fail.
* **`tools/verify_merge_m4534.py`** writes `artifacts/m4534/merge_verification.json` and a
  three-panel figure; the numbers below come from that file, not from this text.
* **Version** 0.23.0 (`pyproject.toml`, `tensorcad.__version__`, `tests/test_public_surface.py`),
  and four pre-existing lint debts on `main` (trailing newlines, an unused import, an unused
  local, a bare `zip`) fixed so `ruff check .` is clean.

## 3. Numbers on the merged tree

| check | measured |
|---|---|
| barge, three routes, worst `\|dGZ\|` over 0-45 deg | `7.915e-14 m` |
| barge vs wall-sided closed form, worst relative | `2.360e-14` |
| `41x13` hull at `z=-0.3 m`, station route vs `gz`/Newton, worst `\|dGZ\|` | `7.299e-07 m` |
| station-area volume vs closed-surface volume (same cap) | `0.0` relative |
| 40 vs 80 station intervals, `GZ_max` relative | `4.804e-04` |
| negative control: `KG + 5 mm` moves `GZ(20 deg)` | `1.710e-03 m` (2343x the inter-route gap) |

The barge agreement is round-off because a box is a polyhedron and each of the three routes is
exact on it; the hull gap (`7.3e-07 m` on arms of order
`1e-2 m`, i.e. about `1e-5` relative) is **not** round-off and is not claimed to be --- the two
routes discretise the body differently, a station grid against a triangle mesh, and the gap is
four orders of magnitude inside the 15 % band that counts as unacceptable.

The example's own claims, re-run on this tree (`examples/righting/run_example.py`, exits
non-zero if one breaks; table lifted programmatically from `docs/righting/results.json`):

| claim | measured | limit | unit | result |
|---|---|---|---|---|
| barge GZ matches the wall-sided closed form (0-45 deg) | `5.770e-15` | `1e-09` | relative | PASS |
| barge dGZ/dphi at 0 recovers KB + BM - KG | `5.820e-08` | `1e-06` | relative | PASS |
| barge deck-edge immersion occurs at the classic angle | `0.000e+00` | `1e-09` | deg | PASS |
| station-area volume equals the closed-surface volume | `0.000e+00` | `1e-09` | relative | PASS |
| hull dGZ/dphi at 0 recovers KB + I_T/V - KG | `5.546e-05` | `2e-03` | relative | PASS |
| hull curve is grid-converged (40 vs 80 intervals) | `2.159e-04` | `5e-03` | relative | PASS |
| hull curve vanishing angle is grid-converged (40 vs 80 intervals) | `2.852e-04` | `1e-02` | relative | PASS |
| zero-freeboard hull is rejected rather than answered | `0.000e+00` | `0e+00` | bool | PASS |

## 4. Scope carried over, unchanged

Free-surface corrections, tanks, flooding, damage stability, trim/sinkage coupling and any
statutory criterion (IS Code, SOLAS) are outside the module. Sections are assumed
single-valued in `z` per `y` (no overhangs, no immersed deck area), the mesh must carry
freeboard above the upright waterline (a capped body is rejected, not answered), and past
`deck_immersion_onset_deg` the curve describes a watertight box rather than a ship. A `GZ`
computed here is geometry, never compliance.
