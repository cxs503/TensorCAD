# PR #41 cure: merging `agent/m3213` into `main` as a semantic union

**Mission** m3845 · **Repo** TensorCAD · **Branch of record** `agent/m3213` (PR #41 was
`mergeable=false`) · **Merge base** `df31744` · **main merged in** `2ff80b2` (Merge #39)
· **Branch head going in** `292d794` ("m3213: matured mission deliverables", 14 files,
+2543/−4) · **Result** `3ff58b1`.

## 1. Why the PR went dirty

Both sides filled the same 0.21.0 release with *additive* work, so every conflict is an
"append at the same anchor" conflict rather than a disagreement:

| side | what it added in 0.21.0 |
|---|---|
| `main` (12 commits past the base) | waterplane second moments on the `Hydrostatics` record, the transom-stern hull family (`tensorcad/transom.py`), PR #39's hull routes |
| `agent/m3213` (1 commit) | finite-heel righting arms (`tensorcad/heel.py`), `GZCurve`, docs/example/tests, `__version__` bump |

`pyproject.toml` merged clean: both sides bumped `0.20.0 → 0.21.0` identically.
`README.md` merged clean with both sections intact (checked line-by-line, see §3).

## 2. Conflict ledger and per-file resolution

| # | file | main side | branch side | resolution |
|---|---|---|---|---|
| 1 | `CHANGELOG.md` | new `0.21.0` / `### Added` section (inertia + transom entries) | a *second* `0.21.0` section without an `### Added` heading (heel entries) | one `### Added` heading, main's four paragraphs verbatim, then the branch's four paragraphs verbatim; nothing reworded, nothing reordered inside a side |
| 2 | `tensorcad/__init__.py` | `"MAX_TRANSOM_FRACTION"` in the constant block; `"generate_transom_hull"` in the callables block; `from .transom import (...)` | `"GZCurve"` at the same anchor; `"gz_curve"`, `"heeled_displacement"`, `"righting_arm"` at the same anchor; `from .heel import (...)` | `__all__` = exact union, each name in its case-insensitive alphabetical slot (`GZCurve` before `MAX_TRANSOM_FRACTION`; `generate_transom_hull` before `generate_wigley` before `gz_curve`; `righting_arm` after `resolve_resolution`). Both `from .heel import` and `from .transom import` blocks were already present from the automatic merge and were not touched |
| 3 | `tests/test_public_surface.py` | `"transom"` added to the layered-module tuple | `"heel"` added to the same tuple | tuple carries both: `… "voxel", "export", "heel", "wigley", "transom", "hull", …`. The `EXPECTED_API` set itself auto-merged to the correct 130 names |

No file needed a behavioural decision (no opposing implementations of one function), so no
side had to be dropped — which is exactly what a mechanical overlay would have done
silently: `git checkout --ours/--theirs` on these three files would have deleted 13 or 10
exports respectively without any test necessarily noticing.

## 3. Verification

### 3.1 Machine-checkable union audit

`artifacts/merge_cure_pr41/verify_merge_union.py` (re-runnable, exit code is the verdict):

```
base   = df31744: 107 names
main   = origin/main: +13 names -> 120
branch = origin/agent/m3213: +10 names -> 117
merged = worktree: 130 names (union 130)
union adds: main-only 13, branch-only 10, both 0
added lines preserved: main 3865/3865, branch 2282/2282

PASS: semantic union -- no marker left, no export lost, no added line dropped
```

It checks three things: (i) no `<<<<<<<`/`>>>>>>>`/`=======` conflict marker anywhere in
the worktree's text files; (ii) `__all__` equals *exactly* the union of base ∪ main ∪
branch (130 names would show up as both a loss and an invention); (iii) every
non-blank line that either side added relative to the base still occurs in the worktree
(3865 for main, 2282 for the branch, 0 missing).

**Negative controls** (a check that cannot fail is not evidence): deleting
`"righting_arm"` from `__all__` → `exit 1, "__all__ lost 1 names: ['righting_arm']"`;
planting a conflict marker in a tracked markdown file → `exit 1, conflict marker in
./docs/…:1`. Restoring the file returns `exit 0`. The first version of the marker scan
false-positived on docstring underlines (`=======` in `hydrostatics_analytic.py`), so the
scan now requires the angle-bracket markers (which git labels) and only then treats a
bare `=======` line as a separator — a real trap, recorded in the script comment.

### 3.2 Test parity against the main baseline

Both runs: `pytest -q -n 12 --ignore=tests/test_cli_neural.py
--ignore=tests/test_neural_pipeline.py` (the two ignored modules need the optional
`safetensors` dependency, absent in this image on *both* sides; see §3.3), clean checkout
of `origin/main` extracted by `git archive` for the baseline.

| tree | collected | passed | failed | skipped | wall clock |
|---|---|---|---|---|---|
| `main` @ `2ff80b2` (baseline) | 831 | 818 | 2 | 11 | 285.7 s |
| merged @ `3ff58b1` | 861 | **848** | 2 | 11 | 122.8 s |
| Δ (merged − main) | **+30** | **+30** | 0 | 0 | — |

The 30 extra tests are exactly `tests/test_heel.py` (30 tests, the branch's own suite);
no test that main has is skipped or newly failing in the merged tree. The merged tree was
also run sequentially (`848 passed, 2 failed, 11 skipped in 395.9 s`), i.e. the result does
not depend on the runner.

The two failures are identical in both trees and are environment-only:

```
FAILED tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol
FAILED tests/test_geometry_field.py::test_the_learned_field_gradient_is_autograd_and_points_outward
E   ModuleNotFoundError: No module named 'safetensors'
```

Run in isolation on the untouched baseline checkout: `2 failed, 47 passed` — same two
tests. They are pre-existing on `main`, not introduced by this merge.

### 3.3 Both halves still compute after the merge

The suite is necessary but not sufficient: a merge can leave a module importable yet
wrong. So the two acceptance scripts of the two merged halves were re-run in the merged
tree (post-merge worktree, `PYTHONPATH=.`):

* `examples/gz_curve/run_example.py` (branch half, exits non-zero if any claim fails):
  **12/12 PASS** — barge GZ vs wall-sided `1.036e-15` (gate `1e-9`), waterline residual
  `9.474e-16` m, `GM` from waterplane moments vs `KB+BM−KG` `0.0`, odd symmetry
  `2.220e-16`, shallow barge range of stability `45.95°`, sphere refinement `4.018x` per
  doubling.
* `examples/transom_stern/run_example.py` (main half): runs clean, all four truncations
  `watertight=True manifold=True euler=2` (`tau = 0, 0.05, 0.2, 0.35`).

`docs/merge_cure_pr41.png` (generated by `artifacts/merge_cure_pr41/make_figures.py`,
numbers printed by that script — nothing hand-typed here):

| panel | quantity | measured | reference | error |
|---|---|---|---|---|
| (a) heel half | worst \|GZ − wall-sided\| on a 20×6×2 m barge, `KG=2.0`, `GM=0.5 m`, φ ≤ 33.69° | — | `sin φ (GM + BM tan²φ/2)` | **3.886e-16 m** |
| (b) main half | transverse waterplane `I_T` of a faceted circle, `n = 8…128` | `51.68 → 63.566` m⁴ | `πR⁴/4 = 63.617` m⁴ | `1.876e-1 → 8.029e-4`, ratio 3.73/3.93/3.98/4.00 per doubling (`1/n²`) |
| (b) control | `I_T` of a rectangular box, any tessellation | `360.0` m⁴ | `L B³/12 = 360.0` m⁴ | **0.0** |

Panel (a) also shows where the wall-sided formula *stops* being a reference: past the
deck edge (33.69°) the merged code departs from it, which is the branch's documented
behaviour, not a merge artefact.

## 4. Notes, traps, limits

* `git checkout --ours/--theirs` would have been the fast resolution and would have
  silently deleted 10 exports (`heel`) or 13 (`transom`). The union was chosen and then
  *proved* by the line-level audit, because the case's previous mechanical-overlay
  incident is exactly this shape.
* `__all__` keeps the branch's internal ordering choices (`HeeledSection` before
  `HeelError`, `SubmergedBody` before `StabilityError`) — `__all__` order carries no
  semantics here (`tests/test_public_surface.py` compares sets), and re-sorting the list
  would have been an unrelated diff.
* Running the two example scripts regenerates `docs/gz_curve/{results.json,report.md}`
  and writes `outputs/`; those were reverted/removed so the cure diff contains the
  author's published evidence unchanged.
* The example scripts are not `pip install`ed in this image, so they are run with
  `PYTHONPATH=.`; without it they fail at `from tensorcad import __version__` (pre-existing
  packaging expectation, unrelated to the merge).
* Scope: this cure proves *merge integrity* (nothing lost, nothing invented, both halves
  compute, no test regression). It adds no new physics or validation; the accuracy claims
  belong to the two feature sets' own reports.

## 5. Error judgement (§ key quantities vs baseline, per-item "<3%?")

| quantity | baseline | measured | relative error | < 3 %? |
|---|---|---|---|---|
| barge `GZ` vs wall-sided, φ ≤ 33.69° | wall-sided formula | merged `gz_curve` | 1.0e-15 (worst abs `3.89e-16 m`) | **yes** |
| box transverse `I_T` vs `L B³/12` | closed form 360.0 m⁴ | merged `compute_hydrostatics` | 0.0 | **yes** |
| circle transverse `I_T` vs `πR⁴/4` (n=128) | closed form 63.617 m⁴ | merged `compute_hydrostatics` | 8.03e-4 | **yes** |
| merged suite vs main suite | 818 passed / 2 failed / 11 skipped | 848 passed / 2 failed / 11 skipped | 0 regressions (Δ perf: +30 tests) | **yes** |
| export surface | main 120 ∪ branch 117 | 130 = 107 + 13 + 10 | 0 lost, 0 invented | **yes** |

No quantity in this task lacks a baseline: every claim above is against a closed form, the
main checkout, or a set-theoretic identity.

## 6. SHAs of record

```
old agent/m3213 head : 292d794  m3213: matured mission deliverables
merge result         : 3ff58b1  Merge main into agent/m3213 (cure PR #41): keep both 0.21.0 feature sets
```