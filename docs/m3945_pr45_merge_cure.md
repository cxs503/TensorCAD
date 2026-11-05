# PR #45 cure — semantic union of `main` and `agent/m3391` (m3945)

Repo `TensorCAD`. Pull request #45 (`agent/m3391`, mission m3391: large-angle stability /
righting-arm curve) was reported `mergeable=false` against `main`. This document records
the state before the cure, the conflict inventory with the per-file resolution, and the
numbers that show **both** sides still work on the merged tree.

## 1. State before the cure

| ref | SHA | note |
|---|---|---|
| `origin/main` (merge target) | `2ff80b2` | “Merge #39 (m2925) via merge-governor” |
| `origin/agent/m3391` (PR #45 head, old) | `7115535` | one commit on top of `e604f02` |
| merge base | `e604f02` | “Merge #43 (m3313)”; branch 1 commit ahead, `main` 10 commits ahead |

`git merge origin/main` into the PR branch conflicted in exactly one file
(`CHANGELOG.md`); six more files were touched by both sides and auto-merged, so each was
diffed against both parents by hand before the result was accepted — an auto-merge can be
textually clean and still drop one side's edits ('additive union' is not something git
guarantees).

What each side brought since `e604f02`:

* **`agent/m3391`** — `tensorcad/gz.py` (heeled hydrostatics, `righting_arm_curve`,
  free-surface correction, `with_freeboard`), `docs/gz.md`, `examples/gz/`, `tests/test_gz.py`,
  14 new public names, version bump `0.20.0 → 0.21.0`.
* **`main`** — waterplane second moments on the `Hydrostatics` record, the transom-stern
  hull family (`tensorcad/transom.py`), `docs/transom_stern.md`, `docs/asymmetric_hull*`,
  `tools/`, new tests, version bump `0.20.0 → 0.21.0`.

## 2. Conflict inventory and per-file resolution

| file | kind | resolution (union — nothing dropped) |
|---|---|---|
| `CHANGELOG.md` | **textual**, 1 hunk spanning both releases | both sides prepended a `## 0.21.0` section for *different* features. Kept one release section: `main`'s entry verbatim (`### Added`, waterplane inertias + transom family), then the branch's entry demoted to `### Large-angle stability: the righting-arm curve, and two exact benchmarks` — its prose and all six bullets unchanged. Single `## 0.21.0` heading, `## 0.20.0` header untouched below. Verified line-by-line against both parents by `artifacts/m3945/verify_merge.py` (0 missing lines). |
| `pyproject.toml` | auto (both bumped) | both sides wrote `version = "0.21.0"` — *byte-identical*, so nothing to choose; the merged file equals both parents. No third version number invented. |
| `tensorcad/__init__.py` | auto (both added imports + `__all__` entries) | verified by AST against both parents: every name of either parent is exported exactly once. Merged surface = 134 names = exact set union of the parents' `__all__` (`artifacts/m3945/merged_tree_evidence.py`, check 1). |
| `tests/test_public_surface.py` | auto (both extended `EXPECTED_API`) | diffed against `main`: the only change is `agent/m3391`'s 14 names, inserted as one block ahead of `main`'s inertias/transom names — `main`'s assertions untouched, none duplicated. |
| `tests/test_hull.py`, `tensorcad/cli.py`, `tensorcad/hull.py`, `tensorcad/hydrostatics.py` | auto (`main` only) | branch does not touch them; kept `main` side bit for bit (`git diff origin/main -- <these>` empty). |
| `README.md` | auto (`main` only) | kept `main`'s wording (the branch did not edit the README). |
| `docs/gz.md`, `examples/gz/*`, `tests/test_gz.py`, `tensorcad/gz.py` | branch only | kept verbatim, i.e. added to the tree. |

No file needed a semantic rewrite: the two feature sets are disjoint in code and meet only
in the changelog and the export list, where a union is the correct answer. That is a
property of *this* pair of PRs, not something assumed — §3 tests the joined tree.

## 3. What the joined tree computes

`artifacts/m3945/merged_tree_evidence.py` (run it; exits non-zero if a claim breaks)
re-measures both sides on the merged tree against closed forms:

* **branch side** — box barge `10 × 2 × 2 m`, draft `0.8 m`, `KG = 0.5 m`, `GZ(φ)` over
  0–35° against `GZ = sin φ (GM₀ + ½BM₀ tan²φ)`; then heeled past the form's 38.66° limit to
  prove the agreement is not vacuous (it must *stop* agreeing — and does);
* **main side** — waterplane second moments of the same barge against `L B³/12`, `B L³/12`;
  the transom-stern family at τ = 0, 0.1, 0.2, 0.3 against `analytic_transom_hydrostatics`.

| φ [deg] | GZ merged tree [m] | GZ closed form [m] | ΔGZ [m] |
|---|---|---|---|
| 0 | -0.000000000000 | 0.000000000000 | -2.602e-18 |
| 5 | 0.027738300466 | 0.027738300466 | +2.082e-17 |
| 10 | 0.056113366873 | 0.056113366873 | +3.469e-17 |
| 15 | 0.085830691654 | 0.085830691654 | -6.939e-17 |
| 20 | 0.117745730766 | 0.117745730766 | -5.551e-17 |
| 25 | 0.152973972824 | 0.152973972824 | +5.551e-17 |
| 30 | 0.193055555556 | 0.193055555556 | +5.551e-17 |
| 35 | 0.240219857585 | 0.240219857585 | -2.776e-17 |

| φ [deg] | GZ merged tree [m] | wall-sided form [m] | difference [m] |
|---|---|---|---|
| 40 | 0.297092478 | 0.297836736 | -0.0007 |
| 45 | 0.353553391 | 0.371231060 | -0.0177 |
| 60 | 0.489957146 | 0.815507255 | -0.3256 |

| τ | V mesh [m³] | V closed form [m³] | ΔV/V | ΔA_wp/A_wp | ΔLCB/LCB | Cb mesh / form |
|---|---|---|---|---|---|---|
| 0 | 7.446554 | 7.512698 | -0.8804 % | -0.1749 % | +0.0000 % | 0.465410 / 0.469544 |
| 0.1 | 7.974741 | 8.039094 | -0.8005 % | -0.0944 % | -0.0627 % | 0.498464 / 0.502443 |
| 0.2 | 8.216624 | 8.281656 | -0.7852 % | -0.0791 % | -0.0651 % | 0.513539 / 0.517604 |
| 0.3 | 8.209587 | 8.273675 | -0.7746 % | -0.0683 % | -0.0653 % | 0.513199 / 0.517105 |

| quantity | merged tree [m⁴] | closed form [m⁴] | relative |
|---|---|---|---|
| I_T = L B³/12 | 6.666666667 | 6.666666667 | -1.33e-16 |
| I_L = B L³/12 | 166.666667 | 166.666667 | +5.12e-16 |

Reading of the tables:

* `GZ` reproduces the wall-sided closed form to **6.9e-17 m** worst case (round-off; the
  barge is a polygonal mesh, not an approximation to one), with the displacement residual
  held at `2.2e-16` relative at every heel. Past 38.66° the two diverge as they must
  (−0.7 mm at 40°, −17.7 mm at 45°, −326 mm at 60°), which is the non-vacuity control.
* The two waterplane inertias are exact to **1.3e-16 / 5.1e-16** relative on the box, so
  `main`'s flux-identity route survived the merge unchanged.
* The transom family is within **0.88 %** of its analytic reference on the default 41×13
  grid (τ = 0 reproduces the number `main` documents, −0.88 %), and stays watertight for
  every τ; the error *decreases* as τ grows, as expected because the plate replaces the
  pointed tip that the inscribed chords under-integrate.

## 4. Test evidence

Both trees ran the same command, same interpreter, same machine, serial:

```
python3 -m pytest -q --ignore=tests/test_cli_neural.py --ignore=tests/test_neural_pipeline.py
```

| tree | commit | passed | failed | skipped | wall clock |
|---|---|---|---|---|---|
| `main` (baseline worktree, clean) | `2ff80b2` | 818 | 2 | 11 | 471.6 s |
| merged (`main` ∪ `agent/m3391`) | see §5 | 844 | 2 | 11 | 458.1 s |

The delta is **+26 passed, 0 new failures, 0 new skips** — exactly `tests/test_gz.py`'s 26
tests, the branch's whole test file. The same two failures occur on both trees, and on
`main` they are unrelated to the merge: `tests/test_geometry_field.py::test_the_learned_
revision_speaks_the_same_protocol` and `::test_the_learned_field_gradient_is_autograd_and_
points_outward` fail with `ModuleNotFoundError: No module named 'safetensors'` — the same
missing optional dependency that makes `tests/test_cli_neural.py` and
`tests/test_neural_pipeline.py` uncollectable, which is why both are excluded from the run
above (they are excluded identically on the baseline, so the comparison is fair).

Logs: `runs/m3945/pytest_merged.log`, `runs/m3945/pytest_main_baseline.log`.

## 5. Result

* Old PR head `agent/m3391` = `7115535`; new head = the merge commit below (`git log -1`),
  i.e. `origin/main@2ff80b2` plus the resolution described in §2.
* `git diff origin/main --stat` on the merged tree touches only the branch's own files plus
  `CHANGELOG.md`, `tensorcad/__init__.py`, `tests/test_public_surface.py` — no `main`-only
  file is modified, so the merge is additive with respect to `main`.
* No force-push, no history rewrite: the merge commit has both `7115535` and `2ff80b2` as
  parents, so `main` can be fast-forwarded onto it by the merge governor.

## 6. Notes and limits

* The version number is left at `0.21.0` because both sides chose it independently and the
  merged release does contain both feature sets; bumping again would have been an invented
  number with no release behind it.
* This cure verifies *hydrostatics and stability numerics*. It says nothing about whether
  the righting-arm curve is fit for a statutory stability check: `docs/gz.md` is explicit
  that damaged / free-surface-open cases and the IMO weather criterion are outside the
  module, and the merge changes none of that.
* `safetensors` is not installed in this container, so the four neural test modules cannot
  run here at all; the parity claim in §4 is a parity claim, not a statement that the
  neural path is healthy.