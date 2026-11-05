# PR #40 cure — semantic union of `main` and `agent/m3169` (m3923)

Repository `TensorCAD`. Pull request #40 (head branch `agent/m3169`) was reported
`mergeable=false` against `main`: both sides had moved after the last common commit, and
two of the four files touched by both sides conflicted. This document is the merge
record — the state before the cure, the per-file resolution, the machine-checked
preservation audit, and the test evidence that neither side's capability was lost.

## 1. State

| ref | SHA | note |
|---|---|---|
| `origin/main` | `2ff80b2` | “Merge #39 (m2925) via merge-governor” (PR #39 landed) |
| `origin/agent/m3169` (PR #40 head, **old**) | `fcf15ad` | “m3169: matured mission deliverables” — 1 commit |
| merge base | `df31744` | PR #38 merge; `main` is 12 commits ahead, the branch 1 |
| cured branch head (**new**) | `2e5e934` | “Merge main (2ff80b2, PR #39) into agent/m3169 (cure PR #40)” |

The branch brings the **righting-arm curve** (`tensorcad/righting_arm.py`, 674 lines,
41 tests, one example, three figures); `main` brings the **transom-stern hull family**,
the **waterplane second moments** on the `Hydrostatics` record, the asymmetric hull and
the curvature rewrite. The union is therefore a plain superset at module level — the
only genuinely shared contracts are the export list, the export-set test, the CHANGELOG
and `pyproject.toml`.

## 2. Conflict inventory and per-file resolution

`git merge --no-ff origin/main` from `fcf15ad` produced exactly two conflicted files.
No resolution below was produced by taking one side wholesale: each keeps every line
both sides added.

| file | conflict | resolution (union — nothing dropped) |
|---|---|---|
| `CHANGELOG.md` | both sides prepended bullets under `## 0.21.0`; the branch wrote them as a bare list, `main` under a `### Added` heading | one `### Added` section: `main`'s heading + its transom/inertia/Wigley entries verbatim, followed by the branch's four righting-arm bullets verbatim. No entry reworded, none removed, no duplicate heading. |
| `tests/test_public_surface.py` | both edited the layered-module tuple in `test_layered_modules_are_all_importable` | union of the two module lists: `"stability", "righting_arm", "bonjean",` (branch) **and** `"wigley", "transom", "hull", …` (`main`). Assertion strength unchanged (it still imports every module). |
| `tensorcad/__init__.py` | auto-merged (both appended exports) | verified, not assumed: the branch adds 10 names, `main` adds 13, and the merged `__all__` (130 names) contains **all 23** — checked by parsing both sides' `__all__` and diffing against the merged file (0 missing). |
| `pyproject.toml` | auto-merged | both sides bumped the version to `0.21.0`; the merged file has a single `version = "0.21.0"`, consistent with the merged CHANGELOG heading. |
| 43 other paths | one side only | byte-identical to that side's blob (verified by `git hash-object` against `rev-parse <ref>:<path>`, 0 violations). |

## 3. Preservation audit (machine-checked, not asserted)

`tools/verify_merge_m3169.py` recomputes the audit from git plumbing, so the claim does
not rest on the merge having been done carefully:

| quantity | value |
|---|---|
| paths touched by either side | 47 |
| files missing from the merged tree | **0** |
| deletions declared by either side | **0** |
| paths touched by both sides | 4 (`CHANGELOG.md`, `pyproject.toml`, `tensorcad/__init__.py`, `tests/test_public_surface.py`) |
| paths touched by one side only | 43, blob-identity violations **0** |
| added lines checked for presence in the merged file (both sides) | 6155 |
| added lines lost | **0** |

This is the exact failure mode the mission warns about (a mechanical overlay silently
dropping one side): the audit is per *added line*, so a content loss would have to hide
inside a reworded line to survive it — and no line was reworded.

## 4. Test evidence

Same command on both trees (`python3 -m pytest -q`, the two neural modules that cannot
import in this environment — `tests/test_cli_neural.py` (`safetensors` missing) and
`tests/test_neural_pipeline.py` — excluded identically on both sides):

| tree | passed | failed | skipped | wall |
|---|---|---|---|---|
| `main` `2ff80b2` (baseline) | 818 | 2 | 11 | 472 s |
| merged `2e5e934` | **859** | 2 | 11 | 823 s |

The two failures are the **same two tests on both trees**, both environment-related
(`tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol`,
`…::test_the_learned_field_gradient_is_autograd_and_points_outward`). The merged tree is
therefore not below the baseline: +41 passing tests, an identical failure set.

Merge-critical modules on the merged tree, run separately:
`tests/test_righting_arm.py tests/test_transom.py tests/test_public_surface.py` →
**143 passed in 40 s**. `tests/test_merge_m3169.py` (new, 4 tests) pins the composition
described in §5.

## 5. What the union actually does — and one honest finding

Measured by `tools/verify_merge_m3169.py` (`artifacts/m3923/merge_verification.json`,
figure `artifacts/m3923/merge_verification.png`):

1. **The branch's route is alive.** A wall-sided barge (`L=10, B=2, T=0.8`, freeboard
   0.4, `KG=0.55` above keel) reproduces `GZ = sin φ (GM + ½ BM tan²φ)` with
   **max |error| = 1.2e-14 m** inside the formula's validity window — the binding limit
   there is deck-edge immersion at `atan(2·freeboard/B) = 21.8°`, not the 38.7° wall-sided
   window, so only heels below 21.8° are compared. A floating sphere reproduces
   `GM sin φ` to **1.2e-4 m** at resolution 32.
2. **`main`'s route is alive.** The transom family matches its closed-form taper
   integrals: **max |Δ∇/∇| = 0.88 %** over `τ = 0, 0.2, 0.4` at the default 41×13 grid
   (inscribed chords, so the error is one-signed). The `Hydrostatics` record's waterplane
   second moments reproduce `L B³/12` and `B L³/12` to **1.7e-16** relative on two
   rectangular barges, and agree with `stability.waterplane_moments` to the last bit.
3. **The two routes compose.** On a transom hull the heeled integrator
   (`heeled_displacement`) and the hydrostatics record (`compute_hydrostatics`) give the
   same displacement (**6.9e-13** relative) and the same buoyancy-centre height
   (**2.2e-13 m**); the solved displacement is held at every heel (**residual 6.9e-13**)
   and `kg` enters exactly as a lever, `dGZ/dkg = −sin φ` (**1.4e-17 m**).
4. **Finding (reported, not smoothed over).** On a *fine* hull the small-angle slope from
   the heeled solver is `KB − KG`, while the record predicts `KB + I_T/∇ − KG`; the
   difference is **exactly `BM`** (parametric hull: slope 0.48301 vs 0.92162, absolute
   gap 0.43862 m = `BM`; transom hull `τ=0.2`: 0.48761 vs 0.95752, gap 0.46991 = `BM`).
   The barge of §5.1 shows the solver is right when the waterline does not coincide with
   a tessellation ring; the parametric hull is capped exactly at the design waterline, and
   m3169's own documentation already records a `waterplane_moments` failure for a
   waterline that lands on a vertex ring. The likely cause is therefore that degenerate
   cap/ring coincidence, **but this is a hypothesis, not a verified diagnosis**. The merge
   touched neither module, so the effect predates the union and is a cross-module
   question for the two owners, not a conflict-resolution defect — hence verdict
   `WITHHELD` on "GM from the heeled curve equals GM from the record on a fine hull",
   with both numbers published above rather than a weakened assertion.

## 6. Publication

Branch `agent/m3169` old head `fcf15ad` → new head `2e5e934`, pushed so PR #40 becomes
mergeable; the merge into `main` is the merge-governor's step, not this mission's.
Reproduce the audit with:

```
python3 tools/verify_merge_m3169.py     # exits non-zero if any claim above breaks
python3 -m pytest -q tests/test_merge_m3169.py
```