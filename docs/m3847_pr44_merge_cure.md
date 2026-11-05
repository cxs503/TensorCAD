# PR #44 cure — semantic union of `main` and `agent/m3347` (m3847)

Repo `TensorCAD`. Pull request #44 (`agent/m3347`) was reported `mergeable=false`
against `main`. This document records the conflict inventory, the per-file resolution,
the test evidence, and the cross-feature numbers that show the merged tree still runs
both sides' physics. Every number below is produced by `tools/report_m3847.py`
(inventory + experiments + figure) and `pytest`, and stored in
`artifacts/m3847/{merge_inventory.json,numbers.json,pytest_main.log,
pytest_merged_full.log}` — nothing in the tables is typed by hand.

## 1. State before the cure

| ref | SHA | note |
|---|---|---|
| `origin/main` | `2ff80b2` | “Merge #39 (m2925) via merge-governor” |
| `origin/agent/m3347` (PR #44 head, **old**) | `9f1572f` | “m3347: matured mission deliverables” |
| merge base | `e604f02` | “Merge #43 (m3313) via merge-governor” |

`main` had moved 10 commits past the merge base (transom-stern hull family, waterplane
second moments on the `Hydrostatics` record, asymmetric fore/aft hull exponents,
curvature rewrite); the branch had 1 commit (righting-arm / `GZ` curve, wall-sided
`freeboard=` topsides). No commit on either side is dropped: the cure is a true merge,
not a rebase or an overlay, so both sides' history stays reachable from the PR head.

**New SHA after the cure (merge commit): `5fe2912`.**

## 2. Conflict inventory and per-file resolution

41 paths were touched by one side or both since `e604f02`: 31 by `main` only, 6 by the
branch only, 4 by both. `git merge` raised **one** textual conflict
(`CHANGELOG.md`); the other three shared files auto-merged and were then checked for
semantic loss by symbol-level comparison against *both* parents (`merge_inventory.json`,
`single_sided_content_lost: []`).

| file | conflict | resolution (union, nothing dropped) |
|---|---|---|
| `CHANGELOG.md` | both sides opened the same `## 0.21.0 — 2026-11-04` section (branch: 3 bullets on `righting_arm_curve`/`freeboard=`; main: `### Added` with 6 bullets on waterplane inertia + transom family) | kept main's `### Added` heading and **both bullet sets verbatim** — branch's GZ/freeboard/verification bullets first, then main's six. Verified line-by-line: every non-empty line of both parents' section text is present in the merged file, and no conflict marker survives. |
| `tensorcad/__init__.py` | auto (both extended the export surface) | union: **125** exported names = 112 (branch) ∪ 120 (main), 107 shared. Nothing from either `__all__` is missing, no duplicates; `GZCurve/GZPoint/SCHEMA_GZ/DEFAULT_HEEL_ANGLES/righting_arm_curve` **and** `generate_transom_hull/TransomModel/MAX_TRANSOM_FRACTION` all import from the package root. |
| `tests/test_public_surface.py` | auto (both extended `EXPECTED_API`) | union: the 125-name set is asserted once each; `transom` also joined the layered-module import list. |
| `pyproject.toml` | auto (both bumped the version) | both blobs were the *same* `0.21.0` bump, so the merged file is byte-identical to both — one version header, consistent with the merged CHANGELOG. |
| the other 37 paths | none | single-sided and taken byte-for-byte from the side that changed them (verified in `merge_inventory.json`; `main`'s 31 files are byte-identical to `2ff80b2`, the branch's 6 to `9f1572f`). |

A mechanical overlay (take one side for the shared files) would have silently dropped,
depending on the direction, either the whole GZ/`freeboard` feature or the whole
transom/inertia work — including their CHANGELOG entries and their export-surface
assertions, which is exactly the failure mode the cure was asked to avoid.

## 3. Cross-feature verification on the merged tree

The seam of this merge is real: the branch's GZ solver was developed against a `main`
that had neither the transom hull family nor the waterplane-inertia record, and it
consumes both (`freeboard=` topsides, waterplane area/inertia). So the check is not
“does each side still pass its own tests” but “does the branch's solver run on main's
hull families and agree with independent references”.

### 3.1 Branch solver vs the classical wall-sided formula (box `2 × 0.8` m, freeboard 1.2 m)

`KG` is chosen so the metacentre sits half a `BM` above `G`. Wall-sided window:
`tan φ < 2T/B` ⇒ φ < 38.66°.

| quantity | value |
|---|---|
| `GM` from the curve | **0.233333 m** |
| `GM` from the closed form `B²/12T − (KG − T/2)` | **0.233333 m** (relative 2e-16) |
| worst relative error inside the window (0…30°) | **1.11e-08** |
| first relative error just outside the window (45°) | **5.66 %** (asserted departure, not tolerance) |

### 3.2 Branch solver on main's hull families

| body | `GM` (`initial_stability`) | `GM` (curve) | rel. error | `GZ_max` | heel of `GZ_max` | vanishing angle | displacement volume spread over 0…90° |
|---|---|---|---|---|---|---|---|
| transom hull `τ=0.20, p=0.70`, `KG=0.8`, `freeboard=0.3` | 0.157521898 m | 0.157521898 m | **0.0** | 0.039198 m | 20° | 41.44° | **0.0 m³** (exactly conserved) |
| Wigley + 0.3 m topsides, `KG=0.264` | 0.264014 m | 0.264014 m | **0.0** | 0.277458 m | 90° | none | — |

Displacement being bit-identical at every heel is the branch's own invariant: the plane
solve re-establishes the unheeled volume to 1e-11 relative, and it survives main's transom
plate and the freeboard extension unchanged.

### 3.3 Main's transom hydrostatics against its closed form (same tree)

| quantity | mesh | closed form | rel. error |
|---|---|---|---|
| displacement volume | 8.216624 m³ | 8.281656 m³ | **0.785 %** |
| waterplane area | 15.515828 m² | 15.528105 m² | **0.079 %** |
| LCB (`x`) | 4.374557 m | 4.377408 m | **0.065 %** |
| block coefficient | 0.513539 | 0.517604 | **0.785 %** |

These are the reference grid (`41 × 13`) errors and are monotone in the grid
(`docs/transom_stern.md`); the point here is only that they are unchanged by the merge.

## 4. Test evidence

Same command on both trees (`pytest tests -q --ignore` the two neural modules, whose
`safetensors` dependency is absent in this image and which fail to collect on `main`
too):

| tree | passed | failed | skipped | log |
|---|---|---|---|---|
| `main` = `2ff80b2` (baseline) | 818 | 2 | 11 | `artifacts/m3847/pytest_main.log` |
| merged cure tree = `5fe2912` | **879** | 2 | 11 | `artifacts/m3847/pytest_merged_full.log` |

The delta is **+61**, exactly `tests/test_gz_curve.py` — the branch's own suite — and
nothing main had regressed. The two failures are identical on both trees and are
pre-existing (`tests/test_geometry_field.py`, learned-geometry protocol/autograd tests
that need optional packages); they are not merge-caused and are not touched here. A
chunked run of the same suite on the cure tree
(`artifacts/m3847/pytest_merged.log`) adds up to the same 879/2/11, so the total does not
depend on the worker count.

The branch's self-checking example also passes on the merged tree:
`PYTHONPATH=. python3 examples/gz_curve/run_example.py` reproduces
`GM 0.321622996 m == initial_stability 0.321622996 m` and the 0.99 % one-degree secant
figure, and `PYTHONPATH=. python3 examples/transom_stern/run_example.py` (main's) prints
watertight/edge-manifold/Euler-2 for every `τ` — both written as assert-and-exit-non-zero
scripts.

## 5. Notes, limits, honest gaps

* **Push not performed from this sandbox.** GitHub was unreachable from the worker
  (`git fetch`/`ls-remote` returned nothing; `git push` is blocked by the task policy),
  so `origin/agent/m3347` still points at the old `9f1572f`. The cured merge commit
  `5fe2912` is on this branch and is also shipped as an incremental bundle
  `artifacts/m3847/agent_m3347_cured.bundle` whose prerequisite is `main` `2ff80b2`:
  `git fetch artifacts/m3847/agent_m3347_cured.bundle m3347-cured:agent/m3347` reconstitutes
  the branch without network. The one command still owed is
  `git push origin 5fe2912:agent/m3347` (or the bundle fetch above, then push).
* `GZ(0)` on the transom hull is **−4.0e-05 m** rather than exactly 0 — a 0.026 %
  offset on `GM`, from the clipped-surface `Y`-moment at exactly zero heel. It is the
  same discretisation residue the branch documents; it is not a merge artefact (the
  branch's box case returns exactly 0 on a symmetric section).
* The transom case above is deliberately marginal (`GM = 0.158 m`, vanishing at 41°).
  Its signed `dynamic_stability_area` over 0…90° is therefore **negative**
  (−0.0551 m·rad): the area aft of the vanishing angle is counted with its sign. That is
  a property of the chosen `KG`, not of the merge; the report quotes it rather than
  quoting a one-sided positive area.
* Only one textual conflict existed. That makes this cure cheap, but it also means the
  risk was in the *auto*-merged files (export surface, CHANGELOG section, version) —
  which is why they were verified by comparison against both parents rather than trusted.

## 6. Error verdict

| key quantity | reference | error | < 3 % ? |
|---|---|---|---|
| box `GM` (curve) | closed form `B²/12T − BG` | 2.2e-16 | **yes** |
| box `GZ` inside wall-sided window | wall-sided formula | 1.11e-08 | **yes** |
| transom `GM` (curve) | `initial_stability` on the same mesh | 0.0 | **yes** |
| Wigley+topsides `GM` (curve) | `initial_stability` on the same mesh | 0.0 | **yes** |
| transom displacement volume | closed-form incomplete-taper integral | 0.785 % | **yes** |
| transom waterplane area | closed form | 0.079 % | **yes** |
| test counts | `main` baseline 818 passed | +61 passed, 0 new failures | **yes** |
| `origin/agent/m3347` updated | required by the task | **not done** (no network) | **no** |