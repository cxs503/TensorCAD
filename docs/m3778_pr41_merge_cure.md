# PR #41 (`agent/m3213`) merge cure — conflict ledger and evidence

Mission 3778. Target: `TensorCAD` PR #41, branch `agent/m3213`, mergeable=false
against `main`.

| item | value |
|---|---|
| merge target (main side) | `2ff80b2` — *Merge #39 (m2925) via merge-governor* |
| branch before cure | `5f46b63` — *m3213: righting arms at finite heel (`tensorcad.heel`, GZ curve)* |
| merge base | `df31744` — *Merge pull request #38 from cxs503/agent/m2877* |
| branch after cure | `2051f08` — *Merge main into agent/m3213 (cure PR #41): keep heel and transom surfaces* |
| conflict files | 3 (`CHANGELOG.md`, `tensorcad/__init__.py`, `tests/test_public_surface.py`) |
| files auto-merged | 5 of the 5 files both sides touched; 2 of them audited by hand (`pyproject.toml`, `README.md`) |
| branch adds to main | 14 files, +2543 / −4 lines |
| main adds since merge base | 38 files (transom-stern family, waterplane inertias, asymmetric hull docs, m2873/m3768 evidence) |
| merge probes | 22/22 pass (`tools/verify_merge_pr41.py`) |

## Why the branch was dirty

`main` moved 38 files after `df31744`; the branch had one commit touching 14 files.
Both sides edited the same five names-only / release-metadata files, so git emitted
three content conflicts — and, more dangerously, two clean textual auto-merges
(`pyproject.toml`, `README.md`) that still had to be read, because "no conflict
marker" is not "no lost content".

## Per-file resolution

### `tensorcad/__init__.py` (conflict, 2 hunks) — **union**

The `__all__` list is the public contract and `tests/test_public_surface.py` asserts
set equality with a literal. The branch added the heel surface, main added the
transom surface; both landed in the same alphabetical slots, so git could only pick
one side per hunk.

* hunk 1 (uppercase block): kept **both** `"GZCurve"` (branch) and
  `"MAX_TRANSOM_FRACTION"` (main), in alphabetical order.
* hunk 2 (lowercase block): kept `"generate_transom_hull"` (main) and
  `"gz_curve"`, `"heeled_displacement"` (branch); `"righting_arm"` was re-inserted
  in its alphabetical slot after `"resolve_resolution"`.

Result: `__all__` has **130 names = |117 branch ∪ 120 main|**, no duplicates, no
missing names, and every name resolves:

```
merged n=130 union n=130  missing=[]  extra=[]  duplicates=[]  unimportable=[]
```

The two new module imports (`from .heel import …`, `from .transom import …`) were in
different regions and auto-merged; both were verified present.

### `CHANGELOG.md` (conflict, 1 hunk) — **union**

Both sides wrote their entries under the same `## 0.21.0 — 2026-11-04` heading. The
branch wrote flat bullets; main wrote an `### Added` block. Resolution keeps the
branch's "Righting arms at finite heel" bullets **verbatim**, followed by main's
`### Added` block **verbatim** (waterplane second moments, transom-stern family,
hydrostatic references, example). No bullet was reworded or dropped.

### `tests/test_public_surface.py` (conflict, 1 hunk) — **union**

The layered-module import list: branch appended `"heel"`, main inserted `"transom"`
before `"hull"`. Resolution lists both:
`… "hydrostatics_analytic", "voxel", "export", "heel",` / `"wigley", "transom", "hull", …`.
The surrounding `EXPECTED_API` literal (the big set the test compares against)
auto-merged; the union of its names equals the merged `__all__` — which is exactly
what the test asserts, and it passes.

### `pyproject.toml` (auto-merged) — audited, nothing lost

Both sides bumped `0.20.0 → 0.21.0` with the **same** line, so git had no conflict to
report and nothing to lose; the merged file carries the single `0.21.0` bump.

### `README.md` (auto-merged) — audited, nothing lost

Merged README = main's README **plus** the branch's "Righting arms at finite heel
(0.21)" section (verified by `git diff origin/main -- README.md`, whose only hunk is
that added section). Residual known inconsistency, **pre-existing on main and not
introduced by this merge**: main's own "what this is not" paragraph still lists "the
righting-arm curve" as out of scope while `tensorcad.heel` now provides it. Left
verbatim on purpose (deleting a parent's sentence during a merge is exactly the
silent-loss failure mode this cure is guarding against); recorded here as doc debt.

## Anti-silent-loss evidence

`tools/verify_merge_pr41.py` compares the merged tree against **both** parents, file
by file, and writes `artifacts/m3778/pr41_merge_evidence.json`:

| check | result |
|---|---|
| `base_is_ancestor_of_branch` (`2ff80b2` ⊆ `2051f08`) | PASS |
| `no_conflict_markers` (`git grep -E '^(<<<<<<< \|=======$\|>>>>>>> )'`) | PASS |
| `__all___equals_parent_union` (130 names, no dup) | PASS |
| `all_public_names_importable` | PASS |
| `contract_lists_both_new_modules` (`heel`, `transom`) | PASS |
| branch-side probes (heel.py, test_heel.py, docs/gz_curve.md, example, README, CHANGELOG, `__init__`) | 7/7 PASS |
| main-side probes (transom.py, test_transom.py, test_hydrostatics_inertia.py, docs/transom_stern.md, docs/asymmetric_hull.md, CHANGELOG ×2, `__init__`) | 9/9 PASS |
| `pyproject_version_matches_main_bump` | PASS |
| `nothing_deleted_from_main` (name-status accounting both ways) | PASS |

`git diff --name-status 2ff80b2 2051f08` contains **no `D` rows**: the cure deletes
nothing from `main`; the merged tree is main ∪ the branch's 14 files.

## Tests (local, this sandbox)

The machine is shared and was under load average ≈3500 during this mission, so a
full-suite run could not be completed inside the sandbox's per-command budget
(120 s); two full-suite attempts were started and abandoned, and their logs are
**not** used as evidence (they were also polluted by another mission's pytest
writing the same `/tmp` filenames). What was run and is reported is the **core set
that the merge can actually affect**, compared file-by-file against a clean
checkout of `main` (`2ff80b2`, `git status` clean, 36 test files) with the same
interpreter and the same two exclusions (`tests/test_cli_neural.py`,
`tests/test_neural_pipeline.py` — both fail to import on both trees because the
optional `safetensors` extra is not installed here).

| tree | command | result |
|---|---|---|
| merged `2051f08` | `pytest -q tests/test_public_surface.py tests/test_heel.py tests/test_transom.py` | **132 passed** in 42.77 s |
| `main` `2ff80b2` | `pytest -q tests/test_public_surface.py tests/test_transom.py` | **102 passed** in 25.88 s |
| merged `2051f08` | `pytest -q tests/test_hull.py tests/test_hydrostatics.py tests/test_stability.py` | **71 passed** in 23.40 s |
| `main` `2ff80b2` | same three files | **71 passed** in 14.33 s |

Reading: the 30 branch-only tests (`tests/test_heel.py`) plus the shared 102 pass on
the merged tree (132 = 102 + 30), and the shared core modules are green and
byte-identical in outcome to `main`'s baseline (71 = 71). No test that passes on
`main` fails after the merge in any set that was run.

## Functional evidence from the merged tree

`tools/m3778_pr41_figures.py` re-executes the branch's own physics **from the
merged commit**: righting arms of the barge `L=20 m, B=6 m, T=2 m, depth=4 m,
KG=2 m` (deck edge 33.7°) against the wall-sided closed form
`GZ = sinφ (GM + BM tan²φ/2)`, `KB = 1.0 m`, `BM = B²/12T = 1.5 m`:

* worst relative error over φ = 2…30°: **4.94e-15** (≈ round-off; ≪ 3 %),
* `GZ_max = 0.375000 m` over the sampled range,
* figures: `artifacts/m3778/pr41_gz_after_merge.png`,
  `artifacts/m3778/pr41_content_preservation.png`;
  numbers: `artifacts/m3778/pr41_gz_curve_after_merge.json`.

## Reproduce

```bash
git worktree add --detach /tmp/pr41-m3213 agent/m3213
cd /tmp/pr41-m3213 && git merge origin/main          # 3 conflicts
# resolve as documented above, then:
python3 tools/verify_merge_pr41.py --base 2ff80b2 --branch HEAD \
    --json artifacts/m3778/pr41_merge_evidence.json
python3 tools/m3778_pr41_figures.py --tensorcad-root . \
    --evidence artifacts/m3778/pr41_merge_evidence.json --out artifacts/m3778
python3 -m pytest -q tests/test_public_surface.py tests/test_heel.py tests/test_transom.py
```

## Limitations

* Full-suite green is **not** claimed: only the core set above was executed to
  completion (see the load note). merge_governor's next cycle should re-run the full
  suite on a quiet machine.
* The stale "righting-arm curve is out of scope" sentence on `README.md` was left as
  main wrote it (deliberate, see above).