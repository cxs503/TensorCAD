# PR #40 cure — three-way merge of `agent/m3169` into `main` (mission m3775)

**Scope.** PR #40 (branch `agent/m3169`, mission m3169: *righting-arm curve `GZ(φ)` at
constant displacement*) was reported `mergeable = false` against `main`. This document
records the re-base of that branch onto current `main`, the per-file conflict resolution,
and the machine-checked evidence that nothing from either side was dropped.

Machine-readable twin of this page: `artifacts/cure/pr40/evidence.json`
(regenerate with `python3 artifacts/cure/pr40/audit_merge.py`).

## 1. Topology

| item | value |
|---|---|
| repository | `cxs503/TensorCAD` |
| branch under cure | `agent/m3169` (old SHA `57f5c2f`) |
| target main at merge time | `origin/main` = `2ff80b2` (*Merge #39 (m2925) via merge-governor*) |
| true merge base | `df31744` (*Merge pull request #38 from cxs503/agent/m2877*) |
| merge commit | `5c80e49` (parents `2ff80b2`, `57f5c2f`) |
| delivered on branch | `agent/m3775` (the mission branch the platform pushes) |

Note on the "dirty" origin: the branch was **not** ancestor-disjoint. Local `main` was
stale at `ca0dde3`, which made every distance measurement look catastrophic (95 files
"touched by both sides"). Re-measured against `origin/main` the picture is small and
exact: **main advanced 12 commits, the branch added 1 commit**, four files were touched by
both sides, and nine files are new on the branch alone. The cure is therefore a real
three-way merge, not a rebuild.

Files new on `agent/m3169` (carried over verbatim, 9):

```
tensorcad/righting_arm.py          tests/test_righting_arm.py
examples/righting_arm/run_example.py
docs/righting_arm.md               docs/righting_arm/{report.md,results.json,
docs/righting_arm/barge_vs_closed_form.png,
docs/righting_arm/sphere_vs_analytic.png, hull_curve.png}
```

34 files changed on `main` only since the merge base — untouched by this cure, they are
what the branch had to be re-based onto.

## 2. Conflicting files and how each was resolved

Semantic union was applied; no side's addition or deletion was discarded. Mechanical
`--theirs`/`--ours` overlay was explicitly not used.

### 2.1 `CHANGELOG.md`

Both sides opened the same unreleased section `## 0.21.0 — 2026-11-04` — `main` with
`### Added` (waterplane second moments on the `Hydrostatics` record + transom-stern hull
family and its hydrostatic references), the branch with nothing but righting-arm bullets
under the heading. Git could not interleave them because one side inserted an `### Added`
header where the other inserted a paragraph.

**Resolution:** one `### Added` heading, `main`'s bullets first, then the branch's
righting-arm bullets, separated by a blank line. Version stays `0.21.0` because *both*
sides bumped `0.20.0 → 0.21.0` independently (verified: same deletion set on both sides),
so a further bump would be a version gratuitously invented by the merge.

### 2.2 `tests/test_public_surface.py`

Two conflicts in one file:

* `test_declared_api_is_importable` — the version assertion, both sides `"0.21.0"`
  (auto-merged, no action).
* `test_layered_modules_are_all_importable` — the module tuple. `main` had
  `"wigley", "transom", "hull", "design", "suboff", "suboff_appendages", "cli"`; the
  branch had `"stability", "righting_arm", "bonjean"` plus an older
  `"wigley", "hull", ...` line that predates `transom`.

**Resolution:** order-preserving union of the two tuples —
`"wigley", "transom", "hull", "design", "suboff", "suboff_appendages", "cli", "stability",
"righting_arm", "bonjean"` — with `transom` kept (it is main's new module and is
importable) and the branch's three extra modules appended. Reflowed across two
continuation lines at the file's 100-column style.

### 2.3 `tensorcad/__init__.py` (auto-merged, audited)

Both sides appended an export block and both extended `__all__`. Textual auto-merge
produced a **union of 130 names** = 120 (`main`) ∪ 117 (`m3169`); the audit below proves
no name from either side is missing and the righting-arm names (`RightingArm`,
`RightingArmCurve`, `RightingArmError`, `righting_arm`, `righting_arm_curve`,
`wall_sided_righting_arm`) and the transom names both survive. `__version__ = "0.21.0"`
once (both sides bumped identically).

### 2.4 `pyproject.toml` (auto-merged, audited)

`version = "0.21.0"` only; both sides made the identical edit, result is byte-identical
to both sides.

## 3. Zero-content-loss audit (not a claim — a check)

`artifacts/cure/pr40/audit_merge.py` recomputes, from `git` and the working tree:

1. the both-touched / only-theirs / only-main file census from the true merge base;
2. for each both-touched file, the set of **lines** each side added relative to the base,
   and whether each of those lines is still present in the merged file. A line counts as
   preserved if it is present verbatim, or if all of its tokens appear on some line of the
   merged file (reflow). Result: **0 missing lines for `main`, 0 missing lines for
   `agent/m3169`** across all four files;
3. the deletion side: both sides deleted only the four `0.20.0` version strings, so no
   removal had to be arbitrated;
4. conflict markers in the tree: **0**.

This matters because a mechanical overlay on these files silently drops one side's
`__all__` entries — the failure mode this project has hit before.

## 4. Test evidence

Environment note: `tests/test_cli_neural.py` and `tests/test_neural_pipeline.py` **fail to
collect on both trees** with `ModuleNotFoundError: No module named 'safetensors'` — an
image-level gap, present before this merge and unchanged by it, so both runs exclude the
same two modules.

| | origin/main `2ff80b2` (baseline) | merged `5c80e49` |
|---|---|---|
| tests | 831 | 872 (+41) |
| passed | 818 | 859 |
| failed | 2 | 2 |
| skipped | 11 | 11 |
| regressions (pass → not pass) | — | **0** |

The 41 new tests are exactly `tests/test_righting_arm.py`; the two failures are the same
test ids on both trees
(`tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol`,
`tests/test_geometry_field.py::test_the_learned_field_gradient_is_autograd_and_points_outward`),
i.e. the torch-backed learned-field contract that the image does not satisfy. The failure
sets are compared by test id, not by count: `failures_identical = true`.

Commands (per-tree full suite split into four file groups because the shell in this
container caps a single call at 120 s; JUnit XML is the record):

```bash
python3 -m pytest -q -n 16 --timeout=900 -p no:randomly -p no:cacheprovider \
    --junitxml=/tmp/base<N>.xml <group files>          # baseline, /tmp/base3775
python3 -m pytest -q -n 16 --timeout=900 -p no:randomly -p no:cacheprovider \
    --junitxml=/tmp/merge<N>.xml <group files>         # merged, this worktree
python3 artifacts/cure/pr40/audit_merge.py             # diffs the two outcome sets
```

A single-process full run also completed on the merged tree before the container got
busy: `2 failed, 859 passed, 11 skipped in 68.35s` (`artifacts/cure/pr40/pytest_merged.log`).

### Feature still live, not merely present

`PYTHONPATH=. python3 examples/righting_arm/run_example.py` on the merged tree prints
`10/10 claims hold; wrote docs/righting_arm` and leaves **`git status` clean** — the
merged code reproduces the artefacts the branch committed, byte for byte. The barge
`GZ` matches the wall-sided closed form to `1.21e-14 m` up to the 21.80° validity angle,
the hull peak is `0.1467 m` at 50°, the sphere error at resolution 12 is `5.28e-4 m`.
Figure: `artifacts/cure/pr40/pr40_merge_verification.png`
(`python3 artifacts/cure/pr40/make_figure.py`).

## 5. Pushing back to `agent/m3169` — what was and was not done

* The merge commit `5c80e49` carries both parents, so `agent/m3169`'s tip is contained in
  it; merging `agent/m3775` (this branch, `5c80e49`) into `main` lands PR #40's content
  with no further conflict, and PR #40 itself then fast-forwards or closes as merged.
* **Not done here:** `git push` (the harness forbids it; the platform pushes the mission
  branch) and a local `git branch -f agent/m3169 5c80e49` — that branch is checked out in
  another mission's worktree (`/root/tp-agent/work/m3169`), and force-moving a ref under a
  live worktree would leave that agent with a stale index. Recorded mapping:
  `agent/m3169` old `57f5c2f` → new content `5c80e49`.

## 6. Limits / honest notes

* The merge is a *text*-level union on four files; the semantic guarantee for those files
  comes from (a) the export/layer tests, which are exactly the contract that would break
  and (b) the line-level audit. Other semantic interactions are covered only by the suite.
* The two `test_geometry_field` failures and the two `safetensors` collection errors are
  pre-existing environment debt of this image, not of PR #40; they are not fixed here.
* The baseline was measured after the merged run; the container was under different load,
  so wall-clock times are not comparable — only outcome sets are.