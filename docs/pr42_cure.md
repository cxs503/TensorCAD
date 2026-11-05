# PR #42 cure — `agent/m3300` (large-angle righting-arm curve) re-based onto current `main`

Mission m3779. The pull request `agent/m3300` → `main` was reported `mergeable=false`,
because the branch was cut from `main` at `ca0dde3` (merge base `df31744`), while `main`
advanced 12 commits (merge-governor merges #38 … #39 and the CI-debt clearance
`dc9b66c`: `ruff 50 → 0`, `mypy 7 → 0`).

This note records the three-way merge, the per-file resolution, and the evidence that
neither side lost content.

## 1. Topology

| ref | commit | note |
|---|---|---|
| merge base (`git merge-base main agent/m3300` on the true main) | `df31744` | `Merge pull request #38 from cxs503/agent/m2877` |
| `main` (as published, `origin/main`) | `2ff80b2` | +12 commits over the base |
| `agent/m3300` (old tip) | `cc8985b` | +1 commit over the base |
| local `main` branch in this worktree | `ca0dde3` | stale (41 commits behind `origin/main`); **not** used as the merge target |

Note: the local `main` ref predates 41 upstream commits, so a naive
`git merge agent/m3300 main` would have "succeeded" against the wrong parent and
re-introduced 41 reverted commits. The merge was therefore done against **`origin/main`
(`2ff80b2`)**, which is what the governor compares against.

## 2. What each side changed

`agent/m3300` (`cc8985b`, 13 files, +3818/−4): `tensorcad/gz.py`, `tests/test_gz.py`,
`docs/gz.md`, `docs/gz/{report.md,results.json,3 × .png}`,
`examples/gz/run_example.py`, plus edits to `CHANGELOG.md`, `pyproject.toml`,
`tensorcad/__init__.py`, `tests/test_public_surface.py`.

`main` (12 commits, 38 files): transom-stern family (`tensorcad/transom.py`,
`docs/transom_stern.md`, `tests/test_transom.py`), waterplane second moments
(`tensorcad/hydrostatics.py`, `tensorcad/hydrostatics_analytic.py`,
`tests/test_hydrostatics_inertia.py`), asymmetric hull, CI debt clearance, plus edits to
`CHANGELOG.md`, `pyproject.toml`, `tensorcad/__init__.py`, `tests/test_public_surface.py`.

File-set overlap (the only places a conflict is possible):

```
CHANGELOG.md  pyproject.toml  tensorcad/__init__.py  tests/test_public_surface.py
```

`pyproject.toml` auto-merged: both sides bumped `version` `0.20.0 → 0.21.0` to the same
value, so the merge is a no-op there.

## 3. Conflicts and how each was resolved (semantic union, no side dropped)

`git merge --no-ff --no-commit agent/m3300` on `origin/main` reported 3 content conflicts
(`CHANGELOG.md`, `tensorcad/__init__.py`, `tests/test_public_surface.py`). All three are
*additive* conflicts — both sides append to a shared list — so the resolution is a union,
never a pick-a-side.

### 3.1 `CHANGELOG.md`
Both sides opened a `## 0.21.0 — 2026-11-04` section: `main` with a `### Added` heading
carrying the waterplane-inertia / transom bullets, `agent/m3300` with the four GZ bullets
written without a heading (its file had no `### Added` sub-header).

Resolution: one `### Added` section — `main`'s heading and four bullets first, then
`agent/m3300`'s four GZ bullets appended verbatim under the same heading. Both bullet sets
are kept in full; the only added characters are the shared heading (which already existed
on the main side). Ordering inside the section is the merge-base order of each side.

### 3.2 `tensorcad/__init__.py` — `__all__`
```
HEAD:    "sinc_power_coefficients",
m3300:   "solve_waterline", "station_offsets", "submerged_properties",
```
Resolution: keep all four (`__all__` is alphabetically sorted here, and the union happens to
be in sorted order):
`sinc_power_coefficients`, `solve_waterline`, `station_offsets`, `submerged_properties`.
The `from .gz import (...)` / `from .transom import (...)` blocks landed in *different*
places in the file and were auto-merged by git without conflict; they are both present
(`__init__.py:70` and `:141`). Verified: `len(tensorcad.__all__) == 133` and
`[n for n in tensorcad.__all__ if not hasattr(tensorcad, n)] == []`.

### 3.3 `tests/test_public_surface.py` — the layered-module import list
```
HEAD:  "wigley", "transom", "hull", "design", "suboff", "suboff_appendages", "cli",
m3300: "wigley", "hull", "design", "suboff", "suboff_appendages", "gz", "cli",
```
Resolution: union into a single line, keeping both new modules and each side's relative
order:
`"wigley", "transom", "hull", "design", "suboff", "suboff_appendages", "gz", "cli",`
(`transom` from main, `gz` from the PR). This is the test that would have caught a silent
drop of either module.

Also preserved from main, on a line adjacent to the conflict:
`"wigley_reference_parameters",    # transom-stern family (0.21.0)` — main's trailing
comment survives.

## 4. Silent-loss audit (the failure mode this cure exists to avoid)

A mechanical overlay (`git checkout --theirs` / `-X ours`) would have dropped one side's
bullets and exports with a green test run. Instead the merged tree was diffed against
**both** parents, file by file:

* for every file either side touched, `git diff <parent> -- <file>` must contain **no
  deletion** of the other side's lines;
* result: only the three conflict files differ from either parent (they are supersets), and
  there are **no missing files**: all 10 new files from `agent/m3300` and all 38 files from
  `main` are present in the merge.

## 5. Gate parity — the merge must not lower `main`'s bar

`agent/m3300` was cut before `main`'s CI-debt clearance (`dc9b66c`, `mypy 7 → 0`), so the
merge, taken verbatim, regressed `mypy` by 5 errors — all in the PR's own `tensorcad/gz.py`.
Because "not lower than the main baseline" applies to the type gate as much as to pytest,
the 5 were fixed **typing-only** (2 lines, behaviour-neutral, zero numeric change):

* `_submerged(mesh, unit: tuple[float, float, float], …)` → `unit: Sequence[float]`
  (`_submerged` only indexes `unit[0..2]`, exactly like `_plane_point`/`_clip_triangle`,
  which already take `Sequence[float]`) — clears the three `arg-type` errors at the
  `unit = tuple(component / length …)` call sites;
* `polygon` in `_submerged` declared `Sequence[tuple[float, float, float]]` before the face
  loop, because the variable holds either a 3-tuple (fully wet face) or a clipped *list*;
  the declared type also makes the `len(polygon) < 3` guard reachable again (clears the
  `assignment` + `unreachable` pair).

## 6. Evidence

Full suite, `python3 -m pytest tests -q` with
`--ignore=tests/test_cli_neural.py --ignore=tests/test_neural_pipeline.py`
(both modules fail to *collect* in this image: `ModuleNotFoundError: No module named
'safetensors'` — pre-existing, identical on `origin/main`):

| tree | result | wall |
|---|---|---|
| `origin/main` `2ff80b2` (pristine `git archive`) | **2 failed, 818 passed, 11 skipped** | 646 s |
| merge result (`agent/m3779`) | **2 failed, 863 passed, 11 skipped** | 668 s |

`863 − 818 = 45` = exactly the 45 cases of `tests/test_gz.py`, the PR's own new tests: the
merge adds coverage and removes none. The 2 failures are the same two tests on both trees
(`tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol`,
`…::test_the_learned_field_gradient_is_autograd_and_points_outward`), both caused by the
missing `safetensors` wheel — an environment gap, not a regression.

Per-file sweep (`runs/cure_pr42/per_file_merged.log`, 35 modules run one at a time):
every module `rc=0` except `test_geometry_field.py` (`2 failed / 47 passed`, the same
safetensors gap). Notable: `test_gz.py` 45 passed, `test_transom.py` 96 passed,
`test_benchmark.py` 32 passed + 1 skipped, `test_public_surface.py` 6 passed.

Static gates on the merged tree: `ruff check tensorcad tests` → *All checks passed*;
`mypy` → *Success: no issues found in 33 source files* (parity with main's 0/0).

## 7. Push

`agent/m3300`: `cc8985b` → see `runs/cure_pr42/push_record.txt` for the new SHA and the
exact `git push` outcome. The merge commit has first parent `origin/main` `2ff80b2` and
second parent `agent/m3300` `cc8985b`, so the branch tip is a fast-forward of the old tip
and the merge-governor can re-evaluate `mergeable` without a force-push.

## 8. Reproduce

```sh
git worktree add wt origin/main && cd wt
git merge --no-ff --no-commit agent/m3300          # 3 content conflicts
# resolve per §3 (union), then:
python3 -m pytest tests -q --ignore=tests/test_cli_neural.py \
        --ignore=tests/test_neural_pipeline.py     # 2 failed, 863 passed, 11 skipped
python3 -m ruff check tensorcad tests && python3 -m mypy
```