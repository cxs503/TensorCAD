# PR #47 cure (m3917): 3-way merge agent/m3431 -> main, no silent loss

## 1. Situation

TensorCAD PR #47 (`agent/m3431`, mid 3431) delivered 0.21.0 **large-angle righting arms**
(`tensorcad/righting.py`, GZ curves solved by Newton on the plane offset) on top of merge
base `e604f02`. Meanwhile main advanced 10 commits and itself published 0.21.0 with the
**waterplane second moments** and the **transom-stern hull family**. Both sides therefore
bumped `pyproject.toml` to 0.21.0 and both appended to the same four shared files, so
GitHub reported `mergeable=false` / `dirty`.

`main` API surface 120 names, `agent/m3431` 122 names,
union 135.

## 2. Conflict set and per-file resolution

| file | conflict | resolution |
|---|---|---|
| `tensorcad/__init__.py` | 2 hunks, both in `__all__` (righting names vs transom/sinc names) | union of both name sets in approximate alphabetical order; both import blocks kept (`from .righting import (`, `from .transom import (`) |
| `CHANGELOG.md` | one `## 0.21.0` section per side | one section holds both: main's `### Added` block followed by m3431's righting bullets |
| `README.md` | same "what it does **not** provide" list rewritten by both | every item of both lists kept; the two statements that the merged tree falsifies (no transom family; no GZ / righting-arm curve) re-stated as three families + GZ exists while the *criteria* do not |
| `tests/test_public_surface.py` | auto-merged names, but the layered-module tuple only had `transom` | tuple extended to `"transom", "righting"` (superset, nothing dropped) |
| `pyproject.toml` | none (both wrote 0.21.0) | kept 0.21.0 -- one release, both feature sets |

## 3. No-silent-loss audit (`tools/verify_merge_m3917.py`)

Line-set audit over all 42 files touched by either side: every line *added*
since the merge base by either parent must survive in the merged blob.

* main lines lost: **0**; m3431 lines lost: **0**
* extended (superset) lines: 1 (`test_public_surface.py` module tuple)
* re-stated lines (README, contradictory factual claim): 3
* `__all__` merged = 135 = |main ∪ branch| = 135 (no name missing on either side)
* `EXPECTED_API == __all__`: True

Machine artifact: `artifacts/m3917/merge_audit.json` (VERDICT True).

## 4. Numerical re-verification of the merged tree (not just "tests are green")

Box 6x3x4 m, T=2.0 m, KG=1.0 m
(GM = KB + BM - KG = 0.375 m). `righting_arm` (3-D tetrahedron clip on the merged
tree) vs `analytic_box_righting_arm` (independent 2-D section clip):

| heel (deg) | GZ merged tree (m) | GZ analytic (m) | abs err (m) |
|---|---|---|---|
| 0 | 0.000000000 | 0.000000000 | 0.00e+00 |
| 5 | 0.032808487 | 0.032808487 | 8.33e-17 |
| 10 | 0.066130366 | 0.066130366 | 1.39e-16 |
| 15 | 0.100541337 | 0.100541337 | 1.67e-16 |
| 20 | 0.136752971 | 0.136752971 | 2.22e-16 |
| 25 | 0.175712219 | 0.175712219 | 1.67e-16 |
| 30 | 0.218750000 | 0.218750000 | 3.05e-16 |
| 35 | 0.267819751 | 0.267819751 | 4.44e-16 |
| 40 | 0.325903947 | 0.325903947 | 3.89e-16 |

max abs err **4.44e-16 m**, max rel err **2.54e-15** -> far below the
3 % acceptance bar (it is at round-off, i.e. the merge changed no arithmetic).

## 5. Test evidence

| run | tree | tests |
|---|---|---|
| `pytest tests/test_public_surface.py test_transom.py test_stability.py test_hydrostatics_inertia.py` | merged `agent/m3431` | 137 passed |
| `pytest tests/test_hull.py test_wigley.py test_righting.py` | merged `agent/m3431` | 114 passed |
| `pytest` 10-file core group (baseline) | `origin/main` 2ff80b2 | 310 passed |

Zero failures on the merged tree; the merged runs cover the same files as the 310-test
baseline group plus `test_righting.py` (40 tests, new in this PR). `test_cli_neural.py` and
`test_neural_pipeline.py` are excluded on **both** trees: `safetensors` (the `neural` extra)
is not installed in this image -- a pre-existing environment limit, identical on main.

## 6. Reproduce

```
git merge origin/main            # 3 conflicts
python3 tools/resolve_merge_m3917.py    # semantic-union resolution
python3 tools/verify_merge_m3917.py c699a0d   # no-silent-loss audit -> PASS
python3 tools/report_merge_m3917.py     # GZ re-verification + docs/m3917/results.json + figure
```

## 7. Limitations

* `pyproject.toml` stays at 0.21.0: both parents already claimed that version for their own
  feature set, so the merge does not invent a third number. A reviewer may prefer 0.22.0.
* The audit is line-set based: it proves no line disappeared, not that prose still reads
  agreeably. The three re-stated README lines are printed by the audit for review.
* Only the righting feature was numerically re-verified here; transom/hydrostatics inertia
  keep their own PR-level evidence (`tests/test_transom.py`, `tests/test_hydrostatics_inertia.py`).
