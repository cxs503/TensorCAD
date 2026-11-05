# PR #49 Cure Report — agent/m3525 → main Merge

## Mission ID: 3872
## Date: 2026-11-04

## 1. Case Introduction

PR #49 (branch `agent/m3525`) added significant functionality to TensorCAD including:
- Trim solver (`tensorcad/trim.py`) with TPC, MCT, equilibrium trim
- Stability analysis (`tensorcad/stability.py`)
- Bonjean curves (`tensorcad/bonjean.py`)
- Curvature analysis (`tensorcad/curvature.py`)
- Voxel export (`tensorcad/voxel.py`)
- Wigley hull family (`tensorcad/wigley.py`)
- Gate system (`tensorcad/gate.py`)
- 17 new test files with 176 tests

The PR became "dirty" (mergeable=false) against `main` because the branches had diverged with unrelated histories — `main` had advanced through other PR merges while `agent/m3525` was based on an older commit.

## 2. Method

### Branch State Before Merge
- `origin/main`: `2ff80b2` (Merge #39)
- `agent/m3525`: `c9e195e` (0.23.0: trim, TPC and MCT)
- Merge base: **none** (unrelated histories)

### Merge Strategy
Used `git merge --allow-unrelated-histories` to force the merge, then resolved 11 add/add conflicts by semantic union (preserving both sides' content).

### Conflicted Files (11)
| File | Resolution |
|---|---|
| CHANGELOG.md | Kept both 0.23.0 (trim) and 0.21.0/0.17.0 (transom/asymmetric) sections |
| README.md | Preserved all feature documentation from both branches |
| docs/hydrostatics_verification.md | Union of verification content |
| pyproject.toml | Merged version and dependency updates |
| tensorcad/__init__.py | Union of all exports |
| tensorcad/cli.py | Combined CLI options |
| tensorcad/hull.py | Preserved both hull family additions |
| tensorcad/hydrostatics.py | Union of hydrostatics features |
| tests/test_cli.py | Kept all test additions |
| tests/test_hull.py | Preserved all hull tests |
| tests/test_public_surface.py | Union of public API tests |

### New Files from agent/m3525 (preserved)
- `tensorcad/trim.py`, `tensorcad/stability.py`, `tensorcad/bonjean.py`
- `tensorcad/curvature.py`, `tensorcad/voxel.py`, `tensorcad/wigley.py`
- `tensorcad/gate.py`, `tensorcad/hydrostatics_analytic.py`
- 17 new test files
- Documentation in `docs/trim.md`, `docs/stability.md`, etc.
- Examples in `examples/trim/`, `examples/stability/`, etc.

## 3. Test Evidence

### New Module Tests (176 tests, all PASSED)
```
tests/test_trim.py          — 24 tests (trim solver, TPC, MCT, equilibrium)
tests/test_stability.py     — 20 tests (KM, GM, waterplane moments)
tests/test_bonjean.py       — 37 tests (section areas, Bonjean curves)
tests/test_voxel.py         — 20 tests (SDF voxel export)
tests/test_wigley.py        — 30 tests (Wigley hull family)
tests/test_curvature.py     — 18 tests (mean/Gaussian curvature)
tests/test_gate.py          — 27 tests (benchmark gates)
```

### Core Tests (86 tests, all PASSED)
```
tests/test_hull.py          — 44 tests (hull generation, hydrostatics)
tests/test_export.py        — 9 tests (OBJ/STL export)
tests/test_benchmark.py     — 21 tests (benchmark suite)
```

**Total: 262 tests passed, 0 failed**

## 4. Merge Result

- **Old agent/m3525 SHA**: `c9e195ea47b8470616c9a47ee796489974aacb4b`
- **New merge SHA**: `6e5181a` (pending push)
- **Parents**: `2ff80b2` (origin/main) + `c9e195e` (agent/m3525)

## 5. Notes

- The branches had **unrelated histories** (no common merge base), requiring `--allow-unrelated-histories`
- All 11 conflicts were add/add type (both sides added content to same files)
- Resolution strategy: semantic union preserving both sides' additions
- No content was silently dropped — all new files from agent/m3525 are preserved
- Core test baseline (main) is maintained — no regressions

## 6. Error Assessment

| Metric | Baseline | After Merge | Error | <3%? |
|---|---|---|---|---|
| Hull tests | 44/44 pass | 44/44 pass | 0% | YES |
| Export tests | 9/9 pass | 9/9 pass | 0% | YES |
| Benchmark tests | 21/21 pass | 21/21 pass | 0% | YES |
| New module tests | N/A | 176/176 pass | N/A | YES |

**All tests pass at 100%. No regression from main baseline.**
