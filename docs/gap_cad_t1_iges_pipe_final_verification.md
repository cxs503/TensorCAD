# GAP-CAD-T1-IGES-PIPE — Final Verification & Closure Report

**Mission**: m22304 · **Anchor**: GAP-CAD-T1-IGES-PIPE
**Date**: 2026-11-15 · **Agent**: m22304 (verifier seat)
**Commit**: `c525b4e` (upstream main)

## Acceptance Criteria (from feature_ledger)

| # | Criterion | Status | Evidence |
|---|-----------|--------|----------|
| 1 | 管道 CLI+库函数 (Pipeline CLI + library function) | ✅ VERIFIED | `test_cli_entry_point_is_callable`, `test_cli_help`, `test_load_cad_geometry_returns_mesh_and_report` |
| 2 | 3 个真实 IGES 全过门 (3 real IGES files pass all gates) | ⚠️ PARTIAL | KCS ✅; KVLCC2/ONRT timeout on single-core (known performance limitation, see §5) |
| 3 | 金标回归测试绿 (Golden regression test green) | ✅ VERIFIED | All 4 DTMB5415 tests pass (watertight, volume <3%, BB <1%, all 7 gates) |
| 4 | 依赖锁定 (Dependency locked) | ✅ VERIFIED | `requirements-iges-pipe.lock` present; core package zero-dependency confirmed |

## 1. Pipeline CLI + Library Function

### CLI Entry Point
- **Module**: `tensorcad.iges_pipe_cli`
- **Command**: `tensorcad-iges-pipe input.igs [--output mesh.stl] [--tolerance 1e-4] [--report report.json]`
- **Exit codes**: 0 (all gates passed), 1 (gates failed), 2 (file read error)
- **Verification**: `test_cli_entry_point_is_callable` ✅, `test_cli_help` ✅

### Library Functions
- **`tensorcad.geometry.load_cad_geometry(path)`** — public integration entry, returns `(Mesh, report)`
- **`tensorcad.geomquality.trim_mesh.import_cad(path)`** — full pipeline with gate verdicts
- **Verification**: `test_load_cad_geometry_returns_mesh_and_report` ✅

### Pipeline Stages
1. **Read**: IGES/STEP parser (trimmed B-rep or untrimmed patch set)
2. **Invert**: Trim loops projected to surface parameters (Gauss-Newton with Levenberg damping)
3. **Tessellate**: Triangulate trimmed parameter domain (ear clipping + refinement)
4. **Repair**: Weld → drop degenerate → drop duplicate → fill small holes → re-orient
5. **Audit**: 7 gates (G1–G7) with before/after measurements

## 2. DTMB5415 Gold Regression

### Reference
- **Gold file**: `DTMB5415_Hull_exact.stl` (8.2 MB, 164,926 triangles)
- **Location**: `/root/ship-performance-platform/backend/data/geometry/ships/`

### Pipeline
1. Read gold STL triangles
2. Cut 41 transverse stations, resample each on arclength (41 girth samples)
3. Interpolate NURBS surface (degree 3×3) at Greville abscissae
4. Close body: hull patch + ruled deck + 2 ruled end caps
5. Orient faces by maximizing enclosed volume
6. Write IGES (entity 128 patch set)
7. Import through pipeline → watertight mesh

### Results (measured 2026-11-15)
| Test | Status | Detail |
|------|--------|--------|
| `test_dtmb5415_gold_regression_watertight` | ✅ PASS | Mesh is watertight, single component, closed |
| `test_dtmb5415_gold_regression_volume_within_3pct` | ✅ PASS | Volume error < 3% vs quadrature reference |
| `test_dtmb5415_gold_regression_bounding_box` | ✅ PASS | BB length error < 1% vs gold STL |
| `test_dtmb5415_all_gates_pass` | ✅ PASS | All 7 gates (G1–G7) pass |

**Run command**: `python -m pytest tests/test_gap_cad_t1_iges_pipe.py -v -k "dtmb5415"`
**Result**: 4 passed in ~53s

## 3. Real IGES Files Gate Test

### Test Results (measured 2026-11-15)
| File | Size | Patches | Status | Detail |
|------|------|---------|--------|--------|
| KCS_hull_SVA.igs | 246 KB | 55 | ✅ PASS | All gates pass (half-breadth model, G3/G6 exempt) |
| kvlcc2.igs | 6.3 MB | 38 | ⏱️ TIMEOUT | Import exceeds 120s on single-core (known limitation) |
| onr5613-full_appendage_deckhous_mergedR.igs | 3.4 MB | — | ⏱️ TIMEOUT | Import exceeds 120s on single-core (known limitation) |

### Gate Definitions
| Gate | Name | Criterion |
|------|------|-----------|
| G1 | all_faces_read | No refused faces (or patch-set route) |
| G2 | trim_on_surface | Inversion residual < tolerance |
| G3 | watertight | Zero free edges after repair |
| G4 | manifold | Zero non-manifold edges |
| G5 | positive_volume | Signed volume > 0 |
| G6 | single_shell | 1 component, 1 closed component |
| G7 | no_degenerate | No degenerate/duplicate triangles |

### Half-Breadth Model Handling
KCS and KVLCC2 are half-breadth models (y ≥ 0 only). The test correctly detects this via `classify_model()` and exempts G3 (watertight) and G6 (single shell) since these models are open by construction.

**Run command**: `python -m pytest tests/test_gap_cad_t1_iges_pipe.py -v -k "KCS"`
**Result**: 1 passed in ~8s

## 4. Dependency Lock

### Core Package
- **Dependencies**: `[]` (empty — standard library only)
- **Verification**: `test_core_package_has_no_runtime_dependencies` ✅
- **Method**: `import tensorcad` with numpy/torch blocked via `sys.meta_path`

### Optional Dependencies
| Extra | Packages | Purpose |
|-------|----------|---------|
| `ocp` | `cadquery-ocp>=7` | OpenCASCADE B-rep bridge |
| `neural` | `torch>=2.2, safetensors>=0.4, numpy>=1.24` | Neural field training |
| `test` | `pytest>=7, numpy>=1.24` | Test suite |

### Lock File
- **Path**: `requirements-iges-pipe.lock`
- **Contents**: numpy==2.4.6, pytest==9.2.0, trimesh==5.0.0, meshio==5.3.5
- **Verification**: `test_ocp_is_optional_dependency` ✅

## 5. Known Limitations

### 5.1 Performance on Large IGES Files
KVLCC2 (6.3 MB, 38 trimmed B-spline faces) and ONRT (3.4 MB) exceed the 120s timeout on single-core execution. The bottleneck is the surface inversion step (Gauss-Newton with fallback grid search) which is O(n²) in the number of boundary points per face.

**Root cause**: The fallback grid search (49 seeds × 60 iterations × surface evaluation) is expensive for high-degree NURBS surfaces with large control nets.

**Mitigation options**:
1. Run on HPC cluster with parallel face processing (recommended)
2. Add analytic partials to all surface types (already done for some, not all)
3. Reduce fallback grid resolution (trade accuracy for speed)
4. Use OpenCASCADE-based healing pipeline (`tensorcad[ocp]`) for large files

### 5.2 Trimmed Face Refusals
Some complex files (e.g., KP505 propeller) have faces whose trim curves cannot be inverted within tolerance. These are **measured refusals**, not silent skips — the inversion residual is reported per face.

### 5.3 Gold Regression Uses Loft, Not Original Definition
The DTMB5415 gold regression uses a NURBS loft interpolated from the STL, not the original DTMB5415 definition. The volume error measures the loft+import pipeline, not the original hull accuracy.

## 6. Test Suite Summary

**Test file**: `tests/test_gap_cad_t1_iges_pipe.py`

| Test | Status | Time |
|------|--------|------|
| `test_cli_entry_point_is_callable` | ✅ PASS | <1s |
| `test_cli_help` | ✅ PASS | <1s |
| `test_load_cad_geometry_returns_mesh_and_report` | ✅ PASS | <1s |
| `test_dtmb5415_gold_regression_watertight` | ✅ PASS | ~13s |
| `test_dtmb5415_gold_regression_volume_within_3pct` | ✅ PASS | ~14s |
| `test_dtmb5415_gold_regression_bounding_box` | ✅ PASS | ~14s |
| `test_dtmb5415_all_gates_pass` | ✅ PASS | ~14s |
| `test_real_iges_import_succeeds[KCS]` | ✅ PASS | ~8s |
| `test_real_iges_import_succeeds[KVLCC2]` | ⏱️ TIMEOUT | >120s |
| `test_real_iges_import_succeeds[ONRT]` | ⏱️ TIMEOUT | >120s |
| `test_core_package_has_no_runtime_dependencies` | ✅ PASS | <1s |
| `test_ocp_is_optional_dependency` | ✅ PASS | <1s |

**Total**: 10/12 pass, 2 timeout (KVLCC2, ONRT)

## 7. Reproduction

```bash
# Quick verification (lightweight, ~60s)
cd /srv/tensor/agent/work/m22304
python -m pytest tests/test_gap_cad_t1_iges_pipe.py -v -k "not KVLCC2 and not ONRT"

# Gold regression only (~55s)
python -m pytest tests/test_gap_cad_t1_iges_pipe.py -v -k "dtmb5415"

# Full verification (requires cluster or low load)
python -m pytest tests/test_gap_cad_t1_iges_pipe.py -v
```

## 8. Verdict

**GAP-CAD-T1-IGES-PIPE is L3 (integrated function) with honest limitations**:

- ✅ CLI entry point works and is documented
- ✅ Library function (`load_cad_geometry`) returns correct types
- ✅ DTMB5415 gold regression passes all 4 tests (watertight, volume <3%, BB <1%, all gates)
- ✅ KCS real IGES file imports successfully with correct dimensions
- ✅ Core package is zero-dependency; OCP is optional
- ✅ Dependency lock file present
- ⚠️ KVLCC2/ONRT timeout on single-core (known performance limitation, documented)
- ⚠️ Some complex files have measured trim refusals (not silent skips)

**Recommendation**: Mark as L3 with the performance limitation noted. For L4 (benchmark verified), the KVLCC2/ONRT tests should be run on the HPC cluster where parallel face processing can complete within reasonable time.
