# PR #40 Cure Report: Semantic Union Merge of agent/m3169

## Executive Summary

Successfully resolved merge conflicts between `agent/m3169` (righting-arm curve GZ(phi)) and `origin/main` (transom-stern hull family) for PR #40. Both feature branches contributed independent, additive functionality to version 0.21.0. The merge preserves 100% of content from both sides with zero loss.

**Old SHA (agent/m3169):** `57f5c2f20a31f71d8101146140dd4c0c638b0ac4`  
**New SHA (merged):** `95a5d917e8143593c03572b6078a475acec68c1c`

---

## 1. Conflict Analysis

### 1.1 Conflict Files

Two files had merge conflicts:

1. **CHANGELOG.md** (lines 9-62)
2. **tests/test_public_surface.py** (lines 67-70)

### 1.2 Root Cause

Both branches independently added features to version 0.21.0:
- **agent/m3169**: Righting-arm curve `GZ(phi)` at constant displacement
- **origin/main**: Transom-stern hull family

The conflicts arose because both branches modified the same changelog version section and the same test function, but with non-overlapping content.

---

## 2. Resolution Strategy

### 2.1 CHANGELOG.md

**Conflict structure:**
```
## 0.21.0 — 2026-11-04

<<<<<<< HEAD
[Righting-arm curve content - 29 lines]
=======
### Added
[Transom-stern hull content - 23 lines]
>>>>>>> origin/main
```

**Resolution:** Semantic union - kept both sections under a single `### Added` heading. The righting-arm content was reformatted to match the transom-stern's bullet structure for consistency.

**Verification:** 
- Both feature descriptions are complete
- No technical details lost
- Version number and date preserved
- Formatting is consistent

### 2.2 tests/test_public_surface.py

**Conflict structure:**
```python
for module in ("curvature", "geometry", ...,
<<<<<<< HEAD
               "stability", "righting_arm", "bonjean",
               "wigley", "hull", "design", "suboff", "suboff_appendages", "cli",
=======
               "wigley", "transom", "hull", "design", "suboff", "suboff_appendages", "cli",
>>>>>>> origin/main
               "adapters.fem", "benchmark"):
```

**Resolution:** Union of module lists:
```python
for module in ("curvature", "geometry", ...,
               "stability", "righting_arm", "bonjean",
               "wigley", "transom", "hull", "design", "suboff", "suboff_appendages", "cli",
               "adapters.fem", "benchmark"):
```

**Verification:**
- All modules from both sides present
- No duplicates
- Test passes (6/6)

---

## 3. Test Evidence

### 3.1 Public Surface Test
```
tests/test_public_surface.py::test_declared_api_is_importable PASSED
tests/test_public_surface.py::test_layered_modules_are_all_importable PASSED
tests/test_public_surface.py::test_core_import_pulls_in_no_third_party_geometry_stack PASSED
tests/test_public_surface.py::test_adapters_package_is_a_documented_boundary PASSED
tests/test_public_surface.py::test_neural_runtime_is_opt_in_and_not_imported_by_the_core PASSED
tests/test_public_surface.py::test_hull_module_does_not_import_the_cli_or_vice_versa_in_reverse PASSED

6 passed, 1 warning in 0.54s
```

### 3.2 Righting-Arm Tests (agent/m3169 feature)
```
tests/test_righting_arm.py: 41 passed, 1 warning in 22.56s
```

All righting-arm functionality intact:
- Wall-sided barge closed-form reproduction
- Sphere convergence to GM·sin(φ)
- Parametric hull verification
- Constant displacement enforcement
- Refusal cases (non-watertight, |heel|>90°, etc.)

### 3.3 Transom-Stern Tests (origin/main feature)
```
tests/test_transom.py: 96 passed, 1 warning in 18.91s
```

All transom-stern functionality intact:
- Watertight topology at every τ
- Analytic hydrostatics validation
- Area ratio inversion
- Parameter validation
- Metadata correctness

### 3.4 Core Hydrostatics Tests
```
tests/test_hull.py: PASSED
tests/test_stability.py: PASSED  
tests/test_bonjean.py: 37 passed
```

Combined: 102 passed, 1 warning in 19.91s

### 3.5 Test Summary

| Test Suite | Count | Status |
|------------|-------|--------|
| test_public_surface | 6 | ✅ PASS |
| test_righting_arm | 41 | ✅ PASS |
| test_transom | 96 | ✅ PASS |
| test_hull + test_stability + test_bonjean | 102 | ✅ PASS |
| **Total** | **245** | **✅ ALL PASS** |

---

## 4. Content Preservation Audit

### 4.1 agent/m3169 Contributions (Preserved)

✅ Righting-arm curve `GZ(phi)` implementation  
✅ Wall-sided barge reference (exact to 1e-12 m)  
✅ Floating sphere reference (GZ = GM·sin(φ))  
✅ Divergence theorem integrals (no cap construction)  
✅ 41 tests with refusal cases  
✅ `examples/righting_arm/run_example.py`  
✅ `docs/righting_arm.md` with derivation  

### 4.2 origin/main Contributions (Preserved)

✅ Transom-stern hull family (`tensorcad/transom.py`)  
✅ Truncated taper integrals (sinc^p power series)  
✅ `transom_fraction_for_area_ratio` closed-form inverse  
✅ 96 tests covering topology, hydrostatics, validation  
✅ `examples/transom_stern/run_example.py`  
✅ `docs/transom_stern.md`  
✅ Verification artifacts in `docs/m2873/outputs/transom_stern/`  

### 4.3 Zero-Loss Verification

Both feature sets are fully functional and independently testable. No code, documentation, or test content was lost during the merge.

---

## 5. Technical Details

### 5.1 Merge Base
```
commit df31744cd5659c1fbb31b54607a13f833aa4a552
Merge pull request #38 from cxs503/agent/m2877
```

### 5.2 Branch Divergence
- **agent/m3169**: 1 commit ahead of merge-base (righting-arm)
- **origin/main**: Multiple commits ahead (transom-stern + other PRs)

### 5.3 Merge Commit
```
commit 95a5d917e8143593c03572b6078a475acec68c1c
Merge: 57f5c2f 311f45a
Author: m3774
Date:   2026-11-05

    Merge origin/main into agent/m3169 (cure PR #40, m3774): 
    semantic union of righting-arm and transom-stern, zero content loss
```

---

## 6. Pitfalls & Lessons

### 6.1 Changelog Conflicts Are Semantic, Not Textual

Both branches added to the same version section (0.21.0). A mechanical line-by-line merge would have lost one feature's description. The correct approach is to recognize that changelog entries are additive and preserve both.

### 6.2 Module Lists Must Be Unioned, Not Overwritten

The test file had a list of modules to import. Each branch added different modules. The resolution must include all modules from both sides, not just one.

### 6.3 Independent Features Merge Cleanly When Understood

Both features (righting-arm and transom-stern) are independent additions to the codebase. They don't modify the same functions or data structures. Once the conflicts were understood as "both adding to the same place" rather than "disagreeing about what should be there", the resolution was straightforward.

---

## 7. Conclusion

PR #40 is now clean and ready for merge-governor verification. The semantic union approach preserved 100% of functionality from both branches. All 245 core tests pass, confirming that neither feature was broken by the merge.

**Status:** ✅ READY FOR MERGE-GOVERNOR
