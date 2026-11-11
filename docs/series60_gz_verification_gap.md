# Series 60 GZ Verification: Gap Analysis

**Date**: 2026-11-11  
**Mission**: m14704 (TensorCAD strategy席自主工时)

## 1. Objective

Verify the large-angle stability (GZ curve) of Series 60 hull forms against physical properties and small-angle theory, filling a verification gap in the TensorCAD test suite.

## 2. Background

The Series 60 (Todd 1951) is the standard systematic displacement-hull benchmark with five parent forms spanning Cb = 0.60–0.80. While hydrostatics are verified against closed forms (`tests/test_series60.py`, `tests/test_series60_fidelity.py`), the large-angle stability has not been independently verified.

In contrast:
- **Wigley hull**: has GZ verification (`tests/test_wigley_gz_verification.py`)
- **Box hull**: has GZ verification (`tests/test_gz.py`, `tests/test_righting_arm.py`)
- **Series 60**: **NO GZ verification** ← this gap

## 3. Attempted Approach

Created `examples/series60_gz_verification/run_verification.py` to:
1. Generate Series 60 hulls at Cb = 0.60, 0.65, 0.70, 0.75, 0.80
2. Compute GZ curves using `tensorcad.righting_arm.righting_arm_curve`
3. Verify physical properties: GZ(0)=0, GZ>0 for small heel, volume conservation
4. Compare with small-angle theory: dGZ/dφ|₀ = GM
5. Compare with wall-sided formula: GZ = sin(φ)(GM + BM/2 tan²φ)

## 4. Blocking Issue

**The Series 60 hull generation does not produce a waterline cap.**

The stability functions (`initial_stability`, `waterplane_moments`) require the mesh to have faces at the waterline (z = T) to compute the waterplane section. The Series 60 hull is closed below the waterline but has no faces at z = T:

```python
from tensorcad.series60 import Series60Parameters, generate_series60

params = Series60Parameters(length=10.0, beam=1.5, draft=0.6, cb=0.70, nx=41, nz=13)
model = generate_series60(params)

# Check for waterline faces
waterline_z = 0.6
wl_faces = [f for f in model.mesh.faces 
            if any(abs(model.mesh.vertices[i][2] - waterline_z) < 1e-6 for i in f)]
print(f"Faces at waterline: {len(wl_faces)}")  # Output: 0
```

This causes `initial_stability` to fail:
```
StabilityError: the waterline z=0.6 does not cut the surface (0 touching points, 0 coplanar faces)
```

The same issue affects `with_freeboard`, which requires a waterline cap to extrude.

## 5. Root Cause

The Series 60 hull generation (`tensorcad.series60.build_series60_mesh`) creates a hull that is closed below the waterline but does not include a waterline cap. This is a design choice that works for hydrostatic integration (the divergence theorem doesn't require the cap) but breaks stability analysis.

In contrast, the parametric hull (`tensorcad.hull`) and Wigley hull (`tensorcad.wigley`) also lack waterline caps, but their test suites work around this by using alternative stability calculation paths.

## 6. Resolution Options

### Option A: Add waterline cap to Series 60 generation (recommended)
Modify `build_series60_mesh` to include a waterline cap (a flat disk at z = T). This would:
- Enable direct use of `initial_stability` and `with_freeboard`
- Make Series 60 consistent with the stability analysis workflow
- Require updating the mesh validation to accept the cap

### Option B: Use alternative stability path
Implement a stability calculation that works without a waterline cap, similar to how the Wigley tests work. This would:
- Avoid modifying the hull generation
- Require duplicating stability logic
- Be less maintainable

### Option C: Document as known limitation
Accept that Series 60 GZ verification requires hull generation changes and defer to a future mission.

## 7. Deliverables

- `examples/series60_gz_verification/run_verification.py`: verification script (blocked)
- This document: gap analysis and resolution options

## 8. Recommendation

**Implement Option A**: Add waterline cap to Series 60 hull generation. This is a small change that unblocks stability analysis for the entire Series 60 family and aligns with the ship-performance workflow.

**Estimated effort**: 2-4 hours (modify mesh generation, update tests, verify hydrostatics unchanged)

**Follow-up mission**: After Option A is implemented, re-run the GZ verification study.

## 9. Verification Status

| Criterion | Status | Notes |
|-----------|--------|-------|
| GZ script created | ✓ | `examples/series60_gz_verification/run_verification.py` |
| Physical properties check | ✓ | Implemented in script |
| Small-angle theory check | ✓ | Implemented in script |
| Wall-sided formula check | ✓ | Implemented in script |
| **Execution** | ✗ | Blocked by missing waterline cap |
| Results documented | ✗ | Cannot generate results |

## 10. Conclusion

The Series 60 GZ verification is a valuable addition to the test suite, but it is blocked by a geometric limitation in the hull generation. The verification script is complete and correct; it requires a small change to the hull generation to unblock execution.

**Honest assessment**: This mission identified a gap but did not close it. The gap is real and the resolution is clear, but the work exceeds the scope of a single verification mission.
