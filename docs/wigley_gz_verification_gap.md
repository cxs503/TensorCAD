# Wigley GZ verification gap

> **STATUS: CLOSED** (2026-11-09, TensorCAD). Both open items in §5 are fixed and verified:
> the waterplane cut now accepts the freeboard seam (`_plane_segments` traces a plane-lying
> edge once) and `initial_stability` measures the immersed body instead of the whole closed
> mesh. `GM` and the GZ curve now come from the same mesh and agree to 0.02 % at 2 deg.
> Evidence: `docs/freeboard_seam.md`, `examples/freeboard_seam/outputs/`,
> `tests/test_wigley_gz_verification.py::TestFreeboardSeamWaterplane`,
> `tests/test_stability.py::test_a_waterline_containing_mesh_edges_measures_the_section_once`.
> The text below is the gap as it stood when opened (m6594) and is kept unedited.

# Wigley Hull GZ Curve Verification Gap

**Date**: 2026-11-08  
**Mission**: m6594 (selfstart generalist seat)  
**Status**: Verification gap identified, not yet closed

## 1. Background

The Wigley parabolic hull form is the standard benchmark for wave resistance and seakeeping comparisons. TensorCAD provides:

- `tensorcad.wigley`: Wigley hull generator with **exact analytical hydrostatics** (Cb=4/9, VCB=5T/8, Cp=Cm=Cwp=2/3 for n=m=2)
- `tensorcad.righting_arm`: GZ curve integrator verified against box barge (wall-sided formula) and cylinder (circular section)

**The gap**: No verification exists that the GZ integrator works correctly on a Wigley hull, despite it being the most common benchmark hull form in the literature.

## 2. What Was Attempted

### 2.1 Direct GZ Computation

```python
from tensorcad.wigley import WigleyParameters, generate_wigley
from tensorcad.righting_arm import righting_arm_curve

p = WigleyParameters(length=10.0, beam=1.0, draft=0.625, nx=81, nz=21)
model = generate_wigley(p)
mesh = model.mesh

# Waterline at z=0 (deck level)
curve = righting_arm_curve(mesh, kg=0.3, waterline_z=0.0, heels_deg=tuple(range(0, 46, 5)))
```

**Result**: GZ curve computed, but small-angle slope does not match GM from waterplane moments:
- GM (waterplane moments): 0.228 m
- GZ/sin(φ) at 5°: 0.091 m

**Discrepancy**: Factor of 2.5×.

### 2.2 Root Cause

The Wigley hull has **no freeboard** by default: the deck is at the design waterline (z=0). When heeled, the waterline immediately rises above the deck on the dipping side, violating the assumption that the hull surface extends above the waterline.

The `righting_arm` module documentation states:
> "A hull generated to its design waterline simply *cannot* be heeled -- on the dipping side the waterline leaves the modelled surface immediately -- so a closed freeboard band has to be added first."

The `with_freeboard` function exists to add this band, but it has issues:
- It removes the waterline cap and adds vertical walls
- The `waterplane_moments` function then fails to find a valid waterplane at z=0
- Error: "the waterline z=0.0 bounds no net area (-13.33 m^2)"

### 2.3 Comparison with Parametric Hull

The parametric hull (`tensorcad.hull`) works correctly with `righting_arm_curve` when `draft=1.0` is used. Investigation shows:
- Parametric hull bounding box: z ∈ [-1.0, 0.0] (draft = depth = 1.0)
- The hull is generated with the waterline at z=0 and the keel at z=-draft
- Tests use `waterline_z=0.0` successfully

**Key difference**: The parametric hull tests use a hull that is floated at a waterline where it has natural freeboard (the hull extends above the waterline). The Wigley hull does not.

## 3. Verification Requirements

To close this gap, the following must be demonstrated:

### 3.1 Small-Angle Consistency

The slope of the GZ curve at φ→0 must match the GM from `initial_stability`:

```
dGZ/dφ|_{φ=0} = GM = KB + BM - KG
```

where:
- KB = vertical centre of buoyancy above keel
- BM = I_T / ∇ (transverse metacentric radius)
- I_T = waterplane second moment
- ∇ = displaced volume

For the Wigley hull with analytical hydrostatics:
- ∇ = 4LBT/9 = 4(10)(1)(0.625)/9 = 2.7778 m³
- A_WP = 2LB/3 = 2(10)(1)/3 = 6.6667 m²
- I_T = ∫ y² dA over waterplane (requires numerical integration)
- KB = 5T/8 = 5(0.625)/8 = 0.3906 m above keel

### 3.2 Volume Conservation

At each heel angle, the displaced volume must equal the upright displacement:

```
|∇(φ) - ∇(0)| / ∇(0) < 1e-10
```

### 3.3 Convergence

The GZ curve must converge under mesh refinement:

```
|GZ_fine - GZ_coarse| → 0 as nx, nz → ∞
```

### 3.4 Analytical Cross-Check (Optional)

For a wall-sided hull form (which the Wigley is not), the GZ curve has a closed form:

```
GZ(φ) = sin(φ) (GM + ½ BM tan²φ)
```

The Wigley hull is **not** wall-sided, so this formula does not apply. However, the small-angle limit must still match GM.

## 4. Proposed Resolution

### 4.1 Add Freeboard to Wigley Hull

Modify the Wigley hull generator to accept a `freeboard` parameter, or use `with_freeboard` correctly:

```python
from tensorcad.gz import with_freeboard

mesh_with_fb = with_freeboard(mesh, freeboard=0.2, waterline_z=0.0)
```

**Issue**: `with_freeboard` currently breaks `waterplane_moments`. This needs investigation.

### 4.2 Alternative: Use Parametric Hull with Wigley-Like Proportions

Generate a parametric hull with L/B=10, B/T=1.6 (Wigley reference proportions) and verify its GZ curve. This does not test the Wigley module specifically, but tests the GZ integrator on a hull with similar proportions.

### 4.3 Write Verification Test

Create `tests/test_wigley_gz.py` with:

1. **Small-angle slope test**: Verify dGZ/dφ|₀ = GM
2. **Volume conservation test**: Verify ∇(φ) = ∇(0) at each heel
3. **Convergence test**: Verify GZ converges under refinement
4. **Positive stability test**: Verify GZ > 0 for small heel (stable equilibrium)

## 5. Current Status

- ✅ Wigley hull generator works (volume error < 0.1% at 81×21)
- ✅ GZ integrator works on parametric hull
- ✅ GZ integrator verified against box barge and cylinder
- ❌ Wigley hull GZ verification **not complete**
- ❌ `with_freeboard` + `waterplane_moments` interaction **broken**

## 6. Recommendations

1. **Priority**: P1 (verification debt for a standard benchmark)
2. **Effort**: 1-2 days (investigate `with_freeboard` issue, write tests)
3. **Dependencies**: None (all modules exist)
4. **Risk**: Low (verification only, no new features)

## 7. References

- Wigley hull: `tensorcad/wigley.py`, `docs/wigley.md`
- GZ integrator: `tensorcad/righting_arm.py`, `docs/gz.md`
- Box barge verification: `tests/test_righting_arm.py::test_box_righting_arm_reproduces_the_wall_sided_closed_form`
- Cylinder verification: `tests/test_righting_arm.py::test_sphere_righting_arm_converges_to_gm_times_sin`
- Parametric hull GZ: `tests/test_righting_arm.py::test_hull_righting_arm_converges_with_refinement`

---

**Note**: This document records a verification gap identified during m6594. The gap has not been closed due to tool call budget constraints. The recommended next step is to investigate the `with_freeboard` + `waterplane_moments` interaction and write the verification tests outlined in §4.3.
