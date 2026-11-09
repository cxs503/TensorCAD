# TensorCAD Verification Status Summary

**Mission**: m7564 (selfstart generalist seat)  
**Date**: 2026-11-09  
**Purpose**: Comprehensive audit of verification coverage and remaining gaps

## Executive Summary

TensorCAD has achieved comprehensive verification coverage across all major hull families and hydrostatic modules. All 66 test files pass (1500+ test cases), with verification documentation anchored to analytical solutions, published literature, and independent numerical routes. This document summarizes the verification state and identifies remaining gaps for future work.

## 1. Verified Modules (All Tests Passing)

### 1.1 Hull Families

| Module | Verification Basis | Key Results | Status |
|--------|-------------------|-------------|--------|
| `tensorcad.hull` | Parametric hull with analytical volume | Volume converges to analytical solution (0.14% at 161×41) | ✅ Complete |
| `tensorcad.wigley` | Exact analytical hydrostatics (Cb=4/9, VCB=5T/8) | All coefficients match to machine precision | ✅ Complete |
| `tensorcad.series60` | Published parent forms (Todd 1951) | Cwp matches published values to <0.01% with fidelity fix | ✅ Complete |
| `tensorcad.suboff` | DARPA offset table (Groves et al. 1989) | Volume 0.699211 m³ vs literature 0.699211 m³ (0.000%) | ✅ Complete |
| `tensorcad.planing` | Closed-form hydrostatics for constant-deadrise hull | Wetted area and volume match closed forms | ✅ Complete |
| `tensorcad.asymmetric_hull` | Port/starboard asymmetric loft | Watertight, manifold, correct dimensions | ✅ Complete |

### 1.2 Hydrostatics & Stability

| Module | Verification Basis | Key Results | Status |
|--------|-------------------|-------------|--------|
| `tensorcad.hydrostatics` | Divergence theorem on mesh | Box barge: exact to round-off | ✅ Complete |
| `tensorcad.stability` | Waterplane moments (Green's theorem) | Box/ hemisphere anchors match closed forms | ✅ Complete |
| `tensorcad.righting_arm` | Wall-sided barge formula, circular cylinder | GZ(φ) matches closed forms to <1e-15 m | ✅ Complete |
| `tensorcad.gz` | Large-angle stability with heeled displacement | Volume conservation <1e-13 relative error | ✅ Complete |
| `tensorcad.trim` | TPC/MCT from two independent routes | Box hull: exact; parametric hull: <0.6% | ✅ Complete |
| `tensorcad.damage` | Lost buoyancy method | Compartment flooding matches analytical | ✅ Complete |
| `tensorcad.free_surface` | Exact heel-dependent tank correction | Matches first-order at small heel | ✅ Complete |
| `tensorcad.bonjean` | Section properties along hull | Simpson integration exact for linear areas | ✅ Complete |

### 1.3 Geometry & Mesh Quality

| Module | Verification Basis | Key Results | Status |
|--------|-------------------|-------------|--------|
| `tensorcad.intersect` | Exact rational arithmetic oracle | 0 logic errors in 4000 pairs; resolution contract verified | ✅ Complete |
| `tensorcad.nurbs` | Rational B-spline theory (Piegl & Tiller) | Cylinder area exact to 1e-15; sphere volume converges as O(h²) | ✅ Complete |
| `tensorcad.voxel` | SDF voxelization | SUBOFF volume 0.698982 m³ vs literature 0.699211 m³ (-0.033%) | ✅ Complete |
| `tensorcad.curvature` | Differential geometry | Principal curvatures match analytical for primitives | ✅ Complete |
| `tensorcad.mesh_conditioning` | Mesh quality metrics | Sliver detection and reporting | ✅ Complete |

### 1.4 Specialized Modules

| Module | Verification Basis | Key Results | Status |
|--------|-------------------|-------------|--------|
| `tensorcad.icechannel` | Spherical trigonometry | Corridor distortion <0.5% at 1061 km scale | ✅ Complete |
| `tensorcad.pc6` | PC6 hull particulars | L/B/Cb gate verified | ✅ Complete |
| `tensorcad.statutory` | Intact stability criteria | Criteria evaluation and reporting | ✅ Complete |

## 2. Verification Methodology

### 2.1 Multi-Route Verification

Critical quantities are computed by two or more independent routes:
- **Volume**: divergence theorem on mesh vs. analytical closed forms vs. tanh-sinh quadrature
- **Waterplane moments**: Green's theorem on waterplane loops vs. station sum
- **GZ curve**: heeled displacement integration vs. wall-sided formula (box barge)
- **TPC/MCT**: waterplane area/moments vs. finite difference of displacement

### 2.2 Analytical Anchors

Every module with closed-form solutions is verified against them:
- **Wigley**: Cb=4/9, Cp=Cm=Cwp=2/3, VCB=5T/8 (exact rationals)
- **Series 60**: Published parent forms (Todd 1951) with Cwp fidelity fix
- **SUBOFF**: DARPA offset table (Groves et al. 1989)
- **Box barge**: Wall-sided formula exact to round-off
- **Sphere**: Area=4πR², Volume=4/3πR³

### 2.3 Convergence Studies

Mesh integrals are verified to converge under refinement:
- **Parametric hull**: Volume error 5.29% → 0.14% as nx goes 11 → 161
- **Wigley**: All coefficients converge to analytical values
- **NURBS sphere**: Volume converges as O(h²)
- **Righting arm**: GZ(φ) converges with mesh refinement

## 3. Remaining Gaps & Future Work

### 3.1 High Priority

1. **OCP B-rep Backend**
   - **Gap**: No industrial-standard B-rep geometry (STEP/IGES import/export)
   - **Impact**: Cannot interface with CAD systems or perform boolean operations
   - **Effort**: Large (requires OCP integration, multi-patch NURBS, trimmed surfaces)
   - **Recommendation**: Dedicated mission with CAD kernel expertise

2. **Multi-Patch NURBS Stitching**
   - **Gap**: Single-patch NURBS only; no trimmed patches or multi-patch models
   - **Impact**: Cannot represent complex hull forms (bulbous bow, transom corners)
   - **Effort**: Medium-large (requires surface-surface intersection, reparameterization)
   - **Recommendation**: Build on existing single-patch infrastructure

### 3.2 Medium Priority

3. **Mesh Smoothing/Simplification**
   - **Gap**: No mesh decimation or smoothing algorithms
   - **Impact**: Cannot reduce mesh size for downstream solvers
   - **Effort**: Medium (requires quadric error metrics, edge collapse)
   - **Recommendation**: Integrate existing library (trimesh/pyvista) or implement minimal version

4. **Resistance Prediction**
   - **Gap**: No resistance modules (ITTC-57, Savitsky, Holtrop-Mennen)
   - **Impact**: Cannot predict hull resistance from geometry
   - **Effort**: Medium (empirical formulas, friction line correlations)
   - **Recommendation**: Implement ITTC-57 + form factor method first

5. **Series 60 Station Tables**
   - **Gap**: No tabulated station offsets for Series 60 parents
   - **Impact**: Cannot compare against published station tables
   - **Effort**: Small (extract stations from parametric form, compare with literature)
   - **Recommendation**: Low-hanging fruit for additional verification

### 3.3 Low Priority

6. **Neural Field Shape Family Expansion**
   - **Gap**: Neural fields limited to sphere + dent edits
   - **Impact**: Cannot train neural representations of hull forms
   - **Effort**: Large (requires training pipeline, shape priors)
   - **Recommendation**: Defer until use case is clear

7. **Fairness/Surface Lightness Checks**
   - **Gap**: No curvature continuity metrics beyond basic curvature
   - **Impact**: Cannot assess design quality (fairness)
   - **Effort**: Medium (requires curvature analysis, visualization)
   - **Recommendation**: Build on existing curvature module

## 4. Test Suite Statistics

| Metric | Value |
|--------|-------|
| Total test files | 66 |
| Total test cases | ~1500+ |
| Pass rate | 100% |
| Wall time (full suite) | ~600s |
| Modules with verification docs | 20+ |
| Analytical anchors | 15+ |
| Multi-route verifications | 8+ |

## 5. Verification Discipline

### 5.1 Evidence Chain

Every verification follows the pattern:
1. **Code** (`tensorcad/*.py`) → **Example** (`examples/*/run_example.py`) → **Artifact** (`artifacts/*/metrics.json`) → **Documentation** (`docs/*.md`) → **Test** (`tests/test_*.py`)
2. Tests recompute values rather than reading artifacts, so merge gates cannot turn them red
3. All numbers are program-injected (zero hand-copying)

### 5.2 Error Reporting

- **Honest FAIL**: Defects are reported, not covered
- **Resolution contracts**: Predicates state their resolution explicitly (e.g., `intersect.RELATIVE_TOLERANCE`)
- **Convergence studies**: Mesh integrals report error vs. refinement
- **Multi-route residuals**: Independent routes are compared, residual is reported

### 5.3 Limitations Stated

Every module docstring states what it is **not**:
- No free-surface correction in stability (caller must subtract tank losses)
- No damage stability in righting_arm (separate module)
- No resistance prediction in planing (geometry only)
- No class-rule compliance in statutory (criteria evaluation only)

## 6. Recommendations for Next Selfstart

1. **OCP B-rep backend**: Highest impact, largest effort. Requires dedicated mission with CAD kernel expertise.
2. **Mesh smoothing/simplification**: Medium impact, medium effort. Good candidate for a focused mission.
3. **Resistance prediction**: High impact for ship design, medium effort. ITTC-57 + form factor is a good starting point.
4. **Series 60 station tables**: Low effort, additional verification. Good warm-up task.

## 7. Conclusion

TensorCAD has achieved a mature verification state with comprehensive coverage of hull families, hydrostatics, stability, and geometry quality. All major modules are verified against analytical solutions, published literature, and independent numerical routes. The remaining gaps are primarily in industrial interoperability (B-rep/STEP) and downstream applications (resistance, structural analysis), which are natural next steps for the ecosystem.

The verification discipline is strong: multi-route checks, analytical anchors, convergence studies, and honest error reporting. The codebase is ready for integration with downstream solvers (TensorBEM, TensorFEM, TensorLBM) and ship design tools (TensorShipDesign).

---

**Verification audit performed by**: m7564 (selfstart generalist seat)  
**Test suite status**: All 66 test files passing  
**Documentation status**: 20+ verification documents anchored to analytical solutions  
**Recommendation**: TensorCAD is production-ready for hull geometry and hydrostatics; B-rep backend is the next major milestone.
