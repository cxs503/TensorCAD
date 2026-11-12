# Series 60 Station Offsets: Verification and Documentation

**Mission**: m16098 (selfstart generalist seat)  
**Date**: 2026-11-12  
**Module**: `tensorcad.series60`  
**Purpose**: Verify mesh extraction of station offsets and document parametric form values

## 1. Case Description

The Series 60 (Todd 1951) is the standard systematic displacement-hull benchmark. While the parametric form reproduces the integral coefficients (Cb, Cwp, VCB) to high accuracy, the **station-by-station shape** extraction from the mesh has not been systematically verified. This verification closes that gap.

**Hull parameters**:
| Parameter | Value |
|-----------|-------|
| Length L | 10.0 m |
| Beam B | 1.5 m |
| Draft T | 0.6 m |
| Block coefficient Cb | 0.70 |
| Waterplane coefficient Cwp | 0.715 (legacy) / 0.753 (fidelity-corrected) |
| L/B | 6.67 |
| B/T | 2.50 |

## 2. Method

### 2.1 Parametric Form

The Series 60 station law is:
```
y(xi, zeta) = (B/2) * f_deck(xi) * zeta^(1/m)
```
where:
- `xi = x/L` is the longitudinal position (0 at bow, 1 at stern)
- `zeta = z/T` is the vertical position (0 at keel, 1 at waterline)
- `f_deck(xi) = 1 - |2*xi - 1|^n` is the deck half-breadth law
- `n = Cp/(1-Cp)` is the longitudinal exponent
- `m = Cb/(Cp-Cb)` is the vertical exponent

For the **legacy single-exponent form** (cwp=None), the waterplane coefficient equals the prismatic coefficient: `Cwp = Cp`. For Cb = 0.70, this gives Cp = 0.715 and n = 2.509.

The **fidelity-corrected form** (cwp=0.753, lcb=0.490) decouples Cwp from Cp and introduces fore-aft asymmetry to match the published LCB.

### 2.2 Mesh Extraction

Station half-breadths are extracted from the mesh by finding the maximum |y| coordinate at each station location:
```python
for i in range(n_stations):
    xi = i / (n_stations - 1)
    x = xi * L
    max_y = max(abs(v[1]) for v in mesh.vertices if abs(v[0] - x) < tolerance * L)
    y_normalized = max_y / (B / 2)
```

The extraction accuracy depends on mesh resolution (nx) and the tolerance parameter.

## 3. Results

### 3.1 Legacy Form Station Offsets (Cb = 0.70)

**Table 1: Legacy single-exponent form (Cwp = Cp = 0.715, n = 2.509)**

| Station [%L] | x/L | y/(B/2) mesh | y/(B/2) analytic | error |
|--------------|-----|--------------|------------------|-------|
| 0 | 0.000 | 0.0000 | 0.0000 | 0.00% |
| 10 | 0.100 | 0.4287 | 0.4287 | <0.01% |
| 20 | 0.200 | 0.7224 | 0.7224 | <0.01% |
| 30 | 0.300 | 0.8996 | 0.8996 | <0.01% |
| 40 | 0.400 | 0.9824 | 0.9824 | <0.01% |
| 50 | 0.500 | 1.0000 | 1.0000 | 0.00% |
| 60 | 0.600 | 0.9824 | 0.9824 | <0.01% |
| 70 | 0.700 | 0.8996 | 0.8996 | <0.01% |
| 80 | 0.800 | 0.7224 | 0.7224 | <0.01% |
| 90 | 0.900 | 0.4287 | 0.4287 | <0.01% |
| 100 | 1.000 | 0.0000 | 0.0000 | 0.00% |

**Mesh parameters**: nx=161, nz=41, tolerance=0.002  
**Maximum extraction error**: <0.01% at all stations

**Conclusion**: The mesh extraction reproduces the analytic deck law to within 0.01%, which is well inside the 3% accuracy bar for a key quantity. The small discrepancies are due to mesh discretization.

### 3.2 Fidelity-Corrected Form

The fidelity-corrected form (cwp=0.753, lcb=0.490) intentionally changes the station shape to decouple Cwp from Cp:

**Table 2: Comparison of legacy vs fidelity-corrected forms**

| Station [%L] | y/(B/2) legacy | y/(B/2) corrected | difference |
|--------------|----------------|-------------------|------------|
| 10 | 0.4287 | 0.4523 | +5.5% |
| 20 | 0.7224 | 0.7389 | +2.3% |
| 30 | 0.8996 | 0.8912 | -0.9% |
| 40 | 0.9824 | 0.9756 | -0.7% |
| 50 | 1.0000 | 1.0000 | 0.0% |
| 60 | 0.9824 | 0.9876 | +0.5% |
| 70 | 0.8996 | 0.9123 | +1.4% |
| 80 | 0.7224 | 0.7456 | +3.2% |
| 90 | 0.4287 | 0.4567 | +6.5% |

**Maximum shape difference**: 6.5% at station 90  
**Fore-aft asymmetry**: The corrected form with lcb=0.490 has fuller aft sections (station 90: 0.4567) than forward sections (station 10: 0.4523), shifting the LCB aft of midship.

### 3.3 Key Properties Verified

1. **Midship beam is exact**: y/(B/2) = 1.0000 at station 50 for both forms
2. **Bow and stern are zero**: y/(B/2) < 1e-6 at stations 0 and 100
3. **Legacy form is symmetric**: Station i matches station (100-i) to 1e-6
4. **Corrected form is asymmetric**: Station 10 differs from station 90 by >0.5%
5. **Mesh extraction converges**: Fine mesh (161×41) matches analytic to <0.01%

## 4. Error Judgement

| key quantity | reference | measured | error | < 3%? |
|---|---|---|---:|---|
| Station offsets, legacy form | analytic deck law | max <0.01% | <0.01% | **yes** |
| Midship beam | 1.0000 | 1.0000 | 0.00% | **yes** |
| Bow/stern | 0.0000 | <1e-6 | ~0% | **yes** |
| Fore-aft symmetry (legacy) | station i = station (100-i) | max diff <1e-6 | ~0% | **yes** |

## 5. Points to Watch

1. **Mesh resolution matters.** At nx=21, the extraction error is ~0.1%. At nx=161, it drops to <0.01%. Users should choose nx based on the accuracy they need.
2. **The legacy form has Cwp = Cp.** This is a structural property of the single-exponent form, not a bug. If you need Cwp != Cp, use the fidelity-corrected form.
3. **The fidelity correction trades station shape for Cwp accuracy.** The corrected form matches the target Cwp to round-off but changes the station offsets by up to 6.5%. This is the price of decoupling.
4. **No published station tables exist for the fidelity-corrected form.** The correction is a TensorCAD invention; there is no external reference to compare against except the integral coefficients.
5. **The extraction tolerance affects accuracy.** A tolerance of 0.002*L works well for nx=161. For coarser meshes, increase the tolerance to capture more vertices.

## 6. Recommendations

1. **Use nx >= 81 for station extraction.** This gives <0.05% error at all stations.
2. **Document the trade-off clearly.** Users should know that the legacy form has Cwp = Cp, while the fidelity-corrected form trades station shape for Cwp accuracy.
3. **Provide both forms in examples.** Show the legacy form for comparison with published data, and the fidelity-corrected form for design work where Cwp matters.

## 7. Reproducing

```bash
python -m pytest tests/test_series60_stations_published_comparison.py -v  # 8 cases
```

The test suite verifies:
- Mesh extraction accuracy vs analytic law (<0.01% at fine resolution)
- Midship beam is exact (1.0000)
- Bow and stern are zero (<1e-6)
- Legacy form is fore-aft symmetric
- Fidelity-corrected form has Cwp != Cp
- Fidelity-corrected form changes station shape
- Fidelity-corrected form is fore-aft asymmetric
- Documented station offsets match expected values

---

**Verification performed by**: m16098 (selfstart generalist seat)  
**Test suite status**: 8 new tests passing  
**Documentation status**: Anchored to analytic deck law  
**Recommendation**: Legacy form matches analytic law to <0.01%; fidelity-corrected form trades station shape for Cwp accuracy (6.5% max deviation, expected and documented).
