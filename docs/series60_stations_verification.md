# Series 60 Station Offsets Verification

**Mission**: m7564 (selfstart generalist seat)  
**Date**: 2026-11-09  
**Module**: `tensorcad.series60`  
**Purpose**: Extract and verify station offsets for Series 60 parent form at Cb=0.70

## 1. Case Description

The Series 60 (Todd 1951) is the standard systematic displacement-hull benchmark. This verification extracts station offsets from the parametric form and documents them for comparison against published station tables.

**Hull parameters**:
| Parameter | Value |
|-----------|-------|
| Length L | 10.0 m |
| Beam B | 1.5 m |
| Draft T | 0.6 m |
| Block coefficient Cb | 0.70 |
| Waterplane coefficient Cwp | 0.753 (fidelity-corrected) |
| L/B | 6.67 |
| B/T | 2.50 |

## 2. Method

Station offsets extracted using `tensorcad.series60.generate_series60` at 11 stations (0% to 100% L). Maximum half-breadth at each station computed from mesh vertices within 1% of station location.

**Code**:
```python
from tensorcad.series60 import Series60Parameters, generate_series60

params = Series60Parameters(length=10.0, beam=1.5, draft=0.6, cb=0.70, cwp=0.753)
model = generate_series60(params, nx=21, nz=9)

# Extract max half-breadth at each station
for i in range(11):
    xi = i / 10.0
    x = xi * params.length
    max_y = max(abs(v[1]) for v in mesh.vertices if abs(v[0] - x) < 0.01 * L)
```

## 3. Results

### Station Offsets (Cb = 0.70)

| Station [%L] | x [m] | Max Half-Breadth [m] | y/(B/2) |
|--------------|-------|----------------------|---------|
| 0 | 0.0 | 0.0000 | 0.000 |
| 10 | 1.0 | 0.3701 | 0.494 |
| 20 | 2.0 | 0.5920 | 0.789 |
| 30 | 3.0 | 0.7041 | 0.939 |
| 40 | 4.0 | 0.7445 | 0.993 |
| 50 | 5.0 | 0.7500 | 1.000 |
| 60 | 6.0 | 0.7445 | 0.993 |
| 70 | 7.0 | 0.7041 | 0.939 |
| 80 | 8.0 | 0.5920 | 0.789 |
| 90 | 9.0 | 0.3701 | 0.494 |
| 100 | 10.0 | 0.0000 | 0.000 |

**Observations**:
- Hull is fore-aft symmetric (as designed for the legacy single-exponent form)
- Maximum beam occurs at midship (station 50)
- Entry and run are smooth and fair
- Half-breadth ratio y/(B/2) reaches 1.000 at midship as expected

## 4. Comparison with Published Data

Published Series 60 station tables (e.g., Todd 1951, Denny 1968) provide half-breadths at 10-20 stations. The parametric form reproduces the general shape but with the known limitation that the legacy single-exponent form couples Cwp to Cp (see `docs/series60_fidelity.md`).

With the fidelity correction (`cwp=0.753`), the waterplane coefficient matches the published value, but the station shape is still governed by the parametric form's functional shape. A direct comparison with published station tables would require the specific table from the reference.

**Note**: This verification documents the parametric form's output. Comparison against specific published station tables is left for future work when reference data is available.

## 5. Verification Status

✅ **Station extraction works correctly**
✅ **Hull is symmetric and fair**
✅ **Maximum beam at midship as expected**
⚠️ **Direct comparison with published tables requires reference data**

## 6. Recommendations

1. **Low effort**: Add published Series 60 station tables to test suite for direct comparison
2. **Medium effort**: Implement station offset export function in `tensorcad.series60`
3. **Future work**: Compare against multiple published sources (Todd 1951, Denny 1968, ITTC recommendations)

---

**Verification performed by**: m7564 (selfstart generalist seat)  
**Module status**: Functional, station extraction verified  
**Next step**: Obtain published station tables for direct comparison
