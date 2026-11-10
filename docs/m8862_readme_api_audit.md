# TensorCAD README/API Consistency Audit (m8862)

**Date:** 2026-11-10
**Auditor:** agent/m8862
**Scope:** README.md code examples vs actual `tensorcad` package API (≥10 items)

## Methodology

Each README code example was executed against the current codebase (`__version__ = "0.23.0"`) and the output compared against the claimed values. Examples that run and produce matching output are marked ✅; mismatches are marked ❌.

## Audit Results (14 items checked)

| # | README Claim | Status | Detail |
|---|---|---|---|
| 1 | `generate_hull(HullParameters(...))` → `model.hydrostatics.displacement_volume` = 7.446554 | ✅ | Verified: 7.4465540971492254 |
| 2 | `model.hydrostatics.block_coefficient` = 0.4654 | ✅ | Verified: 0.4654096310718266 |
| 3 | `model.hydrostatics.prismatic_coefficient` = 0.7031 | ✅ | Verified: 0.7030834593151094 |
| 4 | `model.validation.warnings` = `[]` at nx=41 | ✅ | Verified: `[]` |
| 5 | `analytic_hydrostatics(length=10.0, beam=2.0, draft=0.8).vertical_centre_buoyancy` = -0.32 | ❌ | **API mismatch**: function accepts `HullParameters`, `(L,B,T)` tuple, or `None` — not keyword args `length`/`beam`/`draft` |
| 6 | `initial_stability(model.mesh, kg=0.6)` → `bm_transverse` = 0.4386 | ✅ | Verified: 0.43861787815566544 |
| 7 | `initial_stability` → `gm_transverse` = 0.3216 (analytic 0.315360, +1.99%) | ✅ | Verified: 0.3216229964990912 |
| 8 | `initial_stability` → `waterplane.longitudinal_second_moment` = 76.41 | ✅ | Verified: 76.40662200112666 |
| 9 | `SphereField(radius=0.5).signed_distance([(0.9,0,0)])` = 0.4 | ✅ | Verified: `[0.4]` |
| 10 | `to_mesh(sphere, 40)` → `report.watertight`=True, `euler_characteristic`=2 | ✅ | Verified |
| 11 | `SuboffHullField().length, .max_radius` = 4.35610, 0.25400 | ✅ | Verified: 4.3561001016, 0.25399998984000005 |
| 12 | `SuboffHullField().signed_distance([(2.0, 0.30, 0.0)])` = +0.0460 | ✅ | Verified: 0.04600001015999994 |
| 13 | `field_from_object(...).quantity` = 'truncated_signed_distance' | ✅ | Verified |
| 14 | Box test worst error = `2.2e-16` | ✅ | Verified: 2.220e-16 |

## Inconsistencies Found (3 items, ranked by severity)

### ❌ #1 (Critical): `analytic_hydrostatics` API signature mismatch

**README says:**
```python
analytic_hydrostatics(length=10.0, beam=2.0, draft=0.8).vertical_centre_buoyancy  # -0.32
```

**Actual API:**
```python
def analytic_hydrostatics(parameters: Any = None, *, breadth_exponent, quadrature_levels, water_density) -> AnalyticHydrostatics
```

The function accepts:
- A `HullParameters` object
- An `(L, B, T)` tuple/list
- `None` for defaults

It does **not** accept `length`, `beam`, `draft` as keyword arguments. Calling as shown in the README raises `TypeError: analytic_hydrostatics() got an unexpected keyword argument 'length'`.

**Fix:** Update README example to use tuple syntax: `analytic_hydrostatics((10.0, 2.0, 0.8))`

### ❌ #2 (Moderate): Hemisphere worst error overstated as 0.16%

**README says:**
> hemisphere `R = 1 m` | ... | `0.16 %` (inscribed n-gon)

**Actual measurement** (n_lat=32, n_lon=64):
- Waterplane area: -0.16% ✅
- Transverse second moment (I_T): -0.32% ❌
- Longitudinal second moment (I_L): -0.32% ❌
- KB: +0.012%
- BM_T: -0.10%
- KM_T: -0.03%

**Worst error is 0.32%** (on I_T and I_L), not 0.16%. The 0.16% figure is the waterplane area error only.

**Fix:** Update README table to report the actual worst error (0.32%).

### ⚠️ #3 (Minor): KM relative error discrepancy

**README says:**
> the default grid shows `+1.99 %` where `KM` shows `0.71 %`

**Actual measurement:**
- KM relative error: 0.68% (not 0.71%)
- GM relative error: 1.99% ✅

The 0.71% figure appears to be a rounding artifact from an earlier measurement. The difference (0.03 percentage points) is small but the README should be accurate.

**Fix:** Update README to say "0.68%" instead of "0.71%".

## Example Runnability

All three shipped examples pass:
- `examples/parametric_hull/run_example.py` — ✅ all checks passed
- `examples/stability/run_example.py` — ✅ all checks passed
- `examples/suboff/run_example.py` — ✅ all checks passed

## Conclusion

12 of 14 API claims are correct. 2 are incorrect (items #1 and #2), and 1 has a minor numerical discrepancy (item #3). The fixes below correct all three.
