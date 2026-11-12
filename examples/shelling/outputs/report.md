# Shelling: exact erosion, measured wall, loud refusals

## 1. Box (GAP-CAD-011 acceptance criterion)

| t/L | inner volume [m^3] | analytic (1-2t)^3 [m^3] | relative error | wall [m] |
|---:|---:|---:|---:|---:|
| 0.01 | 0.941191999999972 | 0.941192000000000 | 2.97e-14 | 0.010000000000 |
| 0.02 | 0.884735999999974 | 0.884736000000000 | 2.97e-14 | 0.020000000000 |
| 0.03 | 0.830583999999975 | 0.830584000000000 | 2.99e-14 | 0.030000000000 |
| 0.04 | 0.778687999999977 | 0.778688000000000 | 3.01e-14 | 0.040000000000 |
| 0.05 | 0.728999999999978 | 0.729000000000000 | 3.05e-14 | 0.050000000000 |
| 0.06 | 0.681471999999980 | 0.681472000000000 | 3.00e-14 | 0.060000000000 |
| 0.07 | 0.636055999999981 | 0.636056000000000 | 3.02e-14 | 0.070000000000 |
| 0.08 | 0.592703999999982 | 0.592704000000000 | 2.96e-14 | 0.080000000000 |
| 0.09 | 0.551367999999984 | 0.551368000000000 | 2.98e-14 | 0.090000000000 |
| 0.10 | 0.511999999999985 | 0.512000000000000 | 3.04e-14 | 0.100000000000 |

Worst relative error over the sweep: **3.05e-14** (threshold 2%).

## 2. Pentagon prism vs independent 2-D clipping

| t [m] | inner volume [m^3] | clipped reference [m^3] | relative error | eroded corners |
|---:|---:|---:|---:|---:|
| 0.02 | 1.100325623474289 | 1.100325623474316 | 2.46e-14 | 5 |
| 0.05 | 0.923302866069921 | 0.923302866069944 | 2.49e-14 | 5 |
| 0.10 | 0.671989350204373 | 0.671989350204390 | 2.48e-14 | 5 |
| 0.15 | 0.470598690123310 | 0.470598690123322 | 2.45e-14 | 5 |

Worst relative error: **2.49e-14**.

## 3. Sphere, t = 0.2 m

| meridians x parallels | faces | shell volume [m^3] | analytic [m^3] | rel. error | observed order |
|---|---:|---:|---:|---:|---:|
| 12 x 6 | 120 | 1.905196 | 2.044130 | 6.797% |  |
| 24 x 12 | 528 | 2.008892 | 2.044130 | 1.724% | 1.98 |
| 48 x 24 | 2208 | 2.035288 | 2.044130 | 0.433% | 1.99 |

## 4. The wall measurement, falsified on purpose

| offset applied [m] | measured wall [m] | relative error |
|---:|---:|---:|
| 0.020 | 0.020000000000 | 1.02e-14 |
| 0.050 | 0.050000000000 | 9.99e-15 |
| 0.070 | 0.070000000000 | 1.01e-14 |
| 0.100 | 0.100000000000 | 9.99e-15 |

## 5. Refusals

| case | refused | reason |
|---|---|---|
| plate_0.1_shelled_at_0.08 | True | the inward offset of 0.08 m inverts the solid: the inner signed volume is -0.042335999999998736 m^3 (outer 0.09999999999999999 m^3), so the thickness exceeds half the local wall somewhere. Reduce the thickness or split the body |
| taper_shelled_at_0.10 | True | the inward offset folds the surface: min(dot(n_inner, n_outer)) = -1.0 over the faces, so some offset faces are wound backwards |

Hollow export: watertight True, material 0.488000000000 m^3, cavity 0.512000000000 m^3, signed volume of the exported surface 0.488000000000 m^3.

![shelling](shelling.png)

## 6. Error judgement

| key quantity | reference | error | below 3 %? |
|---|---|---|---|
| box inner volume, t/L 0.01-0.10 | (L-2t)^3 | 3.05e-14 | yes |
| pentagon prism inner volume | 2-D half-plane clipping | 2.49e-14 | yes |
| sphere shell volume, 48 x 24 | 4/3 pi (R^3-(R-t)^3) | 0.43% | yes |
| sphere shell volume, 12 x 6 | 4/3 pi (R^3-(R-t)^3) | 6.80% | no (faceting, converges at order ~2) |

No quantisation of the mesh, no smoothing of the inner surface and no tolerance tuning is applied in any of the numbers above.
