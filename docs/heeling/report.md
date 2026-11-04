# Large-angle stability (GZ) — measured evidence

TensorCAD 0.20.0; wall-sided barge L=10.0 m, B=2.0 m, D=1.5 m, T=0.8 m, KG=0.5 m.

## Wall-sided barge vs the closed form

| heel [deg] | GZ measured [m] | GZ closed form [m] | difference [m] | waterplane [m^2] |
|---:|---:|---:|---:|---:|
| 2 | 0.011060373653 | 0.011060373653 | 1.39e-17 | 20.012190886 |
| 5 | 0.027738300466 | 0.027738300466 | 3.47e-17 | 20.076396751 |
| 10 | 0.056113366873 | 0.056113366873 | -4.86e-17 | 20.308532238 |
| 15 | 0.085830691654 | 0.085830691654 | -4.16e-17 | 20.705523608 |
| 20 | 0.117745730766 | 0.117745730766 | -5.55e-17 | 21.283555450 |
| 25 | 0.152973972824 | 0.152973972824 | 0.00e+00 | 22.067558379 |
| 30 | 0.193055555556 | 0.193055555556 | 0.00e+00 | 23.094010768 |

## Initial slope vs GM

| step [deg] | dGZ/dphi | GM [m] | difference [m] |
|---:|---:|---:|---:|
| 0.8 | 0.316696997083 | 0.316666666667 | 3.033e-05 |
| 0.4 | 0.316674248510 | 0.316666666667 | 7.582e-06 |
| 0.2 | 0.316668562080 | 0.316666666667 | 1.895e-06 |
| 0.1 | 0.316667140517 | 0.316666666667 | 4.739e-07 |

## Curved hull (Wigley at 70 % draught)

| heel [deg] | KN [m] | GZ [m] | displacement [m^3] |
|---:|---:|---:|---:|
| 0 | -0.000000000 | -0.000000000 | 1.560945095486 |
| 5 | -0.014036854 | 0.023004337 | 1.560945095486 |
| 10 | -0.027860565 | 0.045939910 | 1.560945095486 |
| 15 | -0.041225454 | 0.068772640 | 1.560945095486 |
| 20 | -0.053917384 | 0.091441177 | 1.560945095486 |

## Mesh refinement at 10 deg of heel

| grid | faces | GZ(10 deg) [m] | difference from 161x41 [m] | relative |
|---|---:|---:|---:|---:|
| 21x9 | 678 | 0.045833897 | -1.60e-04 | 3.48e-03 |
| 41x13 | 1998 | 0.045939910 | -5.39e-05 | 1.17e-03 |
| 81x25 | 7838 | 0.045983883 | -9.90e-06 | 2.15e-04 |
| 161x41 | 25918 | 0.045993787 | 0.00e+00 | 0.00e+00 |

## Claims

- PASS — barge GZ reproduces the wall-sided closed form at every angle below deck-edge immersion: worst relative difference 1.255e-15 over 7 angles
- PASS — past deck-edge immersion the measured arm falls below the wall-sided claim: at 44.99 deg: measured 0.320276 < wall-sided 0.371097
- PASS — at 90 degrees KN equals half the depth, so GZ = D/2 - KG: KN 0.750000000000 vs D/2 0.750000000000
- PASS — dGZ/dphi at zero reproduces GM on the barge: worst |slope - GM| = 3.033e-05 over steps (0.8, 0.4, 0.2, 0.1)
- PASS — the slope difference is second order in the step (a central difference of a smooth curve): observed orders 2.00, 2.00, 2.00
- PASS — the clipped waterplane area matches the Green-theorem waterplane area: 6.051770833 vs 6.051770833
- PASS — dGZ/dphi at zero reproduces GM on the curved hull: slope 0.263342159 vs GM 0.263340808
- PASS — every heel angle displaces the same volume (a GZ curve is constant displacement): worst relative volume drift 1.423e-16
- PASS — the curved hull's righting arm grows monotonically over the sampled range: GZ -0.000000, 0.023004, 0.045940, 0.068773, 0.091441
- PASS — the curved-hull GZ converges as the mesh refines: relative difference vs the 161x41 grid: 21x9 3.48e-03, 41x13 1.17e-03, 81x25 2.15e-04
