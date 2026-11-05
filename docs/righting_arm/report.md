# Righting-arm curve: measured report

TensorCAD 0.21.0, written by `examples/righting_arm/run_example.py`.

## Claims

| claim | result | detail |
|---|---|---|
| barge displaces B*T*L at the waterline | PASS | 16.000000000 m^3 vs 16.000000000 m^3 |
| barge GZ equals the wall-sided closed form inside its validity window | PASS | worst |GZ - formula| = 1.21e-14 m up to 21.80 deg |
| the wall-sided form is no longer exact past the deck-edge angle | PASS | at 40 deg the formula is -0.0955 m out |
| sphere GZ reproduces GM sin(phi) to better than 1 mm at the finest mesh | PASS | worst |error| = 6.51e-05 m at resolution 32 |
| sphere error falls with refinement | PASS | 12: 5.28e-04, 16: 2.67e-04, 24: 1.13e-04, 32: 6.51e-05 |
| hull small-angle slope equals the waterplane-moment GM | PASS | worst relative difference 6.92e-08 over 4 meshes |
| hull GZ(30 deg) converges as the mesh is refined | PASS | GZ(30) = 0.101859, 0.097583, 0.093945, 0.092795 |
| hull curve has an interior maximum and a positive arm at 30 deg | PASS | peak 0.146677 m at 50 deg; area 0.091849 m rad |
| constant displacement is held along the curve | PASS | worst relative volume residual 2.10e-12 |
| deck-edge immersion is found and reported, not assumed | PASS | the closed-solid arm is a ship's arm only up to 35 deg heel |

## Barge inside the wall-sided window

| heel [deg] | GZ mesh [m] | GZ closed form [m] | difference [m] |
|---|---|---|---|
| 0 | -0.000000000 | 0.000000000 | -1.39e-17 |
| 2 | 0.009315399 | 0.009315399 | +1.89e-16 |
| 5 | 0.023380513 | 0.023380513 | +1.22e-15 |
| 10 | 0.047430958 | 0.047430958 | +7.08e-16 |
| 15 | 0.072889739 | 0.072889739 | -5.36e-15 |
| 20 | 0.100644724 | 0.100644724 | -1.21e-14 |

## Barge beyond it

| heel [deg] | GZ mesh [m] | GZ closed form [m] | formula error [m] | solved waterline [m] |
|---|---|---|---|---|
| 25 | 0.128539596 | 0.131843060 | -0.0033 | +0.0023 |
| 30 | 0.147825274 | 0.168055556 | -0.0202 | +0.0141 |
| 35 | 0.160989961 | 0.211541036 | -0.0506 | +0.0342 |
| 40 | 0.170155631 | 0.265697355 | -0.0955 | +0.0616 |

## Sphere

| resolution | worst |GZ - GM sin(phi)| [m] |
|---|---|
| 12 | 5.283e-04 |
| 16 | 2.666e-04 |
| 24 | 1.133e-04 |
| 32 | 6.514e-05 |

## Hull ladder

| mesh | faces | displacement [m^3] | GM (waterplane) [m] | slope [m] | rel. diff | GZ(30) [m] |
|---|---|---|---|---|---|---|
| 21x9 | 678 | 3.486834 | 0.134729 | 0.134729 | +6.18e-08 | 0.101859 |
| 41x13 | 1998 | 3.552393 | 0.126776 | 0.126776 | +6.52e-08 | 0.097583 |
| 81x25 | 7838 | 3.597980 | 0.120299 | 0.120299 | +6.88e-08 | 0.093945 |
| 161x41 | 25918 | 3.611760 | 0.118228 | 0.118228 | +6.92e-08 | 0.092795 |

## Hull curve (41x13, waterline -0.47 m, KG 0.55 m)

| heel [deg] | GZ [m] | waterline [m] | y_B [m] | volume residual |
|---|---|---|---|---|
| 0 | +0.000000 | -0.47000 | +0.00000 | +5.8e-13 |
| 5 | +0.011207 | -0.47133 | +0.01295 | -3.2e-13 |
| 10 | +0.023159 | -0.47534 | +0.02663 | -3.5e-13 |
| 15 | +0.036700 | -0.48210 | +0.04188 | -3.7e-13 |
| 20 | +0.052849 | -0.49174 | +0.05969 | -1.0e-12 |
| 25 | +0.072913 | -0.50447 | +0.08137 | -3.4e-13 |
| 30 | +0.097583 | -0.52016 | +0.10758 | -2.1e-12 |
| 35 | +0.118554 | -0.53542 | +0.13003 | +6.0e-13 |
| 40 | +0.133146 | -0.54865 | +0.14600 | +6.9e-14 |
| 45 | +0.142167 | -0.55969 | +0.15631 | +1.7e-12 |
| 50 | +0.146677 | -0.56856 | +0.16200 | +6.9e-13 |
| 60 | +0.144216 | -0.58091 | +0.16154 | -1.8e-12 |

Peak 0.146677 m at 50 deg; dynamical stability 0.091849 m rad; deck-edge immersion from the solved waterlines: 35 deg.

## Figures

- `barge_vs_closed_form.png`
- `sphere_vs_analytic.png`
- `hull_curve.png`
