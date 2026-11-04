# Generated tables (do not edit by hand)

Written by `examples/gz/run_example.py`. Geometry: barge `10.0 x 2.0 x 0.8 m` with `0.4 m` freeboard and `KG = 0.6 m`; cylinder `R = 1.0 m`, `L = 4.0 m`, `G` `0.3 R` below the axis; hull `10.0 x 2.0 x 0.8 m` at waterline `z = -0.11 m` (0.11 m freeboard) and `KG = 0.7 m`.

## Barge against the wall-sided formula

| heel [deg] | GZ measured [m] | GZ wall-sided [m] | relative error |
|---:|---:|---:|---:|
| 1 | 0.003782463 | 0.003782463 | -7.96e-14 |
| 5 | 0.019022726 | 0.019022726 | -7.57e-14 |
| 10 | 0.038748549 | 0.038748549 | +5.93e-14 |
| 15 | 0.059948787 | 0.059948787 | +8.90e-14 |
| 20 | 0.083543716 | 0.083543716 | +1.24e-13 |

Deck-edge immersion is at `tan(phi) = 2 f / B = 21.80 deg`; past it the formula has no
standing, and the measured arm falls below it:

| heel [deg] | GZ measured [m] | GZ wall-sided [m] | relative error |
|---:|---:|---:|---:|
| 25 | 0.107408683 | 0.110712147 | -2.98e-02 |
| 30 | 0.122825274 | 0.143055556 | -1.41e-01 |
| 40 | 0.138016251 | 0.233557975 | -4.09e-01 |

## Barge righting-arm curve

| quantity | value |
|---|---:|
| maximum arm | 0.138016 m at 40 deg |
| vanishing angle | 90.000 deg |
| area under the curve to 30 deg | 0.032296 m*rad |
| area under the curve to 40 deg | 0.055059 m*rad |
| GM from the slope at 5 deg | 0.218261 m |
| worst displacement residual | 9.2e-13 (relative) |

## Circular cylinder against the exact arm

| heel [deg] | GZ exact [m] | GZ 48-gon [m] |
|---:|---:|---:|
| 0 | 0.000000000 | 0.000000000 |
| 10 | 0.052094453 | 0.051700619 |
| 20 | 0.102606043 | 0.101895792 |
| 30 | 0.150000000 | 0.148929462 |
| 45 | 0.212132034 | 0.210618064 |
| 60 | 0.259807621 | 0.257953394 |
| 75 | 0.289777748 | 0.287709627 |
| 90 | 0.300000000 | 0.297858923 |
| 105 | 0.289777748 | 0.287709627 |
| 120 | 0.259807621 | 0.257953394 |

| segments | worst relative error over 5..120 deg |
|---:|---:|
| 24 | 3.53e-02 |
| 48 | 7.14e-03 |
| 96 | 1.89e-03 |

## Parametric hull with freeboard

| grid | faces | displacement [m3] | BM [m] | GM (waterplane) [m] | GM (curve slope) [m] | relative error | GZ(10 deg) [m] | vanishing angle [deg] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 41x13 | 1998 | 5.954490 | 0.439028 | 0.156214 | 0.156315 | +6.50e-04 | 0.022845 | 20.86 |
| 81x25 | 7838 | 5.995478 | 0.436746 | 0.151946 | 0.152038 | +6.07e-04 | 0.022106 | 20.49 |

## Claims

| claim | holds | detail |
|---|---|---|
| wall-sided barge matches the closed form to 1e-9 | yes | worst 1.24e-13 |
| past deck-edge immersion the wall-sided formula over-predicts | yes | deck edge at 21.80 deg |
| barge displacement is held to 1e-12 at every station | yes | 9.2e-13 |
| barge curve peaks before it vanishes | yes | 0.1380 m at 40 deg |
| barge vanishing angle is 90 deg | yes | 90.00000000000003 |
| cylinder reproduces the exact arm at every heel | yes | 1.89e-03 at 96 segments |
| cylinder error is second order in the segment size | yes | 3.53e-02 -> 1.89e-03 |
| hull GM from the curve slope matches the waterplane integral | yes | worst 6.50e-04 |
| hull arm is grid-converged at 10 deg | yes | 0.02285 -> 0.02211 m |
| a body at exactly full submergence is refused | yes | displacement volume 7.374345349674991 m^3 is the whole body volume 7.37434534967 |
