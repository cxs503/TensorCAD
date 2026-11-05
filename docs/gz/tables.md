# Righting-arm tables

Rendered from `docs/gz/results.json` by `examples/gz/publish_evidence.py`; the JSON
is written by `examples/gz/run_example.py`. Every number below is a formatting of a
value in that file -- `tests/test_gz_evidence.py` re-renders this page and compares
it byte for byte, so a hand-typed figure cannot survive CI.

## 1. Cases

| body | geometry | reference used |
|---|---|---|
| box barge | L 10.0 m, B 2.0 m, D 2.0 m, T 0.8 m, KG 0.5 m | wall-sided trapezoid, valid to 38.66 deg |
| circular log (regular 64-gon) | R 1.0 m, half submerged, KG 0.5 m | exact `(R - KG) sin(phi)`, every heel |
| parametric hull, 41x13 mesh | L 10.0 m, B 2.0 m, T 0.8 m, freeboard 1.0 m, KG 0.50 m | `GM` from `tensorcad.stability.initial_stability` (waterplane moments) |

## 2. Box barge against the wall-sided closed form (0-38 deg)

`GM0 = 0.3167 m`, `BM0 = 0.4167 m`.

| heel [deg] | GZ integrated [m] | wall-sided [m] | abs. error [m] |
|---:|---:|---:|---:|
| 0 | 0.000000000 | 0.000000000 | 2.60e-18 |
| 5 | 0.027738300 | 0.027738300 | 2.08e-17 |
| 10 | 0.056113367 | 0.056113367 | 3.47e-17 |
| 15 | 0.085830692 | 0.085830692 | 6.94e-17 |
| 20 | 0.117745731 | 0.117745731 | 5.55e-17 |
| 25 | 0.152973973 | 0.152973973 | 5.55e-17 |
| 30 | 0.193055556 | 0.193055556 | 5.55e-17 |
| 35 | 0.240219858 | 0.240219858 | 2.78e-17 |

Worst absolute error inside the validity window: **6.94e-17 m**
(round-off: the integrator and the closed form are both exact for this polyhedron).

## 3. Negative control: past the wall-sided limit

The formula is invalid above the limit, so the agreement above must *stop*; this
table is the proof that it can fail.

| heel [deg] | GZ integrated [m] | wall-sided [m] | abs. error [m] |
|---:|---:|---:|---:|
| 45 | 0.353553391 | 0.371231060 | 1.77e-02 |
| 55 | 0.456144383 | 0.607470652 | 1.51e-01 |
| 70 | 0.523110651 | 1.775359736 | 1.25e+00 |

## 4. Circular log against the exact `(R - KG) sin(phi)`

| heel [deg] | GZ integrated [m] | exact [m] | abs. error [m] |
|---:|---:|---:|---:|
| 0.0000 | 0.000000000 | 0.000000000 | 1.81e-17 |
| 2.8125 | 0.024533837 | 0.024533837 | 1.39e-17 |
| 5.6250 | 0.049008570 | 0.049008570 | 1.39e-17 |
| 8.4375 | 0.073365237 | 0.073365237 | 2.78e-17 |
| 11.2500 | 0.097545161 | 0.097545161 | 5.55e-17 |
| 14.0625 | 0.121490090 | 0.121490090 | 0.00e+00 |
| 16.8750 | 0.145142339 | 0.145142339 | 5.55e-17 |
| 19.6875 | 0.168444927 | 0.168444927 | 0.00e+00 |
| 22.5000 | 0.191341716 | 0.191341716 | 0.00e+00 |
| 25.3125 | 0.213777547 | 0.213777547 | 5.55e-17 |
| 28.1250 | 0.235698368 | 0.235698368 | 1.11e-16 |
| 30.9375 | 0.257051372 | 0.257051372 | 0.00e+00 |
| 33.7500 | 0.277785117 | 0.277785117 | 0.00e+00 |
| 36.5625 | 0.297849652 | 0.297849652 | 0.00e+00 |
| 39.3750 | 0.317196642 | 0.317196642 | 1.11e-16 |
| 42.1875 | 0.335779477 | 0.335779477 | 0.00e+00 |
| 45.0000 | 0.353553391 | 0.353553391 | 1.11e-16 |
| 47.8125 | 0.370475563 | 0.370475563 | 0.00e+00 |
| 50.6250 | 0.386505227 | 0.386505227 | 0.00e+00 |
| 53.4375 | 0.401603766 | 0.401603766 | 0.00e+00 |
| 56.2500 | 0.415734806 | 0.415734806 | 0.00e+00 |
| 59.0625 | 0.428864305 | 0.428864305 | 1.11e-16 |
| 61.8750 | 0.440960632 | 0.440960632 | 0.00e+00 |
| 64.6875 | 0.451994647 | 0.451994647 | 0.00e+00 |
| 67.5000 | 0.461939766 | 0.461939766 | 0.00e+00 |
| 70.3125 | 0.470772033 | 0.470772033 | 0.00e+00 |
| 73.1250 | 0.478470168 | 0.478470168 | 0.00e+00 |
| 75.9375 | 0.485015627 | 0.485015627 | 1.11e-16 |
| 78.7500 | 0.490392640 | 0.490392640 | 1.11e-16 |
| 81.5625 | 0.494588255 | 0.494588255 | 0.00e+00 |
| 84.3750 | 0.497592363 | 0.497592363 | 1.11e-16 |
| 87.1875 | 0.499397728 | 0.499397728 | 0.00e+00 |
| 90.0000 | 0.500000000 | 0.500000000 | 0.00e+00 |

Worst absolute error over the whole sweep: **1.11e-16 m**.

## 5. Hull: the curve, and its slope at zero

| quantity | value |
|---|---:|
| `GZ'(0)` from the curve (Richardson) [m/rad] | 0.421622996 |
| `GM` from `initial_stability` (waterplane moments) [m] | 0.421622996 |
| relative difference | 1.77e-09 |
| `KB` (upright) [m] | 0.483005 |
| `BM` transverse (upright) [m] | 0.438618 |
| maximum `GZ` [m] | 0.558859 |
| heel at maximum `GZ` [deg] | 90.0 |
| area under the curve to 30 deg [m·rad] | 0.052100 |
| displacement held at [m³] | 7.446554 |
| worst displacement residual [m³] | 8.59e-14 |

| heel [deg] | GZ [m] | KN [m] | waterline offset [m] | V [m³] | volume residual [m³] |
|---:|---:|---:|---:|---:|---:|
| 0.0 | 0.000000000 | 0.000000 | 0.000000 | 7.446554 | 0.00e+00 |
| 5.0 | 0.035521387 | 0.079099 | -0.003871 | 7.446554 | 0.00e+00 |
| 10.0 | 0.068975985 | 0.155800 | -0.015429 | 7.446554 | 0.00e+00 |
| 15.0 | 0.100967261 | 0.230377 | -0.034589 | 7.446554 | 0.00e+00 |
| 20.0 | 0.132034360 | 0.303044 | -0.061243 | 7.446554 | 2.00e-15 |
| 25.0 | 0.162714694 | 0.374024 | -0.095265 | 7.446554 | 8.60e-14 |
| 30.0 | 0.193605846 | 0.443606 | -0.136519 | 7.446554 | 0.00e+00 |
| 35.0 | 0.225437698 | 0.512226 | -0.184869 | 7.446554 | 0.00e+00 |
| 40.0 | 0.259171095 | 0.580565 | -0.240201 | 7.446554 | -1.00e-15 |
| 45.0 | 0.296144005 | 0.649697 | -0.302449 | 7.446554 | 0.00e+00 |
| 50.0 | 0.338035171 | 0.721057 | -0.371536 | 7.446554 | 0.00e+00 |
| 55.0 | 0.380909342 | 0.790485 | -0.445080 | 7.446554 | 0.00e+00 |
| 60.0 | 0.420731907 | 0.853745 | -0.520761 | 7.446554 | 0.00e+00 |
| 65.0 | 0.456602750 | 0.909757 | -0.597711 | 7.446554 | 1.00e-15 |
| 70.0 | 0.487984230 | 0.957831 | -0.675547 | 7.446554 | 0.00e+00 |
| 75.0 | 0.515140832 | 0.998104 | -0.753765 | 7.446554 | 0.00e+00 |
| 80.0 | 0.537894651 | 1.030299 | -0.832217 | 7.446554 | 0.00e+00 |
| 85.0 | 0.552362122 | 1.050459 | -0.912653 | 7.446554 | 0.00e+00 |
| 90.0 | 0.558859016 | 1.058859 | -0.994745 | 7.446554 | 0.00e+00 |

## 6. Free surface and upright hydrostatics

| quantity | value |
|---|---:|
| tank | 2.0 m x 0.5 m rectangle |
| `i_x` [m⁴] | 0.020833 |
| flat loss on every point [m] | 0.002798 |
| peak `GZ`, solid [m] | 0.558859 |
| peak `GZ`, with free surface [m] | 0.556061 |
| upright displacement [m³] | 7.446554 |
| upright `KB` [m] | 0.483005 |
| waterplane area [m²] | 14.061669 |
| wetted surface area [m²] | 22.838010 |

## 7. Claim ledger

- verdict: **all claims hold**
- box barge: worst 6.94e-17 m < 1e-12 m (pass)
- log: worst 1.11e-16 m < 1e-9 m (pass)
- hull `GZ'(0)` vs `GM`: 1.77e-09 relative < 1e-6 (pass)
- displacement held at every point: worst 8.59e-14 m³ < 1e-10 (pass)
- negative control: wall-sided formula off by 1.25e+00 m at 70 deg > 1e-2 m (pass)
