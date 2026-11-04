# Large-angle stability: GZ from the heeled displacement

Barge 10.0 x 2.0 x 1.2 m, draft 0.8 m (`KG = 0.78 m`); deep barge depth 2.4 m (`KG = 0.3` / `0.9 m`); sphere `R = 1 m`, `KG = 0.5 m`.

## Numbers

| item | value |
|---|---:|
| deck_edge_deg | 21.801409486351808 |
| wall_sided_limit_deg | 38.659808254090095 |
| angle_of_maximum_deg | 30.0 |
| maximum_arm_m | 0.03282527428325002 |
| vanishing_angle_deg | 49.45212071731248 |
| dynamical_stability_m_rad | 0.015181390916075203 |
| loll_angle_deg | 32.3115332374232 |
| expected_loll_deg | 32.31153323742385 |

## Claims

| claim | measured | limit | holds |
|---|---:|---|:--:|
| barge: mesh arm vs independent 2-D section integration | 1.4405e-14 m | < 1e-12 m | yes |
| deep barge: mesh arm vs wall-sided closed form | 6.1448e-16  | < 1e-12 relative | yes |
| heeled displacement balance (all angles, both barges) | 6.8390e-14  | < 1e-13 relative | yes |
| GM from the curve at 0.25 deg, minus the BM/2*tan^2 term | 7.9936e-15 m | < 1e-12 m | yes |
| halving the heel quarters the GM deficit (second order) | 4.0002e+00  | > 3.9 | yes |
| angle of loll vs tan^2(phi) = -2 GM / BM | 6.5370e-13 deg | < 1e-6 deg | yes |
| GZ at the bisected angle of vanishing stability | 3.6082e-16 m | < 1e-12 m | yes |
| dynamical stability 0-30 deg vs closed-form integral | 1.3166e-05  | < 1e-3 relative | yes |
| sphere (20^3) worst GZ error vs (KM - KG) sin(phi) | 2.0404e-03  | < 5e-3 | yes |
| sphere: error ratio between 12^3 and 20^3 | 3.0276e+00  | > 2 | yes |
| dGZ/dKG = -sin(phi) | 2.7756e-17 m | < 1e-13 m | yes |
| a mesh capped at the design waterline is refused | 1.0000e+00  | == 1 | yes |

## Curve, sampled

| heel [deg] | mesh GZ [m] | section reference [m] | difference [m] |
|---:|---:|---:|---:|
| 0 | -0.000000000 | 0.000000000 | 9.25e-18 |
| 5 | 0.003334692 | 0.003334692 | 2.08e-17 |
| 10 | 0.007491877 | 0.007491877 | 0.00e+00 |
| 15 | 0.013361359 | 0.013361359 | 0.00e+00 |
| 20 | 0.021980091 | 0.021980091 | 1.39e-17 |
| 25 | 0.031337396 | 0.031337396 | 1.11e-16 |
| 30 | 0.032825274 | 0.032825274 | 1.44e-14 |
| 35 | 0.029067380 | 0.029067380 | 2.78e-17 |
| 40 | 0.022314481 | 0.022314481 | 8.33e-17 |
| 45 | 0.012963624 | 0.012963624 | 5.55e-17 |
| 50 | -0.001842411 | -0.001842411 | 2.78e-17 |
| 55 | -0.020531963 | -0.020531963 | 2.78e-17 |
| 60 | -0.041717906 | -0.041717906 | 5.55e-17 |
| 65 | -0.064433580 | -0.064433580 | 2.78e-17 |
| 70 | -0.087958670 | -0.087958670 | 2.78e-17 |
| 75 | -0.111723554 | -0.111723554 | 8.33e-17 |
| 80 | -0.135253958 | -0.135253958 | 1.11e-16 |
| 85 | -0.158137168 | -0.158137168 | 5.55e-17 |
| 90 | -0.180000000 | -0.180000000 | 8.33e-17 |
