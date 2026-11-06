### Area quadrature against the closed forms

| curve | interval [deg] | computed [m rad] | closed form [m rad] | rel. error |
|---|---|---|---|---|
| GZ = 0.5 sin(phi) | 0-30 | 0.066987298108 | 0.066987298108 | 4.14e-16 |
| GZ = 0.5 sin(phi) | 0-40 | 0.116977778441 | 0.116977778441 | 1.19e-16 |
| GZ = 0.5 sin(phi) | 30-40 | 0.049990480333 | 0.049990480333 | 6.94e-16 |
| wall-sided barge (B=6 m, T=3 m, KG=2 m) | 0-30 | 0.077350269190 | 0.077350269190 | 1.79e-16 |
| wall-sided barge (B=6 m, T=3 m, KG=2 m) | 0-40 | 0.152703644666 | 0.152703644666 | 1.82e-16 |
| wall-sided barge (B=6 m, T=3 m, KG=2 m) | 30-40 | 0.075353375477 | 0.075353375477 | 1.84e-16 |

### GM0 read from the slope at upright

| curve | computed GM0 [m] | reference [m] | rel. error |
|---|---|---|---|
| GZ = 0.5 sin(phi) | 0.50000000 | 0.50000000 | 0.00e+00 |
| wall-sided barge | 0.49999988 | 0.50000000 | 2.48e-07 |

### Table coarsening: error and reported uncertainty

| step [deg] | samples | area 0-30 [m rad] | error [m rad] | reported uncertainty [m rad] | covered |
|---|---|---|---|---|---|
| 20 | 4 | 0.077589369 | +2.39e-04 | 2.24e-03 | True |
| 10 | 6 | 0.077356601 | +6.33e-06 | 1.16e-03 | True |
| 5 | 11 | 0.077350634 | +3.65e-07 | 2.92e-04 | True |
| 2.5 | 21 | 0.077350292 | +2.23e-08 | 7.33e-05 | True |
| 1 | 51 | 0.077350270 | +5.68e-10 | 1.17e-05 | True |

### Verdicts

* **compliant analytic (GZ = 0.5 sin phi)** -> `pass` (2.2.1a=pass, 2.2.1b=pass, 2.2.1c=pass, 2.2.2=pass, 2.2.3=pass, 2.2.4=pass)
  areas 0.06698729810778067 / 0.11697777844051098 / 0.04999048033273031 m rad, GZ(30) = 0.24999999999999997, GM0 = 0.5
* **wall-sided barge, analytic curve** -> `pass` (2.2.1a=pass, 2.2.1b=pass, 2.2.1c=pass, 2.2.2=pass, 2.2.3=pass, 2.2.4=pass)
  areas 0.07735026918962577 / 0.1527036446661393 / 0.07535337547651354 m rad, GZ(30) = 0.33333333333333326, GM0 = 0.49999987617078295
* **wall-sided barge, mesh -> righting_arm_curve** -> `pass` (2.2.1a=pass, 2.2.1b=pass, 2.2.1c=pass, 2.2.2=pass, 2.2.3=pass, 2.2.4=pass)
  areas 0.07735660053515858 / 0.15358214813019114 / 0.07622554759503256 m rad, GZ(30) = 0.33333333333333337, GM0 = 0.4986484141785438
* **tender control (GZ = 0.1 sin phi)** -> `fail` (2.2.1a=fail, 2.2.1b=fail, 2.2.1c=fail, 2.2.2=fail, 2.2.3=pass, 2.2.4=fail)
  areas 0.013397459621556133 / 0.023395555688102195 / 0.009998096066546062 m rad, GZ(30) = 0.049999999999999996, GM0 = 0.1

### Cross-check against `tensorcad.intact_stability`

| quantity | `tensorcad.statutory` | `tensorcad.intact_stability` | closed form |
|---|---|---|---|
| area 0-30 deg [m rad] | 0.077350269192 | 0.077351002569 | 0.077350269190 |
| area 0-40 deg [m rad] | 0.152703644670 | 0.152705345159 | 0.152703644666 |
| verdict, compliant 0.5 sin(phi) | `pass` (failed none) | `passes IS Code 2008 part A 2.2.1-2.2.6` (failed none) | - |
| verdict, tender 0.1 sin(phi) | `fail` (failed 2.2.1a, 2.2.1b, 2.2.1c, 2.2.2, 2.2.4) | `fails IS Code 2008 part A 2.2.1, 2.2.2, 2.2.3, 2.2.4, 2.2.6` (failed 2.2.1, 2.2.2, 2.2.3, 2.2.4, 2.2.6) | - |

### Mesh route: GZ at the sampled heels against the closed form

| heel [deg] | mesh GZ [m] | closed form [m] | |error| [m] |
|---|---|---|---|
| 0 | 0.000000000000000 | 0.000000000000000 | 0.00e+00 |
| 10 | 0.089523554302420 | 0.089523554302420 | 8.33e-17 |
| 20 | 0.193664516574472 | 0.193664516574472 | 2.78e-17 |
| 30 | 0.333333333333333 | 0.333333333333333 | 1.11e-16 |
| 40 | 0.547683387507425 | 0.547683387507424 | 8.88e-16 |
