# Surface--surface intersection: measured against analytic references

Acceptance (GAP-CAD-002): relative arc-length error < 3%, residual < 1e-09 L.

| study | measured [m] | reference [m] | rel. error | branches | closed |
|---|---|---|---|---|---|
| plane through the cylinder axis (two meridians) | 8.000000000 | 8.000000000 | 0.000e+00 | 2 | 0 |
| inclined plane cuts an ellipse at 30 deg | 6.777877070 | 6.777917613 | 5.982e-06 | 1 | 1 |
| sphere patch cut by a plane | 5.441337392 | 5.441398093 | 1.116e-05 | 1 | 1 |
| cylinder cut by a concentric sphere | 6.283115216 | 6.283185307 | 1.116e-05 | 2 | 2 |
| NURBS rational cylinder cut by a plane | 6.283113706 | 6.283185307 | 1.140e-05 | 1 | 1 |
| perpendicular equal cylinders (Viviani) | 15.280707677 | 15.280791156 | 5.463e-06 | 2 | 2 |
| plane patch cut by a cylinder (truncated) | 1.178089830 | 1.178097245 | 6.294e-06 | 1 | 0 |

## Convergence (sphere cap vs plane, refine_levels=2)

| u samples | length [m] | rel. error | order to next row |
|---|---|---|---|
| 24 | 5.440426927 | 1.785e-04 | 2.00 |
| 48 | 5.441155292 | 4.462e-05 | 2.00 |
| 96 | 5.441337392 | 1.116e-05 | 2.00 |
| 192 | 5.441382917 | 2.789e-06 | - |

## Residuals on both surfaces (Viviani)

- reported max residual: 6.661e-16 m
- recomputed independently: 6.661e-16 m over 338 points
- max closure defect: 0.000e+00 m
