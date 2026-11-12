# patch_network: numbers injected from `outputs/results.json`

## Box 2 x 3 x 4 m (6 bilinear patches, exact tessellation)

| quantity | measured | analytic | relative error |
| --- | --- | --- | --- |
| volume [m3] | 24 | 24 | +0.000e+00 |
| area [m2] | 52 | 52 | +0.000e+00 |
| watertight | True | - | - |
| seam gap max [m] | 0.000e+00 | 0 | - |

## Torus R = 2 m, r = 0.5 m (16 patches, 32 seams)

| samples | triangles | volume [m3] | volume err | area [m2] | area err | order (vol) |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 2048 | 9.742651548 | -1.286e-02 | 39.255706568 | -5.641e-03 | - |
| 16 | 8192 | 9.837731861 | -3.229e-03 | 39.422606552 | -1.414e-03 | 1.99 |
| 32 | 32768 | 9.861627839 | -8.082e-04 | 39.464456486 | -3.536e-04 | 2.00 |
| 64 | 131072 | 9.867609733 | -2.021e-04 | 39.474926802 | -8.842e-05 | 2.00 |

analytic: V = 2 pi^2 R r^2 = 9.869604401 m3, A = 4 pi^2 R r = 39.478417604 m2
seam gap max = 6.123e-16 m, seam dihedral max = 1.479e-06 deg

## Double cone R = 1 m, H = 2 x 1 m (8 patches, 2 apex poles)

| samples | triangles | volume [m3] | volume err | order (vol) |
| --- | --- | --- | --- | --- |
| 8 | 960 | 2.080881389 | -6.452e-03 | - |
| 16 | 3968 | 2.091010586 | -1.616e-03 | 2.00 |
| 32 | 16128 | 2.093548592 | -4.042e-04 | 2.00 |
| 64 | 65024 | 2.094183451 | -1.011e-04 | 2.00 |

analytic: V = 2 pi R^2 H / 3 = 2.094395102 m3, A = 2 pi R sqrt(R^2 + H^2) = 8.885765876 m2
pole sides collapsed = 8, open sides = 0

## Negative controls

```json
{
  "moved_patch": {
    "offset_m": 0.001,
    "sides_left_unclaimed_by_derive": 8,
    "watertight_with_declared_seams": true,
    "measured_gap_m": 0.0009999999999998899,
    "volume_error": 0.00028124999999927053
  },
  "missing_patch": {
    "open_sides": 4,
    "watertight": false,
    "boundary_edges": 16
  },
  "unclaimed_but_coincident_seam": {
    "open_sides": 2,
    "watertight": true
  },
  "reversed_flag": {
    "gap_declared_correctly_m": 0.0,
    "gap_flag_omitted_m": 1.0,
    "edge_length_m": 1.0
  },
  "sample_mismatch": {
    "refused": true,
    "message": "seam 0 joins a side sampled with 8 nodes to one sampled with 4; a u1/v0 seam needs nu == nv (got 4 x 8)"
  }
}
```
