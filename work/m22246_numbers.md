| body | groups | verts | tets | boundary tri | V file [m³] | V tets [m³] | rel | watertight |
|---|---|---|---|---|---|---|---|---|
| box 2x1x0.6 | bottom top fore aft hull | 140 | 432 | 216 | 1.200000000000 | 1.200000000000 | 0 | yes |
| Wigley body L10xB1xT0.625 | bottom hull | 1696 | 8280 | 1134 | 2.762123390480 | 2.762123390480 | 0 | yes |

Readers of the same file (each an independent implementation):

| body | this module | meshio | FEM handoff (`solver_interface`) | cells seen by handoff |
|---|---|---|---|---|
| box 2x1x0.6 | 0 | 0 | 0 | 432 tet / 216 tri |
| Wigley body L10xB1xT0.625 | 0 | 0 | 0 | 8280 tet / 1134 tri |

Physical-group areas against the closed-form / quadrature anchors:

| body | patch | area [m²] | anchor [m²] | rel |
|---|---|---|---|---|
| box 2x1x0.6 | `bottom` | 2.000000000 | 2.000000000 | 8.882e-16 |
| box 2x1x0.6 | `top` | 2.000000000 | 2.000000000 | 8.882e-16 |
| box 2x1x0.6 | `fore` | 0.600000000 | 0.600000000 | 1.850e-16 |
| box 2x1x0.6 | `aft` | 0.600000000 | 0.600000000 | 1.850e-16 |
| box 2x1x0.6 | `hull` | 2.400000000 | 2.400000000 | 1.850e-16 |
| Wigley body L10xB1xT0.625 | `bottom` | 6.655092593 | 6.666666667 | 1.736e-03 |
| Wigley body L10xB1xT0.625 | `hull` | 14.839761791 | 14.879052246 | 2.641e-03 |

`all_pass` of the run that produced these numbers: **True** (criteria in the JSON: `gmsh_volume_rel` < 1e-12, `gmsh_area_rel` < 0.01).
