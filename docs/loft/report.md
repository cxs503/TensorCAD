# Exact section skinning -- measured results

| claim | holds |
|---|---|
| `cylinder_radial_deviation_lt_1e-13` | yes |
| `cylinder_area_lt_1e-12` | yes |
| `cone_deviation_lt_1e-12` | yes |
| `cone_area_lt_1e-9` | yes |
| `parabola_deviation_lt_1e-13` | yes |
| `cubic_skin_is_c2` | yes |
| `cubic_jump_lt_1e-12` | yes |
| `quadratic_skin_jumps` | yes |
| `ruled_skin_jumps_in_first_derivative` | yes |
| `mixed_degree_defect_lt_1e-13` | yes |
| `mixed_degree_transverse_degree_2` | yes |
| `rotating_sections_defect_lt_1e-13` | yes |

Cylinder: max radial deviation 3.331e-16 m, area rel. error 2.120e-15.
Cone: max deviation 4.441e-16 m, area rel. error 2.426e-15.
Parabola family: max |z - y^2| = 5.551e-16 m.
Continuity (relative, one-sided at interior knots): degree_u=1: order 0, degree_u=2: order 1, degree_u=3: order 2

