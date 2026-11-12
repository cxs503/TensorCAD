| fixture | what it checks | status | DOF | rank / rows | iterations | max residual [m] | max position error [m] |
|---|---|---|---:|---:|---:|---:|---:|
| `fixed_point` | one point pinned at (2.5, -1.25) | fully_constrained | 0 | 2 / 2 | 2 | 0.0e+00 | 0.0e+00 |
| `distance_circle` | point held 5 m from a pinned point (one equation, two parameters) | under_constrained | 1 | 3 / 3 | 3 | 8.9e-16 | 0.0e+00 |
| `horizontal_arm` | horizontal 3 m arm off a pinned point | fully_constrained | 0 | 4 / 4 | 3 | 1.8e-17 | 1.8e-17 |
| `sss_triangle` | triangle from three side lengths (4, 3, 2 m), apex from the law of cosines | fully_constrained | 0 | 6 / 6 | 5 | 4.4e-16 | 4.4e-16 |
| `rectangle` | 4 x 3 rectangle from horizontal/vertical pairs and two dimensions | fully_constrained | 0 | 8 / 8 | 3 | 3.8e-18 | 7.1e-18 |
| `angle_arm` | 2 m arm at 30 degrees to a pinned x reference line | fully_constrained | 0 | 6 / 6 | 5 | 4.4e-16 | 2.5e-16 |
| `perpendicular_chain` | vertical 3 m arm followed by a 4 m arm perpendicular to it | fully_constrained | 0 | 6 / 6 | 4 | 2.7e-18 | 3.9e-18 |
| `point_on_line` | point on the 3-4-5 line at 2.5 m from its start | fully_constrained | 0 | 6 / 6 | 3 | 1.8e-16 | 2.2e-16 |
| `point_line_offset` | point 1 m to the left of the 3-4-5 line, x anchored | fully_constrained | 0 | 6 / 6 | 2 | 7.8e-16 | 2.6e-15 |
| `circle_between_lines` | circle tangent to two parallel lines 2 m apart (r = h / 2) | fully_constrained | 0 | 11 / 11 | 2 | 1.1e-16 | 1.1e-16 |
| `point_on_circle` | point on a radius-2 circle with its ordinate anchored | fully_constrained | 0 | 5 / 5 | 3 | 2.0e-13 | 2.4e-13 |
| `under_constrained_line` | a bare horizontal line: 4 parameters, one equation | under_constrained | 3 | 1 / 1 | 2 | 2.8e-17 | 0.0e+00 |
| `inconsistent` | point pinned at (3, 0) and asked to be 5 m from the origin | inconsistent | 0 | 4 / 5 | 4 | 6.7e-01 | 0.0e+00 |
| `redundant` | the same 5 m dimension stated twice over two pinned points | fully_constrained | 0 | 4 / 6 | 3 | 8.0e-17 | 1.0e-16 |

Perturbed-restart study: **325 / 325** solves converged (25 deterministic restarts per solvable fixture, perturbation ±0.25 m), median 3 iterations, worst position drift 3.9e-12 m.

| drawing scale | DOF | rank / rows | max residual [m] | relative position error |
|---:|---:|---:|---:|---:|
| 0.001 | 0 | 8 / 8 | 3.8e-21 | 1.8e-18 |
| 1 | 0 | 8 / 8 | 3.8e-18 | 1.8e-18 |
| 1000 | 0 | 8 / 8 | 3.8e-15 | 1.8e-18 |
