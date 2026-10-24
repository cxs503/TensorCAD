# Geometry representation benchmark

Generated from `results.json` by `tensorcad.benchmark`; every number below was
produced by the run that wrote this file.

- command: `python examples/benchmarks/run_benchmark.py --output docs/benchmarks/full`
- tensorcad 0.9.0, python 3.11.2, numpy 2.4.6, torch 2.14.1+cu130, ocp None
- wall clock: 6.22 s
- **verdict: pass**

## case_a_sphere_dent

*A local edit with a closed-form reference: how much of the error is the edit, and how much is the representation?*

Reference: sphere R = 0.5 m with a Gaussian dent of depth 0.06 m and sigma 0.22 rad on +x; 512 stratified rays, dent < 1 sigma, protected > 2.5 sigma + 1.5 cells

- the mesh row moves only vertices inside 2.5 sigma and does not re-project: the error it shows is interpolation between vertices, not a defect of meshes
- every row is probed along the same stratified ray set: a third of the rays lie inside the edit cone so a small edit is resolved, the rest are uniform, so no row sees an easier sample
- an unweighted mean over rays is not a volume average: the rows are compared with each other, not with a surface integral

### `analytic_field` — reference: exact closed form (measured, 0.06 s)

| metric | value | unit | status | gate |
|---|---|---|---|---|
| `reference_radius_error_max_m` | 0 | m | measured | pass |
| `volume_m3` | 0.519374 | m^3 | measured | — |
| `surface_area_m2` | 3.12963 | m^2 | measured | — |
| `triangulated_volume_error` | 6.671e-04 | relative | measured | — |
| `triangulated_area_error` | 3.347e-04 | relative | measured | — |
| `watertight` | 1 | bool | measured | pass |

> the surface radius is a formula, so this row's surface error is the harness's own; its useful numbers are the volume/area of the triangulation a solver would receive

### `mesh_local_edit` — control: explicit triangles, vertices moved in place (measured, 2.97 s)

| metric | value | unit | status | gate |
|---|---|---|---|---|
| `dent_depth_measured_mean_m` | 0.047199 | m | measured | — |
| `dent_depth_mean_error_m` | 2.255e-06 | m | measured | pass |
| `dent_depth_axis_error_m` | 1.276e-04 | m | measured | pass |
| `dent_depth_peak_m` | 0.0597846 | m | measured | — |
| `dent_surface_error_mean_m` | 4.546e-05 | m | measured | pass |
| `dent_surface_error_max_m` | 1.498e-04 | m | measured | — |
| `protected_displacement_mean_m` | 0 | m | measured | pass |
| `protected_displacement_max_m` | 0 | m | measured | — |
| `protected_total_error_mean_m` | 1.264e-04 | m | measured | — |
| `surface_error_mean_m` | 1.024e-04 | m | measured | — |
| `surface_error_max_m` | 0.00122204 | m | measured | — |
| `transition_band_displacement_m` | 0.00195981 | m | measured | — |
| `locality_ratio` | 0 | relative | measured | pass |
| `volume_error` | 3.833e-04 | relative | measured | pass |
| `discretisation_floor_m` | 1.904e-04 | m | measured | — |
| `ray_coverage` | 1 | fraction | measured | pass |
| `watertight` | 1 | bool | measured | pass |
| `euler_characteristic` | 2 | count | measured | pass |

> 2241 of 12642 vertices moved; every other vertex is bit-identical, so the protected-region leak is exactly zero wherever no moved vertex is reachable, and the edit region pays for linear interpolation instead
> no re-projection, no fairing, no smoothing: this is what a plain mesh edit costs
> protected_displacement_* is measured against the unedited mesh on the same ray; protected_total_error_* is measured against the analytic surface and includes the chord sag of the triangulation

### `neural_field` — learned representation: versioned .ngo object, queried as an implicit field (not_supported, 0.00 s)

| metric | value | unit | status | gate |
|---|---|---|---|---|
| `dent_surface_error_mean_m` | — | m | not_supported | — |
| `protected_displacement_mean_m` | — | m | not_supported | — |
| `protected_total_error_mean_m` | — | m | not_supported | — |
| `locality_ratio` | — | relative | not_supported | — |

> this row is deliberately not simulated: a learned field's numbers must come from a trained artifact, and an untrained stand-in would be a mock result

### `ocp_step` — kernel path: faceted solid exported as STEP (not_supported, 0.00 s)

| metric | value | unit | status | gate |
|---|---|---|---|---|
| `step_volume_error` | — | relative | not_supported | — |

> a faceted solid is not a feature B-rep; when this row runs, it reports the facet count next to the round-trip volume so the distinction is measurable

#### gates

| gate | result | detail |
|---|---|---|
| `analytic_row_is_the_reference` | pass | the harness reproduces the closed form to 1e-12 m |
| `mesh_edit_hits_the_target_depth` | pass | axis depth error 1.276e-04 m (mean over the cone 0.047199 m) |
| `mesh_edit_is_local` | pass | locality ratio 0.000e+00 |
| `every_ray_hits_every_surface` | pass | 1.0000 coverage |
| `the_depth_threshold_is_the_same_for_every_measured_row` | pass | 1 % of the edit depth, applied unchanged to every representation that measured it (the analytic row is the reference and is not a candidate) |
| `no_representation_went_unmeasured_silently` | pass | every row is measured, imported or explicitly not_supported |

## case_b_hull_handoff

*The engineering case: can a target become a verified hull, and does the artifact a solver reads describe the same surface?*

Reference: L = 12.0 m, B = 2.4 m, T = 0.9 m; analytic volume 12.170571118 m^3

- the mesh volume here is a discretisation of an analytic form, so its error is the grid's: the refinement row measures how much

### `parametric_hull` — engineering case: three design variables, analytic volume reference (measured, 1.78 s)

| metric | value | unit | status | gate |
|---|---|---|---|---|
| `volume_error` | 0.00880431 | relative | measured | pass |
| `design_loop_volume_error` | 2.220e-16 | relative | measured | pass |
| `design_loop_solved` | 1 | bool | measured | pass |
| `design_loop_evaluations` | 6 | count | measured | — |
| `design_loop_seconds` | 1.76485 | s | measured | — |
| `watertight` | 1 | bool | measured | pass |
| `euler_characteristic` | 2 | count | measured | pass |
| `export_node_error_m` | 0 | m | measured | pass |
| `export_area_error` | 6.833e-16 | relative | measured | pass |
| `surface_area_m2` | 51.9923 | m^2 | measured | — |
| `grid_refinement_volume_change` | 0.0057745 | relative | measured | pass |

> the design loop converged from 85 % of the target dimensions in 6 evaluations (converged)
> analytic volume 12.170571 m^3 vs mesh 12.063418 m^3 at (41, 13); refining to (81, 25) moves it by 0.577 %

### `ocp_step` — kernel path: faceted solid exported as STEP (not_supported, 0.00 s)

| metric | value | unit | status | gate |
|---|---|---|---|---|
| `step_volume_error` | — | relative | not_supported | — |

#### gates

| gate | result | detail |
|---|---|---|
| `hull_mesh_matches_the_analytic_volume` | pass | 0.880 % error at (41, 13) |
| `design_loop_reaches_the_target` | pass | converged, 6 evaluations |
| `export_round_trip_is_exact` | pass | worst node displacement 0.0e+00 m, area 6.8e-16 |

