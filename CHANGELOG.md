# Changelog

All notable changes to TensorCAD are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
semantic versioning.

## 0.21.0 — 2026-11-04

### Added
- **Waterplane second moments on the `Hydrostatics` record** (`tensorcad/hydrostatics.py`):
  `compute_hydrostatics` now reports `transverse_waterplane_inertia` (`∫ y² dA`),
  `longitudinal_waterplane_inertia` (`∫ (x−LCF)² dA`) and `keel_z` next to the waterplane
  area. They come out of the flux identity the module already uses for the waterplane area
  (`∫_cap f dA = −∫_wetted f n_z dA` with `f` independent of `z`), so no cap has to exist
  and no waterline polygon has to be reconstructed; the per-face second moments are exact
  (`x`, `y` linear over a triangle). An open or capless surface reports `nan` rather than a
  number the shell cannot support, and `to_dict` carries all three quantities.
- `tensorcad.hull.analytic_hydrostatics` reports the same two inertias from the closed forms
  (`hydrostatics_analytic.waterplane_inertia_transverse` — gamma identity — and
  `waterplane_inertia_longitudinal` — tanh-sinh quadrature with its own convergence
  evidence), so the mesh record and the analytic reference speak the same quantities.
- `tests/test_hydrostatics_inertia.py` — the record against three independent references: a
  vertex-by-vertex rectangular barge (`L·B³/12`, `B·L³/12`, `keel_z`) exact to `1e-12` on
  three tessellations, a circular waterplane converging to `πR⁴/4` from below at the
  `O(1/n²)` rate, and the hull mesh against the closed forms at three sizes.
- **Transom-stern hull family** (`tensorcad/transom.py`): the aft body is truncated at a
  station `tau` and closed by a planar transom plate, so the parametric hull can carry
  the stern real displacement vessels have instead of a knife edge. Two parameters —
  `transom_fraction` (`tau`) and `breadth_exponent` (`p`, the sine-taper exponent) —
  turn it into a design space: `tau = 0` reproduces the pointed hull exactly (identical
  faces, LOA-based frame with the transom at `x = 0`), and
  `transom_fraction_for_area_ratio` inverts `A_t/A_m = sin(pi*tau)^p` in closed form so a
  loop can dial an immersed transom area directly. Topology is unchanged by the cut:
  watertight, edge-manifold, `V-E+F = 2`, outward, no zero-area faces for every `tau`
  tested. See `docs/transom_stern.md`.
- **Hydrostatic references that survive the truncation.** The incomplete taper integrals
  `int_tau^1 sin(pi*u)^p du` are evaluated from a `sinc^p` power series (the `p`-singular
  factor integrates exactly, the even remainder is a series in `s^2`), giving
  `V`, `A_wp`, `LCB`, `Cb`, `Cp`, `Cwp` and the immersed transom area in closed form.
  Validated three ways: against `Gamma((1+p)/2)/(sqrt(pi)*Gamma(1+p/2))` at `tau = 0`
  (2.2e-16), against the `M(0)/J(0) = 1/2` symmetry (1e-14), and against
  `tanh_sinh_quadrature` for `tau > 0` (1.1e-15). Mesh integration then lands within
  0.045–0.88 % of the reference depending on grid, monotonically converging, and the
  transom plate area converges to its analytic section (−0.71 % at `nz = 13`, −0.043 %
  at `nz = 81`).
- `examples/transom_stern/run_example.py`: reference ladder, 4 truncations × 4 grids
  convergence table, design table and three figures.
- **Large-angle stability: the righting-arm curve** (`tensorcad/gz.py`). `tensorcad.stability`
  stops at `GM`, the slope of `GZ(phi)` at the upright; this module is the integration it
  pointed at. A hull heeled by `phi` floats on the half-space `p . (0, -sin phi, cos phi) <= d`;
  the offset `d` is solved by bisection so the **displacement is held** (not the draft), the
  submerged volume and centroid come from an exact clip of the faceted surface (apex of the
  tetrahedral fan in the cutting plane, so the cap is never built), and the arm is the
  transverse lever `GZ = (B - G) . (0, cos phi, sin phi)` about a `KG` measured from the keel.
- **Two independent routes.** `righting_arm` / `righting_arm_curve` integrate the 3-D heeled
  displacement; `cross_curve_righting_arm` / `cross_curve_righting_arm_curve` integrate the
  classical cross curves from station cuts (`bonjean.section_loops` + inclined waterline). They
  share only the mesh, and they agree to 1.2e-3 relative on the parametric hull.
- **Verified, not asserted** (`docs/gz.md`, `docs/gz/report.md`, 7 self-checking claims in
  `examples/gz/run_example.py`, 45 cases in `tests/test_gz.py`): a wall-sided box reproduces
  `sin(phi)(GM + BM tan²(phi)/2)` to 2.3e-15 relative inside the closed form's validity window
  (phi < 38.66 deg for that box); the curve's slope at the origin equals `KB + BM - KG` to
  2.3e-5; `GZ(-phi) = -GZ(phi)` to 8.6e-16 m; the heeled displacement stays on target to
  1.1e-15; a KG ladder lowers `GZ_max` monotonically and goes negative past `KM`.
- **Summary numbers on the curve**: `maximum_righting_arm`, `angle_of_maximum`,
  `range_of_stability` (interpolated crossing) and `area_under_curve(limit)` — the dynamical
  stability `integral GZ dphi` that statutory criteria consume. No free-surface correction, no
  wind arm, no damage case and no rule pass/fail: the curve is the input to those, not a
  verdict. Also exported: `submerged_properties`, `solve_waterline`, `keel_height`,
  `station_offsets`, `wall_sided_righting_arm`, `HeelPoint`, `RightingArmCurve`, `GZError`.

## 0.20.0 — 2026-11-03

- **A second hull family, and the one benchmarks use** (`tensorcad/wigley.py`): the
  parabolic Wigley `y = (B/2)(1-|2x/L-1|^n)(1-(w/T)^m)` lofted with the same
  ring-per-station topology as `tensorcad.hull`, so it is watertight, outward, Euler-2 and
  capped at the design waterline without a welding pass. `generate_wigley`,
  `wigley_analytic_hydrostatics`, `wigley_analytic_relative_errors`,
  `wigley_convergence_study` and `export_wigley` are exported; the shape exponents `n, m`
  are exposed because the closed forms hold for any `n, m >= 1`.
- **Exact hydrostatics, not a fitted reference**: `Cb = nm/((n+1)(m+1))` (= `4/9`),
  `Cp = Cwp = n/(n+1)` (= `2/3`), `Cm = m/(m+1)` (= `2/3`), `VCB = T - T(m+1)/(2(m+2))`
  (= `5T/8` above the keel), `V = Cb L B T`; only the wetted area needs quadrature, and it
  is a composite Gauss-Legendre rule on the analytic flank that reports its own
  convergence (4e-16 at `n = m = 2`).
- **Verification debt closed for the family**: on the ITTC reference proportions
  (`L/B = 10`, `B/T = 1.6`) the default 41×13 grid integrates volume 0.236 % low,
  waterplane 0.063 %, wetted area 0.020 %, and 161×41 comes down to 0.0195 %; observed
  order 2.17 / 2.58 / 2.01, every signed error negative (inscribed chords). The same 41×13
  grid is 0.88 % on the `sin^0.7 sqrt` family -- the parabolic section has a finite keel
  angle, and that is asserted as a test.
  Evidence: `docs/wigley.md`, `docs/wigley/{results.json,report.md,plot.png}`,
  `examples/wigley/run_example.py` (exits non-zero if any of its 10 claims breaks),
  `tests/test_wigley.py` (30 tests).


## 0.17.0 — 2026-11-03

An asymmetric hull, so that the longitudinal centres stop being an identity. The family
was mirror-symmetric, which made `LCB = LCF = length/2` a consequence of the symmetry: the
check validated the frame and the sign of the `x` moments and could not detect an
antisymmetric bias, because the geometry had none.

### Added
- **`HullParameters.fore_exponent` / `aft_exponent`** — the sine-taper exponent per half,
  both defaulting to the shipped `BREADTH_EXPONENT`, so the default surface is the previous
  one bit for bit (same vertices, faces and hydrostatics). Midship belongs to both halves
  and `sin(pi/2) = 1` for every exponent, so the halves meet at the requested beam without
  a step and the maximum section stays at `x = L/2`.
- **`hull.taper_integral(exponent)`**, **`hull.taper_integral_pair(fore, aft)`** — the
  closed form `Gamma((p+1)/2)/(sqrt(pi) Gamma(p/2+1))` of `int_0^1 sin(pi x)^p dx`; the
  shipped exponent returns `TAPER_INTEGRAL` verbatim so no published digit moves.
- **Analytic references now take the split**: `analytic_volume`,
  `analytic_waterplane_area`, `analytic_longitudinal_centre`, `analytic_wetted_surface_area`
  and `analytic_hydrostatics` accept `fore_exponent`/`aft_exponent`. Volume, waterplane
  area, section area, VCB and both coefficients stay closed form; LCB and LCF are
  `L*(m(p_aft) + h(p_fore) - m(p_fore))/(h(p_aft) + h(p_fore))`.
- **`hull._half_first_moment`** — `int_0^(1/2) x sin(pi x)^p dx`, the one quantity of this
  family without an elementary form, by Gauss-Legendre after `x = t^3/2`;
  **`hull._half_first_moment_exp_sinh`** recomputes it with an exponential substitution and
  a different rule. They agree to `5e-15` relative.
- `tests/test_asymmetric_hull.py` — 35 tests: bit-for-bit back-compatibility, the taper
  integral against Simpson and the half moment across two schemes, the involution on the
  reference and on the mesh, monotonicity and continuity in the exponent, mesh vs reference
  on a four-grid ladder, closedness of asymmetric meshes, parameter validation, metadata,
  refinement plumbing and the CLI.
- `examples/asymmetric_hull/run_example.py` → `docs/asymmetric_hull/{results.json,report.md,plot.png}`,
  exiting 1 if any of its 13 claims stops holding; `docs/asymmetric_hull.md`.
- CLI: `--fore-exponent`, `--aft-exponent`.

### Measured (`docs/asymmetric_hull/`)
- mesh LCB error for `(p_f, p_a) = (0.5, 0.9)`: `-7.79 mm` at 41x13, `-2.92` at 81x25,
  `-1.08` at 161x49, `-0.39` at 321x97 — one sign throughout, so a discretisation bias and
  not noise; `0.011 %` of the length at 161x49.
- the involution `LCB(p_a, p_f) = L - LCB(p_f, p_a)` holds to `1e-12 m` on the reference and
  `1e-6 m` on the mesh, while the mesh volume is unchanged to `1e-15` under the swap.
- `LCB == LCF` on an asymmetric mesh to `1e-9` — on the symmetric family that is one
  identity, here it is two different sums that have to agree.
- volume error `-0.95 %` at 41x13 and `-0.12 %` at 161x49 for `(0.5, 0.9)`; the analytic
  volume for the symmetric default is still `7.512698221250894 m^3` exactly.


## 0.16.0 — Bonjean curves: a station's immersed area against the waterline

`section_area_at` could only report the area of a station the mesh had a vertex
ring on, at the full depth of that ring; nothing in the package could answer "at
station x, how much area is under waterline z", which is the input to damage
stability and to the buoyancy curve of a still-water bending calculation.

- `tensorcad/bonjean.py`: `section_loops` (plane cut of a closed mesh, segments
  matched by shared edge key, plane nudged 1e-9 off vertex rings), `section_properties`
  (clipped area, perimeter, centroid), `bonjean_curve`, `bonjean_table` with
  `integrated_volume` (Simpson when the stations are even and odd in number,
  trapezoid otherwise), and `analytic_section_area` — the closed form
  `(2/3) B sin(pi x/L)^0.7 zeta^(3/2) / sqrt(draft)` for the parametric family.
- Verification: bit-exact on a barge; −5.31 % worst on the default 41x13 grid at the
  shallowest waterline, converging as 5.31 / 2.54 / 0.92 / 0.33 % for nz = 13/21/41/81;
  the trapezoid integral of the station areas equals the surface-flux displacement to
  1.1e-16 relative (the loft is linear in x, so the identity is exact), leaving the
  −0.16 % residual as surface sampling only.
- `docs/bonjean.md`, `docs/bonjean/{results.json,plot.png}`,
  `examples/bonjean/run_example.py`, `tests/test_bonjean.py` (37 tests).
- Additive only: no existing key, function or number changed.

## 0.14.0 — 2026-11-03

### Added
- **`tensorcad.stability` — initial stability (KM, GM) from the waterplane moments.**
  `waterplane_moments(mesh, waterline_z)` cuts the surface with the plane, chains the
  segments into closed loops on the mesh's own vertex/edge identity, and integrates the
  area and its first/second moments by Green's theorem — exact for the polygonal
