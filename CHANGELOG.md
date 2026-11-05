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
- **Righting arms at finite heel** (`tensorcad.stability.righting_arm_curve`, `GZCurve`,
  `GZPoint`, `SCHEMA_GZ`, `DEFAULT_HEEL_ANGLES`): the integration the initial-stability solver
  deliberately left out. At each heel the body is rotated about its longitudinal axis and the
  waterplane height is solved so that the clipped solid displaces the same volume, giving
  `GZ(phi) = Y_G - Y_B` with `Y_G = -(keel_z + KG) sin(phi)`; volume, `Y` and `Z` moments come from
  exact divergence-theorem integrals over the waterplane-clipped triangles (Sutherland-Hodgman
  clipping, fan triangulation, no cap mesh needed), and the plane solve is Illinois regula falsi
  converging to `1e-11` of the displacement. Curve features (`gz_max`, `angle_of_gz_max_deg`,
  `vanishing_angle_deg`, `dynamic_stability_area`, per-point `deck_immersed`) are reported with the
  grid they were read from.
- **Wall-sided topsides** (`freeboard=`): `tensorcad.hull` and `tensorcad.wigley` produce wetted
  surfaces closed at the design waterline, so a body capped there cannot heel at fixed displacement.
  `freeboard=<h>` extends the design waterplane vertically — stated in `basis`, refused when the
  mesh already carries topsides, and refused for multi-loop waterplanes rather than guessing.
- **Verification**: on a `10 x 2 x 0.5 + 1.5` m box the measured arm matches an independent 2-D
  section clipping to **1.49e-14 m over 0..180 deg**; the classical wall-sided formula is recovered
  to **4.4e-16 relative inside `tan(phi) < 2T/B = 0.5`** and departs by 2.1 % one step outside it
  (asserted, not tolerated); on the `10 x 2 x 0.8 m` parametric hull the curve's `GM` equals
  `initial_stability`'s to the last digit (`0.321622996 m`) while the 1-degree secant slope is
  0.99 % off `GM`, the `O(phi^2)` truncation.
  Evidence: `docs/gz_curve.md`, `docs/gz_curve/{results.json,gz_curves.png}`,
  `examples/gz_curve/run_example.py` (exits non-zero if any of its claims breaks),
  `tests/test_gz_curve.py` (61 tests).

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
  waterplane the mesh has. Sec
