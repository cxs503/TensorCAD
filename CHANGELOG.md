# Changelog

All notable changes to TensorCAD are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
semantic versioning.

## 0.22.0 — 2026-11-05

### Added
- **A second, independent route to `GZ(phi)`** (`tensorcad/righting_arm.py`), carried over
  from PR #40 (`agent/m3169`) and redeveloped on the current `main`: the heeled section is
  clipped face by face and integrated in the body frame (divergence theorem with no cap,
  each contribution exact for a planar triangle), and the equilibrium waterline is found by
  **bracket-and-bisect** on the monotone `V(h)` instead of `tensorcad/gz.py`'s Newton
  iteration on the projected wetted area. It adds `heeled_displacement`, `righting_arm`
  (one heel, with `lcg`), `RightingArmCurve.peak/vanishing_crossed/area_under_curve`,
  `wall_sided_validity_deg` and per-point volume residuals.
- **The overlap between the two routes is now the check.** `tools/verify_merge_m4038.py`
  puts both on the same meshes and records the agreement: `1.2e-13 m` on a wall-sided barge
  (`L=10, B=2, T=0.8 m`), `1.6e-13 m` on the `41x13` hull and `6.9e-14 m` on `81x25`
  (floated at `z=-0.47 m`), `2.8e-16 m` on a sphere of resolution 32 -- round-off, because
  both integrals are exact for the polygonal hull and only the waterline solver differs.
  A **negative control** is part of the record: shifting `gz`'s displacement target by 10 %
  moves `GZ(20 deg)` by 3.2 mm, so the comparison is capable of failing.
  Evidence: `artifacts/m4038/{merge_verification.json,merge_verification.png}` and
  `docs/m4038_pr40_merge_cure.md`.
- **Each route against its own reference on the merged tree.** The barge reproduces the
  wall-sided closed form to `1.2e-14 m` inside its validity window (21.80 deg here) and the
  formula is 9.6 cm out at 40 deg; the sphere sits within `1.3e-04 m` of `GM sin(phi)` at
  resolution 32; the hull's small-angle slope matches `KB + I_T/nabla - KG` from
  `tensorcad.stability` (a plane-section integration of the *same* mesh) to `6.8e-10`
  relative on three grids; the volume residual stays below `2.1e-12`.
- `tests/test_righting_arm.py` (46 tests: PR #40's 41 plus 5 cross-route/alias tests),
  `examples/righting_arm/run_example.py` (10 claims, exits non-zero if one breaks, writes
  `docs/righting_arm/{results.json,report.md}` and three figures), `docs/righting_arm.md`
  and `docs/righting_arm/`.

### Changed
- `tensorcad/__init__.py` exports the new route **without shadowing `gz`**: the three
  colliding names are published as `bisection_righting_arm_curve`,
  `BisectionRightingArmCurve` and `barge_wall_sided_righting_arm`, while
  `righting_arm_curve`, `RightingArmCurve` and `wall_sided_righting_arm` keep pointing at
  `tensorcad/gz.py`. `tests/test_righting_arm.py` pins both directions.
- `README.md`: the "what it does not provide" sentence still listed the righting-arm curve
  and free-surface stability as absent, and claimed no hull family has a transom; corrected
  to what the tree has done since 0.21.0 (and to what it still does not do: damage and
  flooded-compartment stability, free-trim equilibrium, class-rule criteria).

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


`stability.py` stopped at `GM` and said so; the backlog tracked 大倾角 GZ as missing. This
release adds the curve itself, by integrating the **heeled displacement** instead of a
waterplane moment.

- `tensorcad/gz.py`: `heeled_hydrostatics` (rotate about the keel line, cut at the waterline,
  solve the waterline that restores the displacement by Newton with the exact derivative
  `dV/dc = A(c)`; two to five iterations to a `1e-13` relative volume residual),
  `righting_arm_curve` (`GZ = KN − KG sin φ` with `KG` above the keel, matching
  `initial_stability`'s datum), `max_gz` / `angle_of_max_gz_deg` / `range_of_stability_deg`
  (interpolated) / `area_under_curve` (dynamical stability, m·rad, truncatable).
- **Verified against two closed forms**, both exact for the polygonal mesh: the box barge
  wall-sided trapezoid `GZ = sin φ (GM₀ + ½BM₀ tan²φ)` — *derived* in the module docstring,
  not quoted — agrees to `6.9e-17 m` over 0–35° on a 10×2×2 m barge, and provably stops
  agreeing past its 38.66° limit (4.0e-04 m at 45°, 1.25 m at 70°); and the half-submerged
  log `GZ = (R − KG) sin φ`, exact at all 33 symmetry-aligned heels from 0 to 90°
  (`1.1e-16 m`). Third check: `KB` (0.4830051 m) and `dGZ/dφ` at zero (0.421622996 m/rad vs
  `GM` 0.421622996 m) both reproduce `initial_stability`, which integrates the waterplane
  loops by Green's theorem instead.
- `with_freeboard`: a hull generated to its design waterline cannot be heeled (the waterline
  leaves the modelled surface on the dipping side immediately), so this removes the waterline
  cap, extrudes the boundary ring with vertical walls and closes it with a fan — watertight,
  outward, Euler-2, upright displacement unchanged to round-off.
- `free_surface_correction` / `rectangular_tank_inertia`: the flat `Σi_x/V` reduction a slack
  tank takes off every point of the curve (2.0 × 0.5 m tank in the 7.45 m³ example: 2.80 mm).
- Evidence: `docs/gz.md`, `examples/gz/run_example.py` (writes `results.json`, `report.md`,
  `plot.png`; exits non-zero if a claim breaks), `tests/test_gz.py` (26 tests; the suite
  pins the public surface, which grew by 14 names).

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
  waterplane the mesh has. Sections are referred to their own centroid
  (`I_T = ∫y²dA − A·y_c²`), so asymmetric and heeled waterplanes integrate correctly;
  nested loops are detected by point-in-polygon parity and subtract.
  `initial_stability(mesh, kg=...)` returns `KB`, `BM_T`, `BM_L`, `KM` and `GM` with `KG`
  as a required loading input (no silent `KG = 0`). Free-surface corrections, `GZ(φ)`,
  trim and class rules are explicitly *not* implemented, and `basis` says `solid GM`.
- **`hydrostatics_analytic.waterplane_inertia_transverse` / `waterplane_inertia_longitudinal`** —
  the closed form `I_T = (B³L/12)·taper_integral(3e)` (note the `3e`, not `e`) and the
  quadrature `I_L = 2BL³∫₀^{1/2} u² cos(πu)^e du` with its self-convergence estimate.
- `examples/stability/run_example.py` writes `docs/stability/{results.json,report.md}` plus
  three figures and exits non-zero if one of its ten claims stops holding;
  `tests/test_stability.py` adds 19 tests; `docs/stability.md` derives the formulas, reports
  the ladder and names the limits.

### Measured
- **Box** (`10 × 2 × 0.8 m`, `KG = 0.55 m`): every waterplane and metacentric quantity
  reproduces its closed form to **2.2e-16** — one ulp — so the section integration is exact
  and the hull's residual error is the mesh's.
- **Hemisphere** (`R = 1 m`, 32×64): `KB = 0.625075 m` vs `5R/8`, `BM = 0.374624 m` vs
  `3R/8`, `KM = 0.999699 m` vs `R` (**−0.030 %**): two independent integrations summing to
  the sphere's radius.
- **Hull** (`10 × 2 × 0.8 m`, `KG = 0.6 m`): `I_T` `−0.139 %` / `−0.0086 %`, `I_L`
  `−0.604 %` / `−0.064 %`, `BM_T` `+0.748 %` / `+0.127 %`, `GM_T` `+1.99 %` / `+0.35 %` at
  `41×13` / `161×41`. Fitted order of `I_T` in the grid step: **2.01**. `I_L`'s reference is
  rule-independent: tanh-sinh vs Composite Simpson under `u = 1/2 − s²` agree to **7.8e-16**.
- **`GM` amplifies.** `rel(GM) ≈ (KM/GM)·rel(KM) ≈ 2.87·rel(KM)` at `KG = 0.6 m`, so the
  coarsest rung (`21×9`) leaves the 3 % band at `+3.47 %` while `I_T` is still inside it at
  `−0.56 %`; `GM = 0` falls at `KG = KM = 0.915360 m`, bracketed by the `GM(KG)` ladder.

## 0.13.0 — 2026-11-03

- **Closed-form hydrostatics reference** (`tensorcad/hydrostatics_analytic.py`): the
  coefficients, centres, areas and wetted surface of the analytic hull surface are now
  available as exact formulas (`analytic_hydrostatics`, `taper_integral`) instead of
  only as mesh integrals. `Cb = 2I/3`, `Cp = I`, `Cm = 2/3`, `Cw = I`, `LCB = LCF = L/2`,
  `VCB = -0.4T`, `V = B·L·T·I·2/3`, with `I = ∫₀¹ sin(πx)^0.7 dx` from the gamma
  function. The wetted area needs quadrature and uses a double-exponential rule
  (`tanh_sinh_quadrature`) that handles the `ζ^-1/2` keel edge and the `ξ^-0.3` ends;
  it reports its own convergence (`5e-13` relative).
- **Verification debt closed**: every hydrostatic quantity the mesh integrator reports
  now has an independent reference. Worst error is 0.939 % on the default 41×13 grid,
  0.169 % on 161×41 and 0.061 % on 321×81; the wetted-area reference is independently
  confirmed by a Richardson extrapolation of the mesh series to +8.8e-6.
  See `docs/hydrostatics_verification.md`; the self-checking ladder lives in
  `examples/hydrostatics_verification/run_verification.py`.
- **Release hygiene**: `tensorcad.__version__` had drifted to "0.9.0" while
  `pyproject.toml` said 0.10.0. Both now read 0.13.0, and the public-surface test pins
  the same string in one place.

## 0.12.0 — 2026-11-03

### Added
- **`tensorcad.curvature`** — the fairness check the backlog had left open (P3). Two
  independent estimates, because the question has two audiences:
  - **`principal_curvatures(field, points)`** reads the principal curvatures off the
    analytic field through the shape operator `M = P (Hess f) P / |grad f|`. `P` kills
    the normal direction, so the estimate is invariant to whether the field is a true
    signed distance or only an implicit function — which is what makes it work on the
    ellipsoid and on the learned `.ngo` fields, not just the exact primitives. The
    seven gradients per probe (point plus three axis pairs) go out in one batched call;
    measuring that was worth an order of magnitude on the AFF-8 min-union. Each sample
    is Richardson-extrapolated from `h` and `h/2` and carries its own `step_error`.
  - **`vertex_mean_curvature(mesh)`** is the cotangent-Laplacian estimate a consumer of
    the mesh can recompute, with the standard cotangent clamp against slivers.
  - **`assess_fairness(mesh, field=...)`** returns the audit: the dihedral distribution
    over interior edges, the crease count and length fraction above a threshold, the
    mesh curvature summary, and — when a field is supplied — the curvature jump across
    each interior edge, which is then a property of the geometry rather than of the
    triangulation.
- `tests/test_curvature.py` — 20 tests: the four exact primitives against their closed
  forms, the invariance to `|grad f| != 1`, input rejection, the mesh estimator's
  convergence and its integral conservation, twelve creases on a cube, a sphere with
  none, and the stdlib-only import.
- `examples/fairness/run_example.py` — writes `docs/fairness/{results.json,report.md,plot.png}`
  and exits 1 if one of its 15 claims stops holding.
- `docs/fairness.md` — the derivation, what each estimate is good for, and the limits.

### Measured
- **The analytic estimate is exact**: sphere `1/R, 1/R`, cylinder `1/R, 0`, plane
  `0, 0` (exactly zero), ellipsoid at all three poles — worst relative error
  `1.33e-14`, reported step error below
  `1.3e-06 1/m`.
- **The mesh estimate is exact in the integral and needs a conditioned mesh
  pointwise.** On a raw marching-tet sphere: pointwise median error
  `48.9 %`, but the
  area-weighted mean `1.000490 1/m`
  against the analytic `1.0` (`0.049 %`).
  The cotangent operator is a discrete divergence, so slivers cannot move the integral;
  they own the pointwise field. After `taubin_smooth` + `project_to_field` the pointwise
  median is `1.42 %`.
- **AFF-8** (22060 triangles): p50 dihedral
  `0.049°`, p95
  `11.94°`, min
  `0°` (the NACA trailing edges),
  `1171` crease edges over 15° =
  `3.04 %` of the interior edge length.
- **The parametric hull family is not G1 at its ends.** `b(x) = B/2 sin(pi x/L)^0.7` gives
  `db/dx ~ x^-0.3000` (fit over
  `x/L = 1e-2 .. 1e-7`, expected `-0.3`):
  infinite waterline slope at the stem and stern, so curvature diverges there. The keel
  line has the same problem because the section profile is `sqrt(z/T)`. Both are
  properties of the family, not of a mesher.
- **A crease threshold on a coarse mesh measures the mesh.** On the hull, p95 dihedral
  falls `26.80° -> 3.59°`
  and the crease length fraction `7.55 % -> 1.44 %`
  between `41x13` and `161x41`, while the crease *count* rises
  `232 -> 628`. Use the distribution,
  not the count.

## 0.11.0 — 2026-11-03

The sampled hull stops being its own reference. `analytic_hydrostatics` returns the
hydrostatics of the analytic surface the generator samples as the same
`Hydrostatics` record the mesh path returns, so displacement, waterplane area,
midship section area, LCB, LCF, VCB and the form coefficients get a closed-form
baseline and the mesh numbers get an error rather than a refinement trend.

### Added
- **`tensorcad.hull.analytic_hydrostatics(length, beam, draft, water_density=1025)`**
  — the surface factorises as `(B/2)*sin(pi*x/L)^0.7*sqrt((z+T)/T)`, so every
  volume integral splits and has a closed form: volume
  `L*B*T*TAPER_INTEGRAL*2/3`, waterplane area `L*B*TAPER_INTEGRAL`, midship section
  area `2/3*B*T`, LCB and LCF `L/2`, VCB `-2/5*T`, `Cb = 2/3*TAPER_INTEGRAL` and
  `Cp = TAPER_INTEGRAL`. `basis` names which values came from a formula and which
  from a quadrature.
- **`analytic_wetted_surface_area`** — the one quantity without a closed form,
  integrated with an 80-point Gauss-Legendre rule after `w = u^2` (removes the keel
  `1/sqrt(w)`) and `psi = (pi/2)*tau^3` (removes the bow/stern `x^(p-1)` tip slope)
  put both endpoints inside the rule's convergence radius. Refining to 120 nodes
  moves the value by `< 1e-9` relative. Legendre nodes are generated by Newton
  iteration, so the core stays standard-library only.
- **`analytic_relative_errors(parameters)`** — signed relative error per quantity
  for a sampled grid, plus `outside_discretisation_estimate`, the quantities whose
  true error exceeds the coarse/refined relative change that
  `estimate_discretisation_error` reports.
- **`tests/test_analytic_hydrostatics.py`** — 21 tests: the closed forms against
  independent Simpson quadratures, the wetted-area Gauss rule against a Simpson
  product rule and against node refinement, mesh-vs-closed-form error at four
  resolutions with the convergence direction and sign pinned, and input rejection.
- **`examples/analytic_hydrostatics/run_analytic_check.py`** — writes
  `docs/analytic_hydrostatics/{results.json,report.md,convergence.png}`, exits 1 if
  a claim stops holding.
- **`docs/hull_analytic_reference.md`** — the derivation, how the reference is
  itself checked, and what it does not cover.

### Measured (default grid `41x13`, all within 3 %; `161x41` within 0.2 %)
- volume `-0.880 %` / `-0.136 %`, block coefficient the same by construction;
- waterplane area `-0.175 %` / `-0.017 %`, prismatic coefficient the same;
- wetted surface area `-0.326 %` / `-0.074 %` (reference `22.9127295 m²`);
- midship section area `-0.707 %` / `-0.119 %`, the largest default-grid error;
- VCB `+0.939 %` / `+0.169 %` about its closed form `-2/5*T = -0.32 m`: the
  trapezoidal section rule is high near the keel, where the `sqrt` profile is
  steepest, so this is the one quantity that over-reports;
- LCB and LCF land on `L/2` to `1e-9 m` at every grid tested;
- the refinement estimator is optimistic at the default grid: the true error is
  `1.53x` its reported relative change for volume, `1.77x` for wetted area and
  `1.44x` for waterplane area, so it is an indicator and not a bound.

### Fixed
- The declared and runtime versions had drifted: `pyproject.toml` said `0.10.0`
  while `tensorcad.__version__` still said `0.9.0` and `tests/test_public_surface.py`
  asserted the old string. All three now read `0.11.0`, and the 0.10.0 voxel-export
  and benchmark work that shipped without a changelog entry is noted here rather
  than left silent.

## 0.9.0 — 2026-10-24

Representations stop being an argument and become a measurement. One harness runs an
analytic formula, a mesh edit, a trained neural field and an optional kernel export through
the same rays, the same closed form and the same thresholds, and writes `results.json` plus
a generated `report.md` and `plot.svg`.

### Added
- **`tensorcad.benchmark`** — `Metric` (value, unit, reference, tolerance, status ∈
  measured/imported/not_supported/not_run/failed, and a gate verdict; a missing value can
  never be `0` and never passes a gate), `RepresentationResult`, `Gate`, `CaseResult`,
  `BenchmarkResult` (`passed`, `failures`, `write_json`, `write_markdown`) and
  `run_benchmark(quick=…)`. Standard library only; the neural row imports `tensorcad.neural`
  behind an importability check.
- **Case A, `case_a_sphere_dent`** — a sphere R = 0.5 m with a Gaussian dent of depth
  0.06 m, σ = 0.22 rad, whose radius is a closed form; exact volume and area by Simpson
  quadrature, a Möller–Trumbore ray caster, and a stratified ray set (a third of the rays
  inside the edit cone so a 0.22 rad feature is actually resolved).
- **Case B, `case_b_hull_handoff`** — a displacement target → damped Newton → hull → gmsh
  round trip, measured against `hull.analytic_volume` rather than against its own mesh.
- **`examples/benchmarks/run_benchmark.py`** — writes `results.json`, `report.md` and
  `plot.svg` (standard-library SVG, surface error against wall clock), and exits 1 when any
  gate fails, so it is usable as a CI step.
- **Contract coverage** — the row IDs are the ones the merged contract defines
  (`analytic_field`, `mesh`, `neural_field`, `ocp_step`, `parametric_hull`), P95 alongside
  mean/max, query throughput and artifact bytes, and training steps/time/seed imported from
  the bundle. Peak memory is not measured; contract Case B (AFF-8 round trip) and Case C
  (real solver run) are not implemented and the document says so.

### Measured (full sampling, 512 rays, grid 80×160, in `docs/benchmarks/`)
- mesh edit: cone-mean depth error **2.3e-6 m**, dent surface error mean **4.5e-5 m**, leak
  outside the edit support **exactly 0 m**, all-ray surface error 1.0e-4 m (P95 1.8e-4 m);
- trained neural field (`.ngo`, weights `db1df71e…`, 221,620 bytes): cone-mean depth error
  **1.9e-3 m** (3.2 % of the edit depth), dent surface error mean 2.1e-3 m, imported
  parent-relative leak 1.3e-3 m (P95 3.3e-3 m), all-ray surface error 3.4e-3 m (P95 9.2e-3 m);
- the same 1 %-of-depth gate applied to both: the mesh row passes, the neural row fails, and
  the run exits 1. No looser tolerance was invented for the learned field;
- the shipped grid is 80×160 because its measured chord-sag floor (1.9e-4 m) sits below the
  6.0e-4 m gate; 40×80 happens to pass with a floor *above* the gate, which is luck;
- the bundle's own claim (2.32e-3 m) and this harness's re-measurement (2.12e-3 m) agree,
  which is why imported numbers are kept beside measured ones instead of replacing them;
- case B: 0.880 % volume error at (41, 13), residual 2.2e-16 in 6 evaluations, gmsh node
  round trip exact (0 m) and area to 6.8e-16, refinement +0.577 %.

### Fixed
- `examples/neural_geometry/edit_sphere_dent.py` wrote a bundle that verified but could not
  be evaluated: `state.json` lacked `field_state.radius` and `coordinates` lacked `bounds`,
  so `field_from_object` raised a format error. Both are now written, and re-running the
  example reproduces the same weights hash bit-for-bit.
- `at_most` / `at_least` / `equals` now record the limit they compared against, so a
  published `results.json` states the threshold that produced each verdict.

## 0.8.0 — 2026-10-24

The first release whose input is a *requirement* rather than a parameter set. 0.7.0 could
measure and hand off a surface; 0.8.0 finds the surface that meets a brief and refuses to
deliver one that does not.

### Added
- **`tensorcad.design`** — `DesignTargets` (equalities `displacement_m3`,
  `beam_to_draft`, `length_to_beam`, `draft_m`; budget `wetted_area_m2_max`),
  `DesignProblem` (damped Newton on `log` of the free dimensions, central-difference
  Jacobian, backtracking line search, bounds from the hull model's own specs) and
  `DesignSolution` (`status`, `solved`, `failed_checks`, `convergence_order`,
  `sensitivity`, `manifest`). Every iteration is kept with its hydrostatics, its
  `tensorcad.geometry.handoff/1` record and its full check list. Standard library only.
- **Measured**: recovering L = 12 m, B = 2.4 m, T = 0.9 m from its own measured
  displacement (12.063418 m³), B/T and L/B converges in 6 evaluations and 1.1 s to
  **1e-7 relative** on all three dimensions; observed contraction order **1.940**; the
  finite-difference sensitivities come out at +1.000000 / −1.000000 for a ratio and
  +1.000000 for displacement against beam and draft; a replay reproduces the mesh digest
  `814f4770b0189c58`; re-opening the design for a 15 % heavier brief converges from the
  previous answer in 5 evaluations.
- **The verdict is a conjunction.** A 200 m³ target at L = 10 m drives the residual to
  1.13e-12 and is still refused: at L/B = 0.965 the triangles fail the minimum-angle,
  aspect-ratio and sliver checks, so `solved` is `False`, `failed_checks` names all
  three, and the manifest reports `solver.status: "not_ready"`.
- **Grid sensitivity is part of the record**: the same design measures 12.063418 m³ at
  (41, 13) and 12.133078 m³ at (81, 25) — +0.577 % — so a reader can see how much of the
  answer is the sampling grid.
- **`write_solver_manifest` / `write_design_json`** — the manifest a solver reads
  (`tensorcad.design/1`: targets, every iteration, handoff, reference quantities,
  verification) and the cheap per-iteration digest an agent polls between steps.
- **`examples/design_loop/run_example.py`** — 21 acceptance checks, writes OBJ, gmsh 2.2,
  the FEM JSON and the solver manifest; in CI.

### Notes
- Displacement and wetted area are deliberately not both offered as equalities: on this
  family they both scale with beam times draft, so the pair is nearly rank deficient. The
  loop stalls at a residual of 4.0e-3 and cannot remove it. Wetted area is a budget
  instead, and the docs record the measurement.
- `solver.executed_here` is `false` in every manifest: this repository measures and hands
  off geometry, it does not run a flow or a structure solve.
- 369 tests (355 pass without the `neural` extra; 14 skip), 93.2 % line coverage, strict
  `mypy` clean, `ruff` clean, seven examples in CI.

### Still missing
- A hydrodynamic objective (resistance, stability, seakeeping), multiple hull families,
  and topology search. The loop matches the stated targets and proves the geometry; it
  does not judge whether the ship is good.

## 0.7.0 — 2026-10-24

Mesh conditioning and the two solver exits. The 0.6.0 release closed with "mesh
smoothing/decimation: 1.5–3.2 % of triangles are thinner than 2°; repairing them is the
next piece of work"; this release repairs them, measures the shape drift that repair
causes, and hands the result to TensorFEM (and to any tool that reads gmsh or STEP)
without adding a mandatory dependency.

### Added
- **`tensorcad.mesh_conditioning`** — `taubin_smooth` (λ 0.5 then μ −0.53, a band-pass so
  the body does not shrink), `project_to_field` (one Newton step per vertex along the
  field gradient, with a per-vertex step limit and a 4-step backtracking line search),
  `decimate` (quadric-error collapse under the link, orientation, duplicate-face and
  aspect-ratio guards, with every rejected collapse counted), `measure_mesh` and
  `condition_mesh`, which refuses to return a surface whose watertightness, manifolding,
  Euler characteristic, volume (0.5 %) or area (1 %) moved beyond its declared tolerance.
- **Measured on AFF-8 at `(150, 18, 26)`** (65,812 triangles, 1.664 % slivers):
  Taubin ×30 → min angle 0.0080° → 8.3202°, no slivers, volume −0.127 % of the analytic
  0.70754064 m³ in 0.4 s; + projection ×30 (53.6 s) → 3.6590°, aspect 17.6, volume
  −0.117 %, area −0.198 %, worst displacement 18.4 mm, worst residual `|f|` 0.57 mm;
  + decimation to 45,000 / 30,000 / 20,000 faces in 5.4 / 9.7 / 11.6 s → min angle
  2.70 / 2.56 / 2.40°, zero slivers, volume drift −0.12 / −0.14 / −0.21 %.
- **The decimate-first order is reported, not recommended.** `condition_mesh` at the same
  budgets reaches 2.21° (45 k) and 0.149° with 28 slivers (30 k); at 45 k faces it costs
  40 s against 5 s smoothing-first. Both pass the shape gates; the example prints both
  tables.
- **`tensorcad.adapters.fem`** — solver-neutral handoff needing no extra: node/element
  tables, boundary patches assigned by nearest analytic component, gmsh 2.2 ASCII with one
  physical tag per patch, and a JSON manifest. Round trip: worst node displacement **0.0 m**,
  worst per-patch area difference 2.2e−14 m², digest identical to the source surface.
- **`tensorcad.adapters.ocp`** — optional STEP bridge (`pip install "tensorcad[ocp]"`):
  sew → shell → solid → STEP, with the free-edge count, shell closure, `BRepCheck_Analyzer`
  verdict, kernel volume, written size and re-imported volume all reported. Measured on
  14,764 AFF-8 triangles: 12.8 s sewing, 4.4 s writing, 0 free edges, kernel volume
  0.694708267 m³ against a mesh signed volume of 0.694708267 m³ (**+5.3e−15 relative**),
  re-import −3.1e−13 relative, 37.14 MB. Detached islands are kept and counted
  (`shells`/`solids`) rather than dropped.

### Notes
- The STEP file is a **faceted solid**: one planar face per triangle (14,764 faces). That
  is what a solver-facing exchange file needs; a B-rep feature model is issue #4.
  `write_step` raises `OcpUnavailable` with the pip command when the extra is absent — it
  never writes a mesh file with a `.step` extension.
- 357 tests in the default lane (336 without the `neural` extra), plus 10 OCP tests that
  run when the wheels are installed; 93.2 % line coverage overall, with
  `tensorcad/adapters/ocp.py` at 32 % in the default lane because its tests skip there.
  Strict `mypy` clean on 22 files (OCP is an untyped import behind a documented
  override), `ruff` clean.

### Still missing
- Mesh *refinement* or remeshing to a target size field: conditioning only coarsens.
- A B-rep/STEP feature model, booleans, fillets, sweeps.

## 0.6.0 — 2026-10-24

The DARPA SUBOFF appendages (AFF-8): the sail and the four cruciform stern fins, on top
of the 0.5.0 bare hull. Both are published closed-form predicates, so they are
implemented as formulas — no boolean engine, no CAD kernel, no CSG tree. The three
solids are unioned as `min` of their fields, which is exact for the sign and the zero
set. `examples/suboff_appendages/run_example.py` enforces every claim below;
`docs/suboff_appendages.md` records the numbers.

### Added
- **`tensorcad.suboff_appendages.SuboffSailField`** — the published three-segment sail:
  entrance polynomial, parallel middle, exit polynomial, plus the elliptical cap
  (`Zmax = 0.109375 ft` half width, body top `1.507813 ft`, axial
  `[3.032986, 4.241319] ft` of the 14.291667 ft hull). The section is one radial
  coordinate `ρ = |y|` below the body top and `√(y² + (2(z − YTMP))²)` above it, so the
  cap and the body need no union. `height` (the cap apex, 0.47625 m) is a property.
- **`tensorcad.suboff_appendages.SuboffFinsField`** — the swept NACA 4-digit cruciform:
  chord `c(r) = −0.466308 r + 0.88859` ft, thickness
  `(0.2969, 0.126, 0.3516, 0.2852, 0.1045)`, span `r ∈ [0.075, 0.825] ft`, trailing edge
  at `x = 13.146284 ft`. The planform is a convex quadrilateral in `(x, r)`, whose
  distance is computed exactly (four edge distances), not as a gap.
- **`SuboffAppendagesField` and `SuboffFullField`** — the sail ∪ fins, and the whole
  AFF-8 model hull ∪ sail ∪ fins. `sail_scale` / `fin_scale` reproduce TensorLBM's
  similarity axes (sail about its footprint centre, centreplane and deck plane; fins
  about their trailing edge and root radius), so a scaled appendage stays attached.
- **Verified against the source twice.** The section values match TensorLBM's
  transcription exactly (`0.00e+00 ft` on the sail, `6.9e-18 ft` on the NACA section);
  containment matches the actual torch predicates on 18 hand-picked labels and on
  20,000 random points at each of three scales with **0 mismatches**, against a
  verbatim in-repo transcription of those predicates that CI can run.
- **An independent volume reference.** Quadrature over the published sections gives
  sail 0.00916043 m³ (46.2 % of it outside the hull), fins 0.00785718 m³, and the full
  model 0.70754064 m³ (+1.191 % over the bare hull). The exported
  `(436, 52, 76)` mesh lands at 0.707120 m³ (−0.059 %) with 569,056 triangles,
  6.350646 m² of surface, watertight, edge-manifold, genus 0.
- **The magnitude is measured, not claimed.** The cap and the thickness surface have no
  closed-form nearest point, so the fields declare `quantity = "implicit_function"` and
  report the first-order normalised gap. Along the analytic normals: sail 0.0012 mm
  median error at a 0.5 mm offset, fins 0.00001 mm; `|∇f|` is 0.996 (sail) and 1.000
  (fins) at the median, and exact on the flat faces.

### Notes
- The fins alone are a **genus-1** closed shell (`V − E + F = 0`) and that is the
  published geometry, not a meshing bug: the blades start at `r = 0.075 ft` and touch
  near the root, leaving the axis column empty. The hull fills it, so the full model is
  genus 0. Both are asserted.
- Throughput: 77 k queries/s on a 20,000-point batch for the full model (~1.6× the
  bare hull's cost).
- 314 tests (295 without the `neural` extra), 95.9 % line coverage; `suboff_appendages`
  99 %. Strict `mypy` clean on 19 files, `ruff` clean.

### Still missing
- Mesh smoothing/decimation: 1.5–3.2 % of the triangles are thinner than 2°, on the
  NACA knife edges. Reported as warnings; repairing them is the next piece of work.
- The AFF-2/3/4/7 appendage layouts, explicit intersection curves, and therefore STEP.

## 0.5.0 — 2026-10-24

The first real hull form on the geometry core: the DARPA SUBOFF bare hull (AFF-1) as an
analytic signed-distance field, with its offset-table provenance, an analytic reference,
anisotropic meshing, and a measured export. `examples/suboff/run_example.py` enforces
every claim below; `docs/suboff.md` records the numbers.

### Added
- **`tensorcad.suboff.SuboffHullField`** — the four-segment DARPA offset table (bow
  polynomial, parallel midbody, sextic stern taper, elliptical stern cap) in its
  published feet, carried as a true signed distance in metres. `L = 4.356100 m`,
  `R = 0.254000 m`, `L/D = 8.575`, bow at `x = 0`, stern at `x = length`.
  `|s|` matches a brute-force 400,001-sample meridian to `2.4e-08 m` (0.09 ppm of R);
  the gradient is the exact unit normal (`1.000000`, agreeing with central differences
  to `1.7e-08`).
- **Provenance in the record, not in a document.** `parameters()` / `provenance()`
  carry the source (Groves, Huang & Chang 1989, DTRC/SHD-1298-01), the transcription
  (TensorLBM `suboff_cad.suboff_radius_profile`), the unit convention, the axial origin
  and the profile description. The profile agrees with that transcription to
  `0.00e+00` at every station tested.
- **`suboff_reference_geometry()`** — analytic volume, wetted area, midship area and
  prismatic coefficient by quadrature, reported at two resolutions with their
  difference (volume 0.000 ppm, area 0.005 ppm between 10⁵ and 10⁶ samples).
  `V = 0.699211 m³`, `A = 5.988264 m²`, `Cp = 0.79194`.
- **Anisotropic export.** `to_mesh` and `extract_surface` now accept three cell counts
  per axis, and `sample_grid` marches a non-cubic grid. A 8.575:1 hull at
  `(218, 26, 26)` (147,368 cells) reaches −0.138 % volume / −0.066 % wetted-area
  error with 4.3 % slivers, where an isotropic `64³` (262,144 cells) reaches
  −0.089 % / −0.080 % with 10.0 % slivers and 64 stations along the hull.

### Fixed
- **`test_public_surface` poisoned the interpreter.** It deleted `numpy` from
  `sys.modules` and reloaded `tensorcad` to prove the core needs no third-party stack,
  which left torch's numpy bridge broken for every later test in the session — a
  real order-dependent failure, observed when the file was run before
  `test_geometry_field.py`. The check now runs in a fresh interpreter, like the
  neighbouring opt-in tests.
- **A scale-dependent "degenerate face" criterion.** 0.4.0 counted faces with an area
  below 10⁻¹² m², which is a different geometric claim on a 0.3 m cylinder and a
  4.4 m hull. Tests and benchmarks now use area below 10⁻⁶ of the longest edge
  squared. The corrected 0.4.0 record: the cylinder at 40 cells emitted 8,108 faces
  below 10⁻¹² m² (12.4 %), 2,936 of them with an aspect ratio above 10⁶ — not, as
  0.4.0's notes said, "exactly zero area".

### Notes
- 281 tests (262 without the `neural` extra), 95 % line coverage.
- TensorLBM's `r_over_l` default is `1/(2·8.57) = 0.058343` (`R = 0.254148 m`); the
  offset table gives `0.058310` (`R = 0.254000 m`). This module uses the table, which
  makes its wetted area 0.060 % smaller than TensorLBM's published value. At
  TensorLBM's radius the two implementations agree to `+0.00000 %`, so the whole
  difference is the radius convention.
- Mesh quality is reported, not repaired: at `(218, 26, 26)` 4.3 % of triangles are
  thinner than 2°, concentrated at the two tips where the meridian slope is infinite.
  There is no smoothing or decimation yet.

## 0.4.0 — 2026-10-24

The continuous geometry core (CGC v0.1): analytic and learned geometry now answer the
same questions, so they can be interchanged, compared and handed to a solver without
each consumer learning a second API. Everything here is measured by
`examples/geometry_core/run_benchmark.py`, which fails the build if a claim stops
being true; see `docs/geometry-core.md` for the protocol and the numbers.

### Added
- **One protocol for every representation.** `tensorcad.geometry_field.GeometryField`
  defines batch `signed_distance`, `gradient`, `bounds`, editable `state`/`with_state`,
  `parameters`, `provenance`, `metadata`, `contains`, `numerical_gradient` and an
  optional `numpy_callable` fast path. Declared per field, not assumed: `kind`,
  `quantity`, `inside_sign`, `surface_isovalue`, `units`, `coordinate_frame`.
- **Four exact primitives** — `PlaneField`, `SphereField`, `CylinderField` (capped,
  arbitrary axis), `EllipsoidField`. Measured against independently recomputed closed
  forms, the difference is `0.00e+00` in double precision. Gradients are analytic and
  match central differences to `2.7e-11 … 2.0e-10`. `EllipsoidField` declares
  `quantity="implicit_function"`: its zero set is exact, its magnitude is not a distance.
- **Edits that provably do not drift.** `with_state` returns a new field and leaves the
  original untouched. Measured: sphere radius 0.5 → 0.65 moves the surface 0.150000 m in
  all 24 probe directions with the centre unchanged; a cylinder height edit leaves the
  barrel value bit-identical; an ellipsoid axis edit leaves the other two axes alone;
  unknown controls raise `FieldError` instead of being ignored.
- **`to_mesh(field, resolution)`** — field → closed mesh → the existing hull acceptance
  pipeline. At 40 cells per axis: sphere −0.129 %, capped cylinder −0.167 %, ellipsoid
  −0.134 % volume error; all watertight, edge-manifold, `V − E + F = 2`, no hard errors.
- **`inspect_mesh_quality`** — smallest interior angle, worst aspect ratio and sliver
  share as *warnings*: a closed solid can still be a poor discretisation.
- **`SampledField.handoff()` and `mesh_digest()`** — a self-describing record (schema,
  mesh SHA-256 over float64 bits, units, frame, counts, topology, volume, area, bounds,
  validation, quality) for the TensorFEM handoff.
- **`surface_point(field, direction)`** — ray-bisected surface position, so
  representations are compared by where their surfaces are rather than by field values
  that are not comparable.
- **`NeuralGeometryField` / `field_from_object`** — a published `.ngo` revision behind
  the same protocol, with `quantity="truncated_signed_distance"`, autograd gradients and
  conditioned edits. Measured on the committed bundle: surface 3.98 mm mean / 9.33 mm
  worst from the analytic sphere away from the dent; dent axis 441.73 mm against a
  440.00 mm target; conditioned radius 0.45/0.55/0.60 reproduced within 0.33/0.08/0.01 %.

### Fixed
- **Squashed triangles from grid-vertex coincidence.** A radius-0.3 cylinder in a ±0.3
  box is exactly 20 cells across, so the level set passed through grid vertices and
  marching tetrahedra emitted 8,108 faces with an area below 10⁻¹² m² (12.4 % of
  65,144), 2,936 of them with an aspect ratio above 10⁶, and a 22.4 % sliver share.
  `to_mesh` now offsets the sampling box by an irrational fraction of a cell per axis:
  every count → 0, sliver share → 2.3 %, volume error unchanged at −0.17 %.
- **Finite differences are the wrong derivative for a float32 field.** At `step=1e-6`
  a central difference of the network is dominated by rounding (0.023 disagreement with
  autograd); at `1e-4` it agrees to `3e-4`. `NeuralGeometryField.gradient` now uses
  torch autograd, and `numerical_gradient`'s docstring states the step requirement.

### Notes
- The learned field's `|∇f|` ranges 0.65–1.39, i.e. it is not a metric: distance queries
  must read the field value and surface queries the zero set. Beyond the truncation band
  the outward cosine drops from 0.983 to 0.721 — the truncation caveat, measured.
- The bundle's recorded `dent_peak_measured_m = 0.052690` is parent-relative; re-measuring
  against the analytic sphere gives 0.058268 m. Both are correct and both are printed
  side by side; the difference is the parent's own 5.58 mm axis error.
- 262 tests (245 without the `neural` extra), 95 % line / 91 % branch coverage.

## 0.3.1 — 2026-10-24

The review of 0.3.0 found four defects that could crash a legal model, corrupt an
active bundle or delete an export. They are fixed here, with regression tests that
fail on the old behaviour; nothing was re-thresholded to make a gate pass.

### Fixed
- **Convergence estimate crashed on a legal model.** `estimate_discretisation_error`
  guarded the requested grid (`nx*nz`), then built the *refined* parameters, which
  could exceed the `4001 x 513` limits (`nx=2002, nz=3` raised
  `nx must be <= 4001 stations, got 4003`). It now reports `skipped` with the
  refined grid in the reason. The refined parameter set also keeps
  `water_density` instead of silently falling back to 1025.
- **Half-written export on failure.** `export_hull` wrote the OBJ and CSV before
  the metadata, so a failure during reporting left a directory with a mesh and no
  record of the parameters that produced it. The metadata document is now built
  first and all three files are written afterwards.
- **`tensorcad-hull` had a half-open error boundary.** Only `generate_hull` was
  inside a `try`; export and metadata failures escaped as tracebacks. Every phase
  now reports `stage` + `error` (JSON on stdout under `--json`, human text on
  stderr) and exits 2. The CLI also computes the metadata document once instead of
  generating the refined grid twice.
- **`tensorcad-neural verify` required torch.** The `build_field` import sat
  outside the error boundary, so a torch-free environment got a `ModuleNotFoundError`
  traceback even though the manifest, hash and history checks need no neural
  extra. Those checks now run everywhere, `field_query` reports `skipped` when the
  extra is absent, and `verify` gained derived-artifact integrity checks
  (presence, SHA-256, size) as warnings.
- **Derived meshes could be written into the live bundle before publishing.**
  `apply_dent_edit` extracted straight into `bundle/derived/` and then committed,
  so a rejected or failed commit left files behind that no manifest described.
  Extraction now happens in a scratch directory and travels inside the staged
  bundle.
- **Extraction was invisible and destructible.** `extract_object` wrote
  `derived/surface.obj|stl` without registering them, and the next commit -- which
  swaps the whole bundle directory -- deleted them silently. Extraction into the
  object now publishes its own revision (`lineage.operation = "extract"`, weights
  untouched) and registers the records; an edit carries the parent's registered
  artifacts forward, marked with `source_revision`. `extract --output DIR` remains
  an unregistered export, which is now stated in the docstring.
- **`parse_manifest` did not validate `derived[]`.** Readers accepted an escaping
  `../../../../etc/passwd` path, a non-hex digest or a negative size, so `info`
  would print them. Records are now validated exactly like the writer builds them.
  `ObjectStore.history()` also rejects `NaN`/`Infinity` instead of accepting the
  tokens this project refuses to write.
- **`max_section_area` was quadratic on non-grid meshes.** It re-scanned every
  vertex per station, which is invisible on a lofted hull and ruinous on an
  extracted implicit surface: 14.1 s of a 14.1 s `compute_hydrostatics` at 64^3.
  One bucketing pass makes it `O(V + stations)`; the same call now takes 0.34 s and
  returns the identical numbers (a test pins it against the old definition).

### Added
- Coarse-grid warning: `generate_hull` compares the mesh-integrated volume against
  the analytic volume of the generating functions and warns above 2 % (`nx=3, nz=3`
  is -35.7 % and used to report `warnings == []`). `analytic_volume` is public, and
  the benchmark example now imports it instead of duplicating the formula.
- `tensorcad/py.typed`, so the strict typing this package is developed under is
  actually visible to PEP 561 consumers.
- Regression tests for every item above, including a parametrised overflow matrix
  and a "derived files survive an edit" test.

### Changed
- Test suite: 194 tests without the neural extra, 213 with it (was 175 / 189);
  96 % branch coverage.
- Corpus bundles: `verify` reports recorded-but-absent derived meshes as warnings.
  The committed experiment bundle intentionally omits the 20 MB exports, so it now
  shows two warnings -- which is the honest description of what is in the repo.

## 0.3.0 — 2026-10-24

Geometry can now be stored as an object rather than a mesh: a learned continuous
field plus the metadata, versioning and transactions needed to trust it. STL/OBJ
become derived exports instead of the source of truth.

### Added
- `tensorcad.neural_geometry` (standard library only) — the NGO v0.2 bundle format:
  manifest with lineage, field semantics, state, training, validation, provenance
  and derived-artifact records; SHA-256 verification; path-escape rejection;
  `history.json`; `ObjectStore` with stage → validate → commit/rollback, atomic
  directory swaps, optimistic concurrency (`expected_parent`), and JSON-safe
  history records for failed attempts.
- `tensorcad.neural` (optional `neural` extra) — `ConditionalTSDF` field runtime,
  base training, replay-protected local edits, measurement of dent peak / dent-region
  error / protected-region mean, P95, P99, max, and a marching-tetrahedra
  `extract_surface` that emits an already-outward-oriented, watertight mesh.
- `tensorcad.cli_neural` / `tensorcad-neural` — `fit`, `edit`, `extract`, `info`,
  `verify`; every command prints measured numbers, supports `--json`, and exits
  non-zero when a gate fails (a rolled-back edit exits 1, not 0).
- `write_stl` (binary STL) for derived meshes, alongside the existing OBJ writer.
- `examples/neural_geometry/run_object_pipeline.py` — the end-to-end benchmark:
  fit → edit → extract → verify → deliberate rejection, with assertions.
- `schemas/neural-geometry-manifest-v0.2.schema.json` and
  `docs/neural-geometry-format.md` describing what the format stores, what it
  replaces, and what it does not.
- Coverage gate (`fail_under = 90`, currently 96 %) and a torch-free `neural`
  CLI help check in the default test job.

### Changed
- Manifest v0.2 keeps the 0.1 top-level file keys, so older bundles still load;
  unknown fields are preserved rather than dropped.
- Protected-region gates now include P99 and max, not only the mean, and failed
  edits report `failed_checks` instead of a single generic reason.
- The dent metrics are split into `dent_peak_*` (on the axis) and
  `dent_region_error_mean_m` (inside the cone), so a peak target is never compared
  against a region average.
- The dent cone is always sampled with a deterministic stratified grid, so the
  dent metrics do not silently become infinite when few random rays land in it.
- Torch↔NumPy conversion falls back to list conversion when the installed wheels
  have mismatched ABIs (torch 2.2 with NumPy 2.x) rather than failing a run.
- Test suite: 189 tests (from 142), 96 % branch coverage on the measured modules.

### Known limits
- The only fitted shape family is the sphere (with radial state and local dents);
  ellipsoids, cylinders and hulls are not wired into the object layer yet.
- Field accuracy is about 3 % of the radius at 0.5 m; extraction at 96³ takes
  several seconds on CPU.

## 0.2.0 — 2026-10-24

The exported surface becomes a solid a solver can actually use, and the checks
that matter are now enforced rather than documented.

### Fixed
- **Inverted winding.** The port half of the hull was lofted with the same face
  order as the starboard half, so its normals pointed inward and the mesh's
  signed volume cancelled to exactly zero. Solved structurally: one closed
  section ring per station, plus a generic `orient_outward` pass, instead of two
  mirrored sheets.
- **Tolerance welding artefacts.** Merging the two mirrored sheets by a floating
  point tolerance produced coincident duplicates along the keel and one
  non-manifold edge. Vertices are now shared by exact position at build time.
- **Open, unclosed hull.** The surface ended at the waterline with no deck face;
  it is now capped, closed, edge-manifold, and has Euler characteristic V-E+F = 2.
- **Non-exact waterline plane.** `z = -draft + draft * (nz-1)/(nz-1)` left the
  waterline vertices about 1e-16 off the plane, breaking coplanarity assumptions;
  the mapping is now `z = -draft * (1 - fraction)`, which is exact.
- **Public API.** `export_hull` was not importable from `tensorcad` although it was
  the documented entry point; it is re-exported with the rest of the API.
- **OBJ precision silently altered the model.** Files were written with 9
  significant digits; export now uses 17 so `read_obj(write_obj(mesh))` returns the
  identical vertex table (a test locks this in).

### Added
- `tensorcad.geometry` / `tensorcad.mesh` / `tensorcad.parameters` /
  `tensorcad.validation` / `tensorcad.export` / `tensorcad.hydrostatics` /
  `tensorcad.adapters` — the module boundaries the README promised in 0.1 are now
  real, with generic mesh algorithms kept out of the ship-specific model.
- Validity checks that a geometry foundation has to provide: watertightness, edge
  and vertex manifoldness, outward orientation, zero-area faces, Euler
  characteristic, volume-bounded-by-bbox, and realised-versus-requested L/B/T.
  `generate_hull` raises instead of returning a model that fails them.
- Hydrostatics integrated from the mesh with the divergence theorem: displaced
  volume, displacement mass, wetted surface, waterplane area, LCB, LCF, VCB, Cb,
  Cp, and maximum section area — all reported in `metadata.json`.
- Measured discretisation evidence: `estimate_discretisation_error` regenerates on
  a refined grid and records the relative change (default grid: 0.88 % under the
  analytic volume, 0.14 % at nx=161/nz=41).
- Dimensionless ratio checks (`L/B`, `B/T`, `T/L`) reported as warnings, and hard
  bounds on `nx`/`nz` so an accidental `1e8` fails fast instead of allocating.
- Explicit `nx` parity reporting: with an even station count no station lands on
  midship, so the realised beam falls short; that is now a named warning instead of
  a silent ~1 % dimensional error.
- CLI: `--density`, `--json`, `--no-export`, `--no-convergence`, non-zero exit
  codes with the full validation list on stderr.
- `examples/parametric_hull/run_example.py`: a runnable benchmark that verifies
  closure, orientation, dimensional fidelity, hydrostatics against the analytic
  value, and a resolution study — and exits non-zero if a claim fails.
- `LICENSE` (Apache-2.0 text was declared in `pyproject.toml` but absent),
  `CHANGELOG.md`, `.gitignore`, ruff + mypy (strict) + coverage (90 % gate, currently 99 %)
  configuration, and a Python 3.13 CI matrix entry with a byte-identical
  reproducibility check.
- 129 tests, up from 8.

### Changed
- `HullModel` holds a `Mesh` (`vertices`/`faces` remain as accessors), plus
  `validation` and `hydrostatics`; `metadata()` no longer re-runs validation on every
  call.
- Face indices in the API are 0-based; the OBJ writer applies the 1-based offset.
- `validate_hull(model)` keeps its signature and returns the list of error strings.
- `metadata.json` schema version is `2.0`.

## 0.1.0 — 2026-10-23

- Initial package: deterministic parametric sampled hull surface, OBJ/CSV/JSON
  export, CLI entry point, structural validation, and CI.

## Unreleased

Export becomes a verified step instead of a dump.

### Added
- **`tensorcad.gate`** — `gate_mesh` / `enforce_gate` / `MeshGateReport` / `MeshGateError`:
  a hard-fail gate over closed, edge- and vertex-manifold, single-component,
  `V - E + F = 2`, outward-oriented, finite, zero-area-free meshes, plus optional
  reference-volume, facet-count and sliver limits. Every reference is external to the mesh
  (analytic displacement, requested component count, genus-0 topology).
- **`tensorcad.hull.resolve_resolution`** — `ResolutionPlan`: the coarsest grid whose
  integrated displacement is within a requested relative error of `hull.analytic_volume`,
  with an honest `accepted = False` when the generator caps or a triangle budget is reached.
- **`export_hull(..., stl=True, target_volume_error=..., gate=...)`** — resolves the grid,
  regenerates, gates, writes `hull.stl`, reads it back and re-runs the gate on the recovered
  soup; the plan, gate evidence and read-back volume error land in `metadata.json`.
- **`tensorcad.export.read_stl`** — binary STL reader that rebuilds shared vertices by exact
  float32 welding, so an exported soup can be checked for watertightness; ASCII STL and
  truncated/non-finite files are rejected with explicit messages.
- **`examples/hull_stl_export/run_example.py`** — resolution ladder against the analytic
  volume, a rejected 178-triangle export, a resolution-controlled STL export and its
  read-back gate; exits non-zero if any claim fails.
- **CLI** — `--stl`, `--gate`, `--target-volume-error`.

### Measured (reference hull, `L = 10 m`, `B = 2 m`, `T = 0.8 m`; analytic 7.512698 m³)
- 11×5 / 178 triangles: −5.290 % (gate rejects at the default 2 % tolerance);
- 41×13 / 1998 triangles: −0.880 %; 161×41 / 25 918 triangles: −0.136 %;
- `target_volume_error = 0.002` from an 11×5 request resolves to 161×65 / 41 278 triangles,
  −0.0759 %;
- STL round trip: watertight, `V - E + F = 2`, volume error **2.4e-8** relative (float32).

### Fixed
- The STL header no longer claims every STL came from a neural object.
### Added — waterplane second moments, and the transom-stern hull family
### Added — large-angle stability: the righting-arm curve, and two exact benchmarks
### Large-angle stability: the righting-arm curve, and two exact benchmarks
## 0.23.0 — 2026-11-04
- **Trim, and the two constants that price it** (`tensorcad/trim.py`): `TPC`
  (`rho A_wp / 1e5`), `MCT` (`rho_t I_L / (100 L)`), `LCB - LCF`, a
  `hydrostatic_table` over drafts, and `equilibrium_trim` — the angle a hull takes
  when `weight` is shifted `shift` along it. `hydrostatic_particulars` now reports
  volume, LCB, KB, KM and the centimetre constants from one waterline;
  `trimmed_hydrostatics` / `waterline_for_displacement` float a hull at a trim angle
  with the displaced volume held by a safeguarded secant.
- **No linearisation anywhere in the solve**: at a fixed station "below the tilted
  plane" is a half plane, so the station clip is reused at a per-station waterline.
  On the `20x4x5` m box that reproduces the wallsided closed form at 8 degrees of
  trim to 5.4e-15 relative, and the two centimetre routes (`w d / MCT` against
  `100 L tan(theta)`) agree to 2.0e-16.
- **Two routes, one answer**: section areas integrated in `x` (Simpson) against the
  plane section's own moments. All nine box particulars reproduce their closed forms
  to <=5.1e-16; on the parametric hull `TPC` equals `dV/dz` to 1.3e-3 and `MCT` at
  the design waterline equals the analytic `I_L` route to 6.0e-3; the trim solve
  agrees with the small-angle answer to 4.0e-3 at 0.25 deg.
- `_SectionCutter` indexes stations to the faces that can straddle them: a sweep
  drops from 75 ms to 37 ms on the default hull while reproducing
  `bonjean.section_properties` to 1e-12 (asserted).
- Honest limits, in `docs/trim.md`: a hull capped at its own waterline cannot be
  trimmed at that displacement (`TrimError`, not a wrong angle), ring waterlines are
  nudged by 1e-6 of the mesh scale and flagged `waterline_nudged`, and a
  `station_count` equal to the mesh's `nx` puts every plane on a vertex ring and
  costs 0.4 % on the change of trim.
- `docs/trim.md`, `docs/trim/{results.json,trim.png}`,
  `examples/trim/run_example.py` (exits non-zero if any of its 20 claims breaks),
  `tests/test_trim.py` (24 tests, incl. the indexed cutter against the reference
  route, the refusals, and the stdlib-only guarantee).
  segments into closed loops on the mesh's own vertex/edge identity, and
