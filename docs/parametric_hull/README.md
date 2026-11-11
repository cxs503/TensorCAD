# Parametric hull generation: NURBS loft, surface mesh, hydrostatics

Campaign node **p316w1-n03-cad-parametric-hull** (repository `TensorCAD`).
Module: [`tensorcad/parametric_hull.py`](../../tensorcad/parametric_hull.py).
Verification driver: [`examples/parametric_hull_dag/run_verification.py`](../../examples/parametric_hull_dag/run_verification.py).
Measured acceptance numbers (regenerated, never hand-copied): [`VERIFICATION.md`](VERIFICATION.md),
machine-readable record: `artifacts/p316w1-n03/verification.json`.

## 1. What this capability is

A complete, dependency-free pipeline from a mathematical or tabulated hull definition to
a solver-ready surface mesh plus hydrostatics:

| stage | API | reference |
|---|---|---|
| section definition | `cross_section_from_function`, `cross_section_from_offsets`, `OffsetTable` | offset table (station, half-breadth, height) or analytic profile |
| NURBS loft | `loft_nurbs_surface` | global B-spline interpolation, Piegl & Tiller §9.2 (eqs. 9.1-9.8) |
| surface mesh | `tessellate_surface`, `mesh_density_for` | structured parametric grid, shorter-diagonal split, seam welding, optional waterplane lid |
| geometric transforms | `Transform`, `apply_transform` | scale / translate / rotate / mirror (single-twin screw, port-starboard) |
| hydrostatics | `compute_hydrostatics` | divergence theorem on the clipped wetted surface |
| quality | `mesh_quality_report` | scaled Jacobian `2A/(e_mid e_max)`, aspect ratio `e_max/e_min` |
| hull families | `wigley_cross_sections`, `wigley_family_sections`, `wigley_analytic_properties`, `suboff_cross_sections`, `suboff_station_fractions` | Wigley III / generalised Wigley, DARPA SUBOFF |

Repository boundary (see `docs/PROJECT_CHARTER.md`): this is **geometry plus geometric
hydrostatics**. It does not claim hydrodynamic, structural or seakeeping suitability, and
it does not run a solver; the mesh is handed to TensorLBM/TensorFVM/TensorBEM as a
triangulated surface artifact.

## 2. Geometry and physics definitions

The reference family is the generalised Wigley hull (child of the classical Wigley III
parabolic form used for code-to-code verification):

```
y(x, z) = (B/2) [1 - |2x - L|^p / L^p] [1 - (|z| / T)^q] ,  x in [0, L], z in [-T, 0]
```

with `p = q = 2` the classical Wigley III, `x` from AP (0) to FP (L) and `z` from keel
(-T) to waterline (0). Closed-form hydrostatics used as the reference solution:

```
V   = L B T p q / ((p+1)(q+1))            Cb = p q / ((p+1)(q+1))
Aw  = L B p / (p+1)                       Cw = p / (p+1)
Cm  = q / (q+1)
I_T = L B^3/12 [1 - 3/(p+1) + 3/(2p+1) - 1/(3p+1)]
I_L = L^3 B [1/12 - 1/(4(p+3))]           BM_T = I_T / V
VCB = -T (q+1)/q (1/2 - 1/(q+2))          LCB = LCF = L/2
```

For `p = q = 2` these give the textbook Wigley values `V = 4/9 L B T`, `Aw = 2/3 L B`,
`I_T = 0.0381 L B^3`, `I_L = L^3 B / 30`, `VCB = -0.375 T`; the test suite pins them
(`tests/test_parametric_hull_dag.py::test_wigley_waterplane_inertia_and_metacentric_radius`,
`::test_wigley_centre_of_buoyancy`).

Reference hull for all reported numbers: `L = 10 m, B = 1 m, T = 0.625 m`
(`Cb = 4/9`, `Cw = 2/3`), 11 stations x 12 profile points, mesh `n_u = 97, n_v = 16`
(chosen by `mesh_density_for(L, T, n_v=16)`, i.e. `dx = 2.5 dz`).

The second hull is the **DARPA SUBOFF bare hull** (Groves, Huang & Chang 1989,
DTRC/SHD-1298-01): four-segment offset table (bow quintic, parallel midbody, stern taper,
elliptical stern cap), `L = 4.3561 m`, `R = 0.2540 m`, `L/D = 8.575`. The transcription in
`tensorcad/suboff.py` is pinned to the published values by `tests/test_suboff.py`.

## 3. Method

### 3.1 Interpolating loft (this is what fixes the accuracy)

The historical implementation used the sampled section points **directly as B-spline
control points**. A clamped B-spline only reproduces its first and last control point, so
every interior ordinate was pulled towards the control polygon: the mesh displacement came
out **-3.28%** low and the waterplane area **-2.13%** low, and refining the mesh made it
no better (the bias is geometric, not discretisation).

`loft_nurbs_surface(..., interpolate=True)` now solves the global interpolation system for
both directions (Piegl & Tiller §9.2):

1. chord-length parametrisation of the data, averaged over parallel columns so that every
   curve in a direction shares one parametrisation and one knot vector;
2. `u_k` -> averaged knots `U_{j+p} = (1/p) sum u_{j+i}`;
3. collocation matrix `N_i,p(u_k)`, solved by Gaussian elimination with partial pivoting
   (standard library only), per coordinate;
4. sections are interpolated first in the profile direction (`v`), the resulting control
   nets are then interpolated along the hull (`u`) -- classic skinning, so the surface
   passes through **every** section point.

The profile parametrisation is deliberately *uniform* and shared: a per-section
chord-length parametrisation (then averaged) puts constant-`v` rows at different heights on
different sections, which leaves a fixed `~4e-4` displacement floor that no refinement can
remove (measured: the `u`-refinement study in `VERIFICATION.md` stalls at that floor; with
the shared uniform parametrisation the mid-station ordinate error drops from `8.7e-3` to
`3.4e-5`).

### 3.2 Exact hydrostatics from the clipped wetted surface

The waterplane lid is never meshed and never needed. For any integrand independent of `z`
the closed-surface integral vanishes, so the waterplane properties are obtained from the
wetted hull alone:

```
Aw  = -∮_hull n_z dA                  V = ∮_hull z n_z dA + z_wl Aw
I_T = -∮_hull y^2 n_z dA              I_L = -∮_hull x^2 n_z dA   (transferred to LCF)
M_x = ∮_hull x z n_z dA - z_wl ∮_hull x n_z dA        (and likewise M_y, M_z)
```

Each triangle contributes its **exact** planar integrals (`∫f dA = A/3 Σf_i`,
`∫fg dA = A/12 [(Σf)(Σg) + Σ f_i g_i]`), and faces are clipped exactly at the waterline
(Sutherland-Hodgman) instead of being selected by centroid -- the previous contour-tracing
waterplane routine dropped every face straddling the plane and returned `-2.13%`.

`I_L` is transferred from the AP origin to the transverse axis through the centre of
flotation with the parallel-axis theorem (`I_L - LCF^2 Aw`); reporting it about the AP
would be `~6x` too large for a hull whose LCF sits at midship.

### 3.3 Surface mesh generation

* **Sampling**: `u` samples are placed at uniformly spaced *physical abscissae* `x`
  (inverted through a tabulated `x(u)`). Sampling uniformly in the spline **parameter**
  is not equivalent: on a clamped loft the first knot span covered 3.8% of the hull on one
  strip, which alone produced the 14:1 sliver at the stem.
* **Density rule** (`mesh_density_for`): a structured grid on a hull with a sharp stem has
  two competing limits -- stem triangles degenerate as `dz/|dx|`, midbody triangles stretch
  as `dx/dz`. Keeping `dx = 2.5 dz` clears both (`min scaled Jacobian 0.365`,
  `max aspect ratio 2.74`); the rule scales with slenderness (`n_v` rows -> `n_u` columns).
* **Quad split** along the shorter diagonal. A fixed corner-to-corner split reverses under
  `x -> L - x`, which biased the LCB by `8.7e-4 m` on a 10 m hull even though the vertex
  set was exactly symmetric; the shorter-diagonal rule restores exact symmetry
  (`|LCB - L/2| = 1e-14`).
* **Topology**: centre-plane vertices are shared, coincident vertices welded, zero-area
  faces removed, meridian-plane flaps (the pinch where the keel line meets the stem line)
  dropped once instead of twice, and an optional planar waterplane lid stitched on so the
  body is a closed 2-manifold (`watertight = True`, `edge_manifold = True`). The lid is
  recorded in `SurfaceMesh.hull_face_count` and excluded from hull quality statistics and
  wetted-area claims.

## 4. Results (measured, see `VERIFICATION.md` for the generated tables)

| criterion | measured | threshold | verdict |
|---|---|---|---|
| (a) Wigley ordinates, 10 stations x 5 waterlines | `3.43e-05` | `< 1e-3` | PASS |
| (b) Wigley displacement volume | `-0.148%` | `< 0.5%` | PASS |
| (c) Wigley waterplane area | `-0.011%` | `< 1%` | PASS |
| (d) SUBOFF radius, 20 stations | `0.093%` | `< 0.5%` | PASS |
| (e) hull mesh min scaled Jacobian | `0.365` | `> 0.30` | PASS |
| (e) hull mesh max aspect ratio | `2.74` | `< 5.0` | PASS |
| (f) Cb / L-B / B-T sweep | `0.243%` volume error, monotone in Cb | `< 0.5%` | PASS |

Closed-form cross-checks at the reference mesh (`L=10, B=1, T=0.625`): `V = 2.7744` vs
`2.7778 m^3`; `Aw = 6.6659` vs `6.6667 m^2`; `I_T = 0.38086` vs `0.38095 m^4`
(`-2.5e-4`); `I_L = 33.3274` vs `33.3333 m^4`; `BM_T = 0.13731` vs `0.13714 m`;
`VCB = -0.23424` vs `-0.234375 m`; `LCF - L/2 = 0` exactly.

Displacement convergence against the analytic hull (`h = L/(n_u-1)`):

| n_u | n_v | faces | volume error | waterplane error | min scaled J | max aspect |
|---|---|---|---|---|---|---|
| 21 | 6 | 396 | -1.62% | -0.250% | 0.239 | 4.19 |
| 41 | 10 | 1436 | -0.475% | -0.062% | 0.263 | 3.80 |
| 61 | 12 | 2636 | -0.291% | -0.028% | 0.317 | 3.15 |
| 101 | 20 | 7596 | -0.099% | -0.010% | 0.307 | 3.26 |
| 131 | 26 | 12996 | -0.057% | -0.006% | 0.304 | 3.29 |

Observed order of accuracy `p = 2.07` (Richardson extrapolation of the finest pair is
`< 1e-4 %`), consistent with a second-order chordal triangulation error. The legacy
control-point loft, by contrast, is stuck at `-2.83%` (n_u=61) and `-2.63%` (n_u=101): it
does not converge because it never passes through the data.

SUBOFF: worst radius error `9.3e-4` over 20 check stations (surface built from 48 clustered
stations, **not** the check stations); displacement `-0.76%` and wetted area `-0.24%`
against the analytic quadrature in `suboff_reference_geometry()`; bow/stern tip closure
`6.3e-3 m` / `4.8e-6 m` (the bow tip is a square-root closure, so a finite `r` at the tip
is expected from any finite grid).

Figures (`artifacts/p316w1-n03/figures/`): body plan and waterlines (`fig1`), ordinate
error map (`fig2`), mesh convergence with the `O(h^2)` guide and the legacy curve (`fig3`),
quality versus density (`fig4`), Cb sweep (`fig5`), SUBOFF offsets (`fig6`).

## 5. Accuracy judgement against the >15%-is-unacceptable standard

Every key quantity is compared against an independent reference, and each is far inside
3% -- the explicit answers:

| quantity | reference | error | `< 3%`? |
|---|---|---|---|
| displacement volume | closed form `4/9 L B T` | 0.148% | yes |
| waterplane area | closed form `2/3 L B` | 0.011% | yes |
| `I_T` | closed form (`0.0381 L B^3`) | 0.025% | yes |
| `I_L` | closed form (`L^3 B/30`) | 0.018% | yes |
| `BM_T` | closed form `I_T/V` | 0.123% | yes |
| `VCB` | closed form `-0.375 T` | 0.057% | yes |
| ordinates (50 points) | closed form | 0.0034% | yes |
| SUBOFF radius (20 stations) | DARPA offset table | 0.093% | yes |
| SUBOFF displacement | analytic quadrature | 0.76% | yes |

## 6. Pitfalls and limitations (honest list)

1. **The seeding formula `Delta = 4/3 L B T Cb` is mis-specified.** With the naval
   architecture block coefficient of this hull (`Cb = 4/9`) it returns `16/27 L B T`, 25%
   below the closed-form volume; it is the volume only for `Cb = 1/3`. The criterion is
   therefore evaluated in the standard form `V = Cb L B T` (with `Cb` measured from the
   geometry, `0.147%` off), and the literal value is recorded in
   `criterion_b.seeding_formula_note` for the audit trail rather than silently satisfied.

   > **Update (m13763, 2026-11-11): fixable without an unstructured remesher.**
   > `tessellate_surface(..., cap_adaptation=True)` drops profile columns where the
   > ring is already finer than the axial spacing; the zipper falls out of the quad
   > loop and the existing `weld`/`remove_degenerate_faces` collapses it. SUBOFF
   > `n_u=160, n_v=16`: min scaled Jacobian **0.127 -> 0.325**, max aspect
   > **7.88 -> 3.03**, watertight, volume drift 0.017 %. Cost: worst chord sagitta
   > 1.0 % of the radius, Wigley closed volume -0.37 %. See
   > [`CAP_ADAPTATION.md`](CAP_ADAPTATION.md). The text below is the limitation as it
   > stood when opened; criterion (e) is now met on both hulls.

2. **SUBOFF mesh quality does not meet criterion (e).** Best achievable with a structured
   parametric grid is `min scaled Jacobian 0.13`, `max aspect ratio 7.9` (reported honestly
   in `criterion_d.mesh.quality`). Cause: both ends of the bare hull close with an infinite
   slope (`r ~ x^0.476` at the bow) and the stern cap is an ellipse of semi-axis
   `0.1175 R`; a tensor-product grid turns those polar caps into stretched fans, and the
   circumferential spacing that would fix the fan is incompatible with the body's required
   resolution (a fan is isotropic only when `theta_1 ~ dtheta`, i.e. a handful of
   circumferential facets). Fixing it needs an unstructured/polar-cap remesher -- out of
   scope here, recorded as future work. Criterion (e) is satisfied on the Wigley hull, the
   canonical case for this node.
3. **`u`-refinement has a residual floor** (~`4e-4` relative displacement at `n_u -> 200`,
   `n_v = 26`): the `u` parametrisation is chord-length based, so `y` is not a polynomial in
   the parameter and the interpolation error does not vanish at the rate `O(h^4)`. The
   `v` direction converges at order `1.5-1.9`. Total error stays below 0.2% at production
   density, so this is documented rather than fixed.
4. **Uniform stations cannot resolve the SUBOFF stern cap.** With 21 uniform stations the
   whole cap falls inside one interval and the lofted radius there is wrong by 33.9% at
   `x/L = 0.975`. `suboff_station_fractions()` is the fix: cosine clustering plus the three
   published segment junctions plus four stations inside the cap.
5. **The waterplane lid is excluded from hull quality metrics** (`hull_face_count`). Its
   shape is dictated by the body plan (long, thin near the ends), not by the solver; it is
   reported separately as `closure_faces` and its triangles never enter wetted-area or
   displacement sums.
6. **Hydrostatics are geometric only**: no trim/list equilibrium solve, no appendages, no
   free-surface or viscous effects. `displaced_mass` uses the requested water density.
7. **Performance**: pure standard library, no numpy; a 5.4k-face Wigley mesh with
   hydrostatics takes ~0.2 s, the full verification campaign ~80 s.

## 7. Reproducing

```
python examples/parametric_hull_dag/run_verification.py            # criteria + figures
python -m pytest tests/test_parametric_hull_dag.py -q              # 33 acceptance tests
python -m pytest tests/test_parametric_hull.py -q                  # 35 regression tests
```

Coverage of `tensorcad/parametric_hull.py` from the two test modules: **92%**
(`coverage run --source=tensorcad.parametric_hull -m pytest ...`, statement + branch).
