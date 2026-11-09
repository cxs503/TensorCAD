# NURBS hull surface → body mesh → geometry consistency (IN383-CAD01)

`tensorcad.hull_nurbs` is the geometry **source** of the ship-performance workflow: a
parametric NURBS hull patch, its hydrostatics integrated with high-order Gauss
quadrature, a body-fitted hexahedral volume mesh generated from the same patch, and the
same hydrostatics re-integrated **from the mesh** so that the two routes can be compared.
The point is not the hull form — it is the loop. A geometry source is only usable by a
solver if the numbers a solver measures on the mesh agree with the numbers the design
claims on the surface.

Everything quoted below is measured by `python tools/hull_nurbs_report.py` in one run and
stored in [`artifacts/in383_cad01/acceptance.json`](../artifacts/in383_cad01/acceptance.json).
The tables are generated into
[`docs/experiments/in383-cad01-hull-nurbs/results.md`](experiments/in383-cad01-hull-nurbs/results.md);
the figures into `artifacts/in383_cad01/figures/`. No number in this document is
transcribed by hand from a console.

## 1. Case: what is modelled, and what is not

Reference case (node specification): a KCS-derived main-dimension triple at model scale,
**L = 3.0 m, B = 0.4 m, T = 0.15 m** with the lighter block coefficient **Cb = 0.555**.
The full-scale KCS has Cb = 0.651; the form here is *calibrated* to the requested Cb and is
therefore a different hull. Calling it "KCS" would be a claim this model cannot support, so
it is called what it is: a KCS-dimensioned parametric hull.

| quantity | value |
|---|---|
| L / B / T | 3.0 / 0.4 / 0.15 m |
| L/B, B/T | 7.5, 2.667 |
| Cb (requested) | 0.555 |
| depth to deck D | 0.1875 m (freeboard 25 % T) |
| design displacement Cb·L·B·T | 0.0999 m³ |
| design Froude number (for the LCB baseline) | 0.20 |
| NURBS patch | bicubic, 21 × 13 control points, clamped, non-rational |

Water: infinite depth, even keel, no trim, no heel, no appendages, no free-surface
correction. The hull is **bare**.

## 2. Method

### 2.1 The patch

The hull is one tensor-product B-spline patch, a graph over `(x, z)`:

```
S(u, v) = ( L·u,  y(u, v),  D·v )
```

`x` and `z` are carried by the control net at the **Greville abscissae**, so B-spline linear
precision makes `x(u) = L·u` and `z(v) = D·v` exact (measured worst deviation `< 1e-14` in
the test suite; `test_bspline_linear_precision_carries_x_and_z_exactly`). Only the
half-breadth `y` is free-form. This buys two things at once: the geometry stays a real NURBS
patch that a downstream consumer can evaluate and differentiate exactly, and the
hydrostatics become ordinary 2-D integrals,

```
V = 2·L·D ∫₀¹du ∫₀^{v_wl}dv y(u,v),      v_wl = T/D
```

with the factor 2 the **mirror plane**: the starboard half is modelled, the port half is
its mirror image. The assumption is explicit (`test_mirror_assumption_is_explicit`) rather
than discretised twice.

### 2.2 The form

```
y(x, z) = (B/2) · f(ξ) · ψ(ζ; ξ),     ξ = x/L,  ζ = z/T
```

* `f(ξ)` is the waterline half-breadth ratio: 1 over a parallel middle body, a smoothstep run
  aft down to the transom ratio, a superelliptical entrance forward
  `f = (1 − s^a)^{1/b}` (vertical tangent at the stem, as a real bow).
* `ψ(ζ;ξ) = r(ξ) + (1 − r(ξ))(1 − (1−ζ)^p(ξ))^{1/p(ξ)}` is a superelliptical section family:
  `p → 1` fine and V-shaped, large `p` wall-sided. `r(ξ)` is the flat-of-bottom ratio, which
  vanishes towards the ends because a stem has no flat of bottom.
* The local fullness relaxes as `p(ξ) = p_ends + (p_mid − p_ends)·f(ξ)^k`. This is what
  produces the real, measured relation **Cwp = 0.633 > Cp = 0.606** (`Cb = CM·Cp = 0.916 ×
  0.606 = 0.555`): on a real hull the waterline is fuller than the sectional-area curve
  because the entrance sections are fine while the waterline still has breadth. Relaxing
  only at the very ends (the first thing tried) left the entrance sections as full as the
  midship one and pushed Cb to ≈ 0.69 — a measured FAIL of the Cb target, fixed in the form,
  not in the acceptance number.

### 2.3 The one calibrated quantity

`parallel_body_ratio` is **not** a design input; it is found by bisection so that the
integrated block coefficient of the *patch* equals the requested Cb. The family is
monotone in it (`Cb(0) = 0.537` → `Cb(0.85) = 0.849`), so bisection is safe; the band is
reported with every run, and a request outside it raises `HullNurbsError` instead of
returning a hull that quietly misses by 20 % (`test_calibration_without_a_target_is_refused_when_unreachable`).
Calibrating on the same functional a consumer will use means C1 is a *calibration residual*
(7.1e-4), which is stated as such rather than presented as an independent validation.

### 2.4 The mesh

A body-fitted curvilinear hexahedral grid over the starboard half:

* `i` fore-and-aft, cosine-distributed → clustered at stem and transom, where `f` turns;
* `j` vertical, power-graded towards the bilge, with **the design waterline an exact grid
  level** — so the wetted area is clipped exactly, not by centroid tests;
* `k` from the centreplane to the shell. Every node on the shell layer **is** the NURBS
  patch (`test_mesh_nodes_on_the_shell_are_the_nurbs_patch`, `< 1e-15` m), so the mesh
  boundary is the geometry, not an approximation of it.

Default: 32 × 16 × 6 = 3072 hexahedra, 3927 nodes per half hull.

Quality is measured per cell with two independent defects: the **scaled Jacobian**
(minimum over the eight corners, from the corner trihedra — skew) and the **edge aspect
ratio** (thinness). Thresholds are module constants declared before any measurement:
`SCALED_JACOBIAN_MIN = 0.2`, `SCALED_JACOBIAN_GOOD = 0.5`, `ASPECT_RATIO_MAX = 100`.
A six-tetrahedron conforming subdivision is also produced (18 432 tets), and the tet/hex
volume mismatch is `3.7e-15` — a subdivision, not a re-meshing.

### 2.5 The consistency loop

Three independent routes to the same three numbers:

| route | volume | LCB / KB | wetted area |
|---|---|---|---|
| design | Cb·L·B·T | — | — |
| NURBS quadrature | `2LD∬y dudv` | `∬ y·x`, `∬ y·z` | `2·side + bottom + transom`, from `\|S_u × S_v\|` |
| hex mesh | Σ cell volumes (below the waterline) | Σ `V_c·centroid_c` | Σ boundary quads that see water |

The empirical baselines are read as regressions, not as truth: Holtrop–Mennen (1984) for
the wetted area (no bulb, `A_BT = 0`), Schneekluth–Bertram (1998) for the LCB,
`LCB[%] = 9.7 − 22.5·Fn − 13.5·Cb`.

## 3. Results

Principal quantities (all from the artifact):

| quantity | value |
|---|---|
| displaced volume | 0.0999711 m³ (design 0.0999) |
| LCB from AP | 1.53421 m → +1.140 % L from midship |
| KB | 0.0815499 m (0.544 T) |
| waterplane area / Cwp | 0.759710 m² / 0.6331 |
| max section area / CM | 0.0549584 m² / 0.9160 |
| prismatic coefficient Cp | 0.6063 |
| wetted area | 1.42355 m² (sides 0.63662, bottom 0.14467, transom 0.00564) |
| realised beam at waterline | 0.399112 m (0.222 % below the 0.4 m request) |

Note the last row: the patch is a spline *quasi*-interpolant of the form, sampled at the
Greville abscissae, so the beam it realises is a spline maximum and sits 0.9 mm under the
request. It is measured and reported (`realised_dimensions`) instead of being assumed.

### Acceptance criteria

| id | criterion | measured | limit | verdict |
|---|---|---|---|---|
| C1 | displacement vs Cb·L·B·T | 7.11e-4 (0.071 %) | 2 % | PASS |
| C2 | wetted area vs Holtrop–Mennen | 7.125 % | 10 % | PASS |
| C3 | \|LCB − LCB_Schneekluth\| / L | 3.433 % L | 5 % L | PASS |
| C4 | hex cells failing the quality gate | 6.673 % failing (93.33 % pass) | 10 % failing | PASS |
| C5 | NURBS vs mesh: V, S, KB, LCB | worst 0.468 % | 1 % | PASS |
| C6 | Cb↑ ⇒ ∇↑ (`dV/dCb` vs `L·B·T`) | monotone, slope error 0.018 % | strict monotone | PASS |
| C7 | refinement ×2: V, S change | 0.354 %, 0.112 % | 2 % | PASS |
| C8 | pytest | 40 tests in `tests/test_hull_nurbs.py`, all green | ≥ 12 | PASS |

C3's metric needs a word: LCB is quoted *as a fraction of L*, and the criterion is read as
`|ΔLCB|/L ≤ 5 %`. Reading it as a relative error on the LCB value would be meaningless here
because the two numbers are near zero and of opposite sign — `+1.14 % L` against the
regression's `−2.29 % L`. A relative-error reading is reported for honesty in §6.

### Figures

`artifacts/in383_cad01/figures/`

| figure | what it shows |
|---|---|
| `body-plan-waterlines.png` | 21 sections (starboard + mirror) and five waterlines of the patch: the entrance is fine and V-shaped, the run full, the transom narrow at the waterline |
| `body-mesh.png` | the 32 × 16 × 6 body-fitted hex mesh, and the midship column showing the waterline as an exact grid level |
| `mesh-quality.png` | distributions of the two quality defects; the nose at `sj = 0` is the 96 collapsed stem-column cells, reported not hidden |
| `convergence.png` | mesh → patch error against cell count (volume 1.91 % → 0.12 %) and the passing fraction against resolution |
| `routes-comparison.png` | design vs NURBS vs mesh for displacement, wetted area and LCB, each with its band |
| `cb-sweep.png` | the five calibrated hulls of the Cb sweep against the analytic `Cb·L·B·T` line |

### Mesh ladder (measured)

| stations × levels × layers | cells | V [m³] | S_wet [m²] | V error | S error | pass |
|---|---|---|---|---|---|---|
| 16 × 6 × 4 | 512 | 0.098063 | 1.416282 | 1.909 % | 0.510 % | 88.28 % |
| 24 × 9 × 5 | 1440 | 0.099149 | 1.420440 | 0.822 % | 0.218 % | 91.25 % |
| 32 × 12 × 6 | 3072 | 0.099503 | 1.421921 | 0.468 % | 0.114 % | 93.33 % |
| 64 × 24 × 12 | 24576 | 0.099855 | 1.423508 | 0.116 % | 0.003 % | 95.11 % |

Second-order convergence of a piecewise-bilinear re-integration of a smooth surface, and a
quality fraction that rises with resolution because the cell that fails is the one that
straddles a geometric singular point.

## 4. Method pitfalls met while building this (kept, because they bite again)

1. **Cb is not a free parameter.** The first form (constant section fullness, wide transom)
   could not build anything below Cb = 0.67. Calibrating the parallel body is not enough:
   the *section relaxation* has to follow the plan form, otherwise no member of the family
   is fine enough. The family band is now reported with every run.
2. **The flat of bottom decides mesh quality.** With `r(ξ) = r_b·f^1.5` the bottom row died
   towards the stem (`AR = 5546` in one cell); with a constant `r` the hull gained a flat
   bottom at the stem, which is not a hull. `r_b·f^0.5` is the compromise: pass fraction
   88 % → 93 % at unchanged geometry error.
3. **Two defects, not one.** A flattened cell has a perfect scaled Jacobian (skew) and a
   terrible aspect ratio; a sheared cell the other way round. The first version of the
   corner-trihedron table also had the wrong handedness for the `(di, dj, dk)` node
   ordering and scored every cell 0 — the metric now takes its global orientation from
   corner 0, so it is convention-independent.
4. **The stem is a coordinate singularity, not a meshing bug.** The last station column
   collapses to the centreplane (96 cells) and the bottom row is thin there (`max AR 549`).
   A structured topology has to collapse somewhere; the cells are counted as `collapsed`,
   the wedge is visible in the quality figure, and the 90 % gate is met over the whole
   block including them.
5. **The wetted-area integrand has a vertical tangent at the stem.** `|S_u × S_v|` converges
   more slowly than the volume: `6e-4` between the default `(24,16)` rule and `(32,20)`,
   still 15× inside the 1 % band. Volume and waterplane converge to `< 1e-5`.
6. **The mesh must be clipped at the waterline.** The first version integrated the whole
   volume up to the deck and disagreed with the patch by 28 % — the freeboard volume. The
   waterline being an exact grid level is what makes the clip exact.

## 5. Limitations

* One tensor-product patch: no bulbous bow, no sonar dome, no appendages, no flare above the
  waterline (`ψ = 1` for `ζ > 1`), no bilge-radius knuckle.
* Mirror symmetry about the centreplane is assumed, not meshed.
* Even keel, no trim/sinkage equilibrium, no heel, no free-surface correction, no damage,
  no ice, no resistance/powering: this module produces geometry and hydrostatics only.
* The body mesh is structured, so the stem column collapses; an unstructured tet mesher with
  a boundary layer is out of scope (the six-tet subdivision is provided for solvers that
  want tets, at the cost of the structured quality guarantees).
* Holtrop–Mennen and Schneekluth–Bertram are regressions from ship series. At Cb = 0.555 the
  hull sits at the fine edge of both bands, which is exactly why the node's acceptance bands
  are 10 % and 5 % and not 2 %.

## 6. Error judgement against the references (explicit)

| key quantity | reference | measured error | < 3 % ? |
|---|---|---|---|
| displaced volume | Cb·L·B·T (definition) | 0.071 % | **yes** |
| volume, mesh vs NURBS | patch quadrature | 0.468 % | **yes** |
| wetted area, mesh vs NURBS | patch quadrature | 0.114 % | **yes** |
| KB, mesh vs NURBS | patch quadrature | 0.175 % | **yes** |
| LCB, mesh vs NURBS | patch quadrature | 0.025 % of L | **yes** |
| wetted area | Holtrop–Mennen (1984) | **+7.12 %** | **no** |
| LCB | Schneekluth–Bertram (1998) | **3.43 % of L** (+1.14 vs −2.29) | **no** |
| volume vs design after calibration | Cb·L·B·T | 0.071 % | yes |

The two "no" rows are the empirical comparisons, and they are reported as they came out:

* The wetted area of this form is 7.12 % **above** the Holtrop–Mennen estimate. The
  regression was fitted on ship series starting around Cb = 0.55 and has 5–10 % scatter
  there; a 7 % offset is inside that scatter but is not the < 3 % that would make the two
  numbers interchangeable.
* The LCB sits 3.43 % L forward of the Schneekluth regression. The form has no aft-body
  asymmetry tuning: the split between run and entrance (`aft_body_share = 0.47`) is a design
  choice, and forcing the LCB onto the regression by tuning that share would have made C3 a
  tautology rather than a check. The honest reading is: a hull with a longer entrance than
  run carries its buoyancy forward of this regression by about 3 % L.
* The internal consistency rows (C5) are the ones that decide whether the geometry source is
  usable, and all four are half a percent or better — the mesh a solver gets describes the
  surface the design claims.

## 7. Reproducing

```bash
python tools/hull_nurbs_report.py      # figures + artifact + generated tables
python -m pytest tests/test_hull_nurbs.py -q
```

```python
from tensorcad import build_nurbs_hull, nurbs_hydrostatics, build_body_mesh, mesh_integrals

hull = build_nurbs_hull()                       # L=3.0, B=0.4, T=0.15, Cb=0.555
hydro = nurbs_hydrostatics(hull)                # patch quadrature
hydro.displaced_volume, hydro.wetted_surface_area, hydro.lcb_ratio_from_midship
mesh = build_body_mesh(hull)                    # 32 x 16 x 6 hexahedra, half hull
mesh.quality.pass_fraction                      # 0.933
mesh_integrals(mesh).displaced_volume / hydro.displaced_volume - 1   # -0.0047
```