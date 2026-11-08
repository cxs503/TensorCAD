# Closed NURBS patch with singular pole rings: the sphere anchor

**Mission** m6439 / TensorCAD. **Module** `tensorcad/nurbs.py`
(`sphere_surface`, `NURBSSurface.singular_boundary_point`, the pole-welding branch of
`NURBSSurface.tessellate`; plus a one-line contract fix in
`tensorcad/geometry.py::MeshBuilder.add_quad`). **Artifacts**
`artifacts/nurbs_sphere/` — `metrics.json`, `fig1_pole_shell.png`,
`fig2_convergence.png`, `fig3_pole_physics.png`, `sphere_solid.stl`. **Tests**
`tests/test_nurbs_pole_surface.py` (27 cases).

This closes the follow-up that `docs/m6410_selfstart_verification.md` left open: the
surface service had an anchor for an *open* tube and for a closed cylinder whose caps
are triangle fans added **by hand**, but no anchor for a closed patch that pinches
**by itself**. That is the singular case, and it is where the tessellator had a real
defect (§5.2).

---

## 1. Case: why a sphere, and what it stresses that a cylinder does not

A hull surface must eventually be closed — a solver needs a body, not a sheet. The
cheap way to close a patch is to bolt caps on (`cylinder_solid_mesh` does exactly
that). The structural way is to let the parametrisation pinch: one row of the control
net collapses to a single point, and the patch becomes a closed surface of its own.
Every hull form introduces that situation — a stem, a transom corner, a bulb tip, a
skeg root — so the tessellator has to survive it.

The minimal object that has nothing in it *except* that difficulty is the sphere: two
singular pole rings, no boundary at all, a closed formula for both area and volume,
and a curvature that makes any normal error visible.

| property | value | why this value |
|---|---|---|
| object | rational biquadratic sphere, **one** patch | no cap is added by hand |
| radius | `R = 12 m` | PC6 bilge radius `B/2 = 12 m`, same scale as the tube anchor |
| degree | `(2, 2)` | exact arcs, not an approximation |
| control net | `7 x 5` points | 3 x 120° circle arcs, 2 x 90° meridian arcs |
| weights | `w_i^circle * w_j^meridian`, circle `1, 1/2, 1, ...`, meridian `1, sqrt(1/2), 1, sqrt(1/2), 1` | the quotient separates (§2) |
| poles | `S(u, v_min) = (0, 0, -R)`, `S(u, v_max) = (0, 0, +R)` | singular control rows |
| analytic anchors | `A = 4 pi R^2 = 1809.5573684677208 m^2`, `V = 4/3 pi R^3 = 7238.229473870882 m^3` | closed form |

## 2. Method

**Why the patch is exact.** Put the unit circle in the `xy` plane with control points
`(X_i, Y_i)` and weights `w_i`, and a meridian semicircle in a half-plane with control
points `(r_j, z_j)` and weights `m_j`. Build the net `P_ij = R (r_j X_i, r_j Y_i, z_j)`
with weights `w_ij = w_i m_j`. Because the weights separate, both sums factor:

```
S(u,v) = [ sum_ij w_i m_j N_i(u) N_j(v) R (r_j X_i, r_j Y_i, z_j) ] / [ (sum_i w_i N_i)(sum_j m_j N_j) ]
       = R ( r(v) X(u), r(v) Y(u), z(v) )
```

where `r(v), z(v)` is the exact meridian and `(X(u), Y(u))` the exact unit circle. So
`|S|^2 = R^2 (r^2 (X^2+Y^2) + z^2) = R^2 (r^2 + z^2) = R^2` identically — the patch is
the sphere, and the measured radial error `5.33e-15 m` at `R = 12 m` is the round-off
of that identity, not a fitting residual.

**Pole rows are a control-net fact, not a tolerance.** `singular_boundary_point(axis,
upper=)` returns the point when *all* control points of a boundary row coincide, else
`None`. That row forces the boundary curve to that point, so `tessellate` keys the
whole row on one vertex instead of on the per-sample evaluation. The evaluation at
`v = v_min` is not bit-for-bit constant: the rational quotient `(a·b)/b` rounds twice,
so the 49 samples of the pole row produce **3** distinct positions spanning **1 ulp**
(`fig3`, left). Exact-position welding — the module's deliberate "no tolerance"
policy — then leaves a ragged ring; keying the row is a modelling fact and keeps the
policy intact.

**Degenerate quads.** A pole quad `(a, b, c, d)` has `a == b`, so only `(a, c, d)`
survives: one triangle instead of two, counted in `Tessellation.degenerate_skipped`
(`2 x u_samples` per tessellation). `MeshBuilder.add_quad` now **returns the number of
triangles it actually appended**; before this mission it returned `2` while appending
one, which mis-indexed every later quad (§5.2).

**Audit chain.** `tensorcad/nurbs.py` → `examples/nurbs_sphere/run_example.py` →
`artifacts/nurbs_sphere/metrics.json` → this document → `tests/test_nurbs_pole_surface.py`.
Every number below is printed by the second step and re-asserted by the last one; the
tests recompute their values instead of reading the artifact, so a merge gate that
filters an artifact out cannot turn them red.

## 3. Results

### 3.1 Analytic anchors

| quantity | computed | analytic | relative error |
|---|---|---|---|
| `max | |S| - R |` over 17 x 17 parameter grid | — | 0 | `5.33e-15 m` absolute |
| parametric area, 24 Gauss points/span | `1809.5573684677186 m^2` | `1809.5573684677208 m^2` | `-1.22e-15` |
| arc-length-free check: volume from the closed shell | see 3.3 | `7238.229473870882 m^3` | `-7.94e-04` at `192 x 64` |

### 3.2 Quadrature floor (parametric area, per-knot-span Gauss–Legendre)

| points per span | 4 | 6 | 8 | 12 | 16 | 24 | 32 |
|---|---|---|---|---|---|---|---|
| relative error | `-7.50e-05` | `-2.92e-07` | `-1.27e-09` | `-2.78e-14` | `+1.33e-15` | `-1.22e-15` | `-3.77e-15` |

The integrand contains `sqrt(r'^2 + z'^2)`, so unlike the cylinder it is **not** a
rational function and the quadrature is not algebraically exact — it converges
spectrally to the double-precision floor. 16 points per span already sits on it; the
default 24 is kept for margin.

### 3.3 Tessellated closed shell

`V_exact = 7238.229473870882 m^3`, `A_exact = 1809.5573684677208 m^2`,
`vertices_expected = u_samples (v_samples - 1) + 2`.

| grid | triangles | vertices | collapsed quads | `V` error | `A` error | max normal | watertight | outward |
|---|---|---|---|---|---|---|---|---|
| 12 x 4 | 72 | 38 / 38 | 24 | `-1.855e-01` | `-9.791e-02` | 1.728° | yes | yes |
| 24 x 8 | 336 | 170 / 170 | 48 | `-4.965e-02` | `-2.517e-02` | 0.424° | yes | yes |
| 48 x 16 | 1440 | 722 / 722 | 96 | `-1.263e-02` | `-6.337e-03` | 0.109° | yes | yes |
| 96 x 32 | 5952 | 2978 / 2978 | 192 | `-3.172e-03` | `-1.587e-03` | 0.027° | yes | yes |
| 192 x 64 | 24192 | 12098 / 12098 | 384 | `-7.938e-04` | `-3.970e-04` | 0.007° | yes | yes |

Every doubling of the grid divides both errors by ≈4 (volume `3.74, 3.93, 3.98, 4.00`,
area `3.89, 3.97, 3.99, 4.00`), i.e. `O(h^2)`; the sign is negative throughout because the facets are inscribed. At
`48 x 16` the STL round-trip keeps the shell watertight and the volume to
`-1.263e-02` (source `-1.263e-02`), the difference being the written precision.

### 3.4 Negative controls

| control | expected | measured |
|---|---|---|
| meridian reversed (`+z` pole first) | still watertight, every face inward, negative volume | watertight `True`, volume `-6878.82 m^3`, `flipped_faces = 336 / 336`, radius error `1.78e-15 m` |
| `sphere_surface(0.0)` | `NURBSError` | raised |
| `singular_boundary_point("w")` | `NURBSError` | raised |

## 4. Figures

* `fig1_pole_shell.png` — the 48 x 16 closed shell, and the outward face normals on
  both pole rings (the ring where a naive weld would leave a hole).
* `fig2_convergence.png` — volume and area error against `1/u_samples` on log–log
  with the `O(h^2)` reference slope; and the parametric-area quadrature floor
  against the `2^-52` round-off line.
* `fig3_pole_physics.png` — left: the south-pole samples' `z` offset in ulps (the
  3-valued wobble that defeats exact-position welding); right: signed volume for the
  two pole orders — watertight in both, only the sign changes.

## 5. Pitfalls (each one cost a debugging round)

1. **The pole ring does not evaluate bit-for-bit.** `(a·b)/b` rounds twice, so 49
   samples of `S(u, v_min)` give 3 distinct points within 1 ulp. Before the fix the
   `48 x 16` shell was **not** watertight (12 inward faces, a spurious 90.4° normal):
   the sliver triangles had garbage normals. The `24 x 8` grid happened to weld, which
   is exactly how such a defect survives a spot check — *the failure depends on the
   density*. Fixed by keying singular rows on the control point.
2. **`add_quad` lied about how many faces it emitted.** It returned `2` for a pole
   quad but appended only `(a, c, d)`. `NURBSSurface.tessellate` records
   `faces_per_quad` from that return value, so every quad after the first pole row was
   shifted by one; `tessellation_report` then compared the wrong facet and finally
   raised `IndexError: tuple index out of range` (`383` indices for `336` faces). The
   cylinder anchor never noticed because it has no degenerate quad. Fixed at the
   contract level (`MeshBuilder.add_quad` returns what it appended) and pinned by
   `test_a_pole_quad_reports_the_face_it_actually_emitted` plus the
   `faces_per_quad` alignment test at five densities.
3. **The pole order decides the sign, silently.** A north→south meridian is the *same
   sphere* (radius error `1.8e-15 m`, area exact) and tessellates to the *same*
   watertight shell — with every face inward and `V < 0`. Watertightness is therefore
   not evidence of usable orientation. `sphere_surface` runs south→north and the
   reversed case is a pinned test, not a comment.
4. **`orientation_flips` is a walk-order statistic, not a winding test.** It compares
   each face with its predecessor in row-major order, so on a patch closed in `u` it
   jumps from the south-pole row to the north-pole row and reports exactly
   `u_samples - 1` "flips" (23 at `24 x 8`) with a perfectly consistent orientation.
   The topology-aware invariant — every directed edge appears once and has exactly one
   opposite twin — is what the test asserts.
5. **Inscribed facets always underestimate.** Volume and area converge from below;
   a test that only checks `|error| < tol` would accept a systematically wrong sign.
   The convergence test checks the `O(h^2)` ratio instead.
6. **Quadrature is not exact here.** Copying the cylinder's "24 points is exact to
   round-off" claim to the sphere would be wrong reasoning: the integrand is not
   rational, and the observed floor is reached at 12–16 points (§3.2), not derived.

## 6. Error judgement (the explicit verdict)

| key quantity | computed | reference | relative error | threshold | verdict |
|---|---|---|---|---|---|
| `|S(u,v)|` (17 x 17 grid, `R = 12 m`) | — | `R` | `5.33e-15 m` abs. (`4.4e-16` rel.) | `< 1e-12` | **pass** |
| parametric area (24 pts/span) | `1809.5573684677186 m^2` | `1809.5573684677208 m^2` | `-1.22e-15` | `< 1e-12` | **pass** |
| closed-shell volume at `192 x 64` | `7232.4838 m^3` | `7238.2295 m^3` | `-7.94e-04` | `< 1e-3` | **pass** |
| closed-shell area at `192 x 64` | `1808.839 m^2` | `1809.557 m^2` | `-3.97e-04` | `< 1e-2` | **pass** |
| STL round-trip volume (`48 x 16`) | `7146.79 m^3` | `7238.23 m^3` | `-1.263e-02` | `< 3e-2` | **pass** (discretisation of that grid) |
| pole weld | 2 pole vertices | 2 | 0 | exact | **pass** |
| reversed-pole orientation | `flipped_faces = 336/336`, `V < 0` | rejected as inward | — | must be rejected | **pass** |
| max face normal deviation (`96 x 32`) | `0.027°` | `0°` | — | `< 1°` | **pass** |

**Is every key quantity below 3 %?** Yes. The parametric quantities are exact to
round-off (`1e-15`); the *tessellated* quantities are discretisation errors reported
at the grid they were measured on, and they converge as `O(h^2)` with measured ratios
`3.74, 3.93, 3.98, 4.00`. At `192 x 64` the volume error is `0.079 %`. The only row
above `3e-3` is the STL round-trip measured at the coarse `48 x 16` grid (`1.26 %`),
which is a property of the grid and not of the round-trip: the same grid in memory
gives `1.2631 %` and the recovered file `1.2632 %`. No threshold was tuned to make a
row pass.

### Machine-checkable criteria of this run

`artifacts/nurbs_sphere/metrics.json` carries them under `criteria`, all `true`:

| criterion | value |
|---|---|
| `sphere_exact_within_1e-12` | pass |
| `parametric_area_within_1e-12` | pass |
| `poles_weld_to_two_vertices` | pass |
| `vertex_count_closed_form_holds` | pass (5 densities) |
| `faces_per_quad_aligned_at_every_density` | pass |
| `no_zero_area_faces_at_every_density` | pass |
| `watertight_at_every_density` | pass |
| `outward_at_every_density` | pass |
| `volume_converges_o_h2` | pass |
| `area_converges_o_h2` | pass |
| `volume_error_below_threshold_at_192x64` | pass |
| `stl_round_trip_watertight` | pass |
| `reversed_poles_rejected_as_inward` | pass |
| `bad_inputs_rejected` | pass |

## 7. Limitations (what this does **not** claim)

* **One patch is not a hull.** The sphere is the minimal closed patch, not a ship: a
  real hull needs several patches stitched along trimmed seams, which the type still
  does not support (`docs/nurbs-surface.md` §7 keeps that limitation).
* **Only boundary rows are checked for singularity.** A control net that pinches on an
  interior row (a degenerate middle ring) is not detected; the sphere's poles are both
  boundaries. This is a named gap, not an oversight: interior degeneracy would need a
  different test (the row's *evaluated* curve collapsing, not its control polygon).
* **Welding is exact, keyed on the control point.** If a caller builds an equivalent
  sphere by hand with non-coincident pole control points (`(0, 0, -R + 1e-9)`), the
  rows are not singular, the welding does not happen, and the shell is not watertight —
  correctly so, because such a net *is* a slightly different surface.
* **Volume accuracy is bounded by the grid, not by the quadrature.** `192 x 64` is
  `24192` triangles for `0.08 %`; a finer grid is a wall-clock decision, not a
  correctness one.

## 8. Reproduce

```bash
python3 examples/nurbs_sphere/run_example.py      # writes artifacts/nurbs_sphere/*
python3 -m pytest tests/test_nurbs_pole_surface.py -q   # 27 tests, 2 rejection tests
python3 -m pytest tests/test_nurbs.py tests/test_mesh.py -q   # regression, unchanged
```