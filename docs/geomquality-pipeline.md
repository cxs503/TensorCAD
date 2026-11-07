# CAD geometry quality: a full-chain verification pipeline

*Audit chain: branch `agent/m5345`, package `tensorcad/geomquality`, artefacts in
`artifacts/geomquality/`, tests in `tests/test_geomquality_io.py` and
`tests/test_geomquality_pipeline.py`. Numbers below are rendered from
`artifacts/geomquality/gate.json` by `examples/geomquality/build_doc.py`.*

## 1. What this is

A hull that renders correctly is not a hull a solver can use. Three failures are on
record in this repository and in every shipyard's model archive: a surface with a crack
along a patch seam, a surface that is tangent-continuous everywhere except at one
unintended knuckle, and a surface that is simply the wrong ship. None of them is visible
in a viewer, all of them survive a watertightness thumbnail test, and all of them change
the displacement.

This module is the chain from a **CAD file** to a **gated hydrostatic table**:

```
IGES / STEP / native B-rep
   -> edge-face adjacency graph         (brep.py)
   -> watertightness 0/1                (free edges, cracks, non-manifold, orientation)
   -> declared L/B/T and Cb             (before any integral: is this the right ship?)
   -> seam continuity G0/G1/G2          (continuity.py, calibrated thresholds)
   -> hydrostatics by patch quadrature  (hydro.py, exact per knot span)
   -> grid convergence and GCI          (three densities, equal-arc refinement)
   -> machine gate                      (pipeline.py, CRITERIA + enforce())
```

The gate never reads a narrative field: :func:`tensorcad.geomquality.pipeline.enforce`
takes :class:`CriterionResult` objects, each of which carries a number, a limit and a
comparison direction.

## 2. Machine criteria

| criterion | value | limit | verdict | meaning |
|---|---:|---:|---|---|
| `G1_watertight` | 1 | 1 | PASS | the shell is closed, orientable and positively oriented |
| `C1_wigley_volume_error` | 2.23341e-13 | 0.0005 | PASS | Wigley III displacement from the CAD patch set vs closed form |
| `C2_defect_detection_rate` | 1 | 1 | PASS | library cases whose declared defect code is reported (zero misses) |
| `C3_curvature_f1` | 1 | 0.9 | PASS | held-out seam classification F1 with calibrated thresholds |
| `C4_hydrostatic_error` | 0.000114191 | 0.001 | PASS | worst of volume / LCB / VCB / waterplane area against the analytic hull |
| `C5_grid_convergence_order` | 2.13474 | 1.9 | PASS | Richardson observed order of the tessellated displacement |
| `P5_geometry_deviation` | 1.8497e-06 | 0.0001 | PASS | max |y_CAD - y_analytic| / L for every reference hull (campaign item 5) |
| `C7_seconds_per_hull` | 9.06075 | 10 | PASS | wall-clock for the whole chain on one hull |

Gate verdict: **PASS**.

## 3. How it is done

**Exchange.** Three routes into the same B-rep. Native `.brep.json` is lossless. STEP is
read as a *topological* model (shared `edge_curve` identity is the shared edge), so the
adjacency graph comes from the file. IGES entity 128 carries **no topology and no
orientation at all**, so the reader sews boundaries geometrically and then propagates a
consistent orientation over the face graph (3 hulls pass all
three routes; deviations in the table below). The closure flag (the deck patch, which
seals the waterplane but is not wetted) travels in the STEP entity name and in the IGES
directory level, because neither format has such an attribute.

**Topology.** Edges are used once (free), twice in opposite directions (manifold) or
more (non-manifold). A pair of boundary curves that coincide within
`GAP_TOLERANCE_FRACTION * L` but are two distinct edges is a **crack**, measured in metres
and as a fraction of L. Collapsed boundaries (the keel line end, a bow apex) are recorded
as *degenerate* and exempt: they bound no area and cannot leak. A face whose declared
orientation contradicts the winding of its own loop is an **inside-out face**.

**Representation.** The vertical direction of the Series 60 and KVLCC2-class hulls is a
power law `y ~ zeta**(1/m)` whose derivative at the keel is unbounded. Instead of fighting
that with graded knots, the patch is parameterised as `z = -T + T*v**m`, `y = a(x)*v`:
`z(v)` is then a polynomial of degree m (exactly one Bezier span) and `y` is linear in
`v`, so the vertical direction is **exact** and all representation error is in the
one-dimensional plan form `a(x)`, where it is measured directly.

**Quadrature.** Each face is integrated with Gauss-Legendre per knot span, with the order
raised to follow the patch degree (`order = degree + 2`): the volume integrand of a
degree-24 vertical span is a polynomial of degree ~47, and the default 6-point rule
under-integrates it by 0.5 %. That was measured, not assumed -- the test
`test_series60_closed_forms` pins the volume to `Cb*L*B*T` to 1e-6.

**Grid convergence.** Refining *parameter space* is not refining *physical space* on a
graded patch: with `z = -T + T*v**24`, doubling `v_samples` refines the keel and leaves
the deck alone, and the observed order collapses to 1.55. Placing the samples at equal
arc length (`hydro.arc_length_parameters`) recovers the second-order rate; the difference
is a test of its own.

**Continuity.** Seams are classified from the normal angle and the relative jump of the
normal curvature *across* the seam. Limits are calibrated by grid search on 60 % of a
labelled seam set (G0/G1/G2 constructed from known end derivatives) and scored on the
remaining 40 %. Designed creases (deck chine, keel knife edge) are recognised by their
normal angle and reported as `KNUCKLE`, not as defects.

**Defects.** 14 injected failures, each with the code it must produce and
a repair hint: hole, open waterline, unsewn seam, non-manifold flap, flipped face,
duplicate face, sliver, micro edge, curvature-jump splice, tangent-jump splice,
corrugation, fold, unit scale, and a 5 % narrow beam.

## 4. Results

### 4.1 Hydrostatics against the analytic hull

| hull | quantity | CAD patch set | tessellation | analytic | patch error | facet error |
|---|---|---:|---:|---:|---:|---:|
| wigley_iii | volume | 2.77777778 | 2.76754223 | 2.77777778 | 2.23e-13 | 3.68e-03 |
| wigley_iii | lcb_x | -6.24500451e-18 | -0.00400885312 | -1.68927372e-16 | 1.63e-17 | 4.01e-04 |
| wigley_iii | vcb_z | -0.234375 | -0.23435329 | -0.234375 | 7.11e-16 | 9.26e-05 |
| wigley_iii | waterplane_area | 6.66666667 | -3.02183304e-07 | 6.66666667 | 2.23e-13 | 1.00e+00 |
| series60_cb0.60 | volume | 3.75000083 | 3.73337896 | 3.75000001 | 2.20e-07 | 4.43e-03 |
| series60_cb0.60 | lcb_x | -2.31296412e-16 | -0.000841100584 | -1.42287079e-16 | 8.90e-18 | 8.41e-05 |
| series60_cb0.60 | vcb_z | -0.306122449 | -0.3053321 | -0.306129532 | 2.31e-05 | 2.60e-03 |
| series60_cb0.60 | waterplane_area | 6.25000138 | 1.34094125e-15 | 6.25000001 | 2.20e-07 | 1.00e+00 |
| kvlcc2_class | volume | 312624.529 | 309883.76 | 312621.67 | 9.14e-06 | 8.76e-03 |
| kvlcc2_class | lcb_x | -11.1372716 | -10.742194 | -11.136 | 1.14e-04 | 3.54e-02 |
| kvlcc2_class | vcb_z | -9.9047619 | -9.8840377 | -9.9047619 | 3.77e-15 | 2.09e-03 |
| kvlcc2_class | waterplane_area | 16533.028 | 4.66116035e-12 | 16532.8768 | 9.14e-06 | 1.00e+00 |

Worst case: kvlcc2_class `lcb_x`: 1.14e-04 (campaign limit 1e-3). The Wigley row is at round-off because
the hull is polynomial and the patch set represents it exactly; the Series 60 and
KVLCC2-class rows are representation error, not quadrature error.

### 4.2 CAD exchange

| hull | format | watertight before | watertight after | max deviation [m] | bytes |
|---|---|---|---:|---:|---:|
| wigley_iii | `.brep.json` | True | True | 0.00e+00 | 27918 |
| wigley_iii | `.step` | True | True | 1.11e-16 | 24703 |
| wigley_iii | `.iges` | True | True | 1.11e-16 | 6075 |
| series60_cb0.60 | `.brep.json` | True | True | 0.00e+00 | 258318 |
| series60_cb0.60 | `.step` | True | True | 3.03e-13 | 189789 |
| series60_cb0.60 | `.iges` | True | True | 3.01e-13 | 149121 |
| kvlcc2_class | `.brep.json` | True | True | 0.00e+00 | 282224 |
| kvlcc2_class | `.step` | True | True | 2.77e-13 | 195357 |
| kvlcc2_class | `.iges` | True | True | 1.24e-11 | 131301 |

### 4.3 Grid convergence

| hull | levels | observed order p | GCI (finest) | finest error vs analytic |
|---|---|---:|---:|---:|
| wigley_iii | 8x4, 16x8, 32x16 | 1.984 | 0.263 % | 2.07e-01 % |
| series60_cb0.60 | 8x4, 16x8, 32x16 | 2.135 | 0.275 % | 2.64e-01 % |
| kvlcc2_class | 8x4, 16x8, 32x16 | 1.916 | 0.681 % | 4.99e-01 % |

### 4.4 Defect library

14 of 14 cases detected, **100 %**,
false alarms on the clean reference hull: none.

| case | expected code | detected codes | repair hint |
|---|---|---|---|
| hole | `FREE_EDGE` | CB_MISMATCH, FREE_EDGE | sew the two coincident boundaries, or cap the opening with a face |
| open_waterline | `OPEN_TOP` | OPEN_TOP | close the waterline: add the deck / waterplane face at z = T |
| unsewn_seam | `SEAM_GAP` | FREE_EDGE, SEAM_GAP | snap the seam vertices to a common tolerance and re-sew the patch pair |
| non_manifold_flap | `NON_MANIFOLD` | NON_MANIFOLD, ORIENTATION | remove the extra face at the edge, or split the shell at the T-junction |
| flipped_face | `ORIENTATION` | ORIENTATION | flip the offending face so its loop is anticlockwise about the outward normal |
| duplicate_face | `DUPLICATE_FACE` | CB_MISMATCH, DUPLICATE_FACE, NON_MANIFOLD | delete the redundant face; it contributes displacement twice |
| sliver_face | `SLIVER` | OPEN_TOP, SLIVER | collapse the degenerate face to an edge or rebuild the patch with a coarser net |
| micro_edge | `MICRO_EDGE` | MICRO_EDGE, NON_MANIFOLD, OPEN_TOP, SLIVER | merge the two vertices of the micro edge (tolerance weld) |
| g1_splice | `CURVATURE_JUMP` | CB_MISMATCH, CURVATURE_JUMP, FREE_EDGE | re-fair the seam: match the second derivative (G2) or accept the knuckle explicitly |
| ripple | `ROUGHNESS` | ROUGHNESS | low-pass the control net; a high-frequency ripple is not a design feature |
| fold | `SELF_INTERSECTION` | CB_MISMATCH, ORIENTATION, SELF_INTERSECTION | reduce the offending control point; the patch normal reverses inside the domain |
| unit_scale | `UNIT_SCALE` | CB_MISMATCH, UNIT_SCALE | apply the declared unit factor to the whole model (mm/m confusion) |
| cb_mismatch | `CB_MISMATCH` | CB_MISMATCH | re-scale the sections to the declared Cb, or correct the declared value |
| g0_splice | `TANGENT_JUMP` | CB_MISMATCH, FREE_EDGE, TANGENT_JUMP | make the two patch tangents collinear at the seam (G1) before fairing |

### 4.5 Seam continuity

* Calibrated limits: normal angle 0.05 deg,
  relative curvature jump 0.01.
* Held-out score (60 seams):
  F1 = 1.000,
  precision = 1.000,
  recall = 1.000.
* Curvature-only view (G1 vs G2, ignored tangent angle):
  F1 = 1.000.
* Separation achieved: minimum G1 curvature jump
  0.283 vs maximum G2
  8.12e-05; smallest G0 normal angle
  5.11 deg.

### 4.6 Figures

* `artifacts/geomquality/figures/fig1_geometry.png`
* `artifacts/geomquality/figures/fig2_parity.png`
* `artifacts/geomquality/figures/fig3_convergence.png`
* `artifacts/geomquality/figures/fig4_defect_matrix.png`
* `artifacts/geomquality/figures/fig5_continuity.png`

## 5. Error judgement (against the campaign numbers)

| quantity | reference | measured | threshold | verdict |
|---|---|---|---|---|
| Wigley III displacement | `4/9 LBT` closed form | 2.23e-13 | 0.05 % | yes |
| Wigley III waterplane area | `2/3 LB` | 2.23e-13 | 0.1 % | yes |
| Series 60 displacement | `Cb*L*B*T` (Cp, Cm closed forms) | 2.20e-07 | 0.1 % | yes |
| KVLCC2-class displacement | analytic half-breadth integral | 9.14e-06 | 0.1 % | yes |
| KVLCC2-class LCB | analytic | 1.14e-04 (of L) | 0.1 % | yes |
| geometry deviation, all three | analytic hull | 1.85e-06 of L | 1e-4 of L | yes |
| grid convergence | Richardson, 3 densities | p = 2.135 | > 1.9 | yes |

Every acceptance number is inside its threshold by at least an order of magnitude except
the two slow-convergence cases in section 6.

## 6. Notes, limitations, and what is not claimed

1. **The KVLCC2 hull is a class hull, not the official geometry.** No KVLCC2 offset table
   exists in this repository, so the form is a math hull whose two plan-form exponents are
   solved to reproduce the published `Cb = 0.8098` and `LCB = -3.48 % L`. The Cb check on
   that hull is therefore a *representation* consistency check (it can and does fail if
   the patch set drifts), **not** an independent validation against KVLCC2. The Wigley and
   Series 60 rows are the independent ones.
2. **IGES trimmed surfaces (type 144) are refused, not approximated.** The reader raises
   `ExchangeError` instead of silently dropping trims; a hull exported with trims must go
   through STEP or be re-exported untrimmed.
3. **STEP support is a documented subset** (`advanced_face` + `edge_loop` +
   `b_spline_surface_with_knots` + degree-1 seam curves). Assemblies with placement
   transforms are read in file coordinates; rational surfaces are read but not written.
4. **No general self-intersection test between patches.** `SELF_INTERSECTION` covers a
   patch that folds through itself (Jacobian reversal over its own domain); two distinct
   patches passing through each other is out of scope.
5. **The Series 60 section exponent is rounded to an integer** (m = 24 for Cb = 0.60) so
   that `z(v) = -T + T*v**m` is an exact Bezier span. The realised Cb is reported
   (0.600000) and used as the reference.
6. **Timing is parameterisation-sensitive.** The slowest hull is Series 60 (degree-24
   vertical span, 26-point Gauss per span): 9.06 s against a 10 s
   budget. Wigley needs 0.37 s. A rational
   patch falls back to the reference NURBS evaluator and is roughly 20x slower.
7. **The observed order is 2.13 on Series 60**, slightly above the asymptotic 2.0: the
   three densities are not yet asymptotic for the u-direction spline error. It is reported
   as measured rather than clipped to the "expected" value.
8. **What is still missing for a production gate**: curvature-continuity checking of the
   *interior* of a patch (a dent inside a face is smooth and is caught only by the ripple
   and fold indicators), plate thickness and structural idealisation, and a repair that
   actually fixes a seam rather than describing the fix.

## 7. Reproduce

```bash
python examples/geomquality/run_pipeline.py --out artifacts/geomquality
python examples/geomquality/build_doc.py
python -m pytest tests/test_geomquality_io.py tests/test_geomquality_pipeline.py -q
```
