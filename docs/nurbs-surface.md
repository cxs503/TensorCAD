# NURBS parametric surface service: exact arcs, exact area, verified tessellation

**Mission** IN16-K5 / TensorCAD. **Module** `tensorcad/nurbs.py` (standard library
only). **Artifacts** `artifacts/nurbs_surface/` — `metrics.json`,
`pc6_nurbs_block.json`, `cylinder_solid.stl`, `fig1..fig4`. **Schema** of the
parametric-surface block `tensorcad.nurbs.surface/1.0`, additive to P1's
`tensorcad.pc6.geom/1.0`.

## 1. Case: why a parametric surface, and what it has to be good for

The geometry service could already publish a *sampled* hull, an *implicit* field and a
*triangle mesh*. None of those three carries the surface's own analytic structure: a
circle sampled at 96 points is a 96-gon whose area is wrong by 1e-3, a tessellated
patch has no derivative at the facet edges, and "the surface is smooth" is an opinion
about the sampling density rather than a property of the data. A parametric (NURBS)
surface is the fourth representation, and it is the one that can be *audited*:

| representation | area | derivatives | continuity |
|---|---|---|---|
| sampled mesh | integral over facets (discretisation error) | piecewise constant normal | not applicable |
| implicit field | needs a re-meshing step | gradient of the field | field-dependent |
| **NURBS patch** | integral of `\|S_u x S_v\|` (exact/quadrature) | exact `S^(k,l)` | property of the knot vector |

The minimal surface the campaign asks for is a **rational quadratic arc revolved into a
cylinder** of radius 12 m and height 11 m — the scale of the PC6 hull's bilge region
(R = B/2 = 12 m, H = T = 11 m) — with three things proven about it:

1. the area equals `2 pi r h` to round-off (analytic anchor, criterion `< 1e-6`);
2. the parameter derivatives are continuous across the interior knots to the order the
   data promises, *measured* from both sides;
3. the tessellation used for STL export is watertight, every face points outward, and
   its normals converge to the analytic ones.

## 2. Method

**Representation.** Rational tensor-product B-spline with a rectangular control net.
The circle is three exact 120-degree rational quadratic arcs (7 control points,
weights `1, 1/2, 1, 1/2, 1, 1/2, 1`, interior knot multiplicity 2); the cylinder is that
circle times a linear direction, and because both `v`-rows share the same weights the
quotient leaves `z` untouched — the surface *is* the cylinder, it is not fitted to it.

**Evaluation.** Basis functions by Cox–de Boor (Piegl & Tiller A2.2); derivatives by
forming the **derivative control net** (A3.4) and then unwinding the rational quotient
with the standard recursion (A3.6). Building the derivative net instead of hand-rolling
a basis-derivative table is a deliberate choice: this module's first implementation used
a transcription of A2.3 and produced `|C'| = 0` at the arc endpoints, i.e. an arc length
of 0.689 instead of 2.0. The derivative net has no such failure mode, and the test suite
checks it against central differences of the evaluations.

**Area.** Gauss–Legendre quadrature applied *per knot span* in both directions on
`|S_u x S_v|`. A 24-point rule per span is exact to round-off for these integrands
(they are rational functions of low degree after the substitution), which is why the
parametric area matches the closed form at 1e-15 rather than at 1e-6.

**Continuity.** At each interior knot the derivative is evaluated twice, once from the
span on the left and once from the span on the right; the difference is the jump. The
B-spline bookkeeping declares `p - m` (degree minus multiplicity); the measurement says
what is actually delivered. A derivative of order `p - m + 1` does not exist as a
function across the knot (it is a distribution), so the measurement is clamped to what
is representable and the clamp is reported (`max_measurable_order`).

**Tessellation and export.** `tessellate(nu, nv)` samples the parameter domain, emits
quads wound so that their normal is `S_u x S_v`, and shares vertices by *exact*
position. The revolution seam therefore welds itself (`S(0, v) == S(1, v)` to the last
bit), so watertightness is a property of the data rather than of a welding tolerance.
`cylinder_solid_mesh` adds two triangle-fan caps whose rim vertices *are* the lateral
vertices, giving a closed solid that `tensorcad.export.write_stl` writes and
`read_stl` recovers as a closed shell.

**Toolchain.** Python 3.11, standard library only (`math`, `dataclasses`); matplotlib
is used *only* by the example script to draw the figures and is not imported by the
package (the stdlib-only promise is tested in `tests/test_public_surface.py`).

## 3. Results

### 3.1 Analytic anchors

| quantity | computed | analytic | relative error |
|---|---|---|---|
| quarter-circle arc, r = 2 m, sweep 1 rad | 2.000000000000 m | 2.000000000000 m | 0 |
| full circle, r = 12 m (3 rational arcs) | 75.398223686 m | 75.398223686 m | -1.88e-16 |
| max radius deviation, arc | 4.44e-16 m | 0 m | round-off |
| max radius deviation, circle | 5.33e-15 m | 0 m | round-off |

| radius [m] | height [m] | parametric area [m^2] | 2 pi r h [m^2] | relative error |
|---|---|---|---|---|
| 12 | 11 | 829.380460548 | 829.380460548 | +3.02e-15 |
| 1 | 1 | 6.283185307 | 6.283185307 | +7.07e-16 |
| 0.5 | 40 | 125.663706144 | 125.663706144 | +1.47e-15 |

The lateral area error is `+3.02e-15` at
R = 12 m, H = 11 m, i.e. five orders inside
the `< 1e-6` criterion. It is *quadrature* error only: the geometry is exact, as
figure 4 shows (deviation from `sqrt(x^2 + y^2) - r` is round-off everywhere).

### 3.2 Quadrature convergence

| GL order per knot span | area [m^2] | relative error |
|---|---|---|
| 2 | 823.070543756731 | +7.61e-03 |
| 3 | 829.844786915221 | +5.60e-04 |
| 4 | 829.346680800628 | +4.07e-05 |
| 6 | 829.380283870709 | +2.13e-07 |
| 8 | 829.380459629879 | +1.11e-09 |
| 12 | 829.380460547681 | +3.21e-14 |
| 16 | 829.380460547705 | +3.29e-15 |
| 24 | 829.380460547708 | +1.37e-16 |

### 3.3 Tessellation: topology, normals, area

| grid | triangles | boundary edges | area error | max normal angle [deg] | inward faces | orientation flips |
|---|---|---|---|---|---|---|
| 12 x 8 | 192 | 24 | -1.16e-02 | 0.7243 | 0 | 0 |
| 24 x 8 | 384 | 48 | -2.91e-03 | 0.1909 | 0 | 0 |
| 48 x 8 | 768 | 96 | -7.29e-04 | 0.0483 | 0 | 0 |
| 96 x 8 | 1536 | 192 | -1.82e-04 | 0.0121 | 0 | 0 |
| 192 x 8 | 3072 | 384 | -4.56e-05 | 0.0030 | 0 | 0 |

The lateral surface is an **open tube**: it has exactly `2 nu` boundary edges (the two
rims) and no others, which is the statement that the seam welded. All faces point
outward (`inward = 0`), the winding is globally consistent (`flips = 0`), and the
maximum angle between a facet normal and the analytic normal at the quad's parametric
centre falls by a factor ~4 per doubling of the grid (O(h^2)):
0.7243 deg at
12 samples to
0.0030 deg at
192.

### 3.4 Closed solid and STL round-trip

Reference volume `pi r^2 h = 4976.282763 m^3`.

| u samples | triangles | volume [m^3] | volume error | area error | watertight |
|---|---|---|---|---|---|
| 96 | 1344 | 4972.6545 | -7.29e-04 | -4.68e-04 | yes |
| 192 | 2688 | 4975.3755 | -1.82e-04 | -1.17e-04 | yes |
| 384 | 5376 | 4976.0559 | -4.56e-05 | -2.92e-05 | yes |

STL round-trip (binary float32, re-welded by exact position):
watertight **True**, one shell, volume error `-4.56e-05`.
The remaining error is the *inscribed* cap polygon, not the format: it halves by four
per doubling, exactly as the table above shows.

### 3.5 Parameter-derivative continuity

| curve / surface | interior knot | multiplicity | declared | measured | first-derivative jump |
|---|---|---|---|---|---|
| three-arc circle, r = 12 m | 0.333333 | 2 | C0 | C1 | +1.59e-14 |
| three-arc circle, r = 12 m | 0.666667 | 2 | C0 | C1 | +3.78e-14 |
| cylinder surface (u direction) | 0.333333 | 2 | C0 | C1 | +1.59e-14 |
| cubic control curve | 0.4 | 1 | C2 | C2 | 0 |

The interesting result is the first two rows: the knot bookkeeping only promises `C0`
at the arc joins, yet the parameter derivative is continuous to 4e-14. The reason is
symmetry — each rational quadratic arc has the same parametric speed at both of its
ends, so the join matches in value *and* in slope. The **positive control** is the
cubic row: there the third derivative jumps by
`495.117` (against a scale of `342.155`), so the
measurement does detect a genuine discontinuity; it is not a test that always passes.

### 3.6 Negative controls (the ticket's two rejection criteria, plus four more)

| control | expected rejection reason | outcome | message |
|---|---|---|---|
| `non_monotone_knots` | knot vector must be non-decreasing | rejected | knot vector must be non-decreasing: knots[4]=0.4 < knots[3]=0.6 |
| `zero_weight` | weights must be strictly positive | rejected | weights[1] = 0.0: rational B-spline weights must be strictly positive (zero divides by z... |
| `negative_weight` | weights must be strictly positive | rejected | weights[0] = -0.5: rational B-spline weights must be strictly positive (zero divides by ... |
| `parameter_outside_domain` | parameter outside [U[p], U[n+1]] is an error, not a silent clamp | rejected | parameter 1.5 is outside the domain [0.0, 1.0] |
| `unclamped_end` | a non-clamped end would miss its end control point | rejected | knot vector is not clamped at the start: knots[0]=-0.5 != knots[p]=0.0 |
| `knot_multiplicity_too_high` | multiplicity above degree + 1 | rejected | a degree-2 spline with 3 control points needs 6 knots, got 11 |

### 3.7 PC6 contract compatibility

P1's geometry contract (`tensorcad.pc6.geom/1.0`, `tensorcad/pc6_contract.py` in
`main`) is extended **additively**: the block carries the same `schema` tag, adds one
key (`nurbs_surface`), and `tensorcad.pc6_contract.load_contract` reads only
`schema` / `parameters` / `artifacts` / `geometry` / `acceptance`, so no existing check
is modified or relaxed. The test suite proves the additivity by re-parsing a contract
with the block appended through P1's **own** `_parameters_from_contract` and asserting
the parameters are unchanged.

Geometry-first gate (recomputed from P1's closed form, never copied):

| quantity | value |
|---|---|
| length L [m] | 200 |
| beam B [m] | 24 |
| draft T [m] | 11 |
| displacement volume [m^3] | 41629.884 |
| block coefficient C_b | 0.788445 |
| prismatic coefficient C_p | 0.875063 |
| midship coefficient C_m | 0.901014 |

A declared value that deviates by more than 1e-3 is **refused** (tampered C_b at -10 %
=> `True`, message
"declared block_coefficient = 0.709600300367423 deviates 1.000e-01 from the recomputed 0.78844477"), because a surface block whose
stated C_b disagrees with its own parameters is how two solvers end up with two hulls.

## 4. Figures

* `fig1_cylinder_normals.png` — the tessellated cylinder with its outward face normals.
* `fig2_convergence.png` — parametric area versus quadrature order (round-off), and
  tessellated area versus `1/nu` with the `O(h^2)` reference slope.
* `fig3_continuity.png` — one-sided derivative jumps per order for the arc joins and
  for the cubic control curve, plus `|C'(u)|` of the three-arc circle showing the
  symmetry that makes the joins `C1`.
* `fig4_exactness.png` — `sqrt(x^2 + y^2) - r` over the whole parameter domain: the
  arc is exact, the remaining error is round-off.

## 5. Pitfalls found (and what they cost)

1. **A transcribed derivative table is a trap.** The first implementation of the basis
   derivatives was a literal transcription of Piegl & Tiller A2.3; it returned `-1.5`
   instead of `-1.0` for `N'_0(0.5)` and zero at the endpoints, giving an arc length of
   0.689 m for a 2 m arc. Replacing it with the derivative **control net** (A3.4) fixed
   it, and the central-difference test keeps it fixed.
2. **`face_normal` returns an unnormalised cross product.** Comparing it with a unit
   reference normal produced a spurious 24.9 degree "normal error" that looked like a
   geometry bug. Normalising the facet normal before the dot product turned the error
   into the expected O(h^2) decay.
3. **Quad index versus vertex stride.** `Tessellation.faces_of_quad` maps a quad to its
   triangles with the *quad* stride (`nv`), while the vertex grid has `nv + 1` columns;
   mixing the two silently compares the wrong facet. The confusion cost one debugging
   round because the resulting angle error looked like geometry.
4. **A C0 knot is not a differentiable point.** Asking for the second derivative of the
   three-arc circle is asking for a distribution. It now raises
   `NURBSError("zero-length knot interval ...")` instead of returning a number, and the
   continuity measurement clamps itself to the representable order.
5. **A coarse grid hides two different errors.** At `nu = 12` the tessellated area is
   off by -1.2 % and the normal by 0.72 deg; both are discretisation, while the
   parametric area is already exact. Reporting only "the area error" would conflate them,
   so the report always states *which* area it means.

## 6. Error judgement (the explicit verdict)

| key quantity | computed | reference | relative error | threshold | verdict |
|---|---|---|---|---|---|
| cylinder lateral area (R=12, H=11) | 829.380460548 m^2 | 829.380460548 m^2 | +3.02e-15 | < 1e-6 | **pass, 1e-9 x margin** |
| arc length (r = 2 m, 1 rad) | 2.000000000000 m | 2 m | 0 | < 1e-6 | **pass** |
| full circle length (r = 12 m) | 75.398223686 m | 75.398223686 m | -1.88e-16 | < 1e-6 | **pass** |
| circle first-derivative jump at the joins | 3.777e-14 m | 0 (C1) | rel. to \|C'\| = 62.354: 6.06e-16 | < 1e-9 | **pass** |
| closed-solid volume (u = 384) | 4976.0559 m^3 | 4976.2828 m^3 | -4.56e-05 | < 1e-3 | **pass (inscribed facets)** |
| STL round-trip volume | 4976.0559 m^3 | 4976.2828 m^3 | -4.56e-05 | < 1e-3 | **pass** |
| tessellation max normal angle | 0.0483 deg | 0 deg | — | < 1 deg | **pass** |

**Is every key quantity below 3 %?** Yes, with the caveat stated per row: the
parametric quantities are exact to round-off (1e-15), and the tessellated quantities
are *discretisation* errors that converge as `O(h^2)` and are reported at the grid they
were measured on. No quantity is above 3 %, and none of the pass margins depends on a
tolerance chosen to make it pass.

### Machine-checkable criteria of this run

| criterion | value |
|---|---|
| `all_negative_controls_rejected` | **pass** |
| `continuity_declared_holds` | **pass** |
| `cylinder_area_error_below_1e-6` | **pass** |
| `pc6_gate_passed` | **pass** |
| `pc6_tamper_rejected` | **pass** |
| `solid_volume_error_below_1e-3` | **pass** |
| `solid_watertight` | **pass** |
| `stl_round_trip_watertight` | **pass** |
| `tessellation_max_normal_angle_deg_below_1` | **pass** |
| `tessellation_no_inward_faces` | **pass** |

## 7. Limitations, and what is deliberately not done

* **No trimmed or multi-patch surface.** The type is one rectangular tensor-product
  patch. A bulbous bow or a reentrant stern needs trimming, which is out of scope; the
  PC6 hull is still lofted by `tensorcad.pc6` from analytic stations.
* **The PC6 block carries the representation and the gate, not a fitted hull surface.**
  The cylinder is the *demonstration surface* of the representation format and the
  L/B/Cb gate; re-fitting the PC6 body as one NURBS patch is a separate piece of work
  and is not claimed here.
* **Derivatives above the declared continuity order raise** rather than return a
  distribution (documented in the module docstring and pinned by a test).
* **Watertightness is claimed for the closed solid, not for the lateral tube.** The
  tube is open at its two rims by construction, and the report says so explicitly
  instead of quietly testing only the closed case.

## 8. Reproduce

```bash
python3 examples/nurbs_surface/run_example.py     # writes artifacts/nurbs_surface/*
python3 examples/nurbs_surface/build_doc.py       # regenerates this document
python3 -m pytest tests/test_nurbs.py -q          # 26 tests, 6 of them rejection tests
```

Audit chain: `tensorcad/nurbs.py` -> `examples/nurbs_surface/run_example.py` ->
`artifacts/nurbs_surface/metrics.json` -> this document -> `tests/test_nurbs.py`.
Every number above is printed by the first step and re-asserted by the last one.
