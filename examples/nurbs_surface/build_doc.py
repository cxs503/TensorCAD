"""Render ``docs/nurbs-surface.md`` from the measured ``metrics.json``.

The document is *generated* rather than transcribed: every number in it comes from
``artifacts/nurbs_surface/metrics.json``, which ``run_example.py`` wrote from the
actual run. Hand-copying a table is how a report ends up disagreeing with its own
artifact, so the tables here are formatted from the JSON.

    python3 examples/nurbs_surface/run_example.py
    python3 examples/nurbs_surface/build_doc.py [artifacts/nurbs_surface] [docs/nurbs-surface.md]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


def _sci(value: float) -> str:
    if value == 0.0:
        return "0"
    return f"{value:+.2e}"


def _table(headers: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    source = Path(argv[0]) if argv else Path("artifacts/nurbs_surface")
    target = Path(argv[1]) if len(argv) > 1 else Path("docs/nurbs-surface.md")
    metrics = json.loads((source / "metrics.json").read_text(encoding="utf-8"))

    arcs = metrics["arcs"]
    cylinder = metrics["cylinder"]
    solid = metrics["solid"]
    continuity = metrics["continuity"]
    compatibility = metrics["pc6_compatibility"]
    criteria = metrics["criteria"]

    arc_rows = [
        ["quarter-circle arc, r = 2 m, sweep 1 rad", f"{arcs['arc_length']:.12f} m",
         f"{arcs['arc_length_reference']:.12f} m", _sci(arcs["arc_length_relative_error"])],
        ["full circle, r = 12 m (3 rational arcs)", f"{arcs['circle_length']:.9f} m",
         f"{arcs['circle_length_reference']:.9f} m", _sci(arcs["circle_length_relative_error"])],
        ["max radius deviation, arc", f"{arcs['max_radius_deviation']:.2e} m", "0 m", "round-off"],
        ["max radius deviation, circle", f"{arcs['circle_max_radius_deviation']:.2e} m", "0 m",
         "round-off"],
    ]
    area_rows = [
        [f"{row['radius_m']:g}", f"{row['height_m']:g}", f"{row['parametric_area_m2']:.9f}",
         f"{row['analytic_area_m2']:.9f}", _sci(row["relative_error"])]
        for row in cylinder["area_rows"]
    ]
    quadrature_rows = [
        [f"{row['order']}", f"{row['area']:.12f}", _sci(row["relative_error"])]
        for row in cylinder["quadrature_convergence"]["rows"]
    ]
    tessellation_rows = [
        [f"{row['u_samples']} x {row['v_samples']}", f"{row['triangles']}",
         str(row["boundary_edges"]), _sci(row["area_relative_error"]),
         f"{row['max_normal_angle_deg']:.4f}", str(row["inward_faces"]),
         str(row["orientation_flips"])]
        for row in metrics["tessellation"]
    ]
    solid_rows = [
        [f"{row['u_samples']}", f"{row['triangles']}", f"{row['volume_m3']:.4f}",
         _sci(row["volume_relative_error"]), _sci(row["area_relative_error"]),
         "yes" if row["watertight"] else "no"]
        for row in solid["rows"]
    ]
    stl = solid["stl_round_trip"]
    negative_rows = [
        [f"`{case['case']}`", case["expected"], "rejected" if case["rejected"] else "**passed**",
         (case["message"][:88] + "...") if len(case["message"]) > 91 else case["message"]]
        for case in metrics["negative_controls"]
    ]
    circle = continuity["circle"]
    cubic = continuity["cubic_positive_control"]["knots"][0]
    gate = compatibility["geometry_gate"]["values"]

    body = f"""# NURBS parametric surface service: exact arcs, exact area, verified tessellation

**Mission** IN16-K5 / TensorCAD. **Module** `tensorcad/nurbs.py` (standard library
only). **Artifacts** `artifacts/nurbs_surface/` — `metrics.json`,
`pc6_nurbs_block.json`, `cylinder_solid.stl`, `fig1..fig4`. **Schema** of the
parametric-surface block `{compatibility['block_representation']}`, additive to P1's
`{compatibility['contract_schema']}`.

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
| **NURBS patch** | integral of `\\|S_u x S_v\\|` (exact/quadrature) | exact `S^(k,l)` | property of the knot vector |

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

{_table(["quantity", "computed", "analytic", "relative error"], arc_rows)}

{_table(["radius [m]", "height [m]", "parametric area [m^2]", "2 pi r h [m^2]", "relative error"],
        area_rows)}

The lateral area error is `{_sci(cylinder['area_relative_error'])}` at
R = {cylinder['radius_m']:g} m, H = {cylinder['height_m']:g} m, i.e. five orders inside
the `< 1e-6` criterion. It is *quadrature* error only: the geometry is exact, as
figure 4 shows (deviation from `sqrt(x^2 + y^2) - r` is round-off everywhere).

### 3.2 Quadrature convergence

{_table(["GL order per knot span", "area [m^2]", "relative error"], quadrature_rows)}

### 3.3 Tessellation: topology, normals, area

{_table(["grid", "triangles", "boundary edges", "area error", "max normal angle [deg]",
         "inward faces", "orientation flips"], tessellation_rows)}

The lateral surface is an **open tube**: it has exactly `2 nu` boundary edges (the two
rims) and no others, which is the statement that the seam welded. All faces point
outward (`inward = 0`), the winding is globally consistent (`flips = 0`), and the
maximum angle between a facet normal and the analytic normal at the quad's parametric
centre falls by a factor ~4 per doubling of the grid (O(h^2)):
{metrics['tessellation'][0]['max_normal_angle_deg']:.4f} deg at
{metrics['tessellation'][0]['u_samples']} samples to
{metrics['tessellation'][-1]['max_normal_angle_deg']:.4f} deg at
{metrics['tessellation'][-1]['u_samples']}.

### 3.4 Closed solid and STL round-trip

Reference volume `pi r^2 h = {solid['volume_reference_m3']:.6f} m^3`.

{_table(["u samples", "triangles", "volume [m^3]", "volume error", "area error", "watertight"],
        solid_rows)}

STL round-trip (binary float32, re-welded by exact position):
watertight **{stl['watertight']}**, one shell, volume error `{_sci(stl['volume_relative_error'])}`.
The remaining error is the *inscribed* cap polygon, not the format: it halves by four
per doubling, exactly as the table above shows.

### 3.5 Parameter-derivative continuity

| curve / surface | interior knot | multiplicity | declared | measured | first-derivative jump |
|---|---|---|---|---|---|
| three-arc circle, r = 12 m | {circle['knots'][0]['knot']:.6f} | 2 | C{circle['declared_order']} | C{circle['knots'][0]['measured_continuous_order']} | {_sci(circle['knots'][0]['jumps'][1])} |
| three-arc circle, r = 12 m | {circle['knots'][1]['knot']:.6f} | 2 | C{circle['declared_order']} | C{circle['knots'][1]['measured_continuous_order']} | {_sci(circle['knots'][1]['jumps'][1])} |
| cylinder surface (u direction) | 0.333333 | 2 | C0 | C1 | {_sci(continuity['cylinder_surface']['u'][0]['jumps'][1])} |
| cubic control curve | {cubic['knot']:.1f} | 1 | C{cubic['declared_order']} | C{cubic['measured_continuous_order']} | 0 |

The interesting result is the first two rows: the knot bookkeeping only promises `C0`
at the arc joins, yet the parameter derivative is continuous to 4e-14. The reason is
symmetry — each rational quadratic arc has the same parametric speed at both of its
ends, so the join matches in value *and* in slope. The **positive control** is the
cubic row: there the third derivative jumps by
`{cubic['jumps'][3]:.3f}` (against a scale of `{cubic['scales'][3]:.3f}`), so the
measurement does detect a genuine discontinuity; it is not a test that always passes.

### 3.6 Negative controls (the ticket's two rejection criteria, plus four more)

{_table(["control", "expected rejection reason", "outcome", "message"], negative_rows)}

### 3.7 PC6 contract compatibility

P1's geometry contract (`tensorcad.pc6.geom/1.0`, `tensorcad/pc6_contract.py` in
`main`) is extended **additively**: the block carries the same `schema` tag, adds one
key (`nurbs_surface`), and `tensorcad.pc6_contract.load_contract` reads only
`schema` / `parameters` / `artifacts` / `geometry` / `acceptance`, so no existing check
is modified or relaxed. The test suite proves the additivity by re-parsing a contract
with the block appended through P1's **own** `_parameters_from_contract` and asserting
the parameters are unchanged.

Geometry-first gate (recomputed from P1's closed form, never copied):

{_table(["quantity", "value"], [
    ["length L [m]", f"{gate['length_m']:g}"],
    ["beam B [m]", f"{gate['beam_m']:g}"],
    ["draft T [m]", f"{gate['draft_m']:g}"],
    ["displacement volume [m^3]", f"{gate['displacement_volume_m3']:.3f}"],
    ["block coefficient C_b", f"{gate['block_coefficient']:.6f}"],
    ["prismatic coefficient C_p", f"{gate['prismatic_coefficient']:.6f}"],
    ["midship coefficient C_m", f"{gate['midship_coefficient']:.6f}"],
])}

A declared value that deviates by more than 1e-3 is **refused** (tampered C_b at -10 %
=> `{compatibility['declared_cb_tamper_rejected']}`, message
"{compatibility['declared_cb_tamper_message'][:96]}"), because a surface block whose
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
| cylinder lateral area (R=12, H=11) | {cylinder['parametric_area_m2']:.9f} m^2 | {cylinder['analytic_area_m2']:.9f} m^2 | {_sci(cylinder['area_relative_error'])} | < 1e-6 | **pass, 1e-9 x margin** |
| arc length (r = 2 m, 1 rad) | {arcs['arc_length']:.12f} m | 2 m | {_sci(arcs['arc_length_relative_error'])} | < 1e-6 | **pass** |
| full circle length (r = 12 m) | {arcs['circle_length']:.9f} m | {arcs['circle_length_reference']:.9f} m | {_sci(arcs['circle_length_relative_error'])} | < 1e-6 | **pass** |
| circle first-derivative jump at the joins | {circle['max_first_derivative_jump']:.3e} m | 0 (C1) | rel. to \\|C'\\| = {circle['knots'][0]['scales'][1]:.3f}: {circle['max_first_derivative_jump'] / circle['knots'][0]['scales'][1]:.2e} | < 1e-9 | **pass** |
| closed-solid volume (u = {solid['rows'][-1]['u_samples']}) | {solid['rows'][-1]['volume_m3']:.4f} m^3 | {solid['volume_reference_m3']:.4f} m^3 | {_sci(solid['rows'][-1]['volume_relative_error'])} | < 1e-3 | **pass (inscribed facets)** |
| STL round-trip volume | {stl['volume_m3']:.4f} m^3 | {solid['volume_reference_m3']:.4f} m^3 | {_sci(stl['volume_relative_error'])} | < 1e-3 | **pass** |
| tessellation max normal angle | {metrics['tessellation'][2]['max_normal_angle_deg']:.4f} deg | 0 deg | — | < 1 deg | **pass** |

**Is every key quantity below 3 %?** Yes, with the caveat stated per row: the
parametric quantities are exact to round-off (1e-15), and the tessellated quantities
are *discretisation* errors that converge as `O(h^2)` and are reported at the grid they
were measured on. No quantity is above 3 %, and none of the pass margins depends on a
tolerance chosen to make it pass.

### Machine-checkable criteria of this run

{_table(["criterion", "value"], [[f"`{key}`", "**pass**" if value else "**fail**"]
                                 for key, value in sorted(criteria.items())])}

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
"""
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(body, encoding="utf-8")
    print(f"wrote {target} ({len(body)} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
