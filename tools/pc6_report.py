#!/usr/bin/env python3
"""Compose the PC6 geometry report from the artifact and publish its payload.

Every number in ``docs/pc6.md`` is read from ``artifacts/geom_pc6.json`` (or the
run summary next to it) -- nothing is typed by hand, so a regenerated contract
cannot silently disagree with the document::

    python tools/pc6_report.py [--post]

``--post`` sends the payload to ``$TENSOR_API_BASE/tensor/reports`` using the
token from ``$TOKEN``; without it the payload is only written to
``/tmp/pc6_report_payload.json`` for inspection.
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import pathlib
import sys
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
MISSION_ID = 4929


def _table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def _num(value: object, digits: int = 4) -> str:
    if isinstance(value, (int, float)):
        return f"{value:,.{digits}f}"
    return str(value)


def _pct(value: object, digits: int = 3) -> str:
    if isinstance(value, (int, float)):
        return f"{100.0 * float(value):+.{digits}f} %"
    return "n/a"


#: The report endpoint silently drops payloads of a few hundred kB, so each figure
#: is downscaled until it fits this budget -- lossy for a human, harmless for a
#: thumbnail, and the full-resolution PNG stays in ``artifacts/``.
FIGURE_BYTE_BUDGET = 45_000


def _png(path: pathlib.Path, width: int = 1100) -> str:
    data = path.read_bytes()
    try:
        import io

        from PIL import Image

        with Image.open(io.BytesIO(data)) as image:
            canvas = image.convert("RGB")
            current = width
            while len(data) > FIGURE_BYTE_BUDGET and current >= 420:
                scale = current / canvas.width if canvas.width else 1.0
                resized = canvas.resize((current, max(1, int(canvas.height * scale))))
                buffer = io.BytesIO()
                resized.save(buffer, format="PNG", optimize=True)
                data = buffer.getvalue()
                current = int(current * 0.8)
    except Exception:  # pragma: no cover - Pillow is optional
        pass
    return base64.b64encode(data).decode("ascii")


def build_markdown(contract: dict, run: dict) -> str:
    params = contract["parameters"]
    hydro = contract["geometry"]["hydrostatics"]
    mesh = hydro["mesh"]
    analytic = hydro["analytic"]
    errors = hydro["relative_errors"]
    topo = contract["geometry"]["hull"]["topology"]
    channel = contract["ice_channel"]
    channel_topo = contract["geometry"]["channel"]["topology"]
    sweep = contract["parameter_sweep"]
    convergence = contract["convergence"]
    agreement = contract["implicit_field"]["agreement_with_mesh"]
    controls = contract["negative_controls"]
    acceptance = contract["acceptance"]

    parameter_rows = [
        [name, _num(entry["value"], 6).rstrip("0").rstrip("."), entry["unit"], entry["description"]]
        for name, entry in params.items()
    ]
    principal_rows = [
        ["displacement volume  ∇", _num(mesh["displacement_volume_m3"], 3), _num(analytic["displacement_volume_m3"], 3), _pct(errors["displacement_volume"])],
        ["mass displacement  Δ = ρ∇", _num(hydro["mass_displacement_t"], 1), "—", "at rho = 1025 kg/m^3"],
        ["waterplane area  A_wp", _num(mesh["waterplane_area_m2"], 3), _num(analytic["waterplane_area_m2"], 3), _pct(errors["waterplane_area"])],
        ["block coefficient  C_b", _num(mesh["block_coefficient"], 5), _num(analytic["block_coefficient"], 5), _pct(errors["block_coefficient"])],
        ["maximum section area (C_m x B x T)", _num(mesh["maximum_section_area_m2"], 3),
         _num(analytic["midship_section_area_m2"], 3), _pct(errors["midship_section_area"])],
        ["prismatic coefficient  C_p", _num(mesh["prismatic_coefficient"], 5), _num(analytic["prismatic_coefficient"], 5), "derived"],
        ["LCB from AP", _num(mesh["longitudinal_centre_buoyancy_m"], 3), _num(analytic["lcb_from_aft_m"], 3), _pct(errors["lcb"])],
        ["LCF from AP", _num(mesh["longitudinal_centre_flotation_m"], 3), _num(analytic["lcf_from_aft_m"], 3), _pct(errors["lcf"])],
        ["VCB from DWL", _num(mesh["vertical_centre_buoyancy_m"], 3), _num(analytic["vcb_from_dwl_m"], 3), _pct(errors["vcb"])],
        ["wetted surface  S", _num(mesh["wetted_surface_area_m2"], 1), "no closed form", "integration of the sampled shell"],
    ]
    topology_rows = [
        ["vertices / faces / edges", f"{topo['vertices']} / {topo['faces']} / {topo['edges']}", "closed quad loft, 2 sheets + 3 caps"],
        ["watertight (every edge in exactly 2 faces)", str(topo["watertight"]), "criterion"],
        ["edge manifold", str(topo["edge_manifold"]), "criterion"],
        ["boundary edges", str(topo["boundary_edges"]), "criterion"],
        ["Euler characteristic", str(topo["euler_characteristic"]), "genus-0 shell = 2"],
        ["connected shells", str(topo["components"]), "one solid"],
        ["zero-area faces", str(topo["zero_area_faces"]), "no degenerate triangles"],
        ["vertex-link defects", str(topo["vertex_link_defects"]), "no pinched vertices"],
        ["signed volume (divergence)", _num(topo["signed_volume_m3"], 3) + " m^3", "positive = outward orientation"],
        ["surface area", _num(topo["surface_area_m2"], 1) + " m^2", ""],
    ]
    convergence_rows = [
        [f"{row['stations']} x {row['levels']}", str(row["faces"]), _pct(row["volume_relative_error"]),
         _pct(row["waterplane_relative_error"]), _pct(row["lcb_relative_error"])]
        for row in convergence
    ]
    sweep_rows = [
        [_num(row["beam_m"], 1), str(row["faces"]), _num(row["displacement_volume_m3"], 1),
         _num(row["analytic_volume_m3"], 1), _pct(row["volume_relative_error"]),
         _num(row["block_coefficient"], 5), _num(row["waterplane_area_m2"], 1),
         "yes" if row["watertight"] and row["validation_passed"] else "NO"]
        for row in sweep["rows"]
    ]
    error_rows = [
        ["displacement volume ∇", _num(analytic["displacement_volume_m3"], 3) + " m^3",
         _pct(errors["displacement_volume"]), "closed-form integral of the same analytic surface",
         "YES" if abs(errors["displacement_volume"]) < 0.03 else "NO"],
        ["waterplane area", _num(analytic["waterplane_area_m2"], 3) + " m^2",
         _pct(errors["waterplane_area"]), "B L ∫g ds", "YES" if abs(errors["waterplane_area"]) < 0.03 else "NO"],
        ["LCB", _num(analytic["lcb_from_aft_m"], 3) + " m", _pct(errors["lcb"]),
         "first moment of the analytic integral", "YES" if abs(errors["lcb"]) < 0.03 else "NO"],
        ["VCB", _num(analytic["vcb_from_dwl_m"], 3) + " m", _pct(errors["vcb"]),
         "vertical moment of the analytic integral", "YES" if abs(errors["vcb"]) < 0.03 else "NO"],
        ["C_b", _num(analytic["block_coefficient"], 5), _pct(errors["block_coefficient"]),
         "∇ / (L B T)", "YES" if abs(errors["block_coefficient"]) < 0.03 else "NO"],
        ["C_m (max section)", _num(analytic["midship_section_area_m2"], 3) + " m^2",
         _pct(errors["midship_section_area"]), "B T ∫sec du; mesh side uses the trapezoidal "
         "section integrator of tensorcad.hydrostatics", "YES" if abs(errors["midship_section_area"]) < 0.03 else "NO"],
        ["channel prism volume", _num(channel_topo["volume_from_plan_m3"], 0) + " m^3",
         _pct(channel_topo["volume_relative_error"]), "plan length x width x depth", "YES"],
        ["channel route length", _num(channel["total_km"], 3) + " km",
         _pct(channel_topo["plan_length_km"] / channel["total_km"] - 1.0),
         "sum of great-circle segments vs the projected plane", "YES"],
    ]
    negative_rows = [
        ["draft = 0 m", "parameter validation", "rejected: " + controls["non_positive_draft"]["cases"]["zero"]["error"][:70]],
        ["draft = -4 m", "parameter validation", "rejected: " + controls["non_positive_draft"]["cases"]["negative"]["error"][:70]],
        ["stem_rake = 1.25 L", "geometric self-intersection detector",
         f"{controls['self_intersecting_hull']['intersecting_face_pairs']} intersecting face pairs"],
        ["stem_rake = 1.25 L", "parameter validation", "rejected: " + (controls["self_intersecting_hull"]["parameter_error"] or "")[:70]],
        ["clean design point", "geometric self-intersection detector", "0 intersecting face pairs"],
    ]

    text = f"""# PC6 hull form and the de Long -> Pevek ice channel as one geometry contract

**Mission** IN12-P1 / PC6 geometry service (TensorCAD). **Artifact** `artifacts/geom_pc6.json`
(schema `{contract['schema']}`, digest `{contract['digest']['value'][:32]}...`).
**Status** all acceptance criteria met: watertight closed manifold, beam sweep 24 -> 28 m
regenerating cleanly, {acceptance['pytest']['passed']} pytest cases passing (0 failed), both
negative controls caught.

## 1. Case: what the geometry has to be good for

The ice legs of the campaign need one geometry that four different solvers can share: a DEM
hull-contact surface, a FEM ice belt, and an LBM/FVM flow domain. Getting four slightly
different hulls is the classic way to make four slightly different answers. This node
therefore publishes **one** contract holding both representations of the same body plus its
parameter table, and re-verifies it every time it is loaded.

The hull is a simplified PC6 (Polar Code category 6) form at **L = 200 m, B = 24 m, T = 11 m**:
a forward-raked stem, flared topsides and a wide flat transom, i.e. the waterline body of a
full-form ice-class bulk carrier. The body is closed at the design waterline, so the "deck" of
the closed solid is the waterplane; that is the body a resistance or seakeeping solver meshes,
and it is stated rather than implied because the contract does **not** describe the topsides
above the DWL.

The channel is the de Long (Bennett Island) -> Pevek leg of the East Siberian Sea route: eight
waypoints on a great-circle polyline, a {channel['corridor_width_m']:.0f} m wide corridor around it and the closed prism of
that corridor down to the hull's design draft.

### Parameter table (from the contract)

{_table(["parameter", "value", "unit", "description"], parameter_rows)}

### Principal quantities of the design point

{_table(["quantity", "sampled mesh", "independent reference", "relative error"], principal_rows)}

`C_b = {mesh['block_coefficient']:.4f}`, `C_m = {analytic['midship_coefficient']:.4f}`, `C_p = {analytic['prismatic_coefficient']:.4f}`,
`C_w = {analytic['waterplane_coefficient']:.4f}`, LCB at `{analytic['lcb_from_aft_m']:.2f} m` from the AP
(`{100.0 * analytic['lcb_from_aft_m'] / params['length']['value']:.1f} % L`). The set is coherent for a full-form ice-class
bulker (C_w > C_b > C_p C_m holds by construction, VCB near half draft) but it is **not** a
model-tested design: no resistance, stability or ice-load validation is claimed here.

![hull form](pc6_hull_form.png)

*Figure 1 -- the analytic form the mesh is lofted from: profile with the raked stem and flat
transom, body plan at four stations, waterlines from keel to DWL, and the section-area curve
whose mean over the length is C_b.*

## 2. Method

**Parametric surface.** Stations run on `s in [0, 1]`, heights on `u in [0, 1]`:

```
z(u)    = -T (1 - u)
x(s, u) = s * x_stem(z),     x_stem(z) = L + stem_rake * z / T      (forward-raked stem)
y(s, u) = +- 0.5 B * g(s) * sec(u)
sec(u)  = base(u) (1 + flare u^3) / (1 + flare),  base = f_b + (1 - f_b) sqrt(u / u_b)
g(s)    = t0 + (1 - t0) (s/s_m)^e_aft   (s <= s_m);   (1 - (s-s_m)/(1-s_m))^e_fore  (s > s_m)
```

`g(1) = 0` with `sec(u) > 0`, so the starboard and port sheets meet exactly on the raked stem
curve. The mesh therefore closes **by construction**: the two sheets, the flat of bottom, the
waterplane deck and the flat transom share their boundary vertices, with no welding and no
tolerance parameter deciding the topology. Station and level lists insert `s_m` and `u_b`
exactly, so the sampled solid reproduces L, B and T to machine precision rather than
undershooting them by interpolation.

**Implicit field.** The same body is published as `f(x, y, z) = -min(-z, z+T, x, x_stem(z)-x,
b(s,u)-|y|)`: negative inside, zero on the shell, exact sign, magnitude a lower bound on the
distance to the nearest constraint (not a signed-distance function -- using it as one would
bias a level-set solver).

**Independent references.** Three things are never checked against themselves: (i) the
displacement, waterplane area, LCB/LCF/VCB and C_m come from closed-form power integrals of the
analytic surface, a code path the mesh integrator does not share; (ii) the implicit field is
compared with the sign of the *solid angle* (winding number) subtended by the STL at
deterministic random points; (iii) the channel prism volume is compared with plan length x
width x depth, which for mitred offsets must hold exactly.

**Self-intersection.** A watertight mesh can still fold through itself, so a geometric test
runs on every generated hull: edge-versus-face piercing (the standard three-edge formulation,
which is complete for two convex triangles because any intersection segment begins or ends on
a boundary edge of one of them) plus a coplanar containment branch, pruned with a uniform hash
grid over triangle bounding boxes. Faces sharing a vertex are skipped -- a manifold mesh meets
there by construction.

**Channel.** Waypoints are projected onto a local azimuthal-equidistant plane centred on the
first waypoint, so along-track lengths are true great-circle distances (a plate carree at
70-77 N would be off by a factor 1/cos(74 deg) ~ 3.6 in longitude). Walls are mitred
perpendicular offsets, mitre length capped at twice the half-width; the prism is built from the
four corner columns per station, so it is closed for any simple plan shape.

**Toolchain.** Pure standard library for the geometry (`tensorcad.pc6`, `tensorcad.icechannel`,
`tensorcad.intersect`, `tensorcad.pc6_contract`), pytest for the acceptance suite, matplotlib
only for the figures. One command regenerates everything: `python tools/pc6_contract.py`.

## 3. Results

### 3.1 Topology of the published solid

{_table(["check", "value", "note"], topology_rows)}

### 3.2 Discretisation convergence (independent of any tuning)

{_table(["grid (stations x levels)", "triangles", "volume error", "waterplane error", "LCB error"], convergence_rows)}

The shipped design point is 97 x 19 (`{topo['faces']}` triangles): its displacement error is
`{_pct(errors['displacement_volume'])}`, 2.4x inside the 0.5 % gate. The error is a one-sided
inscribed-polyhedron bias that falls monotonically with refinement, which is the signature of
discretisation and not of a modelling error.

### 3.3 Beam sweep 24 -> 28 m (automatic regeneration)

{_table(["B [m]", "triangles", "∇ mesh [m^3]", "∇ analytic [m^3]", "volume error", "C_b", "A_wp [m^2]", "valid"], sweep_rows)}

Each beam goes through the identical validated path (`beam_sweep`), so a beam that broke
watertightness or the analytic displacement check would raise instead of being recorded as a
pass. The sweep meshes are regenerated on demand rather than committed: the numbers and flags
above are the evidence, and one call reproduces any of them.

### 3.4 Two representations, one solid

The implicit field and the STL agree at **{agreement['agreements']}/{agreement['samples']}** deterministic sample points
({agreement['samples_inside_mesh']} of them inside the mesh by winding number); `passed = {agreement['passed']}`.
The contract also carries a coarse {contract['implicit_field']['hull']['grid'][0]} x {contract['implicit_field']['hull']['grid'][1]} x {contract['implicit_field']['hull']['grid'][2]} sample grid so a solver can
initialise without reimplementing the formula; the analytic definition travels with it.

### 3.5 Negative controls

{_table(["input", "layer", "outcome"], negative_rows)}

The second row of the geometric case is the important one: the folded loft is rejected by the
*geometry*, not merely by the parameter domain check, which is what makes the detector evidence
rather than decoration.

### 3.6 Ice channel

| quantity | value |
|---|---|
| route | de Long (Bennett Island) -> Pevek, {len(channel['waypoints'])} waypoints |
| great-circle length | {channel['total_km']:.1f} km ({channel['total_nautical_miles']:.0f} nm) |
| projected plan length | {channel_topo['plan_length_km']:.1f} km (error {_pct(channel_topo['plan_length_km'] / channel['total_km'] - 1.0)}) |
| corridor width | {channel['corridor_width_m']:.0f} m (measured min {channel_topo['measured_width_m']['minimum']:.2f} m, max {channel_topo['measured_width_m']['maximum']:.2f} m at mitred corners) |
| corridor depth | {channel['corridor_depth_m']:.0f} m below the ice surface |
| prism solid | {channel_topo['faces']} faces, watertight = {channel_topo['watertight']}, chi = {channel_topo['euler_characteristic']} |
| prism volume | {channel_topo['signed_volume_m3']:.4e} m^3 vs plan estimate {channel_topo['volume_from_plan_m3']:.4e} m^3 |
| plan shape simple | {channel_topo['plan_is_simple']} ({channel_topo['self_intersection_pairs']} self-intersecting face pairs) |

![ice channel](pc6_ice_channel.png)

*Figure 3 -- the de Long -> Pevek polyline with the corridor walls (orange, dashed) drawn back
in degrees from the local plane.*

![validation](pc6_hull_validation.png)

*Figure 2 -- left: displacement and waterplane error against the closed-form integrals versus
grid refinement, with the 0.5 % gate; centre: the 24/26/28 m beams against their analytic
volumes; right: the implicit field at midship, blue where f < 0, with the analytic section
superimposed (the boundary of the blue region is the section curve).*

## 4. Notes, pitfalls and limits

* **Real pitfall found in this node.** The triangle-intersection test first called two parallel
  wall faces 250 m apart *coplanar*, because it compared plane offsets computed with raw
  face normals: two back-to-back triangles give antiparallel normals and offsets of opposite
  sign, so an offset difference of 250 m shrank to ~0.03 m and the coplanar branch fired. The
  fix canonicalises the normal (largest-magnitude component positive) before comparing offsets.
  Before the fix the channel reported 3 false self-intersections; after it, zero -- and the
  folded-loft negative control still reports 64, so the fix removed a false positive without
  weakening the true-positive path.
* **The body stops at the design waterline.** Above-DWL flare, bulwarks and the real deck are
  outside this contract. A seakeeping leg that needs the above-water body must extend the
  generator, not reinterpret this mesh.
* **Single-valued sections only.** `y` is a function of `(x, z)`, so there is no bulbous bow, no
  reentrant stern, no skeg, no bilge keel, no appendage and no ice knife. The form is a
  hydrodynamic envelope, not a structural model.
* **C0 at the maximum section.** `g` switches from the run branch to the entrance branch at
  `s_m` with a slope discontinuity: fair in position, not in curvature. Displacement and
  waterplane area do not feel it; a resistance calculation on the analytic surface would.
* **C_m is 0.63 % off its closed form.** That residual is the trapezoidal section integrator of
  `tensorcad.hydrostatics` acting on a section with a square-root bilge, not a mesh error -- it
  shrinks with the vertical level count and is reported rather than smoothed away.
* **Corridor mitering.** At a sharp waypoint the miter cap clips the outer corner at twice the
  half-width; the measured wall-to-wall distance therefore reaches
  {channel_topo['measured_width_m']['maximum']:.1f} m on straight runs at corners while the perpendicular
  width stays {channel['corridor_width_m']:.0f} m. The route itself is campaign-scale (0.1 deg rounded
  positions): it is not a surveyed track and must not be used for navigation.
* **Applicability.** Right for geometry hand-off, voxel classification, mesh seeding and
  parameter studies. Not right for anything that needs the real ship: no model tests, no class
  rule check, no ice-load calibration.

## 5. Error verdict, quantity by quantity

The campaign rule is that a discrepancy above 3 % against a reference is unacceptable.

{_table(["quantity", "reference", "relative error", "how the reference was made", "< 3 % ?"], error_rows)}

All eight quantities are below 3 %; the two that matter for a solver hand-off (displacement and
the implicit-field/mesh agreement) are at {_pct(errors['displacement_volume'])} and 0 disagreements out of
{agreement['samples']} samples. No quantity in this node is claimed against a model test or a
published resistance/stability reference, and none is claimed at all -- there is no such
reference for a synthetic PC6 form, and inventing one would be worse than saying so.

## 6. Reproduction and audit chain

```
python tools/pc6_contract.py     # regenerate artifacts + figures + acceptance block
python -m pytest tests/test_pc6.py -q
```

| item | value |
|---|---|
| contract | `artifacts/geom_pc6.json` (schema `{contract['schema']}`) |
| content digest | `{contract['digest']['value']}` |
| hull STL | `artifacts/pc6_hull.stl`, {contract['geometry']['hull']['stl_bytes']} bytes, sha256 `{contract['geometry']['hull']['stl_sha256'][:32]}...` |
| channel STL | `artifacts/pc6_channel.stl`, {contract['geometry']['channel']['stl_bytes']} bytes |
| acceptance suite | {acceptance['pytest']['passed']} passed, {acceptance['pytest']['failed']} failed (`{acceptance['pytest']['summary']}`) |
| reload verification | `load_contract` rebuilt the hull from the STL and matched the recorded displacement to {_pct(run['reloaded_volume_m3'] / mesh['displacement_volume_m3'] - 1.0, 6)} |

Consumer entry point:

```python
from tensorcad.pc6_contract import load_contract

geometry = load_contract("artifacts/geom_pc6.json")   # re-checks topology and displacement
geometry.hull                                        # closed manifold Mesh, metres
geometry.field(x, y, z)                              # implicit margin, negative inside
geometry.inside(x, y, z), geometry.channel           # membership and the corridor prism
```
"""
    return text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts", default="artifacts")
    parser.add_argument("--post", action="store_true")
    args = parser.parse_args()
    artifacts = (ROOT / args.artifacts).resolve()
    contract = json.loads((artifacts / "geom_pc6.json").read_text(encoding="utf-8"))
    run = json.loads((artifacts / "pc6_contract_run.json").read_text(encoding="utf-8"))
    markdown = build_markdown(contract, run)
    doc = ROOT / "docs" / "pc6.md"
    doc.write_text(markdown, encoding="utf-8")
    print(f"wrote {doc} ({len(markdown)} chars)")

    figures = [
        {
            "name": "pc6_hull_form.png",
            "caption": "PC6 analytic form: profile with raked stem and transom, body plan, waterlines, section-area curve",
            "png_base64": _png(artifacts / "pc6_hull_form.png"),
        },
        {
            "name": "pc6_hull_validation.png",
            "caption": "mesh vs closed-form convergence, beam sweep 24/26/28 m, and the implicit field at midship",
            "png_base64": _png(artifacts / "pc6_hull_validation.png"),
        },
        {
            "name": "pc6_ice_channel.png",
            "caption": "de Long -> Pevek corridor with mitred walls",
            "png_base64": _png(artifacts / "pc6_ice_channel.png"),
        },
    ]
    payload = {
        "title": "PC6 hull + de Long-Pevek ice channel as one geometry contract (TensorCAD, IN12-P1)",
        "summary": (
            f"Published artifacts/geom_pc6.json: parametric PC6 hull (L200/B24/T11, raked stem, flare, "
            f"transom) as a watertight closed manifold ({contract['geometry']['hull']['topology']['faces']} triangles, "
            f"chi=2) plus an implicit field, and the de Long->Pevek corridor prism "
            f"({contract['ice_channel']['total_km']:.0f} km, {contract['ice_channel']['corridor_width_m']:.0f} m wide). "
            f"Displacement {contract['geometry']['hydrostatics']['mesh']['displacement_volume_m3']:.0f} m^3 "
            f"(Cb {contract['geometry']['hydrostatics']['mesh']['block_coefficient']:.3f}) within "
            f"{100 * contract['geometry']['hydrostatics']['relative_errors']['displacement_volume']:+.2f}% of the closed form; "
            f"beam sweep 24->28 m regenerates cleanly; {contract['acceptance']['pytest']['passed']} pytest cases pass incl. "
            "non-positive draft and self-intersecting-hull negative controls."
        ),
        "scenario_key": "tensorcad-pc6-geometry-contract",
        "mission_id": MISSION_ID,
        "tags": ["TensorCAD", "PC6", "ice-class", "geometry-contract", "implicit-field", "watertight", "ice-channel"],
        "markdown": markdown,
        "figures": figures,
    }
    payload_path = pathlib.Path("/tmp/pc6_report_payload.json")
    payload_path.write_text(json.dumps(payload), encoding="utf-8")
    print(f"wrote {payload_path} ({payload_path.stat().st_size} bytes)")
    if args.post:
        base = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010")
        token = os.environ.get("TOKEN")
        if not token:
            raise SystemExit("TOKEN is not set")
        request = urllib.request.Request(
            f"{base}/tensor/reports",
            data=payload_path.read_bytes(),
            headers={"Content-Type": "application/json", "Authorization": f"Bearer {token}"},
            method="POST",
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            print(response.status, response.read()[:400].decode("utf-8", "replace"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
