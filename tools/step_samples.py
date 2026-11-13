"""Generate small AP214 analytic B-rep solids with closed-form nominal volumes.

The files this writes are *fixtures* with an honest provenance: their geometry is
defined here in closed form (plate with a through hole, stepped shaft with a
conical transition, bushing, L-bracket), so the nominal volume is a formula rather
than a number copied from anywhere, and the acceptance check "mesh volume vs CAD
nominal < 1 %" has an independent reference. They are written as genuine
``MANIFOLD_SOLID_BREP`` / ``ADVANCED_FACE`` AP214 files -- analytic PLANE /
CYLINDRICAL_SURFACE / CONICAL_SURFACE entities trimmed by ``EDGE_CURVE`` loops --
which is exactly the entity family :mod:`tensorcad.step_brep` reads and
:mod:`tensorcad.exchange` refuses.

Run::

    python tools/step_samples.py --out artifacts/gap-cad-b1-step-reader/samples

The script writes ``nominal_volumes.json`` next to the STEP files so a report can
quote measured numbers from an artifact instead of from a transcript.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

SCHEMA = "AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }"


def _f(value: float) -> str:
    text = f"{value:.12g}"
    if "e" in text or "E" in text:
        return text
    if "." not in text:
        text += "."
    return text


def _add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _mul(a, k):
    return (a[0] * k, a[1] * k, a[2] * k)


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _unit(a):
    length = math.sqrt(_dot(a, a))
    return (a[0] / length, a[1] / length, a[2] / length)


class SolidWriter:
    """Minimal AP214 writer: enough product structure to look like a real part."""

    def __init__(self, part_name: str):
        self.part_name = part_name
        self.entities: list[str] = []
        self._vertices: dict[tuple[float, float, float], int] = {}
        self.length_unit = self.raw("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))")
        self.angle_unit = self.raw("(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))")
        self.solid_angle_unit = self.raw("(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())")
        self.uncertainty = self.raw(
            f"UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#{self.length_unit},"
            "'distance_accuracy_value','confusion accuracy')"
        )
        self.context = self.raw(
            "(GEOMETRIC_REPRESENTATION_CONTEXT(3)"
            f"GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{self.uncertainty}))"
            f"GLOBAL_UNIT_ASSIGNED_CONTEXT((#{self.length_unit},#{self.angle_unit},#{self.solid_angle_unit}))"
            "REPRESENTATION_CONTEXT('',''))"
        )

    # -- primitive emitters --------------------------------------------------
    def raw(self, text: str) -> int:
        self.entities.append(text)
        return len(self.entities)

    def point(self, p) -> int:
        return self.raw(f"CARTESIAN_POINT('',({_f(p[0])},{_f(p[1])},{_f(p[2])}))")

    def direction(self, d) -> int:
        d = _unit(d)
        return self.raw(f"DIRECTION('',({_f(d[0])},{_f(d[1])},{_f(d[2])}))")

    def vector(self, d, magnitude: float = 1.0) -> int:
        return self.raw(f"VECTOR('',#{self.direction(d)},{_f(magnitude)})")

    def axis2(self, origin, axis, ref_direction) -> int:
        return self.raw(
            f"AXIS2_PLACEMENT_3D('',#{self.point(origin)},#{self.direction(axis)},"
            f"#{self.direction(ref_direction)})"
        )

    def vertex(self, p) -> int:
        key = (round(p[0], 9), round(p[1], 9), round(p[2], 9))
        cached = self._vertices.get(key)
        if cached is not None:
            return cached
        created = self.raw(f"VERTEX_POINT('',#{self.point(p)})")
        self._vertices[key] = created
        return created

    def line_edge(self, start_point, end_point) -> int:
        line = self.raw(f"LINE('',#{self.point(start_point)},#{self.vector(_sub(end_point, start_point))})")
        start = self.vertex(start_point)
        end = self.vertex(end_point)
        return self.raw(f"EDGE_CURVE('',#{start},#{end},#{line},.T.)")

    def circle_edge(self, center, normal, ref_direction, radius: float) -> int:
        axis = self.axis2(center, normal, ref_direction)
        circle = self.raw(f"CIRCLE('',#{axis},{_f(radius)})")
        start = self.vertex(_add(center, _mul(_unit(ref_direction), radius)))
        return self.raw(f"EDGE_CURVE('',#{start},#{start},#{circle},.T.)")

    def oriented(self, edge: int, orientation: bool = True) -> int:
        return self.raw(f"ORIENTED_EDGE('',*,*,#{edge},{'.T.' if orientation else '.F.'})")

    def loop(self, oriented_edges) -> int:
        refs = ",".join(f"#{edge}" for edge in oriented_edges)
        return self.raw(f"EDGE_LOOP('',({refs}))")

    # -- faces ---------------------------------------------------------------
    def planar_face(self, normal, x_dir, loop3d, holes=()) -> int:
        normal = _unit(normal)
        axis = self.axis2(loop3d[0], normal, x_dir)
        surface = self.raw(f"PLANE('',#{axis})")
        oriented = [
            self.oriented(self.line_edge(loop3d[index], loop3d[(index + 1) % len(loop3d)]))
            for index in range(len(loop3d))
        ]
        bounds = [self.raw(f"FACE_OUTER_BOUND('',#{self.loop(oriented)},.T.)")]
        for center, radius in holes:
            edge = self.circle_edge(center, normal, _unit(x_dir), radius)
            bounds.append(self.raw(f"FACE_BOUND('',#{self.loop([self.oriented(edge)])},.T.)"))
        refs = ",".join(f"#{bound}" for bound in bounds)
        return self.raw(f"ADVANCED_FACE('',({refs}),#{surface},.T.)")

    def disc_face(self, center, normal, x_dir, radius: float) -> int:
        normal = _unit(normal)
        axis = self.axis2(center, normal, x_dir)
        surface = self.raw(f"PLANE('',#{axis})")
        edge = self.circle_edge(center, normal, x_dir, radius)
        bound = self.raw(f"FACE_OUTER_BOUND('',#{self.loop([self.oriented(edge)])},.T.)")
        return self.raw(f"ADVANCED_FACE('',(#{bound}),#{surface},.T.)")

    def annulus_face(self, center, normal, x_dir, outer_radius: float, inner_radius: float) -> int:
        normal = _unit(normal)
        axis = self.axis2(center, normal, x_dir)
        surface = self.raw(f"PLANE('',#{axis})")
        outer = self.circle_edge(center, normal, x_dir, outer_radius)
        inner = self.circle_edge(center, normal, x_dir, inner_radius)
        outer_bound = self.raw(f"FACE_OUTER_BOUND('',#{self.loop([self.oriented(outer)])},.T.)")
        inner_bound = self.raw(f"FACE_BOUND('',#{self.loop([self.oriented(inner)])},.T.)")
        return self.raw(f"ADVANCED_FACE('',(#{outer_bound},#{inner_bound}),#{surface},.T.)")

    def _revolved_face(self, surface_text: str, bottom_center, top_center, axis_dir, x_dir, r_lo: float, r_hi: float, same_sense: bool) -> int:
        """A full 2*pi revolved face: two closed rings joined by a seam edge used twice.

        That is how a valid EDGE_LOOP closes around a periodic surface: the ring at
        the bottom is traversed forward, the seam up, the top ring backward, the same
        seam down. Emitting two independent closed rings instead would describe a
        disconnected wire whose parameter trace cancels to zero area.
        """
        surface = self.raw(surface_text)
        low_start = _add(bottom_center, _mul(_unit(x_dir), r_lo))
        high_start = _add(top_center, _mul(_unit(x_dir), r_hi))
        low_edge = self.circle_edge(bottom_center, axis_dir, x_dir, r_lo)
        high_edge = self.circle_edge(top_center, axis_dir, x_dir, r_hi)
        seam = self.line_edge(low_start, high_start)
        oriented = [
            self.oriented(low_edge, True),
            self.oriented(seam, True),
            self.oriented(high_edge, False),
            self.oriented(seam, False),
        ]
        bound = self.raw(f"FACE_OUTER_BOUND('',#{self.loop(oriented)},.T.)")
        return self.raw(f"ADVANCED_FACE('',(#{bound}),#{surface},{'.T.' if same_sense else '.F.'})")

    def cylindrical_face(self, origin, axis_dir, x_dir, radius: float, v_lo: float, v_hi: float, same_sense=True) -> int:
        axis = self.axis2(origin, axis_dir, x_dir)
        unit_axis = _unit(axis_dir)
        bottom_center = _add(origin, _mul(unit_axis, v_lo))
        top_center = _add(origin, _mul(unit_axis, v_hi))
        return self._revolved_face(
            f"CYLINDRICAL_SURFACE('',#{axis},{_f(radius)})",
            bottom_center, top_center, axis_dir, x_dir, radius, radius, same_sense,
        )

    def conical_face(self, origin, axis_dir, x_dir, r_lo: float, r_hi: float, height: float, same_sense=True) -> int:
        semi_angle = math.atan2(r_hi - r_lo, height)
        axis = self.axis2(origin, axis_dir, x_dir)
        high_center = _add(origin, _mul(_unit(axis_dir), height))
        return self._revolved_face(
            f"CONICAL_SURFACE('',#{axis},{_f(r_lo)},{_f(semi_angle)})",
            origin, high_center, axis_dir, x_dir, r_lo, r_hi, same_sense,
        )

    # -- assembly ------------------------------------------------------------
    def solid(self, faces) -> int:
        refs = ",".join(f"#{face}" for face in faces)
        shell = self.raw(f"CLOSED_SHELL('',({refs}))")
        return self.raw(f"MANIFOLD_SOLID_BREP('{self.part_name}',#{shell})")

    def text(self, solid: int) -> str:
        application_context = self.raw("APPLICATION_CONTEXT('automotive design')")
        protocol = self.raw(
            "APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,"
            f"#{application_context})"
        )
        product_context = self.raw(f"PRODUCT_CONTEXT('',#{application_context},'mechanical')")
        product = self.raw(f"PRODUCT('{self.part_name}','{self.part_name}','',(#{product_context}))")
        self.raw(f"PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#{product}))")
        formation = self.raw(
            f"PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE('','',#{product},.NOT_KNOWN.)"
        )
        definition_context = self.raw(f"PRODUCT_DEFINITION_CONTEXT('part definition',#{protocol},'design')")
        definition = self.raw(f"PRODUCT_DEFINITION('design','',#{formation},#{definition_context})")
        shape = self.raw(f"PRODUCT_DEFINITION_SHAPE('','',#{definition})")
        representation = self.raw(
            f"ADVANCED_BREP_SHAPE_REPRESENTATION('{self.part_name}',(#{solid},#{self.context}),#{self.context})"
        )
        self.raw(f"SHAPE_DEFINITION_REPRESENTATION(#{shape},#{representation})")
        header = [
            "ISO-10303-21;",
            "HEADER;",
            "FILE_DESCRIPTION(('analytic B-rep sample for tensorcad.step_brep'),'2;1');",
            f"FILE_NAME('{self.part_name}.step','2026-11-13T00:00:00',('TensorCAD mission m17349'),"
            "('TensorCAD'),'tensorcad.tools.step_samples','TensorCAD','');",
            f"FILE_SCHEMA(('{SCHEMA}'));",
            "ENDSEC;",
            "DATA;",
        ]
        body = [f"#{index}={text};" for index, text in enumerate(self.entities, start=1)]
        return "\n".join(header + body + ["ENDSEC;", "END-ISO-10303-21;", ""])


# --------------------------------------------------------------------------- #
# the parts
# --------------------------------------------------------------------------- #
def plate_with_hole(out: Path) -> dict:
    """200 x 120 x 20 mm plate with a 40 mm through hole (plane + cylinder, one hole)."""
    x0, x1, y0, y1, z0, z1, radius = -100.0, 100.0, -60.0, 60.0, 0.0, 20.0, 20.0
    writer = SolidWriter("plate_with_hole")
    faces = [
        writer.planar_face(
            (0, 0, 1), (1, 0, 0),
            [(x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)],
            holes=[((0.0, 0.0, z1), radius)],
        ),
        writer.planar_face(
            (0, 0, -1), (1, 0, 0),
            [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0)],
            holes=[((0.0, 0.0, z0), radius)],
        ),
        writer.planar_face((-1, 0, 0), (0, 1, 0), [(x0, y0, z0), (x0, y1, z0), (x0, y1, z1), (x0, y0, z1)]),
        writer.planar_face((1, 0, 0), (0, 1, 0), [(x1, y0, z0), (x1, y1, z0), (x1, y1, z1), (x1, y0, z1)]),
        writer.planar_face((0, -1, 0), (1, 0, 0), [(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)]),
        writer.planar_face((0, 1, 0), (1, 0, 0), [(x0, y1, z0), (x1, y1, z0), (x1, y1, z1), (x0, y1, z1)]),
        writer.cylindrical_face((0.0, 0.0, 0.0), (0, 0, 1), (1, 0, 0), radius, z0, z1, same_sense=False),
    ]
    solid = writer.solid(faces)
    path = out / "plate_with_hole.step"
    path.write_text(writer.text(solid), encoding="utf-8")
    volume = (x1 - x0) * (y1 - y0) * (z1 - z0) - math.pi * radius**2 * (z1 - z0)
    return {
        "file": path.name,
        "nominal_volume_mm3": volume,
        "formula": "(200*120*20) - pi*20^2*20",
        "geometry": "plate 200x120x20 mm, through hole d40 mm",
        "entities": "PLANE + CYLINDRICAL_SURFACE, two faces with a FACE_BOUND hole",
    }


def stepped_shaft(out: Path) -> dict:
    """Cylinder d60 x 80, cone d60 -> d30 over 50, cylinder d30 x 40 (plane+cylinder+cone)."""
    r1, r2, h1, hc, h2 = 30.0, 15.0, 80.0, 50.0, 40.0
    writer = SolidWriter("stepped_shaft")
    faces = [
        writer.disc_face((0, 0, 0), (0, 0, -1), (1, 0, 0), r1),
        writer.disc_face((0, 0, h1 + hc + h2), (0, 0, 1), (1, 0, 0), r2),
        writer.cylindrical_face((0, 0, 0), (0, 0, 1), (1, 0, 0), r1, 0.0, h1),
        writer.conical_face((0, 0, h1), (0, 0, 1), (1, 0, 0), r1, r2, hc),
        writer.cylindrical_face((0, 0, h1 + hc), (0, 0, 1), (1, 0, 0), r2, 0.0, h2),
    ]
    solid = writer.solid(faces)
    path = out / "stepped_shaft.step"
    path.write_text(writer.text(solid), encoding="utf-8")
    volume = (
        math.pi * r1**2 * h1
        + math.pi * hc / 3.0 * (r1**2 + r1 * r2 + r2**2)
        + math.pi * r2**2 * h2
    )
    return {
        "file": path.name,
        "nominal_volume_mm3": volume,
        "formula": "pi*30^2*80 + pi*50/3*(30^2+30*15+15^2) + pi*15^2*40",
        "geometry": "d60x80 cylinder + d60->d30 cone (h=50) + d30x40 cylinder",
        "entities": "PLANE + CYLINDRICAL_SURFACE + CONICAL_SURFACE",
    }


def bushing(out: Path) -> dict:
    """Hollow cylinder d80 outside, d50 bore, 60 long (annular planar faces with holes)."""
    router, rinner, height = 40.0, 25.0, 60.0
    writer = SolidWriter("bushing")
    faces = [
        writer.cylindrical_face((0, 0, 0), (0, 0, 1), (1, 0, 0), router, 0.0, height),
        writer.cylindrical_face((0, 0, 0), (0, 0, 1), (1, 0, 0), rinner, 0.0, height, same_sense=False),
        writer.annulus_face((0, 0, height), (0, 0, 1), (1, 0, 0), router, rinner),
        writer.annulus_face((0, 0, 0), (0, 0, -1), (1, 0, 0), router, rinner),
    ]
    solid = writer.solid(faces)
    path = out / "bushing.step"
    path.write_text(writer.text(solid), encoding="utf-8")
    volume = math.pi * (router**2 - rinner**2) * height
    return {
        "file": path.name,
        "nominal_volume_mm3": volume,
        "formula": "pi*(40^2-25^2)*60",
        "geometry": "hollow cylinder d80/d50, length 60 mm",
        "entities": "PLANE with holes + two CYLINDRICAL_SURFACE (one same_sense=.F. bore)",
    }


def l_bracket(out: Path) -> dict:
    """L profile 60 x 50, leg width 20, extruded 10 mm: a non-convex planar trim."""
    profile = [(0.0, 0.0), (60.0, 0.0), (60.0, 20.0), (20.0, 20.0), (20.0, 50.0), (0.0, 50.0)]
    depth = 10.0
    writer = SolidWriter("l_bracket")
    top = [(x, y, depth) for x, y in profile]
    bottom = [(x, y, 0.0) for x, y in profile]
    faces = [
        writer.planar_face((0, 0, 1), (1, 0, 0), top),
        writer.planar_face((0, 0, -1), (1, 0, 0), bottom),
    ]
    for index in range(len(profile)):
        p0 = profile[index]
        p1 = profile[(index + 1) % len(profile)]
        edge = _sub((p1[0], p1[1], 0.0), (p0[0], p0[1], 0.0))
        normal = _unit((edge[1], -edge[0], 0.0))
        faces.append(
            writer.planar_face(
                normal,
                _unit(edge),
                [(p0[0], p0[1], 0.0), (p1[0], p1[1], 0.0), (p1[0], p1[1], depth), (p0[0], p0[1], depth)],
            )
        )
    solid = writer.solid(faces)
    path = out / "l_bracket.step"
    path.write_text(writer.text(solid), encoding="utf-8")
    area = 0.0
    for index in range(len(profile)):
        x1, y1 = profile[index]
        x2, y2 = profile[(index + 1) % len(profile)]
        area += x1 * y2 - x2 * y1
    area = abs(area) / 2.0
    return {
        "file": path.name,
        "nominal_volume_mm3": area * depth,
        "formula": "shoelace(profile 60x50, leg 20) * 10",
        "geometry": "L profile 60x50, leg width 20 mm, extruded 10 mm",
        "entities": "PLANE only; two non-convex (6-vertex) trims exercise ear clipping",
    }


BUILDERS = (plate_with_hole, stepped_shaft, bushing, l_bracket)


def main() -> int:
    parser = argparse.ArgumentParser(description="Write AP214 analytic B-rep sample solids")
    parser.add_argument("--out", default="artifacts/gap-cad-b1-step-reader/samples")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    records = []
    for builder in BUILDERS:
        record = builder(out)
        record["file_bytes"] = (out / record["file"]).stat().st_size
        records.append(record)
    payload = {
        "note": "nominal volumes are closed-form for the geometry defined in tools/step_samples.py",
        "units": "millimetres in the file (SI_UNIT(.MILLI.,.METRE.)); volumes converted to m^3 by the reader",
        "samples": records,
    }
    (out / "nominal_volumes.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())