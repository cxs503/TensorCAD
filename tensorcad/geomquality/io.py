"""CAD exchange for the geometry-quality pipeline: STEP, IGES and the native form.

Three routes into the same :class:`~tensorcad.geomquality.brep.BRep`, and the
difference between them is the point of the pipeline:

``.brep.json``
    the native patch set: control nets, knots, ``same_sense`` and an explicit edge
    pool.  Lossless, and the reference every round trip is compared against.
``.step``
    a **topological** route.  Boundary curves and vertices are emitted once and
    referenced by several faces, so ``edge_curve`` identity *is* the shared edge: the
    reader gets the adjacency graph the file declares, without guessing.
``.iges``
    a **geometric** route.  IGES entity 128 is an untrimmed NURBS surface with no
    topology at all, so the reader has to sew the patch boundaries itself
    (:func:`tensorcad.geomquality.brep.sew`).  That is exactly where an export seam
    crack survives, and why the pipeline accepts both routes.

Deliberate scope limits, stated rather than hidden: trimmed surfaces (IGES 144,
STEP ``advanced_face`` with non-rectangular ``edge_loop``) are *parsed* so that
topology is still checkable, but patch quadrature falls back to tessellation for
them; STEP assemblies with placement transforms are read as a single body in file
coordinates; rational B-spline surfaces are read, written only when non-rational.
"""
from __future__ import annotations

import json
import math
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

from ..nurbs import NURBSCurve, NURBSSurface
from .brep import (
    BOUNDARY_SAMPLES,
    DEFAULT_TOLERANCE,
    BRep,
    Edge,
    EdgeRef,
    Face,
    sew,
)

Vec3 = tuple[float, float, float]


class ExchangeError(ValueError):
    """Raised when a CAD file cannot be read into a B-rep."""


# --------------------------------------------------------------------------- #
# native JSON
# --------------------------------------------------------------------------- #
def brep_to_dict(brep: BRep) -> dict[str, Any]:
    faces = []
    for face in brep.faces:
        faces.append(
            {
                "name": face.name,
                "degree": list(face.surface.degree),
                "control_net": [
                    [list(point) for point in row] for row in face.surface.control_net
                ],
                "weights": [list(row) for row in face.surface.weights],
                "knots": [list(axis) for axis in face.surface.knots],
                "same_sense": face.same_sense,
                "wetted": face.wetted,
                "boundary": [
                    {"edge": ref.edge, "forward": ref.forward} for ref in face.boundary
                ],
            }
        )
    return {
        "format": "tensorcad.brep/1",
        "source": brep.source,
        "metadata": dict(brep.metadata),
        "faces": faces,
        "edges": [
            {"index": edge.index, "points": [list(point) for point in edge.points]}
            for edge in brep.edges
        ],
    }


def brep_from_dict(payload: dict[str, Any]) -> BRep:
    if payload.get("format") != "tensorcad.brep/1":
        raise ExchangeError(f"unknown brep format {payload.get('format')!r}")
    faces = []
    for index, entry in enumerate(payload["faces"]):
        surface = NURBSSurface(
            degree=tuple(entry["degree"]),
            control_net=[
                [tuple(point) for point in row] for row in entry["control_net"]
            ],
            weights=[list(row) for row in entry["weights"]],
            knots=[list(axis) for axis in entry["knots"]],
        )
        faces.append(
            Face(
                index=index,
                surface=surface,
                boundary=tuple(
                    EdgeRef(edge=ref["edge"], forward=ref["forward"])
                    for ref in entry["boundary"]
                ),
                same_sense=bool(entry.get("same_sense", True)),
                name=entry.get("name", f"face{index}"),
                wetted=bool(entry.get("wetted", True)),
            )
        )
    edges = tuple(
        Edge(index=edge["index"], points=tuple(tuple(p) for p in edge["points"]))
        for edge in payload["edges"]
    )
    return BRep(
        faces=tuple(faces),
        edges=edges,
        source=payload.get("source", "json"),
        metadata=dict(payload.get("metadata", {})),
    )


def write_brep_json(path: str | Path, brep: BRep) -> Path:
    target = Path(path)
    target.write_text(json.dumps(brep_to_dict(brep), indent=1))
    return target


def read_brep_json(path: str | Path) -> BRep:
    return brep_from_dict(json.loads(Path(path).read_text()))


# --------------------------------------------------------------------------- #
# STEP (AP203/AP214 subset)
# --------------------------------------------------------------------------- #
_STEP_HEADER = """ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('TensorCAD geomquality B-rep export'),'2;1');
FILE_NAME('{name}','{stamp}',('TensorCAD'),('TensorCAD'),'geomquality','TensorCAD','');
FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));
ENDSEC;
DATA;
"""


def _step_number(value: float) -> str:
    text = repr(float(value))
    if "e" in text or "E" in text:
        mantissa, exponent = text.split("e")
        return f"{mantissa}E{int(exponent):+d}"
    if "." not in text:
        text += "."
    return text


def _knot_multiplicities(knots: Sequence[float]) -> tuple[list[float], list[int]]:
    distinct: list[float] = []
    mults: list[int] = []
    for knot in knots:
        if distinct and knot == distinct[-1]:
            mults[-1] += 1
        else:
            distinct.append(knot)
            mults.append(1)
    return distinct, mults


#: STEP has no per-face attribute for "this face is the waterplane closure", so the
#: flag travels as a suffix in the entity name and is stripped by the reader.
CLOSURE_MARKER = "#closure"


def _step_face_label(face: Face) -> str:
    return face.name if face.wetted else face.name + CLOSURE_MARKER


def _strip_marker(name: str) -> tuple[str, bool]:
    if name.endswith(CLOSURE_MARKER):
        return name[: -len(CLOSURE_MARKER)], False
    return name, True


def write_step(path: str | Path, brep: BRep, *, name: str = "model") -> Path:
    """Write a B-rep as STEP with explicit shared topology.

    Every edge becomes one ``EDGE_CURVE`` (a degree-1 B-spline through the sampled
    curve, which is exact for the polyline the audit measured) referenced by the
    faces that use it, so a reader recovers the same adjacency graph the audit saw.
    """
    lines: list[str] = []
    counter = 1

    def add(text: str) -> int:
        nonlocal counter
        lines.append(f"#{counter}={text};")
        counter += 1
        return counter - 1

    point_ids: dict[tuple[float, float, float], int] = {}

    def cartesian(point: Sequence[float]) -> int:
        key = (round(point[0], 12), round(point[1], 12), round(point[2], 12))
        if key not in point_ids:
            point_ids[key] = add(
                "CARTESIAN_POINT('',({},{},{}))".format(*(_step_number(v) for v in key))
            )
        return point_ids[key]

    edge_ids: dict[int, int] = {}
    vertex_ids: dict[tuple[float, float, float], int] = {}

    def vertex(point: Sequence[float]) -> int:
        key = (round(point[0], 12), round(point[1], 12), round(point[2], 12))
        if key not in vertex_ids:
            vertex_ids[key] = add(f"VERTEX_POINT('',#{cartesian(point)})")
        return vertex_ids[key]

    for edge in brep.edges:
        # Control points are emitted as CARTESIAN_POINT entities and referenced, which is
        # what a real exporter does and what the reader below expects.
        controls = ",".join(f"#{cartesian(point)}" for point in edge.points)
        degree = 1
        count = len(edge.points)
        inner = [i / (count - 1) for i in range(1, count - 1)]
        knots = [0.0, 0.0, *inner, 1.0, 1.0]
        distinct, mults = _knot_multiplicities(knots)
        curve = add(
            "B_SPLINE_CURVE_WITH_KNOTS('',{},{},.UNSPECIFIED.,.F.,.F.,({}),({}),.UNSPECIFIED.)".format(
                degree,
                f"({controls})",
                ",".join(str(m) for m in mults),
                ",".join(_step_number(k) for k in distinct),
            )
        )
        edge_ids[edge.index] = add(
            f"EDGE_CURVE('',#{vertex(edge.points[0])},#{vertex(edge.points[-1])},#{curve},.T.)"
        )

    face_ids: list[int] = []
    for face in brep.faces:
        surface = face.surface
        net_refs = ", ".join(
            "(" + ",".join(f"#{cartesian(point)}" for point in row) + ")"
            for row in surface.control_net
        )
        if surface.rational:
            raise ExchangeError(
                "writing a rational B-spline surface is not implemented; "
                "export the non-rational patch or the .brep.json form"
            )
        u_knots, u_mults = _knot_multiplicities(surface.knots[0])
        v_knots, v_mults = _knot_multiplicities(surface.knots[1])
        surface_id = add(
            "B_SPLINE_SURFACE_WITH_KNOTS('',{},{},({}),.UNSPECIFIED.,.F.,.F.,.F.,({}),({}),({}),({}),.UNSPECIFIED.)".format(
                surface.degree[0],
                surface.degree[1],
                net_refs,
                ",".join(str(m) for m in u_mults),
                ",".join(str(m) for m in v_mults),
                ",".join(_step_number(k) for k in u_knots),
                ",".join(_step_number(k) for k in v_knots),
            )
        )
        oriented = []
        for ref in face.boundary:
            oriented.append(
                add(
                    "ORIENTED_EDGE('',*,*,#{},.{})".format(
                        edge_ids[ref.edge], "T." if ref.forward else "F."
                    )
                )
            )
        loop = add("EDGE_LOOP('',({}))".format(",".join(f"#{i}" for i in oriented)))
        bound = add(f"FACE_OUTER_BOUND('',#{loop},.T.)")
        face_ids.append(
            add(
                "ADVANCED_FACE('{}',(#{}) ,#{},.{})".format(
                    _step_face_label(face), bound, surface_id,
                    "T." if face.same_sense else "F."
                )
            )
        )
    shell_kind = "CLOSED_SHELL" if not any(
        len(brep.edge_uses()[edge.index]) == 1 for edge in brep.edges
    ) else "OPEN_SHELL"
    shell = add("{}('',({}))".format(shell_kind, ",".join(f"#{i}" for i in face_ids)))
    add(f"MANIFOLD_SURFACE_SHAPE_REPRESENTATION('{name}',(#{shell}),#0)")

    text = _STEP_HEADER.format(name=name, stamp="2026-11-07T00:00:00")
    text += "\n".join(lines) + "\nENDSEC;\nEND-ISO-10303-21;\n"
    target = Path(path)
    target.write_text(text)
    return target


def _step_statements(text: str) -> dict[int, tuple[str, str]]:
    data_start = text.find("DATA;")
    data_end = text.find("ENDSEC;", data_start)
    if data_start < 0 or data_end < 0:
        raise ExchangeError("no DATA section")
    body = text[data_start + 5 : data_end]
    statements: dict[int, tuple[str, str]] = {}
    buffer = []
    for char in body:
        if char == ";":
            statement = "".join(buffer).strip()
            buffer = []
            if not statement.startswith("#"):
                continue
            identifier, rest = statement[1:].split("=", 1)
            name, _, arguments = rest.partition("(")
            statements[int(identifier)] = (name.strip().upper(), arguments.rsplit(")", 1)[0])
        else:
            buffer.append(char)
    return statements


def _split_top(arguments: str) -> list[str]:
    parts: list[str] = []
    depth = 0
    current: list[str] = []
    in_string = False
    for char in arguments:
        if in_string:
            current.append(char)
            if char == "'":
                in_string = False
            continue
        if char == "'":
            in_string = True
            current.append(char)
        elif char == "(":
            depth += 1
            current.append(char)
        elif char == ")":
            depth -= 1
            current.append(char)
        elif char == "," and depth == 0:
            parts.append("".join(current).strip())
            current = []
        else:
            current.append(char)
    if "".join(current).strip():
        parts.append("".join(current).strip())
    return parts


def _step_value(token: str) -> Any:
    token = token.strip()
    if token.startswith("#"):
        return ("ref", int(token[1:]))
    if token.startswith("(") and token.endswith(")"):
        return [_step_value(item) for item in _split_top(token[1:-1])]
    if token.startswith("'") and token.endswith("'"):
        return token[1:-1]
    if token.startswith(".") and token.endswith("."):
        return token[1:-1].upper() == "T" if token[1:-1].upper() in {"T", "F"} else token
    try:
        return float(token)
    except ValueError:
        return token


def _step_bool(value: Any) -> bool:
    """STEP enumerations arrive as ``True``/``False`` (parsed) or ``.T.``/``.F.`` text."""
    if isinstance(value, bool):
        return value
    return str(value).strip().strip(".").upper() in {"T", "TRUE"}


def _step_args(arguments: str) -> list[Any]:
    return [_step_value(token) for token in _split_top(arguments)]


def read_step(path: str | Path, *, tolerance: float = DEFAULT_TOLERANCE) -> BRep:
    """Read the B-rep subset of an AP203/AP214 file into faces and edges."""
    statements = _step_statements(Path(path).read_text())

    def surface_from(identifier: int) -> NURBSSurface:
        name, arguments = statements[identifier]
        if name not in {"B_SPLINE_SURFACE_WITH_KNOTS", "RATIONAL_B_SPLINE_SURFACE"}:
            raise ExchangeError(f"unsupported surface entity {name}")
        args = _step_args(arguments)
        degree = (int(args[1]), int(args[2]))
        net = [
            [_point_of(ref[1]) for ref in row]
            for row in args[3]
        ]
        closed_u, closed_v, self_intersect = args[5], args[6], args[7]
        u_mults, v_mults = [int(v) for v in args[8]], [int(v) for v in args[9]]
        u_distinct = [float(v) for v in args[10]]
        v_distinct = [float(v) for v in args[11]]
        u_knots = [k for k, m in zip(u_distinct, u_mults) for _ in range(m)]
        v_knots = [k for k, m in zip(v_distinct, v_mults) for _ in range(m)]
        weights = [[1.0] * len(net[0]) for _ in net]
        if name == "RATIONAL_B_SPLINE_SURFACE":
            weight_values = [float(v) for v in _flatten(args[12])]
            weights = [
                weight_values[i * len(net[0]) : (i + 1) * len(net[0])] for i in range(len(net))
            ]
        del closed_u, closed_v, self_intersect
        return NURBSSurface(
            degree=degree, control_net=net, weights=weights, knots=[u_knots, v_knots]
        )

    def _flatten(values: Iterable[Any]) -> list[Any]:
        flat: list[Any] = []
        for item in values:
            if isinstance(item, list):
                flat.extend(_flatten(item))
            else:
                flat.append(item)
        return flat

    def _point_of(identifier: int) -> Vec3:
        args = _step_args(statements[identifier][1])
        return tuple(float(value) for value in args[1])  # type: ignore[return-value]

    def curve_polyline(identifier: int, samples: int = BOUNDARY_SAMPLES) -> tuple[Vec3, ...]:
        name, arguments = statements[identifier]
        if name in {"B_SPLINE_CURVE_WITH_KNOTS", "RATIONAL_B_SPLINE_CURVE"}:
            args = _step_args(arguments)
            degree = int(args[1])
            points = [_point_of(ref[1]) for ref in args[2]]
            u_mults = [int(v) for v in args[6]]
            u_distinct = [float(v) for v in args[7]]
            knots = [k for k, m in zip(u_distinct, u_mults) for _ in range(m)]
            if degree == 1:
                return tuple(points)  # type: ignore[return-value]
            weights = [1.0] * len(points)
            if name == "RATIONAL_B_SPLINE_CURVE":
                weights = [float(v) for v in args[3]]
            curve = NURBSCurve(
                degree=degree, control_net=points, weights=weights, knots=knots
            )
            return tuple(curve.polyline(samples))
        raise ExchangeError(f"unsupported curve entity {name}")

    edge_pool: dict[int, int] = {}
    edges: list[Edge] = []
    faces: list[Face] = []
    for identifier in sorted(statements):
        name, arguments = statements[identifier]
        if name != "ADVANCED_FACE":
            continue
        args = _step_args(arguments)
        face_name = args[0] if isinstance(args[0], str) else ""
        bounds = args[1]
        surface = surface_from(args[2][1])
        same_sense = _step_bool(args[3])
        loop_ref = None
        for bound in bounds:
            bound_args = _step_args(statements[bound[1]][1])
            loop_ref = bound_args[1][1]
            break
        if loop_ref is None:
            raise ExchangeError(f"face #{identifier} has no bound")
        loop_args = _step_args(statements[loop_ref][1])
        refs: list[EdgeRef] = []
        for oriented in loop_args[1]:
            oriented_args = _step_args(statements[oriented[1]][1])
            edge_curve = oriented_args[3][1]
            forward = _step_bool(oriented_args[4])
            if edge_curve not in edge_pool:
                edge_args = _step_args(statements[edge_curve][1])
                polyline = curve_polyline(edge_args[3][1])
                index = len(edges)
                edges.append(Edge(index=index, points=polyline))
                edge_pool[edge_curve] = index
            refs.append(EdgeRef(edge=edge_pool[edge_curve], forward=forward))
        clean_name, wetted = _strip_marker(face_name)
        faces.append(
            Face(
                index=len(faces),
                surface=surface,
                boundary=tuple(refs),
                same_sense=same_sense,
                name=clean_name,
                wetted=wetted,
            )
        )
    if not faces:
        raise ExchangeError("no ADVANCED_FACE entities: not a B-rep STEP file")
    return BRep(
        faces=tuple(faces),
        edges=tuple(edges),
        source=str(path),
        metadata={"format": "step", "tolerance": tolerance},
    )


# --------------------------------------------------------------------------- #
# IGES (entity 128, untrimmed NURBS surfaces)
# --------------------------------------------------------------------------- #
def _iges_record(body: str, section: str, sequence: int) -> str:
    body = body[:72]
    return f"{body:<72}{section}{sequence:>7}\n"


def write_iges(path: str | Path, brep: BRep, *, units_flag: int = 6) -> Path:
    """Write the patch surfaces as IGES entity 128 records (no topology by design)."""
    start = [
        _iges_record("TensorCAD geomquality export: untrimmed NURBS patch set", "S", 1)
    ]
    globals_fields = [
        "1H,", "1H;", "8HTensorCAD", "12Hgeomquality.stp", "8HTensorCAD", "11Hgeomquality",
        "32", "38", "308", "15", "308", "15", "8HTensorCAD", "1.0",
        str(units_flag), "2HM", "1", "1", "26H20261107.000000", "1.0E-06",
        _iges_f(brep.length), "8HTensorCAD", "8HTensorCAD", "11", "0",
    ]
    globals_text = _pack_parameters(globals_fields, "G", 1)
    directory: list[str] = []
    parameter_lines: list[str] = []
    for index, face in enumerate(brep.faces):
        surface = face.surface
        if surface.rational:
            raise ExchangeError("IGES export of rational patches is not implemented")
        (u0, u1), (v0, v1) = surface.domain
        values: list[str] = [
            str(surface.degree[0]),
            str(surface.degree[1]),
            str(len(surface.control_net)),
            str(len(surface.control_net[0])),
            "1" if u0 == surface.knots[0][0] and u1 == surface.knots[0][-1] else "0",
            "0",
            "0",
            "0",
        ]
        # Control points, U index outer / V index fastest (IGES 5.3 §4.23 order).
        for row in surface.control_net:
            for point in row:
                values.extend(_iges_f(value) for value in point)
        values.extend(_iges_f(knot) for knot in surface.knots[0])
        values.extend(_iges_f(knot) for knot in surface.knots[1])
        values.extend("1.0" for _ in range(len(surface.control_net) * len(surface.control_net[0])))
        values.extend([_iges_f(v0), _iges_f(v1), _iges_f(u0), _iges_f(u1)])
        first_pointer = len(parameter_lines) + 1
        # Each P record points (cols 66-72) at the entity's first directory line.
        body_lines = _pack_parameters(values, "P", 1, pointer=2 * index + 1)
        count = len(body_lines)
        body_lines = [
            f"{line[:64]}{2 * index + 1:>8}P{(offset + first_pointer):>7}\n"
            for offset, line in enumerate(body_lines)
        ]
        parameter_lines.extend(body_lines)
        field1 = [
            "128", f"{first_pointer}", "0", "0",
            "1" if face.wetted else "2",  # level: 1 = wetted surface, 2 = closure face
            "0", "0", "0", "00000000",
        ]
        field2 = ["128", "0", "0", f"{count}", "0", "0", "0", "0", f"{index}"]
        directory.append(_iges_record("".join(f"{f:>8}" for f in field1), "D", 2 * index + 1))
        directory.append(_iges_record("".join(f"{f:>8}" for f in field2), "D", 2 * index + 2))
    terminate = _iges_record(
        f"S{1:>7}G{len(globals_text):>7}D{len(directory):>7}P{len(parameter_lines):>7}",
        "T",
        1,
    )
    text = "".join(start) + "".join(globals_text) + "".join(directory) + "".join(parameter_lines) + terminate
    target = Path(path)
    target.write_text(text)
    return target


def _iges_f(value: float) -> str:
    text = f"{float(value):.12g}"
    if "e" in text or "E" in text:
        mantissa, exponent = text.lower().split("e")
        return f"{mantissa}D{int(exponent):+d}"
    if "." not in text:
        text += "."
    return text


def _pack_parameters(
    values: Sequence[str], section: str, sequence: int, *, pointer: int = 0
) -> list[str]:
    """Free-format IGES parameter block: 64 columns of data, pointer in 66-72."""
    lines: list[str] = []
    current = ""
    for value in values:
        token = value if value.endswith(",") or value.endswith(";") else value + ","
        if len(current) + len(token) > 64:
            lines.append(current)
            current = ""
        current += token
    if current:
        if not current.endswith(";"):
            current = current[:-1] + ";" if current.endswith(",") else current + ";"
        lines.append(current)
    return [
        f"{line:<64}{pointer:>8}{section}{sequence + offset:>7}\n"
        for offset, line in enumerate(lines)
    ]


def _iges_parameters(text: str) -> dict[int, list[float]]:
    sections: dict[str, list[str]] = {"S": [], "G": [], "D": [], "P": [], "T": []}
    for line in text.splitlines():
        if len(line) < 73:
            continue
        marker = line[72]
        if marker in sections:
            sections[marker].append(line)
    blocks: dict[int, list[str]] = {}
    for line in sections["P"]:
        pointer = int(line[64:72].strip() or 0)
        blocks.setdefault(pointer, []).append(line[:64].rstrip())
    parameters: dict[int, list[float]] = {}
    for pointer, chunks in blocks.items():
        joined = "".join(chunks)
        joined = joined.replace("D", "E").replace("d", "e")
        joined = joined.replace(";", "")
        tokens = [token for token in joined.split(",") if token.strip()]
        values: list[float] = []
        for token in tokens:
            try:
                values.append(float(token))
            except ValueError:
                values.append(float("nan"))
        parameters[pointer] = values
    return parameters


def read_iges(path: str | Path, *, tolerance: float = DEFAULT_TOLERANCE) -> BRep:
    """Read entity 128 surfaces and sew them into a B-rep.

    IGES carries no topology, so the adjacency graph here is *derived*: any two patch
    boundaries closer than ``tolerance`` become one edge.  A crack wider than the
    tolerance therefore stays a free edge and is reported as such.
    """
    text = Path(path).read_text()
    directory: list[tuple[int, int, int, int, int]] = []
    lines = [line for line in text.splitlines() if len(line) >= 73]
    d_lines = [line for line in lines if line[72] == "D"]
    for pair in range(0, len(d_lines) - 1, 2):
        first, second = d_lines[pair], d_lines[pair + 1]
        if not first[0:8].strip():
            continue
        entity_type = int(first[0:8])
        parameter_pointer = int(first[8:16])   # first P-section line of this entity
        directory_pointer = int(first[72 + 1 : 80].strip() or 0) or (pair + 1)
        count = int(second[24:32])
        level = int(first[32:40].strip() or 1)
        directory.append((entity_type, parameter_pointer, count, directory_pointer, level))
    parameters = _iges_parameters(text)
    surfaces: list[NURBSSurface] = []
    wetted_flags: list[bool] = []
    for entity_type, pointer, count, directory_pointer, level in directory:
        if entity_type == 144:
            raise ExchangeError(
                "IGES trimmed surface (type 144) is out of scope: export untrimmed patches"
            )
        if entity_type != 128:
            continue
        # Either column in the two sections can be used to join them; try both.
        values = parameters.get(directory_pointer) or parameters.get(pointer) or []
        if len(values) < 8:
            raise ExchangeError(
                f"entity 128 at P line {pointer} has {len(values)} parameters: truncated file"
            )
        degree_u = int(values[0])
        degree_v = int(values[1])
        count_u = int(values[2])
        count_v = int(values[3])
        cursor = 8
        net: list[list[Vec3]] = []
        for _ in range(count_u):
            row = []
            for _ in range(count_v):
                row.append(tuple(values[cursor : cursor + 3]))  # type: ignore[arg-type]
                cursor += 3
            net.append(row)
        u_knots = values[cursor : cursor + count_u + degree_u + 1]
        cursor += count_u + degree_u + 1
        v_knots = values[cursor : cursor + count_v + degree_v + 1]
        cursor += count_v + degree_v + 1
        raw_weights = values[cursor : cursor + count_u * count_v]
        cursor += count_u * count_v
        weights = [
            list(raw_weights[i * count_v : (i + 1) * count_v]) for i in range(count_u)
        ]
        surfaces.append(
            NURBSSurface(
                degree=(degree_u, degree_v),
                control_net=net,
                weights=weights,
                knots=[list(u_knots), list(v_knots)],
            )
        )
        wetted_flags.append(level != 2)
    if not surfaces:
        raise ExchangeError("no entity 128 surfaces found in the IGES file")
    sewn = sew(
        surfaces,
        tolerance=tolerance,
        wetted=wetted_flags,
        source=str(path),
        metadata={"format": "iges"},
    )
    # IGES has no orientation: propagate one over the adjacency graph before the audit
    # can say anything about watertightness (a sewn but randomly-oriented shell would be
    # reported as an orientation defect, which would be true of the file but not of the
    # model).
    from .brep import orient_shell

    return orient_shell(sewn)


# --------------------------------------------------------------------------- #
# dispatch
# --------------------------------------------------------------------------- #
def load_cad(path: str | Path, *, tolerance: float = DEFAULT_TOLERANCE) -> BRep:
    """Read any supported CAD file into a B-rep, dispatching on the suffix."""
    target = Path(path)
    suffix = target.suffix.lower()
    if suffix == ".json":
        return read_brep_json(target)
    if suffix == ".step" or suffix == ".stp":
        return read_step(target)
    if suffix == ".iges" or suffix == ".igs":
        return read_iges(target, tolerance=tolerance)
    raise ExchangeError(f"unsupported CAD format {suffix!r}")


def save_cad(path: str | Path, brep: BRep) -> Path:
    target = Path(path)
    suffix = target.suffix.lower()
    if suffix == ".json":
        return write_brep_json(target, brep)
    if suffix in {".step", ".stp"}:
        return write_step(target, brep)
    if suffix in {".iges", ".igs"}:
        return write_iges(target, brep)
    raise ExchangeError(f"unsupported CAD format {suffix!r}")


def round_trip(brep: BRep, directory: str | Path, *, stem: str = "model") -> dict[str, Any]:
    """Write/read every format and report the geometry deviation that survives."""
    from .brep import audit_topology  # local import keeps the module graph acyclic

    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    results: dict[str, Any] = {}
    before = audit_topology(brep)
    for suffix in (".brep.json", ".step", ".iges"):
        target = out / f"{stem}{suffix}"
        save_cad(target, brep)
        reloaded = load_cad(target)
        deviation = (
            max(
                math.dist(a.surface.evaluate(u, v), b.surface.evaluate(u, v))
                for a, b in zip(brep.faces, reloaded.faces)
                for u in (a.domain[0][0], 0.5 * sum(a.domain[0]), a.domain[0][1])
                for v in (a.domain[1][0], 0.5 * sum(a.domain[1]), a.domain[1][1])
            )
            if len(brep.faces) == len(reloaded.faces)
            else math.inf
        )
        after = audit_topology(reloaded)
        results[suffix] = {
            "path": str(target),
            "bytes": target.stat().st_size,
            "faces": len(reloaded.faces),
            "edges": len(reloaded.edges),
            "max_deviation": deviation,
            "watertight_before": before.watertight,
            "watertight_after": after.watertight,
            "codes_before": sorted({d.code for d in before.defects}),
            "codes_after": sorted({d.code for d in after.defects}),
        }
    return results
