"""Geometry exchange seam for the quality audit pipeline.

The audit pipeline in :mod:`tensorcad.geometry_audit` judges a *triangle mesh*.  A
shipyard hands it something else: a STEP or IGES file.  This module is the seam,
and like :mod:`tensorcad.adapters.ocp` it is deliberately narrow and honest about
what it can and cannot do:

===========================================  ==============================  ============
route                                        status here                     evidence
===========================================  ==============================  ============
faceted STEP (``FACETED_BREP``/``CLOSED_SHELL``)  written and read here      REAL
neutral ``.patchset.json`` (patch table + facet triangles)  read here         REAL
analytic B-rep STEP (``ADVANCED_FACE``, ``CYLINDRICAL_SURFACE`` ...)  refused  NOT_SUPPORTED
IGES                                                      refused             NOT_SUPPORTED
===========================================  ==============================  ============

What "refused" means: :class:`ExchangeError` is raised and the message names the
entities that block the import.  It never guesses a tessellation.  A silently
mis-tessellated B-rep is the worst possible input to a watertightness audit,
because the audit would then report defects of the tessellator and not of the
model.  Tessellating a real B-rep needs OpenCASCADE (``pip install
"tensorcad[ocp]"``) and that route lives in :mod:`tensorcad.adapters.ocp`.

The faceted STEP writer emits AP203-style ``FACETED_BREP`` solids: every facet is
a planar ``FACE`` bounded by a single ``POLY_LOOP`` over ``VERTEX_POINT``s.  That
is exactly what :func:`tensorcad.adapters.ocp.write_step` produces by sewing a
triangle soup, so the two writers describe the same kind of artifact.  No
third-party CAD kernel is available in this environment, so the file is validated
by an internal round trip (:func:`round_trip_facets`) and **not** by a foreign
reader; the module says so rather than implying interoperability it has not
tested.

Units: STEP files carry a millimetre ``LENGTH_UNIT`` by convention, so
:func:`write_step_facets` multiplies by ``length_scale`` (default 1000) on the way
out and :func:`read_step_facets` divides it back out.  The scale travels with the
data instead of being assumed.
"""
from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

from .geometry import Mesh

__all__ = [
    "ExchangeError",
    "ExchangeUnavailable",
    "PatchGeometry",
    "read_step_facets",
    "read_patchset",
    "load_exchange",
    "round_trip_facets",
    "write_step_facets",
    "write_patchset",
    "STEP_FACETS_SCHEMA",
    "PATCHSET_SCHEMA",
]

STEP_FACETS_SCHEMA = "tensorcad.exchange.step_facets/1"
PATCHSET_SCHEMA = "tensorcad.geometry_audit.patchset/1"

#: STEP is written in millimetres by convention.
DEFAULT_LENGTH_SCALE = 1000.0

#: Entities that mean "this is an analytic B-rep, not a faceted solid". Their
#: presence is what makes a refusal more useful than a guess.
_BREP_ENTITIES = (
    "ADVANCED_FACE",
    "B_SPLINE_SURFACE",
    "B_SPLINE_SURFACE_WITH_KNOTS",
    "CYLINDRICAL_SURFACE",
    "CONICAL_SURFACE",
    "SPHERICAL_SURFACE",
    "TOROIDAL_SURFACE",
    "SURFACE_OF_REVOLUTION",
    "SURFACE_OF_LINEAR_EXTRUSION",
    "RATIONAL_B_SPLINE_SURFACE",
)


class ExchangeError(ValueError):
    """Raised when an exchange file cannot be turned into a triangle mesh."""


class ExchangeUnavailable(RuntimeError):
    """Raised when a format needs an optional dependency that is not installed."""


@dataclass(frozen=True)
class PatchGeometry:
    """A mesh plus the patch (face-group) each triangle came from.

    The audit needs the grouping because a *seam gap* between two B-rep faces is a
    defect that a single welded triangle soup cannot express: the triangles are
    fine, the faces simply do not meet.  ``patch_ids[i]`` is the patch of
    ``mesh.faces[i]``.
    """

    mesh: Mesh
    patch_ids: tuple[str, ...]
    source: str = "unknown"
    #: Anything the importer measured about the file (entity counts, scale, ...).
    metadata: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        if len(self.patch_ids) != len(self.mesh.faces):
            raise ExchangeError(
                f"patch_ids has {len(self.patch_ids)} entries for {len(self.mesh.faces)} faces"
            )

    def patch_of_face(self, index: int) -> str:
        return self.patch_ids[index]

    def patch_names(self) -> tuple[str, ...]:
        seen: list[str] = []
        for name in self.patch_ids:
            if name not in seen:
                seen.append(name)
        return tuple(seen)


# --------------------------------------------------------------------------- #
# patchset JSON: the neutral, dependency-free carrier used by the audit tests
# --------------------------------------------------------------------------- #
def write_patchset(path: str | Path, geometry: PatchGeometry) -> Path:
    """Write a patch table plus facets as JSON (round-trippable, human readable)."""
    target = Path(path)
    payload = {
        "schema": PATCHSET_SCHEMA,
        "units": geometry.mesh.units,
        "coordinate_frame": geometry.mesh.coordinate_frame,
        "source": geometry.source,
        "metadata": geometry.metadata or {},
        "vertices": [list(vertex) for vertex in geometry.mesh.vertices],
        "faces": [list(face) for face in geometry.mesh.faces],
        "patch_ids": list(geometry.patch_ids),
    }
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    return target


def read_patchset(path: str | Path) -> PatchGeometry:
    """Read a ``.patchset.json`` file written by :func:`write_patchset`."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    schema = data.get("schema")
    if schema != PATCHSET_SCHEMA:
        raise ExchangeError(f"unexpected patchset schema {schema!r}, want {PATCHSET_SCHEMA!r}")
    vertices = tuple(tuple(float(value) for value in vertex) for vertex in data["vertices"])
    faces = tuple(tuple(int(index) for index in face) for face in data["faces"])
    mesh = Mesh(
        vertices=vertices,
        faces=faces,
        units=str(data.get("units", "m")),
        coordinate_frame=str(data.get("coordinate_frame", "x aft->fore, y port, z up")),
    )
    return PatchGeometry(
        mesh=mesh,
        patch_ids=tuple(str(name) for name in data["patch_ids"]),
        source=str(data.get("source", str(path))),
        metadata=dict(data.get("metadata") or {}),
    )


# --------------------------------------------------------------------------- #
# faceted STEP (AP203 FACETED_BREP over POLY_LOOP)
# --------------------------------------------------------------------------- #
def write_step_facets(
    path: str | Path,
    mesh: Mesh,
    *,
    length_scale: float = DEFAULT_LENGTH_SCALE,
    name: str = "tensorcad_facets",
) -> Path:
    """Write ``mesh`` as an AP203 faceted solid (``FACETED_BREP`` over ``POLY_LOOP``).

    The file is a plain faceted representation: no analytic surfaces, no sewing,
    one planar face per triangle.  Coordinates are multiplied by ``length_scale``
    and the millimetre ``LENGTH_UNIT`` context says so.
    """
    if not mesh.faces:
        raise ExchangeError("cannot write a STEP file for a mesh with no faces")
    if length_scale <= 0:
        raise ExchangeError(f"length_scale must be positive, got {length_scale!r}")

    lines: list[str] = []
    counter = 0

    def emit(body: str) -> str:
        nonlocal counter
        counter += 1
        entity = f"#{counter}"
        lines.append(f"{entity} = {body};")
        return entity

    def number(value: float) -> str:
        # STEP wants a real (``1.``, ``-31.25``); appending a second decimal point to
        # ``-31.25`` produces ``-31.25.``, which reads back as an invalid float.
        text = f"{value * length_scale:.17g}"
        if "." not in text and "e" not in text.lower():
            text += "."
        return text

    unit = emit("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) )")
    angle = emit("( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )")
    solid_angle = emit("( NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT() )")
    uncertainty = emit(f"UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-06),{unit},'distance_accuracy_value','')")
    context = emit(
        "( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT"
        f"(({uncertainty})) GLOBAL_UNIT_ASSIGNED_CONTEXT(({unit},{angle},{solid_angle}))"
        " REPRESENTATION_CONTEXT('','faceted audit export') )"
    )

    point_ids: list[str] = []
    vertex_ids: list[str] = []
    for x, y, z in mesh.vertices:
        point_ids.append(emit(f"CARTESIAN_POINT('',({number(x)},{number(y)},{number(z)}))"))
        vertex_ids.append(emit(f"VERTEX_POINT('',{point_ids[-1]})"))

    face_ids: list[str] = []
    for a, b, c in mesh.faces:
        loop = emit(f"POLY_LOOP('',({vertex_ids[a]},{vertex_ids[b]},{vertex_ids[c]}))")
        bound = emit(f"FACE_OUTER_BOUND('',{loop},.T.)")
        face_ids.append(emit(f"FACE('',({bound}))"))

    shell = emit(f"CLOSED_SHELL('',({','.join(face_ids)}))")
    solid = emit(f"FACETED_BREP('{name}',{shell})")
    emit(
        "ADVANCED_BREP_SHAPE_REPRESENTATION"
        f"('{name}',({solid},{context}))"
    )

    header = [
        "ISO-10303-21;",
        "HEADER;",
        "FILE_DESCRIPTION(('faceted triangle mesh for geometry audit'),'2;1');",
        f"FILE_NAME('{Path(path).name}','',(''),(''),'tensorcad.exchange','tensorcad','');",
        "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));",
        "ENDSEC;",
        "DATA;",
    ]
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text("\n".join(header + lines + ["ENDSEC;", "END-ISO-10303-21;", ""]), encoding="utf-8")
    return target


_STATEMENT = re.compile(r"#(\d+)\s*=\s*([A-Z0-9_]+)?\s*\((.*)\)\s*$", re.DOTALL)
_REF_LIST = re.compile(r"\(\s*((?:#\d+\s*,?\s*)+)\)")


def _statements(text: str) -> dict[str, tuple[str, str]]:
    """Map entity id -> (entity type, raw argument string)."""
    data_start = text.find("DATA;")
    data_end = text.rfind("ENDSEC;")
    body = text[data_start + len("DATA;") : data_end] if data_start >= 0 else text
    result: dict[str, tuple[str, str]] = {}
    for chunk in body.split(";"):
        statement = chunk.strip()
        if not statement:
            continue
        if not statement.startswith("#"):
            continue
        match = _STATEMENT.match(statement)
        if match is None:
            raise ExchangeError(f"cannot parse STEP statement {statement[:80]!r}")
        # A complex entity (``#1 = ( A() B() )``) has no single keyword; keeping it
        # under the empty type means the readers below simply ignore it.
        result[match.group(1)] = ((match.group(2) or ""), match.group(3))
    return result


def _refs(arguments: str) -> list[str]:
    """All entity references in an argument string, in file order.

    Both shapes occur: an entity reference list ``('',(#3,#5))`` and a single
    reference ``('',#1)``. Scanning for ``#n`` finds both without a grammar.
    """
    return re.findall(r"#(\d+)", arguments)


def _floats(arguments: str) -> tuple[float, float, float]:
    start = arguments.find("(")
    inner = arguments[start + 1 : arguments.rfind(")")]
    values = [part.strip() for part in inner.split(",") if part.strip()]
    if len(values) != 3:
        raise ExchangeError(f"expected a 3-tuple of numbers, got {arguments!r}")
    try:
        return (float(values[0]), float(values[1]), float(values[2]))
    except ValueError as error:  # pragma: no cover - defensive
        raise ExchangeError(f"non-numeric coordinate in {arguments!r}") from error


def read_step_facets(path: str | Path, *, length_scale: float | None = None) -> PatchGeometry:
    """Read a faceted STEP file (``FACETED_BREP`` over ``POLY_LOOP``) into a mesh.

    ``length_scale`` defaults to the millimetre convention
    (:data:`DEFAULT_LENGTH_SCALE`), i.e. file coordinates are divided by 1000 to
    give metres.  Pass ``1.0`` for a file already written in metres.

    Analytic B-rep files are refused with :class:`ExchangeError` listing the
    entities that block the import rather than being tessellated by guesswork.
    """
    source = Path(path)
    text = source.read_text(encoding="utf-8", errors="replace")
    if "ISO-10303-21" not in text:
        raise ExchangeError(f"{source} is not an ISO-10303-21 (STEP) file")
    statements = _statements(text)
    types = {entity_type for entity_type, _ in statements.values()}

    patches = [name for name in types if name in ("FACE", "ADVANCED_FACE")]
    unsupported = sorted(types.intersection(_BREP_ENTITIES) - {"FACE"})
    if not patches:
        raise ExchangeError(
            f"{source} contains no face entities; found {len(statements)} entities of "
            f"{len(types)} types"
        )
    if unsupported:
        raise ExchangeError(
            "STEP file is an analytic B-rep, not a faceted solid; refusing to guess a "
            f"tessellation. Blocking entities: {', '.join(unsupported)}. "
            "Install the optional OCP bridge (pip install 'tensorcad[ocp]') to tessellate it."
        )
    scale = DEFAULT_LENGTH_SCALE if length_scale is None else float(length_scale)
    if scale == 0.0:
        raise ExchangeError("length_scale must be non-zero")

    vertex_by_point: dict[str, int] = {}
    vertices: list[tuple[float, float, float]] = []
    for entity_id, (entity_type, arguments) in statements.items():
        if entity_type != "VERTEX_POINT":
            continue
        refs = _refs(arguments)
        if not refs or refs[0] not in statements:
            raise ExchangeError(f"VERTEX_POINT {entity_id} has no CARTESIAN_POINT")
        point_type, point_arguments = statements[refs[0]]
        if point_type != "CARTESIAN_POINT":
            raise ExchangeError(f"VERTEX_POINT {entity_id} points at {point_type}")
        x, y, z = _floats(point_arguments)
        vertex_by_point[entity_id] = len(vertices)
        vertices.append((x / scale, y / scale, z / scale))
    if not vertices:
        raise ExchangeError(f"{source} has no VERTEX_POINT entities")

    loop_faces: dict[str, tuple[int, ...]] = {}
    for entity_id, (entity_type, arguments) in statements.items():
        if entity_type != "POLY_LOOP":
            continue
        refs = _refs(arguments)
        try:
            loop_faces[entity_id] = tuple(vertex_by_point[ref] for ref in refs)
        except KeyError as error:  # pragma: no cover - malformed file
            raise ExchangeError(f"POLY_LOOP {entity_id} references unknown vertex {error}") from error

    faces: list[tuple[int, int, int]] = []
    patch_ids: list[str] = []
    bounds: dict[str, list[str]] = {}
    for entity_id, (entity_type, arguments) in statements.items():
        if entity_type.endswith("BOUND"):
            bounds[entity_id] = _refs(arguments)
    skipped = 0
    for entity_id, (entity_type, arguments) in statements.items():
        if entity_type not in ("FACE", "ADVANCED_FACE"):
            continue
        # FACE -> FACE_OUTER_BOUND -> POLY_LOOP; tolerate a FACE that points at a
        # loop directly (some writers do).
        face_bounds = [ref for ref in _refs(arguments) if ref in bounds]
        loops = [loop for bound in face_bounds for loop in bounds[bound] if loop in loop_faces]
        if not loops:
            loops = [ref for ref in _refs(arguments) if ref in loop_faces]
        polygon: list[int] = []
        for loop in loops:
            polygon.extend(loop_faces.get(loop, ()))
        if len(polygon) < 3:
            skipped += 1
            continue
        for index in range(1, len(polygon) - 1):
            faces.append((polygon[0], polygon[index], polygon[index + 1]))
            patch_ids.append(entity_id)
    if not faces:
        raise ExchangeError(f"{source} yielded no triangles from {len(patches)} face entities")

    mesh = Mesh(vertices=tuple(vertices), faces=tuple(faces))
    return PatchGeometry(
        mesh=mesh,
        patch_ids=tuple(patch_ids),
        source=str(source),
        metadata={
            "schema": STEP_FACETS_SCHEMA,
            "length_scale": scale,
            "entities": len(statements),
            "entity_types": sorted(types),
            "vertex_points": len(vertices),
            "faces": len(faces),
            "skipped_faces": skipped,
            "validation": "internal round trip only; no third-party CAD kernel available",
        },
    )


def round_trip_facets(mesh: Mesh, *, directory: str | Path | None = None) -> dict[str, float]:
    """Write ``mesh`` to faceted STEP and read it back; report the geometry drift."""
    import tempfile

    with tempfile.TemporaryDirectory(dir=directory) as scratch:
        target = Path(scratch) / "round_trip.step"
        write_step_facets(target, mesh)
        back = read_step_facets(target)
    if back.mesh.faces != mesh.faces:
        raise ExchangeError("faceted STEP round trip changed the face table")
    deviation = 0.0
    for original, restored in zip(mesh.vertices, back.mesh.vertices):
        for a, b in zip(original, restored):
            deviation = max(deviation, abs(a - b))
    return {
        "vertex_deviation_m": deviation,
        "vertex_count": float(len(mesh.vertices)),
        "triangle_count": float(len(mesh.faces)),
        "length_scale": float(back.metadata["length_scale"]),
    }


# --------------------------------------------------------------------------- #
# dispatch
# --------------------------------------------------------------------------- #
def load_exchange(path: str | Path, *, length_scale: float | None = None) -> PatchGeometry:
    """Load a geometry exchange file, dispatching on the file name.

    ``.json``/``.patchset`` -> :func:`read_patchset`; ``.step``/``.stp`` ->
    :func:`read_step_facets`; ``.igs``/``.iges`` -> :class:`ExchangeUnavailable`
    (IGES needs OpenCASCADE, which is not installed here).
    """
    source = Path(path)
    name = source.name.lower()
    if name.endswith(".patchset.json") or name.endswith(".json"):
        return read_patchset(source)
    if name.endswith(".step") or name.endswith(".stp"):
        return read_step_facets(source, length_scale=length_scale)
    if name.endswith(".igs") or name.endswith(".iges"):
        raise ExchangeUnavailable(
            "IGES import needs OpenCASCADE (pip install 'tensorcad[ocp]'); IGES also has no "
            "closed-shell guarantee, so a gap audit would need the B-rep, not a tessellation"
        )
    raise ExchangeError(
        f"unrecognised exchange extension {source.suffix!r}; supported: .step, .stp, .json"
    )


def _finite(values: Iterable[float]) -> bool:
    return all(math.isfinite(value) for value in values)


def supported_routes() -> dict[str, str]:
    """Evidence level per exchange route, for audit reports."""
    return {
        "step_facets_read": "REAL",
        "step_facets_write": "REAL",
        "patchset_read": "REAL",
        "step_brep_read": "NOT_SUPPORTED (needs OCP bridge)",
        "iges_read": "NOT_SUPPORTED (needs OCP bridge)",
        "step_file_external_validation": "WITHHELD (no third-party CAD kernel installed)",
    }