"""Independent verification of the geometry-exchange seam (:mod:`tensorcad.exchange`).

The seam exists to take what a shipyard hands over -- a faceted STEP or a neutral
patch table -- and give the quality audit a mesh it can judge. Two failure modes are
worse than a refusal, because they are silent:

* **a wrong scale** (a millimetre label on metre numbers, or the reverse), and
* **a wrong triangulation** (a fan across a hole, or across a concave loop).

Both were present before the exchange verification mission: the writer labelled every
file ``SI_UNIT(.MILLI.,.METRE.)`` whatever ``length_scale`` it used, the reader ignored
the label and always divided by 1000, and a loop of more than three vertices was fanned
from its first vertex without a convexity test. The tests below assert the *fixed*
contract and re-derive the failure magnitudes, so a regression is red rather than
quietly wrong.

Structure of this file:

1. an independent structural parse of the written STEP (framing, entity graph,
   reference closure, entity counts, unit token) -- written against the file text, not
   against the writer's own bookkeeping;
2. scale honesty: writer declaration vs numbers, reader from declaration, override;
3. duplicate-vertex welding (per-corner ``VERTEX_POINT`` writers);
4. refusals: analytic B-rep, holes, non-convex loops, truncation, dangling refs,
   zero/negative scales, empty mesh, unrepresentable scale, IGES;
5. round-trip fidelity at three scales, patchset round trip, dispatch.
"""
from __future__ import annotations

import math
import re

import pytest

from tensorcad.exchange import (
    DEFAULT_LENGTH_SCALE,
    ExchangeError,
    ExchangeUnavailable,
    PatchGeometry,
    load_exchange,
    read_patchset,
    read_step_facets,
    round_trip_facets,
    supported_routes,
    write_patchset,
    write_step_facets,
)
from tensorcad.geometry import Mesh
from tensorcad.mesh import is_watertight

# --------------------------------------------------------------------------- #
# fixtures: a closed unit cube (12 triangles) and helpers to hand-build STEP files
# --------------------------------------------------------------------------- #
CUBE_VERTICES = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
)
CUBE_FACES = (
    (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
    (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
    (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
)


def cube() -> Mesh:
    return Mesh(vertices=CUBE_VERTICES, faces=CUBE_FACES)


def triangle_area(triangle, vertices) -> float:
    a, b, c = (vertices[index] for index in triangle)
    u = [b[i] - a[i] for i in range(3)]
    w = [c[i] - a[i] for i in range(3)]
    cross = (
        u[1] * w[2] - u[2] * w[1],
        u[2] * w[0] - u[0] * w[2],
        u[0] * w[1] - u[1] * w[0],
    )
    return 0.5 * math.sqrt(sum(value * value for value in cross))


def polygon_area(points) -> float:
    """Shoelace area of a planar polygon in the xy-plane (absolute value)."""
    total = 0.0
    count = len(points)
    for index in range(count):
        x1, y1 = points[index][0], points[index][1]
        x2, y2 = points[(index + 1) % count][0], points[(index + 1) % count][1]
        total += x1 * y2 - x2 * y1
    return 0.5 * abs(total)


def fan_area(points, start: int = 0) -> float:
    """Area the *unguarded* fan triangulation would report, for the failure record."""
    ordered = list(points[start:]) + list(points[:start])
    total = 0.0
    for index in range(1, len(ordered) - 1):
        total += triangle_area((0, index, index + 1), ordered)
    return total


HEADER = (
    "ISO-10303-21;",
    "HEADER;",
    "FILE_DESCRIPTION(('hand-built fixture'),'2;1');",
    "FILE_NAME('fixture.step','',(''),(''),'tests','tensorcad','');",
    "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));",
    "ENDSEC;",
    "DATA;",
)


class StepWriter:
    """Minimal hand-written STEP emitter, so fixtures are not produced by the code under test."""

    def __init__(self, length_unit: str = "SI_UNIT($,.METRE.)") -> None:
        self.statements: list[str] = []
        self.counter = 0
        self.emit(f"( LENGTH_UNIT() NAMED_UNIT(*) {length_unit} )")
        self.emit("( GEOMETRIC_REPRESENTATION_CONTEXT(3) REPRESENTATION_CONTEXT('','t') )")

    def emit(self, body: str) -> str:
        self.counter += 1
        entity = f"#{self.counter}"
        self.statements.append(f"{entity} = {body};")
        return entity

    def point(self, x: float, y: float, z: float) -> str:
        return self.emit(f"CARTESIAN_POINT('',({x:.17g},{y:.17g},{z:.17g}))")

    def vertex(self, x: float, y: float, z: float) -> str:
        return self.emit(f"VERTEX_POINT('',{self.point(x, y, z)})")

    def poly_loop(self, vertices) -> str:
        return self.emit(f"POLY_LOOP('',({','.join(vertices)}))")

    def face(self, loops, *, outer_first: bool = True) -> str:
        bounds = []
        for index, loop in enumerate(loops):
            kind = "FACE_OUTER_BOUND" if (index == 0 and outer_first) else "FACE_BOUND"
            bounds.append(self.emit(f"{kind}('',{loop},.T.)"))
        return self.emit(f"FACE('',({','.join(bounds)}))")

    def solid(self, faces, name: str = "fixture", context: str | None = None) -> str:
        shell = self.emit(f"CLOSED_SHELL('',({','.join(faces)}))")
        brep = self.emit(f"FACETED_BREP('{name}',{shell})")
        context = context or "#2"
        self.emit(f"ADVANCED_BREP_SHAPE_REPRESENTATION('{name}',({brep},{context}))")
        return brep

    def text(self, *, terminate: bool = True) -> str:
        tail = ["ENDSEC;", "END-ISO-10303-21;", ""] if terminate else []
        return "\n".join(list(HEADER) + self.statements + tail)


def write(path, text: str):
    path.write_text(text, encoding="utf-8")
    return path


# --------------------------------------------------------------------------- #
# 1. independent structural parse of the written file
# --------------------------------------------------------------------------- #
ENTITY = re.compile(r"^#(\d+)\s*=\s*(.+)$", re.DOTALL)


def parse_entities(text: str) -> dict[int, str]:
    """Entity id -> body, read straight off the file with a separate grammar."""
    start = text.index("DATA;") + len("DATA;")
    stop = text.rindex("ENDSEC;")
    entities: dict[int, str] = {}
    for chunk in text[start:stop].split(";"):
        statement = chunk.strip()
        if not statement:
            continue
        match = ENTITY.match(statement)
        assert match is not None, f"unparsable statement {statement[:60]!r}"
        entities[int(match.group(1))] = match.group(2).strip()
    return entities


def keyword(body: str) -> str:
    match = re.match(r"\s*(?:\(\s*)?([A-Z_]+)", body)
    return match.group(1) if match else ""


def test_the_written_file_is_structurally_closed(tmp_path):
    """Every reference resolves, ids are contiguous, and the framing is complete."""
    path = write_step_facets(tmp_path / "cube.step", cube())
    text = path.read_text()
    assert text.startswith("ISO-10303-21;")
    assert text.rstrip().endswith("END-ISO-10303-21;")
    assert text.count("HEADER;") == 1 and text.count("ENDSEC;") == 2

    entities = parse_entities(text)
    assert sorted(entities) == list(range(1, len(entities) + 1)), "entity ids must be 1..n"
    for entity_id, body in entities.items():
        for ref in re.findall(r"#(\d+)", body):
            assert int(ref) in entities, f"#{entity_id} references missing #{ref}"


def test_the_entity_counts_are_the_mesh_counts(tmp_path):
    """8 points, 8 vertices, one loop+bound+face per triangle, one shell and one brep."""
    mesh = cube()
    path = write_step_facets(tmp_path / "cube.step", mesh)
    bodies = parse_entities(path.read_text())
    counts: dict[str, int] = {}
    for body in bodies.values():
        counts[keyword(body)] = counts.get(keyword(body), 0) + 1
    assert counts["CARTESIAN_POINT"] == len(mesh.vertices) == 8
    assert counts["VERTEX_POINT"] == 8
    assert counts["POLY_LOOP"] == len(mesh.faces) == 12
    assert counts["FACE_OUTER_BOUND"] == 12
    assert counts["FACE"] == 12
    assert counts["CLOSED_SHELL"] == 1
    assert counts["FACETED_BREP"] == 1
    assert counts["ADVANCED_BREP_SHAPE_REPRESENTATION"] == 1


def test_every_loop_is_a_triangle_of_distinct_vertices(tmp_path):
    """The writer's claim: one planar face per triangle, no n-gons, no repeated corner."""
    path = write_step_facets(tmp_path / "cube.step", cube())
    bodies = parse_entities(path.read_text())
    loops = [body for body in bodies.values() if keyword(body) == "POLY_LOOP"]
    for body in loops:
        refs = re.findall(r"#(\d+)", body)
        assert len(refs) == 3
        assert len(set(refs)) == 3


def test_the_closed_shell_references_every_face_exactly_once(tmp_path):
    path = write_step_facets(tmp_path / "cube.step", cube())
    bodies = parse_entities(path.read_text())
    faces = {entity for entity, body in bodies.items() if keyword(body) == "FACE"}
    shell_id = next(entity for entity, body in bodies.items() if keyword(body) == "CLOSED_SHELL")
    referenced = [int(ref) for ref in re.findall(r"#(\d+)", bodies[shell_id])]
    assert sorted(referenced) == sorted(faces)
    assert len(referenced) == 12


def test_the_declared_length_unit_matches_the_numbers(tmp_path):
    """The label is only honest if the coordinates are in that unit. Check both scales."""
    for scale, token in ((1000.0, ".MILLI."), (1.0, "$"), (1e6, ".MICRO.")):
        path = write_step_facets(tmp_path / f"cube_{scale}.step", cube(), length_scale=scale)
        text = path.read_text()
        assert f"SI_UNIT({token},.METRE.)" in text
        bodies = parse_entities(text)
        coords = [
            float(value)
            for body in bodies.values()
            if keyword(body) == "CARTESIAN_POINT"
            for value in re.findall(r"[-+0-9.eE]+", body.split("(", 1)[1].rsplit(")", 1)[0])
        ]
        assert max(abs(value) for value in coords) == pytest.approx(scale)


# --------------------------------------------------------------------------- #
# 2. scale honesty, writer -> file -> reader
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("scale", [1.0, 1000.0, 1e6])
def test_a_round_trip_at_any_declared_scale_recovers_the_coordinates(tmp_path, scale):
    """The reader is given no scale, so this measures the declaration path end to end."""
    report = round_trip_facets(cube(), directory=tmp_path, length_scale=scale)
    assert report["used_declared_unit"] == 1.0
    assert report["length_scale"] == pytest.approx(scale)
    assert report["triangle_count"] == 12.0
    assert report["vertex_deviation_m"] <= 1e-15


def test_the_old_failure_was_1000x(tmp_path):
    """A metre file read with the old hard-coded millimetre divisor was 0.001x too small.

    The number is re-derived here from the current reader with the scale forced to the
    old assumption, so the record is reproducible rather than quoted from history.
    """
    path = write_step_facets(tmp_path / "metre.step", cube(), length_scale=1.0)
    honest = read_step_facets(path)
    assumed_mm = read_step_facets(path, length_scale=DEFAULT_LENGTH_SCALE)
    assert max(abs(value) for vertex in honest.mesh.vertices for value in vertex) == pytest.approx(1.0)
    assert max(abs(value) for vertex in assumed_mm.mesh.vertices for value in vertex) == pytest.approx(1e-3)
    assert honest.metadata["length_scale_source"].startswith("declared")
    assert assumed_mm.metadata["length_scale_source"] == "caller override"


def test_a_metre_declaring_third_party_file_is_read_in_metres(tmp_path):
    """Hand-built (not produced by this writer): SI_UNIT($,.METRE.), numbers in metres."""
    writer = StepWriter("SI_UNIT($,.METRE.)")
    faces = []
    for triangle in CUBE_FACES:
        vertices = [writer.vertex(*CUBE_VERTICES[index]) for index in triangle]
        faces.append(writer.face([writer.poly_loop(vertices)]))
    writer.solid(faces)
    path = write(tmp_path / "third_party_metre.step", writer.text())

    geometry = read_step_facets(path)
    assert geometry.metadata["length_scale"] == pytest.approx(1.0)
    assert geometry.metadata["length_scale_source"] == "declared SI_UNIT($,.METRE.)"
    assert max(abs(value) for vertex in geometry.mesh.vertices for value in vertex) == pytest.approx(1.0)
    assert len(geometry.mesh.faces) == 12

    override = read_step_facets(path, length_scale=1.0)
    assert override.metadata["length_scale_source"] == "caller override"


def test_a_file_without_a_length_unit_falls_back_and_says_so(tmp_path):
    writer = StepWriter("SI_UNIT(.MILLI.,.METRE.)")
    writer.statements = [s for s in writer.statements if "LENGTH_UNIT" not in s]
    vertices = [writer.vertex(x, y, z) for (x, y, z) in CUBE_VERTICES]
    faces = []
    for triangle in CUBE_FACES:
        faces.append(writer.face([writer.poly_loop([vertices[index] for index in triangle])]))
    writer.solid(faces)
    path = write(tmp_path / "no_unit.step", writer.text())
    geometry = read_step_facets(path)
    assert geometry.metadata["length_scale"] == pytest.approx(DEFAULT_LENGTH_SCALE)
    assert geometry.metadata["length_scale_source"].startswith("assumed millimetre convention")


def test_an_unrepresentable_scale_is_refused_rather_than_mislabelled(tmp_path):
    """25.4 file units per metre has no SI prefix; writing it would lie about the unit."""
    with pytest.raises(ExchangeError, match="no SI length prefix"):
        write_step_facets(tmp_path / "inch.step", cube(), length_scale=25.4)


def test_a_quote_in_the_solid_name_is_refused(tmp_path):
    with pytest.raises(ExchangeError, match="quote"):
        write_step_facets(tmp_path / "q.step", cube(), name="o'brien")


# --------------------------------------------------------------------------- #
# 3. duplicate-vertex welding
# --------------------------------------------------------------------------- #
def per_corner_step() -> str:
    """A closed cube written the way a per-corner writer does it: 36 VERTEX_POINTs."""
    writer = StepWriter("SI_UNIT(.MILLI.,.METRE.)")
    faces = []
    for triangle in CUBE_FACES:
        vertices = [writer.vertex(*(1000.0 * c for c in CUBE_VERTICES[index])) for index in triangle]
        faces.append(writer.face([writer.poly_loop(vertices)]))
    writer.solid(faces)
    return writer.text()


def test_per_corner_vertices_weld_into_a_watertight_solid(tmp_path):
    """Without welding, a closed solid reads back open and the audit blames the model."""
    path = write(tmp_path / "dup.step", per_corner_step())
    welded = read_step_facets(path)
    assert len(welded.mesh.vertices) == 8
    assert welded.metadata["welded_vertex_rows"] == 28
    assert len(welded.mesh.faces) == 12
    assert is_watertight(welded.mesh) is True

    raw = read_step_facets(path, weld_duplicates=False)
    assert len(raw.mesh.vertices) == 36
    assert raw.metadata["welded_vertex_rows"] == 0
    assert is_watertight(raw.mesh) is False


def test_welding_does_not_move_any_coordinate(tmp_path):
    path = write(tmp_path / "dup.step", per_corner_step())
    geometry = read_step_facets(path)
    assert sorted(geometry.mesh.vertices) == sorted(CUBE_VERTICES)
    # welding renumbers vertices in order of first appearance, so the face table is
    # index-remapped; the triangles (as coordinates) must be identical as a set.
    def as_coordinate_triangles(mesh):
        return sorted(
            tuple(sorted(mesh.vertices[index] for index in triangle))
            for triangle in mesh.faces
        )

    assert as_coordinate_triangles(geometry.mesh) == as_coordinate_triangles(cube())


def test_patch_ids_line_up_with_the_faces_they_came_from(tmp_path):
    path = write_step_facets(tmp_path / "cube.step", cube())
    geometry = read_step_facets(path)
    assert len(geometry.patch_ids) == len(geometry.mesh.faces) == 12
    bodies = parse_entities(path.read_text())
    face_ids = {str(entity) for entity, body in bodies.items() if keyword(body) == "FACE"}
    assert set(geometry.patch_ids) <= face_ids
    assert len(set(geometry.patch_ids)) == 12  # one triangle per emitted FACE


# --------------------------------------------------------------------------- #
# 4. refusals: never a silent wrong tessellation
# --------------------------------------------------------------------------- #
def test_an_analytic_brep_is_refused_by_name(tmp_path):
    writer = StepWriter()
    writer.solid([])
    text = writer.text().replace(
        "DATA;",
        "DATA;\n"
        + f"#{writer.counter + 1} = ADVANCED_FACE('',(#2),#2,.T.);\n"
        + f"#{writer.counter + 2} = CYLINDRICAL_SURFACE('',#2,5.);",
    )
    path = write(tmp_path / "brep.step", text)
    with pytest.raises(ExchangeError, match="analytic B-rep"):
        read_step_facets(path)


def test_a_face_with_a_hole_is_refused_not_fanned(tmp_path):
    """A fan over (outer + inner) invents the hole's area: 6.25 vs the true 3.0."""
    writer = StepWriter("SI_UNIT($,.METRE.)")
    outer_points = ((0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0))
    inner_points = ((0.5, 0.5, 0), (1.5, 0.5, 0), (1.5, 1.5, 0), (0.5, 1.5, 0))
    outer = writer.poly_loop([writer.vertex(*point) for point in outer_points])
    inner = writer.poly_loop([writer.vertex(*point) for point in inner_points])
    writer.solid([writer.face([outer, inner])])
    path = write(tmp_path / "hole.step", writer.text())

    # what the unguarded concatenate-and-fan used to report
    concatenated = list(outer_points) + list(inner_points)
    assert round(fan_area(concatenated), 6) == 6.25
    assert round(polygon_area(outer_points) - polygon_area(inner_points), 6) == 3.0

    with pytest.raises(ExchangeError, match="hole"):
        read_step_facets(path)


def test_a_non_convex_loop_is_refused_not_fanned(tmp_path):
    """U-shaped loop: the fan reports 11.0 where the polygon is 7.0 (+57 %)."""
    writer = StepWriter("SI_UNIT($,.METRE.)")
    points = ((0, 0, 0), (3, 0, 0), (3, 3, 0), (2, 3, 0), (2, 1, 0), (1, 1, 0), (1, 3, 0), (0, 3, 0))
    loop = writer.poly_loop([writer.vertex(*point) for point in points])
    writer.solid([writer.face([loop])])
    path = write(tmp_path / "u.step", writer.text())

    assert round(fan_area(list(points)), 6) == 11.0
    assert round(polygon_area(list(points)), 6) == 7.0
    with pytest.raises(ExchangeError, match="not convex"):
        read_step_facets(path)


def test_a_convex_quad_is_still_accepted_and_exact(tmp_path):
    """A fan over a convex quad is exact, so it must keep working (no over-refusal)."""
    writer = StepWriter("SI_UNIT($,.METRE.)")
    points = ((0, 0, 0), (2, 0, 0), (2, 1, 0), (0, 1, 0))
    loop = writer.poly_loop([writer.vertex(*point) for point in points])
    writer.solid([writer.face([loop])])
    path = write(tmp_path / "quad.step", writer.text())
    geometry = read_step_facets(path)
    fan = sum(triangle_area(triangle, geometry.mesh.vertices) for triangle in geometry.mesh.faces)
    assert len(geometry.mesh.faces) == 2
    assert fan == pytest.approx(2.0, rel=1e-15)


def test_a_truncated_file_is_refused(tmp_path):
    """Dropping the closing statements must not read as a smaller, valid body."""
    text = write_step_facets(tmp_path / "cube.step", cube()).read_text()
    truncated = "\n".join(line for line in text.splitlines() if not line.startswith(("ENDSEC", "END-ISO")))
    path = write(tmp_path / "truncated.step", truncated)
    with pytest.raises(ExchangeError, match="truncated"):
        read_step_facets(path)


def test_a_dangling_reference_is_refused(tmp_path):
    writer = StepWriter()
    vertices = [writer.vertex(x, y, z) for (x, y, z) in CUBE_VERTICES]
    faces = [
        writer.face([writer.poly_loop([vertices[index] for index in triangle])])
        for triangle in CUBE_FACES
    ]
    writer.solid(faces)
    bodies = writer.text()
    # point the last FACE_OUTER_BOUND at an entity that does not exist
    broken = re.sub(r"(FACE_OUTER_BOUND\('',)#\d+(,\.T\.\));$", r"\g<1>#9999\g<2>", bodies, count=1, flags=re.MULTILINE)
    assert broken != bodies
    path = write(tmp_path / "dangling.step", broken)
    with pytest.raises(ExchangeError, match="does not exist"):
        read_step_facets(path)


def test_negative_controls(tmp_path):
    with pytest.raises(ExchangeError, match="no faces"):
        write_step_facets(tmp_path / "empty.step", Mesh(vertices=(), faces=()))
    with pytest.raises(ExchangeError, match="positive"):
        write_step_facets(tmp_path / "zero.step", cube(), length_scale=0.0)
    with pytest.raises(ExchangeError, match="positive"):
        write_step_facets(tmp_path / "neg.step", cube(), length_scale=-1.0)
    good = write_step_facets(tmp_path / "cube.step", cube())
    with pytest.raises(ExchangeError, match="non-zero"):
        read_step_facets(good, length_scale=0.0)
    path = write(tmp_path / "plain.txt", "not a step file\n")
    with pytest.raises(ExchangeError, match="not an ISO-10303-21"):
        read_step_facets(path)


# --------------------------------------------------------------------------- #
# 5. patchset round trip, dispatch, evidence table
# --------------------------------------------------------------------------- #
def test_patchset_round_trip_preserves_mesh_and_patches(tmp_path):
    geometry = PatchGeometry(
        mesh=cube(),
        patch_ids=("hull",) * 12,
        source="fixture",
        metadata={"note": "cube"},
    )
    path = write_patchset(tmp_path / "cube.patchset.json", geometry)
    back = read_patchset(path)
    assert back.mesh.vertices == geometry.mesh.vertices
    assert back.mesh.faces == geometry.mesh.faces
    assert back.patch_ids == geometry.patch_ids
    assert back.metadata["note"] == "cube"


def test_load_exchange_dispatches_on_the_extension(tmp_path):
    step = write_step_facets(tmp_path / "cube.step", cube())
    assert len(load_exchange(step).mesh.faces) == 12

    patchset = write_patchset(
        tmp_path / "cube.json", PatchGeometry(mesh=cube(), patch_ids=("hull",) * 12)
    )
    assert len(load_exchange(patchset).mesh.faces) == 12

    iges = write(tmp_path / "hull.iges", "start\n")
    with pytest.raises(ExchangeUnavailable, match="OpenCASCADE"):
        load_exchange(iges)
    with pytest.raises(ExchangeError, match="unrecognised exchange extension"):
        load_exchange(write(tmp_path / "hull.dwg", "x"))


def test_the_evidence_table_keeps_its_honest_levels():
    """Analytic B-rep reads are REAL for the supported table; external validation is WITHHELD."""
    routes = supported_routes()
    assert routes["step_facets_write"] == "REAL"
    assert routes["step_facets_read"] == "REAL"
    assert routes["patchset_read"] == "REAL"
    # tensorcad.step_brep reads analytic trims; anything outside its table is refused
    # by name, which is why the level is REAL *for that table* and not a blanket claim.
    assert routes["step_brep_read"].startswith("REAL for")
    assert "refused by name" in routes["step_brep_read"]
    assert routes["iges_read"].startswith("NOT_SUPPORTED")
    assert routes["step_file_external_validation"].startswith("WITHHELD")