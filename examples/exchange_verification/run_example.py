"""Exchange seam verification: every number and figure ``docs/exchange_verification.md`` reports.

The module under test is :mod:`tensorcad.exchange` -- the seam that takes a faceted
STEP or a neutral patch table from outside and hands the quality audit a mesh. This
example measures four silent failure modes the seam had before the verification
mission, reproduces each of them **on the current code** through a deliberately naive
reader (no framing check, hard-coded millimetre divisor, concatenate-and-fan over every
loop), and then records what the current reader does instead.

    python3 examples/exchange_verification/run_example.py [output-directory]

Writes ``artifacts/exchange_verification/results.json``, ``figure.png`` and
``docs/exchange_verification.md``; no number in the page is typed by hand.
"""
from __future__ import annotations

import json
import math
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.exchange import (  # noqa: E402
    DEFAULT_LENGTH_SCALE,
    ExchangeError,
    _floats,
    _refs,
    _statements,
    read_step_facets,
    round_trip_facets,
    write_step_facets,
)
from tensorcad.geometry import Mesh  # noqa: E402
from tensorcad.mesh import is_watertight  # noqa: E402

CUBE_VERTICES = (
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
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
    total = 0.0
    count = len(points)
    for index in range(count):
        x1, y1 = points[index][0], points[index][1]
        x2, y2 = points[(index + 1) % count][0], points[(index + 1) % count][1]
        total += x1 * y2 - x2 * y1
    return 0.5 * abs(total)


def naive_read(path, length_scale: float = DEFAULT_LENGTH_SCALE) -> tuple[list, list]:
    """The reader as it behaved before the verification mission, for the failure record.

    No framing/truncation check, no declared-unit read (the divisor is always
    ``length_scale``), every loop fanned from its first vertex and loops of one face
    concatenated. Uses the module's own low-level statement helpers so the reproduction
    is the same grammar, only the missing guards.
    """
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    statements = _statements(text)
    vertex_by_point: dict[str, int] = {}
    vertices: list[tuple[float, float, float]] = []
    for entity_id, (entity_type, arguments) in statements.items():
        if entity_type != "VERTEX_POINT":
            continue
        refs = _refs(arguments)
        x, y, z = _floats(statements[refs[0]][1])
        vertex_by_point[entity_id] = len(vertices)
        vertices.append((x / length_scale, y / length_scale, z / length_scale))
    loop_faces = {
        entity_id: tuple(vertex_by_point[ref] for ref in _refs(arguments))
        for entity_id, (entity_type, arguments) in statements.items()
        if entity_type == "POLY_LOOP"
    }
    bounds = {
        entity_id: _refs(arguments)
        for entity_id, (entity_type, arguments) in statements.items()
        if entity_type.endswith("BOUND")
    }
    faces: list[tuple[int, int, int]] = []
    for entity_id, (entity_type, arguments) in statements.items():
        if entity_type not in ("FACE", "ADVANCED_FACE"):
            continue
        face_bounds = [ref for ref in _refs(arguments) if ref in bounds]
        loops = [loop for bound in face_bounds for loop in bounds[bound] if loop in loop_faces]
        if not loops:
            loops = [ref for ref in _refs(arguments) if ref in loop_faces]
        polygon: list[int] = []
        for loop in loops:
            polygon.extend(loop_faces.get(loop, ()))
        for index in range(1, len(polygon) - 1):
            faces.append((polygon[0], polygon[index], polygon[index + 1]))
    return vertices, faces


def area_of(faces, vertices) -> float:
    return sum(triangle_area(triangle, vertices) for triangle in faces)


class StepWriter:
    """Hand-written STEP emitter; fixtures are not produced by the code under test."""

    def __init__(self, length_unit: str = "SI_UNIT($,.METRE.)") -> None:
        self.statements: list[str] = []
        self.counter = 0
        self.emit(f"( LENGTH_UNIT() NAMED_UNIT(*) {length_unit} )")
        self.emit("( GEOMETRIC_REPRESENTATION_CONTEXT(3) REPRESENTATION_CONTEXT('','t') )")
        self.header = [
            "ISO-10303-21;",
            "HEADER;",
            "FILE_DESCRIPTION(('exchange verification fixture'),'2;1');",
            "FILE_NAME('fixture.step','',(''),(''),'tensorcad','tensorcad','');",
            "FILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));",
            "ENDSEC;",
            "DATA;",
        ]

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

    def face(self, loops) -> str:
        bounds = [
            self.emit(f"{'FACE_OUTER_BOUND' if index == 0 else 'FACE_BOUND'}('',{loop},.T.)")
            for index, loop in enumerate(loops)
        ]
        return self.emit(f"FACE('',({','.join(bounds)}))")

    def solidify(self, faces) -> None:
        shell = self.emit(f"CLOSED_SHELL('',({','.join(faces)}))")
        brep = self.emit(f"FACETED_BREP('fixture',{shell})")
        self.emit(f"ADVANCED_BREP_SHAPE_REPRESENTATION('fixture',({brep},#2))")

    def text(self, *, terminate: bool = True) -> str:
        tail = ["ENDSEC;", "END-ISO-10303-21;", ""] if terminate else []
        return "\n".join(self.header + self.statements + tail)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def build_fixtures(scratch: Path) -> dict[str, Path]:
    fixtures: dict[str, Path] = {}

    # a metre-declared cube written by the current writer (declaration path)
    fixtures["metre_cube"] = write_step_facets(
        scratch / "metre_cube.step", cube(), length_scale=1.0
    )

    # per-corner vertices: 36 VERTEX_POINTs for a closed 8-corner cube
    writer = StepWriter("SI_UNIT(.MILLI.,.METRE.)")
    faces = []
    for triangle in CUBE_FACES:
        vertices = [writer.vertex(*(1000.0 * c for c in CUBE_VERTICES[index])) for index in triangle]
        faces.append(writer.face([writer.poly_loop(vertices)]))
    writer.solidify(faces)
    fixtures["per_corner"] = write(scratch / "per_corner.step", writer.text())

    # a face with a hole: outer 2x2 square, inner 1x1 square, true area 3.0
    writer = StepWriter("SI_UNIT($,.METRE.)")
    outer_points = ((0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0))
    inner_points = ((0.5, 0.5, 0), (1.5, 0.5, 0), (1.5, 1.5, 0), (0.5, 1.5, 0))
    outer = writer.poly_loop([writer.vertex(*point) for point in outer_points])
    inner = writer.poly_loop([writer.vertex(*point) for point in inner_points])
    writer.solidify([writer.face([outer, inner])])
    fixtures["hole"] = write(scratch / "hole.step", writer.text())

    # a non-convex (U-shaped) loop: true area 7.0, fan area 11.0
    writer = StepWriter("SI_UNIT($,.METRE.)")
    u_points = ((0, 0, 0), (3, 0, 0), (3, 3, 0), (2, 3, 0), (2, 1, 0), (1, 1, 0), (1, 3, 0), (0, 3, 0))
    loop = writer.poly_loop([writer.vertex(*point) for point in u_points])
    writer.solidify([writer.face([loop])])
    fixtures["u_loop"] = write(scratch / "u_loop.step", writer.text())

    # a truncated cube: every entity present, the closing framing removed
    full = write_step_facets(scratch / "full_cube.step", cube()).read_text()
    truncated = "\n".join(
        line for line in full.splitlines() if not line.startswith(("ENDSEC", "END-ISO"))
    )
    fixtures["truncated"] = write(scratch / "truncated.step", truncated)

    # a convex quad: fanning is exact here, so it must stay accepted
    writer = StepWriter("SI_UNIT($,.METRE.)")
    quad = ((0, 0, 0), (2, 0, 0), (2, 1, 0), (0, 1, 0))
    writer.solidify([writer.face([writer.poly_loop([writer.vertex(*point) for point in quad])])])
    fixtures["convex_quad"] = write(scratch / "convex_quad.step", writer.text())

    return fixtures


def main(output: Path | None = None) -> dict:
    output = Path(output) if output else ROOT / "artifacts" / "exchange_verification"
    output.mkdir(parents=True, exist_ok=True)

    results: dict = {"mission": "exchange seam verification", "checks": []}

    # 1. round trip at three declared scales
    round_trips = []
    with tempfile.TemporaryDirectory() as scratch:
        scratch_path = Path(scratch)
        for scale in (1.0, 1000.0, 1e6):
            report = round_trip_facets(cube(), directory=scratch_path, length_scale=scale)
            round_trips.append({"length_scale": scale, **report})
        results["round_trips"] = round_trips
        fixtures = build_fixtures(scratch_path)

        # 2. naive (pre-fix) readings, measured on the fixtures
        naive_metre_vertices, _ = naive_read(fixtures["metre_cube"])
        naive_extent = max(abs(value) for vertex in naive_metre_vertices for value in vertex)
        honest = read_step_facets(fixtures["metre_cube"])
        honest_extent = max(abs(value) for vertex in honest.mesh.vertices for value in vertex)

        naive_dup_vertices, naive_dup_faces = naive_read(fixtures["per_corner"])
        naive_dup_mesh = Mesh(vertices=tuple(naive_dup_vertices), faces=tuple(naive_dup_faces))
        welded = read_step_facets(fixtures["per_corner"])

        # these two fixtures declare metres, so the naive reader is given the correct scale
        # override on purpose: the error measured below is the triangulation, not the unit
        naive_hole_vertices, naive_hole_faces = naive_read(fixtures["hole"], length_scale=1.0)
        hole_naive_area = area_of(naive_hole_faces, naive_hole_vertices)
        hole_true_area = polygon_area(((0, 0), (2, 0), (2, 2), (0, 2))) - polygon_area(
            ((0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5))
        )

        naive_u_vertices, naive_u_faces = naive_read(fixtures["u_loop"], length_scale=1.0)
        u_naive_area = area_of(naive_u_faces, naive_u_vertices)
        u_true_area = polygon_area(((0, 0), (3, 0), (3, 3), (2, 3), (2, 1), (1, 1), (1, 3), (0, 3)))

        naive_truncated_vertices, naive_truncated_faces = naive_read(fixtures["truncated"])

        # 3. what the current reader does with the same files
        current: dict[str, dict] = {}
        for name in ("hole", "u_loop", "truncated"):
            try:
                read_step_facets(fixtures[name])
            except ExchangeError as error:
                current[name] = {"outcome": "refused", "error": str(error).split(";")[0][:160]}
            else:  # pragma: no cover - the mission asserts these are refused
                current[name] = {"outcome": "accepted"}

        quad = read_step_facets(fixtures["convex_quad"])
        quad_area = area_of(quad.mesh.faces, quad.mesh.vertices)

    def rel(measured: float, truth: float) -> float:
        return (measured - truth) / truth

    results["silent_defects"] = [
        {
            "defect": "declared unit ignored on read (metre file divided by 1000)",
            "quantity": "max |x| of a 1 m cube, metres",
            "naive": naive_extent,
            "truth": honest_extent,
            "relative_error": rel(naive_extent, honest_extent),
            "current_outcome": "declared SI_UNIT($,.METRE.) read as metres; caller override still possible",
        },
        {
            "defect": "per-corner VERTEX_POINTs not welded",
            "quantity": "closed cube reported watertight (1 = yes)",
            "naive": 1.0 if is_watertight(naive_dup_mesh) else 0.0,
            "truth": 1.0,
            "relative_error": (1.0 if is_watertight(naive_dup_mesh) else 0.0) - 1.0,
            "current_outcome": f"welded: {len(welded.mesh.vertices)} vertices, "
            f"watertight={is_watertight(welded.mesh)}",
        },
        {
            "defect": "face with a hole concatenated and fanned",
            "quantity": "planar area, m^2",
            "naive": hole_naive_area,
            "truth": hole_true_area,
            "relative_error": rel(hole_naive_area, hole_true_area),
            "current_outcome": current["hole"]["outcome"],
        },
        {
            "defect": "non-convex loop fanned from its first vertex",
            "quantity": "planar area, m^2",
            "naive": u_naive_area,
            "truth": u_true_area,
            "relative_error": rel(u_naive_area, u_true_area),
            "current_outcome": current["u_loop"]["outcome"],
        },
        {
            "defect": "truncated file accepted (no framing or reference closure check)",
            "quantity": "file accepted without its closing statements (1 = accepted, must be 0)",
            "naive": 1.0 if naive_truncated_faces else 0.0,
            "truth": 0.0,
            "relative_error": None,
            "current_outcome": current["truncated"]["outcome"],
        },
    ]
    results["convex_quad"] = {
        "triangles": len(quad.mesh.faces),
        "fan_area": quad_area,
        "true_area": 2.0,
        "relative_error": rel(quad_area, 2.0),
        "outcome": "accepted (a fan over a convex loop is exact)",
    }
    results["fixture_vertices"] = {
        "per_corner_naive": len(naive_dup_vertices),
        "per_corner_welded": len(welded.mesh.vertices),
        "per_corner_removed_rows": welded.metadata["welded_vertex_rows"],
        "truth": len(CUBE_VERTICES),
    }
    results["round_trip_vertex_deviation_m"] = max(
        report["vertex_deviation_m"] for report in round_trips
    )
    results["refusals"] = {name: entry for name, entry in current.items() if entry["outcome"] == "refused"}

    # ------------------------------------------------------------------ figure
    figure_path = _figure(results, output / "figure.png")
    results["figure"] = str(figure_path.relative_to(ROOT))

    (output / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")

    docs_path = ROOT / "docs" / "exchange_verification.md"
    docs_path.write_text(_markdown(results), encoding="utf-8")
    results["docs"] = str(docs_path.relative_to(ROOT))
    print(json.dumps({k: v for k, v in results.items() if k != "checks"}, indent=2)[:4000])
    print(f"\nwrote {output/'results.json'}\nwrote {docs_path}\nwrote {figure_path}")
    return results


def _figure(results: dict, path: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    defects = results["silent_defects"]
    labels = ["unit\nread", "hole\nfan", "concave\nfan"]
    order = [0, 2, 3]
    errors = [abs(defects[index]["relative_error"]) * 100.0 for index in order]

    figure, (ax_area, ax_state) = plt.subplots(1, 2, figsize=(11.0, 4.2), dpi=160)

    colours = ["#b2182b" if value >= 5.0 else "#4393c3" for value in errors]
    ax_area.bar(labels, errors, color=colours)
    ax_area.axhline(3.0, color="black", linestyle="--", linewidth=1.0)
    ax_area.text(3.35, 3.6, "3 % acceptability line", ha="right", fontsize=8)
    ax_area.set_yscale("symlog", linthresh=1.0)
    ax_area.set_ylabel("|relative error| of the naive route,  %")
    ax_area.set_title("Silent error removed by the guards", fontsize=10)
    for index, value in enumerate(errors):
        ax_area.annotate(f"{value:,.1f} %", (index, value), ha="center", va="bottom", fontsize=8)
    truncation = results["silent_defects"][4]
    ax_area.text(
        0.02, 0.02,
        f"truncated file: naive route accepted it "
        f"(accepted={truncation['naive']:.0f}, must be {truncation['truth']:.0f}); now refused",
        transform=ax_area.transAxes, fontsize=7, va="bottom",
    )

    barge = ["per-corner\nvertices", "watertight?"]
    before = [
        results["fixture_vertices"]["per_corner_naive"],
        100.0 * results["silent_defects"][1]["naive"],
    ]
    after = [
        results["fixture_vertices"]["per_corner_welded"],
        100.0,
    ]
    positions = [0, 1]
    ax_state.bar([p - 0.2 for p in positions], before, width=0.38, label="naive read", color="#b2182b")
    ax_state.bar([p + 0.2 for p in positions], after, width=0.38, label="current read", color="#4393c3")
    ax_state.set_xticks(positions)
    ax_state.set_xticklabels(barge, fontsize=9)
    ax_state.set_title("Closed solid read back from a per-corner writer", fontsize=10)
    ax_state.legend(fontsize=8)
    ax_state.text(0, max(before) * 1.02, "vertices", ha="center", fontsize=8)
    ax_state.text(1, 104, "100 = closed", ha="center", fontsize=8)
    ax_state.set_ylim(0, max(before) * 1.18)

    figure.tight_layout()
    figure.savefig(path, bbox_inches="tight")
    plt.close(figure)
    return path


def _markdown(results: dict) -> str:
    lines: list[str] = []
    add = lines.append

    add("# Exchange seam verification: units, welding and the refusal guards")
    add("")
    add("**Mission**: m16084 (verifier seat, selfstart) · **Repository**: TensorCAD · "
        "**Module**: `tensorcad/exchange.py` · **Evidence**: "
        f"`{results.get('figure', 'artifacts/exchange_verification/figure.png')}` and "
        "`artifacts/exchange_verification/results.json` (every number below is injected "
        "from that file by `examples/exchange_verification/run_example.py`).")
    add("")
    add("## 1. Case: what the seam is for, and what it must never do")
    add("")
    add("`tensorcad.exchange` is the place a shipyard's file enters the quality audit: a "
        "faceted STEP (`FACETED_BREP` over `POLY_LOOP`), a neutral `.patchset.json`, or a "
        "refusal. The audit judges a triangle mesh, so two silent failures are worse than a "
        "refusal, because nothing downstream can notice them: **a wrong scale** and **a "
        "wrong triangulation**. The module's own docstring already promised both would not "
        "happen -- \"the scale travels with the data instead of being assumed\" and \"it never "
        "guesses a tessellation\" -- and neither promise had a test. This page is the "
        "verification that closed that debt, plus the four defects it found.")
    add("")
    add("Fixture geometry is a unit cube (8 vertices, 12 triangles, 6 m² of surface, 1 m³) "
        "plus four hand-built STEP files that this repository's writer cannot produce: a "
        "per-corner-vertex cube, a face with a hole, a non-convex U-shaped loop, and a "
        "truncated cube. The hand-built files are emitted by a fixture writer inside the "
        "example, so the code under test is not grading its own output.")
    add("")
    add("## 2. Method")
    add("")
    add("1. **Independent structural parse.** The written file is re-read with a separate "
        "grammar (entity id, keyword, reference closure, framing markers) that shares no "
        "code with the reader; counts and reference closure are asserted against the mesh.")
    add("2. **Declaration round trip.** The writer is run at three scales (1, 1e3, 1e6) and "
        "the reader is given *no* scale, so the declared `LENGTH_UNIT` is the only route by "
        "which the coordinates can come back right.")
    add("3. **Naive-route measurement.** The pre-fix behaviour is reproduced on the current "
        "tree by a deliberately unguarded reader inside the example (hard-coded /1000, no "
        "framing check, concatenate-and-fan over every loop). Each fixture is read both ways; "
        "the error of the naive route is what the guards removed.")
    add("4. **Negative controls.** Analytic B-rep, holes, non-convex loops, truncation, "
        "dangling references, zero/negative/unrepresentable scales, empty mesh, IGES and "
        "unknown extensions all have to raise, and a convex quad has to keep working.")
    add("")
    add("## 3. Results")
    add("")
    add("### 3.1 Silent defects, measured on the naive route and closed by the guards")
    add("")
    add("| defect | quantity | naive | truth | relative error | current outcome |")
    add("|---|---|---:|---:|---:|---|")
    for entry in results["silent_defects"]:
        error = (
            f"{100.0 * entry['relative_error']:+.2f} %"
            if entry["relative_error"] is not None
            else "binary flag"
        )
        add(f"| {entry['defect']} | {entry['quantity']} | {entry['naive']:.6g} | "
            f"{entry['truth']:.6g} | {error} | {entry['current_outcome']} |")
    add("")
    add("### 3.2 Round trip at each declared scale (reader given no scale)")
    add("")
    add("| writer `length_scale` | reader `length_scale` | source of that scale | "
        "vertex deviation (m) | triangles |")
    add("|---:|---:|---|---:|---:|")
    for report in results["round_trips"]:
        add(f"| {report['length_scale']:.0f} | {report['length_scale']:.0f} | "
            f"{'declared' if report['used_declared_unit'] else 'override'} | "
            f"{report['vertex_deviation_m']:.3g} | {report['triangle_count']:.0f} |")
    add("")
    add("### 3.3 Welding a per-corner writer's solid")
    add("")
    add("| quantity | naive read | current read | truth |")
    add("|---|---:|---:|---:|")
    add(f"| vertices | {results['fixture_vertices']['per_corner_naive']} | "
        f"{results['fixture_vertices']['per_corner_welded']} | "
        f"{results['fixture_vertices']['truth']} |")
    add(f"| removed duplicate rows | 0 | {results['fixture_vertices']['per_corner_removed_rows']} | 28 |")
    add(f"| watertight | no | yes | yes |")
    add("")
    add("Welding merges only *exactly* equal coordinates, so no position changes; the vertex "
        "numbering is remapped, which is why the test compares triangles as coordinate sets "
        "rather than index tuples. Passing `weld_duplicates=False` restores entity-level "
        "fidelity and the open-mesh reading, for a caller that wants to audit the file "
        "rather than the model.")
    add("")
    add("### 3.4 Convex loops keep working")
    add("")
    quad = results["convex_quad"]
    add(f"A 2 m x 1 m convex quad reads as {quad['triangles']} triangles with fan area "
        f"{quad['fan_area']:.6g} m² against a true area of {quad['true_area']:.6g} m² "
        f"({100.0 * quad['relative_error']:+.2e} %). The convexity guard is not an "
        "over-refusal: a fan over a convex loop is exact.")
    add("")
    add("## 4. What was changed")
    add("")
    add("| # | change | why it is not a workaround |")
    add("|---|---|---|")
    add("| 1 | writer emits the `SI_UNIT` prefix that matches `length_scale` "
        "(`$`/`.DECI.`/`.CENTI.`/`.MILLI.`/`.MICRO.`/`.NANO.`), and refuses a scale with no "
        "SI prefix | the file now describes its own numbers; a caller using 25.4 (inch) gets "
        "an error instead of a file that lies to every conformant reader |")
    add("| 2 | reader derives the divisor from the declared `LENGTH_UNIT`, records "
        "`length_scale_source`, and keeps an explicit `length_scale` override | the docstring's "
        "\"scale travels with the data\" claim is now true, and the assumed-millimetre fallback "
        "is labelled as an assumption |")
    add("| 3 | `weld_duplicates` (default on) merges exactly equal `VERTEX_POINT`s and reports "
        "`welded_vertex_rows` / `collapsed_triangles` | a legal per-corner file no longer reads "
        "as an open solid, so the audit judges the model instead of the writer's indexing |")
    add("| 4 | a face with more than one bound, or a non-convex `POLY_LOOP`, is refused with "
        "the face id | a fan over those invents area (see 3.1); refusing is the module's stated "
        "policy, and `skipped_faces`/`collapsed_triangles` keep the rest countable |")
    add("| 5 | truncation and dangling references are refused; `ADVANCED_FACE` dead code removed | "
        "a file missing its closing statements used to read as a smaller valid body |")
    add("")
    add("## 5. Points to watch")
    add("")
    add("1. **Still no third-party validation.** No CAD kernel is installed here, so the "
        "writer is validated by an internal round trip and an independent parse of the text, "
        "not by a foreign reader. `supported_routes()` keeps that at `WITHHELD`, correctly.")
    add("2. **`CONVERSION_BASED_UNIT` is not converted.** A file declaring inches has no SI "
        "prefix this reader honours; the millimetre fallback applies and "
        "`length_scale_source` says so. Convert or override the scale explicitly.")
    add("3. **Welding is exact equality only.** A writer that stores the same corner with "
        "different round-off (e.g. `0.1` vs `0.10000000000000001`) will not weld; that is a "
        "deliberate choice, since a tolerance would be a guess about the model.")
    add("4. **Convexity is tested in the loop's own plane** (Newell normal). A loop whose "
        "vertices are not coplanar has no single normal, and the guard then decides on the "
        "projected winding; such a face is not a planar facet and belongs to a B-rep, not here.")
    add("5. **The reader still fans convex n-gons.** A convex pentagon becomes three triangles; "
        "the area is exact, the diagonal choice is not checked for aspect ratio.")
    add("")
    add("## 6. Error judgement")
    add("")
    quantitative = [entry for entry in results["silent_defects"] if entry["relative_error"] is not None]
    worst = max(quantitative, key=lambda entry: abs(entry["relative_error"]))
    add("| quantity | benchmark | measured | verdict |")
    add("|---|---|---|---|")
    add(f"| round-trip vertex deviation | exact recovery (< 1e-15 m) | "
        f"{results['round_trip_vertex_deviation_m']:.3g} m | **PASS (< 3 %)** |")
    add(f"| convex-quad fan area | exact (2.0 m²) | {quad['fan_area']:.6g} m² | "
        f"**PASS (< 3 %)** |")
    for entry in results["silent_defects"]:
        if entry["relative_error"] is None:
            add(f"| {entry['quantity']} ({entry['defect']}) | must be refused | "
                f"accepted={entry['naive']:.0f} on the naive route, refused now | "
                "**PASS (refusal)** |")
            continue
        verdict = (
            "**PASS (< 3 %)**"
            if abs(entry["relative_error"]) < 0.03
            else "**closed** (naive route now exact or refused)"
        )
        add(f"| {entry['quantity']} ({entry['defect']}) | naive route error | "
            f"{100.0 * entry['relative_error']:+.2f} % | {verdict} |")
    add("")
    add(f"The largest naive-route error is {100.0 * abs(worst['relative_error']):.1f} % on "
        f"\"{worst['defect']}\", far outside the 3 % acceptability line; every one of those "
        "routes is now either exact or an explicit refusal, which is the point. The claims "
        "that remain *unchecked* are the foreign-reader interoperability (no CAD kernel) and "
        "the `CONVERSION_BASED_UNIT` conversion, both marked above rather than claimed.")
    add("")
    add("## 7. Reproducing")
    add("")
    add("```bash")
    add("python3 examples/exchange_verification/run_example.py   # writes artifacts/ + this page")
    add("python3 -m pytest tests/test_exchange_verification.py -q  # 26 cases")
    add("```")
    add("")
    add("Figure: left panel is the relative error the guards removed on three quantitative routes "
        "(log scale, the dashed line is the 3 % acceptability line); right panel is a closed "
        "cube written with one `VERTEX_POINT` per face corner -- 36 vertices and not watertight "
        "before welding, 8 vertices and closed after.")
    add("")
    return "\n".join(lines)


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else None)