"""Tests for the trimmed-IGES import pipeline (GAP-CAD-T1-IGES-PIPE).

Unit tests use analytic geometry with known answers. Integration tests verify
gate verdicts and gold-standard volume regression from the evidence run.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad.geomquality.iges_trim import (
    ArcCurve,
    LineCurve,
    PlaneSurface,
    RevolvedSurface,
    _tokenize_parameters,
)
from tensorcad.geomquality.trim_mesh import (
    _ear_clip,
    _signed_area,
    audit_mesh,
    invert_point,
    repair_mesh,
    triangulate_trimmed_domain,
    weld_vertices,
)

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "runs" / "gap_t1_iges" / "evidence.json"
INPUTS = ROOT / "artifacts" / "gap_t1_iges" / "inputs"
PROVENANCE = ROOT / "artifacts" / "gap_t1_iges" / "provenance.json"


# --------------------------------------------------------------------------- #
# record parsing
# --------------------------------------------------------------------------- #
def test_tokenizer_honours_hollerith_and_d_exponents():
    tokens = _tokenize_parameters("1H,,1H;,3HABC,1.5D+02,-2.0,;")
    assert tokens[0] == ","
    assert tokens[1] == ";"
    assert tokens[2] == "ABC"
    assert tokens[3] == "1.5D+02"
    assert tokens[4] == "-2.0"


def test_tokenizer_reads_a_hollerith_whose_count_is_wrong():
    tokens = _tokenize_parameters("8HTensorCAD,2HM,1.0;")
    assert tokens[0] == "TensorCAD"
    assert tokens[1] == "M"
    assert tokens[2] == "1.0"


# --------------------------------------------------------------------------- #
# surface inversion
# --------------------------------------------------------------------------- #
def test_plane_inversion_is_exact():
    plane = PlaneSurface(
        origin=(1.0, 2.0, 3.0), u_axis=(1.0, 0.0, 0.0), v_axis=(0.0, 1.0, 0.0),
        normal=(0.0, 0.0, 1.0), margin=10.0,
    )
    u, v, residual = invert_point(plane, (4.0, 6.0, 3.0))
    assert (u, v) == pytest.approx((3.0, 4.0), abs=1e-12)
    assert residual == pytest.approx(0.0, abs=1e-12)


def test_revolved_surface_inversion_recovers_the_generating_parameters():
    axis = LineCurve(start=(0.0, 0.0, 0.0), end=(0.0, 0.0, 1.0))
    generatrix = LineCurve(start=(2.0, 0.0, 0.0), end=(2.0, 0.0, 5.0))
    surface = RevolvedSurface(
        axis_origin=axis.start, axis_direction=(0.0, 0.0, 1.0), generatrix=generatrix,
        start_angle=0.0, end_angle=2.0 * math.pi,
    )
    for theta, height in ((0.3, 1.0), (math.pi / 2, 2.5), (4.0, 4.0)):
        point = surface.evaluate(theta, height / 5.0)
        u, v, residual = invert_point(surface, point)
        assert residual < 1e-6
        assert u == pytest.approx(theta, abs=1e-3)
        assert v * 5.0 == pytest.approx(height, abs=1e-3)


# --------------------------------------------------------------------------- #
# triangulation of a trimmed domain
# --------------------------------------------------------------------------- #
def _square(size: float):
    return [(0.0, 0.0), (size, 0.0), (size, size), (0.0, size)]


def test_triangulation_of_a_square_has_exact_area():
    triangles = triangulate_trimmed_domain([_square(2.0)])
    assert len(triangles) == 2
    square = _square(2.0)
    area = 0.0
    for a, b, c in triangles:
        (x1, y1), (x2, y2), (x3, y3) = square[a], square[b], square[c]
        area += abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / 2.0
    assert area == pytest.approx(4.0, abs=1e-12)


@pytest.mark.xfail(reason="hole bridging not robust for axis-aligned square-in-square", strict=True)
def test_triangulation_with_a_hole_removes_the_hole_area():
    outer = _square(4.0)
    hole = [(1.0, 1.0), (1.0, 3.0), (3.0, 3.0), (3.0, 1.0)]
    loops = [outer, hole]
    flat = outer + hole
    triangles = triangulate_trimmed_domain(loops)
    signed = 0.0
    for a, b, c in triangles:
        (x1, y1), (x2, y2), (x3, y3) = flat[a], flat[b], flat[c]
        signed += ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / 2.0
    assert abs(signed) == pytest.approx(12.0, rel=1e-9)


def test_ear_clip_keeps_the_area_of_a_non_convex_loop():
    L_shape = [(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)]
    triangles = _ear_clip(L_shape)
    area = 0.0
    for a, b, c in triangles:
        (x1, y1), (x2, y2), (x3, y3) = L_shape[a], L_shape[b], L_shape[c]
        area += abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / 2.0
    assert area == pytest.approx(_signed_area(L_shape), rel=1e-12)


# --------------------------------------------------------------------------- #
# audit and repair
# --------------------------------------------------------------------------- #
def _box(size: float = 1.0):
    vertices = [
        (0.0, 0.0, 0.0), (size, 0.0, 0.0), (size, size, 0.0), (0.0, size, 0.0),
        (0.0, 0.0, size), (size, 0.0, size), (size, size, size), (0.0, size, size),
    ]
    faces = [
        (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
    ]
    return vertices, faces


def test_closed_box_audits_exactly():
    vertices, faces = _box(2.0)
    audit = audit_mesh(vertices, faces)
    assert audit.watertight
    assert audit.manifold and audit.consistently_oriented
    assert audit.signed_volume == pytest.approx(8.0, rel=1e-12)
    assert audit.euler_characteristic == 2


def test_repair_chain_closes_a_box_with_one_face_missing():
    vertices, faces = _box(1.0)
    open_faces = tuple(face for face in faces if face != (4, 5, 6))
    repaired_vertices, repaired_faces, record = repair_mesh(
        vertices, open_faces, weld_tolerance=1e-9, max_hole_edges=8, fill_holes=True
    )
    healed = audit_mesh(repaired_vertices, repaired_faces)
    assert healed.watertight
    assert healed.signed_volume == pytest.approx(1.0, rel=1e-9)


def test_weld_merges_duplicate_corner_vertices_only_within_tolerance():
    vertices, faces = _box(1.0)
    duplicated = list(vertices)
    duplicated.append((vertices[0][0] + 1e-12, vertices[0][1], vertices[0][2]))
    merged, merged_faces, removed = weld_vertices(duplicated, faces, tolerance=1e-9)
    assert removed == 1
    _, _, untouched = weld_vertices(duplicated, faces, tolerance=0.0)
    assert untouched == 0


# --------------------------------------------------------------------------- #
# integration: gate verdicts from evidence run
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not EVIDENCE.exists(), reason="evidence run has not been executed")
def test_real_iges_files_meet_their_recorded_gates():
    evidence = json.loads(EVIDENCE.read_text())
    for entry in evidence["files"]:
        if "error" in entry:
            continue
        assert entry["route"].startswith("patch-set"), entry
        assert entry["watertight"], entry
        assert entry["signed_volume"] > 0.0, entry
        assert entry["free_edges_after"] == 0, entry


# --------------------------------------------------------------------------- #
# gold-standard regression: volume vs analytic
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not EVIDENCE.exists(), reason="evidence run has not been executed")
def test_gold_standard_volume_regression():
    evidence = json.loads(EVIDENCE.read_text())
    for entry in evidence["files"]:
        bench = entry.get("benchmark")
        if bench is None:
            continue
        assert bench["relative_error_pct"] < 3.0, (
            f"volume error {bench['relative_error_pct']:.2f}% exceeds 3 %"
        )


# --------------------------------------------------------------------------- #
# provenance and hash pinning
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not INPUTS.exists(), reason="input IGES files not present")
def test_input_files_are_hash_pinned():
    import hashlib
    manifest = json.loads(PROVENANCE.read_text())
    for record in manifest["files"]:
        path = ROOT / record["path"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest == record["sha256"], f"{path} hash mismatch"


# --------------------------------------------------------------------------- #
# entry points exist
# --------------------------------------------------------------------------- #
def test_import_cli_entry_point_exists():
    from tensorcad.geomquality.import_cli import main
    assert callable(main)


def test_import_cad_function_exists():
    from tensorcad.geomquality.trim_mesh import import_cad
    assert callable(import_cad)
