"""Opening classification and planar-cut capping (GAP-CAD-T1-IGES-PIPE).

The claim under test is narrow and decidable: a boundary loop whose vertices lie on a
plane to round-off is the file's own cut plane and may be capped with no invented
surface; a loop that wanders off any plane is a defect gap and must be left open and
reported.  The box references below have closed-form volumes, so a sealed mesh is
compared against arithmetic, not against another mesher.
"""
from __future__ import annotations

import math
from pathlib import Path

import pytest

from tensorcad.geomquality.openings import (
    DEFAULT_PLANAR_RATIO,
    boundary_loops,
    classify_openings,
    seal_planar_openings,
)
from tensorcad.geomquality.trim_mesh import import_cad, seal_stage

ASSET_ROOT = Path("/srv/tensor/agent/assets/geometry/hulls")
KCS = ASSET_ROOT / "KCS" / "KCS_hull_SVA.igs"


# --------------------------------------------------------------------------- #
# fixtures: closed-form boxes
# --------------------------------------------------------------------------- #
def _box(width: float, height: float, depth: float) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]:
    """Axis-aligned box as 12 triangles, corners at the origin and (w, h, d)."""
    corners = [
        (0.0, 0.0, 0.0), (width, 0.0, 0.0), (width, height, 0.0), (0.0, height, 0.0),
        (0.0, 0.0, depth), (width, 0.0, depth), (width, height, depth), (0.0, height, depth),
    ]
    quads = [
        (0, 3, 2, 1),  # z = 0
        (4, 5, 6, 7),  # z = depth
        (0, 1, 5, 4),  # y = 0
        (1, 2, 6, 5),  # x = width
        (2, 3, 7, 6),  # y = height
        (3, 0, 4, 7),  # x = 0
    ]
    faces: list[tuple[int, int, int]] = []
    for a, b, c, d in quads:
        faces.append((a, b, c))
        faces.append((a, c, d))
    return corners, faces


def _volume(vertices, faces) -> float:
    total = 0.0
    for a, b, c in faces:
        p, q, r = vertices[a], vertices[b], vertices[c]
        total += (p[0] * (q[1] * r[2] - q[2] * r[1])
                  - p[1] * (q[0] * r[2] - q[2] * r[0])
                  + p[2] * (q[0] * r[1] - q[1] * r[0])) / 6.0
    return total


def _open_box(width=2.0, height=3.0, depth=4.0):
    """The same box with its top face (y = height) removed, in original winding."""
    vertices, faces = _box(width, height, depth)
    kept = [f for f in faces if not ({2, 3, 7, 6} >= set(f))]
    return vertices, kept


# --------------------------------------------------------------------------- #
# 1. loop walking
# --------------------------------------------------------------------------- #
def test_boundary_loops_of_open_box_is_one_quad():
    vertices, faces = _open_box()
    loops = boundary_loops(vertices, faces)
    assert len(loops) == 1
    assert sorted(loops[0]) == [2, 3, 6, 7]


def test_boundary_loops_are_deterministic():
    vertices, faces = _open_box()
    first = boundary_loops(vertices, faces)
    second = boundary_loops(vertices, faces)
    assert first == second


def test_closed_box_has_no_boundary():
    vertices, faces = _box(2.0, 3.0, 4.0)
    assert boundary_loops(vertices, faces) == []
    report = classify_openings(vertices, faces)
    assert report.free_edges == 0
    assert report.verdict == "closed"


# --------------------------------------------------------------------------- #
# 2. classification: planar cut vs defect gap
# --------------------------------------------------------------------------- #
def test_open_box_loop_is_a_planar_cut():
    vertices, faces = _open_box()
    report = classify_openings(vertices, faces, planar_tolerance=1e-9)
    assert report.free_edges == 4
    assert len(report.loops) == 1
    loop = report.loops[0]
    assert loop.kind == "planar_cut"
    assert loop.planarity_residual == pytest.approx(0.0, abs=1e-12)
    assert sorted(abs(value) for value in loop.normal) == pytest.approx([0.0, 0.0, 1.0], abs=1e-9)
    assert report.verdict == "closable_by_planar_caps"


def test_displacing_one_loop_vertex_turns_it_into_a_defect_gap():
    vertices, faces = _open_box()
    vertices = list(vertices)
    vertices[6] = (vertices[6][0], vertices[6][1] + 0.5, vertices[6][2])
    report = classify_openings(vertices, faces)
    assert report.loops[0].kind == "defect_gap"
    # The least-squares plane tilts, so the worst vertex sits ~0.13 m off it --
    # two orders of magnitude above round-off, which is what makes it a gap.
    assert report.loops[0].planarity_residual > 0.1
    assert report.verdict == "defect_gaps_present"


def test_default_tolerance_scales_with_the_model():
    """A millimetre model and a metre model are judged by the same relative standard."""
    small_vertices, small_faces = _open_box(0.002, 0.003, 0.004)
    report = classify_openings(small_vertices, small_faces)
    assert report.planar_tolerance == pytest.approx(DEFAULT_PLANAR_RATIO * 0.004)
    assert report.loops[0].kind == "planar_cut"


# --------------------------------------------------------------------------- #
# 3. capping: exact volume, no invented surface
# --------------------------------------------------------------------------- #
def test_sealing_the_open_box_restores_the_analytic_volume():
    width, height, depth = 2.0, 3.0, 4.0
    vertices, faces = _open_box(width, height, depth)
    sealed_vertices, sealed_faces, record = seal_planar_openings(vertices, faces)
    assert record["loops_capped"] == 1
    assert record["triangles_added"] == 2
    assert record["capped_area"] == pytest.approx(width * depth, rel=1e-12)
    assert boundary_loops(sealed_vertices, sealed_faces) == []
    assert abs(_volume(sealed_vertices, sealed_faces)) == pytest.approx(
        width * height * depth, rel=1e-12)


def test_sealing_leaves_a_defect_gap_open_and_says_so():
    vertices, faces = _open_box()
    vertices = list(vertices)
    vertices[6] = (vertices[6][0], vertices[6][1] + 0.5, vertices[6][2])
    sealed_vertices, sealed_faces, record = seal_planar_openings(vertices, faces)
    assert record["loops_capped"] == 0
    assert record["triangles_added"] == 0
    assert len(sealed_faces) == len(faces)
    assert record["classification"]["verdict"] == "defect_gaps_present"


def test_surface_vertices_are_never_moved_by_sealing():
    vertices, faces = _open_box()
    sealed_vertices, _, _ = seal_planar_openings(vertices, faces)
    assert sealed_vertices == [tuple(float(v) for v in vertex) for vertex in vertices]


# --------------------------------------------------------------------------- #
# 4. pipeline stage and public import report
# --------------------------------------------------------------------------- #
def test_seal_stage_measures_without_touching_the_mesh_by_default():
    vertices, faces = _open_box()
    out_vertices, out_faces, record = seal_stage(vertices, faces, seal=False)
    assert record["applied"] is False
    assert out_faces == [tuple(face) for face in faces]
    assert record["analysis_before"]["verdict"] == "closable_by_planar_caps"


def test_seal_stage_caps_when_asked_and_reseals_idempotently():
    vertices, faces = _open_box()
    sealed_vertices, sealed_faces, record = seal_stage(vertices, faces, seal=True)
    assert record["applied"] is True
    assert record["analysis_after"]["verdict"] == "closed"
    again_vertices, again_faces, again = seal_stage(sealed_vertices, sealed_faces, seal=True)
    assert again["applied"] is False
    assert len(again_faces) == len(sealed_faces)


@pytest.mark.skipif(not KCS.exists(), reason="KCS publisher IGES not present")
def test_real_kcs_reports_a_mixed_boundary():
    """The publisher patch set is neither closed nor fully cappable -- and says which."""
    result = import_cad(KCS, seal_openings=False)
    record = result.report["openings"]
    analysis = record["analysis_before"]
    assert analysis["free_edges"] > 0
    assert analysis["verdict"] in ("mixed_planar_cuts_and_defect_gaps",
                                   "defect_gaps_present")
    # Some of the boundary really is the file's own cut plane (zero residual loops).
    assert analysis["planar_edges"] > 0
    planar = [loop for loop in analysis["loops"] if loop["kind"] == "planar_cut"]
    assert min(loop["planarity_residual"] for loop in planar) < 1e-9
    # The gates must still be honest: sealing planar cuts cannot close defect gaps.
    assert result.report["gates"]["G3_watertight"]["passed"] is False
