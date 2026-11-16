"""Open-boundary accounting: the four causes are separated on cases with a known answer.

Each fixture is built so that exactly one class can be correct, and the assertion is on
the class *and* on the length in metres -- a count-only test would pass on a
classification that mixes a 5 mm knuckle chord with a 20 m deck chord.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad.geomquality.open_boundary import (
    CLASS_OVERLAP,
    CLASS_PLANE_PREFIX,
    CLASS_SEAM,
    CLASS_UNEXPLAINED,
    classify_open_boundary,
    closure_verdict,
    point_triangle_distance,
)
from tensorcad.geomquality.trim_mesh import import_patch_set

ASSETS = Path("/root/tp-agent/assets/geometry/hulls")


def _quad(origin, edge_u, edge_v):
    """Two triangles for the parallelogram ``origin + s*edge_u + t*edge_v``."""
    ox, oy, oz = origin
    ux, uy, uz = edge_u
    vx, vy, vz = edge_v
    vertices = [
        (ox, oy, oz),
        (ox + ux, oy + uy, oz + uz),
        (ox + ux + vx, oy + uy + vy, oz + uz + vz),
        (ox + vx, oy + vy, oz + vz),
    ]
    return vertices, [(0, 1, 2), (0, 2, 3)]


def test_point_triangle_distance_matches_analytic_value():
    a, b, c = (0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0)
    # Directly above the centroid: the distance is the height.
    assert point_triangle_distance((2 / 3, 2 / 3, 3.0), a, b, c) == pytest.approx(3.0)
    # Outside the triangle and beyond a vertex: distance to that vertex.
    assert point_triangle_distance((-1.0, -1.0, 0.0), a, b, c) == pytest.approx(math.sqrt(2.0))
    # Degenerate (zero area) triangle falls back to its sides, not to NaN or a division
    # by zero.
    assert point_triangle_distance((0.0, 0.0, 1.0), a, (2.0, 0.0, 0.0),
                                   (4.0, 0.0, 0.0)) == pytest.approx(1.0)


def test_closed_tetrahedron_has_no_open_boundary():
    vertices = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
    faces = [(0, 2, 1), (0, 1, 3), (1, 2, 3), (2, 0, 3)]
    audit = classify_open_boundary(vertices, faces, tolerance=1e-3)
    assert audit["free_edges"] == 0
    assert audit["open_length"] == 0.0
    assert all(entry["edges"] == 0 for entry in audit["classes"].values())


def test_differently_sampled_shared_curve_is_classified_as_seam():
    """Patch A samples the shared line once, patch B twice: a T-junction seam.

    Nothing is stacked, so no chord may be an overlap: the shared line has to come out
    as ``seam`` and the rest as ``unexplained`` (the two patches have an outer boundary).
    """
    vertices = [
        (0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (2.0, 1.0, 0.0), (0.0, 1.0, 0.0),   # A
        (2.0, 0.5, 0.0), (3.0, 0.0, 0.0), (3.0, 1.0, 0.0),                     # B
    ]
    faces = [(0, 1, 2), (0, 2, 3),          # patch A: long chord 1-2 along x = 2
             (1, 5, 6), (1, 6, 4), (4, 6, 2)]  # patch B: chords 1-4 and 4-2
    audit = classify_open_boundary(vertices, faces, tolerance=1e-3)
    classes = audit["classes"]
    # Exactly the three chords on the shared line are a seam: A's long chord 1-2 (1 m)
    # and B's two short chords 1-4 and 4-2 (0.5 m each).  The outer boundary of both
    # patches is an unexplained opening, and no chord is stacked on a face interior.
    assert classes[CLASS_SEAM]["edges"] == 3
    assert classes[CLASS_SEAM]["length"] == pytest.approx(2.0)
    assert classes[CLASS_OVERLAP]["edges"] == 0
    assert classes[CLASS_UNEXPLAINED]["edges"] == 6
    assert audit["open_length"] == pytest.approx(
        sum(entry["length"] for entry in classes.values()))
    assert audit["family"]["seam"]["length"] == pytest.approx(
        classes[CLASS_SEAM]["length"])


def test_boundary_lying_inside_another_patch_is_classified_as_overlap():
    """Untrimmed patches stacked on each other: the inner patch's boundary lies on the
    outer patch's *interior*, which stitching cannot close (it needs a trim)."""
    outer_vertices, outer_faces = _quad((0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (0.0, 2.0, 0.0))
    inner_vertices, inner_faces = _quad((1.0, 0.5, 0.0), (2.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    offset = len(outer_vertices)
    vertices = outer_vertices + inner_vertices
    faces = outer_faces + [tuple(index + offset for index in face) for face in inner_faces]
    audit = classify_open_boundary(vertices, faces, tolerance=1e-6)
    classes = audit["classes"]
    # The inner patch contributes 4 chords (2+2+1+1 = 6 m) that sit inside the outer one.
    assert classes[CLASS_OVERLAP]["edges"] == 4
    assert classes[CLASS_OVERLAP]["length"] == pytest.approx(6.0)
    # The outer patch's boundary is explained by nothing: it is an unexplained opening.
    assert classes[CLASS_UNEXPLAINED]["edges"] == 4
    assert classes[CLASS_SEAM]["edges"] == 0
    witnesses = audit["unexplained_worst"]
    assert witnesses and math.isclose(witnesses[0]["length"], 4.0)


def test_declared_plane_explains_a_half_hull_cut():
    """A boundary on a declared plane is an opening, not a defect -- and the plane wins
    over a coincident seam chord (ordering is part of the contract)."""
    vertices, faces = _quad((0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (0.0, 2.0, 0.0))
    plane = ("y0", ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)), 1e-9)
    audit = classify_open_boundary(vertices, faces, tolerance=1e-6, planes=[plane])
    plane_key = CLASS_PLANE_PREFIX + "y0"
    assert audit["classes"][plane_key]["edges"] == 1          # chord 0-1 runs on y = 0
    assert audit["classes"][plane_key]["length"] == pytest.approx(4.0)
    assert audit["family"]["declared_plane"]["length"] == pytest.approx(4.0)
    assert audit["classes"][CLASS_UNEXPLAINED]["edges"] == 3  # x = 0, x = 4, y = 2


def test_plane_declaration_requires_every_sample_to_lie_on_the_plane():
    """A chord that merely crosses the plane is not an opening on it."""
    vertices = [(0.0, -1.0, 0.0), (4.0, 1.0, 0.0), (4.0, 1.0, 2.0), (0.0, -1.0, 2.0)]
    faces = [(0, 1, 2), (0, 2, 3)]
    plane = ("y0", ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)), 1e-9)
    audit = classify_open_boundary(vertices, faces, tolerance=1e-6, planes=[plane])
    assert audit["family"]["declared_plane"]["edges"] == 0


def test_rejects_non_positive_tolerance_and_zero_samples():
    vertices, faces = _quad((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    with pytest.raises(ValueError):
        classify_open_boundary(vertices, faces, tolerance=0.0)
    with pytest.raises(ValueError):
        classify_open_boundary(vertices, faces, tolerance=1e-3, samples=0)


def test_closure_verdict_is_decided_by_the_budgets():
    audit = {"family": {"unexplained": {"length": 0.4}, "overlap": {"length": 1.5}}}
    assert closure_verdict(audit, max_unexplained_length=0.5,
                           max_overlap_length=2.0)["passed"] is True
    assert closure_verdict(audit, max_unexplained_length=0.1,
                           max_overlap_length=2.0)["passed"] is False
    assert closure_verdict(audit, max_unexplained_length=0.5,
                           max_overlap_length=1.0)["passed"] is False


@pytest.mark.skipif(not (ASSETS / "JBC" / "jbc-boss.igs").exists(),
                    reason="publisher asset registry not mounted")
def test_real_publisher_file_reports_a_closure_budget():
    """A real IGES: every class is reported, the lengths add up, and the unexplained
    witnesses point at coordinates instead of only at a count."""
    source = ASSETS / "JBC" / "jbc-boss.igs"
    result = import_patch_set(source, closure_tolerance=5e-3,
                              closure_planes=[("y0", ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                                               1e-6)])
    audit = result.report["open_boundary"]
    assert audit["free_edges"] > 0
    assert audit["tolerance"] == pytest.approx(5e-3)
    for name in ("seam", "overlap", "declared_plane", "unexplained"):
        assert name in audit["family"]
    lengths = sum(entry["length"] for entry in audit["classes"].values())
    edges = sum(entry["edges"] for entry in audit["classes"].values())
    assert lengths == pytest.approx(audit["open_length"])
    assert edges == audit["free_edges"]
    for witness in audit["unexplained_worst"]:
        assert len(witness["from"]) == 3 and len(witness["to"]) == 3
    # The report has to survive the JSON round trip the CLI writes.
    assert json.loads(json.dumps(result.report, default=str))["open_boundary"]["free_edges"] \
        == audit["free_edges"]


def test_closure_audit_is_off_by_default():
    """The measurement is opt-in: an import without ``closure_tolerance`` must not pay
    for it and must not gain a report key that other code could mistake for a result."""
    source = ASSETS / "JBC" / "jbc-boss.igs"
    if not source.exists():
        pytest.skip("publisher asset registry not mounted")
    result = import_patch_set(source)
    assert result.report["open_boundary"] == {}