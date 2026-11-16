"""Symmetry closure of a half-hull import: analytic verification.

Publisher ship IGES files are commonly the starboard half of the hull, cut on the centre
plane ``y = 0``.  That boundary is a modelling cut, not a leak, so the watertight gate can
never pass on such a file however good the seam stage is -- measured on
``work/assets/kcs.igs`` 58.6 % of the surviving open-boundary length lies on the centre
plane (``artifacts/iges_seam/kcs_census_002.json``).  ``geomquality.symmetry`` mirrors the
half across that plane.

These tests verify the stage against **closed-form volumes**, not against the stage's own
bookkeeping: a half cylinder of radius 1 and height 2 has volume ``pi``, so its mirror must
be watertight with volume ``2*pi``; a half sphere of radius 1 must close to ``4*pi/3``.
The tolerances quoted are the tessellation's, measured, not assumed.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tensorcad.geomquality.symmetry import (  # noqa: E402
    detect_symmetry_plane,
    mirror_across_plane,
    symmetry_stage,
)
from tensorcad.geomquality.trim_mesh import audit_mesh  # noqa: E402


# --------------------------------------------------------------------------- #
# fixtures: analytic half solids, boundary on y = 0
# --------------------------------------------------------------------------- #
def half_cylinder(radius: float = 1.0, height: float = 2.0, segments: int = 96):
    """Solid half cylinder, ``y >= 0``, axis along +z, open only on ``y = 0``."""
    vertices: list[tuple[float, float, float]] = []
    for level in (0.0, height):
        for index in range(segments + 1):
            theta = math.pi * index / segments
            vertices.append((radius * math.cos(theta), radius * math.sin(theta), level))
    faces: list[tuple[int, int, int]] = []
    stride = segments + 1
    for index in range(segments):
        a, b = index, index + 1
        c, d = stride + index, stride + index + 1
        faces.append((a, b, d))
        faces.append((a, d, c))
    # half-disc caps, fanned from the centre point which itself lies on y = 0
    for level_offset, sign in ((0, -1), (stride, 1)):
        centre = len(vertices)
        vertices.append((0.0, 0.0, vertices[level_offset][2]))
        for index in range(segments):
            a, b = level_offset + index, level_offset + index + 1
            faces.append((centre, b, a) if sign < 0 else (centre, a, b))
    return vertices, faces


def half_sphere(radius: float = 1.0, meridians: int = 96, parallels: int = 48):
    """Solid half sphere, ``y >= 0``, open only on the great circle in ``y = 0``."""
    vertices: list[tuple[float, float, float]] = []
    for ring in range(1, parallels):
        phi = math.pi * ring / parallels
        for column in range(meridians + 1):
            theta = math.pi * column / meridians
            rho = radius * math.sin(phi)
            vertices.append((rho * math.cos(theta), rho * math.sin(theta),
                             radius * math.cos(phi)))
    north = len(vertices)
    vertices.append((0.0, 0.0, radius))
    south = len(vertices)
    vertices.append((0.0, 0.0, -radius))
    stride = meridians + 1
    faces: list[tuple[int, int, int]] = []
    for column in range(meridians):
        faces.append((north, column, column + 1))
    for ring in range(parallels - 2):
        upper, lower = ring * stride, (ring + 1) * stride
        for column in range(meridians):
            a, b = upper + column, upper + column + 1
            c, d = lower + column, lower + column + 1
            faces.append((a, d, b))
            faces.append((a, c, d))
    last = (parallels - 2) * stride
    for column in range(meridians):
        faces.append((south, last + column + 1, last + column))
    return vertices, faces


def _volume(vertices, faces) -> float:
    total = 0.0
    for a, b, c in faces:
        p, q, r = vertices[a], vertices[b], vertices[c]
        total += (p[0] * (q[1] * r[2] - q[2] * r[1])
                  - p[1] * (q[0] * r[2] - q[2] * r[0])
                  + p[2] * (q[0] * r[1] - q[1] * r[0])) / 6.0
    return total


# --------------------------------------------------------------------------- #
# detection
# --------------------------------------------------------------------------- #
def test_half_cylinder_detects_the_centre_plane():
    vertices, faces = half_cylinder()
    plane = detect_symmetry_plane(vertices, faces)
    assert plane is not None
    assert plane["axis_name"] == "y"
    assert plane["value"] == pytest.approx(0.0, abs=1e-12)
    assert plane["coverage"] > 0.99


def test_closed_mesh_offers_no_symmetry_plane():
    """A mesh with no open boundary must not be mirrored -- that would double the solid."""
    vertices, faces = half_sphere()
    closed_vertices, closed_faces, record = symmetry_stage(
        vertices, faces, enabled=True, weld_tolerance=1e-9, min_coverage=0.0)
    assert record["applied"] is True          # the half sphere does close
    assert detect_symmetry_plane(closed_vertices, closed_faces) is None
    full = symmetry_stage(closed_vertices, closed_faces, enabled=True, weld_tolerance=1e-9)
    assert full[2]["applied"] is False
    assert full[2]["detected"] is None


def test_a_mesh_that_crosses_its_candidate_plane_is_refused():
    """A full cylinder has its open boundary on two end planes and material on both sides.

    Reflecting it would produce a doubled solid, so the one-sidedness guard has to refuse
    the plane even though the end rims sit exactly on a bounding-box face.
    """
    vertices, faces = half_cylinder()
    # close the y = 0 cut with a degenerate wall so the mesh straddles the plane
    doubled = [(x, y, z) for (x, y, z) in vertices] + [(x, -y, z) for (x, y, z) in vertices]
    offset = len(vertices)
    doubled_faces = list(faces) + [(a + offset, c + offset, b + offset) for (a, b, c) in faces]
    plane = detect_symmetry_plane(doubled, doubled_faces)
    # the only surviving open boundary is the pair of end rims; y = 0 carries none of it
    assert plane is None or plane["axis_name"] != "y"


# --------------------------------------------------------------------------- #
# closure against closed-form volumes
# --------------------------------------------------------------------------- #
def test_half_cylinder_mirrors_to_a_watertight_cylinder():
    vertices, faces = half_cylinder(radius=1.0, height=2.0, segments=96)
    before = audit_mesh(vertices, faces)
    assert before.watertight is False

    out_vertices, out_faces, record = symmetry_stage(
        vertices, faces, enabled=True, weld_tolerance=1e-9)
    assert record["applied"] is True
    after = audit_mesh(out_vertices, out_faces)

    assert after.free_edge_count == 0
    assert after.watertight is True
    assert after.component_count == 1 and after.closed_component_count == 1
    assert after.consistently_oriented is True
    exact = math.pi * 1.0 ** 2 * 2.0
    # a 96-segment chordal approximation of a circle under-estimates the area by
    # (1 - sin(2*pi/96)/(2*pi/96)) ~ 0.07 %, which is the tessellation's error, not the
    # mirror's: the mirrored half is exactly twice the half it came from
    assert _volume(out_vertices, out_faces) == pytest.approx(exact, rel=2e-3)
    assert record["free_edges_after"] == 0
    assert record["triangles_mirrored"] == len(faces)


def test_half_sphere_mirrors_to_a_watertight_sphere():
    vertices, faces = half_sphere(radius=1.0)
    out_vertices, out_faces, record = symmetry_stage(
        vertices, faces, enabled=True, weld_tolerance=1e-9)
    assert record["applied"] is True
    after = audit_mesh(out_vertices, out_faces)
    assert after.watertight is True
    assert after.free_edge_count == 0
    exact = 4.0 / 3.0 * math.pi
    assert _volume(out_vertices, out_faces) == pytest.approx(exact, rel=3e-3)


def test_measurement_only_does_not_change_the_mesh():
    vertices, faces = half_cylinder(segments=32)
    out_vertices, out_faces, record = symmetry_stage(
        vertices, faces, enabled=False, weld_tolerance=1e-9)
    assert record["applied"] is False
    assert record["detected"]["axis_name"] == "y"
    assert out_vertices == list(vertices)
    assert out_faces == [tuple(face) for face in faces]


def test_mirror_is_exactly_twice_the_half_it_came_from():
    """No loft: the mirrored volume must be 2x the half's own signed volume."""
    vertices, faces = half_cylinder(segments=48)
    out_vertices, out_faces, _record = mirror_across_plane(
        vertices, faces, plane=(1, 0.0), weld_tolerance=1e-9)
    half = _volume(vertices, faces)
    assert _volume(out_vertices, out_faces) == pytest.approx(2.0 * half, rel=1e-9)