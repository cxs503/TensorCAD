"""Tests for the validity checks that gate an export."""
from __future__ import annotations

import json

import pytest

from tensorcad import HullParameters, build_mesh, generate_hull
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.validation import (
    Check,
    DimensionTarget,
    ValidationReport,
    check_dimensions,
    inspect_geometry,
    merge,
    vertex_link_defects,
)


def box(
    origin: tuple[float, float, float], size: float, builder: MeshBuilder | None = None
) -> MeshBuilder:
    """Append an axis-aligned box's faces to ``builder`` so corners can be shared."""
    target = builder or MeshBuilder()
    x0, y0, z0 = origin
    corners = [
        (x0, y0, z0), (x0 + size, y0, z0), (x0 + size, y0 + size, z0), (x0, y0 + size, z0),
        (x0, y0, z0 + size), (x0 + size, y0, z0 + size),
        (x0 + size, y0 + size, z0 + size), (x0, y0 + size, z0 + size),
    ]
    ids = [target.add_vertex(point) for point in corners]
    for face in (
        (0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
        (2, 3, 7), (2, 7, 6), (1, 2, 6), (1, 6, 5), (0, 4, 7), (0, 7, 3),
    ):
        target.add_triangle(*[ids[index] for index in face])
    return target


def test_degenerate_face_vertices_are_tolerated_by_the_link_walk() -> None:
    sliver = Mesh(vertices=((0, 0, 0), (1, 0, 0), (1, 1, 0)), faces=((0, 1, 1), (0, 1, 2)))  # type: ignore[arg-type]
    defects = vertex_link_defects(sliver)
    assert len(defects) == 3, "every vertex of an open sliver fails the closed-ring test"


def test_generated_hull_passes_every_check() -> None:
    model = generate_hull(HullParameters(nx=21, nz=9))
    report = model.validation
    assert report.passed and not report.errors and not report.warnings
    assert len(report.checks) >= 20
    json.dumps(report.to_dict())


def test_inverted_orientation_is_an_error() -> None:
    """The defect this release fixes: winding that cancels the volume to zero."""
    from tensorcad import mesh as meshlib

    mesh, _ = build_mesh(HullParameters(nx=9, nz=5))
    assert meshlib.total_signed_volume(mesh) > 0
    inverted = Mesh(vertices=mesh.vertices, faces=tuple(face[::-1] for face in mesh.faces))
    assert meshlib.total_signed_volume(inverted) < 0

    checks = {check.name: check for check in inspect_geometry(inverted)}
    assert not checks["outward_orientation"].passed
    assert checks["outward_orientation"].severity == "error"
    assert checks["watertight"].passed, "inversion does not open the surface"
    report = merge(inspect_geometry(inverted))
    assert not report.passed and any("outward_orientation" in error for error in report.errors)

    # Mirrored half-sheets, the original bug: identical magnitudes, opposite signs.
    half = len(mesh.faces) // 2
    mirrored = Mesh(
        vertices=mesh.vertices,
        faces=mesh.faces[:half] + tuple(face[::-1] for face in mesh.faces[half:]),
    )
    assert abs(meshlib.total_signed_volume(mirrored)) < abs(meshlib.total_signed_volume(mesh))
    mirrored_checks = inspect_geometry(mirrored)
    assert not any(c.passed for c in mirrored_checks if c.name == "outward_orientation")


def test_open_surface_is_an_error() -> None:
    builder = box((0.0, 0.0, 0.0), 1.0)
    mesh = builder.build()
    opened = Mesh(vertices=mesh.vertices, faces=mesh.faces[:-2])
    failed = {check.name for check in inspect_geometry(opened) if not check.passed}
    assert {"watertight", "euler_characteristic_genus_zero"} <= failed


def test_non_manifold_edge_is_an_error() -> None:
    builder = box((0.0, 0.0, 0.0), 1.0)
    mesh = builder.build()
    doubled = Mesh(vertices=mesh.vertices, faces=(*mesh.faces, mesh.faces[0]))
    failed = {check.name for check in inspect_geometry(doubled) if not check.passed}
    assert "edge_manifold" in failed


def test_degenerate_faces_are_an_error() -> None:
    sliver = Mesh(vertices=((0, 0, 0), (1, 0, 0), (0, 1, 0)), faces=((0, 1, 2), (0, 0, 1)))  # type: ignore[arg-type]
    failed = {check.name for check in inspect_geometry(sliver) if not check.passed}
    assert "no_degenerate_faces" in failed
    assert "watertight" in failed


def test_two_separate_lumps_pass_closure_but_fail_the_characteristic() -> None:
    builder = box((0.0, 0.0, 0.0), 1.0)
    box((1.5, 1.5, 1.5), 1.0, builder)
    mesh = builder.build()
    checks = {check.name: check for check in inspect_geometry(mesh)}
    assert checks["watertight"].passed, "each lump is individually closed"
    assert not checks["euler_characteristic_genus_zero"].passed, "two components give V-E+F = 4"
    assert checks["vertex_manifold"].passed


def test_pinch_vertices_are_reported_as_a_warning() -> None:
    """Two cubes touching at one corner are closed but not vertex-manifold."""
    builder = box((0.0, 0.0, 0.0), 1.0)
    box((1.0, 1.0, 1.0), 1.0, builder)
    mesh = builder.build()
    defects = vertex_link_defects(mesh)
    assert defects == (6,), "the shared corner is where two rings meet"
    checks = {check.name: check for check in inspect_geometry(mesh)}
    assert not checks["vertex_manifold"].passed
    assert checks["vertex_manifold"].severity == "warning"
    report = merge([check for check in inspect_geometry(mesh) if check.name == "vertex_manifold"])
    assert report.passed, "a documented geometric degeneracy is a warning, never an error"
    assert any("vertex_manifold" in warning for warning in report.warnings)
    assert report.errors == []


def test_volume_beyond_the_bounding_box_is_caught() -> None:
    builder = box((0.0, 0.0, 0.0), 1.0)
    box((0.25, 0.25, 0.25), 0.5, builder)
    mesh = builder.build()
    checks = {check.name: check for check in inspect_geometry(mesh)}
    assert not checks["volume_within_bounding_box"].passed
    assert "self-intersection" in checks["volume_within_bounding_box"].detail


def test_a_beam_larger_than_the_mesh_has_is_an_error_not_an_allowance() -> None:
    """The sampling allowance covers a shortfall, never an over-wide hull."""
    model = generate_hull(HullParameters(beam=2.0, nx=11, nz=5))
    checks = {check.name: check for check in check_dimensions(model.mesh, DimensionTarget(10.0, 1.0, 0.8))}
    assert not checks["beam_matches_request"].passed


def test_dimension_checks_report_realised_versus_requested() -> None:
    model = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13))
    ok = check_dimensions(model.mesh, DimensionTarget(length=10.0, beam=2.0, draft=0.8))
    assert all(check.passed for check in ok)
    wrong = check_dimensions(model.mesh, DimensionTarget(length=11.0, beam=3.0, draft=1.0))
    failed = {check.name for check in wrong if not check.passed}
    assert failed == {
        "length_matches_request",
        "draft_matches_request",
        "beam_matches_request",
    }
    detail = next(check.detail for check in wrong if check.name == "length_matches_request")
    assert "10" in detail and "11" in detail


def test_beam_allowance_covers_sampling_but_not_a_wrong_shape() -> None:
    """A coarse grid legitimately misses the taper maximum; a wrong beam does not."""
    model = generate_hull(HullParameters(beam=2.0, nx=10, nz=5))
    realised = model.mesh.bounding_box().sizes[1]
    assert realised < 2.0
    checks = {check.name: check for check in check_dimensions(model.mesh, DimensionTarget(10.0, 2.0, 0.8))}
    assert checks["beam_matches_request"].passed, "the sampling deficit is allowed"
    assert "realised beam" in checks["beam_matches_request"].detail
    assert "grid deficit" in checks["beam_matches_request"].detail


def test_report_records_pass_and_fail_states() -> None:
    report = ValidationReport(
        (
            Check("ok_check", True, "nothing to say"),
            Check("broken", False, "bad"),
            Check("soft", False, "hm", severity="warning"),
        )
    )
    assert report.errors == ["broken: bad"]
    assert report.warnings == ["soft: hm"]
    assert not report.passed
    data = report.to_dict()
    assert data["passed"] is False
    assert {entry["name"] for entry in data["checks"]} == {"ok_check", "broken", "soft"}
    assert ValidationReport().passed, "no checks means nothing failed"
    assert ValidationReport().to_dict()["errors"] == []


def test_empty_mesh_is_caught_before_geometry_checks_run() -> None:
    checks = inspect_geometry(Mesh(vertices=(), faces=()))
    failed = {check.name for check in checks if not check.passed}
    assert failed == {"non_empty_mesh"}, "no vertices is vacuously finite, but not usable"
    with pytest.raises(ValueError, match="empty"):
        Mesh(vertices=(), faces=()).bounding_box()
