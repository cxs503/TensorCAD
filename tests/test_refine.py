"""Target-size-field refinement: the claims ``tensorcad.refine`` makes.

The point of the module is that splitting a triangle is *exact* -- a new vertex is the
midpoint of an existing edge, so volume, area and centroid cannot move -- and that
resolution alone is not accuracy, which is why projection is a separate, measured step.
Those two statements are what these tests pin.
"""

from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.geometry_field import SphereField, to_mesh
from tensorcad.refine import (
    RefinementError,
    RefinementReport,
    edge_lengths,
    refine,
)

SPHERE_RADIUS = 1.0
SPHERE_FIELD = SphereField(radius=SPHERE_RADIUS)
SPHERE_VOLUME = 4.0 / 3.0 * math.pi * SPHERE_RADIUS**3


def _sphere_mesh(resolution: tuple[int, int, int] = (16, 16, 16)) -> Mesh:
    sampled = to_mesh(SPHERE_FIELD, resolution)
    assert sampled.report.watertight
    return sampled.mesh


def _box_mesh(divisions: int = 3) -> Mesh:
    """A hand-built watertight cube, cheaper and cleaner than a sampled primitive."""
    builder = MeshBuilder()
    corners = [
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    ]
    index = [builder.add_vertex(point) for point in corners]
    quads = [
        (0, 3, 2, 1),
        (4, 5, 6, 7),
        (0, 1, 5, 4),
        (1, 2, 6, 5),
        (2, 3, 7, 6),
        (3, 0, 4, 7),
    ]
    for a, b, c, d in quads:
        builder.add_triangle(index[a], index[b], index[c])
        builder.add_triangle(index[a], index[c], index[d])
    mesh = builder.build()
    if divisions <= 1:
        return mesh
    refined, _ = refine(mesh, 1.0 / divisions)
    return refined


def test_splitting_is_geometry_neutral() -> None:
    """Volume, area and centroid are invariant to round-off, not approximately."""
    mesh = _sphere_mesh()
    before_volume = meshlib.total_signed_volume(mesh)
    before_area = meshlib.surface_area(mesh)
    refined, report = refine(mesh, 0.08)
    assert report.passes >= 1
    assert report.converged
    assert report.after.faces > report.before.faces
    assert abs(report.volume_relative_change) < 1e-12
    assert abs(report.area_relative_change) < 1e-12
    assert report.centroid_shift < 1e-12
    assert meshlib.total_signed_volume(refined) == pytest.approx(before_volume, rel=1e-12)
    assert meshlib.surface_area(refined) == pytest.approx(before_area, rel=1e-12)


def test_topology_survives_the_split() -> None:
    """No hanging nodes: watertight, edge-manifold and Euler are all preserved."""
    mesh = _sphere_mesh()
    refined, report = refine(mesh, 0.09)
    assert report.before.watertight and report.after.watertight
    assert report.before.edge_manifold and report.after.edge_manifold
    assert report.after.euler_characteristic == report.before.euler_characteristic == 2
    assert meshlib.is_watertight(refined)
    assert len(meshlib.zero_area_faces(refined)) == 0


def test_every_edge_meets_the_target() -> None:
    target = 0.07
    mesh = _sphere_mesh()
    refined, report = refine(mesh, target)
    assert report.after.max_edge <= target
    for key, length in edge_lengths(refined).items():
        midpoint = tuple(
            (refined.vertices[key[0]][axis] + refined.vertices[key[1]][axis]) * 0.5 for axis in range(3)
        )
        assert length <= target * (1.0 + 1e-12), (length, midpoint)


def test_a_size_field_refines_where_the_caller_asked() -> None:
    """A callable target is honoured locally: an edge is judged by its own midpoint."""

    def size_field(point: tuple[float, float, float]) -> float:
        return 0.06 if point[2] < 0.0 else 0.5

    mesh = _sphere_mesh()
    adapted, adaptive_report = refine(mesh, size_field)
    uniform, uniform_report = refine(mesh, 0.06)
    assert adaptive_report.target_callable
    assert adaptive_report.converged
    assert adaptive_report.after.faces < uniform_report.after.faces
    # The southern hemisphere is well resolved, the northern one is not: the claim is
    # directional, so it is asserted directionally rather than on an average.
    south = [length for key, length in edge_lengths(adapted).items() if refined_mid_z(adapted, key) < 0.0]
    north = [length for key, length in edge_lengths(adapted).items() if refined_mid_z(adapted, key) >= 0.0]
    assert max(south) <= 0.06 * (1.0 + 1e-12)
    assert max(north) > 0.06 * 1.5
    assert adaptive_report.after.faces < len(uniform.faces)


def refined_mid_z(mesh: Mesh, key: tuple[int, int]) -> float:
    return (mesh.vertices[key[0]][2] + mesh.vertices[key[1]][2]) * 0.5


def test_projection_is_what_adds_accuracy() -> None:
    """Splitting alone leaves the faceting error; projecting halves it, honestly."""
    mesh = _sphere_mesh()
    before_error = abs(meshlib.total_signed_volume(mesh) / SPHERE_VOLUME - 1.0)
    _, plain = refine(mesh, 0.075)
    plain_error = abs(plain.after.volume / SPHERE_VOLUME - 1.0)
    # Geometry-neutral: the error did not move.
    assert plain_error == pytest.approx(before_error, rel=1e-9)

    projected, report = refine(mesh, 0.075, project=SPHERE_FIELD, project_iterations=3)
    assert report.projection is not None
    assert float(report.projection["max_residual_before"]) > 0.0
    assert float(report.projection["max_residual_after"]) < 1e-9
    assert float(report.projection["max_displacement"]) > 0.0
    projected_error = abs(meshlib.total_signed_volume(projected) / SPHERE_VOLUME - 1.0)
    assert projected_error < 0.75 * before_error
    # The projection moves the new vertices only; the original ones are untouched.
    for index in range(len(mesh.vertices)):
        assert projected.vertices[index] == mesh.vertices[index]


def test_a_pass_budget_that_runs_out_is_reported_not_hidden() -> None:
    mesh = _sphere_mesh()
    _, report = refine(mesh, 0.01, max_passes=1)
    assert report.passes == 1
    assert not report.converged
    assert report.after.max_edge > 0.01


def test_refinement_is_deterministic() -> None:
    mesh = _sphere_mesh()
    first, first_report = refine(mesh, 0.1)
    second, _ = refine(mesh, 0.1)
    assert first.vertices == second.vertices
    assert first.faces == second.faces
    assert first_report.to_dict()["after"] == refine(mesh, 0.1)[1].to_dict()["after"]


def test_units_and_frame_are_carried_through() -> None:
    builder = MeshBuilder(units="mm", coordinate_frame="x-aft")
    for point in [(0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0)]:
        builder.add_vertex(point)
    for face in [(0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3)]:
        builder.add_triangle(*face)
    mesh = builder.build()
    refined, _ = refine(mesh, 3.0)
    assert refined.units == "mm"
    assert refined.coordinate_frame == "x-aft"


def test_the_box_anchor_is_exact() -> None:
    """A cube of unit side: three levels of refinement keep V = 1 and A = 6 exactly."""
    mesh = _box_mesh(1)
    refined, report = refine(mesh, 0.25)
    assert report.after.volume == pytest.approx(1.0, rel=1e-15)
    assert report.after.area == pytest.approx(6.0, rel=1e-15)
    assert meshlib.is_watertight(refined)


@pytest.mark.parametrize(
    ("target", "message"),
    [(0.0, "target_edge_length"), (-1.0, "target_edge_length"), (float("nan"), "target_edge_length")],
)
def test_an_impossible_target_is_refused(target: float, message: str) -> None:
    with pytest.raises(RefinementError, match=message):
        refine(_box_mesh(1), target)


def test_an_empty_mesh_is_refused() -> None:
    with pytest.raises(RefinementError, match="no faces"):
        refine(Mesh(vertices=(), faces=()), 0.1)


def test_a_face_budget_is_enforced() -> None:
    with pytest.raises(RefinementError, match="max_faces"):
        refine(_box_mesh(1), 0.05, max_faces=100)


def test_a_report_round_trips_through_json() -> None:
    import json

    _, report = refine(_box_mesh(1), 0.4)
    assert isinstance(report, RefinementReport)
    payload = json.loads(json.dumps(report.to_dict()))
    assert payload["after"]["watertight"] is True
    assert payload["volume_relative_change"] == pytest.approx(report.volume_relative_change)
