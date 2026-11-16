"""Mesh-backed SDF field: exactness, truncation semantics, voxel accuracy."""
from __future__ import annotations

import math
import random

import pytest

from tensorcad.geometry_field import FieldError, SphereField, to_mesh
from tensorcad.hull import HullParameters, analytic_volume, generate_hull
from tensorcad.mesh_field import MeshSDFField, voxelize_mesh
from tensorcad.spatial_index import TriangleGrid, brute_force_surface_distance


def _hull_mesh(nx: int = 41, nz: int = 13):
    return generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=nx, nz=nz)).mesh


def _sphere_mesh(resolution: int = 24):
    field = SphereField(center=(0.0, 0.0, 0.0), radius=1.0)
    return field, to_mesh(field, resolution).mesh


def test_nearest_surface_matches_brute_force_exactly() -> None:
    mesh = _hull_mesh()
    grid = TriangleGrid(mesh)
    random.seed(11)
    points = [
        (random.uniform(-3.0, 13.0), random.uniform(-2.0, 2.0), random.uniform(-3.0, 2.0))
        for _ in range(60)
    ]
    for point in points:
        hit = grid.nearest_surface(point)
        assert hit is not None
        assert hit.distance == pytest.approx(brute_force_surface_distance(mesh, point), rel=1e-12, abs=1e-12)


def test_nearest_surface_is_exact_outside_and_below_the_box() -> None:
    mesh = _hull_mesh()
    grid = TriangleGrid(mesh)
    for point in [(500.0, 0.0, 0.0), (5.0, 0.0, -5.0), (-50.0, 40.0, -3.0)]:
        hit = grid.nearest_surface(point)
        assert hit is not None
        assert hit.distance == pytest.approx(brute_force_surface_distance(mesh, point), rel=1e-12)


def test_band_truncates_instead_of_searching_the_interior() -> None:
    mesh = _hull_mesh()
    grid = TriangleGrid(mesh)
    interior = (5.0, 0.0, -0.4)
    exact = grid.nearest_surface(interior)
    truncated = grid.nearest_surface(interior, band=0.05)
    assert exact is not None and truncated is not None
    assert exact.distance > 0.05
    assert truncated.distance == pytest.approx(0.05)
    assert truncated.triangle == -1  # marker: "nothing within the band"
    # A near-surface point is still exact under the same band.
    near = exact.point
    outward = (near[0], near[1], near[2] - 0.01)
    hit = grid.nearest_surface(outward, band=0.05)
    assert hit is not None and hit.triangle >= 0
    assert hit.distance == pytest.approx(brute_force_surface_distance(mesh, outward), rel=1e-9)


def test_band_must_be_positive() -> None:
    grid = TriangleGrid(_hull_mesh())
    with pytest.raises(ValueError):
        grid.nearest_surface((1.0, 0.0, 0.0), band=0.0)


def test_mesh_field_reproduces_the_analytic_sphere_away_from_the_seam() -> None:
    field, mesh = _sphere_mesh(24)
    mesh_field = MeshSDFField(mesh)
    random.seed(5)
    worst = 0.0
    for _ in range(120):
        point = tuple(random.uniform(-1.5, 1.5) for _ in range(3))
        analytic = field.signed_distance([point])[0]
        meshed = mesh_field.signed_distance([point])[0]
        worst = max(worst, abs(analytic - meshed))
    # Discretisation of a 24-cell sphere is ~1 % of the radius; the SDF itself is
    # exact, so the gap is the mesh's, not the query's.
    assert worst < 0.03


def test_mesh_field_sign_inside_negative_outside_positive() -> None:
    field, mesh = _sphere_mesh(24)
    mesh_field = MeshSDFField(mesh)
    inside = mesh_field.signed_distance([(0.0, 0.0, 0.0)])[0]
    outside = mesh_field.signed_distance([(2.0, 0.0, 0.0)])[0]
    assert inside < 0.0 < outside
    assert mesh_field.contains([(0.0, 0.0, 0.0)]) == [True]
    assert mesh_field.contains([(2.0, 0.0, 0.0)]) == [False]


def test_mesh_field_metadata_and_bounds() -> None:
    field, mesh = _sphere_mesh(16)
    mesh_field = MeshSDFField(mesh)
    box = mesh_field.bounds()
    assert box is not None
    assert box.sizes[0] == pytest.approx(2.0, abs=0.2)
    params = mesh_field.parameters()
    assert params["triangles"] == len(mesh.faces)
    assert params["closed"] is True
    metadata = mesh_field.metadata()
    assert metadata["kind"] == "mesh_sdf"
    assert metadata["quantity"] == "signed_distance"
    with pytest.raises(FieldError):
        mesh_field.with_state(radius=2.0)


def test_mesh_field_rejects_empty_mesh() -> None:
    from tensorcad.geometry import Mesh
    with pytest.raises(FieldError):
        MeshSDFField(Mesh(vertices=(), faces=()))


def test_voxelised_mesh_volume_beats_the_3pct_bar() -> None:
    _field, mesh = _sphere_mesh(24)
    voxels = voxelize_mesh(mesh, (32, 32, 32))
    exact = 4.0 / 3.0 * math.pi
    rel = (voxels.volume() - exact) / exact
    assert abs(rel) < 0.03


def test_hull_mesh_voxel_volume_against_analytic_hull() -> None:
    mesh = _hull_mesh(101, 33)
    voxels = voxelize_mesh(mesh, (32, 16, 8))
    exact = analytic_volume(10.0, 2.0, 0.8)
    rel = (voxels.volume() - exact) / exact
    # At 32x16x8 the spacing is 0.31/0.13/0.10 m, so the discretisation of a 0.8 m
    # draft is the error source, not the field; the committed artifact
    # (benchmarks/perf_scale/artifacts/ladder_voxelize_mesh.json) records the exact
    # value and the refinement trend. The test pins sign convention + partial-fill
    # agreement with the analytic volume at a coarse-but-finite band.
    assert abs(rel) < 0.05
