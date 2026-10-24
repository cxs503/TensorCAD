"""Mesh conditioning: the invariants it promises, on meshes small enough to check by hand.

The tests are split in two: exact, hand-checkable cases on synthetic meshes (a cube, a
single triangle strip, a deliberately illegal collapse) and one integration case on the
SUBOFF union, where the numbers that matter are the ones a solver would ask for. The
integration thresholds are the measured ones from `examples/suboff_conditioning`, and
they are one-sided where the measurement is: a smoother that suddenly improves a mesh
by more than this has changed the shape, not the mesh.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.geometry import Mesh
from tensorcad.geometry_field import SphereField, to_mesh
from tensorcad.mesh import (
    edge_face_counts,
    is_edge_manifold,
    is_watertight,
    surface_area,
    total_signed_volume,
)
from tensorcad.mesh_conditioning import (
    ConditionedMesh,
    MeshConditioningError,
    condition_mesh,
    decimate,
    measure_mesh,
    project_to_field,
    taubin_smooth,
    vertex_normals,
)
from tensorcad.validation import inspect_mesh_quality, measure_quality


def _cube() -> Mesh:
    """A unit cube as 12 triangles, consistently wound, one shell."""
    vertices = (
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (1.0, 1.0, 1.0),
        (0.0, 1.0, 1.0),
    )
    faces = (
        (0, 2, 1), (0, 3, 2),  # bottom (z = 0), wound so the normal points -z
        (4, 5, 6), (4, 6, 7),  # top
        (0, 1, 5), (0, 5, 4),  # front (y = 0)
        (1, 2, 6), (1, 6, 5),  # right (x = 1)
        (2, 3, 7), (2, 7, 6),  # back (y = 1)
        (3, 0, 4), (3, 4, 7),  # left (x = 0)
    )
    return Mesh(vertices=vertices, faces=faces, units="m", coordinate_frame="unit")


def _uv_sphere(rings: int = 12, segments: int = 24, radius: float = 1.0) -> Mesh:
    """A UV sphere with shared poles: closed, manifold, genus 0."""
    vertices: list[tuple[float, float, float]] = [(0.0, 0.0, radius)]
    for ring in range(1, rings):
        phi = math.pi * ring / rings
        for segment in range(segments):
            theta = 2.0 * math.pi * segment / segments
            vertices.append(
                (
                    radius * math.sin(phi) * math.cos(theta),
                    radius * math.sin(phi) * math.sin(theta),
                    radius * math.cos(phi),
                )
            )
    vertices.append((0.0, 0.0, -radius))
    bottom = len(vertices) - 1

    def index(ring: int, segment: int) -> int:
        return 1 + (ring - 1) * segments + (segment % segments)

    faces: list[tuple[int, int, int]] = []
    for segment in range(segments):
        # Wound outward: the pole is at +z, so the top cap runs anticlockwise seen
        # from outside, and total_signed_volume() is positive.
        faces.append((0, index(1, segment), index(1, segment + 1)))
    for ring in range(1, rings - 1):
        for segment in range(segments):
            a, b = index(ring, segment), index(ring, segment + 1)
            c, d = index(ring + 1, segment + 1), index(ring + 1, segment)
            faces.append((a, c, b))
            faces.append((a, d, c))
    for segment in range(segments):
        faces.append((bottom, index(rings - 1, segment + 1), index(rings - 1, segment)))
    return Mesh(vertices=tuple(vertices), faces=tuple(faces), units="m")


# --------------------------------------------------------------------------- #
# Measurement
# --------------------------------------------------------------------------- #


def test_measure_mesh_reads_the_same_numbers_as_the_quality_checks() -> None:
    mesh = _cube()
    measured = measure_mesh(mesh)
    quality = measure_quality(mesh)
    assert measured.triangles == 12
    assert measured.watertight and measured.edge_manifold
    assert measured.zero_area_faces == 0
    assert measured.signed_volume == pytest.approx(1.0, rel=1e-12)
    assert measured.surface_area == pytest.approx(6.0, rel=1e-12)
    assert measured.min_angle_degrees == pytest.approx(quality.min_angle_degrees, rel=1e-12)
    assert measured.max_aspect_ratio == pytest.approx(quality.max_aspect_ratio, rel=1e-12)
    assert measured.sliver_share == pytest.approx(quality.sliver_share, rel=1e-12)
    assert measured.passed


def test_the_shared_quality_measure_is_what_the_checks_report() -> None:
    """One definition of \"good triangle\": the report and the conditioning agree."""
    mesh = to_mesh(SphereField(), 16, margin=2).mesh
    measured = measure_quality(mesh)
    checks = {check.name: check for check in inspect_mesh_quality(mesh)}
    assert f"{measured.min_angle_degrees:.3f}" in checks["mesh_min_angle"].detail
    assert f"{measured.max_aspect_ratio:.3g}" in checks["mesh_aspect_ratio"].detail
    assert f"{measured.sliver_count}/{measured.triangle_count}" in checks["mesh_sliver_share"].detail


# --------------------------------------------------------------------------- #
# Smoothing
# --------------------------------------------------------------------------- #


def _grid(size: int = 8) -> Mesh:
    """A flat triangulated patch in the z = 0 plane."""
    vertices = tuple(
        (column / size, row / size, 0.0) for row in range(size + 1) for column in range(size + 1)
    )
    faces: list[tuple[int, int, int]] = []
    for row in range(size):
        for column in range(size):
            a = row * (size + 1) + column
            b = a + 1
            c = a + (size + 1)
            d = c + 1
            faces.append((a, b, d))
            faces.append((a, d, c))
    return Mesh(vertices=vertices, faces=tuple(faces))


def test_smoothing_a_plane_keeps_it_flat() -> None:
    """A Laplacian is tangential on a plane, so the lam/mu pair cannot leave it."""
    mesh = _grid()
    smoothed, moved = taubin_smooth(mesh, iterations=8)
    assert moved > 0.0  # the boundary does move in-plane
    assert all(abs(vertex[2]) <= 1e-15 for vertex in smoothed.vertices)
    # It slides *within* the plane: the one-sided diagonals of this triangulation make
    # the in-plane Laplacian asymmetric, so the drift along the diagonal is expected --
    # the invariant is that z stays exactly zero.
    assert max(abs(vertex[0]) for vertex in smoothed.vertices) <= 1.0 + 1e-12
    assert max(abs(vertex[1]) for vertex in smoothed.vertices) <= 1.0 + 1e-12


def test_smoothing_is_shrink_free_on_a_sphere() -> None:
    """A plain Laplacian shrinks a sphere; the Taubin pair does not."""
    mesh = _uv_sphere()
    radius = sum(vertex[2] ** 2 for vertex in mesh.vertices)
    plain = mesh
    numpy = pytest.importorskip("numpy")
    points = numpy.asarray(mesh.vertices)
    from tensorcad.mesh_conditioning import _neighbour_csr

    indptr, indices = _neighbour_csr(len(points), mesh.faces)
    counts = numpy.diff(indptr)
    for _ in range(20):
        average = numpy.add.reduceat(points[indices], indptr[:-1], axis=0) / counts[:, None]
        points = points + 0.5 * (average - points)
    laplacian_volume = abs(total_signed_volume(Mesh(vertices=tuple(map(tuple, points)), faces=plain.faces)))
    smoothed, moved = taubin_smooth(mesh, iterations=20)
    taubin_volume = abs(total_signed_volume(smoothed))
    reference = abs(total_signed_volume(mesh))
    assert moved > 0.0
    assert abs(taubin_volume / reference - 1.0) < 0.35 * abs(laplacian_volume / reference - 1.0)
    assert radius > 0.0


def test_smoothing_keeps_the_topology_and_the_units() -> None:
    mesh = _uv_sphere()
    smoothed, _ = taubin_smooth(mesh, iterations=10)
    assert is_watertight(smoothed) and is_edge_manifold(smoothed)
    assert smoothed.faces == mesh.faces
    assert smoothed.units == mesh.units
    assert smoothed.coordinate_frame == mesh.coordinate_frame


def test_smoothing_validates_its_parameters() -> None:
    mesh = _cube()
    with pytest.raises(MeshConditioningError):
        taubin_smooth(mesh, iterations=-1)
    with pytest.raises(MeshConditioningError):
        taubin_smooth(mesh, lam=0.0)
    with pytest.raises(MeshConditioningError):
        taubin_smooth(mesh, mu=0.5)
    assert taubin_smooth(mesh, iterations=0) == (mesh, 0.0)


# --------------------------------------------------------------------------- #
# Projection
# --------------------------------------------------------------------------- #


def test_the_vertex_normal_is_the_area_weighted_one() -> None:
    numpy = pytest.importorskip("numpy")
    cube = _cube()
    normals = vertex_normals(cube)
    assert numpy.allclose(numpy.linalg.norm(normals, axis=1), 1.0)
    # Every cube normal points away from the centre. The corner at the origin is *not*
    # (-1,-1,-1)/sqrt(3): the twelve-triangle cube gives that corner two triangles on
    # the bottom, two on the front and one on the left, so the area-weighted average is
    # (-1,-2,-2)/3 -- which is the point of weighting by area rather than by face count.
    centre = numpy.array([0.5, 0.5, 0.5])
    offsets = numpy.asarray(cube.vertices) - centre
    assert numpy.all(numpy.einsum("ij,ij->i", normals, offsets) > 0.0)
    assert normals[0] == pytest.approx(numpy.array([-1.0, -2.0, -2.0]) / 3.0, abs=1e-12)
    # A flat patch has one normal, exactly.
    flat = vertex_normals(_grid(4))
    assert numpy.allclose(flat, numpy.array([0.0, 0.0, 1.0]))


def test_projection_puts_a_perturbed_sphere_back_on_its_surface() -> None:
    numpy = pytest.importorskip("numpy")
    field = SphereField(radius=1.0)
    mesh = _uv_sphere()
    rng = numpy.random.default_rng(4)
    noisy = numpy.asarray(mesh.vertices) * (1.0 + 0.02 * rng.standard_normal(len(mesh.vertices)))[:, None]
    perturbed = Mesh(vertices=tuple(map(tuple, noisy)), faces=mesh.faces)
    projected, moved = project_to_field(perturbed, field, iterations=3)
    radii = numpy.linalg.norm(numpy.asarray(projected.vertices), axis=1)
    assert numpy.abs(radii - 1.0).max() < 2e-3, "projection did not reach the surface"
    assert 0.0 < moved < 0.1
    # The residual is monotone: more iterations never make it worse.
    two, _ = project_to_field(perturbed, field, iterations=2)
    residual_two = numpy.abs(numpy.linalg.norm(numpy.asarray(two.vertices), axis=1) - 1.0).max()
    assert numpy.abs(radii - 1.0).max() <= residual_two + 1e-15


def test_projection_is_monotone_and_cannot_be_thrown_off_by_a_thin_feature() -> None:
    """The line search is what stops the SUBOFF fins' 0.1 ft wedge from oscillating."""
    numpy = pytest.importorskip("numpy")
    field = SphereField(radius=1.0)
    mesh = _uv_sphere(rings=8, segments=16)
    far = numpy.asarray(mesh.vertices) * 1.6  # every vertex 60 % too far out
    perturbed = Mesh(vertices=tuple(map(tuple, far)), faces=mesh.faces)
    projected, moved = project_to_field(perturbed, field, iterations=1)
    before = numpy.abs(numpy.linalg.norm(far, axis=1) - 1.0)
    after = numpy.abs(numpy.linalg.norm(numpy.asarray(projected.vertices), axis=1) - 1.0)
    assert numpy.all(after <= before + 1e-12), "a projection step made some vertex worse"
    assert moved <= 0.6 + 1e-9  # clamped to the local edge length, not to the extent


def test_projection_validates_its_parameters() -> None:
    mesh = _cube()
    field = SphereField()
    with pytest.raises(MeshConditioningError):
        project_to_field(mesh, field, iterations=-1)
    with pytest.raises(MeshConditioningError):
        project_to_field(mesh, field, max_step=0.0)
    with pytest.raises(MeshConditioningError):
        project_to_field(mesh, field, method="magic")
    assert project_to_field(mesh, field, iterations=0) == (mesh, 0.0)
    assert project_to_field(Mesh(vertices=(), faces=()), field) == (
        Mesh(vertices=(), faces=()),
        0.0,
    )


def test_the_gradient_method_agrees_with_the_normal_method_on_a_smooth_field() -> None:
    numpy = pytest.importorskip("numpy")
    field = SphereField(radius=1.0)
    mesh = _uv_sphere()
    perturbed = Mesh(
        vertices=tuple((x * 1.01, y * 1.01, z * 1.01) for x, y, z in mesh.vertices),
        faces=mesh.faces,
    )
    by_normal, moved_normal = project_to_field(perturbed, field, iterations=2, method="normal")
    by_gradient, moved_gradient = project_to_field(
        perturbed, field, iterations=2, method="gradient"
    )
    normal_radii = numpy.linalg.norm(numpy.asarray(by_normal.vertices), axis=1)
    gradient_radii = numpy.linalg.norm(numpy.asarray(by_gradient.vertices), axis=1)
    assert numpy.abs(normal_radii - gradient_radii).max() < 1e-6
    assert moved_normal == pytest.approx(moved_gradient, rel=0.05)


# --------------------------------------------------------------------------- #
# Decimation
# --------------------------------------------------------------------------- #


def test_decimation_preserves_the_topology_of_a_sphere() -> None:
    mesh = _uv_sphere(rings=16, segments=32)
    before = len(mesh.faces)
    decimated, stats = decimate(mesh, target_faces=400)
    assert len(decimated.faces) <= before
    assert len(decimated.faces) >= 400
    assert stats["collapsed"] > 0
    assert is_watertight(decimated) and is_edge_manifold(decimated)
    # A genus-0 closed surface stays genus 0: V - E + F == 2.
    edges = len(edge_face_counts(decimated))
    euler = len(decimated.vertices) - edges + len(decimated.faces)
    assert euler == 2
    # The volume of a sphere survives a collapse-based coarsening to within a percent.
    assert total_signed_volume(decimated) == pytest.approx(
        total_signed_volume(mesh), rel=0.02
    )


def test_decimation_removes_short_edges_and_leaves_long_ones() -> None:
    """The sliver rule: only edges under the threshold may collapse."""
    mesh = _uv_sphere(rings=16, segments=32)
    lengths = sorted(
        math.dist(mesh.vertices[a], mesh.vertices[b]) for a, b in edge_face_counts(mesh)
    )
    threshold = 0.5 * lengths[len(lengths) // 2]
    decimated, stats = decimate(mesh, max_edge_length=threshold)
    assert stats["collapsed"] > 0
    shortest_before = lengths[0]
    lengths_after = sorted(
        math.dist(decimated.vertices[a], decimated.vertices[b])
        for a, b in edge_face_counts(decimated)
    )
    assert lengths_after[0] > shortest_before
    assert is_watertight(decimated)


def test_decimation_refuses_a_collapse_that_would_make_a_needle() -> None:
    """The quality gate is what stops QEM from trading one bad triangle for another."""
    mesh = _uv_sphere(rings=10, segments=20)
    lenient, lenient_stats = decimate(mesh, target_faces=60)
    strict, strict_stats = decimate(mesh, target_faces=60, max_aspect_ratio=2.0)
    assert strict_stats["rejected_quality"] > lenient_stats["rejected_quality"]
    assert len(strict.faces) >= len(lenient.faces)
    assert is_watertight(strict) and is_watertight(lenient)


def test_decimation_into_a_field_keeps_the_result_on_the_surface() -> None:
    numpy = pytest.importorskip("numpy")
    field = SphereField(radius=1.0)
    mesh = _uv_sphere(rings=12, segments=24)
    decimated, stats = decimate(mesh, target_faces=200, field=field)
    radii = numpy.linalg.norm(numpy.asarray(decimated.vertices), axis=1)
    assert stats["projected"] == len(decimated.vertices)
    assert numpy.abs(radii - 1.0).max() < 5e-3


def test_decimation_validates_its_parameters() -> None:
    mesh = _cube()
    with pytest.raises(MeshConditioningError):
        decimate(mesh)
    with pytest.raises(MeshConditioningError):
        decimate(mesh, target_faces=2)
    with pytest.raises(MeshConditioningError):
        decimate(mesh, max_edge_length=0.0)
    with pytest.raises(MeshConditioningError):
        decimate(mesh, target_faces=6, length_weight=-1.0)
    empty, stats = decimate(Mesh(vertices=(), faces=()), target_faces=10)
    assert empty.faces == () and stats["collapsed"] == 0


def test_decimation_stops_at_the_target() -> None:
    """It used to stop when the queue drained: a stale-version bug, pinned here."""
    mesh = _uv_sphere(rings=18, segments=36)
    decimated, _ = decimate(mesh, target_faces=500)
    assert 500 <= len(decimated.faces) <= 520


# --------------------------------------------------------------------------- #
# The whole step
# --------------------------------------------------------------------------- #


def test_conditioning_a_sphere_makes_it_more_spherical() -> None:
    """The mesh is an inscribed polyhedron; relaxing it moves *towards* the sphere.

    That is the useful direction to check: a smoother that shrank the shape would move
    away from 4/3 pi r^3, and one that merely rearranged triangles in the tangent plane
    would not move at all.
    """
    exact = 4.0 / 3.0 * math.pi
    result = condition_mesh(
        _uv_sphere(rings=14, segments=28),
        smooth_iterations=10,
        volume_tolerance=0.05,
        area_tolerance=0.05,
    )
    assert isinstance(result, ConditionedMesh)
    assert result.passed
    assert result.after.triangles == result.before.triangles
    assert result.after.sliver_share <= result.before.sliver_share
    before_error = abs(result.before.signed_volume - exact)
    after_error = abs(result.after.signed_volume - exact)
    assert after_error <= before_error, "smoothing moved the sphere the wrong way"
    assert result.after.surface_area == pytest.approx(
        result.before.surface_area, rel=0.05
    )
    summary = result.summary()
    assert summary["schema"] == "tensorcad.mesh.conditioning/1"
    assert summary["passed"] is True
    assert summary["before"]["triangles"] == summary["after"]["triangles"]


def test_conditioning_refuses_to_hand_over_a_shape_that_moved() -> None:
    """Decimating a coarse sphere to 20 faces has to change the sphere; say so."""
    mesh = _uv_sphere(rings=8, segments=16)
    with pytest.raises(MeshConditioningError) as excinfo:
        condition_mesh(mesh, target_faces=20, smooth_iterations=2, volume_tolerance=1e-6)
    assert "volume changed" in str(excinfo.value)


def test_conditioning_with_nothing_to_do_is_a_no_op() -> None:
    result = condition_mesh(_cube(), smooth_iterations=0)
    assert result.after.watertight
    assert result.after.triangles == 12
    assert result.mesh.vertices == _cube().vertices
    assert result.max_projection_move == 0.0
    assert result.collapsed_edges == 0


# --------------------------------------------------------------------------- #
# Integration: the SUBOFF union, where the numbers are the point
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def suboff_union() -> object:
    from tensorcad.suboff_appendages import SuboffFullField

    return to_mesh(SuboffFullField(), (109, 13, 19), margin=2).mesh


def test_smoothing_the_suoffer_union_removes_the_slivers(suboff_union: Mesh) -> None:
    """The claim this module exists for, measured on the mesh the sampler produced."""
    before = measure_mesh(suboff_union)
    assert before.sliver_share > 0.01, "the input must still be the sliver-heavy mesh"
    smoothed, moved = taubin_smooth(suboff_union, iterations=20)
    after = measure_mesh(smoothed)
    assert after.watertight and after.edge_manifold
    assert after.triangles == before.triangles
    assert after.sliver_share < 0.001, f"{100 * after.sliver_share:.3f} % slivers left"
    assert after.min_angle_degrees > 10.0 * before.min_angle_degrees
    assert after.max_aspect_ratio < 0.05 * before.max_aspect_ratio
    # ... and it did not do it by moving the surface: the largest displacement is a
    # few cells, and the volume and area are within a few parts in a thousand.
    assert moved < 0.05 * max(suboff_union.bounding_box().sizes)
    assert after.signed_volume == pytest.approx(before.signed_volume, rel=0.005)
    assert after.surface_area == pytest.approx(before.surface_area, rel=0.02)
    # The reference is the analytic quadrature volume of the same model, not the
    # sampler's own output: smoothing must move *towards* it.
    assert abs(after.signed_volume - 0.70754064) < abs(before.signed_volume - 0.70754064)


def test_the_conditioned_union_passes_every_quality_gate(suboff_union: Mesh) -> None:
    conditioned = condition_mesh(suboff_union, smooth_iterations=20)
    checks = {check.name: check for check in inspect_mesh_quality(conditioned.mesh)}
    assert checks["mesh_sliver_share"].passed, checks["mesh_sliver_share"].detail
    assert conditioned.after.watertight
    assert conditioned.after.triangles == conditioned.before.triangles
    edges = len(edge_face_counts(conditioned.mesh))
    euler = (
        len(conditioned.mesh.vertices) - edges + len(conditioned.mesh.faces)
    )
    assert euler == 2
    # The summary is what the handoff record carries.
    summary = conditioned.summary()
    assert summary["after"]["sliver_share"] < 0.001
    assert summary["max_smoothing_move"] > 0.0
    assert summary["collapsed_edges"] == 0
    assert surface_area(conditioned.mesh) > 0.0
