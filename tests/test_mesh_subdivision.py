"""Uniform 1-to-4 subdivision: the invariants that make it a valid large-mesh source.

The point of these tests is not "the code runs" but "refining must not change the
solid": every assertion below is a property an imported-mesh substitute has to
keep, otherwise the 10^7-triangle benchmark would measure a *different* body
than the 10^6 rung and the ladder would be meaningless.
"""
from __future__ import annotations

import pytest

from tensorcad.geometry import Mesh
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.mesh import (
    edge_face_counts,
    face_normal,
    is_edge_manifold,
    is_watertight,
    surface_area,
    total_signed_volume,
    triangle_area,
)
from tensorcad.mesh_subdivision import (
    expected_face_count,
    subdivide_mesh,
    subdivision_edges,
)


def _tetrahedron() -> Mesh:
    vertices = (
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
    )
    faces = ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3))
    return Mesh(vertices=vertices, faces=faces)


def _octahedron() -> Mesh:
    vertices = (
        (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0), (0.0, -1.0, 0.0),
        (0.0, 0.0, 1.0), (0.0, 0.0, -1.0),
    )
    faces = (
        (0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
        (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5),
    )
    return Mesh(vertices=vertices, faces=faces)


def test_face_count_is_four_to_the_levels():
    mesh = _octahedron()
    assert expected_face_count(len(mesh.faces), 0) == 8
    assert expected_face_count(len(mesh.faces), 1) == 32
    assert expected_face_count(len(mesh.faces), 3) == 8 * 64
    for levels in (0, 1, 2):
        assert len(subdivide_mesh(mesh, levels).faces) == expected_face_count(8, levels)


def test_vertex_count_matches_v_plus_e():
    """Midpoints are welded per undirected edge, so V' = V + E exactly."""
    mesh = _octahedron()
    edges = len(edge_face_counts(mesh))
    refined = subdivide_mesh(mesh)
    assert len(refined.vertices) == len(mesh.vertices) + edges
    assert len(subdivision_edges(mesh)) == edges


def test_shared_edge_gets_one_vertex():
    """Two faces sharing an edge must reference the *same* midpoint index."""
    mesh = _tetrahedron()
    refined = subdivide_mesh(mesh)
    # Every refined vertex must be used by at least one face (no orphan midpoint).
    used = {v for face in refined.faces for v in face}
    assert used == set(range(len(refined.vertices)))
    # And no two distinct midpoint vertices may coincide in space.
    assert len(set(refined.vertices)) == len(refined.vertices)


def test_closed_manifold_euler_preserved():
    for mesh in (_tetrahedron(), _octahedron(), generate_hull(HullParameters(nx=41, nz=13))):
        refined = subdivide_mesh(mesh, 2)
        assert is_watertight(refined)
        assert is_edge_manifold(refined)
        euler = len(refined.vertices) - len(edge_face_counts(refined)) + len(refined.faces)
        base_euler = len(mesh.vertices) - len(edge_face_counts(mesh)) + len(mesh.faces)
        assert euler == base_euler == 2


def test_volume_and_area_are_preserved():
    """Each planar parent is tiled exactly, so the integrals must not move."""
    mesh = generate_hull(HullParameters(nx=41, nz=13))
    refined = subdivide_mesh(mesh, 2)
    v0, v1 = total_signed_volume(mesh), total_signed_volume(refined)
    a0, a1 = surface_area(mesh), surface_area(refined)
    assert abs(v1 - v0) / abs(v0) < 1e-12
    assert abs(a1 - a0) / a0 < 1e-12


def test_children_of_one_parent_tile_it():
    mesh = generate_hull(HullParameters(nx=21, nz=9))
    refined = subdivide_mesh(mesh, 1)
    # Face order is parent-major, four children per parent.
    for index in range(0, len(refined.faces), 4):
        parent = mesh.faces[index // 4]
        parent_area = triangle_area(mesh, parent)
        child_area = sum(triangle_area(refined, f) for f in refined.faces[index:index + 4])
        assert abs(child_area - parent_area) <= 1e-12 * max(parent_area, 1.0)


def test_child_winding_matches_parent():
    mesh = generate_hull(HullParameters(nx=21, nz=9))
    refined = subdivide_mesh(mesh, 1)
    for index in range(0, len(refined.faces), 4):
        nx, ny, nz = face_normal(mesh, mesh.faces[index // 4])
        for child in refined.faces[index:index + 4]:
            cx, cy, cz = face_normal(refined, child)
            dot = cx * nx + cy * ny + cz * nz
            assert dot > 0.0


def test_input_mesh_is_not_mutated():
    mesh = _octahedron()
    before_vertices, before_faces = mesh.vertices, mesh.faces
    subdivide_mesh(mesh, 2)
    assert mesh.vertices is before_vertices
    assert mesh.faces is before_faces
    assert len(mesh.faces) == 8


def test_levels_zero_is_a_copy():
    mesh = _tetrahedron()
    copy = subdivide_mesh(mesh, 0)
    assert copy.faces == mesh.faces
    assert copy.vertices == mesh.vertices
    assert copy is not mesh


def test_negative_levels_rejected():
    with pytest.raises(ValueError):
        subdivide_mesh(_tetrahedron(), -1)
    with pytest.raises(ValueError):
        expected_face_count(10, -1)


def test_empty_mesh_rejected():
    with pytest.raises(ValueError):
        subdivide_mesh(Mesh(vertices=(), faces=()), 1)