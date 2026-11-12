"""Tests for the repair path of ``tensorcad.mesh_repair`` (GAP-CAD-009).

The fixtures are analytic: a box (volume ``l*w*h`` exactly, any triangulation), a UV sphere
(volume ``4/3 pi r^3`` analytically, its tessellation error is what the fixture already has)
and the parametric hull (whose analytic volume the hull audit records). Defects are then
*injected* -- facets deleted, a facet wound backwards, a facet duplicated, a sheet turned
inside out -- and the repair has to put the measured numbers back.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.geometry import Mesh
from tensorcad.hull import HullParameters, analytic_volume, generate_hull
from tensorcad.mesh_repair import (
    DEFAULT_MAX_HOLE_EDGES,
    MeshRepairError,
    fill_holes,
    manifoldize,
    measure_mesh,
    orient_consistently,
    repair_mesh,
)

SIZE = (1.0, 1.0, 1.0)


# --------------------------------------------------------------------------- fixtures


def box_mesh(n: int = 2, size=SIZE) -> Mesh:
    """A closed box, every face split into ``n x n`` quads, wound outward."""
    lx, ly, lz = size
    centre = (0.5 * lx, 0.5 * ly, 0.5 * lz)
    vertices: list[tuple[float, float, float]] = []
    lookup: dict[tuple[int, int, int], int] = {}

    def vid(i, j, k):
        if (i, j, k) not in lookup:
            lookup[(i, j, k)] = len(vertices)
            vertices.append((lx * i / n, ly * j / n, lz * k / n))
        return lookup[(i, j, k)]

    faces: list[tuple[int, int, int]] = []
    for axis in range(3):
        u, v = [a for a in range(3) if a != axis]
        for side in (0, 1):
            for a in range(n):
                for b in range(n):
                    quad = []
                    for du, dv in ((0, 0), (1, 0), (1, 1), (0, 1)):
                        index = [0, 0, 0]
                        index[axis] = side * n
                        index[u] = a + du
                        index[v] = b + dv
                        quad.append(vid(*index))
                    for tri in ((quad[0], quad[1], quad[2]), (quad[0], quad[2], quad[3])):
                        faces.append(_outward(vertices, tri, centre))
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


def _outward(vertices, tri, centre) -> tuple[int, int, int]:
    """Wind ``tri`` so that its normal points away from ``centre``."""
    ax, ay, az = vertices[tri[0]]
    bx, by, bz = vertices[tri[1]]
    cx, cy, cz = vertices[tri[2]]
    ux, uy, uz = bx - ax, by - ay, bz - az
    vx, vy, vz = cx - ax, cy - ay, cz - az
    nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
    mx = (ax + bx + cx) / 3.0 - centre[0]
    my = (ay + by + cy) / 3.0 - centre[1]
    mz = (az + bz + cz) / 3.0 - centre[2]
    if nx * mx + ny * my + nz * mz < 0.0:
        return (tri[0], tri[2], tri[1])
    return tri


def sphere_mesh(radius: float = 0.5, rings: int = 12, segments: int = 24) -> Mesh:
    vertices = [(0.0, 0.0, radius)]
    for ring in range(1, rings):
        theta = math.pi * ring / rings
        for segment in range(segments):
            phi = 2.0 * math.pi * segment / segments
            vertices.append(
                (
                    radius * math.sin(theta) * math.cos(phi),
                    radius * math.sin(theta) * math.sin(phi),
                    radius * math.cos(theta),
                )
            )
    vertices.append((0.0, 0.0, -radius))
    south = len(vertices) - 1

    def index(ring, segment):
        return 1 + (ring - 1) * segments + (segment % segments)

    faces: list[tuple[int, int, int]] = []
    centre = (0.0, 0.0, 0.0)
    for segment in range(segments):
        faces.append(_outward(vertices, (0, index(1, segment), index(1, segment + 1)), centre))
    for ring in range(1, rings - 1):
        for segment in range(segments):
            a = index(ring, segment)
            b = index(ring, segment + 1)
            c = index(ring + 1, segment + 1)
            d = index(ring + 1, segment)
            faces.append(_outward(vertices, (a, b, c), centre))
            faces.append(_outward(vertices, (a, c, d), centre))
    for segment in range(segments):
        faces.append(
            _outward(vertices, (south, index(rings - 1, segment + 1), index(rings - 1, segment)), centre)
        )
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


def drop_faces(mesh: Mesh, indices) -> Mesh:
    doomed = set(indices)
    return Mesh(
        vertices=mesh.vertices,
        faces=tuple(face for index, face in enumerate(mesh.faces) if index not in doomed),
    )


# --------------------------------------------------------------------------- clean mesh


def test_a_clean_box_is_left_alone():
    box = box_mesh(2)
    metrics = measure_mesh(box)
    assert metrics.closed and metrics.manifold and metrics.consistently_wound
    assert metrics.signed_volume == pytest.approx(1.0)
    assert metrics.surface_area == pytest.approx(6.0)

    result = repair_mesh(box)
    assert result.accepted
    assert result.mesh.faces == box.faces
    assert "already closed" in result.reason
    assert all(step.changed == 0 for step in result.steps)


# --------------------------------------------------------------------------- hole filling


def test_filling_a_triangle_hole_restores_the_volume_exactly():
    box = box_mesh(2)
    holed = drop_faces(box, [0])
    before = measure_mesh(holed)
    assert before.open_edges == 3 and not before.closed

    fill = fill_holes(holed)
    assert fill.filled == 1
    hole = dict(fill.holes[0])
    assert hole["loop_vertices"] == 3 and hole["triangles"] == 1
    assert hole["planar_deviation_m"] == pytest.approx(0.0, abs=1e-15)

    after = measure_mesh(fill.mesh)
    assert after.closed and after.manifold and after.consistently_wound
    assert after.faces == len(box.faces)
    # the fill lies in the face's own plane, so the volume is the analytic 1 m3
    assert after.signed_volume == pytest.approx(1.0, abs=1e-12)


def test_filling_a_quad_hole_uses_two_facets_and_keeps_the_volume():
    box = box_mesh(2)
    holed = drop_faces(box, [0, 1])  # the two triangles of one quad
    assert measure_mesh(holed).open_edges == 4
    fill = fill_holes(holed)
    assert fill.filled == 1
    assert fill.holes[0]["loop_vertices"] == 4
    assert fill.holes[0]["triangles"] == 2
    assert measure_mesh(fill.mesh).signed_volume == pytest.approx(1.0, abs=1e-12)


def test_two_separate_holes_are_both_closed():
    box = box_mesh(3)
    holed = drop_faces(box, [0, 54])  # two facets on different faces of the box
    fill = fill_holes(holed)
    assert fill.filled == 2
    after = measure_mesh(fill.mesh)
    assert after.closed and after.signed_volume == pytest.approx(1.0, abs=1e-12)


def test_a_hole_in_a_curved_surface_is_closed_and_the_flat_fill_is_measured():
    sphere = sphere_mesh()
    analytic = 4.0 / 3.0 * math.pi * 0.5**3
    tessellated = measure_mesh(sphere).signed_volume
    holed = drop_faces(sphere, [3, 4])
    fill = fill_holes(holed)
    assert fill.filled == 1
    # the patch is planar; the deviation from the sphere is reported, not zero
    assert fill.holes[0]["planar_deviation_m"] > 0.0
    after = measure_mesh(fill.mesh)
    assert after.closed
    # closing the hole by a chord patch can only remove volume, and by less than the
    # tessellation error the fixture already carries
    assert after.signed_volume <= tessellated + 1e-12
    assert abs(after.signed_volume - tessellated) / tessellated < 1e-3
    assert abs(after.signed_volume - analytic) / analytic < 0.03  # fixture tessellation error is 2.8 %


def test_a_hull_with_a_hole_is_repaired_within_the_surface_error():
    parameters = HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13)
    hull = generate_hull(parameters).mesh
    reference = analytic_volume(10.0, 2.0, 0.8)
    holed = drop_faces(hull, [7, 8, 9])
    result = repair_mesh(holed, reference_volume=reference)
    assert result.accepted, result.reason
    assert result.after.closed and result.after.manifold
    assert result.after.self_intersections == 0
    assert result.volume_relative_change is not None
    assert result.volume_relative_change < 0.01
    # the repaired body keeps the surface discretisation error the hull audit records
    assert abs(result.after.signed_volume - reference) / reference < 0.01


# --------------------------------------------------------------------------- winding


def test_a_single_backwards_facet_is_turned_round():
    box = box_mesh(2)
    broken = Mesh(
        vertices=box.vertices,
        faces=tuple([(box.faces[0][0], box.faces[0][2], box.faces[0][1])] + list(box.faces[1:])),
    )
    assert measure_mesh(broken).winding_conflicts >= 1
    fixed, step = orient_consistently(broken)
    assert step.changed == 1
    metrics = measure_mesh(fixed)
    assert metrics.winding_conflicts == 0
    assert metrics.signed_volume == pytest.approx(1.0, abs=1e-12)


def test_an_inside_out_shell_is_turned_outward():
    box = box_mesh(2)
    inward = Mesh(vertices=box.vertices, faces=tuple((f[0], f[2], f[1]) for f in box.faces))
    assert measure_mesh(inward).signed_volume == pytest.approx(-1.0)
    fixed, _ = orient_consistently(inward)
    assert measure_mesh(fixed).signed_volume == pytest.approx(1.0, abs=1e-12)


def test_orientation_repair_is_idempotent():
    box = box_mesh(2)
    inward = Mesh(vertices=box.vertices, faces=tuple((f[0], f[2], f[1]) for f in box.faces))
    once, _ = orient_consistently(inward)
    twice, step = orient_consistently(once)
    assert step.changed == 0
    assert twice.faces == once.faces


# --------------------------------------------------------------------------- non-manifold


def test_a_duplicated_facet_is_dropped_and_the_volume_survives():
    box = box_mesh(2)
    duplicated = Mesh(vertices=box.vertices, faces=box.faces + (box.faces[0],))
    assert measure_mesh(duplicated).nonmanifold_edges == 3
    result = repair_mesh(duplicated, reference_volume=1.0)
    assert result.accepted, result.reason
    assert result.after.faces == len(box.faces)
    assert result.after.nonmanifold_edges == 0
    assert result.after.signed_volume == pytest.approx(1.0, abs=1e-12)


def _two_boxes(shift) -> Mesh:
    """Two unit boxes, the second translated by ``shift`` so they touch on a corner/edge."""
    from tensorcad import mesh as meshlib

    box = box_mesh(1)
    moved = tuple((x + shift[0], y + shift[1], z + shift[2]) for x, y, z in box.vertices)
    both = Mesh(
        vertices=box.vertices + moved,
        faces=box.faces + tuple(tuple(v + len(box.vertices) for v in face) for face in box.faces),
    )
    return meshlib.weld(both, tolerance=1e-9)


def test_two_boxes_joined_at_one_edge_are_split_into_two_manifold_sheets():
    welded = _two_boxes((1.0, 1.0, 0.0))  # shares the edge x=1, y=1, z in [0, 1]
    assert measure_mesh(welded).nonmanifold_edges > 0
    split, step = manifoldize(welded)
    assert step.changed > 0
    healed = measure_mesh(split)
    assert healed.nonmanifold_edges == 0
    assert healed.closed
    assert healed.components == 2
    assert healed.signed_volume == pytest.approx(2.0, abs=1e-9)


def test_a_bow_tie_vertex_is_split_without_moving_geometry():
    welded = _two_boxes((1.0, 1.0, 1.0))  # shares exactly one corner
    split, step = manifoldize(welded)
    healed = measure_mesh(split)
    assert healed.nonmanifold_edges == 0
    assert healed.closed
    assert healed.components == 2
    assert healed.signed_volume == pytest.approx(2.0, abs=1e-9)


# --------------------------------------------------------------------------- refusals


def test_a_self_intersecting_mesh_is_refused_rather_than_patched():
    sphere = sphere_mesh()
    vertices = list(sphere.vertices)
    # pull the north pole through the south pole: opposite facets now pierce each other
    vertices[0] = (0.0, 0.0, -0.9)
    broken = Mesh(vertices=tuple(vertices), faces=sphere.faces)
    assert measure_mesh(broken).self_intersections > 0
    result = repair_mesh(broken)
    assert not result.accepted
    assert result.mesh.faces == broken.faces
    assert "self-intersecting" in result.reason


def test_a_hole_larger_than_the_cap_is_refused_not_filled():
    box = box_mesh(6)
    holed = drop_faces(box, range(12))  # a loop well above a 6-edge cap
    result = repair_mesh(holed, reference_volume=1.0, max_hole_edges=6)
    assert not result.accepted
    assert "could not be filled" in result.reason or "above the 6 cap" in result.reason


def test_the_default_cap_is_recorded_and_finite():
    assert DEFAULT_MAX_HOLE_EDGES >= 3


def test_a_repair_that_would_move_the_volume_beyond_tolerance_is_refused():
    sphere = sphere_mesh()
    holed = drop_faces(sphere, list(range(0, 24)))  # the whole north cap
    exact = 4.0 / 3.0 * math.pi * 0.5**3
    result = repair_mesh(holed, reference_volume=exact, volume_tolerance=1e-9)
    assert not result.accepted
    assert "volume moved" in result.reason


def test_a_mesh_with_no_faces_is_refused():
    empty = Mesh(vertices=((0.0, 0.0, 0.0),), faces=())
    with pytest.raises(MeshRepairError):
        fill_holes(empty)


def test_negative_tolerances_are_refused():
    with pytest.raises(MeshRepairError):
        repair_mesh(box_mesh(1), volume_tolerance=-0.1)


# --------------------------------------------------------------------------- reporting


def test_the_result_serialises_to_json():
    import json

    result = repair_mesh(drop_faces(box_mesh(2), [0]), reference_volume=1.0)
    payload = json.dumps(result.to_dict(), allow_nan=False)
    decoded = json.loads(payload)
    assert decoded["accepted"] is True
    assert decoded["after"]["closed"] is True
    assert decoded["holes"]["filled"] == 1
    assert decoded["volume_relative_change"] < 1e-9


def test_a_clean_repair_says_so_in_its_notes():
    result = repair_mesh(box_mesh(2), reference_area=6.0)
    assert any("reference volume" in note for note in result.notes)
    assert result.area_relative_change == pytest.approx(0.0, abs=1e-12)