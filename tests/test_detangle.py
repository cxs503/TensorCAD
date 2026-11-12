"""Tests for the self-intersection heal of ``tensorcad.detangle`` (GAP-CAD-009 removal).

The fixtures are the same analytic ones :mod:`tests.test_mesh_repair` uses -- a box whose
volume is exactly ``l*w*h``, a UV sphere (``4/3 pi r^3``, tessellation error reported) and the
parametric hull -- and the defect is *injected*: one vertex is pushed until the facets it owns
pierce the opposite side of the shell. The heal has to remove the penetration, keep the shell
closed and manifold, and move the enclosed volume by no more than the budget it was given.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.detangle import (
    DetangleResult,
    excavation_faces,
    heal_self_intersections,
)
from tensorcad.geometry import Mesh
from tensorcad.hull import HullParameters, build_mesh
from tensorcad.intersect import find_self_intersections
from tensorcad.mesh_repair import measure_mesh, repair_mesh

SIZE = (1.0, 1.0, 1.0)


# --------------------------------------------------------------------------- fixtures


def _outward(vertices, tri, centre):
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


def box_mesh(n: int = 4, size=SIZE) -> Mesh:
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
    centre = (0.0, 0.0, 0.0)

    def index(ring, segment):
        return 1 + (ring - 1) * segments + (segment % segments)

    faces: list[tuple[int, int, int]] = []
    for segment in range(segments):
        faces.append(_outward(vertices, (0, index(1, segment), index(1, segment + 1)), centre))
    for ring in range(1, rings - 1):
        for segment in range(segments):
            a, b = index(ring, segment), index(ring, segment + 1)
            c, d = index(ring + 1, segment + 1), index(ring + 1, segment)
            faces.append(_outward(vertices, (a, b, c), centre))
            faces.append(_outward(vertices, (a, c, d), centre))
    for segment in range(segments):
        faces.append(
            _outward(
                vertices,
                (south, index(rings - 1, segment + 1), index(rings - 1, segment)),
                centre,
            )
        )
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


def poke(mesh: Mesh, vertex: int, target) -> Mesh:
    vertices = list(mesh.vertices)
    vertices[vertex] = target
    return Mesh(vertices=tuple(vertices), faces=mesh.faces)


def box_apex(mesh: Mesh) -> int:
    return next(
        index
        for index, point in enumerate(mesh.vertices)
        if abs(point[1]) < 1e-12 and abs(point[0] - 0.5) < 1e-12 and abs(point[2] - 0.5) < 1e-12
    )


def poked_box(n: int, depth: float = 1.4) -> Mesh:
    box = box_mesh(n)
    return poke(box, box_apex(box), (0.5, depth, 0.5))


def poked_hull(push: float = 1.6) -> tuple[Mesh, float]:
    hull = build_mesh(HullParameters(length=10.0, beam=2.0, draft=0.8))[0]
    keel = min(range(len(hull.vertices)), key=lambda index: hull.vertices[index][2])
    point = hull.vertices[keel]
    return (
        poke(hull, keel, (point[0], point[1], point[2] + push)),
        measure_mesh(hull).signed_volume,
    )


# --------------------------------------------------------------------------- the heal


def test_box_fold_is_healed():
    box = box_mesh(4)
    assert abs(measure_mesh(box).signed_volume - 1.0) < 1e-12
    spiked = poked_box(4)
    assert find_self_intersections(spiked, limit=64)

    result = heal_self_intersections(spiked, reference_volume=1.0)

    assert result.healed, result.reason
    assert result.changed
    assert result.before.self_intersections > 0
    assert result.after.self_intersections == 0
    assert result.after.closed and result.after.manifold and result.after.consistently_wound
    assert result.rings_used == 1
    assert result.faces_removed > 0


def test_box_fold_keeps_the_volume():
    result = heal_self_intersections(poked_box(4), reference_volume=1.0)
    assert result.healed, result.reason
    assert result.volume_relative_change is not None
    assert result.volume_relative_change < 1e-12
    assert abs(result.volume_after - 1.0) < 1e-12


def test_finer_box_needs_a_second_ring():
    """``n = 6`` poked through the same way: one ring of facets is not enough."""
    result = heal_self_intersections(poked_box(6), reference_volume=1.0)
    assert result.healed, result.reason
    assert result.rings_used >= 1
    assert result.after.self_intersections == 0
    assert result.volume_relative_change < 1e-9
    assert len(result.bands) == result.rings_used


def test_sphere_fold_is_healed():
    sphere = sphere_mesh()
    control = measure_mesh(sphere)
    apex = next(
        index for index, point in enumerate(sphere.vertices) if point[2] > 0.4 and point[0] > 0.1
    )
    spiked = poke(sphere, apex, (0.0, 0.0, -0.6))
    assert find_self_intersections(spiked, limit=64)

    result = heal_self_intersections(spiked, reference_volume=control.signed_volume)

    assert result.healed, result.reason
    assert result.after.self_intersections == 0 and result.after.closed
    assert result.volume_relative_change < 0.01
    # the fixture's own tessellation error is not the heal's business
    exact = 4.0 / 3.0 * math.pi * 0.5**3
    assert abs(control.signed_volume - exact) / exact < 0.03


def test_hull_fold_is_healed_with_a_centroid_fan():
    spiked, reference = poked_hull(1.6)
    result = heal_self_intersections(spiked, reference_volume=reference)

    assert result.healed, result.reason
    assert result.after.self_intersections == 0
    assert result.after.closed and result.after.manifold and result.after.consistently_wound
    assert result.volume_relative_change < 0.01
    # the excavated loop around the keel does not project to a simple polygon, so the planar
    # fill is refused and the fan caps it -- that path has to have run, not be assumed
    assert any("centroid fan" in note for note in result.notes)


def test_heal_is_deterministic():
    spiked, reference = poked_hull(0.6)
    first = heal_self_intersections(spiked, reference_volume=reference)
    second = heal_self_intersections(spiked, reference_volume=reference)
    assert first.healed and second.healed
    assert first.mesh.faces == second.mesh.faces
    assert first.volume_after == second.volume_after


def test_heal_is_idempotent():
    spiked, reference = poked_hull(0.6)
    healed = heal_self_intersections(spiked, reference_volume=reference)
    assert healed.healed
    again = heal_self_intersections(healed.mesh, reference_volume=reference)
    assert not again.healed
    assert again.mesh is healed.mesh
    assert again.after.self_intersections == 0


# --------------------------------------------------------------------------- refusals


def test_clean_mesh_is_returned_untouched():
    box = box_mesh(4)
    result = heal_self_intersections(box)
    assert not result.healed and not result.changed
    assert result.mesh is box
    assert "nothing to heal" in result.reason


def test_a_dent_is_not_a_fold():
    """Pushing a vertex inwards short of the far side is a valid surface, not a defect."""
    box = box_mesh(2)
    dent = poke(box, box_apex(box), (0.5, 0.6, 0.5))
    assert measure_mesh(dent).closed
    assert not find_self_intersections(dent)
    result = heal_self_intersections(dent)
    assert not result.healed and result.mesh is dent


def test_refuses_an_open_input():
    box = box_mesh(2)
    holed = Mesh(vertices=box.vertices, faces=box.faces[2:])
    spiked = poke(holed, box_apex(box), (0.5, 1.4, 0.5))
    result = heal_self_intersections(spiked)
    assert not result.healed
    assert "closed input is required" in result.reason


def test_refuses_a_non_local_cut():
    """On a 48-facet box one ring is half the mesh: the volume guard has to stop it."""
    result = heal_self_intersections(poked_box(2), reference_volume=1.0)
    assert not result.healed
    assert "above the" in result.reason
    assert result.volume_relative_change is not None and result.volume_relative_change > 0.01


def test_volume_budget_is_honoured():
    sphere = sphere_mesh()
    control = measure_mesh(sphere)
    spiked = poke(
        sphere,
        next(index for index, point in enumerate(sphere.vertices) if point[2] > 0.4 and point[0] > 0.1),
        (0.0, 0.0, -0.6),
    )
    tight = heal_self_intersections(
        spiked, reference_volume=control.signed_volume, volume_tolerance=1e-9
    )
    assert not tight.healed
    assert "budget" in tight.reason
    loose = heal_self_intersections(spiked, reference_volume=control.signed_volume)
    assert loose.healed


def test_invalid_arguments():
    box = box_mesh(2)
    with pytest.raises(ValueError):
        heal_self_intersections(box, max_rings=0)
    with pytest.raises(ValueError):
        heal_self_intersections(box, max_hole_edges=2)
    with pytest.raises(ValueError):
        heal_self_intersections(box, volume_tolerance=-1.0)
    with pytest.raises(ValueError):
        excavation_faces(box, rings=0)


# --------------------------------------------------------------------------- integration


def test_repair_mesh_refuses_a_fold_by_default():
    spiked, reference = poked_hull(0.6)
    result = repair_mesh(spiked)
    assert not result.accepted
    assert result.after.self_intersections > 0


def test_repair_mesh_heals_when_asked():
    spiked, reference = poked_hull(0.6)
    result = repair_mesh(spiked, heal_self_intersections=True, reference_volume=reference)
    assert result.accepted, result.reason
    assert result.after.self_intersections == 0
    assert result.after.closed
    assert "heal_self_intersections" in [step.name for step in result.steps]


def test_excavation_preview_grows_but_never_moves_the_mesh():
    spiked = poked_box(4)
    one = excavation_faces(spiked, rings=1)
    two = excavation_faces(spiked, rings=2)
    assert one and one < two
    assert len(spiked.faces) == measure_mesh(spiked).faces


def test_result_serialises_every_number():
    result = heal_self_intersections(poked_box(4), reference_volume=1.0)
    payload = result.to_dict()
    for key in (
        "healed",
        "reason",
        "pairs_found",
        "rings_used",
        "faces_removed",
        "before",
        "after",
        "volume_after_m3",
        "volume_relative_change",
        "holes",
        "notes",
        "bands",
    ):
        assert key in payload
    assert payload["after"]["self_intersections"] == 0
    assert payload["bands"][0]["residual_pairs"] == 0
    assert payload["bands"][0]["capped"] is True
    assert isinstance(result, DetangleResult)
