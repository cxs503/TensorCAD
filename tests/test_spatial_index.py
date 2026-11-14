"""Tests for the uniform-grid spatial index (``tensorcad.spatial_index``).

The contract that matters is *answer equivalence*: the acceleration structure
must return exactly what a brute-force Möller--Trumbore loop returns, including
the misses.  Everything else (memory layout, resolution heuristic, oversized
fallback) is an implementation detail that must not perturb the geometry.
"""
from __future__ import annotations

import math
import random

import pytest

from tensorcad.geometry import Mesh
from tensorcad.mesh import is_watertight, total_signed_volume
from tensorcad.spatial_index import (
    RayHit,
    TriangleGrid,
    brute_force_all,
    brute_force_nearest,
)


# ---------------------------------------------------------------------------
# Fixtures: a closed unit cube and a coarse hull
# ---------------------------------------------------------------------------

def _cube(size: float = 2.0) -> Mesh:
    h = size / 2.0
    vertices = (
        (-h, -h, -h), (h, -h, -h), (h, h, -h), (-h, h, -h),
        (-h, -h, h), (h, -h, h), (h, h, h), (-h, h, h),
    )
    faces = (
        (0, 2, 1), (0, 3, 2),   # -z
        (4, 5, 6), (4, 6, 7),   # +z
        (0, 1, 5), (0, 5, 4),   # -y
        (2, 3, 7), (2, 7, 6),   # +y
        (1, 2, 6), (1, 6, 5),   # +x
        (0, 4, 7), (0, 7, 3),   # -x
    )
    return Mesh(vertices=vertices, faces=faces)


@pytest.fixture(scope="module")
def hull_mesh() -> Mesh:
    from tensorcad.hull import HullParameters, generate_hull

    return generate_hull(
        HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13)
    ).mesh


@pytest.fixture(scope="module")
def cube_mesh() -> Mesh:
    return _cube()


def _random_directions(n: int, seed: int = 11):
    rng = random.Random(seed)
    out = []
    for _ in range(n):
        v = (rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(-1, 1))
        norm = math.sqrt(sum(c * c for c in v))
        if norm < 1e-6:
            continue
        out.append((v[0] / norm, v[1] / norm, v[2] / norm))
    return out


def _assert_same_hit(grid: TriangleGrid, mesh: Mesh, origin, direction) -> None:
    fast = grid.cast_ray(origin, direction)
    slow = brute_force_nearest(mesh, origin, direction)
    if slow is None:
        assert fast is None, f"grid reported {fast} where brute force found no hit"
        return
    assert fast is not None, f"grid missed a hit at t={slow}"
    assert fast.distance == pytest.approx(slow, rel=1e-9, abs=1e-12)


# ---------------------------------------------------------------------------
# Answer equivalence
# ---------------------------------------------------------------------------

def test_matches_brute_force_radial_rays(hull_mesh: Mesh) -> None:
    """Radial rays from the hull origin agree with the O(N) reference."""
    grid = TriangleGrid(hull_mesh)
    from tensorcad.perf_benchmark import _fibonacci_directions

    for direction in _fibonacci_directions(120):
        _assert_same_hit(grid, hull_mesh, (0.0, 0.0, 0.0), direction)


def test_matches_brute_force_arbitrary_rays(hull_mesh: Mesh) -> None:
    """Random origins/directions, including rays that miss entirely."""
    grid = TriangleGrid(hull_mesh)
    rng = random.Random(4242)
    checked_miss = 0
    for direction in _random_directions(90):
        origin = (
            rng.uniform(-14.0, 24.0),
            rng.uniform(-6.0, 6.0),
            rng.uniform(-6.0, 6.0),
        )
        slow = brute_force_nearest(hull_mesh, origin, direction)
        if slow is None:
            checked_miss += 1
        _assert_same_hit(grid, hull_mesh, origin, direction)
    assert checked_miss > 0, "test should exercise misses as well as hits"


def test_matches_brute_force_cube(cube_mesh: Mesh) -> None:
    """Axis-aligned and diagonal rays on a box: exact slab distances."""
    grid = TriangleGrid(cube_mesh, resolution=(4, 4, 4))
    # From outside, +x ray hits the face at x=-1 -> t = 4.
    hit = grid.cast_ray((-5.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert hit is not None
    assert hit.distance == pytest.approx(4.0, abs=1e-12)
    for direction in [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),
                      (1.0, 1.0, 0.0), (1.0, -1.0, 1.0)]:
        _assert_same_hit(grid, cube_mesh, (-3.0, -3.0, -3.0), direction)


def test_all_hits_match_brute_force(hull_mesh: Mesh) -> None:
    """Every crossing along a ray matches the reference, in order."""
    grid = TriangleGrid(hull_mesh)
    for direction in _random_directions(15, seed=7):
        fast = [hit.distance for hit in grid.cast_all((0.0, 0.0, 0.0), direction)]
        slow = list(brute_force_all(hull_mesh, (0.0, 0.0, 0.0), direction))
        # Reference may report the same crossing twice (shared edge); the grid
        # collapses those, so compare the deduplicated sequences.
        deduped: list[float] = []
        for t in slow:
            if deduped and abs(t - deduped[-1]) <= 1e-9 * max(1.0, abs(t)):
                continue
            deduped.append(t)
        assert len(fast) == len(deduped)
        for a, b in zip(fast, deduped):
            assert a == pytest.approx(b, rel=1e-9, abs=1e-12)


# ---------------------------------------------------------------------------
# Point location (inside / outside)
# ---------------------------------------------------------------------------

def test_contains_closed_hull(hull_mesh: Mesh) -> None:
    grid = TriangleGrid(hull_mesh)
    assert grid.contains((5.0, 0.0, -0.4)) is True      # amidships, inside
    assert grid.contains((5.0, 0.0, 1.0)) is False      # above the deck
    assert grid.contains((15.0, 0.0, -0.4)) is False    # ahead of the stem
    assert grid.contains((-2.0, 0.0, -0.4)) is False    # abaft the transom
    assert grid.contains((5.0, 3.0, -0.4)) is False     # outboard


def test_contains_cube(cube_mesh: Mesh) -> None:
    grid = TriangleGrid(cube_mesh)
    assert grid.contains((0.0, 0.0, 0.0)) is True
    assert grid.contains((0.9, 0.9, 0.9)) is True
    assert grid.contains((1.1, 0.0, 0.0)) is False
    assert grid.contains((0.0, 0.0, -2.0)) is False


def test_contains_agrees_with_volume_sign(hull_mesh: Mesh) -> None:
    """A point at the centroid of a closed outward-oriented mesh is inside."""
    volume = total_signed_volume(hull_mesh)
    assert volume > 0.0  # outward orientation
    grid = TriangleGrid(hull_mesh)
    cx = sum(v[0] for v in hull_mesh.vertices) / len(hull_mesh.vertices)
    cz = sum(v[2] for v in hull_mesh.vertices) / len(hull_mesh.vertices)
    assert grid.contains((cx, 0.0, cz * 0.5)) is True


# ---------------------------------------------------------------------------
# Structure: oversized fallback, stats, robustness
# ---------------------------------------------------------------------------

def test_oversized_triangle_fallback_is_exact(hull_mesh: Mesh) -> None:
    """A forced fallback list must not change any answer."""
    normal = TriangleGrid(hull_mesh)
    forced = TriangleGrid(hull_mesh, resolution=(2, 2, 2), oversized_limit=0)
    assert forced.n_oversized == normal.n_triangles  # everything on the fallback
    for direction in _random_directions(40, seed=3):
        origin = (5.0, 0.0, -0.4)
        a = forced.nearest_distance(origin, direction)
        b = normal.nearest_distance(origin, direction)
        if b is None:
            assert a is None
        else:
            assert a == pytest.approx(b, rel=1e-9, abs=1e-12)


def test_stats_are_consistent(hull_mesh: Mesh) -> None:
    grid = TriangleGrid(hull_mesh)
    stats = grid.stats()
    assert stats["triangles"] == len(hull_mesh.faces)
    assert stats["cells"] == grid.n_cells
    assert math.prod(stats["resolution"]) == stats["cells"]
    assert stats["entries"] >= stats["triangles"]
    assert stats["occupied_cells"] <= stats["cells"]
    assert stats["max_triangles_per_cell"] >= 1
    assert stats["index_bytes"] > 0
    assert grid.build_seconds >= 0.0


def test_resolution_override_and_clamping(cube_mesh: Mesh) -> None:
    grid = TriangleGrid(cube_mesh, resolution=(8, 2, 5))
    assert grid.resolution == (8, 2, 5)
    tiny = TriangleGrid(cube_mesh, resolution=(0, -3, 2))
    assert tiny.resolution == (1, 1, 2)


def test_max_distance_culls_far_hits(cube_mesh: Mesh) -> None:
    grid = TriangleGrid(cube_mesh)
    assert grid.cast_ray((-5.0, 0.0, 0.0), (1.0, 0.0, 0.0), max_distance=3.0) is None
    hit = grid.cast_ray((-5.0, 0.0, 0.0), (1.0, 0.0, 0.0), max_distance=4.5)
    assert hit is not None and hit.distance == pytest.approx(4.0, abs=1e-12)


def test_ray_far_outside_grid_misses(cube_mesh: Mesh) -> None:
    grid = TriangleGrid(cube_mesh)
    assert grid.cast_ray((100.0, 100.0, 100.0), (0.0, 0.0, 1.0)) is None
    assert grid.cast_ray((100.0, 0.0, 0.0), (-1.0, 0.0, 0.0)) is not None


def test_zero_direction_raises(cube_mesh: Mesh) -> None:
    grid = TriangleGrid(cube_mesh)
    with pytest.raises(ValueError):
        grid.cast_ray((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    with pytest.raises(ValueError):
        grid.cast_all((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))


def test_empty_mesh_raises() -> None:
    with pytest.raises(ValueError):
        TriangleGrid(Mesh(vertices=(), faces=()))


def test_answers_are_deterministic(hull_mesh: Mesh) -> None:
    a = TriangleGrid(hull_mesh)
    b = TriangleGrid(hull_mesh)
    directions = _random_directions(30, seed=99)
    hits_a = a.query_many((5.0, 0.0, -0.4), directions)
    hits_b = b.query_many((5.0, 0.0, -0.4), directions)
    assert [None if h is None else h.to_dict() for h in hits_a] == [
        None if h is None else h.to_dict() for h in hits_b
    ]


def test_hit_payload(hull_mesh: Mesh) -> None:
    grid = TriangleGrid(hull_mesh)
    hit = grid.cast_ray((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert isinstance(hit, RayHit)
    assert hit.point == pytest.approx((hit.distance, 0.0, 0.0), abs=1e-9)
    assert 0 <= hit.triangle < len(hull_mesh.faces)
    # The reported triangle must actually contain the reported point.
    face = hull_mesh.faces[hit.triangle]
    verts = [hull_mesh.vertices[i] for i in face]
    for axis in range(3):
        lo = min(v[axis] for v in verts)
        hi = max(v[axis] for v in verts)
        assert lo - 1e-9 <= hit.point[axis] <= hi + 1e-9


def test_index_does_not_mutate_mesh(hull_mesh: Mesh) -> None:
    before = (len(hull_mesh.vertices), len(hull_mesh.faces), is_watertight(hull_mesh))
    TriangleGrid(hull_mesh)
    after = (len(hull_mesh.vertices), len(hull_mesh.faces), is_watertight(hull_mesh))
    assert before == after