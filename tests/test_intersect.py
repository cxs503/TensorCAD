"""Verification suite for :mod:`tensorcad.intersect`.

The module is the geometric gate behind the PC6 / ice-channel services: a
watertight mesh that folds through itself is rejected on the strength of
:func:`tensorcad.intersect.find_self_intersections`. Until now the primitive was
covered by two hand-picked pairs inside ``tests/test_pc6.py``.

Three independent lines of evidence are used here:

* **analytic cases** -- configurations whose answer is known by construction
  (crossing, disjoint, touching at a vertex or along an edge, coplanar overlap,
  containment, parallel-offset planes, degenerate triangles);
* **differential testing against an exact oracle** -- ``tests/intersect_oracle.py``
  solves the same question with rational arithmetic through a completely
  different route (plane/plane line clipping and 2-D polygon clipping instead of
  edge piercing), so an agreement is not a shared rounding artefact;
* **mesh-level controls** -- a clean closed solid must report nothing, and the
  uniform-hash pruning must return exactly what an all-pairs scan returns.

The suite found one real defect: coplanar pairs were decided by "is a vertex of
one triangle inside the other", which misses overlaps whose corners are all edge
crossings (a Star-of-David configuration, and also two crossing strips). See
``TestCoplanarRegression``.
"""
from __future__ import annotations

import itertools
import random
import sys
from fractions import Fraction
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from intersect_oracle import (  # noqa: E402
    triangles_intersect_clip3d,
    triangles_intersect_exact,
)

from tensorcad.geometry import Mesh  # noqa: E402
from tensorcad.intersect import (  # noqa: E402
    RELATIVE_TOLERANCE,
    find_self_intersections,
    triangles_intersect,
)

#: Lattice with exact binary coordinates: every product and sum the module
#: evaluates on these inputs is exact in double precision, so a disagreement with
#: the rational oracle is a logic error and never a rounding ambiguity.
DYADIC_GRID = (-1.0, -0.5, 0.0, 0.5, 1.0)

FLAT = ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0))
PIERCING = ((0.5, 0.5, -1.0), (0.5, 0.5, 1.0), (1.5, 1.5, 0.5))


def random_triangle(rng: random.Random, grid=DYADIC_GRID):
    return tuple(tuple(rng.choice(grid) for _ in range(3)) for _ in range(3))


class TestAnalyticCases:
    """Configurations with an answer that needs no computation to state."""

    def test_perpendicular_piercing_pair_intersects(self):
        # A vertical blade through a horizontal plate, crossing at (0.5, 0.5, 0).
        assert triangles_intersect(FLAT, PIERCING)

    def test_mirrored_piercing_pair_intersects(self):
        assert triangles_intersect(PIERCING, FLAT)

    def test_far_apart_pair_is_disjoint(self):
        far = ((5.0, 5.0, 0.0), (6.0, 5.0, 0.0), (5.0, 6.0, 0.0))
        assert not triangles_intersect(FLAT, far)
        assert not triangles_intersect(far, FLAT)

    def test_gap_just_above_tolerance_is_disjoint(self):
        # A blade whose lower tip stops 1 mm above the plate: no shared point.
        blade = ((0.5, 0.5, 1e-3), (0.5, 0.5, 2.0), (1.5, 1.5, 1.5))
        assert not triangles_intersect(FLAT, blade)

    def test_touching_at_a_shared_vertex_intersects(self):
        # Closed sets: meeting at a single point counts.
        blade = ((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (1.0, 1.0, 1.0))
        assert triangles_intersect(FLAT, blade)

    def test_touching_along_an_edge_intersects(self):
        # The blade edge (0,0,0)-(2,0,0) is exactly the plate edge, lifted in y.
        blade = ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 1.0, 1.0))
        assert triangles_intersect(FLAT, blade)

    def test_vertex_tip_exactly_on_face_interior_intersects(self):
        tip = ((0.5, 0.5, 0.0), (0.5, 0.5, 1.0), (1.5, 1.5, 1.0))
        assert triangles_intersect(FLAT, tip)

    def test_parallel_offset_planes_are_disjoint(self):
        lifted = tuple((x, y, z + 0.5) for x, y, z in FLAT)
        assert not triangles_intersect(FLAT, lifted)

    def test_degenerate_triangle_carries_no_intersection(self):
        sliver = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0))  # collinear
        assert not triangles_intersect(FLAT, sliver)
        assert not triangles_intersect(sliver, FLAT)

    def test_argument_order_does_not_matter(self):
        rng = random.Random(11)
        for _ in range(300):
            first = random_triangle(rng)
            second = random_triangle(rng)
            assert triangles_intersect(first, second) == triangles_intersect(second, first)

    def test_contract_rejects_non_triangles(self):
        with pytest.raises(ValueError):
            triangles_intersect(((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), FLAT)


class TestCoplanarRegression:
    """Overlaps whose corners are all edge crossings (the defect this suite found).

    Before the fix the coplanar branch only asked whether a vertex of one triangle
    lay inside the other. In both cases below every vertex is outside, so the
    answer was a false negative: two coplanar shells that fold through each other
    passed the gate.
    """

    #: Two coplanar triangles in the plane x + y = 0 -- a plate and a paddle that
    #: cross in a lens whose corners are edge intersections only.
    LENS_A = ((1.0, -1.0, 0.0), (-0.5, 0.5, 0.5), (1.0, -1.0, -0.5))
    LENS_B = ((0.5, -0.5, -0.5), (0.0, 0.0, 0.5), (1.0, -1.0, 0.5))

    def test_coplanar_overlap_without_contained_vertex_intersects(self):
        assert triangles_intersect(self.LENS_A, self.LENS_B)
        assert triangles_intersect(self.LENS_B, self.LENS_A)

    def test_no_vertex_of_either_triangle_is_inside_the_other(self):
        # Documents why the old test could not work: both containment answers False.
        from tensorcad.intersect import _point_in_triangle

        assert not any(_point_in_triangle(p, self.LENS_B, 1e-9) for p in self.LENS_A)
        assert not any(_point_in_triangle(p, self.LENS_A, 1e-9) for p in self.LENS_B)

    def test_star_of_david_pair_intersects(self):
        # Two congruent coplanar triangles rotated by 180 degrees about the
        # centroid of a hexagram: the overlap is a hexagon, no vertex contained.
        up = ((0.0, 1.0, 0.0), (-1.0, -0.5, 0.0), (1.0, -0.5, 0.0))
        down = ((0.0, -1.0, 0.0), (1.0, 0.5, 0.0), (-1.0, 0.5, 0.0))
        assert triangles_intersect(up, down)

    def test_coplanar_disjoint_pair_stays_disjoint(self):
        assert not triangles_intersect(self.LENS_A, ((3.0, -3.0, 0.0), (4.0, -4.0, 0.0), (3.0, -3.0, 1.0)))

    def test_coplanar_containment_intersects(self):
        outer = ((0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (0.0, 10.0, 0.0))
        inner = ((1.0, 1.0, 0.0), (2.0, 1.0, 0.0), (1.0, 2.0, 0.0))
        assert triangles_intersect(outer, inner)
        assert triangles_intersect(inner, outer)

    def test_coplanar_patches_are_reported_at_mesh_level(self):
        # Two coplanar overlapping plates with no shared vertices: the mesh-level
        # gate must report the pair, not return an empty tuple.
        mesh = Mesh(
            vertices=(
                (0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0),
                (0.5, 0.2, 0.0), (2.4, 0.9, 0.0), (1.6, 2.1, 0.0),
            ),
            faces=((0, 1, 2), (3, 4, 5)),
        )
        assert find_self_intersections(mesh) == ((0, 1),)


class TestExactDifferential:
    """Randomised comparison against the rational-arithmetic oracle."""

    def test_dyadic_pairs_match_the_exact_oracle(self):
        rng = random.Random(20261107)
        mismatches = []
        for _ in range(2000):
            first = random_triangle(rng)
            second = random_triangle(rng)
            if triangles_intersect(first, second) != triangles_intersect_exact(first, second):
                mismatches.append((first, second))
        assert not mismatches, f"{len(mismatches)} disagreements, first: {mismatches[0]}"

    def test_three_exact_routes_agree_on_dyadic_pairs(self):
        """The two exact routes must also agree with each other, or neither is a reference."""
        rng = random.Random(4242)
        disagreements = []
        for _ in range(800):
            first = random_triangle(rng)
            second = random_triangle(rng)
            line_route = triangles_intersect_exact(first, second)
            clip_route = triangles_intersect_clip3d(first, second)
            if line_route != clip_route:
                disagreements.append((first, second))
        assert not disagreements, f"{len(disagreements)} route disagreements: {disagreements[0]}"

    def test_thirds_lattice_matches_halfspace_clipping(self):
        """Non-dyadic inputs: the library is compared with 3-D half-space clipping.

        Third-based coordinates are not representable in binary, so the float
        implementation cannot be held to the rational answer bit for bit; the
        honest statement is that it must agree with an exact route on every pair
        whose geometry is unambiguous. The pairs scored here are the ones where
        both exact routes agree with each other.
        """
        grid = (-1.0, -1.0 / 3.0, 0.0, 1.0 / 3.0, 1.0)
        rng = random.Random(7)
        mismatches = []
        for _ in range(600):
            first = random_triangle(rng, grid)
            second = random_triangle(rng, grid)
            line_route = triangles_intersect_exact(first, second)
            clip_route = triangles_intersect_clip3d(first, second)
            if line_route != clip_route:
                continue  # near-degenerate: the exact routes disagree, no verdict
            if triangles_intersect(first, second) != clip_route:
                mismatches.append((first, second))
        assert not mismatches, f"{len(mismatches)} disagreements, first: {mismatches[0]}"


class TestInvariance:
    """The verdict must not depend on where the part sits or how big it is."""

    CASES = (
        (FLAT, PIERCING),
        (TestCoplanarRegression.LENS_A, TestCoplanarRegression.LENS_B),
        (FLAT, ((5.0, 5.0, 0.0), (6.0, 5.0, 0.0), (5.0, 6.0, 0.0))),
    )

    @pytest.mark.parametrize("scale", [2.0**-10, 2.0**-4, 2.0**4, 2.0**10])
    def test_verdict_is_scale_invariant(self, scale):
        for first, second in self.CASES:
            scaled_first = tuple(tuple(c * scale for c in p) for p in first)
            scaled_second = tuple(tuple(c * scale for c in p) for p in second)
            assert triangles_intersect(scaled_first, scaled_second) == triangles_intersect(first, second)

    @pytest.mark.parametrize("offset", [(1000.0, 0.0, 0.0), (-250.0, 125.0, 40.0)])
    def test_verdict_survives_translation(self, offset):
        for first, second in self.CASES:
            moved_first = tuple(tuple(p[i] + offset[i] for i in range(3)) for p in first)
            moved_second = tuple(tuple(p[i] + offset[i] for i in range(3)) for p in second)
            assert triangles_intersect(moved_first, moved_second) == triangles_intersect(first, second)

    def test_axis_permutation_does_not_change_the_verdict(self):
        for trail in itertools.permutations(range(3)):
            permuted_first = tuple(tuple(p[i] for i in trail) for p in FLAT)
            permuted_second = tuple(tuple(p[i] for i in trail) for p in PIERCING)
            assert triangles_intersect(permuted_first, permuted_second)


def _octahedron(scale: float = 1.0) -> Mesh:
    v = (
        (scale, 0.0, 0.0), (-scale, 0.0, 0.0),
        (0.0, scale, 0.0), (0.0, -scale, 0.0),
        (0.0, 0.0, scale), (0.0, 0.0, -scale),
    )
    f = (
        (0, 2, 4), (2, 1, 4), (1, 3, 4), (3, 0, 4),
        (2, 0, 5), (1, 2, 5), (3, 1, 5), (0, 3, 5),
    )
    return Mesh(vertices=v, faces=f)


def _crossing_sheets() -> Mesh:
    """Two triangles that pierce each other, sharing no vertex index."""
    return Mesh(
        vertices=((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (0.0, 2.0, 0.0),
                  (0.5, 0.5, -1.0), (0.5, 0.5, 1.0), (1.5, 1.5, 0.5)),
        faces=((0, 1, 2), (3, 4, 5)),
    )


class TestSelfIntersectionDetection:
    def test_clean_closed_solid_reports_nothing(self):
        assert find_self_intersections(_octahedron()) == ()

    def test_solid_scales_down_to_a_millimetre_part(self):
        assert find_self_intersections(_octahedron(scale=1e-3)) == ()

    def test_crossing_sheets_are_reported_once(self):
        assert find_self_intersections(_crossing_sheets()) == ((0, 1),)

    def test_adjacent_faces_sharing_a_vertex_are_not_reported(self):
        # The octahedron is full of faces that touch at vertices and edges.
        assert find_self_intersections(_octahedron()) == ()

    def test_single_face_mesh_is_empty(self):
        assert find_self_intersections(Mesh(vertices=FLAT, faces=((0, 1, 2),))) == ()

    def test_limit_is_validated(self):
        with pytest.raises(ValueError):
            find_self_intersections(_octahedron(), limit=0)

    @pytest.mark.parametrize("cells", [4, 12, 24, 48])
    def test_pruning_does_not_change_the_answer(self, cells):
        """The uniform hash grid is an accelerator, never a filter."""
        mesh = _crossing_sheets()
        expected = _brute_force_pairs(mesh)
        assert find_self_intersections(mesh, cells_per_axis=cells) == expected

    def test_pruning_matches_brute_force_on_a_warped_grid(self):
        mesh = _warped_grid(3, 3)
        assert find_self_intersections(mesh, cells_per_axis=8) == _brute_force_pairs(mesh)

    def test_limit_truncates_but_stays_sorted(self):
        mesh = _crossing_sheets()
        pairs = find_self_intersections(mesh, limit=1)
        assert pairs == ((0, 1),)


def _brute_force_pairs(mesh: Mesh) -> tuple[tuple[int, int], ...]:
    found = []
    for i, j in itertools.combinations(range(len(mesh.faces)), 2):
        if set(mesh.faces[i]) & set(mesh.faces[j]):
            continue
        if triangles_intersect(
            [mesh.vertices[k] for k in mesh.faces[i]],
            [mesh.vertices[k] for k in mesh.faces[j]],
        ):
            found.append((i, j))
    return tuple(sorted(found))


def _warped_grid(rows: int, columns: int) -> Mesh:
    """A triangular lattice with a deterministic alternation in z.

    Alternating the height of every other interior vertex makes neighbouring
    triangles tilt through one another at some vertices and stay clear at others,
    which gives the pruning test a mixture of intersecting and clean pairs.
    """
    vertices = []
    for r in range(rows + 1):
        for c in range(columns + 1):
            z = 0.0
            if 0 < r < rows and 0 < c < columns and (r + c) % 2 == 0:
                z = 0.6
            vertices.append((float(c), float(r), z))
    stride = columns + 1
    faces = []
    for r in range(rows):
        for c in range(columns):
            a = r * stride + c
            b = a + 1
            d = a + stride
            e = d + 1
            faces.append((a, b, e))
            faces.append((a, e, d))
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))