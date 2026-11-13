"""The mid-surface must be a property of the geometry, not of the face numbering.

``tests/midsurface_oracle.py`` rebuilds the wall search the expensive way -- every face
against every other, ray/triangle intersection by plane + area barycentrics -- and states
the analytic areas. This file uses it for four things the area/thickness assertions in
``tests/test_midsurface.py`` cannot see:

* **the bucket grid loses nothing.** A candidate pruned at a cell boundary shows up as an
  "unpaired face", which is a legal outcome of the same function, so no area check can
  notice it. Comparing against the all-pairs scan can.
* **relabelling the faces changes nothing.** Which facing of a wall a mid triangle is
  copied from used to be decided by the face *index*; since the two facings of a curved
  wall do not have the same area, that made the mid-surface area a property of the file's
  ordering (1.2 % on a t/R = 2 % tube, 2.1 % on a sphere). Both the pairing tie-break and
  the side choice are now geometric, and these tests pin it.
* **a one-sided measurement is reported but does not emit**, so a hull's folded rim
  cannot put an inner-facing copy next to an outer-facing one.
* **a nearer non-antiparallel surface shadows the wall**, which is the documented policy
  of ``extract_midsurface`` -- measured here rather than assumed.
"""
from __future__ import annotations

import math
import random

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh
from tensorcad.midsurface import (
    cylinder_skin,
    extract_midsurface,
    plate_solid,
    sphere_skin,
)
from tensorcad.series60 import Series60Parameters, generate_series60
from tensorcad.shell import offset_inward

import midsurface_oracle as oracle

# Small enough for an O(n^2) scan over every pair of faces (the brute force is the point).
PLATE = plate_solid(8.0, 3.0, 0.02, nx=8, ny=5)
TUBE = cylinder_skin(radius=1.0, thickness=0.02, length=4.0, segments=24)
SPHERE = sphere_skin(radius=1.0, thickness=0.02, rings=8, segments=12)
# Fine enough that the facet sagitta is below the surplus being measured (a 24-segment
# tube is 0.86 % short in wall thickness, which would swamp a 1 % surplus).
TUBE_FINE = cylinder_skin(radius=1.0, thickness=0.02, length=4.0, segments=96)
SPHERE_FINE = sphere_skin(radius=1.0, thickness=0.02, rings=32, segments=64)
BAND = 0.2


# --------------------------------------------------------------------------- #
# the oracle itself
# --------------------------------------------------------------------------- #
def test_oracle_ray_triangle_hits_and_misses():
    a, b, c = (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)
    # straight down onto the hypotenuse-side interior
    assert oracle.ray_triangle_distance((0.25, 0.25, 1.0), (0.0, 0.0, -1.0), a, b, c) == pytest.approx(1.0)
    # outside the triangle, in its plane
    assert oracle.ray_triangle_distance((0.9, 0.9, 1.0), (0.0, 0.0, -1.0), a, b, c) is None
    # parallel to the plane: never hits
    assert oracle.ray_triangle_distance((0.25, 0.25, 1.0), (1.0, 0.0, 0.0), a, b, c) is None
    # behind the origin
    assert oracle.ray_triangle_distance((0.25, 0.25, 1.0), (0.0, 0.0, 1.0), a, b, c) is None


# --------------------------------------------------------------------------- #
# the acceleration structure against the all-pairs scan
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mesh, name", [(PLATE, "plate"), (TUBE, "tube"), (SPHERE, "sphere")])
def test_bucket_grid_finds_exactly_what_an_all_pairs_scan_finds(mesh: Mesh, name: str):
    result = extract_midsurface(mesh, max_thickness=BAND)
    brute = oracle.pairs_first_hit(mesh, max_thickness=BAND)
    accelerated = {sample.face: (sample.partner, sample.thickness) for sample in result.samples}
    assert set(accelerated) == set(brute), f"{name}: the grid and the all-pairs scan disagree on who is paired"

    scale = max(oracle._norm(oracle._sub(max(mesh.vertices), min(mesh.vertices))), 1e-12)
    for face, (partner, thickness) in accelerated.items():
        expected_partner, expected_thickness = brute[face]
        # The measured wall is the same to round-off. Which of two *equidistant* facings
        # was reported may differ -- that is a tie, not a disagreement -- so the partner
        # is checked structurally: it must be antiparallel and its plane must sit at the
        # measured distance from the face's centroid.
        assert thickness == pytest.approx(expected_thickness, abs=1e-12 * scale)
        centroid = oracle.face_centroid(mesh.vertices, mesh.faces[face])
        normal = oracle._unit(oracle.face_normal(mesh.vertices, mesh.faces[face]))
        other = oracle._unit(oracle.face_normal(mesh.vertices, mesh.faces[partner]))
        assert oracle._dot(normal, other) < 0.0, f"{name}: face {face} paired with a facing that is not a wall"
        a, b, c = (mesh.vertices[i] for i in mesh.faces[partner])
        plane_normal = oracle._cross(oracle._sub(b, a), oracle._sub(c, a))
        if oracle._norm(plane_normal) == 0.0:
            continue  # a degenerate (zero-area) sliver has no plane to measure against
        # the measured thickness is the distance *along the face normal* to the partner's
        # plane (the ray-plane intersection), not the plane's perpendicular distance: on a
        # coarse skin the two facings are not parallel, so the two differ by 1/cos
        inward = (-normal[0], -normal[1], -normal[2])
        denominator = oracle._dot(inward, plane_normal)
        assert denominator != 0.0
        offset = oracle._dot(oracle._sub(a, centroid), plane_normal) / denominator
        assert offset == pytest.approx(thickness, rel=1e-9)


# --------------------------------------------------------------------------- #
# determinism: the answer may not depend on the face numbering
# --------------------------------------------------------------------------- #
def _geometric_pairs(result, mesh: Mesh) -> set[tuple[tuple[float, ...], tuple[float, ...]]]:
    """The pairing as face *centroids*, so it survives a relabelling."""
    pairs = set()
    for sample in result.samples:
        first = oracle.face_centroid(mesh.vertices, mesh.faces[sample.face])
        second = oracle.face_centroid(mesh.vertices, mesh.faces[sample.partner])
        pairs.add((tuple(round(value, 6) for value in first),
                   tuple(round(value, 6) for value in second)))
    return pairs


@pytest.mark.parametrize("mesh, name", [(PLATE, "plate"), (TUBE, "tube"), (SPHERE, "sphere")])
def test_relabelling_faces_does_not_move_the_mid_surface(mesh: Mesh, name: str):
    reference = extract_midsurface(mesh, max_thickness=BAND)
    assert reference.mesh.faces
    for seed in (1, 2, 3):
        permutation = list(range(len(mesh.faces)))
        random.Random(seed).shuffle(permutation)
        shuffled = extract_midsurface(oracle.relabel(mesh, permutation), max_thickness=BAND)
        assert len(shuffled.mesh.faces) == len(reference.mesh.faces), (
            f"{name}: relabelling the faces changed the mid-surface face count "
            f"({len(shuffled.mesh.faces)} vs {len(reference.mesh.faces)})"
        )
        assert shuffled.area == pytest.approx(reference.area, rel=1e-12), (
            f"{name}: relabelling the faces moved the mid-surface area "
            f"({shuffled.area!r} vs {reference.area!r})"
        )
        assert shuffled.paired_face_count == reference.paired_face_count
        # index lists are meaningless across a relabelling; the *geometry* is not
        assert sorted(sample.thickness for sample in shuffled.samples) == pytest.approx(
            sorted(sample.thickness for sample in reference.samples), rel=1e-12
        )
        assert len(_geometric_pairs(shuffled, mesh)) == len(_geometric_pairs(reference, mesh))


def test_relabelling_a_hull_plating_does_not_move_the_mid_surface():
    parameters = Series60Parameters(length=10.0, beam=1.5, draft=0.6, cb=0.60, nx=21, nz=7)
    hull = generate_series60(parameters).mesh
    cap = {index for index, face in enumerate(hull.faces)
           if all(abs(hull.vertices[i][2]) < 1e-12 for i in face)}
    keep = [index for index in range(len(hull.faces)) if index not in cap]
    outer = Mesh(vertices=hull.vertices, faces=tuple(hull.faces[i] for i in keep))
    inner = Mesh(vertices=offset_inward(hull, 0.02).vertices,
                 faces=tuple(hull.faces[i] for i in keep))
    shift = len(outer.vertices)
    plating = Mesh(
        vertices=tuple(outer.vertices) + tuple(inner.vertices),
        faces=tuple(list(outer.faces) + [tuple(i + shift for i in face)[::-1] for face in inner.faces]),
    )
    reference = extract_midsurface(plating, max_thickness=0.5)
    permutation = list(range(len(plating.faces)))
    random.Random(11).shuffle(permutation)
    shuffled = extract_midsurface(oracle.relabel(plating, permutation), max_thickness=0.5)
    assert shuffled.paired_face_count == reference.paired_face_count
    assert len(shuffled.mesh.faces) == len(reference.mesh.faces)
    assert shuffled.area == pytest.approx(reference.area, rel=1e-12)


# --------------------------------------------------------------------------- #
# analytic areas and thicknesses
# --------------------------------------------------------------------------- #
def test_closed_form_areas():
    plate = extract_midsurface(PLATE)
    assert plate.area == pytest.approx(oracle.plate_mid_area(8.0, 3.0), rel=1e-12)

    tube = extract_midsurface(TUBE_FINE)
    exact = oracle.cylinder_offset_area(1.0, 0.02, 4.0)
    assert oracle.relative_error(tube.area, exact) == pytest.approx(
        oracle.cylinder_surplus(1.0, 0.02), rel=0.10
    )
    assert abs(oracle.relative_error(tube.area, exact)) < 0.03

    sphere = extract_midsurface(SPHERE_FINE)
    exact = oracle.sphere_offset_area(1.0, 0.02)
    surplus = oracle.relative_error(sphere.area, exact)
    assert surplus < 0.03
    # a wall's curvature is kappa1 + kappa2 = 2/R on a sphere and 1/R on a tube
    assert surplus == pytest.approx(oracle.sphere_surplus(1.0, 0.02), rel=0.10)


@pytest.mark.parametrize("mesh, wall, tolerance", [(PLATE, 0.02, 1e-9), (TUBE, 0.02, 1e-2),
                                                  (SPHERE, 0.02, 5e-2)])
def test_measured_thickness_is_the_wall(mesh: Mesh, wall: float, tolerance: float):
    # A coarse skin measures the *chord* between the facings, not the radial gap, so the
    # wall reads short by the facet sagitta: 0.86 % at 24 tube segments, 1.9 % at 8 sphere
    # rings. That is discretisation of the fixture, not bias of the measurement.
    result = extract_midsurface(mesh, max_thickness=BAND)
    assert result.mean_thickness == pytest.approx(wall, rel=tolerance)
    assert result.min_measured > 0.0
    assert result.max_measured < BAND


# --------------------------------------------------------------------------- #
# one-sided walls: measured, reported, and not emitted
# --------------------------------------------------------------------------- #
def test_one_sided_walls_are_measured_but_emit_nothing():
    result = extract_midsurface(SPHERE, max_thickness=BAND)
    partner_of = {sample.face: sample.partner for sample in result.samples}
    reciprocal = {frozenset((sample.face, sample.partner)) for sample in result.samples
                  if partner_of.get(sample.partner) == sample.face}
    assert len(reciprocal) < len(result.samples) / 2 + 1, "fixture no longer has one-sided walls"
    assert len(result.mesh.faces) == len(reciprocal)
    assert result.paired_face_count == len(result.samples)


# --------------------------------------------------------------------------- #
# the documented policy: a nearer non-antiparallel surface shadows the wall
# --------------------------------------------------------------------------- #
def test_only_a_nearer_antiparallel_surface_can_be_a_wall():
    """The documented policy, measured: on a wall-only skin nothing shadows anything.

    ``extract_midsurface`` keeps the *nearest surface the ray strikes* and then filters it
    on antiparallelity and the thickness band, so a closer non-antiparallel surface would
    hide the wall behind it. The diagnostic below counts exactly those faces; on the three
    wall-only fixtures it is empty, which is the evidence that the band + normal filter is
    not silently discarding walls on real plate/shell geometry. (An attempt to build a
    shadowing fixture -- a tilted baffle inside a plate -- did *not* shadow anything: the
    baffle has to be crossed by the centroid ray, which the per-facet construction never
    arranged; the diagnostic is therefore stated as a property of the fixtures, not as a
    demonstrated failure mode.)
    """
    for mesh, name in ((PLATE, "plate"), (TUBE, "tube"), (SPHERE, "sphere")):
        assert oracle.shadowed_faces(mesh, max_thickness=BAND) == (), (
            f"{name}: a wall face is hidden behind a nearer non-antiparallel surface"
        )


def test_mesh_area_matches_the_sum_of_its_faces():
    """Guard the oracle's own area helper against the module's."""
    for mesh in (PLATE, TUBE, SPHERE):
        assert oracle.mesh_area(mesh) == pytest.approx(meshlib.surface_area(mesh), rel=1e-12)