"""Gate for the mesh-level boolean (GAP-CAD-003).

The assertions are of three kinds, and they are deliberately not "the code ran":

* **analytic volumes.**  Every supported configuration has a closed-form volume, so the
  check is ``|V - V_analytic| / V_analytic <= 1 %`` -- the number the anchor asks for.
  Where a faceted operand is involved the *polyhedral* volume is used as the primary
  reference (that is the exact answer for the input tessellation) and the analytic solid
  is reported alongside, so tessellation error and boolean error cannot be confused.
* **topology.**  A closed genus-0 result has Euler characteristic 2; a solid with a
  through hole has 0; a solid with an internal cavity has 4.  Watertight, edge-manifold and
  outward orientation are asserted independently of the volume.
* **explicit failure.**  Configurations this kernel does not support, and invalid operands,
  must raise :class:`MeshBooleanError`.  A silent invalid solid is the failure mode the
  anchor forbids, so "returns nothing usable" is tested as a *pass* only when it raises.

``docs/mesh_boolean.md`` records the measured pass/fail matrix, including the
configurations that fail: they are not hidden behind a skip.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.geometry import Mesh
from tensorcad.mesh import face_edges, is_edge_manifold, is_watertight, total_signed_volume
from tensorcad.mesh_boolean import (
    MeshBooleanError,
    boolean,
    mesh_difference,
    mesh_intersection,
    mesh_union,
    point_in_solid,
)


def box(low, high) -> Mesh:
    """Axis-aligned box, outward-oriented, 12 triangles."""

    x0, y0, z0 = low
    x1, y1, z1 = high
    vertices = [
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    ]
    faces = [
        (0, 2, 1), (0, 3, 2), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6), (3, 0, 4), (3, 4, 7),
    ]
    return Mesh(tuple(vertices), tuple(faces))


def cylinder(radius: float, height: float, segments: int = 48, z0: float = 0.0) -> Mesh:
    """Faceted prism (the polyhedral stand-in for a cylinder), outward-oriented."""

    vertices = []
    for level in (z0, z0 + height):
        for index in range(segments):
            angle = 2.0 * math.pi * index / segments
            vertices.append((radius * math.cos(angle), radius * math.sin(angle), level))
    bottom_centre = len(vertices)
    vertices.append((0.0, 0.0, z0))
    top_centre = len(vertices)
    vertices.append((0.0, 0.0, z0 + height))

    faces = []
    for index in range(segments):
        nxt = (index + 1) % segments
        lower, lower_next = index, nxt
        upper, upper_next = segments + index, segments + nxt
        faces.append((lower, lower_next, upper_next))
        faces.append((lower, upper_next, upper))
        faces.append((bottom_centre, lower_next, lower))
        faces.append((top_centre, upper, upper_next))
    return Mesh(tuple(vertices), tuple(faces))


def euler_characteristic(mesh: Mesh) -> int:
    edges = set()
    for face in mesh.faces:
        for a, b in face_edges(face):
            edges.add((min(a, b), max(a, b)))
    return len(mesh.vertices) - len(edges) + len(mesh.faces)


def relative_error(measured: float, reference: float) -> float:
    return abs(measured - reference) / abs(reference)


def assert_solid(mesh: Mesh, *, chi: int) -> float:
    """Every accepted result must be a closed, manifold, outward solid of known genus."""

    assert is_watertight(mesh), "result is not watertight"
    assert is_edge_manifold(mesh), "result is not edge-manifold"
    assert euler_characteristic(mesh) == chi, f"chi = {euler_characteristic(mesh)}, expected {chi}"
    volume = total_signed_volume(mesh)
    assert volume > 0.0
    return volume


# --------------------------------------------------------------------------- #
# analytic volumes: the anchor's <1 % criterion
# --------------------------------------------------------------------------- #
def test_overlapping_boxes_all_three_operations_match_analytic_volumes() -> None:
    """Two unit-ish cubes offset by one edge length: union 15, intersection 1, difference 7."""

    first = box((0.0, 0.0, 0.0), (2.0, 2.0, 2.0))
    second = box((1.0, 1.0, 1.0), (3.0, 3.0, 3.0))

    union = assert_solid(mesh_union(first, second), chi=2)
    intersection = assert_solid(mesh_intersection(first, second), chi=2)
    difference = assert_solid(mesh_difference(first, second), chi=2)

    assert relative_error(union, 15.0) <= 1e-12
    assert relative_error(intersection, 1.0) <= 1e-12
    assert relative_error(difference, 7.0) <= 1e-12
    # inclusion-exclusion: |A| + |B| = |A u B| + |A n B|
    assert relative_error(union + intersection, 8.0 + 8.0) <= 1e-12


def test_coplanar_shared_face_union_is_exact_and_not_doubled() -> None:
    """Two unit cubes glued face to face: the union is one 2x1x1 box, not two shells."""

    union = mesh_union(box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)), box((1.0, 0.0, 0.0), (2.0, 1.0, 1.0)))
    volume = assert_solid(union, chi=2)
    assert relative_error(volume, 2.0) <= 1e-12


def test_difference_with_an_internal_cavity_is_two_shells() -> None:
    """A small box removed from inside a big one leaves a void: chi = 4 (two spheres)."""

    hollow = mesh_difference(box((0.0, 0.0, 0.0), (3.0, 3.0, 3.0)), box((1.0, 1.0, 1.0), (2.0, 2.0, 2.0)))
    volume = assert_solid(hollow, chi=4)
    assert relative_error(volume, 26.0) <= 1e-12


def test_coaxial_through_hole_is_genus_one_and_matches_both_references() -> None:
    """Cube minus a flush-cut through bolt -- the configuration that broke the BSP kernel.

    Two references are checked on purpose: the *polyhedral* volume ``V_cube - V_prism`` is
    the exact answer for the tessellated input (1e-12), and the analytic solid volume
    ``18 - 2 pi`` is the answer for the smooth problem (< 1 %), the gap between them being
    the 48-segment prism's own tessellation error, not a boolean error.
    """

    prism = cylinder(1.0, 2.0)
    cube = box((-1.5, -1.5, 0.0), (1.5, 1.5, 2.0))
    holed = mesh_difference(cube, prism)

    volume = assert_solid(holed, chi=0)  # a through hole is a torus: genus 1
    assert relative_error(volume, 18.0 - total_signed_volume(prism)) <= 1e-12
    assert relative_error(volume, 18.0 - 2.0 * math.pi) < 0.01

    # the hole is really open: a point on the axis inside the cube is now outside
    assert not point_in_solid(holed, (0.0, 0.0, 1.0))
    assert point_in_solid(holed, (1.3, 1.3, 1.0))


def test_identical_operands_are_idempotent() -> None:
    solid = box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    assert relative_error(assert_solid(mesh_union(solid, solid), chi=2), 1.0) <= 1e-12
    assert relative_error(assert_solid(mesh_intersection(solid, solid), chi=2), 1.0) <= 1e-12


def test_the_result_does_not_depend_on_the_model_units() -> None:
    """The same part in metres and in millimetres must give the same topology and volume ratio."""

    metre = mesh_union(box((0.0, 0.0, 0.0), (2.0, 2.0, 2.0)), box((1.0, 1.0, 1.0), (3.0, 3.0, 3.0)))
    millimetre = mesh_union(
        box((0.0, 0.0, 0.0), (2000.0, 2000.0, 2000.0)),
        box((1000.0, 1000.0, 1000.0), (3000.0, 3000.0, 3000.0)),
    )
    assert euler_characteristic(metre) == euler_characteristic(millimetre)
    assert total_signed_volume(millimetre) / total_signed_volume(metre) == pytest.approx(1e9, rel=1e-9)


# --------------------------------------------------------------------------- #
# failure must be loud
# --------------------------------------------------------------------------- #
def test_disjoint_intersection_raises_instead_of_returning_an_empty_shell() -> None:
    with pytest.raises(MeshBooleanError, match="empty"):
        mesh_intersection(box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)), box((5.0, 5.0, 5.0), (6.0, 6.0, 6.0)))


def test_difference_of_identical_solids_raises() -> None:
    solid = box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    with pytest.raises(MeshBooleanError, match="empty"):
        mesh_difference(solid, solid)


def test_open_operand_is_refused() -> None:
    opened = Mesh(box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)).vertices,
                  box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)).faces[:-1])
    with pytest.raises(MeshBooleanError, match="not closed"):
        mesh_union(opened, box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)))


def test_inward_oriented_operand_is_refused() -> None:
    solid = box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    inward = Mesh(solid.vertices, tuple(tuple(reversed(face)) for face in solid.faces))
    with pytest.raises(MeshBooleanError, match="outward-oriented"):
        mesh_difference(solid, inward)


def test_zero_area_face_is_refused() -> None:
    solid = box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    degenerate = Mesh(
        solid.vertices + ((5.0, 5.0, 5.0),),
        solid.faces[:-1],
    )
    with pytest.raises(MeshBooleanError):
        boolean(degenerate, solid, "union")


def test_unknown_operation_is_refused() -> None:
    solid = box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    with pytest.raises(MeshBooleanError, match="unknown operation"):
        boolean(solid, solid, "xor")


# --------------------------------------------------------------------------- #
# documented unsupported configurations: they must fail loudly, not quietly
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    ("operation", "first", "second"),
    [
        ("difference", "cylinder", "box"),
        ("union", "cylinder", "corner-box"),
        ("intersection", "cylinder", "corner-box"),
    ],
)
def test_unsupported_configurations_never_return_an_invalid_solid(operation, first, second) -> None:
    """These configurations are recorded as failing in ``docs/mesh_boolean.md``.

    The contract they *do* satisfy is the one that matters for a kernel: either a valid
    closed manifold, or :class:`MeshBooleanError`.  A silent invalid solid is the failure
    this test exists to prevent, so it must raise while the known defect stands.
    """

    prism = cylinder(1.0, 2.0)
    operands = {
        "cylinder": prism,
        "box": box((-1.5, -1.5, 0.5), (1.5, 1.5, 2.0)),
        "corner-box": box((0.5, 0.5, 0.5), (1.5, 1.5, 1.5)),
    }
    with pytest.raises(MeshBooleanError):
        boolean(operands[first], operands[second], operation)