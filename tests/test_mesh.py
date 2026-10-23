"""Tests for the generic mesh algorithms and the geometry contracts."""
from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import BoundingBox, Mesh, MeshBuilder

# A unit cube, outward-oriented, split into triangles. Faces 4..5 are the y=0
# side, 6..7 the y=1 side, and so on.
CUBE_VERTICES = (
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0),
    (1.0, 0.0, 1.0),
    (1.0, 1.0, 1.0),
    (0.0, 1.0, 1.0),
)
CUBE_FACES = (
    (0, 3, 2), (0, 2, 1),
    (4, 5, 6), (4, 6, 7),
    (0, 1, 5), (0, 5, 4),
    (2, 3, 7), (2, 7, 6),
    (1, 2, 6), (1, 6, 5),
    (0, 4, 7), (0, 7, 3),
)


def cube(*, invert: tuple[int, ...] = (), mirror_from: int | None = None) -> Mesh:
    faces = list(CUBE_FACES)
    if mirror_from is not None:
        # Emulate the old bug: a second sheet built by negating y with the same
        # winding, which inverts every normal on that half.
        faces = faces[:mirror_from] + [tuple(reversed(face)) for face in faces[mirror_from:]]
    faces = [tuple(reversed(face)) if index in invert else face for index, face in enumerate(faces)]
    return Mesh(vertices=CUBE_VERTICES, faces=tuple(faces))  # type: ignore[arg-type]


def test_clean_cube_is_closed_and_manifold() -> None:
    mesh = cube()
    assert meshlib.is_watertight(mesh)
    assert meshlib.is_edge_manifold(mesh)
    assert not meshlib.boundary_edges(mesh)
    assert meshlib.total_signed_volume(mesh) == pytest.approx(1.0)
    assert meshlib.surface_area(mesh) == pytest.approx(6.0)
    assert len(meshlib.edge_face_counts(mesh)) == 18  # 12 edges + 6 face diagonals


def test_inverted_faces_are_detected_and_repaired() -> None:
    broken = cube(invert=(5, 9))
    assert meshlib.total_signed_volume(broken) == pytest.approx(2.0 / 3.0)
    result = meshlib.orient_outward(broken)
    assert result.flipped_faces == 2
    assert meshlib.total_signed_volume(result.mesh) == pytest.approx(1.0)
    assert meshlib.is_watertight(result.mesh)


def test_mirrored_patch_is_reoriented_without_changing_topology() -> None:
    broken = cube(mirror_from=6)
    assert abs(meshlib.total_signed_volume(broken)) < 1.0
    result = meshlib.orient_outward(broken)
    assert meshlib.total_signed_volume(result.mesh) == pytest.approx(1.0)
    assert result.components == 1
    assert len(result.mesh.faces) == len(broken.faces)


def test_opening_a_face_creates_boundary_edges() -> None:
    opened = Mesh(vertices=CUBE_VERTICES, faces=CUBE_FACES[2:])  # drop both top triangles
    assert len(meshlib.boundary_edges(opened)) == 4  # the shared diagonal is simply gone
    assert not meshlib.is_watertight(opened)


def test_degenerate_and_zero_area_faces_are_found_and_droppable() -> None:
    sliver = Mesh(
        vertices=((0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 0, 0)),
        faces=((0, 1, 2), (0, 2, 3), (1, 1, 2)),  # type: ignore[arg-type]
    )
    assert meshlib.zero_area_faces(sliver) == (2,)  # the face that repeats a vertex
    cleaned = meshlib.remove_degenerate_faces(sliver)
    assert len(cleaned.faces) == 2
    assert len(cleaned.vertices) == 4  # coincident-but-distinct points stay distinct
    with pytest.raises(ValueError, match="out-of-range"):
        Mesh(vertices=((0, 0, 0),), faces=((0, 1, 2),))  # type: ignore[arg-type]


def test_zero_area_face_area_reports() -> None:
    collapsed = Mesh(vertices=((0, 0, 0), (1, 0, 0), (1, 0, 0)), faces=((0, 1, 2),))  # type: ignore[arg-type]
    assert meshlib.triangle_areas(collapsed) == (0.0,)
    assert meshlib.zero_area_faces(collapsed) == (0,)


def test_weld_merges_coincident_imported_vertices() -> None:
    doubled = Mesh(
        vertices=((0, 0, 0), (1, 0, 0), (0, 1, 0), (-0.0, 0.0, -0.0), (1.0, 0.0, 1e-12)),
        faces=((0, 1, 2), (3, 4, 2)),  # type: ignore[arg-type]
    )
    welded = meshlib.weld(doubled, tolerance=1e-9)
    assert len(welded.vertices) == 3
    assert welded.faces == ((0, 1, 2), (0, 1, 2))
    with pytest.raises(ValueError, match="non-negative"):
        meshlib.weld(doubled, tolerance=-1.0)


def test_fully_inverted_cube_is_flipped_as_one_component() -> None:
    inverted = cube(invert=tuple(range(12)))
    assert meshlib.total_signed_volume(inverted) == pytest.approx(-1.0)
    result = meshlib.orient_outward(inverted)
    assert result.flipped_faces == 12
    assert result.components == 1
    assert meshlib.total_signed_volume(result.mesh) == pytest.approx(1.0)


def test_orienting_an_already_correct_mesh_is_a_no_op() -> None:
    result = meshlib.orient_outward(cube())
    assert result.flipped_faces == 0
    assert result.mesh.faces == CUBE_FACES


def test_degenerate_face_does_not_break_the_orientation_walk() -> None:
    weird = Mesh(vertices=((0, 0, 0), (1, 0, 0), (1, 1, 0)), faces=((0, 1, 1), (0, 1, 2)))  # type: ignore[arg-type]
    result = meshlib.orient_outward(weird)
    assert len(result.mesh.faces) == 2
    assert not meshlib.is_watertight(weird), "an empty or open mesh is never called closed"
    assert not meshlib.is_watertight(Mesh(vertices=((0, 0, 0),), faces=()))


def test_weld_uses_neighbouring_cells_and_tolerates_zero_tolerance() -> None:
    # 0.9 and 1.1 are within tolerance 1.0 but land in different cells.
    straddling = Mesh(
        vertices=((0.9, 0.0, 0.0), (1.1, 0.0, 0.0), (0.0, 5.0, 0.0)),
        faces=((0, 1, 2),),  # type: ignore[arg-type]
    )
    assert len(meshlib.weld(straddling, tolerance=1.0).vertices) == 2
    assert len(meshlib.weld(straddling, tolerance=0.0).vertices) == 3, "exact matching only"
    with pytest.raises(ValueError, match="finite"):
        meshlib.weld(Mesh(vertices=((float("inf"), 0, 0),), faces=()))


def test_remove_degenerate_faces_is_a_no_op_when_there_is_nothing_to_do() -> None:
    mesh = cube()
    assert meshlib.remove_degenerate_faces(mesh) is mesh, "the same object, not a copy"


def test_mesh_accessors_and_validation() -> None:
    mesh = cube()
    assert len(mesh) == 12
    assert next(iter(mesh.iter_vertices())) == CUBE_VERTICES[0]
    with pytest.raises(ValueError, match="not a triangle"):
        Mesh(vertices=CUBE_VERTICES, faces=((0, 1),))  # type: ignore[arg-type]


def test_face_normal_and_centroid_follow_winding() -> None:
    mesh = cube()
    top = meshlib.face_normal(mesh, CUBE_FACES[2])
    assert top[2] > 0, "the +z face must have a positive z normal"
    assert meshlib.face_centroid(mesh, (0, 1, 2)) == pytest.approx((2.0 / 3.0, 1.0 / 3.0, 0.0))


class TestMeshBuilder:
    def test_dedups_by_exact_position_including_signed_zero(self) -> None:
        builder = MeshBuilder()
        first = builder.add_vertex((0.0, 0.0, 0.0))
        assert builder.add_vertex((-0.0, 0.0, -0.0)) == first
        assert builder.add_vertex((0.0, 0.0, 1e-16)) != first  # no tolerance guessing
        assert builder.vertex_count == 2

    def test_quad_with_a_pinched_corner_becomes_one_triangle(self) -> None:
        builder = MeshBuilder()
        a = builder.add_vertex((0, 0, 0))
        b = builder.add_vertex((1, 0, 0))
        c = builder.add_vertex((1, 1, 0))
        assert builder.add_quad(a, b, c, a) == 1
        assert builder.add_quad(a, b, c, d=a) == 1
        assert len(builder.build().faces) == 2

    def test_quad_with_a_collapsed_diagonal_is_dropped(self) -> None:
        builder = MeshBuilder()
        a, b, c = builder.add_vertex((0, 0, 0)), builder.add_vertex((1, 0, 0)), builder.add_vertex((0, 1, 0))
        assert builder.add_quad(a, b, a, c) == 0
        assert builder.degenerate_skipped == 1
        assert builder.build().faces == ()

    def test_rejects_non_finite_and_misshapen_points(self) -> None:
        builder = MeshBuilder()
        with pytest.raises(ValueError, match="not finite"):
            builder.add_vertex((math.nan, 0.0, 0.0))
        with pytest.raises(ValueError, match="exactly 3 coordinates"):
            builder.add_vertex((0.0, 0.0))

    def test_mesh_carries_units_and_frame_into_every_artifact(self) -> None:
        builder = MeshBuilder(units="mm", coordinate_frame="x=fore, z=up")
        builder.add_vertex((0, 0, 0))
        builder.add_vertex((1, 0, 0))
        builder.add_vertex((0, 1, 0))
        mesh = builder.build()
        assert mesh.units == "mm"
        assert mesh.coordinate_frame == "x=fore, z=up"
        assert mesh.metadata()["bounding_box"]["size"] == (1.0, 1.0, 0.0)

    def test_obj_face_indices_are_one_based(self) -> None:
        mesh = cube()
        assert next(iter(mesh.to_obj_faces())) == (1, 4, 3)


def test_bounding_box_and_metadata() -> None:
    box = BoundingBox((-1.0, -2.0, -3.0), (1.0, 2.0, 3.0))
    assert box.sizes == (2.0, 4.0, 6.0)
    assert box.to_dict()["size"] == (2.0, 4.0, 6.0)
    with pytest.raises(ValueError, match="empty mesh"):
        Mesh(vertices=(), faces=()).bounding_box()
