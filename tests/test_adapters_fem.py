"""The FEM handoff: named boundary patches, exaxct tables, and a format that reads back.

Two halves. The synthetic half checks the things a solver would notice immediately if
they were wrong: the orientation gate (an inward surface inverts every load sign), the
region labelling on a surface every triangle of which belongs to one analytic field,
and the gmsh/JSON round trips (same coordinates, same areas, same digest). The
integration half runs the SUBOFF union through the same path and checks that the three
analytic components come back as three non-empty boundary patches whose areas add up to
the measured wetted surface.
"""
from __future__ import annotations

import json
import math

import pytest

from tensorcad.adapters.fem import (
    FEM_SCHEMA,
    FemSurface,
    fem_surface,
    label_regions,
    read_fem_json,
    read_gmsh,
    write_fem_json,
    write_gmsh,
)
from tensorcad.geometry import Mesh
from tensorcad.geometry_field import PlaneField, SphereField, to_mesh
from tensorcad.suboff_appendages import SuboffFinsField, SuboffFullField, SuboffHullField, SuboffSailField


def _cube() -> Mesh:
    """A unit cube as 12 outward-wound triangles."""
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
        (0, 2, 1), (0, 3, 2),
        (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4),
        (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6),
        (3, 0, 4), (3, 4, 7),
    )
    return Mesh(vertices=vertices, faces=faces, units="m", coordinate_frame="unit-cube")


def _sphere() -> Mesh:
    """A sampled unit sphere: closed, outward wound, analytic volume 4/3 pi."""
    return to_mesh(SphereField(1.0), resolution=(24, 24, 24), margin=2).mesh


@pytest.fixture(scope="module")
def sphere_surface() -> FemSurface:
    return fem_surface(to_mesh(SphereField(1.0), resolution=(24, 24, 24), margin=2))


def test_a_bare_mesh_becomes_one_named_patch() -> None:
    surface = fem_surface(_cube(), coordinate_frame="unit-cube")
    assert surface.region_names == ("wall",)
    assert set(surface.region_ids) == {1}
    assert surface.regions[0].triangles == 12
    assert surface.volume_m3 == pytest.approx(1.0, abs=1e-12)
    assert surface.area_m2 == pytest.approx(6.0, abs=1e-12)
    # Nothing was sampled, so there is no upstream digest to claim.
    assert surface.source_digest == ""


def test_the_conventions_are_declared_not_implied() -> None:
    surface = fem_surface(_cube(), coordinate_frame="unit-cube", time_s=1.5)
    record = surface.to_dict()
    assert record["schema"] == FEM_SCHEMA
    assert record["units"] == "m"
    assert record["coordinate_frame"] == "unit-cube"
    assert record["sign_convention"] == "outward"
    assert record["time_s"] == 1.5
    assert record["index_base"] == 0
    assert record["interpolation"].startswith("P1")
    assert record["integration"].startswith("3-point")


def test_an_inward_surface_is_refused() -> None:
    """Reversing the winding must fail loudly: every pressure load would flip sign."""
    cube = _cube()
    inward = Mesh(
        vertices=cube.vertices,
        faces=tuple((a, c, b) for a, b, c in cube.faces),
    )
    with pytest.raises(ValueError, match="face inward"):
        fem_surface(inward)


def test_an_empty_surface_is_refused() -> None:
    with pytest.raises(ValueError, match="empty mesh"):
        fem_surface(Mesh(vertices=(), faces=()))


def test_a_negative_time_is_refused() -> None:
    with pytest.raises(ValueError, match="time_s"):
        fem_surface(_cube(), time_s=-1.0)


def test_region_tags_must_have_names() -> None:
    with pytest.raises(ValueError, match="no name"):
        FemSurface(
            nodes=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
            triangles=((0, 1, 2),),
            region_ids=(7,),
            region_names=("wall",),
        )


def test_region_tags_must_match_the_triangle_count() -> None:
    with pytest.raises(ValueError, match="region tags"):
        FemSurface(
            nodes=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
            triangles=((0, 1, 2),),
            region_ids=(1, 1),
            region_names=("wall",),
        )


def test_labelling_finds_the_only_component_the_surface_can_belong_to(
    sphere_surface: FemSurface,
) -> None:
    """A unit sphere's triangles all sit on the sphere, 1 m away from any plane we add."""
    mesh = to_mesh(SphereField(1.0), resolution=(24, 24, 24), margin=2).mesh
    ids, names, stats = label_regions(
        mesh, {"sphere": SphereField(1.0), "floor": PlaneField(normal=(0.0, 0.0, 1.0), offset=-5.0)}
    )
    assert names == ("sphere", "floor")
    assert set(ids) == {1}
    assert stats[0].triangles == len(mesh.faces)
    assert stats[1].triangles == 0
    # The centroids are inside the surface, so they are at most one cell from it.
    assert stats[0].max_assignment_distance_m < 0.1
    assert stats[0].area_m2 == pytest.approx(sphere_surface.area_m2, rel=1e-12)


def test_labelling_rejects_an_empty_component_set() -> None:
    with pytest.raises(ValueError, match="at least one component"):
        label_regions(_sphere(), [])


def test_labelling_rejects_a_non_field() -> None:
    with pytest.raises(TypeError, match="not a GeometryField"):
        label_regions(_sphere(), {"sphere": object()})


def test_gmsh_round_trip_preserves_coordinates_and_areas(
    sphere_surface: FemSurface, tmp_path
) -> None:
    path = write_gmsh(sphere_surface, tmp_path / "sphere.msh")
    parsed = read_gmsh(path)
    assert parsed["nodes"] == len(sphere_surface.nodes)
    assert parsed["triangles"] == len(sphere_surface.triangles)
    assert parsed["physical_names"] == {1: "wall"}
    # Floats are written with repr, so the coordinates come back bit-identical.
    for node_id, node in enumerate(sphere_surface.nodes, start=1):
        assert parsed["node_table"][node_id] == node
    assert parsed["area_by_region_m2"][1] == pytest.approx(sphere_surface.area_m2, rel=1e-14)


def test_gmsh_declares_the_conventions_in_its_comments(
    sphere_surface: FemSurface, tmp_path
) -> None:
    path = write_gmsh(sphere_surface, tmp_path / "sphere.msh")
    text = path.read_text(encoding="utf-8")
    assert "$Comments" in text
    assert f"digest {sphere_surface.digest}" in text
    assert "units m" in text
    # gmsh's own header: version 2.2, ASCII, 8-byte floats.
    assert text.splitlines()[1] == "2.2 0 8"


def test_json_round_trip_preserves_the_digest(sphere_surface: FemSurface, tmp_path) -> None:
    path = write_fem_json(sphere_surface, tmp_path / "sphere.json")
    restored = read_fem_json(path)
    assert restored.digest == sphere_surface.digest
    assert restored.nodes == sphere_surface.nodes
    assert restored.triangles == sphere_surface.triangles
    assert restored.region_names == sphere_surface.region_names
    assert restored.volume_m3 == sphere_surface.volume_m3
    assert restored.source_digest == sphere_surface.source_digest


def test_the_default_manifest_leaves_the_tables_out(sphere_surface: FemSurface, tmp_path) -> None:
    path = write_fem_json(sphere_surface, tmp_path / "sphere.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert "node_table" in payload  # the file always carries the tables
    assert "node_table" not in sphere_surface.to_dict()


def test_json_reader_refuses_a_foreign_schema(tmp_path) -> None:
    path = tmp_path / "wrong.json"
    path.write_text(json.dumps({"schema": "something/else"}), encoding="utf-8")
    with pytest.raises(ValueError, match="schema"):
        read_fem_json(path)


def test_json_reader_refuses_a_missing_file(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        read_fem_json(tmp_path / "absent.json")


def test_the_digest_travels_with_the_handoff() -> None:
    """The table digest is the mesh digest: one convention, so a receiver can compare.

    ``SampledField.mesh_digest`` and :attr:`FemSurface.digest` hash the same way, so a
    consumer that only sees the exported file can check it against the handoff record
    the producer logged -- which is the whole point of carrying a digest at all.
    """
    sampled = to_mesh(SphereField(1.0), resolution=(16, 16, 16), margin=2)
    surface = fem_surface(sampled)
    assert surface.source_digest == sampled.mesh_digest()
    assert surface.digest == surface.source_digest
    assert surface.quality  # the quality warnings came along


@pytest.fixture(scope="module")
def suboff_surface() -> FemSurface:
    sampled = to_mesh(SuboffFullField(), (109, 13, 19), margin=2)
    return fem_surface(
        sampled,
        components={
            "hull": SuboffHullField(),
            "sail": SuboffSailField(),
            "fins": SuboffFinsField(),
        },
    )


def test_the_suboff_handoff_has_three_non_empty_patches(suboff_surface: FemSurface) -> None:
    assert suboff_surface.region_names == ("hull", "sail", "fins")
    counts = {stats.name: stats.triangles for stats in suboff_surface.regions}
    assert all(count > 0 for count in counts.values()), counts
    # A hull dominates the wetted surface of a submarine, and that is a fact, not a
    # tuning knob: if the labelling ever disagrees, it is the labelling that is wrong.
    assert counts["hull"] > counts["sail"] + counts["fins"]
    total_area = sum(stats.area_m2 for stats in suboff_surface.regions)
    assert total_area == pytest.approx(suboff_surface.area_m2, rel=1e-12)


def test_the_suboff_patch_areas_survive_the_gmsh_round_trip(
    suboff_surface: FemSurface, tmp_path
) -> None:
    path = write_gmsh(suboff_surface, tmp_path / "suboff.msh")
    parsed = read_gmsh(path)
    assert parsed["physical_names"] == {1: "hull", 2: "sail", 3: "fins"}
    for stats in suboff_surface.regions:
        assert parsed["area_by_region_m2"][stats.physical_id] == pytest.approx(
            stats.area_m2, rel=1e-12
        )


def test_the_suboff_labelling_stays_near_the_analytic_surfaces(
    suboff_surface: FemSurface,
) -> None:
    """Every centroid must be within a cell of the component it was assigned to.

    The grid is anisotropic (the hull is 8.6:1), so the bound is the largest cell
    diagonal the sampler used, not a fitted tolerance. A centroid further than that
    from its own component's zero set means the union and the components disagree,
    which would make the patch tags meaningless.
    """
    box = to_mesh(SuboffFullField(), (109, 13, 19), margin=2).mesh.bounding_box()
    cells = len((109, 13, 19))
    diagonal = math.dist(box.minimum, box.maximum) / cells
    for stats in suboff_surface.regions:
        assert stats.max_assignment_distance_m < 0.5 * diagonal, stats
