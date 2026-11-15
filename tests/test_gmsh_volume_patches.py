"""gmsh 2.2 volume export with named physical groups.

The SU2 marker path (``tests/test_volume_mesh_su2_markers.py``) proved that the body grid
leaves the package with its wall *named*.  This file pins the same contract for gmsh, the
format the open-source solver stack reads, and it insists on three separate agreements
rather than one:

1. **Self round trip is exact, not approximate** -- the same vertices (bit for bit, because
   coordinates are written with ``repr``), the same tetrahedra, the same triangles per patch.
2. **A closed-form anchor** -- a box has sides ``L*B``, ``B*H``, ``2*L*H``; those areas and
   its volume are known without any mesh, so agreement is a measurement, not self-comparison.
3. **Two foreign readers** -- meshio's gmsh reader and this package's own FEM-handoff reader
   (``tensorcad.solver_interface.read_gmsh_volume``) both recover the volume, the cell counts
   and the physical names.  A file only its own writer can read would prove nothing.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from tensorcad import volume_mesh as vmlib
from tensorcad.meshio_exchange import (
    GMSH_SOLID_GROUP,
    GMSH_VOLUME_VERSION,
    PATCH_NAMES,
    MeshioFormatError,
    read_gmsh_volume_groups,
    tetrahedral_volume,
    volume_mesh_patches,
    volume_mesh_to_tetrahedral,
    write_gmsh_volume_mesh,
)

pytest.importorskip("meshio")


# ---------------------------------------------------------------------------------------
# Bodies
# ---------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def box() -> vmlib.VolumeMesh:
    """2 x 1 x 0.6 box: every side has a closed form, so it is a zero-error control."""
    grid = vmlib.box_section_grid(length=2.0, beam=1.0, draft=0.6, nx=5, nz=4)
    return vmlib.build_volume_mesh(grid, first_layer=0.05, layers=3)


@pytest.fixture(scope="module")
def hull() -> vmlib.VolumeMesh:
    """Wigley-type body: breadth vanishes at both ends and at the keel (collapsed columns)."""
    length, beam, draft = 10.0, 1.0, 0.625
    xs = tuple(-length / 2 + i * (length / 15) for i in range(16))
    zs = tuple(-draft + j * (draft / 5) for j in range(6))

    def breadth(i: int, j: int) -> float:
        x, z = xs[i], zs[j]
        return max(0.0, (beam / 2) * (1 - (2 * x / length) ** 2) * (1 - (z / draft) ** 2))

    return vmlib.build_volume_mesh(
        vmlib.SectionGrid(xs=xs, zs=zs, breadth=breadth), first_layer=0.02, layers=3
    )


def _relative(measured: float, reference: float) -> float:
    return abs(measured - reference) / abs(reference)


# ---------------------------------------------------------------------------------------
# 1. Self round trip: exact
# ---------------------------------------------------------------------------------------


def test_box_file_parses_to_the_same_geometry(box, tmp_path: Path) -> None:
    source = volume_mesh_patches(box)
    path = write_gmsh_volume_mesh(tmp_path / "box.msh", box, patches=source)
    back = read_gmsh_volume_groups(path)

    assert back.vertices == source.vertices, "coordinates must survive repr/float exactly"
    assert back.tets == source.tets
    assert back.patches == source.patches
    assert back.names() == source.names()


def test_header_states_the_version_and_the_source_units(box, tmp_path: Path) -> None:
    path = write_gmsh_volume_mesh(tmp_path / "box.msh", box)
    text = path.read_text(encoding="utf-8")
    assert text.splitlines()[1].split()[0] == GMSH_VOLUME_VERSION
    # ``VolumeMesh`` carries no units field, so the bridge declares metres and the ship
    # frame (the defaults of ``volume_mesh_to_tetrahedral``) rather than inventing one.
    back = read_gmsh_volume_groups(path)
    assert back.units == "m"
    assert "ship frame" in back.coordinate_frame


def test_write_is_deterministic(box, tmp_path: Path) -> None:
    first = write_gmsh_volume_mesh(tmp_path / "a.msh", box)
    second = write_gmsh_volume_mesh(tmp_path / "b.msh", box)
    assert first.read_bytes() == second.read_bytes()


def test_patch_names_are_the_ones_volume_mesh_patches_makes(box, tmp_path: Path) -> None:
    back = read_gmsh_volume_groups(write_gmsh_volume_mesh(tmp_path / "box.msh", box))
    assert set(back.names()) <= set(PATCH_NAMES)
    assert back.names() == ("bottom", "top", "fore", "aft", "hull")


# ---------------------------------------------------------------------------------------
# 2. Closed-form anchors
# ---------------------------------------------------------------------------------------


def test_box_volume_and_side_areas_match_closed_form(box, tmp_path: Path) -> None:
    length, beam, draft = 2.0, 1.0, 0.6
    source = volume_mesh_patches(box)
    back = read_gmsh_volume_groups(write_gmsh_volume_mesh(tmp_path / "box.msh", box))

    assert _relative(back.tet_volume(), length * beam * draft) < 1e-12
    assert _relative(back.skin_volume(), length * beam * draft) < 1e-12
    areas = back.areas()
    assert _relative(areas["bottom"], length * beam) < 1e-12
    assert _relative(areas["top"], length * beam) < 1e-12
    assert _relative(areas["fore"], beam * draft) < 1e-12
    assert _relative(areas["aft"], beam * draft) < 1e-12
    assert _relative(areas["hull"], 2 * length * draft) < 1e-12
    # The file's own skin must bound the file's own tetrahedra, not merely the source's.
    assert back.is_watertight()
    for name, faces in back.patches:
        assert len(faces) == source.counts()[name]


def test_hull_skin_is_watertight_and_encloses_its_tetrahedra(hull, tmp_path: Path) -> None:
    source = volume_mesh_patches(hull)
    back = read_gmsh_volume_groups(write_gmsh_volume_mesh(tmp_path / "hull.msh", hull))

    assert back.is_watertight(), "a slit face would be used in one direction only"
    assert _relative(back.skin_volume(), back.tet_volume()) < 1e-12
    assert _relative(back.tet_volume(), source.tet_volume()) < 1e-12
    # The collapsed end stations contribute zero-volume cells; they are dropped, not exported.
    assert source.degenerate_cells > 0
    assert len(back.tets) == len(volume_mesh_to_tetrahedral(hull).tets)
    # Every patch the grid plane naming produced survives the file, with its triangles.
    assert back.counts() == source.counts()
    for name in back.names():
        assert back.area(name) > 0.0
    assert all(area > 0.0 for area in back.areas().values())


# ---------------------------------------------------------------------------------------
# 3. Foreign readers
# ---------------------------------------------------------------------------------------


def test_meshio_reads_the_same_volume_and_the_same_tags(hull, tmp_path: Path) -> None:
    import meshio

    source = volume_mesh_patches(hull)
    path = write_gmsh_volume_mesh(tmp_path / "hull.msh", hull)
    mio = meshio.read(str(path))

    blocks = {cell.type: cell.data for cell in mio.cells}
    assert set(blocks) == {"triangle", "tetra"}, "meshio must see one boundary and one volume block"

    points = [tuple(float(value) for value in row) for row in mio.points]
    volume = sum(tetrahedral_volume(points, tuple(int(i) for i in tet)) for tet in blocks["tetra"])
    assert _relative(volume, source.tet_volume()) < 1e-12

    # meshio keeps the tags as integers plus a name table; both must agree with our naming.
    names = {name: tuple(int(value) for value in data) for name, data in mio.field_data.items()}
    assert names[GMSH_SOLID_GROUP] == (1, 3)
    for tag, (name, _faces) in enumerate(source.patches, start=2):
        assert names[name] == (tag, 2)

    physical = np.asarray(mio.cell_data["gmsh:physical"][0])
    assert set(physical.tolist()) == set(range(2, len(source.patches) + 2))
    per_tag_area: dict[int, float] = {}
    for element, tag in zip(blocks["triangle"], physical.tolist(), strict=True):
        face = tuple(int(i) for i in element)
        per_tag_area[tag] = per_tag_area.get(tag, 0.0) + source_area(points, face)
    for tag, (name, _faces) in enumerate(source.patches, start=2):
        assert _relative(per_tag_area[tag], source.area(name)) < 1e-12


def source_area(points, face) -> float:
    """Triangle area from plain coordinates (independent of the module's helper)."""
    p0, p1, p2 = (np.asarray(points[index]) for index in face)
    return 0.5 * float(np.linalg.norm(np.cross(p1 - p0, p2 - p0)))


def test_the_packages_fem_handoff_reads_the_file(hull, tmp_path: Path) -> None:
    """The gmsh bundle consumer (``solver_interface``) must accept what we export."""
    from tensorcad.solver_interface import read_gmsh_volume

    source = volume_mesh_patches(hull)
    path = write_gmsh_volume_mesh(tmp_path / "hull.msh", hull)
    bundle = read_gmsh_volume(path)

    assert len(bundle["tets"]) == len(source.tets)
    assert len(bundle["triangles"]) == len(source.skin())
    volume = sum(tetrahedral_volume(bundle["nodes"], tet) for tet in bundle["tets"])
    assert _relative(volume, source.tet_volume()) < 1e-12
    assert set(source.names()) <= set(bundle["physical_names"].values())
    assert bundle["physical_names"][1] == GMSH_SOLID_GROUP


def test_a_foreign_writer_is_read_with_its_names(box, tmp_path: Path) -> None:
    """A gmsh 2.2 file this package did *not* write, with the same patch names."""
    import meshio

    source = volume_mesh_patches(box)
    triangles = [face for _name, faces in source.patches for face in faces]
    physical_tri = np.array(
        [tag for tag, (_name, faces) in enumerate(source.patches, start=2) for _ in faces]
    )
    foreign = meshio.Mesh(
        np.array(source.vertices, dtype=float),
        [("triangle", np.array(triangles)), ("tetra", np.array(source.tets))],
        cell_data={
            "gmsh:physical": [physical_tri, np.ones(len(source.tets), dtype=int)],
        },
        field_data={
            **{name: np.array([tag, 2]) for tag, (name, _f) in enumerate(source.patches, start=2)},
            GMSH_SOLID_GROUP: np.array([1, 3]),
        },
    )
    path = tmp_path / "foreign.msh"
    meshio.write(str(path), foreign, file_format="gmsh22", binary=False)

    back = read_gmsh_volume_groups(path)
    assert back.names() == source.names()
    assert _relative(back.tet_volume(), source.tet_volume()) < 1e-12
    assert back.counts() == source.counts()


# ---------------------------------------------------------------------------------------
# 4. Refusals: named, not silent
# ---------------------------------------------------------------------------------------


def test_a_solid_name_that_collides_with_a_patch_is_refused(box, tmp_path: Path) -> None:
    with pytest.raises(MeshioFormatError, match="collides with a patch name"):
        write_gmsh_volume_mesh(tmp_path / "box.msh", box, solid_name="hull")


def test_reader_refuses_a_surface_only_file(tmp_path: Path) -> None:
    path = tmp_path / "surface.msh"
    path.write_text(
        "\n".join(
            [
                "$MeshFormat",
                "2.2 0 8",
                "$EndMeshFormat",
                "$PhysicalNames",
                "1",
                '2 1 "wall"',
                "$EndPhysicalNames",
                "$Nodes",
                "3",
                "1 0.0 0.0 0.0",
                "2 1.0 0.0 0.0",
                "3 0.0 1.0 0.0",
                "$EndNodes",
                "$Elements",
                "1",
                "1 2 2 1 1 1 2 3",
                "$EndElements",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(MeshioFormatError, match="no tetrahedra"):
        read_gmsh_volume_groups(path)


def test_reader_refuses_a_binary_file(hull, tmp_path: Path) -> None:
    import meshio

    source = volume_mesh_patches(hull)
    triangles = [face for _name, faces in source.patches for face in faces]
    binary = meshio.Mesh(
        np.array(source.vertices, dtype=float),
        [("triangle", np.array(triangles)), ("tetra", np.array(source.tets))],
    )
    path = tmp_path / "binary.msh"
    meshio.write(str(path), binary, file_format="gmsh22", binary=True)
    with pytest.raises(MeshioFormatError, match="not gmsh ASCII"):
        read_gmsh_volume_groups(path)


def test_reader_refuses_another_gmsh_version(tmp_path: Path) -> None:
    path = tmp_path / "v4.msh"
    path.write_text("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n", encoding="utf-8")
    with pytest.raises(MeshioFormatError, match=r"only 2\.x ASCII"):
        read_gmsh_volume_groups(path)


def test_reader_refuses_an_unmarked_element(tmp_path: Path) -> None:
    """An element with zero tags belongs to no patch; a solver could not bind it."""
    path = tmp_path / "untagged.msh"
    path.write_text(
        "\n".join(
            [
                "$MeshFormat",
                "2.2 0 8",
                "$EndMeshFormat",
                "$Nodes",
                "4",
                "1 0.0 0.0 0.0",
                "2 1.0 0.0 0.0",
                "3 0.0 1.0 0.0",
                "4 0.0 0.0 1.0",
                "$EndNodes",
                "$Elements",
                "1",
                "1 4 0 1 2 3 4",
                "$EndElements",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(MeshioFormatError, match="no physical tag"):
        read_gmsh_volume_groups(path)


def test_reader_refuses_quad_boundary_cells(tmp_path: Path) -> None:
    path = tmp_path / "quads.msh"
    path.write_text(
        "\n".join(
            [
                "$MeshFormat",
                "2.2 0 8",
                "$EndMeshFormat",
                "$PhysicalNames",
                "2",
                '2 2 "bottom"',
                '3 1 "solid"',
                "$EndPhysicalNames",
                "$Nodes",
                "5",
                "1 0.0 0.0 0.0",
                "2 1.0 0.0 0.0",
                "3 1.0 1.0 0.0",
                "4 0.0 1.0 0.0",
                "5 0.0 0.0 1.0",
                "$EndNodes",
                "$Elements",
                "2",
                "1 3 2 2 2 1 2 3 4",
                "2 4 2 3 3 1 2 3 5",
                "$EndElements",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    with pytest.raises(MeshioFormatError, match="cannot be expressed as a tetrahedron"):
        read_gmsh_volume_groups(path)


def test_a_single_tag_element_still_reads(box, tmp_path: Path) -> None:
    """gmsh allows one tag (physical only); the reader must not depend on the second."""
    source = volume_mesh_patches(box)
    written = write_gmsh_volume_mesh(tmp_path / "box.msh", box).read_text(encoding="utf-8")
    rewritten: list[str] = []
    in_elements = False
    for line in written.splitlines():
        stripped = line.strip()
        if stripped == "$Elements":
            in_elements = True
            rewritten.append(line)
            continue
        if stripped == "$EndElements":
            in_elements = False
            rewritten.append(line)
            continue
        tokens = line.split()
        if in_elements and len(tokens) > 6 and tokens[1] in {"2", "4"}:
            element_id, kind, physical, nodes = tokens[0], tokens[1], tokens[3], tokens[5:]
            rewritten.append(" ".join([element_id, kind, "1", physical, *nodes]))
        else:
            rewritten.append(line)
    path = tmp_path / "onetag.msh"
    path.write_text("\n".join(rewritten) + "\n", encoding="utf-8")

    back = read_gmsh_volume_groups(path)
    assert back.tets == source.tets
    assert back.counts() == source.counts()
