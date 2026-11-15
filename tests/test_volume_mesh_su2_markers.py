"""The SU2 file a solver gets must carry a *named*, closed skin -- not just tetrahedra.

Why these tests exist
=====================
``volume_mesh_to_tetrahedral`` closed the format-matrix gap at the level of *cells*: the
body grid becomes tetrahedra and SU2/CGNS accept them.  What it did not do is give the file
a boundary.  ``write_meshio(..., file_format="su2")`` on that mesh writes ``NMARK= 0``
(measured, not assumed -- see ``test_unmarked_su2_has_no_markers``), because meshio's SU2
writer emits markers only for triangle blocks carrying integer cell data and a
:class:`TetrahedralMesh` has no boundary.  A solver config that says
``MARKER_HEATFLUX= (hull, 0.0)`` then has nothing to bind to.

Naming the skin is where two defects of the grid have to be dealt with, which is why the
assertions here are identities rather than smoke tests:

* **The skin must be closed.**  The builder keeps a distinct node for every
  ``(station, level, breadth)`` index even where the breadth collapses, so dropping the
  zero-volume cells exposes the *interior* faces of the surviving cells.  On the Wigley body
  below that is 1050 of the 2184 triangles that look like boundary: coincident pairs of a
  zero-thickness slit.  ``volume_mesh_patches`` welds the coincident points, the pairs cancel,
  and the result is checked with the divergence theorem and an edge census -- the skin must
  enclose exactly the volume the tetrahedra measure, and every edge must appear once in each
  direction.
* **The areas must be the geometry's.**  A box has a closed-form surface area, per side; a
  Wigley body has a wetted-surface integral that a quadrature over the *breadth function*
  computes independently of any mesh.  Both are compared below.

The SU2 file itself is read back three ways: by the module's own parser (tags as names), by
meshio's reader (an implementation that knows nothing about us), and by re-checking the
element counts against the mesh that produced them.
"""
from __future__ import annotations

import math
import warnings

import pytest

from tensorcad import volume_mesh as vmlib
from tensorcad.meshio_exchange import (
    PATCH_NAMES,
    MeshioFormatError,
    read_su2_markers,
    volume_mesh_patches,
    volume_mesh_to_tetrahedral,
    write_meshio,
    write_su2_volume_mesh,
)

REL = 1e-12
BOX = {"length": 2.0, "beam": 1.0, "draft": 0.6}
HULL = {"length": 10.0, "beam": 1.0, "draft": 0.625}


# ---------------------------------------------------------------------------------------
# Fixtures: one analytic box, one real (collapsed-column) hull body
# ---------------------------------------------------------------------------------------

def box_volume_mesh() -> vmlib.VolumeMesh:
    grid = vmlib.box_section_grid(length=BOX["length"], beam=BOX["beam"], draft=BOX["draft"],
                                  nx=5, nz=4)
    return vmlib.build_volume_mesh(grid, first_layer=0.05, growth=1.2, layers=3)


def hull_section_grid(
    length: float = HULL["length"], beam: float = HULL["beam"], draft: float = HULL["draft"],
    nx: int = 25, nz: int = 9,
) -> vmlib.SectionGrid:
    """Wigley-type half breadth ``(B/2)(1-(2x/L)^2)(1-(z/T)^2)``: collapses at both ends."""
    xs = tuple(-length / 2 + i * (length / (nx - 1)) for i in range(nx))
    zs = tuple(j * (draft / (nz - 1)) for j in range(nz))

    def breadth(i: int, j: int) -> float:
        x, z = xs[i], zs[j]
        return max(0.0, (beam / 2) * (1 - (2 * x / length) ** 2) * (1 - (z / draft) ** 2))

    return vmlib.SectionGrid(xs=xs, zs=zs, breadth=breadth)


@pytest.fixture(scope="module")
def box_mesh() -> vmlib.VolumeMesh:
    return box_volume_mesh()


@pytest.fixture(scope="module")
def hull_mesh() -> vmlib.VolumeMesh:
    return vmlib.build_volume_mesh(hull_section_grid(), first_layer=0.004, growth=1.25, layers=4)


@pytest.fixture(scope="module")
def box_patches(box_mesh: vmlib.VolumeMesh):
    return volume_mesh_patches(box_mesh)


@pytest.fixture(scope="module")
def hull_patches(hull_mesh: vmlib.VolumeMesh):
    return volume_mesh_patches(hull_mesh)


# ---------------------------------------------------------------------------------------
# Geometry helpers used by the assertions (independent of the module under test)
# ---------------------------------------------------------------------------------------

def triangle_area(p: tuple[float, float, float], q: tuple[float, float, float],
                  r: tuple[float, float, float]) -> float:
    ux, uy, uz = q[0] - p[0], q[1] - p[1], q[2] - p[2]
    vx, vy, vz = r[0] - p[0], r[1] - p[1], r[2] - p[2]
    cx, cy, cz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
    return 0.5 * math.sqrt(cx * cx + cy * cy + cz * cz)


def triangle_normal(p, q, r) -> tuple[float, float, float]:
    ux, uy, uz = q[0] - p[0], q[1] - p[1], q[2] - p[2]
    vx, vy, vz = r[0] - p[0], r[1] - p[1], r[2] - p[2]
    cx, cy, cz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
    norm = math.sqrt(cx * cx + cy * cy + cz * cz)
    return (cx / norm, cy / norm, cz / norm)


def raw_boundary_faces(tets, vertices) -> list[tuple[int, int, int]]:
    """Faces used by exactly one tetrahedron of the *unwelded*, un-welded kept split."""
    del vertices
    counts: dict[tuple[int, int, int], int] = {}
    for tet in tets:
        for a in range(4):
            for b in range(a + 1, 4):
                for c in range(b + 1, 4):
                    face = tuple(sorted((tet[a], tet[b], tet[c])))
                    counts[face] = counts.get(face, 0) + 1
    return [face for face, count in counts.items() if count == 1]


def wetted_area_quadrature(length: float, beam: float, draft: float,
                           nx: int = 400, nz: int = 200) -> float:
    """``2 \\iint sqrt(1 + b_x^2 + b_z^2) dx dz`` by the midpoint rule.

    The integrand comes from the *breadth function*, not from the mesh: this is the
    independent anchor the mesh's ``hull`` patch area is compared against.
    """
    total = 0.0
    dx, dz = length / nx, draft / nz
    for i in range(nx):
        x = -length / 2 + (i + 0.5) * dx
        for j in range(nz):
            z = (j + 0.5) * dz
            bx = (beam / 2) * (-8 * x / length**2) * (1 - (z / draft) ** 2)
            bz = (beam / 2) * (1 - (2 * x / length) ** 2) * (-2 * z / draft**2)
            total += math.sqrt(1.0 + bx * bx + bz * bz)
    return 2.0 * total * dx * dz


# ---------------------------------------------------------------------------------------
# The gap, measured: an unmarked SU2 file has no markers
# ---------------------------------------------------------------------------------------

def test_unmarked_su2_has_no_markers(tmp_path) -> None:
    """The pre-existing export path writes ``NMARK= 0`` -- this is the state being fixed."""
    mesh = volume_mesh_to_tetrahedral(box_volume_mesh())
    path = write_meshio(tmp_path / "plain.su2", mesh, file_format="su2")
    lines = path.read_text(encoding="utf-8").splitlines()
    assert "NMARK= 0" in lines
    assert not any(line.startswith("MARKER_TAG=") for line in lines)


# ---------------------------------------------------------------------------------------
# The skin: closed, complete, and of the right measure
# ---------------------------------------------------------------------------------------

@pytest.mark.parametrize("fixture", ["box_patches", "hull_patches"])
def test_skin_is_watertight(fixture: str, request: pytest.FixtureRequest) -> None:
    """Every edge of the skin is used once in each direction: no holes, no folds."""
    patches = request.getfixturevalue(fixture)
    counts = patches.edge_counts()
    for (tail, head), count in counts.items():
        assert count == 1, f"edge {tail}->{head} used {count} times"
        assert counts.get((head, tail)) == 1, f"edge {tail}->{head} has no return edge"
    assert patches.is_watertight()


@pytest.mark.parametrize("fixture", ["box_patches", "hull_patches"])
def test_skin_encloses_the_tetrahedra(fixture: str, request: pytest.FixtureRequest) -> None:
    """Divergence theorem on the skin == the tetrahedra's volume == the grid's volume."""
    patches = request.getfixturevalue(fixture)
    mesh = request.getfixturevalue(fixture.replace("_patches", "_mesh"))
    reference = mesh.volume()
    assert patches.tet_volume() == pytest.approx(reference, rel=REL)
    assert patches.skin_volume() == pytest.approx(reference, rel=REL)


def test_patches_partition_the_skin(box_patches, hull_patches) -> None:
    """Every skin triangle belongs to one patch; the patches are the whole skin, once."""
    for patches in (box_patches, hull_patches):
        assert patches.names() == tuple(
            name for name in PATCH_NAMES if name in patches.names()
        ), "patch order must follow PATCH_NAMES"
        seen: set[tuple[int, int, int]] = set()
        total = 0
        for name, faces in patches.patches:
            assert faces, f"patch {name} is empty and should not have been emitted"
            for face in faces:
                key = tuple(sorted(face))
                assert key not in seen, f"triangle {key} is in two patches"
                seen.add(key)
                assert len(set(face)) == 3
                assert max(face) < len(patches.vertices)
            total += len(faces)
        assert total == len(seen) == len(patches.skin())
        assert sum(patches.counts().values()) == len(patches.skin())


def test_every_output_point_is_used(box_patches, hull_patches) -> None:
    """No orphan points: the file's ``NPOIN`` is what the cells actually reference."""
    for patches in (box_patches, hull_patches):
        used = {vertex for tet in patches.tets for vertex in tet}
        assert used == set(range(len(patches.vertices)))


def test_coincident_points_are_welded(box_mesh, hull_mesh, box_patches, hull_patches) -> None:
    """The box has no coincident nodes; the collapsed hull body has hundreds."""
    assert box_patches.merged_points == 0
    assert hull_patches.merged_points > 0
    # The weld removes exactly one point per duplicate coordinate, and what is exported is
    # the welded set minus points no kept tetrahedron references.
    distinct = len({tuple(vertex) for vertex in hull_mesh.vertices})
    assert distinct == len(hull_mesh.vertices) - hull_patches.merged_points
    assert len(hull_patches.vertices) <= distinct
    assert len(box_patches.vertices) == len(box_mesh.vertices)


def test_slit_faces_cancel_exactly(hull_mesh, hull_patches) -> None:
    """The 1050 zero-thickness slit triangles are pairs, and the skin is what is left.

    Recomputed here from the raw split, without the module's welding: group the raw
    boundary triangles by their rounded vertex coordinates; every coordinate signature that
    occurs twice is a slit pair (both copies are dropped), and what occurs once is skin.
    """
    kept = volume_mesh_to_tetrahedral(hull_mesh)
    raw = raw_boundary_faces(kept.tets, hull_mesh.vertices)
    signature: dict[tuple, int] = {}
    for face in raw:
        key = tuple(sorted(tuple(round(c, 12) for c in hull_mesh.vertices[v]) for v in face))
        signature[key] = signature.get(key, 0) + 1
    pairs = sum(1 for count in signature.values() if count == 2)
    singles = sum(1 for count in signature.values() if count == 1)
    assert pairs > 0, "this fixture is supposed to exercise the collapsed columns"
    assert singles == len(hull_patches.skin())
    assert len(raw) - 2 * pairs == len(hull_patches.skin())


# ---------------------------------------------------------------------------------------
# The names mean something geometric
# ---------------------------------------------------------------------------------------

def test_patch_faces_lie_on_their_planes(box_mesh, hull_mesh, box_patches, hull_patches) -> None:
    """``bottom``/``top``/``fore``/``aft`` are the grid's own planes; ``hull`` is off them."""
    for mesh, patches in ((box_mesh, box_patches), (hull_mesh, hull_patches)):
        xs, zs = mesh.station_lines, mesh.level_lines
        for name, faces in patches.patches:
            for face in faces:
                points = [patches.vertices[index] for index in face]
                ys, zs_of, xs_of = ([p[1] for p in points], [p[2] for p in points],
                                    [p[0] for p in points])
                tol = 1e-9 * max(1.0, max(xs) - min(xs))
                on = {
                    "symmetry": all(abs(y) <= tol for y in ys),
                    "bottom": all(abs(z - zs[0]) <= tol for z in zs_of),
                    "top": all(abs(z - zs[-1]) <= tol for z in zs_of),
                    "fore": all(abs(x - xs[0]) <= tol for x in xs_of),
                    "aft": all(abs(x - xs[-1]) <= tol for x in xs_of),
                }
                if name == "hull":
                    assert not any(on.values()), f"hull triangle {face} lies on a grid plane"
                else:
                    assert on[name], f"{name} triangle {face} is off its plane ({on})"


def test_box_patch_areas_are_analytic(box_patches) -> None:
    """A box's sides have closed-form areas -- a zero-error control for the patch split."""
    length, beam, depth = BOX["length"], BOX["beam"], BOX["draft"]
    areas = box_patches.areas()
    assert areas["bottom"] == pytest.approx(length * beam, rel=REL)
    assert areas["top"] == pytest.approx(length * beam, rel=REL)
    assert areas["fore"] == pytest.approx(beam * depth, rel=REL)
    assert areas["aft"] == pytest.approx(beam * depth, rel=REL)
    assert areas["hull"] == pytest.approx(2.0 * length * depth, rel=REL)
    total = 2.0 * (length * beam + length * depth + beam * depth)
    assert sum(areas.values()) == pytest.approx(total, rel=REL)
    assert "symmetry" not in areas, "a full-breadth box grid has no centreplane skin"


def test_box_normals_point_outward(box_patches) -> None:
    """Winding, not just coordinates: each side's normal points away from the body."""
    expected = {
        "bottom": (0.0, 0.0, -1.0),
        "top": (0.0, 0.0, 1.0),
        "fore": (-1.0, 0.0, 0.0),
        "aft": (1.0, 0.0, 0.0),
    }
    for name, direction in expected.items():
        for face in box_patches.faces(name):
            normal = triangle_normal(*(box_patches.vertices[i] for i in face))
            assert sum(n * d for n, d in zip(normal, direction, strict=True)) == pytest.approx(1.0, abs=1e-12)
    for face in box_patches.faces("hull"):
        points = [box_patches.vertices[i] for i in face]
        normal = triangle_normal(*points)
        side = 1.0 if sum(p[1] for p in points) > 0 else -1.0
        assert side * normal[1] == pytest.approx(1.0, abs=1e-12)


def test_hull_area_matches_quadrature_of_the_breadth_function(hull_patches) -> None:
    """Wetted area vs an integral over ``b(x, z)`` -- no mesh involved in the reference."""
    continuum = wetted_area_quadrature(HULL["length"], HULL["beam"], HULL["draft"])
    mesh_area = hull_patches.area("hull")
    error = abs(mesh_area - continuum) / continuum
    assert error < 0.01, f"hull patch area {mesh_area} vs quadrature {continuum} ({error:.2%})"
    # The bottom patch is 2*integral(b(x, 0) dx) = B*L*2/3 for this breadth function; the
    # mesh integrates the piecewise-linear interpolant, which is slightly short of it.
    bottom = hull_patches.area("bottom")
    analytic = HULL["beam"] * HULL["length"] * 2.0 / 3.0
    assert abs(bottom - analytic) / analytic < 0.01, f"bottom {bottom} vs {analytic}"


def test_hull_patch_areas_sum_to_the_skin_area(hull_patches) -> None:
    """The per-patch areas are a partition of one surface, so they add up to its measure."""
    total = sum(
        triangle_area(*(hull_patches.vertices[i] for i in face)) for face in hull_patches.skin()
    )
    assert sum(hull_patches.areas().values()) == pytest.approx(total, rel=1e-15)


# ---------------------------------------------------------------------------------------
# The file
# ---------------------------------------------------------------------------------------

def test_su2_file_layout(tmp_path, box_mesh, box_patches) -> None:
    """Header counts, element types and marker order as SU2 defines them."""
    path = write_su2_volume_mesh(tmp_path / "box.su2", box_mesh, patches=box_patches)
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "NDIME= 3"
    assert lines[1] == f"NPOIN= {len(box_patches.vertices)}"
    assert f"NELEM= {len(box_patches.tets)}" in lines
    assert f"NMARK= {len(box_patches.patches)}" in lines

    counts = box_patches.counts()
    for name, count in counts.items():
        assert f"MARKER_TAG= {name}" in lines
        assert f"MARKER_ELEMS= {count}" in lines
    assert sum(counts.values()) == len(box_patches.skin())

    # Every element line starts with its SU2 type code and references only existing points.
    volume_lines = 0
    boundary_lines = 0
    marker = False
    for line in lines:
        if line.startswith("MARKER_TAG="):
            marker = True
        parts = line.split()
        if len(parts) == 5 and not marker:
            assert parts[0] == "10", f"volume element type must be 10 (tetra): {line}"
            volume_lines += 1
            assert all(0 <= int(i) < len(box_patches.vertices) for i in parts[1:])
        elif len(parts) == 4 and marker and parts[0].isdigit():
            assert parts[0] == "5", f"marker element type must be 5 (triangle): {line}"
            boundary_lines += 1
            assert all(0 <= int(i) < len(box_patches.vertices) for i in parts[1:])
    assert volume_lines == len(box_patches.tets)
    assert boundary_lines == len(box_patches.skin())


def test_own_reader_recovers_the_patches(tmp_path, box_mesh, box_patches) -> None:
    """Tags survive as *names*: this is what a solver config binds to."""
    path = write_su2_volume_mesh(tmp_path / "box.su2", box_mesh, patches=box_patches)
    read_back = read_su2_markers(path)
    assert tuple(read_back) == box_patches.names()
    for name, faces in box_patches.patches:
        assert set(read_back[name]) == set(faces)


def test_writer_reuses_supplied_patches(tmp_path, hull_mesh, hull_patches) -> None:
    """Passing the patches in must produce the same file as computing them again."""
    computed = write_su2_volume_mesh(tmp_path / "a.su2", hull_mesh)
    supplied = write_su2_volume_mesh(tmp_path / "b.su2", hull_mesh, patches=hull_patches)
    assert computed.read_text() == supplied.read_text()
    assert volume_mesh_patches(hull_mesh).counts() == hull_patches.counts()


def test_meshio_reads_the_marked_file(tmp_path, hull_mesh, hull_patches) -> None:
    """An implementation that knows nothing about us must still recover the mesh.

    meshio cannot keep string tags (it renumbers them, warning as it goes), so the check is
    made on what it *can* be trusted with: the tetrahedra, their volume, and the boundary
    triangles as geometry.
    """
    meshio = pytest.importorskip("meshio")
    path = write_su2_volume_mesh(tmp_path / "hull.su2", hull_mesh)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        read_back = meshio.read(str(path), file_format="su2")
    blocks = {block.type: block for block in read_back.cells}
    assert len(blocks["tetra"].data) == len(hull_patches.tets)
    assert len(blocks["triangle"].data) == len(hull_patches.skin())

    points = [tuple(map(float, point)) for point in read_back.points]
    volume = 0.0
    for tet in blocks["tetra"].data:
        p0, p1, p2, p3 = (points[index] for index in tet)
        u = [p1[i] - p0[i] for i in range(3)]
        v = [p2[i] - p0[i] for i in range(3)]
        w = [p3[i] - p0[i] for i in range(3)]
        det = (
            u[0] * (v[1] * w[2] - v[2] * w[1])
            - u[1] * (v[0] * w[2] - v[2] * w[0])
            + u[2] * (v[0] * w[1] - v[1] * w[0])
        )
        volume += det / 6.0
    assert volume == pytest.approx(hull_patches.tet_volume(), rel=1e-12)

    def signature(triangle: tuple[int, int, int], vertices) -> tuple:
        return tuple(sorted(tuple(round(c, 12) for c in vertices[index]) for index in triangle))

    ours = {signature(face, hull_patches.vertices) for face in hull_patches.skin()}
    theirs = {
        signature(tuple(int(i) for i in face), points) for face in blocks["triangle"].data
    }
    assert theirs == ours


def test_split_tets_false_is_refused_with_the_remedy() -> None:
    """An unsplit grid cannot be exported: the refusal has to say what to do instead."""
    grid = vmlib.box_section_grid(length=2.0, beam=1.0, draft=0.6, nx=4, nz=3)
    unsplit = vmlib.build_volume_mesh(grid, first_layer=0.05, layers=2, split_tets=False)
    with pytest.raises(MeshioFormatError, match="split_tets=True"):
        volume_mesh_patches(unsplit)


def test_real_hull_generator_body_is_exportable(tmp_path) -> None:
    """The package's own hull generator, not just a synthetic breadth function."""
    from tensorcad.hull import HullParameters

    mesh = vmlib.volume_mesh_from_hull(HullParameters(nx=13, nz=7), first_layer=0.05,
                                       growth=1.2, layers=1)
    patches = volume_mesh_patches(mesh)
    assert patches.merged_points > 0
    assert patches.is_watertight()
    assert patches.tet_volume() == pytest.approx(mesh.volume(), rel=REL)
    assert patches.skin_volume() == pytest.approx(mesh.volume(), rel=REL)
    path = write_su2_volume_mesh(tmp_path / "hull_body.su2", mesh, patches=patches)
    assert set(read_su2_markers(path)) == set(patches.names())
    assert patches.counts()["hull"] > 0
