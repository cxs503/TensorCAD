"""meshio format matrix: round-trip fidelity for 6 formats × 4 geometries = 24 combos.

Verifies that exporting a mesh through meshio and reading it back preserves:

- **Volume** (divergence-theorem volume, relative error < 1e-3; exact formats
  < 1e-12)
- **Topology** (vertex count and cell count exact for every format; exact vertex
  positions and cell indices for the lossless formats OBJ/VTK/Gmsh/SU2/CGNS)
- **Orientation** (sign of the volume unchanged → outward normals / positive
  tetrahedra preserved)

Matrix layout
=============
Each geometry is exercised in the topology its formats speak:

* surface formats (STL, OBJ, VTK, Gmsh) get the triangle surface mesh;
* volume formats (SU2, CGNS) get the tetrahedralisation of the same geometry
  (:func:`tensorcad.meshio_exchange.tetrahedralise_surface`).

6 formats × 4 geometries = 24 combinations, all green.  Before this change
SU2 and CGNS were declared NOT_SUPPORTED and only 12 of the 18 combinations
were actually exercised.

The fourth geometry is a **real ship hull** -- ``tensorcad.wigley``'s own
validated Wigley generator (L=10 m, B=1 m, T=0.625 m, 1001 vertices / 1998
triangles, closed by a deck cap).  It is in the matrix because the previous
version's volume path was only ever demonstrated on convex toy bodies (cube,
sphere, capsule), while the module docstring claimed ship hulls work "in
practice".  That claim was false against this very hull: where the section ring
collapses to zero breadth at the stem and the post, the tip facets lie exactly
in the centreplane ``y = 0``, so every off-plane apex *inverts* them and the
strict-positivity rule refused the hull as "not star-shaped".  Certification by
volume identity (no inverted cell + tets reproduce the surface volume to 1e-9)
accepts it: the two edge-on facets span zero-volume cells, are dropped, and
counted in ``TetrahedralMesh.edge_on_facets``.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.meshio_exchange import (
    FORMATS_ALL,
    FORMATS_SURFACE,
    FORMATS_VOLUME,
    MeshioFormatError,
    TetrahedralisationError,
    TetrahedralMesh,
    mesh_to_meshio,
    meshio_to_mesh,
    meshio_to_tetrahedral,
    round_trip_meshio,
    surface_signed_volume,
    tetrahedral_volume,
    tetrahedralise_surface,
    write_meshio,
)

# ---------------------------------------------------------------------------
# Test geometries
# ---------------------------------------------------------------------------

def _make_cube() -> Mesh:
    """Unit cube, outward-oriented, 8 vertices, 12 triangles.  Volume = 1."""
    vertices = (
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    )
    faces = (
        (0, 3, 2), (0, 2, 1),  # bottom z=0
        (4, 5, 6), (4, 6, 7),  # top    z=1
        (0, 1, 5), (0, 5, 4),  # front  y=0
        (2, 3, 7), (2, 7, 6),  # back   y=1
        (1, 2, 6), (1, 6, 5),  # right  x=1
        (0, 4, 7), (0, 7, 3),  # left   x=0
    )
    return Mesh(vertices=vertices, faces=faces, units="m")


def _make_sphere(n_lat: int = 24, n_lon: int = 48) -> Mesh:
    """Unit sphere: pole fans + quad bands, oriented outward, watertight.

    Analytic volume = 4/3*pi = 4.188790...; the flat-faceted mesh is slightly
    smaller, and the fixture asserts that discretisation error is < 1 %.

    Note: the previous fixture used a full lat/lon *grid including the poles*,
    which left a ring of boundary edges (32 unpaired edges) and an inward
    normal -- a mesh whose "volume preserved" check was self-consistent but
    physically meaningless.  Poles are single vertices here.
    """
    builder = MeshBuilder(units="m")
    north = builder.add_vertex((0.0, 0.0, 1.0))
    rings: list[list[int]] = []
    for i_lat in range(1, n_lat):
        phi = math.pi * i_lat / n_lat
        ring = []
        for i_lon in range(n_lon):
            theta = 2 * math.pi * i_lon / n_lon
            ring.append(
                builder.add_vertex(
                    (math.sin(phi) * math.cos(theta), math.sin(phi) * math.sin(theta), math.cos(phi))
                )
            )
        rings.append(ring)
    south = builder.add_vertex((0.0, 0.0, -1.0))

    for i_lon in range(n_lon):
        j = (i_lon + 1) % n_lon
        builder.add_triangle(north, rings[0][i_lon], rings[0][j])
        builder.add_triangle(south, rings[-1][j], rings[-1][i_lon])
    for i_lat in range(len(rings) - 1):
        upper, lower = rings[i_lat], rings[i_lat + 1]
        for i_lon in range(n_lon):
            j = (i_lon + 1) % n_lon
            builder.add_quad(upper[i_lon], lower[i_lon], lower[j], upper[j])

    return meshlib.orient_outward(builder.build()).mesh


def _make_capsule(n_theta: int = 32, n_axial: int = 12) -> Mesh:
    """Ship-like body of revolution: cylinder + hemispherical ends.

    A submarine hull without the appendages, i.e. a *closed, star-shaped* body
    with an analytic volume: V = pi*r^2*L + 4/3*pi*r^3 with r = 0.5, L = 4.

    Why this replaced the hand-rolled Wigley shelf
    ----------------------------------------------
    Two fixtures were probed against the invariants before being trusted:

    * the old latitude/longitude sphere included the poles as full rings, which
      left 32 unpaired boundary edges and an inward normal (its "volume
      preserved" check compared the mesh with itself, so it passed anyway);
    * the old Wigley shelf was non-manifold at its knife edges (42 duplicated
      directed edges, 40 unpaired) and integrated to -0.034 m^3 -- less than
      half of the Wigley volume it was labelled with.

    A body of revolution keeps this fixture free of the trap found in the
    package's own hull output: a tapering hull puts its stem/post tip facets
    exactly in the centreplane, so an off-plane apex inverts them and an
    on-plane apex sees them edge-on.  That case is exercised deliberately by
    :func:`_make_wigley_hull` (see :func:`test_hull_tetrahedralisation_is_certified`).
    """
    radius, length = 0.5, 4.0
    builder = MeshBuilder(units="m")
    ring_size = 2 * n_axial + 2  # cylinder rings plus one ring per hemisphere step
    xs = []
    # aft hemisphere: x from -L/2-r to -L/2, radius sqrt(r^2-(x+L/2+r)^2)
    for i in range(n_axial + 1):
        phi = (math.pi / 2) * i / n_axial  # 0 at aft pole .. pi/2 at shoulder
        x = -length / 2 - radius * math.cos(phi)
        xs.append((x, radius * math.sin(phi)))
    # cylinder
    xs.append((length / 2, radius))
    # fore hemisphere (i = n_axial is the shoulder already added above)
    for i in range(n_axial - 1, -1, -1):
        phi = (math.pi / 2) * i / n_axial
        x = length / 2 + radius * math.cos(phi)
        xs.append((x, radius * math.sin(phi)))
    assert len(xs) == ring_size

    rings: list[list[int]] = []
    for x, r in xs:
        if r == 0.0:  # pole: single vertex on the axis
            rings.append([builder.add_vertex((x, 0.0, 0.0))])
            continue
        ring = []
        for k in range(n_theta):
            theta = 2 * math.pi * k / n_theta
            ring.append(builder.add_vertex((x, r * math.cos(theta), r * math.sin(theta))))
        rings.append(ring)

    def connect(r0: list[int], r1: list[int]) -> None:
        if len(r0) == 1:
            for k in range(n_theta):
                builder.add_triangle(r0[0], r1[k], r1[(k + 1) % n_theta])
        elif len(r1) == 1:
            for k in range(n_theta):
                builder.add_triangle(r1[0], r0[(k + 1) % n_theta], r0[k])
        else:
            for k in range(n_theta):
                j = (k + 1) % n_theta
                builder.add_quad(r0[k], r0[j], r1[j], r1[k])

    for i in range(len(rings) - 1):
        connect(rings[i], rings[i + 1])

    return meshlib.orient_outward(builder.build()).mesh


def _make_wigley_hull():
    """Real benchmark hull: ``tensorcad.wigley``'s closed, deck-capped Wigley.

    Not hand-rolled on purpose -- the failing fixtures in this file's history
    were hand-rolled hulls.  The package generator is the one that carries its
    own validation report and closed-form hydrostatics, so the volume reference
    below is an independent analytic number, not a self-comparison:

        V = L B T n m / ((n+1)(m+1)) = 4 L B T / 9 at n = m = 2 = 25/9 m^3
    """
    from tensorcad import wigley

    return wigley.generate_wigley().mesh


#: Analytic volumes: cube = 1, sphere = 4/3*pi, capsule = pi*r^2*L + 4/3*pi*r^3,
#: wigley hull = 4*L*B*T/9 (closed form, see :func:`tensorcad.wigley.analytic_volume`).
ANALYTIC_VOLUME = {
    "cube": 1.0,
    "sphere": 4.0 / 3.0 * math.pi,
    "capsule": math.pi * 0.5**2 * 4.0 + 4.0 / 3.0 * math.pi * 0.5**3,
    "wigley_hull": 4.0 * 10.0 * 1.0 * 0.625 / 9.0,
}


_TET_CACHE: dict[str, TetrahedralMesh] = {}


@dataclass(frozen=True)
class GeometryCase:
    name: str
    mesh: Mesh

    @property
    def tetra(self) -> TetrahedralMesh:
        cached = _TET_CACHE.get(self.name)
        if cached is None:
            cached = tetrahedralise_surface(self.mesh)
            _TET_CACHE[self.name] = cached
        return cached


GEOMETRIES: tuple[GeometryCase, ...] = (
    GeometryCase("cube", _make_cube()),
    GeometryCase("sphere", _make_sphere()),
    GeometryCase("capsule", _make_capsule()),
    GeometryCase("wigley_hull", _make_wigley_hull()),
)


def _volume(mesh: Mesh | TetrahedralMesh) -> float:
    """Enclosed volume of either representation."""
    if isinstance(mesh, TetrahedralMesh):
        return mesh.volume()
    return surface_signed_volume(mesh)


def _relative_error(a: float, b: float) -> float:
    scale = max(abs(a), abs(b), 1e-300)
    return abs(a - b) / scale


def _matrix_mesh(geom: GeometryCase, fmt: str) -> Mesh | TetrahedralMesh:
    return geom.tetra if fmt in FORMATS_VOLUME else geom.mesh


def _cell_count(mesh: Mesh | TetrahedralMesh) -> int:
    return len(mesh.tets) if isinstance(mesh, TetrahedralMesh) else len(mesh.faces)


# ---------------------------------------------------------------------------
# Tetrahedralisation (the capability that unlocks SU2/CGNS)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("geom", GEOMETRIES, ids=lambda g: g.name)
def test_tetrahedralise_conserves_volume_exactly(geom: GeometryCase) -> None:
    """Star tets must sum to the surface's signed volume (same closed boundary)."""
    tet = tetrahedralise_surface(geom.mesh)
    # Every surface facet yields one tet, except facets the apex sees exactly
    # edge-on; those are dropped (zero measure) and reported, never silent.
    assert len(tet.tets) == len(geom.mesh.faces) - tet.edge_on_facets
    assert len(tet.vertices) == len(geom.mesh.vertices) + 1
    convex = {"cube", "sphere", "capsule"}
    if geom.name in convex:
        assert tet.edge_on_facets == 0, f"{geom.name}: convex body, no edge-on facets"
    surf_vol = surface_signed_volume(geom.mesh)
    assert surf_vol > 0.0, f"{geom.name}: surface must be outward oriented"
    assert tet.signed_volume() == pytest.approx(surf_vol, rel=1e-12), (
        f"{geom.name}: tet volume {tet.signed_volume()} vs surface {surf_vol}"
    )


@pytest.mark.parametrize("geom", GEOMETRIES, ids=lambda g: g.name)
def test_tetrahedra_all_positive(geom: GeometryCase) -> None:
    """Every tetrahedron must have strictly positive volume (no inverted cells)."""
    tet = geom.tetra
    volumes = [tetrahedral_volume(tet.vertices, t) for t in tet.tets]
    assert min(volumes) > 0.0, f"{geom.name}: min tet volume {min(volumes)}"


def test_tetrahedralise_rejects_open_surface() -> None:
    """A surface with a hole must be refused, not silently filled."""
    cube = _make_cube()
    open_mesh = Mesh(vertices=cube.vertices, faces=cube.faces[:-1], units=cube.units)
    with pytest.raises(TetrahedralisationError, match="not closed"):
        tetrahedralise_surface(open_mesh)


def test_tetrahedralise_flips_inward_surface() -> None:
    """An inward-oriented surface is flipped globally, volume stays positive."""
    cube = _make_cube()
    inward = Mesh(vertices=cube.vertices, faces=tuple((a, c, b) for a, b, c in cube.faces))
    assert surface_signed_volume(inward) < 0.0
    tet = tetrahedralise_surface(inward)
    assert tet.signed_volume() == pytest.approx(1.0, rel=1e-12)


def test_tetrahedralise_refuses_non_star_shaped_body() -> None:
    """A body whose centroid cannot see the surface must raise, not approximate."""
    # Two unit cubes joined into an L with the apex per-vertex mean inside one lobe
    # is hard to build; instead force an exterior apex, which sees some faces back-face-on.
    cube = _make_cube()
    with pytest.raises(TetrahedralisationError, match="not star-shaped"):
        tetrahedralise_surface(cube, interior_point=(2.0, 2.0, 2.0))


def test_geometry_gate_discrete_volume_vs_analytic() -> None:
    """Geometry gate: every fixture must be watertight and within 2 % of analytic.

    This is the check the old matrix was missing -- it verified that a mesh
    round-tripped to itself, never that the mesh was the body it claimed to be
    (the old sphere was open and inward-facing, the old Wigley reference was
    non-manifold and 2x off its own volume).
    """
    for geom in GEOMETRIES:
        assert meshlib.is_watertight(geom.mesh), f"{geom.name}: surface is not closed"
        assert meshlib.is_edge_manifold(geom.mesh), f"{geom.name}: non-manifold edges"
        measured = surface_signed_volume(geom.mesh)
        analytic = ANALYTIC_VOLUME[geom.name]
        assert measured > 0.0, f"{geom.name}: inward-oriented surface"
        rel = abs(measured - analytic) / analytic
        assert rel < 2e-2, (
            f"{geom.name}: discrete volume {measured:.6f} vs analytic {analytic:.6f} "
            f"(rel {rel:.2%}, limit 2%)"
        )


# ---------------------------------------------------------------------------
# The 18-combination matrix
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("geom", GEOMETRIES, ids=lambda g: g.name)
@pytest.mark.parametrize("fmt", FORMATS_ALL)
def test_round_trip_volume_preserved(fmt: str, geom: GeometryCase) -> None:
    """Volume must be preserved for all 6 formats × 3 geometries."""
    original = _matrix_mesh(geom, fmt)
    rt = round_trip_meshio(original, file_format=fmt)
    vol_orig = _volume(original)
    vol_rt = _volume(rt)
    rel_err = _relative_error(vol_rt, vol_orig)
    limit = 1e-3 if fmt == "stl" else 1e-12
    assert rel_err < limit, (
        f"{fmt}/{geom.name}: volume {vol_rt} vs {vol_orig} "
        f"(rel_err={rel_err:.2e}, limit={limit:g})"
    )


@pytest.mark.parametrize("geom", GEOMETRIES, ids=lambda g: g.name)
@pytest.mark.parametrize("fmt", FORMATS_ALL)
def test_round_trip_topology_preserved(fmt: str, geom: GeometryCase) -> None:
    """Vertex count and cell count must be preserved exactly."""
    original = _matrix_mesh(geom, fmt)
    rt = round_trip_meshio(original, file_format=fmt)
    assert len(rt.vertices) == len(original.vertices), (
        f"{fmt}/{geom.name}: vertex count {len(rt.vertices)} vs {len(original.vertices)}"
    )
    assert _cell_count(rt) == _cell_count(original), (
        f"{fmt}/{geom.name}: cell count {_cell_count(rt)} vs {_cell_count(original)}"
    )
    assert isinstance(rt, TetrahedralMesh) == isinstance(original, TetrahedralMesh)


@pytest.mark.parametrize("geom", GEOMETRIES, ids=lambda g: g.name)
@pytest.mark.parametrize("fmt", FORMATS_ALL)
def test_round_trip_orientation_preserved(fmt: str, geom: GeometryCase) -> None:
    """Signed volume sign must be unchanged (outward normals / positive tets)."""
    original = _matrix_mesh(geom, fmt)
    rt = round_trip_meshio(original, file_format=fmt)
    vol_orig = _volume(original)
    vol_rt = _volume(rt)
    assert math.copysign(1, vol_rt) == math.copysign(1, vol_orig), (
        f"{fmt}/{geom.name}: orientation flipped "
        f"(vol_orig={vol_orig}, vol_rt={vol_rt})"
    )


@pytest.mark.parametrize("geom", GEOMETRIES, ids=lambda g: g.name)
@pytest.mark.parametrize("fmt", FORMATS_ALL)
def test_exact_topology_for_lossless_formats(fmt: str, geom: GeometryCase) -> None:
    """OBJ/VTK/Gmsh/SU2/CGNS must preserve exact positions and cell indices."""
    if fmt == "stl":
        pytest.skip("STL is a triangle soup: it reorders vertices and faces by design")
    original = _matrix_mesh(geom, fmt)
    rt = round_trip_meshio(original, file_format=fmt)
    for i, (v_orig, v_rt) in enumerate(zip(original.vertices, rt.vertices, strict=True)):
        assert v_orig == pytest.approx(v_rt, abs=1e-12), (
            f"{fmt}/{geom.name}: vertex {i} mismatch {v_orig} vs {v_rt}"
        )
    orig_cells = original.tets if isinstance(original, TetrahedralMesh) else original.faces
    rt_cells = rt.tets if isinstance(rt, TetrahedralMesh) else rt.faces
    for i, (c_orig, c_rt) in enumerate(zip(orig_cells, rt_cells, strict=True)):
        assert c_orig == c_rt, f"{fmt}/{geom.name}: cell {i} mismatch {c_orig} vs {c_rt}"


# ---------------------------------------------------------------------------
# Topology / format mismatch handling
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("geom", GEOMETRIES, ids=lambda g: g.name)
@pytest.mark.parametrize("fmt", FORMATS_VOLUME)
def test_volume_formats_reject_surface_mesh(fmt: str, geom: GeometryCase) -> None:
    """SU2/CGNS must refuse a surface mesh and point at tetrahedralise_surface."""
    with pytest.raises(MeshioFormatError, match="volume cells") as excinfo:
        round_trip_meshio(geom.mesh, file_format=fmt)
    assert "tetrahedralise_surface" in str(excinfo.value)


@pytest.mark.parametrize("geom", GEOMETRIES, ids=lambda g: g.name)
@pytest.mark.parametrize("fmt", FORMATS_SURFACE)
def test_surface_formats_reject_tetrahedral_mesh(fmt: str, geom: GeometryCase) -> None:
    """Surface formats cannot carry tetrahedra; the mismatch must be explicit."""
    with pytest.raises(MeshioFormatError, match="surface format"):
        round_trip_meshio(geom.tetra, file_format=fmt)


def test_unknown_format_rejected(tmp_path: Path) -> None:
    with pytest.raises(MeshioFormatError, match="unknown format"):
        round_trip_meshio(_make_cube(), file_format="exodus")


def test_write_meshio_su2_file_is_volume_mesh(tmp_path: Path) -> None:
    """The SU2 file we emit must be NDIME=3 with tetra (type 10) elements."""
    cube = _make_cube()
    path = write_meshio(tmp_path / "cube.su2", tetrahedralise_surface(cube), file_format="su2")
    lines = path.read_text().splitlines()
    assert "NDIME= 3" in lines
    nel_index = next(i for i, ln in enumerate(lines) if ln.startswith("NELEM="))
    nel = int(lines[nel_index].split("=")[1])
    assert nel == len(cube.faces)
    elements = [ln.split() for ln in lines[nel_index + 1:] if len(ln.split()) == 5]
    assert len(elements) == nel
    types = {int(fields[0]) for fields in elements}
    assert types == {10}, f"SU2 element types written: {types} (10 = tetra)"


def test_write_meshio_cgns_file_is_volume_mesh(tmp_path: Path) -> None:
    """The CGNS file we emit must carry tetra connectivity and an element range.

    Ground truth is read straight out of the HDF5 file, not through our own
    reader: `ElementConnectivity` must hold 4*n 1-based indices and
    `ElementRange` must be [1, n].
    """
    h5py = pytest.importorskip("h5py")
    cube = _make_cube()
    tet = tetrahedralise_surface(cube)
    path = write_meshio(tmp_path / "cube.cgns", tet, file_format="cgns")
    with h5py.File(path, "r") as handle:
        elements = handle["Base"]["Zone1"]["GridElements"]
        connectivity = np.asarray(elements["ElementConnectivity"][" data"])
        element_range = np.asarray(elements["ElementRange"][" data"])
        coords = {
            axis: np.asarray(handle["Base"]["Zone1"]["GridCoordinates"][f"Coordinate{axis}"][" data"])
            for axis in "XYZ"
        }
    assert connectivity.shape == (4 * len(tet.tets),)
    assert connectivity.min() == 1 and connectivity.max() == len(tet.vertices)
    assert element_range.tolist() == [1, len(tet.tets)]
    for axis_index, axis in enumerate("XYZ"):
        np.testing.assert_allclose(coords[axis], np.array(tet.vertices)[:, axis_index])


# ---------------------------------------------------------------------------
# meshio conversion functions
# ---------------------------------------------------------------------------

def test_mesh_to_meshio_preserves_data() -> None:
    mesh = _make_cube()
    mio = mesh_to_meshio(mesh)
    assert mio.points.shape == (8, 3)
    assert len(mio.cells) == 1
    assert mio.cells[0].type == "triangle"
    assert mio.cells[0].data.shape == (12, 3)


def test_meshio_to_mesh_roundtrip() -> None:
    mesh = _make_cube()
    mio = mesh_to_meshio(mesh)
    recovered = meshio_to_mesh(mio, units=mesh.units, coordinate_frame=mesh.coordinate_frame)
    assert recovered.vertices == mesh.vertices
    assert recovered.faces == mesh.faces
    assert recovered.units == mesh.units


def test_mesh_to_meshio_tetra() -> None:
    tet = tetrahedralise_surface(_make_cube())
    mio = mesh_to_meshio(tet)
    assert mio.points.shape == (len(tet.vertices), 3)
    assert mio.cells[0].type == "tetra"
    assert mio.cells[0].data.shape == (12, 4)


def test_meshio_to_tetrahedral_rejects_surface_only() -> None:
    mio = mesh_to_meshio(_make_cube())
    with pytest.raises(MeshioFormatError, match="no tetra cells"):
        meshio_to_tetrahedral(mio)


# ---------------------------------------------------------------------------
# Real hull: the case that falsified the "ship hulls work in practice" claim
# ---------------------------------------------------------------------------

def test_hull_tetrahedralisation_is_certified() -> None:
    """The Wigley hull tetrahedralises, with the exact number of edge-on facets.

    Regression guard for the failure that motivated this work: the hull's
    stem/post tip facets lie in the centreplane, so the strict-positivity rule
    raised ``TetrahedralisationError`` ("best candidate ... still gives a
    worst-case tetrahedron volume 0.000e+00").  Two facets are seen edge-on from
    a centreplane apex and exactly two are dropped; if the generator or the
    tessellation changes, this count changes with it.
    """
    hull = _make_wigley_hull()
    tet = tetrahedralise_surface(hull)
    assert tet.edge_on_facets == 2, f"edge-on facets: {tet.edge_on_facets} (expected 2)"
    assert len(tet.tets) == len(hull.faces) - 2
    surf_vol = surface_signed_volume(hull)
    assert tet.signed_volume() == pytest.approx(surf_vol, rel=1e-12)


def test_hull_dropped_facets_lie_in_the_centreplane() -> None:
    """The dropped cells are exactly the facets whose three vertices have y == 0.

    This binds the "edge-on" claim to geometry instead of trusting the count:
    the dropped facets must be the centreplane facets, and every centreplane
    facet must be dropped (a zero-area cell cannot be tiled by a finite apex).
    """
    hull = _make_wigley_hull()
    verts = np.array(hull.vertices)
    centreplane = {
        index
        for index, face in enumerate(hull.faces)
        if all(verts[v][1] == 0.0 for v in face)
    }
    assert len(centreplane) == 2, f"{len(centreplane)} centreplane facets (expected 2)"
    tet = tetrahedralise_surface(hull)
    kept = {tuple(sorted(face)) for *_apex, a, b, c in tet.tets for face in [(a, b, c)]}
    dropped = [
        face for index, face in enumerate(hull.faces) if tuple(sorted(face)) not in kept
    ]
    assert len(dropped) == 2
    assert all(index in centreplane for index, face in enumerate(hull.faces)
               if tuple(sorted(face)) not in kept)


def test_hull_off_centreplane_apex_is_refused() -> None:
    """Off the centreplane the hull really is not star-shaped -- refuse, do not fudge.

    Documented counter-example to any "just nudge the apex" fix: lifting the
    apex 5 % of the beam off the centreplane inverts three of the tip cells, so
    the certificate must fail and say so.
    """
    hull = _make_wigley_hull()
    bbox_centre = (np.array(hull.vertices).min(axis=0) + np.array(hull.vertices).max(axis=0)) / 2
    with pytest.raises(TetrahedralisationError, match="not star-shaped"):
        tetrahedralise_surface(hull, interior_point=(bbox_centre[0], 0.05, bbox_centre[2]))


def test_hull_volume_mesh_boundary_slit_is_quantified() -> None:
    """Known limitation, pinned: the hull volume mesh boundary is open along 2 slits.

    Dropping the zero-volume cells keeps the *volume* identity exact, but the
    two edge-on facets are skin patches of the hull, so the boundary of the
    delivered volume mesh is the hull surface minus those two patches: 2 skin
    facets missing and 6 rim edges unpaired (each slit is bounded by three
    edges).  This is stated rather than hidden -- an SU2/CGNS solver receiving
    the mesh sees an open wall there, which is why *volume-conservation* is REAL
    for the hull while *solver-ready watertight volume mesh* is WITHHELD.
    If a real volume mesher is wired in later, this test is the one to invert.
    """
    hull = _make_wigley_hull()
    tet = tetrahedralise_surface(hull)
    apex = len(tet.vertices) - 1
    counts: dict[tuple[int, ...], int] = {}
    for a, b, c, d in tet.tets:
        for face in ((a, b, c), (b, a, d), (c, b, d), (a, c, d)):
            key = tuple(sorted(face))
            counts[key] = counts.get(key, 0) + 1
    boundary = [face for face, n in counts.items() if n == 1]
    outer = [face for face in boundary if apex not in face]
    assert len(outer) == len(hull.faces) - 2, (
        f"{len(outer)} outer boundary faces vs {len(hull.faces) - 2} expected"
    )
    rim: dict[tuple[int, int], int] = {}
    for face in outer:
        for edge in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            key = (min(edge), max(edge))
            rim[key] = rim.get(key, 0) + 1
    open_edges = [edge for edge, n in rim.items() if n != 2]
    assert len(open_edges) == 6, f"{len(open_edges)} unpaired rim edges (expected 6)"
    # the slit area is small but not zero -- quote it, do not round it away
    verts = np.array(hull.vertices)
    dropped = [
        face for face in hull.faces if tuple(sorted(face)) not in {tuple(f) for f in outer}
    ]
    slit_area = sum(
        0.5 * float(np.linalg.norm(np.cross(verts[b] - verts[a], verts[c] - verts[a])))
        for a, b, c in dropped
    )
    assert 1e-3 < slit_area < 5e-2, f"slit area {slit_area:.5f} m^2 outside expected band"


def _make_torus(major_radius: float = 2.0, minor_radius: float = 0.5,
                n_major: int = 64, n_minor: int = 32) -> Mesh:
    """Closed torus: *not* star-shaped from any point (the certification must refuse it)."""
    builder = MeshBuilder(units="m")
    rings: list[list[int]] = []
    for i in range(n_major):
        u = 2.0 * math.pi * i / n_major
        ring = []
        for j in range(n_minor):
            v = 2.0 * math.pi * j / n_minor
            ring.append(
                builder.add_vertex(
                    (
                        (major_radius + minor_radius * math.cos(v)) * math.cos(u),
                        (major_radius + minor_radius * math.cos(v)) * math.sin(u),
                        minor_radius * math.sin(v),
                    )
                )
            )
        rings.append(ring)
    for i in range(n_major):
        a, b = rings[i], rings[(i + 1) % n_major]
        for j in range(n_minor):
            k = (j + 1) % n_minor
            builder.add_quad(a[j], b[j], b[k], a[k])
    mesh = builder.build()
    return meshlib.orient_outward(mesh).mesh if surface_signed_volume(mesh) < 0 else mesh


def test_certification_refuses_non_star_body_torus() -> None:
    """The volume-identity certificate must not become a rubber stamp.

    A torus is the canonical non-star body: from *every* interior point some
    tetrahedra are inverted (the hole is behind the apex), so the certificate
    fails and the error names the inverted cells.  This guards the relaxation
    that made the hull work -- zero-volume cells are tolerated only together
    with a volume identity that a torus cannot satisfy.
    """
    torus = _make_torus()
    assert meshlib.is_watertight(torus)
    analytic = 2.0 * math.pi**2 * 2.0 * 0.5**2
    assert abs(surface_signed_volume(torus) - analytic) / analytic < 0.02
    with pytest.raises(TetrahedralisationError, match="not star-shaped"):
        tetrahedralise_surface(torus)


def test_hull_hull_form_gate_volume_vs_analytic() -> None:
    """Displaced volume from the hull's triangles vs the closed form (gate: < 1 %).

    The 41 x 13 station/level grid under-resolves the ends, so the sampled solid
    is slightly smaller than the analytic ``4LBT/9``; the gate is the geometry
    validity check the matrix needs before any round-trip claim is meaningful.
    """
    hull = _make_wigley_hull()
    measured = surface_signed_volume(hull)
    analytic = ANALYTIC_VOLUME["wigley_hull"]
    rel = abs(measured - analytic) / analytic
    assert rel < 1e-2, f"hull volume {measured:.6f} vs analytic {analytic:.6f} ({rel:.2%})"


# ---------------------------------------------------------------------------
# Summary table generation (for artifacts)
# ---------------------------------------------------------------------------

def generate_format_matrix_table() -> str:
    """Generate a markdown table of all 24 combinations, measured not asserted."""
    lines = [
        "# meshio Format Matrix — Round-Trip Fidelity (6 formats × 4 geometries)",
        "",
        "| Format | Geometry | Topology sent | Volume orig | Volume round-trip | rel. err | "
        "Topology | Orientation | Status |",
        "|--------|----------|---------------|-------------|-------------------|----------|"
        "----------|-------------|--------|",
    ]
    for fmt in FORMATS_ALL:
        for geom in GEOMETRIES:
            original = _matrix_mesh(geom, fmt)
            rt = round_trip_meshio(original, file_format=fmt)
            vol_orig = _volume(original)
            vol_rt = _volume(rt)
            rel_err = _relative_error(vol_rt, vol_orig)
            vol_ok = rel_err < (1e-3 if fmt == "stl" else 1e-12)
            topo_ok = len(rt.vertices) == len(original.vertices) and _cell_count(rt) == _cell_count(original)
            orient_ok = math.copysign(1, vol_rt) == math.copysign(1, vol_orig)
            status = "PASS" if (vol_ok and topo_ok and orient_ok) else "FAIL"
            kind = "tetra" if isinstance(original, TetrahedralMesh) else "triangle"
            lines.append(
                f"| {fmt} | {geom.name} | {kind} "
                f"({len(original.vertices)}v/{_cell_count(original)}c) | "
                f"{vol_orig:.9f} | {vol_rt:.9f} | {rel_err:.2e} | "
                f"{'ok' if topo_ok else 'MISMATCH'} | {'ok' if orient_ok else 'FLIPPED'} | "
                f"{status} |"
            )

    lines += [
        "",
        "## Notes",
        "- STL reorders vertices and faces (triangle soup) but preserves volume and orientation.",
        "- `wigley_hull` is `tensorcad.wigley`'s deck-capped benchmark hull (L=10, B=1, "
        "T=0.625); its volume reference is the closed form 4LBT/9 = 25/9 m^3.",
        "- For the hull the SU2/CGNS rows carry 1996 tets: 2 surface facets lie exactly in "
        "the centreplane and are seen edge-on by the apex, so their zero-volume tets are "
        "dropped (`edge_on_facets=2`). Dropping zero-measure cells keeps the tiling exact; "
        "`test_hull_tetrahedralisation_is_certified` pins the volume identity.",
        "- OBJ, VTK, Gmsh preserve exact vertex positions and triangle indices.",
        "- SU2 and CGNS carry **volume** cells (tetra); the matrix feeds them "
        "`tetrahedralise_surface(mesh)`, whose tets tile the body exactly, so the volume "
        "invariant is the same divergence-theorem number as for the surface formats.",
        "- CGNS caveat: meshio's writer stores a single uniform tetra block and omits full "
        "CGNS/HDF5 metadata; round-trip through meshio is REAL, third-party CGNS reader "
        "conformance is WITHHELD (no cgnslib on this machine).",
        "",
    ]
    return "\n".join(lines)


if __name__ == "__main__":  # pragma: no cover
    print(generate_format_matrix_table())
