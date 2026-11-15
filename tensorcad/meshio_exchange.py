"""meshio-based format exchange for TensorCAD meshes.

Provides round-trip conversion between TensorCAD geometry (surface triangles or
tetrahedra) and meshio for the six formats required by the format-matrix audit:
STL, OBJ, VTK, SU2, CGNS, Gmsh.

Why topology, not just "format", decides the matrix
===================================================
The earlier version of this module declared SU2 and CGNS *NOT_SUPPORTED* and
stopped there, with a comment that "their meshio writers only emit volume cell
types".  Half of that is true and half of it hid a real capability: SU2 and CGNS
are **volume-cell formats** -- their native element is the tetrahedron -- so a
surface triangle mesh cannot be expressed in them, but the *geometry* can, once
it is tetrahedralised.  :func:`tetrahedralise_surface` does that (star
tessellation from an interior point, refused for non star-shaped inputs), after
which both formats round-trip with full fidelity.  Claims in the old docstring
about "tetrahedralise the surface first" were unbacked: nothing in the package
actually offered it and no test exercised the volume path.

Format capability matrix
========================
======  ===================  ==============================  ============
Format  surface triangles    tetrahedra                      verification
======  ===================  ==============================  ============
STL     REAL                 -                               round trip
OBJ     REAL                 -                               round trip
VTK     REAL                 -                               round trip
Gmsh    REAL                 -                               round trip
SU2     NOT_SUPPORTED (*)    REAL                            round trip
CGNS    NOT_SUPPORTED (*)    REAL                            round trip
======  ===================  ==============================  ============

The matrix is about meshio's surface/volume split, and it stays as it is: ``gmsh`` is
listed as a *surface* format because that is what the 18-combination gate measures through
meshio.  A body grid still reaches gmsh -- through
:func:`write_gmsh_volume_mesh`/:func:`read_gmsh_volume_groups`, which write msh 2.2 ASCII
with one named physical group per patch (the same patches :func:`volume_mesh_patches`
derives for SU2).  That is a separate entry point rather than a new member of the matrix
because reclassifying a format would silently change what the existing gate measures.

Import direction (files TensorCAD did not write)
------------------------------------------------
The table above is about what TensorCAD can *emit*.  Reading is a different
problem: files arriving from a CAD kernel, Salome/Pointwise or a structured
panel generator are usually **quad** or **polygon** cages, often in several
cell blocks (several Gmsh volume regions, several CGNS zones), and meshes
written in one block often come back split in two.  :func:`meshio_to_mesh`
therefore consumes *every* surface block and triangulates quads (diagonal
0-2) and polygons (fan from vertex 0), preserving winding, and
:func:`meshio_to_tetrahedral` concatenates *every* tetra block instead of
taking the first.  Before that, a six-quad cube -- the simplest closed body
there is -- was rejected as "contains no triangle cells", and a file whose
first triangle block was a boundary patch silently lost the rest of its
geometry.  ``tests/test_meshio_import_cells.py`` gates both directions.

Two producers of the tetrahedra
--------------------------------
:func:`tetrahedralise_surface` star-tessellates a *closed surface*; it is exact for
star-shaped bodies and refuses the rest.  :func:`volume_mesh_to_tetrahedral` takes the
body-fitted grid :mod:`tensorcad.volume_mesh` builds (a real hull, boundary layer grading
and all) and hands over its Kuhn-split tetrahedra.  That second path is the one that works
for the package's own hull generator, whose rings collapse to zero breadth at the stem and
post and are therefore refused by star tessellation.

(*) SU2/CGNS surface triangles are not supported *by meshio* in this environment,
but they are **not** unsupported by the formats: SU2 carries triangles as element
type 5 and CGNS carries ``TRI_3``/``QUAD_4``, and :mod:`tensorcad.mesh_formats`
writes and reads both.  The limitation is upstream and reproducible (meshio 5.3.5):
its SU2 writer raises ``TypeError: cannot unpack non-iterable CellBlock object``
for surface cells, and its CGNS reader raises ``KeyError: object ' data' doesn't
exist`` on its own triangle output -- tetrahedra work through both.  So
:func:`round_trip_meshio` raises a ``MeshioFormatError`` pointing the caller at
:func:`tetrahedralise_surface` (volume path) or at
``tensorcad.mesh_formats.write``/``read`` (surface path) instead of silently
having nothing to offer.  ``tests/test_mesh_formats_foreign_intake.py`` pins the
surface path; this docstring previously claimed the formats themselves could not
carry a surface, which was wrong.

A caveat that is stated rather than hidden: the CGNS path goes through meshio's
CGNS writer, which stores a *single* uniform element block (tetra only) and does
not emit the full CGNS/HDF5 metadata (CGNSLibraryVersion, ElementType, ...).
Round-trip fidelity through meshio is REAL; conformance as a CGNS file for a
third-party CGNS reader is **WITHHELD** -- no ``cgnslib`` is available on this
machine to check it independently.
"""
from __future__ import annotations

import math
import tempfile
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from . import mesh as meshlib
from .geometry import Mesh

__all__ = [
    "FORMATS_ALL",
    "FORMATS_SURFACE",
    "FORMATS_VOLUME",
    "GMSH_SOLID_GROUP",
    "GMSH_VOLUME_VERSION",
    "PATCH_NAMES",
    "SURFACE_CELL_TYPES",
    "VOLUME_CELL_TYPES",
    "MeshioFormatError",
    "TetrahedralMesh",
    "TetrahedralisationError",
    "VolumePatches",
    "mesh_to_meshio",
    "meshio_to_mesh",
    "meshio_to_tetrahedral",
    "read_gmsh_volume_groups",
    "read_meshio",
    "read_su2_markers",
    "round_trip_meshio",
    "surface_signed_volume",
    "tetrahedral_volume",
    "tetrahedralise_surface",
    "volume_mesh_patches",
    "volume_mesh_to_tetrahedral",
    "write_gmsh_volume_mesh",
    "write_meshio",
    "write_su2_volume_mesh",
]

Vec3 = tuple[float, float, float]
Tet = tuple[int, int, int, int]

#: Formats that support surface triangle meshes.
FORMATS_SURFACE: tuple[str, ...] = ("stl", "obj", "vtk", "gmsh")

#: Formats that only support volume cells (require tetrahedralisation first).
FORMATS_VOLUME: tuple[str, ...] = ("su2", "cgns")

#: Every format in the matrix.
FORMATS_ALL: tuple[str, ...] = FORMATS_SURFACE + FORMATS_VOLUME

#: Extension map for each format.
_EXT: dict[str, str] = {
    "stl": ".stl",
    "obj": ".obj",
    "vtk": ".vtk",
    "su2": ".su2",
    "cgns": ".cgns",
    "gmsh": ".msh",
}


class MeshioFormatError(ValueError):
    """Raised when a format cannot handle the given mesh type."""


class TetrahedralisationError(ValueError):
    """Raised when a surface mesh cannot be tetrahedralised by star tessellation."""


# ---------------------------------------------------------------------------
# Geometry helpers (no meshio involved)
# ---------------------------------------------------------------------------

def surface_signed_volume(mesh: Mesh) -> float:
    """Signed volume enclosed by a closed triangle mesh (divergence theorem).

    Positive for outward-oriented (counter-clockwise seen from outside) faces.
    """
    total = 0.0
    verts = mesh.vertices
    for a, b, c in mesh.faces:
        ax, ay, az = verts[a]
        bx, by, bz = verts[b]
        cx, cy, cz = verts[c]
        total += (
            ax * (by * cz - bz * cy)
            - ay * (bx * cz - bz * cx)
            + az * (bx * cy - by * cx)
        )
    return total / 6.0


def _tet_signed_volume_from_points(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3) -> float:
    a = np.subtract(p1, p0)
    b = np.subtract(p2, p0)
    c = np.subtract(p3, p0)
    return float(np.dot(a, np.cross(b, c)) / 6.0)


def tetrahedral_volume(vertices: Sequence[Vec3], tet: Sequence[int]) -> float:
    """Signed volume of one tetrahedron given by four vertex indices."""
    p0, p1, p2, p3 = (vertices[i] for i in tet)
    return _tet_signed_volume_from_points(
        tuple(map(float, p0)), tuple(map(float, p1)), tuple(map(float, p2)), tuple(map(float, p3))
    )


# ---------------------------------------------------------------------------
# Tetrahedral mesh container
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TetrahedralMesh:
    """Linear tetrahedral volume mesh: vertex positions plus 4-tuples.

    ``degenerate_cells`` counts zero-volume tetrahedra that were present in the
    source topology and dropped here: :func:`volume_mesh_to_tetrahedral` hits
    them on every real hull body, where the keel line and the two end stations
    collapse to zero breadth, so the cells of those columns have no measure.
    Dropping them is exact (a zero-volume cell contributes nothing to the
    volume), and the count is carried so it is visible instead of silent.

    ``tets`` holds 0-based indices into ``vertices``.  All tetrahedra are
    strictly positively oriented (positive signed volume); the star
    tessellation produced by :func:`tetrahedralise_surface` guarantees it.

    ``edge_on_facets`` records how many surface facets the tessellation had to
    drop because the star apex sees them exactly edge-on (the tetrahedron
    spanned by the apex and such a facet has zero volume).  This is not an
    approximation: a zero-volume cell carries no measure, so dropping it leaves
    the tiling of the body exact -- :func:`tetrahedralise_surface` only reports
    a non-zero count after it has certified that the remaining tets reproduce
    the surface volume.  Real hull meshes hit this case: where a Wigley-type
    ring collapses to zero breadth, the tip facets lie exactly in the
    centreplane and are seen edge-on from any centreplane apex.
    """

    vertices: tuple[Vec3, ...]
    tets: tuple[Tet, ...]
    units: str = "m"
    coordinate_frame: str = "unspecified right-handed frame"
    edge_on_facets: int = 0
    degenerate_cells: int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "vertices", tuple(tuple(float(c) for c in v) for v in self.vertices))
        object.__setattr__(self, "tets", tuple(tuple(int(i) for i in t) for t in self.tets))
        for index, tet in enumerate(self.tets):
            if len(tet) != 4:
                raise ValueError(f"tetrahedron {index} does not have 4 nodes")
            if len(set(tet)) != 4:
                raise ValueError(f"tetrahedron {index} repeats a node: {tet}")
            for vertex in tet:
                if not 0 <= vertex < len(self.vertices):
                    raise ValueError(
                        f"tetrahedron {index} references out-of-range vertex {vertex}"
                    )

    def __len__(self) -> int:
        return len(self.tets)

    def volume(self) -> float:
        """Total volume (sum of absolute tetrahedron volumes)."""
        return sum(abs(tetrahedral_volume(self.vertices, tet)) for tet in self.tets)

    def signed_volume(self) -> float:
        """Total signed volume; equals :meth:`volume` when all tets are positive."""
        return sum(tetrahedral_volume(self.vertices, tet) for tet in self.tets)


# ---------------------------------------------------------------------------
# Surface -> tetrahedra (star tessellation)
# ---------------------------------------------------------------------------

def _closed_oriented_faces(mesh: Mesh) -> tuple[tuple[int, int, int], ...]:
    """Validate that ``mesh`` is a closed, consistently oriented manifold.

    Two independent checks, both of them *invariants* rather than smoke tests:

    * :func:`tensorcad.mesh.edge_face_counts` -- every undirected edge must be
      shared by exactly two faces (a boundary edge means a hole, three means
      non-manifold);
    * directed winding -- every directed edge must appear exactly once, so that
      neighbouring faces traverse their shared edge in opposite directions and
      the divergence-theorem volume is meaningful.
    """
    counts = meshlib.edge_face_counts(mesh)
    bad = [edge for edge, count in counts.items() if count != 2]
    if bad:
        raise TetrahedralisationError(
            f"surface is not closed: {len(bad)} of {len(counts)} undirected edges are "
            f"shared by {counts[bad[0]]} face(s) instead of 2 (first edge {bad[0]})"
        )
    directed: dict[tuple[int, int], int] = {}
    for a, b, c in mesh.faces:
        for edge in ((a, b), (b, c), (c, a)):
            directed[edge] = directed.get(edge, 0) + 1
    duplicated = [e for e, count in directed.items() if count > 1]
    if duplicated:
        raise TetrahedralisationError(
            "surface is not consistently oriented: directed edge "
            f"{duplicated[0]} is traversed {directed[duplicated[0]]} times"
        )
    return tuple(mesh.faces)


def _bbox(mesh: Mesh) -> tuple[Vec3, Vec3]:
    verts = np.array(mesh.vertices, dtype=np.float64)
    return tuple(verts.min(axis=0)), tuple(verts.max(axis=0))


def _interior_candidates(mesh: Mesh) -> list[Vec3]:
    """Deterministic candidate interior points, most conservative first.

    The vertex mean is a poor apex on a ship hull: the deck cap and the hull
    sides sit in the *same* planes as the mean (this repo's hull mesh has the
    mean at z = -2e-18, i.e. exactly in the deck plane), which makes some
    tetrahedra exactly degenerate.  So a small set of candidates is offered and
    the one with the best worst-case tetrahedron wins.
    """
    verts = np.array(mesh.vertices, dtype=np.float64)
    lo, hi = verts.min(axis=0), verts.max(axis=0)
    extent = hi - lo
    mean = verts.mean(axis=0)
    centre = (lo + hi) / 2.0

    areas = np.array([meshlib.triangle_area(mesh, f) for f in mesh.faces])
    face_centroids = np.array(
        [verts[list(f)].mean(axis=0) for f in mesh.faces], dtype=np.float64
    )
    total_area = float(areas.sum())
    weighted = (
        (face_centroids * areas[:, None]).sum(axis=0) / total_area
        if total_area > 0.0
        else mean
    )

    seeds = [weighted, mean, centre]
    candidates = [tuple(float(c) for c in seed) for seed in seeds]
    # Nudge off any plane the seeds happen to lie in (deck, centreplane, ends).
    for seed in seeds:
        for axis in range(3):
            step = 0.05 * float(extent[axis])
            for sign in (+1.0, -1.0):
                point = np.array(seed, dtype=np.float64).copy()
                point[axis] += sign * step
                candidates.append(tuple(float(c) for c in point))
    return candidates


def tetrahedralise_surface(mesh: Mesh, *, interior_point: Vec3 | None = None) -> TetrahedralMesh:
    """Fill a closed surface mesh with tetrahedra (star tessellation).

    A single interior vertex is inserted and every surface triangle is joined to
    it, giving ``len(mesh.faces)`` tetrahedra.  This is exact (no approximation
    of the geometry) and reverses cleanly: the tets' boundary *is* the input
    surface, so their total volume equals the surface's signed volume.

    Valid for bodies that are star-shaped with respect to the interior point --
    which covers convex bodies (cube, sphere) and, in practice, ship hulls whose
    centroid sees the whole shell.  A non star-shaped input is **refused** with
    :class:`TetrahedralisationError` rather than approximated.

    When ``interior_point`` is None an automatic search runs over a small
    deterministic candidate set (area-weighted surface centroid, vertex mean,
    bounding-box centre, and +/- 5 % nudges along each axis).  Acceptance is a
    *certificate*, not a heuristic: a candidate is accepted when (a) no
    tetrahedron is inverted and (b) the tets it spans reproduce the surface
    volume to ``1e-9`` relative.  Among accepted candidates the one that leaves
    the fewest edge-on facets wins.

    Facets the apex sees exactly edge-on span zero-volume tetrahedra.  Those
    cells are dropped and counted in ``TetrahedralMesh.edge_on_facets``; the
    certificate above is what makes the drop exact rather than a silent loss.
    This matters for real hulls, not for toy bodies: where a Wigley-type section
    ring collapses to zero breadth (see :func:`tensorcad.wigley.build_mesh`) the
    resulting tip facets lie exactly in the centreplane, so *every* off-plane
    apex inverts them while an on-plane apex sees them edge-on.  Requiring
    strictly positive cells (the previous rule) refused such a hull outright and
    falsely reported it as non star-shaped.

    Parameters
    ----------
    mesh :
        Closed, consistently oriented triangle mesh.
    interior_point :
        Explicit star centre.  When given, it is used as-is and the result is
        validated -- no silent fallback to another point.
    """
    faces = _closed_oriented_faces(mesh)
    verts = [tuple(float(c) for c in v) for v in mesh.vertices]
    volume = surface_signed_volume(mesh)
    if volume < 0.0:
        # Inward-oriented input: flipping the whole surface is exact.
        faces = tuple((a, c, b) for a, b, c in faces)
        volume = -volume
    if volume <= 0.0:
        raise TetrahedralisationError("surface encloses no volume (degenerate mesh)")

    lo, hi = _bbox(mesh)
    extent = np.subtract(hi, lo)
    scale = float(np.linalg.norm(extent))
    # Volume below this is numerically coplanar: such a cell is not emitted.
    tol = 1e-12 * (scale**3 if scale > 0.0 else 1.0)
    # A tessellation is certified exact when the emitted tets reproduce the
    # surface volume to this relative tolerance.  Together with "no negative
    # tet" this is what distinguishes a genuinely star-shaped body that has
    # edge-on facets (certified) from a non star-shaped one (refused).
    rel_volume_tol = 1e-9

    faces_array = np.array(faces, dtype=np.int64)
    verts_array = np.array(verts, dtype=np.float64)

    def tet_volumes(apex: Vec3) -> np.ndarray:
        """Signed volumes of (apex, face) tetrahedra, vectorised over faces."""
        origin = np.array(apex, dtype=np.float64)
        pa = verts_array[faces_array[:, 0]] - origin
        pb = verts_array[faces_array[:, 1]] - origin
        pc = verts_array[faces_array[:, 2]] - origin
        triple = np.einsum("ij,ij->i", pa, np.cross(pb, pc))
        return triple / 6.0

    def assess(apex: Vec3) -> tuple[bool, float, int, int, float]:
        """(certified, worst volume, n_negative, n_edge_on, relative volume error)."""
        volumes = tet_volumes(apex)
        n_negative = int(np.count_nonzero(volumes < -tol))
        n_edge_on = int(np.count_nonzero(np.abs(volumes) <= tol))
        rel_error = abs(float(volumes.sum()) - volume) / volume
        certified = n_negative == 0 and rel_error <= rel_volume_tol
        return certified, float(volumes.min()), n_negative, n_edge_on, rel_error

    if interior_point is not None:
        apex = tuple(float(c) for c in interior_point)
        certified, worst, n_negative, n_edge_on, rel_error = assess(apex)
        if not certified:
            raise TetrahedralisationError(
                f"the body is not star-shaped with respect to {apex!r}: "
                f"{n_negative} of {len(faces)} tetrahedra are inverted (volume < -{tol:.3e}), "
                f"worst volume {worst:.3e}, and the tiling reproduces the surface volume "
                f"only to {rel_error:.3e} (needs <= {rel_volume_tol:.3e}). Pass a different "
                f"interior_point= or use a volume mesher."
            )
    else:
        best: tuple[tuple[int, float], Vec3, int] | None = None
        rejected: tuple[float, Vec3, int, float] | None = None
        for candidate in _interior_candidates(mesh):
            certified, worst, n_negative, n_edge_on, rel_error = assess(candidate)
            if certified:
                # Prefer strictly positive tets (fewer edge-on facets), then the
                # largest worst-case cell.
                key = (n_edge_on, -worst)
                if best is None or key < best[0]:
                    best = (key, candidate, n_edge_on)
            elif rejected is None or n_negative < rejected[2]:
                rejected = (worst, candidate, n_negative, rel_error)
        if best is None:
            assert rejected is not None
            worst, apex_bad, n_negative, rel_error = rejected
            raise TetrahedralisationError(
                f"no interior point sees the whole surface: best candidate {apex_bad!r} still "
                f"leaves {n_negative} inverted tetrahedra (worst volume {worst:.3e}) and "
                f"reproduces the surface volume only to {rel_error:.3e} "
                f"(needs <= {rel_volume_tol:.3e}). The body is not star-shaped; use a volume "
                f"mesher instead."
            )
        _, apex, n_edge_on = best

    apex_index = len(verts)
    verts.append(apex)
    volumes = tet_volumes(apex)
    kept = tuple(
        (apex_index, a, b, c)
        for (a, b, c), cell_volume in zip(faces, volumes, strict=True)
        if cell_volume > tol
    )
    return TetrahedralMesh(
        vertices=tuple(verts),
        tets=kept,
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
        edge_on_facets=len(faces) - len(kept),
    )


# ---------------------------------------------------------------------------
# meshio conversion
# ---------------------------------------------------------------------------

def _cell_blocks(cells: Iterable[Any]) -> list[Any]:
    """Normalise ``mesh.cells`` to a list of cell blocks.

    meshio 5 returns ``CellBlock`` objects; older versions returned
    ``(type, data)`` tuples.  Accept both so the module is not pinned to one
    meshio release (5.3.5 is what this repo tests against).
    """
    blocks = []
    for block in cells:
        if hasattr(block, "type") and hasattr(block, "data"):
            blocks.append(block)
        else:  # pragma: no cover - legacy meshio (<5)
            cell_type, data = block
            import meshio

            blocks.append(meshio.CellBlock(cell_type, np.asarray(data)))
    return blocks


def mesh_to_meshio(mesh: Mesh | TetrahedralMesh) -> Any:
    """Convert a surface or tetrahedral TensorCAD mesh to a ``meshio.Mesh``."""
    import meshio

    if isinstance(mesh, TetrahedralMesh):
        points = np.array(mesh.vertices, dtype=np.float64)
        cells = np.array(mesh.tets, dtype=np.int64)
        return meshio.Mesh(points, [("tetra", cells)])

    points = np.array(mesh.vertices, dtype=np.float64)
    cells = np.array(mesh.faces, dtype=np.int64)
    return meshio.Mesh(points, [("triangle", cells)])


#: Cell types that carry a 2-D surface facet, and how many triangles each splits into.
SURFACE_CELL_TYPES: tuple[str, ...] = ("triangle", "quad", "polygon")

#: Cell types that carry a 3-D volume cell.
VOLUME_CELL_TYPES: tuple[str, ...] = ("tetra", "hexahedron", "wedge", "pyramid")


def _triangulate_block(block_type: str, data: Any) -> list[tuple[int, int, int]]:
    """Split one meshio cell block of a surface type into 0-based triangles.

    ``triangle`` passes through; ``quad`` splits along the ``0-2`` diagonal and
    ``polygon`` is fanned from its first vertex.  Both splits preserve the
    winding order of the parent facet, so an outward-oriented quad cage stays
    outward-oriented as triangles (verified by the signed volume in the tests).
    Triangles with a repeated index are dropped: they are topologically invalid
    (zero area, self-edge) and cannot come from a conforming writer.
    """
    if block_type == "triangle" or block_type == "quad" or block_type == "polygon":
        rows = data
    else:  # pragma: no cover - guarded by the caller
        raise MeshioFormatError(f"cell type '{block_type}' is not a surface facet")

    triangles: list[tuple[int, int, int]] = []
    for row in rows:
        idx = [int(i) for i in row]
        if block_type == "triangle":
            if len(idx) != 3:
                raise MeshioFormatError(
                    f"triangle block carries a cell with {len(idx)} indices"
                )
            if len(set(idx)) == 3:
                triangles.append((idx[0], idx[1], idx[2]))
            continue
        if block_type == "quad":
            if len(idx) != 4:
                raise MeshioFormatError(f"quad block carries a cell with {len(idx)} indices")
            for tri in ((idx[0], idx[1], idx[2]), (idx[0], idx[2], idx[3])):
                if len(set(tri)) == 3:
                    triangles.append(tri)
            continue
        # polygon: fan from vertex 0
        if len(idx) < 3:
            continue
        for k in range(1, len(idx) - 1):
            tri = (idx[0], idx[k], idx[k + 1])
            if len(set(tri)) == 3:
                triangles.append(tri)
    return triangles


def _present_cell_types(mio: Any) -> str:
    blocks = _cell_blocks(mio.cells)
    return ", ".join(sorted({block.type for block in blocks})) or "none"


def meshio_to_mesh(
    mio: Any,
    *,
    units: str = "m",
    coordinate_frame: str = "unspecified right-handed frame",
) -> Mesh:
    """Convert a ``meshio.Mesh`` back to a :class:`tensorcad.geometry.Mesh`.

    Every surface cell block is consumed, not just the first one: ``triangle``
    blocks are concatenated and ``quad``/``polygon`` blocks are triangulated
    (:func:`_triangulate_block`, winding preserved).  One-sided and vertex
    blocks (boundary ``line`` sets, isolated points) are not geometry and are
    skipped; 3-D volume blocks are *also* skipped, so importing a volume file
    through this function returns its boundary surface -- use
    :func:`meshio_to_tetrahedral` when the cells themselves are wanted.
    """
    vertices = tuple(tuple(float(c) for c in row) for row in np.asarray(mio.points))
    faces: list[tuple[int, int, int]] = []
    for block in _cell_blocks(mio.cells):
        if block.type in SURFACE_CELL_TYPES:
            faces.extend(_triangulate_block(block.type, np.asarray(block.data)))

    if not faces:
        present = _present_cell_types(mio)
        if any(t in present.split(", ") for t in VOLUME_CELL_TYPES):
            raise MeshioFormatError(
                f"meshio mesh contains no surface facets (cell types present: {present}); "
                "the file is a volume mesh -- use meshio_to_tetrahedral() for its cells"
            )
        raise MeshioFormatError(
            f"meshio mesh contains no surface facets (cell types present: {present})"
        )

    max_index = max(max(face) for face in faces)
    if max_index >= len(vertices):
        raise MeshioFormatError(
            f"cell block references vertex {max_index} but the file only has {len(vertices)} points"
        )
    return Mesh(
        vertices=vertices,
        faces=tuple(faces),
        units=units,
        coordinate_frame=coordinate_frame,
    )


def meshio_to_tetrahedral(
    mio: Any,
    *,
    units: str = "m",
    coordinate_frame: str = "unspecified right-handed frame",
) -> TetrahedralMesh:
    """Convert a ``meshio.Mesh`` carrying tetra cells to a ``TetrahedralMesh``.

    *All* tetra blocks are concatenated (multi-region/multi-zone files: a CGNS
    or Gmsh volume written as several tetra zones used to lose every zone but
    the first, silently).  Non-tetra volume cells (hexahedra, wedges, pyramids)
    are not convertible without splitting, so they are reported as an error
    rather than dropped.
    """
    vertices = tuple(tuple(float(c) for c in row) for row in np.asarray(mio.points))
    tets: list[Tet] = []
    other_volume: set[str] = set()
    for block in _cell_blocks(mio.cells):
        if block.type == "tetra":
            for row in np.asarray(block.data):
                idx = tuple(int(i) for i in row)
                if len(idx) != 4:
                    raise MeshioFormatError(
                        f"tetra block carries a cell with {len(idx)} indices"
                    )
                if len(set(idx)) == 4:
                    tets.append(idx)  # type: ignore[arg-type]
        elif block.type in VOLUME_CELL_TYPES:
            other_volume.add(block.type)

    if not tets:
        present = _present_cell_types(mio)
        detail = f" (cell types present: {present})"
        if other_volume:
            detail += (
                "; non-tetra volume cells "
                + ", ".join(sorted(other_volume))
                + " would have to be split into tetrahedra first"
            )
        raise MeshioFormatError("meshio mesh contains no tetra cells" + detail)

    if any(max(tet) >= len(vertices) for tet in tets):
        raise MeshioFormatError(
            "tetra block references a vertex outside the point array"
        )
    return TetrahedralMesh(
        vertices=vertices, tets=tets, units=units, coordinate_frame=coordinate_frame
    )



# ---------------------------------------------------------------------------
# Body-fitted volume mesh (tensorcad.volume_mesh) -> tetrahedra
# ---------------------------------------------------------------------------

def volume_mesh_to_tetrahedral(
    volume_mesh: Any,
    *,
    units: str = "m",
    coordinate_frame: str = "ship frame: x aft-to-fore, y to starboard, z up",
    tolerance: float | None = None,
) -> TetrahedralMesh:
    """Hand a body-fitted :class:`tensorcad.volume_mesh.VolumeMesh` to the volume formats.

    Why this exists
    ---------------
    :func:`tetrahedralise_surface` star-tessellates a *surface*, which is exact only for
    star-shaped bodies and is refused for the package's own hull generator (its rings
    collapse to zero breadth at the stem and post).  The documented remedy was "use
    :mod:`tensorcad.volume_mesh` and feed the tetrahedra instead" -- but nothing performed
    that hand-off, so the remedy was unbacked exactly like the earlier claim it replaced.
    This function is the missing link: a hull body grid (hexahedra, Kuhn-split into six
    tetrahedra each) becomes a :class:`TetrahedralMesh`, which SU2 and CGNS accept.

    What is checked, and not assumed
    --------------------------------
    * Every tetrahedron is classified on its *signed* volume.  A negative one is a real
      orientation defect (the solver would see an inverted cell), so it raises
      :class:`MeshioFormatError` naming the offending cell -- it is never reordered here.
    * Zero-volume cells are dropped and *counted* in
      :attr:`TetrahedralMesh.degenerate_cells`.  Structured hull grids always contain
      them (collapsed columns at the keel and the end stations); dropping is exact.
    * Conservation is asserted, not hoped for: the kept tetrahedra must reproduce
      ``volume_mesh.volume()`` to ``tolerance``, so a silent loss of measure is a raise.

    Parameters
    ----------
    volume_mesh :
        A :class:`tensorcad.volume_mesh.VolumeMesh` built with ``split_tets=True``.
    units, coordinate_frame :
        Carried onto the returned :class:`TetrahedralMesh`.
    tolerance :
        Absolute volume tolerance for the degenerate/conservation tests.  Defaults to
        ``1e-12 * volume()`` so it scales with the body.
    """
    from .volume_mesh import VolumeMesh

    if not isinstance(volume_mesh, VolumeMesh):
        raise TypeError(
            f"volume_mesh_to_tetrahedral expects a VolumeMesh, got {type(volume_mesh).__name__}; "
            "surface meshes belong to tetrahedralise_surface()"
        )
    if not volume_mesh.tets:
        raise MeshioFormatError(
            "volume mesh carries no tetrahedra (it was built with split_tets=False); "
            "rebuild with tensorcad.volume_mesh.build_volume_mesh(..., split_tets=True) "
            "before exporting to " + ", ".join(FORMATS_VOLUME)
        )

    vertices = tuple(tuple(float(c) for c in v) for v in volume_mesh.vertices)
    reference = volume_mesh.volume()
    tol = 1e-12 * abs(reference) if tolerance is None else float(tolerance)

    kept: list[Tet] = []
    degenerate = 0
    inverted: list[tuple[int, float]] = []
    for index, tet in enumerate(volume_mesh.tets):
        signed = tetrahedral_volume(vertices, tuple(int(i) for i in tet))
        if signed < -tol:
            inverted.append((index, signed))
        elif signed <= tol:
            degenerate += 1
        else:
            kept.append(tuple(int(i) for i in tet))

    if inverted:
        index, signed = inverted[0]
        raise MeshioFormatError(
            f"{len(inverted)} of {len(volume_mesh.tets)} tetrahedra are inverted "
            f"(cell {index} has signed volume {signed:.6e}); the source volume mesh is "
            "inconsistent and must be fixed at the builder, not reordered here"
        )

    kept_volume = sum(tetrahedral_volume(vertices, tet) for tet in kept)
    if abs(kept_volume - reference) > tol * max(1.0, len(kept)):
        raise MeshioFormatError(
            f"dropping {degenerate} zero-volume cells lost measure: kept tetrahedra sum to "
            f"{kept_volume:.16e} but the volume mesh measures {reference:.16e} "
            f"(difference {kept_volume - reference:.3e} exceeds {tol:.3e})"
        )

    return TetrahedralMesh(
        vertices=vertices,
        tets=tuple(kept),
        units=units,
        coordinate_frame=coordinate_frame,
        degenerate_cells=degenerate,
    )


# ---------------------------------------------------------------------------
# Boundary patches: the outer skin of the body grid, named for a solver
# ---------------------------------------------------------------------------
#
# A tetrahedral block on its own is not a usable solver mesh: the formats that carry volume
# cells expect the outer skin to be *named*, and SU2 in particular does nothing with a mesh
# whose ``NMARK`` is zero.  ``write_meshio(..., file_format="su2")`` writes exactly that --
# ``NMARK= 0`` -- because meshio's SU2 writer only emits markers for triangle blocks that
# carry integer cell data, and :class:`TetrahedralMesh` has no boundary at all.
#
# Naming the skin is therefore less a formatting chore than the place where the grid's two
# real defects have to be dealt with.
#
# **Coincident points.**  :mod:`tensorcad.volume_mesh` keeps a *distinct* node for every
# ``(station, level, breadth)`` index even where the breadth collapses, on purpose (see its
# docstring).  The cells of those columns have zero volume and
# :func:`volume_mesh_to_tetrahedral` drops them exactly, but dropping them exposes the
# *interior* faces of the surviving cells as if they were skin.  On the Wigley body used by
# ``examples/meshio_volume_bridge/run_example.py``, 1050 of the 2184 triangles that look like
# boundary are such zero-thickness slit faces -- each has an exact coincident twin with the
# opposite winding.  Handing those to a solver is not an approximation, it is a defect: a
# doubly covered, non-manifold wall.  :func:`volume_mesh_patches` welds coincident points, so
# a slit pair collapses onto one triangle and cancels, and what remains is a closed surface.
#
# **Anonymity.**  An unnamed skin cannot be referred to from a solver config.  Patches are
# assigned from the grid's own station and level planes, so a name means what the geometry
# means: ``hull`` is the body surface proper, ``fore``/``aft`` the two end stations,
# ``bottom``/``top`` the first and last level lines, ``symmetry`` the centreplane.
#
# The result is asserted, not assumed: the skin must be watertight (every edge used exactly
# twice, in opposite directions) and its divergence-theorem volume must reproduce the volume
# of the tetrahedra it was taken from.  A patch set that satisfied neither would still write
# a syntactically valid SU2 file -- which is exactly why the checks are here and not in a
# solver's error message.

#: Patch names, in the order they are tested.  The first plane a triangle lies on names it;
#: anything left over is ``hull``.  ``symmetry`` comes first because a triangle on the
#: centreplane is a symmetry boundary wherever else it also lies, and ``hull`` last because
#: the body surface is the general case, not a plane.
PATCH_NAMES: tuple[str, ...] = ("symmetry", "bottom", "top", "fore", "aft", "hull")

Tri = tuple[int, int, int]


def _oriented_tet_faces(tet: Tet) -> tuple[Tri, Tri, Tri, Tri]:
    """The four faces of a *positively oriented* tetrahedron, wound outward.

    For a tet whose signed volume is positive, reversing the cyclic order of three of the
    four index triples (and keeping the first) points every normal away from the fourth
    node; the derivation is the reason the ordering looks arbitrary and is not.
    """
    a, b, c, d = tet
    return ((a, c, b), (a, b, d), (b, c, d), (a, d, c))


def _sorted_face(face: Tri) -> Tri:
    """The triangle's three indices in ascending order, as a 3-tuple."""
    a, b, c = sorted(face)
    return (a, b, c)


def _winding(face: Tri) -> int:
    """``+1``/``-1`` for the two cyclic orientations of a triangle, independent of rotation.

    The face is compared with its own sorted order: rotating the triple to start at its
    smallest index and then asking whether the next two indices run upward gives the same
    answer for all three rotations of one orientation and the opposite for its reverse --
    which is exactly the property the slit cancellation relies on.
    """
    low = min(face)
    for start in range(3):
        if face[start] == low:
            following = (face[(start + 1) % 3], face[(start + 2) % 3])
            return 1 if following[0] < following[1] else -1
    raise AssertionError("unreachable: min() of a 3-tuple is one of its entries")


def _weld(vertices: Sequence[Vec3]) -> tuple[tuple[Vec3, ...], tuple[int, ...], int]:
    """Merge exactly coincident points, keeping first-use order.

    ``0.0`` and ``-0.0`` compare and hash equal, so the centreplane nodes the builder
    normalises to zero weld onto each other as well.
    """
    index: dict[Vec3, int] = {}
    welded: list[Vec3] = []
    remap: list[int] = []
    for vertex in vertices:
        key = (float(vertex[0]), float(vertex[1]), float(vertex[2]))
        found = index.get(key)
        if found is None:
            found = len(welded)
            index[key] = found
            welded.append(key)
        remap.append(found)
    return tuple(welded), tuple(remap), len(vertices) - len(welded)


def _skin_from_tets(
    tets: Sequence[Tet],
) -> tuple[tuple[Tri, ...], int]:
    """Outer skin of a positively oriented tetrahedral block, with the slit cancelled.

    A face that occurs once is skin.  A face that occurs twice must be the same triangle
    seen from both sides (the zero-thickness slit between two coincident columns); both
    copies are cancelled.  Anything else -- a face used by three cells, or twice with the
    *same* winding -- is a non-manifold or inverted tiling and raises.
    """
    entries: dict[Tri, list[tuple[Tri, int]]] = {}
    for tet in tets:
        for face in _oriented_tet_faces(tet):
            entries.setdefault(_sorted_face(face), []).append((face, _winding(face)))
    skin: list[Tri] = []
    cancelled = 0
    for rotated, seen in entries.items():
        if len(seen) == 1:
            skin.append(seen[0][0])
        elif len(seen) == 2:
            if seen[0][1] == seen[1][1]:
                raise MeshioFormatError(
                    f"triangle {rotated} is used twice with the same winding: two "
                    "coincident cells face the same way, so the block is folded onto "
                    "itself (a zero-thickness sheet), not a tiling of a body"
                )
            cancelled += 1
        else:
            raise MeshioFormatError(
                f"triangle {rotated} is shared by {len(seen)} tetrahedra; a surface that "
                "is not edge-to-edge is not a skin"
            )
    return tuple(sorted(skin)), cancelled


def _on_plane(values: Sequence[float], plane: float, tolerance: float) -> bool:
    return all(abs(value - plane) <= tolerance for value in values)


@dataclass(frozen=True)
class VolumePatches:
    """A welded tetrahedral body mesh together with its named outer skin.

    ``vertices`` holds only the points the kept tetrahedra actually use (coincident points
    merged, unreferenced points dropped), ``tets`` index into it, and ``patches`` is a tuple
    of ``(name, faces)`` in :data:`PATCH_NAMES` order with outward-wound triangles indexing
    the same ``vertices``.

    ``merged_points`` counts the coincident points the weld removed and ``cancelled_faces``
    the slit triangles that cancelled -- both are properties of the source grid that a
    reader of the exported file can no longer see, so they are carried rather than hidden.
    """

    vertices: tuple[Vec3, ...]
    tets: tuple[Tet, ...]
    patches: tuple[tuple[str, tuple[Tri, ...]], ...]
    merged_points: int = 0
    cancelled_faces: int = 0
    degenerate_cells: int = 0
    units: str = "m"
    coordinate_frame: str = "unspecified right-handed frame"

    def names(self) -> tuple[str, ...]:
        """Non-empty patch names, in priority order."""
        return tuple(name for name, _ in self.patches)

    def faces(self, name: str) -> tuple[Tri, ...]:
        """Outward-wound triangles of one patch."""
        for patch_name, faces in self.patches:
            if patch_name == name:
                return faces
        raise KeyError(f"no patch named {name!r}; this mesh has {self.names()}")

    def counts(self) -> dict[str, int]:
        """Triangle count per patch."""
        return {name: len(faces) for name, faces in self.patches}

    def area(self, name: str) -> float:
        """Total area of one patch, in square metres."""
        return sum(_triangle_area(self.vertices, face) for face in self.faces(name))

    def areas(self) -> dict[str, float]:
        """Total area per patch, in square metres."""
        return {name: self.area(name) for name, _ in self.patches}

    def skin(self) -> tuple[Tri, ...]:
        """Every skin triangle, all patches together."""
        return tuple(face for _, faces in self.patches for face in faces)

    def skin_volume(self) -> float:
        """Volume enclosed by the oriented skin (divergence theorem)."""
        return sum(_triangle_solid_angle_volume(self.vertices, face) for face in self.skin())

    def tet_volume(self) -> float:
        """Total volume of the tetrahedra, by their signed volumes."""
        return sum(tetrahedral_volume(self.vertices, tet) for tet in self.tets)

    def edge_counts(self) -> dict[tuple[int, int], int]:
        """How often each undirected *directed-consistent* edge is used by the skin.

        A closed, consistently oriented surface uses every undirected edge exactly twice,
        once in each direction; :meth:`is_watertight` is the assertion of that.
        """
        counts: dict[tuple[int, int], int] = {}
        for a, b, c in self.skin():
            for tail, head in ((a, b), (b, c), (c, a)):
                counts[(tail, head)] = counts.get((tail, head), 0) + 1
        return counts

    def is_watertight(self) -> bool:
        """True when every edge is used once in each direction on the skin."""
        counts = self.edge_counts()
        return all(
            counts.get(edge, 0) == 1 and counts.get((edge[1], edge[0]), 0) == 1
            for edge in counts
        )


def _triangle_area(vertices: Sequence[Vec3], face: Tri) -> float:
    p0, p1, p2 = (vertices[index] for index in face)
    ux, uy, uz = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
    vx, vy, vz = p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]
    cx = uy * vz - uz * vy
    cy = uz * vx - ux * vz
    cz = ux * vy - uy * vx
    return 0.5 * math.sqrt(cx * cx + cy * cy + cz * cz)


def _triangle_solid_angle_volume(vertices: Sequence[Vec3], face: Tri) -> float:
    """One term of the divergence-theorem volume of a closed oriented surface."""
    ax, ay, az = vertices[face[0]]
    bx, by, bz = vertices[face[1]]
    cx, cy, cz = vertices[face[2]]
    return (
        ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)
    ) / 6.0


def volume_mesh_patches(
    volume_mesh: Any,
    *,
    tetrahedral: TetrahedralMesh | None = None,
    tolerance: float | None = None,
) -> VolumePatches:
    """Weld the body grid's tetrahedra and name their outer skin.

    Parameters
    ----------
    volume_mesh :
        A :class:`tensorcad.volume_mesh.VolumeMesh` with ``split_tets=True`` (the default).
    tetrahedral :
        The already-converted tetrahedra, to avoid converting twice.  ``None`` calls
        :func:`volume_mesh_to_tetrahedral`.
    tolerance :
        Absolute distance under which a point counts as lying on a grid plane.  ``None``
        uses ``1e-12`` of the mesh extent, which is far above the round-off of the builder's
        own coordinates (they *are* ``grid.xs``/``grid.zs``) and far below any real geometry.

    Raises
    ------
    MeshioFormatError
        If the skin is not watertight, if a patch is empty, or if the enclosed volume does
        not reproduce the tetrahedra's.
    """
    tet_mesh = volume_mesh_to_tetrahedral(volume_mesh) if tetrahedral is None else tetrahedral
    welded, remap, merged = _weld(tet_mesh.vertices)

    kept: list[Tet] = []
    for index, tet in enumerate(tet_mesh.tets):
        mapped = (remap[tet[0]], remap[tet[1]], remap[tet[2]], remap[tet[3]])
        if len(set(mapped)) != 4:
            raise MeshioFormatError(
                f"tetrahedron {index} has two nodes that coincide after welding "
                f"({tet} -> {mapped}); a positive-volume cell cannot, so the source "
                "topology is inconsistent"
            )
        kept.append(mapped)

    skin, cancelled = _skin_from_tets(kept)
    if not skin:
        raise MeshioFormatError("the tetrahedral block has no skin: it is empty")

    used = sorted({vertex for tet in kept for vertex in tet})
    compact = {old: new for new, old in enumerate(used)}
    vertices = tuple(welded[old] for old in used)
    tets = tuple(
        (compact[a], compact[b], compact[c], compact[d]) for a, b, c, d in kept
    )

    xs = tuple(float(x) for x in getattr(volume_mesh, "station_lines", ()))
    zs = tuple(float(z) for z in getattr(volume_mesh, "level_lines", ()))
    if not xs or not zs:
        raise MeshioFormatError(
            "the volume mesh carries no station/level lines, so its skin cannot be named"
        )
    span = max(
        max(coordinate[axis] for coordinate in vertices) - min(coordinate[axis] for coordinate in vertices)
        for axis in range(3)
    )
    tol = (1e-12 * max(1.0, span)) if tolerance is None else float(tolerance)

    bucketed: dict[str, list[Tri]] = {name: [] for name in PATCH_NAMES}
    for face in skin:
        points = [vertices[compact[index]] for index in face]
        if _on_plane([point[1] for point in points], 0.0, tol):
            name = "symmetry"
        elif _on_plane([point[2] for point in points], zs[0], tol):
            name = "bottom"
        elif _on_plane([point[2] for point in points], zs[-1], tol):
            name = "top"
        elif _on_plane([point[0] for point in points], xs[0], tol):
            name = "fore"
        elif _on_plane([point[0] for point in points], xs[-1], tol):
            name = "aft"
        else:
            name = "hull"
        bucketed[name].append((compact[face[0]], compact[face[1]], compact[face[2]]))

    patches = tuple((name, tuple(faces)) for name, faces in bucketed.items() if faces)
    result = VolumePatches(
        vertices=vertices,
        tets=tets,
        patches=patches,
        merged_points=merged,
        cancelled_faces=cancelled,
        degenerate_cells=tet_mesh.degenerate_cells,
        units=tet_mesh.units,
        coordinate_frame=tet_mesh.coordinate_frame,
    )

    if not result.is_watertight():
        raise MeshioFormatError(
            "the named skin is not watertight: some edge is not used exactly once in each "
            "direction, so the marker set would leave a hole in the body"
        )
    reference = tet_mesh.volume()
    if abs(result.tet_volume() - reference) > max(abs(reference), 1.0) * 1e-12:
        raise MeshioFormatError(
            f"welding changed the tetrahedra's volume: {result.tet_volume():.16e} versus "
            f"{reference:.16e}"
        )
    if abs(result.skin_volume() - reference) > max(abs(reference), 1.0) * 1e-12:
        raise MeshioFormatError(
            f"the skin encloses {result.skin_volume():.16e} but the tetrahedra measure "
            f"{reference:.16e}: the patch set does not bound the body"
        )
    return result


def write_su2_volume_mesh(
    path: str | Path,
    volume_mesh: Any,
    *,
    patches: VolumePatches | None = None,
) -> Path:
    """Write a body-fitted volume grid as an SU2 ASCII mesh *with named markers*.

    The file is written here rather than through meshio because meshio's SU2 writer only
    accepts integer marker tags: it replaces a string tag with its own running index, and a
    solver config that says ``MARKER_HEATFLUX= (hull, 0.0)`` cannot then find its wall.
    Writing string tags is four lines of ASCII, so the format is emitted directly and read
    back by :func:`read_su2_markers` (and, independently, by meshio's reader, which accepts
    the file -- see ``tests/test_volume_mesh_su2_markers.py``).

    Returns
    -------
    Path
        The file written.  Call :func:`volume_mesh_patches` (or pass ``patches=`` to obtain
        the same object twice-free) for the counts and areas that went into it.
    """
    skin = volume_mesh_patches(volume_mesh) if patches is None else patches
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)

    lines = ["NDIME= 3", f"NPOIN= {len(skin.vertices)}"]
    lines += [f"{x:.16e} {y:.16e} {z:.16e}" for x, y, z in skin.vertices]
    lines.append(f"NELEM= {len(skin.tets)}")
    lines += [f"10 {a} {b} {c} {d}" for a, b, c, d in skin.tets]
    lines.append(f"NMARK= {len(skin.patches)}")
    for name, faces in skin.patches:
        lines.append(f"MARKER_TAG= {name}")
        lines.append(f"MARKER_ELEMS= {len(faces)}")
        lines += [f"5 {a} {b} {c}" for a, b, c in faces]
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def read_su2_markers(path: str | Path) -> dict[str, tuple[Tri, ...]]:
    """Read the marker sections of an SU2 ASCII mesh, keeping the tags as names.

    Only what the writer above emits is parsed (plus the coordinate and element blocks, which
    are skipped): this is the reader that closes the loop on a file whose tags meshio would
    flatten to indices.  Boundary element types other than triangles are accepted and
    discarded, so a hand-edited file still reads.
    """
    text = Path(path).read_text(encoding="utf-8").splitlines()
    markers: dict[str, tuple[Tri, ...]] = {}
    tag: str | None = None
    index = 0
    while index < len(text):
        line = text[index].strip()
        index += 1
        if not line or line.startswith("%"):
            continue
        if line.startswith(("NPOIN=", "NELEM=")):
            index += int(line.split("=", 1)[1].split()[0])
        elif line.startswith("MARKER_TAG="):
            tag = line.split("=", 1)[1].strip()
        elif line.startswith("MARKER_ELEMS="):
            count = int(line.split("=", 1)[1])
            if tag is None:
                raise MeshioFormatError("MARKER_ELEMS before any MARKER_TAG")
            faces: list[Tri] = []
            for _ in range(count):
                entry = text[index].split()
                index += 1
                if int(entry[0]) == 5:
                    faces.append((int(entry[1]), int(entry[2]), int(entry[3])))
            if tag in markers:
                raise MeshioFormatError(f"marker {tag!r} appears twice")
            markers[tag] = tuple(faces)
            tag = None
    return markers


#: gmsh mesh format written by :func:`write_gmsh_volume_mesh`.  Version 2.2 is the last
#: ASCII version that carries per-element physical tags in one flat ``$Elements`` block, and
#: it is what ``gmsh`` itself, ``dolfin-convert``, ``gmshio`` and TensorCAD's own FEM handoff
#: (``tensorcad.solver_interface.read_gmsh_volume``) all read.
GMSH_VOLUME_VERSION = "2.2"

#: Physical name given to the tetrahedral volume group.
GMSH_SOLID_GROUP = "solid"


def write_gmsh_volume_mesh(
    path: str | Path,
    volume_mesh: Any,
    *,
    patches: VolumePatches | None = None,
    solid_name: str = GMSH_SOLID_GROUP,
) -> Path:
    """Write a body-fitted volume grid as gmsh 2.2 ASCII *with named physical groups*.

    This is the gmsh counterpart of :func:`write_su2_volume_mesh`: the same named skin, written
    into the format the open-source solver stack actually consumes.  It is a separate entry
    point rather than a new member of :data:`FORMATS_ALL` on purpose -- ``FORMATS_SURFACE``
    classifies ``gmsh`` as a *surface* format (that is what the 18-combination matrix of
    ``tests/test_meshio_format_matrix.py`` measures, through meshio), and reclassifying it
    would silently change the meaning of that matrix.  A volume grid leaves through *this*
    function, and comes back through :func:`read_gmsh_volume_groups`.

    What lands in the file:

    * one ``$PhysicalNames`` entry per patch -- ``hull``, ``symmetry``, ``bottom``, ``top``,
      ``fore``, ``aft`` as :func:`volume_mesh_patches` names them, dimension 2 -- plus the
      tetrahedral volume group ``solid`` at dimension 3 (tag 1, patches from tag 2 on);
    * the welded points, ``repr``-formatted, so a read returns the *same doubles* the grid
      builder produced (a ``%.6f`` would move a node by up to a micrometre and turn a
      zero-error control into a 1e-6 error);
    * every skin triangle (gmsh element type 2) carrying its physical tag, then every
      tetrahedron (type 4), both with the two mandatory tags (physical, elementary) set to the
      physical tag -- gmsh's reader requires the elementary entity to exist, and reusing the
      physical tag keeps the file readable by tools that only understand one of the two.

    Parameters
    ----------
    path :
        Output ``.msh`` path.
    volume_mesh :
        A :class:`tensorcad.volume_mesh.VolumeMesh` built with ``split_tets=True``.
    patches :
        The already-computed :func:`volume_mesh_patches` result, to avoid deriving it twice.
    solid_name :
        Physical name of the volume group.  Must differ from every patch name.

    Returns
    -------
    Path
        The file written.
    """
    skin = volume_mesh_patches(volume_mesh) if patches is None else patches
    if solid_name in skin.names():
        raise MeshioFormatError(
            f"solid name {solid_name!r} collides with a patch name; the physical groups "
            f"would be indistinguishable (patches: {skin.names()})"
        )
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)

    # (name, dimension, physical tag); the volume group is tag 1, patches follow.
    groups = [(solid_name, 3, 1)] + [
        (name, 2, index + 2) for index, (name, _) in enumerate(skin.patches)
    ]
    lines = [
        "$MeshFormat",
        f"{GMSH_VOLUME_VERSION} 0 8",
        "$EndMeshFormat",
        "$Comments",
        "generated by TensorCAD tensorcad.meshio_exchange.write_gmsh_volume_mesh",
        f"units {skin.units}",
        f"coordinate_frame {skin.coordinate_frame}",
        "$EndComments",
        "$PhysicalNames",
        str(len(groups)),
    ]
    lines += [f'{dimension} {tag} "{name}"' for name, dimension, tag in groups]
    lines.append("$EndPhysicalNames")
    lines.append("$Nodes")
    lines.append(str(len(skin.vertices)))
    lines += [
        f"{index} {x!r} {y!r} {z!r}"
        for index, (x, y, z) in enumerate(skin.vertices, start=1)
    ]
    lines.append("$EndNodes")

    elements: list[str] = []
    for name, _dimension, tag in groups[1:]:
        for a, b, c in skin.faces(name):
            elements.append(f"{len(elements) + 1} 2 2 {tag} {tag} {a + 1} {b + 1} {c + 1}")
    for a, b, c, d in skin.tets:
        elements.append(f"{len(elements) + 1} 4 2 1 1 {a + 1} {b + 1} {c + 1} {d + 1}")
    lines.append("$Elements")
    lines.append(str(len(elements)))
    lines += elements
    lines.append("$EndElements")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def _gmsh_block(text: str, name: str, path: str | Path) -> list[str]:
    """Lines of one ``$name`` block of a gmsh ASCII file, stripped, blanks removed."""
    marker, end = f"${name}", f"$End{name}"
    if marker not in text:
        raise MeshioFormatError(f"{marker} is missing from {path}")
    body = text.split(marker, 1)[1].split(end, 1)[0]
    return [line.strip() for line in body.splitlines() if line.strip()]


def read_gmsh_volume_groups(path: str | Path) -> VolumePatches:
    """Read a gmsh 2.2 ASCII volume mesh, keeping the physical groups as patch names.

    An independent parser, not a call into meshio: the point of writing the format by hand is
    that reading it back with the same library would only prove self-consistency.  Two
    readers that agree (this one and meshio's, see
    ``tests/test_gmsh_volume_patches.py``) plus the volume an independent writer measures is
    the evidence the file is what it claims to be.

    Returns
    -------
    VolumePatches
        ``vertices``/``tets``/``patches`` exactly as the file states them, with ``units`` and
        ``coordinate_frame`` taken from the ``$Comments`` block when present.  ``merged_points``
        and ``cancelled_faces`` are **zero**: they describe the hex-grid conversion that
        produced the file, which a reader cannot reconstruct, and reporting a number there
        would be inventing one.

    Raises
    ------
    MeshioFormatError
        On a version other than 2.2, on a file with no tetrahedra (that is a surface mesh --
        use :func:`read_meshio`), or on a boundary element type this reader cannot express as
        a triangle.
    """
    raw = Path(path).read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise MeshioFormatError(
            f"{path} is not gmsh ASCII text (byte {exc.start} is not UTF-8): a binary gmsh "
            "file cannot carry the per-element tags this reader groups by name; write it "
            "with binary=False, or with write_gmsh_volume_mesh"
        ) from exc

    header = _gmsh_block(text, "MeshFormat", path)
    if not header:
        raise MeshioFormatError(f"empty $MeshFormat block in {path}")
    version = header[0].split()[0]
    if not version.startswith("2."):
        raise MeshioFormatError(
            f"{path} is gmsh {version}, but only 2.x ASCII carries per-element physical tags "
            "in one $Elements block; re-export with gmsh -format msh2"
        )

    units, coordinate_frame = "m", "unspecified right-handed frame"
    if "$Comments" in text:
        for line in _gmsh_block(text, "Comments", path):
            key, _, value = line.partition(" ")
            if key == "units" and value:
                units = value
            elif key == "coordinate_frame" and value:
                coordinate_frame = value

    names: dict[int, tuple[str, int]] = {}
    if "$PhysicalNames" in text:
        entries = _gmsh_block(text, "PhysicalNames", path)
        count = int(entries[0])
        for line in entries[1 : count + 1]:
            parts = line.replace('"', " ").split()
            if len(parts) < 3:
                raise MeshioFormatError(f"malformed $PhysicalNames entry {line!r} in {path}")
            names[int(parts[1])] = (parts[2], int(parts[0]))

    nodes = _gmsh_block(text, "Nodes", path)
    count = int(nodes[0])
    if len(nodes) - 1 < count:
        raise MeshioFormatError(f"{path} declares {count} nodes but carries {len(nodes) - 1}")
    vertices: list[Vec3] = []
    for line in nodes[1 : count + 1]:
        parts = line.split()
        vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))

    elements = _gmsh_block(text, "Elements", path)
    nelem = int(elements[0])
    if len(elements) - 1 < nelem:
        raise MeshioFormatError(f"{path} declares {nelem} elements but carries {len(elements) - 1}")
    tets: list[Tet] = []
    bucketed: dict[str, list[Tri]] = {}
    for line in elements[1 : nelem + 1]:
        parts = line.split()
        kind, tags = int(parts[1]), int(parts[2])
        if tags < 1:
            raise MeshioFormatError(
                f"element {parts[0]} in {path} carries no physical tag, so it belongs to no "
                "patch; a solver file needs every boundary element marked"
            )
        physical = int(parts[3])
        connectivity = [int(value) - 1 for value in parts[3 + tags :]]
        name, _dimension = names.get(physical, (f"group_{physical}", 0))
        if kind == 4:
            if len(connectivity) != 4:
                raise MeshioFormatError(
                    f"tetrahedron {parts[0]} in {path} has {len(connectivity)} nodes"
                )
            tets.append((connectivity[0], connectivity[1], connectivity[2], connectivity[3]))
        elif kind == 2:
            if len(connectivity) != 3:
                raise MeshioFormatError(
                    f"triangle {parts[0]} in {path} has {len(connectivity)} nodes"
                )
            bucketed.setdefault(name, []).append(
                (connectivity[0], connectivity[1], connectivity[2])
            )
        else:
            raise MeshioFormatError(
                f"element type {kind} (element {parts[0]}) in {path} cannot be expressed as a "
                "tetrahedron or a boundary triangle; only types 2 and 4 are read here"
            )

    if not tets:
        raise MeshioFormatError(
            f"{path} carries no tetrahedra: it is a surface mesh. Read it with "
            "read_meshio(file_format='gmsh') instead"
        )
    patches = tuple((name, tuple(faces)) for name, faces in bucketed.items())
    if not patches:
        raise MeshioFormatError(
            f"{path} carries no marked boundary triangles; a solver cannot bind a wall to it"
        )
    return VolumePatches(
        vertices=tuple(vertices),
        tets=tuple(tets),
        patches=patches,
        units=units,
        coordinate_frame=coordinate_frame,
    )


def _check_topology_supported(mesh: Mesh | TetrahedralMesh, file_format: str) -> None:
    """Raise :class:`MeshioFormatError` if the format cannot carry this topology."""
    if file_format not in _EXT:
        raise MeshioFormatError(
            f"unknown format '{file_format}'; expected one of {FORMATS_ALL}"
        )
    if file_format in FORMATS_VOLUME and not isinstance(mesh, TetrahedralMesh):
        raise MeshioFormatError(
            f"format '{file_format}' only supports volume cells (tetrahedra); "
            f"a surface triangle mesh cannot be written to it. "
            f"Call tetrahedralise_surface(mesh) first, or use one of {FORMATS_SURFACE}"
        )
    if file_format in FORMATS_SURFACE and isinstance(mesh, TetrahedralMesh):
        raise MeshioFormatError(
            f"format '{file_format}' is a surface format and only supports triangle "
            f"cells; use one of {FORMATS_VOLUME} or export the boundary surface instead"
        )


def write_meshio(path: str | Path, mesh: Mesh | TetrahedralMesh, *, file_format: str) -> Path:
    """Write a TensorCAD mesh to a file using meshio.

    Parameters
    ----------
    path :
        Output file path.
    mesh :
        Surface triangle mesh (STL/OBJ/VTK/Gmsh) or tetrahedral mesh
        (SU2/CGNS).
    file_format :
        One of ``stl``, ``obj``, ``vtk``, ``gmsh``, ``su2``, ``cgns``.

    Raises
    ------
    MeshioFormatError
        If the format cannot carry the mesh's topology.
    """
    _check_topology_supported(mesh, file_format)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    if file_format == "cgns":
        # meshio's CGNS driver writes no SIDS metadata at all (no
        # CGNSLibraryVersion, no ElementType, no node labels), so use the
        # conforming exporter instead; it keeps the node names meshio reads.
        from .mesh_cgns import write_cgns

        return write_cgns(target, mesh)
    import meshio

    mio = mesh_to_meshio(mesh)
    meshio.write(str(target), mio, file_format=file_format)
    return target


def read_meshio(
    path: str | Path,
    *,
    file_format: str,
    units: str = "m",
    coordinate_frame: str = "unspecified right-handed frame",
) -> Mesh | TetrahedralMesh:
    """Read a meshio file into a TensorCAD mesh.

    Volume formats (SU2/CGNS) come back as :class:`TetrahedralMesh`; surface
    formats as :class:`tensorcad.geometry.Mesh`.
    """
    import meshio

    if file_format not in _EXT:
        raise MeshioFormatError(
            f"unknown format '{file_format}'; expected one of {FORMATS_ALL}"
        )
    if file_format == "cgns":
        # Prefer the spec reader: it understands both the SIDS layout and
        # meshio's, and meshio's reader refuses anything but tetrahedra.
        from .mesh_cgns import read_cgns

        try:
            return read_cgns(path, units=units, coordinate_frame=coordinate_frame)
        except Exception:  # pragma: no cover - fall back to meshio's reader
            pass
    mio = meshio.read(str(path), file_format=file_format)
    if file_format in FORMATS_VOLUME:
        return meshio_to_tetrahedral(mio, units=units, coordinate_frame=coordinate_frame)
    return meshio_to_mesh(mio, units=units, coordinate_frame=coordinate_frame)


def round_trip_meshio(
    mesh: Mesh | TetrahedralMesh,
    *,
    file_format: str,
    tmpdir: str | Path | None = None,
) -> Mesh | TetrahedralMesh:
    """Write then read back a mesh via meshio, returning the round-tripped mesh.

    Returns the same kind of object as the input (``Mesh`` for surface formats,
    ``TetrahedralMesh`` for SU2/CGNS).

    Parameters
    ----------
    mesh :
        Original mesh (surface triangles or tetrahedra).
    file_format :
        Target format (stl/obj/vtk/gmsh/su2/cgns).
    tmpdir :
        Directory for temporary files.  If None, uses the system temp dir.
    """
    import meshio

    _check_topology_supported(mesh, file_format)
    ext = _EXT[file_format]
    with tempfile.TemporaryDirectory(dir=tmpdir) as td:
        path = Path(td) / f"roundtrip{ext}"
        write_meshio(path, mesh, file_format=file_format)
        back = meshio.read(str(path), file_format=file_format)
        if file_format in FORMATS_VOLUME:
            return meshio_to_tetrahedral(
                back, units=mesh.units, coordinate_frame=mesh.coordinate_frame
            )
        return meshio_to_mesh(back, units=mesh.units, coordinate_frame=mesh.coordinate_frame)
