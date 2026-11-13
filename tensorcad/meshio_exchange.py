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

(*) SU2/CGNS surface triangles remain genuinely unsupported -- the formats only
carry volume cells -- but :func:`round_trip_meshio` now raises a
``MeshioFormatError`` telling the caller to call
:func:`tetrahedralise_surface` instead of silently having nothing to offer.

A caveat that is stated rather than hidden: the CGNS path goes through meshio's
CGNS writer, which stores a *single* uniform element block (tetra only) and does
not emit the full CGNS/HDF5 metadata (CGNSLibraryVersion, ElementType, ...).
Round-trip fidelity through meshio is REAL; conformance as a CGNS file for a
third-party CGNS reader is **WITHHELD** -- no ``cgnslib`` is available on this
machine to check it independently.
"""
from __future__ import annotations

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
    "MeshioFormatError",
    "TetrahedralMesh",
    "TetrahedralisationError",
    "mesh_to_meshio",
    "meshio_to_mesh",
    "meshio_to_tetrahedral",
    "read_meshio",
    "round_trip_meshio",
    "surface_signed_volume",
    "tetrahedral_volume",
    "tetrahedralise_surface",
    "write_meshio",
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

    ``tets`` holds 0-based indices into ``vertices``.  All tetrahedra are
    expected to be positively oriented (positive signed volume); the star
    tessellation produced by :func:`tetrahedralise_surface` guarantees it.
    """

    vertices: tuple[Vec3, ...]
    tets: tuple[Tet, ...]
    units: str = "m"
    coordinate_frame: str = "unspecified right-handed frame"

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
    bounding-box centre, and +/- 5 % nudges along each axis); the candidate with
    the largest minimum tetrahedron volume wins.  Degenerate (coplanar) cells are
    rejected, never emitted: a solver receiving a zero-volume tetrahedron fails
    much later and much less clearly than one that never got the mesh.

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
    # Volume below this is numerically coplanar: reject the apex instead of
    # handing a zero-volume cell to a solver.
    tol = 1e-12 * (scale**3 if scale > 0.0 else 1.0)

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

    def min_tet_volume(apex: Vec3) -> float:
        return float(tet_volumes(apex).min())

    if interior_point is not None:
        apex = tuple(float(c) for c in interior_point)
        volumes = tet_volumes(apex)
        worst = float(volumes.min())
        if worst <= tol:
            bad = [int(i) for i in np.flatnonzero(volumes <= tol)]
            raise TetrahedralisationError(
                f"{len(bad)} of {len(faces)} tetrahedra are non-positive or coplanar "
                f"(worst volume {worst:.3e} <= {tol:.3e}, first face {bad[0]}); the body "
                f"is not star-shaped with respect to {apex!r}. Pass a different "
                f"interior_point= or use a volume mesher."
            )
    else:
        best: tuple[float, Vec3] | None = None
        for candidate in _interior_candidates(mesh):
            worst = min_tet_volume(candidate)
            if best is None or worst > best[0]:
                best = (worst, candidate)
        assert best is not None
        worst, apex = best
        if worst <= tol:
            raise TetrahedralisationError(
                f"no interior point sees the whole surface: best candidate {apex!r} still "
                f"gives a worst-case tetrahedron volume {worst:.3e} (needs > {tol:.3e}). "
                f"The body is not star-shaped; use a volume mesher instead."
            )

    apex_index = len(verts)
    verts.append(apex)
    tets: Tet = tuple((apex_index, a, b, c) for a, b, c in faces)

    return TetrahedralMesh(
        vertices=tuple(verts),
        tets=tets,
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
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


def _find_block(mio: Any, cell_type: str) -> np.ndarray | None:
    for block in _cell_blocks(mio.cells):
        if block.type == cell_type:
            return np.asarray(block.data)
    return None


def meshio_to_mesh(
    mio: Any,
    *,
    units: str = "m",
    coordinate_frame: str = "unspecified right-handed frame",
) -> Mesh:
    """Convert a ``meshio.Mesh`` back to a :class:`tensorcad.geometry.Mesh`.

    Only the first triangle cell block is used; other cell types are ignored.
    """
    triangle_data = _find_block(mio, "triangle")
    if triangle_data is None:
        raise MeshioFormatError("meshio mesh contains no triangle cells")

    vertices = tuple(tuple(float(c) for c in row) for row in np.asarray(mio.points))
    faces = tuple(tuple(int(idx) for idx in row) for row in triangle_data)
    return Mesh(vertices=vertices, faces=faces, units=units, coordinate_frame=coordinate_frame)


def meshio_to_tetrahedral(
    mio: Any,
    *,
    units: str = "m",
    coordinate_frame: str = "unspecified right-handed frame",
) -> TetrahedralMesh:
    """Convert a ``meshio.Mesh`` carrying a tetra cell block to a TetrahedralMesh."""
    tetra_data = _find_block(mio, "tetra")
    if tetra_data is None:
        found = ", ".join(sorted({b.type for b in _cell_blocks(mio.cells)})) or "none"
        raise MeshioFormatError(
            f"meshio mesh contains no tetra cells (cell types present: {found})"
        )

    vertices = tuple(tuple(float(c) for c in row) for row in np.asarray(mio.points))
    tets = tuple(tuple(int(idx) for idx in row) for row in tetra_data)
    return TetrahedralMesh(
        vertices=vertices, tets=tets, units=units, coordinate_frame=coordinate_frame
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
    import meshio

    _check_topology_supported(mesh, file_format)
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
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
        io = mesh_to_meshio(mesh)
        meshio.write(str(path), io, file_format=file_format)
        back = meshio.read(str(path), file_format=file_format)
        if file_format in FORMATS_VOLUME:
            return meshio_to_tetrahedral(
                back, units=mesh.units, coordinate_frame=mesh.coordinate_frame
            )
        return meshio_to_mesh(back, units=mesh.units, coordinate_frame=mesh.coordinate_frame)
