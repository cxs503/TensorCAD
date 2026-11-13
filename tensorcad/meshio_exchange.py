"""meshio-based format exchange for TensorCAD meshes.

Provides round-trip conversion between :class:`tensorcad.geometry.Mesh` and
meshio for the six formats required by the format-matrix audit:
STL, OBJ, VTK, SU2, CGNS, Gmsh.

Format capability matrix (surface triangle meshes)
===================================================
======  ===================  =============================================
Format  surface triangles    notes
======  ===================  =============================================
STL     REAL                 reorders vertices/faces; volume preserved
OBJ     REAL                 exact round-trip (points + topology)
VTK     REAL                 exact round-trip (points + topology)
SU2     NOT_SUPPORTED        only volume cells (tetra/hex/wedge/pyramid)
CGNS    NOT_SUPPORTED        only tetrahedral volume cells
Gmsh    REAL                 exact round-trip (points + topology)
======  ===================  =============================================

SU2 and CGNS are documented as NOT_SUPPORTED for *surface* meshes because
their meshio writers only emit volume cell types.  A caller that needs to
export to these formats must first tetrahedralise the surface (e.g. via
:mod:`tensorcad.volume_mesh` or an external mesher).
"""
from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Any

import numpy as np

from .geometry import Mesh

__all__ = [
    "FORMATS_SURFACE",
    "FORMATS_VOLUME",
    "MeshioFormatError",
    "mesh_to_meshio",
    "meshio_to_mesh",
    "write_meshio",
    "read_meshio",
    "round_trip_meshio",
]

#: Formats that support surface triangle meshes.
FORMATS_SURFACE: tuple[str, ...] = ("stl", "obj", "vtk", "gmsh")

#: Formats that only support volume cells (require tetrahedralisation first).
FORMATS_VOLUME: tuple[str, ...] = ("su2", "cgns")

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


def mesh_to_meshio(mesh: Mesh) -> Any:
    """Convert a :class:`tensorcad.geometry.Mesh` to a ``meshio.Mesh``.

    Returns a meshio Mesh object with points (N×3 float64) and a single
    triangle cell block.
    """
    import meshio

    points = np.array(mesh.vertices, dtype=np.float64)
    cells = np.array(mesh.faces, dtype=np.int64)
    return meshio.Mesh(points, [("triangle", cells)])


def meshio_to_mesh(mio: Any, *, units: str = "m", coordinate_frame: str = "unspecified right-handed frame") -> Mesh:
    """Convert a ``meshio.Mesh`` back to a :class:`tensorcad.geometry.Mesh`.

    Only the first triangle cell block is used; other cell types are ignored.
    """
    # Find triangle cells
    triangle_data = None
    for block in mio.cells:
        if block.type == "triangle":
            triangle_data = block.data
            break
    if triangle_data is None:
        raise MeshioFormatError("meshio mesh contains no triangle cells")

    vertices = tuple(tuple(float(c) for c in row) for row in mio.points)
    faces = tuple(tuple(int(idx) for idx in row) for row in triangle_data)
    return Mesh(vertices=vertices, faces=faces, units=units, coordinate_frame=coordinate_frame)


def write_meshio(path: str | Path, mesh: Mesh, *, file_format: str) -> Path:
    """Write a TensorCAD mesh to a file using meshio.

    Parameters
    ----------
    path :
        Output file path.
    mesh :
        The surface mesh to write.
    file_format :
        One of ``stl``, ``obj``, ``vtk``, ``gmsh``, ``su2``, ``cgns``.

    Raises
    ------
    MeshioFormatError
        If the format does not support surface triangle meshes (SU2, CGNS).
    """
    import meshio

    if file_format in FORMATS_VOLUME:
        raise MeshioFormatError(
            f"format '{file_format}' only supports volume cells; "
            f"tetrahedralise the surface first or use one of {FORMATS_SURFACE}"
        )

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    mio = mesh_to_meshio(mesh)
    meshio.write(str(target), mio, file_format=file_format)
    return target


def read_meshio(path: str | Path, *, file_format: str, units: str = "m", coordinate_frame: str = "unspecified right-handed frame") -> Mesh:
    """Read a meshio file into a TensorCAD Mesh.

    Parameters
    ----------
    path :
        Input file path.
    file_format :
        Format hint (passed to meshio.read).
    units :
        Units to attach to the resulting Mesh.
    coordinate_frame :
        Coordinate frame description to attach to the resulting Mesh.
    """
    import meshio

    mio = meshio.read(str(path), file_format=file_format)
    return meshio_to_mesh(mio, units=units, coordinate_frame=coordinate_frame)


def round_trip_meshio(mesh: Mesh, *, file_format: str, tmpdir: str | Path | None = None) -> Mesh:
    """Write then read back a mesh via meshio, returning the round-tripped mesh.

    Parameters
    ----------
    mesh :
        Original surface mesh.
    file_format :
        Target format (stl/obj/vtk/gmsh).
    tmpdir :
        Directory for temporary files.  If None, uses the system temp dir.

    Returns
    -------
    Mesh
        The mesh after write→read round trip.
    """
    import meshio

    if file_format in FORMATS_VOLUME:
        raise MeshioFormatError(
            f"format '{file_format}' only supports volume cells; "
            f"tetrahedralise the surface first or use one of {FORMATS_SURFACE}"
        )

    ext = _EXT[file_format]
    with tempfile.TemporaryDirectory(dir=tmpdir) as td:
        path = Path(td) / f"roundtrip{ext}"
        mio = mesh_to_meshio(mesh)
        meshio.write(str(path), mio, file_format=file_format)
        mio2 = meshio.read(str(path), file_format=file_format)
        return meshio_to_mesh(mio2, units=mesh.units, coordinate_frame=mesh.coordinate_frame)
