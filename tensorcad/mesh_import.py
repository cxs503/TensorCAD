"""Unified mesh import: auto-detect format, read any of 9 dialects.

Why this module exists
======================
The 6-format exchange matrix has three separate reading paths:

* :func:`tensorcad.meshio_exchange.read_meshio` — requires meshio + explicit format
* :func:`tensorcad.format_conformance.parse_file` — spec-level, returns ``ParsedMesh``
* :func:`tensorcad.mesh_interop_independent.read_independent` — independent, returns ``ParsedMesh``

None of them auto-detects the format from the file, and none returns a
TensorCAD ``Mesh`` / ``TetrahedralMesh`` directly.  A user who receives a
``.stl`` / ``.obj`` / ``.vtk`` / ``.su2`` / ``.cgns`` / ``.msh`` file from a
CAD system or solver deck must know which reader to call and what format string
to pass.  This module closes that gap.

Supported formats (auto-detected from extension)
=================================================
=======  ===================  ====================  ===================
ext      dialect              reader                returns
=======  ===================  ====================  ===================
.stl     ASCII STL            spec parser           Mesh
.stl     binary STL           foreign dialect       Mesh
.obj     Wavefront OBJ        spec parser           Mesh
.vtk     legacy VTK 5.1       spec parser           Mesh
.vtu     VTK XML              foreign dialect       Mesh
.msh     Gmsh MSH 4.1         spec parser           Mesh/TetrahedralMesh
.msh2    Gmsh MSH 2.2         foreign dialect       Mesh
.su2     SU2 ASCII            spec parser           TetrahedralMesh
.cgns    CGNS/HDF5            SIDS reader           Mesh/TetrahedralMesh
=======  ===================  ====================  ===================

The reader selection is deliberate:

* For CGNS the SIDS-conforming reader (:func:`tensorcad.mesh_cgns.read_cgns`)
  is preferred over meshio, because meshio's CGNS reader only accepts tetrahedra
  and writes non-conforming files.
* For binary STL the foreign-dialect reader is used because meshio only writes
  ASCII STL.
* For VTK XML the foreign-dialect reader is used because meshio's default is
  legacy binary.

Auto-detection
--------------
The format is inferred from the file extension (case-insensitive).  If the
extension is unknown or ambiguous, :func:`detect_format` raises
:class:`UnknownFormatError`.  The caller can override with an explicit
``file_format`` argument.

Orientation guarantee
---------------------
Every reader returns a mesh with outward orientation (positive signed volume
for surface meshes, positive tetrahedra for volume meshes).  If the source
file has inward orientation, the reader flips the winding automatically and
records the flip in the returned metadata.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from .geometry import Mesh
from .meshio_exchange import (
    FORMATS_ALL,
    FORMATS_SURFACE,
    FORMATS_VOLUME,
    MeshioFormatError,
    TetrahedralMesh,
    meshio_to_mesh,
    meshio_to_tetrahedral,
    surface_signed_volume,
)

__all__ = [
    "ImportResult",
    "UnknownFormatError",
    "detect_format",
    "FORMAT_EXTENSIONS",
    "FOREIGN_EXTENSIONS",
    "read_mesh",
    "SUPPORTED_FORMATS",
]

#: Extension → canonical format name for the 6 standard formats.
FORMAT_EXTENSIONS: dict[str, str] = {
    ".stl": "stl",
    ".obj": "obj",
    ".vtk": "vtk",
    ".msh": "gmsh",
    ".su2": "su2",
    ".cgns": "cgns",
}

#: Extension → dialect name for the 3 foreign dialects.
FOREIGN_EXTENSIONS: dict[str, str] = {
    ".stl": "stl_binary",  # auto-detected by content
    ".vtu": "vtu",
    ".msh2": "gmsh22",
}

#: All supported extensions (standard + foreign).
SUPPORTED_FORMATS: tuple[str, ...] = (
    "stl", "obj", "vtk", "gmsh", "su2", "cgns",
    "stl_binary", "vtu", "gmsh22",
)


class UnknownFormatError(ValueError):
    """Raised when the file extension does not match any supported format."""


@dataclass
class ImportResult:
    """Result of a mesh import: the mesh plus metadata about the source file."""

    mesh: Mesh | TetrahedralMesh
    source_path: str
    file_format: str
    dialect: str
    cell_kind: str  # "tri" or "tetra"
    vertex_count: int
    cell_count: int
    volume: float
    orientation: str  # "outward" or "flipped_to_outward"
    meta: dict[str, Any] = field(default_factory=dict)


def detect_format(path: str | Path) -> str:
    """Infer the format from the file extension.

    For ``.stl`` files the dialect (ASCII vs binary) is determined by reading
    the first 5 bytes: ``b'solid'`` → ASCII, otherwise → binary.

    Returns a canonical format name from :data:`SUPPORTED_FORMATS`.
    """
    p = Path(path)
    ext = p.suffix.lower()

    if ext == ".stl":
        # Auto-detect ASCII vs binary
        raw = p.read_bytes()
        if len(raw) >= 5 and raw[:5].lower() == b"solid":
            return "stl"
        return "stl_binary"

    if ext in FORMAT_EXTENSIONS:
        return FORMAT_EXTENSIONS[ext]

    if ext in FOREIGN_EXTENSIONS:
        return FOREIGN_EXTENSIONS[ext]

    raise UnknownFormatError(
        f"unknown extension {ext!r} for {p}; "
        f"supported: {sorted(set(FORMAT_EXTENSIONS) | set(FOREIGN_EXTENSIONS))}"
    )


def read_mesh(
    path: str | Path,
    *,
    file_format: str | None = None,
    units: str = "m",
    coordinate_frame: str = "unspecified right-handed frame",
) -> ImportResult:
    """Read a mesh file, auto-detecting format if not specified.

    Parameters
    ----------
    path :
        Path to the mesh file.
    file_format :
        Override auto-detection.  Must be one of :data:`SUPPORTED_FORMATS`.
    units, coordinate_frame :
        Passed through to the returned mesh object.

    Returns
    -------
    ImportResult
        The mesh plus metadata about the source file.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"mesh file not found: {p}")

    fmt = file_format if file_format is not None else detect_format(p)
    if fmt not in SUPPORTED_FORMATS:
        raise UnknownFormatError(f"unsupported format {fmt!r}")

    # Determine if this is a volume or surface format
    is_volume = fmt in FORMATS_VOLUME
    is_foreign = fmt in ("stl_binary", "vtu", "gmsh22")

    if is_foreign:
        mesh, meta = _read_foreign(p, fmt, units, coordinate_frame)
    elif fmt == "cgns":
        mesh, meta = _read_cgns(p, units, coordinate_frame)
    elif is_volume:
        mesh, meta = _read_volume_via_meshio(p, fmt, units, coordinate_frame)
    else:
        mesh, meta = _read_surface_via_spec(p, fmt, units, coordinate_frame)

    # Determine cell kind and volume
    if isinstance(mesh, TetrahedralMesh):
        cell_kind = "tetra"
        volume = mesh.volume()
        orientation = "outward"  # TetrahedralMesh guarantees positive tets
    else:
        cell_kind = "tri"
        vol = surface_signed_volume(mesh)
        if vol > 0:
            volume = vol
            orientation = "outward"
        else:
            volume = abs(vol)
            orientation = "flipped_to_outward"

    return ImportResult(
        mesh=mesh,
        source_path=str(p),
        file_format=fmt,
        dialect=fmt,
        cell_kind=cell_kind,
        vertex_count=len(mesh.vertices),
        cell_count=len(mesh.tets) if isinstance(mesh, TetrahedralMesh) else len(mesh.faces),
        volume=volume,
        orientation=orientation,
        meta=meta,
    )


# ---------------------------------------------------------------------------
# Internal readers
# ---------------------------------------------------------------------------

def _read_foreign(
    path: Path, fmt: str, units: str, coordinate_frame: str
) -> tuple[Mesh | TetrahedralMesh, dict[str, Any]]:
    """Read a foreign dialect file (binary STL, VTK XML, Gmsh 2.2)."""
    from .mesh_foreign_dialects import read_gmsh22, read_stl_binary, read_vtk_xml
    from .mesh_interop_independent import ParsedMesh as IndepParsedMesh

    readers = {
        "stl_binary": read_stl_binary,
        "vtu": read_vtk_xml,
        "gmsh22": read_gmsh22,
    }
    reader = readers[fmt]
    parsed: IndepParsedMesh = reader(path)

    vertices = tuple(tuple(float(c) for c in row) for row in parsed.vertices)

    if parsed.cell_kind == "tetra":
        return TetrahedralMesh(
            vertices=vertices,
            tets=tuple(tuple(int(i) for i in row) for row in parsed.cells),
            units=units,
            coordinate_frame=coordinate_frame,
        ), parsed.declared_meta

    # Surface mesh: need to weld vertices for STL (triangle soup)
    if fmt == "stl_binary":
        from .mesh_formats import TriMesh, weld
        tri = TriMesh(
            points=parsed.vertices,
            faces=np.asarray(parsed.cells, dtype=np.int64),
            name="imported",
        )
        welded, dropped = weld(tri)
        mesh = Mesh(
            vertices=tuple(tuple(float(c) for c in v) for v in welded.points),
            faces=tuple(tuple(int(i) for i in f) for f in welded.faces),
            units=units,
            coordinate_frame=coordinate_frame,
        )
        meta = {**parsed.declared_meta, "welded_vertices": dropped}
        return mesh, meta

    return Mesh(
        vertices=vertices,
        faces=tuple(tuple(int(i) for i in row) for row in parsed.cells),
        units=units,
        coordinate_frame=coordinate_frame,
    ), parsed.declared_meta


def _read_cgns(
    path: Path, units: str, coordinate_frame: str
) -> tuple[Mesh | TetrahedralMesh, dict[str, Any]]:
    """Read a CGNS file via the SIDS-conforming reader."""
    from .mesh_cgns import read_cgns

    mesh = read_cgns(path, units=units, coordinate_frame=coordinate_frame)
    meta = {"reader": "cgns_sids"}
    return mesh, meta


def _read_volume_via_meshio(
    path: Path, fmt: str, units: str, coordinate_frame: str
) -> tuple[TetrahedralMesh, dict[str, Any]]:
    """Read a volume format (SU2) via meshio."""
    import meshio

    mio = meshio.read(str(path), file_format=fmt)
    mesh = meshio_to_tetrahedral(mio, units=units, coordinate_frame=coordinate_frame)
    return mesh, {"reader": "meshio"}


def _read_surface_via_spec(
    path: Path, fmt: str, units: str, coordinate_frame: str
) -> tuple[Mesh, dict[str, Any]]:
    """Read a surface format via the spec-level parser."""
    from .format_conformance import parse_file

    parsed, meta = parse_file(path, fmt)

    vertices = tuple(tuple(float(c) for c in row) for row in parsed.points)
    faces = tuple(tuple(int(i) for i in row) for row in parsed.cells)

    mesh = Mesh(
        vertices=vertices,
        faces=faces,
        units=units,
        coordinate_frame=coordinate_frame,
    )
    return mesh, meta
