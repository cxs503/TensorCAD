"""CGNS export in the SIDS HDF5 mapping (and a reader that accepts it).

Why this module exists
======================
The CGNS leg of :mod:`tensorcad.meshio_exchange` is delegated to
``meshio.cgns``.  That driver is a placeholder: it writes no
``CGNSLibraryVersion``, no node labels/types, no ``CellDimension`` /
``PhysicalDimension`` / ``ZoneType``, and no ``ElementType`` / ``ElementSize``
on the elements node -- and it stores connectivity 1-based while omitting the
attribute that tells a reader the width.  A file like that is only readable by
meshio itself, which is why the 18-cell matrix could go green on a file that no
CGNS-conforming tool would accept (see
:func:`tensorcad.format_conformance.cgns_issues`, which reports 14 violations
for a ``meshio``-written file).

What this module writes
=======================
The layout of CGNS 3.x/4.x over HDF5 (SIDS):

* ``/CGNSLibraryVersion`` (``R4``, ``CNST``);
* ``/Base`` labelled ``CGNSBase_t`` with ``CellDimension`` and
  ``PhysicalDimension``;
* ``/Base/Zone1`` labelled ``Zone_t`` with a ``ZoneType`` of ``Unstructured``
  and a ``GridCoordinates_t`` holding ``CoordinateX|Y|Z`` as ``DataArray_t``;
* ``/Base/Zone1/GridElements`` labelled ``Elements_t`` carrying ``ElementType``
  (``TRI_3`` = 5 / ``TETRA_4`` = 10) and ``ElementSize``, with 1-based
  ``ElementRange`` (``IndexRange_t``) and ``ElementConnectivity`` children.

Every node carries the three attributes a SIDS reader walks the tree with:
``label``, ``type`` and ``flags``.

The zone and elements node names stay ``Zone1`` / ``GridElements`` on purpose:
that keeps the file readable by ``meshio`` as well, so switching the CGNS
exporter to this module does not cost compatibility -- it only *adds* the
metadata a conforming reader needs.  Conformance is structural: the checks in
:func:`tensorcad.format_conformance.cgns_issues` verify the node/attribute
table, not that a particular third-party CGNS library accepts the file.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .format_conformance import CGNS_ELEMENT_NODES, parse_cgns

__all__ = [
    "CGNS_LIBRARY_VERSION",
    "read_cgns",
    "write_cgns",
]

#: Version stamped into ``/CGNSLibraryVersion``.
CGNS_LIBRARY_VERSION = 4.0

#: SIDS ``ElementType_t`` codes this exporter can write.
_ELEMENT_TYPE_CODE = {"triangle": 5, "tetra": 10}


def _node(group: Any, name: str, label: str, node_type: str, *, constant: bool = False) -> Any:
    """Create one SIDS node (an HDF5 group plus its ``label``/``type``/``flags``)."""
    node = group.create_group(name)
    node.attrs["label"] = np.bytes_(label)
    node.attrs["type"] = np.bytes_(node_type)
    node.attrs["flags"] = np.array([0, 0], dtype=np.int32)
    if constant:
        node.attrs["name_type"] = np.bytes_("CNST")
    return node


def _data_array(
    parent: Any,
    name: str,
    values: Any,
    *,
    label: str | None = None,
    dtype: str | None = None,
    constant: bool = False,
) -> Any:
    """Create a ``DataArray_t`` child holding ``values`` in the SIDS ``" data"`` slot."""
    node = _node(parent, name, label or name, "DataArray_t", constant=constant)
    array = np.asarray(values, dtype=dtype)
    node.create_dataset(
        " data",
        data=array,
        compression="gzip",
        compression_opts=4,
    )
    return node


def write_cgns(
    path: str | Path,
    mesh: Any,
    *,
    base_name: str = "Base",
    zone_name: str = "Zone1",
    elements_name: str = "GridElements",
) -> Path:
    """Write ``mesh`` as a CGNS/SIDS file over HDF5.

    Parameters
    ----------
    path :
        Destination ``.cgns`` file.
    mesh :
        A surface :class:`tensorcad.geometry.Mesh` (written as ``TRI_3``) or a
        :class:`tensorcad.meshio_exchange.TetrahedralMesh` (``TETRA_4``).  Any
        object exposing ``vertices`` plus either ``faces`` or ``tets`` works.
    base_name, zone_name, elements_name :
        Node names.  The defaults are the ones ``meshio`` looks for, so the
        output stays readable by meshio as well.

    Returns
    -------
    Path
        The file that was written.
    """
    import h5py

    from .geometry import Mesh
    from .meshio_exchange import TetrahedralMesh

    if isinstance(mesh, TetrahedralMesh):
        kind, cells = "tetra", mesh.tets
    elif isinstance(mesh, Mesh):
        kind, cells = "triangle", mesh.faces
        # A surface CGNS file needs no tet topology; triangles are TRI_3.
    else:
        if hasattr(mesh, "tets"):
            kind, cells = "tetra", mesh.tets
        elif hasattr(mesh, "faces"):
            kind, cells = "triangle", mesh.faces
        else:
            raise TypeError(
                "write_cgns needs a mesh with 'tets' or 'faces'; "
                f"got {type(mesh).__name__}"
            )

    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    conn = np.asarray(cells, dtype=np.int32)
    if not len(conn):
        raise ValueError("write_cgns: mesh has no cells")
    per_element = conn.shape[1]
    expected = CGNS_ELEMENT_NODES[_ELEMENT_TYPE_CODE[kind]][1]
    if per_element != expected:
        raise ValueError(
            f"write_cgns: {kind} cells must have {expected} nodes, mesh has {per_element}"
        )
    if conn.min() < 0 or conn.max() >= len(vertices):
        raise ValueError("write_cgns: cell indices out of range")

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(str(target), "w") as handle:
        version = handle.create_dataset(
            "CGNSLibraryVersion", data=np.array([CGNS_LIBRARY_VERSION], dtype=np.float32)
        )
        version.attrs["label"] = np.bytes_("CGNSLibraryVersion_t")
        version.attrs["type"] = np.bytes_("R4")
        version.attrs["flags"] = np.array([0, 0], dtype=np.int32)
        version.attrs["name_type"] = np.bytes_("CNST")

        base = _node(handle, base_name, "CGNSBase_t", "CGNSBase_t")
        _data_array(base, "CellDimension", [3, 3, 0], label="CellDimension", dtype=np.int32)
        _data_array(base, "PhysicalDimension", [3, 3, 0], label="PhysicalDimension", dtype=np.int32)

        zone = _node(base, zone_name, "Zone_t", "Zone_t")
        zone_type = _node(zone, "ZoneType", "ZoneType", "ZoneType_t", constant=True)
        zone_type.create_dataset(" data", data=np.array([b"Unstructured"], dtype="S13"))

        coords = _node(zone, "GridCoordinates", "GridCoordinates_t", "GridCoordinates_t")
        for axis, column in zip("XYZ", range(3)):
            _data_array(
                coords,
                f"Coordinate{axis}",
                vertices[:, column],
                label=f"Coordinate{axis}",
                dtype=np.float64,
            )

        elements = _node(zone, elements_name, "Elements_t", "Elements_t")
        elements.attrs["ElementType"] = np.int32(_ELEMENT_TYPE_CODE[kind])
        elements.attrs["ElementSize"] = np.int32(per_element)
        _data_array(
            elements,
            "ElementRange",
            np.array([1, len(conn)], dtype=np.int64),
            label="IndexRange_t",
            dtype=np.int64,
        )
        _data_array(
            elements,
            "ElementConnectivity",
            (conn.reshape(-1) + 1).astype(np.int32),
            label="ElementConnectivity",
            dtype=np.int32,
        )

    return target


def read_cgns(
    path: str | Path,
    *,
    units: str = "m",
    coordinate_frame: str = "unspecified right-handed frame",
) -> Any:
    """Read a CGNS file written by :func:`write_cgns` or by any SIDS writer.

    Accepts both the SIDS layout and meshio's (``" data"`` child with no
    ``ElementType``); the parsing itself lives in
    :func:`tensorcad.format_conformance.parse_cgns`, which is written from the
    format specification rather than on top of meshio.

    Returns a :class:`tensorcad.geometry.Mesh` for triangle data and a
    :class:`tensorcad.meshio_exchange.TetrahedralMesh` for tetrahedra.
    """
    from .geometry import Mesh
    from .meshio_exchange import TetrahedralMesh

    parsed, _meta = parse_cgns(path)
    vertices = tuple(tuple(float(c) for c in row) for row in parsed.points)
    if parsed.cell_kind == "tetra":
        return TetrahedralMesh(
            vertices=vertices,
            tets=tuple(tuple(int(i) for i in row) for row in parsed.cells),
            units=units,
            coordinate_frame=coordinate_frame,
        )
    if parsed.cell_kind == "triangle":
        return Mesh(
            vertices=vertices,
            faces=tuple(tuple(int(i) for i in row) for row in parsed.cells),
            units=units,
            coordinate_frame=coordinate_frame,
        )
    raise ValueError(f"read_cgns: unsupported cell kind '{parsed.cell_kind}'")