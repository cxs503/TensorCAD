"""Foreign mesh dialects: binary STL, Gmsh 2.2 ASCII, VTK XML (``.vtu``).

Why this module exists
======================
The 6-format matrix (:mod:`tensorcad.meshio_exchange`) writes every file with
meshio's *default* driver for that format and reads it back with
:mod:`tensorcad.mesh_interop_independent`.  The dialects meshio's defaults emit
are only a subset of what the six formats look like in the wild:

============  =========================  ==========================
format        dialect this package emits  dialects exchanged in practice
============  =========================  ==========================
STL           ASCII (``solid``)          **binary** (what CAD/CATIA/
                                         SolidWorks/Fusion/FreeCAD and
                                         every printing/solver pipeline
                                         actually writes: 80-byte header,
                                         uint32 facet count, float32)
Gmsh          MSH 4.1 binary             **MSH 2.2 ASCII** (still the
                                         default of most solver decks and
                                         of Gmsh's own ``-format msh2``)
VTK           legacy 5.1 binary          **XML ``.vtu``** (ParaView's
                                         default for unstructured grids)
============  =========================  ==========================

Reading only the emitted dialect means a hull exported by a real CAD system
cannot be ingested at all -- the independent reader raises rather than accept a
file that the rest of the world considers valid STL.  This module adds readers
for those three dialects, plus spec-level writers so the matrix can produce
them without depending on a third-party library (the third-party consumers are
exercised separately, see :func:`tensorcad.mesh_interop_thirdparty.trimesh_probe`).

Every reader parses the bytes from the format specification alone: no meshio, no
``h5py``, no trimesh.  Anything outside the documented dialect raises
:class:`tensorcad.mesh_interop_independent.IndependentReadError` -- an
unparsed file must never be reported as verified.

Precision note (matters for the acceptance table)
-------------------------------------------------
Binary STL stores positions as IEEE-754 **float32** (the format has no float64
option), so a binary-STL round trip cannot be exact: the volume error is
bounded by ~``3 * 2**-24`` relative, i.e. ``~2e-7``.  The 18-combination matrix
used to carry a blanket ``1e-3`` STL tolerance described as "a leftover from an
assumption that STL stores float32".  Here the bound is derived instead of
guessed: :data:`STL_BINARY_RELATIVE_TOLERANCE` and the campaign records the
measured error, so a loss beyond float32 would still be visible.
"""
from __future__ import annotations

import base64
import struct
import zlib
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

from .mesh_interop_independent import (
    IndependentReadError,
    ParsedMesh,
    _cells_array,
)

__all__ = [
    "STL_BINARY_RELATIVE_TOLERANCE",
    "read_stl_binary",
    "read_gmsh22",
    "read_vtk_xml",
    "write_stl_binary",
    "write_gmsh22",
    "write_vtu_ascii",
    "foreign_identity",
    "FOREIGN_DIALECTS",
]

#: relative volume tolerance justified by binary STL's float32 positions:
#: ``3 * eps32 = 3 * 1.19e-7``, doubled for the accumulation of 3 coordinates
#: per vertex over the divergence-theorem sum.
STL_BINARY_RELATIVE_TOLERANCE = 1e-6

FOREIGN_DIALECTS = ("stl_binary", "gmsh22", "vtu")

_STL_RECORD = np.dtype(
    [
        ("normal", "<f4", (3,)),
        ("vertices", "<f4", (3, 3)),
        ("attribute", "<u2"),
    ]
)


# ---------------------------------------------------------------------------
# geometry helper shared by the writers
# ---------------------------------------------------------------------------

def _mesh_arrays(mesh: Any) -> tuple[np.ndarray, tuple[tuple[int, ...], ...], str]:
    """Return ``(vertices, cells, cell_kind)`` from a TensorCAD mesh object."""
    vertices = np.asarray(getattr(mesh, "vertices"), dtype=np.float64)
    if hasattr(mesh, "tets"):
        cells = tuple(tuple(int(i) for i in t) for t in mesh.tets)
        kind = "tetra"
    elif hasattr(mesh, "faces"):
        cells = tuple(tuple(int(i) for i in f) for f in mesh.faces)
        kind = "tri"
    else:  # pragma: no cover - defensive
        raise TypeError("mesh has neither .tets nor .faces")
    if vertices.ndim != 2 or vertices.shape[1] != 3:
        raise ValueError(f"vertices must be (n, 3), got {vertices.shape}")
    widths = {len(c) for c in cells}
    if len(widths) != 1:
        raise ValueError(f"mixed cell widths {sorted(widths)} are not supported")
    return vertices, cells, kind


def _check_width(kind: str, width: int, path: Path) -> None:
    expected = {"tri": 3, "tetra": 4}[kind]
    if width != expected:
        raise IndependentReadError(
            f"{path}: {width}-node elements are not {kind} ({expected} nodes)"
        )


# ---------------------------------------------------------------------------
# binary STL
# ---------------------------------------------------------------------------

def read_stl_binary(path: str | Path) -> ParsedMesh:
    """Parse a binary STL (80-byte header, uint32 LE count, 50-byte facets).

    The stored per-facet normal is kept: it is redundant with the winding and a
    writer that emits the two inconsistently produces a file that renders one
    way and meshes the other, so the caller can compare them.
    """
    p = Path(path)
    raw = p.read_bytes()
    if len(raw) < 84:
        raise IndependentReadError(f"{p}: too short for a binary STL ({len(raw)} bytes)")
    head = raw[:80].rstrip(b"\x00").lstrip()
    if head[:5].lower() == b"solid" and b"facet" in raw[:1024]:
        raise IndependentReadError(
            f"{p}: this is an ASCII STL, not binary; use read_stl_ascii"
        )
    (count,) = struct.unpack("<I", raw[80:84])
    available = (len(raw) - 84) // _STL_RECORD.itemsize
    if count == 0:
        raise IndependentReadError(f"{p}: binary STL declares 0 facets")
    if count > available:
        raise IndependentReadError(
            f"{p}: binary STL declares {count} facets but only {available} are present"
        )
    records = np.frombuffer(raw, dtype=_STL_RECORD, count=count, offset=84)
    points = records["vertices"].reshape(-1, 3).astype(np.float64)
    cells = [(3 * i, 3 * i + 1, 3 * i + 2) for i in range(count)]
    return ParsedMesh(
        file_format="stl_binary",
        vertices=points,
        cells=_cells_array(cells),
        cell_kind="tri",
        declared_meta={
            "facet_count": int(count),
            "normals": [tuple(float(x) for x in n) for n in records["normal"]],
            "header": raw[:80].rstrip(b"\x00").decode("ascii", "replace"),
            "position_dtype": "float32",
            "decimal_precision": 7,
        },
    )


def write_stl_binary(
    path: str | Path,
    mesh: Any,
    *,
    header: str = "tensorcad binary STL",
) -> Path:
    """Write ``mesh`` as a binary STL from the specification.

    Positions are cast to float32 because that is the only precision the format
    has; the cast is deliberate and measured, never silent (see
    :data:`STL_BINARY_RELATIVE_TOLERANCE`).

    The facet normal is computed **after** the float32 cast, from the
    coordinates that are actually stored.  Computing it from the float64
    positions (which is what ``meshio``'s writer does) leaves the stored normal
    describing a triangle that is not in the file: for geometry far from the
    origin the two disagree by up to ~0.45 deg at ``x ~ 1e4 m``, which is
    visible to renderers and to readers that trust the normal instead of the
    winding.  See ``docs/m23297_stl_encoding_verification.md``.
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    vertices, cells, kind = _mesh_arrays(mesh)
    if kind != "tri":
        raise ValueError("binary STL carries triangles only")
    _check_width(kind, len(cells[0]), p)
    tri32 = vertices[np.asarray(cells, dtype=np.int64)].astype(np.float32)
    tri = tri32.astype(np.float64)  # normals from the *quantised* geometry
    normals = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    lengths = np.linalg.norm(normals, axis=1)
    lengths[lengths == 0.0] = 1.0
    normals = normals / lengths[:, None]
    records = np.zeros(len(cells), dtype=_STL_RECORD)
    records["normal"] = normals.astype(np.float32)
    records["vertices"] = tri32
    raw_header = header.encode("ascii")[:80].ljust(80, b"\x00")
    with p.open("wb") as fh:
        fh.write(raw_header)
        fh.write(struct.pack("<I", len(cells)))
        fh.write(records.tobytes())
    return p


# ---------------------------------------------------------------------------
# Gmsh MSH 2.2 (ASCII)
# ---------------------------------------------------------------------------

_GMSH22_TYPE_WIDTH = {2: 3, 4: 4}
_GMSH22_TYPE_NAME = {2: "triangle", 4: "tetra"}
_GMSH22_KIND_TYPE = {"tri": 2, "tetra": 4}


def read_gmsh22(path: str | Path) -> ParsedMesh:
    """Parse an ASCII MSH 2.2 file (``$Nodes`` / ``$Elements`` lists).

    This is the dialect Gmsh writes with ``-format msh2`` and the one most
    solver decks ship, so it is the common case for third-party files.
    """
    p = Path(path)
    text = p.read_text(errors="replace")
    if not text.startswith("$MeshFormat"):
        raise IndependentReadError(f"{p}: no $MeshFormat section")
    lines = [ln.strip() for ln in text.splitlines()]
    try:
        version_line = lines[lines.index("$MeshFormat") + 1]
    except (ValueError, IndexError) as exc:  # pragma: no cover - defensive
        raise IndependentReadError(f"{p}: truncated $MeshFormat section") from exc
    parts = version_line.split()
    if not parts or not parts[0].startswith("2."):
        raise IndependentReadError(
            f"{p}: MSH version {version_line!r} is not 2.x "
            "(MSH 4.x is read by read_gmsh41)"
        )
    if len(parts) > 1 and parts[1] not in ("0", "1"):
        raise IndependentReadError(f"{p}: unknown MSH file-type field {parts[1]!r}")
    if len(parts) > 1 and parts[1] == "1":
        raise IndependentReadError(f"{p}: MSH {parts[0]} binary is not parsed here")

    def _section(name: str) -> list[str]:
        start = f"${name}"
        end = f"$End{name}"
        try:
            i = lines.index(start) + 1
            j = lines.index(end, i)
        except ValueError as exc:
            raise IndependentReadError(f"{p}: missing {start} section") from exc
        return [ln for ln in lines[i:j] if ln]

    node_lines = _section("Nodes")
    if not node_lines:
        raise IndependentReadError(f"{p}: empty $Nodes section")
    declared_nodes = int(node_lines[0].split()[0])
    coords_by_tag: dict[int, tuple[float, float, float]] = {}
    for line in node_lines[1:]:
        fields = line.split()
        if len(fields) < 4:
            raise IndependentReadError(f"{p}: short node record {line!r}")
        coords_by_tag[int(fields[0])] = (float(fields[1]), float(fields[2]), float(fields[3]))
    if declared_nodes != len(coords_by_tag):
        raise IndependentReadError(
            f"{p}: declared {declared_nodes} nodes, found {len(coords_by_tag)}"
        )

    element_lines = _section("Elements")
    if not element_lines:
        raise IndependentReadError(f"{p}: empty $Elements section")
    declared_elements = int(element_lines[0].split()[0])
    cells: list[tuple[int, ...]] = []
    widths: set[int] = set()
    seen_tags: list[int] = []
    skipped = 0
    for line in element_lines[1:]:
        fields = line.split()
        if len(fields) < 3:
            raise IndependentReadError(f"{p}: short element record {line!r}")
        tag, etype, ntags = int(fields[0]), int(fields[1]), int(fields[2])
        seen_tags.append(tag)
        width = _GMSH22_TYPE_WIDTH.get(etype)
        if width is None:
            # MSH 2.2 legitimately carries points/lines/quads alongside the
            # bulk cells; skipping them is not a loss of the bulk topology, but
            # it is counted so the report can state what was ignored.
            skipped += 1
            continue
        nodes = fields[3 + ntags :]
        if len(nodes) != width:
            raise IndependentReadError(
                f"{p}: element {tag} declares type {etype} but has {len(nodes)} nodes"
            )
        try:
            cells.append(tuple(int(n) for n in nodes))
        except ValueError as exc:
            raise IndependentReadError(f"{p}: non-integer node tag in {line!r}") from exc
        widths.add(width)
    if not cells:
        raise IndependentReadError(f"{p}: no triangle/tetra elements")
    if len(widths) != 1:
        raise IndependentReadError(f"{p}: mixed bulk cell widths {sorted(widths)}")
    width = widths.pop()
    kind = {3: "tri", 4: "tetra"}[width]
    if len(seen_tags) != declared_elements:
        raise IndependentReadError(
            f"{p}: declared {declared_elements} elements, found {len(seen_tags)}"
        )
    if seen_tags != list(range(1, len(seen_tags) + 1)):
        raise IndependentReadError(f"{p}: element tags do not run 1..N")
    tags = sorted(coords_by_tag)
    index = {tag: i for i, tag in enumerate(tags)}
    try:
        cells_idx = [tuple(index[t] for t in c) for c in cells]
    except KeyError as exc:
        raise IndependentReadError(f"{p}: element references unknown node tag {exc}") from exc
    return ParsedMesh(
        file_format="gmsh22",
        vertices=np.asarray([coords_by_tag[t] for t in tags], dtype=np.float64),
        cells=_cells_array(cells_idx),
        cell_kind=kind,
        declared_meta={
            "version": parts[0],
            "element_tag_range": [min(seen_tags), max(seen_tags)],
            "node_count_declared": declared_nodes,
            "element_count_declared": declared_elements,
            "cells_skipped_non_bulk": skipped,
        },
    )


def write_gmsh22(path: str | Path, mesh: Any, *, physical_tag: int = 1) -> Path:
    """Write ASCII MSH 2.2 (``$Nodes``/``$Elements``, two tags per element)."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    vertices, cells, kind = _mesh_arrays(mesh)
    etype = _GMSH22_KIND_TYPE[kind]
    with p.open("w") as fh:
        fh.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        fh.write(f"$Nodes\n{len(vertices)}\n")
        for i, (x, y, z) in enumerate(vertices, start=1):
            fh.write(f"{i} {x:.17g} {y:.17g} {z:.17g}\n")
        fh.write("$EndNodes\n")
        fh.write(f"$Elements\n{len(cells)}\n")
        for i, cell in enumerate(cells, start=1):
            nodes = " ".join(str(int(n) + 1) for n in cell)
            fh.write(f"{i} {etype} 2 {physical_tag} {physical_tag} {nodes}\n")
        fh.write("$EndElements\n")
    return p


# ---------------------------------------------------------------------------
# VTK XML UnstructuredGrid (.vtu)
# ---------------------------------------------------------------------------

_VTU_CELL_KIND = {5: "tri", 10: "tetra"}
_VTU_KIND_CELL = {v: k for k, v in _VTU_CELL_KIND.items()}
#: names an XML writer may use for the zlib compressor (root attribute is the
#: VTK class name, the per-array attribute is the short name)
_VTU_ZLIB_NAMES = frozenset(("zlib", "vtkzlib", "vtkzlibdatacompressor"))

_VTU_DTYPE = {
    "Float32": np.float32,
    "Float64": np.float64,
    "Int32": np.int32,
    "Int64": np.int64,
    "UInt8": np.uint8,
    "UInt32": np.uint32,
    "UInt64": np.uint64,
}


def _b64_chars(num_bytes: int) -> int:
    """Number of base64 characters an ``num_bytes``-byte chunk encodes to."""
    return -(-num_bytes // 3) * 4


def _decode_compressed(
    text: str, path: Path, *, header_dtype: str = "UInt32"
) -> bytes:
    """Decode VTK compressed inline binary (``vtkZLibDataCompressor``).

    Layout, read off the bytes rather than assumed (VTK writes it this way and
    so does meshio 5.x): the header is one base64 chunk of
    ``3 + nr_blocks`` unsigned integers -- ``nr_blocks``, ``block_size``,
    ``last_block_size``, then the **compressed** size of every block -- and the
    blocks themselves follow as a *single* base64 chunk, each compressed
    independently with zlib.

    The first version of this function guessed "base64 chunk per block, with a
    size prefix per block" and returned a 2-value array for an 8-point mesh
    instead of failing.  That silent short read is why every length here is
    cross-checked against the header.
    """
    item = {"uint32": 4, "uint64": 8}[(header_dtype or "UInt32").lower()]
    fmt = "<I" if item == 4 else "<Q"
    if len(text) < _b64_chars(item):
        raise IndependentReadError(f"{path}: compressed VTU array shorter than one header item")
    (num_blocks,) = struct.unpack(fmt, base64.b64decode(text[: _b64_chars(item)])[:item])
    if num_blocks == 0:
        raise IndependentReadError(f"{path}: compressed VTU array declares 0 blocks")
    header_bytes = item * (3 + int(num_blocks))
    header_chars = _b64_chars(header_bytes)
    if len(text) < header_chars:
        raise IndependentReadError(
            f"{path}: compressed VTU header needs {header_chars} base64 chars, {len(text)} present"
        )
    header = np.frombuffer(base64.b64decode(text[:header_chars]), dtype=fmt)
    num_blocks, block_size, last_block_size = (int(x) for x in header[:3])
    sizes = [int(x) for x in header[3:]]
    if len(sizes) != num_blocks:
        raise IndependentReadError(
            f"{path}: {num_blocks} blocks declared but {len(sizes)} sizes in the header"
        )
    payload = base64.b64decode(text[header_chars:])
    total = sum(sizes)
    if total != len(payload):
        raise IndependentReadError(
            f"{path}: compressed VTU blocks declare {total} bytes, {len(payload)} present"
        )
    out = bytearray()
    offset = 0
    for index, size in enumerate(sizes):
        try:
            out.extend(zlib.decompress(payload[offset : offset + size]))
        except zlib.error as exc:
            raise IndependentReadError(
                f"{path}: compressed VTU block {index} is not a zlib stream ({exc})"
            ) from exc
        offset += size
    expected = block_size * (num_blocks - 1) + last_block_size
    if len(out) != expected:
        raise IndependentReadError(
            f"{path}: compressed VTU blocks decompressed to {len(out)} bytes, "
            f"header declares {expected}"
        )
    return bytes(out)


def _vtu_array(
    node: Any,
    path: Path,
    *,
    root_compressor: str = "",
    root_header_type: str = "UInt32",
) -> np.ndarray:
    """Decode one ``DataArray`` (ascii, inline binary, or compressed binary)."""
    dtype = _VTU_DTYPE.get(node.get("type", "Float64"))
    if dtype is None:
        raise IndependentReadError(f"{path}: unsupported VTU array type {node.get('type')!r}")
    fmt = (node.get("format") or "ascii").lower()
    text = (node.text or "").strip()
    if fmt == "ascii":
        if not text:
            raise IndependentReadError(f"{path}: empty ascii VTU array {node.get('Name')!r}")
        return np.fromstring(text, sep=" ").astype(dtype)  # type: ignore[deprecated]
    if fmt == "appended":
        raise IndependentReadError(
            f"{path}: VTU appended data (offset={node.get('offset')!r}) is not parsed here"
        )
    if fmt != "binary":
        raise IndependentReadError(f"{path}: unknown VTU encoding {fmt!r}")
    compressor = (node.get("compressor") or root_compressor or "").lower()
    if compressor and compressor not in _VTU_ZLIB_NAMES:
        raise IndependentReadError(f"{path}: unsupported VTU compressor {compressor!r}")
    if compressor:
        raw = _decode_compressed(text, path, header_dtype=root_header_type)
    else:
        raw = base64.b64decode(text)
        itemsize = np.dtype(dtype).itemsize
        # VTK inline binary may carry a length header; trust it only when the
        # declared count matches what is actually there, and never fall through
        # to a short array (that is how a wrong decode looks like a small mesh).
        for header_len, fmt_char in ((4, "<I"), (8, "<Q")):
            if len(raw) > header_len and (len(raw) - header_len) % itemsize == 0:
                (declared,) = struct.unpack_from(fmt_char, raw, 0)
                if declared == len(raw) - header_len:
                    raw = raw[header_len:]
                    break
    if not raw or len(raw) % np.dtype(dtype).itemsize:
        raise IndependentReadError(
            f"{path}: VTU array {node.get('Name')!r} payload is {len(raw)} bytes, "
            f"not a whole number of {np.dtype(dtype).name}"
        )
    return np.frombuffer(raw, dtype=dtype).copy()


def read_vtk_xml(path: str | Path) -> ParsedMesh:
    """Parse a VTK XML ``UnstructuredGrid`` (``.vtu``), ascii or inline binary."""
    from xml.etree import ElementTree

    p = Path(path)
    try:
        root = ElementTree.fromstring(p.read_bytes())
    except ElementTree.ParseError as exc:
        raise IndependentReadError(f"{p}: not well-formed XML ({exc})") from exc
    if root.tag != "VTKFile":
        raise IndependentReadError(f"{p}: root element is {root.tag!r}, not 'VTKFile'")
    if root.get("type") != "UnstructuredGrid":
        raise IndependentReadError(f"{p}: VTKFile type {root.get('type')!r} is not UnstructuredGrid")
    piece = root.find(".//Piece")
    if piece is None:
        raise IndependentReadError(f"{p}: no <Piece>")

    root_compressor = (root.get("compressor") or "").lower()
    root_header_type = root.get("header_type") or "UInt32"

    def _named(parent: Any, name: str) -> np.ndarray:
        for node in parent.findall(".//DataArray"):
            if node.get("Name") == name:
                return _vtu_array(
                    node, p, root_compressor=root_compressor, root_header_type=root_header_type
                )
        raise IndependentReadError(f"{p}: no DataArray named {name!r}")

    n_points = int(piece.get("NumberOfPoints", 0))
    n_cells = int(piece.get("NumberOfCells", 0))
    points_node = piece.find("Points")
    if points_node is None:
        raise IndependentReadError(f"{p}: no <Points>")
    points = _vtu_array(
        points_node.find("DataArray"),
        p,
        root_compressor=root_compressor,
        root_header_type=root_header_type,
    )
    if points.size != 3 * n_points:
        raise IndependentReadError(
            f"{p}: {points.size} point coordinates decoded for {n_points} points"
        )
    points = np.asarray(points, dtype=np.float64).reshape(n_points, 3)
    cells_node = piece.find("Cells")
    if cells_node is None:
        raise IndependentReadError(f"{p}: no <Cells>")
    connectivity = np.asarray(_named(cells_node, "connectivity"), dtype=np.int64)
    offsets = np.asarray(_named(cells_node, "offsets"), dtype=np.int64)
    types = np.asarray(_named(cells_node, "types")).astype(np.int64)
    if types.size != n_cells or offsets.size != n_cells:
        raise IndependentReadError(
            f"{p}: {n_cells} cells declared, {types.size} types / {offsets.size} offsets"
        )
    kinds = {_VTU_CELL_KIND.get(int(t)) for t in types}
    if None in kinds or len(kinds) != 1:
        raise IndependentReadError(f"{p}: cell types {sorted(set(types.tolist()))} are not one of tri/tetra")
    kind = kinds.pop()
    if offsets[-1] != connectivity.size:
        raise IndependentReadError(
            f"{p}: last offset {int(offsets[-1])} != connectivity length {connectivity.size}"
        )
    starts = np.concatenate(([0], offsets[:-1]))
    widths = offsets - starts
    _check_width(kind, int(widths[0]), p)
    if not np.all(widths == widths[0]):
        raise IndependentReadError(f"{p}: mixed cell widths {sorted(set(widths.tolist()))}")
    cells = [tuple(int(i) for i in connectivity[s:e]) for s, e in zip(starts, offsets)]
    return ParsedMesh(
        file_format="vtu",
        vertices=points,
        cells=_cells_array(cells),
        cell_kind=kind,
        declared_meta={
            "version": root.get("version"),
            "byte_order": root.get("byte_order"),
            "encoding": (piece.find(".//DataArray") is not None)
            and (piece.find(".//DataArray").get("format") or "ascii"),
            "vtk_cell_type": _VTU_KIND_CELL[kind],
            "point_count_declared": n_points,
            "cell_count_declared": n_cells,
        },
    )


def write_vtu_ascii(path: str | Path, mesh: Any) -> Path:
    """Write an ASCII VTK XML ``UnstructuredGrid`` from the specification."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    vertices, cells, kind = _mesh_arrays(mesh)
    offsets = np.cumsum([len(c) for c in cells])
    connectivity = np.concatenate([np.asarray(c, dtype=np.int64) for c in cells])
    types = np.full(len(cells), _VTU_KIND_CELL[kind], dtype=np.uint8)

    def _fmt(arr: np.ndarray) -> str:
        return " ".join(f"{x:.17g}" for x in arr.reshape(-1))

    with p.open("w") as fh:
        fh.write('<?xml version="1.0"?>\n')
        fh.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        fh.write("  <UnstructuredGrid>\n")
        fh.write(f'    <Piece NumberOfPoints="{len(vertices)}" NumberOfCells="{len(cells)}">\n')
        fh.write("      <Points>\n")
        fh.write(
            '        <DataArray type="Float64" Name="Points" NumberOfComponents="3" '
            f'format="ascii">{_fmt(vertices)}</DataArray>\n'
        )
        fh.write("      </Points>\n")
        fh.write("      <Cells>\n")
        fh.write(
            '        <DataArray type="Int64" Name="connectivity" format="ascii">'
            f"{_fmt(connectivity)}</DataArray>\n"
        )
        fh.write(
            '        <DataArray type="Int64" Name="offsets" format="ascii">'
            f"{_fmt(offsets)}</DataArray>\n"
        )
        fh.write(
            '        <DataArray type="UInt8" Name="types" format="ascii">'
            f"{_fmt(types)}</DataArray>\n"
        )
        fh.write("      </Cells>\n")
        fh.write("    </Piece>\n")
        fh.write("  </UnstructuredGrid>\n")
        fh.write("</VTKFile>\n")
    return p


# ---------------------------------------------------------------------------
# byte-level identity facts for the foreign dialects
# ---------------------------------------------------------------------------

def foreign_identity(path: str | Path, *, dialect: str) -> dict[str, Any]:
    """Facts about the bytes that a conforming reader of ``dialect`` relies on."""
    p = Path(path)
    raw = p.read_bytes()
    if dialect == "stl_binary":
        head = raw[:80]
        count = struct.unpack("<I", raw[80:84])[0] if len(raw) >= 84 else None
        return {
            "header_len": len(head),
            "header_is_ascii": head[:5].lower() != b"solid",
            "declared_facets": count,
            "bytes_expected": None if count is None else 84 + 50 * count,
            "bytes_actual": len(raw),
            "position_width_bits": 32,
        }
    if dialect == "gmsh22":
        return {
            "has_meshformat": raw.startswith(b"$MeshFormat"),
            "version_line": raw.split(b"\n", 2)[1].decode("ascii", "replace")
            if raw.startswith(b"$MeshFormat")
            else None,
            "has_end_nodes": b"$EndNodes" in raw,
            "has_end_elements": b"$EndElements" in raw,
            "is_text": all(b < 128 for b in raw[:4096]),
        }
    if dialect == "vtu":
        head = raw[:200]
        return {
            "xml_declaration": head.lstrip().startswith(b"<?xml"),
            "vtkfile_tag": b"<VTKFile" in head,
            "type_declared": b'type="UnstructuredGrid"' in head,
            "is_text": all(b < 128 for b in raw[:4096]),
        }
    raise IndependentReadError(f"no foreign identity reader for dialect {dialect!r}")

# ---------------------------------------------------------------------------
# The foreign-dialect campaign (single source of truth for the second table)
# ---------------------------------------------------------------------------

#: dialects the campaign writes, with their producers and precision class
DIALECT_PRODUCERS: tuple[tuple[str, str, str, str], ...] = (
    # (dialect, producer, value_precision, producer_note)
    ("stl_binary", "tensorcad", "float32", "write_stl_binary (this module)"),
    ("stl_binary", "trimesh", "float32", "trimesh binary STL exporter"),
    ("gmsh22", "tensorcad", "float64", "write_gmsh22 (this module)"),
    ("gmsh22", "meshio", "float64", "meshio 'gmsh22' ASCII driver"),
    ("vtu", "tensorcad", "float64", "write_vtu_ascii (this module)"),
    ("vtu", "meshio", "float64", "meshio 'vtu' driver (inline base64 binary)"),
)

_DIALECT_EXT = {"stl_binary": "stl", "gmsh22": "msh", "vtu": "vtu"}


def _produce(dialect: str, producer: str, path: Path, vertices: np.ndarray, cells, kind: str) -> Path:
    """Write one file with the requested producer."""
    if producer == "tensorcad":
        payload = _SurfaceLike(vertices, cells, kind)
        writer = {
            "stl_binary": write_stl_binary,
            "gmsh22": write_gmsh22,
            "vtu": write_vtu_ascii,
        }[dialect]
        return writer(path, payload)
    if producer == "trimesh":
        from .mesh_interop_thirdparty import trimesh_write

        return trimesh_write(path, vertices, cells, file_type="stl")
    if producer == "meshio":
        import meshio

        vtk_type = {"tri": ("triangle", 3), "tetra": ("tetra", 4)}[kind][0]
        m = meshio.Mesh(
            np.asarray(vertices, dtype=np.float64),
            [(vtk_type, np.asarray(cells, dtype=np.int64))],
        )
        kwargs = {"binary": False} if dialect == "gmsh22" else {}
        meshio.write(str(path), m, file_format=dialect if dialect != "vtu" else "vtu", **kwargs)
        return path
    raise ValueError(f"unknown producer {producer!r}")


class _SurfaceLike:
    """Minimal mesh object the spec writers accept (``.vertices`` + cells)."""

    def __init__(self, vertices: np.ndarray, cells, kind: str):
        self.vertices = np.asarray(vertices, dtype=np.float64)
        if kind == "tetra":
            self.tets = cells
        else:
            self.faces = cells


def run_foreign_matrix(
    geometries: Sequence[str] = ("cube", "capsule", "wigley_hull"),
    *,
    workdir: str | Path | None = None,
    include_trimesh: bool = True,
) -> list[dict[str, Any]]:
    """Write/read every (geometry, dialect, producer) and measure the exchange.

    Invariants per row: volume vs the reference mesh, cell count, per-cell
    winding signature.  ``float32`` dialects are measured against a
    float32-quantised reference so their tolerance reflects the format's
    precision instead of hiding a genuine loss.
    """
    import tempfile

    from .mesh_interop_independent import (
        ANALYTIC_VOLUME_3,
        campaign_mesh,
        file_identity,
        oriented_signature,
        read_independent,
        surface_signed_volume,
        tetrahedral_volume,
    )
    from .meshio_exchange import tetrahedralise_surface

    base = Path(workdir) if workdir is not None else Path(tempfile.mkdtemp(prefix="foreigndialects-"))
    base.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for geom_name in geometries:
        surface = campaign_mesh(geom_name)
        for dialect, producer, precision, note in DIALECT_PRODUCERS:
            if producer == "trimesh" and not include_trimesh:
                continue
            mesh = surface
            kind = "tri"
            cells_src = surface.faces
            vol_ref = surface_signed_volume(
                np.asarray(surface.vertices, dtype=np.float64), cells_src
            )
            rows.append(
                _measure(
                    base, geom_name, dialect, producer, precision, note,
                    mesh, cells_src, kind, vol_ref,
                )
            )
        # volume-cell dialects: gmsh22/vtu must also carry tetrahedra
        tet = tetrahedralise_surface(surface)
        vol_tet = tetrahedral_volume(
            np.asarray(tet.vertices, dtype=np.float64), tet.tets
        )
        for dialect in ("gmsh22", "vtu"):
            rows.append(
                _measure(
                    base, geom_name, dialect, "tensorcad", "float64",
                    "tetra cells (this module)", tet, tet.tets, "tetra", vol_tet,
                )
            )
    return rows


def _measure(
    base: Path,
    geom_name: str,
    dialect: str,
    producer: str,
    precision: str,
    note: str,
    mesh: Any,
    cells_src,
    kind: str,
    vol_ref: float,
) -> dict[str, Any]:
    from .mesh_interop_independent import (
        ANALYTIC_VOLUME_3,
        file_identity,
        oriented_signature,
        read_independent,
    )

    ext = _DIALECT_EXT[dialect]
    tag = f"{geom_name}.{producer}.{dialect}"
    path = base / f"{tag}.{ext}"
    reference_vertices = np.asarray(mesh.vertices, dtype=np.float64)
    _produce(dialect, producer, path, reference_vertices, cells_src, kind)
    parsed = read_independent(path, file_format=dialect)
    vol_ind = parsed.volume()
    rel = abs(vol_ind - vol_ref) / max(abs(vol_ref), 1e-300)

    # float32 dialects are compared against the float32 image of the reference:
    # otherwise the signature check would report the format's precision as a
    # geometry change.
    if precision == "float32":
        quant = reference_vertices.astype(np.float32).astype(np.float64)
        tolerance = STL_BINARY_RELATIVE_TOLERANCE
    else:
        quant = reference_vertices
        tolerance = 1e-12
    sig_ref = oriented_signature(quant, cells_src, cell_kind=kind)
    sig_ind = oriented_signature(parsed.vertices, parsed.cells, cell_kind=parsed.cell_kind)
    row: dict[str, Any] = {
        "geometry": geom_name,
        "dialect": dialect,
        "producer": producer,
        "producer_note": note,
        "precision": precision,
        "cell_kind": parsed.cell_kind,
        "cells": len(parsed.cells),
        "cells_reference": len(cells_src),
        "vertices": int(parsed.vertices.shape[0]),
        "vertices_reference": len(reference_vertices),
        "volume_reference": vol_ref,
        "volume_independent": vol_ind,
        "volume_relative_error": rel,
        "analytic_volume": ANALYTIC_VOLUME_3.get(geom_name),
        "cell_count_preserved": len(parsed.cells) == len(cells_src),
        "winding_preserved": sig_ref == sig_ind,
        "identity": file_identity(path, file_format=dialect),
        "declared_meta": parsed.declared_meta,
        "bytes": int(Path(path).stat().st_size),
        "volume_tolerance": tolerance,
    }
    if dialect == "stl_binary":
        from .mesh_interop_independent import stl_normal_consistency

        try:
            row["stl_normals"] = stl_normal_consistency(parsed)
        except Exception:  # pragma: no cover - defensive
            row["stl_normals"] = None
    row["ok"] = bool(
        row["cell_count_preserved"] and row["winding_preserved"] and rel <= tolerance
    )
    return row


def trimesh_consumer_rows(
    geometries: Sequence[str] = ("cube", "capsule", "wigley_hull"),
    *,
    workdir: str | Path | None = None,
) -> list[dict[str, Any]]:
    """TensorCAD-written files loaded by trimesh, volume compared to reference.

    This is the direction the 18-combination matrix could not evidence at all:
    a genuinely foreign implementation doing the reading.
    """
    import tempfile

    from .mesh_interop_independent import (
        ANALYTIC_VOLUME_3,
        campaign_mesh,
        surface_signed_volume,
    )
    from .mesh_interop_thirdparty import trimesh_available, trimesh_probe

    if not trimesh_available():  # pragma: no cover - environment dependent
        return []
    base = Path(workdir) if workdir is not None else Path(tempfile.mkdtemp(prefix="trimeshrows-"))
    base.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for geom_name in geometries:
        surface = campaign_mesh(geom_name)
        vertices = np.asarray(surface.vertices, dtype=np.float64)
        cells = surface.faces
        vol_ref = surface_signed_volume(vertices, cells)
        for dialect in ("stl", "stl_binary", "obj"):
            tag = f"{geom_name}.consumer.{dialect}"
            ext = "stl" if dialect.startswith("stl") else "obj"
            path = base / f"{tag}.{ext}"
            if dialect == "stl":
                from .meshio_exchange import write_meshio

                write_meshio(path, _SurfaceLike(vertices, cells, "tri"), file_format="stl")
            elif dialect == "stl_binary":
                write_stl_binary(path, _SurfaceLike(vertices, cells, "tri"))
            else:
                from .meshio_exchange import write_meshio

                write_meshio(path, _SurfaceLike(vertices, cells, "tri"), file_format="obj")
            probe = trimesh_probe(path)
            rel = abs(probe["volume"] - vol_ref) / max(abs(vol_ref), 1e-300)
            tolerance = STL_BINARY_RELATIVE_TOLERANCE if dialect == "stl_binary" else 1e-12
            row = {
                "geometry": geom_name,
                "dialect": dialect,
                "consumer": "trimesh",
                "consumer_version": probe["tool_version"],
                "volume_reference": vol_ref,
                "volume_thirdparty": probe["volume"],
                "volume_relative_error": rel,
                "analytic_volume": ANALYTIC_VOLUME_3.get(geom_name),
                "thirdparty_watertight": probe["watertight"],
                "thirdparty_winding_consistent": probe["winding_consistent"],
                "thirdparty_faces": probe["faces"],
                "thirdparty_euler": probe["euler_number"],
                "volume_tolerance": tolerance,
            }
            row["ok"] = bool(
                rel <= tolerance
                and probe["watertight"]
                and probe["winding_consistent"]
                and probe["faces"] == len(cells)
            )
            rows.append(row)
    return rows


def foreign_matrix_markdown(rows: Sequence[dict[str, Any]]) -> str:
    """Render :func:`run_foreign_matrix` + consumer rows as one table."""
    lines = [
        "| geometry | dialect | cells | producer | consumer | cells | V [m^3] | V ref [m^3] | rel. err | tol | winding | count | ok |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            "| {g} | {d} | {k} | {p} | {c} | {n}/{nr} | {v:.12f} | {vr:.12f} | {rel:.2e} | "
            "{tol:.0e} | {w} | {cc} | {ok} |".format(
                g=r["geometry"],
                d=r["dialect"],
                k=r.get("cell_kind", "tri"),
                p=r.get("producer", "-"),
                c=r.get("consumer", "independent reader"),
                n=r.get("cells", r.get("thirdparty_faces")),
                nr=r.get("cells_reference", r.get("thirdparty_faces")),
                v=r.get("volume_independent", r.get("volume_thirdparty")),
                vr=r["volume_reference"],
                rel=r["volume_relative_error"],
                tol=r["volume_tolerance"],
                w="PASS" if r.get("winding_preserved", r.get("thirdparty_winding_consistent")) else "FAIL",
                cc="PASS" if r.get("cell_count_preserved", True) else "FAIL",
                ok="PASS" if r["ok"] else "FAIL",
            )
        )
    return "\n".join(lines) + "\n"
