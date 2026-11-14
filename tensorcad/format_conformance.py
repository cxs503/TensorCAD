"""Spec-level readers for the six exchange formats, written from the format
specifications and *not* on top of ``meshio`` or of the TensorCAD reader set.

Why this module exists
======================
``tensorcad.meshio_exchange`` proves its matrix by writing a file with meshio and
reading it back with meshio.  That is a self round trip: a writer and a reader
that share a misunderstanding of a format agree with each other and the matrix
still goes green.  The meshio CGNS leg is exactly that case -- both directions
are a ``TODO`` in meshio 5.3.5 (``meshio/cgns/_cgns.py`` stores coordinates in a
child dataset named ``" data"`` and writes no ``ElementType``), so the cell
passes while the file is not a legal CGNS file at all.

This module is the second opinion.  It parses the bytes actually on disk and
reports, for each format:

* points / cells as the file itself declares them (counts taken from the file
  headers, never from the mesh the caller started from),
* the divergence-theorem volume of that geometry,
* whether the surface is closed and consistently oriented,
* a list of format-spec violations (``issues``).

The mesh representation is deliberately dumb (``points`` + integer ``cells``) so
it cannot inherit TensorCAD's conventions.

Formats
-------
========  =====================================================================
``stl``   binary STL (80-byte header, ``uint32`` count, 50 bytes/triangle) and
          ASCII STL
``obj``   Wavefront OBJ (``v`` lines, 1-based ``f`` lines, n-gons fanned)
``vtk``   legacy VTK 5.x unstructured grid / polydata, ASCII and BINARY
          (big-endian data)
``gmsh``  Gmsh MSH 4.1, ASCII and binary (data size 8)
``su2``   SU2 native (``NDIME`` / ``NPOIN`` / ``NELEM`` / ``NMARK``)
``cgns``  CGNS over HDF5 -- both the SIDS layout and meshio's layout
========  =====================================================================

Public API
----------
:class:`ParsedMesh`
    ``points``, ``cells``, ``cell_kind`` plus ``volume()``, ``open_edges()``,
    ``orientation()`` and ``summary()``.
:func:`parse_file`
    Dispatch on ``fmt``; returns ``(ParsedMesh, meta)``.
:func:`cgns_issues`
    CGNS/SIDS conformance check (node labels and types, ``CGNSLibraryVersion``,
    ``ElementType``/``ElementSize``/``ElementRange`` consistency).
"""
from __future__ import annotations

import struct
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

__all__ = [
    "CGNS_ELEMENT_NODES",
    "FormatParseError",
    "GMSH_TYPES",
    "ParsedMesh",
    "SU2_TYPES",
    "cgns_issues",
    "parse_cgns",
    "parse_file",
    "parse_gmsh",
    "parse_obj",
    "parse_stl",
    "parse_su2",
    "parse_vtk",
]


class FormatParseError(ValueError):
    """Raised when a file does not match the syntax of the declared format."""


# ---------------------------------------------------------------------------
# Mesh container + geometric invariants
# ---------------------------------------------------------------------------


@dataclass
class ParsedMesh:
    """A mesh exactly as a file declares it: points plus index cells."""

    points: np.ndarray
    cells: Any
    cell_kind: str  # "triangle" | "polygon" | "tetra"
    meta: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.points = np.asarray(self.points, dtype=np.float64).reshape(-1, 3)
        rows = [np.asarray(c, dtype=np.int64).reshape(-1) for c in self.cells]
        self.cells = np.asarray(rows, dtype=np.int64) if rows else np.zeros((0, 0), np.int64)
        if self.cells.size and (
            int(self.cells.min()) < 0 or int(self.cells.max()) >= len(self.points)
        ):
            raise FormatParseError(
                f"{self.cell_kind} cell index out of range: max {int(self.cells.max())} "
                f"vs {len(self.points)} points"
            )

    # -- derived quantities -------------------------------------------------
    def triangles(self) -> np.ndarray:
        """Facets as an ``(m, 3)`` index array (polygons fan-triangulated)."""
        if self.cell_kind == "triangle":
            return self.cells
        if self.cell_kind == "tetra":
            raise FormatParseError("a tetrahedral mesh has no facet list")
        if self.cell_kind == "polygon":
            tris = [
                (int(row[0]), int(row[k]), int(row[k + 1]))
                for row in self.cells
                for k in range(1, len(row) - 1)
            ]
            return np.asarray(tris, dtype=np.int64).reshape(-1, 3)
        raise FormatParseError(f"unknown cell kind '{self.cell_kind}'")

    def volume(self) -> float:
        """Signed volume: divergence theorem on facets, tet sum on cells."""
        if self.cell_kind == "tetra":
            p = self.points[self.cells]
            a, b, c, d = p[:, 0], p[:, 1], p[:, 2], p[:, 3]
            return float(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a)).sum() / 6.0)
        tris = self.triangles()
        if not len(tris):
            return 0.0
        p = self.points[tris]
        return float(np.einsum("ij,ij->i", p[:, 0], np.cross(p[:, 1], p[:, 2])).sum() / 6.0)

    def area(self) -> float:
        tris = self.triangles()
        if not len(tris):
            return 0.0
        a, b, c = self.points[tris[:, 0]], self.points[tris[:, 1]], self.points[tris[:, 2]]
        return float(0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum())

    def orientation(self) -> str:
        """``"outward"`` (positive signed volume) or ``"inward"``."""
        return "outward" if self.volume() > 0 else "inward"

    def open_edges(self) -> int:
        """Number of undirected edges with valence != 2 (0 means closed)."""
        tris = self.triangles()
        if not len(tris):
            return 0
        edges = np.concatenate([tris[:, [0, 1]], tris[:, [1, 2]], tris[:, [2, 0]]])
        _, counts = np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)
        return int((counts != 2).sum())

    def degenerate_facets(self) -> int:
        """Facets with a repeated index or zero area."""
        tris = self.triangles()
        if not len(tris):
            return 0
        repeated = (np.diff(np.sort(tris, axis=1), axis=1) == 0).any(axis=1)
        a, b, c = self.points[tris[:, 0]], self.points[tris[:, 1]], self.points[tris[:, 2]]
        zero_area = np.linalg.norm(np.cross(b - a, c - a), axis=1) <= 1e-15
        return int((repeated | zero_area).sum())

    def summary(self) -> dict[str, Any]:
        data: dict[str, Any] = {
            "cell_kind": self.cell_kind,
            "points": int(len(self.points)),
            "cells": int(len(self.cells)),
            "volume": self.volume(),
            "area": self.area(),
            "orientation": self.orientation(),
        }
        if self.cell_kind != "tetra":
            data["open_edges"] = self.open_edges()
            data["degenerate_facets"] = self.degenerate_facets()
        data.update(self.meta)
        return data

# ---------------------------------------------------------------------------
# STL
# ---------------------------------------------------------------------------


def parse_stl(path: str | Path) -> tuple[ParsedMesh, dict[str, Any]]:
    """Parse binary or ASCII STL.

    Binary layout: 80-byte header, ``uint32`` triangle count, then 50 bytes per
    facet (3 float32 normal, 9 float32 vertex coordinates, uint16 attribute
    byte count).  The triangle count is cross-checked against the file size --
    a mismatch is a spec violation, not something to guess around.
    """
    data = Path(path).read_bytes()
    expected = None
    if len(data) >= 84:
        (count,) = struct.unpack("<I", data[80:84])
        expected = 84 + 50 * int(count)
    is_binary = expected is not None and expected == len(data)
    if not is_binary and data.lstrip()[:5].lower() != b"solid":
        raise FormatParseError(
            f"STL: size {len(data)} matches neither a binary STL (84+50n={expected}) "
            "nor an ASCII 'solid' header"
        )

    issues: list[str] = []
    if is_binary:
        count = int(struct.unpack("<I", data[80:84])[0])
        if count == 0:
            raise FormatParseError("STL: binary file declares zero triangles")
        rows = np.frombuffer(data, dtype=np.uint8, count=50 * count, offset=84).reshape(count, 50)
        floats = rows[:, :48].copy().view("<f4").reshape(count, 12).astype(np.float64)
        normals = floats[:, 0:3]
        verts = floats[:, 3:12].reshape(count, 3, 3)
        attr = rows[:, 48:50].copy().view("<u2").reshape(count)
        meta: dict[str, Any] = {
            "encoding": "binary",
            "declared_triangles": count,
            "file_bytes": len(data),
            "expected_bytes": expected,
            "bytes_ok": expected == len(data),
            "nonzero_attribute_bytes": int((attr != 0).sum()),
        }
        if meta["nonzero_attribute_bytes"]:
            issues.append(
                f"STL: {meta['nonzero_attribute_bytes']} facets set the 2-byte attribute field"
            )
    else:
        tokens = data.split()
        normals_list: list[list[float]] = []
        verts_list: list[list[list[float]]] = []
        current: list[list[float]] = []
        for i, tok in enumerate(tokens):
            if tok == b"normal":
                normals_list.append(
                    [float(tokens[i + 1]), float(tokens[i + 2]), float(tokens[i + 3])]
                )
            elif tok == b"vertex":
                current.append([float(tokens[i + 1]), float(tokens[i + 2]), float(tokens[i + 3])])
                if len(current) == 3:
                    verts_list.append(current)
                    current = []
        if not verts_list:
            raise FormatParseError("STL: ASCII file contains no facets")
        count = len(verts_list)
        normals = np.asarray(normals_list[:count], dtype=np.float64).reshape(count, 3)
        verts = np.asarray(verts_list, dtype=np.float64)
        meta = {"encoding": "ascii", "declared_triangles": count, "file_bytes": len(data)}

    flat = verts.reshape(-1, 3)
    points, inverse = np.unique(flat, axis=0, return_inverse=True)
    cells = np.asarray(inverse, dtype=np.int64).reshape(count, 3)
    mesh = ParsedMesh(points, cells, "triangle", meta)

    # STL stores one normal per facet; its sign must follow the vertex winding.
    tris = mesh.points[mesh.cells]
    geo = np.cross(tris[:, 1] - tris[:, 0], tris[:, 2] - tris[:, 0])
    good = np.linalg.norm(geo, axis=1) > 0
    agree = np.einsum("ij,ij->i", normals[good], geo[good]) > 0
    meta["normal_present"] = int((np.linalg.norm(normals, axis=1) > 0).sum())
    meta["normal_consistent"] = int(agree.sum())
    meta["nondegenerate_facets"] = int(good.sum())
    if meta["normal_present"] and meta["normal_consistent"] != int(good.sum()):
        issues.append(
            f"STL: {int(good.sum()) - meta['normal_consistent']} facet normals disagree with "
            "the vertex winding order"
        )
    meta["points_deduped"] = int(len(points))
    meta["issues"] = issues
    return mesh, meta


# ---------------------------------------------------------------------------
# OBJ
# ---------------------------------------------------------------------------


def parse_obj(path: str | Path) -> tuple[ParsedMesh, dict[str, Any]]:
    """Parse a Wavefront OBJ: ``v`` lines and 1-based ``f`` lines."""
    points: list[list[float]] = []
    faces: list[list[int]] = []
    issues: list[str] = []
    counts = {"v": 0, "vt": 0, "vn": 0, "f": 0, "other": 0}
    for lineno, raw in enumerate(Path(path).read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        tag = parts[0]
        if tag == "v":
            if len(parts) < 4:
                raise FormatParseError(f"OBJ line {lineno}: 'v' needs 3 coordinates")
            points.append([float(parts[1]), float(parts[2]), float(parts[3])])
            counts["v"] += 1
        elif tag in ("vt", "vn"):
            counts[tag] += 1
        elif tag == "f":
            idx: list[int] = []
            for ref in parts[1:]:
                head = ref.split("/")[0]
                if not head:
                    raise FormatParseError(f"OBJ line {lineno}: facet reference '{ref}' has no vertex")
                i = int(head)
                if i == 0:
                    raise FormatParseError(f"OBJ line {lineno}: facet index 0 is invalid (OBJ is 1-based)")
                idx.append(i - 1 if i > 0 else len(points) + i)
            if len(idx) < 3:
                raise FormatParseError(f"OBJ line {lineno}: facet with {len(idx)} vertices")
            faces.append(idx)
            counts["f"] += 1
        else:
            counts["other"] += 1
    if not faces:
        raise FormatParseError("OBJ: no facets")
    width = {len(f) for f in faces}
    meta: dict[str, Any] = {"counts": counts, "facet_widths": sorted(width)}
    kind = "triangle" if width == {3} else "polygon"
    if kind == "polygon":
        issues.append(f"OBJ: mixed/n-gon facets {sorted(width)} were fan-triangulated")
    meta["issues"] = issues
    return ParsedMesh(points, faces, kind, meta), meta


# ---------------------------------------------------------------------------
# Legacy VTK
# ---------------------------------------------------------------------------


class _Cursor:
    """Byte cursor with the primitives a legacy VTK file needs."""

    def __init__(self, data: bytes) -> None:
        self.data = data
        self.pos = 0

    def line(self) -> str:
        end = self.data.find(b"\n", self.pos)
        if end < 0:
            raise FormatParseError(f"VTK: truncated header at byte {self.pos}")
        text = self.data[self.pos : end].decode("latin-1").strip()
        self.pos = end + 1
        return text

    def peek_line(self) -> str:
        end = self.data.find(b"\n", self.pos)
        chunk = self.data[self.pos : end if end >= 0 else len(self.data)]
        return chunk.decode("latin-1").strip()

    def block(self, dtype: str, count: int) -> np.ndarray:
        size = np.dtype(dtype).itemsize * count
        if self.pos + size > len(self.data):
            raise FormatParseError(
                f"VTK: wanted {size} bytes at offset {self.pos}, file has {len(self.data)}"
            )
        chunk = self.data[self.pos : self.pos + size]
        self.pos += size
        return np.frombuffer(chunk, dtype=dtype, count=count).astype(np.float64)

    def skip_separator(self) -> None:
        while self.pos < len(self.data) and self.data[self.pos : self.pos + 1] in (b"\n", b"\r"):
            self.pos += 1

    def rest_ascii(self, count: int) -> list[str]:
        tokens = self.data[self.pos :].split()[:count]
        self.pos = len(self.data)
        return [t.decode("latin-1") for t in tokens]


#: legacy VTK binary data is big-endian.
_VTK_DTYPE = {
    "double": ">f8",
    "float": ">f4",
    "int": ">i4",
    "long": ">i8",
    "vtktypeint64": ">i8",
    "vtktypeint32": ">i4",
    "vtkidtype": ">i8",
    "char": "i1",
    "unsigned_char": "u1",
}


def parse_vtk(path: str | Path) -> tuple[ParsedMesh, dict[str, Any]]:
    """Parse a legacy VTK unstructured grid (ASCII or BINARY, file version <= 5.1)."""
    data = Path(path).read_bytes()
    cur = _Cursor(data)
    header = cur.line()
    if not header.startswith("# vtk DataFile Version"):
        raise FormatParseError(f"VTK: bad magic line '{header}'")
    version = header.split("Version")[-1].strip()
    cur.line()  # free-form title
    encoding = cur.line().upper()
    if encoding not in ("ASCII", "BINARY"):
        raise FormatParseError(f"VTK: unknown encoding '{encoding}'")
    dataset = cur.line().upper()
    if not dataset.startswith("DATASET"):
        raise FormatParseError(f"VTK: expected DATASET, got '{dataset}'")
    dtype_name = dataset.split()[1].upper()
    if dtype_name not in ("UNSTRUCTURED_GRID", "POLYDATA"):
        raise FormatParseError(f"VTK: unsupported DATASET {dtype_name}")

    points: np.ndarray | None = None
    conn: list[int] = []
    offsets: list[int] = []
    cell_types: list[int] = []
    meta: dict[str, Any] = {"encoding": encoding, "vtk_version": version, "dataset": dtype_name}
    issues: list[str] = []

    while cur.pos < len(data):
        text = cur.line()
        if not text:
            continue
        parts = text.split()
        key = parts[0].upper()
        if key == "POINTS":
            n, kind = int(parts[1]), parts[2]
            if encoding == "ASCII":
                arr = np.asarray([float(t) for t in cur.rest_ascii(3 * n)], dtype=np.float64)
            else:
                arr = cur.block(_VTK_DTYPE[kind], 3 * n)
                cur.skip_separator()
            points = arr.reshape(n, 3)
        elif key == "CELLS":
            n_cells, total = int(parts[1]), int(parts[2])
            if encoding == "ASCII":
                vals = [int(t) for t in cur.rest_ascii(total)]
                pos = 0
                offsets = [0]
                while pos < len(vals):
                    k = vals[pos]
                    conn.extend(vals[pos + 1 : pos + 1 + k])
                    pos += 1 + k
                    offsets.append(pos)
            else:
                head = cur.peek_line().split()
                if head and head[0].upper() == "OFFSETS":
                    cur.line()
                    # In the VTK 5.x OFFSETS/CONNECTIVITY layout the number after
                    # CELLS is already the offset count (cells + 1).
                    offsets = [int(v) for v in cur.block(_VTK_DTYPE[head[1]], n_cells)]
                    cur.skip_separator()
                    head2 = cur.peek_line().split()
                    if not head2 or head2[0].upper() != "CONNECTIVITY":
                        raise FormatParseError("VTK: CELLS without a CONNECTIVITY block")
                    cur.line()
                    conn = [int(v) for v in cur.block(_VTK_DTYPE[head2[1]], total)]
                    cur.skip_separator()
                else:  # VTK <= 5.0: size-prefixed connectivity
                    raw = [int(v) for v in cur.block(">i4", total)]
                    cur.skip_separator()
                    pos = 0
                    offsets = [0]
                    while pos < len(raw):
                        k = raw[pos]
                        conn.extend(raw[pos + 1 : pos + 1 + k])
                        pos += 1 + k
                        offsets.append(pos)
        elif key in ("CONNECTIVITY", "OFFSETS"):
            raise FormatParseError(f"VTK: stray '{key}' block")
        elif key == "CELL_TYPES":
            n = int(parts[1])
            if encoding == "ASCII":
                cell_types = [int(t) for t in cur.rest_ascii(n)]
            else:
                cell_types = [int(v) for v in cur.block(">i4", n)]
                cur.skip_separator()
        elif key == "POLYGONS" and encoding == "ASCII":
            vals = [int(t) for t in cur.rest_ascii(int(parts[2]))]
            pos = 0
            offsets = [0]
            while pos < len(vals):
                k = vals[pos]
                conn.extend(vals[pos + 1 : pos + 1 + k])
                pos += 1 + k
                offsets.append(pos)
            cell_types = cell_types or [7] * int(parts[1])
        elif key in ("POINT_DATA", "CELL_DATA", "VERTICES", "LINES", "POLYGONS"):
            issues.append(f"VTK: '{key}' section not consumed by this audit parser")
            break
        else:
            if encoding == "ASCII":
                continue
            raise FormatParseError(f"VTK: unexpected header '{text}'")

    if points is None:
        raise FormatParseError("VTK: no POINTS block")
    if not offsets:
        raise FormatParseError("VTK: no cell connectivity")
    faces = [conn[offsets[i] : offsets[i + 1]] for i in range(len(offsets) - 1)]

    surface_types = {5}  # VTK_TRIANGLE
    polygonal = {7}  # VTK_POLYGON
    if cell_types and set(cell_types) - surface_types - polygonal:
        keep = [f for f, t in zip(faces, cell_types) if t in surface_types | polygonal]
        meta["dropped_cells"] = len(faces) - len(keep)
        issues.append(f"VTK: dropped {meta['dropped_cells']} non-surface cells")
        faces = keep
    if not faces:
        raise FormatParseError("VTK: no surface cells")
    widths = {len(f) for f in faces}
    kind = "triangle" if widths == {3} else "polygon"
    if kind == "polygon":
        issues.append(f"VTK: facet widths {sorted(widths)} were fan-triangulated")
    meta["n_cell_types"] = len(cell_types) if cell_types else 0
    meta["issues"] = issues
    return ParsedMesh(points, faces, kind, meta), meta


# ---------------------------------------------------------------------------
# Gmsh MSH 4.1
# ---------------------------------------------------------------------------

#: Gmsh element type -> (dimension, nodes per element, name).
GMSH_TYPES: dict[int, tuple[int, int, str]] = {
    1: (0, 2, "line"),
    2: (2, 3, "triangle"),
    3: (2, 4, "quad"),
    4: (3, 4, "tetra"),
    5: (3, 8, "hexahedron"),
    15: (0, 1, "point"),
}


def parse_gmsh(path: str | Path) -> tuple[ParsedMesh, dict[str, Any]]:
    """Parse Gmsh MSH 4.1, ASCII or binary (data size 8).

    In binary 4.1 the section counters are ``size_t`` (8 bytes when the data size
    line says ``8``), element headers are 3 x int32 + size_t, and node
    coordinates follow the node tags of each block.
    """
    data = Path(path).read_bytes()
    if not data.startswith(b"$MeshFormat"):
        raise FormatParseError("Gmsh: file does not start with $MeshFormat")
    end = data.index(b"$EndMeshFormat")
    head = data[len(b"$MeshFormat") : end].split()
    version, binary, data_size = head[0].decode(), int(head[1]), int(head[2])
    if not version.startswith("4.1"):
        raise FormatParseError(f"Gmsh: this parser implements 4.1, file says {version}")
    pos = end + len(b"$EndMeshFormat")
    while pos < len(data) and data[pos : pos + 1] in (b"\n", b"\r"):
        pos += 1
    is_binary = bool(binary)
    size_code = "Q" if data_size == 8 else "I"
    ssize = 8 if data_size == 8 else 4
    meta: dict[str, Any] = {
        "msh_version": version,
        "encoding": "binary" if is_binary else "ascii",
        "data_size": data_size,
    }
    issues: list[str] = []
    node_coords: dict[int, np.ndarray] = {}
    elements: list[tuple[int, list[int]]] = []

    while pos < len(data):
        if data[pos : pos + 1] != b"$":
            pos += 1
            continue
        line_end = data.index(b"\n", pos)
        section = data[pos + 1 : line_end].decode().strip()
        pos = line_end + 1
        stop = data.index(f"$End{section}".encode(), pos)
        if section == "Nodes":
            if is_binary:
                nblocks, nnodes, _tmin, _tmax = struct.unpack(
                    f"<4{size_code}", data[pos : pos + 4 * ssize]
                )
                pos += 4 * ssize
                for _ in range(int(nblocks)):
                    _dim, _tag, _parametric, n = struct.unpack("<iiiQ", data[pos : pos + 20])
                    pos += 20
                    tags = struct.unpack(f"<{int(n)}{size_code}", data[pos : pos + int(n) * ssize])
                    pos += int(n) * ssize
                    xyz = np.frombuffer(data, dtype="<f8", count=3 * int(n), offset=pos).reshape(
                        int(n), 3
                    )
                    pos += 24 * int(n)
                    for tag, row in zip(tags, xyz):
                        node_coords[int(tag)] = np.array(row, dtype=np.float64)
                meta["node_blocks"] = int(nblocks)
                meta["declared_nodes"] = int(nnodes)
            else:
                body = data[pos:stop].split(b"\n")
                nblocks, nnodes = int(body[0].split()[0]), int(body[0].split()[1])
                i = 1
                for _ in range(nblocks):
                    _dim, _tag, _par, n = (int(x) for x in body[i].split())
                    i += 1
                    tags = [int(x) for x in body[i].split()]
                    i += 1
                    for j in range(n):
                        val = [float(x) for x in body[i + j].split()]
                        node_coords[tags[j]] = np.asarray(val[:3], dtype=np.float64)
                    i += n
                meta["declared_nodes"] = nnodes
        elif section == "Elements":
            if is_binary:
                nblocks, nelem, _tmin, _tmax = struct.unpack(
                    f"<4{size_code}", data[pos : pos + 4 * ssize]
                )
                pos += 4 * ssize
                counts: dict[str, int] = {}
                for _ in range(int(nblocks)):
                    _dim, _tag, etype, n = struct.unpack("<iiiQ", data[pos : pos + 20])
                    pos += 20
                    if etype not in GMSH_TYPES:
                        raise FormatParseError(f"Gmsh: unknown element type {etype}")
                    per = GMSH_TYPES[etype][1]
                    for _k in range(int(n)):
                        row = struct.unpack(
                            f"<{per + 1}{size_code}", data[pos : pos + (per + 1) * ssize]
                        )
                        pos += (per + 1) * ssize
                        elements.append((etype, [int(v) for v in row[1:]]))
                    name = GMSH_TYPES[etype][2]
                    counts[name] = counts.get(name, 0) + int(n)
                meta["element_blocks"] = int(nblocks)
                meta["declared_elements"] = int(nelem)
                meta["element_types"] = counts
            else:
                body = data[pos:stop].split(b"\n")
                nblocks, nelem = int(body[0].split()[0]), int(body[0].split()[1])
                i, counts = 1, {}
                for _ in range(nblocks):
                    _dim, _tag, etype, n = (int(x) for x in body[i].split())
                    i += 1
                    for j in range(n):
                        row = [int(x) for x in body[i + j].split()]
                        elements.append((etype, row[1:]))
                    i += n
                    name = GMSH_TYPES[etype][2]
                    counts[name] = counts.get(name, 0) + n
                meta["declared_elements"] = nelem
                meta["element_types"] = counts
        pos = stop + len(f"$End{section}")
        while pos < len(data) and data[pos : pos + 1] in (b"\n", b"\r"):
            pos += 1

    if not node_coords:
        raise FormatParseError("Gmsh: no nodes")
    tags_sorted = sorted(node_coords)
    index = {tag: i for i, tag in enumerate(tags_sorted)}
    points = np.asarray([node_coords[t] for t in tags_sorted], dtype=np.float64)

    surface = [(t, row) for t, row in elements if t in (2, 3)]
    volumes = [(t, row) for t, row in elements if t in (4, 5)]
    if volumes:
        tet_rows = [row for t, row in volumes if t == 4]
        if len(tet_rows) != len(volumes):
            issues.append("Gmsh: non-tetra volume cells present (only tets are read)")
        if not tet_rows:
            raise FormatParseError("Gmsh: volume cells are not tetrahedra")
        cells = [[index[i] for i in row] for row in tet_rows]
        kind = "tetra"
    elif surface:
        cells = [[index[i] for i in row] for _t, row in surface]
        widths = {len(c) for c in cells}
        kind = "triangle" if widths == {3} else "polygon"
        if kind == "polygon":
            issues.append(f"Gmsh: facet widths {sorted(widths)} were fan-triangulated")
    else:
        raise FormatParseError(
            f"Gmsh: no surface/volume elements (types {meta.get('element_types')})"
        )
    meta["issues"] = issues
    return ParsedMesh(points, cells, kind, meta), meta


# ---------------------------------------------------------------------------
# SU2
# ---------------------------------------------------------------------------

#: SU2 element type code -> (name, nodes per element).
SU2_TYPES: dict[int, tuple[str, int]] = {
    3: ("line", 2),
    5: ("triangle", 3),
    9: ("quad", 4),
    10: ("tetra", 4),
    12: ("hexahedron", 8),
    13: ("wedge", 6),
    14: ("pyramid", 5),
}


def parse_su2(path: str | Path) -> tuple[ParsedMesh, dict[str, Any]]:
    """Parse an SU2 native mesh (``NDIME`` / ``NPOIN`` / ``NELEM`` / ``NMARK``)."""
    lines = Path(path).read_text().splitlines()
    meta: dict[str, Any] = {"encoding": "ascii"}
    issues: list[str] = []
    i = 0
    ndime = npoin = nelem = None
    points: list[list[float]] = []
    cells: list[list[int]] = []
    type_count: dict[str, int] = {}
    while i < len(lines):
        raw = lines[i].strip()
        i += 1
        if not raw:
            continue
        if raw.startswith("NDIME"):
            ndime = int(raw.split("=")[1])
        elif raw.startswith("NPOIN"):
            npoin = int(raw.split("=")[1])
            for _ in range(npoin):
                vals = lines[i].split()
                i += 1
                xyz = [float(v) for v in vals[:ndime]]
                while len(xyz) < 3:
                    xyz.append(0.0)
                points.append(xyz)
        elif raw.startswith("NELEM"):
            nelem = int(raw.split("=")[1])
            for _ in range(nelem):
                vals = [int(v) for v in lines[i].split()]
                i += 1
                etype, nodes = vals[0], vals[1:]
                if etype not in SU2_TYPES:
                    raise FormatParseError(f"SU2: unknown element type {etype}")
                name, per = SU2_TYPES[etype]
                if len(nodes) != per:
                    raise FormatParseError(
                        f"SU2: element type {etype} ({name}) carries {len(nodes)} nodes, "
                        f"expected {per}"
                    )
                cells.append(nodes)
                type_count[name] = type_count.get(name, 0) + 1
        elif raw.startswith("NMARK"):
            nmark = int(raw.split("=")[1])
            meta["nmark"] = nmark
            for _ in range(nmark):
                tagline = lines[i].strip()
                i += 1
                if not tagline.upper().startswith("MARKER_TAG"):
                    raise FormatParseError(f"SU2: expected MARKER_TAG, got '{tagline}'")
                meta.setdefault("marker_tags", []).append(tagline.split("=")[1].strip())
                n_elem = int(lines[i].split("=")[1])
                i += 1
                for _k in range(n_elem):
                    vals = [int(v) for v in lines[i].split()]
                    i += 1
                    if vals[0] not in SU2_TYPES:
                        raise FormatParseError(f"SU2: unknown marker element type {vals[0]}")
                    key = "marker_" + SU2_TYPES[vals[0]][0]
                    type_count[key] = type_count.get(key, 0) + 1
        else:
            raise FormatParseError(f"SU2: unexpected line '{raw}'")
    if ndime is None or npoin is None or nelem is None:
        raise FormatParseError("SU2: missing NDIME/NPOIN/NELEM")
    meta.update(
        {
            "ndime": ndime,
            "declared_points": npoin,
            "declared_elements": nelem,
            "element_types": type_count,
        }
    )
    if ndime != 3:
        issues.append(f"SU2: NDIME={ndime} is planar; volume is not defined")
    if len(points) != npoin:
        issues.append(f"SU2: read {len(points)} points but NPOIN={npoin}")
    if len(cells) != nelem:
        issues.append(f"SU2: read {len(cells)} elements but NELEM={nelem}")
    if "nmark" not in meta:
        issues.append("SU2: file has no NMARK block (the format requires one, even if empty)")
    width = {len(c) for c in cells}
    if width == {4}:
        kind = "tetra"
    elif width == {3}:
        kind = "triangle"
    else:
        raise FormatParseError(f"SU2: mixed cell widths {sorted(width)}")
    meta["issues"] = issues
    return ParsedMesh(points, cells, kind, meta), meta


# ---------------------------------------------------------------------------
# CGNS (HDF5)
# ---------------------------------------------------------------------------

#: CGNS SIDS ``ElementType_t`` code -> (name, nodes per element).
CGNS_ELEMENT_NODES: dict[int, tuple[str, int]] = {
    2: ("point", 1),  # NODE
    5: ("triangle", 3),  # TRI_3
    7: ("quad", 4),  # QUAD_4
    10: ("tetra", 4),  # TETRA_4
    11: ("tetra", 10),  # TETRA_10
    17: ("hexahedron", 8),  # HEXA_8
}


def _attr(obj: Any, name: str) -> Any:
    """Read one HDF5 attribute, normalising bytes / ``np.bytes_``."""
    value = obj.attrs.get(name) if hasattr(obj, "attrs") else None
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("latin-1")
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def _data(node: Any) -> np.ndarray | None:
    """Fetch a CGNS data node's array.

    The SIDS HDF5 mapping stores the array in a child dataset literally named
    ``" data"`` (leading space) inside the node group; some writers (meshio) make
    the same layout by hand, others store the node itself as a dataset.
    """
    for child in (" data", "data"):
        if child in node:
            return np.asarray(node[child][()])
    try:
        return np.asarray(node[()])
    except Exception:  # pragma: no cover - group without data
        return None


def _find_zone(base: Any) -> tuple[str, Any] | tuple[None, None]:
    """Locate the zone group inside a base node.

    A SIDS base node also holds ``CellDimension`` and ``PhysicalDimension``
    groups, so "first child group" is not the zone.  Prefer the ``Zone_t``
    label, then a name containing ``Zone``, then any group that looks like a
    zone (holds coordinates or elements).
    """
    children = [(k, base[k]) for k in base if hasattr(base[k], "keys")]
    for name, node in children:
        if _attr(node, "label") == "Zone_t":
            return name, node
    for name, node in children:
        if "zone" in name.lower():
            return name, node
    for name, node in children:
        keys = set(node.keys())
        if any("GridCoordinates" in k or "Element" in k for k in keys):
            return name, node
    return None, None


def parse_cgns(path: str | Path) -> tuple[ParsedMesh, dict[str, Any]]:
    """Parse a CGNS file over HDF5 (SIDS layout, and meshio's layout)."""
    import h5py

    meta: dict[str, Any] = {"encoding": "hdf5"}
    issues: list[str] = []
    with h5py.File(str(path), "r") as handle:
        base_names = [k for k in handle if k != "CGNSLibraryVersion"]
        if not base_names:
            raise FormatParseError("CGNS: file contains no base node")
        base = handle[base_names[0]]
        zone_name, zone = _find_zone(base)
        if zone is None:
            raise FormatParseError("CGNS: base node contains no zone")
        meta["base"] = base_names[0]
        meta["zone"] = zone_name

        coords_name = next(
            (k for k in zone if "GridCoordinates" in k or k == "Coordinates"), None
        )
        if coords_name is None:
            raise FormatParseError("CGNS: zone has no GridCoordinates/Coordinates node")
        coords = zone[coords_name]
        axes = sorted(
            (k for k in coords if k.lower().startswith("coordinate")), key=lambda k: k.lower()
        )
        if len(axes) < 2:
            raise FormatParseError(f"CGNS: fewer than 2 coordinate arrays in {coords_name}")
        columns = []
        for axis in axes[:3]:
            arr = _data(coords[axis])
            if arr is None:
                raise FormatParseError(f"CGNS: {axis} carries no data")
            columns.append(np.asarray(arr, dtype=np.float64).reshape(-1))
        points = np.column_stack(columns)
        if points.shape[1] < 3:
            points = np.column_stack([points, np.zeros(len(points))])
        npoints = len(points)

        elem_names = [k for k in zone if "Element" in k]
        if not elem_names:
            raise FormatParseError("CGNS: zone has no elements node")
        rows: list[list[int]] = []
        kinds: set[str] = set()
        for elem_name in elem_names:
            node = zone[elem_name]
            etype = _attr(node, "ElementType")
            esize = _attr(node, "ElementSize")
            rng_name = next((k for k in node if "ElementRange" in k), None)
            conn_name = next((k for k in node if "ElementConnectivity" in k), None)
            if conn_name is None:
                raise FormatParseError(f"CGNS: {elem_name} has no ElementConnectivity")
            rng = None
            if rng_name is not None:
                raw_rng = _data(node[rng_name])
                if raw_rng is not None:
                    rng = [int(v) for v in np.asarray(raw_rng).reshape(-1)[:2]]
                    meta.setdefault("element_ranges", {})[elem_name] = rng
            conn = np.asarray(_data(node[conn_name]), dtype=np.int64).reshape(-1)
            if conn.size and conn.min() >= 1:
                base_index = 1  # SIDS connectivity is 1-based
            else:
                base_index = 0
                issues.append(f"CGNS: {conn_name} is 0-based (SIDS is 1-based)")
            if etype is not None:
                if int(etype) not in CGNS_ELEMENT_NODES:
                    raise FormatParseError(f"CGNS: unsupported ElementType {int(etype)}")
                name, per = CGNS_ELEMENT_NODES[int(etype)]
                if esize is not None and int(esize) != per:
                    issues.append(
                        f"CGNS: {elem_name} ElementSize={int(esize)} contradicts "
                        f"ElementType={int(etype)} ({per} nodes)"
                    )
            else:
                if not rng:
                    raise FormatParseError(
                        f"CGNS: {elem_name} has neither ElementType nor ElementRange; "
                        "cannot determine the connectivity width"
                    )
                n_elem = rng[1] - rng[0] + 1
                if n_elem <= 0 or conn.size % n_elem:
                    raise FormatParseError(
                        f"CGNS: {elem_name} connectivity of {conn.size} entries does not "
                        f"divide into ElementRange {rng}"
                    )
                per = conn.size // n_elem
                name = {3: "triangle", 4: "tetra", 8: "hexahedron"}.get(per, f"nodes{per}")
                issues.append(
                    f"CGNS: {elem_name} has no ElementType; width {per} inferred from "
                    f"{conn.size} entries over ElementRange {rng}"
                )
            if conn.size % per:
                raise FormatParseError(
                    f"CGNS: {conn_name} has {conn.size} entries, not a multiple of {per}"
                )
            if base_index == 1 and conn.size and conn.max() > npoints:
                raise FormatParseError(
                    f"CGNS: connectivity references node {int(conn.max())} but zone has "
                    f"{npoints} points"
                )
            for row in conn.reshape(-1, per) - base_index:
                rows.append([int(v) for v in row])
            kinds.add(name)
        if not rows:
            raise FormatParseError("CGNS: no element connectivity")
        if len(kinds) > 1:
            raise FormatParseError(f"CGNS: mixed element kinds {sorted(kinds)}")
        kind = kinds.pop()
        if kind not in ("tetra", "triangle"):
            raise FormatParseError(f"CGNS: unsupported element kind '{kind}'")

    meta["issues"] = issues
    return ParsedMesh(points, rows, kind, meta), meta


# ---------------------------------------------------------------------------
# CGNS SIDS conformance
# ---------------------------------------------------------------------------


def cgns_issues(path: str | Path) -> list[str]:
    """Return the CGNS/SIDS conformance violations of ``path`` (empty == clean).

    Checks the parts of the SIDS HDF5 mapping a mesher is responsible for:

    1. ``/CGNSLibraryVersion`` exists at the root;
    2. the base node is labelled ``CGNSBase_t`` and declares ``CellDimension``
       and ``PhysicalDimension``;
    3. the zone is labelled ``Zone_t`` and declares ``ZoneType``;
    4. coordinates live under ``GridCoordinates_t/CoordinateX|Y|Z`` as
       ``DataArray_t`` data nodes in the SIDS ``" data"`` layout;
    5. elements are grouped in ``Elements_t`` nodes that carry ``ElementType``
       and ``ElementSize``, each with its own 1-based ``ElementRange``
       (``IndexRange_t``) and an ``ElementConnectivity`` whose length equals
       ``ElementSize`` x number of elements and whose indices lie in
       ``[1, npoints]``.

    The check is structural: it reads the attribute/node table of the file, it
    cannot prove that a particular third-party CGNS library accepts the file.
    """
    import h5py

    problems: list[str] = []
    with h5py.File(str(path), "r") as handle:
        if "CGNSLibraryVersion" not in handle:
            problems.append("missing /CGNSLibraryVersion (SIDS requires it)")
        elif _attr(handle["CGNSLibraryVersion"], "type") is None:
            problems.append("/CGNSLibraryVersion carries no type attribute")
        base_names = [k for k in handle if k != "CGNSLibraryVersion"]
        if len(base_names) != 1:
            problems.append(f"expected exactly one base node at /, found {base_names}")
        if not base_names:
            return problems
        base = handle[base_names[0]]
        if _attr(base, "label") != "CGNSBase_t":
            problems.append(
                f"/{base_names[0]} label is {_attr(base, 'label')!r}, expected 'CGNSBase_t'"
            )
        for required in ("CellDimension", "PhysicalDimension"):
            if required not in base:
                problems.append(f"base node has no {required}")
        _zone_name, zone = _find_zone(base)
        if zone is None:
            problems.append("base node has no zone")
            return problems
        if _attr(zone, "label") != "Zone_t":
            problems.append(f"zone label is {_attr(zone, 'label')!r}, expected 'Zone_t'")
        if "ZoneType" not in zone:
            problems.append("zone has no ZoneType node")
        coords_name = next((k for k in zone if "GridCoordinates" in k), None)
        npoints = None
        if coords_name is None:
            problems.append("zone has no GridCoordinates_t node")
        else:
            coords = zone[coords_name]
            if _attr(coords, "label") != "GridCoordinates_t":
                problems.append(
                    f"{coords_name} label is {_attr(coords, 'label')!r}, "
                    "expected 'GridCoordinates_t'"
                )
            for axis in ("CoordinateX", "CoordinateY", "CoordinateZ"):
                if axis not in coords:
                    problems.append(f"GridCoordinates has no {axis}")
                    continue
                node = coords[axis]
                if _attr(node, "type") != "DataArray_t":
                    problems.append(
                        f"{axis} type is {_attr(node, 'type')!r}, expected 'DataArray_t'"
                    )
                if " data" not in node:
                    problems.append(f"{axis} has no ' data' dataset (SIDS data layout)")
                arr = _data(node)
                if arr is None:
                    problems.append(f"{axis} carries no data")
                    continue
                if np.asarray(arr).dtype.kind != "f":
                    problems.append(
                        f"{axis} data dtype {np.asarray(arr).dtype} is not floating point"
                    )
                npoints = len(np.asarray(arr).reshape(-1))
        elem_names = [k for k in zone if "Element" in k]
        if not elem_names:
            problems.append("zone has no Elements_t node")
            return problems
        for elem_name in elem_names:
            node = zone[elem_name]
            if _attr(node, "label") != "Elements_t":
                problems.append(
                    f"{elem_name} label is {_attr(node, 'label')!r}, expected 'Elements_t'"
                )
            etype = _attr(node, "ElementType")
            esize = _attr(node, "ElementSize")
            if etype is None:
                problems.append(
                    f"{elem_name} carries no ElementType attribute -- SIDS requires it and "
                    "without it a conforming reader cannot size the connectivity"
                )
            if esize is None:
                problems.append(f"{elem_name} carries no ElementSize attribute")
            rng_name = next((k for k in node if "ElementRange" in k), None)
            conn_name = next((k for k in node if "ElementConnectivity" in k), None)
            if rng_name is None:
                problems.append(f"{elem_name} has no ElementRange child")
            if conn_name is None:
                problems.append(f"{elem_name} has no ElementConnectivity child")
            if rng_name is None or conn_name is None:
                continue
            raw_rng = _data(node[rng_name])
            raw_conn = _data(node[conn_name])
            if raw_rng is None or raw_conn is None:
                problems.append(f"{elem_name} ElementRange/ElementConnectivity carry no data")
                continue
            rng = np.asarray(raw_rng).reshape(-1)
            conn = np.asarray(raw_conn).reshape(-1)
            if rng.size < 2:
                problems.append(f"{elem_name}/ElementRange does not hold [start, end]")
                continue
            start, stop = int(rng[0]), int(rng[1])
            n_elem = stop - start + 1
            if start < 1:
                problems.append(f"{elem_name}/ElementRange starts at {start}; SIDS is 1-based")
            if esize is not None and int(esize) > 0 and conn.size != int(esize) * n_elem:
                problems.append(
                    f"{elem_name}: ElementSize={int(esize)} x {n_elem} elements != "
                    f"{conn.size} connectivity entries"
                )
            if conn.size and int(conn.min()) < 1:
                problems.append(
                    f"{elem_name}/ElementConnectivity is 0-based; SIDS is 1-based"
                )
            if npoints is not None and conn.size and int(conn.max()) > npoints:
                problems.append(
                    f"{elem_name}/ElementConnectivity references node {int(conn.max())} but "
                    f"the zone has {npoints} points"
                )
            if _attr(node[rng_name], "label") != "IndexRange_t":
                problems.append(
                    f"{elem_name}/{rng_name} label is {_attr(node[rng_name], 'label')!r}, "
                    "expected 'IndexRange_t'"
                )
    return problems


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

_PARSERS = {
    "stl": parse_stl,
    "obj": parse_obj,
    "vtk": parse_vtk,
    "gmsh": parse_gmsh,
    "su2": parse_su2,
    "cgns": parse_cgns,
}


def parse_file(path: str | Path, fmt: str) -> tuple[ParsedMesh, dict[str, Any]]:
    """Parse ``path`` with the spec parser for ``fmt``."""
    try:
        parser = _PARSERS[fmt]
    except KeyError:
        raise FormatParseError(
            f"no spec parser for format '{fmt}'; available: {sorted(_PARSERS)}"
        ) from None
    return parser(path)
