"""Independent (non-meshio) readers for the 6-format mesh exchange matrix.

Why this module exists
======================
``tests/test_meshio_format_matrix.py`` proves round-trip fidelity by writing
with meshio and reading back *with meshio*.  That is a self-consistency check,
not an interoperability check: a writer and a reader that share a private
convention (a wrong cell-type code, a swapped axis, a reordered connectivity
block, a file that is not really the format its extension claims) agree with
each other and the matrix still looks green.  This module closes that gap by
parsing the on-disk bytes from the format specifications alone -- no meshio,
no ``tensorcad.mesh_cgns`` -- so a green table means *another tool can read the
file and gets the same geometry*.

Supported dialects (those actually produced by
:func:`tensorcad.meshio_exchange.write_meshio`):

===========  =========================================================
``stl``      ASCII STL (``solid``/``facet normal``/``outer loop``)
``obj``      Wavefront OBJ, ``v`` / ``f`` records, 1-based indices
``vtk``      legacy VTK 5.1, ``BINARY``; ``OFFSETS``/``CONNECTIVITY``
``gmsh``     MSH 4.1, binary, ``$Nodes``/``$Elements`` entity blocks
``su2``      SU2 ASCII (``NDIME``/``NPOIN``/``NELEM``, VTK type codes)
``cgns``     CGNS/HDF5 read through ``h5py`` SIDS node paths
===========  =========================================================

Anything else raises :class:`IndependentReadError` explicitly rather than
guessing -- an unreachable format must never be reported as verified.

Orientation invariant
=====================
The package's own matrix checks only ``sign(volume)``.  A closed surface can
keep the sign of its total volume while individual facets flip (two opposite
flips cancel).  :func:`oriented_signature` therefore compares *per-cell*
winding: every cell is keyed by its sorted vertex coordinates and carries the
sign of its orientation, so a mesh whose facets/tets were re-wound -- or whose
vertex numbering changed, as STL legitimately does -- is detected.
"""
from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

__all__ = [
    "IndependentReadError",
    "ParsedMesh",
    "read_stl_ascii",
    "read_obj",
    "read_vtk_legacy",
    "read_gmsh41",
    "read_su2",
    "read_cgns_hdf5",
    "read_independent",
    "surface_signed_volume",
    "tetrahedral_volume",
    "oriented_signature",
    "file_identity",
]


class IndependentReadError(Exception):
    """The file cannot be parsed as the format its name claims."""


@dataclass(frozen=True)
class ParsedMesh:
    """A mesh as recovered from the file bytes by an independent parser."""

    file_format: str
    vertices: np.ndarray            # (n, 3) float64
    cells: tuple[tuple[int, ...], ...]
    cell_kind: str                  # "triangle" | "tetra"
    declared_meta: dict[str, Any] = field(default_factory=dict)

    def volume(self) -> float:
        if self.cell_kind == "tri":
            return surface_signed_volume(self.vertices, self.cells)
        return tetrahedral_volume(self.vertices, self.cells)


def _sids_data(node: Any) -> Any:
    """Return the SIDS array of a CGNS node.

    CGNS/SIDS keeps the values in a child node whose name is ``" data"`` (a
    leading blank is part of the name), so a plain ``node["data"]`` lookup
    fails on conforming files while a raw ``node[...]`` works on datasets.
    """
    if hasattr(node, "keys"):
        for key in node.keys():
            if key.strip() == "data":
                return node[key][...]
        raise IndependentReadError(f"CGNS node {node.name!r} has no ' data' child")
    return node[...]


def _cells_array(cells: Iterable[Sequence[int]]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(int(i) for i in c) for c in cells)


# ---------------------------------------------------------------------------
# Geometry invariants (independent of the writer's own volume helpers)
# ---------------------------------------------------------------------------

def surface_signed_volume(vertices: np.ndarray, triangles: Iterable[Sequence[int]]) -> float:
    """Divergence-theorem volume ``sum (a . (b x c)) / 6`` of a closed surface."""
    v = np.asarray(vertices, dtype=np.float64)
    tris = np.asarray(list(triangles), dtype=np.int64)
    if tris.size == 0:
        return 0.0
    a, b, c = v[tris[:, 0]], v[tris[:, 1]], v[tris[:, 2]]
    return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)


def tetrahedral_volume(vertices: np.ndarray, tets: Iterable[Sequence[int]]) -> float:
    """``sum det(b-a, c-a, d-a) / 6`` over positively and negatively signed cells."""
    v = np.asarray(vertices, dtype=np.float64)
    t = np.asarray(list(tets), dtype=np.int64)
    if t.size == 0:
        return 0.0
    a = v[t[:, 0]]
    m = np.stack((v[t[:, 1]] - a, v[t[:, 2]] - a, v[t[:, 3]] - a), axis=1)
    return float(np.linalg.det(m).sum() / 6.0)


def _key_coords(points: np.ndarray) -> str:
    """Order-independent coordinate key; 12-decimal rounding (all formats here
    are ASCII or float64, so no legitimate loss needs smoothing over)."""
    unique = sorted(
        (round(float(p[0]), 12), round(float(p[1]), 12), round(float(p[2]), 12))
        for p in points
    )
    return "|".join(f"{x!r},{y!r},{z!r}" for x, y, z in unique)


def oriented_signature(
    vertices: np.ndarray,
    cells: Iterable[Sequence[int]],
    *,
    cell_kind: str,
) -> Counter:
    """Per-cell winding signature: ``Counter{(coordinate key, sign): count}``.

    The key ignores vertex numbering (STL is a triangle soup), the sign records
    which way the cell is wound, so two signatures compare equal exactly when
    every cell survived with the same geometry *and* the same orientation.
    """
    v = np.asarray(vertices, dtype=np.float64)
    signature: Counter = Counter()
    for cell in cells:
        idx = [int(i) for i in cell]
        pts = v[idx]
        key = _key_coords(pts)
        if cell_kind == "tri":
            normal = np.cross(pts[1] - pts[0], pts[2] - pts[0])
            norm = float(np.linalg.norm(normal))
            if norm == 0.0:
                sign = 0
            else:
                # quantised unit normal: a re-wound facet gets the opposite sign
                sign = tuple(np.round(normal / norm, 9).tolist())
        else:
            det = float(np.linalg.det(np.stack((pts[1] - pts[0], pts[2] - pts[0], pts[3] - pts[0]))))
            sign = 0 if det == 0.0 else (1 if det > 0 else -1)
        signature[(key, sign)] += 1
    return signature


# ---------------------------------------------------------------------------
# Format identity (magic numbers / required records)
# ---------------------------------------------------------------------------

def file_identity(path: str | Path, *, file_format: str) -> dict[str, Any]:
    """Facts about the bytes that a format-conforming reader relies on."""
    p = Path(path)
    raw = p.read_bytes()
    if file_format == "cgns":
        return {
            "hdf5_magic": raw[:8] == b"\x89HDF\r\n\x1a\n",
            "has_cgns_library_version": b"CGNSLibraryVersion" in raw,
        }
    head = raw[:512]
    if file_format == "stl":
        return {
            "starts_solid": head.lstrip().lower().startswith(b"solid"),
            "has_facet": b"facet normal" in head or b"facet" in raw[:4096],
            "ends_endsolid": b"endsolid" in raw[-256:].lower(),
        }
    if file_format == "obj":
        return {
            "has_v": b"\nv " in raw or raw.startswith(b"v "),
            "has_f": b"\nf " in raw or raw.startswith(b"f "),
        }
    if file_format == "vtk":
        first = raw.split(b"\n", 1)[0]
        return {
            "header": first.decode("ascii", "replace"),
            "is_legacy_vtk": first.startswith(b"# vtk DataFile Version"),
            "has_unstructured_grid": b"DATASET UNSTRUCTURED_GRID" in raw,
        }
    if file_format == "gmsh":
        version_line = b""
        if raw.startswith(b"$MeshFormat"):
            version_line = raw.split(b"\n", 2)[1]
        return {
            "has_meshformat": raw.startswith(b"$MeshFormat"),
            "version_line": version_line.decode("ascii", "replace"),
            "has_end_elements": b"$EndElements" in raw,
        }
    if file_format == "su2":
        return {
            "has_ndime": raw.startswith(b"NDIME="),
            "has_npoin": b"NPOIN=" in raw,
            "has_nelem": b"NELEM=" in raw,
        }
    raise IndependentReadError(f"no independent reader for format '{file_format}'")


# ---------------------------------------------------------------------------
# STL (ASCII)
# ---------------------------------------------------------------------------

def read_stl_ascii(path: str | Path) -> ParsedMesh:
    p = Path(path)
    raw = p.read_bytes()
    if p.read_bytes()[:5].lower() != b"solid":
        raise IndependentReadError(f"{p}: not an ASCII STL (no 'solid' header)")
    text = raw.decode("ascii")
    vertices: list[tuple[float, float, float]] = []
    cells: list[tuple[int, int, int]] = []
    normals: list[tuple[float, float, float]] = []
    current: list[int] = []
    stated: tuple[float, float, float] | None = None
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "facet" and len(parts) >= 5 and parts[1] == "normal":
            stated = (float(parts[2]), float(parts[3]), float(parts[4]))
            current = []
        elif parts[0] == "vertex" and len(parts) >= 4:
            if stated is None:
                raise IndependentReadError(f"{p}: 'vertex' outside a facet block")
            vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
            current.append(len(vertices) - 1)
        elif parts[0] == "endfacet":
            if len(current) != 3:
                raise IndependentReadError(
                    f"{p}: facet with {len(current)} vertices (expected 3)"
                )
            cells.append((current[0], current[1], current[2]))
            assert stated is not None
            normals.append(stated)
            stated = None
    if not cells:
        raise IndependentReadError(f"{p}: no facets found")
    # a degenerate STL duplicates vertices per facet; that is legal but then the
    # soup volume is still well defined, so no de-duplication is needed here.
    return ParsedMesh(
        file_format="stl",
        vertices=np.asarray(vertices, dtype=np.float64),
        cells=_cells_array(cells),
        cell_kind="tri",
        declared_meta={"facet_count": len(cells), "normals": normals},
    )


# ---------------------------------------------------------------------------
# OBJ
# ---------------------------------------------------------------------------

def read_obj(path: str | Path) -> ParsedMesh:
    p = Path(path)
    vertices: list[tuple[float, float, float]] = []
    cells: list[tuple[int, ...]] = []
    for line in p.read_text().splitlines():
        parts = line.split()
        if not parts or parts[0].startswith("#"):
            continue
        if parts[0] == "v":
            if len(parts) < 4:
                raise IndependentReadError(f"{p}: short vertex record {line!r}")
            vertices.append((float(parts[1]), float(parts[2]), float(parts[3])))
        elif parts[0] == "f":
            nodes = []
            for token in parts[1:]:
                # face records may carry v/vt/vn; the vertex index is the first field
                nodes.append(int(token.split("/")[0]))
            if len(nodes) != 3:
                raise IndependentReadError(
                    f"{p}: face with {len(nodes)} nodes is not a triangle: {line!r}"
                )
            if any(i == 0 for i in nodes):
                raise IndependentReadError(f"{p}: OBJ indices are 1-based, found 0")
            cells.append(tuple(i - 1 for i in nodes))
    if not vertices or not cells:
        raise IndependentReadError(f"{p}: OBJ has no v/f records")
    return ParsedMesh(
        file_format="obj",
        vertices=np.asarray(vertices, dtype=np.float64),
        cells=_cells_array(cells),
        cell_kind="tri",
        declared_meta={"vertex_records": len(vertices)},
    )


# ---------------------------------------------------------------------------
# legacy VTK (ASCII or BINARY)
# ---------------------------------------------------------------------------

_VTK_DTYPE = {
    "double": np.float64,
    "float": np.float32,
    "vtktypeint32": np.int32,
    "vtktypeint64": np.int64,
    "int": np.int32,
    "long": np.int64,
}


def _vtk_line(raw: bytes, pos: int) -> tuple[str, int]:
    end = raw.index(b"\n", pos)
    return raw[pos:end].decode("ascii"), end + 1


def _vtk_array(raw: bytes, pos: int, count: int, dtype_name: str) -> tuple[np.ndarray, int]:
    try:
        dtype = _VTK_DTYPE[dtype_name.lower()]
    except KeyError as exc:  # pragma: no cover - unsupported dialect
        raise IndependentReadError(f"unsupported VTK array type '{dtype_name}'") from exc
    # legacy VTK binary arrays are big-endian (spec; meshio forces ">i8" too)
    dtype = np.dtype(dtype).newbyteorder(">")
    nbytes = count * np.dtype(dtype).itemsize
    end = pos + nbytes
    if raw[end:end + 1] == b"\n":  # binary blocks are newline-terminated
        end += 1
    return np.frombuffer(raw, dtype=dtype, count=count, offset=pos), end


def read_vtk_legacy(path: str | Path) -> ParsedMesh:
    p = Path(path)
    raw = p.read_bytes()
    header, pos = _vtk_line(raw, 0)
    if not header.startswith("# vtk DataFile Version"):
        raise IndependentReadError(f"{p}: {header!r} is not a legacy VTK header")
    _title, pos = _vtk_line(raw, pos)
    encoding, pos = _vtk_line(raw, pos)
    encoding = encoding.strip().upper()
    if encoding not in ("ASCII", "BINARY"):
        raise IndependentReadError(f"{p}: unknown VTK encoding {encoding!r}")
    dataset, pos = _vtk_line(raw, pos)
    if "UNSTRUCTURED_GRID" not in dataset.upper():
        raise IndependentReadError(f"{p}: dataset {dataset!r} is not UNSTRUCTURED_GRID")

    if encoding == "BINARY":
        points_at = raw.find(b"POINTS ", pos)
        if points_at < 0:
            raise IndependentReadError(f"{p}: no POINTS section")
        points_hdr, pos = _vtk_line(raw, points_at)
        _, n_pts, points_type = points_hdr.split()
        n_pts = int(n_pts)
        pts, pos = _vtk_array(raw, pos, 3 * n_pts, points_type)

        cells_at = raw.find(b"CELLS ", pos)
        if cells_at < 0:
            raise IndependentReadError(f"{p}: no CELLS section")
        cells_hdr, pos = _vtk_line(raw, cells_at)
        _cells_parts = cells_hdr.split()
        # VTK 5.1 (as ParaView writes it) declares the number of *offsets*
        # (= cells + 1) first, then the connectivity length; meshio follows.
        n_offsets, n_conn = int(_cells_parts[1]), int(_cells_parts[2])
        n_cells = n_offsets - 1
        if raw[pos:].startswith(b"OFFSETS"):
            offsets_hdr, pos = _vtk_line(raw, pos)
            offsets, pos = _vtk_array(raw, pos, n_offsets, offsets_hdr.split()[1])
            conn_hdr, pos = _vtk_line(raw, pos)
            if not conn_hdr.upper().startswith("CONNECTIVITY"):
                raise IndependentReadError(f"{p}: expected CONNECTIVITY, got {conn_hdr!r}")
            conn, pos = _vtk_array(raw, pos, n_conn, conn_hdr.split()[1])
        else:  # VTK < 5.1: connectivity then offsets, both in the CELLS type
            offsets, pos = _vtk_array(raw, pos, n_cells + 1, "vtktypeint32")
            conn, pos = _vtk_array(raw, pos, n_conn, "vtktypeint32")
    else:
        tokens = raw.decode("ascii").split()
        pos_idx = tokens.index("POINTS") + 1
        n_pts, points_type = int(tokens[pos_idx]), tokens[pos_idx + 1]
        if points_type not in ("double", "float"):
            raise IndependentReadError(f"{p}: unsupported ASCII POINTS type {points_type}")
        flat = [float(x) for x in tokens[pos_idx + 2:pos_idx + 2 + 3 * n_pts]]
        pts = np.asarray(flat, dtype=np.float64)
        c_idx = tokens.index("CELLS") + 1
        n_offsets, n_conn = int(tokens[c_idx]), int(tokens[c_idx + 1])
        n_cells = n_offsets - 1
        values = [int(x) for x in tokens[c_idx + 2:c_idx + 2 + n_offsets + n_conn]]
        offsets = values[:n_offsets]
        conn = np.asarray(values[n_cells + 1:], dtype=np.int64)

    cell_types = None
    types_at = raw.find(b"CELL_TYPES ", pos)
    if types_at >= 0:
        types_hdr, types_pos = _vtk_line(raw, types_at)
        n_types = int(types_hdr.split()[1])
        if encoding == "BINARY":
            cell_types, _ = _vtk_array(raw, types_pos, n_types, "vtktypeint32")
        else:
            tail = raw.decode("ascii")[types_pos - types_at:].split()
            cell_types = np.asarray([int(x) for x in tail[1:1 + n_types]], dtype=np.int32)

    cells = []
    for i in range(n_cells):
        lo, hi = int(offsets[i]), int(offsets[i + 1])
        cells.append(tuple(int(x) for x in conn[lo:hi]))
    if not cells:
        raise IndependentReadError(f"{p}: no cells")
    width = len(cells[0])
    if any(len(c) != width for c in cells):
        raise IndependentReadError(f"{p}: mixed cell widths")
    kind = {3: "tri", 4: "tetra"}.get(width)
    if kind is None:
        raise IndependentReadError(f"{p}: cells of width {width} are not tri/tetra")
    if cell_types is not None:
        expected = 5 if kind == "tri" else 10
        bad = sorted({int(t) for t in cell_types if int(t) != expected})
        if bad:
            raise IndependentReadError(
                f"{p}: CELL_TYPES {bad} disagree with {width}-node cells (expected {expected})"
            )
    return ParsedMesh(
        file_format="vtk",
        vertices=np.asarray(pts, dtype=np.float64).reshape(-1, 3),
        cells=_cells_array(cells),
        cell_kind=kind,
        declared_meta={
            "encoding": encoding,
            "cell_types": None if cell_types is None else sorted({int(t) for t in cell_types}),
            "points_type": points_type,
        },
    )


# ---------------------------------------------------------------------------
# Gmsh MSH 4.1 (binary)
# ---------------------------------------------------------------------------

_GMSH_NODES_PER_TYPE = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8, 10: 4}


def read_gmsh41(path: str | Path) -> ParsedMesh:
    p = Path(path)
    raw = p.read_bytes()
    if not raw.startswith(b"$MeshFormat"):
        raise IndependentReadError(f"{p}: no $MeshFormat section")
    version_line = raw.split(b"\n", 2)[1].decode("ascii")
    fields = version_line.split()
    version, file_type = fields[0], int(fields[1])
    size_t = np.int64 if (len(fields) > 2 and int(fields[2]) == 8) else np.int32
    if not version.startswith("4."):
        raise IndependentReadError(f"{p}: MSH {version} not parsed here (4.x only)")
    if file_type != 1:
        raise IndependentReadError(f"{p}: MSH 4.x ASCII not parsed here (binary only)")

    nodes_at = raw.index(b"$Nodes\n") + len(b"$Nodes\n")
    n_blocks = int(np.frombuffer(raw, dtype=size_t, count=1, offset=nodes_at)[0])
    pos = nodes_at + np.dtype(size_t).itemsize
    n_nodes, _min_tag, _max_tag = np.frombuffer(raw, dtype=size_t, count=3, offset=pos)
    pos += 3 * np.dtype(size_t).itemsize
    coords_by_tag: dict[int, tuple[float, float, float]] = {}
    for _ in range(int(n_blocks)):
        _dim, _tag, _param = np.frombuffer(raw, dtype=np.int32, count=3, offset=pos)
        pos += 3 * np.dtype(np.int32).itemsize
        count = int(np.frombuffer(raw, dtype=size_t, count=1, offset=pos)[0])
        pos += np.dtype(size_t).itemsize
        tags = np.frombuffer(raw, dtype=size_t, count=count, offset=pos)
        pos += count * np.dtype(size_t).itemsize
        xyz = np.frombuffer(raw, dtype=np.float64, count=3 * count, offset=pos).reshape(-1, 3)
        pos += 3 * count * np.dtype(np.float64).itemsize
        for tag, xyz_row in zip(tags, xyz):
            coords_by_tag[int(tag)] = (float(xyz_row[0]), float(xyz_row[1]), float(xyz_row[2]))

    elements_at = raw.index(b"$Elements\n") + len(b"$Elements\n")
    e_blocks = int(np.frombuffer(raw, dtype=size_t, count=1, offset=elements_at)[0])
    pos = elements_at + np.dtype(size_t).itemsize
    pos += 3 * np.dtype(size_t).itemsize  # n_elements, min_tag, max_tag
    cells: list[tuple[int, ...]] = []
    element_types: set[int] = set()
    element_tags: set[int] = set()
    for _ in range(int(e_blocks)):
        _dim, _tag, etype = (int(v) for v in np.frombuffer(raw, dtype=np.int32, count=3, offset=pos))
        pos += 3 * np.dtype(np.int32).itemsize
        count = int(np.frombuffer(raw, dtype=size_t, count=1, offset=pos)[0])
        pos += np.dtype(size_t).itemsize
        per = _GMSH_NODES_PER_TYPE.get(etype)
        if per is None:
            raise IndependentReadError(f"{p}: gmsh element type {etype} not handled")
        # MSH 4.1 element record: elementTag(size_t) then nodeTag(size_t) x per
        conn = np.frombuffer(
            raw, dtype=size_t, count=count * (1 + per), offset=pos
        ).reshape(count, 1 + per)
        pos += count * (1 + per) * np.dtype(size_t).itemsize
        element_tags.update(int(v) for v in conn[:, 0])
        element_types.add(etype)
        for row in conn:
            cells.append(tuple(int(x) for x in row[1:]))

    if not cells:
        raise IndependentReadError(f"{p}: no elements")
    width = len(cells[0])
    kind = {3: "tri", 4: "tetra"}.get(width)
    if kind is None:
        raise IndependentReadError(f"{p}: {width}-node elements are not tri/tetra")
    tags = sorted(coords_by_tag)
    index = {tag: i for i, tag in enumerate(tags)}
    vertices = np.asarray([coords_by_tag[t] for t in tags], dtype=np.float64)
    try:
        cells = [tuple(index[t] for t in c) for c in cells]
    except KeyError as exc:
        raise IndependentReadError(f"{p}: element references unknown node tag {exc}") from exc
    return ParsedMesh(
        file_format="gmsh",
        vertices=vertices,
        cells=_cells_array(cells),
        cell_kind=kind,
        declared_meta={
            "version": version,
            "element_types": sorted(element_types),
            "element_tag_range": [min(element_tags), max(element_tags)] if element_tags else None,
            "node_count_declared": int(n_nodes),
        },
    )


# ---------------------------------------------------------------------------
# SU2 (ASCII)
# ---------------------------------------------------------------------------

_SU2_NODES = {5: 3, 9: 4, 10: 4, 12: 8, 3: 2}


def read_su2(path: str | Path) -> ParsedMesh:
    p = Path(path)
    lines = [ln.strip() for ln in p.read_text().splitlines()]
    n_dime = n_poin = None
    vertices: list[tuple[float, ...]] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("NDIME="):
            n_dime = int(line.split("=")[1])
            i += 1
        elif line.startswith("NPOIN="):
            n_poin = int(line.split("=")[1])
            i += 1
            for _ in range(n_poin):
                values = [float(x) for x in lines[i].split()]
                if len(values) < n_dime:
                    raise IndependentReadError(f"{p}: short point record {lines[i]!r}")
                vertices.append(tuple(values[:n_dime]))
                i += 1
        elif line.startswith("NELEM="):
            n_elem = int(line.split("=")[1])
            i += 1
            cells: list[tuple[int, ...]] = []
            element_types: set[int] = set()
            for _ in range(n_elem):
                parts = [int(x) for x in lines[i].split()]
                etype = parts[0]
                per = _SU2_NODES.get(etype)
                if per is None:
                    raise IndependentReadError(f"{p}: SU2 element type {etype} not handled")
                if len(parts) < 1 + per:
                    raise IndependentReadError(f"{p}: short element record {lines[i]!r}")
                element_types.add(etype)
                cells.append(tuple(parts[1:1 + per]))
                i += 1
            return _su2_result(p, vertices, cells, element_types, n_dime)
        else:
            i += 1
    raise IndependentReadError(f"{p}: no NELEM section")


def _su2_result(
    p: Path,
    vertices: list[tuple[float, ...]],
    cells: list[tuple[int, ...]],
    element_types: set[int],
    n_dime: int | None,
) -> ParsedMesh:
    if not cells:
        raise IndependentReadError(f"{p}: no elements")
    width = len(cells[0])
    kind = {3: "tri", 4: "tetra"}.get(width)
    if kind is None:
        raise IndependentReadError(f"{p}: {width}-node elements are not tri/tetra")
    v = np.zeros((len(vertices), 3), dtype=np.float64)
    for i, row in enumerate(vertices):
        v[i, :len(row)] = row
    return ParsedMesh(
        file_format="su2",
        vertices=v,
        cells=_cells_array(cells),
        cell_kind=kind,
        declared_meta={"ndime": n_dime, "element_types": sorted(element_types)},
    )


# ---------------------------------------------------------------------------
# CGNS (HDF5, read straight off the SIDS node paths)
# ---------------------------------------------------------------------------

def read_cgns_hdf5(path: str | Path) -> ParsedMesh:
    p = Path(path)
    raw = p.read_bytes()
    if raw[:8] != b"\x89HDF\r\n\x1a\n":
        raise IndependentReadError(f"{p}: not an HDF5/CGNS file (bad magic)")
    if b"CGNSLibraryVersion" not in raw:
        raise IndependentReadError(f"{p}: no CGNSLibraryVersion node; not a CGNS file")
    try:
        import h5py
    except ImportError as exc:  # pragma: no cover - h5py is a test dependency
        raise IndependentReadError("h5py is required to read CGNS independently") from exc

    with h5py.File(p, "r") as f:
        bases = [k for k in f if not k.startswith("CGNS")]
        if not bases:
            raise IndependentReadError(f"{p}: no CGNS base node")
        base = f[bases[0]]
        zones = [k for k in base if k.startswith("Zone")]
        if not zones:
            raise IndependentReadError(f"{p}: base {bases[0]!r} has no zone")
        zone = base[zones[0]]
        grid = zone.get("GridCoordinates")
        if grid is None:
            raise IndependentReadError(f"{p}: zone {zones[0]!r} has no GridCoordinates")

        def coord(axis: str) -> np.ndarray:
            node = grid.get(axis)
            if node is None:
                raise IndependentReadError(f"{p}: missing {axis}")
            return np.asarray(_sids_data(node), dtype=np.float64)

        x, y, z = coord("CoordinateX"), coord("CoordinateY"), coord("CoordinateZ")
        if not (len(x) == len(y) == len(z)):
            raise IndependentReadError(f"{p}: coordinate arrays differ in length")
        vertices = np.stack((x, y, z), axis=1)

        elems = zone.get("GridElements")
        if elems is None:
            raise IndependentReadError(f"{p}: zone {zones[0]!r} has no GridElements")
        conn_node = elems.get("ElementConnectivity")
        conn = np.asarray(
            _sids_data(conn_node), dtype=np.int64,
        ).ravel()
        rng_node = elems.get("ElementRange")
        rng = np.asarray(
            _sids_data(rng_node), dtype=np.int64,
        ).ravel()
        n_elems = int(rng[1] - rng[0] + 1) if rng.size >= 2 else 0
        if n_elems == 0:
            raise IndependentReadError(f"{p}: empty ElementRange")
        if conn.size % n_elems != 0:
            raise IndependentReadError(
                f"{p}: connectivity {conn.size} not divisible by {n_elems} elements"
            )
        per = conn.size // n_elems
        if per not in (3, 4):
            raise IndependentReadError(f"{p}: {per} nodes per element is not tri/tetra")
        cells = tuple(tuple(int(i) - 1 for i in row) for row in conn.reshape(n_elems, per))
        cell_dim = int(np.asarray(_sids_data(base["CellDimension"]), dtype=np.int64).ravel()[0])
        phys_dim = int(np.asarray(_sids_data(base["PhysicalDimension"]), dtype=np.int64).ravel()[0])
        lib_version = float(np.asarray(f["CGNSLibraryVersion"][...]).ravel()[0])

    kind = "tri" if per == 3 else "tetra"
    if (kind == "tetra") != (cell_dim == 3):
        raise IndependentReadError(
            f"{p}: CellDimension {cell_dim} disagrees with {per}-node elements"
        )
    if any(int(i) < 0 or int(i) >= len(vertices) for row in cells for i in row):
        raise IndependentReadError(f"{p}: connectivity indexes outside the point array")
    return ParsedMesh(
        file_format="cgns",
        vertices=vertices,
        cells=_cells_array(cells),
        cell_kind=kind,
        declared_meta={
            "zone": zones[0],
            "cell_dimension": cell_dim,
            "physical_dimension": phys_dim,
            "cgns_library_version": lib_version,
        },
    )


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------

_READERS = {
    "stl": read_stl_ascii,
    "obj": read_obj,
    "vtk": read_vtk_legacy,
    "gmsh": read_gmsh41,
    "su2": read_su2,
    "cgns": read_cgns_hdf5,
}


def read_independent(path: str | Path, *, file_format: str) -> ParsedMesh:
    """Parse ``path`` with the independent reader for ``file_format``."""
    try:
        reader = _READERS[file_format]
    except KeyError as exc:
        raise IndependentReadError(
            f"no independent reader for '{file_format}'; supported: {sorted(_READERS)}"
        ) from exc
    return reader(path)

# ---------------------------------------------------------------------------
# The 6-format x 3-geometry campaign (single source of truth for the matrix)
# ---------------------------------------------------------------------------

def _mesh_builder_cube():
    """Unit cube, outward oriented, 12 triangles, 8 vertices, V = 1."""
    from .geometry import MeshBuilder

    b = MeshBuilder(units="m")
    points = (
        (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
    )
    idx = [b.add_vertex(p) for p in points]
    faces = (
        (0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4), (2, 3, 7), (2, 7, 6),
        (1, 2, 6), (1, 6, 5), (0, 4, 7), (0, 7, 3),
    )
    for f in faces:
        b.add_triangle(idx[f[0]], idx[f[1]], idx[f[2]])
    return b.build()


def _mesh_capsule(n_theta: int = 32, n_axial: int = 12):
    """Closed body of revolution (cylinder + two hemispheres), V = pi r^2 L + 4/3 pi r^3.

    Same fixture the package matrix uses: a closed, star-shaped, *curved* body,
    i.e. the interesting case for faceting/vertex-order losses.
    """
    import math

    from . import mesh as meshlib
    from .geometry import MeshBuilder

    radius, length = 0.5, 4.0
    b = MeshBuilder(units="m")
    sections: list[tuple[float, float]] = []
    for i in range(n_axial + 1):
        phi = (math.pi / 2) * i / n_axial
        sections.append((-length / 2 - radius * math.cos(phi), radius * math.sin(phi)))
    sections.append((length / 2, radius))
    for i in range(n_axial - 1, -1, -1):
        phi = (math.pi / 2) * i / n_axial
        sections.append((length / 2 + radius * math.cos(phi), radius * math.sin(phi)))

    rings: list[list[int]] = []
    for x, r in sections:
        if r == 0.0:
            rings.append([b.add_vertex((x, 0.0, 0.0))])
            continue
        rings.append([
            b.add_vertex((x, r * math.cos(2 * math.pi * k / n_theta),
                          r * math.sin(2 * math.pi * k / n_theta)))
            for k in range(n_theta)
        ])

    for r0, r1 in zip(rings, rings[1:]):
        if len(r0) == 1:
            for k in range(n_theta):
                b.add_triangle(r0[0], r1[k], r1[(k + 1) % n_theta])
        elif len(r1) == 1:
            for k in range(n_theta):
                b.add_triangle(r1[0], r0[(k + 1) % n_theta], r0[k])
        else:
            for k in range(n_theta):
                j = (k + 1) % n_theta
                b.add_quad(r0[k], r0[j], r1[j], r1[k])
    return meshlib.orient_outward(b.build()).mesh


def _mesh_wigley_hull():
    """The package's own validated Wigley hull (real ship geometry, 1001 pts)."""
    from . import wigley

    return wigley.generate_wigley().mesh


#: Analytic volumes of the three campaign fixtures.
ANALYTIC_VOLUME_3: dict[str, float] = {
    "cube": 1.0,
    "capsule": math.pi * 0.5 ** 2 * 4.0 + 4.0 / 3.0 * math.pi * 0.5 ** 3,
    "wigley_hull": 4.0 * 10.0 * 1.0 * 0.625 / 9.0,
}

_FIXTURES = {
    "cube": _mesh_builder_cube,
    "capsule": _mesh_capsule,
    "wigley_hull": _mesh_wigley_hull,
}

_MESH_CACHE: dict[str, Any] = {}


def campaign_mesh(name: str):
    """Surface mesh of one campaign fixture (cached)."""
    if name not in _MESH_CACHE:
        _MESH_CACHE[name] = _FIXTURES[name]()
    return _MESH_CACHE[name]


def stl_normal_consistency(parsed: ParsedMesh) -> dict[str, float]:
    """STL carries a redundant per-facet normal; does it agree with the winding?

    Independent readers may trust either the stored normal or the vertex order.
    A writer that emits the two inconsistently produces a file that renders one
    way and meshes the other, so this is checked rather than assumed.
    """
    v = np.asarray(parsed.vertices, dtype=np.float64)
    stated = np.asarray(parsed.declared_meta["normals"], dtype=np.float64)
    computed = np.array([
        np.cross(v[c[1]] - v[c[0]], v[c[2]] - v[c[0]])
        for c in parsed.cells
    ])
    norms = np.linalg.norm(computed, axis=1)
    norms[norms == 0.0] = 1.0
    unit = computed / norms[:, None]
    cos = np.einsum("ij,ij->i", unit, stated)
    return {
        "min_cos": float(cos.min()),
        "max_abs_deviation": float(np.abs(cos - 1.0).max()),
        "facet_count": int(len(parsed.cells)),
    }


def run_format_matrix(
    *,
    geometries: Sequence[str] = ("cube", "capsule", "wigley_hull"),
    formats: Sequence[str] | None = None,
    workdir: str | Path | None = None,
    volume_tolerance: float = 1e-9,
) -> list[dict[str, Any]]:
    """Write every geometry in every format, read it back with the independent
    parser and record volume / cell-count / per-cell winding evidence.

    Returns one row per (geometry, format) pair; ``row["ok"]`` is the conjunction
    of every invariant.  The table is the artifact the ledger acceptance asks
    for, computed from the file bytes rather than from the writer's own reader.
    """
    import tempfile

    from .meshio_exchange import FORMATS_ALL, tetrahedralise_surface, write_meshio

    formats = tuple(FORMATS_ALL if formats is None else formats)
    base = Path(workdir) if workdir is not None else Path(tempfile.mkdtemp(prefix="fmtmatrix-"))
    base.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    for geom_name in geometries:
        surface = campaign_mesh(geom_name)
        for fmt in formats:
            volume_format = fmt in ("su2", "cgns")
            mesh = tetrahedralise_surface(surface) if volume_format else surface
            kind = "tetra" if volume_format else "tri"
            cells_src = mesh.tets if volume_format else mesh.faces
            path = write_meshio(base / f"{geom_name}.{fmt}", mesh, file_format=fmt)
            parsed = read_independent(path, file_format=fmt)
            vol_ref = (
                tetrahedral_volume(np.asarray(mesh.vertices, dtype=np.float64), cells_src)
                if volume_format
                else surface_signed_volume(np.asarray(mesh.vertices, dtype=np.float64), cells_src)
            )
            vol_ind = parsed.volume()
            rel = abs(vol_ind - vol_ref) / max(abs(vol_ref), 1e-300)
            sig_ref = oriented_signature(
                np.asarray(mesh.vertices, dtype=np.float64), cells_src, cell_kind=kind
            )
            sig_ind = oriented_signature(parsed.vertices, parsed.cells, cell_kind=parsed.cell_kind)
            identity = file_identity(path, file_format=fmt)
            row: dict[str, Any] = {
                "geometry": geom_name,
                "format": fmt,
                "cell_kind": parsed.cell_kind,
                "cells": len(parsed.cells),
                "cells_reference": len(cells_src),
                "vertices": int(parsed.vertices.shape[0]),
                "vertices_reference": len(mesh.vertices),
                "volume_reference": vol_ref,
                "volume_independent": vol_ind,
                "volume_relative_error": rel,
                "analytic_volume": ANALYTIC_VOLUME_3.get(geom_name),
                "cell_count_preserved": len(parsed.cells) == len(cells_src),
                "vertex_count_preserved": int(parsed.vertices.shape[0]) == len(mesh.vertices),
                "winding_preserved": sig_ref == sig_ind,
                "identity": identity,
                "declared_meta": parsed.declared_meta,
                "bytes": int(Path(path).stat().st_size),
                "volume_tolerance": volume_tolerance,
            }
            if fmt == "stl":
                row["stl_normals"] = stl_normal_consistency(parsed)
            row["ok"] = bool(
                row["cell_count_preserved"]
                and row["winding_preserved"]
                and rel <= volume_tolerance
            )
            rows.append(row)
    return rows


def matrix_markdown(rows: Sequence[dict[str, Any]]) -> str:
    """Render :func:`run_format_matrix` output as the acceptance table."""
    lines = [
        "| geometry | format | cells | extracted V [m^3] | reference V [m^3] | rel. err | V analytic | winding | count | ok |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            "| {geometry} | {format} | {cells}/{cells_reference} | {volume_independent:.12f} | "
            "{volume_reference:.12f} | {volume_relative_error:.2e} | {analytic:.9f} | "
            "{winding} | {count} | {ok} |".format(
                geometry=r["geometry"],
                format=r["format"],
                cells=r["cells"],
                cells_reference=r["cells_reference"],
                volume_independent=r["volume_independent"],
                volume_reference=r["volume_reference"],
                volume_relative_error=r["volume_relative_error"],
                analytic=r["analytic_volume"] or float("nan"),
                winding="PASS" if r["winding_preserved"] else "FAIL",
                count="PASS" if r["cell_count_preserved"] else "FAIL",
                ok="PASS" if r["ok"] else "FAIL",
            )
        )
    return "\n".join(lines) + "\n"
