"""Native six-format surface-mesh exchange for TensorCAD (the product writer set).

One mesh, six on-disk grammars, one comparison contract::

    MESH_FORMATS = ("stl", "obj", "vtk", "su2", "cgns", "gmsh")

    m = icosphere(radius=1.0, subdivisions=2)
    write(m, "cgns", "out/sphere.cgns")     # -> Path
    back = read("out/sphere.cgns", "cgns")  # -> TriMesh
    report = compare_meshes(m, back)        # -> dict, `pass` is the machine gate

Design decisions that are not arbitrary, each bought with a measured failure:

* **VTK defaults to ``UNSTRUCTURED_GRID`` + ``CELLS`` + ``CELL_TYPES=5``.**
  Legacy ``POLYDATA`` is legal VTK and ParaView reads it, but ``meshio``'s legacy
  reader only accepts ``UNSTRUCTURED_GRID``/``STRUCTURED_*``
  (``meshio/vtk/_vtk_42.py``), so a ``POLYDATA`` file is invisible to the Python
  mesh ecosystem.  ``POLYDATA`` stays available behind ``dataset=``.
* **SU2 defaults to ``nelem_mode="volume"``** (triangles inside ``NELEM``): legal
  SU2 and readable by ``meshio``.  ``nelem_mode="boundary"`` emits the faithful
  boundary-only subset (``NELEM=0`` + one ``MARKER_ELEMS`` block), also legal, but
  ``meshio/su2/_su2.py:113`` raises ``StopIteration`` on it -- a reader defect,
  documented rather than worked around.
* **CGNS is written in the SIDS-to-HDF5 dialect and is checked by
  :func:`sids_violations`**, including the root ``CGNSLibraryVersion_t`` node.
  That node is the reason this module exists: the earlier reference dialect had
  name/label/type/order but no version node, so calling it SIDS-conformant was an
  overstatement.
* **Precision is float64 everywhere except STL**, which the format itself fixes to
  float32.  STL is also a facet soup (no topology): weld before asking topological
  questions, and keep models in a local frame -- float32 carries ~7 significant
  digits, so a 1000 m coordinate offset already costs ~4e-5 m.

Nothing here is ship-specific: it is geometry-exchange plumbing.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

__all__ = [
    "MESH_FORMATS",
    "TriMesh",
    "weld",
    "cube",
    "icosphere",
    "torus",
    "write",
    "read",
    "sids_violations",
    "compare_meshes",
    "CGNS_RULES",
    "WRITERS",
    "READERS",
]

#: Canonical format names, in report order.
MESH_FORMATS: Tuple[str, ...] = ("stl", "obj", "vtk", "su2", "cgns", "gmsh")

_ALIASES: Dict[str, str] = {
    "msh": "gmsh",
    "msh2": "gmsh",
    "stlb": "stl",
    "wavefront": "obj",
}

#: SU2 v7 element type id for a 3-node triangle.
SU2_TRIANGLE = 5
#: Gmsh MSH 2.2 element type id for a 3-node triangle.
GMSH_TRIANGLE = 2
#: CGNS library version written into ``CGNSLibraryVersion_t`` (R4).
CGNS_LIBRARY_VERSION = 4.3


def canonical_format(name: str) -> str:
    """Map a user spelling onto :data:`MESH_FORMATS` or raise ``ValueError``."""
    key = str(name).strip().lower().lstrip(".")
    key = _ALIASES.get(key, key)
    if key not in MESH_FORMATS:
        raise ValueError(f"unknown mesh format {name!r}; expected one of {MESH_FORMATS}")
    return key


# ---------------------------------------------------------------------------
# mesh container
# ---------------------------------------------------------------------------
@dataclass
class TriMesh:
    """Triangle mesh with explicit connectivity (float64 points, int64 faces)."""

    points: np.ndarray
    faces: np.ndarray
    name: str = "mesh"

    def __post_init__(self) -> None:
        self.points = np.ascontiguousarray(np.asarray(self.points, dtype=float))
        self.faces = np.ascontiguousarray(np.asarray(self.faces, dtype=np.int64))
        if self.points.ndim != 2 or self.points.shape[1] != 3:
            raise ValueError(f"points must be (n,3), got {self.points.shape}")
        if self.faces.ndim != 2 or self.faces.shape[1] != 3:
            raise ValueError(f"faces must be (m,3), got {self.faces.shape}")
        if self.faces.size and (self.faces.min() < 0 or self.faces.max() >= len(self.points)):
            raise ValueError("faces index out of range")

    def volume(self) -> float:
        """Enclosed volume via the divergence theorem (0 for an open surface)."""
        if not len(self.faces):
            return 0.0
        tri = self.points[self.faces]
        return float(np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0)

    def area(self) -> float:
        tri = self.points[self.faces]
        return float(0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1).sum())

    def euler(self) -> int:
        edges = np.sort(
            np.vstack([self.faces[:, [0, 1]], self.faces[:, [1, 2]], self.faces[:, [2, 0]]]), axis=1
        )
        uniq, counts = np.unique(edges, axis=0, return_counts=True)
        used = int(np.unique(self.faces).size)
        return int(used - uniq.shape[0] + len(self.faces))

    def boundary_edges(self) -> int:
        """Edges not shared by exactly two triangles (0 for a closed surface)."""
        edges = np.sort(
            np.vstack([self.faces[:, [0, 1]], self.faces[:, [1, 2]], self.faces[:, [2, 0]]]), axis=1
        )
        _uniq, counts = np.unique(edges, axis=0, return_counts=True)
        return int((counts != 2).sum())


def weld(mesh: TriMesh, decimals: int = 12) -> Tuple[TriMesh, int]:
    """Recover connectivity from a corner soup by rounding-invariant welding."""
    key = np.round(mesh.points, decimals)
    uniq, inverse = np.unique(key, axis=0, return_inverse=True)
    faces = np.asarray(inverse).reshape(-1)[mesh.faces]
    good = (faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2]) & (faces[:, 0] != faces[:, 2])
    return TriMesh(uniq, faces[good], mesh.name + "/welded"), int(len(mesh.points) - len(uniq))


# ---------------------------------------------------------------------------
# analytic bodies (closed; genus 0 and genus 1) with exact volumes
# ---------------------------------------------------------------------------
def cube(edge: float = 2.0) -> TriMesh:
    """Axis-aligned cube; volume is exactly ``edge**3`` (analytic oracle)."""
    h = edge / 2.0
    corners = np.array([[x, y, z] for x in (-h, h) for y in (-h, h) for z in (-h, h)], dtype=float)
    quads = [
        (0, 1, 3, 2),
        (4, 6, 7, 5),
        (0, 4, 5, 1),
        (2, 3, 7, 6),
        (0, 2, 6, 4),
        (1, 5, 7, 3),
    ]
    tris = []
    for a, b, c, d in quads:
        tris.append((a, b, c))
        tris.append((a, c, d))
    return TriMesh(corners, np.array(tris, dtype=np.int64), "cube")


def icosphere(radius: float = 1.0, subdivisions: int = 2) -> TriMesh:
    """Icosphere of exact radius; volume converges to ``4/3*pi*r**3``."""
    t = (1.0 + 5.0 ** 0.5) / 2.0
    verts = np.array(
        [
            [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
            [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
            [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1],
        ],
        dtype=float,
    )
    faces = np.array(
        [
            [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
            [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
            [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
            [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
        ],
        dtype=np.int64,
    )
    for _ in range(int(subdivisions)):
        cache: Dict[Tuple[int, int], int] = {}
        verts = list(verts.tolist())

        def mid(u: int, v: int) -> int:
            key = (u, v) if u < v else (v, u)
            i = cache.get(key)
            if i is None:
                a, b = verts[key[0]], verts[key[1]]
                verts.append([(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0, (a[2] + b[2]) / 2.0])
                i = len(verts) - 1
                cache[key] = i
            return i

        new_faces = []
        for a, b, c in faces.tolist():
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            new_faces += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        verts = np.array(verts, dtype=float)
        faces = np.array(new_faces, dtype=np.int64)
    verts = radius * verts / np.linalg.norm(verts, axis=1, keepdims=True)
    return TriMesh(verts, faces, f"icosphere{subdivisions}")


def torus(major: float = 1.0, minor: float = 0.3, nu: int = 64, nv: int = 32) -> TriMesh:
    """Torus (genus 1); volume is exactly ``2*pi**2*major*minor**2``."""
    u = np.linspace(0.0, 2 * np.pi, nu, endpoint=False)
    v = np.linspace(0.0, 2 * np.pi, nv, endpoint=False)
    uu, vv = np.meshgrid(u, v, indexing="ij")
    pts = np.stack(
        [
            (major + minor * np.cos(vv)) * np.cos(uu),
            (major + minor * np.cos(vv)) * np.sin(uu),
            minor * np.sin(vv),
        ],
        axis=-1,
    ).reshape(-1, 3)
    idx = np.arange(nu * nv).reshape(nu, nv)
    faces = []
    for i in range(nu):
        for j in range(nv):
            a, b = idx[i, j], idx[(i + 1) % nu, j]
            c, d = idx[(i + 1) % nu, (j + 1) % nv], idx[i, (j + 1) % nv]
            faces += [(a, b, c), (a, c, d)]
    return TriMesh(pts, np.array(faces, dtype=np.int64), f"torus{nu}x{nv}")

# ---------------------------------------------------------------------------
# writers
# ---------------------------------------------------------------------------
def write_stl(mesh: TriMesh, path, header: str = "TensorCAD native binary STL") -> Path:
    """Binary STL (float32 by format definition): 80-byte header, uint32 count, facets."""
    tri = mesh.points[mesh.faces]
    nrm = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(nrm, axis=1, keepdims=True)
    nrm = np.divide(nrm, ln, out=np.zeros_like(nrm), where=ln > 0)
    with open(path, "wb") as fh:
        fh.write(header.encode("ascii", "replace")[:80].ljust(80, b" "))
        fh.write(struct.pack("<I", len(mesh.faces)))
        for i in range(len(mesh.faces)):
            fh.write(struct.pack("<3f", *nrm[i].astype(np.float32)))
            fh.write(struct.pack("<9f", *tri[i].astype(np.float32).ravel()))
            fh.write(struct.pack("<H", 0))
    return Path(path)


def write_obj(mesh: TriMesh, path) -> Path:
    """Wavefront OBJ, minimal subset (``v``/``f``), 17 significant digits."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# TensorCAD native OBJ {mesh.name}\n")
        for x, y, z in mesh.points:
            fh.write(f"v {x:.17g} {y:.17g} {z:.17g}\n")
        for a, b, c in mesh.faces + 1:
            fh.write(f"f {a} {b} {c}\n")
    return Path(path)


def write_vtk(mesh: TriMesh, path, dataset: str = "UNSTRUCTURED_GRID") -> Path:
    """VTK legacy ASCII 3.0; ``UNSTRUCTURED_GRID`` by default (ecosystem-readable)."""
    dataset = dataset.upper()
    if dataset not in ("UNSTRUCTURED_GRID", "POLYDATA"):
        raise ValueError(dataset)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# vtk DataFile Version 3.0\n")
        fh.write(f"TensorCAD native {mesh.name} {dataset}\n")
        fh.write("ASCII\n")
        fh.write(f"DATASET {dataset}\n")
        fh.write(f"POINTS {len(mesh.points)} double\n")
        for x, y, z in mesh.points:
            fh.write(f"{x:.17g} {y:.17g} {z:.17g}\n")
        if dataset == "POLYDATA":
            fh.write(f"POLYGONS {len(mesh.faces)} {4 * len(mesh.faces)}\n")
        else:
            fh.write(f"CELLS {len(mesh.faces)} {4 * len(mesh.faces)}\n")
        for a, b, c in mesh.faces:
            fh.write(f"3 {a} {b} {c}\n")
        if dataset == "UNSTRUCTURED_GRID":
            fh.write(f"CELL_TYPES {len(mesh.faces)}\n")
            fh.write("".join("5\n" for _ in mesh.faces))
    return Path(path)


def write_su2(mesh: TriMesh, path, marker: str = "WALL", nelem_mode: str = "volume") -> Path:
    """SU2 ASCII mesh; ``nelem_mode='volume'`` (default) keeps triangles in ``NELEM``."""
    if nelem_mode not in ("volume", "boundary"):
        raise ValueError(nelem_mode)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"% TensorCAD native SU2 {mesh.name}\n")
        fh.write("NDIME= 3\n")
        fh.write(f"NPOIN= {len(mesh.points)}\n")
        for x, y, z in mesh.points:
            fh.write(f"{x:.17g} {y:.17g} {z:.17g}\n")
        if nelem_mode == "volume":
            fh.write(f"NELEM= {len(mesh.faces)}\n")
            for a, b, c in mesh.faces:
                fh.write(f"{SU2_TRIANGLE} {a} {b} {c} 0\n")
            fh.write("NMARK= 0\n")
        else:
            fh.write("NELEM= 0\n")
            fh.write("NMARK= 1\n")
            fh.write(f"MARKER_TAG= {marker}\n")
            fh.write(f"MARKER_ELEMS= {len(mesh.faces)}\n")
            for a, b, c in mesh.faces:
                fh.write(f"{SU2_TRIANGLE} {a} {b} {c}\n")
    return Path(path)


def _h5py():
    try:
        import h5py
    except ImportError as exc:  # pragma: no cover - environment guard
        raise RuntimeError("CGNS support needs h5py (`pip install h5py`)") from exc
    return h5py


def _sids_attrs(node, props: Dict[str, object]) -> None:
    """SIDS-to-HDF5 node attributes, written with the spec's leading space.

    A CGNS node in the HDF5 mapping is a group carrying ``name``/``label``/
    ``type``/``order``; the leading space keeps them apart from HDF5's own
    attributes (SIDS-to-HDF5, "Basic Node").  :func:`sids_violations` accepts the
    space-prefixed and the bare spelling, since both occur in the wild.
    """
    for key, value in props.items():
        node.attrs[" " + key] = np.bytes_(value) if isinstance(value, str) else np.int32(value)


def _sids_node(parent, name: str, label: str, dtype: str = "MT", order: int = 1, data=None):
    grp = parent.create_group(name)
    _sids_attrs(grp, {"name": name, "label": label, "type": dtype, "order": order})
    if data is not None:
        ds = grp.create_dataset(" data", data=data)
        ds.attrs[" type"] = np.bytes_(dtype)
    return grp


def write_cgns(mesh: TriMesh, path, library_version: float = CGNS_LIBRARY_VERSION) -> Path:
    """CGNS/HDF5, SIDS-to-HDF5 dialect, single unstructured zone, ``TRI_3``.

    Node tree::

        /                    reserved root attributes version/format
        CGNSLibraryVersion   label CGNSLibraryVersion_t, R4 = 4.3
        Base (CGNSBase_t)    I4 [CellDim, PhysDim] = [3, 3]
          Zone (Zone_t)      I4 [nvertex, ncell, 0]
            ZoneType         C1 "Unstructured"
            GridCoordinates  CoordinateX/Y/Z (DataArray_t, R8)
            Elements         I4 [TRI_3, 0]
              ElementRange       I4 [1, ncell]
              ElementConnectivity I4 1-based, ncell x 3
    """
    h5py = _h5py()
    nv, nf = len(mesh.points), len(mesh.faces)
    with h5py.File(path, "w") as f:
        _sids_attrs(f, {"version": "3.4", "format": "0.9"})  # reserved root attributes
        _sids_node(
            f, "CGNSLibraryVersion", "CGNSLibraryVersion_t", "R4", 1,
            np.array([library_version], dtype=np.float32),
        )
        base = _sids_node(f, "Base", "CGNSBase_t", "I4", 1, np.array([3, 3], dtype=np.int32))
        zone = _sids_node(base, "Zone", "Zone_t", "I4", 1, np.array([nv, nf, 0], dtype=np.int32))
        zt = _sids_node(zone, "ZoneType", "ZoneType_t", "C1", 1)
        ds = zt.create_dataset(" data", data=np.bytes_("Unstructured"))
        ds.attrs[" type"] = np.bytes_("C1")
        gc = _sids_node(zone, "GridCoordinates", "GridCoordinates_t", "MT", 1)
        for k, axis in enumerate("XYZ"):
            _sids_node(
                gc, f"Coordinate{axis}", "DataArray_t", "R8", k + 1,
                np.ascontiguousarray(mesh.points[:, k]),
            )
        el = _sids_node(zone, "Elements", "Elements_t", "I4", 2,
                        np.array([SU2_TRIANGLE, 0], dtype=np.int32))
        _sids_node(el, "ElementRange", "IndexRange_t", "I4", 1, np.array([1, nf], dtype=np.int32))
        _sids_node(el, "ElementConnectivity", "DataArray_t", "I4", 2,
                   (mesh.faces + 1).astype(np.int32).ravel())
    return Path(path)


def write_gmsh(mesh: TriMesh, path) -> Path:
    """Gmsh MSH 2.2 ASCII (widest reader support): element type 2, 1-based, 0 tags."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        fh.write(f"$Nodes\n{len(mesh.points)}\n")
        for i, (x, y, z) in enumerate(mesh.points, 1):
            fh.write(f"{i} {x:.17g} {y:.17g} {z:.17g}\n")
        fh.write("$EndNodes\n")
        fh.write(f"$Elements\n{len(mesh.faces)}\n")
        for i, (a, b, c) in enumerate(mesh.faces + 1, 1):
            fh.write(f"{i} {GMSH_TRIANGLE} 0 {a} {b} {c}\n")
        fh.write("$EndElements\n")
    return Path(path)


WRITERS = {
    "stl": write_stl,
    "obj": write_obj,
    "vtk": write_vtk,
    "su2": write_su2,
    "cgns": write_cgns,
    "gmsh": write_gmsh,
}


def write(mesh: TriMesh, fmt: str, path, **kwargs) -> Path:
    """Write ``mesh`` as ``fmt`` to ``path`` (parent directories are created)."""
    key = canonical_format(fmt)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return WRITERS[key](mesh, path, **kwargs)


# ---------------------------------------------------------------------------
# readers: our own implementations, so a round trip exercises the grammar
# rather than a third-party library
# ---------------------------------------------------------------------------
_STL_REC = np.dtype([("n", "<3f4"), ("v", "<9f4"), ("attr", "<u2")])


def read_stl(path) -> TriMesh:
    data = Path(path).read_bytes()
    n = struct.unpack("<I", data[80:84])[0]
    if len(data) < 84 + 50 * n:
        raise ValueError(f"STL: truncated ({len(data)} bytes for {n} facets)")
    rec = np.frombuffer(data[84 : 84 + 50 * n], dtype=_STL_REC)
    pts = rec["v"].reshape(-1, 3).astype(float)
    faces = np.arange(len(pts), dtype=np.int64).reshape(-1, 3)
    return TriMesh(pts, faces, "stl")


def read_obj(path) -> TriMesh:
    pts, faces = [], []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if not parts or parts[0].startswith("#"):
            continue
        if parts[0] == "v":
            pts.append([float(v) for v in parts[1:4]])
        elif parts[0] == "f":
            idx = [int(tok.split("/")[0]) for tok in parts[1:]]
            for k in range(1, len(idx) - 1):  # fan-triangulate polygons
                faces.append([idx[0] - 1, idx[k] - 1, idx[k + 1] - 1])
    if not pts:
        raise ValueError("OBJ: no 'v' records")
    return TriMesh(np.array(pts, dtype=float), np.array(faces, dtype=np.int64), "obj")


def read_vtk(path) -> TriMesh:
    tokens = Path(path).read_text(encoding="utf-8").split()
    if len(tokens) < 6 or tokens[0] != "#":
        raise ValueError("VTK: not a legacy ASCII file")
    dataset = tokens[tokens.index("DATASET") + 1]
    i = tokens.index("POINTS")
    nv = int(tokens[i + 1])
    body = np.array(tokens[i + 3 : i + 3 + 3 * nv], dtype=float).reshape(nv, 3)
    if dataset == "POLYDATA":
        i = tokens.index("POLYGONS")
    elif dataset == "UNSTRUCTURED_GRID":
        i = tokens.index("CELLS")
    else:  # pragma: no cover - only the two dialects are emitted
        raise ValueError(f"VTK: unsupported DATASET {dataset}")
    ncells = int(tokens[i + 1])
    faces: List[List[int]] = []
    j = i + 3
    for _ in range(ncells):
        k = int(tokens[j])
        corners = [int(t) for t in tokens[j + 1 : j + 1 + k]]
        for t in range(1, k - 1):
            faces.append([corners[0], corners[t], corners[t + 1]])
        j += 1 + k
    return TriMesh(body, np.array(faces, dtype=np.int64), "vtk")


def read_su2(path) -> TriMesh:
    lines = [ln for ln in Path(path).read_text(encoding="utf-8").splitlines() if not ln.lstrip().startswith("%")]
    pts, faces = [], []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("NPOIN="):
            npoin = int(line.split("=")[1])
            pts = [[float(v) for v in lines[i + 1 + k].split()[:3]] for k in range(npoin)]
            i += 1 + npoin
            continue
        if line.startswith("NELEM="):
            nelem = int(line.split("=")[1])
            for k in range(nelem):
                row = lines[i + 1 + k].split()
                if int(row[0]) != SU2_TRIANGLE:
                    raise ValueError(f"SU2: element type {row[0]} is not a triangle")
                faces.append([int(v) for v in row[1:4]])
            i += 1 + nelem
            continue
        if line.startswith("MARKER_ELEMS="):
            nmark = int(line.split("=")[1])
            for k in range(nmark):
                row = lines[i + 1 + k].split()
                if int(row[0]) != SU2_TRIANGLE:
                    raise ValueError(f"SU2: marker element type {row[0]} is not a triangle")
                faces.append([int(v) for v in row[1:4]])
            i += 1 + nmark
            continue
        i += 1
    if not pts:
        raise ValueError("SU2: NPOIN block missing")
    return TriMesh(np.array(pts, dtype=float), np.array(faces, dtype=np.int64), "su2")


#: nodes per element for the MSH 2.2 element types a surface mesh can contain
_GMSH_NNODE = {1: 2, 2: 3, 3: 4, 15: 1}


def read_gmsh(path) -> TriMesh:
    tokens = Path(path).read_text(encoding="utf-8").split()
    if "$MeshFormat" not in tokens:
        raise ValueError("Gmsh: $MeshFormat missing")
    version = tokens[tokens.index("$MeshFormat") + 1]
    if not version.startswith("2."):
        raise ValueError(f"Gmsh: only MSH 2.x is read, got {version}")
    i = tokens.index("$Nodes")
    nv = int(tokens[i + 1])
    block = tokens[i + 2 : i + 2 + 4 * nv]  # per node: tag x y z
    tags = [int(t) for t in block[0::4]]
    coords = np.array([v for t in block[1::4] for v in (t,)] + block[2::4] + block[3::4], dtype=float)
    coords = np.array(block, dtype=float).reshape(nv, 4)[:, 1:]
    lookup = {tag: k for k, tag in enumerate(tags)}
    i = tokens.index("$Elements")
    ne = int(tokens[i + 1])
    j = i + 2
    faces = []
    for _ in range(ne):
        etype = int(tokens[j + 1])
        ntag = int(tokens[j + 2])
        nn = _GMSH_NNODE.get(etype, 0)
        nodes = tokens[j + 3 + ntag : j + 3 + ntag + nn]
        if etype == GMSH_TRIANGLE:
            faces.append([lookup[int(t)] for t in nodes])
        j += 3 + ntag + nn
    return TriMesh(coords, np.array(faces, dtype=np.int64), "gmsh")


# ---------------------------------------------------------------------------
# CGNS: SIDS conformance checker + reader
# ---------------------------------------------------------------------------
def _attr(node, key: str):
    for spelling in (" " + key, key):
        if spelling in node.attrs:
            return node.attrs[spelling]
    return None


def _text(value) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    if hasattr(value, "tobytes") and not isinstance(value, (int, float)):
        return value.tobytes().decode("utf-8", "replace").rstrip("\x00")
    return str(value)


def _label(node) -> str:
    return _text(_attr(node, "label"))


def _named(node, label: str) -> List:
    return [g for g in node.values() if hasattr(g, "values") and _label(g) == label]


def _data(node):
    return node[" data"][()] if " data" in node else None


#: The SIDS-to-HDF5 rules this writer is held to.  :func:`sids_violations` checks
#: them one by one and returns ``[]`` for a conformant file.
CGNS_RULES: Tuple[str, ...] = (
    "root carries the reserved 'version' and 'format' attributes",
    "root holds a CGNSLibraryVersion_t node whose R4 data is a CGNS 4.x version",
    "a node labelled CGNSBase_t is the database root; every node carries name/label/type/order",
    "Zone_t carries data [nvertex, ncell, 0] and a ZoneType_t child = 'Unstructured'",
    "GridCoordinates_t holds CoordinateX/Y/Z DataArray_t of length nvertex",
    "Elements_t carries data [TRI_3, 0]; ElementRange [1,ncell]; ElementConnectivity 1-based 3*ncell",
)


def _cgns_zone(root):
    bases = _named(root, "CGNSBase_t")
    if not bases:
        raise ValueError("CGNS: no node labelled CGNSBase_t (database root)")
    zones = _named(bases[0], "Zone_t")
    if not zones:
        raise ValueError("CGNS: no node labelled Zone_t")
    return zones[0]


def sids_violations(path) -> List[str]:
    """Check a CGNS/HDF5 file against :data:`CGNS_RULES`; ``[]`` = conformant."""
    h5py = _h5py()
    out: List[str] = []
    with h5py.File(path, "r") as f:
        if _attr(f, "version") is None:
            out.append("root: reserved attribute 'version' missing")
        if _attr(f, "format") is None:
            out.append("root: reserved attribute 'format' missing")

        lib = _named(f, "CGNSLibraryVersion_t")
        if not lib:
            out.append("root: CGNSLibraryVersion_t node missing")
        else:
            raw = _data(lib[0])
            val = np.asarray(raw, dtype=float).ravel() if raw is not None else np.array([])
            if val.size != 1 or not (4.0 <= float(val[0]) <= 5.0):
                out.append(f"root: CGNSLibraryVersion_t data {val.tolist()} is not one CGNS 4.x version")
            if _text(_attr(lib[0], "type")) != "R4":
                out.append(f"root: CGNSLibraryVersion_t type {_text(_attr(lib[0], 'type'))!r} != 'R4'")

        try:
            zone = _cgns_zone(f)
        except ValueError as exc:
            return out + [str(exc)]

        def check_attrs(node) -> None:
            where = node.name.rsplit("/", 1)[-1]
            for req in ("name", "label", "type", "order"):
                if _attr(node, req) is None:
                    out.append(f"{where}: required attribute {req!r} missing")

        for grp in f.values():
            if hasattr(grp, "values"):
                check_attrs(grp)
        for grp in _named(f, "CGNSBase_t")[0].values():
            if hasattr(grp, "values"):
                check_attrs(grp)

        zsize = _data(zone)
        zsize = np.asarray(zsize, dtype=np.int64).ravel() if zsize is not None else np.zeros(0)
        if zsize.size != 3:
            out.append(f"Zone: data {zsize.tolist()} is not [nvertex, ncell, 0]")
            return out
        nv, ncell = int(zsize[0]), int(zsize[1])
        if zsize[2] != 0:
            out.append(f"Zone: VertexSize[2] = {int(zsize[2])} != 0")

        ztypes = _named(zone, "ZoneType_t")
        if not ztypes:
            out.append("Zone: ZoneType_t child missing")
        else:
            zname = _text(_data(ztypes[0]))
            if zname != "Unstructured":
                out.append(f"Zone: ZoneType {zname!r} != 'Unstructured'")

        gcs = _named(zone, "GridCoordinates_t")
        if not gcs:
            out.append("Zone: GridCoordinates_t child missing")
        else:
            kids = {g.name.rsplit("/", 1)[-1]: g for g in gcs[0].values() if hasattr(g, "values")}
            for axis in "XYZ":
                node = kids.get(f"Coordinate{axis}")
                if node is None:
                    out.append(f"GridCoordinates: Coordinate{axis} missing")
                    continue
                if _label(node) != "DataArray_t":
                    out.append(f"GridCoordinates: Coordinate{axis} label {_label(node)!r} != 'DataArray_t'")
                data = _data(node)
                n = -1 if data is None else int(np.asarray(data).size)
                if n != nv:
                    out.append(f"GridCoordinates: Coordinate{axis} size {n} != nvertex {nv}")

        els = _named(zone, "Elements_t")
        if not els:
            out.append("Zone: Elements_t child missing")
            return out
        el = els[0]
        etype = _data(el)
        etype = np.asarray(etype, dtype=np.int64).ravel() if etype is not None else np.zeros(0)
        if etype.size != 2 or int(etype[0]) != SU2_TRIANGLE or int(etype[1]) != 0:
            out.append(f"Elements: data {etype.tolist()} is not [TRI_3(=5), 0]")
        kids = {g.name.rsplit("/", 1)[-1]: g for g in el.values() if hasattr(g, "values")}
        if "ElementRange" not in kids:
            out.append("Elements: ElementRange missing")
        else:
            rng = _data(kids["ElementRange"])
            rng = [int(v) for v in np.asarray(rng, dtype=np.int64).ravel()] if rng is not None else []
            if rng != [1, ncell]:
                out.append(f"Elements: ElementRange {rng} != [1, {ncell}]")
        if "ElementConnectivity" not in kids:
            out.append("Elements: ElementConnectivity missing")
        else:
            raw = _data(kids["ElementConnectivity"])
            conn = np.asarray(raw, dtype=np.int64).ravel() if raw is not None else np.zeros(0)
            if conn.size != 3 * ncell:
                out.append(f"Elements: ElementConnectivity size {conn.size} != 3*ncell = {3 * ncell}")
            elif conn.size and (conn.min() < 1 or conn.max() > nv):
                out.append(f"Elements: ElementConnectivity not 1-based within [1, {nv}]")
    return out


def read_cgns(path) -> TriMesh:
    """Read the SIDS dialect; a non-conformant file raises instead of being laundered."""
    violations = sids_violations(path)
    if violations:
        raise ValueError("CGNS: non-conformant file -> " + "; ".join(violations))
    h5py = _h5py()
    with h5py.File(path, "r") as f:
        zone = _cgns_zone(f)
        nv = int(np.asarray(_data(zone), dtype=np.int64).ravel()[0])
        gc = _named(zone, "GridCoordinates_t")[0]
        kids = {g.name.rsplit("/", 1)[-1]: g for g in gc.values() if hasattr(g, "values")}
        coords = [np.asarray(_data(kids[f"Coordinate{axis}"]), dtype=float).ravel() for axis in "XYZ"]
        el = _named(zone, "Elements_t")[0]
        ekids = {g.name.rsplit("/", 1)[-1]: g for g in el.values() if hasattr(g, "values")}
        faces = np.asarray(_data(ekids["ElementConnectivity"]), dtype=np.int64).ravel().reshape(-1, 3) - 1
    return TriMesh(np.stack(coords, axis=-1), faces, "cgns")


READERS = {
    "stl": read_stl,
    "obj": read_obj,
    "vtk": read_vtk,
    "su2": read_su2,
    "cgns": read_cgns,
    "gmsh": read_gmsh,
}


def read(path, fmt: str) -> TriMesh:
    """Read ``path`` as ``fmt`` with this module's own reader."""
    return READERS[canonical_format(fmt)](Path(path))


# ---------------------------------------------------------------------------
# comparison contract
# ---------------------------------------------------------------------------
def _bijection(ref: np.ndarray, got: np.ndarray, tol: float) -> Optional[np.ndarray]:
    """Map ``got`` vertices onto ``ref`` indices, or ``None`` when impossible."""
    if len(ref) != len(got):
        return None
    mapping = np.full(len(got), -1, dtype=np.int64)
    used = np.zeros(len(ref), dtype=bool)
    by_key: Dict[tuple, int] = {}
    for i, p in enumerate(ref):
        by_key.setdefault(tuple(np.round(p, 12)), i)
    if tol <= 0.0:
        for j, p in enumerate(got):
            i = by_key.get(tuple(np.round(p, 12)))
            if i is None or used[i]:
                return None
            used[i] = True
            mapping[j] = i
        return mapping
    for j, p in enumerate(got):
        i = by_key.get(tuple(np.round(p, 12)))
        if i is None or used[i]:
            live = np.flatnonzero(~used)
            if live.size == 0:
                return None
            dist = np.linalg.norm(ref[live] - p, axis=1)
            if float(dist.min()) > tol:
                return None
            i = int(live[int(np.argmin(dist))])
            if used[i]:
                return None
        used[i] = True
        mapping[j] = i
    return mapping


def _oriented_keys(points: np.ndarray, faces: np.ndarray, mapping=None) -> List[Tuple[int, int, int]]:
    f = faces if mapping is None else mapping[faces]
    out = []
    for a, b, c in np.asarray(f, dtype=np.int64):
        t = (int(a), int(b), int(c))
        i = t.index(min(t))
        out.append(t[i:] + t[:i])  # cyclic rotation keeps orientation
    return out


def compare_meshes(ref: TriMesh, got: TriMesh, tol: float = 0.0) -> dict:
    """Geometry + topology + orientation comparison with a vertex bijection.

    ``tol`` is an absolute length; ``0`` demands an exact bijection (right for the
    float64 formats).  STL is a facet soup, so when the face count matches but the
    vertex count does not, ``got`` is welded first and ``welded`` reports it.
    ``pass`` is the machine gate: bijection, equal face count, equal face multiset,
    equal *oriented* face multiset, equal volume sign and a closed surface.
    """
    welded_from = int(len(got.points))
    if len(got.points) != len(ref.points) and len(got.faces) == len(ref.faces):
        got, _merged = weld(got)
    v_ref, v_got = ref.volume(), got.volume()
    denom = abs(v_ref) if abs(v_ref) > 0 else 1.0
    out: dict = {
        "vertices": [int(len(ref.points)), int(len(got.points))],
        "vertices_read": [int(len(ref.points)), welded_from],
        "faces": [int(len(ref.faces)), int(len(got.faces))],
        "welded": bool(welded_from != len(got.points)),
        "volume_ref": v_ref,
        "volume_got": v_got,
        "volume_rel_err": abs(v_got - v_ref) / denom,
        "euler": [ref.euler(), got.euler()],
        "boundary_edges": [ref.boundary_edges(), got.boundary_edges()],
    }
    mapping = _bijection(ref.points, got.points, tol)
    out["bijection"] = mapping is not None
    if mapping is not None:
        key_ref = set(_oriented_keys(ref.points, ref.faces))
        key_got = _oriented_keys(got.points, got.faces, mapping)
        out["face_multiset_equal"] = sorted(
            tuple(sorted(f)) for f in ref.faces.tolist()
        ) == sorted(tuple(sorted(mapping[f].tolist())) for f in got.faces.tolist())
        out["oriented_face_multiset_equal"] = key_ref == set(key_got)
        out["flipped_faces"] = int(len(set(key_got) - key_ref))
        out["unmatched_faces"] = int(len(key_ref - set(key_got)))
    else:
        out["face_multiset_equal"] = False
        out["oriented_face_multiset_equal"] = False
        out["flipped_faces"] = -1
        out["unmatched_faces"] = -1
    out["volume_sign_equal"] = bool(v_ref * v_got == 0 or np.sign(v_ref) == np.sign(v_got))
    out["closed"] = out["boundary_edges"][1] == 0
    out["pass"] = bool(
        out["bijection"]
        and out["faces"][0] == out["faces"][1]
        and out["face_multiset_equal"]
        and out["oriented_face_multiset_equal"]
        and out["volume_sign_equal"]
        and out["closed"]
    )
    return out


def roundtrip(mesh: TriMesh, fmt: str, path, tol: float = 0.0, **kwargs) -> dict:
    """Write then read one format and return the comparison report (+ the path)."""
    fmt = canonical_format(fmt)
    write(mesh, fmt, path, **kwargs)
    got = read(path, fmt)
    report = compare_meshes(mesh, got, tol=tol)
    report["format"] = fmt
    report["path"] = str(path)
    report["bytes"] = int(Path(path).stat().st_size)
    return report
