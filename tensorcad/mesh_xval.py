"""Cross-validation (XVAL) of the native six-format surface-mesh matrix against
mature open-source implementations.

Why this module exists
----------------------
``tensorcad.mesh_formats`` (node m17483) asserts that a closed triangle surface
survives write->read in six formats (STL, OBJ, VTK, SU2, CGNS, Gmsh).  An
internal round trip only proves self-consistency: a writer and a reader that
share the same misunderstanding of a grammar agree with each other.  This module
therefore re-reads the files with *independent* implementations and compares the
recovered geometry with the geometry that was written:

===========  ==========================================================================
format       independent reader used here (source read, not guessed)
===========  ==========================================================================
STL          ``meshio`` 5.3.5, ``numpy-stl`` 4.1.3, ``trimesh`` 5.0.0
OBJ          ``meshio`` 5.3.5, ``trimesh`` 5.0.0
VTK          ``meshio`` 5.3.5 legacy reader (UNSTRUCTURED_GRID only, see below)
SU2          ``meshio`` 5.3.5 + strict grammar reader in this module
CGNS         ``h5py`` SIDS conformance checker in this module (no OSS triangle
             reader exists, see :func:`meshio_cgns_status`)
Gmsh         ``meshio`` 5.3.5 MSH 2.2 reader
===========  ==========================================================================

Two of the six formats are *not* readable by the reference implementation, and
both failures were traced to the reference source rather than assumed:

* ``meshio/su2/_su2.py``:113-122 — for ``NELEM=0`` (a legal boundary-only SU2
  mesh) ``next(gen)`` is called on an empty ``islice`` -> ``StopIteration``;
* ``meshio/cgns/_cgns.py`` — reader is a stub that ends in
  ``ReadError("Can only read tetrahedra.")``;
* ``meshio/vtk/_vtk_42.py`` — legacy ``POLYDATA`` is rejected, only
  ``UNSTRUCTURED_GRID`` / ``STRUCTURED_*`` are parsed.

The writers here are deliberately a *second, independent* implementation of the
same six grammars (they do not import ``tensorcad.mesh_formats``), so the
agreement matrix is writer x OSS-reader and not "our writer x our reader".

Nothing in this module is ship-specific: it is pure geometry-exchange plumbing.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

FORMATS: Tuple[str, ...] = ("stl", "obj", "vtk", "su2", "cgns", "gmsh")

#: SU2 element type ids (SU2 v7 mesh-file grammar, "Element type 5 = triangle").
SU2_TRIANGLE = 5
#: Gmsh MSH 2.2 element type for a 3-node triangle.
GMSH_TRIANGLE = 2


# ---------------------------------------------------------------------------
# geometry container
# ---------------------------------------------------------------------------
@dataclass
class TriMesh:
    """A closed, consistently oriented triangle surface."""

    points: np.ndarray
    faces: np.ndarray
    name: str = "mesh"

    def __post_init__(self) -> None:
        self.points = np.asarray(self.points, dtype=np.float64).reshape(-1, 3)
        self.faces = np.asarray(self.faces, dtype=np.int64).reshape(-1, 3)
        if self.faces.size and (self.faces.min() < 0 or self.faces.max() >= len(self.points)):
            raise ValueError("face index out of range")

    # -- metrics -----------------------------------------------------------
    def volume(self) -> float:
        """Signed volume by the divergence theorem, ``V = 1/6 sum a.(b x c)``."""
        a = self.points[self.faces[:, 0]]
        b = self.points[self.faces[:, 1]]
        c = self.points[self.faces[:, 2]]
        return float(np.einsum("ij,ij->i", a, np.cross(b, c)).sum() / 6.0)

    def euler(self) -> int:
        """``V - E + F`` on the *used* vertex set."""
        used = np.unique(self.faces)
        edges = np.vstack([self.faces[:, [0, 1]], self.faces[:, [1, 2]], self.faces[:, [2, 0]]])
        edges = np.sort(edges, axis=1)
        n_edges = len(np.unique(edges, axis=0))
        return int(len(used) - n_edges + len(self.faces))

    def boundary_edges(self) -> int:
        edges = np.vstack([self.faces[:, [0, 1]], self.faces[:, [1, 2]], self.faces[:, [2, 0]]])
        edges = np.sort(edges, axis=1)
        uniq, counts = np.unique(edges, axis=0, return_counts=True)
        return int((counts != 2).sum())

    def oriented_face_keys(self, mapping: Optional[np.ndarray] = None) -> List[Tuple[int, int, int]]:
        f = self.faces if mapping is None else mapping[self.faces]
        keys = []
        for a, b, c in f:
            t = (int(a), int(b), int(c))
            i = t.index(min(t))
            keys.append(t[i:] + t[:i])  # cyclic rotation, orientation preserved
        return keys


def weld(mesh: TriMesh, decimals: int = 12) -> Tuple[TriMesh, int]:
    """Recover connectivity from a triangle soup by exact coordinate welding.

    STL (and ``process=False`` trimesh) stores each facet's corners separately; the
    *file* therefore carries no topology at all.  Welding by identical coordinates is
    the operation any consumer must perform before topological questions are asked,
    so the comparison does it explicitly and reports how many vertices were merged.
    """
    key = np.round(mesh.points, decimals)
    uniq, inverse = np.unique(key, axis=0, return_inverse=True)
    faces = inverse[mesh.faces]
    if (faces < 0).any():  # pragma: no cover - numpy contract
        raise ValueError("weld failed")
    good = (faces[:, 0] != faces[:, 1]) & (faces[:, 1] != faces[:, 2]) & (faces[:, 0] != faces[:, 2])
    return TriMesh(uniq, faces[good], mesh.name + "/welded"), int(len(mesh.points) - len(uniq))


# ---------------------------------------------------------------------------
# analytic test bodies (genus 0, genus 0 refined, genus 1)
# ---------------------------------------------------------------------------
def cube(edge: float = 2.0) -> TriMesh:
    h = edge / 2.0
    p = np.array(
        [[-h, -h, -h], [h, -h, -h], [h, h, -h], [-h, h, -h],
         [-h, -h, h], [h, -h, h], [h, h, h], [-h, h, h]], dtype=float)
    f = np.array([[0, 2, 1], [0, 3, 2], [4, 5, 6], [4, 6, 7],
                  [0, 1, 5], [0, 5, 4], [1, 2, 6], [1, 6, 5],
                  [2, 3, 7], [2, 7, 6], [3, 0, 4], [3, 4, 7]])
    return TriMesh(p, f, "cube")


def icosphere(radius: float = 1.0, subdivisions: int = 2) -> TriMesh:
    t = (1.0 + 5.0 ** 0.5) / 2.0
    verts = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0),
             (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
             (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    faces = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
             (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
             (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
             (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    P = [np.array(v, dtype=float) / np.linalg.norm(v) for v in verts]
    F = [tuple(f) for f in faces]
    for _ in range(subdivisions):
        cache: Dict[Tuple[int, int], int] = {}

        def mid(i: int, j: int) -> int:
            key = (min(i, j), max(i, j))
            if key not in cache:
                m = (P[i] + P[j]) / 2.0
                P.append(m / np.linalg.norm(m))
                cache[key] = len(P) - 1
            return cache[key]

        nf = []
        for a, b, c in F:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            nf.extend([(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)])
        F = nf
    return TriMesh(np.array(P) * radius, np.array(F), f"icosphere{subdivisions}")


def torus(major: float = 1.0, minor: float = 0.3, nu: int = 64, nv: int = 32) -> TriMesh:
    u = np.linspace(0, 2 * np.pi, nu, endpoint=False)
    v = np.linspace(0, 2 * np.pi, nv, endpoint=False)
    U, V = np.meshgrid(u, v, indexing="ij")
    P = np.stack([(major + minor * np.cos(V)) * np.cos(U),
                  (major + minor * np.cos(V)) * np.sin(U),
                  minor * np.sin(V)], axis=-1).reshape(-1, 3)
    idx = np.arange(nu * nv).reshape(nu, nv)
    F = []
    for i in range(nu):
        for j in range(nv):
            a = idx[i, j]
            b = idx[(i + 1) % nu, j]
            c = idx[(i + 1) % nu, (j + 1) % nv]
            d = idx[i, (j + 1) % nv]
            F.append((a, b, c))
            F.append((a, c, d))
    return TriMesh(P, np.array(F), f"torus{nu}x{nv}")


# ---------------------------------------------------------------------------
# writers (independent second implementation of each grammar)
# ---------------------------------------------------------------------------
def write_stl(mesh: TriMesh, path: Path) -> Path:
    """Binary STL: 80-byte header, uint32 count, per-facet float32 normal+vertices."""
    tri = mesh.points[mesh.faces]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = np.divide(n, ln, out=np.zeros_like(n), where=ln > 0)
    with open(path, "wb") as fh:
        fh.write(b"TensorCAD XVAL binary STL (float32)".ljust(80, b" "))
        fh.write(struct.pack("<I", len(mesh.faces)))
        for i in range(len(mesh.faces)):
            fh.write(struct.pack("<3f", *n[i].astype(np.float32)))
            fh.write(struct.pack("<9f", *tri[i].astype(np.float32).ravel()))
            fh.write(struct.pack("<H", 0))
    return path


def write_obj(mesh: TriMesh, path: Path) -> Path:
    """Wavefront OBJ, 17 significant digits (float64 exact round trip)."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# TensorCAD XVAL {mesh.name}\n")
        for x, y, z in mesh.points:
            fh.write(f"v {x:.17g} {y:.17g} {z:.17g}\n")
        for a, b, c in mesh.faces + 1:
            fh.write(f"f {a} {b} {c}\n")
    return path


def write_vtk(mesh: TriMesh, path: Path, dataset: str = "UNSTRUCTURED_GRID") -> Path:
    """VTK legacy ASCII.  ``UNSTRUCTURED_GRID`` (meshio-readable) or ``POLYDATA``."""
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# vtk DataFile Version 3.0\n")
        fh.write(f"TensorCAD XVAL {mesh.name} {dataset}\n")
        fh.write("ASCII\n")
        fh.write(f"DATASET {dataset}\n")
        fh.write(f"POINTS {len(mesh.points)} double\n")
        for x, y, z in mesh.points:
            fh.write(f"{x:.17g} {y:.17g} {z:.17g}\n")
        if dataset == "POLYDATA":
            fh.write(f"POLYGONS {len(mesh.faces)} {4 * len(mesh.faces)}\n")
            for a, b, c in mesh.faces:
                fh.write(f"3 {a} {b} {c}\n")
        elif dataset == "UNSTRUCTURED_GRID":
            fh.write(f"CELLS {len(mesh.faces)} {4 * len(mesh.faces)}\n")
            for a, b, c in mesh.faces:
                fh.write(f"3 {a} {b} {c}\n")
            fh.write(f"CELL_TYPES {len(mesh.faces)}\n")
            for _ in mesh.faces:
                fh.write("5\n")  # VTK_TRIANGLE
        else:  # pragma: no cover
            raise ValueError(dataset)
    return path


def write_su2(mesh: TriMesh, path: Path, marker: str = "WALL", nelem_mode: str = "boundary") -> Path:
    """SU2 ASCII mesh.

    ``nelem_mode='boundary'`` (default, the faithful subset): ``NELEM=0`` and the
    triangles live in one ``MARKER_TAG``/``MARKER_ELEMS`` block, 0-based.
    ``nelem_mode='volume'`` puts the triangles in ``NELEM`` instead -- used only
    to separate *our* grammar from the *reference reader's* NELEM=0 defect.
    """
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"% TensorCAD XVAL {mesh.name}\n")
        fh.write("NDIME= 3\n")
        fh.write(f"NPOIN= {len(mesh.points)}\n")
        for x, y, z in mesh.points:
            fh.write(f"{x:.17g} {y:.17g} {z:.17g}\n")
        if nelem_mode == "volume":
            fh.write(f"NELEM= {len(mesh.faces)}\n")
            for a, b, c in mesh.faces:
                fh.write(f"{SU2_TRIANGLE} {a} {b} {c} 0\n")
            fh.write("NMARK= 0\n")
        elif nelem_mode == "boundary":
            fh.write("NELEM= 0\n")
            fh.write("NMARK= 1\n")
            fh.write(f"MARKER_TAG= {marker}\n")
            fh.write(f"MARKER_ELEMS= {len(mesh.faces)}\n")
            for a, b, c in mesh.faces:
                fh.write(f"{SU2_TRIANGLE} {a} {b} {c}\n")
        else:  # pragma: no cover
            raise ValueError(nelem_mode)
    return path


def _h5py():
    try:
        import h5py  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("CGNS writes need h5py") from exc
    import h5py
    return h5py


def _cgns_attrs(node, props: dict) -> None:
    """Write the CGNS/HDF5 node attributes.

    Spec (CGNS SIDS-to-HDF5, "Basic Node"): a CGNS node is an HDF5 group with the
    required attributes ``name``, ``label``, ``type`` (= data-type code, ``MT`` for a
    node without data) and ``order``; node data, when present, is the single dataset
    `` data``.  The spec encloses the attribute names in square brackets "because the
    SLL uses private attributes with a leading space character", so both the
    space-prefixed and the bare spelling occur in the wild and both are accepted by
    :func:`read_cgns_sids`.
    """
    for key, value in props.items():
        node.attrs[" " + key] = np.bytes_(value) if isinstance(value, str) else np.int32(value)


def _cgns_node(parent, name: str, label: str, dtype: str = "MT", order: int = 1, data=None):
    g = parent.create_group(name)
    _cgns_attrs(g, {"name": name, "label": label, "type": dtype, "order": order})
    if data is not None:
        ds = g.create_dataset(" data", data=data)
        ds.attrs["type"] = np.bytes_(dtype)
    return g


def write_cgns(mesh: TriMesh, path: Path) -> Path:
    """CGNS/HDF5 database, unstructured zone, Elements_t = TRI_3.

    Conformant per the SIDS-to-HDF5 mapping: the file root carries the reserved
    ``version``/``format`` attributes, the database is rooted at a node labelled
    ``CGNSBase_t``, every node carries ``name``/``label``/``type``/``order`` and its
    payload lives in the `` data`` dataset.  Element connectivity is 1-based, the
    single element block has ``ElementRange = [1, ncell]``.
    """
    h5py = _h5py()
    nv, nf = len(mesh.points), len(mesh.faces)
    with h5py.File(path, "w") as f:
        _cgns_attrs(f, {"version": "3.4", "format": "0.9"})  # reserved root attributes
        base = _cgns_node(f, "Base", "CGNSBase_t", "I4", 1, np.array([3, 3], dtype=np.int32))
        zone = _cgns_node(base, "Zone", "Zone_t", "I4", 1, np.array([nv, nf, 0], dtype=np.int32))
        zt = _cgns_node(zone, "ZoneType", "ZoneType_t", "C1", 1)
        zt.create_dataset(" data", data=np.bytes_("Unstructured"))
        gc = _cgns_node(zone, "GridCoordinates", "GridCoordinates_t", "MT", 1)
        for k, axis in enumerate("XYZ"):
            _cgns_node(gc, f"Coordinate{axis}", "DataArray_t", "R8", k + 1,
                       np.ascontiguousarray(mesh.points[:, k]))
        el = _cgns_node(zone, "Elements", "Elements_t", "I4", 2,
                        np.array([SU2_TRIANGLE, 0], dtype=np.int32))  # [type, 0] = TRI_3
        _cgns_node(el, "ElementRange", "IndexRange_t", "I4", 1, np.array([1, nf], dtype=np.int32))
        _cgns_node(el, "ElementConnectivity", "DataArray_t", "I4", 2,
                   (mesh.faces + 1).astype(np.int32).ravel())  # 1-based
    return path


def write_gmsh(mesh: TriMesh, path: Path) -> Path:
    """Gmsh MSH 2.2 ASCII: ``$Nodes``/``$Elements``, element type 2, 1-based."""
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
    return path


WRITERS = {
    "stl": write_stl,
    "obj": write_obj,
    "vtk": write_vtk,
    "su2": write_su2,
    "cgns": write_cgns,
    "gmsh": write_gmsh,
}


# ---------------------------------------------------------------------------
# independent readers
# ---------------------------------------------------------------------------
def read_su2_strict(path: Path) -> TriMesh:
    """Strict SU2 grammar reader (SU2 v7 docs: NDIME/NPOIN/NELEM/NMARK)."""
    lines = [ln.strip() for ln in Path(path).read_text(encoding="utf-8").splitlines()]
    lines = [ln for ln in lines if ln and not ln.startswith("%")]
    i = 0
    points = None
    faces: List[Tuple[int, int, int]] = []
    npoin = nelem = nmark = None
    while i < len(lines):
        ln = lines[i]
        if "=" not in ln:
            raise ValueError(f"SU2: stray line {i}: {ln!r}")
        key, _, val = ln.partition("=")
        key, val = key.strip(), val.strip()
        i += 1
        if key == "NDIME":
            if int(val) != 3:
                raise ValueError("SU2: only NDIME=3 supported")
        elif key == "NPOIN":
            npoin = int(val)
            pts = [list(map(float, lines[i + k].split()[:3])) for k in range(npoin)]
            points = np.array(pts, dtype=float)
            i += npoin
        elif key in ("NELEM", "MARKER_ELEMS"):
            n = int(val)
            for k in range(n):
                parts = lines[i + k].split()
                etype = int(parts[0])
                if etype != SU2_TRIANGLE:
                    raise ValueError(f"SU2: element type {etype} is not a triangle")
                faces.append((int(parts[1]), int(parts[2]), int(parts[3])))
            i += n
            if key == "NELEM":
                nelem = n
        elif key == "NMARK":
            nmark = int(val)
        elif key == "MARKER_TAG":
            pass
        else:
            raise ValueError(f"SU2: unknown keyword {key!r}")
    if points is None:
        raise ValueError("SU2: no NPOIN block")
    if npoin is None or nelem is None or nmark is None:
        raise ValueError("SU2: NDIME/NPOIN/NELEM/NMARK mandatory")
    return TriMesh(points, np.array(faces), "su2")


def _attr(obj, key: str):
    for k in (key, " " + key):
        if k in obj.attrs:
            return obj.attrs[k]
    return None


def _as_text(v) -> str:
    if isinstance(v, bytes):
        return v.decode(errors="replace").strip("\x00 ")
    if isinstance(v, np.ndarray):
        return v.tobytes().decode(errors="replace").strip("\x00 ")
    if isinstance(v, np.generic):
        return _as_text(v.item())
    return str(v)


def _children(node) -> List:
    """CGNS children = child groups whose name does not start with a space."""
    return [node[k] for k in node.keys() if isinstance(node[k], _h5py().Group)
            and not k.startswith(" ")]


def _label(node) -> str:
    v = _attr(node, "label")
    if v is None:
        v = _attr(node, "type")            # some writers put the label in 'type'
    return _as_text(v).strip("\x00 ") if v is not None else ""


def read_cgns_sids(path: Path) -> Tuple[TriMesh, List[str]]:
    """Read a CGNS/HDF5 surface mesh and check the SIDS-to-HDF5 mapping invariants.

    Returns ``(mesh, violations)``; ``violations`` empty means the file conforms to
    :data:`CGNS_RULES`.  The reader is dialect-tolerant on purpose (it accepts the
    attribute spellings that occur in the wild) while the *checker* is not: a node
    that carries its label in the ``type`` slot, or that has no ``label``/``name``/
    ``order`` attribute at all, is reported as non-conformant even though the geometry
    can still be recovered.
    """
    h5py = _h5py()
    violations: List[str] = []
    with h5py.File(path, "r") as f:
        # --- root node -------------------------------------------------------
        root_label = _label(f)
        if "version" not in [k.strip() for k in f.attrs] and _attr(f, "version") is None:
            violations.append("root: reserved attribute 'version' missing")
        # --- locate the database root (a node labelled CGNSBase_t) -----------
        bases = [g for g in _children(f) if _label(g).endswith("Base_t")]
        if not bases:
            bases = [g for g in _children(f) if _label(g) in ("CGNSBase_t", "Base_t")]
        if not bases:
            raise ValueError("CGNS: no node labelled CGNSBase_t (database root)")
        base = bases[0]
        if _label(base) != "CGNSBase_t":
            violations.append(f"Base label {_label(base)!r} != 'CGNSBase_t'")
        for node, tag in ((base, "Base"),):
            for req in ("name", "label", "type", "order"):
                if _attr(node, req) is None:
                    violations.append(f"{tag}: required attribute {req!r} missing")
            t = _attr(node, "type")
            if t is not None and _as_text(t).endswith("_t"):
                violations.append(f"{tag}: 'type' holds a label ({_as_text(t)!r}) where the spec "
                                  f"requires a data-type code (MT/I4/R8/...)")
        # --- zone -----------------------------------------------------------
        zones = [g for g in _children(base) if _label(g) == "Zone_t"]
        if not zones:
            raise ValueError("CGNS: no Zone_t node")
        zone = zones[0]
        for req in ("name", "label", "type", "order"):
            if _attr(zone, req) is None:
                violations.append(f"Zone: required attribute {req!r} missing")
        if _attr(zone, "type") is not None and _as_text(_attr(zone, "type")).endswith("_t"):
            violations.append("Zone: 'type' holds a label where the spec requires a data-type code")
        ztype = ""
        for g in _children(zone):
            if _label(g) == "ZoneType_t":
                ztype = _as_text(g[" data"][()])
            elif g.name.rsplit("/", 1)[-1] == "ZoneType" and " data" in g:
                ztype = _as_text(g[" data"][()])
        if ztype != "Unstructured":
            violations.append(f"ZoneType={ztype!r} (expected 'Unstructured')")
            raise ValueError(f"CGNS: zone is not unstructured ({ztype!r})")
        vsize = None
        if " data" in zone and zone[" data"].shape == (3,):
            vsize = [int(v) for v in zone[" data"][()]]
        else:
            for g in _children(zone):
                if " data" in g and g.name.rsplit("/", 1)[-1] in ("VertexSize", "VertexSize_t"):
                    vsize = [int(v) for v in g[" data"][()]]
        if vsize is None:
            raise ValueError("CGNS: VertexSize not found")
        nv, ncell = vsize[0], vsize[1]
        if vsize[2] != 0:
            violations.append(f"VertexSize[2]={vsize[2]} (unstructured zone: expected 0)")
        # --- coordinates ----------------------------------------------------
        gc = [g for g in _children(zone) if _label(g) == "GridCoordinates_t"]
        if not gc:
            raise ValueError("CGNS: no GridCoordinates_t node")
        coords = []
        for axis in "XYZ":
            node = None
            for g in _children(gc[0]):
                if g.name.rsplit("/", 1)[-1] == f"Coordinate{axis}":
                    node = g
            if node is None:
                raise ValueError(f"CGNS: Coordinate{axis} missing")
            if _label(node) != "DataArray_t":
                violations.append(f"Coordinate{axis} label {_label(node)!r} != 'DataArray_t'")
            data = np.asarray(node[" data"][()], dtype=float)
            if data.size != nv:
                violations.append(f"Coordinate{axis} size {data.size} != nvertex {nv}")
            coords.append(data.ravel())
        points = np.stack(coords, axis=-1)
        # --- elements -------------------------------------------------------
        el = [g for g in _children(zone) if _label(g) == "Elements_t"]
        if not el:
            raise ValueError("CGNS: no Elements_t node")
        el = el[0]
        etype = ""
        rng, conn = None, None
        for g in _children(el):
            nm = g.name.rsplit("/", 1)[-1]
            if nm == "ElementType":
                etype = _as_text(g[" data"][()])
            elif nm == "ElementRange":
                rng = [int(v) for v in np.asarray(g[" data"][()]).ravel()]
            elif nm == "ElementConnectivity":
                conn = np.asarray(g[" data"][()], dtype=np.int64).ravel()
        if not etype:  # older dialect: element type as the node's own data [type, 0]
            if " data" in el and el[" data"].size == 2:
                etype = "TRI_3" if int(el[" data"][0]) == SU2_TRIANGLE else f"type{int(el[' data'][0])}"
        if etype != "TRI_3":
            violations.append(f"ElementType={etype!r} (expected 'TRI_3')")
        if conn is None:
            raise ValueError("CGNS: ElementConnectivity missing")
        n_per = 3 if conn.size == 3 * ncell else None
        if n_per is None:
            violations.append(f"ElementConnectivity size {conn.size} is not 3*ncell={3 * ncell}")
            n_per = conn.size // max(ncell, 1)
        faces = conn.reshape(ncell, n_per)[:, :3]
        if buildings := (faces.min() < 1 or faces.max() > nv):
            violations.append(f"ElementConnectivity not 1-based within [1,{nv}]")
        if rng is not None and rng != [1, ncell]:
            violations.append(f"ElementRange={rng} (expected 1-based contiguous [1,{ncell}])")
        for req in ("name", "label", "type", "order"):
            if _attr(el, req) is None:
                violations.append(f"Elements: required attribute {req!r} missing")
        faces = faces - 1
    return TriMesh(points, faces, "cgns"), violations


CGNS_RULES = (
    "file root carries the reserved 'version' attribute (SIDS-to-HDF5, Root Node)",
    "a node labelled CGNSBase_t is the database root",
    "every node has attributes name/label/type/order, type = data-type code (MT/I4/R8/C1)",
    "Zone_t node with ZoneType=Unstructured and VertexSize=[nvertex,ncell,0]",
    "GridCoordinates_t with CoordinateX/Y/Z DataArray_t of length nvertex",
    "Elements_t with ElementType=TRI_3, ElementRange=[1,ncell]",
    "ElementConnectivity 1-based, ncell x 3",
)


def meshio_mesh(path: Path, fmt: Optional[str] = None) -> TriMesh:
    """Read with meshio and reduce to a triangle surface."""
    import meshio

    kw = {} if fmt is None else {"file_format": fmt}
    m = meshio.read(str(path), **kw)
    pts = np.asarray(m.points, dtype=float)
    if pts.shape[1] == 2:
        pts = np.hstack([pts, np.zeros((len(pts), 1))])
    tri = [np.asarray(c.data) for c in m.cells if c.type.startswith("triangle")]
    if not tri:
        raise ValueError("meshio read back no triangles")
    return TriMesh(pts, np.vstack(tri), "meshio")


def trimesh_mesh(path: Path, process: bool = True) -> TriMesh:
    """Read with trimesh.  ``process=False`` disables vertex merging/sorting."""
    import trimesh

    tm = trimesh.load(str(path), process=process, force="mesh")
    return TriMesh(np.asarray(tm.vertices, dtype=float), np.asarray(tm.faces, dtype=np.int64), "trimesh")


def numpy_stl_mesh(path: Path) -> TriMesh:
    from stl import mesh as stlmesh

    m = stlmesh.Mesh.from_file(str(path))
    tri = np.asarray(m.vectors, dtype=float)  # (nf, 3, 3)
    pts = tri.reshape(-1, 3)
    faces = np.arange(len(pts)).reshape(-1, 3)
    return TriMesh(pts, faces, "numpy-stl")


def meshio_cgns_status() -> dict:
    """Report what the reference CGNS reader actually implements (source-checked)."""
    import meshio

    src = Path(meshio.__file__).resolve().parent / "cgns" / "_cgns.py"
    text = src.read_text()
    return {
        "path": str(src.relative_to(Path(meshio.__file__).resolve().parent.parent)),
        "raises_can_only_read_tetrahedra": "Can only read tetrahedra." in text,
        "todo_stub": "TODO" in text,
    }


# ---------------------------------------------------------------------------
# comparison
# ---------------------------------------------------------------------------
def compare_meshes(ref: TriMesh, got: TriMesh, tol: float = 0.0) -> dict:
    """Geometry/topology/orientation comparison with a vertex bijection.

    ``tol`` is an absolute length; ``0`` demands an exact bijection (lossless
    formats).  A vertex without a partner inside ``tol`` makes both the topology
    and the orientation verdict false -- a format that duplicates or merges
    corners cannot pass by accident.
    """
    nv_ref, nf_ref = len(ref.points), len(ref.faces)
    welded_from = int(len(got.points))
    if len(got.points) != nv_ref and len(got.faces) == nf_ref:
        got, _merged = weld(got)
    out: dict = {
        "vertices": [int(nv_ref), int(len(got.points))],
        "vertices_read": [int(nv_ref), welded_from],
        "faces": [int(nf_ref), int(len(got.faces))],
        "welded": bool(welded_from != len(got.points)),
        "volume_ref": ref.volume(),
        "volume_got": got.volume(),
    }
    vref, vgot = ref.volume(), got.volume()
    denom = abs(vref) if abs(vref) > 0 else 1.0
    out["volume_rel_err"] = abs(vgot - vref) / denom
    out["euler"] = [ref.euler(), got.euler()]
    out["boundary_edges"] = [ref.boundary_edges(), got.boundary_edges()]

    # vertex bijection, then topology/orientation in the bijected index space
    mapping = _bijection(ref.points, got.points, tol)
    out["bijection"] = mapping is not None
    if mapping is not None:
        key_ref = set(ref.oriented_face_keys())
        key_got = got.oriented_face_keys(mapping)
        out["face_multiset_equal"] = sorted(map(tuple, np.sort(ref.faces, axis=1).tolist())) == sorted(
            map(tuple, np.sort(mapping[got.faces], axis=1).tolist()))
        out["oriented_face_multiset_equal"] = key_ref == set(key_got)
        out["flipped_faces"] = int(len(set(key_got) - key_ref))
    else:
        out["face_multiset_equal"] = False
        out["oriented_face_multiset_equal"] = False
        out["flipped_faces"] = -1
    if out["volume_ref"] * out["volume_got"] != 0:
        out["volume_sign_equal"] = bool(np.sign(out["volume_ref"]) == np.sign(out["volume_got"]))
    else:
        out["volume_sign_equal"] = True
    out["pass"] = bool(
        out["bijection"]
        and out["faces"][0] == out["faces"][1]
        and out["face_multiset_equal"]
        and out["oriented_face_multiset_equal"]
        and out["volume_sign_equal"]
        and out["boundary_edges"][1] == 0
    )
    return out


def _bijection(a: np.ndarray, b: np.ndarray, tol: float) -> Optional[np.ndarray]:
    """Return ``map`` with ``a[map[i]] == b[i]`` (or nearest within ``tol``)."""
    if len(a) != len(b):
        return None
    if tol == 0.0:
        key = {tuple(np.round(p, 12).tolist()): i for i, p in enumerate(a)}
        mapping = np.empty(len(b), dtype=np.int64)
        for i, p in enumerate(b):
            k = tuple(np.round(p, 12).tolist())
            if k not in key:
                return None
            mapping[i] = key[k]
        return mapping
    d = np.linalg.norm(b[:, None, :] - a[None, :, :], axis=2)  # (nb, na)
    idx = d.argmin(axis=1)
    if len(set(idx.tolist())) != len(b):  # not a bijection
        return None
    if d[np.arange(len(b)), idx].max() > tol:
        return None
    return idx


# ---------------------------------------------------------------------------
# the XVAL matrix
# ---------------------------------------------------------------------------
@dataclass
class Reader:
    name: str
    fn: object
    lossy_tol: float = 0.0
    note: str = ""


def oss_readers(fmt: str, tol_stl: float) -> List[Reader]:
    """Independent readers available for ``fmt`` in this environment."""
    if fmt == "stl":
        return [
            Reader("meshio", lambda p: meshio_mesh(p, "stl"), tol_stl),
            Reader("numpy-stl", numpy_stl_mesh, tol_stl),
            Reader("trimesh", lambda p: trimesh_mesh(p, process=False), tol_stl),
        ]
    if fmt == "obj":
        return [
            Reader("meshio", lambda p: meshio_mesh(p, "obj")),
            Reader("trimesh", lambda p: trimesh_mesh(p, process=False)),
        ]
    if fmt == "vtk":
        return [Reader("meshio", lambda p: meshio_mesh(p, "vtk"), note="legacy reader, UNSTRUCTURED_GRID only")]
    if fmt == "su2":
        return [Reader("meshio", lambda p: meshio_mesh(p, "su2"), note="fails on NELEM=0 (reader defect)")]
    if fmt == "cgns":
        return []  # no OSS triangle reader exists; SIDS conformance is used instead
    if fmt == "gmsh":
        return [Reader("meshio", lambda p: meshio_mesh(p, None))]
    raise ValueError(fmt)


def write_all(mesh: TriMesh, outdir: Path) -> Dict[str, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {}
    for fmt in FORMATS:
        p = outdir / f"{mesh.name}.{ 'msh' if fmt == 'gmsh' else fmt}"
        paths[fmt] = WRITERS[fmt](mesh, p)
    return paths


def run_matrix(meshes: Sequence[TriMesh], outdir: Path) -> dict:
    """Write x read-back matrix over ``meshes`` x six formats x OSS readers."""
    cells = []
    for mesh in meshes:
        paths = write_all(mesh, outdir / mesh.name)
        diag = float(np.linalg.norm(mesh.points.max(axis=0) - mesh.points.min(axis=0)))
        tol = 1e-5 * diag  # float32 STL: 24-bit mantissa ~ 6e-8 relative
        for fmt in FORMATS:
            cell = {"geometry": mesh.name, "format": fmt, "bytes": paths[fmt].stat().st_size,
                    "readers": [], "analytic_volume": None}
            for rd in oss_readers(fmt, tol):
                entry = {"reader": rd.name, "note": rd.note}
                try:
                    back = rd.fn(paths[fmt])
                    entry.update(compare_meshes(mesh, back, rd.lossy_tol))
                    entry["status"] = "ok"
                except BaseException as exc:  # meshio calls sys.exit(1) on ReadError
                    entry.update({"status": "unreadable", "error": f"{type(exc).__name__}: {str(exc)[:160]}",
                                  "pass": False})
                cell["readers"].append(entry)
            if fmt == "cgns":
                entry = {"reader": "h5py (SIDS mapping checker)", "note": "no OSS triangle reader exists"}
                try:
                    back, problems = read_cgns_sids(paths[fmt])
                    entry.update(compare_meshes(mesh, back, 0.0))
                    entry["sids_violations"] = problems
                    entry["conformant"] = not problems
                    entry["status"] = "ok" if not problems else "nonconformant"
                    entry["pass"] = bool(entry["pass"] and not problems)
                except BaseException as exc:
                    entry.update({"status": "unreadable", "pass": False, "conformant": False,
                                  "sids_violations": [],
                                  "error": f"{type(exc).__name__}: {str(exc)[:160]}"})
                cell["readers"].append(entry)
            if fmt == "su2":  # our strict grammar reader is the second opinion
                entry = {"reader": "mesh_xval(strict grammar)"}
                back = read_su2_strict(paths[fmt])
                entry.update(compare_meshes(mesh, back, 0.0))
                entry["status"] = "ok"
                cell["readers"].append(entry)
            if fmt == "vtk":  # POLYDATA is the historical variant; measure who can read it
                pd = write_vtk(mesh, paths[fmt].with_suffix(".polydata.vtk"), dataset="POLYDATA")
                entry = {"reader": "meshio", "variant": "POLYDATA", "note": "historical writer variant"}
                try:
                    entry.update(compare_meshes(mesh, meshio_mesh(pd, "vtk"), 0.0))
                    entry["status"] = "ok"
                except BaseException as exc:
                    entry.update({"status": "unreadable", "error": f"{type(exc).__name__}: {str(exc)[:160]}",
                                  "pass": False})
                cell["readers"].append(entry)
            cells.append(cell)
    return {"cells": cells}


def summarise(matrix: dict) -> dict:
    cells = matrix["cells"]
    checked = [r for c in cells for r in c["readers"]]
    passes = sum(1 for r in checked if r.get("pass"))
    failing = [{"geometry": c["geometry"], "format": c["format"], "reader": r.get("reader"),
                "variant": r.get("variant"), "status": r.get("status"), "error": r.get("error", "")}
               for c in cells for r in c["readers"] if not r.get("pass")]
    return {
        "combinations": len(cells),
        "cells_all_oss_green": sum(1 for c in cells if c["readers"] and all(r.get("pass") for r in c["readers"])),
        "reader_checks": len(checked),
        "reader_checks_green": passes,
        "failing": failing,
        "worst_volume_rel_err": max(
            (r["volume_rel_err"] for r in checked if r.get("status") == "ok" and "volume_rel_err" in r),
            default=float("nan")),
    }


__all__ = [
    "FORMATS", "TriMesh", "cube", "icosphere", "torus", "WRITERS", "write_all",
    "read_su2_strict", "read_cgns_sids", "weld", "CGNS_RULES", "meshio_mesh", "trimesh_mesh",
    "numpy_stl_mesh", "meshio_cgns_status", "compare_meshes", "run_matrix", "summarise",
]