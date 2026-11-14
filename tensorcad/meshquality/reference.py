"""Reference implementations and external oracles for the mesh-quality metrics.

The point of this module is that :mod:`tensorcad.meshquality.metrics` is never its own
judge. Four independent references are provided:

``reference_triangles`` / ``reference_tets`` / ``reference_hexahedra`` /
``reference_prisms`` / ``reference_pyramids``
    Second implementations of the same quantities by *different arithmetic*: Heron's
    formula instead of a cross product, ``atan2`` instead of ``acos``, Cayley-Menger
    instead of a triple product, and central finite differences of the isoparametric
    map instead of its analytic derivative. They are not a proof, but they catch the
    class of bug this library can have (a wrong normalisation, a transposed Jacobian,
    a corner table that is off by one), because those do not survive two derivations.

``icn_from_singular_values``
    The same inverse condition number written through the singular values of the
    Jacobian, ``d (prod sigma)^(2/d) / sum sigma^2``, instead of through its
    determinant and Frobenius norm. Requires numpy; returns ``None`` without it.

External tools, when present, are used directly:

* **Gmsh** -- ``external_availability()["gmsh"]``; :func:`gmsh_element_qualities` calls
  ``gmsh.model.mesh.getElementQualities(tags, "minSICN")``, which is the number the
  ``gmsh_icn`` metric predicts.
* **trimesh** -- a live third-party implementation of triangle angles, area and aspect
  ratio; :func:`trimesh_triangle_metrics` wraps it.
* **meshmagick** -- its ``print_quality()`` delegates to ``vtk.vtkMeshQuality``, so the
  Verdict definitions mirrored in the metric module are its definitions. A runtime
  adapter is not shipped because meshmagick is GPL-3.0 and its quality path needs VTK;
  the benchmark records availability rather than silently skipping.

Nothing here is imported by the package core at import time.
"""
from __future__ import annotations

import importlib.util
import math
from typing import Any, Mapping, Sequence

__all__ = [
    "external_availability",
    "gmsh_element_qualities",
    "icn_from_singular_values",
    "reference_hexahedra",
    "reference_non_orthogonality",
    "reference_prisms",
    "reference_pyramids",
    "reference_skewness",
    "reference_tets",
    "reference_triangles",
    "trimesh_topology",
    "trimesh_triangle_metrics",
]

Vec3 = tuple[float, float, float]


def _sub(a: Sequence[float], b: Sequence[float]) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _norm(a: Vec3) -> float:
    return math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def _linspace_six(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3) -> tuple[float, ...]:
    return (
        _norm(_sub(p1, p0)),
        _norm(_sub(p2, p1)),
        _norm(_sub(p3, p2)),
        _norm(_sub(p0, p3)),
        _norm(_sub(p0, p2)),
        _norm(_sub(p1, p3)),
    )


# --------------------------------------------------------------------------------------
# independent second implementations
# --------------------------------------------------------------------------------------
def reference_triangles(
    points: Sequence[Sequence[float]], triangles: Sequence[Sequence[int]]
) -> dict[str, tuple[float, ...]]:
    """Triangle metrics via Heron and ``atan2`` -- deliberately unlike the library path."""
    out: dict[str, list[float]] = {
        "area": [],
        "min_angle": [],
        "max_angle": [],
        "aspect_ratio": [],
        "gamma": [],
    }
    for face in triangles:
        p = [tuple(float(c) for c in points[i]) for i in face]
        sides = (_norm(_sub(p[1], p[0])), _norm(_sub(p[2], p[1])), _norm(_sub(p[0], p[2])))
        semiperimeter = 0.5 * sum(sides)
        # Heron with the numerically stable product form.
        area = math.sqrt(
            max(
                0.0,
                semiperimeter
                * (semiperimeter - sides[0])
                * (semiperimeter - sides[1])
                * (semiperimeter - sides[2]),
            )
        )
        angles = []
        for k in range(3):
            a = _sub(p[(k + 1) % 3], p[k])
            b = _sub(p[(k + 2) % 3], p[k])
            angles.append(math.degrees(math.atan2(_norm(_cross(tuple(a), tuple(b))),
                                                  a[0] * b[0] + a[1] * b[1] + a[2] * b[2])))
        # R = a / (2 sin A) with ``a`` the side opposite the angle -- a different path to
        # the circumradius than the library's a*b*c/(4A).
        circumradius = (
            sides[1] / (2.0 * math.sin(math.radians(angles[0]))) if area > 0.0 else math.inf
        )
        inradius = area / semiperimeter if semiperimeter > 0.0 else 0.0
        out["area"].append(area)
        out["min_angle"].append(min(angles))
        out["max_angle"].append(max(angles))
        out["aspect_ratio"].append(circumradius / (2.0 * inradius) if inradius > 0.0 else math.inf)
        out["gamma"].append(2.0 * inradius / circumradius if circumradius not in (0.0, math.inf) else 0.0)
    return {name: tuple(values) for name, values in out.items()}


def reference_skewness(
    points: Sequence[Sequence[float]], triangles: Sequence[Sequence[int]]
) -> tuple[float, ...]:
    """Equiangle skewness of triangles via the ``atan2`` angles of :func:`reference_triangles`.

    Same definition, different angles: no code path is shared with
    :mod:`tensorcad.meshquality.metrics`, so a wrong angle set or a swapped min/max in the
    library shows up as a mismatch rather than as a silently equal number.
    """
    angles = reference_triangles(points, triangles)
    out: list[float] = []
    for low, high in zip(angles["min_angle"], angles["max_angle"]):
        out.append(max((high - 60.0) / 120.0, (60.0 - low) / 60.0))
    return tuple(out)


def _cayley_menger(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3) -> float:
    """Tetrahedron volume from the 5x5 Cayley-Menger determinant (a different algebra)."""
    d = [[0.0] * 5 for _ in range(5)]
    pts = (p0, p1, p2, p3)
    for i in range(4):
        for j in range(4):
            v = _sub(pts[i], pts[j])
            d[i][j] = v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
    for i in range(4):
        d[i][4] = d[4][i] = 1.0
    # Bareiss-free Gaussian elimination is enough at this size and keeps the code short.
    matrix = [row[:] for row in d]
    sign = 1.0
    for column in range(5):
        pivot = max(range(column, 5), key=lambda r: abs(matrix[r][column]))
        if matrix[pivot][column] == 0.0:
            return 0.0
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
            sign = -sign
        for row in range(column + 1, 5):
            factor = matrix[row][column] / matrix[column][column]
            for k in range(column, 5):
                matrix[row][k] -= factor * matrix[column][k]
    det = sign
    for i in range(5):
        det *= matrix[i][i]
    return math.sqrt(max(0.0, det / 288.0))


def reference_tets(
    points: Sequence[Sequence[float]], tets: Sequence[Sequence[int]]
) -> dict[str, tuple[float, ...]]:
    """Tetrahedron volume (Cayley-Menger), eta (edge sum) and dihedral angles (face normals)."""
    out: dict[str, list[float]] = {"volume": [], "eta": [], "min_dihedral": [], "max_dihedral": []}
    for tet in tets:
        p = [tuple(float(c) for c in points[i]) for i in tet]
        volume = _cayley_menger(p[0], p[1], p[2], p[3])
        lengths = [
            _norm(_sub(p[j], p[i]))
            for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
        ]
        total = sum(length * length for length in lengths)
        out["volume"].append(volume)
        out["eta"].append(12.0 * (3.0 * volume) ** (2.0 / 3.0) / total if total > 0.0 else 0.0)
        # Dihedral angles from the inward normals of the two faces sharing an edge.
        faces = ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1))
        normals = []
        for face in faces:
            a, b, c = (p[i] for i in face)
            normals.append(_cross(_sub(b, a), _sub(c, a)))
        interior = []
        for edge in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
            adjacent = [k for k, face in enumerate(faces) if edge[0] in face and edge[1] in face]
            if len(adjacent) != 2:
                continue
            n0, n1 = normals[adjacent[0]], normals[adjacent[1]]
            scale = _norm(n0) * _norm(n1)
            if scale == 0.0:  # a collapsed cell has no dihedral angle to measure
                continue
            cosine = (n0[0] * n1[0] + n0[1] * n1[1] + n0[2] * n1[2]) / scale
            interior.append(180.0 - math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
        interior = [a for a in interior if not math.isnan(a)]
        out["min_dihedral"].append(min(interior) if interior else math.nan)
        out["max_dihedral"].append(max(interior) if interior else math.nan)
    return {name: tuple(values) for name, values in out.items()}


#: Corner signs for the reference cube, repeated here so the reference code does not
#: import the tables it is meant to check.
_CORNER_SIGNS = (
    (-1.0, -1.0, -1.0),
    (1.0, -1.0, -1.0),
    (1.0, 1.0, -1.0),
    (-1.0, 1.0, -1.0),
    (-1.0, -1.0, 1.0),
    (1.0, -1.0, 1.0),
    (1.0, 1.0, 1.0),
    (-1.0, 1.0, 1.0),
)


def _trilinear(points: Sequence[Vec3], xi: float, eta: float, zeta: float) -> Vec3:
    total = [0.0, 0.0, 0.0]
    for sign, p in zip(_CORNER_SIGNS, points):
        weight = 0.125
        weight *= 1.0 + sign[0] * xi
        weight *= 1.0 + sign[1] * eta
        weight *= 1.0 + sign[2] * zeta
        for k in range(3):
            total[k] += weight * p[k]
    return (total[0], total[1], total[2])


def reference_hexahedra(
    points: Sequence[Sequence[float]], hexes: Sequence[Sequence[int]], *, step: float = 1e-5
) -> dict[str, tuple[float, ...]]:
    """Corners Jacobian by central finite differences of the trilinear map.

    Independent of the analytic shape-function derivatives in the library: it evaluates
    the map and differentiates numerically, so an error in the derivative table shows up
    as a finite difference. Accuracy is ~``step**2`` (1e-10 at the default step), which
    is why the comparison tolerance is 1e-7 and not 1e-12.
    """
    out: dict[str, list[float]] = {"scaled_jacobian": [], "min_sicn": []}
    for hex_nodes in hexes:
        p = [tuple(float(c) for c in points[i]) for i in hex_nodes]
        worst_sj = math.inf
        worst_icn = math.inf
        for xi, eta, zeta in _CORNER_SIGNS:
            columns = []
            for axis in range(3):
                plus = [xi, eta, zeta]
                minus = [xi, eta, zeta]
                plus[axis] += step
                minus[axis] -= step
                a = _trilinear(p, plus[0], plus[1], plus[2])
                b = _trilinear(p, minus[0], minus[1], minus[2])
                columns.append(
                    ((a[0] - b[0]) / (2 * step), (a[1] - b[1]) / (2 * step), (a[2] - b[2]) / (2 * step))
                )
            # J[row][column] with rows x, y, z.
            jac = [[columns[column][row] for column in range(3)] for row in range(3)]
            det = (
                jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1])
                - jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0])
                + jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0])
            )
            norms = [math.sqrt(sum(jac[r][c] ** 2 for r in range(3))) for c in range(3)]
            fro2 = sum(jac[r][c] ** 2 for r in range(3) for c in range(3))
            sj = det / (norms[0] * norms[1] * norms[2]) if all(n > 0.0 for n in norms) else 0.0
            icn = 3.0 * abs(det) ** (2.0 / 3.0) / fro2 if fro2 > 0.0 else 0.0
            worst_sj = min(worst_sj, sj)
            worst_icn = min(worst_icn, icn if det >= 0.0 else -icn)
        out["scaled_jacobian"].append(worst_sj)
        out["min_sicn"].append(worst_icn)
    return {name: tuple(values) for name, values in out.items()}


def icn_from_singular_values(jacobian: Sequence[Sequence[float]]) -> float | None:
    """``d (prod sigma)^(2/d) / sum sigma^2`` -- the same ICN through the singular values."""
    try:
        import numpy as np
    except ImportError:  # pragma: no cover - numpy is a test/dev dependency
        return None
    matrix = np.asarray(jacobian, dtype=float)
    singular = np.linalg.svd(matrix, compute_uv=False)
    if not np.any(singular):
        return 0.0
    dim = matrix.shape[0]
    return float(dim * float(math.prod(singular)) ** (2.0 / dim) / float((singular**2).sum()))


# --------------------------------------------------------------------------------------
# external oracles
# --------------------------------------------------------------------------------------
def _version(module: str) -> str | None:
    spec = importlib.util.find_spec(module)
    if spec is None:
        return None
    try:
        loaded = importlib.import_module(module)
    except Exception:  # pragma: no cover - a broken install is reported, not hidden
        return "import-failed"
    return str(getattr(loaded, "__version__", "unknown"))


def external_availability() -> dict[str, str | None]:
    """Version of each external oracle, or ``None`` when it is not importable."""
    return {
        "gmsh": _version("gmsh"),
        "meshmagick": _version("meshmagick"),
        "trimesh": _version("trimesh"),
    }


def trimesh_triangle_metrics(
    points: Sequence[Sequence[float]], triangles: Sequence[Sequence[int]]
) -> dict[str, tuple[float, ...]] | None:
    """Triangle area and angles from ``trimesh`` (third party).

    trimesh 5.0 does not export a Verdict aspect ratio (``trimesh.triangles`` has
    ``area`` and ``angles`` but no ``aspect_ratio``), so the comparison covers the two
    primitives and the library's aspect ratio is checked against the analytic values and
    :func:`reference_triangles` instead of against trimesh.
    """
    try:
        import numpy as np
        import trimesh.triangles as tri
    except ImportError:  # pragma: no cover - exercised only where trimesh is missing
        return None
    vertices = np.asarray(points, dtype=float)
    faces = np.asarray(triangles, dtype=int)
    tris = vertices[faces]
    angles = np.degrees(tri.angles(tris))
    return {
        "area": tuple(float(v) for v in tri.area(tris)),
        "min_angle": tuple(float(v) for v in angles.min(axis=1)),
        "max_angle": tuple(float(v) for v in angles.max(axis=1)),
    }


#: Gmsh element type numbers, for the adapter below.
_GMSH_TYPES = {"tri": 2, "quad": 3, "tet": 4, "hex": 5}


def gmsh_element_qualities(
    points: Sequence[Sequence[float]],
    elements: Mapping[str, Sequence[Sequence[int]]],
    measure: str = "minSICN",
) -> dict[str, tuple[float, ...]] | None:
    """Run the real Gmsh on the mesh and return its own quality numbers.

    Returns ``None`` when Gmsh is not importable (it ships no wheel for this
    architecture), so callers can record ``not_run`` instead of silently passing.
    The tags are 1-based and match the input order.
    """
    try:
        import gmsh
    except ImportError:  # pragma: no cover - depends on the host
        return None
    points_per_kind = {"tri": 3, "quad": 4, "tet": 4, "hex": 8}
    result: dict[str, tuple[float, ...]] = {}
    gmsh.initialize()
    try:
        gmsh.model.add("tensorcad_meshquality")
        tags = [gmsh.model.addDiscreteEntity(0, i + 1, []) for i in range(len(points))]
        for tag, xyz in zip(tags, points):
            gmsh.model.mesh.setCoordinates([(0, tag)], [list(xyz)])
        offset = 0
        for kind, cells in elements.items():
            if kind not in _GMSH_TYPES:
                continue
            flat: list[int] = []
            entities: list[int] = []
            for cell in cells:
                entities.append(offset)
                offset += 1
                flat.extend(tags[i] for i in cell)
            gmsh.model.mesh.addElements(0, _GMSH_TYPES[kind], entities, [flat])
            qualities = gmsh.model.mesh.getElementQualities(entities, measure)
            result[kind] = tuple(float(v) for v in qualities)
    finally:
        gmsh.finalize()
    return result


def trimesh_topology(
    points: Sequence[Sequence[float]], triangles: Sequence[Sequence[int]]
) -> dict[str, Any] | None:
    """Closure, winding, Euler number, shell count and volume from ``trimesh``.

    A live third-party implementation of exactly the checks
    :mod:`tensorcad.meshquality.topology` makes, so the library is not its own judge:
    ``is_watertight`` for closure, ``is_winding_consistent`` for orientation,
    ``euler_number`` for the topology, ``body_count`` for shells and ``volume`` for the
    enclosed volume. Returns ``None`` when trimesh is not importable, which the benchmark
    records instead of silently skipping. ``process=False`` is deliberate: trimesh's
    merge/repair passes would hide the defects being compared.
    """
    try:
        import trimesh  # noqa: PLC0415 -- optional dependency, imported on demand
    except Exception:  # noqa: BLE001 - any import failure means "no oracle here"
        return None
    mesh = trimesh.Trimesh(
        vertices=[[float(c) for c in p] for p in points],
        faces=[[int(i) for i in face] for face in triangles],
        process=False,
        validate=False,
    )
    return {
        "version": getattr(trimesh, "__version__", None),
        "watertight": bool(mesh.is_watertight),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "euler_number": int(mesh.euler_number),
        "body_count": int(mesh.body_count),
        "volume": float(mesh.volume),
        "area": float(mesh.area),
    }


def reference_non_orthogonality(
    points: Sequence[Sequence[float]], cells: Sequence[Sequence[int]], kind: str
) -> tuple[float, ...]:
    """Internal-face non-orthogonality by a different arithmetic route than the library.

    :func:`tensorcad.meshquality.topology.non_orthogonality` computes the angle with
    ``acos`` of a normalised dot product and orients each face normal by comparing it with
    the centre-to-centre vector. This route instead builds the normal by summing the
    triangles of a fan around the face centroid and takes the angle with ``atan2(|u x v|,
    u . v)``, which stays accurate for nearly parallel vectors where ``acos`` loses half
    its digits. Both must agree to about 1e-9 deg on a well-conditioned face.
    """
    from .topology import CELL_FACES  # local import: reference must not import at module load

    table = CELL_FACES[kind]
    pts = [tuple(float(c) for c in p) for p in points]
    owners: dict[tuple[int, ...], list[int]] = {}
    order: dict[tuple[int, ...], tuple[int, ...]] = {}
    for index, cell in enumerate(cells):
        for nodes in table:
            face = tuple(int(cell[i]) for i in nodes)
            key = tuple(sorted(face))
            owners.setdefault(key, []).append(index)
            order.setdefault(key, face)
    centres = [
        tuple(sum(pts[i][axis] for i in cell) / len(cell) for axis in range(3)) for cell in cells
    ]

    def fan_normal(face: tuple[int, ...]) -> tuple[float, float, float]:
        centre = tuple(sum(pts[i][axis] for i in face) / len(face) for axis in range(3))
        total = [0.0, 0.0, 0.0]
        for position in range(len(face)):
            a = pts[face[position]]
            b = pts[face[(position + 1) % len(face)]]
            for axis in range(3):
                total[axis] += (a[axis] - centre[axis] + b[axis] - centre[axis]) * 0.0
            # triangle (centre, a, b) contributes cross(a-centre, b-centre)
            u = [a[i] - centre[i] for i in range(3)]
            v = [b[i] - centre[i] for i in range(3)]
            total[0] += u[1] * v[2] - u[2] * v[1]
            total[1] += u[2] * v[0] - u[0] * v[2]
            total[2] += u[0] * v[1] - u[1] * v[0]
        return (total[0], total[1], total[2])

    angles: list[float] = []
    for key, sharing in owners.items():
        if len(sharing) != 2:
            continue
        low, high = sorted(sharing)
        normal = fan_normal(order[key])
        d = [centres[high][i] - centres[low][i] for i in range(3)]
        if sum(normal[i] * d[i] for i in range(3)) < 0.0:
            normal = [-normal[i] for i in range(3)]
        cross = [
            normal[1] * d[2] - normal[2] * d[1],
            normal[2] * d[0] - normal[0] * d[2],
            normal[0] * d[1] - normal[1] * d[0],
        ]
        magnitude = math.sqrt(sum(c * c for c in cross))
        dot = sum(normal[i] * d[i] for i in range(3))
        angles.append(math.degrees(math.atan2(magnitude, dot)))
    return tuple(angles)


# --------------------------------------------------------------------------------------
# prisms and pyramids: second implementations
# --------------------------------------------------------------------------------------
# Same quantities as :func:`tensorcad.meshquality.metrics.prism_metrics` /
# ``pyramid_metrics``, by *different arithmetic*: the Jacobian is a central finite
# difference of the isoparametric map instead of its analytic derivative, volumes come
# from Cayley-Menger determinants of a conforming tet split instead of triple products,
# and angles come from ``atan2(|u x w|, u.w)`` instead of ``acos`` of the normalised dot
# product. What survives both derivations is a property of the shape; what does not is a
# bug in one of them.


def _prism_map(p: Sequence[Vec3], xi: float, eta: float, zeta: float) -> Vec3:
    """Physical point of the prism map at ``(xi, eta, zeta)``."""
    one_m_z = 1.0 - zeta
    w = (
        (1.0 - xi - eta) * one_m_z,
        xi * one_m_z,
        eta * one_m_z,
        (1.0 - xi - eta) * zeta,
        xi * zeta,
        eta * zeta,
    )
    return (
        sum(weight * point[0] for weight, point in zip(w, p)),
        sum(weight * point[1] for weight, point in zip(w, p)),
        sum(weight * point[2] for weight, point in zip(w, p)),
    )


def _pyramid_map(p: Sequence[Vec3], xi: float, eta: float, zeta: float) -> Vec3:
    """Physical point of the pyramid map at ``(xi, eta, zeta)``."""
    one_m_z = 1.0 - zeta
    w = (
        one_m_z * (1.0 - xi) * (1.0 - eta),
        one_m_z * xi * (1.0 - eta),
        one_m_z * xi * eta,
        one_m_z * (1.0 - xi) * eta,
        zeta,
    )
    return (
        sum(weight * point[0] for weight, point in zip(w, p)),
        sum(weight * point[1] for weight, point in zip(w, p)),
        sum(weight * point[2] for weight, point in zip(w, p)),
    )


def _fd_jacobian(shape, p: Sequence[Vec3], xi: float, eta: float, zeta: float, h: float = 1.0e-5) -> list[list[float]]:
    """Central finite differences of the isoparametric map in ``(xi, eta, zeta)``."""
    jac = [[0.0, 0.0, 0.0] for _ in range(3)]
    for column, (d_xi, d_eta, d_zeta) in enumerate(
        ((h, 0.0, 0.0), (0.0, h, 0.0), (0.0, 0.0, h))
    ):
        plus = shape(p, xi + d_xi, eta + d_eta, zeta + d_zeta)
        minus = shape(p, xi - d_xi, eta - d_eta, zeta - d_zeta)
        for row in range(3):
            jac[row][column] = (plus[row] - minus[row]) / (2.0 * h)
    return jac


def _icn_fd(jac: Sequence[Sequence[float]]) -> float:
    """``3 |det|^(2/3) / ||J||_F^2`` of a finite-difference Jacobian, signed by ``det``."""
    det = (
        jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1])
        - jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0])
        + jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0])
    )
    fro2 = sum(value * value for row in jac for value in row)
    if fro2 == 0.0:
        return 0.0
    value = 3.0 * abs(det) ** (2.0 / 3.0) / fro2
    return value if det >= 0.0 else -value


def _sj_fd(jac: Sequence[Sequence[float]]) -> float:
    """``det J / prod|column|`` of a finite-difference Jacobian."""
    det = (
        jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1])
        - jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0])
        + jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0])
    )
    cols = [math.sqrt(sum(jac[row][column] ** 2 for row in range(3))) for column in range(3)]
    if not all(c > 0.0 for c in cols):
        return 0.0
    return det / (cols[0] * cols[1] * cols[2])


def _atan2_angle(u: Vec3, w: Vec3) -> float:
    """Angle between two vectors in degrees via ``atan2`` (not ``acos``)."""
    cross = _cross(u, w)
    return math.degrees(math.atan2(_norm(cross), u[0] * w[0] + u[1] * w[1] + u[2] * w[2]))


def _atan2_triangle_angles(a: Vec3, b: Vec3, c: Vec3) -> tuple[float, float, float]:
    return (
        _atan2_angle(_sub(b, a), _sub(c, a)),
        _atan2_angle(_sub(a, b), _sub(c, b)),
        _atan2_angle(_sub(a, c), _sub(b, c)),
    )


def _atan2_quad_angles(p: Sequence[Vec3]) -> tuple[float, float, float, float]:
    return tuple(  # type: ignore[return-value]
        _atan2_angle(_sub(p[(i - 1) % 4], p[i]), _sub(p[(i + 1) % 4], p[i])) for i in range(4)
    )


def _equiangle(angles: Sequence[float], ideal: float) -> float:
    valid = [a for a in angles if not math.isnan(a)]
    if not valid:
        return math.nan
    low, high = min(valid), max(valid)
    return max((high - ideal) / (180.0 - ideal), (ideal - low) / ideal)


def _relative_icn(shape, reference: Sequence[Vec3], p: Sequence[Vec3], point: Sequence[float]) -> float:
    """Finite-difference ICN of the physical cell divided by the reference element's."""
    ideal = _icn_fd(_fd_jacobian(shape, reference, *point))
    return _icn_fd(_fd_jacobian(shape, p, *point)) / ideal


#: Sample points, edges and face tables, duplicated here on purpose: the reference route
#: must not import the tables it is meant to check.
_REF_PRISM_CORNERS = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0),
                      (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (0.0, 1.0, 1.0))
_REF_PRISM_EDGES = ((0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5))
_REF_PRISM_TRIS = ((0, 1, 2), (3, 4, 5))
_REF_PRISM_QUADS = ((0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5))
_REF_PRISM_TETS = ((0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5))
_REF_PYRAMID_POINTS = ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
                       (0.5, 0.5, 0.0))
_REF_PYRAMID_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4))
_REF_PYRAMID_BASE = (0, 1, 2, 3)
_REF_PYRAMID_SIDES = ((0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4))
_REF_PYRAMID_TETS = ((0, 1, 2, 4), (0, 2, 3, 4))


def _ideal_prism_points() -> tuple[Vec3, ...]:
    h = math.sqrt(3.0) / 2.0
    return ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, h, 0.0),
            (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (0.5, h, 1.0))


def _ideal_pyramid_points() -> tuple[Vec3, ...]:
    return ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0), (0.5, 0.5, 1.0))


def _cayley_menger_split(p: Sequence[Vec3], tets: Sequence[Sequence[int]]) -> float:
    total = 0.0
    for a, b, c, d in tets:
        total += _cayley_menger(p[a], p[b], p[c], p[d])
    return abs(total)


def reference_prisms(
    points: Sequence[Vec3], prisms: Sequence[Sequence[int]]
) -> dict[str, tuple[float, ...]]:
    """Second route to the prism metric set (FD Jacobian, Cayley-Menger, atan2)."""
    ideal = _ideal_prism_points()
    out: dict[str, list[float]] = {
        "volume": [],
        "min_edge": [],
        "max_edge": [],
        "edge_ratio": [],
        "min_angle": [],
        "max_angle": [],
        "skewness": [],
        "scaled_jacobian": [],
        "min_sicn": [],
    }
    for prism in prisms:
        p = [points[i] for i in prism]
        out["volume"].append(_cayley_menger_split(p, _REF_PRISM_TETS))
        lengths = [_norm(_sub(p[j], p[i])) for i, j in _REF_PRISM_EDGES]
        out["min_edge"].append(min(lengths))
        out["max_edge"].append(max(lengths))
        out["edge_ratio"].append(max(lengths) / min(lengths))
        tri = [a for face in _REF_PRISM_TRIS for a in _atan2_triangle_angles(*[p[i] for i in face])]
        quad = [a for face in _REF_PRISM_QUADS for a in _atan2_quad_angles([p[i] for i in face])]
        out["min_angle"].append(min(tri + quad))
        out["max_angle"].append(max(tri + quad))
        out["skewness"].append(max(_equiangle(tri, 60.0), _equiangle(quad, 90.0)))
        sjs: list[float] = []
        sicns: list[float] = []
        for point in _REF_PRISM_CORNERS:
            jac = _fd_jacobian(_prism_map, p, *point)
            sjs.append(_sj_fd(jac))
            sicns.append(_relative_icn(_prism_map, ideal, p, point))
        out["scaled_jacobian"].append(min(sjs))
        out["min_sicn"].append(min(sicns))
    return {name: tuple(values) for name, values in out.items()}


def reference_pyramids(
    points: Sequence[Vec3], pyramids: Sequence[Sequence[int]]
) -> dict[str, tuple[float, ...]]:
    """Second route to the pyramid metric set (the apex is excluded, as in the library)."""
    ideal = _ideal_pyramid_points()
    out: dict[str, list[float]] = {
        "volume": [],
        "min_edge": [],
        "max_edge": [],
        "edge_ratio": [],
        "min_angle": [],
        "max_angle": [],
        "skewness": [],
        "scaled_jacobian": [],
        "min_sicn": [],
    }
    for pyramid in pyramids:
        p = [points[i] for i in pyramid]
        out["volume"].append(_cayley_menger_split(p, _REF_PYRAMID_TETS))
        lengths = [_norm(_sub(p[j], p[i])) for i, j in _REF_PYRAMID_EDGES]
        out["min_edge"].append(min(lengths))
        out["max_edge"].append(max(lengths))
        out["edge_ratio"].append(max(lengths) / min(lengths))
        base = list(_atan2_quad_angles([p[i] for i in _REF_PYRAMID_BASE]))
        side = [a for face in _REF_PYRAMID_SIDES for a in _atan2_triangle_angles(*[p[i] for i in face])]
        out["min_angle"].append(min(base + side))
        out["max_angle"].append(max(base + side))
        out["skewness"].append(max(_equiangle(base, 90.0), _equiangle(side, 60.0)))
        sjs: list[float] = []
        sicns: list[float] = []
        for point in _REF_PYRAMID_POINTS:
            jac = _fd_jacobian(_pyramid_map, p, *point)
            sjs.append(_sj_fd(jac))
            sicns.append(_relative_icn(_pyramid_map, ideal, p, point))
        out["scaled_jacobian"].append(min(sjs))
        out["min_sicn"].append(min(sicns))
    return {name: tuple(values) for name, values in out.items()}