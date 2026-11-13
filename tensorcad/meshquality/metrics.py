"""Element quality metrics for the four cell shapes a solver actually sees.

Every metric here is a number, not a verdict, and every one carries an explicit
normalisation so that two implementations can be compared without guessing what the
"1.0" means. The definitions come from two published sets:

* **Gmsh** (``src/mesh/qualityMeasures.cpp``, ``src/mesh/qualityMeasuresJacobian.cpp``,
  ``live-clones/gmsh`` mirror, blob sha ``dcc7972`` / ``ab2b34f``). ``gamma``, ``eta``
  and ``gmsh_icn`` are independent transcriptions of Gmsh's own expressions; the
  ``gmsh_icn`` unit is Gmsh's reference element (the right-corner tetrahedron), not the
  regular one. No Gmsh code was copied -- the formulas are arithmetic and the
  provenance is recorded so a reader can check them line by line.
* **Verdict / VTK** (``vtkMeshQuality``: aspect ratio, edge ratio, scaled Jacobian,
  min/max angle). These are the metrics ``meshmagick.print_quality()`` reports,
  because meshmagick delegates to ``vtk.vtkMeshQuality``.

Two conventions are worth stating before the code, because most disputes about mesh
quality are really disputes about these two sentences:

``scaled_jacobian``
    shape *only*: the element Jacobian with its column norms divided out. A cell that
    is thin but orthogonal scores 1.0, which is what a boundary layer wants. A regular
    tetrahedron scores ``1/sqrt(2)``; a cube, a square and an equilateral triangle
    score 1.0.

``min_sicn``
    shape *and* anisotropy: ``d * |det J|**(2/d) / ||J||_F**2`` sampled at the element
    corners, normalised against the *ideal* element of that type (equilateral
    triangle, square, **regular** tetrahedron, cube), so an ideal element scores 1.0.
    Negative means inverted. Gmsh's own ``minSICN`` normalises against its reference
    tetrahedron (a right-corner one); :func:`tet_metrics` therefore also exposes
    ``gmsh_icn``, the same quantity in Gmsh's unit, so a comparison against Gmsh is an
    identity rather than a calibration.

Only the standard library is used: the package core must import with no extras.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

__all__ = [
    "ELEMENT_KINDS",
    "METRIC_NAMES",
    "ElementQuality",
    "hex_metrics",
    "quad_metrics",
    "tet_metrics",
    "tri_metrics",
]

Vec3 = tuple[float, float, float]

#: Cell shapes with a metric set. ``tri``/``quad`` are surface cells, ``tet``/``hex``
#: volume cells.
ELEMENT_KINDS = ("tri", "quad", "tet", "hex")

#: Metric name -> what one value means. Kept as data so a gate can be written without
#: importing the implementation, and so the documentation cannot drift from the code.
METRIC_NAMES: dict[str, str] = {
    "area": "element area (m^2, absolute)",
    "volume": "element volume (m^3, absolute; the sign is reported by min_sicn)",
    "min_edge": "shortest edge (m)",
    "max_edge": "longest edge (m)",
    "edge_ratio": "max_edge / min_edge (Verdict edge ratio, >= 1)",
    "min_angle": "smallest interior/face angle (deg)",
    "max_angle": "largest interior/face angle (deg)",
    "min_dihedral": "smallest interior dihedral angle (deg)",
    "max_dihedral": "largest interior dihedral angle (deg)",
    "aspect_ratio": "triangle circumradius / (2 * inradius) (Verdict; 1 = equilateral)",
    "gamma": "Gmsh gamma (2 r/R for triangles, 3 r/R for tetrahedra; 1 = ideal)",
    "eta": "Gmsh eta (1 - |60 - a|/60 for triangles, 12(3V)^(2/3)/sum l^2 for "
           "tetrahedra, 1 - max|90 - a|/90 for quadrilaterals; 1 = ideal)",
    "scaled_jacobian": "min corner det(J)/prod|column| (shape only; 1 = orthogonal, "
                       "regular tetrahedron = 1/sqrt(2), <= 0 = inverted/flat)",
    "min_sicn": "min corner d |det J|^(2/d) / ||J||_F^2 normalised to the ideal element "
                "(<= 0 = inverted/flat)",
    "gmsh_icn": "Gmsh minSICN unit: same expression normalised to Gmsh's reference "
                "tetrahedron (1 = right-corner tetrahedron, regular tetrahedron = 0.7937)",
}

_NAN = math.nan


# --------------------------------------------------------------------------------------
# small vector helpers (tuples, no numpy: the core has no dependencies)
# --------------------------------------------------------------------------------------
def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: Vec3) -> float:
    return math.sqrt(_dot(a, a))


def _unit(a: Vec3) -> Vec3:
    n = _norm(a)
    if n == 0.0:
        return (0.0, 0.0, 0.0)
    return (a[0] / n, a[1] / n, a[2] / n)


def _scale(a: Vec3, s: float) -> Vec3:
    return (a[0] * s, a[1] * s, a[2] * s)


def _angle_deg(u: Vec3, v: Vec3) -> float:
    """Angle between two vectors in degrees, clamped so round-off cannot raise."""
    nu, nv = _norm(u), _norm(v)
    if nu == 0.0 or nv == 0.0:
        return _NAN
    c = _dot(u, v) / (nu * nv)
    return math.degrees(math.acos(max(-1.0, min(1.0, c))))


# --------------------------------------------------------------------------------------
# container
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class ElementQuality:
    """Per-element metric values of one cell kind, in input order."""

    kind: str
    metrics: dict[str, tuple[float, ...]]

    def __getitem__(self, metric: str) -> tuple[float, ...]:
        return self.metrics[metric]

    def __len__(self) -> int:
        return len(next(iter(self.metrics.values()))) if self.metrics else 0


# --------------------------------------------------------------------------------------
# triangles
# --------------------------------------------------------------------------------------
def _triangle_angles(p0: Vec3, p1: Vec3, p2: Vec3) -> tuple[float, float, float]:
    return (
        _angle_deg(_sub(p1, p0), _sub(p2, p0)),
        _angle_deg(_sub(p2, p1), _sub(p0, p1)),
        _angle_deg(_sub(p0, p2), _sub(p1, p2)),
    )


def _triangle_gamma(area: float, sides: tuple[float, float, float]) -> float:
    """``2 r / R`` -- algebraically the same as Gmsh's sine expression, 1 for equilateral."""
    a, b, c = sides
    if area <= 0.0 or a <= 0.0 or b <= 0.0 or c <= 0.0:
        return 0.0
    circumradius = a * b * c / (4.0 * area)
    inradius = area / (0.5 * (a + b + c))
    return 2.0 * inradius / circumradius


def _triangle_sicn(p0: Vec3, p1: Vec3, p2: Vec3, signed_area: float) -> float:
    """``2 det(M)/||M||_F^2`` with ``M`` measured against the equilateral triangle.

    ``M = J_phys . J_ideal^-1``; writing ``p2`` in the frame ``x = p1 - p0``,
    ``z = (p1-p0) x (p2-p0)`` gives the closed form below, so no 2x2 solve is needed.
    """
    base = _sub(p1, p0)
    length = _norm(base)
    normal = _cross(base, _sub(p2, p0))
    if length == 0.0 or _norm(normal) == 0.0:
        return 0.0
    x_hat = _scale(base, 1.0 / length)
    z_hat = _unit(normal)
    y_hat = _cross(z_hat, x_hat)
    v = _sub(p2, p0)
    x2, y2 = _dot(v, x_hat), _dot(v, y_hat)
    det_m = 4.0 * signed_area / math.sqrt(3.0)
    # ||M||_F^2 of [[L, (-L + 2 x2)/sqrt(3)], [0, 2 y2/sqrt(3)]]
    row1 = (-length + 2.0 * x2) / math.sqrt(3.0)
    row2 = 2.0 * y2 / math.sqrt(3.0)
    fro2 = length * length + row1 * row1 + row2 * row2
    if fro2 == 0.0:
        return 0.0
    return 2.0 * det_m / fro2


def tri_metrics(points: Sequence[Vec3], triangles: Sequence[Sequence[int]]) -> ElementQuality:
    """Area, angles, Verdict aspect ratio and Gmsh gamma/eta/SICN of triangles."""
    area: list[float] = []
    min_edge: list[float] = []
    max_edge: list[float] = []
    edge_ratio: list[float] = []
    min_angle: list[float] = []
    max_angle: list[float] = []
    aspect: list[float] = []
    gamma: list[float] = []
    eta: list[float] = []
    sj: list[float] = []
    sicn: list[float] = []
    for face in triangles:
        p0, p1, p2 = (points[i] for i in face)
        cr = _cross(_sub(p1, p0), _sub(p2, p0))
        magnitude = _norm(cr)
        signed_area = 0.5 * _dot(cr, _unit(cr)) if magnitude else 0.0
        abs_area = 0.5 * magnitude
        sides = (_norm(_sub(p1, p0)), _norm(_sub(p2, p1)), _norm(_sub(p0, p2)))
        shortest, longest = min(sides), max(sides)
        angles = _triangle_angles(p0, p1, p2)
        g = _triangle_gamma(abs_area, sides)
        area.append(abs_area)
        min_edge.append(shortest)
        max_edge.append(longest)
        edge_ratio.append(longest / shortest if shortest > 0.0 else math.inf)
        min_angle.append(min(angles))
        max_angle.append(max(angles))
        aspect.append(1.0 / g if g > 0.0 else math.inf)
        gamma.append(g)
        eta.append(1.0 - abs(60.0 - min(angles)) / 60.0)
        # Gmsh NCJ == Verdict scaled Jacobian for triangles: (2/sqrt3) sin(angle),
        # signed by the winding so a flipped triangle is negative.
        sign = 1.0 if signed_area > 0.0 else -1.0
        sj.append(sign * (2.0 / math.sqrt(3.0)) * math.sin(math.radians(min(angles))))
        sicn.append(_triangle_sicn(p0, p1, p2, signed_area))
    return ElementQuality(
        "tri",
        {
            "area": tuple(area),
            "min_edge": tuple(min_edge),
            "max_edge": tuple(max_edge),
            "edge_ratio": tuple(edge_ratio),
            "min_angle": tuple(min_angle),
            "max_angle": tuple(max_angle),
            "aspect_ratio": tuple(aspect),
            "gamma": tuple(gamma),
            "eta": tuple(eta),
            "scaled_jacobian": tuple(sj),
            "min_sicn": tuple(sicn),
        },
    )


# --------------------------------------------------------------------------------------
# quadrilaterals and hexahedra: bilinear / trilinear corner Jacobians
# --------------------------------------------------------------------------------------
_QUAD_CORNERS = ((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0))
#: VTK/Gmsh hexahedron corner signs, in the node order (0-3 bottom, 4-7 top).
_HEX_CORNERS = (
    (-1.0, -1.0, -1.0),
    (1.0, -1.0, -1.0),
    (1.0, 1.0, -1.0),
    (-1.0, 1.0, -1.0),
    (-1.0, -1.0, 1.0),
    (1.0, -1.0, 1.0),
    (1.0, 1.0, 1.0),
    (-1.0, 1.0, 1.0),
)
_HEX_FACES = ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))


def _quad_jacobian(p: Sequence[Vec3], xi: float, eta: float) -> tuple[float, float, float, float]:
    """2x2 Jacobian of the bilinear quad at ``(xi, eta)`` as ``(dx/dxi, dy/dxi, dx/deta, dy/deta)``."""
    signs = ((-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0))
    dx_dxi = dy_dxi = dx_deta = dy_deta = 0.0
    for k, (sx, sy) in enumerate(signs):
        x, y = p[k][0], p[k][1]
        dx_dxi += 0.25 * sx * (1.0 + sy * eta) * x
        dy_dxi += 0.25 * sx * (1.0 + sy * eta) * y
        dx_deta += 0.25 * sy * (1.0 + sx * xi) * x
        dy_deta += 0.25 * sy * (1.0 + sx * xi) * y
    return dx_dxi, dy_dxi, dx_deta, dy_deta


def _to_plane(points: Sequence[Vec3], normal: Vec3) -> list[Vec3]:
    """Project element points onto their own best-fit plane (2D coords, z = 0)."""
    n = _unit(normal)
    if n == (0.0, 0.0, 0.0):
        return [(p[0], p[1], 0.0) for p in points]
    helper: Vec3 = (1.0, 0.0, 0.0) if abs(n[0]) < 0.9 else (0.0, 1.0, 0.0)
    x_hat = _unit(_cross(helper, n))
    y_hat = _cross(n, x_hat)
    origin = points[0]
    return [(_dot(_sub(p, origin), x_hat), _dot(_sub(p, origin), y_hat), 0.0) for p in points]


def _quad_area(p: Sequence[Vec3]) -> float:
    total = 0.0
    for i in range(4):
        total += _norm(_cross(p[i], p[(i + 1) % 4])) * 0.5
    return total


def _quad_angles(p: Sequence[Vec3]) -> tuple[float, float, float, float]:
    return tuple(  # type: ignore[return-value]
        _angle_deg(_sub(p[(i - 1) % 4], p[i]), _sub(p[(i + 1) % 4], p[i])) for i in range(4)
    )


def quad_metrics(points: Sequence[Vec3], quads: Sequence[Sequence[int]]) -> ElementQuality:
    """Metric set of bilinear quadrilaterals (reference square; an ideal square scores 1)."""
    area: list[float] = []
    min_edge: list[float] = []
    max_edge: list[float] = []
    edge_ratio: list[float] = []
    min_angle: list[float] = []
    max_angle: list[float] = []
    eta: list[float] = []
    sj: list[float] = []
    sicn: list[float] = []
    for quad in quads:
        raw = [points[i] for i in quad]
        normal = (0.0, 0.0, 0.0)
        for i in range(4):  # Newell normal, robust for non-planar quads
            a, b = raw[i], raw[(i + 1) % 4]
            normal = _add(
                normal,
                (
                    (a[1] - b[1]) * (a[2] + b[2]),
                    (a[2] - b[2]) * (a[0] + b[0]),
                    (a[0] - b[0]) * (a[1] + b[1]),
                ),
            )
        p = _to_plane(raw, normal)
        sides = [_norm(_sub(p[i], p[(i + 1) % 4])) for i in range(4)]
        area.append(_quad_area(p))
        min_edge.append(min(sides))
        max_edge.append(max(sides))
        edge_ratio.append(max(sides) / min(sides) if min(sides) > 0.0 else math.inf)
        angles = _quad_angles(p)
        min_angle.append(min(angles))
        max_angle.append(max(angles))
        eta.append(1.0 - max(abs(90.0 - a) for a in angles) / 90.0)
        worst_sj = math.inf
        worst_sicn = math.inf
        for xi, eta_c in _QUAD_CORNERS:
            dx_dxi, dy_dxi, dx_deta, dy_deta = _quad_jacobian(p, xi, eta_c)
            det = dx_dxi * dy_deta - dx_deta * dy_dxi
            fro2 = dx_dxi * dx_dxi + dy_dxi * dy_dxi + dx_deta * dx_deta + dy_deta * dy_deta
            col0 = math.hypot(dx_dxi, dy_dxi)
            col1 = math.hypot(dx_deta, dy_deta)
            sj_corner = det / (col0 * col1) if col0 > 0.0 and col1 > 0.0 else 0.0
            sicn_corner = 2.0 * det / fro2 if fro2 > 0.0 else 0.0
            worst_sj = min(worst_sj, sj_corner)
            worst_sicn = min(worst_sicn, sicn_corner)
        sj.append(worst_sj)
        sicn.append(worst_sicn)
    return ElementQuality(
        "quad",
        {
            "area": tuple(area),
            "min_edge": tuple(min_edge),
            "max_edge": tuple(max_edge),
            "edge_ratio": tuple(edge_ratio),
            "min_angle": tuple(min_angle),
            "max_angle": tuple(max_angle),
            "eta": tuple(eta),
            "scaled_jacobian": tuple(sj),
            "min_sicn": tuple(sicn),
        },
    )


def _hex_jacobian(points: Sequence[Vec3], xi: float, eta: float, zeta: float) -> list[list[float]]:
    """3x3 Jacobian of the trilinear hexahedron at ``(xi, eta, zeta)`` (rows d/dxi, d/deta, d/dzeta)."""
    jac = [[0.0, 0.0, 0.0] for _ in range(3)]
    for k, (sx, sy, sz) in enumerate(_HEX_CORNERS):
        x, y, z = points[k]
        d_xi = 0.125 * sx * (1.0 + sy * eta) * (1.0 + sz * zeta)
        d_eta = 0.125 * sy * (1.0 + sx * xi) * (1.0 + sz * zeta)
        d_zeta = 0.125 * sz * (1.0 + sx * xi) * (1.0 + sy * eta)
        for row, value in ((0, x), (1, y), (2, z)):
            jac[row][0] += value * d_xi
            jac[row][1] += value * d_eta
            jac[row][2] += value * d_zeta
    return jac


def _det3(jac: Sequence[Sequence[float]]) -> float:
    return (
        jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1])
        - jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0])
        + jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0])
    )


def _frobenius2(jac: Sequence[Sequence[float]]) -> float:
    return sum(value * value for row in jac for value in row)


def _corner_quality(jac: Sequence[Sequence[float]]) -> tuple[float, float]:
    """``(scaled_jacobian, min_sicn_normalised_to_the_ideal_element)`` at one corner."""
    det = _det3(jac)
    fro2 = _frobenius2(jac)
    cols = [math.sqrt(sum(jac[row][column] ** 2 for row in range(3))) for column in range(3)]
    sj = det / (cols[0] * cols[1] * cols[2]) if all(c > 0.0 for c in cols) else 0.0
    if fro2 == 0.0:
        return sj, 0.0
    # 3 |det|^(2/3) / ||J||_F^2, signed so an inverted corner is negative. The reference
    # cube is similar to the unit cube, so this is already the normalised value.
    magnitude = abs(det) ** (2.0 / 3.0)
    sicn = 3.0 * magnitude / fro2
    return sj, sicn if det >= 0.0 else -sicn


def _hex_volume(points: Sequence[Vec3]) -> float:
    """Volume by the divergence theorem over the six (possibly warped) faces."""
    total = 0.0
    for face in _HEX_FACES:
        a, b, c, d = (points[i] for i in face)
        for tri in ((a, b, c), (a, c, d)):
            total += _dot(tri[0], _cross(tri[1], tri[2])) / 6.0
    return abs(total)


_HEX_EDGES = (
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7),
)


def hex_metrics(points: Sequence[Vec3], hexes: Sequence[Sequence[int]]) -> ElementQuality:
    """Metric set of trilinear hexahedra (reference cube; an ideal cube scores 1)."""
    volume: list[float] = []
    min_edge: list[float] = []
    max_edge: list[float] = []
    edge_ratio: list[float] = []
    min_angle: list[float] = []
    max_angle: list[float] = []
    sj: list[float] = []
    sicn: list[float] = []
    for hex_nodes in hexes:
        p = [points[i] for i in hex_nodes]
        lengths = [_norm(_sub(p[b], p[a])) for a, b in _HEX_EDGES]
        volume.append(_hex_volume(p))
        min_edge.append(min(lengths))
        max_edge.append(max(lengths))
        edge_ratio.append(max(lengths) / min(lengths) if min(lengths) > 0.0 else math.inf)
        face_angles: list[float] = []
        for face in _HEX_FACES:
            face_angles.extend(_quad_angles([p[i] for i in face]))
        valid = [a for a in face_angles if not math.isnan(a)]
        min_angle.append(min(valid) if valid else _NAN)
        max_angle.append(max(valid) if valid else _NAN)
        worst_sj = math.inf
        worst_sicn = math.inf
        for xi, eta, zeta in _HEX_CORNERS:
            sj_corner, sicn_corner = _corner_quality(_hex_jacobian(p, xi, eta, zeta))
            worst_sj = min(worst_sj, sj_corner)
            worst_sicn = min(worst_sicn, sicn_corner)
        sj.append(worst_sj)
        sicn.append(worst_sicn)
    return ElementQuality(
        "hex",
        {
            "volume": tuple(volume),
            "min_edge": tuple(min_edge),
            "max_edge": tuple(max_edge),
            "edge_ratio": tuple(edge_ratio),
            "min_angle": tuple(min_angle),
            "max_angle": tuple(max_angle),
            "scaled_jacobian": tuple(sj),
            "min_sicn": tuple(sicn),
        },
    )


# --------------------------------------------------------------------------------------
# tetrahedra
# --------------------------------------------------------------------------------------
#: Corner -> the other three vertices in an order that keeps det(J) positive.
_TET_CORNERS = ((0, (1, 2, 3)), (1, (0, 3, 2)), (2, (0, 1, 3)), (3, (0, 2, 1)))
#: Edge pairs (i, j) whose squared lengths multiply in Gmsh's circumradius expression.
_TET_OPPOSITE = ((1, 0, 3, 2), (2, 0, 3, 1), (3, 0, 2, 1))


def _tet_signed_volume(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3) -> float:
    return _dot(_sub(p1, p0), _cross(_sub(p2, p0), _sub(p3, p0))) / 6.0


def _inv3(m: Sequence[Sequence[float]]) -> list[list[float]]:
    det = _det3(m)
    if det == 0.0:
        return [[0.0, 0.0, 0.0] for _ in range(3)]
    return [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det,
        ],
    ]


def _regular_tet() -> tuple[Vec3, ...]:
    """Unit-edge regular tetrahedron, vertices 0..3 (same ordering as the cell)."""
    h = math.sqrt(3.0) / 2.0
    z = math.sqrt(2.0 / 3.0)
    return ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, h, 0.0), (0.5, h / 3.0, z))


def _gmsh_tet() -> tuple[Vec3, ...]:
    """Gmsh's reference tetrahedron: (0,0,0), (1,0,0), (0,1,0), (0,0,1)."""
    return ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))


def _tet_corner_icn(
    p: Sequence[Vec3], reference: Sequence[Vec3], corner: int
) -> float:
    """``3 sign(det M) |det M|^(2/3) / ||M||_F^2`` with ``M`` the physical edge matrix."""
    _, others = _TET_CORNERS[corner]
    edges = [_sub(p[j], p[corner]) for j in others]
    ideal = [_sub(reference[j], reference[corner]) for j in others]
    inv = _inv3([[ideal[c][r] for c in range(3)] for r in range(3)])
    e = [[edges[c][r] for c in range(3)] for r in range(3)]
    m = [[sum(e[r][k] * inv[k][c] for k in range(3)) for c in range(3)] for r in range(3)]
    det = _det3(m)
    fro2 = _frobenius2(m)
    if fro2 == 0.0:
        return 0.0
    magnitude = abs(det) ** (2.0 / 3.0)
    value = 3.0 * magnitude / fro2
    return value if det >= 0.0 else -value


def _tet_dihedral_angles(p: Sequence[Vec3]) -> list[float]:
    """Interior dihedral angles (deg) at the six edges, via the vectors to the opposite vertices."""
    edges = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    out = []
    for i, j in edges:
        opposite = [k for k in range(4) if k not in (i, j)]
        axis = _unit(_sub(p[j], p[i]))
        if axis == (0.0, 0.0, 0.0):
            out.append(_NAN)
            continue
        vecs = []
        for k in opposite:
            v = _sub(p[k], p[i])
            vecs.append(_sub(v, _scale(axis, _dot(v, axis))))
        out.append(_angle_deg(vecs[0], vecs[1]))
    return out


def tet_metrics(points: Sequence[Vec3], tets: Sequence[Sequence[int]]) -> ElementQuality:
    """Volume, dihedral angles, Gmsh gamma/eta, scaled Jacobian and both SICN units."""
    regular = _regular_tet()
    gmsh_ref = _gmsh_tet()
    volume: list[float] = []
    min_edge: list[float] = []
    max_edge: list[float] = []
    edge_ratio: list[float] = []
    min_dihedral: list[float] = []
    max_dihedral: list[float] = []
    gamma: list[float] = []
    eta: list[float] = []
    sj: list[float] = []
    sicn: list[float] = []
    gmsh_icn: list[float] = []
    for tet in tets:
        p = [points[i] for i in tet]
        signed = _tet_signed_volume(*p)  # type: ignore[misc]
        abs_volume = abs(signed)
        lengths = [
            _norm(_sub(p[j], p[i]))
            for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
        ]
        volume.append(abs_volume)
        min_edge.append(min(lengths))
        max_edge.append(max(lengths))
        edge_ratio.append(max(lengths) / min(lengths) if min(lengths) > 0.0 else math.inf)
        dihedrals = _tet_dihedral_angles(p)
        valid = [a for a in dihedrals if not math.isnan(a)]
        min_dihedral.append(min(valid) if valid else _NAN)
        max_dihedral.append(max(valid) if valid else _NAN)

        # Gmsh qmTetrahedron::gamma: rho = 3 * inradius, gamma = rho / R.
        areas = [
            _norm(_cross(_sub(p[1], p[0]), _sub(p[2], p[0]))) * 0.5,
            _norm(_cross(_sub(p[2], p[0]), _sub(p[3], p[0]))) * 0.5,
            _norm(_cross(_sub(p[1], p[0]), _sub(p[3], p[0]))) * 0.5,
            _norm(_cross(_sub(p[2], p[1]), _sub(p[3], p[1]))) * 0.5,
        ]
        if abs_volume == 0.0:
            gamma.append(0.0)
        else:
            a, b, c = _dot(_sub(p[1], p[0]), _sub(p[1], p[0])), _dot(
                _sub(p[2], p[0]), _sub(p[2], p[0])
            ), _dot(_sub(p[3], p[0]), _sub(p[3], p[0]))
            big_a = _dot(_sub(p[3], p[2]), _sub(p[3], p[2]))
            big_b = _dot(_sub(p[3], p[1]), _sub(p[3], p[1]))
            big_c = _dot(_sub(p[2], p[1]), _sub(p[2], p[1]))
            prod_a = math.sqrt(a * big_a)
            prod_b = math.sqrt(b * big_b)
            prod_c = math.sqrt(c * big_c)
            inside = (prod_a + prod_b + prod_c) * (prod_a + prod_b - prod_c) * (
                prod_a - prod_b + prod_c
            ) * (-prod_a + prod_b + prod_c)
            if inside <= 0.0:
                gamma.append(0.0)
            else:
                part_r = math.sqrt(inside) / 24.0
                rho = 9.0 * abs_volume / sum(areas)
                gamma.append(rho * abs_volume / part_r)
        # Gmsh qmTetrahedron::eta.
        sum_l2 = sum(length * length for length in lengths)
        eta.append(12.0 * (3.0 * abs_volume) ** (2.0 / 3.0) / sum_l2 if sum_l2 > 0.0 else 0.0)

        worst_sj = math.inf
        for corner in range(4):
            _, others = _TET_CORNERS[corner]
            edges = [_unit(_sub(p[j], p[corner])) for j in others]
            sj_corner = _dot(edges[0], _cross(edges[1], edges[2]))
            worst_sj = min(worst_sj, sj_corner)
        sj.append(worst_sj)
        sicn.append(min(_tet_corner_icn(p, regular, corner) for corner in range(4)))
        gmsh_icn.append(min(_tet_corner_icn(p, gmsh_ref, corner) for corner in range(4)))
    return ElementQuality(
        "tet",
        {
            "volume": tuple(volume),
            "min_edge": tuple(min_edge),
            "max_edge": tuple(max_edge),
            "edge_ratio": tuple(edge_ratio),
            "min_dihedral": tuple(min_dihedral),
            "max_dihedral": tuple(max_dihedral),
            "gamma": tuple(gamma),
            "eta": tuple(eta),
            "scaled_jacobian": tuple(sj),
            "min_sicn": tuple(sicn),
            "gmsh_icn": tuple(gmsh_icn),
        },
    )