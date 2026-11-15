"""Vectorised (numpy) second implementation of the per-element metrics.

Why this module exists
----------------------
The scalar implementation in :mod:`tensorcad.meshquality.metrics` costs about 260
microseconds per hexahedron, because every corner Jacobian, face angle and edge length is
a Python-level tuple computation. A CFD mesh of 10^6 hexahedra -- the size the solver
consumers of this repository actually hand over -- therefore needs roughly four minutes
just to be measured, which is not a gate, it is a coffee break. The acceptance of
``GAP-CAD-D1-MESH-QUALITY-LIB`` asks for a library usable *as a CI gate*, so the cost has
to come down by the factor that separates "per cell in Python" from "per cell in a vector
instruction".

This module is not a rewrite of the definitions: every expression below is the same
arithmetic as the scalar path, evaluated on arrays of cells instead of one cell at a time,
using the same corner/face/edge tables (imported from :mod:`.metrics`, so the two
implementations cannot drift apart). It is therefore both

* the fast path (``backend="auto"`` picks it above :data:`AUTO_MIN_CELLS`, ``backend="numpy"``
  forces it), and
* a **second implementation** of the same contract, which is what lets the test suite check
  the scalar reference against something other than itself.

How close "the same" is, is measured rather than asserted: see
``docs/mesh_quality_numpy_backend.md``. The residual is summation order inside dot
products; the guards are bit-identical -- an ideal cell still reads exactly 1.0 for
``scaled_jacobian``, a zero-measure cell still reads ``0.0`` with ``edge_ratio == inf``, a
flipped tetrahedron still reads a negative ``min_sicn`` and a degenerate edge still reads
``NaN`` for the angles.

Deliberate limits
-----------------
* ``prism`` and ``pyramid`` are **not** vectorised: they sample a finite-difference
  Jacobian at seven/nine points and are rare in the meshes this library gates. They keep
  the scalar path, and asking for ``backend="numpy"`` on them raises instead of silently
  falling back.
* numpy is imported lazily by :mod:`.metrics`; this module imports it at module level, so
  the dependency-free core never touches it unless a caller asks for it.
"""
from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from .metrics import (
    _HEX_CORNERS,
    _HEX_EDGES,
    _HEX_FACES,
    _QUAD_CORNERS,
    _TET_CORNERS,
    ElementQuality,
)

UNSUPPORTED_KINDS = ("prism", "pyramid")

__all__ = [
    "AUTO_MIN_CELLS",
    "UNSUPPORTED_KINDS",
    "hex_metrics",
    "quad_metrics",
    "tet_metrics",
    "tri_metrics",
]

#: Cell count at or above which ``backend="auto"`` prefers this module.
AUTO_MIN_CELLS = 128

_NAN = float("nan")
_INF = float("inf")


# --------------------------------------------------------------------------------------
# vector primitives (arrays of shape (n, ..., 3))
# --------------------------------------------------------------------------------------
def _norm(v: np.ndarray) -> np.ndarray:
    return np.sqrt((v * v).sum(axis=-1))


def _unit(v: np.ndarray) -> np.ndarray:
    n = _norm(v)
    safe = np.where(n > 0.0, n, 1.0)
    out = v / safe[..., None]
    return np.where((n > 0.0)[..., None], out, 0.0)


def _angle_deg(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    nu = _norm(u)
    nv = _norm(v)
    ok = (nu > 0.0) & (nv > 0.0)
    denom = np.where(ok, nu * nv, 1.0)
    c = np.clip((u * v).sum(axis=-1) / denom, -1.0, 1.0)
    return np.where(ok, np.degrees(np.arccos(c)), _NAN)


def _det3(jac: np.ndarray) -> np.ndarray:
    return (
        jac[..., 0, 0] * (jac[..., 1, 1] * jac[..., 2, 2] - jac[..., 1, 2] * jac[..., 2, 1])
        - jac[..., 0, 1] * (jac[..., 1, 0] * jac[..., 2, 2] - jac[..., 1, 2] * jac[..., 2, 0])
        + jac[..., 0, 2] * (jac[..., 1, 0] * jac[..., 2, 1] - jac[..., 1, 1] * jac[..., 2, 0])
    )


def _frobenius2(jac: np.ndarray) -> np.ndarray:
    return (jac * jac).sum(axis=(-1, -2))


def _angle_set_stats(angles: np.ndarray, ideal: float) -> tuple[np.ndarray, ...]:
    """``(min, max, equiangle_skew)`` over the last axis, NaN-aware.

    Matches the scalar helpers: an angle of ``NaN`` (degenerate edge) is dropped from the
    minimum and the maximum, and a set with no valid angle gives ``NaN`` throughout.
    """
    valid = ~np.isnan(angles)
    count = valid.sum(axis=-1)
    low = np.where(valid, angles, _INF).min(axis=-1)
    high = np.where(valid, angles, -_INF).max(axis=-1)
    none = count == 0
    low = np.where(none, _NAN, low)
    high = np.where(none, _NAN, high)
    skew = np.maximum((high - ideal) / (180.0 - ideal), (ideal - low) / ideal)
    return low, high, np.where(none, _NAN, skew)


def _points(
    points: Sequence[Sequence[float]], cells: Sequence[Sequence[int]], width: int
) -> tuple[np.ndarray, np.ndarray]:
    p = np.asarray(points, dtype=float)
    if p.ndim != 2 or p.shape[1] != 3:
        raise ValueError("points must be an (n, 3) array")
    c = np.asarray(cells, dtype=np.int64)
    if c.size == 0:
        return p, c.reshape(0, width)
    if c.ndim != 2 or c.shape[1] != width:
        raise ValueError(f"cells must be an (m, {width}) array")
    if c.min() < 0 or c.max() >= p.shape[0]:
        raise IndexError("cell connectivity references a point that does not exist")
    return p, c


def _empty(kind: str, names: Sequence[str]) -> ElementQuality:
    return ElementQuality(kind, {name: () for name in names})


def _result(kind: str, values: dict[str, np.ndarray]) -> ElementQuality:
    return ElementQuality(
        kind, {name: tuple(float(v) for v in series) for name, series in values.items()}
    )


def _ratio(longest: np.ndarray, shortest: np.ndarray) -> np.ndarray:
    return np.where(
        shortest > 0.0, longest / np.where(shortest > 0.0, shortest, 1.0), _INF
    )


# --------------------------------------------------------------------------------------
# triangles
# --------------------------------------------------------------------------------------
_TRI_KEYS = (
    "area",
    "min_edge",
    "max_edge",
    "edge_ratio",
    "min_angle",
    "max_angle",
    "skewness",
    "aspect_ratio",
    "gamma",
    "eta",
    "scaled_jacobian",
    "min_sicn",
)


def tri_metrics(
    points: Sequence[Sequence[float]], triangles: Sequence[Sequence[int]]
) -> ElementQuality:
    """Vectorised :func:`tensorcad.meshquality.metrics.tri_metrics`."""
    p_all, cells = _points(points, triangles, 3)
    if cells.shape[0] == 0:
        return _empty("tri", _TRI_KEYS)
    q = p_all[cells]
    p0, p1, p2 = q[:, 0], q[:, 1], q[:, 2]
    base = p1 - p0
    cr = np.cross(base, p2 - p0)
    magnitude = _norm(cr)
    # The scalar path defines signed_area as 0.5 * dot(cr, unit(cr)) == 0.5 * |cr|, so it is
    # never negative: a triangle floating in 3-D carries no orientation without a reference
    # normal. Reproduced here, including that consequence (see docs section 5).
    signed_area = 0.5 * magnitude
    abs_area = 0.5 * magnitude

    sides = np.stack([_norm(base), _norm(p2 - p1), _norm(p0 - p2)], axis=-1)
    shortest = sides.min(axis=-1)
    longest = sides.max(axis=-1)
    angles = np.stack(
        [
            _angle_deg(p1 - p0, p2 - p0),
            _angle_deg(p2 - p1, p0 - p1),
            _angle_deg(p0 - p2, p1 - p2),
        ],
        axis=-1,
    )
    min_angle, max_angle, skew = _angle_set_stats(angles, 60.0)

    a, b, c = sides[:, 0], sides[:, 1], sides[:, 2]
    ok = (abs_area > 0.0) & (a > 0.0) & (b > 0.0) & (c > 0.0)
    safe_area = np.where(abs_area > 0.0, abs_area, 1.0)
    perimeter = np.where(a + b + c > 0.0, a + b + c, 1.0)
    circumradius = np.where(ok, a * b * c / (4.0 * safe_area), 0.0)
    inradius = np.where(ok, abs_area / (0.5 * perimeter), 0.0)
    gamma = np.where(
        ok, 2.0 * inradius / np.where(circumradius > 0.0, circumradius, 1.0), 0.0
    )
    aspect = np.where(gamma > 0.0, 1.0 / np.where(gamma > 0.0, gamma, 1.0), _INF)

    sign = np.where(signed_area > 0.0, 1.0, -1.0)
    sj = sign * (2.0 / math.sqrt(3.0)) * np.sin(np.radians(min_angle))

    # _triangle_sicn, verbatim.
    length = _norm(base)
    normal = np.cross(base, p2 - p0)
    normal_len = _norm(normal)
    live = (length > 0.0) & (normal_len > 0.0)
    x_hat = _unit(base)
    z_hat = _unit(normal)
    y_hat = np.cross(z_hat, x_hat)
    v = p2 - p0
    x2 = (v * x_hat).sum(axis=-1)
    y2 = (v * y_hat).sum(axis=-1)
    det_m = 4.0 * signed_area / math.sqrt(3.0)
    row1 = (-length + 2.0 * x2) / math.sqrt(3.0)
    row2 = 2.0 * y2 / math.sqrt(3.0)
    fro2 = length * length + row1 * row1 + row2 * row2
    sicn = np.where(live & (fro2 > 0.0), 2.0 * det_m / np.where(fro2 > 0.0, fro2, 1.0), 0.0)

    return _result(
        "tri",
        {
            "area": abs_area,
            "min_edge": shortest,
            "max_edge": longest,
            "edge_ratio": _ratio(longest, shortest),
            "min_angle": min_angle,
            "max_angle": max_angle,
            "skewness": skew,
            "aspect_ratio": aspect,
            "gamma": gamma,
            "eta": 1.0 - np.abs(60.0 - min_angle) / 60.0,
            "scaled_jacobian": sj,
            "min_sicn": sicn,
        },
    )


# --------------------------------------------------------------------------------------
# quadrilaterals
# --------------------------------------------------------------------------------------
_QUAD_KEYS = (
    "area",
    "min_edge",
    "max_edge",
    "edge_ratio",
    "min_angle",
    "max_angle",
    "skewness",
    "eta",
    "scaled_jacobian",
    "min_sicn",
)


def _to_plane(raw: np.ndarray, normal: np.ndarray) -> np.ndarray:
    """Vectorised ``_to_plane``: project each quad onto its own Newell plane."""
    n_hat = _unit(normal)
    degenerate = (n_hat == 0.0).all(axis=-1)
    helper = np.where(
        np.abs(n_hat[:, 0:1]) < 0.9,
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
    )
    x_hat = _unit(np.cross(helper, n_hat))
    y_hat = np.cross(n_hat, x_hat)
    rel = raw - raw[:, 0:1, :]
    px = (rel * x_hat[:, None, :]).sum(axis=-1)
    py = (rel * y_hat[:, None, :]).sum(axis=-1)
    projected = np.stack([px, py, np.zeros_like(px)], axis=-1)
    fallback = np.stack([raw[..., 0], raw[..., 1], np.zeros_like(px)], axis=-1)
    return np.where(degenerate[:, None, None], fallback, projected)


def _quad_angles(p: np.ndarray) -> np.ndarray:
    return np.stack(
        [_angle_deg(p[:, (i - 1) % 4] - p[:, i], p[:, (i + 1) % 4] - p[:, i]) for i in range(4)],
        axis=-1,
    )


def _quad_jacobian_corners(
    p: np.ndarray, xi: np.ndarray, eta: np.ndarray
) -> tuple[np.ndarray, ...]:
    sx = np.array([c[0] for c in _QUAD_CORNERS])
    sy = np.array([c[1] for c in _QUAD_CORNERS])
    dx_dxi = np.zeros(p.shape[0])
    dy_dxi = np.zeros(p.shape[0])
    dx_deta = np.zeros(p.shape[0])
    dy_deta = np.zeros(p.shape[0])
    for k in range(4):
        wx = 0.25 * sx[k] * (1.0 + sy[k] * eta)
        wy = 0.25 * sy[k] * (1.0 + sx[k] * xi)
        dx_dxi = dx_dxi + wx * p[:, k, 0]
        dy_dxi = dy_dxi + wx * p[:, k, 1]
        dx_deta = dx_deta + wy * p[:, k, 0]
        dy_deta = dy_deta + wy * p[:, k, 1]
    return dx_dxi, dy_dxi, dx_deta, dy_deta


def quad_metrics(
    points: Sequence[Sequence[float]], quads: Sequence[Sequence[int]]
) -> ElementQuality:
    """Vectorised :func:`tensorcad.meshquality.metrics.quad_metrics`."""
    p_all, cells = _points(points, quads, 4)
    if cells.shape[0] == 0:
        return _empty("quad", _QUAD_KEYS)
    raw = p_all[cells]
    normal = np.zeros((raw.shape[0], 3))
    for i in range(4):
        a = raw[:, i]
        b = raw[:, (i + 1) % 4]
        normal = normal + np.stack(
            [
                (a[:, 1] - b[:, 1]) * (a[:, 2] + b[:, 2]),
                (a[:, 2] - b[:, 2]) * (a[:, 0] + b[:, 0]),
                (a[:, 0] - b[:, 0]) * (a[:, 1] + b[:, 1]),
            ],
            axis=-1,
        )
    p = _to_plane(raw, normal)

    sides = np.stack([_norm(p[:, i] - p[:, (i + 1) % 4]) for i in range(4)], axis=-1)
    shortest = sides.min(axis=-1)
    longest = sides.max(axis=-1)
    area = 0.5 * sum(_norm(np.cross(p[:, i], p[:, (i + 1) % 4])) for i in range(4))
    angles = _quad_angles(p)
    min_angle, max_angle, skew = _angle_set_stats(angles, 90.0)
    eta = 1.0 - np.abs(90.0 - angles).max(axis=-1) / 90.0

    worst_sj = np.full(raw.shape[0], _INF)
    worst_sicn = np.full(raw.shape[0], _INF)
    for xi_corner, eta_corner in _QUAD_CORNERS:
        # NB: not ``eta``/``xi`` -- ``eta`` is the Gmsh metric computed above, and shadowing
        # it here silently returned the corner coordinate as the metric.
        xi = np.full(raw.shape[0], xi_corner)
        eta_corner_values = np.full(raw.shape[0], eta_corner)
        dx_dxi, dy_dxi, dx_deta, dy_deta = _quad_jacobian_corners(p, xi, eta_corner_values)
        det = dx_dxi * dy_deta - dx_deta * dy_dxi
        fro2 = dx_dxi * dx_dxi + dy_dxi * dy_dxi + dx_deta * dx_deta + dy_deta * dy_deta
        col0 = np.hypot(dx_dxi, dy_dxi)
        col1 = np.hypot(dx_deta, dy_deta)
        ok = (col0 > 0.0) & (col1 > 0.0)
        sj_corner = np.where(ok, det / np.where(ok, col0 * col1, 1.0), 0.0)
        sicn_corner = np.where(fro2 > 0.0, 2.0 * det / np.where(fro2 > 0.0, fro2, 1.0), 0.0)
        worst_sj = np.minimum(worst_sj, sj_corner)
        worst_sicn = np.minimum(worst_sicn, sicn_corner)

    return _result(
        "quad",
        {
            "area": area,
            "min_edge": shortest,
            "max_edge": longest,
            "edge_ratio": _ratio(longest, shortest),
            "min_angle": min_angle,
            "max_angle": max_angle,
            "skewness": skew,
            "eta": eta,
            "scaled_jacobian": worst_sj,
            "min_sicn": worst_sicn,
        },
    )


# --------------------------------------------------------------------------------------
# hexahedra
# --------------------------------------------------------------------------------------
_HEX_KEYS = (
    "volume",
    "min_edge",
    "max_edge",
    "edge_ratio",
    "min_angle",
    "max_angle",
    "skewness",
    "scaled_jacobian",
    "min_sicn",
)

_HEX_SIGNS = np.array(_HEX_CORNERS)  # node signs, indexed like the eight nodes
_HEX_TRIANGLES = tuple(
    tri
    for face in _HEX_FACES
    for tri in ((face[0], face[1], face[2]), (face[0], face[2], face[3]))
)


def _corner_quality(
    det: np.ndarray, fro2: np.ndarray, cols: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised ``_corner_quality``: ``(scaled_jacobian, signed min_sicn)`` per corner."""
    ok = (cols > 0.0).all(axis=-1)
    prod = cols[..., 0] * cols[..., 1] * cols[..., 2]
    sj = np.where(ok, det / np.where(ok, prod, 1.0), 0.0)
    live = fro2 != 0.0
    magnitude = np.abs(det) ** (2.0 / 3.0)
    value = np.where(live, 3.0 * magnitude / np.where(live, fro2, 1.0), 0.0)
    return sj, np.where(det >= 0.0, value, -value)


def hex_metrics(
    points: Sequence[Sequence[float]], hexes: Sequence[Sequence[int]]
) -> ElementQuality:
    """Vectorised :func:`tensorcad.meshquality.metrics.hex_metrics`."""
    p_all, cells = _points(points, hexes, 8)
    if cells.shape[0] == 0:
        return _empty("hex", _HEX_KEYS)
    p = p_all[cells]
    n = p.shape[0]

    volume = np.zeros(n)
    for a, b, c in _HEX_TRIANGLES:
        volume = volume + (p[:, a] * np.cross(p[:, b], p[:, c])).sum(axis=-1) / 6.0
    volume = np.abs(volume)

    sides = np.stack([_norm(p[:, j] - p[:, i]) for i, j in _HEX_EDGES], axis=-1)
    shortest = sides.min(axis=-1)
    longest = sides.max(axis=-1)

    face_angles = np.concatenate([_quad_angles(p[:, list(face)]) for face in _HEX_FACES], axis=-1)
    min_angle, max_angle, skew = _angle_set_stats(face_angles, 90.0)

    # Trilinear corner Jacobians: the weight of node k at corner c, broadcast over corners.
    sx, sy, sz = _HEX_SIGNS[:, 0], _HEX_SIGNS[:, 1], _HEX_SIGNS[:, 2]
    signs = np.array(_HEX_CORNERS)
    xi = signs[:, 0][:, None]
    eta = signs[:, 1][:, None]
    zeta = signs[:, 2][:, None]
    w_xi = 0.125 * sx[None, :] * (1.0 + sy[None, :] * eta) * (1.0 + sz[None, :] * zeta)
    w_eta = 0.125 * sy[None, :] * (1.0 + sx[None, :] * xi) * (1.0 + sz[None, :] * zeta)
    w_zeta = 0.125 * sz[None, :] * (1.0 + sx[None, :] * xi) * (1.0 + sy[None, :] * eta)
    jac = np.empty((n, 8, 3, 3))
    for col, w in enumerate((w_xi, w_eta, w_zeta)):
        jac[:, :, :, col] = np.einsum("ck,nkr->ncr", w, p)
    det = _det3(jac)
    fro2 = _frobenius2(jac)
    cols = np.sqrt((jac * jac).sum(axis=-2))
    sj_corner, sicn_corner = _corner_quality(det, fro2, cols)

    return _result(
        "hex",
        {
            "volume": volume,
            "min_edge": shortest,
            "max_edge": longest,
            "edge_ratio": _ratio(longest, shortest),
            "min_angle": min_angle,
            "max_angle": max_angle,
            "skewness": skew,
            "scaled_jacobian": sj_corner.min(axis=-1),
            "min_sicn": sicn_corner.min(axis=-1),
        },
    )


# --------------------------------------------------------------------------------------
# tetrahedra
# --------------------------------------------------------------------------------------
_TET_KEYS = (
    "volume",
    "min_edge",
    "max_edge",
    "edge_ratio",
    "min_dihedral",
    "max_dihedral",
    "skewness",
    "gamma",
    "eta",
    "scaled_jacobian",
    "min_sicn",
    "gmsh_icn",
)

_TET_EDGES = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
_TET_FACES = ((0, 2, 1), (0, 1, 3), (0, 3, 2), (1, 2, 3))


def _regular_tet_points() -> np.ndarray:
    h = math.sqrt(3.0) / 2.0
    z = math.sqrt(2.0 / 3.0)
    return np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, h, 0.0], [0.5, h / 3.0, z]])


def _gmsh_tet_points() -> np.ndarray:
    return np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])


def _inverse_columns(reference: np.ndarray, corner: int) -> np.ndarray:
    """Inverse of the ideal edge matrix at ``corner`` (columns = edges to the other vertices)."""
    _, others = _TET_CORNERS[corner]
    ideal = np.stack([reference[j] - reference[corner] for j in others], axis=1)
    return np.linalg.inv(ideal)


def _tet_corner_icn(p: np.ndarray, reference: np.ndarray, corner: int) -> np.ndarray:
    """Vectorised ``_tet_corner_icn`` for one corner and one reference element."""
    _, others = _TET_CORNERS[corner]
    edges = np.stack([p[:, j] - p[:, corner] for j in others], axis=2)
    m = edges @ _inverse_columns(reference, corner)
    det = _det3(m)
    fro2 = _frobenius2(m)
    live = fro2 != 0.0
    magnitude = np.abs(det) ** (2.0 / 3.0)
    value = np.where(live, 3.0 * magnitude / np.where(live, fro2, 1.0), 0.0)
    return np.where(det >= 0.0, value, -value)


def tet_metrics(
    points: Sequence[Sequence[float]], tets: Sequence[Sequence[int]]
) -> ElementQuality:
    """Vectorised :func:`tensorcad.meshquality.metrics.tet_metrics`."""
    p_all, cells = _points(points, tets, 4)
    if cells.shape[0] == 0:
        return _empty("tet", _TET_KEYS)
    p = p_all[cells]
    n = p.shape[0]

    signed_volume = (p[:, 1] - p[:, 0]) * np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0])
    abs_volume = np.abs(signed_volume.sum(axis=-1) / 6.0)

    sides = np.stack([_norm(p[:, j] - p[:, i]) for i, j in _TET_EDGES], axis=-1)
    shortest = sides.min(axis=-1)
    longest = sides.max(axis=-1)

    face_angles = np.concatenate(
        [
            np.stack(
                [
                    _angle_deg(p[:, j] - p[:, i], p[:, k] - p[:, i]),
                    _angle_deg(p[:, k] - p[:, j], p[:, i] - p[:, j]),
                    _angle_deg(p[:, i] - p[:, k], p[:, j] - p[:, k]),
                ],
                axis=-1,
            )
            for i, j, k in _TET_FACES
        ],
        axis=-1,
    )
    _, _, skew = _angle_set_stats(face_angles, 60.0)

    dihedrals = []
    for i, j in _TET_EDGES:
        opposite = [k for k in range(4) if k not in (i, j)]
        axis = _unit(p[:, j] - p[:, i])
        projections = []
        for k in opposite:
            v = p[:, k] - p[:, i]
            projections.append(v - (v * axis).sum(axis=-1, keepdims=True) * axis)
        dihedrals.append(_angle_deg(projections[0], projections[1]))
    dihedral_set = np.stack(dihedrals, axis=-1)
    min_dihedral, max_dihedral, _ = _angle_set_stats(dihedral_set, 60.0)

    # Gmsh qmTetrahedron::gamma: rho = 3 * inradius, gamma = rho / R.
    areas = (
        0.5 * _norm(np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]))
        + 0.5 * _norm(np.cross(p[:, 2] - p[:, 0], p[:, 3] - p[:, 0]))
        + 0.5 * _norm(np.cross(p[:, 1] - p[:, 0], p[:, 3] - p[:, 0]))
        + 0.5 * _norm(np.cross(p[:, 2] - p[:, 1], p[:, 3] - p[:, 1]))
    )
    a = ((p[:, 1] - p[:, 0]) ** 2).sum(axis=-1)
    b = ((p[:, 2] - p[:, 0]) ** 2).sum(axis=-1)
    c = ((p[:, 3] - p[:, 0]) ** 2).sum(axis=-1)
    big_a = ((p[:, 3] - p[:, 2]) ** 2).sum(axis=-1)
    big_b = ((p[:, 3] - p[:, 1]) ** 2).sum(axis=-1)
    big_c = ((p[:, 2] - p[:, 1]) ** 2).sum(axis=-1)
    prod_a = np.sqrt(a * big_a)
    prod_b = np.sqrt(b * big_b)
    prod_c = np.sqrt(c * big_c)
    inside = (
        (prod_a + prod_b + prod_c)
        * (prod_a + prod_b - prod_c)
        * (prod_a - prod_b + prod_c)
        * (-prod_a + prod_b + prod_c)
    )
    part_r = np.where(inside > 0.0, np.sqrt(np.where(inside > 0.0, inside, 1.0)) / 24.0, 0.0)
    rho = np.where(areas > 0.0, 9.0 * abs_volume / np.where(areas > 0.0, areas, 1.0), 0.0)
    gamma = np.where(
        (abs_volume > 0.0) & (inside > 0.0) & (part_r > 0.0),
        rho * abs_volume / np.where(part_r > 0.0, part_r, 1.0),
        0.0,
    )

    sum_l2 = (sides * sides).sum(axis=-1)
    eta = np.where(
        sum_l2 > 0.0, 12.0 * (3.0 * abs_volume) ** (2.0 / 3.0) / np.where(sum_l2 > 0.0, sum_l2, 1.0), 0.0
    )

    sj_corners = []
    for corner in range(4):
        _, others = _TET_CORNERS[corner]
        edges = [_unit(p[:, j] - p[:, corner]) for j in others]
        sj_corners.append((edges[0] * np.cross(edges[1], edges[2])).sum(axis=-1))
    sj = np.stack(sj_corners, axis=-1).min(axis=-1)

    regular = _regular_tet_points()
    gmsh_ref = _gmsh_tet_points()
    sicn = np.stack([_tet_corner_icn(p, regular, corner) for corner in range(4)], axis=-1).min(axis=-1)
    gmsh_icn = np.stack(
        [_tet_corner_icn(p, gmsh_ref, corner) for corner in range(4)], axis=-1
    ).min(axis=-1)

    return _result(
        "tet",
        {
            "volume": abs_volume,
            "min_edge": shortest,
            "max_edge": longest,
            "edge_ratio": _ratio(longest, shortest),
            "min_dihedral": min_dihedral,
            "max_dihedral": max_dihedral,
            "skewness": skew,
            "gamma": gamma,
            "eta": eta,
            "scaled_jacobian": sj,
            "min_sicn": sicn,
            "gmsh_icn": gmsh_icn,
        },
    )