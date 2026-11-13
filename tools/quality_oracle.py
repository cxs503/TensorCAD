"""Independent recomputation of the canonical mesh-quality measures.

This file exists to *disagree* with :mod:`tensorcad.geomquality.quality` if that module
is wrong.  It is written on purpose along different arithmetic routes, and it never
imports the library it audits:

===================  =====================================  ============================
measure              library route (:mod:`.quality`)        oracle route (here)
===================  =====================================  ============================
area                 cross product ``|u x v|/2``            Heron from the three sides
angles               unit edge vectors at each apex         law of cosines from the sides
aspect_ratio         edge lengths + cross-product area      Heron area (and, exactly,
                                                           ``fractions.Fraction``)
radius_ratio         ``8A^2/(s l_a l_b l_c)``               ``2r/R`` with ``r=A/s``,
                                                           ``R=abc/(4A)``
scaled_jacobian      ``2A/(l_mid l_max)``                   min over the three apices of
                                                           ``2A/(|e_i||e_j|)``
condition number     algebraic invariants of the Jacobian   LAPACK SVD (``numpy``)
===================  =====================================  ============================

Two entry points:

``oracle_metrics(p0, p1, p2)``
    The floating-point route, used on the hull meshes.
``oracle_metrics_exact(p0, p1, p2)``
    The same quantities on meshes whose coordinates are rational, carried through
    :class:`fractions.Fraction` until the single unavoidable square root.  A float
    library that agrees with this to ``1e-12`` is not merely self-consistent: its
    rounding is bounded.

Run it as a script to print the oracle's view of the canonical triangles::

    python tools/quality_oracle.py
"""
from __future__ import annotations

import math
from fractions import Fraction
from typing import Any

import numpy as np

_SQRT3 = math.sqrt(3.0)

METRICS = (
    "area",
    "min_angle_deg",
    "max_angle_deg",
    "edge_ratio",
    "aspect_ratio",
    "aspect_frobenius",
    "gamma",
    "radius_ratio",
    "scaled_jacobian",
    "condition_number_inverse",
)


def _sides(p0, p1, p2):
    return (
        math.dist(p1, p2),  # opposite p0
        math.dist(p2, p0),  # opposite p1
        math.dist(p0, p1),  # opposite p2
    )


def oracle_metrics(p0, p1, p2) -> dict[str, Any]:
    """Every measure of :data:`METRICS`, by the routes listed in the module docstring."""
    x, y, z = _sides(p0, p1, p2)
    s = 0.5 * (x + y + z)
    heron = s * (s - x) * (s - y) * (s - z)
    area = math.sqrt(heron) if heron > 0.0 else 0.0
    if area <= 0.0:
        return {
            "area": area,
            "min_angle_deg": 0.0,
            "max_angle_deg": 180.0,
            "edge_ratio": math.inf,
            "aspect_ratio": math.inf,
            "aspect_frobenius": math.inf,
            "gamma": 0.0,
            "radius_ratio": 0.0,
            "scaled_jacobian": 0.0,
            "condition_number_inverse": 0.0,
            "degenerate": True,
        }
    perimeter = x + y + z
    squares = x * x + y * y + z * z
    angles = []
    for opposite, first, second in ((x, y, z), (y, z, x), (z, x, y)):
        cosine = (first * first + second * second - opposite * opposite) / (
            2.0 * first * second
        )
        angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
    frobenius = squares / (4.0 * _SQRT3 * area)
    inradius = area / s
    circumradius = (x * y * z) / (4.0 * area)
    # Scaled Jacobian: the worst of the three possible base vertices.
    jacobians = []
    for apex, first, second in (
        (p0, p1, p2),
        (p1, p2, p0),
        (p2, p0, p1),
    ):
        e1 = math.dist(apex, first)
        e2 = math.dist(apex, second)
        if e1 > 0.0 and e2 > 0.0:
            jacobians.append(2.0 * area / (e1 * e2))
    u = np.subtract(p1, p0)
    v = np.subtract(p2, p0)
    singular = np.linalg.svd(np.stack((u, v), axis=1), compute_uv=False)
    sigma_max, sigma_min = float(singular[0]), float(singular[-1])
    return {
        "area": area,
        "min_angle_deg": min(angles),
        "max_angle_deg": max(angles),
        "edge_ratio": max(x, y, z) / min(x, y, z),
        "aspect_ratio": max(x, y, z) * perimeter / (4.0 * _SQRT3 * area),
        "aspect_frobenius": frobenius,
        "gamma": 1.0 / frobenius,
        "radius_ratio": 2.0 * inradius / circumradius,
        "scaled_jacobian": min(jacobians),
        "condition_number_inverse": sigma_min / sigma_max if sigma_max > 0.0 else 0.0,
        "degenerate": False,
    }


def _fdist(p, q) -> Fraction:
    dx = Fraction(p[0]) - Fraction(q[0])
    dy = Fraction(p[1]) - Fraction(q[1])
    dz = Fraction(p[2]) - Fraction(q[2])
    return dx * dx + dy * dy + dz * dz  # squared distance, exact


def oracle_metrics_exact(p0, p1, p2) -> dict[str, Any]:
    """Exact-arithmetic route; coordinates must be rational (int, Fraction, decimals)."""
    sq = (_fdist(p0, p1), _fdist(p1, p2), _fdist(p2, p0))
    total = sum(sq)
    # Heron in squared form: 16 A^2 = 2(a^2 b^2 + b^2 c^2 + c^2 a^2) - (a^4 + b^4 + c^4)
    a2, b2, c2 = sq
    sixteen_area_sq = 2 * (a2 * b2 + b2 * c2 + c2 * a2) - (a2 * a2 + b2 * b2 + c2 * c2)
    if sixteen_area_sq <= 0:
        return {"area": 0.0, "degenerate": True}
    area = math.sqrt(float(sixteen_area_sq)) / 4.0
    lengths = [math.sqrt(float(value)) for value in sq]
    x, y, z = lengths
    s = 0.5 * (x + y + z)
    perimeter = x + y + z
    frobenius = float(total) / (4.0 * _SQRT3 * area)
    angles = []
    for opposite, first, second in ((x, y, z), (y, z, x), (z, x, y)):
        cosine = (first * first + second * second - opposite * opposite) / (
            2.0 * first * second
        )
        angles.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
    inradius = area / s
    circumradius = (x * y * z) / (4.0 * area)
    return {
        "area": area,
        "min_angle_deg": min(angles),
        "max_angle_deg": max(angles),
        "edge_ratio": max(lengths) / min(lengths),
        "aspect_ratio": max(lengths) * perimeter / (4.0 * _SQRT3 * area),
        "aspect_frobenius": frobenius,
        "gamma": 1.0 / frobenius,
        "radius_ratio": 2.0 * inradius / circumradius,
        "scaled_jacobian": 2.0 * area / (sorted(lengths)[1] * sorted(lengths)[2]),
        "condition_number_inverse": _exact_sicn(p0, p1, p2),
        "degenerate": False,
    }


def _exact_sicn(p0, p1, p2) -> float:
    """sigma_min/sigma_max of the exact 2x2 Gram matrix, via its exact eigenvalues."""
    u = tuple(Fraction(b) - Fraction(a) for a, b in zip(p0, p1))
    v = tuple(Fraction(b) - Fraction(a) for a, b in zip(p0, p2))
    guu = sum(value * value for value in u)
    gj = sum(value * value for value in v)
    guv = sum(left * right for left, right in zip(u, v))
    # Eigenvalues of [[guu, guv], [guv, gj]] are the squared singular values.
    half_trace = (guu + gj) / 2
    disc = (half_trace * half_trace) - (guu * gj - guv * guv)
    root = math.sqrt(float(disc))
    large = float(half_trace) + root
    small = float(half_trace) - root
    return math.sqrt(small / large) if large > 0.0 else 0.0


if __name__ == "__main__":
    cases = {
        "equilateral": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, _SQRT3 / 2.0, 0.0)),
        "right_isosceles": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        "three_four_five": ((0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (3.0, 4.0, 0.0)),
        "sliver": ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 1e-4, 0.0)),
    }
    for name, triangle in cases.items():
        values = oracle_metrics(*triangle)
        print(name, {key: round(value, 12) for key, value in values.items() if key != "degenerate"})
    exact = oracle_metrics_exact((0.0, 0.0, 0.0), (3.0, 0.0, 0.0), (3.0, 4.0, 0.0))
    print("three_four_five (exact)", {key: round(value, 15) for key, value in exact.items() if key != "degenerate"})