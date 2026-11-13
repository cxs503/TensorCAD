"""An elliptic-section prism with a closed-form righting arm.

:mod:`tensorcad.righting_arm` verifies its heeled-displacement integrator against two
analytic families: the wall-sided box (:mod:`tensorcad.box_gz`) and the circular
section. Both are degenerate: the box is wall-sided (so the metacentric formula
``BM = I/nabla`` is exact) and the circle has ``I/nabla`` equal to ``KB`` by accident of
geometry. Neither exercises the *coupled* ``y``-``z`` curvature term that a real hull
section has, and neither can catch a sign error in the ``(a^2 - b^2)`` factor, because
for both families that factor is zero or trivially positive.

This module adds the missing family: a prism of **elliptic** cross-section with
semi-axes ``a`` (transverse, ``y``) and ``b`` (vertical, ``z``), floating at exactly
half immersion so the upright waterplane always passes through the centre of the
ellipse. For that geometry the immersed section at every heel is an exact half-ellipse
and everything can be written in closed form -- no series, no quadrature, no mesh. It
is therefore an *independent* route to the same number the mesh integrator reports,
for a shape where the two invariants above are both non-trivial.

Closed form
-----------
Let the body be heeled by ``phi`` about its own longitudinal axis. In the body frame the
world-vertical is ``n = (-sin phi, cos phi)`` and the world-horizontal is
``u = (cos phi, sin phi)``. The immersed half-ellipse is the linear image
``diag(a, b)`` of a half unit disc, whose centroid sits ``4/(3 pi)`` from the chord
along ``-m`` with ``m = diag(a, b) n``, so

    B - C = -(4/(3 pi)) * diag(a^2, b^2) n / |m|,     |m| = sqrt(a^2 sin^2 + b^2 cos^2).

The righting arm is ``GZ = KN - KG sin(phi)`` with ``KN`` the horizontal distance from
the keel to the vertical through ``B``:

    KN(phi) = b sin(phi) + (4/(3 pi)) (a^2 - b^2) sin(phi) cos(phi)
              / sqrt(a^2 sin^2(phi) + b^2 cos^2(phi)).

At ``a = b = R`` this collapses to ``KN = R sin(phi)``, i.e. ``GZ = (R - KG) sin(phi)``,
which is the exact answer for a circular caisson -- and notably *not* the wall-sided
``(KB + BM - KG) sin(phi)``, because ``B`` does not rise with heel on a circular
section the way the wall-sided derivation assumes. The module documents and tests that
difference rather than hiding it.

Only the Python standard library is used.
"""

from __future__ import annotations

import math
from typing import Any

from .mesh import Face, Mesh, Vec3

__all__ = [
    "EllipticSectionError",
    "elliptic_gm0",
    "elliptic_kn",
    "elliptic_prism_mesh",
    "elliptic_righting_arm",
    "elliptic_righting_arm_curve",
    "wall_sided_gm0",
]

#: Quadrature-free constant 4/(3 pi), the half-disc centroid factor.
_FOUR_OVER_THREE_PI = 4.0 / (3.0 * math.pi)


class EllipticSectionError(ValueError):
    """Raised when the elliptic-section family is asked for a non-physical input."""


def _check_axes(a: float, b: float) -> tuple[float, float]:
    for name, value in (("a", a), ("b", b)):
        if not math.isfinite(value) or value <= 0.0:
            raise EllipticSectionError(f"{name} must be finite and positive: {value!r}")
    return float(a), float(b)


def _check_heel(heel_deg: float) -> float:
    if not math.isfinite(heel_deg) or abs(heel_deg) >= 180.0:
        raise EllipticSectionError(f"heel_deg must be finite and |heel| < 180: {heel_deg!r}")
    return float(heel_deg)


def elliptic_kn(a: float, b: float, heel_deg: float) -> float:
    """``KN`` -- horizontal distance from the keel to the vertical through ``B``, m."""
    a, b = _check_axes(a, b)
    phi = math.radians(_check_heel(heel_deg))
    sin_phi, cos_phi = math.sin(phi), math.cos(phi)
    scale = math.sqrt(a * a * sin_phi * sin_phi + b * b * cos_phi * cos_phi)
    coupled = _FOUR_OVER_THREE_PI * (a * a - b * b) * sin_phi * cos_phi / scale
    return b * sin_phi + coupled


def elliptic_righting_arm(a: float, b: float, heel_deg: float, *, kg: float) -> float:
    """``GZ(phi)`` for the half-immersed elliptic prism, m.

    ``kg`` is measured from the keel, as everywhere else in this package. The result is
    analytic: the only operations are ``sin``, ``cos`` and one square root.
    """
    if not math.isfinite(kg):
        raise EllipticSectionError(f"kg must be finite: {kg!r}")
    phi = math.radians(_check_heel(heel_deg))
    return elliptic_kn(a, b, heel_deg) - float(kg) * math.sin(phi)


def elliptic_gm0(a: float, b: float, *, kg: float) -> float:
    """Initial metacentric height of the half-immersed elliptic prism, m.

    ``GM0 = dGZ/dphi`` at upright, from the closed form (the coupled term's derivative
    brings down ``1/D(0) = 1/b``, not ``1/a``):
    ``GM0 = b - kg + (4/(3 pi)) (a^2 - b^2) / b``.
    """
    if not math.isfinite(kg):
        raise EllipticSectionError(f"kg must be finite: {kg!r}")
    a, b = _check_axes(a, b)
    return b - float(kg) + _FOUR_OVER_THREE_PI * (a * a - b * b) / b


def wall_sided_gm0(a: float, b: float, *, kg: float) -> float:
    """The *wall-sided* prediction ``KB + BM - KG`` for the same section, m.

    Reported next to :func:`elliptic_gm0` because the two agree only in the
    ``a -> b`` limit's neighbourhood and disagree completely at ``a = b``: for a circle
    this returns ``8R/(3 pi) - KG`` while the exact slope is ``R - KG``. Keeping both
    makes the point measurable instead of rhetorical.
    """
    if not math.isfinite(kg):
        raise EllipticSectionError(f"kg must be finite: {kg!r}")
    a, b = _check_axes(a, b)
    kb = 2.0 * b / math.pi  # half-ellipse centroid above the keel
    waterplane_inertia = (2.0 / 3.0) * a ** 3  # 2a * (2a)^3 / 12 per unit length
    volume = math.pi * a * b  # per unit length
    return kb + waterplane_inertia / volume - float(kg)


def elliptic_righting_arm_curve(
    a: float,
    b: float,
    *,
    kg: float,
    heels_deg: tuple[float, ...] | list[float] = tuple(range(0, 91, 5)),
) -> dict[str, Any]:
    """The analytic ``GZ`` ladder as a plain record, ready to serialise."""
    heels = tuple(float(heel) for heel in heels_deg)
    levers = tuple(elliptic_righting_arm(a, b, heel, kg=kg) for heel in heels)
    peak_heel = heels[max(range(len(levers)), key=lambda index: levers[index])]
    return {
        "a": float(a),
        "b": float(b),
        "kg": float(kg),
        "heels_deg": heels,
        "gz": levers,
        "gm0_analytic": elliptic_gm0(a, b, kg=kg),
        "gm0_wall_sided": wall_sided_gm0(a, b, kg=kg),
        "peak_heel_deg": peak_heel,
        "peak_gz": max(levers),
    }


def elliptic_prism_mesh(
    a: float,
    b: float,
    length: float,
    *,
    ring: int = 192,
    waterline_z: float = 0.0,
) -> Mesh:
    """A watertight prism of elliptic cross-section, centred on ``waterline_z``.

    The side surface is a surface of translation along ``x``, so a single station pair
    is *exact*: the only discretisation is the ``ring`` points around the ellipse. The
    two end caps are fans, closed, with outward winding, which is what
    :func:`tensorcad.righting_arm.righting_arm_curve` requires.
    """
    a, b = _check_axes(a, b)
    if not math.isfinite(length) or length <= 0.0:
        raise EllipticSectionError(f"length must be finite and positive: {length!r}")
    if ring < 8:
        raise EllipticSectionError(f"ring must be at least 8: {ring!r}")
    if not math.isfinite(waterline_z):
        raise EllipticSectionError(f"waterline_z must be finite: {waterline_z!r}")

    half = 0.5 * float(length)
    vertices: list[Vec3] = []
    for x in (-half, half):
        for index in range(ring):
            angle = 2.0 * math.pi * index / ring
            vertices.append((x, a * math.cos(angle), waterline_z + b * math.sin(angle)))
    aft_centre = len(vertices)
    vertices.append((-half, 0.0, waterline_z))
    fore_centre = len(vertices)
    vertices.append((half, 0.0, waterline_z))

    faces: list[Face] = []
    for index in range(ring):
        nxt = (index + 1) % ring
        low, low_next = index, nxt
        high, high_next = ring + index, ring + nxt
        # side: outward winding seen from outside (x from -half to +half)
        faces.append((low, low_next, high_next))
        faces.append((low, high_next, high))
        # aft cap (x = -half): normal towards -x
        faces.append((aft_centre, nxt, index))
        # fore cap (x = +half): normal towards +x
        faces.append((fore_centre, ring + index, ring + nxt))

    return Mesh(
        vertices=tuple(vertices),
        faces=tuple(faces),
        coordinate_frame="right-handed, x aft->fore, z up, waterline at waterline_z",
    )