"""An exact righting arm for a damaged *box* barge, computed without the mesh.

The library route integrates clipped section polygons at Gauss stations, so every
number it produces carries a discretisation error that has to be argued away. A box
hull with a box compartment has no such excuse: the plane
``z = a + tan(heel) y`` is independent of ``x`` when the trim is zero, so the net
buoyancy is a one-dimensional integral over ``y`` of a function that is *piecewise
linear* -- and each piece is a polynomial of degree at most two whose integral is
written down in closed form below.

Given the waterline, everything follows exactly:

* net volume ``V(a) = L INT h_h dy - mu L_c INT h_c dy``,
* its ``y``-, ``z``- and ``x``- moments (the same integral with ``y``, ``z`` and ``x``
  weights),
* and from those, the arm ``GZ = (y_B - y_G) cos(phi) + (KB - KG) sin(phi)``.

``V(a)`` is strictly increasing in ``a``, so the floating attitude is a bisection, and
the answer is exact to the double-precision round-off of the closed forms -- there is
no grid to refine. This is the oracle :mod:`tests.test_damage_gz` holds the library to.

Restriction: the trim is taken as zero. That is not an assumption bolted on for
convenience -- it is the *solution* whenever the flooded compartment and the hull are
both symmetric about midships, which is the case every oracle test uses (and the
library is separately asserted to return ``lcb == x_G`` and a vanishing trim).
"""
from __future__ import annotations

import math

__all__ = ["ExactArm", "exact_damaged_arm", "wall_sided_arm"]


def _clamped_pieces(
    y0: float, y1: float, a: float, c: float, low: float, high: float
) -> list[tuple[float, float, float, float]]:
    """``[y0, y1]`` split where ``clamp(a + c y, low, high) - low`` changes branch.

    Returns ``(ya, yb, alpha, beta)`` pieces with ``h(y) = alpha + beta y`` exactly.
    """
    if y1 <= y0:
        return []
    breaks = {y0, y1}
    if c != 0.0:
        for level in (low, high):
            y = (level - a) / c
            if y0 < y < y1:
                breaks.add(y)
    edges = sorted(breaks)
    pieces: list[tuple[float, float, float, float]] = []
    for index in range(len(edges) - 1):
        ya, yb = edges[index], edges[index + 1]
        if yb <= ya:
            continue
        mid = 0.5 * (ya + yb)
        value = a + c * mid
        if value <= low:
            pieces.append((ya, yb, 0.0, 0.0))
        elif value >= high:
            pieces.append((ya, yb, high - low, 0.0))
        else:
            pieces.append((ya, yb, a - low, c))
    return pieces


def _moments(ya: float, yb: float, alpha: float, beta: float) -> tuple[float, float, float]:
    """``(INT h dy, INT y h dy, INT h^2 dy)`` over one linear piece, in closed form."""
    d1 = yb - ya
    d2 = yb * yb - ya * ya
    d3 = yb * yb * yb - ya * ya * ya
    integral_h = alpha * d1 + beta * 0.5 * d2
    integral_yh = alpha * 0.5 * d2 + beta * d3 / 3.0
    integral_h2 = alpha * alpha * d1 + alpha * beta * d2 + beta * beta * d3 / 3.0
    return integral_h, integral_yh, integral_h2


class ExactArm:
    """One exact point of the damaged box's arm curve."""

    def __init__(
        self,
        heel_deg: float,
        base_z: float,
        volume: float,
        lcb: float,
        tcb: float,
        kb: float,
        righting_arm: float,
    ) -> None:
        self.heel_deg = heel_deg
        self.base_z = base_z
        self.volume = volume
        self.lcb = lcb
        self.tcb = tcb
        self.kb = kb
        self.righting_arm = righting_arm

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return (
            f"ExactArm(heel={self.heel_deg}, z0={self.base_z!r}, V={self.volume!r}, "
            f"GZ={self.righting_arm!r})"
        )


def exact_damaged_arm(
    *,
    heel_deg: float,
    displacement_volume: float,
    length: float,
    beam: float,
    keel_z: float,
    deck_z: float,
    kg: float,
    y_gravity: float = 0.0,
    compartment: tuple[float, float, float, float, float] | None = None,
) -> ExactArm:
    """Exact damaged arm of a box hull, by closed-form one-dimensional integration.

    ``compartment`` is ``(x0, x1, z0, z1, permeability)`` spanning the full breadth, or
    ``None`` for the intact ship. ``kg`` is above the keel and ``y_gravity`` is the
    centre of gravity's transverse position in the body frame.
    """
    if displacement_volume <= 0.0:
        raise ValueError("displacement_volume must be positive")
    phi = math.radians(heel_deg)
    tangent = math.tan(phi)
    y_low, y_high = -0.5 * beam, 0.5 * beam
    comp_x0 = comp_x1 = comp_z0 = comp_z1 = mu = 0.0
    if compartment is not None:
        comp_x0, comp_x1, comp_z0, comp_z1, mu = compartment
        if comp_z1 <= comp_z0:
            raise ValueError("compartment z1 must exceed z0")

    def integrals(base: float) -> tuple[float, float, float, float]:
        """``(net volume, INT x dV, INT y dV, INT z dV)`` at waterline ``base``."""
        h_h = _clamped_pieces(y_low, y_high, base, tangent, keel_z, deck_z)
        volume = x_moment = y_moment = z_moment = 0.0
        for ya, yb, alpha, beta in h_h:
            integral_h, integral_yh, integral_h2 = _moments(ya, yb, alpha, beta)
            volume += length * integral_h
            x_moment += 0.5 * length * length * integral_h
            y_moment += length * integral_yh
            z_moment += length * (keel_z * integral_h + 0.5 * integral_h2)
        if compartment is not None and mu > 0.0:
            h_c = _clamped_pieces(y_low, y_high, base, tangent, comp_z0, comp_z1)
            length_c = comp_x1 - comp_x0
            for ya, yb, alpha, beta in h_c:
                integral_h, integral_yh, integral_h2 = _moments(ya, yb, alpha, beta)
                volume -= mu * length_c * integral_h
                x_moment -= (
                    mu * 0.5 * (comp_x1 * comp_x1 - comp_x0 * comp_x0) * integral_h
                )
                y_moment -= mu * length_c * integral_yh
                z_moment -= mu * length_c * (comp_z0 * integral_h + 0.5 * integral_h2)
        return volume, x_moment, y_moment, z_moment

    low = keel_z - 0.5 * beam * abs(tangent) - 1.0
    high = deck_z + 0.5 * beam * abs(tangent) + 1.0
    for _ in range(400):
        mid = 0.5 * (low + high)
        volume, *_ = integrals(mid)
        if volume < displacement_volume:
            low = mid
        else:
            high = mid
    base = 0.5 * (low + high)
    volume, x_moment, y_moment, z_moment = integrals(base)
    lcb = x_moment / volume
    tcb = y_moment / volume
    kb = z_moment / volume - keel_z
    arm = (tcb - y_gravity) * math.cos(phi) + (kb - kg) * math.sin(phi)
    return ExactArm(
        heel_deg=heel_deg,
        base_z=base,
        volume=volume,
        lcb=lcb,
        tcb=tcb,
        kb=kb,
        righting_arm=arm,
    )


def wall_sided_arm(heel_deg: float, *, gm: float, bm: float) -> float:
    """The classical wall-sided box formula ``GZ = sin(phi) (GM + BM/2 tan^2 phi)``."""
    phi = math.radians(heel_deg)
    return math.sin(phi) * (gm + 0.5 * bm * math.tan(phi) ** 2)