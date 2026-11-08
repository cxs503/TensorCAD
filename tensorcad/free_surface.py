"""Slack tanks: where the liquid's centre of gravity actually goes, at any heel.

The free-surface correction of a stability book is usually written as one number,
``i_x / nabla``: the liquid in a partly filled tank keeps its surface horizontal while the
ship heels, so the liquid's centroid walks sideways by ``(i_x / V_liquid) phi`` and the
righting arm loses ``(i_x / nabla) sin(phi)``. Both of those are *first-order* statements in
``phi``. They are also what a designer wants 99 % of the time, which is why the classical
form is the default in :func:`tensorcad.gz.righting_arm_curve`.

What they are not is exact, and this module computes the exact answer for the one geometry
where nothing has to be approximated: a **prismatic rectangular tank**, heeled about an axis
parallel to its length, holding a fixed volume of liquid. The fluid region is the part of a
rotated rectangle below a horizontal plane, which is a polygon; its area, centroid and
free-surface chord are closed-form once the level is solved for (the level is one scalar
unknown fixed by the liquid volume, and it is found by bisection to round-off). No mesh, no
quadrature, no small-angle substitution.

Two closed forms come out of the same code and are used as independent references in the
tests:

* the upright free-surface inertia is independent of the fill level (a rectangular tank's
  ``i_x = l b^3 / 12`` does not know how deep it is filled), and
* while the free surface still runs wall to wall, its chord is ``b / cos(phi)``, so the
  inertia grows as ``i_x(phi) = i_x(0) / cos^3(phi)`` -- the widest-angled and most easily
  falsified of the classical statements about slack tanks.

The transverse walk of the liquid centroid is what the module actually solves for, and the
inertia law above is a *consequence* of it, not an input: the tests check that the two agree
rather than that the code reproduces its own assumption.

Sign conventions are the ones :mod:`tensorcad.gz` uses. ``phi > 0`` heels the ``+y``
(starboard) side down. A liquid centroid that walks towards ``+y`` in the body frame raises
the combined centre of gravity's world transverse coordinate, so the arm **loses** righting
moment: :func:`free_surface_gz_shift` is negative at positive heel, and it is odd in the
heel, exactly like the classical ``-(i_x/nabla) sin(phi)`` it replaces.

What is *not* modelled: non-prismatic tanks (a real tank has a sloped top, structure,
suction wells and a sounding table of its own), trim on top of heel, fluids of different
densities occupying one tank, and tanks whose free surface reaches the tank's ends in the
longitudinal direction (the model assumes the free surface runs the full length, which holds
for a heel about the longitudinal axis and no trim). Those are stated here rather than
silently approximated.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

__all__ = [
    "FreeSurfaceError",
    "SlackTank",
    "TankFreeSurface",
    "free_surface_gz_shift",
    "free_surface_gz_shift_curve",
    "free_surface_loss_m",
    "tank_free_surface",
]


class FreeSurfaceError(ValueError):
    """Raised when a slack tank, a heel angle or a displacement is not physical."""


@dataclass(frozen=True)
class SlackTank:
    """A rectangular prismatic tank holding ``fill_volume`` of liquid with a free surface.

    The tank's own frame: ``length`` along ``x``, ``breadth`` along ``y``, ``height`` along
    ``z``, base at ``z = 0`` and centreplane at ``y = 0``. Only *differences* of the liquid
    centroid enter the answer, so translating the tank inside the ship changes nothing and
    the origin of the tank frame is arbitrary.

    ``density_ratio`` is ``rho_liquid / rho_displacement`` -- 1.0 for a tank of sea water in
    a ship floating in sea water, 0.85 for fuel oil.
    """

    length: float
    breadth: float
    height: float
    fill_volume: float
    density_ratio: float = 1.0
    name: str = "tank"

    def __post_init__(self) -> None:
        for label, value in (
            ("length", self.length),
            ("breadth", self.breadth),
            ("height", self.height),
            ("fill_volume", self.fill_volume),
            ("density_ratio", self.density_ratio),
        ):
            if not math.isfinite(value) or value <= 0.0:
                raise FreeSurfaceError(f"{label} must be finite and positive: {value!r}")
        if self.fill_volume >= self.capacity:
            raise FreeSurfaceError(
                f"tank {self.name!r} is full ({self.fill_volume!r} of {self.capacity!r} m^3): "
                "a full tank has no free surface, so it has no free-surface correction"
            )

    @property
    def capacity(self) -> float:
        """``length * breadth * height``, m^3."""
        return self.length * self.breadth * self.height

    @property
    def fill_fraction(self) -> float:
        """Liquid volume over tank capacity."""
        return self.fill_volume / self.capacity

    @property
    def upright_inertia(self) -> float:
        """``i_x = l b^3 / 12`` of the upright free surface, m^4. Independent of the fill."""
        return self.length * self.breadth**3 / 12.0

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "length_m": self.length,
            "breadth_m": self.breadth,
            "height_m": self.height,
            "fill_volume_m3": self.fill_volume,
            "fill_fraction": self.fill_fraction,
            "density_ratio": self.density_ratio,
            "upright_inertia_m4": self.upright_inertia,
        }


@dataclass(frozen=True)
class TankFreeSurface:
    """The fluid state of one slack tank at one heel angle.

    ``level`` is the height of the horizontal free surface along the world vertical,
    expressed in the tank's own frame (``-y sin(phi) + z cos(phi)``), so it is not a sounding:
    use ``centroid_z`` and ``centroid_y`` for anything positional.
    """

    heel_angle_deg: float
    volume: float
    level: float
    free_surface_area: float
    free_surface_width: float
    inertia: float
    centroid_y: float
    centroid_z: float
    fill_fraction: float

    def to_dict(self) -> dict[str, object]:
        return {
            "heel_angle_deg": _round(self.heel_angle_deg),
            "volume_m3": _round(self.volume, 15),
            "level_m": _round(self.level, 15),
            "free_surface_area_m2": _round(self.free_surface_area),
            "free_surface_width_m": _round(self.free_surface_width),
            "inertia_m4": _round(self.inertia),
            "centroid_y_m": _round(self.centroid_y),
            "centroid_z_m": _round(self.centroid_z),
            "fill_fraction": _round(self.fill_fraction),
        }


def _round(value: float, digits: int = 12) -> float | None:
    if not math.isfinite(value):
        return None
    return round(value, digits)


def _corners(tank: SlackTank) -> tuple[tuple[float, float], ...]:
    half = 0.5 * tank.breadth
    return ((-half, 0.0), (half, 0.0), (half, tank.height), (-half, tank.height))


def _height(point: tuple[float, float], sin_phi: float, cos_phi: float) -> float:
    """World-vertical coordinate of a body point: ``-y sin(phi) + z cos(phi)``."""
    return -point[0] * sin_phi + point[1] * cos_phi


def _transverse(point: tuple[float, float], sin_phi: float, cos_phi: float) -> float:
    """World-transverse coordinate of a body point: ``y cos(phi) + z sin(phi)``."""
    return point[0] * cos_phi + point[1] * sin_phi


def _clip_below(
    polygon: Sequence[tuple[float, float]], sin_phi: float, cos_phi: float, level: float
) -> list[tuple[float, float]]:
    """Sutherland-Hodgman clip of ``polygon`` to ``-y sin(phi) + z cos(phi) <= level``."""
    out: list[tuple[float, float]] = []
    count = len(polygon)
    for index in range(count):
        first = polygon[index]
        second = polygon[(index + 1) % count]
        d_first = _height(first, sin_phi, cos_phi) - level
        d_second = _height(second, sin_phi, cos_phi) - level
        if d_first <= 0.0:
            out.append(first)
        if (d_first < 0.0 < d_second) or (d_second < 0.0 < d_first):
            t = d_first / (d_first - d_second)
            out.append(
                (
                    first[0] + t * (second[0] - first[0]),
                    first[1] + t * (second[1] - first[1]),
                )
            )
    cleaned: list[tuple[float, float]] = []
    for point in out:
        if not cleaned or point != cleaned[-1]:
            cleaned.append(point)
    if len(cleaned) > 1 and cleaned[0] == cleaned[-1]:
        cleaned.pop()
    return cleaned


def _area_centroid(polygon: Sequence[tuple[float, float]]) -> tuple[float, float, float]:
    """Signed-free area and centroid of a simple polygon (shoelace formulas)."""
    if len(polygon) < 3:
        return 0.0, 0.0, 0.0
    twice_area = 0.0
    y_moment = 0.0
    z_moment = 0.0
    count = len(polygon)
    for index in range(count):
        y0, z0 = polygon[index]
        y1, z1 = polygon[(index + 1) % count]
        cross = y0 * z1 - y1 * z0
        twice_area += cross
        y_moment += (y0 + y1) * cross
        z_moment += (z0 + z1) * cross
    area = 0.5 * twice_area
    if area == 0.0:
        return 0.0, 0.0, 0.0
    return area, y_moment / (6.0 * area), z_moment / (6.0 * area)


def _solve_level(
    tank: SlackTank, sin_phi: float, cos_phi: float
) -> tuple[float, list[tuple[float, float]]]:
    """Bisect the free-surface level until the wetted polygon has the tank's liquid volume."""
    corners = _corners(tank)
    heights = [_height(point, sin_phi, cos_phi) for point in corners]
    low, high = min(heights), max(heights)
    target = tank.fill_volume / tank.length
    span = high - low
    if not math.isfinite(span) or span <= 0.0:
        raise FreeSurfaceError(f"tank {tank.name!r} has a degenerate section")
    polygon: list[tuple[float, float]] = []
    for _ in range(200):
        middle = 0.5 * (low + high)
        polygon = _clip_below(corners, sin_phi, cos_phi, middle)
        area, _, _ = _area_centroid(polygon)
        if area < target:
            low = middle
        else:
            high = middle
        if high - low <= 1e-15 * span:
            break
    level = 0.5 * (low + high)
    polygon = _clip_below(corners, sin_phi, cos_phi, level)
    return level, polygon


def _free_surface_chord(
    polygon: Sequence[tuple[float, float]],
    sin_phi: float,
    cos_phi: float,
    level: float,
    *,
    scale: float,
) -> tuple[list[float], float]:
    """The chord the free surface cuts in a polygon, and its transverse width.

    Only the vertices that lie *on* the plane belong to the free surface: the wetted
    polygon's other corners are metres away from it and must not be allowed to set the
    width. Vertices are accepted within double-precision round-off of the plane, which is
    what the clipping interpolation leaves behind.
    """
    tolerance = 1e-12 * (abs(scale) + abs(level) + 1.0)
    on_plane: list[tuple[float, float]] = []
    for point in polygon:
        if abs(_height(point, sin_phi, cos_phi) - level) <= tolerance:
            on_plane.append(point)
    if len(on_plane) < 2:
        raise FreeSurfaceError(
            "the free surface is a point or empty: the tank has no liquid surface at this heel"
        )
    transverse = [_transverse(point, sin_phi, cos_phi) for point in on_plane]
    return transverse, max(transverse) - min(transverse)


def tank_free_surface(tank: SlackTank, heel_deg: float = 0.0) -> TankFreeSurface:
    """The exact fluid state of ``tank`` heeled ``heel_deg`` degrees about its long axis.

    The liquid volume is held fixed and the free-surface level is solved for; the answer is a
    closed-form polygon computation on top of that one scalar solve, so it carries
    double-precision round-off and nothing else.
    """
    if not isinstance(tank, SlackTank):
        raise FreeSurfaceError(f"expected a SlackTank, got {type(tank).__name__}")
    if not math.isfinite(heel_deg):
        raise FreeSurfaceError(f"heel angle must be finite: {heel_deg!r}")
    heel_rad = math.radians(heel_deg)
    sin_phi, cos_phi = math.sin(heel_rad), math.cos(heel_rad)
    level, polygon = _solve_level(tank, sin_phi, cos_phi)
    area, centroid_y, centroid_z = _area_centroid(polygon)
    if area <= 0.0:
        raise FreeSurfaceError(f"tank {tank.name!r} holds no liquid at {heel_deg!r} degrees")
    width = _free_surface_chord(polygon, sin_phi, cos_phi, level, scale=tank.height + tank.breadth)[1]
    return TankFreeSurface(
        heel_angle_deg=float(heel_deg),
        volume=area * tank.length,
        level=level,
        free_surface_area=width * tank.length,
        free_surface_width=width,
        inertia=tank.length * width**3 / 12.0,
        centroid_y=centroid_y,
        centroid_z=centroid_z,
        fill_fraction=tank.fill_fraction,
    )


def free_surface_gz_shift(tank: SlackTank, heel_deg: float, displaced_volume: float) -> float:
    """The metres a slack ``tank`` takes off ``GZ`` at ``heel_deg``.

    ``displaced_volume`` is the ship's volume of displacement in m^3; the liquid's weight
    fraction is ``density_ratio * fill_volume / displaced_volume``. The result is **negative
    at positive heel** (the arm loses righting moment) and odd in the heel, so it is exactly
    the signed quantity to add to a righting arm that was computed with the liquid frozen.
    """
    if not math.isfinite(displaced_volume) or displaced_volume <= 0.0:
        raise FreeSurfaceError(
            f"displaced volume must be finite and positive: {displaced_volume!r}"
        )
    upright = tank_free_surface(tank, 0.0)
    heeled = tank_free_surface(tank, heel_deg)
    heel_rad = math.radians(heeled.heel_angle_deg)
    sin_phi, cos_phi = math.sin(heel_rad), math.cos(heel_rad)
    # World-transverse walk of the liquid centroid between upright and heeled.
    walk = (heeled.centroid_y - upright.centroid_y) * cos_phi + (
        heeled.centroid_z - upright.centroid_z
    ) * sin_phi
    weight_fraction = tank.density_ratio * tank.fill_volume / displaced_volume
    return -weight_fraction * walk


def free_surface_gz_shift_curve(
    tanks: Sequence[SlackTank], angles_deg: Sequence[float], displaced_volume: float
) -> tuple[float, ...]:
    """Sum of :func:`free_surface_gz_shift` over ``tanks`` at every heel angle."""
    checked = tuple(tanks)
    for tank in checked:
        if not isinstance(tank, SlackTank):
            raise FreeSurfaceError(f"expected a SlackTank, got {type(tank).__name__}")
    return tuple(
        math.fsum(free_surface_gz_shift(tank, angle, displaced_volume) for tank in checked)
        for angle in angles_deg
    )


def free_surface_loss_m(
    tanks: Sequence[SlackTank], displaced_volume: float
) -> float:
    """``sum(rho_liquid i_x) / nabla``: the metres of ``GM`` the slack tanks cost.

    This is the classical number a stability book quotes, and it is the *slope* of
    :func:`free_surface_gz_shift` at zero heel -- not a flat offset to be subtracted from
    every ordinate of the curve.
    """
    if not math.isfinite(displaced_volume) or displaced_volume <= 0.0:
        raise FreeSurfaceError(
            f"displaced volume must be finite and positive: {displaced_volume!r}"
        )
    total = math.fsum(
        tank.density_ratio * tank.upright_inertia for tank in tanks
    )
    if total < 0.0:
        raise FreeSurfaceError("a free-surface inertia cannot be negative")
    return total / displaced_volume