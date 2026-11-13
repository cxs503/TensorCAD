r"""Still-water hull girder loads: buoyancy along the hull, shear force, bending moment.

Every other structural module in this package takes its load as a given:
:mod:`tensorcad.girder` says so in its own docstring ("Not a stress calculation. No
bending moment is applied"), :mod:`tensorcad.shear_flow` says "no load model: no bending
moment", and :mod:`tensorcad.bonjean` describes itself as "the input to still-water shear
force and bending moment". This module is that missing link: it turns the hull geometry
into the buoyancy load per unit length, adds the weight distribution the caller supplies,
and integrates the difference twice with exact piecewise integration.

**The load model.** With ``x`` running from the aft end to the forward end,

    w(x) = p(x) - b(x)          [N/m]     weight minus buoyancy, positive down
    Q(x) = \int_{x0}^{x} w dx'  [N]       shear force, ``dQ/dx = w``
    M(x) = \int_{x0}^{x} Q dx'  [N m]     bending moment, ``dM/dx = Q``

so a section carries the moment of everything abaft it, and hogging appears as a
positive ``M`` with these signs (weight abaft pushing down, buoyancy abaft pushing up).

The convention is stated once and then used everywhere, because the sign of a bending
moment is the single most-reported-wrong number in hull girder work: ``tests/test_girder_loads.py``
checks it against hand-derived textbook cases rather than against itself.

**Why the integration is exact.** ``b(x)`` is known only at the stations that were cut, so
it enters as a piecewise-linear function through those samples -- that is the *definition*
of the load diagram here, not an approximation of it. ``p(x)`` is piecewise constant
(:class:`WeightItem` is a force spread uniformly over an interval). The sum is therefore
piecewise linear, and

    \Delta Q = h (w_a + w_b) / 2 ,      \Delta M = h Q_a + h^2 (2 w_a + w_b) / 6

integrate it *exactly per segment*. Across a jump in ``p`` the segment ending at the jump
uses the left limit of ``w`` and the segment starting there the right limit; ``Q`` and ``M``
are continuous, so nothing is smeared. All the error in a reported load is therefore
mesh error in ``b(x)``, which the convergence study in ``docs/girder_loads.md`` measures
separately -- quadrature error is zero by construction, and the box-barge cases in the
tests come out bit-exact because a wall-sided box has a linear ``b(x)`` too.

**What the two residuals mean.** ``Q(x1)`` is the total unbalanced force and ``M(x1)`` is
the total unbalanced moment about the forward end. A floating body in equilibrium needs
both to vanish: ``Q(x1) = 0`` is force balance (weight = displacement), and ``M(x1) = 0``
is *trim* balance (the weight resultant acts on the same vertical as the buoyancy
resultant). They are returned rather than enforced, because a caller who has just invented
a weight distribution needs to see by how much it is out of balance instead of getting an
exception; a caller who wants them zero can solve for the missing item's position.

**What is deliberately not here.** No wave-induced loads (this is the *still-water* case
the rules add the wave component to), no dynamic slam/whipping, no shear lag, no
transverse strength, and no rule criterion -- IACS UR S11 and its ilk are policy, applied
downstream to the numbers computed here. The structural check that consumes ``M`` is
:mod:`tensorcad.girder` (``sigma = M z / I``) and :mod:`tensorcad.shear_flow` (where
``Q`` travels in the plating); this module supplies their two inputs.

Only the Python standard library is used.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import Any

from .bonjean import section_properties
from .geometry import Mesh

__all__ = [
    "BuoyancyDistribution",
    "GirderLoadError",
    "GRAVITY",
    "HullGirderLoads",
    "SEA_WATER_DENSITY",
    "WeightItem",
    "buoyancy_distribution",
    "hull_girder_loads",
    "still_water_loads",
]

#: Standard gravity, m/s^2.
GRAVITY = 9.80665
#: Density of sea water, kg/m^3.
SEA_WATER_DENSITY = 1025.0


class GirderLoadError(ValueError):
    """Raised when the load diagram, the stations or the weight items are not physical."""


@dataclass(frozen=True)
class WeightItem:
    """A force ``force`` (N, positive down) spread uniformly over ``[x_start, x_end]``.

    The item occupies the half-open interval ``[x_start, x_end)``: a point load at a
    section is the limit of a narrow item, which is how ``tests/test_girder_loads.py``
    reproduces the textbook ``PL/8`` midspan moment of a box barge. ``force`` may be
    negative (a buoyant item, a lifting force), it just may not be infinite.
    """

    x_start: float
    x_end: float
    force: float
    label: str = ""

    def __post_init__(self) -> None:
        values = (self.x_start, self.x_end, self.force)
        for name, value in zip(("x_start", "x_end", "force"), values):
            try:
                number = float(value)
            except (TypeError, ValueError) as error:  # pragma: no cover - defensive
                raise GirderLoadError(f"{name} must be a number, got {value!r}") from error
            if not math.isfinite(number):
                raise GirderLoadError(f"{name} must be finite, got {value!r}")
        if not self.x_end > self.x_start:
            raise GirderLoadError(
                f"{self.label or 'weight item'}: x_end={self.x_end!r} must exceed "
                f"x_start={self.x_start!r}"
            )

    @property
    def density(self) -> float:
        """The item's load per unit length, N/m."""
        return self.force / (self.x_end - self.x_start)

    @property
    def centroid(self) -> float:
        """Where the item's resultant acts, m."""
        return 0.5 * (self.x_start + self.x_end)


@dataclass(frozen=True)
class BuoyancyDistribution:
    """Buoyancy per unit length at the station samples, and what it integrates to.

    ``loads[i]`` is ``density * gravity * areas[i]`` in N/m at ``stations[i]``. The
    entries of ``stations`` are the knots of the piecewise-linear diagram that
    :func:`still_water_loads` integrates, so the tuple is returned with the
    distribution instead of being hidden inside it.
    """

    stations: tuple[float, ...]
    areas: tuple[float, ...]
    loads: tuple[float, ...]
    volume: float
    lcb: float
    waterline_z: float
    trim_angle_deg: float
    reference_x: float
    density: float
    gravity: float
    draft_aft: float = 0.0
    draft_fore: float = 0.0

    @property
    def displacement_force(self) -> float:
        """``density * gravity * volume``, N."""
        return self.density * self.gravity * self.volume

    @property
    def stations_count(self) -> int:
        return len(self.stations)

    def to_dict(self) -> dict[str, Any]:
        return {
            "stations": list(self.stations),
            "areas": list(self.areas),
            "loads": list(self.loads),
            "volume": self.volume,
            "lcb": self.lcb,
            "waterline_z": self.waterline_z,
            "trim_angle_deg": self.trim_angle_deg,
            "reference_x": self.reference_x,
            "density": self.density,
            "gravity": self.gravity,
            "draft_aft": self.draft_aft,
            "draft_fore": self.draft_fore,
            "displacement_force": self.displacement_force,
        }


@dataclass(frozen=True)
class HullGirderLoads:
    """Shear force and bending moment along the hull, with the balances they imply."""

    stations: tuple[float, ...]
    shear: tuple[float, ...]
    moment: tuple[float, ...]
    total_weight: float
    total_buoyancy: float
    lcb_weight: float
    lcb_buoyancy: float
    #: ``Q`` at the forward end -- the unbalanced force, N.
    residual_force: float
    #: ``M`` at the forward end -- the unbalanced moment about ``x1``, N m.
    residual_moment: float
    items: tuple[WeightItem, ...] = ()

    @property
    def length(self) -> float:
        return self.stations[-1] - self.stations[0]

    @property
    def max_shear(self) -> float:
        """The largest ``|Q|``, N."""
        return max(abs(value) for value in self.shear)

    @property
    def max_moment(self) -> float:
        """The largest ``|M|``, N m."""
        return max(abs(value) for value in self.moment)

    def extreme_shear(self) -> tuple[float, float]:
        """``(x, Q)`` at the largest ``|Q|``."""
        index = max(range(len(self.shear)), key=lambda i: abs(self.shear[i]))
        return self.stations[index], self.shear[index]

    def extreme_moment(self) -> tuple[float, float]:
        """``(x, M)`` at the largest ``|M|`` -- the section the modulus is checked against."""
        index = max(range(len(self.moment)), key=lambda i: abs(self.moment[i]))
        return self.stations[index], self.moment[index]

    def to_dict(self) -> dict[str, Any]:
        x_shear, shear = self.extreme_shear()
        x_moment, moment = self.extreme_moment()
        return {
            "stations": list(self.stations),
            "shear": list(self.shear),
            "moment": list(self.moment),
            "total_weight": self.total_weight,
            "total_buoyancy": self.total_buoyancy,
            "lcb_weight": self.lcb_weight,
            "lcb_buoyancy": self.lcb_buoyancy,
            "residual_force": self.residual_force,
            "residual_moment": self.residual_moment,
            "max_shear": self.max_shear,
            "max_moment": self.max_moment,
            "extreme_shear": [x_shear, shear],
            "extreme_moment": [x_moment, moment],
            "items": [
                {
                    "label": item.label,
                    "x_start": item.x_start,
                    "x_end": item.x_end,
                    "force": item.force,
                    "density": item.density,
                    "centroid": item.centroid,
                }
                for item in self.items
            ],
        }


def _as_stations(stations: Sequence[float], loads: Sequence[float]) -> tuple[tuple[float, ...], tuple[float, ...]]:
    xs = tuple(float(value) for value in stations)
    ys = tuple(float(value) for value in loads)
    if len(xs) != len(ys):
        raise GirderLoadError(
            f"stations and loads must have the same length, got {len(xs)} and {len(ys)}"
        )
    if len(xs) < 2:
        raise GirderLoadError(f"at least two stations are needed, got {len(xs)}")
    for name, values in (("stations", xs), ("loads", ys)):
        for value in values:
            if not math.isfinite(value):
                raise GirderLoadError(f"{name} must be finite, got {value!r}")
    for index in range(len(xs) - 1):
        if not xs[index + 1] > xs[index]:
            raise GirderLoadError(
                f"stations must be strictly increasing, broke at index {index}: "
                f"{xs[index]!r} then {xs[index + 1]!r}"
            )
    return xs, ys


def _interp(xs: tuple[float, ...], ys: tuple[float, ...], x: float) -> float:
    """Linear interpolation on ``xs`` with clamped ends (``xs`` strictly increasing)."""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    low, high = 0, len(xs) - 1
    while high - low > 1:
        middle = (low + high) // 2
        if xs[middle] <= x:
            low = middle
        else:
            high = middle
    span = xs[high] - xs[low]
    fraction = (x - xs[low]) / span
    return ys[low] + (ys[high] - ys[low]) * fraction


def _weight_density(items: Iterable[WeightItem], x: float, side: int) -> float:
    """Load per unit length from the items at ``x``; ``side=-1`` left limit, ``+1`` right."""
    total = 0.0
    for item in items:
        if side < 0:
            inside = item.x_start < x <= item.x_end
        else:
            inside = item.x_start <= x < item.x_end
        if inside:
            total += item.density
    return total


def still_water_loads(
    stations: Sequence[float],
    buoyancy_loads: Sequence[float],
    weight_items: Iterable[WeightItem] = (),
) -> HullGirderLoads:
    """Integrate ``w = p - b`` twice over a piecewise-linear buoyancy diagram.

    ``buoyancy_loads`` is buoyancy per unit length in N/m at ``stations`` (the output of
    :func:`buoyancy_distribution`, or any other source). The domain is
    ``[stations[0], stations[-1]]``; weight items are clipped to it and an item that lies
    entirely outside is an error rather than a silent zero.
    """
    xs, ys = _as_stations(stations, buoyancy_loads)
    items = tuple(weight_items)
    for item in items:
        if not isinstance(item, WeightItem):  # pragma: no cover - defensive
            raise GirderLoadError(f"weight_items must be WeightItem, got {type(item).__name__}")
    x0, x1 = xs[0], xs[-1]
    for item in items:
        if item.x_end <= x0 or item.x_start >= x1:
            raise GirderLoadError(
                f"{item.label or 'weight item'} [{item.x_start}, {item.x_end}] lies outside "
                f"the loaded length [{x0}, {x1}]"
            )

    knots = set(xs)
    for item in items:
        if x0 < item.x_start < x1:
            knots.add(item.x_start)
        if x0 < item.x_end < x1:
            knots.add(item.x_end)
    order = sorted(knots)

    # Per-knot load, right limit first; jumps are handled by taking the left limit at the
    # end of a segment and the right limit at the start of the next one.
    right_limit = [
        _weight_density(items, x, +1) - _interp(xs, ys, x) for x in order
    ]
    shear = [0.0]
    moment = [0.0]
    for index in range(len(order) - 1):
        a, b = order[index], order[index + 1]
        h = b - a
        wa = right_limit[index]
        wb = _weight_density(items, b, -1) - _interp(xs, ys, b)
        q_a = shear[-1]
        m_a = moment[-1]
        shear.append(q_a + h * (wa + wb) / 2.0)
        moment.append(m_a + h * q_a + h * h * (2.0 * wa + wb) / 6.0)

    # Loads: exact integrals of the same diagram (the trapezoid rule is exact for a
    # piecewise-linear integrand, which is why these are not approximations).
    total_buoyancy = 0.0
    first_moment_buoyancy = 0.0
    total_weight = sum(item.force for item in items)
    for index in range(len(xs) - 1):
        a, b = xs[index], xs[index + 1]
        h = b - a
        mean = 0.5 * (ys[index] + ys[index + 1])
        total_buoyancy += h * mean
        first_moment_buoyancy += h * (a * ys[index] + b * ys[index + 1]) / 2.0 + (
            h * h * (ys[index + 1] - ys[index]) / 6.0
        )
    first_moment_weight = sum(item.force * item.centroid for item in items)
    lcb_buoyancy = first_moment_buoyancy / total_buoyancy if total_buoyancy else math.nan
    lcb_weight = first_moment_weight / total_weight if total_weight else math.nan

    return HullGirderLoads(
        stations=tuple(order),
        shear=tuple(shear),
        moment=tuple(moment),
        total_weight=total_weight,
        total_buoyancy=total_buoyancy,
        lcb_weight=lcb_weight,
        lcb_buoyancy=lcb_buoyancy,
        residual_force=shear[-1],
        residual_moment=moment[-1],
        items=items,
    )


def buoyancy_distribution(
    mesh: Mesh,
    *,
    waterline_z_at_reference: float,
    trim_angle_deg: float = 0.0,
    reference_x: float | None = None,
    density: float = SEA_WATER_DENSITY,
    gravity: float = GRAVITY,
    station_count: int = 41,
    stations: Sequence[float] | None = None,
) -> BuoyancyDistribution:
    """Buoyancy per unit length along ``mesh`` at a (possibly trimmed) waterline.

    The waterline is a plane: its height at station ``x`` is
    ``waterline_z_at_reference + tan(trim) * (x - reference_x)``, which is exactly how a
    Bonjean table is read at a prescribed trim -- one waterline per station, then the
    areas are integrated in ``x``. ``trim_angle_deg > 0`` puts the bow down.

    ``reference_x`` defaults to mid-length of the mesh's bounding box. Areas come from
    :func:`tensorcad.bonjean.section_properties`, the same clipped cut used by the rest of
    the package, so a section that is dry at the waterline returns zero area rather than
    an error.
    """
    if not math.isfinite(float(waterline_z_at_reference)):
        raise GirderLoadError(
            f"waterline_z_at_reference must be finite, got {waterline_z_at_reference!r}"
        )
    if not math.isfinite(float(trim_angle_deg)):
        raise GirderLoadError(f"trim_angle_deg must be finite, got {trim_angle_deg!r}")
    for name, value in (("density", density), ("gravity", gravity)):
        if not math.isfinite(float(value)) or float(value) <= 0.0:
            raise GirderLoadError(f"{name} must be finite and positive, got {value!r}")

    box = mesh.bounding_box()
    x_min, x_max = float(box.minimum[0]), float(box.maximum[0])
    keel_z = float(box.minimum[2])
    x_ref = 0.5 * (x_min + x_max) if reference_x is None else float(reference_x)

    if stations is None:
        if station_count < 2:
            raise GirderLoadError(f"station_count must be at least 2, got {station_count!r}")
        stations = tuple(
            x_min + (x_max - x_min) * index / (station_count - 1) for index in range(station_count)
        )
    xs = tuple(sorted(float(value) for value in stations))
    if len(xs) < 2:
        raise GirderLoadError(f"at least two stations are needed, got {len(xs)}")
    for index in range(len(xs) - 1):
        if not xs[index + 1] > xs[index]:
            raise GirderLoadError(f"stations must be strictly increasing (index {index})")

    tan_trim = math.tan(math.radians(float(trim_angle_deg)))
    areas: list[float] = []
    loads: list[float] = []
    for station in xs:
        z = float(waterline_z_at_reference) + tan_trim * (station - x_ref)
        section = section_properties(mesh, station, z)
        area = max(0.0, float(section.area))
        areas.append(area)
        loads.append(float(density) * float(gravity) * area)

    total_area = 0.0
    first_moment_area = 0.0
    for index in range(len(xs) - 1):
        a, b = xs[index], xs[index + 1]
        h = b - a
        total_area += h * 0.5 * (areas[index] + areas[index + 1])
        first_moment_area += h * (a * areas[index] + b * areas[index + 1]) / 2.0 + (
            h * h * (areas[index + 1] - areas[index]) / 6.0
        )

    return BuoyancyDistribution(
        stations=xs,
        areas=tuple(areas),
        loads=tuple(loads),
        volume=total_area,
        lcb=(first_moment_area / total_area) if total_area > 0.0 else math.nan,
        waterline_z=float(waterline_z_at_reference),
        trim_angle_deg=float(trim_angle_deg),
        reference_x=x_ref,
        density=float(density),
        gravity=float(gravity),
        draft_aft=(float(waterline_z_at_reference) + tan_trim * (x_min - x_ref)) - keel_z,
        draft_fore=(float(waterline_z_at_reference) + tan_trim * (x_max - x_ref)) - keel_z,
    )


def hull_girder_loads(
    mesh: Mesh,
    weight_items: Iterable[WeightItem] = (),
    *,
    waterline_z_at_reference: float | None = None,
    trim_angle_deg: float = 0.0,
    reference_x: float | None = None,
    density: float = SEA_WATER_DENSITY,
    gravity: float = GRAVITY,
    station_count: int = 41,
    stations: Sequence[float] | None = None,
) -> tuple[HullGirderLoads, BuoyancyDistribution]:
    """Buoyancy distribution plus the shear force and bending moment it produces.

    ``waterline_z_at_reference`` defaults to the top of the mesh (the hull floating with
    its deck at the waterline), which is only sane for a closed body whose top *is* the
    design waterline -- pass it explicitly otherwise, from
    :func:`tensorcad.trim.waterline_for_displacement` at the design displacement. Both
    the loads and the distribution are returned: the residuals mean nothing without the
    waterline they were computed at.
    """
    if waterline_z_at_reference is None:
        waterline_z_at_reference = float(mesh.bounding_box().maximum[2])
    distribution = buoyancy_distribution(
        mesh,
        waterline_z_at_reference=waterline_z_at_reference,
        trim_angle_deg=trim_angle_deg,
        reference_x=reference_x,
        density=density,
        gravity=gravity,
        station_count=station_count,
        stations=stations,
    )
    loads = still_water_loads(distribution.stations, distribution.loads, weight_items)
    return loads, distribution