"""IMO IS Code 2008 intact-stability criteria evaluated on a measured ``GZ`` curve.

``tensorcad.gz`` and ``tensorcad.righting_arm`` produce :math:`GZ(\\varphi)` from the
surface itself (two independent routes, cross-checked in ``docs/m4038_pr40_merge_cure.md``)
and ``tensorcad.stability`` produces the initial slope :math:`GM_0`; both stop exactly
where the statutory question starts, and say so in as many words --

    "**No weather criterion.** Nothing here computes a wind heeling arm, a `GZ` area
    criterion or an IMO/IS-Code check; the areas `area_under_curve` returns are the input
    to such a check."   -- ``docs/gz.md``

This module is that check. It adds no hydrostatics: it is arithmetic on the polyline
``(angles, gz)`` plus the loading numbers the criteria need (``GM_0``, windage area and
arm, displacement), which is why it stays dependency-free and why every criterion it
reports can be recomputed by hand from the curve it was given.

Clauses implemented -- IS Code 2008 (resolution MSC.267(85)), **part A**:

===========  ======================================================================
2.2.1        area under ``GZ`` to 30 deg                      >= 0.055 m rad
2.2.2        area under ``GZ`` to 40 deg or to the flooding angle if less  >= 0.09 m rad
2.2.3        area between 30 deg and 40 deg (or the flooding angle)  >= 0.03 m rad
2.2.4        ``GZ`` at 30 deg                                 >= 0.20 m
2.2.5        angle of maximum ``GZ``                          >= 25 deg
2.2.6        initial metacentric height ``GM_0``             >= 0.15 m
2.3.1/2.3.2  weather criterion (steady wind heel angle, then gust area ``b >= a``)
2.3.3        windward roll amplitude ``theta_1``             <= ``theta_f``
===========  ======================================================================

**What is not computed, and why.** Clause 2.3.4 obtains the windward roll amplitude
:math:`\\theta_1` from an empirical formula,

    theta1 = 109 k X1 X2 sqrt(r s),

whose factors are read off three tables keyed on ``Ak/B``, ``B/d``, ``Cb`` and the roll
period. Those tables are transcribed nowhere in this repository, and a table copied from
memory is worse than a missing one -- so :math:`\\theta_1` is an **input** here
(:func:`weather_criterion`, ``roll_amplitude_deg``), defaulted to 0 (no windward roll,
which is the case the code's own worked examples start from). That input is what clause
2.3.3 needs: *the angle of heel* ``theta_1`` *shall not exceed the flooding angle*
``theta_f``, so 2.3.3 is evaluated whenever the caller supplies ``flooding_angle_deg`` (and
is then also the ``theta_f`` the 2.3.2 cap uses); when it is not supplied, 2.3.3 is
reported as not evaluated rather than assumed satisfied. Damage stability (part B),
grain and timber-deck codes, and the 2.2 criteria *as amended for* special ship types are
outside the module and are not claimed.

**Interpolation.** Everything below is exact for the piecewise-linear interpolant of the
supplied ``(angle, GZ)`` samples: point values at 30 deg and 40 deg are read off a segment,
areas are trapezoidal on it, and the wind-heel crossings are linear solves. Discretisation
error therefore comes only from the angle grid, is O(h^2) on areas and O(h^2) on a lever
(curvature), and is *measured* -- not assumed -- in ``examples/intact_stability`` (a
5 deg grid moves the barge's 30 deg area by 0.11 %, a 1 deg grid by 0.004 %).
"""

from __future__ import annotations

import bisect
import math
from dataclasses import dataclass
from typing import Any, Iterable, Sequence

__all__ = [
    "IS_CODE_2008_GENERAL",
    "Criterion",
    "IntactStabilityError",
    "IntactStabilityReport",
    "WeatherCriterionReport",
    "WindHeelingLevers",
    "general_criteria",
    "gz_area_between",
    "gz_at_angle",
    "steady_wind_heel_angle",
    "weather_criterion",
    "wind_heeling_levers",
]

SCHEMA_INTACT_STABILITY = "tensorcad.intact-stability/1"

#: IS Code 2008 part A 2.2 thresholds. Units are carried in the tuple so a report can
#: quote them without restating the code.
IS_CODE_2008_GENERAL = {
    "2.2.1": (0.055, "m rad", "area under GZ up to 30 deg"),
    "2.2.2": (0.090, "m rad", "area under GZ up to 40 deg or the flooding angle"),
    "2.2.3": (0.030, "m rad", "area under GZ between 30 deg and 40 deg or the flooding angle"),
    "2.2.4": (0.200, "m", "GZ at 30 deg"),
    "2.2.5": (25.0, "deg", "angle of maximum GZ"),
    "2.2.6": (0.150, "m", "initial metacentric height GM_0"),
}

#: ``2.3.1``: steady wind heel shall not exceed this, or 80 % of deck-edge immersion.
WEATHER_MAX_STEADY_HEEL_DEG = 16.0

#: ``2.3``: wind pressure used for the steady heel, N/m^2 (504 Pa).
WIND_PRESSURE_PA = 504.0

#: ``2.3.2``: gust factor applied to the steady wind heeling lever.
WIND_GUST_FACTOR = 1.5

#: ``2.3.2``: ``theta2`` is capped at 50 deg in addition to down-flooding / second intercept.
WEATHER_THETA2_CAP_DEG = 50.0


class IntactStabilityError(ValueError):
    """Raised when a curve cannot support the check that was asked of it."""


def _round(value: float, digits: int = 9) -> float | None:
    return None if value is None else round(value, digits)


def _finite(value: Any, name: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):  # pragma: no cover - defensive
        raise IntactStabilityError(f"{name} must be a number: {value!r}") from None
    if not math.isfinite(number):
        raise IntactStabilityError(f"{name} must be finite: {value!r}")
    return number


def _polyline(curve: Any) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Accept a curve object (``.angles_deg`` / ``.gz``) or a ``(angles, gz)`` pair."""
    if hasattr(curve, "angles_deg") and hasattr(curve, "gz"):
        # tensorcad.gz.RightingArmCurve
        angles: Iterable[Any] = curve.angles_deg
        values: Iterable[Any] = curve.gz
    elif hasattr(curve, "points"):
        # tensorcad.righting_arm.RightingArmCurve, whose points carry the pair
        try:
            angles = [point.heel_deg for point in curve.points]
            values = [point.righting_arm for point in curve.points]
        except AttributeError:  # pragma: no cover - defensive
            raise IntactStabilityError(
                "a curve with .points must yield objects carrying .heel_deg and "
                f".righting_arm; got {type(curve).__name__}"
            ) from None
    elif isinstance(curve, (tuple, list)) and len(curve) == 2:
        angles, values = curve
        try:
            angles, values = tuple(angles), tuple(values)
        except TypeError:
            raise IntactStabilityError(
                "an (angles, gz) pair must hold two sequences; got "
                f"{type(curve[0]).__name__} and {type(curve[1]).__name__}"
            ) from None
    else:
        raise IntactStabilityError(
            "expected a righting-arm curve exposing .angles_deg and .gz (or .points), or "
            f"an (angles, gz) pair; got {type(curve).__name__}"
        )
    angles = tuple(_finite(angle, "heel angle") for angle in angles)
    values = tuple(_finite(value, "GZ") for value in values)
    if len(angles) != len(values):
        raise IntactStabilityError(
            f"heel angles and GZ values differ in length: {len(angles)} vs {len(values)}"
        )
    if len(angles) < 2:
        raise IntactStabilityError("a curve needs at least two heel angles")
    for left, right in zip(angles, angles[1:]):
        if not right > left:
            raise IntactStabilityError(
                f"heel angles must increase strictly: {left!r} then {right!r}"
            )
    return angles, values


def gz_at_angle(angles_deg: Sequence[float], gz: Sequence[float], angle_deg: float) -> float:
    """``GZ`` at ``angle_deg``: exact on the samples, linear in between, never extrapolated."""
    angles, values = _polyline((angles_deg, gz))
    angle = _finite(angle_deg, "angle_deg")
    if angle < angles[0] or angle > angles[-1]:
        raise IntactStabilityError(
            f"angle {angle!r} deg is outside the sampled range "
            f"[{angles[0]!r}, {angles[-1]!r}] deg"
        )
    index = bisect.bisect_left(angles, angle)
    if index < len(angles) and angles[index] == angle:
        return values[index]
    if index == 0:  # pragma: no cover - guarded by the range check above
        return values[0]
    left, right = index - 1, index
    span = angles[right] - angles[left]
    fraction = (angle - angles[left]) / span
    return values[left] + fraction * (values[right] - values[left])


def gz_area_between(
    angles_deg: Sequence[float],
    gz: Sequence[float],
    lower_deg: float,
    upper_deg: float,
) -> float:
    """Dynamical-stability area ``int GZ dphi`` in m*rad between two heel angles.

    Trapezoidal on the supplied grid, with both limits interpolated, so the four areas the
    code quotes (to 30 deg, to 40 deg, between 30 and 40 deg) are consistent to round-off:
    ``area(0, 40) == area(0, 30) + area(30, 40)``.
    """
    angles, values = _polyline((angles_deg, gz))
    lower = _finite(lower_deg, "lower_deg")
    upper = _finite(upper_deg, "upper_deg")
    if upper <= lower:
        return 0.0
    if lower < angles[0] or upper > angles[-1]:
        raise IntactStabilityError(
            f"the requested interval [{lower}, {upper}] deg is not inside the sampled range "
            f"[{angles[0]}, {angles[-1]}] deg"
        )
    total = 0.0
    for index in range(len(angles) - 1):
        left, right = angles[index], angles[index + 1]
        lo, hi = max(left, lower), min(right, upper)
        if hi <= lo:
            continue
        gz_lo = values[index] if lo == left else gz_at_angle(angles, values, lo)
        gz_hi = values[index + 1] if hi == right else gz_at_angle(angles, values, hi)
        total += 0.5 * (gz_lo + gz_hi) * math.radians(hi - lo)
    return total


@dataclass(frozen=True)
class Criterion:
    """One clause: what it demands, what the curve attained, and whether that passes."""

    clause: str
    quantity: str
    required: float
    attained: float
    unit: str
    comparison: str
    satisfied: bool
    margin: float

    #: ``margin`` is the slack in the direction that satisfies the clause, so it is
    #: non-negative exactly when ``satisfied`` is: ``attained - required`` for a "shall be at
    #: least" clause, ``required - attained`` for a "shall not exceed" one (2.3.1).

    def to_dict(self) -> dict[str, object]:
        return {
            "clause": self.clause,
            "quantity": self.quantity,
            "required": _round(self.required),
            "attained": _round(self.attained),
            "unit": self.unit,
            "comparison": self.comparison,
            "satisfied": self.satisfied,
            "margin": _round(self.margin),
        }


def _clause(clause: str, attained: float) -> Criterion:
    """A 2.2 criterion, with its threshold taken from :data:`IS_CODE_2008_GENERAL`."""
    required, unit, quantity = IS_CODE_2008_GENERAL[clause]
    return _criterion(clause, attained, required, unit, quantity)


def _criterion(clause: str, attained: float, required: float, unit: str, quantity: str) -> Criterion:
    """A "shall be at least" clause: 2.2.1-2.2.6 and 2.3.2."""
    return Criterion(
        clause=clause,
        quantity=quantity,
        required=required,
        attained=attained,
        unit=unit,
        comparison="at least",
        satisfied=attained >= required,
        margin=attained - required,
    )


def _limit(clause: str, attained: float, limit: float, unit: str, quantity: str) -> Criterion:
    """A "shall not exceed" clause: 2.3.1, where a small heel is the passing one."""
    return Criterion(
        clause=clause,
        quantity=quantity,
        required=limit,
        attained=attained,
        unit=unit,
        comparison="not more than",
        satisfied=attained <= limit,
        margin=limit - attained,
    )


@dataclass(frozen=True)
class IntactStabilityReport:
    """The 2.2 general criteria for one curve and one loading condition."""

    criteria: tuple[Criterion, ...]
    gm0_m: float
    flooding_angle_deg: float | None
    upper_angle_deg: float
    area_0_30_m_rad: float
    area_0_40_m_rad: float
    area_30_40_m_rad: float
    gz_30_m: float
    max_gz_m: float
    angle_of_max_gz_deg: float
    angle_step_max_deg: float

    @property
    def satisfied(self) -> bool:
        return all(item.satisfied for item in self.criteria)

    @property
    def failures(self) -> tuple[Criterion, ...]:
        return tuple(item for item in self.criteria if not item.satisfied)

    def verdict(self) -> str:
        if self.satisfied:
            return "passes IS Code 2008 part A 2.2.1-2.2.6"
        failed = ", ".join(item.clause for item in self.failures)
        return f"fails IS Code 2008 part A {failed}"

    def to_dict(self) -> dict[str, object]:
        return {
            "schema": SCHEMA_INTACT_STABILITY,
            "code": "IMO IS Code 2008, part A, 2.2",
            "satisfied": self.satisfied,
            "verdict": self.verdict(),
            "gm0_m": _round(self.gm0_m),
            "flooding_angle_deg": _round(self.flooding_angle_deg),
            "upper_angle_deg": _round(self.upper_angle_deg),
            "area_0_30_m_rad": _round(self.area_0_30_m_rad),
            "area_0_40_m_rad": _round(self.area_0_40_m_rad),
            "area_30_40_m_rad": _round(self.area_30_40_m_rad),
            "gz_30_m": _round(self.gz_30_m),
            "max_gz_m": _round(self.max_gz_m),
            "angle_of_max_gz_deg": _round(self.angle_of_max_gz_deg),
            "angle_step_max_deg": _round(self.angle_step_max_deg),
            "criteria": [item.to_dict() for item in self.criteria],
            "failures": [item.clause for item in self.failures],
        }


def general_criteria(
    curve: Any,
    *,
    gm0_m: float,
    flooding_angle_deg: float | None = None,
    deck_edge_angle_deg: float | None = None,
) -> IntactStabilityReport:
    """Check IS Code 2008 part A 2.2.1-2.2.6 against a measured righting-arm curve.

    ``curve`` is any object with ``angles_deg`` and ``gz`` -- both ``tensorcad.gz`` and
    ``tensorcad.righting_arm`` curves qualify -- or a plain ``(angles, gz)`` pair.
    ``gm0_m`` is the initial metacentric height for the same loading condition, normally
    from :func:`tensorcad.stability.initial_stability` (the loader's ``KG`` is the only
    difference between the two halves of the check, so they must be stated together).
    ``flooding_angle_deg`` is the down-flooding angle: when it is below 40 deg it replaces
    40 deg in 2.2.2 and 2.2.3, as the code says. ``deck_edge_angle_deg`` is not used by
    2.2 (it belongs to 2.3.1) and is accepted only to keep one call site for both checks.
    """
    del deck_edge_angle_deg  # documented above: recorded by the weather report, not here
    angles, values = _polyline(curve)
    gm0 = _finite(gm0_m, "gm0_m")

    flooding = None if flooding_angle_deg is None else _finite(flooding_angle_deg, "flooding_angle_deg")
    if flooding is not None and flooding <= 30.0:
        raise IntactStabilityError(
            f"a flooding angle of {flooding!r} deg cannot support 2.2.2/2.2.3 (needs >= 30 deg)"
        )
    upper = 40.0 if flooding is None or flooding >= 40.0 else flooding
    if angles[-1] < upper:
        raise IntactStabilityError(
            f"the curve stops at {angles[-1]!r} deg but 2.2.2 needs {upper!r} deg: sample "
            "further, or pass the flooding angle if it is lower"
        )
    if angles[-1] < 30.0:
        raise IntactStabilityError(
            f"the curve stops at {angles[-1]!r} deg but 2.2.1 and 2.2.4 need 30 deg"
        )

    area_0_30 = gz_area_between(angles, values, 0.0, 30.0)
    area_0_upper = gz_area_between(angles, values, 0.0, upper)
    area_30_upper = gz_area_between(angles, values, 30.0, upper)
    gz_30 = gz_at_angle(angles, values, 30.0)
    max_gz = max(values)
    peak_index = values.index(max_gz)
    angle_of_max = angles[peak_index]

    criteria = (
        _clause("2.2.1", area_0_30),
        _clause("2.2.2", area_0_upper),
        _clause("2.2.3", area_30_upper),
        _clause("2.2.4", gz_30),
        _clause("2.2.5", angle_of_max),
        _clause("2.2.6", gm0),
    )
    step = max(right - left for left, right in zip(angles, angles[1:]))
    return IntactStabilityReport(
        criteria=criteria,
        gm0_m=gm0,
        flooding_angle_deg=flooding,
        upper_angle_deg=upper,
        area_0_30_m_rad=area_0_30,
        area_0_40_m_rad=area_0_upper,
        area_30_40_m_rad=area_30_upper,
        gz_30_m=gz_30,
        max_gz_m=max_gz,
        angle_of_max_gz_deg=angle_of_max,
        angle_step_max_deg=step,
    )


@dataclass(frozen=True)
class WindHeelingLevers:
    """Steady and gust wind heeling levers in metres, and what produced them."""

    steady_m: float
    gust_m: float
    wind_area_m2: float
    wind_arm_m: float
    displacement_t: float
    pressure_pa: float
    gravity: float

    @property
    def gust_factor(self) -> float:
        return self.gust_m / self.steady_m if self.steady_m else WIND_GUST_FACTOR

    def to_dict(self) -> dict[str, object]:
        return {
            "steady_lever_m": _round(self.steady_m),
            "gust_lever_m": _round(self.gust_m),
            "wind_area_m2": _round(self.wind_area_m2),
            "wind_arm_m": _round(self.wind_arm_m),
            "displacement_t": _round(self.displacement_t),
            "pressure_pa": _round(self.pressure_pa),
            "gravity": _round(self.gravity),
        }


def wind_heeling_levers(
    wind_area_m2: float,
    wind_arm_m: float,
    displacement_t: float,
    *,
    pressure_pa: float = WIND_PRESSURE_PA,
    gust_factor: float = WIND_GUST_FACTOR,
    gravity: float = 9.81,
) -> WindHeelingLevers:
    """``lw1 = P A Z / (1000 g Delta)`` and ``lw2 = 1.5 lw1`` (IS Code 2008, 2.3).

    ``A`` is the projected lateral windage area above the waterline, ``Z`` the vertical
    distance from the centre of ``A`` to the centre of the underwater lateral area, and
    ``Delta`` the displacement in **tonnes** -- the 1000 in the denominator is what carries
    the unit, so the lever comes out in metres.
    """
    area = _finite(wind_area_m2, "wind_area_m2")
    arm = _finite(wind_arm_m, "wind_arm_m")
    displacement = _finite(displacement_t, "displacement_t")
    pressure = _finite(pressure_pa, "pressure_pa")
    factor = _finite(gust_factor, "gust_factor")
    g = _finite(gravity, "gravity")
    if area < 0.0 or arm < 0.0:
        raise IntactStabilityError(f"wind area and arm must be non-negative: {area!r}, {arm!r}")
    if displacement <= 0.0:
        raise IntactStabilityError(f"displacement must be positive: {displacement!r}")
    if pressure <= 0.0 or factor <= 0.0 or g <= 0.0:
        raise IntactStabilityError(
            f"pressure, gust factor and gravity must be positive: {pressure!r}, {factor!r}, {g!r}"
        )
    steady = pressure * area * arm / (1000.0 * g * displacement)
    return WindHeelingLevers(
        steady_m=steady,
        gust_m=factor * steady,
        wind_area_m2=area,
        wind_arm_m=arm,
        displacement_t=displacement,
        pressure_pa=pressure,
        gravity=g,
    )


def _first_crossing(
    angles: Sequence[float], values: Sequence[float], lever: float
) -> float | None:
    """Smallest heel where the polyline crosses ``lever`` upwards, or ``None``."""
    for left, right, value_left, value_right in zip(angles, angles[1:], values, values[1:]):
        if value_left < lever <= value_right:
            span = value_right - value_left
            if span <= 0.0:  # pragma: no cover - strictly increasing guard
                return left
            return left + (right - left) * (lever - value_left) / span
    return None


def _second_crossing(
    angles: Sequence[float], values: Sequence[float], lever: float
) -> float | None:
    """First downward crossing of ``lever`` after the peak: the code's ``theta_c``."""
    peak = values.index(max(values))
    for left, right, value_left, value_right in zip(
        angles[peak:], angles[peak + 1 :], values[peak:], values[peak + 1 :]
    ):
        if value_left >= lever > value_right:
            span = value_left - value_right
            if span <= 0.0:  # pragma: no cover
                return left
            return left + (right - left) * (value_left - lever) / span
    return None


def steady_wind_heel_angle(curve: Any, lever_m: float) -> float | None:
    """The heel where a constant heeling lever balances ``GZ`` (the code's ``theta_0``)."""
    angles, values = _polyline(curve)
    return _first_crossing(angles, values, _finite(lever_m, "lever_m"))


@dataclass(frozen=True)
class WeatherCriterionReport:
    """IS Code 2008 part A 2.3: steady heel (2.3.1), gust areas (2.3.2), roll cap (2.3.3).

    ``area_b_m_rad`` and ``area_a_m_rad`` are ``None`` in the one case where the gust areas
    have no span to be measured over -- a supplied ``flooding_angle_deg`` at or below the
    roll amplitude. There 2.3.3 fails outright and 2.3.2 is not evaluated, instead of a
    number being invented for it.
    """

    levers: WindHeelingLevers
    theta0_deg: float | None
    deck_edge_angle_deg: float | None
    flooding_angle_deg: float | None
    steady_heel_limit_deg: float
    roll_amplitude_deg: float
    theta2_deg: float
    theta_c_deg: float | None
    area_b_m_rad: float | None
    area_a_m_rad: float | None
    criteria: tuple[Criterion, ...]

    @property
    def satisfied(self) -> bool:
        return all(item.satisfied for item in self.criteria)

    @property
    def failures(self) -> tuple[Criterion, ...]:
        return tuple(item for item in self.criteria if not item.satisfied)

    @property
    def clauses(self) -> str:
        """The 2.3 clauses this report actually evaluated, as a readable list."""
        if self.flooding_angle_deg is None:
            return "2.3.1 and 2.3.2"
        if self.area_b_m_rad is None:
            return "2.3.1 and 2.3.3"
        return "2.3.1, 2.3.2 and 2.3.3"

    def verdict(self) -> str:
        if self.satisfied:
            return f"passes IS Code 2008 part A {self.clauses} (weather criterion)"
        failed = ", ".join(item.clause for item in self.failures)
        return f"fails IS Code 2008 part A {failed}"

    def to_dict(self) -> dict[str, object]:
        return {
            "schema": SCHEMA_INTACT_STABILITY,
            "code": "IMO IS Code 2008, part A, 2.3",
            "satisfied": self.satisfied,
            "verdict": self.verdict(),
            "levers": self.levers.to_dict(),
            "theta0_deg": _round(self.theta0_deg),
            "deck_edge_angle_deg": _round(self.deck_edge_angle_deg),
            "flooding_angle_deg": _round(self.flooding_angle_deg),
            "clauses": self.clauses,
            "steady_heel_limit_deg": _round(self.steady_heel_limit_deg),
            "roll_amplitude_deg": _round(self.roll_amplitude_deg),
            "theta2_deg": _round(self.theta2_deg),
            "theta_c_deg": _round(self.theta_c_deg),
            "area_b_m_rad": _round(self.area_b_m_rad),
            "area_a_m_rad": _round(self.area_a_m_rad),
            "criteria": [item.to_dict() for item in self.criteria],
            "failures": [item.clause for item in self.failures],
        }


def weather_criterion(
    curve: Any,
    *,
    wind_area_m2: float,
    wind_arm_m: float,
    displacement_t: float,
    roll_amplitude_deg: float = 0.0,
    deck_edge_angle_deg: float | None = None,
    flooding_angle_deg: float | None = None,
    downflooding_angle_deg: float = WEATHER_THETA2_CAP_DEG,
    gust_factor: float = WIND_GUST_FACTOR,
) -> WeatherCriterionReport:
    """IS Code 2008 part A 2.3.1 (steady heel), 2.3.2 (gust areas ``b >= a``), 2.3.3.

    ``theta_0`` is solved from ``GZ(theta_0) = lw1`` and must not exceed 16 deg, nor 80 % of
    the deck-edge immersion angle when that is known. ``theta_1`` (the windward roll
    amplitude) is **supplied**, not derived -- see the module docstring: clause 2.3.4's
    empirical formula needs roll-period tables that are not transcribed here. ``theta_1 = 0``
    is the case the code's worked example starts from.

    ``theta_2`` is the least of the down-flooding angle, 50 deg, the supplied
    ``flooding_angle_deg`` and the second intercept ``theta_c`` of ``GZ`` with the gust lever
    ``lw2``. The gust area ``a = lw2 (theta_2 - theta_1)`` since the lever is constant, and
    ``b`` is the ``GZ`` area over the same span.

    ``flooding_angle_deg`` is clause 2.3.3's ``theta_f``: the windward roll amplitude
    ``theta_1`` shall not exceed it. Supplying it adds that clause *and* puts it into the
    ``theta_2`` minimum, so one number plays the role the code gives ``theta_f`` in both
    places. When it is omitted, 2.3.3 is not evaluated -- the report says so in ``clauses``
    and leaves ``flooding_angle_deg`` ``None`` rather than presuming a value.
    """
    angles, values = _polyline(curve)
    rotation = _finite(roll_amplitude_deg, "roll_amplitude_deg")
    if rotation < 0.0 or rotation >= 90.0:
        raise IntactStabilityError(f"roll_amplitude_deg must be in [0, 90): {rotation!r}")
    sink = _finite(downflooding_angle_deg, "downflooding_angle_deg")
    if sink <= 0.0:
        raise IntactStabilityError(f"downflooding_angle_deg must be positive: {sink!r}")
    levers = wind_heeling_levers(
        wind_area_m2,
        wind_arm_m,
        displacement_t,
        gust_factor=gust_factor,
    )

    flooding = None if flooding_angle_deg is None else _finite(
        flooding_angle_deg, "flooding_angle_deg"
    )
    if flooding is not None and not 0.0 < flooding < 90.0:
        raise IntactStabilityError(f"flooding_angle_deg must be in (0, 90): {flooding!r}")

    deck_edge = None if deck_edge_angle_deg is None else _finite(
        deck_edge_angle_deg, "deck_edge_angle_deg"
    )
    limit = WEATHER_MAX_STEADY_HEEL_DEG
    if deck_edge is not None:
        limit = min(limit, 0.8 * deck_edge)

    theta0 = _first_crossing(angles, values, levers.steady_m)
    criteria: list[Criterion] = []
    heel_quantity = "steady wind heel angle theta0 (limit min(16 deg, 0.8 deck edge))"
    if theta0 is None:
        # No crossing inside the sampled range: the arm never balances the steady lever, so
        # the body would keep heeling past the last sampled angle. Report the last angle with
        # the clause marked failed rather than inventing an angle the curve does not have.
        criteria.append(
            Criterion(
                clause="2.3.1",
                quantity=heel_quantity,
                required=limit,
                attained=float(angles[-1]),
                unit="deg",
                comparison="not more than",
                satisfied=False,
                margin=limit - float(angles[-1]),
            )
        )
    else:
        criteria.append(_limit("2.3.1", theta0, limit, "deg", heel_quantity))

    theta_c = _second_crossing(angles, values, levers.gust_m)
    caps = [sink, WEATHER_THETA2_CAP_DEG]
    if flooding is not None:
        caps.append(flooding)
    candidates = [*caps, float(angles[-1])]
    if theta_c is not None:
        candidates.append(theta_c)
    theta2 = min(candidates)
    if theta2 <= rotation:
        if flooding is not None and rotation >= flooding:
            # 2.3.3's own failure: the windward roll reaches the flooding angle, so there is
            # no span between theta_1 and theta_2 for the gust areas to be measured over.
            # Report the violated clause and leave the areas unmeasured (None) instead of
            # fabricating a zero or a negative-lever integral.
            criteria.append(
                _limit(
                    "2.3.3",
                    rotation,
                    flooding,
                    "deg",
                    "windward roll amplitude theta1 <= flooding angle theta_f",
                )
            )
            return WeatherCriterionReport(
                levers=levers,
                theta0_deg=theta0,
                deck_edge_angle_deg=deck_edge,
                flooding_angle_deg=flooding,
                steady_heel_limit_deg=limit,
                roll_amplitude_deg=rotation,
                theta2_deg=theta2,
                theta_c_deg=theta_c,
                area_b_m_rad=None,
                area_a_m_rad=None,
                criteria=tuple(criteria),
            )
        raise IntactStabilityError(
            f"theta2 ({theta2!r} deg) does not exceed the roll amplitude ({rotation!r} deg): "
            "the gust areas would be measured over a non-positive span"
        )
    if theta_c is None and theta2 >= float(angles[-1]) - 1e-12 and angles[-1] < min(caps):
        raise IntactStabilityError(
            f"the curve stops at {angles[-1]!r} deg without reaching the gust lever again: "
            "sample further to locate theta2"
        )

    area_b = gz_area_between(angles, values, rotation, theta2)
    area_a = levers.gust_m * math.radians(theta2 - rotation)
    criteria.append(
        _criterion("2.3.2", area_b, area_a, "m rad", "area b under GZ >= area a under the gust lever lw2")
    )
    if flooding is not None:
        criteria.append(
            _limit(
                "2.3.3",
                rotation,
                flooding,
                "deg",
                "windward roll amplitude theta1 <= flooding angle theta_f",
            )
        )
    return WeatherCriterionReport(
        levers=levers,
        theta0_deg=theta0,
        deck_edge_angle_deg=deck_edge,
        flooding_angle_deg=flooding,
        steady_heel_limit_deg=limit,
        roll_amplitude_deg=rotation,
        theta2_deg=theta2,
        theta_c_deg=theta_c,
        area_b_m_rad=area_b,
        area_a_m_rad=area_a,
        criteria=tuple(criteria),
    )