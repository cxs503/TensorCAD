"""IMO IS Code 2008 intact-stability criteria, evaluated on a righting-arm curve.

:mod:`tensorcad.stability` produces the small-angle number ``GM``; :mod:`tensorcad.gz` and
:mod:`tensorcad.righting_arm` produce the curve ``GZ(phi)`` by two independent routes. All
three stop short of a rule check on purpose -- their docstrings say so in as many words
("nothing here should be read as a check against a statutory criterion"). This module is
that check, and only that check: it consumes a curve and returns the verdict, with the
numbers it was derived from.

**The rule being applied.** IMO Resolution MSC.267(85), *International Code on Intact
Stability, 2008*, Part A, chapter 2, section 2.2 -- the general intact-stability criteria:

========  ==============================================================  =============
code      requirement                                                       threshold
========  ==============================================================  =============
2.2.1a    area under ``GZ``, 0 to 30 deg                                       0.055 m rad
2.2.1b    area under ``GZ``, 0 to 40 deg or to ``theta_f`` if that is less     0.090 m rad
2.2.1c    area under ``GZ``, 30 deg to 40 deg (or to ``theta_f``)              0.030 m rad
2.2.2     ``GZ`` at 30 deg                                                      0.20 m
2.2.3     heel at which ``GZ`` reaches its maximum                               25 deg
2.2.4     initial metacentric height ``GM0``                                     0.15 m
========  ==============================================================  =============

``theta_f`` is the flooding angle: the first heel at which an opening that cannot be closed
weathertight is immersed. Passing it in is what makes 2.2.1b and 2.2.1c *the rule* rather
than a 40-degree approximation of it.

**How the integrals are obtained, and why the number comes with an error bar.** A rule check
compares a computed number against a threshold quoted to three decimals, so the quadrature
must not decide the verdict. Two curve sources are accepted:

* a **callable** ``GZ(angle_deg)`` -- integrated with adaptive Simpson to a requested
  tolerance, and the accumulated error estimate is returned alongside the area. This is the
  route used when the curve comes from a converged analytic form or from a mesh solve at
  arbitrary heel angles.
* a **table** of ``(angles, levers)`` (or any object exposing ``angles_deg``/``gz``, which
  the ``RightingArmCurve`` in :mod:`tensorcad.gz` does, or ``points``, which
  :mod:`tensorcad.righting_arm` does) -- interpolated with monotone cubic Hermite
  (Fritsch-Carlson slopes, so no overshoot at the peak) and integrated adaptively; the
  reported uncertainty is the spread between that integral and the plain trapezoidal one on
  the same samples. That spread bounds the *trapezoidal* error, not the Hermite one, and it
  is reported as what it is.

The verdict is three-valued. When the margin between the attained value and the threshold is
smaller than the stated numerical uncertainty the check is reported ``indeterminate`` rather
than ``pass``: a 0.20 m requirement met by 0.2000001 m, read off a 5-degree ladder, is not
something this module signs off.

**What is deliberately not here.** The severe wind and rolling (weather) criterion of section
2.3 needs the windage area, its vertical lever and a gust model that are not geometry this
package holds. The ship-type criteria of sections 2.4-2.5 (timber deck cargo, grain, fishing
vessels, containerships, offshore supply vessels, ...) amend 2.2 rather than replace it.
Damage stability is :mod:`tensorcad.gz`'s documented blind spot as well (no flooding, no
permeability). None of those are implemented, and every assessment says so in
:attr:`IntactStabilityAssessment.scope` instead of implying a complete statutory pass.

Only the Python standard library is used.
"""

from __future__ import annotations

import itertools
import math
import sys
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any

SCHEMA_INTACT_STABILITY = "tensorcad.intact_stability/1"

#: Reference document, quoted in every assessment so a reader can find the rule.
REFERENCE = "IMO Resolution MSC.267(85), IS Code 2008, Part A, ch. 2, s. 2.2"

AREA_TO_30_MIN_M_RAD = 0.055
AREA_TO_40_MIN_M_RAD = 0.090
AREA_30_TO_40_MIN_M_RAD = 0.030
GZ_AT_30_MIN_M = 0.20
ANGLE_OF_MAX_GZ_MIN_DEG = 25.0
GM0_MIN_M = 0.15

#: Default adaptive-Simpson tolerance on an area, in m*rad.
DEFAULT_AREA_TOLERANCE = 1e-11

#: The rules of the code, as data.
INTACT_STABILITY_CRITERIA: Mapping[str, Mapping[str, Any]] = {
    "2.2.1a": {
        "requirement": "area under GZ from 0 to 30 deg",
        "threshold": AREA_TO_30_MIN_M_RAD,
        "unit": "m*rad",
    },
    "2.2.1b": {
        "requirement": "area under GZ from 0 to 40 deg, or to the flooding angle if that is less",
        "threshold": AREA_TO_40_MIN_M_RAD,
        "unit": "m*rad",
    },
    "2.2.1c": {
        "requirement": "area under GZ from 30 deg to 40 deg, or to the flooding angle if that is less",
        "threshold": AREA_30_TO_40_MIN_M_RAD,
        "unit": "m*rad",
    },
    "2.2.2": {
        "requirement": "GZ at 30 deg of heel",
        "threshold": GZ_AT_30_MIN_M,
        "unit": "m",
    },
    "2.2.3": {
        "requirement": "heel at which GZ reaches its maximum",
        "threshold": ANGLE_OF_MAX_GZ_MIN_DEG,
        "unit": "deg",
    },
    "2.2.4": {
        "requirement": "initial metacentric height GM0",
        "threshold": GM0_MIN_M,
        "unit": "m",
    },
}


class StatutoryError(ValueError):
    """Raised when a curve cannot be read as a GZ curve at all."""


@dataclass(frozen=True)
class CriterionCheck:
    """One rule, the number it was measured on, and the verdict.

    ``status`` is ``"pass"``, ``"fail"``, ``"indeterminate"`` (the margin sits inside the
    stated numerical uncertainty) or ``"not_assessed"`` (the curve does not reach the angle
    the rule is about, so the honest answer is no answer).
    """

    code: str
    requirement: str
    threshold: float
    unit: str
    attained: float | None
    status: str
    uncertainty: float = 0.0
    note: str = ""

    @property
    def margin(self) -> float | None:
        """``attained - threshold``; positive is compliant. ``None`` when not assessed."""
        if self.attained is None:
            return None
        return self.attained - self.threshold

    @property
    def passed(self) -> bool:
        return self.status == "pass"

    def to_dict(self) -> dict[str, Any]:
        return {
            "code": self.code,
            "requirement": self.requirement,
            "threshold": self.threshold,
            "unit": self.unit,
            "attained": None if self.attained is None else _round(self.attained, 9),
            "margin": None if self.margin is None else _round(self.margin, 9),
            "uncertainty": _round(self.uncertainty, 12),
            "status": self.status,
            "note": self.note,
        }


@dataclass(frozen=True)
class IntactStabilityAssessment:
    """Verdict of section 2.2 on one curve, with every number that produced it."""

    checks: tuple[CriterionCheck, ...]
    reference: str = REFERENCE
    scope: str = (
        "IS Code 2008 Part A s. 2.2 general criteria only; the weather criterion (s. 2.3), "
        "the ship-type criteria (s. 2.4-2.5) and damage stability are not implemented"
    )
    flooding_angle_deg: float | None = None
    gm0_m: float | None = None
    gm0_source: str = "not assessed"
    angle_of_max_gz_deg: float | None = None
    max_gz_m: float | None = None
    vanishing_angle_deg: float | None = None
    curve_upper_deg: float | None = None
    curve_kind: str = "unknown"
    area_to_30_m_rad: float | None = None
    area_to_limit_m_rad: float | None = None
    area_30_to_limit_m_rad: float | None = None
    area_limit_deg: float | None = None
    area_uncertainty_m_rad: float = 0.0

    def check(self, code: str) -> CriterionCheck:
        for item in self.checks:
            if item.code == code:
                return item
        raise StatutoryError(f"no criterion {code!r} in this assessment")

    @property
    def failed(self) -> tuple[CriterionCheck, ...]:
        return tuple(item for item in self.checks if item.status == "fail")

    @property
    def indeterminate(self) -> tuple[CriterionCheck, ...]:
        return tuple(item for item in self.checks if item.status == "indeterminate")

    @property
    def not_assessed(self) -> tuple[CriterionCheck, ...]:
        return tuple(item for item in self.checks if item.status == "not_assessed")

    @property
    def passed(self) -> bool | None:
        """``True`` if every rule passed, ``False`` if any failed, ``None`` otherwise.

        ``None`` is a third answer, not a soft ``False``: it means the curve could not
        answer one of the rules (or a margin sat inside the numerical noise), so the caller
        must widen the curve or tighten the solve before a rule check exists at all.
        """
        if self.failed:
            return False
        if self.not_assessed or self.indeterminate:
            return None
        return True

    def verdict(self) -> str:
        """``"pass"``, ``"fail"`` or ``"undetermined"`` -- the string form of :attr:`passed`."""
        passed = self.passed
        if passed is None:
            return "undetermined"
        return "pass" if passed else "fail"

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_INTACT_STABILITY,
            "reference": self.reference,
            "scope": self.scope,
            "verdict": self.verdict(),
            "checks": [item.to_dict() for item in self.checks],
            "flooding_angle_deg": _round(self.flooding_angle_deg, 9),
            "area_limit_deg": _round(self.area_limit_deg, 9),
            "area_to_30_m_rad": _round(self.area_to_30_m_rad, 12),
            "area_to_limit_m_rad": _round(self.area_to_limit_m_rad, 12),
            "area_30_to_limit_m_rad": _round(self.area_30_to_limit_m_rad, 12),
            "area_uncertainty_m_rad": _round(self.area_uncertainty_m_rad, 12),
            "gm0_m": _round(self.gm0_m, 12),
            "gm0_source": self.gm0_source,
            "angle_of_max_gz_deg": _round(self.angle_of_max_gz_deg, 9),
            "max_gz_m": _round(self.max_gz_m, 12),
            "vanishing_angle_deg": _round(self.vanishing_angle_deg, 9),
            "curve_upper_deg": _round(self.curve_upper_deg, 9),
            "curve_kind": self.curve_kind,
        }


def _round(value: float | None, digits: int = 12) -> float | None:
    if value is None:
        return None
    if not math.isfinite(value):
        return value
    return round(value, digits)


# --------------------------------------------------------------------------------------
# quadrature
# --------------------------------------------------------------------------------------


def _adaptive_simpson(
    function: Callable[[float], float],
    lower: float,
    upper: float,
    tolerance: float,
) -> tuple[float, float]:
    """Adaptive Simpson on ``[lower, upper]``: ``(value, estimated_absolute_error)``."""
    middle = 0.5 * (lower + upper)
    f_lower, f_middle, f_upper = function(lower), function(middle), function(upper)
    whole = (upper - lower) / 6.0 * (f_lower + 4.0 * f_middle + f_upper)
    return _simpson_refine(
        function, lower, upper, f_lower, f_middle, f_upper, whole, tolerance, 40
    )


def _simpson_refine(
    function: Callable[[float], float],
    lower: float,
    upper: float,
    f_lower: float,
    f_middle: float,
    f_upper: float,
    whole: float,
    tolerance: float,
    depth: int,
) -> tuple[float, float]:
    middle = 0.5 * (lower + upper)
    left_middle = 0.5 * (lower + middle)
    right_middle = 0.5 * (middle + upper)
    f_left_middle = function(left_middle)
    f_right_middle = function(right_middle)
    left = (middle - lower) / 6.0 * (f_lower + 4.0 * f_left_middle + f_middle)
    right = (upper - middle) / 6.0 * (f_middle + 4.0 * f_right_middle + f_upper)
    delta = left + right - whole
    if depth <= 0 or abs(delta) <= 15.0 * tolerance:
        return left + right + delta / 15.0, abs(delta) / 15.0
    left_value, left_error = _simpson_refine(
        function, lower, middle, f_lower, f_left_middle, f_middle, left, tolerance / 2.0, depth - 1
    )
    right_value, right_error = _simpson_refine(
        function, middle, upper, f_middle, f_right_middle, f_upper, right, tolerance / 2.0, depth - 1
    )
    return left_value + right_value, left_error + right_error


def _pchip_slopes(angles: Sequence[float], levers: Sequence[float]) -> list[float]:
    """Fritsch-Carlson monotone slopes for cubic Hermite interpolation.

    A plain cubic spline overshoots next to the peak of a GZ curve, which biases exactly the
    numbers a rule check reads (``GZ_max`` and the area either side of it). Limiting the
    slopes keeps the interpolant inside the local data range.
    """
    count = len(angles)
    if count == 2:
        slope = (levers[1] - levers[0]) / (angles[1] - angles[0])
        return [slope, slope]
    widths = [angles[index + 1] - angles[index] for index in range(count - 1)]
    secants = [(levers[index + 1] - levers[index]) / widths[index] for index in range(count - 1)]
    slopes = [0.0] * count
    slopes[0] = secants[0]
    slopes[-1] = secants[-1]
    for index in range(1, count - 1):
        if secants[index - 1] * secants[index] <= 0.0:
            slopes[index] = 0.0
        else:
            weight_left = 2.0 * widths[index] + widths[index - 1]
            weight_right = widths[index] + 2.0 * widths[index - 1]
            slopes[index] = (weight_left + weight_right) / (
                weight_left / secants[index - 1] + weight_right / secants[index]
            )
    return slopes


class _GzSource:
    """A GZ curve behind one interface -- ``at``, ``area``, ``lower_deg``, ``upper_deg``."""

    kind: str

    def at(self, angle_deg: float) -> float:  # pragma: no cover - interface
        raise NotImplementedError

    def area(self, start_deg: float, stop_deg: float) -> tuple[float, float]:  # pragma: no cover
        raise NotImplementedError

    def interpolation_spread(self, angle_deg: float) -> float:
        """Uncertainty of a single ordinate: zero for a curve that is its own interpolant."""
        return 0.0

    def refinement_error(self) -> float:
        """Self-convergence error of :func:`_gm0_from_curve` on this source."""
        return 0.0

    def small_angle_pair(self) -> tuple[float, float]:
        """Two heel angles ``(h, r*h)`` near upright used to read the slope there."""
        return 1.0, 2.0

    @property
    def lower_deg(self) -> float:  # pragma: no cover - interface
        raise NotImplementedError

    @property
    def upper_deg(self) -> float:  # pragma: no cover - interface
        raise NotImplementedError


class _CallableSource(_GzSource):
    """A GZ curve given as ``f(angle_deg) -> lever``."""

    kind = "callable"

    def __init__(
        self,
        function: Callable[[float], float],
        domain: tuple[float, float],
        tolerance: float,
    ) -> None:
        self._function = function
        self._lower, self._upper = float(domain[0]), float(domain[1])
        self._tolerance = tolerance
        if not (math.isfinite(self._lower) and math.isfinite(self._upper)):
            raise StatutoryError(f"heel domain must be finite: {domain!r}")
        if self._upper <= self._lower:
            raise StatutoryError(f"empty heel domain: {domain!r}")

    def at(self, angle_deg: float) -> float:
        value = float(self._function(angle_deg))
        if not math.isfinite(value):
            raise StatutoryError(f"GZ({angle_deg!r} deg) is not finite: {value!r}")
        return value

    @property
    def lower_deg(self) -> float:
        return self._lower

    @property
    def upper_deg(self) -> float:
        return self._upper

    def area(self, start_deg: float, stop_deg: float) -> tuple[float, float]:
        radians_per_degree = math.pi / 180.0
        value, error = _adaptive_simpson(self.at, start_deg, stop_deg, self._tolerance)
        return value * radians_per_degree, error * radians_per_degree


class _TableSource(_GzSource):
    """A GZ curve given as sampled ``(angles, levers)``, read through a monotone spline."""

    kind = "table"

    def __init__(self, angles: Sequence[float], levers: Sequence[float], tolerance: float) -> None:
        heel_angles = tuple(float(angle) for angle in angles)
        arms = tuple(float(lever) for lever in levers)
        if len(heel_angles) != len(arms):
            raise StatutoryError(f"{len(heel_angles)} angles against {len(arms)} levers")
        if len(heel_angles) < 2:
            raise StatutoryError("a GZ table needs at least two samples")
        for angle, lever in zip(heel_angles, arms, strict=True):
            if not math.isfinite(angle) or not math.isfinite(lever):
                raise StatutoryError(f"non-finite sample: ({angle!r}, {lever!r})")
        for previous, current in itertools.pairwise(heel_angles):
            if current <= previous:
                raise StatutoryError(
                    f"heel angles must be strictly increasing: {previous!r} then {current!r}"
                )
        self._angles = heel_angles
        self._levers = arms
        self._slopes = _pchip_slopes(heel_angles, arms)
        self._tolerance = tolerance

    @property
    def lower_deg(self) -> float:
        return self._angles[0]

    @property
    def upper_deg(self) -> float:
        return self._angles[-1]

    def at(self, angle_deg: float) -> float:
        angle = float(angle_deg)
        if angle < self._angles[0] or angle > self._angles[-1]:
            raise StatutoryError(
                f"heel {angle!r} deg is outside the sampled curve "
                f"[{self._angles[0]!r}, {self._angles[-1]!r}]"
            )
        angles, levers, slopes = self._angles, self._levers, self._slopes
        low, high = 0, len(angles) - 1
        while high - low > 1:
            middle = (low + high) // 2
            if angles[middle] <= angle:
                low = middle
            else:
                high = middle
        width = angles[high] - angles[low]
        if width <= 0.0:  # pragma: no cover - guarded by the strictly-increasing check
            return levers[low]
        t = (angle - angles[low]) / width
        t2, t3 = t * t, t * t * t
        h00 = 2.0 * t3 - 3.0 * t2 + 1.0
        h10 = t3 - 2.0 * t2 + t
        h01 = -2.0 * t3 + 3.0 * t2
        h11 = t3 - t2
        return (
            h00 * levers[low]
            + h10 * width * slopes[low]
            + h01 * levers[high]
            + h11 * width * slopes[high]
        )

    def _linear_at(self, angle_deg: float) -> float:
        """Ordinate from plain linear interpolation, the alternative reading of the table."""
        angles, levers = self._angles, self._levers
        for index in range(len(angles) - 1):
            if angles[index] <= angle_deg <= angles[index + 1]:
                width = angles[index + 1] - angles[index]
                if width <= 0.0:  # pragma: no cover - strictly increasing by construction
                    return levers[index]
                weight = (angle_deg - angles[index]) / width
                return levers[index] + weight * (levers[index + 1] - levers[index])
        return levers[-1]

    def interpolation_spread(self, angle_deg: float) -> float:
        """How far the monotone spline and linear interpolation differ at one ordinate.

        On a sampled curve the spline is the better reading but it is still an
        interpolation, and this spread is what that ambiguity is worth at this angle.
        """
        return abs(self.at(angle_deg) - self._linear_at(angle_deg))

    def small_angle_pair(self) -> tuple[float, float]:
        """The first two sampled heels above upright.

        A table is only known at its own samples, so the slope at upright is read there
        rather than from the interpolant: interpolating between 0 and 10 degrees and then
        differencing at 1 and 2 degrees measures the interpolant, not the hull.
        """
        positive = [angle for angle in self._angles if angle > 0.0]
        if len(positive) < 2:
            raise StatutoryError("a table needs two heels above upright to read GM0 from its slope")
        return positive[0], positive[1]

    def refinement_error(self) -> float:
        """``GM0`` read off the full table against ``GM0`` read off every second sample.

        A table has no interior knob to tighten, so the only honest error estimate is a
        self-convergence one: halve the sampling and see how much the slope at upright moves.
        """
        angles, levers = self._angles, self._levers
        if len(angles) < 5:
            return 0.0
        coarse_angles = list(angles[::2])
        coarse_levers = list(levers[::2])
        if coarse_angles[-1] != angles[-1]:
            coarse_angles.append(angles[-1])
            coarse_levers.append(levers[-1])
        if len(coarse_angles) < 3:
            return 0.0
        try:
            fine_value, _ = _gm0_from_curve(self)
            coarse_value, _ = _gm0_from_curve(_TableSource(coarse_angles, coarse_levers, self._tolerance))
        except StatutoryError:  # pragma: no cover - coarsening keeps the 0..2 deg span
            return 0.0
        return abs(fine_value - coarse_value)

    def _trapezoid(self, start_deg: float, stop_deg: float) -> float:
        """Plain trapezoidal area in m*rad, truncating the end intervals at the limits."""
        total = 0.0
        radians_per_degree = math.pi / 180.0
        angles = self._angles
        for index in range(len(angles) - 1):
            left, right = angles[index], angles[index + 1]
            if right <= start_deg or left >= stop_deg:
                continue
            start = max(left, start_deg)
            stop = min(right, stop_deg)
            if stop <= start:
                continue
            total += 0.5 * (self.at(start) + self.at(stop)) * (stop - start) * radians_per_degree
        return total

    def area(self, start_deg: float, stop_deg: float) -> tuple[float, float]:
        radians_per_degree = math.pi / 180.0
        value, error = _adaptive_simpson(self.at, start_deg, stop_deg, self._tolerance)
        spline = value * radians_per_degree
        spread = abs(spline - self._trapezoid(start_deg, stop_deg))
        return spline, spread + error * radians_per_degree


def _as_source(
    curve: Any,
    *,
    heel_domain_deg: tuple[float, float] | None,
    tolerance: float,
) -> _GzSource:
    """Coerce the accepted curve spellings into one :class:`_GzSource`."""
    if isinstance(curve, _GzSource):  # pragma: no cover - internal reuse
        return curve
    if callable(curve) and not isinstance(curve, (Sequence, Mapping)):
        domain = (0.0, 90.0) if heel_domain_deg is None else heel_domain_deg
        return _CallableSource(curve, domain, tolerance)
    angles = getattr(curve, "angles_deg", None)
    levers = getattr(curve, "gz", None)
    if angles is not None and levers is not None:
        return _TableSource(angles, levers, tolerance)
    points = getattr(curve, "points", None)
    if points is not None:
        try:
            heel_angles = [float(point.heel_deg) for point in points]
            arms = [float(point.righting_arm) for point in points]
        except (AttributeError, TypeError) as error:  # pragma: no cover - exotic input
            raise StatutoryError(f"cannot read points as a GZ curve: {error}") from error
        return _TableSource(heel_angles, arms, tolerance)
    if isinstance(curve, (tuple, list)) and len(curve) == 2 and not isinstance(curve[0], (int, float)):
        heel_angles, arms = curve
        return _TableSource(list(heel_angles), list(arms), tolerance)
    raise StatutoryError(
        "a GZ curve is a callable f(angle_deg), an (angles, levers) pair, or an object with "
        f"angles_deg/gz or points; got {type(curve).__name__}"
    )


# --------------------------------------------------------------------------------------
# curve readings a rule needs
# --------------------------------------------------------------------------------------


def _peak(source: _GzSource, ceiling_deg: float | None = None) -> tuple[float, float]:
    """``(angle_deg, gz_m)`` at the maximum of the curve, golden-section on its domain.

    The GZ curve of an intact hull is unimodal up to its vanishing angle, which is the
    assumption a golden-section search needs and the reason this is not a scan: a
    5-degree ladder can miss a peak sitting between two samples by whole degrees, and
    2.2.3 is a criterion on the *angle*.

    ``ceiling_deg`` truncates the search at the angle past which the ship cannot be heeled
    intact (the flooding or vanishing angle); a peak found at that ceiling is a peak "so
    far", not the peak, and the caller reports it as such.
    """
    lower = max(source.lower_deg, 0.0)
    upper = source.upper_deg if ceiling_deg is None else min(source.upper_deg, ceiling_deg)
    if upper <= lower:  # pragma: no cover - guarded by the domain checks
        raise StatutoryError(f"empty peak-search range [{lower!r}, {upper!r}] deg")
    if upper <= lower:  # pragma: no cover - guarded by the domain checks
        raise StatutoryError("cannot find a peak on an empty heel range")
    golden = (math.sqrt(5.0) - 1.0) / 2.0
    left = upper - golden * (upper - lower)
    right = lower + golden * (upper - lower)
    f_left, f_right = source.at(left), source.at(right)
    for _ in range(200):
        if upper - lower <= 1e-10:
            break
        if f_left >= f_right:
            upper, right, f_right = right, left, f_left
            left = upper - golden * (upper - lower)
            f_left = source.at(left)
        else:
            lower, left, f_left = left, right, f_right
            right = lower + golden * (upper - lower)
            f_right = source.at(right)
    angle = 0.5 * (lower + upper)
    # The end of the domain can hold the true maximum when the curve is still rising there;
    # report that end rather than an interior point that is not the maximum.
    for end in (source.lower_deg, source.upper_deg):
        end_value = source.at(end)
        if end_value > source.at(angle):
            angle = end
    return angle, source.at(angle)


def _vanishing_angle(source: _GzSource) -> float | None:
    """First heel past the peak where ``GZ`` returns to zero, or ``None`` if it does not."""
    lower = max(source.lower_deg, 0.0)
    upper = source.upper_deg
    if upper <= lower:
        return None
    samples = max(round((upper - lower) / 0.5), 2)
    previous_angle = lower
    previous_value = source.at(previous_angle)
    for index in range(1, samples + 1):
        angle = lower + (upper - lower) * index / samples
        value = source.at(angle)
        if value <= 0.0 and previous_value > 0.0 and index > 1:
            low, high = previous_angle, angle
            for _ in range(60):
                middle = 0.5 * (low + high)
                if source.at(middle) > 0.0:
                    low = middle
                else:
                    high = middle
            return 0.5 * (low + high)
        previous_angle, previous_value = angle, value
    return None


def _gm0_from_curve(source: _GzSource) -> tuple[float, float]:
    """``GM0`` from the slope of ``GZ`` at upright, Richardson-extrapolated: ``(value, error)``.

    For a symmetric hull form ``GZ(phi) = GM0 sin(phi) + c sin^3(phi) + ...``, so
    ``GZ(phi)/sin(phi) = GM0 + c phi^2 + O(phi^4)`` and two heel angles an octave apart
    cancel the ``phi^2`` term exactly. The residual is reported as the uncertainty rather
    than hidden: this is an estimate read off the curve, not an independent calculation, and
    :func:`tensorcad.stability.initial_stability` is the better source when it exists.
    """
    if source.lower_deg > 0.0:
        raise StatutoryError(
            f"the curve starts at {source.lower_deg!r} deg, so GM0 cannot be read from its "
            "slope at upright; pass gm0 explicitly"
        )
    fine, coarse = source.small_angle_pair()
    if coarse > source.upper_deg:
        raise StatutoryError(
            f"the curve spans [{source.lower_deg!r}, {source.upper_deg!r}] deg, so GM0 cannot "
            f"be read from its slope at upright (needs 0 to {coarse!r} deg); pass gm0 explicitly"
        )

    def ratio(angle: float) -> float:
        sine = math.sin(math.radians(angle))
        return source.at(angle) / sine if sine > 0.0 else source.at(angle) / math.radians(angle)

    fine_ratio = ratio(fine)
    coarse_ratio = ratio(coarse)
    step = coarse / fine
    # Richardson for an O(h^2) leading error: f(h) + (f(h) - f(rh)) / (r^2 - 1).
    value = fine_ratio + (fine_ratio - coarse_ratio) / (step * step - 1.0)
    return value, abs(value - fine_ratio)


# --------------------------------------------------------------------------------------
# public API
# --------------------------------------------------------------------------------------


def area_under_gz_curve(
    curve: Any,
    start_deg: float = 0.0,
    stop_deg: float = 30.0,
    *,
    heel_domain_deg: tuple[float, float] | None = None,
    tolerance: float = DEFAULT_AREA_TOLERANCE,
) -> float:
    """``int GZ dphi`` in m*rad between two heel angles, for any accepted curve spelling."""
    if not math.isfinite(start_deg) or not math.isfinite(stop_deg):
        raise StatutoryError(f"heel limits must be finite: {start_deg!r}, {stop_deg!r}")
    if stop_deg <= start_deg:
        raise StatutoryError(f"empty heel interval: {start_deg!r} to {stop_deg!r}")
    source = _as_source(curve, heel_domain_deg=heel_domain_deg, tolerance=tolerance)
    if start_deg < source.lower_deg or stop_deg > source.upper_deg:
        raise StatutoryError(
            f"heel interval [{start_deg!r}, {stop_deg!r}] is outside the curve "
            f"[{source.lower_deg!r}, {source.upper_deg!r}]"
        )
    return source.area(start_deg, stop_deg)[0]


def gm0_from_gz_curve(
    curve: Any,
    *,
    heel_domain_deg: tuple[float, float] | None = None,
    tolerance: float = DEFAULT_AREA_TOLERANCE,
) -> float:
    """``GM0`` read from the slope of the curve at upright (Richardson on 1 deg and 2 deg)."""
    source = _as_source(curve, heel_domain_deg=heel_domain_deg, tolerance=tolerance)
    return _gm0_from_curve(source)[0]


def _check(
    code: str,
    attained: float | None,
    uncertainty: float,
    *,
    note: str = "",
) -> CriterionCheck:
    spec = INTACT_STABILITY_CRITERIA[code]
    threshold = float(spec["threshold"])
    if attained is None:
        return CriterionCheck(
            code=code,
            requirement=str(spec["requirement"]),
            threshold=threshold,
            unit=str(spec["unit"]),
            attained=None,
            status="not_assessed",
            note=note,
        )
    margin = attained - threshold
    # The rules read "at least", so a value equal to the threshold passes; binary round-off
    # at the 1e-15 level must not turn that into a fail. Real numerical uncertainty -- a
    # table's interpolation spread, a quadrature error -- is a different animal and any
    # margin inside it is reported as undecidable rather than rounded to a verdict.
    floor = 8.0 * sys.float_info.epsilon * max(1.0, abs(attained))
    if uncertainty > floor and abs(margin) <= uncertainty:
        status = "indeterminate"
        note = (note + "; margin inside the numerical uncertainty").strip("; ")
    elif margin + floor >= 0.0:
        status = "pass"
    else:
        status = "fail"
    return CriterionCheck(
        code=code,
        requirement=str(spec["requirement"]),
        threshold=threshold,
        unit=str(spec["unit"]),
        attained=attained,
        status=status,
        uncertainty=uncertainty,
        note=note,
    )


def check_intact_stability(
    curve: Any,
    *,
    gm0: float | None = None,
    flooding_angle_deg: float | None = None,
    heel_domain_deg: tuple[float, float] | None = None,
    tolerance: float = DEFAULT_AREA_TOLERANCE,
) -> IntactStabilityAssessment:
    """Evaluate the IS Code 2008 s. 2.2 general criteria on a GZ curve.

    Parameters
    ----------
    curve:
        The righting-arm curve: a callable ``f(angle_deg) -> GZ``, an
        ``(angles_deg, levers)`` pair, a :class:`tensorcad.gz.RightingArmCurve`, or a
        :class:`tensorcad.righting_arm.RightingArmCurve`.
    gm0:
        Initial metacentric height in metres. When omitted it is read off the curve's slope
        at upright (see :func:`gm0_from_gz_curve`) and ``gm0_source`` says so; when supplied
        from :func:`tensorcad.stability.initial_stability` it is the better number.
    flooding_angle_deg:
        ``theta_f``, the first heel at which a non-weathertight opening is immersed. Without
        it, 2.2.1b and 2.2.1c are evaluated to 40 degrees and the check notes that ``theta_f``
        was not supplied.
    heel_domain_deg:
        Domain of a *callable* curve, default ``(0, 90)``; ignored for tabulated curves,
        whose own sample range is the domain.
    tolerance:
        Adaptive-Simpson tolerance on an area, in m*rad.

    Returns
    -------
    IntactStabilityAssessment
        One :class:`CriterionCheck` per rule of the code, plus the areas, the peak, the
        vanishing angle and the verdict. ``passed`` is ``None`` -- not ``False`` -- when the
        curve cannot answer a rule or a margin sits inside the numerical uncertainty.
    """
    if gm0 is not None and (not math.isfinite(gm0) or gm0 <= 0.0):
        raise StatutoryError(f"gm0 must be finite and positive: {gm0!r}")
    if flooding_angle_deg is not None and (
        not math.isfinite(flooding_angle_deg) or flooding_angle_deg <= 0.0
    ):
        raise StatutoryError(
            f"flooding_angle_deg must be finite and positive: {flooding_angle_deg!r}"
        )
    source = _as_source(curve, heel_domain_deg=heel_domain_deg, tolerance=tolerance)
    upper = source.upper_deg
    if upper < 2.0:
        raise StatutoryError(
            f"a stability check needs a curve reaching at least 2 deg of heel; got {upper!r}"
        )
    notes: list[str] = []
    if flooding_angle_deg is None:
        notes.append("flooding angle not supplied: 40 deg used for the area limits")

    limit = 40.0 if flooding_angle_deg is None else min(40.0, flooding_angle_deg)
    vanishing = _vanishing_angle(source)

    def blocked(needed: float) -> str | None:
        """Why the curve cannot answer a rule about ``needed`` degrees, or ``None``."""
        if flooding_angle_deg is not None and flooding_angle_deg < needed:
            return (
                f"the flooding angle ({flooding_angle_deg:.3f} deg) is below the {needed:.3f} deg "
                "the rule needs: the ship cannot be heeled that far intact"
            )
        if vanishing is not None and vanishing < needed:
            return f"GZ vanishes at {vanishing:.3f} deg, before the {needed:.3f} deg the rule needs"
        if upper < needed:
            return f"the curve stops at {upper:.3f} deg, before the {needed:.3f} deg the rule needs"
        return None

    limit_block = blocked(limit)
    block_30 = blocked(30.0)

    # The peak is read over the range the ship can actually be heeled through, not over the
    # whole sampled curve: a maximum on a stretch past flooding is not a righting arm.
    ceiling = upper
    if flooding_angle_deg is not None:
        ceiling = min(ceiling, flooding_angle_deg)
    if vanishing is not None:
        ceiling = min(ceiling, vanishing)
    peak_angle, peak_value = _peak(source, ceiling)
    at_edge = peak_angle >= ceiling - 1e-9
    peak_block = blocked(ANGLE_OF_MAX_GZ_MIN_DEG)
    if peak_block is not None and at_edge:
        peak_attained: float | None = None
        peak_note = f"{peak_block} (the maximum on the assessable range sits at its edge)"
    else:
        peak_attained = peak_angle
        peak_note = (
            "peak at the edge of the assessable range: the true maximum may lie beyond it"
            if at_edge
            else ""
        )

    # Areas.
    area_30: float | None = None
    area_30_uncertainty = 0.0
    if block_30 is None:
        area_30, area_30_uncertainty = source.area(0.0, 30.0)
    area_limit: float | None = None
    area_limit_uncertainty = 0.0
    if limit_block is None:
        area_limit, area_limit_uncertainty = source.area(0.0, limit)
    area_30_limit: float | None = None
    area_30_limit_uncertainty = 0.0
    if block_30 is None and limit_block is None:
        area_30_limit, area_30_limit_uncertainty = source.area(30.0, limit)

    gz_30: float | None = None
    if block_30 is None:
        gz_30 = source.at(30.0)

    # GM0: supplied, or read off the curve.
    gm0_value: float | None = None
    gm0_uncertainty = 0.0
    gm0_source = "not assessed"
    if gm0 is not None:
        gm0_value = float(gm0)
        gm0_source = "supplied"
    else:
        try:
            gm0_value, gm0_uncertainty = _gm0_from_curve(source)
            gm0_uncertainty += source.refinement_error()
            fine, coarse = source.small_angle_pair()
            gm0_source = (
                "slope of the supplied curve at upright "
                f"(Richardson, {fine:g} deg / {coarse:g} deg)"
            )
        except StatutoryError as error:
            notes.append(f"GM0 not assessed: {error}")

    checks = (
        _check("2.2.1a", area_30, area_30_uncertainty, note=block_30 or ""),
        _check(
            "2.2.1b",
            area_limit,
            area_limit_uncertainty,
            note="; ".join(
                part
                for part in (
                    limit_block or "",
                    notes[0] if flooding_angle_deg is None else "",
                )
                if part
            ),
        ),
        _check(
            "2.2.1c",
            area_30_limit,
            area_30_limit_uncertainty,
            note=block_30 or limit_block or "",
        ),
        _check(
            "2.2.2",
            gz_30,
            0.0 if gz_30 is None else source.interpolation_spread(30.0),
            note=block_30 or "",
        ),
        _check("2.2.3", peak_attained, 0.0, note=peak_note),
        _check("2.2.4", gm0_value, gm0_uncertainty),
    )
    return IntactStabilityAssessment(
        checks=checks,
        flooding_angle_deg=flooding_angle_deg,
        gm0_m=gm0_value,
        gm0_source=gm0_source,
        angle_of_max_gz_deg=peak_angle,
        max_gz_m=peak_value,
        vanishing_angle_deg=vanishing,
        curve_upper_deg=upper,
        curve_kind=source.kind,
        area_to_30_m_rad=area_30,
        area_to_limit_m_rad=area_limit,
        area_30_to_limit_m_rad=area_30_limit,
        area_limit_deg=limit,
        area_uncertainty_m_rad=max(area_30_uncertainty, area_limit_uncertainty, area_30_limit_uncertainty),
    )
