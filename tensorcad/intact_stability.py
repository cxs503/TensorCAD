"""Intact-stability criteria: the numbers a righting-arm curve is *judged* by.

``tensorcad.stability`` gives ``GM``, ``tensorcad.gz`` gives the ``GZ(phi)`` curve and
``docs/gz.md`` says plainly that the areas it returns are "the input to such a check".
This module is that check: it evaluates the general intact-stability criteria of the
IMO **Intact Stability Code 2008, Part A, §2.2** (the successor of resolution
A.749(18) §3.1.2) on a sampled curve and returns a per-criterion verdict with the
achieved and required values side by side, never a bare "pass".

The six general criteria, with the default thresholds of the Code:

===========  ====================================================  ==============
id           requirement                                           default
===========  ====================================================  ==============
2.2.1        area under ``GZ`` to 30 deg >=                        0.055 m*rad
2.2.2        area under ``GZ`` to 40 deg, or to the flooding       0.090 m*rad
             angle ``theta_f`` when that is smaller, >=
2.2.3        area between 30 deg and 40 deg (or ``theta_f``), >=   0.030 m*rad
2.2.4        ``GZ`` at 30 deg >=                                   0.200 m
2.2.5        heel at which ``GZ`` is greatest >=                   25 deg
2.2.6        initial metacentric height ``GM_0`` >=                0.150 m
===========  ====================================================  ==============

The thresholds are **defaults, not law**: :class:`IntactThresholds` carries them so a
caller can pin another Administration's or another Code revision's numbers, and the ids
in the report say which clause each number came from. Nothing here paraphrases the
regulation text, and nothing here certifies anything.

**How an area is measured.** The curve is a table of ``(phi_i, GZ_i)``, so the area
``int GZ dphi`` is the trapezoid rule on that table, with ``GZ`` linearly interpolated
where a criterion boundary (30 deg, 40 deg, ``theta_f``) falls between samples. The
result is therefore exact for the piecewise-linear curve the samples define, and the
only error is the sampling: halving the step quarters it (second order), which
``examples/intact_stability`` measures instead of assuming. Angles are converted to
radians inside: a curve quoted in degrees and integrated in degrees is off by a factor
of 57.3, which is the classic way this number is got wrong.

**The closed form the tests use.** For a wall-sided (box) hull the righting arm is

    GZ(phi) = sin(phi) * (GM_0 + BM_0/2 * tan^2 phi),

and it integrates exactly:

    int_0^theta GZ dphi = GM_0 * (1 - cos theta)
                          + BM_0/2 * (sec theta + cos theta - 2).

That is derived here, not quoted: ``d/dtheta`` of the first term is ``GM_0 sin phi``,
of the second ``BM_0/2 * sin phi tan^2 phi``, whose sum is the arm. It lets the area
*machinery* be checked against a value that has no discretisation in it at all -- a
stronger test than comparing two numerical integrations. It stays valid only until the
bilge emerges (``tan phi = 2T/B``) or the deck edge immerses, which is why the example
barge is deeper than it is wide and the comparison stops at 40 deg.

**What this is not.** No weather criterion (§2.3: wind heeling lever, roll amplitude),
no parametric rolling or other second-generation (Part B) failure modes, no damage
stability, no grain or timber-deck cargo rules, no free-surface correction of the curve
itself (a curve built with ``free_surface_i_x`` already carries the flat reduction), and
no ``theta_f`` computation from a compartment layout -- the flooding angle is an input
the caller brings, because it belongs to the subdivision, not to the hull surface. The
verdict is a numerical comparison against published thresholds and nothing more.

Only the Python standard library is used.
"""

from __future__ import annotations

import itertools
import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

SCHEMA_ASSESSMENT = "tensorcad.stability.intact_assessment/1"
SCHEMA_CRITERION = "tensorcad.stability.intact_criterion/1"

SATISFIED = "satisfied"
NOT_SATISFIED = "not_satisfied"
UNDETERMINED = "undetermined"

#: The two angles every area criterion of the Code is quoted at.
ANGLE_30_DEG = 30.0
ANGLE_40_DEG = 40.0


class IntactStabilityError(ValueError):
    """Raised when a criterion cannot be evaluated honestly from the curve given."""


def _round(value: float | None, digits: int = 12) -> float | None:
    return None if value is None else round(float(value), digits)


@dataclass(frozen=True)
class IntactThresholds:
    """The numbers from IS Code 2008 Part A §2.2, replaceable per Administration."""

    area_to_30_m_rad: float = 0.055
    area_to_40_m_rad: float = 0.090
    area_30_to_40_m_rad: float = 0.030
    gz_at_30_deg_m: float = 0.200
    angle_of_max_gz_deg: float = 25.0
    initial_gm_m: float = 0.150

    def __post_init__(self) -> None:
        for name in (
            "area_to_30_m_rad",
            "area_to_40_m_rad",
            "area_30_to_40_m_rad",
            "gz_at_30_deg_m",
            "angle_of_max_gz_deg",
            "initial_gm_m",
        ):
            value = getattr(self, name)
            if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                raise IntactStabilityError(f"threshold {name} must be a finite number: {value!r}")
            if float(value) < 0.0:
                raise IntactStabilityError(f"threshold {name} must not be negative: {value!r}")

    def to_dict(self) -> dict[str, float]:
        return {
            "area_to_30_m_rad": float(self.area_to_30_m_rad),
            "area_to_40_m_rad": float(self.area_to_40_m_rad),
            "area_30_to_40_m_rad": float(self.area_30_to_40_m_rad),
            "gz_at_30_deg_m": float(self.gz_at_30_deg_m),
            "angle_of_max_gz_deg": float(self.angle_of_max_gz_deg),
            "initial_gm_m": float(self.initial_gm_m),
        }


@dataclass(frozen=True)
class CriterionResult:
    """One clause: what it asks for, what the curve gives, and whether that is enough.

    ``status`` is one of ``satisfied``, ``not_satisfied``, ``undetermined``. The third is
    not a soft failure: it means the curve as sampled does not reach the angle the clause
    is quoted at (or the flooding angle the caller withheld), so a verdict would be a
    guess. ``achieved`` is ``None`` exactly when ``status`` is ``undetermined``.
    """

    criterion: str
    quantity: str
    required: float
    achieved: float | None
    units: str
    status: str
    note: str = ""

    def __post_init__(self) -> None:
        if self.status not in (SATISFIED, NOT_SATISFIED, UNDETERMINED):
            raise IntactStabilityError(f"unknown criterion status: {self.status!r}")
        if (self.achieved is None) != (self.status == UNDETERMINED):
            raise IntactStabilityError(
                f"criterion {self.criterion!r}: achieved must be None exactly when the "
                f"status is {UNDETERMINED!r}, got achieved={self.achieved!r} and "
                f"status={self.status!r}"
            )

    @property
    def satisfied(self) -> bool | None:
        """``True``/``False``, or ``None`` when the curve could not answer the clause."""
        if self.status == UNDETERMINED:
            return None
        return self.status == SATISFIED

    @property
    def margin(self) -> float | None:
        """``achieved - required`` in the clause's own units, or ``None``."""
        return None if self.achieved is None else float(self.achieved) - float(self.required)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_CRITERION,
            "criterion": self.criterion,
            "quantity": self.quantity,
            "required": _round(self.required),
            "achieved": _round(self.achieved),
            "units": self.units,
            "status": self.status,
            "margin": _round(self.margin),
            "note": self.note,
        }


@dataclass(frozen=True)
class IntactStabilityAssessment:
    """All six clauses of IS Code §2.2 evaluated on one righting-arm curve."""

    criteria: tuple[CriterionResult, ...]
    thresholds: IntactThresholds
    angles_deg: tuple[float, ...]
    gz_m: tuple[float, ...]
    sampled_range_deg: float
    flooding_angle_deg: float | None = None
    initial_gm_m: float | None = None
    initial_gm_source: str = "supplied"

    @property
    def all_satisfied(self) -> bool:
        """``True`` only when every clause is positively satisfied by the curve."""
        return all(result.status == SATISFIED for result in self.criteria)

    @property
    def undetermined(self) -> tuple[str, ...]:
        return tuple(result.criterion for result in self.criteria if result.status == UNDETERMINED)

    @property
    def failed(self) -> tuple[str, ...]:
        return tuple(result.criterion for result in self.criteria if result.status == NOT_SATISFIED)

    def criterion(self, name: str) -> CriterionResult:
        for result in self.criteria:
            if result.criterion == name:
                return result
        raise KeyError(name)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA_ASSESSMENT,
            "code": "IMO IS Code 2008 Part A §2.2 (general intact stability criteria)",
            "thresholds": self.thresholds.to_dict(),
            "heel_angles_deg": [_round(angle, 9) for angle in self.angles_deg],
            "gz_m": [_round(value, 9) for value in self.gz_m],
            "sampled_range_deg": _round(self.sampled_range_deg, 9),
            "flooding_angle_deg": _round(self.flooding_angle_deg, 9),
            "initial_gm_m": _round(self.initial_gm_m),
            "initial_gm_source": self.initial_gm_source,
            "criteria": [result.to_dict() for result in self.criteria],
            "all_satisfied": self.all_satisfied,
            "failed": list(self.failed),
            "undetermined": list(self.undetermined),
        }


def _samples(angles_deg: Sequence[float], gz: Sequence[float]) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """Validate a curve table: same length, at least two points, strictly increasing."""
    angles = tuple(float(angle) for angle in angles_deg)
    arms = tuple(float(value) for value in gz)
    if len(angles) != len(arms):
        raise IntactStabilityError(f"angles and gz must have the same length: {len(angles)} vs {len(arms)}")
    if len(angles) < 2:
        raise IntactStabilityError("a curve needs at least two samples to be integrable")
    for index, angle in enumerate(angles):
        if not math.isfinite(angle):
            raise IntactStabilityError(f"heel angle {angle!r} is not finite")
        if index and angle <= angles[index - 1]:
            raise IntactStabilityError(f"heel angles must strictly increase: {angles!r}")
    for value in arms:
        if not math.isfinite(value):
            raise IntactStabilityError(f"righting arm {value!r} is not finite")
    return angles, arms


def interpolate_gz(angles_deg: Sequence[float], gz: Sequence[float], angle_deg: float) -> float:
    """``GZ`` at ``angle_deg``, linearly interpolated between the samples that bracket it."""
    angles, arms = _samples(angles_deg, gz)
    target = float(angle_deg)
    if not math.isfinite(target):
        raise IntactStabilityError(f"heel angle {angle_deg!r} is not finite")
    if target < angles[0] or target > angles[-1]:
        raise IntactStabilityError(
            f"angle {target!r} deg is outside the sampled range [{angles[0]!r}, {angles[-1]!r}] deg"
        )
    for index in range(len(angles) - 1):
        left, right = angles[index], angles[index + 1]
        if target == left:
            return arms[index]
        if target == right:
            return arms[index + 1]
        if left < target < right:
            fraction = (target - left) / (right - left)
            return arms[index] + fraction * (arms[index + 1] - arms[index])
    return arms[-1]  # pragma: no cover - target == angles[-1] handled in the loop


def area_between(angles_deg: Sequence[float], gz: Sequence[float], start_deg: float, end_deg: float) -> float:
    """Dynamical stability ``int GZ dphi`` from ``start_deg`` to ``end_deg``, in m*rad.

    Trapezoidal on the samples, with ``GZ`` interpolated at both limits. Signed: an
    interval where ``GZ < 0`` subtracts, which is what the integral means and what makes
    a curve that loses stability early score low instead of being flattered by an
    absolute value.
    """
    angles, arms = _samples(angles_deg, gz)
    start, end = float(start_deg), float(end_deg)
    if not math.isfinite(start) or not math.isfinite(end):
        raise IntactStabilityError(f"interval limits must be finite: {start!r}, {end!r}")
    if end < start:
        raise IntactStabilityError(f"end angle {end!r} deg is before start {start!r} deg")
    if end == start:
        return 0.0
    if start < angles[0]:
        raise IntactStabilityError(
            f"start angle {start!r} deg is below the first sample {angles[0]!r} deg; "
            "the arm there is unknown, not zero"
        )
    if end > angles[-1]:
        raise IntactStabilityError(
            f"end angle {end!r} deg is beyond the last sample {angles[-1]!r} deg; "
            "extend the curve rather than assuming the arm"
        )
    total = 0.0
    current_angle = start
    current_arm = interpolate_gz(angles, arms, start)
    for index, angle in enumerate(angles):
        if angle <= start or angle >= end:
            continue
        total += 0.5 * (current_arm + arms[index]) * math.radians(angle - current_angle)
        current_angle = angle
        current_arm = arms[index]
    end_arm = interpolate_gz(angles, arms, end)
    total += 0.5 * (current_arm + end_arm) * math.radians(end - current_angle)
    return total


def wall_sided_area(angle_deg: float, *, gm: float, bm: float, start_deg: float = 0.0) -> float:
    """Closed form ``int GZ dphi`` for a wall-sided hull, in m*rad.

        int_a^b GZ dphi = GM (1 - cos b) - GM (1 - cos a)
                          + BM/2 [(sec b + cos b) - (sec a + cos a)]

    with ``GZ(phi) = sin(phi) (GM + BM/2 tan^2 phi)``. Exact, and valid only while the
    section stays a trapezoid: below bilge emergence or deck-edge immersion.
    """
    for name, value in (("gm", gm), ("bm", bm), ("angle_deg", angle_deg), ("start_deg", start_deg)):
        if not math.isfinite(float(value)):
            raise IntactStabilityError(f"{name} must be finite: {value!r}")
    if abs(float(angle_deg)) >= 90.0 or abs(float(start_deg)) >= 90.0:
        raise IntactStabilityError("the wall-sided area needs |heel| < 90 degrees")
    if float(angle_deg) < float(start_deg):
        raise IntactStabilityError("angle_deg must not be below start_deg")

    def primitive(degrees: float) -> float:
        phi = math.radians(degrees)
        return float(gm) * (1.0 - math.cos(phi)) + 0.5 * float(bm) * (
            1.0 / math.cos(phi) + math.cos(phi) - 2.0
        )

    return primitive(float(angle_deg)) - primitive(float(start_deg))


def initial_gm_from_curve(
    angles_deg: Sequence[float], gz: Sequence[float], *, max_angle_deg: float = 10.0
) -> tuple[float, float]:
    """Estimate ``GM_0`` as ``GZ(phi)/sin(phi)`` at the smallest angle <= ``max_angle_deg``.

    Returns ``(gm, angle_deg)``. This is the *slope* of the curve at the origin read off a
    chord, so it carries the curvature bias of the arm,
    ``GM_0 + BM_0/2 tan^2 phi + O(phi^2)``, and it therefore reads high: a wall-sided hull
    at 5 deg is over-read by ``BM_0/2 tan^2 5deg``, which for ``BM_0 = 0.27 m`` is 1 mm.
    Pass ``initial_stability(...).gm`` when it is available; this exists so a curve alone
    is still assessable, with the bias stated instead of hidden.
    """
    angles, arms = _samples(angles_deg, gz)
    candidates = [
        (angle, arm)
        for angle, arm in zip(angles, arms, strict=True)
        if 0.0 < angle <= float(max_angle_deg) and math.sin(math.radians(angle)) > 1e-12
    ]
    if not candidates:
        raise IntactStabilityError(
            f"no sample between 0 and {max_angle_deg!r} deg to read the initial slope from"
        )
    angle, arm = candidates[0]
    return arm / math.sin(math.radians(angle)), angle


def _verdict(achieved: float, required: float, tolerance: float) -> str:
    if achieved >= required - tolerance:
        return SATISFIED
    return NOT_SATISFIED


def assess_intact_stability(
    angles_deg: Sequence[float],
    gz: Sequence[float],
    *,
    initial_gm: float | None = None,
    flooding_angle_deg: float | None = None,
    thresholds: IntactThresholds | None = None,
    tolerance: float = 1e-12,
) -> IntactStabilityAssessment:
    """Evaluate the six general criteria of IS Code 2008 Part A §2.2 on a `GZ` table.

    ``initial_gm`` in metres, above the keel datum the curve was built on; when omitted
    it is read off the curve by :func:`initial_gm_from_curve` and flagged as such.
    ``flooding_angle_deg`` is the heel at which an opening immerses; the Code's 40 deg
    intervals are shortened to it when it is smaller, and it is the caller's to supply.
    """
    angles, arms = _samples(angles_deg, gz)
    limits = IntactThresholds() if thresholds is None else thresholds
    if not isinstance(limits, IntactThresholds):
        raise IntactStabilityError("thresholds must be an IntactThresholds instance")
    tol = float(tolerance)
    if not math.isfinite(tol) or tol < 0.0:
        raise IntactStabilityError(f"tolerance must be finite and non-negative: {tolerance!r}")

    range_deg = angles[-1]
    flooding: float | None = None
    if flooding_angle_deg is not None:
        flooding = float(flooding_angle_deg)
        if not math.isfinite(flooding) or flooding <= 0.0:
            raise IntactStabilityError(
                f"flooding_angle_deg must be finite and positive: {flooding_angle_deg!r}"
            )
    upper = ANGLE_40_DEG if flooding is None else min(ANGLE_40_DEG, flooding)
    area_30_limit = ANGLE_30_DEG if flooding is None else min(ANGLE_30_DEG, flooding)

    starts_at_zero = angles[0] <= 0.0

    def shortfall(needed: float) -> str:
        return (
            f"curve sampled only to {range_deg:.6g} deg; the clause is quoted at "
            f"{needed:.6g} deg -- widen the curve instead of assuming the arm"
        )

    def missing_origin() -> str:
        return (
            f"curve starts at {angles[0]:.6g} deg, so the area from 0 deg is unknown "
            "-- start the curve at 0 deg instead of assuming the arm below it"
        )

    def unknown(needed: float) -> str:
        return shortfall(needed) if range_deg < needed else missing_origin()

    results: list[CriterionResult] = []

    # 2.2.1 -- area under the curve to 30 deg.
    if not starts_at_zero or range_deg < area_30_limit:
        results.append(
            CriterionResult(
                criterion="IS-2008-PartA-2.2.1",
                quantity=f"area under GZ up to {area_30_limit:g} deg",
                required=limits.area_to_30_m_rad,
                achieved=None,
                units="m*rad",
                status=UNDETERMINED,
                note=unknown(area_30_limit),
            )
        )
    else:
        area_30 = area_between(angles, arms, 0.0, area_30_limit)
        note = ""
        if flooding is not None and flooding < ANGLE_30_DEG:
            note = f"truncated at the flooding angle {flooding:g} deg, which is below 30 deg"
        results.append(
            CriterionResult(
                criterion="IS-2008-PartA-2.2.1",
                quantity=f"area under GZ up to {area_30_limit:g} deg",
                required=limits.area_to_30_m_rad,
                achieved=area_30,
                units="m*rad",
                status=_verdict(area_30, limits.area_to_30_m_rad, tol),
                note=note,
            )
        )

    # 2.2.2 -- area under the curve to 40 deg, or to the flooding angle.
    if not starts_at_zero or range_deg < upper:
        results.append(
            CriterionResult(
                criterion="IS-2008-PartA-2.2.2",
                quantity=f"area under GZ up to {upper:g} deg",
                required=limits.area_to_40_m_rad,
                achieved=None,
                units="m*rad",
                status=UNDETERMINED,
                note=unknown(upper),
            )
        )
    else:
        area_40 = area_between(angles, arms, 0.0, upper)
        note = f"truncated at the flooding angle {flooding:g} deg" if flooding else ""
        results.append(
            CriterionResult(
                criterion="IS-2008-PartA-2.2.2",
                quantity=f"area under GZ up to {upper:g} deg",
                required=limits.area_to_40_m_rad,
                achieved=area_40,
                units="m*rad",
                status=_verdict(area_40, limits.area_to_40_m_rad, tol),
                note=note,
            )
        )

    # 2.2.3 -- area between 30 deg and 40 deg (or the flooding angle).
    if upper <= ANGLE_30_DEG:
        results.append(
            CriterionResult(
                criterion="IS-2008-PartA-2.2.3",
                quantity="area between 30 deg and the upper limit",
                required=limits.area_30_to_40_m_rad,
                achieved=None,
                units="m*rad",
                status=UNDETERMINED,
                note=(
                    f"the upper limit ({upper:g} deg, flooding angle) does not reach 30 deg, "
                    "so the interval the clause asks about is empty"
                ),
            )
        )
    elif angles[0] > ANGLE_30_DEG or range_deg < upper:
        results.append(
            CriterionResult(
                criterion="IS-2008-PartA-2.2.3",
                quantity=f"area between 30 and {upper:g} deg",
                required=limits.area_30_to_40_m_rad,
                achieved=None,
                units="m*rad",
                status=UNDETERMINED,
                note=(missing_origin() if angles[0] > ANGLE_30_DEG else shortfall(upper)),
            )
        )
    else:
        area_mid = area_between(angles, arms, ANGLE_30_DEG, upper)
        results.append(
            CriterionResult(
                criterion="IS-2008-PartA-2.2.3",
                quantity=f"area between 30 and {upper:g} deg",
                required=limits.area_30_to_40_m_rad,
                achieved=area_mid,
                units="m*rad",
                status=_verdict(area_mid, limits.area_30_to_40_m_rad, tol),
                note="",
            )
        )

    # 2.2.4 -- GZ at 30 deg.
    if range_deg < ANGLE_30_DEG:
        results.append(
            CriterionResult(
                criterion="IS-2008-PartA-2.2.4",
                quantity="GZ at 30 deg",
                required=limits.gz_at_30_deg_m,
                achieved=None,
                units="m",
                status=UNDETERMINED,
                note=shortfall(ANGLE_30_DEG),
            )
        )
    else:
        arm_30 = interpolate_gz(angles, arms, ANGLE_30_DEG)
        results.append(
            CriterionResult(
                criterion="IS-2008-PartA-2.2.4",
                quantity="GZ at 30 deg",
                required=limits.gz_at_30_deg_m,
                achieved=arm_30,
                units="m",
                status=_verdict(arm_30, limits.gz_at_30_deg_m, tol),
                note="",
            )
        )

    # 2.2.5 -- heel of the greatest arm. The peak is on the sampled grid, so it is
    # reported with the grid step it was read on, and a peak sitting on the last sample
    # is not evidence of a peak at all.
    peak_index = max(range(len(arms)), key=lambda index: arms[index])
    peak_angle = angles[peak_index]
    finest_step = min(b - a for a, b in itertools.pairwise(angles))
    peak_note = f"peak read on the sampled grid (finest step {finest_step:g} deg)"
    if peak_index == len(arms) - 1:
        peak_note = (
            "the greatest sampled arm is the last sample, so the true peak lies at or beyond "
            f"{peak_angle:g} deg; widen the curve before relying on this clause"
        )
    results.append(
        CriterionResult(
            criterion="IS-2008-PartA-2.2.5",
            quantity="heel at the greatest GZ",
            required=limits.angle_of_max_gz_deg,
            achieved=peak_angle,
            units="deg",
            status=_verdict(peak_angle, limits.angle_of_max_gz_deg, tol),
            note=peak_note,
        )
    )

    # 2.2.6 -- initial metacentric height.
    if initial_gm is None:
        try:
            gm_value, gm_angle = initial_gm_from_curve(angles, arms)
        except IntactStabilityError as error:
            results.append(
                CriterionResult(
                    criterion="IS-2008-PartA-2.2.6",
                    quantity="initial metacentric height GM_0",
                    required=limits.initial_gm_m,
                    achieved=None,
                    units="m",
                    status=UNDETERMINED,
                    note=f"no GM supplied and the curve gives none: {error}",
                )
            )
            gm_source = "unavailable"
            gm_reported: float | None = None
        else:
            gm_source = f"curve chord at {gm_angle:g} deg (GZ/sin phi; reads high by BM/2 tan^2 phi)"
            gm_reported = gm_value
            results.append(
                CriterionResult(
                    criterion="IS-2008-PartA-2.2.6",
                    quantity="initial metacentric height GM_0",
                    required=limits.initial_gm_m,
                    achieved=gm_value,
                    units="m",
                    status=_verdict(gm_value, limits.initial_gm_m, tol),
                    note=f"derived from the curve, {gm_source}; prefer initial_stability()",
                )
            )
    else:
        gm_value = float(initial_gm)
        if not math.isfinite(gm_value):
            raise IntactStabilityError(f"initial_gm must be finite: {initial_gm!r}")
        gm_source = "supplied"
        gm_reported = gm_value
        results.append(
            CriterionResult(
                criterion="IS-2008-PartA-2.2.6",
                quantity="initial metacentric height GM_0",
                required=limits.initial_gm_m,
                achieved=gm_value,
                units="m",
                status=_verdict(gm_value, limits.initial_gm_m, tol),
                note="",
            )
        )

    return IntactStabilityAssessment(
        criteria=tuple(results),
        thresholds=limits,
        angles_deg=angles,
        gz_m=arms,
        sampled_range_deg=range_deg,
        flooding_angle_deg=flooding,
        initial_gm_m=gm_reported,
        initial_gm_source=gm_source,
    )


def assess_curve(
    curve: Any,
    *,
    initial_gm: float | None = None,
    flooding_angle_deg: float | None = None,
    thresholds: IntactThresholds | None = None,
    tolerance: float = 1e-12,
) -> IntactStabilityAssessment:
    """:func:`assess_intact_stability` on a ``RightingArmCurve`` (or any angles/gz pair).

    Anything with ``angles_deg`` and ``gz`` attributes works, so both routes to the curve
    -- ``tensorcad.gz.righting_arm_curve`` and its bisection twin -- can be judged
    without this module knowing which one produced the numbers.
    """
    try:
        angles = curve.angles_deg
        arms = curve.gz
    except AttributeError as error:  # pragma: no cover - defensive
        raise IntactStabilityError(
            f"expected an object with angles_deg and gz attributes: {error}"
        ) from error
    return assess_intact_stability(
        angles,
        arms,
        initial_gm=initial_gm,
        flooding_angle_deg=flooding_angle_deg,
        thresholds=thresholds,
        tolerance=tolerance,
    )
