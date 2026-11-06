"""Validation harness for :mod:`tensorcad.intact_stability`'s weather criterion.

``tensorcad.intact_stability`` turns a measured ``GZ`` curve into the IS Code 2008 part A
rule check. ``examples/intact_stability`` already validates the *curve* it is given (mesh
against the wall-sided closed form) and the *verdict flags*; what it does not do is verify
the numbers the weather criterion compares. This script does exactly that:

* the two areas ``2.3.2`` weighs against each other -- ``b = int GZ`` and ``a = lw2 (theta2 -
  theta1)`` -- against the antiderivative of the curve, on both sides of the ``b >= a``
  boundary;
* ``theta_2`` against the least of its four caps (50 deg, down-flooding, ``theta_f``, the
  second intercept ``theta_c``), with ``theta_c`` re-solved by an independent bisection;
* the substitution of the flooding angle ``theta_f`` for the 40 deg limit in ``2.2``,
  including the negative control that the number actually changes;
* the boundary of the "shall not exceed" clauses (``theta_1 = theta_f`` is legal, a hair
  above it is not), which the ``Criterion.comparison`` field exists to get right;
* the residual of the two areas under grid refinement, to show it is quadrature error and
  not a modelling error;
* the same curve through :mod:`tensorcad.statutory`, so the two implementations of section
  2.2 are read against one closed form rather than against each other.

Every number in ``docs/intact_stability.md`` s. 6 comes from here, and the script exits
non-zero if one of its claims stops holding::

    python tools/intact_stability_validation.py [--output docs/intact_stability]

Artefacts: ``intact_stability_validation.json`` (machine-readable evidence) and
``intact_stability_validation.png`` (the figure the document embeds).
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import tensorcad  # noqa: E402 - path set above
from tensorcad import intact_stability as iscode  # noqa: E402 - path set above
from tensorcad.statutory import check_intact_stability  # noqa: E402 - path set above

SCHEMA = "tensorcad.intact-stability-validation/1"

#: The reference body: a wall-sided box barge, B = 6 m, T = 3 m, KG = 2 m, D = 8 m.
BREADTH, DRAFT, DEPTH, LENGTH, KG = 6.0, 3.0, 8.0, 10.0, 2.0
BM0 = BREADTH**2 / (12.0 * DRAFT)
GM0 = 0.5 * DRAFT + BM0 - KG

#: The wind geometry is only a way of dialling in a steady lever ``lw1``; one windage
#: area, one arm and one displacement are fixed and the area is solved for the lever.
WIND_ARM_M = 5.0
DISPLACEMENT_T = 1000.0
GRAVITY = 9.81


# --------------------------------------------------------------------------------------
# closed forms: the wall-sided barge's arm and its antiderivative
# --------------------------------------------------------------------------------------
def wall_sided_gz(angle_deg: float) -> float:
    """``GZ = sin(phi) (GM0 + BM0/2 tan^2 phi)`` -- the box barge below bilge emergence."""
    phi = math.radians(angle_deg)
    tangent = math.tan(phi)
    return math.sin(phi) * (GM0 + 0.5 * BM0 * tangent * tangent)


def wall_sided_area(angle_deg: float) -> float:
    """``int_0^angle GZ dphi = GM0 (1 - cos) + BM0/2 (sec + cos - 2)``."""
    phi = math.radians(angle_deg)
    return GM0 * (1.0 - math.cos(phi)) + 0.5 * BM0 * (
        (1.0 / math.cos(phi) + math.cos(phi)) - 2.0
    )


def turnover_gz(angle_deg: float) -> float:
    """``GZ = 0.3 sin(2 phi)``: peaks at 45 deg, so the gust lever is crossed twice."""
    return 0.3 * math.sin(2.0 * math.radians(angle_deg))


def turnover_area(angle_deg: float) -> float:
    """``int_0^angle 0.3 sin(2 phi) dphi = 0.15 (1 - cos 2 phi)``."""
    return 0.15 * (1.0 - math.cos(2.0 * math.radians(angle_deg)))


def ladder(step_deg: float, limit_deg: float, curve) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """An ``(angles, gz)`` table of ``curve`` from 0 to ``limit_deg`` at ``step_deg``."""
    count = int(round(limit_deg / step_deg))
    if abs(count * step_deg - limit_deg) > 1e-9:
        raise ValueError(f"step {step_deg!r} deg does not divide {limit_deg!r} deg")
    angles = [index * step_deg for index in range(count + 1)]
    return tuple(angles), tuple(curve(angle) for angle in angles)


def relative_error(attained: float, reference: float) -> float:
    return abs(attained - reference) / abs(reference)


# --------------------------------------------------------------------------------------
# running the module: a steady lever is dialled in through the windage area
# --------------------------------------------------------------------------------------
def weather(
    curve,
    steady_lever_m: float,
    *,
    roll_amplitude_deg: float = 0.0,
    flooding_angle_deg: float | None = None,
    downflooding_angle_deg: float = iscode.WEATHER_THETA2_CAP_DEG,
    deck_edge_angle_deg: float | None = None,
) -> iscode.WeatherCriterionReport:
    """``weather_criterion`` with the windage area solved for the requested ``lw1``."""
    area = steady_lever_m * 1000.0 * GRAVITY * DISPLACEMENT_T / (
        iscode.WIND_PRESSURE_PA * WIND_ARM_M
    )
    return iscode.weather_criterion(
        curve,
        wind_area_m2=area,
        wind_arm_m=WIND_ARM_M,
        displacement_t=DISPLACEMENT_T,
        roll_amplitude_deg=roll_amplitude_deg,
        flooding_angle_deg=flooding_angle_deg,
        downflooding_angle_deg=downflooding_angle_deg,
        deck_edge_angle_deg=deck_edge_angle_deg,
    )


def bisect_root(function, low: float, high: float, iterations: int = 200) -> float:
    """Independent root of a monotone-crossing ``function`` on ``[low, high]``."""
    f_low, f_high = function(low), function(high)
    if f_low == 0.0:
        return low
    if f_high == 0.0:
        return high
    if f_low * f_high > 0.0:
        raise ValueError(f"no sign change on [{low}, {high}]: {f_low!r}, {f_high!r}")
    for _ in range(iterations):
        middle = 0.5 * (low + high)
        f_middle = function(middle)
        if f_middle == 0.0:
            return middle
        if f_low * f_middle < 0.0:
            high, f_high = middle, f_middle
        else:
            low, f_low = middle, f_middle
    return 0.5 * (low + high)


# --------------------------------------------------------------------------------------
# the claims
# --------------------------------------------------------------------------------------
def claim_area_limits_with_flooding_angle() -> dict[str, object]:
    """The two areas cut at ``theta_f = 35`` deg: the substituted limit is the closed form's.

    The sibling module labels these ``2.2.2`` and ``2.2.3``; the code itself calls them
    ``2.2.1b`` and ``2.2.1c`` (the numbering divergence documented in ``docs/statutory.md``
    s. 6). Both spellings are given below so a report cannot cite the wrong line.
    """
    curve = ladder(1.0, 50.0, wall_sided_gz)
    report = iscode.general_criteria(curve, gm0_m=GM0, flooding_angle_deg=35.0)
    reference_to_limit = wall_sided_area(35.0)
    reference_30_to_limit = wall_sided_area(35.0) - wall_sided_area(30.0)
    error_to_limit = relative_error(report.area_0_40_m_rad, reference_to_limit)
    error_30_to_limit = relative_error(report.area_30_40_m_rad, reference_30_to_limit)
    return {
        "claim": "the flooding angle replaces 40 deg and both areas cut at theta_f match the "
        "antiderivative (code 2.2.1b/2.2.1c, sibling 2.2.2/2.2.3)",
        "holds": bool(
            report.upper_angle_deg == 35.0
            and error_to_limit < 5e-4
            and error_30_to_limit < 2e-3
        ),
        "measured": f"upper limit {report.upper_angle_deg} deg; A(0,35) = "
        f"{report.area_0_40_m_rad:.9f} vs {reference_to_limit:.9f} "
        f"(rel {error_to_limit:.2e}); A(30,35) = {report.area_30_40_m_rad:.9f} vs "
        f"{reference_30_to_limit:.9f} (rel {error_30_to_limit:.2e}) on a 1 deg grid",
        "expected": "upper limit 35 deg, both areas inside 5e-4 / 2e-3 relative "
        "(1 deg trapezoid on the barge)",
        "error": max(error_to_limit, error_30_to_limit),
        "tolerance": 2e-3,
        "numbers": {
            "upper_limit_deg": report.upper_angle_deg,
            "area_to_limit_m_rad": report.area_0_40_m_rad,
            "reference_to_limit_m_rad": reference_to_limit,
            "error_to_limit": error_to_limit,
            "area_30_to_limit_m_rad": report.area_30_40_m_rad,
            "reference_30_to_limit_m_rad": reference_30_to_limit,
            "error_30_to_limit": error_30_to_limit,
        },
    }


def claim_flooding_angle_is_substituted() -> dict[str, object]:
    """Negative control: without ``theta_f`` the same call stops at 40 deg, so the two differ."""
    curve = ladder(1.0, 50.0, wall_sided_gz)
    with_flooding = iscode.general_criteria(curve, gm0_m=GM0, flooding_angle_deg=35.0)
    without = iscode.general_criteria(curve, gm0_m=GM0)
    difference = abs(with_flooding.area_0_40_m_rad - without.area_0_40_m_rad)
    return {
        "claim": "theta_f is actually used, not accepted and ignored",
        "holds": bool(without.upper_angle_deg == 40.0 and difference > 1e-3),
        "measured": f"A(0,40) = {without.area_0_40_m_rad:.9f} without theta_f vs "
        f"{with_flooding.area_0_40_m_rad:.9f} at theta_f = 35 deg "
        f"(difference {difference:.2e} m rad)",
        "expected": "40 deg used when theta_f is omitted; the two areas differ by more than "
        "1e-3 m rad",
        "error": difference,
        "tolerance": 1e-3,
        "numbers": {
            "area_with_theta_f_m_rad": with_flooding.area_0_40_m_rad,
            "area_without_theta_f_m_rad": without.area_0_40_m_rad,
            "difference_m_rad": difference,
        },
    }


def claim_gust_area_b() -> dict[str, object]:
    """``2.3.2``'s ``b`` over the 50 deg cap against the antiderivative, ``theta_1 = 8`` deg."""
    curve = ladder(1.0, 50.0, wall_sided_gz)
    report = weather(curve, 0.02, roll_amplitude_deg=8.0)
    reference = wall_sided_area(50.0) - wall_sided_area(8.0)
    error = relative_error(report.area_b_m_rad, reference)
    return {
        "claim": "gust area b over [theta_1, theta_2] matches the antiderivative of GZ",
        "holds": bool(report.theta2_deg == 50.0 and error < 1e-3),
        "measured": f"theta_1 = 8 deg, theta_2 = {report.theta2_deg} deg, b = "
        f"{report.area_b_m_rad:.9f} vs A(50) - A(8) = {reference:.9f} (rel {error:.2e})",
        "expected": "theta_2 = 50 deg (the code's cap), b inside 1e-3 relative of the "
        "antiderivative on a 1 deg grid",
        "error": error,
        "tolerance": 1e-3,
        "numbers": {
            "theta1_deg": 8.0,
            "theta2_deg": report.theta2_deg,
            "area_b_m_rad": report.area_b_m_rad,
            "reference_m_rad": reference,
            "relative_error": error,
        },
    }


def claim_gust_area_a() -> dict[str, object]:
    """``a = lw2 (theta_2 - theta_1)`` exactly: a constant lever integrates to a rectangle."""
    curve = ladder(1.0, 50.0, wall_sided_gz)
    report = weather(curve, 0.02, roll_amplitude_deg=8.0)
    reference = report.levers.gust_m * math.radians(report.theta2_deg - 8.0)
    error = relative_error(report.area_a_m_rad, reference)
    return {
        "claim": "gust area a is lw2 times the span (the lever is constant)",
        "holds": bool(error < 1e-12),
        "measured": f"a = {report.area_a_m_rad:.15f} vs lw2 (theta_2 - theta_1) = "
        f"{reference:.15f} (rel {error:.2e}), lw2 = {report.levers.gust_m:.9f} m",
        "expected": "rel < 1e-12: exact geometry, no quadrature involved",
        "error": error,
        "tolerance": 1e-12,
        "numbers": {
            "area_a_m_rad": report.area_a_m_rad,
            "reference_m_rad": reference,
            "lw2_m": report.levers.gust_m,
            "span_deg": report.theta2_deg - 8.0,
            "relative_error": error,
        },
    }


def claim_theta2_is_the_least_cap() -> dict[str, object]:
    """``theta_2`` follows whichever cap binds: 50 deg, down-flooding, or ``theta_f``."""
    curve = ladder(1.0, 50.0, wall_sided_gz)
    rows: list[dict[str, object]] = []
    holds = True
    worst = 0.0
    for name, kwargs, expected, reference_span in (
        ("no cap but the code's 50 deg", {}, 50.0, (0.0, 50.0)),
        (
            "down-flooding at 38 deg",
            {"downflooding_angle_deg": 38.0},
            38.0,
            (0.0, 38.0),
        ),
        ("flooding angle 35 deg", {"flooding_angle_deg": 35.0}, 35.0, (0.0, 35.0)),
    ):
        report = weather(curve, 0.02, **kwargs)
        reference = wall_sided_area(reference_span[1]) - wall_sided_area(reference_span[0])
        error = relative_error(report.area_b_m_rad, reference)
        worst = max(worst, error)
        nominal = report.theta2_deg == expected and error < 2e-3
        holds = holds and nominal
        rows.append(
            {
                "cap": name,
                "theta2_deg": report.theta2_deg,
                "expected_theta2_deg": expected,
                "area_b_m_rad": report.area_b_m_rad,
                "reference_m_rad": reference,
                "relative_error": error,
                "holds": nominal,
            }
        )
    return {
        "claim": "theta_2 is the least of the 50 deg, down-flooding and flooding caps",
        "holds": bool(holds),
        "measured": "; ".join(
            f"{row['cap']}: theta_2 = {row['theta2_deg']} deg, b rel {row['relative_error']:.2e}"
            for row in rows
        ),
        "expected": "theta_2 equals the binding cap in all three cases and b stays inside "
        "2e-3 relative of the antiderivative",
        "error": worst,
        "tolerance": 2e-3,
        "rows": rows,
    }


def claim_second_intercept() -> dict[str, object]:
    """``theta_c`` against an independent bisection, then ``b`` over ``[0, theta_c]``."""
    curve = ladder(0.1, 60.0, turnover_gz)
    steady = 0.298 / iscode.WIND_GUST_FACTOR
    report = weather(curve, steady)
    gust = report.levers.gust_m
    root = bisect_root(lambda angle: turnover_gz(angle) - gust, 45.0, 60.0)
    heel_error = abs(report.theta_c_deg - root)
    reference = turnover_area(root)
    area_error = relative_error(report.area_b_m_rad, reference)
    return {
        "claim": "the second intercept theta_c is solved right and b is measured to it",
        "holds": bool(
            heel_error < 5e-3 and area_error < 1e-3 and report.theta2_deg == report.theta_c_deg
        ),
        "measured": f"theta_c = {report.theta_c_deg:.6f} deg vs bisection {root:.6f} deg "
        f"(|diff| {heel_error:.2e} deg); b = {report.area_b_m_rad:.9f} vs A(theta_c) = "
        f"{reference:.9f} (rel {area_error:.2e}); theta_2 = {report.theta2_deg:.6f} deg",
        "expected": "|theta_c - bisection| < 5e-3 deg, b inside 1e-3 relative, theta_2 = "
        "theta_c because the cap is 50 deg",
        "error": max(heel_error, area_error),
        "tolerance": 5e-3,
        "numbers": {
            "theta_c_deg": report.theta_c_deg,
            "theta_c_bisection_deg": root,
            "theta_c_error_deg": heel_error,
            "theta2_deg": report.theta2_deg,
            "lw2_m": gust,
            "area_b_m_rad": report.area_b_m_rad,
            "reference_m_rad": reference,
            "area_relative_error": area_error,
        },
    }


def claim_gust_boundary() -> dict[str, object]:
    """The ``b >= a`` verdict flips exactly at the lever the arithmetic predicts."""
    curve = ladder(0.5, 50.0, wall_sided_gz)
    reference_b = wall_sided_area(50.0)
    crossing = reference_b / (iscode.WIND_GUST_FACTOR * math.radians(50.0))
    below = weather(curve, crossing * (1.0 - 1e-3))
    above = weather(curve, crossing * (1.0 + 1e-3))
    clause_below = next(item for item in below.criteria if item.clause == "2.3.2")
    clause_above = next(item for item in above.criteria if item.clause == "2.3.2")
    margin_below, margin_above = clause_below.margin, clause_above.margin
    analytic = below.area_b_m_rad - below.area_a_m_rad
    margin_error = abs(margin_below - analytic) / abs(analytic)
    return {
        "claim": "2.3.2 passes below the analytic lever where b = a and fails above it",
        "holds": bool(
            margin_below > 0.0
            and margin_above < 0.0
            and margin_error < 1e-6
            and clause_below.satisfied
            and not clause_above.satisfied
        ),
        "measured": f"b = a at lw1 = {crossing:.9f} m; margin {margin_below:+.3e} at "
        f"-0.1 % (2.3.2 satisfied {clause_below.satisfied}) and {margin_above:+.3e} at "
        f"+0.1 % (satisfied {clause_above.satisfied}); reported margin vs b - a rel "
        f"{margin_error:.1e}",
        "expected": "margin positive below the crossing, negative above, equal to b - a to "
        "1e-6 relative, verdict flags following the sign",
        "error": margin_error,
        "tolerance": 1e-6,
        "numbers": {
            "crossing_lw1_m": crossing,
            "margin_below_m_rad": margin_below,
            "margin_above_m_rad": margin_above,
            "analytic_margin_below_m_rad": analytic,
            "margin_relative_error": margin_error,
        },
    }


def claim_roll_amplitude_boundary() -> dict[str, object]:
    """``2.3.3``: ``theta_1 = theta_f`` is legal ("shall not exceed"), a hair above is not."""
    curve = ladder(1.0, 50.0, wall_sided_gz)
    at_limit = weather(curve, 0.02, roll_amplitude_deg=20.0, flooding_angle_deg=20.0)
    clause = next(item for item in at_limit.criteria if item.clause == "2.3.3")
    over = weather(curve, 0.02, roll_amplitude_deg=20.5, flooding_angle_deg=20.0)
    clause_over = next(item for item in over.criteria if item.clause == "2.3.3")
    under = weather(curve, 0.02, roll_amplitude_deg=8.0, flooding_angle_deg=20.0)
    return {
        "claim": "2.3.3 at equality passes, one step above fails, and the empty span reports "
        "no areas rather than a zero",
        "holds": bool(
            clause.satisfied
            and clause.margin == 0.0
            and not clause_over.satisfied
            and at_limit.area_b_m_rad is None
            and at_limit.area_a_m_rad is None
            and under.area_b_m_rad is not None
        ),
        "measured": f"theta_1 = theta_f = 20 deg: 2.3.3 satisfied {clause.satisfied} "
        f"(margin {clause.margin}); theta_1 = 20.5 deg: satisfied {clause_over.satisfied} "
        f"(margin {clause_over.margin}); areas at equality b = {at_limit.area_b_m_rad}, "
        f"a = {at_limit.area_a_m_rad}; at theta_1 = 8 deg b = {under.area_b_m_rad:.9f}",
        "expected": "equality passes with margin 0, 20.5 deg fails, and the empty "
        "theta_2 - theta_1 span is reported as None (not 0) for both areas",
        "error": abs(clause.margin - 0.0),
        "tolerance": 0.0,
        "numbers": {
            "flooding_angle_deg": 20.0,
            "theta1_at_limit_deg": 20.0,
            "satisfied_at_limit": clause.satisfied,
            "margin_at_limit": clause.margin,
            "theta1_over_deg": 20.5,
            "satisfied_over": clause_over.satisfied,
            "margin_over": clause_over.margin,
            "area_b_at_limit_m_rad": at_limit.area_b_m_rad,
            "area_a_at_limit_m_rad": at_limit.area_a_m_rad,
        },
    }


def claim_cross_module() -> dict[str, object]:
    """The same curve through :mod:`tensorcad.statutory` and through the sibling."""
    curve = ladder(1.0, 50.0, wall_sided_gz)
    sibling = iscode.general_criteria(curve, gm0_m=GM0, flooding_angle_deg=35.0)
    statutory = check_intact_stability(curve, gm0=GM0, flooding_angle_deg=35.0)
    reference = wall_sided_area(35.0)
    statutory_error = relative_error(statutory.area_to_limit_m_rad, reference)
    sibling_error = relative_error(sibling.area_0_40_m_rad, reference)
    spread = relative_error(sibling.area_0_40_m_rad, statutory.area_to_limit_m_rad)
    return {
        "claim": "both implementations of 2.2 land on the same closed form at theta_f = 35",
        "holds": bool(
            statutory.area_limit_deg == 35.0
            and statutory_error < 1e-7
            and sibling_error < 5e-4
            and spread < 5e-4
        ),
        "measured": f"statutory A(0,35) = {statutory.area_to_limit_m_rad:.12f} "
        f"(rel {statutory_error:.2e}, its own uncertainty "
        f"{statutory.area_uncertainty_m_rad:.2e}); sibling = {sibling.area_0_40_m_rad:.9f} "
        f"(rel {sibling_error:.2e}); spread {spread:.2e}",
        "expected": "statutory inside 1e-7 relative (adaptive Simpson on the wall-sided "
        "curve), sibling inside 5e-4 (1 deg trapezoid), spread inside 5e-4",
        "error": max(statutory_error, sibling_error, spread),
        "tolerance": 5e-4,
        "numbers": {
            "area_limit_deg": statutory.area_limit_deg,
            "statutory_m_rad": statutory.area_to_limit_m_rad,
            "statutory_uncertainty_m_rad": statutory.area_uncertainty_m_rad,
            "statutory_relative_error": statutory_error,
            "sibling_m_rad": sibling.area_0_40_m_rad,
            "sibling_relative_error": sibling_error,
            "spread": spread,
            "reference_m_rad": reference,
        },
    }


def claim_grid_convergence() -> dict[str, object]:
    """The residual of ``b`` is quadrature error: it falls like ``h^2``."""
    rows: list[dict[str, object]] = []
    reference = wall_sided_area(50.0) - wall_sided_area(8.0)
    for step in (5.0, 2.5, 1.25, 0.625):
        curve = ladder(step, 50.0, wall_sided_gz)
        report = weather(curve, 0.02, roll_amplitude_deg=8.0)
        rows.append(
            {
                "step_deg": step,
                "samples": len(curve[0]),
                "area_b_m_rad": report.area_b_m_rad,
                "relative_error": relative_error(report.area_b_m_rad, reference),
            }
        )
    ratios = [
        rows[index]["relative_error"] / rows[index + 1]["relative_error"]
        for index in range(len(rows) - 1)
    ]
    tail = ratios[-2:]
    holds = all(3.5 <= ratio <= 4.5 for ratio in tail)
    return {
        "claim": "the gust-area residual falls as the step squared (trapezoid signature)",
        "holds": bool(holds),
        "measured": "errors "
        + ", ".join(f"{row['relative_error']:.2e} at {row['step_deg']} deg" for row in rows)
        + "; successive ratios " + ", ".join(f"{ratio:.2f}" for ratio in ratios),
        "expected": "successive ratios in 3.5-4.5 from 2.5 deg down, i.e. O(h^2), so the "
        "residual is the grid and not the rule",
        "error": max(abs(ratio - 4.0) for ratio in tail),
        "tolerance": 0.5,
        "rows": rows,
        "ratios": ratios,
    }


CHECKS = (
    ("area_limits_with_flooding_angle", claim_area_limits_with_flooding_angle),
    ("flooding_angle_is_substituted", claim_flooding_angle_is_substituted),
    ("gust_area_b", claim_gust_area_b),
    ("gust_area_a", claim_gust_area_a),
    ("theta2_is_the_least_cap", claim_theta2_is_the_least_cap),
    ("second_intercept", claim_second_intercept),
    ("gust_boundary", claim_gust_boundary),
    ("roll_amplitude_boundary", claim_roll_amplitude_boundary),
    ("cross_module", claim_cross_module),
    ("grid_convergence", claim_grid_convergence),
)


def evaluate() -> list[dict[str, object]]:
    """Run every claim and return the records (``holds`` is the machine verdict)."""
    records: list[dict[str, object]] = []
    for name, function in CHECKS:
        record = dict(function())
        record["name"] = name
        records.append(record)
    return records


MARKDOWN_BEGIN = "<!-- BEGIN intact_stability_validation (generated by tools/intact_stability_validation.py) -->"
MARKDOWN_END = "<!-- END intact_stability_validation -->"

SECTION_TITLE = "## 6. An independent check of the numbers the weather criterion compares"


def _sci(value: float) -> str:
    return f"{value:.2e}"


def render_section(payload: dict[str, object]) -> str:
    """Render the document section from the json payload -- every number comes from it."""
    by_name = {record["name"]: record for record in payload["claims"]}  # type: ignore[index]
    limits = by_name["area_limits_with_flooding_angle"]["numbers"]
    gust_b = by_name["gust_area_b"]["numbers"]
    gust_a = by_name["gust_area_a"]["numbers"]
    second = by_name["second_intercept"]["numbers"]
    boundary = by_name["gust_boundary"]["numbers"]
    rolling = by_name["roll_amplitude_boundary"]["numbers"]
    substituted = by_name["flooding_angle_is_substituted"]["numbers"]
    cross = by_name["cross_module"]["numbers"]
    convergence = by_name["grid_convergence"]["rows"]
    ratios = by_name["grid_convergence"]["ratios"]
    caps = by_name["theta2_is_the_least_cap"]["rows"]

    lines: list[str] = [
        MARKDOWN_BEGIN,
        SECTION_TITLE,
        "",
        "`examples/intact_stability` (§1-§5) validates the *curve* the rule is applied to and the",
        "verdict flags. It reports `b` and `a` but does not check them, and it reads `θ_c` off the",
        "same interpolant the module uses. `tools/intact_stability_validation.py` closes that gap:",
        "ten claims, each compared against an antiderivative or an independent bisection,",
        "aggregated under `docs/intact_stability/intact_stability_validation.json` and",
        "`intact_stability_validation.png`, exiting non-zero if one of them stops holding.",
        "`tests/test_intact_stability_validation.py` re-runs it and pins both the json and this",
        "section to a fresh run, so neither can go stale. **Nothing here is typed by hand: this",
        "block is written by `--update-docs`.**",
        "",
        "The reference is the barge's arm and its antiderivative (`B = 6` m, `T = 3` m, `KG = 2` m,",
        "so `GM_0 = 0.5` m and `BM_0 = 1.0` m),",
        "",
        "```",
        "GZ(phi)          = sin(phi) (GM0 + BM0/2 tan^2 phi)",
        "int_0^a GZ dphi  = GM0 (1 - cos a) + BM0/2 (sec a + cos a - 2)",
        "```",
        "",
        "plus a second curve `GZ = 0.3 sin(2 phi)`, which peaks at 45° and so lets `θ_2` be the",
        "*second intercept* rather than one of the code's caps.",
        "",
        "### 6.1 The two gust areas and the flooding-angle areas",
        "",
        "| quantity | the module | reference | relative error |",
        "|---|---|---|---|",
        f"| `b` over `[θ_1, θ_2]`, `θ_1 = {gust_b['theta1_deg']:.0f}°`, "
        f"`θ_2 = {gust_b['theta2_deg']:.0f}°` | {gust_b['area_b_m_rad']:.9f} m rad | "
        f"antiderivative `A(θ_2) − A(θ_1) = {gust_b['reference_m_rad']:.9f}` | "
        f"{_sci(gust_b['relative_error'])} |",
        f"| `a = lw2 (θ_2 − θ_1)` | {gust_a['area_a_m_rad']:.15f} m rad | "
        f"`{gust_a['lw2_m']:.9f} m × {gust_a['span_deg']:.0f}°` in rad = "
        f"{gust_a['reference_m_rad']:.15f} | {_sci(gust_a['relative_error'])} |",
        f"| area to `θ_f = {limits['upper_limit_deg']:.0f}°` (code `2.2.1b`, sibling `2.2.2`) | "
        f"{limits['area_to_limit_m_rad']:.9f} m rad | `A(θ_f) = "
        f"{limits['reference_to_limit_m_rad']:.9f}` | {_sci(limits['error_to_limit'])} |",
        f"| area `[30°, θ_f]` (code `2.2.1c`, sibling `2.2.3`) | "
        f"{limits['area_30_to_limit_m_rad']:.9f} m rad | "
        f"`A(θ_f) − A(30°) = {limits['reference_30_to_limit_m_rad']:.9f}` | "
        f"{_sci(limits['error_30_to_limit'])} |",
        "",
        "`b` is a trapezoidal integral on the 1° ladder, so its residual is the grid and not the",
        "rule — §6.6 shows it falling as `h²`. `a` is a rectangle and the module multiplies it out",
        "rather than integrating a constant, which is why it lands on round-off.",
        "",
        "### 6.2 `θ_2` is the least of its caps",
        "",
        "| binding cap | `θ_2` [deg] | expected | `b` [m rad] | reference [m rad] | rel. error |",
        "|---|---|---|---|---|---|",
    ]
    for row in caps:
        lines.append(
            f"| {row['cap']} | {row['theta2_deg']:.3f} | {row['expected_theta2_deg']:.0f} | "
            f"{row['area_b_m_rad']:.9f} | {row['reference_m_rad']:.9f} | "
            f"{_sci(row['relative_error'])} |"
        )
    lines += [
        "",
        "Three calls on one curve, one per cap: the code's own 50°, a down-flooding angle of 38°,",
        "and a flooding angle of 35°. In each, the reported `θ_2` is the cap that binds and `b` is",
        "measured to *that* angle — which is what makes the cap a rule rather than a truncation.",
        "",
        "### 6.3 The second intercept `θ_c`",
        "",
        f"On `GZ = 0.3 sin(2 phi)` with `lw2 = {second['lw2_m']:.6f}` m the gust lever is crossed twice, so",
        "`θ_2 = min(50°, θ_c) = θ_c`:",
        "",
        "| quantity | the module | independent bisection / antiderivative | difference |",
        "|---|---|---|---|",
        f"| `θ_c` | {second['theta_c_deg']:.6f} deg | {second['theta_c_bisection_deg']:.6f} deg | "
        f"{second['theta_c_error_deg']:.2e} deg |",
        f"| `b` over `[0, θ_c]` | {second['area_b_m_rad']:.9f} m rad | "
        f"{second['reference_m_rad']:.9f} m rad | {_sci(second['area_relative_error'])} relative |",
        "",
        "The bisection runs on the analytic `GZ` itself (200 halvings in the harness), not on the",
        "module's interpolant, so agreement means the module found the crossing of the real curve.",
        "",
        "### 6.4 Where `b ≥ a` flips",
        "",
        "`b` is a constant for a given curve, `A(θ_2) − A(θ_1)`, while `a` grows linearly in `lw1`,",
        "so clause 2.3.2 flips at `lw1 = b / (1.5 θ_2)`. Measured on the module:",
        "",
        "| quantity | value |",
        "|---|---|",
        f"| analytic flip lever | {boundary['crossing_lw1_m']:.9f} m |",
        f"| reported margin at `−0.1 %` of it | {boundary['margin_below_m_rad']:+.3e} m rad (2.3.2 passes) |",
        f"| reported margin at `+0.1 %` of it | {boundary['margin_above_m_rad']:+.3e} m rad (2.3.2 fails) |",
        f"| reported margin vs `b − a` on that call | {_sci(boundary['margin_relative_error'])} relative |",
        "",
        "Both sides are checked, so this tests the comparison *direction* as well as the areas: a",
        "flipped sign or a single `>=` would put the wrong flag on one of the two calls.",
        "",
        "### 6.5 The \"shall not exceed\" boundary",
        "",
        f"Clause 2.3.3 says the windward roll amplitude *shall not exceed* the flooding angle, so",
        f"equality must pass. At `θ_1 = θ_f = {rolling['flooding_angle_deg']:.0f}°`:",
        "",
        f"* 2.3.3 is `satisfied = {rolling['satisfied_at_limit']}` with `margin = "
        f"{rolling['margin_at_limit']:.1f}` deg;",
        f"* `b` and `a` are `{rolling['area_b_at_limit_m_rad']}` and "
        f"`{rolling['area_a_at_limit_m_rad']}` — the span `θ_2 − θ_1` is empty, and the report",
        "  says so instead of quoting a zero;",
        f"* at `θ_1 = {rolling['theta1_over_deg']}°` the clause is `satisfied = "
        f"{rolling['satisfied_over']}` with `margin = {rolling['margin_over']:.1f}` deg.",
        "",
        "### 6.6 The flooding angle in 2.2, and the sibling implementation",
        "",
        f"Omitting `θ_f` measures to 40°; supplying `θ_f = 35°` measures to 35° and moves the number",
        f"by {substituted['difference_m_rad']:.3f} m rad "
        f"({substituted['area_with_theta_f_m_rad']:.6f} against "
        f"{substituted['area_without_theta_f_m_rad']:.6f}) — the substitution is a computation, not a",
        "field that gets echoed. The two area clauses are the code's `2.2.1b` and `2.2.1c`; the",
        "sibling labels them `2.2.2` and `2.2.3`, which is the numbering divergence",
        "`docs/statutory.md` §6 records — the thresholds and verdicts are the same, a report that",
        "quotes a clause number is not. The same curve through `tensorcad.statutory`:",
        "",
        "| implementation | `A(0, 35°)` [m rad] | relative error vs antiderivative |",
        "|---|---|---|",
        f"| `tensorcad.intact_stability` (1° trapezoid) | {cross['sibling_m_rad']:.9f} | "
        f"{_sci(cross['sibling_relative_error'])} |",
        f"| `tensorcad.statutory` (adaptive Simpson) | {cross['statutory_m_rad']:.12f} | "
        f"{_sci(cross['statutory_relative_error'])} |",
        f"| antiderivative | {cross['reference_m_rad']:.12f} | — |",
        "",
        f"They agree to {_sci(cross['spread'])} relative — the sibling's `O(h²)` trapezoidal error —",
        "and the statutory run's own uncertainty",
        f"({_sci(cross['statutory_uncertainty_m_rad'])} m rad) is above its actual error. Both are read",
        "against the closed form, not against each other.",
        "",
        "### 6.7 The residuals are quadrature",
        "",
        "`b` against the antiderivative as the ladder is refined:",
        "",
        "| step [deg] | samples | `b` [m rad] | relative error | ratio to the next rung |",
        "|---|---|---|---|---|",
    ]
    for index, row in enumerate(convergence):
        ratio = "" if index == 0 else f"{ratios[index - 1]:.2f}"
        lines.append(
            f"| {row['step_deg']:.3f} | {row['samples']} | {row['area_b_m_rad']:.9f} | "
            f"{_sci(row['relative_error'])} | {ratio} |"
        )
    lines += [
        "",
        "Successive ratios " + ", ".join(f"{ratio:.2f}" for ratio in ratios) + ": the trapezoid",
        "signature from the first halving already, so the residual is the sampling of the curve and",
        f"not a disagreement about the rule. The 5° rung is inside {_sci(convergence[0]['relative_error'])}",
        "relative, which is why a coarse ladder still returns the right verdict on this barge — but",
        "that is not where the error bar should be read from.",
        "",
        "![validation](intact_stability/intact_stability_validation.png)",
        "",
        "*Top left: the barge's `GZ` with the steady lever `lw1` (`θ_0`) and the gust lever `lw2`",
        "(`θ_2`). Top right: `b` and `a` against `lw1`, crossing at the 2.3.2 boundary. Bottom left:",
        "the residual of `b` under grid refinement against an `h²` reference. Bottom right: the",
        "`0.3 sin(2 phi)` curve, where `θ_2` is `θ_c` found by the rule and re-solved by bisection.*",
        "",
        "### 6.8 Error judgement for this section",
        "",
        "| quantity | reference | measured error | < 3 %? |",
        "|---|---|---|---|",
    ]
    for quantity, reference, measured, value in (
        (
            "`b` over `[8°, 50°]`",
            "antiderivative of the barge's arm",
            f"{_sci(gust_b['relative_error'])} relative",
            gust_b["relative_error"],
        ),
        (
            "`a = lw2 (θ_2 − θ_1)`",
            "the lever's rectangle",
            f"{_sci(gust_a['relative_error'])} relative",
            gust_a["relative_error"],
        ),
        (
            "`θ_c`, second intercept",
            "independent bisection of `GZ = lw2`",
            f"{second['theta_c_error_deg']:.2e} deg",
            second["theta_c_error_deg"],
        ),
        (
            "`b` over `[0, θ_c]`",
            "antiderivative at `θ_c`",
            f"{_sci(second['area_relative_error'])} relative",
            second["area_relative_error"],
        ),
        (
            "area to `θ_f = 35°` (code `2.2.1b`), this module",
            "antiderivative at 35°, 1° grid",
            f"{_sci(limits['error_to_limit'])} relative",
            limits["error_to_limit"],
        ),
        (
            "same area, `tensorcad.statutory`",
            "antiderivative at 35°, adaptive Simpson",
            f"{_sci(cross['statutory_relative_error'])} relative",
            cross["statutory_relative_error"],
        ),
        (
            "spread between the two implementations",
            "one curve, one closed form",
            f"{_sci(cross['spread'])} relative",
            cross["spread"],
        ),
        (
            "2.3.2 margin against `b − a`",
            "the two areas the report itself quotes",
            f"{_sci(boundary['margin_relative_error'])} relative",
            boundary["margin_relative_error"],
        ),
        (
            "coarsest grid tested (5°)",
            "antiderivative of the gust area",
            f"{_sci(convergence[0]['relative_error'])} relative",
            convergence[0]["relative_error"],
        ),
    ):
        verdict = "✅" if value < 0.03 else "**NO**"
        lines.append(f"| {quantity} | {reference} | {measured} | {verdict} |")
    lines += [
        "",
        "**Every quantity in this section is inside 3 % of its reference, most by four to twelve",
        "orders of magnitude.** What is *not* covered, and is named rather than glossed: `θ_1`",
        "(clause 2.3.4's tables are absent from this repository, so it stays an input — its",
        "boundary is checked in §6.5, its value is not derived), the barge's areas above bilge",
        "emergence (the closed form stops applying there, as §2 shows with a negative control), and",
        "any real ship.",
        "",
        MARKDOWN_END,
    ]
    return "\n".join(lines) + "\n"


def splice_markdown(path: pathlib.Path, section: str) -> bool:
    """Replace the generated block in ``path`` with ``section``; return whether it changed."""
    text = path.read_text(encoding="utf-8")
    start = text.find(MARKDOWN_BEGIN)
    end = text.find(MARKDOWN_END)
    if start < 0 or end < 0 or end < start:
        raise SystemExit(f"{path}: generated-block markers not found")
    end += len(MARKDOWN_END)
    updated = text[:start] + section.rstrip("\n") + text[end:]
    if updated == text:
        return False
    path.write_text(updated, encoding="utf-8")
    return True


# --------------------------------------------------------------------------------------
# figure
# --------------------------------------------------------------------------------------
def figure(path: pathlib.Path, records: list[dict[str, object]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    by_name = {record["name"]: record for record in records}
    wall = ladder(1.0, 50.0, wall_sided_gz)
    report = weather(wall, 0.02, roll_amplitude_deg=8.0)
    turning = ladder(0.1, 60.0, turnover_gz)

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0))
    axes = axes.ravel()

    ax = axes[0]
    angles = [angle for angle in wall[0]]
    ax.plot(angles, wall[1], color="#1f4e9c", lw=1.8, label="GZ (wall-sided barge)")
    ax.axhline(report.levers.steady_m, color="#c0392b", ls="--", lw=1.2, label="lw1 (steady)")
    ax.axhline(report.levers.gust_m, color="#e67e22", ls=":", lw=1.4, label="lw2 = 1.5 lw1")
    ax.axvline(report.theta0_deg, color="#c0392b", lw=1.0, alpha=0.6)
    ax.axvline(report.theta2_deg, color="#e67e22", lw=1.0, alpha=0.6)
    ax.annotate(
        f"$\\theta_0$ = {report.theta0_deg:.2f}$^\\circ$",
        (report.theta0_deg, report.levers.steady_m),
        textcoords="offset points",
        xytext=(6, -14),
        fontsize=8,
        color="#c0392b",
    )
    ax.annotate(
        f"$\\theta_2$ = {report.theta2_deg:.1f}$^\\circ$",
        (report.theta2_deg, report.levers.gust_m),
        textcoords="offset points",
        xytext=(-58, 6),
        fontsize=8,
        color="#e67e22",
    )
    ax.set_xlabel("heel angle [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title("2.3.1 / 2.3.2 on the reference barge")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(alpha=0.3)

    ax = axes[1]
    steady = [0.15 + 0.08 * index / 99 for index in range(100)]
    b_values, a_values = [], []
    for lever in steady:
        item = weather(wall, lever)
        b_values.append(item.area_b_m_rad)
        a_values.append(item.area_a_m_rad)
    crossing = wall_sided_area(50.0) / (
        iscode.WIND_GUST_FACTOR * math.radians(50.0)
    )  # the analytic lever where b = a: what claim_gust_boundary checks
    ax.plot(steady, b_values, color="#1f4e9c", lw=1.8, label="b = $\\int$GZ")
    ax.plot(steady, a_values, color="#c0392b", lw=1.8, label="a = lw2 $(\\theta_2-\\theta_1)$")
    ax.axvline(crossing, color="0.35", ls="--", lw=1.0)
    ax.annotate(
        f"b = a at lw1 = {crossing:.4f} m",
        (crossing, max(b_values)),
        textcoords="offset points",
        xytext=(6, -12),
        fontsize=8,
    )
    ax.set_xlabel("steady wind lever lw1 [m]")
    ax.set_ylabel("area [m rad]")
    ax.set_title("2.3.2 boundary is where the two areas cross")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[2]
    rows = by_name["grid_convergence"]["rows"]
    steps = [row["step_deg"] for row in rows]
    errors = [row["relative_error"] for row in rows]
    ax.loglog(steps, errors, "o-", color="#1f4e9c", lw=1.6, label="|b - $\\int$GZ| / $\\int$GZ")
    reference_line = [errors[0] * (step / steps[0]) ** 2 for step in steps]
    ax.loglog(steps, reference_line, "k--", lw=1.0, label="$h^2$ reference")
    ax.set_xlabel("angle step [deg]")
    ax.set_ylabel("relative error of b")
    ax.set_title("gust-area residual is quadrature")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[3]
    ax.plot(turning[0], turning[1], color="#1f4e9c", lw=1.8, label="GZ = 0.3 sin(2$\\phi$)")
    turning_report = weather(turning, 0.298 / iscode.WIND_GUST_FACTOR)
    ax.axhline(
        turning_report.levers.gust_m, color="#e67e22", ls=":", lw=1.4, label="lw2"
    )
    ax.axvline(
        turning_report.theta_c_deg, color="#c0392b", lw=1.0, alpha=0.7, label="$\\theta_c$"
    )
    span = [angle for angle in turning[0] if angle <= turning_report.theta2_deg]
    ax.fill_between(
        span,
        [turnover_gz(angle) for angle in span],
        color="#1f4e9c",
        alpha=0.15,
        label="b",
    )
    ax.fill_between(
        span,
        [turning_report.levers.gust_m] * len(span),
        color="#c0392b",
        alpha=0.12,
        label="a",
    )
    ax.annotate(
        f"$\\theta_c$ = {turning_report.theta_c_deg:.2f}$^\\circ$",
        (turning_report.theta_c_deg, turning_report.levers.gust_m),
        textcoords="offset points",
        xytext=(-70, 8),
        fontsize=8,
        color="#c0392b",
    )
    ax.set_xlabel("heel angle [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title("second intercept found by the rule, re-solved by bisection")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.3)

    fig.suptitle(
        "tensorcad.intact_stability -- weather criterion validated against closed forms",
        fontsize=11,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path(__file__).resolve().parents[1]
        / "docs"
        / "intact_stability",
        help="directory for the json and png artifacts",
    )
    parser.add_argument("--no-figure", action="store_true", help="skip the png")
    parser.add_argument(
        "--update-docs",
        type=pathlib.Path,
        default=None,
        help="splice the generated section into this markdown file (docs/intact_stability.md)",
    )
    options = parser.parse_args(argv)

    started = time.time()
    records = evaluate()
    failures = [record["name"] for record in records if not record["holds"]]
    payload = {
        "schema": SCHEMA,
        "tensorcad_version": tensorcad.__version__,
        "runtime_s": round(time.time() - started, 3),
        "reference": {
            "body": "wall-sided box barge, closed in the bilge; second curve GZ = 0.3 sin(2 phi)",
            "breadth_m": BREADTH,
            "draft_m": DRAFT,
            "depth_m": DEPTH,
            "length_m": LENGTH,
            "kg_m": KG,
            "bm0_m": BM0,
            "gm0_m": GM0,
            "wind_arm_m": WIND_ARM_M,
            "displacement_t": DISPLACEMENT_T,
            "gravity_m_s2": GRAVITY,
            "wind_pressure_pa": iscode.WIND_PRESSURE_PA,
            "gust_factor": iscode.WIND_GUST_FACTOR,
            "theta2_cap_deg": iscode.WEATHER_THETA2_CAP_DEG,
        },
        "claims": records,
        "failures": failures,
        "all_claims_hold": not failures,
    }
    options.output.mkdir(parents=True, exist_ok=True)
    (options.output / "intact_stability_validation.json").write_text(
        json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8"
    )
    if not options.no_figure:
        figure(options.output / "intact_stability_validation.png", records)
    if options.update_docs is not None:
        changed = splice_markdown(options.update_docs, render_section(payload))
        print(f"{options.update_docs}: {'updated' if changed else 'already current'}")

    for record in records:
        print(f"[{'ok' if record['holds'] else 'FAIL'}] {record['name']}: {record['measured']}")
    if failures:
        print(f"FAILED claims: {failures}")
        return 1
    print(f"all {len(records)} claims hold")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())