"""Section 2.2 of the IS Code 2008, checked against closed-form righting arms.

The reference values here are derived, not quoted:

* a curve ``GZ = GM sin(phi)`` integrates to ``GM (1 - cos(phi))``;
* the wall-sided barge ``GZ = sin(phi) (GM0 + BM0/2 tan^2 phi)`` integrates to
  ``GM0 (1 - cos) + BM0/2 (sec + cos - 2)``, which follows from
  ``int sin(phi) tan^2(phi) dphi = sec(phi) + cos(phi)``.

Both are reproduced to round-off, which is what makes the quadrature in
:mod:`tensorcad.statutory` trustworthy on curves that have no closed form. The
non-vacuity of each check is tested by a curve that fails it.
"""

from __future__ import annotations

import json
import math
from types import SimpleNamespace

import pytest

import tensorcad.mesh as meshlib
from tensorcad import (
    INTACT_STABILITY_CRITERIA,
    CriterionCheck,
    IntactStabilityAssessment,
    StatutoryError,
    area_under_gz_curve,
    box_righting_arms_analytic,
    check_intact_stability,
    gm0_from_gz_curve,
    righting_arm_curve,
)
from tensorcad.geometry import Mesh, MeshBuilder

# --- closed forms ---------------------------------------------------------------------

SINE_GM = 0.5


def sine_curve(angle_deg: float) -> float:
    return SINE_GM * math.sin(math.radians(angle_deg))


def sine_area_m_rad(stop_deg: float, gm: float = SINE_GM) -> float:
    """``int_0^stop GM sin(phi) dphi`` in m*rad."""
    return gm * (1.0 - math.cos(math.radians(stop_deg)))


BARGE_BREADTH = 6.0
BARGE_DRAFT = 3.0
BARGE_DEPTH = 8.0
BARGE_KG = 2.0
BARGE_BM = BARGE_BREADTH**2 / (12.0 * BARGE_DRAFT)
BARGE_GM = 0.5 * BARGE_DRAFT + BARGE_BM - BARGE_KG


def wall_sided_curve(angle_deg: float, gm: float = BARGE_GM, bm: float = BARGE_BM) -> float:
    """Classical wall-sided righting arm, exact for a box up to the bilge/deck limit."""
    tangent = math.tan(math.radians(angle_deg))
    return math.sin(math.radians(angle_deg)) * (gm + 0.5 * bm * tangent * tangent)


def wall_sided_area_m_rad(stop_deg: float, gm: float = BARGE_GM, bm: float = BARGE_BM) -> float:
    """``int_0^stop`` of :func:`wall_sided_curve`, in m*rad."""
    radians = math.radians(stop_deg)
    return gm * (1.0 - math.cos(radians)) + 0.5 * bm * (
        (1.0 / math.cos(radians) + math.cos(radians)) - 2.0
    )


def box(minimum: tuple[float, float, float], maximum: tuple[float, float, float]) -> Mesh:
    """Closed, outward-oriented axis-aligned box built from 12 raw triangles."""
    x0, y0, z0 = minimum
    x1, y1, z1 = maximum
    builder = MeshBuilder()
    corners = [
        (x0, y0, z0),
        (x1, y0, z0),
        (x1, y1, z0),
        (x0, y1, z0),
        (x0, y0, z1),
        (x1, y0, z1),
        (x1, y1, z1),
        (x0, y1, z1),
    ]
    ids = [builder.add_vertex(point) for point in corners]
    for face in (
        (0, 3, 2),
        (0, 2, 1),
        (4, 5, 6),
        (4, 6, 7),
        (0, 1, 5),
        (0, 5, 4),
        (2, 3, 7),
        (2, 7, 6),
        (1, 2, 6),
        (1, 6, 5),
        (0, 4, 7),
        (0, 7, 3),
    ):
        builder.add_triangle(*[ids[index] for index in face])
    return meshlib.orient_outward(builder.build()).mesh


# --- quadrature and curve readings ----------------------------------------------------


def test_sine_curve_areas_match_the_closed_form() -> None:
    for start, stop in ((0.0, 30.0), (0.0, 40.0), (30.0, 40.0), (0.0, 60.0)):
        expected = sine_area_m_rad(stop) - sine_area_m_rad(start)
        assert area_under_gz_curve(sine_curve, start, stop) == pytest.approx(expected, abs=1e-13)


def test_wall_sided_curve_areas_match_the_closed_form() -> None:
    """The tan^2 term is what a trapezoidal rule gets wrong; the adaptive rule must not."""
    for stop in (30.0, 40.0):
        assert area_under_gz_curve(wall_sided_curve, 0.0, stop) == pytest.approx(
            wall_sided_area_m_rad(stop), abs=1e-13
        )
    assert area_under_gz_curve(wall_sided_curve, 30.0, 40.0) == pytest.approx(
        wall_sided_area_m_rad(40.0) - wall_sided_area_m_rad(30.0), abs=1e-13
    )


def test_gm0_from_the_slope_is_exact_for_a_sine_curve() -> None:
    assert gm0_from_gz_curve(sine_curve) == pytest.approx(SINE_GM, abs=1e-12)


def test_gm0_from_the_slope_cancels_the_cubic_term() -> None:
    """Richardson on 1 deg / 2 deg removes the phi^2 term of the wall-sided expansion."""
    assert gm0_from_gz_curve(wall_sided_curve) == pytest.approx(BARGE_GM, abs=1e-6)


# --- the six rules --------------------------------------------------------------------


def test_all_six_criteria_pass_on_a_compliant_curve() -> None:
    assessment = check_intact_stability(sine_curve)
    assert assessment.verdict() == "pass"
    assert len(assessment.checks) == len(INTACT_STABILITY_CRITERIA) == 6
    assert assessment.check("2.2.1a").attained == pytest.approx(sine_area_m_rad(30.0), abs=1e-13)
    assert assessment.check("2.2.1b").attained == pytest.approx(sine_area_m_rad(40.0), abs=1e-13)
    assert assessment.check("2.2.1c").attained == pytest.approx(
        sine_area_m_rad(40.0) - sine_area_m_rad(30.0), abs=1e-13
    )
    assert assessment.check("2.2.2").attained == pytest.approx(sine_curve(30.0), abs=1e-15)
    assert assessment.check("2.2.4").attained == pytest.approx(SINE_GM, abs=1e-12)
    # sin rises all the way to 90 deg, so the peak sits at the edge of the domain and the
    # assessment says so instead of pretending the maximum was located.
    assert assessment.angle_of_max_gz_deg == pytest.approx(90.0)
    assert "edge" in assessment.check("2.2.3").note


def test_shape_of_each_threshold_is_the_code() -> None:
    thresholds = {code: spec["threshold"] for code, spec in INTACT_STABILITY_CRITERIA.items()}
    assert thresholds == {
        "2.2.1a": 0.055,
        "2.2.1b": 0.090,
        "2.2.1c": 0.030,
        "2.2.2": 0.20,
        "2.2.3": 25.0,
        "2.2.4": 0.15,
    }


def test_a_marginal_curve_fails_the_rules_it_should_fail() -> None:
    """Negative control: a tender curve must not slip through on round-off."""
    assessment = check_intact_stability(lambda angle: 0.1 * math.sin(math.radians(angle)))
    assert assessment.verdict() == "fail"
    assert {item.code for item in assessment.failed} == {"2.2.1a", "2.2.1b", "2.2.1c", "2.2.2", "2.2.4"}
    assert assessment.not_assessed == ()
    assert assessment.indeterminate == ()


def test_a_value_equal_to_the_threshold_passes() -> None:
    """The rules read "at least": equality is compliance, and round-off must not undo it."""

    def exactly_at_the_limit(angle_deg: float) -> float:
        # GZ(30) = 0.2 exactly by construction, with the area still above 0.055.
        sine = math.sin(math.radians(angle_deg))
        return 0.44 * sine - 0.16 * sine**3

    assessment = check_intact_stability(exactly_at_the_limit)
    assert assessment.check("2.2.2").attained == pytest.approx(0.20, abs=1e-15)
    assert assessment.check("2.2.2").status == "pass"
    assert assessment.verdict() == "pass"


def test_a_threshold_touched_inside_the_interpolation_error_is_indeterminate() -> None:
    """Same curve, read off a ladder that has no sample at 30 deg: the answer is "cannot tell"."""

    def exactly_at_the_limit(angle_deg: float) -> float:
        sine = math.sin(math.radians(angle_deg))
        return 0.44 * sine - 0.16 * sine**3

    angles = [0.0, 10.0, 20.0, 28.0, 36.0, 44.0, 52.0]
    assessment = check_intact_stability((angles, [exactly_at_the_limit(angle) for angle in angles]))
    item = assessment.check("2.2.2")
    assert item.uncertainty > 1e-6
    assert item.status == "indeterminate"
    assert assessment.verdict() == "undetermined"
    assert assessment.passed is None


# --- flooding angle, vanishing angle, short curves ------------------------------------


def test_the_flooding_angle_replaces_the_40_degree_limit() -> None:
    assessment = check_intact_stability(wall_sided_curve, flooding_angle_deg=25.0)
    assert assessment.area_limit_deg == 25.0
    assert assessment.check("2.2.1b").attained == pytest.approx(
        wall_sided_area_m_rad(25.0), abs=1e-12
    )
    # The ship cannot be heeled to 30 deg intact, so 2.2.1a and 2.2.2 have no answer.
    assert assessment.check("2.2.1a").status == "not_assessed"
    assert assessment.check("2.2.2").status == "not_assessed"
    assert "flooding angle" in assessment.check("2.2.2").note
    assert assessment.verdict() == "fail"  # 0.040 m*rad to 25 deg is genuinely below 0.09


def test_without_a_flooding_angle_the_check_says_so() -> None:
    assessment = check_intact_stability(sine_curve)
    assert assessment.flooding_angle_deg is None
    assert assessment.area_limit_deg == 40.0
    assert "flooding angle not supplied" in assessment.check("2.2.1b").note


def test_a_curve_that_vanishes_before_40_deg_is_not_signed_off() -> None:
    def capsizing(angle_deg: float) -> float:
        if angle_deg < 35.0:
            return 0.5 * math.sin(math.radians(angle_deg))
        return -0.01 * (angle_deg - 35.0)

    assessment = check_intact_stability(capsizing)
    assert assessment.vanishing_angle_deg == pytest.approx(35.0, abs=1e-6)
    assert assessment.check("2.2.1b").status == "not_assessed"
    assert assessment.check("2.2.1c").status == "not_assessed"
    assert "vanishes" in assessment.check("2.2.1b").note
    assert assessment.verdict() == "undetermined"


def test_a_short_curve_cannot_answer_the_rules_it_cannot_reach() -> None:
    angles = [0.0, 5.0, 10.0, 15.0, 20.0]
    levers = [0.02 * angle for angle in angles]
    assessment = check_intact_stability((angles, levers))
    for code in ("2.2.1a", "2.2.1b", "2.2.1c", "2.2.2", "2.2.3"):
        assert assessment.check(code).status == "not_assessed", code
    assert "stops at 20.000 deg" in assessment.check("2.2.1a").note
    assert assessment.check("2.2.4").status == "pass"  # the slope at upright is measurable
    assert assessment.verdict() == "undetermined"


# --- tables versus callables ----------------------------------------------------------


def test_a_dense_table_reproduces_the_callable_curve() -> None:
    angles = [float(angle) for angle in range(0, 91, 5)]
    from_callable = check_intact_stability(wall_sided_curve)
    from_table = check_intact_stability((angles, [wall_sided_curve(angle) for angle in angles]))
    assert from_table.verdict() == from_callable.verdict() == "pass"
    assert from_table.area_to_30_m_rad == pytest.approx(from_callable.area_to_30_m_rad, abs=3e-4)
    assert from_table.area_30_to_limit_m_rad == pytest.approx(
        from_callable.area_30_to_limit_m_rad, abs=3e-4
    )
    assert from_table.gm0_m == pytest.approx(BARGE_GM, abs=5e-3)
    assert from_table.curve_kind == "table"


def test_a_coarse_table_reports_a_bigger_error_bars() -> None:
    """The error bar has to grow when the sampling coarsens, or it is decoration."""
    dense_angles = [float(angle) for angle in range(0, 91, 5)]
    coarse_angles = [float(angle) for angle in range(0, 91, 10)]
    dense = check_intact_stability((dense_angles, [wall_sided_curve(a) for a in dense_angles]))
    coarse = check_intact_stability((coarse_angles, [wall_sided_curve(a) for a in coarse_angles]))
    assert coarse.area_uncertainty_m_rad > dense.area_uncertainty_m_rad
    assert coarse.check("2.2.4").uncertainty > dense.check("2.2.4").uncertainty


def test_both_curve_objects_of_the_package_are_accepted() -> None:
    """``gz.RightingArmCurve`` exposes angles_deg/gz, ``righting_arm`` exposes points."""
    angles = [float(angle) for angle in range(0, 91, 5)]
    levers = [sine_curve(angle) for angle in angles]
    like_gz = SimpleNamespace(angles_deg=tuple(angles), gz=tuple(levers))
    like_bisection = SimpleNamespace(
        points=tuple(
            SimpleNamespace(heel_deg=angle, righting_arm=lever)
            for angle, lever in zip(angles, levers, strict=True)
        )
    )
    reference = check_intact_stability(sine_curve).verdict()
    for curve in (like_gz, like_bisection):
        assessment = check_intact_stability(curve)
        assert assessment.verdict() == reference


def test_the_mesh_pipeline_feeds_the_rule_check() -> None:
    """mesh -> ``gz.righting_arm_curve`` -> IS Code verdict, against the closed form."""
    barge = box((0.0, -BARGE_BREADTH / 2, 0.0), (10.0, BARGE_BREADTH / 2, BARGE_DEPTH))
    angles = [0.0, 10.0, 20.0, 30.0, 40.0]
    curve = righting_arm_curve(barge, BARGE_KG, angles, waterline_z=BARGE_DRAFT)
    analytic = box_righting_arms_analytic(
        BARGE_BREADTH, BARGE_DRAFT, BARGE_KG, angles, depth=BARGE_DEPTH
    )
    for gz, expected in zip(curve.gz, analytic, strict=True):
        assert gz == pytest.approx(expected, abs=1e-13)

    assessment = check_intact_stability(curve)
    assert assessment.verdict() == "pass"
    assert assessment.curve_kind == "table"
    assert assessment.area_to_30_m_rad == pytest.approx(wall_sided_area_m_rad(30.0), abs=2e-3)
    assert assessment.area_to_limit_m_rad == pytest.approx(wall_sided_area_m_rad(40.0), abs=3e-3)
    assert assessment.gm0_m == pytest.approx(BARGE_GM, abs=1e-2)


# --- refusals -------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("curve", "message"),
    [
        (([0.0, 30.0, 20.0], [0.0, 0.3, 0.2]), "strictly increasing"),
        (([0.0, 10.0, 10.0], [0.0, 0.2, 0.2]), "strictly increasing"),
        (([0.0], [0.0]), "at least two samples"),
        (([0.0, 10.0, 20.0], [0.0, 0.2]), "angles against"),
        (([0.0, float("nan")], [0.0, 0.2]), "non-finite"),
        (object(), "a GZ curve is a callable"),
    ],
)
def test_unreadable_curves_are_refused(curve: object, message: str) -> None:
    with pytest.raises(StatutoryError, match=message):
        check_intact_stability(curve)


def test_bad_parameters_are_refused() -> None:
    with pytest.raises(StatutoryError, match="gm0 must be finite and positive"):
        check_intact_stability(sine_curve, gm0=0.0)
    with pytest.raises(StatutoryError, match="flooding_angle_deg must be finite and positive"):
        check_intact_stability(sine_curve, flooding_angle_deg=-5.0)
    with pytest.raises(StatutoryError, match="empty heel domain"):
        check_intact_stability(sine_curve, heel_domain_deg=(30.0, 30.0))
    with pytest.raises(StatutoryError, match="empty heel interval"):
        area_under_gz_curve(sine_curve, 30.0, 30.0)
    with pytest.raises(StatutoryError, match="outside the curve"):
        area_under_gz_curve(([0.0, 20.0], [0.0, 0.1]), 0.0, 30.0)
    with pytest.raises(StatutoryError, match="cannot be read from its slope"):
        gm0_from_gz_curve(([5.0, 10.0, 20.0], [0.1, 0.2, 0.3]))


def test_a_non_finite_ordinate_is_refused() -> None:
    def broken(angle_deg: float) -> float:
        return float("nan") if angle_deg > 20.0 else 0.1

    with pytest.raises(StatutoryError, match="not finite"):
        check_intact_stability(broken)


# --- reporting ------------------------------------------------------------------------


def test_the_assessment_serialises_to_json() -> None:
    assessment = check_intact_stability(sine_curve, gm0=0.42, flooding_angle_deg=50.0)
    payload = assessment.to_dict()
    assert payload["schema"] == "tensorcad.intact_stability/1"
    assert payload["verdict"] == "pass"
    assert payload["gm0_source"] == "supplied"
    assert payload["gm0_m"] == pytest.approx(0.42)
    assert len(payload["checks"]) == 6
    assert "weather criterion" in payload["scope"]
    assert json.loads(json.dumps(payload)) == payload


def test_an_explicit_gm0_overrides_the_slope_estimate() -> None:
    assessment = check_intact_stability(sine_curve, gm0=0.16)
    assert assessment.gm0_source == "supplied"
    assert assessment.check("2.2.4").attained == pytest.approx(0.16)
    assert assessment.check("2.2.4").status == "pass"


def test_the_three_valued_verdict_is_not_a_soft_boole() -> None:
    """``passed`` is ``None`` -- not ``False`` -- when the answer is undecidable."""
    assert CriterionCheck("2.2.2", "GZ at 30 deg", 0.2, "m", 0.2000001, "indeterminate").passed is False
    undecided = IntactStabilityAssessment(
        checks=(CriterionCheck("2.2.2", "GZ at 30 deg", 0.2, "m", 0.21, "indeterminate"),)
    )
    assert undecided.passed is None
    assert undecided.verdict() == "undetermined"
    failing = IntactStabilityAssessment(
        checks=(CriterionCheck("2.2.2", "GZ at 30 deg", 0.2, "m", 0.10, "fail"),)
    )
    assert failing.passed is False
    assert failing.verdict() == "fail"


# --- cross-check against tensorcad.intact_stability ------------------------------------
#
# ``tensorcad.intact_stability`` (merged from PR #229/#230) checks the same section of the
# same code by a different route: it reports the areas by plain trapezoidal integration on
# the caller's grid and treats the caller's ``gm0`` as an input, where this module integrates
# a Hermite interpolant / an adaptive quadrature and *derives* GM0 from the curve. Two
# independent routes to one number are this package's standing way of validating a number
# (``docs/m4038_pr40_merge_cure.md`` does the same for ``GZ``), so the two are put on the
# same curves here and the agreement is recorded rather than assumed.

import tensorcad.intact_stability as iscode  # noqa: E402 - grouped with the cross-check

#: statutory code -> the sibling module's flattened clause label, matched on the quantity
#: each rule demands. The sibling renumbers 2.2.1a/b/c as 2.2.1-2.2.3 and shifts the rest up
#: by one; the thresholds are the same numbers either way, which is what this mapping asserts.
CLAUSE_CROSSWALK = {
    "2.2.1a": "2.2.1",
    "2.2.1b": "2.2.2",
    "2.2.1c": "2.2.3",
    "2.2.2": "2.2.4",
    "2.2.3": "2.2.5",
    "2.2.4": "2.2.6",
}


def test_the_two_modules_hold_the_same_thresholds() -> None:
    """Section 2.2 as data, read from both modules: the same six numbers in the same units."""
    assert set(CLAUSE_CROSSWALK) == set(INTACT_STABILITY_CRITERIA)
    for mine, theirs in CLAUSE_CROSSWALK.items():
        threshold = INTACT_STABILITY_CRITERIA[mine]["threshold"]
        required, unit, _ = iscode.IS_CODE_2008_GENERAL[theirs]
        assert threshold == pytest.approx(required, abs=0.0), mine
        assert INTACT_STABILITY_CRITERIA[mine]["unit"].replace("*", " ") == unit, mine


def test_areas_agree_with_the_sibling_within_its_own_integration_error() -> None:
    """Hermite-adaptive vs the sibling's trapezoidal areas on a 0.25 deg ladder.

    Both are compared to the closed form rather than to each other, so the test says *which*
    number is wrong when they disagree: this module's area is at round-off from
    ``GM0(1-cos) + BM0/2(sec+cos-2)``, and the sibling's trapezoidal area carries the O(h^2)
    discretisation error of the 0.25 deg ladder it is given. The two therefore agree only to
    that error, which is the honest claim -- a rule check compared at three decimals wants
    the first number, and the spread is what the sibling's docstring calls the grid's error.
    """
    angles = tuple(i * 0.25 for i in range(0, 4 * 41))
    levers = tuple(wall_sided_curve(angle) for angle in angles)
    for stop in (30.0, 40.0):
        closed = wall_sided_area_m_rad(stop)
        mine = area_under_gz_curve((angles, levers), 0.0, stop)
        theirs = iscode.gz_area_between(angles, levers, 0.0, stop)
        assert mine == pytest.approx(closed, rel=1e-10), (stop, mine, closed)
        # 2e-5 is the O(h^2) error of a 0.25 deg trapezoidal ladder at 40 deg (measured:
        # 1.1e-5); the Hermite/adaptive area above sits ~7 orders below it.
        assert theirs == pytest.approx(closed, rel=2e-5), (stop, theirs, closed)
        assert mine == pytest.approx(theirs, rel=2e-5)


def test_the_verdict_agrees_with_the_sibling_on_pass_and_fail() -> None:
    """A passing and a failing curve are called the same way by both modules."""
    angles = tuple(float(i) for i in range(0, 41))
    gm0 = gm0_from_gz_curve((angles, tuple(sine_curve(a) for a in angles)))

    passing = check_intact_stability((angles, tuple(sine_curve(a) for a in angles)), gm0=gm0)
    sibling = iscode.general_criteria((angles, tuple(sine_curve(a) for a in angles)), gm0_m=gm0)
    assert passing.verdict() == "pass"
    assert sibling.failures == ()

    weak = tuple(0.15 * math.sin(math.radians(a)) for a in angles)  # GZ(30) = 0.075 < 0.20
    failing = check_intact_stability((angles, weak), gm0=0.15)
    sibling_weak = iscode.general_criteria((angles, weak), gm0_m=0.15)
    assert failing.verdict() == "fail"
    # the same failed rules, once the clause labels are crosswalked: statutory 2.2.1a/b/c +
    # 2.2.2 are the sibling's 2.2.1-2.2.4, and a GM0 of 0.15 m = the 2.2.4 threshold, so it
    # is the four area/GZ rules that fail, not the metacentric height.
    assert {item.code for item in failing.failed} == {"2.2.1a", "2.2.1b", "2.2.1c", "2.2.2"}
    assert {CLAUSE_CROSSWALK[item.code] for item in failing.failed} == {
        item.clause for item in sibling_weak.failures
    }
