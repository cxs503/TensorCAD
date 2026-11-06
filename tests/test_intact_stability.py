"""IS Code 2008 part A criteria on a measured righting-arm curve.

The module is arithmetic on the polyline ``(angle, GZ)``, so its reference is a shape whose
righting arm is known in closed form: **the wall-sided box barge**, where
``GZ(phi) = sin(phi) (GM0 + (BM0/2) tan^2(phi))`` and therefore

    int GZ dphi = GM0 (1 - cos phi) + (BM0/2) (sec phi + cos phi - 2)

has an antiderivative in elementary functions. That gives the tests three separate things to
pin, and they are not the same check:

* the **area* the criteria are quoted from** -- against that antiderivative, on a
  refinement ladder (so a wrong quadrature order cannot pass), and with the additivity the
  code's four areas rely on (``A(0,40) = A(0,30) + A(30,40)``);
* the **chain** -- a box mesh through :func:`tensorcad.gz.righting_arm_curve` and through
  :func:`tensorcad.righting_arm.righting_arm_curve` (two independent waterline solvers)
  into the criteria, both landing on the analytic report, which is what makes the numbers in
  ``docs/intact_stability.md`` a measurement rather than a restatement;
* the **comparison's ability to fail** -- negative controls for every clause, because a
  criterion that always returns "satisfied" is indistinguishable from a criterion that
  works. The weather criterion gets the same treatment at the top of the file's second half:
  the steady heel angle is checked against a bisection solve of the closed form, and the
  gust areas against the antiderivative over the same span.
"""

from __future__ import annotations

import json
import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.gz import righting_arm_curve as gz_curve
from tensorcad.intact_stability import (
    IS_CODE_2008_GENERAL,
    Criterion,
    IntactStabilityError,
    general_criteria,
    gz_area_between,
    gz_at_angle,
    steady_wind_heel_angle,
    weather_criterion,
    wind_heeling_levers,
)
from tensorcad.righting_arm import righting_arm_curve as bisection_curve

LENGTH, BREADTH, DRAFT, DEPTH = 10.0, 4.0, 2.0, 4.0
KG = 1.2
BM0 = BREADTH * BREADTH / (12.0 * DRAFT)
KB0 = 0.5 * DRAFT
GM0 = KB0 + BM0 - KG
DISPLACEMENT_T = LENGTH * BREADTH * DRAFT * 1.025  # seawater, 1.025 t/m^3


def box(minimum: tuple[float, float, float], maximum: tuple[float, float, float]) -> Mesh:
    """Closed, outward-oriented axis-aligned box (same construction as ``test_gz``)."""
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


@pytest.fixture(scope="module")
def barge() -> Mesh:
    """Barge with 2 m of freeboard: keel at z = 0, deck at z = 4, floating at z = 2."""
    return box((0.0, -BREADTH / 2, 0.0), (LENGTH, BREADTH / 2, DEPTH))


def wall_sided_gz(angle_deg: float, *, gm: float = GM0, bm: float = BM0) -> float:
    heel = math.radians(angle_deg)
    return math.sin(heel) * (gm + 0.5 * bm * math.tan(heel) ** 2)


def wall_sided_area(lower_deg: float, upper_deg: float, *, gm: float = GM0, bm: float = BM0) -> float:
    """``int GZ dphi`` from the antiderivative, in m*rad."""

    def primitive(angle_deg: float) -> float:
        heel = math.radians(angle_deg)
        return gm * (1.0 - math.cos(heel)) + 0.5 * bm * (1.0 / math.cos(heel) + math.cos(heel) - 2.0)

    return primitive(upper_deg) - primitive(lower_deg)


def ladder(step_deg: float, limit_deg: float = 45.0) -> tuple[float, ...]:
    count = int(round(limit_deg / step_deg))
    return tuple(index * step_deg for index in range(count + 1))


def analytic_curve(step_deg: float, limit_deg: float = 45.0) -> tuple[tuple[float, ...], tuple[float, ...]]:
    angles = ladder(step_deg, limit_deg)
    return angles, tuple(wall_sided_gz(angle) for angle in angles)


# ---------------------------------------------------------------------------
# the closed form the areas are read against
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("step_deg", [5.0, 1.0, 0.25, 0.05])
def test_area_matches_the_antiderivative_and_converges(step_deg: float) -> None:
    angles, values = analytic_curve(step_deg)
    measured = gz_area_between(angles, values, 0.0, 40.0)
    expected = wall_sided_area(0.0, 40.0)
    relative = abs(measured - expected) / expected
    # trapezoid on a curved integrand is O(h^2): 3.3e-3 at 5 deg, 3.3e-7 at 0.05 deg
    if step_deg == 5.0:
        assert relative < 5e-3
    if step_deg == 0.05:
        assert relative < 1e-6
    assert abs(measured - expected) < relative * expected + 1e-15


def test_area_error_falls_as_the_grid_squared() -> None:
    coarse = abs(gz_area_between(*analytic_curve(4.0), 0.0, 40.0) - wall_sided_area(0.0, 40.0))
    fine = abs(gz_area_between(*analytic_curve(2.0), 0.0, 40.0) - wall_sided_area(0.0, 40.0))
    assert coarse > 0.0 and fine > 0.0
    assert 3.4 < coarse / fine < 4.6  # halving h quarters the error


def test_areas_are_additive_across_the_clause_boundaries() -> None:
    angles, values = analytic_curve(0.25)
    total = gz_area_between(angles, values, 0.0, 40.0)
    split = gz_area_between(angles, values, 0.0, 30.0) + gz_area_between(angles, values, 30.0, 40.0)
    assert abs(total - split) <= 1e-15 * max(1.0, abs(total))


def test_area_handles_a_single_segment_and_limit_inside_it() -> None:
    angles, values = (0.0, 40.0), (0.0, 0.2)
    assert gz_area_between(angles, values, 10.0, 20.0) == pytest.approx(
        0.5 * (0.05 + 0.1) * math.radians(10.0)
    )
    assert gz_area_between(angles, values, 20.0, 20.0) == 0.0


def test_area_refuses_to_extrapolate() -> None:
    angles, values = (0.0, 30.0), (0.0, 0.3)
    with pytest.raises(IntactStabilityError, match="not inside the sampled range"):
        gz_area_between(angles, values, 0.0, 40.0)


def test_interpolation_is_exact_on_samples_and_linear_between() -> None:
    angles, values = (0.0, 10.0, 20.0), (0.0, 0.1, 0.25)
    assert gz_at_angle(angles, values, 10.0) == 0.1
    assert gz_at_angle(angles, values, 5.0) == pytest.approx(0.05)
    assert gz_at_angle(angles, values, 15.0) == pytest.approx(0.175)
    assert gz_at_angle(angles, values, 0.0) == 0.0
    assert gz_at_angle(angles, values, 20.0) == 0.25
    with pytest.raises(IntactStabilityError, match="outside the sampled range"):
        gz_at_angle(angles, values, 20.5)


# ---------------------------------------------------------------------------
# the criteria on the analytic curve
# ---------------------------------------------------------------------------
def test_general_criteria_on_the_analytic_curve_pass_and_match_the_closed_form() -> None:
    report = general_criteria(analytic_curve(0.25), gm0_m=GM0)
    assert report.satisfied
    assert "passes" in report.verdict()
    assert report.failures == ()
    assert report.area_0_30_m_rad == pytest.approx(wall_sided_area(0.0, 30.0), rel=1e-5)
    assert report.area_0_40_m_rad == pytest.approx(wall_sided_area(0.0, 40.0), rel=1e-5)
    assert report.area_30_40_m_rad == pytest.approx(wall_sided_area(30.0, 40.0), rel=1e-5)
    assert report.gz_30_m == pytest.approx(wall_sided_gz(30.0), abs=1e-12)
    assert report.max_gz_m == pytest.approx(wall_sided_gz(45.0), abs=1e-12)
    assert report.angle_of_max_gz_deg == pytest.approx(45.0)
    assert report.upper_angle_deg == 40.0
    attained = {item.clause: item.attained for item in report.criteria}
    assert attained["2.2.1"] == pytest.approx(wall_sided_area(0.0, 30.0), rel=1e-5)
    assert attained["2.2.6"] == pytest.approx(GM0)
    assert all(item.margin == pytest.approx(item.attained - item.required) for item in report.criteria)


def test_general_criteria_report_round_trips_through_json() -> None:
    report = general_criteria(analytic_curve(1.0), gm0_m=GM0)
    payload = json.loads(json.dumps(report.to_dict()))
    assert payload["satisfied"] is True
    assert payload["schema"] == "tensorcad.intact-stability/1"
    assert len(payload["criteria"]) == 6
    assert payload["failures"] == []
    assert payload["criteria"][0]["clause"] == "2.2.1"
    assert payload["criteria"][0]["required"] == IS_CODE_2008_GENERAL["2.2.1"][0]


def test_identical_calls_are_identical() -> None:
    angles, values = analytic_curve(0.5)
    first = general_criteria((angles, values), gm0_m=GM0).to_dict()
    second = general_criteria((angles, values), gm0_m=GM0).to_dict()
    assert first == second


# ---------------------------------------------------------------------------
# negative controls: each clause has to be able to fail on its own
# ---------------------------------------------------------------------------
def test_scaling_the_arm_down_flips_areas_and_the_lever_but_not_the_shape_clauses() -> None:
    angles, values = analytic_curve(0.25)
    halved = tuple(0.5 * value for value in values)
    report = general_criteria((angles, halved), gm0_m=0.5 * GM0)
    failed = {item.clause for item in report.failures}
    assert failed == {"2.2.1", "2.2.2", "2.2.4"}  # 2.2.3/2.2.5/2.2.6 survive halving
    assert report.verdict() == "fails IS Code 2008 part A 2.2.1, 2.2.2, 2.2.4"


def test_a_low_gm_fails_only_the_metacentric_height_clause() -> None:
    report = general_criteria(analytic_curve(0.25), gm0_m=0.14)
    assert {item.clause for item in report.failures} == {"2.2.6"}
    criterion = next(item for item in report.criteria if item.clause == "2.2.6")
    assert criterion.attained == pytest.approx(0.14)
    assert criterion.margin == pytest.approx(0.14 - 0.15)


def test_an_early_peak_fails_the_angle_of_maximum_lever() -> None:
    angles = ladder(0.5, 45.0)
    peaked = tuple(0.25 * math.sin(math.radians(4.5 * angle)) for angle in angles)
    report = general_criteria((angles, peaked), gm0_m=0.5)
    criterion = next(item for item in report.criteria if item.clause == "2.2.5")
    assert criterion.attained == pytest.approx(20.0)
    assert criterion.satisfied is False
    assert criterion in report.failures


def test_the_flooding_angle_replaces_forty_degrees_when_it_is_lower() -> None:
    angles, values = analytic_curve(0.25)
    report = general_criteria((angles, values), gm0_m=GM0, flooding_angle_deg=35.0)
    assert report.upper_angle_deg == 35.0
    assert report.area_0_40_m_rad == pytest.approx(wall_sided_area(0.0, 35.0), rel=1e-5)
    assert report.area_30_40_m_rad == pytest.approx(wall_sided_area(30.0, 35.0), rel=1e-5)


def test_a_flooding_angle_above_forty_leaves_forty_in_charge() -> None:
    report = general_criteria(analytic_curve(1.0), gm0_m=GM0, flooding_angle_deg=60.0)
    assert report.upper_angle_deg == 40.0
    assert report.flooding_angle_deg == 60.0


# ---------------------------------------------------------------------------
# malformed input is refused, not silently averaged
# ---------------------------------------------------------------------------
def test_a_curve_that_stops_short_of_the_clauses_is_refused() -> None:
    with pytest.raises(IntactStabilityError, match="2.2.2 needs 40"):
        general_criteria(analytic_curve(1.0, 39.0), gm0_m=GM0)


def test_a_short_curve_is_accepted_when_the_flooding_angle_covers_it() -> None:
    report = general_criteria(analytic_curve(0.5, 35.0), gm0_m=GM0, flooding_angle_deg=35.0)
    assert report.upper_angle_deg == 35.0


@pytest.mark.parametrize(
    ("curve", "message"),
    [
        (((0.0, 10.0, 8.0), (0.0, 0.1, 0.2)), "increase strictly"),
        (((0.0, 10.0), (0.0,)), "differ in length"),
        (((0.0,), (0.0,)), "at least two heel angles"),
        ((object(), None), "an \\(angles, gz\\) pair must hold two sequences"),
    ],
)
def test_malformed_curves_are_named_in_the_error(curve: object, message: str) -> None:
    with pytest.raises(IntactStabilityError, match=message):
        general_criteria(curve, gm0_m=GM0)


def test_a_flooding_angle_below_thirty_degrees_is_refused() -> None:
    with pytest.raises(IntactStabilityError, match="cannot support 2.2.2/2.2.3"):
        general_criteria(analytic_curve(1.0), gm0_m=GM0, flooding_angle_deg=25.0)


def test_non_finite_inputs_are_refused() -> None:
    with pytest.raises(IntactStabilityError, match="must be finite"):
        general_criteria(analytic_curve(1.0), gm0_m=float("nan"))


# ---------------------------------------------------------------------------
# end to end: mesh -> curve -> criteria, through both waterline solvers
# ---------------------------------------------------------------------------
def test_the_barge_mesh_lands_on_the_analytic_report_through_the_newton_route(barge: Mesh) -> None:
    curve = gz_curve(barge, KG, ladder(1.0), waterline_z=DRAFT)
    worst = max(
        abs(measured - wall_sided_gz(angle))
        for angle, measured in zip(curve.angles_deg, curve.gz)
    )
    assert worst < 1e-12  # the polyhedral integral is exact; only round-off is left
    report = general_criteria(curve, gm0_m=GM0)
    assert report.satisfied
    assert report.area_0_30_m_rad == pytest.approx(wall_sided_area(0.0, 30.0), rel=5e-4)
    assert report.area_0_40_m_rad == pytest.approx(wall_sided_area(0.0, 40.0), rel=5e-4)


def test_both_waterline_solvers_produce_the_same_criteria(barge: Mesh) -> None:
    angles = ladder(1.0)
    newton = general_criteria(gz_curve(barge, KG, angles, waterline_z=DRAFT), gm0_m=GM0)
    bisect = general_criteria(
        bisection_curve(barge, kg=KG, heels_deg=angles, waterline_z=DRAFT), gm0_m=GM0
    )
    assert newton.area_0_30_m_rad == pytest.approx(bisect.area_0_30_m_rad, abs=1e-9)
    assert newton.area_0_40_m_rad == pytest.approx(bisect.area_0_40_m_rad, abs=1e-9)
    assert newton.gz_30_m == pytest.approx(bisect.gz_30_m, abs=1e-9)
    assert {item.clause for item in newton.criteria if not item.satisfied} == set()
    assert {item.clause for item in bisect.criteria if not item.satisfied} == set()


def test_gm_from_the_waterplane_moments_closes_the_check(barge: Mesh) -> None:
    """``GM_0`` from ``tensorcad.stability`` and the curve's slope are the same number."""
    report = general_criteria(gz_curve(barge, KG, ladder(1.0), waterline_z=DRAFT), gm0_m=GM0)
    assert report.gm0_m == pytest.approx(KB0 + BM0 - KG)
    criterion = next(item for item in report.criteria if item.clause == "2.2.6")
    assert criterion.satisfied and criterion.margin == pytest.approx(GM0 - 0.15)


# ---------------------------------------------------------------------------
# weather criterion
# ---------------------------------------------------------------------------
def test_wind_levers_reproduce_the_code_formula() -> None:
    levers = wind_heeling_levers(12.0, 2.5, 100.0)
    assert levers.steady_m == pytest.approx(504.0 * 12.0 * 2.5 / (1000.0 * 9.81 * 100.0))
    assert levers.gust_m == pytest.approx(1.5 * levers.steady_m)
    assert levers.gust_factor == pytest.approx(1.5)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ((10.0, 2.0, 0.0), "displacement must be positive"),
        ((10.0, -1.0, 80.0), "must be non-negative"),
    ],
)
def test_wind_levers_refuse_impossible_loadings(kwargs: tuple, message: str) -> None:
    with pytest.raises(IntactStabilityError, match=message):
        wind_heeling_levers(*kwargs)


def _solve_wall_sided(lever: float, *, upper_deg: float = 45.0) -> float:
    """Bisection of ``GZ(phi) = lever`` on the monotone rising branch (test-side reference)."""
    low, high = 0.0, upper_deg
    for _ in range(200):
        middle = 0.5 * (low + high)
        if wall_sided_gz(middle) < lever:
            low = middle
        else:
            high = middle
    return 0.5 * (low + high)


def test_steady_heel_angle_matches_an_independent_bisection() -> None:
    angles, values = analytic_curve(0.05)
    lever = 0.05
    measured = steady_wind_heel_angle((angles, values), lever)
    assert measured is not None
    assert measured == pytest.approx(_solve_wall_sided(lever), abs=1e-3)
    assert wall_sided_gz(measured) == pytest.approx(lever, abs=2e-5)


def test_weather_criterion_passes_a_barge_under_moderate_wind() -> None:
    angles, values = analytic_curve(0.25, 55.0)
    report = weather_criterion(
        (angles, values),
        wind_area_m2=10.0,
        wind_arm_m=2.0,
        displacement_t=DISPLACEMENT_T,
    )
    assert report.satisfied
    assert report.levers.steady_m == pytest.approx(504.0 * 10.0 * 2.0 / (1000.0 * 9.81 * DISPLACEMENT_T))
    assert report.theta0_deg == pytest.approx(_solve_wall_sided(report.levers.steady_m), abs=1e-3)
    assert report.theta0_deg < 16.0
    assert report.theta2_deg == 50.0  # the wall-sided arm never re-crosses the gust lever
    assert report.area_a_m_rad == pytest.approx(report.levers.gust_m * math.radians(50.0))
    assert report.area_b_m_rad == pytest.approx(wall_sided_area(0.0, 50.0), rel=1e-4)
    assert report.area_b_m_rad > report.area_a_m_rad
    heel = next(item for item in report.criteria if item.clause == "2.3.1")
    assert heel.comparison == "not more than" and heel.satisfied is True
    assert heel.margin == pytest.approx(16.0 - report.theta0_deg)


def test_weather_criterion_fails_both_clauses_in_a_storm() -> None:
    angles, values = analytic_curve(0.25, 55.0)
    report = weather_criterion(
        (angles, values),
        wind_area_m2=240.0,
        wind_arm_m=2.0,
        displacement_t=DISPLACEMENT_T,
    )
    assert report.satisfied is False
    assert report.theta0_deg is not None and report.theta0_deg > 16.0
    assert {item.clause for item in report.failures} == {"2.3.1", "2.3.2"}
    assert report.area_a_m_rad > report.area_b_m_rad
    assert "fails" in report.verdict()


def test_the_eighty_percent_deck_edge_rule_can_bind_before_sixteen_degrees() -> None:
    angles, values = analytic_curve(0.25, 55.0)
    lever = 0.1  # GZ(12 deg) on this barge, so theta0 lands between 9.6 and 16 deg
    area = lever * (1000.0 * 9.81 * DISPLACEMENT_T) / (504.0 * 2.0)
    report = weather_criterion(
        (angles, values),
        wind_area_m2=area,
        wind_arm_m=2.0,
        displacement_t=DISPLACEMENT_T,
        deck_edge_angle_deg=12.0,
    )
    assert report.steady_heel_limit_deg == pytest.approx(9.6)  # 0.8 * 12, not 16
    assert report.theta0_deg is not None
    assert report.theta0_deg == pytest.approx(_solve_wall_sided(lever), abs=1e-3)
    assert 9.6 < report.theta0_deg < 16.0
    assert {item.clause for item in report.failures} == {"2.3.1"}
    criterion = next(item for item in report.criteria if item.clause == "2.3.1")
    assert criterion.comparison == "not more than"
    assert criterion.margin == pytest.approx(9.6 - report.theta0_deg)
    assert criterion.satisfied is False


def test_the_gust_areas_are_measured_from_the_windward_roll_amplitude() -> None:
    angles, values = analytic_curve(0.25, 55.0)
    report = weather_criterion(
        (angles, values),
        wind_area_m2=10.0,
        wind_arm_m=2.0,
        displacement_t=DISPLACEMENT_T,
        roll_amplitude_deg=8.0,
    )
    assert report.roll_amplitude_deg == 8.0
    assert report.area_b_m_rad == pytest.approx(wall_sided_area(8.0, 50.0), rel=1e-4)
    assert report.area_a_m_rad == pytest.approx(report.levers.gust_m * math.radians(50.0 - 8.0))


def test_the_second_intercept_caps_theta_two_when_the_arm_comes_back() -> None:
    """A curve that falls back through the gust lever ends the gust area at that crossing."""
    angles = ladder(0.5, 60.0)
    values = tuple(0.5 * math.sin(math.radians(4.5 * angle)) for angle in angles)  # peak at 20 deg, zero at 40
    report = weather_criterion(
        (angles, values),
        wind_area_m2=5.0,
        wind_arm_m=2.0,
        displacement_t=DISPLACEMENT_T,
    )
    assert report.theta_c_deg is not None
    assert report.theta2_deg == pytest.approx(report.theta_c_deg)
    assert report.theta2_deg < 50.0
    assert isinstance(report.criteria[1], Criterion)


def test_a_curve_too_short_to_locate_theta_two_is_refused() -> None:
    angles, values = analytic_curve(0.5, 40.0)
    with pytest.raises(IntactStabilityError, match="sample further to locate theta2"):
        weather_criterion(
            (angles, values),
            wind_area_m2=5.0,
            wind_arm_m=2.0,
            displacement_t=DISPLACEMENT_T,
        )


def test_weather_criterion_refuses_an_impossible_roll_amplitude() -> None:
    angles, values = analytic_curve(0.25, 55.0)
    with pytest.raises(IntactStabilityError, match="roll_amplitude_deg must be in"):
        weather_criterion(
            (angles, values),
            wind_area_m2=10.0,
            wind_arm_m=2.0,
            displacement_t=DISPLACEMENT_T,
            roll_amplitude_deg=90.0,
        )
    with pytest.raises(IntactStabilityError, match="does not exceed the roll amplitude"):
        weather_criterion(
            (angles, values),
            wind_area_m2=10.0,
            wind_arm_m=2.0,
            displacement_t=DISPLACEMENT_T,
            roll_amplitude_deg=50.0,
            downflooding_angle_deg=50.0,
        )


# ---------------------------------------------------------------------------
# clause 2.3.3: the windward roll amplitude shall not exceed the flooding angle
# ---------------------------------------------------------------------------
def test_clause_2_3_3_compares_the_roll_amplitude_with_the_flooding_angle() -> None:
    """2.3.3 is a "shall not exceed" clause like 2.3.1, judged on the supplied theta_1."""
    angles, values = analytic_curve(0.25, 55.0)
    report = weather_criterion(
        (angles, values),
        wind_area_m2=10.0,
        wind_arm_m=2.0,
        displacement_t=DISPLACEMENT_T,
        roll_amplitude_deg=8.0,
        flooding_angle_deg=20.0,
    )
    assert report.flooding_angle_deg == 20.0
    assert report.clauses == "2.3.1, 2.3.2 and 2.3.3"
    assert [item.clause for item in report.criteria] == ["2.3.1", "2.3.2", "2.3.3"]
    roll = report.criteria[2]
    assert roll.comparison == "not more than"
    assert (roll.attained, roll.required) == (8.0, 20.0)
    assert roll.satisfied is True
    assert roll.margin == pytest.approx(12.0)
    assert report.satisfied
    assert report.verdict() == "passes IS Code 2008 part A 2.3.1, 2.3.2 and 2.3.3 (weather criterion)"


def test_the_flooding_angle_caps_theta_two_and_moves_the_gust_areas() -> None:
    """theta_f enters 2.3.2's minimum too, so the gust areas stop at it, not at 50 deg."""
    angles, values = analytic_curve(0.25, 55.0)
    kwargs = {"wind_area_m2": 10.0, "wind_arm_m": 2.0, "displacement_t": DISPLACEMENT_T}
    uncapped = weather_criterion((angles, values), roll_amplitude_deg=8.0, **kwargs)
    capped = weather_criterion((angles, values), roll_amplitude_deg=8.0, flooding_angle_deg=20.0, **kwargs)
    assert uncapped.theta2_deg == 50.0
    assert capped.theta2_deg == pytest.approx(20.0)
    assert capped.area_b_m_rad == pytest.approx(wall_sided_area(8.0, 20.0), rel=1e-5)  # 0.25 deg grid
    assert capped.area_a_m_rad == pytest.approx(capped.levers.gust_m * math.radians(20.0 - 8.0))
    assert capped.area_b_m_rad < uncapped.area_b_m_rad


def test_a_roll_amplitude_reaching_the_flooding_angle_fails_2_3_3() -> None:
    """theta_1 >= theta_f is a failure of 2.3.3, not an exception and not a fake area."""
    angles, values = analytic_curve(0.25, 55.0)
    report = weather_criterion(
        (angles, values),
        wind_area_m2=10.0,
        wind_arm_m=2.0,
        displacement_t=DISPLACEMENT_T,
        roll_amplitude_deg=25.0,
        flooding_angle_deg=20.0,
    )
    assert report.satisfied is False
    assert [item.clause for item in report.failures] == ["2.3.3"]
    assert report.criteria[-1].margin == pytest.approx(-5.0)
    # the gust areas have no span between theta_1 and theta_2: they are unmeasured, not zero
    assert report.area_b_m_rad is None and report.area_a_m_rad is None
    assert report.theta2_deg == pytest.approx(20.0)
    assert report.clauses == "2.3.1 and 2.3.3"
    assert report.verdict() == "fails IS Code 2008 part A 2.3.3"
    payload = json.loads(json.dumps(report.to_dict()))
    assert payload["area_b_m_rad"] is None and payload["area_a_m_rad"] is None
    assert payload["flooding_angle_deg"] == 20.0
    assert payload["clauses"] == "2.3.1 and 2.3.3"


def test_clause_2_3_3_is_not_evaluated_without_a_flooding_angle() -> None:
    """No theta_f supplied: the clause is absent and said to be absent, not passed."""
    angles, values = analytic_curve(0.25, 55.0)
    report = weather_criterion(
        (angles, values),
        wind_area_m2=10.0,
        wind_arm_m=2.0,
        displacement_t=DISPLACEMENT_T,
    )
    assert report.flooding_angle_deg is None
    assert "2.3.3" not in [item.clause for item in report.criteria]
    assert report.clauses == "2.3.1 and 2.3.2"
    assert report.verdict() == "passes IS Code 2008 part A 2.3.1 and 2.3.2 (weather criterion)"
    payload = report.to_dict()
    assert payload["flooding_angle_deg"] is None
    assert "2.3.3" not in json.dumps(payload)


@pytest.mark.parametrize("flooding", [0.0, -5.0, 90.0, 130.0])
def test_a_flooding_angle_outside_zero_to_ninety_is_refused(flooding: float) -> None:
    angles, values = analytic_curve(1.0, 55.0)
    with pytest.raises(IntactStabilityError, match="flooding_angle_deg must be in"):
        weather_criterion(
            (angles, values),
            wind_area_m2=10.0,
            wind_arm_m=2.0,
            displacement_t=DISPLACEMENT_T,
            flooding_angle_deg=flooding,
        )
