"""The criteria are checked against curves whose areas are known in closed form.

Three levels, because "the number looks plausible" is not a check:

1. **Pure table.** The synthetic curve ``GZ(phi) = A sin(2 phi)`` has the exact area
   ``A/2 (cos 2a - cos 2b)`` on any interval, so the trapezoid rule and the criterion
   wiring are compared against arithmetic with no discretisation in it at all.
2. **The closed form of the wall-sided hull.** ``int_0^theta GZ dphi =
   GM (1 - cos theta) + BM/2 (sec theta + cos theta - 2)`` is integrated numerically as
   an independent check of :func:`tensorcad.intact_stability.wall_sided_area`, and then a
   *real box-barge mesh* is heeled, its ``GZ`` curve built by ``tensorcad.gz``, and the
   areas compared against that form. The barge is 1.6 m wide and 2.2 m deep, deeper than
   wide, so that the wall-sided validity limit (bilge emergence at ``tan phi = 2T/B =
   1``, i.e. 45 deg) sits *above* the 40 deg the Code quotes and the comparison is not
   stopping at the code's boundary.
3. **Negative controls.** A ``KG`` shift must move the areas by exactly ``-dKG (1-cos
   theta)``; a curve sampled only to 20 deg must be reported as *unable to answer*, not as
   a pass; and a threshold placed one nanometre above the achieved value must flip the
   verdict. Without these, agreement could be vacuous.
"""

from __future__ import annotations

import json
import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.gz import righting_arm_curve, wall_sided_righting_arm
from tensorcad.intact_stability import (
    ANGLE_30_DEG,
    ANGLE_40_DEG,
    NOT_SATISFIED,
    SATISFIED,
    UNDETERMINED,
    IntactStabilityError,
    IntactThresholds,
    area_between,
    assess_curve,
    assess_intact_stability,
    initial_gm_from_curve,
    interpolate_gz,
    wall_sided_area,
)
from tensorcad.stability import initial_stability

# Box barge: deeper than wide so the wall-sided limit (45 deg) clears the 40 deg clause.
LENGTH = 10.0
BREADTH = 1.6
DRAFT = 0.8
DEPTH = 2.2
KG = 0.5
KB0 = DRAFT / 2.0
BM0 = BREADTH**2 / (12.0 * DRAFT)
GM0 = KB0 + BM0 - KG  # 0.1666... m


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


@pytest.fixture(scope="module")
def barge() -> Mesh:
    """Freeboarded barge: keel at z=0, deck at DEPTH, floating at z = DRAFT."""
    return box((0.0, -BREADTH / 2, 0.0), (LENGTH, BREADTH / 2, DEPTH))


@pytest.fixture(scope="module")
def barge_submerged() -> Mesh:
    """The same barge clipped at the waterline: the body ``compute_hydrostatics`` expects."""
    return box((0.0, -BREADTH / 2, 0.0), (LENGTH, BREADTH / 2, DRAFT))


def sine_table(amplitude: float = 0.5, step_deg: float = 0.5, limit_deg: float = 90.0):
    """``GZ(phi) = A sin(2 phi)`` on a regular grid: an arm with an exact area."""
    count = round(limit_deg / step_deg)
    angles = tuple(index * step_deg for index in range(count + 1))
    arms = tuple(amplitude * math.sin(2.0 * math.radians(angle)) for angle in angles)
    return angles, arms


def sine_area(amplitude: float, start_deg: float, end_deg: float) -> float:
    """Exact ``int GZ dphi`` for ``A sin(2 phi)``, in m*rad."""
    return 0.5 * amplitude * (math.cos(2.0 * math.radians(start_deg)) - math.cos(2.0 * math.radians(end_deg)))


# --------------------------------------------------------------------------------------
# 1. The table integrator, against a curve with an exact area
# --------------------------------------------------------------------------------------


def test_area_matches_the_exact_integral_of_sin_two_phi() -> None:
    """Trapezoid on a 0.1 deg table: exact to the discretisation, which is ~5e-7 m*rad."""
    angles, arms = sine_table(0.5, step_deg=0.1)
    for start, end in ((0.0, 30.0), (0.0, 40.0), (30.0, 40.0), (0.0, 90.0), (12.5, 37.5)):
        measured = area_between(angles, arms, start, end)
        exact = sine_area(0.5, start, end)
        assert measured == pytest.approx(exact, abs=1e-6), (start, end)


def test_the_table_rule_is_second_order_on_a_curve_with_an_exact_area() -> None:
    exact = sine_area(0.5, 0.0, 40.0)
    errors = []
    for step in (1.0, 0.5, 0.25):
        angles, arms = sine_table(0.5, step_deg=step)
        errors.append(abs(area_between(angles, arms, 0.0, 40.0) - exact))
    ratios = [errors[index] / errors[index + 1] for index in range(len(errors) - 1)]
    assert all(3.2 < ratio < 4.8 for ratio in ratios), (errors, ratios)


def test_area_is_additive_across_an_interior_point() -> None:
    angles, arms = sine_table(0.4, step_deg=0.5)
    whole = area_between(angles, arms, 0.0, 40.0)
    parts = area_between(angles, arms, 0.0, 30.0) + area_between(angles, arms, 30.0, 40.0)
    assert whole == pytest.approx(parts, rel=0.0, abs=1e-12)


def test_area_of_an_empty_interval_is_zero() -> None:
    angles, arms = sine_table()
    assert area_between(angles, arms, 30.0, 30.0) == 0.0


def test_area_is_signed_when_the_arm_goes_negative() -> None:
    """An arm that is negative over an interval must subtract, not be taken absolute.

    The table below is hand-integrable: in deg*m the trapezoids are +0.2*30, -0.2*30,
    -0.4*30, -0.2*30, +0.2*60, i.e. -6 deg*m over the whole 0-180 deg range.
    """
    angles = (0.0, 30.0, 60.0, 90.0, 120.0, 180.0)
    arms = (0.4, 0.0, -0.4, -0.4, 0.0, 0.4)
    assert area_between(angles, arms, 0.0, 180.0) == pytest.approx(math.radians(-6.0), abs=1e-12)
    assert area_between(angles, arms, 60.0, 120.0) == pytest.approx(math.radians(-18.0), abs=1e-12)


def test_interpolation_is_exact_at_the_samples_and_linear_between() -> None:
    angles, arms = sine_table(0.5, step_deg=1.0, limit_deg=10.0)
    assert interpolate_gz(angles, arms, 3.0) == pytest.approx(arms[3], abs=0.0)
    midpoint = interpolate_gz(angles, arms, 3.5)
    assert midpoint == pytest.approx(0.5 * (arms[3] + arms[4]), abs=1e-15)


# --------------------------------------------------------------------------------------
# 2. The wall-sided closed form, and a real barge against it
# --------------------------------------------------------------------------------------


def test_wall_sided_area_agrees_with_simpson_on_its_own_arm() -> None:
    gm, bm = 0.35, 0.30
    for limit in (10.0, 25.0, 30.0, 40.0):
        steps = 20000
        step = math.radians(limit) / steps
        total = 0.0
        for index in range(steps + 1):
            phi = index * step
            weight = 1.0 if index in (0, steps) else (4.0 if index % 2 else 2.0)
            total += weight * wall_sided_righting_arm(math.degrees(phi), gm=gm, bm=bm)
        simpson = total * step / 3.0
        assert wall_sided_area(limit, gm=gm, bm=bm) == pytest.approx(simpson, rel=1e-10, abs=0.0)


def test_wall_sided_area_difference_between_two_limits() -> None:
    gm, bm = 0.3, 0.27
    between = wall_sided_area(40.0, gm=gm, bm=bm) - wall_sided_area(30.0, gm=gm, bm=bm)
    assert wall_sided_area(40.0, gm=gm, bm=bm, start_deg=30.0) == pytest.approx(between, abs=1e-15)


def test_barge_areas_reproduce_the_wall_sided_closed_form(barge: Mesh) -> None:
    angles = tuple(0.0 + 2.5 * index for index in range(17))  # 0 .. 40 deg
    curve = righting_arm_curve(barge, angles_deg=angles, kg=KG, waterline_z=DRAFT)
    area_30 = area_between(curve.angles_deg, curve.gz, 0.0, 30.0)
    area_40 = area_between(curve.angles_deg, curve.gz, 0.0, 40.0)
    assert area_30 == pytest.approx(wall_sided_area(30.0, gm=GM0, bm=BM0), abs=2e-4)
    assert area_40 == pytest.approx(wall_sided_area(40.0, gm=GM0, bm=BM0), abs=4e-4)
    # And the arms themselves are the closed form, so an area match is not a coincidence
    # of two errors: check a mid-curve sample too.
    assert interpolate_gz(curve.angles_deg, curve.gz, 20.0) == pytest.approx(
        wall_sided_righting_arm(20.0, gm=GM0, bm=BM0), abs=1e-9
    )


def test_barge_area_error_is_second_order_in_the_angle_step(barge: Mesh) -> None:
    exact = wall_sided_area(30.0, gm=GM0, bm=BM0)
    errors = []
    for step in (5.0, 2.5, 1.25):
        angles = tuple(index * step for index in range(round(30.0 / step) + 1))
        curve = righting_arm_curve(barge, angles_deg=angles, kg=KG, waterline_z=DRAFT)
        errors.append(abs(area_between(curve.angles_deg, curve.gz, 0.0, 30.0) - exact))
    ratios = [errors[index] / errors[index + 1] for index in range(len(errors) - 1)]
    assert all(3.2 < ratio < 4.8 for ratio in ratios), (errors, ratios)


def test_kg_shift_moves_the_area_by_the_closed_form_amount(barge: Mesh) -> None:
    """Negative control: the area cannot agree by accident if it tracks ``KG`` exactly."""
    angles = tuple(0.0 + 0.5 * index for index in range(61))  # 0 .. 30 deg
    base = righting_arm_curve(barge, angles_deg=angles, kg=KG, waterline_z=DRAFT)
    shifted = righting_arm_curve(barge, angles_deg=angles, kg=KG + 0.1, waterline_z=DRAFT)
    delta = area_between(shifted.angles_deg, shifted.gz, 0.0, 30.0) - area_between(
        base.angles_deg, base.gz, 0.0, 30.0
    )
    # dGZ = -0.1 sin(phi) => d(area) = -0.1 (1 - cos 30 deg), independent of BM.
    expected = -0.1 * (1.0 - math.cos(math.radians(30.0)))
    assert delta == pytest.approx(expected, rel=0.0, abs=2e-7)


# --------------------------------------------------------------------------------------
# 3. Criterion wiring: verdicts, truncation, and honest "cannot say"
# --------------------------------------------------------------------------------------


def test_synthetic_curve_satisfies_all_six_clauses_against_closed_form_values() -> None:
    angles, arms = sine_table(0.5, step_deg=0.05)
    report = assess_intact_stability(angles, arms, initial_gm=1.0)  # slope of A sin 2phi
    assert [result.criterion for result in report.criteria] == [
        "IS-2008-PartA-2.2.1",
        "IS-2008-PartA-2.2.2",
        "IS-2008-PartA-2.2.3",
        "IS-2008-PartA-2.2.4",
        "IS-2008-PartA-2.2.5",
        "IS-2008-PartA-2.2.6",
    ]
    assert report.all_satisfied is True
    assert report.failed == () and report.undetermined == ()
    assert report.criterion("IS-2008-PartA-2.2.1").achieved == pytest.approx(
        sine_area(0.5, 0.0, 30.0), abs=1e-7
    )
    assert report.criterion("IS-2008-PartA-2.2.2").achieved == pytest.approx(
        sine_area(0.5, 0.0, 40.0), abs=1e-7
    )
    assert report.criterion("IS-2008-PartA-2.2.3").achieved == pytest.approx(
        sine_area(0.5, 30.0, 40.0), abs=1e-7
    )
    assert report.criterion("IS-2008-PartA-2.2.4").achieved == pytest.approx(
        0.5 * math.sin(math.radians(60.0))
    )
    # A sin(2 phi) peaks at 45 deg, which clears the 25 deg clause on the sampled grid.
    assert report.criterion("IS-2008-PartA-2.2.5").achieved == pytest.approx(45.0)
    assert report.criterion("IS-2008-PartA-2.2.6").achieved == pytest.approx(1.0)


def test_barge_curve_fails_the_area_and_arm_clauses_honestly(barge: Mesh, barge_submerged: Mesh) -> None:
    """A slender box barge has GM 0.167 m: it passes the initial-GM clause and fails 2.2.1."""
    angles = tuple(0.0 + 5.0 * index for index in range(19))  # 0 .. 90 deg
    curve = righting_arm_curve(barge, angles_deg=angles, kg=KG, waterline_z=DRAFT)
    gm = initial_stability(barge_submerged, kg=KG).gm_transverse
    assert gm == pytest.approx(GM0, abs=1e-9)
    assert initial_stability(barge_submerged, kg=KG).bm_transverse == pytest.approx(BM0, abs=1e-12)
    report = assess_curve(curve, initial_gm=gm)
    assert report.criterion("IS-2008-PartA-2.2.1").status == NOT_SATISFIED
    assert report.criterion("IS-2008-PartA-2.2.4").satisfied is False
    assert report.criterion("IS-2008-PartA-2.2.6").satisfied is True
    assert report.all_satisfied is False
    assert "IS-2008-PartA-2.2.1" in report.failed
    # The 0-30 deg area is the closed form up to the sampling error, and is below 0.055.
    assert report.criterion("IS-2008-PartA-2.2.1").achieved == pytest.approx(
        wall_sided_area(30.0, gm=GM0, bm=BM0), abs=2e-4
    )
    assert report.criterion("IS-2008-PartA-2.2.1").margin < 0.0


def test_a_short_curve_answers_nothing_rather_than_passing(barge: Mesh) -> None:
    angles = tuple(0.0 + 5.0 * index for index in range(5))  # 0 .. 20 deg
    curve = righting_arm_curve(barge, angles_deg=angles, kg=KG, waterline_z=DRAFT)
    report = assess_curve(curve, initial_gm=GM0)
    for name in (
        "IS-2008-PartA-2.2.1",
        "IS-2008-PartA-2.2.2",
        "IS-2008-PartA-2.2.3",
        "IS-2008-PartA-2.2.4",
    ):
        result = report.criterion(name)
        assert result.status == UNDETERMINED, name
        assert result.achieved is None
        assert "widen the curve" in result.note
    assert report.all_satisfied is False
    assert set(report.undetermined) >= {"IS-2008-PartA-2.2.1", "IS-2008-PartA-2.2.2"}
    assert report.criterion("IS-2008-PartA-2.2.6").status == SATISFIED  # GM is known


def test_the_flooding_angle_shortens_the_40_deg_intervals() -> None:
    angles, arms = sine_table(0.5, step_deg=0.05)
    report = assess_intact_stability(angles, arms, initial_gm=1.0, flooding_angle_deg=35.0)
    assert report.flooding_angle_deg == 35.0
    assert report.criterion("IS-2008-PartA-2.2.2").achieved == pytest.approx(
        sine_area(0.5, 0.0, 35.0), abs=1e-7
    )
    assert report.criterion("IS-2008-PartA-2.2.3").achieved == pytest.approx(
        sine_area(0.5, 30.0, 35.0), abs=1e-7
    )
    assert "truncated at the flooding angle 35 deg" in report.criterion("IS-2008-PartA-2.2.2").note


def test_a_flooding_angle_below_30_deg_empties_the_middle_clause() -> None:
    angles, arms = sine_table(0.5, step_deg=0.05)
    report = assess_intact_stability(angles, arms, initial_gm=1.0, flooding_angle_deg=25.0)
    middle = report.criterion("IS-2008-PartA-2.2.3")
    assert middle.status == UNDETERMINED
    assert "empty" in middle.note
    first = report.criterion("IS-2008-PartA-2.2.1")
    assert first.achieved == pytest.approx(sine_area(0.5, 0.0, 25.0), abs=1e-7)
    assert "below 30 deg" in first.note


def test_a_threshold_a_nanometre_above_the_achieved_value_flips_the_verdict() -> None:
    angles, arms = sine_table(0.5, step_deg=0.5)
    achieved = area_between(angles, arms, 0.0, 30.0)
    at = assess_intact_stability(
        angles,
        arms,
        initial_gm=1.0,
        thresholds=IntactThresholds(area_to_30_m_rad=achieved),
        tolerance=0.0,
    )
    assert at.criterion("IS-2008-PartA-2.2.1").status == SATISFIED
    above = assess_intact_stability(
        angles,
        arms,
        initial_gm=1.0,
        thresholds=IntactThresholds(area_to_30_m_rad=achieved + 1e-9),
        tolerance=0.0,
    )
    assert above.criterion("IS-2008-PartA-2.2.1").status == NOT_SATISFIED


def test_a_peak_on_the_last_sample_is_flagged_not_trusted() -> None:
    angles, arms = sine_table(0.5, step_deg=1.0, limit_deg=30.0)  # peak (45 deg) is beyond
    report = assess_intact_stability(angles, arms, initial_gm=1.0)
    peak = report.criterion("IS-2008-PartA-2.2.5")
    assert peak.achieved == pytest.approx(30.0)
    assert "beyond" in peak.note


def test_initial_gm_falls_back_to_the_curve_chord_and_reads_high(barge: Mesh) -> None:
    angles = tuple(0.0 + 1.0 * index for index in range(41))
    curve = righting_arm_curve(barge, angles_deg=angles, kg=KG, waterline_z=DRAFT)
    gm, at_angle = initial_gm_from_curve(curve.angles_deg, curve.gz)
    assert at_angle == pytest.approx(1.0)
    assert gm > GM0  # the chord carries the BM tan^2 curvature and therefore reads high
    assert gm - GM0 == pytest.approx(0.5 * BM0 * math.tan(math.radians(1.0)) ** 2, abs=1e-12)
    report = assess_curve(curve)  # no GM supplied at all
    assert report.initial_gm_source.startswith("curve chord at 1 deg")
    assert report.criterion("IS-2008-PartA-2.2.6").status == SATISFIED


def test_gm_is_reported_as_unavailable_when_the_curve_starts_far_from_zero() -> None:
    angles = (20.0, 25.0, 30.0, 35.0, 40.0)
    arms = tuple(0.3 for _ in angles)
    report = assess_intact_stability(angles, arms)
    gm_clause = report.criterion("IS-2008-PartA-2.2.6")
    assert gm_clause.status == UNDETERMINED
    assert report.initial_gm_source == "unavailable"


def test_free_surface_correction_lowers_the_reported_gm(barge: Mesh) -> None:
    angles = tuple(0.0 + 1.0 * index for index in range(41))
    solid = righting_arm_curve(barge, angles_deg=angles, kg=KG, waterline_z=DRAFT)
    slack = righting_arm_curve(barge, angles_deg=angles, kg=KG, waterline_z=DRAFT, free_surface_i_x=0.05)
    gm_solid, _ = initial_gm_from_curve(solid.angles_deg, solid.gz)
    gm_slack, _ = initial_gm_from_curve(slack.angles_deg, slack.gz)
    assert gm_solid - gm_slack == pytest.approx(
        slack.free_surface_correction / math.sin(math.radians(1.0)), rel=1e-9
    )


def test_report_serialises_to_json_with_the_expected_schema() -> None:
    angles, arms = sine_table(0.5, step_deg=0.5)
    report = assess_intact_stability(angles, arms, initial_gm=1.0, flooding_angle_deg=45.0)
    payload = json.loads(json.dumps(report.to_dict()))
    assert payload["schema"] == "tensorcad.stability.intact_assessment/1"
    assert payload["code"].startswith("IMO IS Code 2008")
    assert len(payload["criteria"]) == 6
    assert all(item["schema"] == "tensorcad.stability.intact_criterion/1" for item in payload["criteria"])
    assert payload["thresholds"]["area_to_30_m_rad"] == 0.055
    assert payload["all_satisfied"] is True
    assert payload["flooding_angle_deg"] == 45.0
    assert payload["initial_gm_source"] == "supplied"


# --------------------------------------------------------------------------------------
# 4. Refusals: the module raises rather than inventing a number
# --------------------------------------------------------------------------------------


def test_mismatched_tables_are_refused() -> None:
    with pytest.raises(IntactStabilityError, match="same length"):
        area_between((0.0, 10.0), (0.1,), 0.0, 10.0)


def test_non_increasing_angles_are_refused() -> None:
    with pytest.raises(IntactStabilityError, match="strictly increase"):
        area_between((0.0, 10.0, 10.0), (0.1, 0.2, 0.3), 0.0, 10.0)


def test_an_interval_beyond_the_sampled_range_is_refused() -> None:
    angles, arms = sine_table(0.5, step_deg=1.0, limit_deg=20.0)
    with pytest.raises(IntactStabilityError, match="beyond the last sample"):
        area_between(angles, arms, 0.0, 30.0)
    with pytest.raises(IntactStabilityError, match="below the first sample"):
        area_between(angles, arms, -5.0, 10.0)
    with pytest.raises(IntactStabilityError, match="before start"):
        area_between(angles, arms, 10.0, 5.0)


def test_non_finite_arms_and_angles_are_refused() -> None:
    with pytest.raises(IntactStabilityError, match="not finite"):
        area_between((0.0, math.nan, 20.0), (0.1, 0.2, 0.3), 0.0, 20.0)
    with pytest.raises(IntactStabilityError, match="not finite"):
        area_between((0.0, 10.0, 20.0), (0.1, math.inf, 0.3), 0.0, 20.0)


def test_wall_sided_area_refuses_the_region_where_the_form_does_not_hold() -> None:
    with pytest.raises(IntactStabilityError, match=r"\|heel\| < 90"):
        wall_sided_area(90.0, gm=0.3, bm=0.3)
    with pytest.raises(IntactStabilityError, match="must not be below"):
        wall_sided_area(20.0, gm=0.3, bm=0.3, start_deg=25.0)


def test_negative_thresholds_and_bad_tolerance_are_refused() -> None:
    with pytest.raises(IntactStabilityError, match="must not be negative"):
        IntactThresholds(initial_gm_m=-0.15)
    with pytest.raises(IntactStabilityError, match="tolerance"):
        assess_intact_stability((0.0, 30.0), (0.0, 0.3), initial_gm=0.2, tolerance=-1e-9)
    with pytest.raises(IntactStabilityError, match="flooding_angle_deg"):
        assess_intact_stability((0.0, 30.0), (0.0, 0.3), initial_gm=0.2, flooding_angle_deg=0.0)


def test_a_single_sample_is_not_a_curve() -> None:
    with pytest.raises(IntactStabilityError, match="at least two samples"):
        area_between((0.0,), (0.0,), 0.0, 0.0)


def test_unknown_criterion_lookup_raises() -> None:
    angles, arms = sine_table(0.5, step_deg=0.5)
    report = assess_intact_stability(angles, arms, initial_gm=1.0)
    with pytest.raises(KeyError):
        report.criterion("IS-2008-PartA-9.9.9")


def test_angle_constants_are_the_code_angles() -> None:
    assert ANGLE_30_DEG == 30.0 and ANGLE_40_DEG == 40.0
