"""Slack tanks at heel: exact box benchmarks, closed forms and the fixed curve convention.

Three independent references are used, none of which is another route through the same
code:

1. **Hand-derived trapezoid.** While the free surface runs wall to wall the wetted section
   is a trapezoid with sides ``z1 = h - (b/2) tan(phi)`` and ``z2 = h + (b/2) tan(phi)``
   (``h`` the fill level), so

       y_bar = b^2 tan(phi) / (12 h),   z_bar = (h^2 + (b tan(phi))^2/12) / (2 h)

   which is written out in the test, not called from the module.
2. **The ``1/cos^3`` inertia law.** A chord across a rectangle tilted by ``phi`` is
   ``b / cos(phi)`` wide, so ``i_x(phi) = l (b/cos(phi))^3 / 12``. The module never uses
   this law; the test checks the module reproduces it, and checks it *stops* holding past
   the wall-sided window so the assertion cannot pass vacuously.
3. **The classical first-order correction.** ``dGZ/sin(phi) -> -i_x/nabla`` as ``phi -> 0``,
   with the difference falling as ``phi^2``.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.free_surface import (
    FreeSurfaceError,
    SlackTank,
    free_surface_gz_shift,
    free_surface_gz_shift_curve,
    free_surface_loss_m,
    tank_free_surface,
)
from tensorcad.gz import GzError, rectangular_tank_inertia, righting_arm_curve
from tensorcad.hull import HullParameters, generate_hull

LENGTH = 2.0
BREADTH = 0.5
HEIGHT = 1.0
FILL = 0.5 * LENGTH * BREADTH * HEIGHT  # half full: level 0.5 m
DISPLACEMENT = 7.446554  # the hull example's displacement, m^3


def tank(fill: float = FILL, *, breadth: float = BREADTH, length: float = LENGTH) -> SlackTank:
    return SlackTank(length=length, breadth=breadth, height=HEIGHT, fill_volume=fill)


def trapezoid_reference(t: SlackTank, heel_deg: float) -> tuple[float, float]:
    """Independent closed form for the liquid centroid in the tank frame."""
    level = t.fill_volume / (t.length * t.breadth)
    tangent = math.tan(math.radians(heel_deg))
    y_bar = t.breadth**2 * tangent / (12.0 * level)
    z_bar = (level**2 + (t.breadth * tangent) ** 2 / 12.0) / (2.0 * level)
    return y_bar, z_bar


def wall_sided_window_deg(t: SlackTank) -> float:
    """Heel at which the free surface stops running wall to wall."""
    level = t.fill_volume / (t.length * t.breadth)
    to_bottom = math.degrees(math.atan(2.0 * level / t.breadth))
    to_top = math.degrees(math.atan(2.0 * (t.height - level) / t.breadth))
    return min(to_bottom, to_top)


# --------------------------------------------------------------------------- tank model


def test_tank_rejects_non_physical_input() -> None:
    with pytest.raises(FreeSurfaceError):
        SlackTank(length=0.0, breadth=1.0, height=1.0, fill_volume=0.5)
    with pytest.raises(FreeSurfaceError):
        SlackTank(length=1.0, breadth=1.0, height=1.0, fill_volume=0.0)
    with pytest.raises(FreeSurfaceError):
        SlackTank(length=1.0, breadth=1.0, height=1.0, fill_volume=1.0)  # full: no free surface
    with pytest.raises(FreeSurfaceError):
        SlackTank(length=1.0, breadth=1.0, height=1.0, fill_volume=float("nan"))
    with pytest.raises(FreeSurfaceError):
        SlackTank(length=1.0, breadth=1.0, height=1.0, fill_volume=0.5, density_ratio=0.0)


def test_upright_state_is_the_textbook_rectangle() -> None:
    t = tank()
    state = tank_free_surface(t, 0.0)
    level = FILL / (LENGTH * BREADTH)
    assert state.volume == pytest.approx(FILL, rel=1e-14)
    assert state.free_surface_width == pytest.approx(BREADTH, rel=1e-14)
    assert state.free_surface_area == pytest.approx(LENGTH * BREADTH, rel=1e-14)
    assert state.inertia == pytest.approx(t.upright_inertia, rel=1e-14)
    assert state.inertia == pytest.approx(rectangular_tank_inertia(LENGTH, BREADTH), rel=1e-14)
    assert state.centroid_y == pytest.approx(0.0, abs=1e-18)
    assert state.centroid_z == pytest.approx(level / 2.0, rel=1e-14)
    assert state.level == pytest.approx(level, rel=1e-14)


def test_upright_inertia_does_not_depend_on_the_fill() -> None:
    """A rectangular tank's i_x = l b^3/12 is the same at 5 % and 90 % full."""
    reference = tank_free_surface(tank(0.05 * LENGTH * BREADTH * HEIGHT), 0.0).inertia
    for fraction in (0.2, 0.5, 0.9):
        fill = fraction * LENGTH * BREADTH * HEIGHT
        state = tank_free_surface(tank(fill), 0.0)
        assert state.free_surface_width == pytest.approx(BREADTH, rel=1e-14)
        assert state.inertia == reference
        assert state.free_surface_area == pytest.approx(LENGTH * BREADTH, rel=1e-14)


def test_liquid_volume_is_conserved_at_every_heel() -> None:
    t = tank()
    for heel in (0.0, 1.0, 10.0, 30.0, 45.0, 60.0, 89.0):
        assert tank_free_surface(t, heel).volume == pytest.approx(FILL, rel=1e-14)


def test_wall_sided_window_matches_the_cosine_law() -> None:
    """i_x(phi) = l (b/cos(phi))^3/12 while the free surface runs wall to wall."""
    t = tank()
    window = wall_sided_window_deg(t)
    assert window == pytest.approx(math.degrees(math.atan(2.0 * 0.5 / BREADTH)), rel=1e-12)
    for heel in (0.0, 5.0, 10.0, 20.0, 30.0, 40.0):
        assert heel < window
        state = tank_free_surface(t, heel)
        width = BREADTH / math.cos(math.radians(heel))
        assert state.free_surface_width == pytest.approx(width, rel=1e-13)
        assert state.inertia == pytest.approx(LENGTH * width**3 / 12.0, rel=1e-13)
        assert state.inertia == pytest.approx(t.upright_inertia / math.cos(math.radians(heel)) ** 3, rel=1e-13)


def test_the_cosine_law_fails_past_the_window() -> None:
    """The window test above must be able to fail: past it the chord is cut by a face."""
    t = tank()
    heel = wall_sided_window_deg(t) + 10.0
    state = tank_free_surface(t, heel)
    law = BREADTH / math.cos(math.radians(heel))
    assert abs(state.free_surface_width / law - 1.0) > 0.01


def test_centroid_matches_the_independent_trapezoid_formula() -> None:
    t = tank()
    upright = tank_free_surface(t, 0.0)
    for heel in (2.0, 5.0, 10.0, 20.0, 30.0, 40.0):
        state = tank_free_surface(t, heel)
        y_bar, z_bar = trapezoid_reference(t, heel)
        assert state.centroid_y - upright.centroid_y == pytest.approx(y_bar, rel=1e-13, abs=1e-15)
        assert state.centroid_z - upright.centroid_z == pytest.approx(
            z_bar - upright.centroid_z, rel=1e-13, abs=1e-15
        )


def test_shift_matches_the_independent_trapezoid_formula() -> None:
    """The primary output, against a hand-derived closed form (no shared code path)."""
    t = tank()
    for heel in (2.0, 10.0, 25.0, 40.0):
        y_bar, z_bar = trapezoid_reference(t, heel)
        walk = y_bar * math.cos(math.radians(heel)) + (z_bar - tank_free_surface(t, 0.0).centroid_z) * math.sin(
            math.radians(heel)
        )
        expected = -(t.density_ratio * FILL / DISPLACEMENT) * walk
        assert free_surface_gz_shift(t, heel, DISPLACEMENT) == pytest.approx(expected, rel=1e-13)


def test_small_angle_limit_is_the_classical_ix_over_nabla() -> None:
    t = tank()
    loss = free_surface_loss_m([t], DISPLACEMENT)
    assert loss == pytest.approx(t.upright_inertia / DISPLACEMENT, rel=1e-15)
    for heel in (0.001, 0.01, 0.1):
        slope = -free_surface_gz_shift(t, heel, DISPLACEMENT) / math.sin(math.radians(heel))
        assert abs(slope / loss - 1.0) < 2e-5
    # the defect of the first-order law grows with heel: it is zero at zero heel
    slopes = [
        -free_surface_gz_shift(t, heel, DISPLACEMENT) / math.sin(math.radians(heel))
        for heel in (0.1, 1.0, 10.0)
    ]
    assert slopes[0] == pytest.approx(loss, rel=1e-5)
    assert slopes[1] > slopes[0]
    assert slopes[2] > slopes[1]
    assert slopes[2] / loss - 1.0 > 1e-3


def test_shift_is_odd_in_the_heel() -> None:
    t = tank()
    for heel in (5.0, 20.0, 35.0):
        assert free_surface_gz_shift(t, -heel, DISPLACEMENT) == pytest.approx(
            -free_surface_gz_shift(t, heel, DISPLACEMENT), rel=1e-14
        )


def test_density_ratio_scales_the_shift_and_the_loss() -> None:
    sea_water = tank()
    oil = SlackTank(
        length=LENGTH,
        breadth=BREADTH,
        height=HEIGHT,
        fill_volume=FILL,
        density_ratio=0.85,
    )
    for heel in (0.0, 10.0, 30.0):
        assert free_surface_gz_shift(oil, heel, DISPLACEMENT) == pytest.approx(
            0.85 * free_surface_gz_shift(sea_water, heel, DISPLACEMENT), rel=1e-15
        )
    assert free_surface_loss_m([oil], DISPLACEMENT) == pytest.approx(
        0.85 * free_surface_loss_m([sea_water], DISPLACEMENT), rel=1e-15
    )


def test_equal_athwartship_subdivisions_reduce_the_inertia_by_n_squared() -> None:
    """The classical 1/n^2: splitting a tank across its breadth, not along it."""
    for compartments in (1, 2, 4, 5):
        breadth = BREADTH / compartments
        tanks = [
            SlackTank(length=LENGTH, breadth=breadth, height=HEIGHT, fill_volume=FILL / compartments)
            for _ in range(compartments)
        ]
        total = math.fsum(t.upright_inertia for t in tanks)
        assert total == pytest.approx(tank().upright_inertia / compartments**2, rel=1e-15)
        for heel in (0.0, 10.0, 30.0):
            combined = free_surface_gz_shift_curve(tanks, [heel], DISPLACEMENT)[0]
            assert combined == pytest.approx(
                free_surface_gz_shift(tank(), heel, DISPLACEMENT) / compartments**2, rel=1e-12
            )


def test_shift_curve_sums_the_tanks() -> None:
    tanks = [tank(), SlackTank(length=1.0, breadth=0.4, height=0.8, fill_volume=0.1)]
    angles = (0.0, 10.0, 30.0)
    combined = free_surface_gz_shift_curve(tanks, angles, DISPLACEMENT)
    expected = tuple(
        math.fsum(free_surface_gz_shift(t, angle, DISPLACEMENT) for t in tanks) for angle in angles
    )
    assert combined == pytest.approx(expected, rel=1e-15)


def test_module_refuses_bad_geometry_or_volume() -> None:
    t = tank()
    with pytest.raises(FreeSurfaceError):
        free_surface_gz_shift(t, 10.0, 0.0)
    with pytest.raises(FreeSurfaceError):
        free_surface_gz_shift(t, float("inf"), DISPLACEMENT)
    with pytest.raises(FreeSurfaceError):
        free_surface_loss_m([t], -1.0)
    with pytest.raises(FreeSurfaceError):
        tank_free_surface("not a tank", 10.0)  # type: ignore[arg-type]
    with pytest.raises(FreeSurfaceError):
        free_surface_gz_shift_curve(["not a tank"], [0.0], DISPLACEMENT)  # type: ignore[list-item]


def test_integral_reflected_heels_give_the_same_small_angle_loss() -> None:
    """A symmetric tank: the loss at +/- a small heel is the same, the shift is not."""
    t = tank()
    loss = free_surface_loss_m([t], DISPLACEMENT)
    for heel in (0.1, 0.5):
        for sign in (1.0, -1.0):
            # first order in the heel: the deviation from loss*sin is the O(phi^2) term
            assert free_surface_gz_shift(t, sign * heel, DISPLACEMENT) == pytest.approx(
                -sign * loss * math.sin(math.radians(heel)), rel=1e-4
            )


# ------------------------------------------------------------------- integration in gz


@pytest.fixture(scope="module")
def hull():
    return generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8))


def test_upright_symmetric_ship_has_zero_arm_with_a_slack_tank(hull) -> None:
    """Regression: the free-surface correction must not put an arm on an upright ship."""
    inertia = rectangular_tank_inertia(2.0, 0.5)
    raised = righting_arm_curve(hull.mesh, 0.5, [0.0, 20.0, 40.0], free_surface_i_x=inertia)
    assert raised.gz[0] == pytest.approx(0.0, abs=1e-12)
    exact = righting_arm_curve(
        hull.mesh, 0.5, [0.0, 20.0, 40.0], free_surface_model="exact", slack_tanks=[tank()]
    )
    assert exact.gz[0] == pytest.approx(0.0, abs=1e-12)
    assert raised.free_surface_model == "kg_raise"
    assert exact.free_surface_model == "exact"


def test_kg_raise_model_is_the_first_order_correction(hull) -> None:
    inertia = rectangular_tank_inertia(2.0, 0.5)
    solid = righting_arm_curve(hull.mesh, 0.5, [0.0, 20.0, 40.0])
    raised = righting_arm_curve(hull.mesh, 0.5, [0.0, 20.0, 40.0], free_surface_i_x=inertia)
    loss = inertia / solid.displacement_volume
    assert raised.free_surface_correction == pytest.approx(loss, rel=1e-12)
    for angle, gz_solid, gz_raised in zip(solid.angles_deg, solid.gz, raised.gz, strict=True):
        assert gz_raised == pytest.approx(
            gz_solid - loss * math.sin(math.radians(angle)), abs=1e-12, rel=0
        )


def test_flat_model_is_retained_but_wrong_at_zero_heel(hull) -> None:
    inertia = rectangular_tank_inertia(2.0, 0.5)
    solid = righting_arm_curve(hull.mesh, 0.5, [0.0, 20.0, 40.0])
    flat = righting_arm_curve(
        hull.mesh, 0.5, [0.0, 20.0, 40.0], free_surface_i_x=inertia, free_surface_model="flat"
    )
    loss = inertia / solid.displacement_volume
    for gz_solid, gz_flat in zip(solid.gz, flat.gz, strict=True):
        assert gz_flat == pytest.approx(gz_solid - loss, abs=1e-12, rel=0)
    assert flat.gz[0] == pytest.approx(-loss, abs=1e-12)


def test_exact_model_agrees_with_the_first_order_one_only_at_small_heel(hull) -> None:
    angles = (0.1, 10.0, 30.0)
    first_order = righting_arm_curve(
        hull.mesh, 0.5, angles, free_surface_i_x=rectangular_tank_inertia(2.0, 0.5)
    )
    exact = righting_arm_curve(
        hull.mesh, 0.5, angles, free_surface_model="exact", slack_tanks=[tank()]
    )
    small = exact.gz[0] - first_order.gz[0]
    large = exact.gz[2] - first_order.gz[2]
    assert abs(small) < 1e-9
    # at 30 degrees the exact shift must exceed the first-order one by >5 % of itself:
    # the tank's free surface has widened by 1/cos^3(30) = 1.54
    first_order_shift = -first_order.free_surface_correction * math.sin(math.radians(30.0))
    assert large < 0.05 * first_order_shift


def test_exact_model_rejects_conflicting_inputs(hull) -> None:
    with pytest.raises(GzError):
        righting_arm_curve(hull.mesh, 0.5, [10.0], free_surface_model="exact")
    with pytest.raises(GzError):
        righting_arm_curve(
            hull.mesh,
            0.5,
            [10.0],
            free_surface_model="exact",
            slack_tanks=[tank()],
            free_surface_i_x=1.0,
        )
    with pytest.raises(GzError):
        righting_arm_curve(hull.mesh, 0.5, [10.0], slack_tanks=[tank()])
    with pytest.raises(GzError):
        righting_arm_curve(hull.mesh, 0.5, [10.0], free_surface_model="kludge")