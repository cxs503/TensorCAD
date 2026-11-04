"""Large-angle stability: ``GZ`` against closed forms and one independent integration.

The claim under test is not "the curve looks like a righting-arm curve" but that the
heeled displacement is *exact for the polytope* and that the geometry-to-arm step is
signed correctly. Both are pinned three ways: a closed-form wall-sided body, a section
integration that clips in two dimensions instead of three, and a sphere whose arm is
``(KM - KG) sin(phi)`` at every heel.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.geometry_field import SphereField, to_mesh
from tensorcad.gz import (
    HeeledDisplacement,
    HeelError,
    RightingArm,
    barge_deck_edge_angle,
    barge_righting_arm,
    barge_section,
    gz_curve,
    heeled_displacement,
    keel_z,
    righting_arm,
    upright_displacement,
    up_direction,
    wall_sided_dynamical_stability,
    wall_sided_metacentric_height,
    wall_sided_range,
    wall_sided_righting_arm,
    waterplane_direction,
)
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.stability import initial_stability

BEAM, DRAFT, LENGTH = 2.0, 0.8, 10.0
#: A barge that can actually be heeled: 0.4 m of freeboard above the 0.8 m draft.
DEPTH = 1.2
#: A barge whose deck edge (58.0 deg) is *beyond* the wall-sided limit (38.7 deg), so
#: the rectangular-section closed forms apply over their whole stated range. The
#: shallower barge leaves the wall-sided regime at its deck edge, 21.8 deg.
DEEP_DEPTH = 2.4
#: ``KM`` of the section above: ``T/2 + B^2/(12 T) = 0.4 + 0.41667``.
STABLE_KG = 0.3
UNSTABLE_KG = 0.9


# ---------------------------------------------------------------------------
# reference geometries
# ---------------------------------------------------------------------------
def barge(
    length: float = LENGTH,
    beam: float = BEAM,
    draft: float = DRAFT,
    depth: float = DEPTH,
    *,
    x0: float = 0.0,
) -> Mesh:
    """Closed box, deck at ``z = depth - draft``, waterline at ``z = 0``, keel at ``-draft``."""
    builder = MeshBuilder()
    corner = {}
    for x in (x0, x0 + length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, -draft + depth):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    left, right = x0, x0 + length
    lo, hi = -0.5 * beam, 0.5 * beam
    bottom, top = -draft, -draft + depth
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


def sphere_mesh(resolution: tuple[int, int, int] = (20, 20, 20), radius: float = 1.0) -> Mesh:
    return to_mesh(SphereField(radius=radius), resolution).mesh


def degrees_to_radians(degrees: float) -> float:
    return math.radians(degrees)


# ---------------------------------------------------------------------------
# the frame
# ---------------------------------------------------------------------------
def test_up_direction_is_transverse_and_orthogonal_to_the_waterplane() -> None:
    for degrees in (0.0, 17.0, 45.0, 90.0):
        angle = degrees_to_radians(degrees)
        up = up_direction(angle)
        in_plane = waterplane_direction(angle)
        assert up[0] == 0.0 and in_plane[0] == 0.0
        assert math.isclose(math.hypot(up[1], up[2]), 1.0, abs_tol=1e-15)
        assert math.isclose(sum(a * b for a, b in zip(up, in_plane)), 0.0, abs_tol=1e-15)


def test_positive_heel_puts_starboard_down() -> None:
    # Up tipping towards -y is starboard down: otherwise GZ would come out negated.
    assert up_direction(math.radians(10.0))[1] < 0.0


# ---------------------------------------------------------------------------
# the upright limit is already the hydrostatics of the package
# ---------------------------------------------------------------------------
def test_upright_volume_matches_hydrostatics() -> None:
    mesh = barge(depth=DRAFT)  # capped at the waterline: the integrator's own regime
    reference = compute_hydrostatics(mesh, 0.0)
    volume, _, y_b, z_b = upright_displacement(mesh, 0.0)
    assert math.isclose(volume, reference.displacement_volume, rel_tol=1e-14)
    assert math.isclose(z_b, reference.vertical_centre_buoyancy, rel_tol=1e-14)
    assert abs(y_b) < 1e-15


def test_upright_displacement_ignores_freeboard() -> None:
    # The half-space slice does not care that the barge has a deck 0.4 m above water,
    # which is exactly where it differs from compute_hydrostatics (closed-volume).
    capped = upright_displacement(barge(depth=DRAFT), 0.0)[0]
    decked = upright_displacement(barge(depth=DEPTH), 0.0)[0]
    assert math.isclose(decked, capped, rel_tol=1e-15)
    assert math.isclose(decked, LENGTH * BEAM * DRAFT, rel_tol=1e-15)


def test_upright_righting_arm_is_zero() -> None:
    mesh = barge()
    arm = righting_arm(mesh, 0.0, STABLE_KG)
    assert abs(arm.righting_arm) < 1e-15
    assert math.isnan(arm.gz_over_sin)


# ---------------------------------------------------------------------------
# the heeled displacement itself
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("degrees", [5.0, 20.0, 45.0, 70.0, 85.0])
def test_heeled_volume_is_conserved(degrees: float) -> None:
    point = righting_arm(barge(), degrees_to_radians(degrees), STABLE_KG)
    assert abs(point.relative_volume_error) < 1e-13
    assert point.iterations <= 10


def test_heeled_displacement_reports_the_loading() -> None:
    heeled = heeled_displacement(barge(), degrees_to_radians(30.0))
    assert isinstance(heeled, HeeledDisplacement)
    assert heeled.target_volume == pytest.approx(LENGTH * BEAM * DRAFT, rel=1e-15)
    assert heeled.kb == pytest.approx(heeled.centre_of_buoyancy[2] - heeled.keel_z)
    assert math.isclose(heeled.keel_z, -DRAFT, abs_tol=1e-15)
    assert math.isclose(keel_z(barge()), -DRAFT, abs_tol=1e-15)


def test_a_mesh_that_is_entirely_above_the_water_raises() -> None:
    floating = barge()
    lifted = Mesh(
        vertices=tuple((x, y, z + 5.0) for x, y, z in floating.vertices),
        faces=floating.faces,
    )
    with pytest.raises(HeelError, match="displaces no volume"):
        upright_displacement(lifted, 0.0)


# ---------------------------------------------------------------------------
# closed form 1: the wall-sided body
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("degrees", [2.0, 5.0, 10.0, 20.0, 30.0, 35.0])
def test_box_matches_the_wall_sided_closed_form(degrees: float) -> None:
    angle = degrees_to_radians(degrees)
    measured = righting_arm(barge(depth=DEEP_DEPTH), angle, STABLE_KG).righting_arm
    closed = wall_sided_righting_arm(BEAM, DRAFT, STABLE_KG, angle)
    assert measured == pytest.approx(closed, rel=1e-13, abs=1e-15)


def test_the_wall_sided_formula_is_not_a_small_angle_truncation() -> None:
    # At 30 deg the tan^2 term is 13 % of the arm, so the closed form used as the
    # large-angle reference is not the small-angle one in disguise.
    angle = degrees_to_radians(30.0)
    gm = wall_sided_metacentric_height(BEAM, DRAFT, STABLE_KG)
    bm = BEAM * BEAM / (12.0 * DRAFT)
    full = wall_sided_righting_arm(BEAM, DRAFT, STABLE_KG, angle)
    small = math.sin(angle) * gm
    assert full / small == pytest.approx(1.0 + 0.5 * bm * math.tan(angle) ** 2 / gm, rel=1e-12)
    assert full > 1.1 * small
    assert math.isclose(full, math.sin(angle) * (gm + 0.5 * bm * math.tan(angle) ** 2))


def test_the_wall_sided_formula_reduces_to_gm_sin_at_small_heel() -> None:
    angle = degrees_to_radians(0.01)
    gm = wall_sided_metacentric_height(BEAM, DRAFT, STABLE_KG)
    assert math.isclose(
        wall_sided_righting_arm(BEAM, DRAFT, STABLE_KG, angle) / math.sin(angle),
        gm,
        rel_tol=1e-6,
    )


def test_the_wall_sided_formula_refuses_a_heel_past_its_range() -> None:
    limit = wall_sided_range(BEAM, DRAFT)
    assert math.isclose(math.degrees(limit), math.degrees(math.atan(2 * DRAFT / BEAM)))
    with pytest.raises(HeelError, match="wall-sided range"):
        wall_sided_righting_arm(BEAM, DRAFT, STABLE_KG, limit)


def test_the_metacentric_height_agrees_with_initial_stability() -> None:
    # Two different modules, two different integrations, one number.
    reference = initial_stability(barge(depth=DRAFT), kg=STABLE_KG)
    assert math.isclose(
        reference.gm_transverse,
        wall_sided_metacentric_height(BEAM, DRAFT, STABLE_KG),
        rel_tol=1e-12,
    )


def test_gm_from_the_curve_converges_to_the_metacentric_height() -> None:
    mesh = barge()
    exact = wall_sided_metacentric_height(BEAM, DRAFT, STABLE_KG)
    errors = []
    for degrees in (4.0, 2.0, 1.0, 0.5):
        arm = righting_arm(mesh, degrees_to_radians(degrees), STABLE_KG).righting_arm
        errors.append(abs(arm / math.sin(degrees_to_radians(degrees)) - exact))
    #  GZ/sin - GM ~ (BM/2) tan^2(phi): second order, so each halving quarters it.
    for coarse, fine in zip(errors, errors[1:]):
        assert 3.4 < coarse / fine < 4.1


# ---------------------------------------------------------------------------
# closed form 2: an independent two-dimensional integration of the same barge
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("degrees", [0.0, 10.0, 25.0, 40.0, 55.0, 70.0, 90.0])
def test_barge_matches_the_independent_section_integration(degrees: float) -> None:
    angle = degrees_to_radians(degrees)
    mesh_arm = righting_arm(barge(), angle, STABLE_KG).righting_arm
    section_arm = barge_righting_arm(BEAM, DRAFT, DEPTH, STABLE_KG, angle)
    assert mesh_arm == pytest.approx(section_arm, abs=2e-13)


def test_heeled_volume_of_the_section_is_the_upright_area() -> None:
    for degrees in (0.0, 15.0, 45.0, 80.0):
        section = barge_section(BEAM, DRAFT, DEPTH, degrees_to_radians(degrees))
        assert math.isclose(section.area, BEAM * DRAFT, rel_tol=1e-14)


def test_upright_section_has_the_textbook_centroid() -> None:
    section = barge_section(BEAM, DRAFT, DEPTH, 0.0)
    assert abs(section.y_centroid) < 1e-15
    assert math.isclose(section.z_centroid, -DRAFT + 0.5 * DRAFT, abs_tol=1e-15)


def test_a_section_without_freeboard_is_refused() -> None:
    with pytest.raises(HeelError, match="no freeboard"):
        barge_section(BEAM, DRAFT, DRAFT, degrees_to_radians(10.0))


def test_deck_edge_angle() -> None:
    assert math.isclose(
        math.degrees(barge_deck_edge_angle(BEAM, DRAFT, DEPTH)),
        math.degrees(math.atan(2 * (DEPTH - DRAFT) / BEAM)),
    )


# ---------------------------------------------------------------------------
# the curve: peak, loll, vanishing angle, area
# ---------------------------------------------------------------------------
def test_angle_of_loll_matches_the_wall_sided_closed_form() -> None:
    # GM < 0: the upright ship is unstable and settles at tan^2(phi) = -2 GM / BM.
    gm = wall_sided_metacentric_height(BEAM, DRAFT, UNSTABLE_KG)
    bm = BEAM * BEAM / (12.0 * DRAFT)
    assert gm < 0.0
    expected = math.degrees(math.atan(math.sqrt(-2.0 * gm / bm)))
    curve = gz_curve(
        barge(depth=DEEP_DEPTH), [float(angle) for angle in range(0, 91, 5)], UNSTABLE_KG
    )
    assert curve.loll_angle_deg == pytest.approx(expected, abs=1e-6)
    assert curve.vanishing_angle_deg is None
    assert curve.arms[0] == 0.0
    assert all(arm < 0.0 for arm in curve.arms[1:6])


def test_a_stable_ship_has_no_angle_of_loll() -> None:
    curve = gz_curve(barge(), [float(angle) for angle in range(0, 91, 5)], STABLE_KG)
    assert curve.loll_angle_deg is None


def test_the_vanishing_angle_is_a_zero_of_the_arm() -> None:
    mesh = barge()
    curve = gz_curve(mesh, [float(angle) for angle in range(0, 91, 5)], 0.78)
    assert curve.vanishing_angle_deg is not None
    angle = degrees_to_radians(curve.vanishing_angle_deg)
    assert abs(righting_arm(mesh, angle, 0.78).righting_arm) < 1e-12
    # and it is a genuine crossing, not a tangent
    delta = degrees_to_radians(1.0)
    assert righting_arm(mesh, angle - delta, 0.78).righting_arm > 1e-3
    assert righting_arm(mesh, angle + delta, 0.78).righting_arm < -1e-3


def test_a_very_stable_barge_has_no_vanishing_angle_below_90_degrees() -> None:
    curve = gz_curve(barge(depth=2.4), [float(angle) for angle in range(0, 91, 5)], 0.7)
    assert curve.vanishing_angle_deg is None
    assert curve.angle_of_maximum_deg == 90.0
    assert curve.maximum_arm > 0.0


def test_dynamical_stability_matches_the_closed_form_area() -> None:
    # Inside the wall-sided range the arm integrates exactly, so a quadrature of the
    # measured curve has something to be wrong against.
    exact = wall_sided_dynamical_stability(BEAM, DRAFT, STABLE_KG, degrees_to_radians(30.0))
    coarse = gz_curve(
        barge(depth=DEEP_DEPTH), [float(angle) for angle in range(0, 31, 5)], STABLE_KG
    )
    fine = gz_curve(barge(depth=DEEP_DEPTH), [2.5 * index for index in range(0, 13)], STABLE_KG)
    assert coarse.integration_rule == "composite Simpson"
    assert math.isclose(coarse.dynamical_stability_m_rad, exact, rel_tol=1e-4)
    # Simpson is fourth order: halving the step takes the quadrature error down hard.
    assert abs(fine.dynamical_stability_m_rad - exact) < 0.1 * abs(
        coarse.dynamical_stability_m_rad - exact
    )


def test_dynamical_stability_stops_at_the_vanishing_angle() -> None:
    curve = gz_curve(barge(), [float(angle) for angle in range(0, 91, 5)], 0.78)
    assert curve.integrated_to_deg == pytest.approx(curve.vanishing_angle_deg, abs=1e-12)
    assert curve.dynamical_stability_m_rad > 0.0
    assert curve.integration_rule.endswith("tail")


def test_the_curve_carries_both_integration_rules() -> None:
    odd = gz_curve(barge(), [0.0, 10.0, 20.0], STABLE_KG)
    even = gz_curve(barge(), [0.0, 10.0, 20.0, 30.0], STABLE_KG)
    assert odd.integration_rule == "composite Simpson"
    assert even.integration_rule == "trapezoid"


# ---------------------------------------------------------------------------
# the two linear laws a righting arm must obey
# ---------------------------------------------------------------------------
def test_the_arm_is_linear_in_kg() -> None:
    angle = degrees_to_radians(40.0)
    base = righting_arm(barge(), angle, STABLE_KG).righting_arm
    for shift in (0.05, 0.2):
        shifted = righting_arm(barge(), angle, STABLE_KG + shift).righting_arm
        assert math.isclose(shifted - base, -shift * math.sin(angle), abs_tol=1e-13)


def test_the_arm_is_linear_in_yg() -> None:
    angle = degrees_to_radians(25.0)
    base = righting_arm(barge(), angle, STABLE_KG).righting_arm
    shifted = righting_arm(barge(), angle, STABLE_KG, yg=0.1).righting_arm
    assert math.isclose(shifted - base, -0.1 * math.cos(angle), abs_tol=1e-14)


def test_the_arm_survives_a_translation_of_the_frame() -> None:
    mesh = barge()
    lifted = Mesh(
        vertices=tuple((x, y, z - 1.0) for x, y, z in mesh.vertices),
        faces=mesh.faces,
    )
    assert math.isclose(
        righting_arm(mesh, degrees_to_radians(37.0), STABLE_KG).righting_arm,
        righting_arm(lifted, degrees_to_radians(37.0), STABLE_KG, waterline_z=-1.0).righting_arm,
        abs_tol=1e-14,
    )


# ---------------------------------------------------------------------------
# a second analytic body: the sphere has GZ = (KM - KG) sin(phi) at every heel
# ---------------------------------------------------------------------------
KG_SPHERE = 0.5


@pytest.mark.parametrize("degrees", [15.0, 30.0, 45.0, 60.0, 75.0, 90.0])
def test_sphere_matches_gm_sin_at_large_heel(degrees: float) -> None:
    mesh = sphere_mesh((20, 20, 20), radius=1.0)
    angle = degrees_to_radians(degrees)
    measured = righting_arm(mesh, angle, KG_SPHERE).righting_arm
    exact = (1.0 - KG_SPHERE) * math.sin(angle)
    assert measured == pytest.approx(exact, rel=4e-3)


def test_sphere_error_falls_with_the_tessellation() -> None:
    errors = []
    for resolution in ((12, 12, 12), (20, 20, 20)):
        mesh = sphere_mesh(resolution, radius=1.0)
        worst = 0.0
        for degrees in (15.0, 45.0, 75.0):
            angle = degrees_to_radians(degrees)
            measured = righting_arm(mesh, angle, KG_SPHERE).righting_arm
            exact = (1.0 - KG_SPHERE) * math.sin(angle)
            worst = max(worst, abs(measured - exact) / exact)
        errors.append(worst)
    assert errors[0] > 2.0 * errors[1]
    assert errors[1] < 5e-3


def test_sphere_conserves_its_heeled_volume() -> None:
    mesh = sphere_mesh((12, 12, 12), radius=1.0)
    for degrees in (20.0, 50.0, 80.0):
        point = righting_arm(mesh, degrees_to_radians(degrees), KG_SPHERE)
        assert abs(point.relative_volume_error) < 1e-13


# ---------------------------------------------------------------------------
# the trap: a mesh capped at the waterline cannot be heeled
# ---------------------------------------------------------------------------
def test_a_mesh_capped_at_the_waterline_is_refused_past_zero_heel() -> None:
    with pytest.raises(HeelError, match="freeboard"):
        righting_arm(barge(depth=DRAFT), degrees_to_radians(5.0), STABLE_KG)


def test_an_upright_query_of_a_capped_mesh_is_still_allowed() -> None:
    heeled = heeled_displacement(barge(depth=DRAFT), 0.0)
    assert math.isclose(heeled.displacement_volume, LENGTH * BEAM * DRAFT, rel_tol=1e-15)


# ---------------------------------------------------------------------------
# argument handling and serialisation
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("angle", [math.radians(91.0), math.radians(-91.0)])
def test_heel_past_ninety_degrees_raises(angle: float) -> None:
    with pytest.raises(HeelError, match="past 90"):
        righting_arm(barge(), angle, STABLE_KG)


@pytest.mark.parametrize("angle", [float("nan"), float("inf")])
def test_a_non_finite_heel_raises(angle: float) -> None:
    with pytest.raises(HeelError, match="finite"):
        righting_arm(barge(), angle, STABLE_KG)


def test_bad_curve_arguments_raise() -> None:
    with pytest.raises(HeelError, match="at least one"):
        gz_curve(barge(), [])
    with pytest.raises(HeelError, match="non-negative"):
        gz_curve(barge(), [-5.0, 0.0, 5.0], STABLE_KG)
    with pytest.raises(HeelError, match="KG must be finite"):
        righting_arm(barge(), 0.1, float("nan"))


def test_a_non_positive_target_volume_raises() -> None:
    with pytest.raises(HeelError, match="positive"):
        heeled_displacement(barge(), 0.1, displacement_volume=0.0)


def test_a_mesh_without_faces_raises() -> None:
    with pytest.raises(HeelError, match="no faces"):
        heeled_displacement(Mesh(vertices=(), faces=()), 0.1)


def test_records_serialise_with_their_schemas() -> None:
    angle = degrees_to_radians(20.0)
    arm = righting_arm(barge(), angle, STABLE_KG)
    assert isinstance(arm, RightingArm)
    payload = arm.to_dict()
    assert payload["schema"] == "tensorcad.stability.righting_arm/1"
    assert payload["righting_arm_m"] == pytest.approx(arm.righting_arm, abs=1e-12)
    assert arm.gz_over_sin == pytest.approx(
        arm.righting_arm / math.sin(angle), rel=1e-14
    )
    curve = gz_curve(barge(), [0.0, 10.0, 20.0, 30.0], STABLE_KG)
    curve_payload = curve.to_dict()
    assert curve_payload["schema"] == "tensorcad.stability.gz_curve/1"
    assert len(curve_payload["righting_arms_m"]) == 4
    heeled = heeled_displacement(barge(), angle).to_dict()
    assert heeled["schema"] == "tensorcad.stability.heeled/1"
    assert heeled["iterations"] >= 0