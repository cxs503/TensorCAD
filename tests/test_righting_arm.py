"""Righting arms: the heeled volume integral against closed forms and itself.

Three references, in increasing order of how much they can go wrong:

* a **box**, where the wall-sided righting arm ``GZ = sin(phi)(GM + BM tan^2(phi)/2)`` is
  rational arithmetic and the integration has to reproduce it to round-off;
* the **metacentric identity** ``GZ'(0) = KB + BM - KG``, which is exact for any
  polyhedron: the slope of the arm at upright is the waterplane moment divided by the
  volume, and the module measures the left side with a three-dimensional heeled volume
  integral while ``tensorcad.stability`` measures the right side with a two-dimensional
  plane cut. Two code paths, one number;
* the **degenerate cases**, which are the ones that produce plausible-looking nonsense:
  a body asked to displace more than it encloses, a body whose centre of buoyancy cannot
  move because it is fully submerged, and a curve with a zero that is a loll angle rather
  than a vanishing angle.
"""

from __future__ import annotations

import itertools
import math

import pytest

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.righting_arm import (
    DEFAULT_HEEL_ANGLES_DEG,
    DEFAULT_VOLUME_TOLERANCE,
    SCHEMA_RIGHTING_ARM,
    SCHEMA_RIGHTING_ARM_CURVE,
    RightingArm,
    RightingArmCurve,
    RightingArmError,
    heeled_displacement,
    heeled_equilibrium,
    is_watertight,
    righting_arm,
    righting_arm_curve,
    wall_sided_righting_arm,
    wall_sided_vanishing_angle,
)
from tensorcad.stability import waterplane_moments
from tensorcad.wigley import generate_wigley

BOX_LENGTH = 20.0
BOX_BEAM = 10.0
BOX_DEPTH = 10.0
BOX_DRAFT = 5.0
BOX_KG = 4.0
BOX_GM = BOX_DRAFT / 2.0 + BOX_BEAM**2 / (12.0 * BOX_DRAFT) - BOX_KG
BOX_KG_UNSTABLE = BOX_DRAFT / 2.0 + BOX_BEAM**2 / (12.0 * BOX_DRAFT) + 0.2


def box(length: float, beam: float, depth: float) -> Mesh:
    """Axis-aligned closed box, top face on ``z = 0``, keel on ``z = -depth``."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-depth, 0.0):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    left, right = 0.0, length
    lo, hi = -0.5 * beam, 0.5 * beam
    top, bottom = 0.0, -depth
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


def inverted_box(length: float, beam: float, depth: float) -> Mesh:
    """The same box with every winding reversed, so the signed volume comes out negative."""
    mesh = box(length, beam, depth)
    return Mesh(
        vertices=mesh.vertices,
        faces=tuple((a, c, b) for a, b, c in mesh.faces),
        units=mesh.units,
        coordinate_frame=mesh.coordinate_frame,
    )


@pytest.fixture(scope="module")
def hull_mesh() -> Mesh:
    return generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8)).mesh


@pytest.fixture(scope="module")
def hull_displacement() -> float:
    return 3.75


# ---------------------------------------------------------------------------
# the box, where everything has a closed form
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("heel", [0.0, 2.5, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0])
def test_box_matches_the_wall_sided_closed_form(heel: float) -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    arm = righting_arm(
        mesh, heel_angle_deg=heel, kg=BOX_KG, displacement_volume=BOX_LENGTH * BOX_BEAM * BOX_DRAFT
    )
    closed = wall_sided_righting_arm(
        heel, beam=BOX_BEAM, depth=BOX_DEPTH, draft=BOX_DRAFT, kg=BOX_KG
    )
    assert arm.gz == pytest.approx(closed, abs=1e-12)


def test_box_waterline_does_not_move_because_the_section_is_wall_sided() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    displacement = BOX_LENGTH * BOX_BEAM * BOX_DRAFT
    waterlines = [
        heeled_equilibrium(mesh, heel, displacement).waterline_z for heel in (0.0, 10.0, 30.0)
    ]
    for waterline in waterlines:
        assert waterline == pytest.approx(-BOX_DRAFT, abs=1e-9)


def test_box_centre_of_buoyancy_moves_to_the_low_side() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    equilibrium = heeled_equilibrium(mesh, 20.0, BOX_LENGTH * BOX_BEAM * BOX_DRAFT)
    _, y_b, z_b = equilibrium.centre_of_buoyancy
    bm = BOX_BEAM**2 / (12.0 * BOX_DRAFT)
    assert y_b == pytest.approx(bm * math.tan(math.radians(20.0)), rel=1e-12)
    assert z_b == pytest.approx(
        -BOX_DEPTH + BOX_DRAFT / 2.0 + 0.5 * bm * math.tan(math.radians(20.0)) ** 2, rel=1e-12
    )


def test_negative_gm_box_vanishes_at_the_closed_form_angle() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    curve = righting_arm_curve(
        mesh,
        kg=BOX_KG_UNSTABLE,
        displacement_volume=BOX_LENGTH * BOX_BEAM * BOX_DRAFT,
        angles_deg=(0.0, 5.0, 10.0, 20.0, 30.0),
    )
    analytic = wall_sided_vanishing_angle(
        beam=BOX_BEAM, depth=BOX_DEPTH, draft=BOX_DRAFT, kg=BOX_KG_UNSTABLE
    )
    assert analytic is not None
    assert curve.zero_crossings_deg
    assert curve.zero_crossings_deg[0] == pytest.approx(analytic, abs=1e-6)
    # A curve that starts non-positive has no range of stability to end.
    assert curve.vanishing_angle_deg is None


def test_positive_gm_box_has_no_zero_inside_the_scanned_range() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    curve = righting_arm_curve(
        mesh,
        kg=BOX_KG,
        displacement_volume=BOX_LENGTH * BOX_BEAM * BOX_DRAFT,
        angles_deg=(0.0, 10.0, 20.0, 30.0, 40.0),
    )
    assert curve.zero_crossings_deg == ()
    assert curve.vanishing_angle_deg is None
    assert curve.maximum_gz > 0.0


def test_dynamical_stability_matches_the_integral_of_the_closed_form() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    curve = righting_arm_curve(
        mesh,
        kg=BOX_KG_UNSTABLE,
        displacement_volume=BOX_LENGTH * BOX_BEAM * BOX_DRAFT,
        angles_deg=tuple(0.5 * step for step in range(0, 51)),
    )
    upper = math.radians(25.0)
    gm = BOX_DRAFT / 2.0 + BOX_BEAM**2 / (12.0 * BOX_DRAFT) - BOX_KG_UNSTABLE
    bm = BOX_BEAM**2 / (12.0 * BOX_DRAFT)
    analytic = -gm * (math.cos(upper) - 1.0) + 0.5 * bm * (
        1.0 / math.cos(upper) + math.cos(upper) - 2.0
    )
    assert curve.dynamical_stability(25.0) == pytest.approx(analytic, rel=1e-3)


def test_box_heel_antisymmetry() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    displacement = BOX_LENGTH * BOX_BEAM * BOX_DRAFT
    starboard = righting_arm(mesh, heel_angle_deg=15.0, kg=BOX_KG, displacement_volume=displacement)
    port = righting_arm(mesh, heel_angle_deg=-15.0, kg=BOX_KG, displacement_volume=displacement)
    assert starboard.gz == pytest.approx(-port.gz, abs=1e-12)


# ---------------------------------------------------------------------------
# the metacentric identity, which is exact for any polyhedron
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("kg", [0.35, 0.55, 0.75])
def test_hull_slope_at_upright_is_the_metacentric_height(
    hull_mesh: Mesh, hull_displacement: float, kg: float
) -> None:
    curve = righting_arm_curve(
        hull_mesh, kg=kg, displacement_volume=hull_displacement, angles_deg=(0.05, 0.1, 1.0)
    )
    upright = righting_arm(hull_mesh, heel_angle_deg=0.0, kg=kg, displacement_volume=hull_displacement)
    section = waterplane_moments(hull_mesh, upright.waterline_z)
    identity = (
        upright.centre_of_buoyancy[2]
        - upright.centre_of_gravity_z
        + section.transverse_second_moment / hull_displacement
    )
    assert curve.initial_slope == pytest.approx(identity, rel=1e-6)


def test_wigley_slope_at_upright_is_the_metacentric_height() -> None:
    mesh = generate_wigley().mesh
    displacement = 0.5 * compute_hydrostatics(mesh, 0.0).displacement_volume
    curve = righting_arm_curve(
        mesh, kg=0.35, displacement_volume=displacement, angles_deg=(0.05, 0.1, 2.0)
    )
    upright = righting_arm(mesh, heel_angle_deg=0.0, kg=0.35, displacement_volume=displacement)
    section = waterplane_moments(mesh, upright.waterline_z)
    identity = (
        upright.centre_of_buoyancy[2]
        - upright.centre_of_gravity_z
        + section.transverse_second_moment / displacement
    )
    assert curve.initial_slope == pytest.approx(identity, rel=1e-6)


def test_hull_arm_is_odd_in_the_heel_and_zero_upright(
    hull_mesh: Mesh, hull_displacement: float
) -> None:
    upright = righting_arm(hull_mesh, heel_angle_deg=0.0, kg=0.55, displacement_volume=hull_displacement)
    assert abs(upright.gz) < 1e-12
    starboard = righting_arm(hull_mesh, heel_angle_deg=10.0, kg=0.55, displacement_volume=hull_displacement)
    port = righting_arm(hull_mesh, heel_angle_deg=-10.0, kg=0.55, displacement_volume=hull_displacement)
    assert starboard.gz == pytest.approx(-port.gz, abs=1e-12)
    assert starboard.gz > 0.0


def test_heeled_body_sinks_to_keep_its_displacement(
    hull_mesh: Mesh, hull_displacement: float
) -> None:
    waterlines = [
        righting_arm(hull_mesh, heel_angle_deg=heel, kg=0.55, displacement_volume=hull_displacement)
        .waterline_z
        for heel in (0.0, 5.0, 10.0, 20.0)
    ]
    assert all(left > right for left, right in itertools.pairwise(waterlines))


def test_equilibrium_reproduces_the_requested_displacement(
    hull_mesh: Mesh, hull_displacement: float
) -> None:
    equilibrium = heeled_equilibrium(hull_mesh, 12.5, hull_displacement)
    volume, _ = heeled_displacement(hull_mesh, 12.5, equilibrium.waterline_z)
    assert volume == pytest.approx(hull_displacement, rel=DEFAULT_VOLUME_TOLERANCE * 100)
    assert abs(equilibrium.residual) <= DEFAULT_VOLUME_TOLERANCE * hull_displacement * 1e3


# ---------------------------------------------------------------------------
# the honest limits
# ---------------------------------------------------------------------------
def test_a_body_that_cannot_float_that_heavily_is_refused(hull_mesh: Mesh) -> None:
    enclosed = compute_hydrostatics(hull_mesh, 0.0).displacement_volume
    with pytest.raises(RightingArmError, match="fully submerged"):
        heeled_equilibrium(hull_mesh, 0.0, enclosed)
    with pytest.raises(RightingArmError, match="fully submerged"):
        heeled_equilibrium(hull_mesh, 10.0, enclosed * 1.5)


def test_inverted_orientation_is_refused_rather_than_signed() -> None:
    mesh = inverted_box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    with pytest.raises(RightingArmError, match="inverted"):
        heeled_equilibrium(mesh, 0.0, BOX_LENGTH * BOX_BEAM * BOX_DRAFT)


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_non_finite_inputs_are_refused(hull_mesh: Mesh, bad: float) -> None:
    with pytest.raises(RightingArmError):
        righting_arm(hull_mesh, heel_angle_deg=bad, kg=0.55, displacement_volume=3.0)
    with pytest.raises(RightingArmError):
        righting_arm(hull_mesh, heel_angle_deg=5.0, kg=bad, displacement_volume=3.0)
    with pytest.raises(RightingArmError):
        righting_arm(hull_mesh, heel_angle_deg=5.0, kg=0.55, displacement_volume=bad)


def test_wall_sided_formula_refuses_angles_outside_its_assumptions() -> None:
    with pytest.raises(RightingArmError, match="wall-sided"):
        wall_sided_righting_arm(50.0, beam=BOX_BEAM, depth=BOX_DEPTH, draft=BOX_DRAFT, kg=BOX_KG)
    with pytest.raises(RightingArmError):
        wall_sided_righting_arm(5.0, beam=BOX_BEAM, depth=BOX_DEPTH, draft=BOX_DEPTH, kg=BOX_KG)


def test_vanishing_angle_only_exists_for_negative_gm() -> None:
    assert (
        wall_sided_vanishing_angle(beam=BOX_BEAM, depth=BOX_DEPTH, draft=BOX_DRAFT, kg=BOX_KG)
        is None
    )
    assert wall_sided_vanishing_angle(
        beam=BOX_BEAM, depth=BOX_DEPTH, draft=BOX_DRAFT, kg=BOX_KG_UNSTABLE
    ) == pytest.approx(26.1001, abs=1e-4)


def test_deck_rim_limit_is_reported_not_extrapolated(
    hull_mesh: Mesh, hull_displacement: float
) -> None:
    curve = righting_arm_curve(
        hull_mesh, kg=0.55, displacement_volume=hull_displacement, angles_deg=(0.05, 0.1, 10.0, 20.0)
    )
    limit = curve.physical_limit_deg
    assert limit is not None
    assert 10.0 < limit < 20.0
    assert not curve.at(10.0).deck_submerged
    assert curve.at(20.0).deck_submerged


# ---------------------------------------------------------------------------
# shape of the published objects
# ---------------------------------------------------------------------------
def test_curve_reports_its_landmarks(hull_mesh: Mesh, hull_displacement: float) -> None:
    curve = righting_arm_curve(hull_mesh, kg=0.55, displacement_volume=hull_displacement)
    assert isinstance(curve, RightingArmCurve)
    assert all(isinstance(point, RightingArm) for point in curve.points)
    assert curve.points[0].heel_angle_deg == DEFAULT_HEEL_ANGLES_DEG[0]
    assert curve.maximum_gz == max(point.gz for point in curve.points)
    payload = curve.to_dict()
    assert payload["schema"] == SCHEMA_RIGHTING_ARM_CURVE
    assert payload["points"][1]["basis"]
    assert isinstance(payload["zero_crossings_deg"], list)
    assert "vanishing_angle_deg" in payload


def test_point_dict_carries_the_frame_it_was_measured_in(
    hull_mesh: Mesh, hull_displacement: float
) -> None:
    point = righting_arm(hull_mesh, heel_angle_deg=10.0, kg=0.55, displacement_volume=hull_displacement)
    payload = point.to_dict()
    assert payload["heel_angle_deg"] == 10.0
    assert len(payload["centre_of_buoyancy_m"]) == 3
    assert payload["centre_of_gravity_z_m"] == pytest.approx(
        point.centre_of_gravity_z
    )
    assert payload["centre_of_gravity_z_m"] == pytest.approx(
        min(vertex[2] for vertex in hull_mesh.vertices) + 0.55
    )
    assert SCHEMA_RIGHTING_ARM.startswith("tensorcad.righting_arm")


def test_curve_lookup_and_errors(hull_mesh: Mesh, hull_displacement: float) -> None:
    curve = righting_arm_curve(
        hull_mesh, kg=0.55, displacement_volume=hull_displacement, angles_deg=(0.0, 5.0)
    )
    assert curve.at(5.0).heel_angle_deg == 5.0
    with pytest.raises(KeyError):
        curve.at(7.5)
    with pytest.raises(RightingArmError):
        righting_arm_curve(hull_mesh, kg=0.55, displacement_volume=hull_displacement, angles_deg=())


def test_watertightness_is_the_assumption_the_volume_needs(hull_mesh: Mesh) -> None:
    assert is_watertight(box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH))
    assert is_watertight(hull_mesh)
