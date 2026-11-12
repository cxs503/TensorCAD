"""Tests for the exact prismatic-section righting-arm oracle (:mod:`tensorcad.section_gz`).

The oracle exists to close a specific hole: every exact large-angle reference in the package
before it was either wall-sided (a box) or *centrally symmetric* and cut at half draft (the
log).  These tests therefore check the oracle in three directions at once -- against the mesh
integrator in :mod:`tensorcad.gz` (machine precision), against the closed forms that already
exist (the wall-sided box and the log, so two independent implementations agree), and against
the continuum limit of its own chord refinement, which is the only way a *new* closed form can
be shown not to be a plausible-looking mistake.
"""

from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad import section_gz
from tensorcad.geometry import Mesh
from tensorcad.gz import (
    box_righting_arms_analytic,
    box_wall_sided_limit_deg,
    righting_arm_curve,
)

TRIANGLE = ((-1.0, 0.0), (1.0, 0.0), (0.0, 2.0))
ANGLES = [0.0, 5.0, 10.0, 20.0, 30.0, 45.0, 60.0, 75.0, 90.0]


def _ellipse_curve(angles=ANGLES, *, semi_beam=2.0, semi_height=1.0, draft=1.2, kg=0.8):
    return section_gz.ellipse_prism_righting_arms(
        semi_beam, semi_height, angles, kg=kg, waterline_z=draft, length=4.0
    )


def test_clipped_section_properties_are_exact() -> None:
    """Half of a unit square: area 0.5, centroid (0.5, 0.25) -- no discretisation anywhere."""
    square = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
    area, cy, cz = section_gz.section_properties(square, 0.5)
    assert area == pytest.approx(0.5, rel=0, abs=1e-16)
    assert cy == pytest.approx(0.5, rel=0, abs=1e-16)
    assert cz == pytest.approx(0.25, rel=0, abs=1e-16)


def test_prism_mesh_is_closed_and_has_the_expected_volume() -> None:
    """The mesh the comparison rests on is closed, outward-oriented and the right size."""
    section = TRIANGLE
    length = 4.0
    mesh = section_gz.prism_mesh(section, length)
    assert meshlib.is_watertight(mesh)
    assert meshlib.is_edge_manifold(mesh)
    assert not meshlib.zero_area_faces(mesh, 0.0)
    area = section_gz.section_properties(section, 10.0)[0]
    assert meshlib.total_signed_volume(mesh) == pytest.approx(area * length, rel=1e-14)


def test_polygon_oracle_matches_the_mesh_integrator_at_every_heel() -> None:
    """A V-shaped section: the mesh integration and the exact section agree to round-off.

    This is the case the existing references cannot reach -- the section is neither wall-sided
    nor centrally symmetric, so ``GZ`` is not a sine and the wall-sided formula does not apply.
    """
    mesh = section_gz.prism_mesh(TRIANGLE, 4.0)
    exact = section_gz.polygon_prism_righting_arms(
        TRIANGLE, ANGLES, kg=0.6, waterline_z=1.0, length=4.0
    )
    curve = righting_arm_curve(mesh, 0.6, ANGLES, waterline_z=1.0)
    assert curve.max_volume_error < 1e-12
    for arm, gz in zip(exact.arms, curve.gz, strict=True):
        assert gz == pytest.approx(arm.righting_arm, rel=0, abs=1e-12), f"heel {arm.heel_angle_deg}"
    # and the curve is genuinely non-trivial: it is not a sine, and it changes sign
    assert abs(exact.gz[3]) > 0.01
    assert exact.max_arm() is not None


def test_polygon_oracle_is_exact_across_a_waterline_vertex_crossing() -> None:
    """The waterline sweeping over a section vertex changes the clip topology; still exact.

    ``atan(1/2) = 26.565 deg`` puts the waterline through the V's starboard bilge corner, where
    the number of clipped edges jumps.  Sampling either side of it and *at* it checks that the
    topology change is not smoothed over by a tolerance.
    """
    crossing = math.degrees(math.atan2(1.0, 2.0))
    angles = [crossing - 1e-6, crossing, crossing + 1e-6]
    mesh = section_gz.prism_mesh(TRIANGLE, 3.0)
    exact = section_gz.polygon_prism_righting_arms(
        TRIANGLE, angles, kg=0.5, waterline_z=1.0, length=3.0
    )
    curve = righting_arm_curve(mesh, 0.5, angles, waterline_z=1.0)
    for arm, gz in zip(exact.arms, curve.gz, strict=True):
        assert gz == pytest.approx(arm.righting_arm, rel=0, abs=1e-11)
        assert arm.submerged_area == pytest.approx(exact.displaced_area, rel=1e-13)
    # continuity of the arm across the crossing, from the oracle alone
    assert abs(exact.gz[1] - 0.5 * (exact.gz[0] + exact.gz[2])) < 1e-6


def test_polygon_oracle_reproduces_the_wall_sided_box_and_then_stops() -> None:
    """Below the bilge limit it agrees with the box closed form; above it, it must not."""
    breadth, draft, kg = 4.0, 1.0, 0.6
    rectangle = (
        (-breadth / 2.0, 0.0),
        (breadth / 2.0, 0.0),
        (breadth / 2.0, 2.0),
        (-breadth / 2.0, 2.0),
    )
    limit = box_wall_sided_limit_deg(breadth, draft)
    inside = [0.0, 5.0, 10.0, 20.0, 0.9 * limit]
    outside = [limit + 5.0, limit + 15.0]
    mesh = section_gz.prism_mesh(rectangle, 5.0)
    bm0 = breadth**2 / (12.0 * draft)
    gm0 = draft / 2.0 + bm0 - kg
    analytic_arms = box_righting_arms_analytic(breadth, draft, kg, inside)
    for angle, analytic in zip(inside, analytic_arms, strict=True):
        oracle = section_gz.polygon_prism_righting_arms(
            rectangle, [angle], kg=kg, waterline_z=draft
        ).gz[0]
        assert oracle == pytest.approx(analytic, rel=0, abs=1e-12)
        assert oracle == pytest.approx(
            math.sin(math.radians(angle)) * (gm0 + 0.5 * bm0 * math.tan(math.radians(angle)) ** 2),
            rel=1e-12,
        )
        assert righting_arm_curve(mesh, kg, [angle], waterline_z=draft).gz[0] == pytest.approx(
            analytic, rel=0, abs=1e-12
        )
    for angle in outside:
        oracle = section_gz.polygon_prism_righting_arms(
            rectangle, [angle], kg=kg, waterline_z=draft
        ).gz[0]
        wall_sided = math.sin(math.radians(angle)) * (
            gm0 + 0.5 * bm0 * math.tan(math.radians(angle)) ** 2
        )
        assert abs(oracle - wall_sided) > 1e-3, "the comparison above would be vacuous"


def test_centrally_symmetric_section_at_half_draft_is_the_log_formula() -> None:
    """A hexagon at half submergence must reduce to ``(R - KG) sin(phi)`` on symmetry axes."""
    radius, kg = 1.0, 0.5
    section = section_gz.regular_polygon_section(6, radius)
    step = 60.0
    angles = [0.0, step, 2 * step, 3 * step - step]  # symmetry-aligned heels
    curve = section_gz.polygon_prism_righting_arms(
        section, angles, kg=kg, waterline_z=radius, length=2.0
    )
    for angle, gz in zip(angles, curve.gz, strict=True):
        assert gz == pytest.approx((radius - kg) * math.sin(math.radians(angle)), rel=0, abs=1e-13)


def test_rotated_section_is_exact_off_the_symmetry_axes() -> None:
    """A section turned off the symmetry axes is where a symmetry-only reference passes vacuously.

    The existing log test only samples heels that are symmetry axes of the polygon; here the
    polygon itself is rotated by 17 deg, so no sampled heel lines up and the mesh integrator is
    exercised on cut directions it has never been checked on.
    """
    section = section_gz.regular_polygon_section(6, 1.0, rotation_deg=17.0)
    angles = [0.0, 7.0, 23.0, 41.0, 68.0, 90.0]
    mesh = section_gz.prism_mesh(section, 2.5)
    assert meshlib.is_watertight(mesh)
    exact = section_gz.polygon_prism_righting_arms(
        section, angles, kg=0.7, waterline_z=0.8, length=2.5
    )
    curve = righting_arm_curve(mesh, 0.7, angles, waterline_z=0.8)
    for arm, gz in zip(exact.arms, curve.gz, strict=True):
        assert gz == pytest.approx(arm.righting_arm, rel=0, abs=1e-12), f"heel {arm.heel_angle_deg}"


def test_ellipse_oracle_slope_at_zero_equals_the_waterplane_GM() -> None:
    """``dGZ/dphi (0) = KB + I_T/A - KG`` with KB and ``I_T`` derived in the test, not copied.

    The oracle integrates heeled areas; this checks it against the *waterplane moment* route
    that :mod:`tensorcad.stability` takes, which uses none of the same code.
    """
    semi_beam, semi_height, draft, kg = 2.0, 1.0, 1.2, 0.8
    # submerged segment of the ellipse: affine image of the unit-disk segment at d = (draft-b)/b
    d = (draft - semi_height) / semi_height
    half_width = semi_beam * math.sqrt(1.0 - d * d)
    segment_area_disk = 0.5 * math.pi + math.asin(d) + d * math.sqrt(1.0 - d * d)
    area = semi_beam * semi_height * segment_area_disk
    moment_coefficient = -(2.0 / 3.0) * (1.0 - d * d) ** 1.5
    kb = semi_height + semi_height * moment_coefficient / segment_area_disk
    inertias = (2.0 / 3.0) * half_width**3  # per unit length, about the centreline
    gm = kb + inertias / area - kg

    delta = 1e-6
    curve = _ellipse_curve([0.0, delta], draft=draft, kg=kg)
    slope = (curve.gz[1] - curve.gz[0]) / math.radians(delta)
    assert curve.gz[0] == pytest.approx(0.0, abs=1e-15)
    assert slope == pytest.approx(gm, rel=1e-8)


def test_ellipse_oracle_converges_at_second_order_under_chord_refinement() -> None:
    """An inscribed polygon under-fills the ellipse; refining it must converge to the oracle.

    This is the oracle's own independent check: nothing but the chord error separates the mesh
    integrator on the polygon from the closed form, and a second-order rate is what an inscribed
    chord in a smooth boundary predicts.  A wrong closed form would not converge to it at all.
    """
    semi_beam, semi_height, draft, kg = 2.0, 1.0, 1.2, 0.8
    exact = _ellipse_curve(semi_beam=semi_beam, semi_height=semi_height, draft=draft, kg=kg)
    errors = []
    for sides in (16, 32, 64, 128):
        section = section_gz.ellipse_section(semi_beam, semi_height, sides)
        mesh = section_gz.prism_mesh(section, 4.0)
        curve = righting_arm_curve(mesh, kg, ANGLES, waterline_z=draft)
        errors.append(max(abs(a - b) for a, b in zip(curve.gz, exact.gz, strict=True)))
    for coarse, fine in zip(errors, errors[1:], strict=False):
        assert fine < coarse
        assert math.log2(coarse / fine) == pytest.approx(2.0, abs=0.25)
    assert errors[-1] < 1e-3  # 128 chords: under a millimetre on a 4 m arm
    assert errors[0] < 0.02


def test_ellipse_curve_is_not_a_sine_and_the_wall_sided_formula_is_not_enough() -> None:
    """The whole point of the new oracle: shape matters, and the sine answer is wrong.

    ``GZ / sin(phi)`` is constant for a circular section at every heel; for a 2:1 ellipse it is
    not, and the wall-sided formula (which would be exact for a box with the same ``GM`` and
    ``BM``) misses the exact arm by more than the 3 % budget at 30 deg.
    """
    semi_beam, semi_height, draft, kg = 2.0, 1.0, 1.2, 0.8
    curve = _ellipse_curve(semi_beam=semi_beam, semi_height=semi_height, draft=draft, kg=kg)
    ratios = [
        gz / math.sin(math.radians(angle))
        for angle, gz in zip(curve.angles_deg, curve.gz, strict=True)
        if angle > 0.0
    ]
    assert max(ratios) - min(ratios) > 0.1
    # wall-sided reference built from this section's own upright GM and BM
    d = (draft - semi_height) / semi_height
    segment_area_disk = 0.5 * math.pi + math.asin(d) + d * math.sqrt(1.0 - d * d)
    area = semi_beam * semi_height * segment_area_disk
    half_width = semi_beam * math.sqrt(1.0 - d * d)
    bm = (2.0 / 3.0) * half_width**3 / area
    moment_coefficient = -(2.0 / 3.0) * (1.0 - d * d) ** 1.5
    kb = semi_height + semi_height * moment_coefficient / segment_area_disk
    gm = kb + bm - kg
    angle = 30.0
    wall_sided = math.sin(math.radians(angle)) * (
        gm + 0.5 * bm * math.tan(math.radians(angle)) ** 2
    )
    exact_arm = curve.angle_of(angle).righting_arm
    assert abs(exact_arm - wall_sided) / abs(exact_arm) > 0.03


def test_vanishing_stability_angle_is_where_the_arm_changes_sign() -> None:
    """A high ``KG`` on the 2:1 ellipse gives a real crossing; interpolation must land on it."""
    semi_beam, semi_height, draft, kg = 2.0, 1.0, 1.2, 1.6
    angles = [float(angle) for angle in range(0, 91, 2)]
    curve = _ellipse_curve(angles, semi_beam=semi_beam, semi_height=semi_height, draft=draft, kg=kg)
    vanishing = curve.vanishing_stability_angle_deg()
    assert vanishing is not None
    assert 0.0 < vanishing < 90.0
    # the crossing is reported by linear interpolation on a 2 deg grid, so the arm there is
    # small rather than exactly zero: quote what that residual actually is
    at_crossing = _ellipse_curve(
        [vanishing], semi_beam=semi_beam, semi_height=semi_height, draft=draft, kg=kg
    ).gz[0]
    assert abs(at_crossing) < 1e-3
    assert abs(at_crossing) < 1e-2 * max(curve.gz)
    # and the crossing is bracketed by the sampled points it was interpolated from
    near = _ellipse_curve(
        [vanishing - 1.0, vanishing + 1.0],
        semi_beam=semi_beam,
        semi_height=semi_height,
        draft=draft,
        kg=kg,
    )
    assert near.gz[0] > 0.0 > near.gz[1]


def test_curve_quantities_are_consistent_with_their_definitions() -> None:
    """``max_arm``, the dynamic-stability area and the reported area are self-consistent."""
    curve = _ellipse_curve()
    best = curve.max_arm()
    assert best is not None
    assert best[1] == 45.0  # the sampled peak of this 2:1 section
    assert best[0] == pytest.approx(max(curve.gz))
    # dynamic stability up to 30 deg, with the *actual* (non-uniform) sampled angles, so this
    # checks the trapezoid weights and the degree-to-radian conversion and nothing else
    area_limit = curve.area_under_curve(30.0)
    manual = 0.0
    for left, right in zip(curve.arms, curve.arms[1:], strict=False):
        if right.heel_angle_deg > 30.0:
            break
        manual += (
            0.5
            * (left.righting_arm + right.righting_arm)
            * math.radians(right.heel_angle_deg - left.heel_angle_deg)
        )
    assert area_limit == pytest.approx(manual, rel=1e-14)
    assert area_limit < curve.area_under_curve()
    assert curve.to_dict()["schema"] == section_gz.SCHEMA_SECTION_CURVE
    for arm in curve.arms:
        assert arm.submerged_area == pytest.approx(curve.displaced_area, rel=1e-12)


def test_the_two_oracles_agree_on_the_continuum_limit() -> None:
    """Polygon oracle on a fine ellipse polygon must meet the exact ellipse oracle.

    Two implementations that share no code path -- 2D clipping of a polygon versus the affine
    disk segment -- landing on the same arm is the strongest evidence either is right.
    """
    semi_beam, semi_height, draft, kg = 2.0, 1.0, 1.2, 0.8
    exact = _ellipse_curve(semi_beam=semi_beam, semi_height=semi_height, draft=draft, kg=kg)
    polygon = section_gz.polygon_prism_righting_arms(
        section_gz.ellipse_section(semi_beam, semi_height, 512),
        ANGLES,
        kg=kg,
        waterline_z=draft,
        length=4.0,
    )
    for angle, left, right in zip(ANGLES, exact.gz, polygon.gz, strict=True):
        assert left == pytest.approx(right, rel=0, abs=2e-4), f"heel {angle}"


def test_negative_heel_is_antisymmetric_for_a_symmetric_section() -> None:
    """Port and starboard arms are equal and opposite, in the oracle as in the mesh."""
    starboard = section_gz.polygon_prism_righting_arms(
        TRIANGLE, [10.0, 30.0, 60.0], kg=0.6, waterline_z=1.0, length=4.0
    )
    port = section_gz.polygon_prism_righting_arms(
        TRIANGLE, [-10.0, -30.0, -60.0], kg=0.6, waterline_z=1.0, length=4.0
    )
    for left, right in zip(starboard.gz, port.gz, strict=True):
        assert left == pytest.approx(-right, rel=0, abs=1e-13)


def test_concave_section_is_refused() -> None:
    """A concave section can float at two levels; the oracle refuses rather than guessing."""
    arrow = ((0.0, 0.0), (2.0, 0.0), (1.0, 0.5), (2.0, 2.0), (0.0, 2.0))
    with pytest.raises(section_gz.GzError, match="not convex"):
        section_gz.polygon_prism_righting_arms(arrow, [0.0], kg=0.5, waterline_z=1.0)
    with pytest.raises(section_gz.GzError, match="not convex"):
        section_gz.prism_mesh(arrow, 1.0)


def test_degenerate_inputs_are_refused() -> None:
    """Waterline outside the section, a zero-area section and a bad ellipse all raise."""
    square = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
    with pytest.raises(section_gz.GzError, match="below the section"):
        section_gz.polygon_prism_righting_arms(square, [0.0], kg=0.2, waterline_z=1.5)
    with pytest.raises(section_gz.GzError, match="above the keel"):
        section_gz.polygon_prism_righting_arms(square, [0.0], kg=0.2, waterline_z=0.0)
    with pytest.raises(section_gz.GzError, match="zero area"):
        section_gz.polygon_prism_righting_arms(((0.0, 0.0), (1.0, 0.0), (2.0, 0.0)), [0.0], kg=0.1, waterline_z=0.1)
    with pytest.raises(section_gz.GzError, match="repeat"):
        section_gz.polygon_prism_righting_arms(
            ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (1.0, 0.0)), [0.0], kg=0.1, waterline_z=0.5
        )
    with pytest.raises(section_gz.GzError, match="semi-axes"):
        section_gz.ellipse_prism_righting_arms(0.0, 1.0, [0.0], kg=0.5, waterline_z=0.5)
    with pytest.raises(section_gz.GzError, match="strictly inside the section height"):
        section_gz.ellipse_prism_righting_arms(1.0, 1.0, [0.0], kg=0.5, waterline_z=2.5)
    with pytest.raises(section_gz.GzError, match="exceeds the section area"):
        section_gz.ellipse_prism_righting_arms(
            1.0, 1.0, [0.0], kg=0.5, waterline_z=1.0, displaced_area=4.0
        )
    with pytest.raises(section_gz.GzError, match="within"):
        section_gz.polygon_prism_righting_arms(square, [190.0], kg=0.2, waterline_z=0.5)
    with pytest.raises(section_gz.GzError, match="prism length"):
        section_gz.prism_mesh(square, 0.0)


def test_mesh_kwarg_matches_the_oracle_frame() -> None:
    """A section drawn with its keel above ``z = 0`` is shifted, not rescaled."""
    raised = tuple((y, z + 3.0) for (y, z) in TRIANGLE)
    flat = section_gz.polygon_prism_righting_arms(TRIANGLE, ANGLES, kg=0.6, waterline_z=1.0)
    lifted = section_gz.polygon_prism_righting_arms(raised, ANGLES, kg=0.6, waterline_z=4.0)
    for left, right in zip(flat.gz, lifted.gz, strict=True):
        assert left == pytest.approx(right, rel=0, abs=1e-12)


def test_mesh_path_is_available_for_the_same_section() -> None:
    """The oracle and the mesh path stay usable together (import-level regression)."""
    mesh: Mesh = section_gz.prism_mesh(TRIANGLE, 2.0)
    curve = righting_arm_curve(mesh, 0.6, [0.0, 30.0], waterline_z=1.0)
    assert curve.max_volume_error < 1e-12
    assert isinstance(curve.gz, tuple)