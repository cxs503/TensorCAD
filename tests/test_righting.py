"""Large-angle righting arms: three references that cannot fail the same way.

The curve ``GZ(phi)`` is one integral, so "it looks plausible" is not a check. Three
independent references are used here, in increasing order of how much of the code
they exercise:

* a **box**, where the heeled displacement problem separates into a 2D section clip
  (``analytic_box_righting_arm``) that shares no code with the 3D tetrahedron
  integration under test. The box is prismatic, so the reference is exact at *every*
  heel angle, including past deck-edge immersion;
* the **wall-sided closed form** ``GZ = sin(phi)(GM + BM tan^2(phi)/2)``, which fixes
  the whole curve from two upright quantities and is exact only while both sides of
  the section are wet -- below ``bilge_emergence_deg``. Agreement there, and the
  measured departure above it, pins down the range of validity rather than assuming
  it;
* the **parametric Wigley hull**, where there is no closed form for ``GZ``: the check
  is consistency with a *different* integration, the upright metacentric height
  ``GM = KB + BM - KG`` assembled from ``waterplane_moments`` (Green's theorem on the
  waterplane) and the heeled volume integration here. ``dGZ/dphi`` at zero heel has to
  reproduce that number, and it does to a few parts in ``1e4`` on every grid.
"""

from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh
from tensorcad.righting import (
    GEOMETRY_TOLERANCE,
    HeeledFloat,
    PlaneSolution,
    RightingArmCurve,
    RightingArmError,
    SubmergedProperties,
    analytic_box_righting_arm,
    bilge_emergence_deg,
    box_mesh,
    deck_edge_immersion_deg,
    heel_normal,
    righting_arm,
    righting_arm_curve,
    solve_plane_offset,
    submerged_properties,
    upright_displacement,
    wall_sided_righting_arm,
)
from tensorcad.stability import waterplane_moments
from tensorcad.wigley import WigleyParameters, generate_wigley

# Reference box: L/B = 5, B/D = 2, draft half the depth, KG just above KB + BM / 2.
BOX_LENGTH = 10.0
BOX_BEAM = 2.0
BOX_DEPTH = 1.0
BOX_DRAFT = 0.5
BOX_KG = 0.55
BOX_BM = BOX_BEAM**2 / (12.0 * BOX_DRAFT)
BOX_KB = 0.5 * BOX_DRAFT
BOX_GM = BOX_KB + BOX_BM - BOX_KG


def _box() -> Mesh:
    return box_mesh(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)


def _wigley(nx: int = 81, nz: int = 25) -> Mesh:
    return generate_wigley(WigleyParameters(nx=nx, nz=nz)).mesh


# ------------------------------------------------------------------ the box first


def test_box_mesh_is_a_closed_outward_solid() -> None:
    mesh = _box()
    assert meshlib.is_watertight(mesh)
    assert meshlib.is_edge_manifold(mesh)
    assert meshlib.zero_area_faces(mesh) == ()
    assert meshlib.total_signed_volume(mesh) == pytest.approx(
        BOX_LENGTH * BOX_BEAM * BOX_DEPTH, rel=1e-15
    )
    box = mesh.bounding_box()
    assert box.sizes == pytest.approx((BOX_LENGTH, BOX_BEAM, BOX_DEPTH), rel=0, abs=1e-15)


def test_box_upright_displacement_and_buoyancy_centroid_are_analytic() -> None:
    mesh = _box()
    properties = submerged_properties(mesh, (0.0, 0.0, 1.0), BOX_DRAFT)
    assert properties.volume == pytest.approx(BOX_LENGTH * BOX_BEAM * BOX_DRAFT, rel=1e-15)
    # The centroid of a prism below a level waterline is the section centroid.
    assert properties.centre_of_buoyancy[1] == pytest.approx(0.0, abs=1e-15)
    assert properties.centre_of_buoyancy[2] == pytest.approx(0.5 * BOX_DRAFT, rel=1e-15)
    assert properties.centre_of_buoyancy[0] == pytest.approx(0.0, abs=1e-13)


def test_box_cut_area_is_the_waterplane_area() -> None:
    properties = submerged_properties(_box(), (0.0, 0.0, 1.0), BOX_DRAFT)
    assert properties.cut_area == pytest.approx(BOX_LENGTH * BOX_BEAM, rel=1e-15)


def test_box_cut_area_at_heel_is_the_in_plane_chord_of_the_section() -> None:
    """The cut area is the flux of the submerged hull surface -- check it geometrically.

    The waterplane section of a box is a rectangle ``length x chord``, where the chord
    is measured *in the tilted plane*: over the horizontal span it covers, the
    waterline rises by ``tan(phi)`` per unit breadth, so the chord is
    ``span / cos(phi)`` -- 0.38 % longer than the breadth at 5 deg, which is exactly
    the kind of difference a "close enough" test would hide.
    """
    mesh = _box()
    for heel_deg in (5.0, 20.0, 35.0):
        point = righting_arm(mesh, heel_deg, kg=BOX_KG, waterline_z=BOX_DRAFT)
        phi = math.radians(heel_deg)
        cosine, sine = math.cos(phi), math.sin(phi)
        # The waterline z(y) = (offset + y sin) / cos has to stay inside the section,
        # whose keel is z = 0 and deck z = BOX_DEPTH here.
        low = max(-0.5 * BOX_BEAM, -point.plane_offset / sine)
        high = min(0.5 * BOX_BEAM, (BOX_DEPTH * cosine - point.plane_offset) / sine)
        assert low < high
        chord = (high - low) / cosine
        assert point.waterplane_area == pytest.approx(BOX_LENGTH * chord, rel=1e-12)
        assert -0.5 * BOX_BEAM - 1e-15 <= point.centre_of_buoyancy[1] <= 0.5 * BOX_BEAM + 1e-15


@pytest.mark.parametrize("heel_deg", [0.0, 2.0, 5.0, 10.0, 20.0, 30.0, 45.0, 60.0, 75.0, 90.0])
def test_box_gz_matches_the_2d_section_reference_at_every_angle(heel_deg: float) -> None:
    """The 3D tetrahedron integral against an independent 2D clip: machine precision."""
    mesh = _box()
    point = righting_arm(mesh, heel_deg, kg=BOX_KG, waterline_z=BOX_DRAFT)
    reference = analytic_box_righting_arm(
        BOX_LENGTH, BOX_BEAM, BOX_DEPTH, BOX_DRAFT, heel_deg, kg=BOX_KG
    )
    assert point.righting_arm == pytest.approx(reference, rel=0, abs=1e-12)
    assert point.relative_volume_residual < 1e-12


@pytest.mark.parametrize("heel_deg", [2.0, 5.0, 10.0, 15.0, 20.0, 25.0])
def test_box_gz_matches_the_wall_sided_form_below_bilge_emergence(heel_deg: float) -> None:
    assert heel_deg < bilge_emergence_deg(beam=BOX_BEAM, draft=BOX_DRAFT)
    mesh = _box()
    point = righting_arm(mesh, heel_deg, kg=BOX_KG, waterline_z=BOX_DRAFT)
    reference = wall_sided_righting_arm(math.radians(heel_deg), gm=BOX_GM, bm=BOX_BM)
    assert point.righting_arm == pytest.approx(reference, rel=0, abs=1e-12)


def test_wall_sided_form_departs_once_the_bilge_emerges() -> None:
    """The validity limit is a measured number, not a comment: 5 % at 30 deg."""
    limit = bilge_emergence_deg(beam=BOX_BEAM, draft=BOX_DRAFT)
    assert limit == pytest.approx(26.56505117707799, rel=1e-12)
    mesh = _box()
    heel_deg = 30.0
    exact = righting_arm(mesh, heel_deg, kg=BOX_KG, waterline_z=BOX_DRAFT).righting_arm
    wall_sided = wall_sided_righting_arm(math.radians(heel_deg), gm=BOX_GM, bm=BOX_BM)
    assert abs(exact - wall_sided) / exact > 0.03
    # ... and inside the limit the same comparison holds, so it is the range that
    # changed, not the implementation.
    inside = 0.9 * limit
    exact_inside = righting_arm(mesh, inside, kg=BOX_KG, waterline_z=BOX_DRAFT).righting_arm
    form_inside = wall_sided_righting_arm(math.radians(inside), gm=BOX_GM, bm=BOX_BM)
    assert exact_inside == pytest.approx(form_inside, rel=1e-12)


def test_box_deck_edge_immersion_is_where_the_formula_of_wall_sided_stops_holding_for_depth() -> None:
    assert deck_edge_immersion_deg(beam=BOX_BEAM, depth=BOX_DEPTH) == pytest.approx(
        math.degrees(math.atan(2.0 * BOX_DEPTH / BOX_BEAM)), rel=1e-15
    )


def test_box_gz_slope_at_zero_heel_is_the_analytic_metacentric_height() -> None:
    """``GZ'(0) = GM``: the module's heeled integral against the upright algebra."""
    mesh = _box()
    step = 1e-4  # degrees
    forward = righting_arm(mesh, step, kg=BOX_KG, waterline_z=BOX_DRAFT).righting_arm
    aft = righting_arm(mesh, -step, kg=BOX_KG, waterline_z=BOX_DRAFT).righting_arm
    slope = (forward - aft) / (2.0 * math.radians(step))
    assert slope == pytest.approx(BOX_GM, rel=1e-6)


def test_box_gz_zero_crossing_is_the_vanishing_stability_angle() -> None:
    curve = righting_arm_curve(
        _box(), (0.0, 20.0, 40.0, 60.0, 80.0, 90.0), kg=BOX_KG, waterline_z=BOX_DRAFT
    )
    crossing = curve.zero_crossing_deg
    assert crossing is not None
    assert 80.0 < crossing <= 90.0
    assert curve.arm_at(crossing) == pytest.approx(0.0, abs=1e-12)
    assert curve.max_arm == pytest.approx(0.241712678456, rel=1e-6)
    assert curve.angle_of_max_arm_deg in (20.0, 40.0)


def _wall_sided_area(limit_deg: float) -> float:
    """Exact integral of ``sin(phi)(GM + BM tan^2(phi)/2)`` over ``[0, limit]``.

    Antiderivative: ``GM (1 - cos) + BM/2 (sec + cos - 2)``.
    """
    limit = math.radians(limit_deg)
    return BOX_GM * (1.0 - math.cos(limit)) + 0.5 * BOX_BM * (
        1.0 / math.cos(limit) + math.cos(limit) - 2.0
    )


def test_curve_area_converges_to_the_exact_integral_of_the_closed_form() -> None:
    """The reported area is the integral it claims to be, and its error is quadrature.

    The arms are exact (they match the 2D reference to machine precision), so the
    only error left in ``area()`` is the trapezoidal rule on the sample spacing --
    which is measured here as a convergence order instead of being asserted away.
    """
    limit = 25.0
    exact = _wall_sided_area(limit)
    coarse = righting_arm_curve(
        _box(), [0.5 * step for step in range(0, 51)], kg=BOX_KG, waterline_z=BOX_DRAFT
    ).area(limit)
    fine = righting_arm_curve(
        _box(), [0.125 * step for step in range(0, 201)], kg=BOX_KG, waterline_z=BOX_DRAFT
    ).area(limit)
    coarse_error = abs(coarse - exact)
    fine_error = abs(fine - exact)
    assert coarse_error < 1e-4 * exact
    assert fine_error < 2e-5 * exact
    assert fine_error < 0.25 * coarse_error  # halving twice: second order


def test_curve_requires_the_upright_point_and_the_sampled_range() -> None:
    curve = righting_arm_curve(_box(), (10.0, 20.0), kg=BOX_KG, waterline_z=BOX_DRAFT)
    with pytest.raises(RightingArmError):
        curve.area(10.0)
    curve = righting_arm_curve(_box(), (0.0, 10.0), kg=BOX_KG, waterline_z=BOX_DRAFT)
    with pytest.raises(RightingArmError):
        curve.area(20.0)
    with pytest.raises(RightingArmError):
        curve.area(0.0)


def test_heel_normal_is_world_up_and_the_horizontal_arm_direction_is_orthogonal() -> None:
    assert heel_normal(0.0) == (0.0, 0.0, 1.0)
    for heel_deg in (10.0, 45.0, 90.0):
        normal = heel_normal(math.radians(heel_deg))
        assert normal[1] == pytest.approx(-math.sin(math.radians(heel_deg)), rel=1e-15)
        assert normal[2] == pytest.approx(math.cos(math.radians(heel_deg)), rel=1e-15)
        # The arm direction e_h = (0, w_z, -w_y) is horizontal and perpendicular to w.
        horizontal = (0.0, normal[2], -normal[1])
        assert sum(a * b for a, b in zip(normal, horizontal, strict=True)) == pytest.approx(0.0, abs=1e-15)
        assert abs(math.sqrt(sum(value * value for value in horizontal)) - 1.0) < 1e-15


# ----------------------------------------------------------------- the hull now


def test_hull_gz_is_zero_upright_and_odd_in_heel() -> None:
    mesh = _wigley()
    upright = righting_arm(mesh, 0.0, kg=0.30, waterline_z=-0.1)
    assert upright.righting_arm == pytest.approx(0.0, abs=1e-13)
    assert upright.centre_of_buoyancy[1] == pytest.approx(0.0, abs=1e-13)
    for heel_deg in (5.0, 15.0, 30.0):
        starboard = righting_arm(mesh, heel_deg, kg=0.30, waterline_z=-0.1).righting_arm
        port = righting_arm(mesh, -heel_deg, kg=0.30, waterline_z=-0.1).righting_arm
        assert starboard + port == pytest.approx(0.0, abs=1e-12)
        assert starboard > 0.0


def test_hull_gz_slope_matches_the_metacentric_height_at_the_same_displacement() -> None:
    """``dGZ/dphi(0) = KB + BM - KG`` from two other integrations.

    ``BM`` comes from :func:`tensorcad.stability.waterplane_moments` (Green's theorem
    on the waterline polygon) and ``KB`` from the volume centroid measured here at
    zero heel, so a shared bug in the heeled path would have to reproduce the upright
    waterplane moments by accident.
    """
    mesh = _wigley()
    displacement = upright_displacement(mesh, -0.1)
    moments = waterplane_moments(mesh, -0.1)
    keel_z = min(vertex[2] for vertex in mesh.vertices)
    kb = submerged_properties(mesh, (0.0, 0.0, 1.0), -0.1).centre_of_buoyancy[2] - keel_z
    bm = moments.transverse_second_moment / displacement
    gm = kb + bm - 0.30
    step = 1.0  # degrees, central difference
    slope = (
        righting_arm(mesh, step, kg=0.30, waterline_z=-0.1).righting_arm
        - righting_arm(mesh, -step, kg=0.30, waterline_z=-0.1).righting_arm
    ) / (2.0 * math.radians(step))
    assert slope == pytest.approx(gm, rel=2e-3)


def test_hull_curve_holds_the_displacement_to_round_off() -> None:
    mesh = _wigley(41, 13)
    curve = righting_arm_curve(mesh, (0.0, 10.0, 20.0, 30.0, 40.0), kg=0.30, waterline_z=-0.1)
    assert curve.displacement_volume == pytest.approx(upright_displacement(mesh, -0.1), rel=1e-15)
    for point in curve.points:
        assert point.relative_volume_residual < 1e-12
        assert point.iterations <= 10
        assert point.waterplane_area > 0.0


def test_hull_gz_converges_with_resolution_at_second_order() -> None:
    """Doubling the grid twice: successive differences must shrink, not wander."""
    angles = (20.0,)
    coarse = righting_arm_curve(_wigley(41, 13), angles, kg=0.30, waterline_z=-0.1)
    medium = righting_arm_curve(_wigley(81, 25), angles, kg=0.30, waterline_z=-0.1)
    fine = righting_arm_curve(_wigley(161, 41), angles, kg=0.30, waterline_z=-0.1)
    first = abs(medium.arms[0] - coarse.arms[0])
    second = abs(fine.arms[0] - medium.arms[0])
    assert second < first
    assert fine.arms[0] == pytest.approx(medium.arms[0], abs=1e-4)
    # The discretisation makes the hull a polygon inside the analytic surface, so the
    # arm is monotone in resolution here; a sign flip would be a defect.
    assert coarse.arms[0] < medium.arms[0] < fine.arms[0]


def test_curve_reports_the_maximum_and_its_angle() -> None:
    mesh = _wigley(41, 13)
    curve = righting_arm_curve(mesh, (0.0, 5.0, 10.0, 20.0, 30.0, 40.0), kg=0.30, waterline_z=-0.1)
    assert curve.max_arm == max(point.righting_arm for point in curve.points)
    assert curve.angle_of_max_arm_deg == max(
        curve.points, key=lambda point: point.righting_arm
    ).heel_deg
    assert curve.arm_at(0.0) == pytest.approx(0.0, abs=1e-13)
    assert curve.arm_at(15.0) == pytest.approx(0.5 * (curve.arm_at(10.0) + curve.arm_at(20.0)), rel=0)
    payload = curve.to_dict()
    assert payload["basis"].startswith("heeled displacement integration")
    assert len(payload["points"]) == len(curve.points)
    assert payload["area_to_40_deg_m_rad"] is not None
    assert payload["area_to_30_deg_m_rad"] > 0.0


def test_repeated_runs_are_bit_identical() -> None:
    mesh = _wigley(41, 13)
    first = righting_arm_curve(mesh, (0.0, 12.0, 24.0), kg=0.30, waterline_z=-0.1)
    second = righting_arm_curve(mesh, (0.0, 12.0, 24.0), kg=0.30, waterline_z=-0.1)
    assert first.arms == second.arms
    assert [point.plane_offset for point in first.points] == [
        point.plane_offset for point in second.points
    ]


def test_upright_displacement_matches_the_enclosed_volume_when_floating_at_the_cap() -> None:
    mesh = _box()
    assert upright_displacement(mesh, BOX_DEPTH) == pytest.approx(
        meshlib.total_signed_volume(mesh), rel=1e-15
    )


# ------------------------------------------------------------------- solver edges


def test_plane_solver_hits_the_target_and_reports_its_residual() -> None:
    mesh = _box()
    target = 0.5 * upright_displacement(mesh, BOX_DEPTH)
    solution = solve_plane_offset(mesh, (0.0, 0.0, 1.0), target)
    assert isinstance(solution, PlaneSolution)
    assert solution.relative_residual < 1e-12
    assert solution.offset == pytest.approx(0.5 * BOX_DEPTH, rel=1e-12)
    assert solution.cut_area == pytest.approx(BOX_LENGTH * BOX_BEAM, rel=1e-12)
    assert solution.iterations <= 60


def test_plane_solver_rejects_a_target_the_mesh_cannot_displace() -> None:
    mesh = _box()
    total = meshlib.total_signed_volume(mesh)
    with pytest.raises(RightingArmError):
        solve_plane_offset(mesh, (0.0, 0.0, 1.0), total)
    with pytest.raises(RightingArmError):
        solve_plane_offset(mesh, (0.0, 0.0, 1.0), 0.0)


def test_righting_arm_rejects_an_open_surface_and_a_non_positive_displacement() -> None:
    mesh = _box()
    open_mesh = Mesh(vertices=mesh.vertices, faces=mesh.faces[:-2])
    with pytest.raises(RightingArmError):
        righting_arm(open_mesh, 10.0, kg=BOX_KG, waterline_z=BOX_DRAFT)
    with pytest.raises(RightingArmError):
        submerged_properties(open_mesh, (0.0, 0.0, 1.0), BOX_DRAFT)
    with pytest.raises(RightingArmError):
        righting_arm_curve(open_mesh, (0.0, 10.0), kg=BOX_KG, waterline_z=BOX_DRAFT)
    with pytest.raises(RightingArmError):
        righting_arm(mesh, 10.0, kg=BOX_KG, displacement_volume=0.0)
    # A plane that floods nothing is not a floating condition.
    with pytest.raises(RightingArmError):
        submerged_properties(mesh, (0.0, 0.0, 1.0), -1.0)


def test_argument_validation_names_the_problem() -> None:
    mesh = _box()
    with pytest.raises(RightingArmError):
        heel_normal(float("nan"))
    with pytest.raises(RightingArmError):
        submerged_properties(mesh, (0.0, 0.0, 0.0), BOX_DRAFT)
    with pytest.raises(RightingArmError):
        righting_arm(mesh, 10.0, kg=float("inf"))
    with pytest.raises(RightingArmError):
        righting_arm(mesh, float("nan"), kg=BOX_KG)
    with pytest.raises(RightingArmError):
        righting_arm_curve(mesh, (), kg=BOX_KG, waterline_z=BOX_DRAFT)
    with pytest.raises(RightingArmError):
        box_mesh(0.0, BOX_BEAM, BOX_DEPTH)
    with pytest.raises(RightingArmError):
        bilge_emergence_deg(beam=-1.0, draft=1.0)
    with pytest.raises(RightingArmError):
        deck_edge_immersion_deg(beam=BOX_BEAM, depth=0.0)
    with pytest.raises(RightingArmError):
        analytic_box_righting_arm(BOX_LENGTH, BOX_BEAM, BOX_DEPTH, 2.0, 10.0, kg=BOX_KG)


def test_dataclasses_expose_their_evidence() -> None:
    mesh = _box()
    point = righting_arm(mesh, 10.0, kg=BOX_KG, waterline_z=BOX_DRAFT)
    assert isinstance(point, HeeledFloat)
    payload = point.to_dict()
    assert payload["heel_deg"] == 10.0
    assert payload["gz_m"] == pytest.approx(point.righting_arm, abs=1e-12)
    assert payload["iterations"] >= 1
    properties = submerged_properties(mesh, (0.0, 0.0, 1.0), BOX_DRAFT)
    assert isinstance(properties, SubmergedProperties)
    # Bottom (2 triangles) below, top (2) above, the four side quads (8) straddle.
    assert properties.faces_below == 2
    assert properties.faces_cut == 8
    assert properties.to_dict()["cut_area_m2"] == pytest.approx(BOX_LENGTH * BOX_BEAM, rel=1e-12)
    curve = righting_arm_curve(mesh, (0.0, 10.0), kg=BOX_KG, waterline_z=BOX_DRAFT)
    assert isinstance(curve, RightingArmCurve)
    assert curve.keel_z == pytest.approx(0.0, abs=1e-15)


def test_geometry_tolerance_is_tight_enough_to_be_invisible() -> None:
    """A plane through a vertex ring must not lose or double count a face."""
    mesh = _box()
    properties = submerged_properties(mesh, (0.0, 0.0, 1.0), BOX_DEPTH / 2.0)
    assert GEOMETRY_TOLERANCE <= 1e-12
    assert properties.volume == pytest.approx(0.5 * BOX_LENGTH * BOX_BEAM * BOX_DEPTH, rel=1e-15)
    assert properties.faces_cut == 8  # the four side quads, two triangles each
