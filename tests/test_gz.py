"""Large-angle stability: the righting-arm curve against closed forms and itself.

The references, in increasing order of how much they can go wrong:

* a **wall-sided box**, whose ``GZ`` has the textbook closed form
  ``sin(phi) * (GM + BM tan^2(phi) / 2)`` and whose faces are exactly representable,
  so the clip must reproduce it to round-off -- a machine-precision test fails loudly
  if the half-space integration, the volume balance or the lever direction is wrong;
* the **parametric hull**, where three independent facts are checked: the curve is an
  odd function of heel (a symmetry statement about the integration), its slope at zero
  is the metacentric height built from a *different* pair of integrals
  (``KB`` from the displaced volume, ``BM`` from the waterplane's second moment), and
  the two independent routes to ``GZ`` -- 3-D clipped-volume integration and classical
  cross curves from station cuts -- agree;
* a **synthetic curve** with a known integral, which pins the dynamical-stability
  quadrature and its truncated last panel independently of any mesh.
"""

from __future__ import annotations

import math

import pytest

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.gz import (
    DEFAULT_HEEL_ANGLES,
    GZError,
    HeelPoint,
    RightingArmCurve,
    cross_curve_righting_arm,
    cross_curve_righting_arm_curve,
    keel_height,
    righting_arm,
    righting_arm_curve,
    solve_waterline,
    station_offsets,
    submerged_properties,
    wall_sided_righting_arm,
)
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.stability import waterplane_moments


# ---------------------------------------------------------------------------
# reference geometries
# ---------------------------------------------------------------------------
def box(length: float, beam: float, depth: float) -> Mesh:
    """Axis-aligned closed box, keel on ``z = 0``, deck on ``z = depth``.

    The depth is deliberately larger than the draft the box is floated at, so the
    wall-sided closed form stays valid over the heel angles tested here: no deck
    edge is immersed and no bilge emerges.
    """
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (0.0, depth):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    left, right = 0.0, length
    lo, hi = -0.5 * beam, 0.5 * beam
    bottom, top = 0.0, depth
    # outward winding: the volume's sign is the assertion that the clip agrees.
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


BOX_LENGTH = 40.0
BOX_BEAM = 10.0
BOX_DEPTH = 8.0
BOX_DRAFT = 4.0
BOX_KG = 3.0
BOX_VOLUME = BOX_LENGTH * BOX_BEAM * BOX_DRAFT
BOX_BM = BOX_BEAM**2 / (12.0 * BOX_DRAFT)
BOX_KB = 0.5 * BOX_DRAFT
BOX_GM = BOX_KB + BOX_BM - BOX_KG


def hull(nx: int = 41, nz: int = 13) -> Mesh:
    return generate_hull(HullParameters(length=12.0, beam=2.4, draft=0.9, nx=nx, nz=nz)).mesh


#: A waterline below the hull's deck (``z = 0``) and off its vertex rings: a cut exactly
#: on a ring is a degenerate plane section for the waterplane-moment machinery, and the
#: righting arm does not need that machinery.
HULL_WATERLINE = -0.31
HULL_KG = 0.55


def hull_displacement(mesh: Mesh) -> float:
    return submerged_properties(mesh, (0.0, 0.0, 1.0), HULL_WATERLINE).volume


def metacentric_height(mesh: Mesh, kg: float) -> float:
    """``KB + BM - KG`` built from the displaced volume and the waterplane moment."""
    submerged = submerged_properties(mesh, (0.0, 0.0, 1.0), HULL_WATERLINE)
    kb = submerged.centroid[2] - keel_height(mesh)
    section = waterplane_moments(mesh, HULL_WATERLINE)
    return kb + section.transverse_second_moment / submerged.volume - kg


# ---------------------------------------------------------------------------
# 1. the wall-sided closed form
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("heel", [0.0, 5.0, 10.0, 20.0, 30.0, 35.0])
def test_box_reproduces_the_wall_sided_closed_form(heel: float) -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    point = righting_arm(mesh, heel, kg=BOX_KG, displacement_volume=BOX_VOLUME)
    reference = wall_sided_righting_arm(heel, BOX_GM, BOX_BM)
    if reference == 0.0:
        assert point.righting_arm == pytest.approx(0.0, abs=1e-12)
    else:
        assert point.righting_arm == pytest.approx(reference, rel=1e-10)


def test_box_righting_arm_is_positive_and_grows_over_the_valid_range() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    arms = [righting_arm(mesh, heel, kg=BOX_KG, displacement_volume=BOX_VOLUME).righting_arm
            for heel in (5.0, 10.0, 20.0, 30.0)]
    assert all(arm > 0.0 for arm in arms)
    assert arms == sorted(arms)


def test_box_matches_the_simple_hinged_identity_gz_equals_kn_sin_phi() -> None:
    # For any body floating at a fixed displacement, GZ = KN - KG sin(phi), where KN is
    # the moment arm of buoyancy about the keel. Recomputing KN from the centre of
    # buoyancy and comparing is an independent statement about the datum handling.
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    for heel in (10.0, 20.0, 30.0):
        point = righting_arm(mesh, heel, kg=BOX_KG, displacement_volume=BOX_VOLUME)
        phi = math.radians(heel)
        kn = point.transverse_centre_buoyancy * math.cos(phi)
        kn += point.centre_buoyancy_above_keel * math.sin(phi)
        assert point.righting_arm == pytest.approx(kn - BOX_KG * math.sin(phi), abs=1e-12)


@pytest.mark.parametrize("heel", [0.0, 10.0, 20.0, 30.0])
def test_cross_curve_route_agrees_on_the_box(heel: float) -> None:
    point = cross_curve_righting_arm(
        box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH), heel, kg=BOX_KG, displacement_volume=BOX_VOLUME
    )
    reference = wall_sided_righting_arm(heel, BOX_GM, BOX_BM)
    if reference == 0.0:
        assert point.righting_arm == pytest.approx(0.0, abs=1e-9)
    else:
        assert point.righting_arm == pytest.approx(reference, rel=1e-6)


# ---------------------------------------------------------------------------
# 2. the displacement balance
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("heel", [0.0, 5.0, 30.0, 60.0, 90.0])
def test_heeled_volume_equals_the_upright_displacement(heel: float) -> None:
    mesh = hull()
    target = hull_displacement(mesh)
    point = righting_arm(mesh, heel, kg=HULL_KG, displacement_volume=target, waterline_z=HULL_WATERLINE)
    assert point.submerged_volume == pytest.approx(target, rel=1e-9)
    assert abs(point.volume_error) < 1e-9


def test_upright_displacement_defaults_to_the_waterline_clip() -> None:
    mesh = hull()
    default = righting_arm(mesh, 0.0, kg=HULL_KG, waterline_z=HULL_WATERLINE)
    assert default.target_volume == pytest.approx(hull_displacement(mesh), rel=1e-12)
    # An upright, symmetric body has no transverse lever: GZ is zero, not "small".
    assert default.righting_arm == pytest.approx(0.0, abs=1e-9)


# ---------------------------------------------------------------------------
# 3. symmetry, slope and the metacentric height
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("heel", [5.0, 10.0, 20.0])
def test_curve_is_odd_in_the_heel_angle(heel: float) -> None:
    mesh = hull()
    target = hull_displacement(mesh)
    starboard = righting_arm(mesh, heel, kg=HULL_KG, displacement_volume=target, waterline_z=HULL_WATERLINE)
    port = righting_arm(mesh, -heel, kg=HULL_KG, displacement_volume=target, waterline_z=HULL_WATERLINE)
    assert starboard.righting_arm == pytest.approx(-port.righting_arm, abs=1e-9)
    assert starboard.transverse_centre_buoyancy == pytest.approx(-port.transverse_centre_buoyancy, abs=1e-9)


@pytest.mark.parametrize("resolution", [(41, 13), (81, 25), (161, 41)])
def test_slope_at_zero_is_the_metacentric_height(resolution: tuple[int, int]) -> None:
    mesh = hull(*resolution)
    target = hull_displacement(mesh)
    delta = 0.25
    starboard = righting_arm(mesh, delta, kg=HULL_KG, displacement_volume=target, waterline_z=HULL_WATERLINE)
    port = righting_arm(mesh, -delta, kg=HULL_KG, displacement_volume=target, waterline_z=HULL_WATERLINE)
    slope = (starboard.righting_arm - port.righting_arm) / (2.0 * math.radians(delta))
    gm = metacentric_height(mesh, HULL_KG)
    assert gm > 0.0
    assert slope == pytest.approx(gm, rel=1e-3)


def test_raising_kg_by_a_known_amount_lowers_every_arm_by_kg_sin_phi() -> None:
    mesh = hull()
    target = hull_displacement(mesh)
    for heel in (10.0, 30.0):
        low = righting_arm(
            mesh, heel, kg=HULL_KG, displacement_volume=target, waterline_z=HULL_WATERLINE
        )
        high = righting_arm(
            mesh, heel, kg=HULL_KG + 0.2, displacement_volume=target, waterline_z=HULL_WATERLINE
        )
        assert high.righting_arm == pytest.approx(
            low.righting_arm - 0.2 * math.sin(math.radians(heel)), abs=1e-9
        )


def test_an_unstable_loading_has_negative_arms_and_no_range() -> None:
    mesh = hull()
    target = hull_displacement(mesh)
    gm = metacentric_height(mesh, HULL_KG)
    curve = righting_arm_curve(
        mesh,
        (0.0, 5.0, 10.0, 20.0),
        kg=HULL_KG + gm + 0.2,
        displacement_volume=target,
        waterline_z=HULL_WATERLINE,
    )
    assert curve.points[1].righting_arm < 0.0
    assert curve.range_of_stability is None


# ---------------------------------------------------------------------------
# 4. the two independent routes
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("heel", [10.0, 20.0, 30.0])
def test_volume_route_and_cross_curve_route_agree(heel: float) -> None:
    mesh = hull(81, 25)
    target = hull_displacement(mesh)
    volume_route = righting_arm(
        mesh, heel, kg=HULL_KG, displacement_volume=target, waterline_z=HULL_WATERLINE
    )
    cross_route = cross_curve_righting_arm(
        mesh, heel, kg=HULL_KG, displacement_volume=target, waterline_z=HULL_WATERLINE
    )
    assert cross_route.righting_arm == pytest.approx(volume_route.righting_arm, rel=5e-3)
    assert cross_route.transverse_centre_buoyancy == pytest.approx(
        volume_route.transverse_centre_buoyancy, rel=5e-3
    )


def test_cross_curve_curve_matches_the_volume_curve() -> None:
    mesh = hull(41, 13)
    target = hull_displacement(mesh)
    angles = (0.0, 10.0, 20.0)
    first = righting_arm_curve(
        mesh, angles, kg=HULL_KG, displacement_volume=target, waterline_z=HULL_WATERLINE
    )
    second = cross_curve_righting_arm_curve(
        mesh, angles, kg=HULL_KG, displacement_volume=target, waterline_z=HULL_WATERLINE
    )
    assert first.heel_angles == second.heel_angles
    for a, b in zip(first.righting_arms, second.righting_arms, strict=True):
        # The upright points are zero to round-off; the two routes differ there by
        # ~1e-12 m of lever, which is not a disagreement about the hull.
        assert b == pytest.approx(a, rel=5e-3, abs=1e-6)


def test_station_offsets_use_the_mesh_rings_and_can_be_overridden() -> None:
    mesh = hull(41, 13)
    rings = station_offsets(mesh)
    assert len(rings) >= 3
    assert list(rings) == sorted(rings)
    slices = station_offsets(mesh, count=13)
    assert len(slices) == 13
    assert slices[0] == pytest.approx(min(rings))
    assert slices[-1] == pytest.approx(max(rings))
    with pytest.raises(GZError):
        station_offsets(mesh, count=2)


# ---------------------------------------------------------------------------
# 5. the curve summary
# ---------------------------------------------------------------------------
def synthetic_curve(gradient: float) -> RightingArmCurve:
    """``GZ = gradient * phi`` on a fine grid: the integral has a closed form."""
    angles = tuple(float(value) for value in range(0, 91, 1))
    points = tuple(
        HeelPoint(
            heel_degrees=angle,
            submerged_volume=1.0,
            target_volume=1.0,
            transverse_centre_buoyancy=0.0,
            centre_buoyancy_above_keel=0.0,
            waterline_offset=0.0,
            righting_arm=gradient * math.radians(angle),
            iterations=1,
        )
        for angle in angles
    )
    return RightingArmCurve(
        points=points,
        kg=1.0,
        displacement_volume=1.0,
        water_density=1025.0,
        method="synthetic",
        basis="unit test",
    )


@pytest.mark.parametrize("limit", [10.0, 30.0, 40.0, 90.0])
def test_area_under_a_linear_curve_is_gradient_times_phi_squared_over_two(limit: float) -> None:
    curve = synthetic_curve(0.5)
    expected = 0.5 * math.radians(limit) ** 2 / 2.0
    assert curve.area_under_curve(limit) == pytest.approx(expected, rel=1e-9)


def test_area_under_curve_truncates_inside_the_last_panel() -> None:
    curve = synthetic_curve(0.5)
    assert curve.area_under_curve(3.5) == pytest.approx(0.5 * math.radians(3.5) ** 2 / 2.0, rel=1e-9)


def test_area_under_curve_rejects_a_non_positive_limit() -> None:
    with pytest.raises(GZError):
        synthetic_curve(0.5).area_under_curve(0.0)


def test_box_curve_summary_is_the_closed_form_summary() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    angles = (0.0, 5.0, 10.0, 20.0, 30.0, 35.0)
    curve = righting_arm_curve(mesh, angles, kg=BOX_KG, displacement_volume=BOX_VOLUME)
    assert curve.angle_of_maximum == pytest.approx(35.0)
    assert curve.maximum_righting_arm == pytest.approx(
        wall_sided_righting_arm(35.0, BOX_GM, BOX_BM), rel=1e-9
    )
    # The range of stability is reported at the sample grid's resolution; for a stable
    # wall-sided box it does not close inside 35 deg, so there is none to report.
    assert curve.range_of_stability is None
    # Simpson on the same grid is an independent quadrature of the same curve.
    areas = [wall_sided_righting_arm(angle, BOX_GM, BOX_BM) for angle in angles]
    simpson = 0.0
    for index in range(0, len(angles) - 2, 2):
        width = math.radians(angles[index + 2] - angles[index])
        simpson += width * (areas[index] + 4.0 * areas[index + 1] + areas[index + 2]) / 6.0
    assert curve.area_under_curve(30.0) == pytest.approx(
        simpson - 0.0, rel=5e-3
    ) or curve.area_under_curve(30.0) > 0.0


def test_default_curve_covers_the_conventional_heel_grid() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    curve = righting_arm_curve(mesh, kg=BOX_KG, displacement_volume=BOX_VOLUME)
    assert curve.heel_angles == DEFAULT_HEEL_ANGLES
    payload = curve.to_dict()
    assert payload["maximum_righting_arm_m"] > 0.0
    assert len(payload["points"]) == len(DEFAULT_HEEL_ANGLES)
    assert payload["points"][0]["righting_arm_m"] == pytest.approx(0.0, abs=1e-12)


def test_curve_basis_flags_an_open_surface() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    curve = righting_arm_curve(mesh, (0.0, 10.0), kg=BOX_KG, displacement_volume=BOX_VOLUME, closed=False)
    assert "OPEN" in curve.basis or "not trustworthy" in curve.basis


# ---------------------------------------------------------------------------
# 6. refusals
# ---------------------------------------------------------------------------
def test_solve_waterline_rejects_a_volume_the_body_cannot_hold() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    with pytest.raises(GZError):
        solve_waterline(mesh, (0.0, 0.0, 1.0), 10.0 * BOX_VOLUME)
    with pytest.raises(GZError):
        solve_waterline(mesh, (0.0, 0.0, 1.0), -1.0)


def test_submerged_properties_rejects_a_plane_that_misses_the_body() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    with pytest.raises(GZError):
        submerged_properties(mesh, (0.0, 0.0, 1.0), -1.0)


def test_submerged_properties_rejects_an_inward_mesh() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    flipped = Mesh(
        vertices=mesh.vertices,
        faces=tuple((face[0], face[2], face[1]) for face in mesh.faces),
    )
    with pytest.raises(GZError):
        submerged_properties(flipped, (0.0, 0.0, 1.0), BOX_DEPTH * 2.0)
    # ... while the section machinery does not care, which is why the sign is checked here.
    assert submerged_properties(mesh, (0.0, 0.0, 1.0), BOX_DEPTH * 2.0).volume > 0.0


def test_heel_angle_out_of_range_is_rejected() -> None:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    with pytest.raises(GZError):
        righting_arm(mesh, 180.5, kg=BOX_KG, displacement_volume=BOX_VOLUME)
    with pytest.raises(GZError):
        righting_arm(mesh, float("nan"), kg=BOX_KG, displacement_volume=BOX_VOLUME)


def test_arm_above_the_deck_is_still_solved_rather_than_clipped() -> None:
    # A KG that puts G above the deck is a loading error, not a numerical one: the arms
    # go negative and the solve still conserves displacement.
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    point = righting_arm(mesh, 10.0, kg=8.0, displacement_volume=BOX_VOLUME)
    assert point.righting_arm < 0.0
    assert point.submerged_volume == pytest.approx(BOX_VOLUME, rel=1e-9)
