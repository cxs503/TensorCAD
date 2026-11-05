"""Righting arms: closed forms for a barge and a sphere, and a hull ladder.

The two exact references are deliberately computed by code that shares nothing with
:mod:`tensorcad.righting_arm`: a 2-D polygon clip of the heeled barge section (plane
geometry, no surface integral) and the sphere's ``GM sin(phi)``.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.geometry import MeshBuilder
from tensorcad.geometry_field import SphereField, to_mesh
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.mesh import is_watertight
from tensorcad.righting_arm import (
    DEFAULT_HEELS_DEG,
    HeeledDisplacement,
    RightingArmError,
    heeled_displacement,
    righting_arm,
    righting_arm_curve,
    wall_sided_righting_arm,
    wall_sided_validity_deg,
)
from tensorcad.stability import waterplane_moments

LENGTH, BEAM, DRAFT, FREEBOARD = 10.0, 2.0, 0.8, 0.4
DEPTH = DRAFT + FREEBOARD


def box(
    length: float,
    beam: float,
    draft: float,
    freeboard: float = 0.0,
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
) -> MeshBuilder:
    """Closed box with its bottom at ``-draft`` and its top at ``freeboard``."""
    builder = MeshBuilder()
    corner = {}
    ox, oy, oz = origin
    for x in (ox, ox + length):
        for y in (oy - 0.5 * beam, oy + 0.5 * beam):
            for z in (oz - draft, oz + freeboard):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi = oy - 0.5 * beam, oy + 0.5 * beam
    bottom, top = oz - draft, oz + freeboard
    left, right = ox, ox + length
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder


def barge() -> object:
    return box(LENGTH, BEAM, DRAFT, FREEBOARD).build()


# ---------------------------------------------------------------------------
# independent reference: 2-D polygon clip of the heeled barge section
# ---------------------------------------------------------------------------
def _rotate(point: tuple[float, float], phi: float) -> tuple[float, float]:
    y, z = point
    return (y * math.cos(phi) + z * math.sin(phi), -y * math.sin(phi) + z * math.cos(phi))


def _clip(polygon: list[tuple[float, float]], height: float) -> list[tuple[float, float]]:
    out: list[tuple[float, float]] = []
    for index in range(len(polygon)):
        current = polygon[index]
        following = polygon[(index + 1) % len(polygon)]
        dc = current[1] - height
        df = following[1] - height
        if dc <= 0.0:
            out.append(current)
        if (dc < 0.0 < df) or (df < 0.0 < dc):
            t = dc / (dc - df)
            out.append((current[0] + t * (following[0] - current[0]), height))
    return out


def _polygon_moments(polygon: list[tuple[float, float]]) -> tuple[float, float, float]:
    """``(area, int y dA, int z dA)`` of an implicitly closed polygon (shoelace)."""
    area = 0.0
    first_y = 0.0
    first_z = 0.0
    for index in range(len(polygon)):
        y1, z1 = polygon[index]
        y2, z2 = polygon[(index + 1) % len(polygon)]
        cross = y1 * z2 - y2 * z1
        area += cross
        first_y += (y1 + y2) * cross
        first_z += (z1 + z2) * cross
    return (0.5 * area, first_y / 6.0, first_z / 6.0)


def barge_section(
    *,
    beam: float,
    bottom: float,
    top: float,
    length: float,
    heel_deg: float,
    volume: float,
) -> tuple[float, float, float]:
    """``(waterline, y_B, z_B)`` of the heeled rectangular section, by plane geometry."""
    phi = math.radians(heel_deg)
    half = 0.5 * beam
    section = [
        _rotate(point, phi)
        for point in ((half, bottom), (half, top), (-half, top), (-half, bottom))
    ]
    low = min(point[1] for point in section)
    high = max(point[1] for point in section)
    for _ in range(200):
        middle = 0.5 * (low + high)
        area = _polygon_moments(_clip(section, middle))[0]
        if area * length < volume:
            low = middle
        else:
            high = middle
    height = 0.5 * (low + high)
    area, moment_y, moment_z = _polygon_moments(_clip(section, height))
    return height, moment_y / area, moment_z / area


# ---------------------------------------------------------------------------
# the box: both references are exact
# ---------------------------------------------------------------------------
def test_upright_integral_reproduces_the_closed_mesh_volume() -> None:
    """A mesh capped at its waterline: two independent integrations must agree."""
    capped = box(LENGTH, BEAM, DRAFT).build()
    hydrostatics = compute_hydrostatics(capped, 0.0)
    state = heeled_displacement(capped, heel_deg=0.0, waterline_z=0.0)
    # The solved waterline is the bisection tolerance away from the exwaterline, so
    # the volumes agree to that tolerance rather than bit-exactly.
    assert state.volume == pytest.approx(hydrostatics.displacement_volume, rel=1e-10)
    assert state.centre_of_buoyancy[2] == pytest.approx(
        hydrostatics.vertical_centre_buoyancy, rel=1e-10
    )


def test_upright_displacement_of_a_freeboard_box_is_the_immersed_part() -> None:
    """Freeboard above the waterline must not be counted -- the whole box is deeper."""
    mesh = barge()
    state = heeled_displacement(mesh, heel_deg=0.0, waterline_z=0.0)
    assert state.volume == pytest.approx(LENGTH * BEAM * DRAFT, rel=1e-10)
    # ... while the enclosed solid is the larger number, which is what the
    # hydrostatics module would report for this mesh.
    assert compute_hydrostatics(mesh, 0.0).displacement_volume > state.volume


@pytest.mark.parametrize("heel", [0.0, 1.0, 5.0, 10.0, 15.0, 20.0])
def test_box_righting_arm_reproduces_the_wall_sided_closed_form(heel: float) -> None:
    """Within its validity window the wall-sided formula is exact, so this is a check."""
    result = righting_arm(barge(), heel_deg=heel, kg=0.55)
    expected = wall_sided_righting_arm(breadth=BEAM, draft=DRAFT, kg=0.55, heel_deg=heel)
    assert result.righting_arm == pytest.approx(expected, abs=1e-12)


@pytest.mark.parametrize("heel", [0.0, 5.0, 10.0, 20.0, 25.0, 30.0, 35.0, 40.0])
def test_box_righting_arm_matches_the_two_dimensional_reference(heel: float) -> None:
    """Beyond the wall-sided window this reference still holds: it clips a polygon."""
    mesh = barge()
    state = heeled_displacement(mesh, heel_deg=heel, waterline_z=0.0)
    height, y_b, z_b = barge_section(
        beam=BEAM,
        bottom=-DRAFT,
        top=FREEBOARD,
        length=LENGTH,
        heel_deg=heel,
        volume=state.volume,
    )
    assert state.waterline_z == pytest.approx(height, abs=1e-9)
    assert state.centre_of_buoyancy[1] == pytest.approx(y_b, abs=1e-9)
    assert state.centre_of_buoyancy[2] == pytest.approx(z_b, abs=1e-9)
    # GZ from the reference: the vertical through G is at y_G = (z_keel + KG) sin(phi).
    phi = math.radians(heel)
    expected = y_b - (-DRAFT + 0.55) * math.sin(phi)
    assert righting_arm(mesh, heel_deg=heel, kg=0.55).righting_arm == pytest.approx(
        expected, abs=1e-9
    )


def test_the_wall_sided_formula_is_wrong_past_its_own_validity_angle() -> None:
    """The formula's stated limit is real: beyond it the two references disagree."""
    limit = wall_sided_validity_deg(breadth=BEAM, draft=DRAFT, depth=DEPTH)
    assert limit == pytest.approx(math.degrees(math.atan2(2.0 * FREEBOARD, BEAM)) , abs=1e-12)
    mesh = barge()
    for heel in (limit + 3.0, limit + 8.0):
        state = heeled_displacement(mesh, heel_deg=heel, waterline_z=0.0)
        _, y_b, _ = barge_section(
            beam=BEAM, bottom=-DRAFT, top=FREEBOARD, length=LENGTH,
            heel_deg=heel, volume=state.volume,
        )
        measured = righting_arm(mesh, heel_deg=heel, kg=0.55).righting_arm
        phi = math.radians(heel)
        expected = y_b - (-DRAFT + 0.55) * math.sin(phi)
        formula = wall_sided_righting_arm(breadth=BEAM, draft=DRAFT, kg=0.55, heel_deg=heel)
        assert measured == pytest.approx(expected, abs=1e-9)
        assert abs(measured - formula) > 1e-3


def test_a_heavily_loaded_box_still_reproduces_the_wall_sided_form() -> None:
    """The formula has no stability requirement: KG above KM gives negative arms."""
    mesh = barge()
    for heel in (5.0, 10.0, 20.0):
        result = righting_arm(mesh, heel_deg=heel, kg=1.5)
        expected = wall_sided_righting_arm(breadth=BEAM, draft=DRAFT, kg=1.5, heel_deg=heel)
        assert result.righting_arm == pytest.approx(expected, abs=1e-12)
        assert result.righting_arm < 0.0


def test_gz_vanishes_upright_and_is_antisymmetric() -> None:
    mesh = barge()
    assert righting_arm(mesh, heel_deg=0.0, kg=0.55).righting_arm == pytest.approx(0.0, abs=1e-15)
    for heel in (5.0, 20.0, 35.0):
        port = righting_arm(mesh, heel_deg=heel, kg=0.55).righting_arm
        starboard = righting_arm(mesh, heel_deg=-heel, kg=0.55).righting_arm
        assert port == pytest.approx(-starboard, abs=1e-12)


def test_small_angle_slope_is_the_wall_sided_metacentric_height() -> None:
    mesh = barge()
    heel = 1e-3
    arm = righting_arm(mesh, heel_deg=heel, kg=0.55).righting_arm
    slope = arm / math.sin(math.radians(heel))
    kb = 0.5 * DRAFT
    bm = BEAM * BEAM / (12.0 * DRAFT)
    assert slope == pytest.approx(kb + bm - 0.55, rel=1e-9)


def test_constant_displacement_is_held_at_every_heel() -> None:
    mesh = barge()
    curve = righting_arm_curve(mesh, kg=0.55, heels_deg=(0.0, 10.0, 20.0, 40.0))
    volumes = {round(point.displacement_volume, 9) for point in curve.points}
    assert len(volumes) == 1
    assert volumes == {round(LENGTH * BEAM * DRAFT, 9)}
    assert all(abs(point.volume_residual) < 1e-9 for point in curve.points)


def test_the_integral_is_translation_invariant() -> None:
    """The divergence integrals are taken about a shifted origin; that must not show."""
    offset = (137.0, 0.0, 61.0)
    base = barge()
    moved = box(LENGTH, BEAM, DRAFT, FREEBOARD, origin=offset).build()
    for heel in (0.0, 15.0, 35.0):
        reference = righting_arm(base, heel_deg=heel, kg=0.55)
        translated = righting_arm(moved, heel_deg=heel, kg=0.55, waterline_z=offset[2])
        assert translated.waterline_z == pytest.approx(
            reference.waterline_z + offset[2], abs=1e-9
        )
        assert translated.centre_of_buoyancy[1] == pytest.approx(
            reference.centre_of_buoyancy[1] + offset[1], abs=1e-9
        )
        assert translated.righting_arm == pytest.approx(reference.righting_arm, abs=1e-12)


# ---------------------------------------------------------------------------
# the sphere: GZ = GM sin(phi) at every heel
# ---------------------------------------------------------------------------
def sphere_mesh(resolution: int):  # type: ignore[no-untyped-def]
    field = to_mesh(SphereField(radius=1.0), resolution)
    return field.mesh


@pytest.mark.parametrize("resolution", [12, 16, 24])
def test_sphere_righting_arm_converges_to_gm_times_sin(resolution: int) -> None:
    """A floating sphere's buoyancy centre stays on the vertical: ``GZ = GM sin(phi)``."""
    mesh = sphere_mesh(resolution)
    kg = 0.9  # the centre of the sphere is at kg = 1.0, so GM = 0.1 m
    for heel in (10.0, 30.0, 45.0):
        result = righting_arm(mesh, heel_deg=heel, kg=kg)
        expected = (1.0 - kg) * math.sin(math.radians(heel))
        # a polyhedral sphere's immersed cap centroid is not exactly the sphere's
        assert result.righting_arm == pytest.approx(expected, rel=6e-2)


def test_sphere_with_gravity_at_its_centre_is_neutral() -> None:
    """``GM = 0``: a sphere rolled about its centre has no righting arm at any heel."""
    mesh = sphere_mesh(24)
    for heel in (10.0, 30.0, 60.0):
        assert abs(righting_arm(mesh, heel_deg=heel, kg=1.0).righting_arm) < 1e-3


def test_sphere_error_shrinks_with_refinement() -> None:
    errors = []
    for resolution in (12, 24):
        mesh = sphere_mesh(resolution)
        measured = righting_arm(mesh, heel_deg=45.0, kg=0.9).righting_arm
        errors.append(abs(measured - 0.1 * math.sin(math.radians(45.0))))
    assert errors[1] < 0.5 * errors[0]


# ---------------------------------------------------------------------------
# the parametric hull: no closed-form GZ, so a ladder and the waterplane moments
# ---------------------------------------------------------------------------
HULL_WATERLINE = -0.47


def hull_mesh(nx: int, nz: int):  # type: ignore[no-untyped-def]
    """A hull generated with depth 1 m and floated at -0.47 m, i.e. with freeboard."""
    return generate_hull(
        HullParameters(length=10.0, beam=2.0, draft=1.0, nx=nx, nz=nz)
    ).mesh


def test_hull_small_angle_slope_matches_the_waterplane_moments() -> None:
    """Two independent integrations of the same hull: heeled volume vs waterplane I."""
    mesh = hull_mesh(41, 13)
    state = heeled_displacement(mesh, heel_deg=0.0, waterline_z=HULL_WATERLINE)
    section = waterplane_moments(mesh, HULL_WATERLINE)
    keel = min(vertex[2] for vertex in mesh.vertices)
    gm = (
        state.centre_of_buoyancy[2]
        - keel
        + section.transverse_second_moment / state.volume
        - 0.55
    )
    heel = 1e-2
    slope = righting_arm(
        mesh, heel_deg=heel, kg=0.55, waterline_z=HULL_WATERLINE
    ).righting_arm / math.sin(math.radians(heel))
    assert slope == pytest.approx(gm, rel=1e-6)


def test_hull_righting_arm_converges_with_refinement() -> None:
    arms = {
        (nx, nz): righting_arm(
            hull_mesh(nx, nz), heel_deg=30.0, kg=0.55, waterline_z=HULL_WATERLINE
        ).righting_arm
        for (nx, nz) in ((41, 13), (81, 25), (161, 41))
    }
    coarse = abs(arms[(41, 13)] - arms[(81, 25)])
    fine = abs(arms[(81, 25)] - arms[(161, 41)])
    assert fine < coarse
    # and the arms are physical: positive, and larger than the upright 0 at 30 deg
    assert all(arm > 0.0 for arm in arms.values())


def test_hull_curve_reports_a_peak_and_a_positive_range() -> None:
    curve = righting_arm_curve(
        hull_mesh(41, 13),
        kg=0.55,
        waterline_z=HULL_WATERLINE,
        heels_deg=(0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0),
    )
    peak_heel, peak_arm = curve.peak()
    assert peak_heel > 20.0
    assert peak_arm > 0.05
    assert curve.gz_at(0.0) == pytest.approx(0.0, abs=1e-12)
    assert not curve.vanishing_crossed()
    assert curve.positive_range_deg() == 60.0
    assert curve.area_under_curve() > 0.0


def test_a_mesh_capped_at_its_waterline_has_no_freeboard_and_the_waterline_climbs() -> None:
    """The documented limitation, asserted: a capped hull is not a ship with freeboard.

    ``generate_hull`` closes the hull *on* its design waterline. Heel it and the
    constant-displacement waterline has to climb over that cap, so the solved height
    rises above zero and the arm stops being ``GM sin(phi)``. The code reports the
    waterline, which is what makes this detectable instead of silently wrong.
    """
    capped = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13)).mesh
    upright = heeled_displacement(capped, heel_deg=0.0, waterline_z=0.0)
    heeled = heeled_displacement(capped, heel_deg=20.0, waterline_z=0.0)
    assert upright.waterline_z == pytest.approx(0.0, abs=1e-12)
    assert heeled.waterline_z > 0.2
    # the warmed-over GM of that solid is not the arm's slope
    assert righting_arm(capped, heel_deg=20.0, kg=0.55).righting_arm < 0.0


# ---------------------------------------------------------------------------
# curve helpers
# ---------------------------------------------------------------------------
def test_curve_interpolates_and_extrapolation_is_refused() -> None:
    curve = righting_arm_curve(barge(), kg=0.55, heels_deg=(0.0, 10.0, 20.0))
    assert curve.gz_at(10.0) == pytest.approx(curve.points[1].righting_arm, abs=1e-15)
    middle = curve.gz_at(15.0)
    assert middle == pytest.approx(
        0.5 * (curve.points[1].righting_arm + curve.points[2].righting_arm), abs=1e-15
    )
    with pytest.raises(RightingArmError):
        curve.gz_at(25.0)


def test_curve_area_is_the_trapezoid_of_the_arms() -> None:
    curve = righting_arm_curve(barge(), kg=0.55, heels_deg=(0.0, 10.0, 20.0))
    arms = [point.righting_arm for point in curve.points]
    expected = 0.5 * arms[0] * 0.0 + 0.5 * (arms[0] + arms[1]) * math.radians(10.0) + 0.5 * (
        arms[1] + arms[2]
    ) * math.radians(10.0)
    assert curve.area_under_curve() == pytest.approx(expected, rel=1e-12)


def test_curve_range_reports_whether_a_crossing_was_actually_sampled() -> None:
    """A stable barge's curve never returns to zero on the ladder; a heavy one starts at it."""
    stable = righting_arm_curve(barge(), kg=0.55, heels_deg=(0.0, 10.0, 20.0, 30.0))
    assert not stable.vanishing_crossed()
    assert stable.positive_range_deg() == 30.0
    heavy = righting_arm_curve(barge(), kg=1.75, heels_deg=(0.0, 5.0, 10.0, 15.0, 20.0))
    assert heavy.vanishing_crossed()
    assert heavy.positive_range_deg() == pytest.approx(0.0, abs=0.1)
    assert all(point.righting_arm <= 0.0 for point in heavy.points[1:])


def test_curve_dictionary_is_serialisable_and_matches_the_points() -> None:
    curve = righting_arm_curve(barge(), kg=0.55)
    payload = curve.to_dict()
    assert payload["schema"] == "tensorcad.stability.gzcurve/1"
    assert len(payload["points"]) == len(DEFAULT_HEELS_DEG)
    assert payload["heels_deg"][0] == 0.0
    assert payload["peak"]["righting_arm"] == pytest.approx(curve.peak()[1], abs=1e-9)
    assert json_round_trip(payload) == payload


def json_round_trip(payload: dict) -> dict:
    import json

    return json.loads(json.dumps(payload))


def test_heeled_displacement_dictionary_is_serialisable() -> None:
    state = heeled_displacement(barge(), heel_deg=15.0, waterline_z=0.0)
    assert isinstance(state, HeeledDisplacement)
    payload = state.to_dict()
    assert payload["schema"] == "tensorcad.stability.heeled/1"
    assert payload["heel_deg"] == 15.0
    assert payload["iterations"] > 0


# ---------------------------------------------------------------------------
# refusals: a silently wrong arm is worse than no arm
# ---------------------------------------------------------------------------
def test_heel_beyond_the_transverse_convention_is_refused() -> None:
    with pytest.raises(RightingArmError):
        righting_arm(barge(), heel_deg=95.0, kg=0.55)
    with pytest.raises(RightingArmError):
        righting_arm(barge(), heel_deg=float("nan"), kg=0.55)


def test_inputs_must_be_finite_and_positive(tmp_path) -> None:  # type: ignore[no-untyped-def]
    mesh = barge()
    with pytest.raises(RightingArmError):
        righting_arm(mesh, heel_deg=10.0, kg=float("inf"))
    with pytest.raises(RightingArmError):
        righting_arm(mesh, heel_deg=10.0, kg=0.55, lcg=float("nan"))
    with pytest.raises(RightingArmError):
        righting_arm(mesh, heel_deg=10.0, kg=0.55, waterline_z=float("nan"))
    with pytest.raises(RightingArmError):
        righting_arm(mesh, heel_deg=10.0, kg=0.55, displacement_volume=0.0)
    with pytest.raises(RightingArmError):
        righting_arm(mesh, heel_deg=10.0, kg=0.55, displacement_volume=-3.0)
    with pytest.raises(RightingArmError):
        righting_arm(mesh, heel_deg=10.0, kg=0.55, displacement_volume=1e9)
    with pytest.raises(RightingArmError):
        righting_arm(mesh, heel_deg=10.0, kg=0.55, relative_tolerance=0.0)


def test_an_open_mesh_is_refused() -> None:
    builder = box(LENGTH, BEAM, DRAFT, FREEBOARD)
    mesh = builder.build()
    open_mesh = type(mesh)(vertices=mesh.vertices, faces=mesh.faces[:-1], units=mesh.units)
    assert not is_watertight(open_mesh)
    with pytest.raises(RightingArmError):
        righting_arm(open_mesh, heel_deg=10.0, kg=0.55)
    with pytest.raises(RightingArmError):
        righting_arm(mesh, heel_deg=10.0, kg=0.55, closed=False)


def test_curve_ladder_must_increase_and_have_two_points() -> None:
    with pytest.raises(RightingArmError):
        righting_arm_curve(barge(), kg=0.55, heels_deg=(0.0,))
    with pytest.raises(RightingArmError):
        righting_arm_curve(barge(), kg=0.55, heels_deg=(0.0, 10.0, 10.0))


def test_wall_sided_helpers_validate_their_inputs() -> None:
    with pytest.raises(RightingArmError):
        wall_sided_righting_arm(breadth=0.0, draft=1.0, kg=0.5, heel_deg=10.0)
    with pytest.raises(RightingArmError):
        wall_sided_righting_arm(breadth=2.0, draft=1.0, kg=0.5, heel_deg=120.0)
    with pytest.raises(RightingArmError):
        wall_sided_validity_deg(breadth=2.0, draft=1.0, depth=0.5)
    assert wall_sided_validity_deg(breadth=BEAM, draft=DRAFT) == pytest.approx(
        math.degrees(math.atan2(2.0 * DRAFT, BEAM))
    )
