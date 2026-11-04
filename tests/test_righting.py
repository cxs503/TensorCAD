"""The righting-arm curve, against the one closed form that exists for it.

``GZ`` has no general analytic solution, so the module is checked where it can be:
a rectangular barge, for which the wall-sided formula

    ``GZ = sin(phi) (GM + BM tan^2(phi) / 2)``

is *exact* up to the angle at which the waterplane leaves the vertical sides.
Because a box section is a polygon and the module clips polygons exactly, that
case pins the whole chain -- the tilted half-plane sign convention, the Green
moments of the clipped section, the Simpson integration along the hull, the
constant-displacement bisection and the lever formula -- to round-off, and any
sign error anywhere in it shows up immediately rather than as "a plausible
curve".

The second family of checks is on a real hull, where the reference is not a
formula but a different integration: the slope of ``GZ`` at zero heel must equal
``KB + BM - KG``, with ``BM = I_T / V`` taken from
:func:`tensorcad.stability.waterplane_moments`, which cuts the surface rather
than integrating station areas. That is the independent half of the check; the
curves' own convergence with the station grid is the rest of it.
"""
from __future__ import annotations

import dataclasses
import json
import math

import pytest

import tensorcad as tc
from tensorcad.righting import (
    DEFAULT_HEEL_ANGLES_DEG,
    DEFAULT_STATION_INTERVALS,
    MINIMUM_STATION_INTERVALS,
    RightingArmPoint,
    RightingError,
    righting_arm_curve,
    wall_sided_righting_arm,
)

LENGTH = 10.0
BEAM = 2.0
DRAFT = 1.0
FREEBOARD = 1.0


def box_barge(
    length: float = LENGTH, beam: float = BEAM, draft: float = DRAFT, freeboard: float = FREEBOARD
) -> tc.Mesh:
    """A watertight rectangular barge: draft below ``z = 0``, deck at ``z = freeboard``."""
    half = beam / 2.0
    vertices = (
        (0.0, -half, -draft),
        (length, -half, -draft),
        (length, half, -draft),
        (0.0, half, -draft),
        (0.0, -half, freeboard),
        (length, -half, freeboard),
        (length, half, freeboard),
        (0.0, half, freeboard),
    )
    faces = (
        (0, 2, 1), (0, 3, 2),   # bottom
        (4, 5, 6), (4, 6, 7),   # deck
        (0, 1, 5), (0, 5, 4),   # port
        (1, 2, 6), (1, 6, 5),   # forward
        (2, 3, 7), (2, 7, 6),   # starboard
        (3, 0, 4), (3, 4, 7),   # aft
    )
    return tc.Mesh(vertices=vertices, faces=faces)


def wall_sided_reference(kg: float) -> tuple[float, float]:
    """``(GM, BM)`` of the barge box section, exact: ``KB = T/2``, ``BM = B^2/(12T)``."""
    kb = DRAFT / 2.0
    bm = BEAM**2 / (12.0 * DRAFT)
    return (kb + bm - kg, bm)


def hull_mesh(nx: int = 41, nz: int = 13) -> tc.Mesh:
    return tc.generate_hull(tc.HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=nx, nz=nz)).mesh


# --------------------------------------------------------------------------- barge


def test_barge_mesh_is_a_valid_solid() -> None:
    mesh = box_barge()
    assert tc.mesh.is_watertight(mesh)
    assert tc.mesh.total_signed_volume(mesh) == pytest.approx(LENGTH * BEAM * (DRAFT + FREEBOARD))


def test_barge_matches_the_wall_sided_closed_form() -> None:
    """Below the deck-edge immersion angle the arm is a closed form; match it exactly."""
    kg = 0.8
    gm, bm = wall_sided_reference(kg)
    angles = (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0)
    curve = righting_arm_curve(box_barge(), kg=kg, heel_angles_deg=angles, station_intervals=8)
    # (B/2) tan(45) = freeboard, so 45 degrees is the last angle that is still
    # wall-sided; the formula is exact there and the mesh is a polygon, so the
    # agreement is round-off rather than "a few percent".
    for point in curve.points:
        reference = wall_sided_righting_arm(point.heel_deg, gm=gm, bm=bm)
        assert point.gz == pytest.approx(reference, rel=1e-12, abs=1e-12)
        assert point.converged
        if point.heel_deg < 45.0:
            assert not point.deck_edge_immersed


def test_barge_deck_edge_immerses_at_the_classic_angle() -> None:
    # (B/2) tan(phi) = freeboard  ->  phi = 45 degrees for B = 2, freeboard = 1.
    curve = righting_arm_curve(
        box_barge(), kg=0.8, heel_angles_deg=tuple(float(a) for a in range(0, 91, 5)), station_intervals=8
    )
    assert curve.deck_immersion_onset_deg == pytest.approx(45.0)
    assert curve.to_dict()["deck_immersion_onset_deg"] == pytest.approx(45.0)


def test_barge_arm_keeps_rising_past_the_deck_edge_but_not_by_the_closed_form() -> None:
    """Past deck-edge immersion the closed form is invalid; the integral keeps going."""
    kg = 0.8
    gm, bm = wall_sided_reference(kg)
    curve = righting_arm_curve(
        box_barge(),
        kg=kg,
        heel_angles_deg=tuple(float(a) for a in range(0, 91, 5)),
        station_intervals=8,
    )
    at_50 = curve.gz_at(50.0)
    reference_50 = wall_sided_righting_arm(50.0, gm=gm, bm=bm)
    assert at_50 < reference_50
    assert curve.gz_max == pytest.approx(0.2374, abs=5e-4)
    assert curve.angle_gz_max_deg == pytest.approx(68.5, abs=0.5)


def test_barge_never_vanishes_because_the_deck_is_watertight() -> None:
    curve = righting_arm_curve(
        box_barge(), kg=0.8, heel_angles_deg=tuple(float(a) for a in range(0, 91, 5)), station_intervals=8
    )
    assert all(point.gz > 0.0 for point in curve.points[1:])
    assert math.isnan(curve.angle_vanishing_stability_deg)


# --------------------------------------------------------------------------- limits


def test_zero_heel_arm_is_zero() -> None:
    curve = righting_arm_curve(box_barge(), kg=0.8, heel_angles_deg=(0.0, 5.0), station_intervals=8)
    assert abs(curve.points[0].gz) < 1e-12
    assert curve.points[0].transverse_centre_buoyancy == pytest.approx(0.0, abs=1e-12)


def test_slope_at_zero_heel_is_the_metacentric_height() -> None:
    """Richardson extrapolation of ``GZ/phi`` must recover ``GM`` -- the module's own check."""
    kg = 0.8
    gm, _bm = wall_sided_reference(kg)
    curve = righting_arm_curve(
        box_barge(), kg=kg, heel_angles_deg=(0.0, 0.5, 1.0), station_intervals=8
    )
    one, half = curve.points[2], curve.points[1]
    coarse = one.gz / math.radians(one.heel_deg)
    fine = half.gz / math.radians(half.heel_deg)
    slope = (4.0 * fine - coarse) / 3.0
    assert slope == pytest.approx(gm, rel=1e-6)


def test_displacement_is_conserved_at_every_heel() -> None:
    curve = righting_arm_curve(
        box_barge(), kg=0.8, heel_angles_deg=(0.0, 10.0, 30.0, 60.0, 90.0), station_intervals=8
    )
    for point in curve.points:
        assert point.submerged_volume == pytest.approx(LENGTH * BEAM * DRAFT, rel=1e-10)
    assert curve.displacement_volume == pytest.approx(LENGTH * BEAM * DRAFT, rel=1e-10)


def test_offset_yields_the_arm_by_definition() -> None:
    """GZ = (y_B - y_G) cos(phi) + (z_B - z_G) sin(phi), recomputed from the point."""
    kg = 0.8
    curve = righting_arm_curve(box_barge(), kg=kg, heel_angles_deg=(0.0, 20.0, 40.0), station_intervals=8)
    z_g = curve.keel_z + kg
    for point in curve.points[1:]:
        phi = math.radians(point.heel_deg)
        expected = point.transverse_centre_buoyancy * math.cos(phi) + (
            point.vertical_centre_buoyancy - z_g
        ) * math.sin(phi)
        assert point.gz == pytest.approx(expected, abs=1e-12)


def test_off_centreline_gravity_shifts_the_upright_arm() -> None:
    curve = righting_arm_curve(
        box_barge(), kg=0.8, y_g=0.1, heel_angles_deg=(0.0, 10.0), station_intervals=8
    )
    assert curve.points[0].gz == pytest.approx(-0.1, abs=1e-12)
    symmetric = righting_arm_curve(box_barge(), kg=0.8, heel_angles_deg=(0.0, 10.0), station_intervals=8)
    shift = symmetric.points[1].gz - curve.points[1].gz
    assert shift == pytest.approx(0.1 * math.cos(math.radians(10.0)))


def test_higher_kg_gives_a_smaller_arm() -> None:
    low = righting_arm_curve(box_barge(), kg=0.7, heel_angles_deg=(0.0, 20.0), station_intervals=8)
    high = righting_arm_curve(box_barge(), kg=0.9, heel_angles_deg=(0.0, 20.0), station_intervals=8)
    assert low.points[1].gz - high.points[1].gz == pytest.approx(0.2 * math.sin(math.radians(20.0)))


# --------------------------------------------------------------------------- hull


def test_hull_freeboard_is_required() -> None:
    """A hull capped at its waterline has no freeboard: every heel floods the deck."""
    with pytest.raises(RightingError, match="freeboard"):
        righting_arm_curve(hull_mesh(), kg=0.6, waterline_z=0.0)


def test_hull_displacement_matches_the_upright_surface_integral() -> None:
    """At the waterline the two integrations must agree: volume below the cap is the solid."""
    mesh = hull_mesh()
    hydrostatics = tc.compute_hydrostatics(mesh, 0.0)
    # Float it lower (freeboard 0.3) so the curve is legal, then check the section
    # integral against the closed-surface volume at the *same* waterline.
    curve = righting_arm_curve(
        mesh, kg=0.6, waterline_z=-0.3, heel_angles_deg=(0.0, 5.0), station_intervals=80
    )
    assert curve.displacement_volume < hydrostatics.displacement_volume
    assert curve.freeboard == pytest.approx(0.3)
    # Floating 0.3 m lower removes the widest slices, so the displaced volume drops by
    # more than the 30 % a linear share of the draft would suggest.
    share = curve.displacement_volume / hydrostatics.displacement_volume
    assert share == pytest.approx(0.583, abs=0.01)


def test_hull_slope_matches_independent_waterplane_moments() -> None:
    """The dominated-by-BM slope must match ``KB + I_T/V - KG`` from another module."""
    mesh = hull_mesh()
    kg = 0.75
    waterline = -0.3
    curve = righting_arm_curve(
        mesh, kg=kg, waterline_z=waterline, heel_angles_deg=(0.0, 0.5, 1.0), station_intervals=80
    )
    coarse = curve.points[2].gz / math.radians(1.0)
    fine = curve.points[1].gz / math.radians(0.5)
    slope = (4.0 * fine - coarse) / 3.0

    volume = curve.displacement_volume
    kb = curve.points[0].vertical_centre_buoyancy - curve.keel_z
    moments = tc.waterplane_moments(mesh, waterline)
    bm = moments.transverse_second_moment / volume
    assert slope == pytest.approx(kb + bm - kg, rel=2e-3)


def test_hull_curve_reports_a_vanishing_angle() -> None:
    curve = righting_arm_curve(
        hull_mesh(),
        kg=0.75,
        waterline_z=-0.3,
        heel_angles_deg=(0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0),
        station_intervals=80,
    )
    assert curve.gz_max > 0.0
    assert curve.angle_gz_max_deg == pytest.approx(23.4, abs=0.5)
    assert curve.angle_vanishing_stability_deg == pytest.approx(36.0, abs=0.5)
    assert curve.deck_immersion_onset_deg == pytest.approx(20.0)


def test_hull_deck_immersion_fraction_grows_with_heel() -> None:
    curve = righting_arm_curve(
        hull_mesh(),
        kg=0.75,
        waterline_z=-0.3,
        heel_angles_deg=(0.0, 10.0, 20.0, 40.0, 90.0),
        station_intervals=80,
    )
    fractions = [point.deck_immersion_fraction for point in curve.points]
    assert fractions[0] == 0.0 and fractions[1] == 0.0
    assert fractions[2] > 0.0
    assert fractions == sorted(fractions)
    assert curve.deck_immersion_onset_deg == pytest.approx(20.0)


def test_hull_curve_converges_with_the_station_grid() -> None:
    kwargs = {
        "kg": 0.75,
        "waterline_z": -0.3,
        "heel_angles_deg": (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0),
    }
    coarse = righting_arm_curve(hull_mesh(), station_intervals=40, **kwargs)
    fine = righting_arm_curve(hull_mesh(), station_intervals=80, **kwargs)
    assert fine.gz_max == pytest.approx(coarse.gz_max, rel=5e-3)
    assert fine.angle_vanishing_stability_deg == pytest.approx(
        coarse.angle_vanishing_stability_deg, rel=1e-2
    )
    # Doubling again is exact for this faceted body: each ring span then holds two
    # Simpson intervals, so the section-area integral is exact rather than merely
    # refined.
    finer = righting_arm_curve(hull_mesh(), station_intervals=160, **kwargs)
    assert finer.displacement_volume == pytest.approx(fine.displacement_volume, rel=1e-9)


def test_unstable_body_has_a_negative_arm() -> None:
    """KG above KM: the arm is negative, the curve is not "clipped" to look stable."""
    curve = righting_arm_curve(
        box_barge(), kg=0.9, heel_angles_deg=(0.0, 5.0, 10.0), station_intervals=8
    )
    assert all(point.gz < 0.0 for point in curve.points[1:])
    assert math.isnan(curve.gz_max)


# --------------------------------------------------------------------------- interface


def test_default_heel_angles_span_zero_to_ninety() -> None:
    assert DEFAULT_HEEL_ANGLES_DEG[0] == 0.0
    assert DEFAULT_HEEL_ANGLES_DEG[-1] == 90.0
    assert len(DEFAULT_HEEL_ANGLES_DEG) == 19
    curve = righting_arm_curve(box_barge(), kg=0.8, station_intervals=8)
    assert curve.heel_angles_deg == DEFAULT_HEEL_ANGLES_DEG
    assert curve.station_intervals == 8


def test_odd_station_count_is_rounded_up_for_simpson() -> None:
    curve = righting_arm_curve(
        box_barge(), kg=0.8, heel_angles_deg=(0.0, 5.0), station_intervals=MINIMUM_STATION_INTERVALS + 1
    )
    assert curve.station_intervals == MINIMUM_STATION_INTERVALS + 2


def test_gz_at_interpolates_and_refuses_to_extrapolate() -> None:
    curve = righting_arm_curve(
        box_barge(), kg=0.8, heel_angles_deg=(0.0, 10.0, 20.0), station_intervals=8
    )
    middle = 0.5 * (curve.points[1].gz + curve.points[2].gz)
    assert curve.gz_at(15.0) == pytest.approx(middle)
    with pytest.raises(RightingError, match="outside the sampled range"):
        curve.gz_at(25.0)
    with pytest.raises(RightingError, match="outside the sampled range"):
        curve.gz_at(-1.0)


def test_report_is_json_serialisable_and_rounds() -> None:
    curve = righting_arm_curve(
        box_barge(), kg=0.8, heel_angles_deg=(0.0, 5.0, 10.0), station_intervals=8, gm_transverse=0.3667
    )
    payload = json.dumps(curve.to_dict(), allow_nan=False)
    loaded = json.loads(payload)
    assert loaded["schema"] == "tensorcad.righting.curve/1"
    assert loaded["gm_transverse_input_m"] == pytest.approx(0.3667)
    assert loaded["points"][1]["schema"] == "tensorcad.righting.point/1"
    assert loaded["gz_max_m"] == pytest.approx(curve.gz_max, abs=1e-9)
    assert loaded["displacement_volume_m3"] == pytest.approx(20.0)


def test_point_dataclass_is_frozen() -> None:
    point = RightingArmPoint(
        heel_deg=10.0,
        gz=0.01,
        submerged_volume=20.0,
        offset=0.0,
        sinkage=0.0,
        transverse_centre_buoyancy=0.0,
        vertical_centre_buoyancy=-0.5,
    )
    assert not point.deck_edge_immersed
    with pytest.raises(dataclasses.FrozenInstanceError):
        point.gz = 0.0  # type: ignore[misc]


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"kg": float("nan")}, "kg must be finite"),
        ({"kg": 0.8, "waterline_z": float("inf")}, "waterline_z must be finite"),
        ({"kg": 0.8, "y_g": float("nan")}, "y_g must be finite"),
        ({"kg": 0.8, "heel_angles_deg": ()}, "at least one heel angle"),
        ({"kg": 0.8, "heel_angles_deg": (5.0, 10.0)}, "must include 0 degrees"),
        ({"kg": 0.8, "heel_angles_deg": (0.0, 10.0, 5.0)}, "strictly increasing"),
        ({"kg": 0.8, "heel_angles_deg": (0.0, 200.0)}, "must lie in \\[0, 180\\]"),
        ({"kg": 0.8, "heel_angles_deg": (0.0, float("nan"))}, "must be finite"),
        ({"kg": 0.8, "station_intervals": 2}, "at least"),
    ],
)
def test_rejects_invalid_arguments(kwargs: dict, message: str) -> None:
    with pytest.raises(RightingError, match=message):
        righting_arm_curve(box_barge(), **kwargs)


def test_rejects_a_body_with_no_vertical_extent() -> None:
    flat = tc.Mesh(vertices=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)), faces=((0, 1, 2),))
    with pytest.raises(RightingError):
        righting_arm_curve(flat, kg=0.1, waterline_z=-1.0, heel_angles_deg=(0.0, 5.0))


def test_wall_sided_helper_is_the_documented_closed_form() -> None:
    assert wall_sided_righting_arm(0.0, gm=0.3, bm=0.2) == pytest.approx(0.0, abs=1e-15)
    # GZ / sin(phi) = GM + BM tan^2(phi) / 2
    arm = wall_sided_righting_arm(30.0, gm=0.3, bm=0.2)
    assert arm / math.sin(math.radians(30.0)) == pytest.approx(0.3 + 0.1 * math.tan(math.radians(30.0)) ** 2)
    with pytest.raises(RightingError, match="heel_deg must be finite"):
        wall_sided_righting_arm(float("nan"), gm=0.3, bm=0.2)


def test_default_station_intervals_are_documented() -> None:
    assert DEFAULT_STATION_INTERVALS == 80

# ------------------------------------------------------------------ submerged body


def test_submerged_body_of_a_barge_is_exact() -> None:
    body = tc.submerged_body(box_barge(), heel_deg=0.0, waterline_z=0.0, station_intervals=8)
    assert body.volume == pytest.approx(LENGTH * BEAM * DRAFT, rel=1e-12)
    assert body.vertical_centre_buoyancy == pytest.approx(-DRAFT / 2.0, abs=1e-12)
    assert body.transverse_centre_buoyancy == pytest.approx(0.0, abs=1e-12)
    assert body.to_dict()["station_intervals"] == 8


def test_submerged_body_heeled_ninety_degrees_is_the_starboard_half() -> None:
    body = tc.submerged_body(box_barge(), heel_deg=90.0, waterline_z=0.0, station_intervals=8)
    # The waterplane stands vertical through the centreline: half the breadth, full depth.
    assert body.volume == pytest.approx(LENGTH * (BEAM / 2.0) * (DRAFT + FREEBOARD), rel=1e-10)
    assert body.transverse_centre_buoyancy == pytest.approx(BEAM / 4.0, rel=1e-9)


def test_submerged_body_agrees_with_the_closed_surface_integral() -> None:
    """Two independent integrations of one body: station areas vs tetrahedra."""
    mesh = hull_mesh()
    body = tc.submerged_body(mesh, heel_deg=0.0, waterline_z=0.0, station_intervals=80)
    hydrostatics = tc.compute_hydrostatics(mesh, 0.0)
    assert body.volume == pytest.approx(hydrostatics.displacement_volume, rel=1e-12)
    assert body.vertical_centre_buoyancy == pytest.approx(
        hydrostatics.vertical_centre_buoyancy, abs=1e-12
    )
    assert body.transverse_centre_buoyancy == pytest.approx(0.0, abs=1e-9)


def test_submerged_body_rejects_non_finite_input() -> None:
    with pytest.raises(RightingError, match="heel_deg must be finite"):
        tc.submerged_body(box_barge(), heel_deg=float("nan"))
    with pytest.raises(RightingError, match="waterline_z must be finite"):
        tc.submerged_body(box_barge(), waterline_z=float("inf"))


def test_submerged_body_reports_an_empty_body_without_raising() -> None:
    body = tc.submerged_body(box_barge(), heel_deg=0.0, waterline_z=-5.0, station_intervals=8)
    assert body.volume == 0.0
    assert math.isnan(body.vertical_centre_buoyancy)
    full = tc.submerged_body(box_barge(), heel_deg=0.0, waterline_z=5.0, station_intervals=8)
    assert full.volume == pytest.approx(LENGTH * BEAM * (DRAFT + FREEBOARD), rel=1e-12)
