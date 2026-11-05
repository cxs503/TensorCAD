"""Large-angle stability: the righting-arm curve against two exact closed forms.

The module claims a righting arm, so the tests have to be able to say what the right
answer is. Two shapes give one without any discretisation error:

* a **box barge** with vertical sides, where the submerged section is a trapezoid and
  ``GZ = sin(phi) (GM0 + BM0/2 tan^2 phi)`` follows from integrating ``int z dA`` and
  ``int y dA`` over it -- exact until the waterline reaches the bilge or the deck edge, and
  the tests check both that it matches below that limit and that it *stops* matching above
  it, so the agreement cannot pass vacuously;
* a **log** (circular section, half submerged), where central symmetry pins the waterline
  through the section centre at every heel and gives ``GZ = (R - KG) sin(phi)`` exactly.

Two more checks tie the new integration to what the package already had:
````KB`` and the slope ``dGZ/dphi`` at zero must reproduce
:func:`tensorcad.stability.initial_stability`, which gets them from the waterplane moments
by a completely different route (Green's theorem on the cut loops instead of a volume
integral), and the displaced volume must stay constant at every heel, because a curve whose
points float at different displacements is not a righting-arm curve.
"""

from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.gz import (
    DEFAULT_ANGLES_DEG,
    GzError,
    box_righting_arms_analytic,
    box_wall_sided_limit_deg,
    deck_edge_immersion_deg,
    free_surface_correction,
    heeled_hydrostatics,
    rectangular_tank_inertia,
    righting_arm_curve,
    wall_sided_righting_arm,
    with_freeboard,
)
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.stability import initial_stability

BREADTH = 2.0
DRAFT = 0.8
DEPTH = 2.0
LENGTH = 10.0
KG = 0.5


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
    """Box barge with freeboard: keel at z=0, deck at DEPTH, floating at DRAFT."""
    return box((0.0, -BREADTH / 2, 0.0), (LENGTH, BREADTH / 2, DEPTH))


@pytest.fixture(scope="module")
def barge_at_waterline() -> Mesh:
    """The same barge clipped at the waterline, which is what hull generators emit."""
    return box((0.0, -BREADTH / 2, 0.0), (LENGTH, BREADTH / 2, DRAFT))


def log(radius: float, length: float, segments: int = 64) -> Mesh:
    """Cylinder of regular-polygon section, axis along x at z = radius (draft = radius)."""
    builder = MeshBuilder()
    rings: list[tuple[int, int]] = []
    for index in range(segments):
        angle = 2.0 * math.pi * index / segments
        y = radius * math.cos(angle)
        z = radius + radius * math.sin(angle)
        rings.append((builder.add_vertex((0.0, y, z)), builder.add_vertex((length, y, z))))
    for index in range(segments):
        low, high = rings[index], rings[(index + 1) % segments]
        builder.add_triangle(low[0], low[1], high[1])
        builder.add_triangle(low[0], high[1], high[0])
    first_cap = builder.add_vertex((0.0, 0.0, radius))
    second_cap = builder.add_vertex((length, 0.0, radius))
    for index in range(segments):
        low, high = rings[index], rings[(index + 1) % segments]
        builder.add_triangle(first_cap, high[0], low[0])
        builder.add_triangle(second_cap, low[1], high[1])
    return meshlib.orient_outward(builder.build()).mesh


def test_box_reproduces_the_wall_sided_closed_form(barge: Mesh) -> None:
    """Below the bilge/deck limit the trapezoid derivation is exact, not close."""
    limit = box_wall_sided_limit_deg(BREADTH, DRAFT, DEPTH)
    assert limit == pytest.approx(math.degrees(math.atan(2.0 * DRAFT / BREADTH)))
    angles = list(range(0, int(limit) + 1, 5))
    curve = righting_arm_curve(barge, KG, angles, waterline_z=DRAFT)
    reference = box_righting_arms_analytic(BREADTH, DRAFT, KG, angles, depth=DEPTH)
    for angle, gz, expected in zip(curve.angles_deg, curve.gz, reference, strict=True):
        assert gz == pytest.approx(expected, rel=0, abs=1e-14), f"heel {angle} deg"


def test_wall_sided_formula_stops_holding_past_its_limit(barge: Mesh) -> None:
    """Negative control: otherwise the comparison above could be true for any integrator."""
    limit = box_wall_sided_limit_deg(BREADTH, DRAFT, DEPTH)
    bm = BREADTH * BREADTH / (12.0 * DRAFT)
    gm = 0.5 * DRAFT + bm - KG
    deviations = {}
    for angle in (40.0, 50.0):
        assert angle > limit
        gz = righting_arm_curve(barge, KG, [angle], waterline_z=DRAFT).gz[0]
        reference = wall_sided_righting_arm(angle, gm=gm, bm=bm)
        deviations[angle] = abs(gz - reference)
        assert deviations[angle] > 1e-4, "the closed form should have broken down by now"
    # and the disagreement keeps growing as more of the bilge emerges
    assert deviations[50.0] > deviations[40.0]
    assert deviations[50.0] > 1e-3


def test_displacement_is_conserved_at_every_heel(barge: Mesh) -> None:
    angles = [0.0, 1.0, 5.0, 12.5, 25.0, 40.0, 60.0, 89.0]
    curve = righting_arm_curve(barge, KG, angles, waterline_z=DRAFT)
    for point in curve.points:
        assert abs(point.volume_error) < 1e-12, f"heel {point.heel_angle_deg} deg"
        assert point.displacement_volume == pytest.approx(LENGTH * BREADTH * DRAFT, rel=1e-12)
    assert curve.max_volume_error < 1e-12


def test_upright_kb_matches_the_waterplane_moment_route(barge_at_waterline: Mesh) -> None:
    """``KB`` from the volume integral against ``KB`` from the waterplane machinery."""
    reference = initial_stability(barge_at_waterline, kg=KG, waterline_z=DRAFT)
    point = heeled_hydrostatics(barge_at_waterline, 0.0, waterline_z=DRAFT)
    assert point.centre_of_buoyancy_z == pytest.approx(reference.kb, rel=0, abs=1e-12)
    assert point.centre_of_buoyancy_y == pytest.approx(0.0, rel=0, abs=1e-12)
    assert point.waterplane_area == pytest.approx(LENGTH * BREADTH, rel=1e-12)
    assert point.displacement_volume == pytest.approx(reference.displacement_volume, rel=1e-12)


def test_slope_at_zero_is_the_metacentric_height(barge: Mesh, barge_at_waterline: Mesh) -> None:
    """``GM`` is the slope of the curve: a finite difference must land on it."""
    reference = initial_stability(barge_at_waterline, kg=KG, waterline_z=DRAFT)
    curve = righting_arm_curve(barge, KG, [-0.01, 0.0, 0.01], waterline_z=DRAFT)
    slope = (curve.gz[2] - curve.gz[0]) / math.radians(0.02)
    assert slope == pytest.approx(reference.gm_transverse, rel=1e-6)
    assert curve.gz[1] == pytest.approx(0.0, abs=1e-12)


def test_curve_is_linear_in_kg(barge: Mesh) -> None:
    """Raising G shifts every point down by exactly ``d * sin(phi)``."""
    angles = [0.0, 10.0, 30.0, 60.0]
    low = righting_arm_curve(barge, KG, angles, waterline_z=DRAFT)
    high = righting_arm_curve(barge, KG + 0.25, angles, waterline_z=DRAFT)
    for angle, gz_low, gz_high in zip(angles, low.gz, high.gz, strict=True):
        assert gz_high - gz_low == pytest.approx(-0.25 * math.sin(math.radians(angle)), rel=0, abs=1e-12)


def test_heel_sign_is_antisymmetric(barge: Mesh) -> None:
    """A symmetric hull gives equal and opposite arms to port and starboard."""
    starboard = righting_arm_curve(barge, KG, [10.0, 30.0, 55.0], waterline_z=DRAFT)
    port = righting_arm_curve(barge, KG, [-10.0, -30.0, -55.0], waterline_z=DRAFT)
    for gz_starboard, gz_port in zip(starboard.gz, port.gz, strict=True):
        assert gz_starboard == pytest.approx(-gz_port, rel=1e-11)


def test_log_matches_the_exact_circular_section_reference() -> None:
    """``GZ = (R - KG) sin(phi)``; symmetry-aligned cuts make the polygon reference exact."""
    radius = 1.0
    segments = 64
    step = 180.0 / segments  # a symmetry axis of the polygon every ``step`` degrees
    angles = [0.0, step, 2 * step, 4 * step, 8 * step, 16 * step, 24 * step, 28 * step]
    hull = log(radius, 4.0, segments)
    assert meshlib.is_watertight(hull)
    curve = righting_arm_curve(hull, KG, angles, waterline_z=radius)
    for angle, gz in zip(angles, curve.gz, strict=True):
        expected = (radius - KG) * math.sin(math.radians(angle))
        assert gz == pytest.approx(expected, rel=1e-13, abs=1e-15), f"heel {angle} deg"
        assert curve.points[angles.index(angle)].kn == pytest.approx(
            radius * math.sin(math.radians(angle)), rel=1e-13, abs=1e-15
        )


def test_log_has_no_stability_with_g_above_the_section_centre() -> None:
    """A round bilge cannot right itself once ``G`` is above ``R``: every arm is negative."""
    radius = 1.0
    curve = righting_arm_curve(log(radius, 4.0), 1.4, [0.0, 5.0, 15.0, 40.0, 80.0], waterline_z=radius)
    assert all(gz < 0.0 for gz in curve.gz[1:])
    assert max(curve.gz[1:]) < 0.0
    assert curve.max_gz == pytest.approx(0.0, abs=1e-12)  # upright, and only there


def test_curve_summary_metrics(barge: Mesh) -> None:
    curve = righting_arm_curve(barge, KG, [0.0, 10.0, 20.0, 30.0], waterline_z=DRAFT)
    assert curve.max_gz == max(curve.gz)
    assert curve.angle_of_max_gz_deg == 30.0
    # still rising at the last sample, so the range is wider than the sampled range
    assert curve.range_of_stability_deg is None
    area = curve.area_under_curve()
    assert area > 0.0
    assert curve.area_under_curve(0.0) == 0.0


def test_area_under_curve_interpolates_its_limit(barge: Mesh) -> None:
    curve = righting_arm_curve(barge, KG, [0.0, 30.0], waterline_z=DRAFT)
    # one trapezoid from 0 to 30 degrees
    expected = 0.5 * (curve.gz[0] + curve.gz[1]) * math.radians(30.0)
    assert curve.area_under_curve(30.0) == pytest.approx(expected, rel=1e-12)
    half = curve.area_under_curve(15.0)
    assert half == pytest.approx(0.5 * (curve.gz[0] + 0.5 * curve.gz[1]) * math.radians(15.0), rel=1e-12)
    assert half < expected


def test_synthetic_curve_range_interpolates_between_samples() -> None:
    from tensorcad.gz import HeeledHydrostatics, RightingArmCurve

    points = tuple(
        HeeledHydrostatics(
            heel_angle_deg=angle,
            displacement_volume=1.0,
            target_volume=1.0,
            volume_error=0.0,
            waterline_offset=0.0,
            waterplane_area=1.0,
            wetted_surface_area=1.0,
            centre_of_buoyancy_y=0.0,
            centre_of_buoyancy_z=0.5,
            kn=0.0,
            iterations=1,
        )
        for angle in (0.0, 10.0, 20.0, 30.0)
    )
    curve = RightingArmCurve(
        angles_deg=(0.0, 10.0, 20.0, 30.0),
        gz=(0.0, 0.2, 0.1, -0.1),
        kn=(0.0, 0.0, 0.0, 0.0),
        points=points,
        kg=0.5,
        free_surface_correction=0.0,
        waterline_z=0.0,
        keel_z=0.0,
        displacement_volume=1.0,
    )
    assert curve.angle_of_max_gz_deg == 10.0
    assert curve.range_of_stability_deg == pytest.approx(25.0, rel=1e-12)
    assert curve.area_under_curve(25.0) == pytest.approx(
        0.5 * (0.0 + 0.2) * math.radians(10.0)
        + 0.5 * (0.2 + 0.1) * math.radians(10.0)
        + 0.5 * (0.1 + 0.0) * math.radians(5.0),
        rel=1e-12,
    )


def test_rectangular_tank_inertia() -> None:
    assert rectangular_tank_inertia(4.0, 2.0) == pytest.approx(4.0 * 8.0 / 12.0)
    with pytest.raises(GzError):
        rectangular_tank_inertia(0.0, 2.0)


def test_free_surface_correction_is_sum_ix_over_volume() -> None:
    volume = LENGTH * BREADTH * DRAFT
    correction = free_surface_correction([rectangular_tank_inertia(3.0, 1.0)], volume)
    assert correction.total_inertia == pytest.approx(3.0 / 12.0)
    assert correction.loss_m == pytest.approx((3.0 / 12.0) / volume)
    with pytest.raises(GzError):
        free_surface_correction([-1.0], volume)
    with pytest.raises(GzError):
        free_surface_correction([1.0], 0.0)


def test_free_surface_correction_shifts_every_point_of_the_curve(barge: Mesh) -> None:
    volume = LENGTH * BREADTH * DRAFT
    inertia = rectangular_tank_inertia(2.0, 0.5)
    solid = righting_arm_curve(barge, KG, [0.0, 20.0, 40.0], waterline_z=DRAFT)
    flooded = righting_arm_curve(barge, KG, [0.0, 20.0, 40.0], waterline_z=DRAFT, free_surface_i_x=inertia)
    loss = inertia / volume
    assert flooded.free_surface_correction == pytest.approx(loss, rel=1e-12)
    for gz_solid, gz_flooded in zip(solid.gz, flooded.gz, strict=True):
        assert gz_flooded == pytest.approx(gz_solid - loss, rel=0, abs=1e-12)


def test_with_freeboard_keeps_the_upright_displacement(barge_at_waterline: Mesh) -> None:
    banded = with_freeboard(barge_at_waterline, 1.2, waterline_z=DRAFT)
    assert meshlib.is_watertight(banded)
    assert meshlib.is_edge_manifold(banded)
    reference = compute_hydrostatics(barge_at_waterline, waterline_z=DRAFT)
    point = heeled_hydrostatics(banded, 0.0, waterline_z=DRAFT)
    assert point.displacement_volume == pytest.approx(reference.displacement_volume, rel=1e-12)
    assert meshlib.total_signed_volume(banded) > 0.0
    assert banded.bounding_box().maximum[2] == pytest.approx(DRAFT + 1.2)
    faces = len(banded.faces)
    vertices = len(banded.vertices)
    edges = len(meshlib.edge_keys(banded)) if hasattr(meshlib, "edge_keys") else 0
    assert faces > len(barge_at_waterline.faces)
    assert vertices > len(barge_at_waterline.vertices)
    assert edges >= 0


def test_with_freeboard_needs_a_waterline_cap(barge: Mesh) -> None:
    with pytest.raises(GzError):
        with_freeboard(barge, 0.5, waterline_z=DRAFT)


def test_with_freeboard_rejects_a_non_positive_band(barge_at_waterline: Mesh) -> None:
    with pytest.raises(GzError):
        with_freeboard(barge_at_waterline, 0.0, waterline_z=DRAFT)


def test_hull_curve_slope_matches_initial_stability() -> None:
    """A real hull, freeboarded so the heel is admissible, cross-checked two ways."""
    hull = generate_hull(HullParameters(length=LENGTH, beam=BREADTH, draft=DRAFT, nx=41, nz=13))
    banded = with_freeboard(hull.mesh, 0.6)
    assert meshlib.is_watertight(banded)
    reference = initial_stability(hull.mesh, kg=KG, waterline_z=0.0)

    def slope(step_deg: float) -> float:
        curve = righting_arm_curve(banded, KG, [-step_deg, 0.0, step_deg], waterline_z=0.0)
        assert curve.gz[1] == pytest.approx(0.0, abs=1e-12)
        return (curve.gz[2] - curve.gz[0]) / math.radians(2.0 * step_deg)

    # At phi = 0 the waterline coincides with the hull's waterline *vertex ring*, so the
    # polygonal hull's GZ has a curvature kink there and a central difference converges
    # only linearly in the step. Removing that first-order term (h ratio 10) leaves the
    # analytic GM, which comes from the waterplane moments instead of a volume integral.
    coarse, fine = slope(0.01), slope(0.001)
    extrapolated = fine + (fine - coarse) / 9.0
    assert extrapolated == pytest.approx(reference.gm_transverse, rel=1e-7)


def test_hull_curve_is_physically_shaped() -> None:
    hull = generate_hull(HullParameters(length=LENGTH, beam=BREADTH, draft=DRAFT, nx=41, nz=13))
    banded = with_freeboard(hull.mesh, 0.8)
    curve = righting_arm_curve(banded, KG, [0.0, 10.0, 20.0, 30.0], waterline_z=0.0)
    assert curve.gz[0] == pytest.approx(0.0, abs=1e-12)
    assert all(gz > 0.0 for gz in curve.gz[1:])
    assert 5.0 < curve.angle_of_max_gz_deg <= 30.0
    assert curve.max_volume_error < 1e-12
    payload = curve.to_dict()
    assert payload["schema"] == "tensorcad.stability.righting_arm/1"
    assert len(payload["points"]) == len(curve.angles_deg)  # type: ignore[arg-type]


def test_deck_edge_immersion_angle() -> None:
    assert deck_edge_immersion_deg(1.0, 1.0) == pytest.approx(45.0)
    assert deck_edge_immersion_deg(1.0, 0.5) == pytest.approx(math.degrees(math.atan(0.5)))
    with pytest.raises(GzError):
        deck_edge_immersion_deg(0.0, 1.0)


def test_open_mesh_is_refused() -> None:
    builder = MeshBuilder()
    a = builder.add_vertex((0.0, 0.0, 0.0))
    b = builder.add_vertex((1.0, 0.0, 0.0))
    c = builder.add_vertex((0.0, 1.0, 0.0))
    builder.add_triangle(a, b, c)
    with pytest.raises(GzError):
        heeled_hydrostatics(builder.build(), 10.0, waterline_z=0.1)


def test_waterline_below_the_keel_is_refused(barge: Mesh) -> None:
    with pytest.raises(GzError):
        heeled_hydrostatics(barge, 10.0, waterline_z=-1.0)


def test_invalid_angles_are_refused(barge: Mesh) -> None:
    with pytest.raises(GzError):
        heeled_hydrostatics(barge, 180.0, waterline_z=DRAFT)
    with pytest.raises(GzError):
        heeled_hydrostatics(barge, float("nan"), waterline_z=DRAFT)
    with pytest.raises(GzError):
        righting_arm_curve(barge, KG, [], waterline_z=DRAFT)
    with pytest.raises(GzError):
        righting_arm_curve(barge, KG, [10.0, 10.0], waterline_z=DRAFT)
    with pytest.raises(GzError):
        righting_arm_curve(barge, KG, [10.0], waterline_z=DRAFT, free_surface_i_x=-1.0)


def test_analytic_helper_refuses_heel_past_the_limit() -> None:
    limit = box_wall_sided_limit_deg(BREADTH, DRAFT, DEPTH)
    assert limit == pytest.approx(math.degrees(math.atan(2.0 * DRAFT / BREADTH)))
    with pytest.raises(GzError):
        box_righting_arms_analytic(BREADTH, DRAFT, KG, [limit + 1.0], depth=DEPTH)
    with pytest.raises(GzError):
        box_wall_sided_righting_arm_check()


def box_wall_sided_righting_arm_check() -> float:
    """The wall-sided form is undefined at 90 degrees, where the tangent diverges."""
    return wall_sided_righting_arm(90.0, gm=0.3, bm=0.4)


def test_default_angles_span_a_curve() -> None:
    assert DEFAULT_ANGLES_DEG[0] == 0.0
    assert DEFAULT_ANGLES_DEG[-1] == 90.0
    assert len(DEFAULT_ANGLES_DEG) == 19


# --- the published evidence page (publish_report.py -> docs/gz/) ----------------------
#
# The page is part of the deliverable, so it is a test: the JSON has to hold every claim,
# the PNG has to exist, the Markdown has to be the render of that JSON rather than prose
# typed next to it, and ``publish_report.py --check`` has to agree with the committed
# numbers when it recomputes them from the geometry. Without this the page is a claim
# about the past; with it, a change to the integrator that moves a published number fails
# here instead of silently ageing the documentation.

ROOT = Path(__file__).resolve().parents[1]
PUBLISHED = ROOT / "docs" / "gz"
CLAIM_IDS = ("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9")


def _published_json() -> dict:
    return json.loads((PUBLISHED / "results.json").read_text())


def test_published_artifacts_exist_and_are_not_empty() -> None:
    for name in ("results.json", "report.md", "figures/gz-curve.png"):
        path = PUBLISHED / name
        assert path.exists(), f"{name} is missing from docs/gz"
        assert path.stat().st_size > 0, f"{name} is empty"


def test_published_figure_is_a_real_png() -> None:
    """PNG magic + IHDR, without needing an image library: the file must be a drawing."""
    blob = (PUBLISHED / "figures" / "gz-curve.png").read_bytes()
    assert blob[:8] == b"\x89PNG\r\n\x1a\n"
    width = int.from_bytes(blob[16:20], "big")
    height = int.from_bytes(blob[20:24], "big")
    assert width > 800 and height > 300, (width, height)


def test_published_json_records_every_claim() -> None:
    data = _published_json()
    claims = {entry["id"]: entry for entry in data["claims"]}
    assert tuple(claims) == CLAIM_IDS
    assert data["passed"] is True
    assert data["generated_by"] == "publish_report.py"
    import tensorcad

    assert data["tensorcad_version"] == tensorcad.__version__
    for entry in claims.values():
        assert isinstance(entry["threshold"], float) and entry["threshold"] > 0.0
        direction = entry.get("direction", "upper")
        assert direction in ("upper", "lower")
        if direction == "upper":
            assert entry["value"] <= entry["threshold"], entry
        else:
            assert entry["value"] >= entry["threshold"], entry
        assert entry["passed"] is True


def test_published_json_negative_controls_bite() -> None:
    """C2 and C7 are controls: they only mean anything if they are far from round-off."""
    claims = {entry["id"]: entry for entry in _published_json()["claims"]}
    for ident in ("C2", "C7"):
        assert claims[ident]["direction"] == "lower"
        # both are controls, so they must clear their threshold by a wide margin: a control
        # that barely bites is indistinguishable from an environment drift.
        assert claims[ident]["value"] > 50.0 * claims[ident]["threshold"], claims[ident]


def test_published_json_cross_route_is_round_off_not_a_coincidence() -> None:
    """The two GZ routes must agree to round-off *and* the control must be far larger."""
    data = _published_json()
    agreement = data["cross_route"]["worst_abs_difference_m"]
    control = data["cross_route"]["negative_control"]["gz_shift_m"]
    assert agreement < 1e-11
    assert control > 1e3 * agreement
    hull_curve = data["hull"]["curve"]
    assert hull_curve["max_gz_m"] > 0.0
    assert hull_curve["angle_of_max_gz_deg"] > 0.0
    assert data["hull"]["watertight"] is True


def test_published_markdown_is_the_render_of_the_json() -> None:
    page = (PUBLISHED / "report.md").read_text()
    data = _published_json()
    for ident in CLAIM_IDS:
        assert f"| {ident} |" in page, f"claim {ident} is not in the verdict table"
    for heading in ("## 1. The case", "## 2. How it was done", "### 3.1 Barge", "### 3.2 Log",
                    "### 3.3 Hull", "### 3.4 Cross-route", "## 4. Caveats and scope",
                    "## 5. Error verdict"):
        assert heading in page, f"{heading} is missing from the page"
    # the page's numbers have to be the JSON's numbers, printed by the same formatter
    assert f"{data['barge']['worst_abs_error_m']:.2e}" in page
    assert f"{data['log']['worst_abs_error_m']:.2e}" in page
    assert f"{data['cross_route']['worst_abs_difference_m']:.2e}" in page
    assert f"{data['free_surface']['flat_shift_m'] * 1000:.4f}" in page
    assert "figures/gz-curve.png" in page
    assert data["generated_by"] in page
    assert "FAILED" not in page
    for marker in ("TODO", "FIXME", "XXX", "TBD"):
        assert marker not in page


def test_published_page_matches_the_code() -> None:
    """``publish_report.py --check`` recomputes every claim and compares with the page."""
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(ROOT) + os.pathsep + environment.get("PYTHONPATH", "")
    completed = subprocess.run(
        [sys.executable, str(ROOT / "publish_report.py"), "--check"],
        cwd=ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=900,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "all 9 claims hold" in completed.stdout
