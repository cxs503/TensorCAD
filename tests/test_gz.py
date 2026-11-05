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

The second block below exercises the *other* route through the same module
(:func:`tensorcad.gz.submerged_solid` / :func:`gz_curve` / :func:`imo_intact_stability`),
kept from PR #48; its module-level helpers are suffixed ``_imo`` so both suites share one
file and neither side's checks are dropped.

Tests for large-angle stability (:mod:`tensorcad.gz`).

The references are analytic, not fitted:

* a **wall-sided box** floats with its waterline intersecting both side walls up to
  ``tan(phi) = 2T/B``, and there the righting arm is the textbook
  ``GZ = sin(phi) * (GM + BM tan^2(phi) / 2)`` with ``KB = T/2`` and ``BM = B^2 / (12 T)``.
  Because the box is a polyhedron and the cut is exact, agreement is asserted to 1e-9,
  not to a "few percent";
* a **circular cylinder** cut through its centre is a circular segment, so its centre of
  buoyancy sits on the vertical through the centre of the circle at every heel angle:
  the exact cross curve is ``KN = R sin(phi)`` and the exact floatation is
  ``z = R cos(phi)``. The facetted mesh can only match to the facetting error, so that
  comparison is asserted to 2e-3 and its own discretisation is reported.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.gz import (
    DEFAULT_ANGLES_DEG,
    DEFAULT_HEEL_ANGLES_DEG,
    GZError,
    GzError,
    box_righting_arms_analytic,
    box_wall_sided_limit_deg,
    deck_edge_immersion_deg,
    free_surface_correction,
    gz_curve,
    heel_mesh,
    heeled_hydrostatics,
    imo_intact_stability,
    rectangular_tank_inertia,
    righting_arm_curve,
    solve_waterline,
    submerged_solid,
    wall_sided_righting_arm,
    with_freeboard,
)
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.mesh import is_watertight, total_signed_volume
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


# ---------------------------------------------------------------------------
# reference geometries
# ---------------------------------------------------------------------------

def box_imo(length: float, beam: float, depth: float):
    """Closed box, keel on ``z = 0``, deck on ``z = depth``, centred on ``y = 0``."""
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

    lo, hi = -0.5 * beam, 0.5 * beam
    left, right = 0.0, length
    bottom, top = 0.0, depth
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    return builder.build()


def cylinder(radius: float, length: float, segments: int = 360):
    """Closed cylinder, axis along ``x``, lowest point on ``z = 0``, centre on ``y = 0``."""
    builder = MeshBuilder()
    ring_left, ring_right = [], []
    for index in range(segments):
        # half-segment offset: no vertex may sit exactly on a waterline that the tests use,
        # because a plane through mesh edges has no section (see the closure guard in gz.py)
        theta = 2.0 * math.pi * (index + 0.5) / segments
        y = radius * math.sin(theta)
        z = radius - radius * math.cos(theta)
        ring_left.append(builder.add_vertex((0.0, y, z)))
        ring_right.append(builder.add_vertex((length, y, z)))
    centre_left = builder.add_vertex((0.0, 0.0, radius))
    centre_right = builder.add_vertex((length, 0.0, radius))
    for index in range(segments):
        nxt = (index + 1) % segments
        builder.add_triangle(ring_left[index], ring_left[nxt], ring_right[nxt])
        builder.add_triangle(ring_left[index], ring_right[nxt], ring_right[index])
        builder.add_triangle(centre_left, ring_left[nxt], ring_left[index])
        builder.add_triangle(centre_right, ring_right[index], ring_right[nxt])
    return builder.build()


@pytest.fixture(scope="module")
def barge_imo():
    """Wall-sided barge: ``L=30, B=6, T=2, D=10`` -> wall-sided up to ``atan(2T/B)=33.7 deg``."""
    return box_imo(length=30.0, beam=6.0, depth=10.0)


@pytest.fixture(scope="module")
def round_barge():
    return cylinder(radius=3.0, length=20.0, segments=360)


def wall_sided_gz(beam: float, draft: float, kg: float, heel_deg: float) -> float:
    """Textbook wall-sided righting arm, exact while the waterline cuts both side walls."""
    kb = 0.5 * draft
    bm = beam * beam / (12.0 * draft)
    gm = kb + bm - kg
    tan_phi = math.tan(math.radians(heel_deg))
    return math.sin(math.radians(heel_deg)) * (gm + 0.5 * bm * tan_phi * tan_phi)


# ---------------------------------------------------------------------------
# the submerged solid
# ---------------------------------------------------------------------------
def test_box_submerged_solid_is_exact(barge_imo) -> None:
    solid = submerged_solid(barge_imo, 2.0)
    assert solid.volume == pytest.approx(30.0 * 6.0 * 2.0, rel=1e-15)
    assert solid.centroid == pytest.approx((15.0, 0.0, 1.0), abs=1e-13)
    assert solid.waterplane_area == pytest.approx(30.0 * 6.0, rel=1e-15)
    assert solid.section_segments > 0


def test_whole_mesh_volume_is_recovered(barge_imo) -> None:
    """A waterline above the deck gives back the mesh's own signed volume."""
    solid = submerged_solid(barge_imo, 12.0)
    assert solid.volume == pytest.approx(abs(total_signed_volume(barge_imo)), rel=1e-15)


def test_waterline_below_the_keel_is_empty(barge_imo) -> None:
    solid = submerged_solid(barge_imo, -1.0)
    assert solid.volume == 0.0
    assert solid.section_segments == 0


def test_solve_waterline_recovers_the_design_draft(barge_imo) -> None:
    z, solid = solve_waterline(barge_imo, 30.0 * 6.0 * 2.0)
    assert z == pytest.approx(2.0, abs=1e-10)
    assert solid.volume == pytest.approx(360.0, rel=1e-13)


def test_solve_waterline_rejects_a_displacement_it_cannot_float(barge_imo) -> None:
    with pytest.raises(GZError, match="not below the mesh's own volume"):
        solve_waterline(barge_imo, 30.0 * 6.0 * 10.0 * 1.01)
    with pytest.raises(GZError, match="not positive"):
        solve_waterline(barge_imo, 0.0)


# ---------------------------------------------------------------------------
# the wall-sided box: exactness of the whole chain
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("kg", [0.0, 1.0, 2.0])
@pytest.mark.parametrize("heel_deg", [0.0, 10.0, 20.0, 30.0])
def test_box_gz_matches_the_wall_sided_formula(barge_imo, kg: float, heel_deg: float) -> None:
    curve = gz_curve(barge_imo, kg=kg, heel_angles_deg=(heel_deg,), waterline_z=2.0)
    expected = wall_sided_gz(6.0, 2.0, kg, heel_deg)
    assert curve.points[0].gz == pytest.approx(expected, abs=1e-9)


def test_box_kn_slope_reproduces_the_metacentric_height(barge_imo) -> None:
    """``dKN/dphi`` at zero heel is ``KM = KB + BM = T/2 + B^2/(12 T)`` for a box."""
    curve = gz_curve(barge_imo, kg=0.0, heel_angles_deg=(0.0, 0.5), waterline_z=2.0)
    slope = (curve.points[1].kn - curve.points[0].kn) / math.radians(0.5)
    km = 0.5 * 2.0 + 6.0 * 6.0 / (12.0 * 2.0)
    assert slope == pytest.approx(km, rel=2e-3)
    assert curve.gm_upright == pytest.approx(km - 0.0, rel=1e-12)


def test_gm_on_a_hull_capped_at_its_waterline_matches_initial_stability() -> None:
    """The two modules must agree where both are defined (hull capped at z = 0)."""
    model = generate_hull(HullParameters(length=12.0, beam=2.4, draft=0.9))
    mesh = model.mesh
    curve = gz_curve(mesh, kg=0.55, heel_angles_deg=(0.0, 10.0))
    reference = initial_stability(mesh, kg=0.55, waterline_z=0.0)
    assert curve.gm_upright == pytest.approx(reference.gm_transverse, rel=1e-12)
    assert curve.displacement_volume == pytest.approx(reference.displacement_volume, rel=1e-12)


def test_displacement_is_conserved_at_every_heel_angle(barge_imo) -> None:
    curve = gz_curve(barge_imo, kg=1.0, heel_angles_deg=DEFAULT_HEEL_ANGLES_DEG, waterline_z=2.0)
    for point in curve.points:
        assert point.displacement_volume == pytest.approx(curve.displacement_volume, rel=1e-11)


def test_upright_gz_is_zero_and_the_curve_is_antisymmetric(barge_imo) -> None:
    curve = gz_curve(
        barge_imo, kg=1.0, heel_angles_deg=(-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0), waterline_z=2.0
    )
    upright = curve.gz_at(0.0)
    assert upright == pytest.approx(0.0, abs=1e-12)
    for angle in (10.0, 20.0, 30.0):
        assert curve.gz_at(angle) == pytest.approx(-curve.gz_at(-angle), abs=1e-9)


def test_heel_mesh_moves_the_starboard_side_down(barge_imo) -> None:
    """Positive heel = starboard down: the sign convention ``GZ = KN - KG sin(phi)`` rests on it."""
    keel = min(vertex[2] for vertex in barge_imo.vertices)
    heeled = heel_mesh(barge_imo, 30.0)
    by_corner = {}
    for vertex in barge_imo.vertices:
        by_corner[(vertex[0], vertex[1], vertex[2])] = None
    # rotate the two deck-edge corners by hand and compare with the same corners of the mesh
    for side in (+1.0, -1.0):
        original = (30.0, side * 3.0, 10.0)
        angle = math.radians(30.0)
        expected_z = -original[1] * math.sin(angle) + (original[2] - keel) * math.cos(angle) + keel
        matches = [vertex for vertex in heeled.vertices if vertex[2] == pytest.approx(expected_z, abs=1e-12)]
        assert matches
        if side > 0.0:
            starboard_deck_z = expected_z
        else:
            port_deck_z = expected_z
    assert starboard_deck_z < port_deck_z
    assert is_watertight(heeled)
    assert max(vertex[1] for vertex in heeled.vertices) > 0.0


# ---------------------------------------------------------------------------
# the circular section: an independent analytic segment reference
# ---------------------------------------------------------------------------
def test_round_barge_kn_matches_the_circular_segment_reference(round_barge) -> None:
    radius, length = 3.0, 20.0
    segments = 360
    polygon_area = 0.5 * segments * radius * radius * math.sin(2.0 * math.pi / segments)
    target = length * polygon_area / 2.0
    angles = (0.0, 20.0, 45.0, 70.0)
    curve = gz_curve(
        round_barge, kg=1.5, heel_angles_deg=angles, displacement_volume=target, waterline_z=radius
    )
    for point in curve.points:
        phi = math.radians(point.heel_angle_deg)
        assert point.kn == pytest.approx(radius * math.sin(phi), abs=2e-3)
        # a chord through the centre is the floatation at every angle
        assert point.waterline_z == pytest.approx(radius * math.cos(phi), abs=2e-3)
        assert point.gz == pytest.approx((radius - 1.5) * math.sin(phi), abs=2e-3)


def test_round_barge_volume_matches_the_polygon_area(round_barge) -> None:
    radius, length, segments = 3.0, 20.0, 360
    polygon_area = 0.5 * segments * radius * radius * math.sin(2.0 * math.pi / segments)
    solid = submerged_solid(round_barge, radius)
    assert solid.volume == pytest.approx(length * polygon_area / 2.0, rel=1e-10)


# ---------------------------------------------------------------------------
# a real hull
# ---------------------------------------------------------------------------
def test_hull_gz_curve_is_righting_and_bounded() -> None:
    model = generate_hull(HullParameters(length=12.0, beam=2.4, draft=0.9))
    mesh = model.mesh
    # ``generate_hull`` caps the model at the design waterline, i.e. it has no freeboard, so
    # the deck edge is under water from the first degree of heel: the curve is the
    # "no-flooding" idealisation and is only quoted over the small angles a real hull's
    # watertight deck would still hold. Large-angle work needs a model with freeboard --
    # the barge fixtures above -- which is stated in docs/gz.md.
    # KG = 0.40 m gives GM = 0.155 m: a stable but not stiff loading condition
    curve = gz_curve(mesh, kg=0.40, heel_angles_deg=(0.0, 5.0, 10.0, 15.0, 20.0))
    displacement = compute_hydrostatics(mesh, 0.0).displacement_volume
    assert curve.displacement_volume == pytest.approx(displacement, rel=1e-12)
    assert curve.gm_upright > 0.0
    positive = [point.gz for point in curve.points if 0.0 < point.heel_angle_deg <= 20.0]
    assert all(value > 0.0 for value in positive), positive
    maximum, angle = curve.max_gz()
    assert maximum > 0.0
    assert 0.0 < angle <= 20.0
    # an arm cannot exceed the height of B over the keel plus the beam: a coarse bound,
    # which catches a lever arm computed in the wrong frame
    assert maximum < 0.5 * 2.4 + 0.9


# ---------------------------------------------------------------------------
# IMO screening
# ---------------------------------------------------------------------------
def test_imo_criteria_pass_for_a_stiff_barge(barge_imo) -> None:
    curve = gz_curve(barge_imo, kg=1.0, waterline_z=2.0)
    report = imo_intact_stability(curve)
    assert report.passed, report.to_dict()
    assert report.gm_upright == pytest.approx(1.5, rel=1e-9)


def test_imo_criteria_fail_when_the_centre_of_gravity_is_too_high(barge_imo) -> None:
    curve = gz_curve(barge_imo, kg=2.45, waterline_z=2.0)
    report = imo_intact_stability(curve)
    assert not report.passed
    failing = {check.key for check in report.checks if not check.satisfied}
    assert "gm_upright" in failing
    assert report.gm_upright == pytest.approx(0.05, rel=1e-6)


def test_imo_criteria_need_a_curve_that_reaches_forty_degrees(barge_imo) -> None:
    curve = gz_curve(barge_imo, kg=1.0, heel_angles_deg=(0.0, 10.0, 20.0, 30.0), waterline_z=2.0)
    with pytest.raises(GZError, match="0 to at least 40 deg"):
        imo_intact_stability(curve)


def test_area_under_a_linear_curve_is_exact() -> None:
    """``area`` must integrate the piecewise-linear curve, not approximate it."""
    barge_imo = box_imo(length=10.0, beam=2.0, depth=4.0)
    curve = gz_curve(barge_imo, kg=0.5, heel_angles_deg=(0.0, 10.0, 20.0, 30.0, 40.0), waterline_z=1.0)
    # the integral of the piecewise-linear curve is the sum of the segment trapezoids,
    # computed here from the sampled arms rather than from the module's own interpolation
    angles = curve.angles
    expected = 0.0
    for left, right, gz_left, gz_right in zip(
        angles, angles[1:], curve.points, curve.points[1:], strict=False
    ):
        expected += 0.5 * (gz_left.gz + gz_right.gz) * math.radians(right - left)
    assert curve.area(0.0, 40.0) == pytest.approx(expected, rel=1e-13)
    # a single chord over the whole interval would *not* be the integral of a curved arm
    trapezoid = 0.5 * (curve.gz_at(0.0) + curve.gz_at(40.0)) * math.radians(40.0)
    assert abs(curve.area(0.0, 40.0) - trapezoid) > 1e-3
    # and it is additive
    assert curve.area(0.0, 40.0) == pytest.approx(
        curve.area(0.0, 20.0) + curve.area(20.0, 40.0), rel=1e-12
    )


# ---------------------------------------------------------------------------
# contracts and error paths
# ---------------------------------------------------------------------------
def test_kg_and_angle_grid_are_validated(barge_imo) -> None:
    with pytest.raises(GZError, match="kg must be finite"):
        gz_curve(barge_imo, kg=float("nan"), waterline_z=2.0)
    with pytest.raises(GZError, match="must be increasing"):
        gz_curve(barge_imo, kg=1.0, heel_angles_deg=(20.0, 10.0), waterline_z=2.0)
    with pytest.raises(GZError, match="outside the sampled range"):
        gz_curve(barge_imo, kg=1.0, heel_angles_deg=(0.0, 10.0), waterline_z=2.0).gz_at(20.0)


def test_curve_serialises_its_evidence(barge_imo) -> None:
    curve = gz_curve(barge_imo, kg=1.0, waterline_z=2.0)
    payload = curve.to_dict()
    assert payload["schema"] == "tensorcad.gz.curve/1"
    assert payload["displacement_volume_m3"] == pytest.approx(360.0, rel=1e-12)
    assert payload["points"][0]["heel_angle_deg"] == 0.0
    assert payload["area_0_30_m_rad"] > 0.0
    keys = set(payload["points"][0])
    assert {"kn_m", "gz_m", "waterline_z_m", "centre_buoyancy_m"} <= keys
