"""Tests for large-angle stability (:mod:`tensorcad.gz`).

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

from tensorcad import HullParameters, MeshBuilder, generate_hull
from tensorcad.gz import (
    DEFAULT_HEEL_ANGLES_DEG,
    GZError,
    gz_curve,
    heel_mesh,
    imo_intact_stability,
    solve_waterline,
    submerged_solid,
)
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.mesh import is_watertight, total_signed_volume
from tensorcad.stability import initial_stability


# ---------------------------------------------------------------------------
# reference geometries
# ---------------------------------------------------------------------------
def box(length: float, beam: float, depth: float):
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
def barge():
    """Wall-sided barge: ``L=30, B=6, T=2, D=10`` -> wall-sided up to ``atan(2T/B)=33.7 deg``."""
    return box(length=30.0, beam=6.0, depth=10.0)


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
def test_box_submerged_solid_is_exact(barge) -> None:
    solid = submerged_solid(barge, 2.0)
    assert solid.volume == pytest.approx(30.0 * 6.0 * 2.0, rel=1e-15)
    assert solid.centroid == pytest.approx((15.0, 0.0, 1.0), abs=1e-13)
    assert solid.waterplane_area == pytest.approx(30.0 * 6.0, rel=1e-15)
    assert solid.section_segments > 0


def test_whole_mesh_volume_is_recovered(barge) -> None:
    """A waterline above the deck gives back the mesh's own signed volume."""
    solid = submerged_solid(barge, 12.0)
    assert solid.volume == pytest.approx(abs(total_signed_volume(barge)), rel=1e-15)


def test_waterline_below_the_keel_is_empty(barge) -> None:
    solid = submerged_solid(barge, -1.0)
    assert solid.volume == 0.0
    assert solid.section_segments == 0


def test_solve_waterline_recovers_the_design_draft(barge) -> None:
    z, solid = solve_waterline(barge, 30.0 * 6.0 * 2.0)
    assert z == pytest.approx(2.0, abs=1e-10)
    assert solid.volume == pytest.approx(360.0, rel=1e-13)


def test_solve_waterline_rejects_a_displacement_it_cannot_float(barge) -> None:
    with pytest.raises(GZError, match="not below the mesh's own volume"):
        solve_waterline(barge, 30.0 * 6.0 * 10.0 * 1.01)
    with pytest.raises(GZError, match="not positive"):
        solve_waterline(barge, 0.0)


# ---------------------------------------------------------------------------
# the wall-sided box: exactness of the whole chain
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("kg", [0.0, 1.0, 2.0])
@pytest.mark.parametrize("heel_deg", [0.0, 10.0, 20.0, 30.0])
def test_box_gz_matches_the_wall_sided_formula(barge, kg: float, heel_deg: float) -> None:
    curve = gz_curve(barge, kg=kg, heel_angles_deg=(heel_deg,), waterline_z=2.0)
    expected = wall_sided_gz(6.0, 2.0, kg, heel_deg)
    assert curve.points[0].gz == pytest.approx(expected, abs=1e-9)


def test_box_kn_slope_reproduces_the_metacentric_height(barge) -> None:
    """``dKN/dphi`` at zero heel is ``KM = KB + BM = T/2 + B^2/(12 T)`` for a box."""
    curve = gz_curve(barge, kg=0.0, heel_angles_deg=(0.0, 0.5), waterline_z=2.0)
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


def test_displacement_is_conserved_at_every_heel_angle(barge) -> None:
    curve = gz_curve(barge, kg=1.0, heel_angles_deg=DEFAULT_HEEL_ANGLES_DEG, waterline_z=2.0)
    for point in curve.points:
        assert point.displacement_volume == pytest.approx(curve.displacement_volume, rel=1e-11)


def test_upright_gz_is_zero_and_the_curve_is_antisymmetric(barge) -> None:
    curve = gz_curve(
        barge, kg=1.0, heel_angles_deg=(-30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0), waterline_z=2.0
    )
    upright = curve.gz_at(0.0)
    assert upright == pytest.approx(0.0, abs=1e-12)
    for angle in (10.0, 20.0, 30.0):
        assert curve.gz_at(angle) == pytest.approx(-curve.gz_at(-angle), abs=1e-9)


def test_heel_mesh_moves_the_starboard_side_down(barge) -> None:
    """Positive heel = starboard down: the sign convention ``GZ = KN - KG sin(phi)`` rests on it."""
    keel = min(vertex[2] for vertex in barge.vertices)
    heeled = heel_mesh(barge, 30.0)
    by_corner = {}
    for vertex in barge.vertices:
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
def test_imo_criteria_pass_for_a_stiff_barge(barge) -> None:
    curve = gz_curve(barge, kg=1.0, waterline_z=2.0)
    report = imo_intact_stability(curve)
    assert report.passed, report.to_dict()
    assert report.gm_upright == pytest.approx(1.5, rel=1e-9)


def test_imo_criteria_fail_when_the_centre_of_gravity_is_too_high(barge) -> None:
    curve = gz_curve(barge, kg=2.45, waterline_z=2.0)
    report = imo_intact_stability(curve)
    assert not report.passed
    failing = {check.key for check in report.checks if not check.satisfied}
    assert "gm_upright" in failing
    assert report.gm_upright == pytest.approx(0.05, rel=1e-6)


def test_imo_criteria_need_a_curve_that_reaches_forty_degrees(barge) -> None:
    curve = gz_curve(barge, kg=1.0, heel_angles_deg=(0.0, 10.0, 20.0, 30.0), waterline_z=2.0)
    with pytest.raises(GZError, match="0 to at least 40 deg"):
        imo_intact_stability(curve)


def test_area_under_a_linear_curve_is_exact() -> None:
    """``area`` must integrate the piecewise-linear curve, not approximate it."""
    barge = box(length=10.0, beam=2.0, depth=4.0)
    curve = gz_curve(barge, kg=0.5, heel_angles_deg=(0.0, 10.0, 20.0, 30.0, 40.0), waterline_z=1.0)
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
def test_kg_and_angle_grid_are_validated(barge) -> None:
    with pytest.raises(GZError, match="kg must be finite"):
        gz_curve(barge, kg=float("nan"), waterline_z=2.0)
    with pytest.raises(GZError, match="must be increasing"):
        gz_curve(barge, kg=1.0, heel_angles_deg=(20.0, 10.0), waterline_z=2.0)
    with pytest.raises(GZError, match="outside the sampled range"):
        gz_curve(barge, kg=1.0, heel_angles_deg=(0.0, 10.0), waterline_z=2.0).gz_at(20.0)


def test_curve_serialises_its_evidence(barge) -> None:
    curve = gz_curve(barge, kg=1.0, waterline_z=2.0)
    payload = curve.to_dict()
    assert payload["schema"] == "tensorcad.gz.curve/1"
    assert payload["displacement_volume_m3"] == pytest.approx(360.0, rel=1e-12)
    assert payload["points"][0]["heel_angle_deg"] == 0.0
    assert payload["area_0_30_m_rad"] > 0.0
    keys = set(payload["points"][0])
    assert {"kn_m", "gz_m", "waterline_z_m", "centre_buoyancy_m"} <= keys
