"""Large-angle stability: ``GZ`` against a closed form that is exact for the body tested.

The reference is the wall-sided barge, for which

    GZ(phi) = sin(phi) * (GM + BM tan^2(phi) / 2)

is not an approximation but the answer, at every heel angle up to deck-edge immersion. The
mesh of a parallelepiped is the parallelepiped, so the only thing under test is the
machinery: the clip, the divergence-theorem integrals, the constant-displacement waterline
solve and the arm itself.

The weaker but more interesting check is on the *curved* Wigley hull, where the mesh is not
the geometry: there the curve has to converge as the mesh refines, and its initial slope has
to reproduce ``GM`` assembled from an independent waterplane integral.
"""

from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import MeshBuilder
from tensorcad.heeling import (
    HeelError,
    barge_mesh,
    cross_curve,
    cross_curves,
    gz_curve,
    gz_initial_slope,
    heel_normal,
    submerged_properties,
    wall_sided_righting_arm,
    waterline_offset,
)
from tensorcad.stability import waterplane_moments
from tensorcad.wigley import generate_wigley, reference_parameters

LENGTH, BEAM, DEPTH = 10.0, 2.0, 1.5
DRAFT = 0.8
KG = 0.5
# Closed-form hydrostatics of the barge at the design draught.
KB = DRAFT / 2.0
BM = BEAM**2 / (12.0 * DRAFT)
GM = KB + BM - KG
DISPLACEMENT = LENGTH * BEAM * DRAFT
# The wall-sided form stops being exact when the waterline reaches the deck edge on the
# high side: tan(phi) = (DEPTH - DRAFT) / (BEAM / 2).
DECK_EDGE_DEGREES = math.degrees(math.atan((DEPTH - DRAFT) / (BEAM / 2.0)))


@pytest.fixture(scope="module")
def barge():
    return barge_mesh(LENGTH, BEAM, DEPTH)


@pytest.fixture(scope="module")
def wigley_mesh():
    # Floated at 70 % of its design draught: the mesh is capped at the design waterline, so
    # it carries 30 % of a draught as freeboard and admits a heeled equilibrium.
    parameters = reference_parameters()
    return generate_wigley(parameters).mesh, -0.3 * parameters.draft


# ---------------------------------------------------------------------------
# The submerged-solid kernel
# ---------------------------------------------------------------------------


def test_barge_mesh_is_a_closed_outward_solid(barge):
    assert meshlib.is_watertight(barge)
    assert meshlib.total_signed_volume(barge) == pytest.approx(LENGTH * BEAM * DEPTH, rel=1e-15)
    assert meshlib.boundary_edges(barge) == ()


def test_upright_displacement_matches_closed_form(barge):
    volume, centroid, area = submerged_properties(barge, (0.0, 0.0, 1.0), DRAFT)
    assert volume == pytest.approx(DISPLACEMENT, rel=1e-15)
    assert centroid[0] == pytest.approx(0.0, abs=1e-15)
    assert centroid[1] == pytest.approx(0.0, abs=1e-15)
    assert centroid[2] == pytest.approx(KB, rel=1e-15)
    assert area == pytest.approx(LENGTH * BEAM, rel=1e-14)


def test_waterplane_area_agrees_with_the_green_theorem_waterplane(barge, wigley_mesh):
    """Two independent integrations of the same waterplane: clipping vs loop chaining."""
    assert submerged_properties(barge, (0.0, 0.0, 1.0), DRAFT)[2] == pytest.approx(
        waterplane_moments(barge, DRAFT).area, rel=1e-14
    )
    mesh, waterline = wigley_mesh
    assert submerged_properties(mesh, (0.0, 0.0, 1.0), waterline)[2] == pytest.approx(
        waterplane_moments(mesh, waterline).area, rel=1e-12
    )


def test_a_45_degree_plane_through_the_centre_of_a_box_halves_it(barge):
    """A slice at 45 degrees: the volume is half by symmetry, the centroid by hand.

    In the ``(y, z)`` section the plane ``z = y + DEPTH / 2`` cuts the rectangle out into a
    quadrilateral with corners ``(-0.75, 0)``, ``(1, 0)``, ``(1, 1.5)``, ``(0.75, 1.5)``,
    which is a triangle of area ``1.3125`` with centroid ``(5/12, 1/2)`` plus one of area
    ``0.1875`` with centroid ``(1/3, 1)``; that gives ``c = (0.40625, 0.5625)`` and the
    length is untouched. Volume is half of the body by symmetry, so both numbers are
    exact and neither comes from this module.
    """
    half = math.sqrt(0.5)
    volume, centroid, _ = submerged_properties(barge, (0.0, -half, half), 0.5 * DEPTH * half)
    assert volume == pytest.approx(0.5 * LENGTH * BEAM * DEPTH, rel=1e-14)
    assert centroid[0] == pytest.approx(0.0, abs=1e-14)
    assert centroid[1] == pytest.approx(0.40625, rel=1e-13)
    assert centroid[2] == pytest.approx(0.5625, rel=1e-13)


@pytest.mark.parametrize("dividing, expected", [(0.25, 0.375), (0.5, 0.75), (0.75, 1.125)])
def test_axis_aligned_planes_give_exact_volumes_and_centroids(barge, dividing, expected):
    """Sweeping the plane up the depth: ``V = L B z``, ``c = (0, 0, z/2)``."""
    offset = dividing * DEPTH
    volume, centroid, area = submerged_properties(barge, (0.0, 0.0, 1.0), offset)
    assert volume == pytest.approx(LENGTH * BEAM * offset, rel=1e-14)
    assert centroid[0] == pytest.approx(0.0, abs=1e-14)
    assert centroid[1] == pytest.approx(0.0, abs=1e-14)
    assert centroid[2] == pytest.approx(0.5 * offset, rel=1e-13)
    assert area == pytest.approx(LENGTH * BEAM, rel=1e-14)


def test_waterline_offset_recovers_the_requested_volume(barge):
    for fraction in (0.1, 0.25, 0.5, 0.75, 0.9):
        target = fraction * LENGTH * BEAM * DEPTH
        offset = waterline_offset(barge, (0.0, 0.0, 1.0), target)
        volume, _, _ = submerged_properties(barge, (0.0, 0.0, 1.0), offset)
        assert volume == pytest.approx(target, rel=1e-12)
        # For a box the offset is the draught itself, exactly.
        assert offset == pytest.approx(fraction * DEPTH, rel=1e-10)


def test_waterline_offset_refuses_a_body_with_no_freeboard(barge):
    """A target volume equal to the solid leaves no heeled equilibrium; say so."""
    solid = LENGTH * BEAM * DEPTH
    with pytest.raises(HeelError, match="no freeboard"):
        waterline_offset(barge, (0.0, 0.0, 1.0), solid)


def test_a_hull_capped_at_its_waterline_cannot_heel(barge):
    """The design-waterline-capped Wigley displaces its whole solid: no GZ curve exists."""
    mesh = generate_wigley(reference_parameters()).mesh
    volume = meshlib.total_signed_volume(mesh)
    with pytest.raises(HeelError, match="no freeboard"):
        gz_curve(mesh, (0.0, 5.0), kg=1.0, displacement_volume=volume)


def test_open_surface_is_rejected():
    builder = MeshBuilder()
    a, b, c = builder.add_vertex((0.0, 0.0, 0.0)), builder.add_vertex((1.0, 0.0, 0.0)), builder.add_vertex((0.0, 1.0, 0.0))
    builder.add_triangle(a, b, c)
    with pytest.raises(HeelError, match="closed solid"):
        submerged_properties(builder.build(), (0.0, 0.0, 1.0), 0.0)


# ---------------------------------------------------------------------------
# The waterline normal
# ---------------------------------------------------------------------------


def test_heel_normal_tilts_towards_the_low_side():
    assert heel_normal(0.0) == (0.0, 0.0, 1.0)
    normal = heel_normal(10.0)
    assert normal[0] == 0.0
    assert normal[1] == pytest.approx(-math.sin(math.radians(10.0)))
    assert normal[2] == pytest.approx(math.cos(math.radians(10.0)))
    # At 90 degrees of starboard heel the starboard side is straight down.
    assert heel_normal(90.0) == pytest.approx((0.0, -1.0, 0.0), abs=1e-15)


def test_heel_normal_rejects_impossible_angles():
    with pytest.raises(HeelError, match=r"\[-90, 90\]"):
        heel_normal(95.0)
    with pytest.raises(HeelError, match="finite"):
        heel_normal(float("nan"))


# ---------------------------------------------------------------------------
# GZ against the exact wall-sided closed form
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("heel", [2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0])
def test_gz_matches_the_wall_sided_closed_form(barge, heel):
    curve = gz_curve(barge, (heel,), kg=KG, displacement_volume=DISPLACEMENT)
    expected = wall_sided_righting_arm(math.radians(heel), kb=KB, bm=BM, kg=KG)
    assert curve.gz[0] == pytest.approx(expected, rel=1e-12, abs=1e-15)


def test_wall_sided_form_is_exposed_and_self_consistent():
    assert wall_sided_righting_arm(0.0, kb=KB, bm=BM, kg=KG) == 0.0
    # Its slope at zero is GM.
    step = 1e-6
    slope = (
        wall_sided_righting_arm(step, kb=KB, bm=BM, kg=KG)
        - wall_sided_righting_arm(-step, kb=KB, bm=BM, kg=KG)
    ) / (2.0 * step)
    assert slope == pytest.approx(GM, rel=1e-12)


def test_the_closed_form_stops_holding_once_the_deck_edge_dips(barge):
    """Past deck immersion the barge rights itself less than the wall-sided form claims."""
    assert DECK_EDGE_DEGREES == pytest.approx(34.99, abs=0.02)
    past = DECK_EDGE_DEGREES + 10.0
    curve = gz_curve(barge, (past,), kg=KG, displacement_volume=DISPLACEMENT)
    expected = wall_sided_righting_arm(math.radians(past), kb=KB, bm=BM, kg=KG)
    assert curve.gz[0] < expected - 1e-3
    # Still an increasing righting arm at that angle -- the fall comes later.
    assert curve.gz[0] > 0.0


def test_upright_righting_arm_is_zero_and_the_curve_is_odd(barge):
    curve = gz_curve(barge, (-10.0, 0.0, 10.0), kg=KG, displacement_volume=DISPLACEMENT)
    minus, upright, plus = curve.gz
    assert upright == pytest.approx(0.0, abs=1e-15)
    assert minus == pytest.approx(-plus, rel=1e-12)


# ---------------------------------------------------------------------------
# The headline cross-check: dGZ/dphi at zero is GM
# ---------------------------------------------------------------------------


def test_initial_slope_reproduces_gm_on_the_barge(barge):
    result = gz_initial_slope(barge, kg=KG, waterline_z=DRAFT, step_degrees=0.1)
    assert result["kb"] == pytest.approx(KB, rel=1e-13)
    assert result["bm"] == pytest.approx(BM, rel=1e-13)
    assert result["gm"] == pytest.approx(GM, rel=1e-13)
    assert result["slope"] == pytest.approx(GM, rel=1e-5)


def test_initial_slope_converges_at_second_order(barge):
    coarse = abs(gz_initial_slope(barge, kg=KG, waterline_z=DRAFT, step_degrees=0.4)["difference"])
    fine = abs(gz_initial_slope(barge, kg=KG, waterline_z=DRAFT, step_degrees=0.1)["difference"])
    # A central difference of a smooth odd curve is O(h^2): a factor of four in h, sixteen
    # in the error.
    assert coarse / fine == pytest.approx(16.0, rel=0.05)


def test_initial_slope_reproduces_gm_on_the_wigley(wigley_mesh):
    """Nothing here is wall-sided: the reference GM comes from the waterplane integral."""
    mesh, waterline = wigley_mesh
    result = gz_initial_slope(mesh, kg=0.2, waterline_z=waterline, step_degrees=0.25)
    gm = waterplane_moments(mesh, waterline).transverse_second_moment / submerged_properties(
        mesh, (0.0, 0.0, 1.0), waterline
    )[0] + result["kb"] - 0.2
    assert result["gm"] == pytest.approx(gm, rel=1e-12)
    assert result["slope"] == pytest.approx(result["gm"], rel=1e-4)


# ---------------------------------------------------------------------------
# Constant displacement, and the KN decomposition
# ---------------------------------------------------------------------------


def test_every_heel_angle_displaces_the_same_volume(wigley_mesh):
    mesh, waterline = wigley_mesh
    target = submerged_properties(mesh, (0.0, 0.0, 1.0), waterline)[0]
    curve = gz_curve(mesh, (0.0, 5.0, 10.0, 15.0), kg=0.2, displacement_volume=target)
    assert curve.displacement_volume == pytest.approx(target, rel=1e-15)
    for point in curve.points:
        assert point.displacement_volume == pytest.approx(target, rel=1e-12)


def test_kn_is_gz_plus_kg_sin_phi(barge):
    curve = gz_curve(barge, (0.0, 10.0, 20.0), kg=KG, displacement_volume=DISPLACEMENT)
    for point, arm in zip(curve.points, curve.gz):
        assert arm == pytest.approx(
            point.kn - (curve.keel_z + KG) * math.sin(point.heel_radians), rel=1e-12, abs=1e-15
        )


def test_a_lateral_centre_of_gravity_shifts_the_curve_by_yg_cos_phi(barge):
    plain = gz_curve(barge, (0.0, 10.0, 20.0), kg=KG, displacement_volume=DISPLACEMENT)
    shifted = gz_curve(barge, (0.0, 10.0, 20.0), kg=KG, yg=0.1, displacement_volume=DISPLACEMENT)
    for before, after, point in zip(plain.gz, shifted.gz, plain.points):
        assert after == pytest.approx(before - 0.1 * math.cos(point.heel_radians), rel=1e-12, abs=1e-15)


def _barge_kn(heel_degrees: float, draught: float) -> float:
    """Closed-form ``KN`` of a wall-sided barge heeled with the waterline on both sides.

    The section ``|y| <= b`` between ``z = 0`` and ``z = T + y tan(phi)`` has area
    ``2 b T`` independent of ``phi`` (the linear term integrates out), so the floating
    draught does not change with heel and

        y_B = tan(phi) b^2 / (3 T),   z_B = T/2 + b^2 tan^2(phi) / (6 T),

    from which ``KN = cos(phi) y_B + sin(phi) z_B``. Every line of that is independent of
    the clipping kernel under test.
    """
    phi = math.radians(heel_degrees)
    half_beam, tan = BEAM / 2.0, math.tan(math.radians(heel_degrees))
    y_b = tan * half_beam**2 / (3.0 * draught)
    z_b = draught / 2.0 + half_beam**2 * tan**2 / (6.0 * draught)
    return math.cos(phi) * y_b + math.sin(phi) * z_b


def test_cross_curves_match_the_closed_form_barge(barge):
    heels = (0.0, 10.0, 20.0)
    draughts = (0.5 * DRAFT, DRAFT)
    table = cross_curves(barge, heels, tuple(LENGTH * BEAM * t for t in draughts))
    assert len(table) == 2 and all(len(row) == 3 for row in table)
    for row, draught in zip(table, draughts):
        for value, heel in zip(row, heels):
            assert value == pytest.approx(_barge_kn(heel, draught), rel=1e-12, abs=1e-15)
    # A symmetric barge has no upright KN at all.
    assert table[1][0] == pytest.approx(0.0, abs=1e-15)
    single = cross_curve(barge, 20.0, displacement_volume=DISPLACEMENT)
    assert single.kn == pytest.approx(table[1][2], rel=1e-15)


def test_at_ninety_degrees_the_box_floats_on_its_deck_edge(barge):
    """Exact at the end of the sweep: the waterplane is the ``L x D`` face and ``KN = D/2``.

    Heeled 90 degrees the body's ``z`` axis is horizontal, so the waterplane spans the
    length and the depth and its centroid sits at mid-depth whatever the displacement --
    no matter how much of the side is under water. Hence ``KN(90) = DEPTH / 2`` and
    ``GZ(90) = DEPTH / 2 - (keel + kg)``, exactly.
    """
    curve = gz_curve(barge, (90.0,), kg=0.2, displacement_volume=DISPLACEMENT)
    assert curve.kn[0] == pytest.approx(0.5 * DEPTH, rel=1e-13)
    assert curve.gz[0] == pytest.approx(0.5 * DEPTH - 0.2, rel=1e-12)


def test_curve_extremes_of_a_flat_topped_barge(barge):
    """A rectangular barge keeps righting to 90 degrees: it has no superstructure to dip.

    The measured curve peaks near 70 degrees and is still positive at 90, so
    ``first_vanishing_angle`` is ``None`` -- reported, not invented.
    """
    curve = gz_curve(barge, tuple(range(0, 95, 5)), kg=0.2, displacement_volume=DISPLACEMENT)
    angle, arm = curve.maximum()
    assert 60.0 <= angle <= 80.0
    assert 0.55 < arm < 0.65
    assert curve.first_vanishing_angle() is None
    assert min(curve.gz) >= -1e-15  # round-off at upright; the arm is zero, not negative


def test_gz_at_only_answers_for_sampled_angles(barge):
    curve = gz_curve(barge, (0.0, 10.0), kg=KG, displacement_volume=DISPLACEMENT)
    assert curve.gz_at(10.0) == pytest.approx(curve.gz[1], rel=1e-15)
    with pytest.raises(HeelError, match="not one of the sampled"):
        curve.gz_at(11.0)


# ---------------------------------------------------------------------------
# Discretisation: the curved hull's own error
# ---------------------------------------------------------------------------


def test_wigley_gz_converges_as_the_mesh_refines():
    """The only error left is the chord error of the hull's sides; measure it."""
    parameters = reference_parameters()
    waterline = -0.3 * parameters.draft
    heels = (5.0, 10.0)
    values = {}
    for nx, nz in ((21, 9), (41, 13), (81, 25)):
        mesh = generate_wigley(
            reference_parameters(nx=nx, nz=nz)
        ).mesh
        target = submerged_properties(mesh, (0.0, 0.0, 1.0), waterline)[0]
        values[(nx, nz)] = gz_curve(mesh, heels, kg=0.2, displacement_volume=target).gz
    fine = values[(81, 25)]
    for index in range(len(heels)):
        coarse_error = abs(values[(21, 9)][index] - fine[index])
        medium_error = abs(values[(41, 13)][index] - fine[index])
        assert medium_error < coarse_error
        # The coarse grid is already within 1 % of the fine one: the hull is curved but the
        # heeled waterline is not a quantity that amplifies the chord error.
        assert coarse_error <= 1e-2 * abs(fine[index])
    assert parameters.nx == 41 and parameters.nz == 13


def test_barge_gz_is_exact_so_the_kernel_carries_no_discretisation_error(barge):
    """A plane-sided body has no chord error: agreement with the closed form is round-off."""
    curve = gz_curve(barge, (7.0,), kg=KG, displacement_volume=DISPLACEMENT)
    expected = wall_sided_righting_arm(math.radians(7.0), kb=KB, bm=BM, kg=KG)
    assert abs(curve.gz[0] - expected) <= 1e-15 * expected


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


def test_kg_must_be_finite(barge):
    with pytest.raises(HeelError, match="kg must be finite"):
        gz_curve(barge, (0.0,), kg=float("inf"), displacement_volume=DISPLACEMENT)


def test_displacement_must_be_positive(barge):
    with pytest.raises(HeelError, match="positive"):
        gz_curve(barge, (0.0,), kg=KG, displacement_volume=0.0)


def test_barge_dimensions_are_validated():
    for kwargs in ({"length": 0.0, "beam": 1.0, "depth": 1.0}, {"length": 1.0, "beam": -1.0, "depth": 1.0}):
        with pytest.raises(HeelError):
            barge_mesh(**kwargs)


def test_step_degrees_is_bounded(barge):
    with pytest.raises(HeelError, match="step_degrees"):
        gz_initial_slope(barge, kg=KG, waterline_z=DRAFT, step_degrees=0.0)
    with pytest.raises(HeelError, match="step_degrees"):
        gz_initial_slope(barge, kg=KG, waterline_z=DRAFT, step_degrees=10.0)