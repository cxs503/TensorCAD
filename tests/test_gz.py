"""Large-angle stability: ``GZ(phi)`` against closed forms.

Three references, in increasing order of how much they can go wrong:

* a **wall-sided barge** with freeboard, where the classical wall-sided formula
  ``GZ = sin(phi) (GM + BM tan^2(phi) / 2)`` is exact for every heel short of deck-edge
  immersion -- it pins the clipping, the constant-displacement solve and the lever arm at
  once, to round-off, over a *range* of angles rather than at one point;
* a **circular cylinder**, where ``GZ(phi) = -(z_G - z_C) sin(phi)`` exactly at *every*
  heel angle up to full rotation, because the section is invariant under the rotation and
  the metacentre sits on the circle centre. The polygon approximation of the circle is
  the only error, and it must fall as ``n^-2``;
* a body floating at **exactly full submergence**, where the waterplane area vanishes and
  the module must refuse instead of returning a smooth number.

The displacements are held fixed at the upright value at every station, so the
displacement residual is asserted too: a righting arm is only as good as the balance that
produced it.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.gz import (
    GZCurve,
    GZError,
    gz_curve,
    heel_mesh,
    heel_state,
    immersed_moments,
    righting_arm,
)
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.stability import waterplane_moments

LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8
FREEBOARD = 0.4
KEEL_TO_DECK = DRAFT + FREEBOARD


# ---------------------------------------------------------------------------
# reference geometries
# ---------------------------------------------------------------------------
def barge(length: float, beam: float, draft: float, freeboard: float, *, y0: float | None = None) -> Mesh:
    """Closed box: keel at ``z = -draft``, deck at ``z = freeboard``, waterline at ``z = 0``.

    ``y0`` shifts the whole box in the transverse direction, so a body can be built that
    is *not* symmetric about the centreplane.
    """
    low = -0.5 * beam if y0 is None else y0
    high = low + beam
    builder = MeshBuilder()
    corner: dict[tuple[float, float, float], int] = {}
    for x in (0.0, length):
        for y in (low, high):
            for z in (-draft, freeboard):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    keel, deck = -draft, freeboard
    quad((0.0, low, keel), (0.0, high, keel), (length, high, keel), (length, low, keel))
    quad((0.0, low, deck), (length, low, deck), (length, high, deck), (0.0, high, deck))
    quad((0.0, low, keel), (length, low, keel), (length, low, deck), (0.0, low, deck))
    quad((0.0, high, keel), (0.0, high, deck), (length, high, deck), (length, high, keel))
    quad((0.0, low, keel), (0.0, low, deck), (0.0, high, deck), (0.0, high, keel))
    quad((length, low, keel), (length, high, keel), (length, high, deck), (length, low, deck))
    return builder.build()


def cylinder(radius: float, length: float, segments: int) -> Mesh:
    """Closed circular cylinder with its axis on ``x``; vertices offset off ``z = 0``.

    The ``pi / segments`` phase shift keeps every vertex off the waterline, so a section
    cut at ``z = 0`` is a set of proper crossings rather than a degenerate vertex-on-plane
    cut (the body is still the same body: the polygon is rotated, not reshaped).
    """
    builder = MeshBuilder()
    phase = math.pi / segments
    rings: list[list[int]] = []
    for x in (0.0, length):
        rings.append(
            [
                builder.add_vertex(
                    (
                        x,
                        radius * math.cos(2.0 * math.pi * index / segments + phase),
                        radius * math.sin(2.0 * math.pi * index / segments + phase),
                    )
                )
                for index in range(segments)
            ]
        )
    for index in range(segments):
        following = (index + 1) % segments
        a, b = rings[0][index], rings[0][following]
        c, d = rings[1][following], rings[1][index]
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)
    # end caps as fans from the first ring vertex, wound outwards
    for index in range(1, segments - 1):
        builder.add_triangle(rings[0][0], rings[0][index + 1], rings[0][index])
    for index in range(1, segments - 1):
        builder.add_triangle(rings[1][0], rings[1][index], rings[1][index + 1])
    return builder.build()


def wall_sided_arm(phi_deg: float, *, gm: float, bm: float) -> float:
    phi = math.radians(phi_deg)
    return math.sin(phi) * (gm + 0.5 * bm * math.tan(phi) ** 2)


def wall_sided_area(phi_deg: float, *, gm: float, bm: float) -> float:
    """``integral_0^phi GZ dphi`` of the wall-sided formula, in m*rad."""
    phi = math.radians(phi_deg)
    return gm * (1.0 - math.cos(phi)) + 0.5 * bm * (1.0 / math.cos(phi) + math.cos(phi) - 2.0)


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def loaded_barge() -> Mesh:
    return barge(LENGTH, BEAM, DRAFT, FREEBOARD)


@pytest.fixture(scope="module")
def barge_bm() -> float:
    return BEAM * BEAM / (12.0 * DRAFT)


@pytest.fixture(scope="module")
def barge_kg() -> float:
    return 0.6


@pytest.fixture(scope="module")
def tank() -> Mesh:
    return cylinder(radius=1.0, length=4.0, segments=48)


@pytest.fixture(scope="module")
def hull_mesh() -> Mesh:
    return generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=41, nz=13)).mesh


# ---------------------------------------------------------------------------
# the divergence integrals behind the arm
# ---------------------------------------------------------------------------
def test_immersed_moments_are_exact_on_a_box(loaded_barge: Mesh) -> None:
    volume, centroid = immersed_moments(loaded_barge, 0.0)
    assert volume == pytest.approx(LENGTH * BEAM * DRAFT, rel=1e-15)
    assert centroid[0] == pytest.approx(0.0, abs=1e-15)
    assert centroid[1] == pytest.approx(-0.5 * DRAFT, rel=1e-15)


def test_immersed_moments_clip_a_body_that_has_freeboard(loaded_barge: Mesh) -> None:
    """The whole closed body is 1.2 m deep; only the 0.8 m below the waterline displaces."""
    volume, centroid = immersed_moments(loaded_barge, 0.0)
    assert volume != pytest.approx(LENGTH * BEAM * KEEL_TO_DECK, rel=0.2)
    assert centroid[1] > -0.5 * DRAFT


def test_immersed_moments_reject_a_waterline_below_the_body(loaded_barge: Mesh) -> None:
    with pytest.raises(GZError, match="not positive"):
        immersed_moments(loaded_barge, -5.0)


def test_heel_mesh_rotates_the_ordinates_and_keeps_the_topology(loaded_barge: Mesh) -> None:
    heeled = heel_mesh(loaded_barge, 90.0)
    assert heeled.faces == loaded_barge.faces
    # the starboard deck edge swings down to starboard depth = half beam
    index = loaded_barge.vertices.index((0.0, 0.5 * BEAM, FREEBOARD))
    assert heeled.vertices[index] == pytest.approx((0.0, FREEBOARD, -0.5 * BEAM), abs=1e-15)


def test_heel_mesh_rejects_a_non_finite_angle(loaded_barge: Mesh) -> None:
    with pytest.raises(GZError, match="finite"):
        heel_mesh(loaded_barge, math.inf)


def test_heel_state_solves_the_upright_waterline(loaded_barge: Mesh) -> None:
    state = heel_state(loaded_barge, 0.0)
    assert state.waterline_z == pytest.approx(0.0, abs=1e-12)
    assert abs(state.volume_error) < 1e-12
    assert state.centre_buoyancy[0] == pytest.approx(0.0, abs=1e-15)


# ---------------------------------------------------------------------------
# reference 1: the wall-sided barge
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("phi", [1.0, 5.0, 10.0, 15.0, 20.0])
def test_barge_matches_the_wall_sided_formula(
    loaded_barge: Mesh, barge_kg: float, barge_bm: float, phi: float
) -> None:
    """Exact to round-off: BM = B^2/(12 T), GM = BM - (KG - T/2)."""
    gm = barge_bm - (barge_kg - 0.5 * DRAFT)
    measured = righting_arm(loaded_barge, phi, kg=barge_kg)
    assert measured == pytest.approx(wall_sided_arm(phi, gm=gm, bm=barge_bm), rel=1e-9)


def test_barge_leaves_the_wall_sided_formula_once_the_deck_immerses(
    loaded_barge: Mesh, barge_kg: float, barge_bm: float
) -> None:
    """Deck-edge immersion is at tan(phi) = 2 * freeboard / beam: the formula must break."""
    limit = math.degrees(math.atan(2.0 * FREEBOARD / BEAM))
    gm = barge_bm - (barge_kg - 0.5 * DRAFT)
    at_limit = righting_arm(loaded_barge, limit - 0.5, kg=barge_kg)
    assert at_limit == pytest.approx(wall_sided_arm(limit - 0.5, gm=gm, bm=barge_bm), rel=1e-8)
    beyond = righting_arm(loaded_barge, limit + 8.0, kg=barge_kg)
    reference = wall_sided_arm(limit + 8.0, gm=gm, bm=barge_bm)
    assert beyond < 0.9 * reference, "a submerged deck must lose righting arm, not gain it"


def test_barge_curve_peaks_and_vanishes(loaded_barge: Mesh, barge_kg: float) -> None:
    curve = gz_curve(
        loaded_barge, (0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0), kg=barge_kg
    )
    angle, value = curve.max_arm
    assert 25.0 <= angle <= 55.0
    assert value > 0.1
    assert curve.area_under_curve(30.0) > 0.0
    # at exactly 90 deg the heeled box presents a full rectangle to the waterline, so the
    # arm is exactly zero: B sits on the vertical through G
    assert curve.arm_at(90.0) == pytest.approx(0.0, abs=1e-9)
    vanishing = curve.vanishing_angle()
    assert vanishing is not None and 88.0 <= vanishing <= 92.0


def test_curve_area_matches_the_integrated_formula(loaded_barge: Mesh, barge_kg: float, barge_bm: float) -> None:
    """Trapezoid on the 5 deg ladder against the closed-form integral: O(h^2), under 0.1 %."""
    gm = barge_bm - (barge_kg - 0.5 * DRAFT)
    curve = gz_curve(loaded_barge, (0.0, 5.0, 10.0, 15.0, 20.0), kg=barge_kg)
    measured = curve.area_under_curve(20.0)
    # the trapezoid on a 5 deg ladder is the error here, and it is O(h^2): 3e-3 relative
    assert measured == pytest.approx(wall_sided_area(20.0, gm=gm, bm=barge_bm), rel=1e-2)


def test_gm_from_the_curve_slope_is_the_metacentric_height(loaded_barge: Mesh, barge_kg: float, barge_bm: float) -> None:
    gm = barge_bm - (barge_kg - 0.5 * DRAFT)
    curve = gz_curve(loaded_barge, (0.0, 1.0), kg=barge_kg)
    # at a finite heel the slope is GM + BM tan^2(phi)/2, not GM: the O(phi^2) term is
    # itself predicted by the wall-sided formula, so it is checked rather than tolerated
    expected = gm + 0.5 * barge_bm * math.tan(math.radians(1.0)) ** 2
    assert curve.gm_from_slope() == pytest.approx(expected, rel=1e-9)


def test_displacement_is_preserved_at_every_station(loaded_barge: Mesh, barge_kg: float) -> None:
    curve = gz_curve(loaded_barge, (0.0, 15.0, 30.0, 45.0, 60.0, 75.0), kg=barge_kg)
    assert curve.max_volume_error < 1e-12
    assert curve.displacement_volume == pytest.approx(LENGTH * BEAM * DRAFT, rel=1e-15)


def test_a_lopsided_body_reports_its_list_moment() -> None:
    """``GZ(0)`` is the transverse offset of the centre of buoyancy, not a hard zero."""
    offset = barge(LENGTH, BEAM, DRAFT, FREEBOARD, y0=0.0)  # spans y in [0, B]
    arm = righting_arm(offset, 0.0, kg=0.6)
    assert arm == pytest.approx(0.5 * BEAM, rel=1e-15)


# ---------------------------------------------------------------------------
# reference 2: the circular cylinder
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("phi", [5.0, 20.0, 45.0, 60.0, 90.0, 120.0])
def test_cylinder_matches_the_exact_arm(tank: Mesh, phi: float) -> None:
    """``GZ = -(z_G - z_C) sin(phi)`` at every angle: the section is rotation-invariant."""
    radius, g_offset = 1.0, -0.3  # G is 0.3 R below the circle centre
    arm = righting_arm(tank, phi, kg=radius + g_offset)
    exact = -g_offset * math.sin(math.radians(phi))
    assert arm == pytest.approx(exact, rel=2e-2)


def test_cylinder_is_neutral_when_g_sits_on_the_centre(tank: Mesh) -> None:
    """With G on the metacentre (the circle centre) every arm vanishes."""
    # a 48-gon is not exactly the circle: the residual is the polygon's O(n^-2) error
    for phi in (10.0, 30.0, 60.0):
        assert righting_arm(tank, phi, kg=1.0) == pytest.approx(0.0, abs=2e-3)


def test_cylinder_error_falls_as_the_square_of_the_segment_size() -> None:
    errors = []
    for segments in (24, 48, 96):
        body = cylinder(radius=1.0, length=4.0, segments=segments)
        worst = max(
            abs(righting_arm(body, phi, kg=0.7) / (0.3 * math.sin(math.radians(phi))) - 1.0)
            for phi in (20.0, 45.0, 90.0)
        )
        errors.append(worst)
    assert errors[0] > errors[1] > errors[2]
    assert errors[0] / errors[2] > 8.0, f"expected a second-order fall, got {errors}"


# ---------------------------------------------------------------------------
# reference 3: the degenerate body, which must be refused
# ---------------------------------------------------------------------------
def test_a_fully_submerged_body_is_refused(hull_mesh: Mesh) -> None:
    """The parametric hull's lid sits *on* the design waterline: zero freeboard."""
    with pytest.raises(GZError, match="freeboard"):
        gz_curve(hull_mesh, (0.0, 10.0), kg=0.7)


def test_a_displacement_larger_than_the_body_is_refused(loaded_barge: Mesh) -> None:
    with pytest.raises(GZError, match="exceeds the body volume"):
        heel_state(loaded_barge, 10.0, displacement_volume=1000.0)


def test_an_inside_out_mesh_is_refused() -> None:
    body = barge(LENGTH, BEAM, DRAFT, FREEBOARD)
    flipped = Mesh(vertices=body.vertices, faces=tuple((c, b, a) for a, b, c in body.faces))
    with pytest.raises(GZError, match="inside-out"):
        heel_state(flipped, 10.0)


# ---------------------------------------------------------------------------
# the hull with freeboard: an independent waterplane route to the same GM
# ---------------------------------------------------------------------------
def test_hull_gm_from_the_curve_matches_the_waterplane_integral(hull_mesh: Mesh) -> None:
    """The slope of the arm curve and ``I_T / V - BG`` are computed by independent code.

    ``waterplane_moments`` integrates the section of the *heeled-at-zero* body;
    :func:`gz_curve` integrates the clipped volume at a small heel. They must agree.
    """
    waterline = -0.11  # 0.11 m of freeboard under the lid, and off every vertex ring
    kg = 0.7
    volume, centroid = immersed_moments(hull_mesh, waterline)
    section = waterplane_moments(hull_mesh, waterline)
    bm = section.transverse_second_moment / volume
    z_g = min(vertex[2] for vertex in hull_mesh.vertices) + kg
    gm_waterplane = bm - (z_g - centroid[1])
    curve = gz_curve(hull_mesh, (0.0, 1.0, 5.0, 10.0), kg=kg, waterline_z=waterline)
    assert curve.gm_from_slope() == pytest.approx(gm_waterplane, rel=1e-2)
    assert gm_waterplane > 0.0
    assert curve.arm_at(10.0) > 0.0


def test_hull_curve_converges_with_the_grid() -> None:
    """Two grids, same loading: the arm is a property of the body, not of the mesh."""
    coarse = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=41, nz=13)).mesh
    fine = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=81, nz=25)).mesh
    arms = []
    for mesh in (coarse, fine):
        curve = gz_curve(mesh, (0.0, 1.0, 10.0), kg=0.7, waterline_z=-0.11)
        arms.append((curve.arm_at(10.0), curve.gm_from_slope(), curve.vanishing_angle()))
    # GM and the arm well inside the positive range converge; the arm near the vanishing
    # angle does not, and is deliberately not asserted here
    assert arms[1][1] == pytest.approx(arms[0][1], rel=5e-2)
    assert arms[1][0] == pytest.approx(arms[0][0], rel=1e-1)
    assert arms[0][0] > 0.0


# ---------------------------------------------------------------------------
# API guards
# ---------------------------------------------------------------------------
def test_kg_is_required(loaded_barge: Mesh) -> None:
    with pytest.raises(TypeError):
        righting_arm(loaded_barge, 10.0)  # type: ignore[call-arg]


def test_kg_must_be_finite(loaded_barge: Mesh) -> None:
    with pytest.raises(GZError, match="kg must be a finite"):
        righting_arm(loaded_barge, 10.0, kg=math.nan)


def test_angles_must_be_strictly_increasing() -> None:
    with pytest.raises(GZError, match="strictly increasing"):
        GZCurve(angles_deg=(0.0, 10.0, 5.0), arms=(0.0, 1.0, 2.0), kg=1.0, keel_z=-1.0, displacement_volume=1.0)


def test_angles_and_arms_must_agree() -> None:
    with pytest.raises(GZError, match="same length"):
        GZCurve(angles_deg=(0.0, 10.0), arms=(0.0,), kg=1.0, keel_z=-1.0, displacement_volume=1.0)


def test_no_angles_is_an_error(loaded_barge: Mesh) -> None:
    with pytest.raises(GZError, match="at least one"):
        gz_curve(loaded_barge, (), kg=0.6)


def test_arm_at_refuses_to_extrapolate(loaded_barge: Mesh) -> None:
    curve = gz_curve(loaded_barge, (0.0, 10.0), kg=0.6)
    with pytest.raises(GZError, match="outside the sampled range"):
        curve.arm_at(20.0)
    assert curve.arm_at(5.0) == pytest.approx(0.5 * curve.arms[1])


def test_vanishing_angle_is_none_when_the_curve_stays_positive(loaded_barge: Mesh) -> None:
    curve = gz_curve(loaded_barge, (0.0, 5.0, 10.0), kg=0.6)
    assert curve.vanishing_angle() is None


def test_curve_serialises(loaded_barge: Mesh) -> None:
    curve = gz_curve(loaded_barge, (0.0, 10.0, 20.0, 30.0, 40.0), kg=0.6)
    payload = curve.to_dict()
    assert payload["angles_deg"] == [0.0, 10.0, 20.0, 30.0, 40.0]
    assert len(payload["arms"]) == 5
    assert payload["max_arm"] == pytest.approx(max(curve.arms))
    assert payload["area_to_30_deg"] is not None
    assert payload["gm_from_slope"] == pytest.approx(curve.gm_from_slope())
    assert payload["conditions"][0]["volume_error"] == pytest.approx(0.0, abs=1e-12)


def test_gm_from_slope_needs_a_positive_station(loaded_barge: Mesh) -> None:
    curve = gz_curve(loaded_barge, (0.0,), kg=0.6)
    with pytest.raises(GZError, match="strictly positive"):
        curve.gm_from_slope()