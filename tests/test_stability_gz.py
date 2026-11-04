"""Large-angle stability: righting arms against two exact references.

``GM`` is a waterplane moment and a slope; ``GZ(phi)`` is a heeled displacement
integration. Two bodies have a righting arm in closed form, and neither of them shares
code with the integration that is being tested:

* a **wall-sided box** (vertical sides, deck edge not immersed), where the textbook
  result is

      GZ(phi) = sin(phi) * ( GM + BM * tan^2(phi) / 2 ),   BM = B^2 / (12 T),

  which is *exact* up to the angle at which a deck corner enters the water,
  ``tan(phi) = 2 f / B`` for freeboard ``f``. Past that angle the reference is wrong, not
  the solver, and one test below pins exactly that boundary.
* a **circular section**, whose transverse metacentre is the centre of the circle
  (``KM = R(1 - 4/3pi) + 4R/3pi = R``), so ``GZ(phi) = (R - KG) sin(phi)`` at *every*
  heel angle while the section stays a half-disc. A regular polygon inscribed in that
  circle must approach it from below at second order in the vertex spacing.

Both references are computed here from their formulas, not read back from the module, so
an error in the heeled integration cannot cancel against an error in the reference.
"""

from __future__ import annotations

import itertools
import math

import pytest

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.stability import (
    GZCurve,
    RightingArm,
    StabilityError,
    SubmergedMoments,
    gz_curve,
    initial_stability,
    righting_arm,
    submerged_moments,
)


# ---------------------------------------------------------------------------
# reference bodies
# ---------------------------------------------------------------------------
def box(length: float, beam: float, draft: float, freeboard: float) -> Mesh:
    """Closed box: ``x`` in ``[0, L]``, ``y`` in ``[-B/2, B/2]``, ``z`` in ``[-T, f]``.

    The design waterline is ``z = 0``, so ``draft`` is below it and ``freeboard`` above
    it. Unlike the meshes :func:`tensorcad.hydrostatics.compute_hydrostatics` accepts,
    the side faces straddle the waterline on purpose: a righting arm is computed by
    clipping, and the body a righting arm needs *has* freeboard.
    """
    builder = MeshBuilder()
    corner: dict[tuple[float, float, float], int] = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, freeboard):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    left, right = 0.0, length
    lo, hi = -0.5 * beam, 0.5 * beam
    quad((left, lo, -draft), (left, hi, -draft), (right, hi, -draft), (right, lo, -draft))
    quad((left, lo, freeboard), (right, lo, freeboard), (right, hi, freeboard), (left, hi, freeboard))
    quad((left, lo, -draft), (right, lo, -draft), (right, lo, freeboard), (left, lo, freeboard))
    quad((right, hi, -draft), (left, hi, -draft), (left, hi, freeboard), (right, hi, freeboard))
    quad((right, lo, -draft), (right, hi, -draft), (right, hi, freeboard), (right, lo, freeboard))
    quad((left, hi, -draft), (left, lo, -draft), (left, lo, freeboard), (left, hi, freeboard))
    return builder.build()


def prism(sides: int, radius: float, length: float) -> Mesh:
    """Regular ``sides``-gon prism about the ``x`` axis, waterline through the axis.

    The polygon is offset by half a vertex spacing so that no vertex sits exactly on the
    waterline; a vertex *on* the cut plane leaves the section with zero-length edges,
    which is a degenerate loop no reasonable waterline integrator should be asked for.
    """
    builder = MeshBuilder()
    ring = [
        (radius * math.cos(2.0 * math.pi * (index + 0.5) / sides),
         radius * math.sin(2.0 * math.pi * (index + 0.5) / sides))
        for index in range(sides)
    ]
    near = [builder.add_vertex((0.0, y, z)) for y, z in ring]
    far = [builder.add_vertex((length, y, z)) for y, z in ring]
    for index in range(sides):
        following = (index + 1) % sides
        builder.add_triangle(near[index], near[following], far[following])
        builder.add_triangle(near[index], far[following], far[index])
    for index in range(1, sides - 1):
        builder.add_triangle(near[0], near[index + 1], near[index])
        builder.add_triangle(far[0], far[index], far[index + 1])
    return builder.build()


def subdivide(mesh: Mesh) -> Mesh:
    """Split every triangle into four: the same body, a different tessellation."""
    vertices = list(mesh.vertices)
    faces: list[tuple[int, int, int]] = []
    midpoints: dict[tuple[int, int], int] = {}

    def midpoint(a: int, b: int) -> int:
        key = (min(a, b), max(a, b))
        if key not in midpoints:
            pa, pb = vertices[a], vertices[b]
            vertices.append(tuple(0.5 * (pa[axis] + pb[axis]) for axis in range(3)))  # type: ignore[arg-type]
            midpoints[key] = len(vertices) - 1
        return midpoints[key]

    for a, b, c in mesh.faces:
        ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
        faces.extend(((a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)))
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


def wall_sided_arm(heel_degrees: float, *, gm: float, bm: float) -> float:
    """Textbook righting arm of a wall-sided vessel, valid until a deck corner immerses."""
    phi = math.radians(heel_degrees)
    return math.sin(phi) * (gm + 0.5 * bm * math.tan(phi) ** 2)


# the tall box: sides deep enough that the closed form survives to 60 deg
TALL = {"length": 6.0, "beam": 4.0, "draft": 4.0, "freeboard": 4.0}
TALL_KG = 1.0
TALL_BM = TALL["beam"] ** 2 / (12.0 * TALL["draft"])
TALL_GM = 0.5 * TALL["draft"] + TALL_BM - TALL_KG


# ---------------------------------------------------------------------------
# reference 1: the wall-sided box
# ---------------------------------------------------------------------------
def test_wall_sided_box_matches_the_closed_form_arm() -> None:
    """Everything up to 60 deg, where ``tan(phi) = 1.732 < 2 f / B = 2``."""
    mesh = box(**TALL)
    curve = gz_curve(mesh, kg=TALL_KG, heel_angles_deg=tuple(range(0, 61, 5)))
    worst = max(
        abs(arm.righting_arm - wall_sided_arm(arm.heel_angle_deg, gm=TALL_GM, bm=TALL_BM))
        for arm in curve.arms
    )
    # The clipped integral is exact for a polyhedron, so this is round-off and not a
    # discretisation error: a tolerance with any slack in it would hide a real mistake.
    assert worst < 1e-12
    assert curve.arms[0].righting_arm == pytest.approx(0.0, abs=1e-12)


def test_the_closed_form_is_valid_exactly_until_a_deck_corner_immerses() -> None:
    """The reference must break at ``atan(2 f / B) = 63.43 deg``, and the solver must not.

    This is the check that the two agree for the *right* reason: a solver with a sign
    error would disagree everywhere, and one that silently ignored the waterline would
    keep agreeing past the corner too.
    """
    mesh = box(**TALL)
    limit = math.degrees(math.atan(2.0 * TALL["freeboard"] / TALL["beam"]))
    assert limit == pytest.approx(63.43494882292201, abs=1e-9)

    inside = gz_curve(mesh, kg=TALL_KG, heel_angles_deg=(60.0,)).arms[0]
    outside = gz_curve(mesh, kg=TALL_KG, heel_angles_deg=(70.0,)).arms[0]
    assert abs(inside.righting_arm - wall_sided_arm(60.0, gm=TALL_GM, bm=TALL_BM)) < 1e-12
    assert outside.heel_angle_deg > limit
    assert abs(outside.righting_arm - wall_sided_arm(70.0, gm=TALL_GM, bm=TALL_BM)) > 0.1


def test_small_angle_slope_is_the_analytic_metacentric_height() -> None:
    """``dGZ/dphi`` at zero heel is ``GM``, from an independent formula."""
    mesh = box(**TALL)
    analytic_gm = 0.5 * TALL["draft"] + TALL_BM - TALL_KG
    for heel in (0.25, 0.5):
        arm = gz_curve(mesh, kg=TALL_KG, heel_angles_deg=(heel,)).arms[0]
        slope = arm.righting_arm / math.radians(heel)
        assert slope == pytest.approx(analytic_gm, rel=2e-4)


# ---------------------------------------------------------------------------
# reference 2: the circular section
# ---------------------------------------------------------------------------
def test_circular_section_righting_arm_is_km_minus_kg_times_sine() -> None:
    """``KM = R`` for a circular section, so ``GZ = (R - KG) sin(phi)`` at every angle."""
    radius, kg = 1.0, 0.6
    expected = radius - kg
    errors = []
    for sides in (24, 48, 96, 192):
        curve = gz_curve(prism(sides, radius, 6.0), kg=kg, heel_angles_deg=(10, 20, 30, 45, 60, 75))
        error = max(
            abs(arm.righting_arm / math.sin(math.radians(arm.heel_angle_deg)) - expected) / expected
            for arm in curve.arms
        )
        errors.append(error)
    # An inscribed polygon is smaller than the circle, so the error must fall at second
    # order in the vertex spacing: each doubling of the sides divides it by about four.
    assert errors[-1] < 5e-4
    for coarse, fine in itertools.pairwise(errors):
        assert fine < coarse / 3.0


# ---------------------------------------------------------------------------
# invariants that hold for any body
# ---------------------------------------------------------------------------
def test_displacement_is_held_at_every_heel_angle() -> None:
    """The sinkage is solved for, so the displaced volume must not drift with heel."""
    for mesh, kg in ((box(**TALL), TALL_KG), (prism(96, 1.0, 6.0), 0.6)):
        curve = gz_curve(mesh, kg=kg, heel_angles_deg=(0, 15, 30, 45, 60))
        assert curve.max_volume_residual < 1e-9
        for arm in curve.arms:
            assert arm.displaced_volume == pytest.approx(curve.displacement_volume, rel=1e-9)


def test_upright_body_has_no_righting_arm_and_no_sinkage() -> None:
    mesh = box(**TALL)
    arm = gz_curve(mesh, kg=TALL_KG, heel_angles_deg=(0.0,)).arms[0]
    assert arm.righting_arm == pytest.approx(0.0, abs=1e-12)
    assert arm.sinkage == pytest.approx(0.0, abs=1e-12)


def test_righting_arm_is_antisymmetric_in_the_heel_angle() -> None:
    """Heeling to port is the mirror image, so the arm changes sign, not magnitude."""
    mesh = box(**TALL)
    starboard = gz_curve(mesh, kg=TALL_KG, heel_angles_deg=(20.0, 40.0))
    port = gz_curve(mesh, kg=TALL_KG, heel_angles_deg=(-20.0, -40.0))
    for right, left in zip(starboard.righting_arms, port.righting_arms, strict=False):
        assert left == pytest.approx(-right, abs=1e-12)
        assert right > 0.0


def test_righting_arm_is_independent_of_the_tessellation() -> None:
    """Splitting every triangle in four describes the same body, so the arm cannot move."""
    mesh = box(**TALL)
    coarse = gz_curve(mesh, kg=TALL_KG, heel_angles_deg=(30.0,)).arms[0].righting_arm
    fine = gz_curve(subdivide(mesh), kg=TALL_KG, heel_angles_deg=(30.0,)).arms[0].righting_arm
    assert fine == pytest.approx(coarse, abs=1e-12)


# ---------------------------------------------------------------------------
# the shape of a real righting-arm curve
# ---------------------------------------------------------------------------
# A barge with little freeboard: the deck edge goes under early, the waterplane shrinks,
# and the arm peaks and then disappears -- the qualitative behaviour a GM alone misses.
BARGE = {"length": 4.0, "beam": 2.0, "draft": 1.0, "freeboard": 0.3}
BARGE_ANGLES = tuple(range(0, 121, 10))


def test_curve_has_a_crest_and_a_vanishing_angle() -> None:
    curve = gz_curve(box(**BARGE), kg=0.6, heel_angles_deg=BARGE_ANGLES)
    assert curve.angle_of_max_arm_deg == pytest.approx(40.0)
    assert curve.max_arm.righting_arm > 0.09
    vanishing = curve.vanishing_angle_deg
    assert vanishing is not None
    assert 105.0 < vanishing < 115.0
    # the arm is positive at the crest and negative past the vanishing angle
    assert curve.righting_arms[curve.heel_angles_deg.index(curve.angle_of_max_arm_deg)] > 0.0
    assert curve.righting_arms[-1] < 0.0


def test_a_heavier_load_shortens_the_range_of_stability() -> None:
    """Raising ``KG`` can only lower the curve, so the vanishing angle cannot grow."""
    angles = tuple(range(0, 141, 5))
    vanishings = []
    for kg in (0.6, 0.65, 0.7):
        vanishing = gz_curve(box(**BARGE), kg=kg, heel_angles_deg=angles).vanishing_angle_deg
        assert vanishing is not None
        vanishings.append(vanishing)
    assert vanishings[0] > vanishings[1] > vanishings[2]


def test_a_curve_without_a_positive_arm_has_no_range_of_stability() -> None:
    """``KG`` above ``KM`` leaves no righting arm at all, so the question has no answer."""
    curve = gz_curve(box(length=4.0, beam=1.5, draft=1.0, freeboard=0.4), kg=0.8,
                     heel_angles_deg=tuple(range(0, 91, 10)))
    assert curve.righting_arms[0] == pytest.approx(0.0, abs=1e-15)  # upright, round-off only
    assert all(arm < 0.0 for arm in curve.righting_arms[1:])
    assert curve.vanishing_angle_deg is None


# ---------------------------------------------------------------------------
# the volume integrator itself
# ---------------------------------------------------------------------------
def test_submerged_moments_of_a_half_immersed_box_are_exact() -> None:
    moments = submerged_moments(box(length=8.0, beam=4.0, draft=2.0, freeboard=2.0), 0.0)
    assert isinstance(moments, SubmergedMoments)
    assert moments.volume == pytest.approx(8.0 * 4.0 * 2.0, rel=1e-15)
    assert moments.centroid[0] == pytest.approx(4.0, abs=1e-15)
    assert moments.centroid[1] == pytest.approx(0.0, abs=1e-15)
    assert moments.centroid[2] == pytest.approx(-1.0, abs=1e-15)
    wetted = 8.0 * 4.0 + 2.0 * 8.0 * 2.0 + 2.0 * 4.0 * 2.0
    assert moments.wetted_surface_area == pytest.approx(wetted, rel=1e-15)


def test_submerged_moments_follow_the_waterline() -> None:
    mesh = box(length=8.0, beam=4.0, draft=2.0, freeboard=2.0)
    half = submerged_moments(mesh, 0.0).volume
    quarter = submerged_moments(mesh, -1.0).volume
    assert quarter == pytest.approx(0.5 * half, rel=1e-15)
    assert submerged_moments(mesh, 5.0).volume == pytest.approx(8.0 * 4.0 * 4.0, rel=1e-15)
    assert submerged_moments(mesh, -5.0).volume == pytest.approx(0.0, abs=1e-15)


def test_the_heeled_body_agrees_with_initial_stability_on_a_wetted_mesh() -> None:
    """``righting_arm`` must not touch the waterplane machinery ``GM`` is built on.

    ``initial_stability`` takes the *wetted* surface (no face straddling the waterline),
    while ``righting_arm`` clips. Handing the wetted box to both is the one mesh they can
    share, and it makes the two computations comparable at zero heel: ``KB``, ``BM`` and
    the displaced volume must agree with the analytic box formulas they are checked
    against elsewhere.
    """
    wetted = box(length=8.0, beam=4.0, draft=2.0, freeboard=0.0)
    state = initial_stability(wetted, kg=1.2)
    assert state.kb == pytest.approx(1.0, abs=1e-12)
    assert state.bm_transverse == pytest.approx(16.0 / 24.0, rel=1e-12)
    assert state.displacement_volume == pytest.approx(64.0, rel=1e-12)
    moments = submerged_moments(wetted, 0.0)
    assert moments.volume == pytest.approx(state.displacement_volume, rel=1e-12)
    assert moments.centroid[2] - state.keel_z == pytest.approx(state.kb, abs=1e-12)


# ---------------------------------------------------------------------------
# failure modes
# ---------------------------------------------------------------------------
def test_righting_arm_rejects_a_non_finite_loading() -> None:
    with pytest.raises(StabilityError):
        righting_arm(box(**TALL), kg=float("nan"), heel_degrees=10.0)
    with pytest.raises(StabilityError):
        righting_arm(box(**TALL), kg=1.0, heel_degrees=float("inf"))


def test_righting_arm_rejects_an_impossible_displacement() -> None:
    mesh = box(**TALL)
    with pytest.raises(StabilityError):
        righting_arm(mesh, kg=1.0, heel_degrees=10.0, displacement=-1.0)
    with pytest.raises(StabilityError):
        # larger than the body itself: no waterline can displace it
        righting_arm(mesh, kg=1.0, heel_degrees=10.0, displacement=1e6)


def test_gz_curve_rejects_an_empty_angle_list() -> None:
    with pytest.raises(StabilityError):
        gz_curve(box(**TALL), kg=1.0, heel_angles_deg=())


def test_an_awash_body_has_no_positive_righting_arm() -> None:
    """A hull modelled as its wetted surface plus a waterline cap has zero freeboard.

    Heeling it drowns the deck, so the arm goes negative: the curve is a statement about
    the body it was given, and the body has to *have* freeboard to have a range of
    stability. Documented here because handing `gz_curve` a wetted-surface mesh is the
    easiest mistake to make with it.
    """
    wetted = box(length=8.0, beam=4.0, draft=2.0, freeboard=0.0)
    curve = gz_curve(wetted, kg=1.2, heel_angles_deg=(10.0, 20.0))
    assert all(arm.righting_arm < 0.0 for arm in curve.arms)
    with_freeboard = gz_curve(box(length=8.0, beam=4.0, draft=2.0, freeboard=2.0), kg=1.2,
                              heel_angles_deg=(10.0, 20.0))
    assert all(arm.righting_arm > 0.0 for arm in with_freeboard.arms)


# ---------------------------------------------------------------------------
# serialisation, which the example and the report both consume
# ---------------------------------------------------------------------------
def test_curve_and_arm_serialise_with_their_evidence() -> None:
    curve = gz_curve(box(**TALL), kg=TALL_KG, heel_angles_deg=(0, 30, 60))
    assert isinstance(curve, GZCurve)
    assert isinstance(curve.arms[0], RightingArm)
    payload = curve.to_dict()
    assert payload["schema"] == "tensorcad.stability.gz/1"
    assert payload["heel_angles_deg"] == (0.0, 30.0, 60.0)
    assert payload["max_volume_residual"] < 1e-9
    assert len(payload["arms"]) == 3
    assert payload["arms"][2]["volume_residual"] < 1e-9
    assert payload["vanishing_angle_deg"] is None  # positive to 60 deg, not sampled past it
    assert curve.sinkages[0] == pytest.approx(0.0, abs=1e-12)
