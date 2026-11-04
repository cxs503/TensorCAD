"""Tests for :mod:`tensorcad.righting_arm` -- large-angle righting arms.

The references are:
* the **wall-sided formula** ``GZ = sin(phi)(GM + BM tan(phi)**2 / 2)``, exact for a
  rectangular barge while the waterline still spans the full beam;
* an independent **two-dimensional** computation for the same barge: the section is a
  rectangle, so the heeled problem reduces exactly to clipping a 2D polygon, which
  reaches the lever through the shoelace formula instead of through a surface flux.
  This is the reference used beyond the wall-sided regime, where the closed form has
  no claim;
* the **waterplane-moment** route to ``GM`` (``tensorcad.stability``) for the slope of a
  hull's curve at zero heel, which is a genuinely different integration.
"""

from __future__ import annotations

import math

import pytest

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.righting_arm import (
    GRAVITY,
    ImmersedVolume,
    RightingArm,
    RightingArmCurve,
    RightingArmError,
    box_reference_arm,
    heeled_mesh,
    immersed_volume,
    righting_arm,
    righting_arms,
    wall_sided_righting_arm,
    wall_sided_validity,
    wall_sided_vanishing_angle,
    zero_crossings,
)
from tensorcad.stability import waterplane_moments

LENGTH, BEAM, DRAFT, DEPTH = 10.0, 4.0, 2.0, 6.0


# ---------------------------------------------------------------------------
# reference geometries
# ---------------------------------------------------------------------------
def box(length: float, beam: float, depth: float, draft: float) -> Mesh:
    """Closed rectangular prism: keel on ``z = -draft``, top on ``z = depth - draft``."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, depth - draft):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi = -0.5 * beam, 0.5 * beam
    bottom, top = -draft, depth - draft
    quad((0.0, lo, bottom), (0.0, hi, bottom), (length, hi, bottom), (length, lo, bottom))
    quad((0.0, lo, top), (length, lo, top), (length, hi, top), (0.0, hi, top))
    quad((0.0, lo, bottom), (length, lo, bottom), (length, lo, top), (0.0, lo, top))
    quad((0.0, hi, bottom), (0.0, hi, top), (length, hi, top), (length, hi, bottom))
    quad((0.0, lo, bottom), (0.0, lo, top), (0.0, hi, top), (0.0, hi, bottom))
    quad((length, lo, bottom), (length, hi, bottom), (length, hi, top), (length, lo, top))
    return builder.build()


def flared_hull(nx: int = 41, nz: int = 13) -> Mesh:
    """A hull with topsides: generated deeper (draft 1.6) and floated at ``z = -0.9``.

    The shipped generators cap the surface at the design waterline, which leaves no
    freeboard and therefore no large-angle curve at all; generating the deeper body and
    floating it lower is the documented way to get one.
    """
    parameters = HullParameters(length=10.0, beam=2.0, draft=1.6, nx=nx, nz=nz)
    return generate_hull(parameters).mesh


# ---------------------------------------------------------------------------
# independent two-dimensional reference for a rectangular prism
# ---------------------------------------------------------------------------
def _clip_below(polygon, waterline):
    clipped = []
    count = len(polygon)
    for index in range(count):
        p = polygon[index]
        q = polygon[(index + 1) % count]
        dp = p[1] - waterline
        dq = q[1] - waterline
        if dp <= 0.0:
            clipped.append(p)
        if (dp < 0.0 < dq) or (dq < 0.0 < dp):
            t = dp / (dp - dq)
            clipped.append((p[0] + t * (q[0] - p[0]), waterline))
    return clipped


def _polygon_area_and_centroid(polygon) -> tuple[float, float, float]:
    area2 = cx = cy = 0.0
    count = len(polygon)
    for index in range(count):
        x1, y1 = polygon[index]
        x2, y2 = polygon[(index + 1) % count]
        cross = x1 * y2 - x2 * y1
        area2 += cross
        cx += (x1 + x2) * cross
        cy += (y1 + y2) * cross
    area = 0.5 * area2
    return abs(area), cx / (6.0 * area), cy / (6.0 * area)


def section_reference_arm(heel_deg: float, *, beam: float, draft: float, depth: float, kg: float) -> float:
    """``GZ`` of a rectangular prism by clipping its section in 2D (independent route)."""
    angle = math.radians(heel_deg)
    cos, sin = math.cos(angle), math.sin(angle)
    corners = [
        (-0.5 * beam, -draft),
        (0.5 * beam, -draft),
        (0.5 * beam, depth - draft),
        (-0.5 * beam, depth - draft),
    ]
    rotated = [(y * cos + z * sin, -y * sin + z * cos) for y, z in corners]
    target = beam * draft
    lower = min(point[1] for point in rotated)
    upper = max(point[1] for point in rotated)
    polygon = rotated
    for _ in range(200):
        waterline = 0.5 * (lower + upper)
        clipped = _clip_below(rotated, waterline)
        if len(clipped) < 3:
            break
        area, y_bar, _ = _polygon_area_and_centroid(clipped)
        polygon = clipped
        if abs(area - target) <= 1e-15 * target:
            break
        if area > target:
            upper = waterline
        else:
            lower = waterline
    _, y_bar, _ = _polygon_area_and_centroid(polygon)
    return y_bar - (-draft + kg) * sin


# ---------------------------------------------------------------------------
# immersed volume
# ---------------------------------------------------------------------------
def test_immersed_volume_of_a_box_matches_the_analytic_slab() -> None:
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    for waterline in (-1.5, -1.0, 0.0, 1.0):
        immersed = immersed_volume(mesh, waterline)
        expected = LENGTH * BEAM * (waterline + DRAFT)
        assert immersed.volume == pytest.approx(expected, rel=1e-14)
        # centroid of the slab, exactly
        assert immersed.centre_of_buoyancy[0] == pytest.approx(0.5 * LENGTH, rel=1e-14)
        assert immersed.centre_of_buoyancy[1] == pytest.approx(0.0, abs=1e-15)
        assert immersed.centre_of_buoyancy[2] == pytest.approx(0.5 * (waterline - DRAFT), rel=1e-12)


def test_immersed_volume_agrees_with_compute_hydrostatics_on_a_cut_hull() -> None:
    """On a surface capped at its waterline the two integrators must give the same volume.

    ``compute_hydrostatics`` sums signed tetrahedra of the *whole* closed surface, which
    is only the submerged volume when the mesh is pre-cut at the waterline; this module
    integrates only what is below the plane. On the shipped hull both must agree.
    """
    from tensorcad.hydrostatics import compute_hydrostatics

    mesh = generate_hull().mesh
    waterline = mesh.bounding_box().maximum[2]
    ours = immersed_volume(mesh, waterline)
    theirs = compute_hydrostatics(mesh, waterline, closed=True)
    assert ours.volume == pytest.approx(theirs.displacement_volume, rel=1e-12)
    # The centroid of the immersed volume is the centre of buoyancy.
    assert ours.centre_of_buoyancy[2] == pytest.approx(theirs.vertical_centre_buoyancy, abs=1e-12)
    assert ours.centre_of_buoyancy[0] == pytest.approx(theirs.longitudinal_centre_buoyancy, abs=1e-9)


def test_immersed_volume_rejects_an_open_surface() -> None:
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    holed = Mesh(vertices=mesh.vertices, faces=mesh.faces[:-1])
    with pytest.raises(RightingArmError, match="not closed"):
        immersed_volume(holed, 0.0)


def test_immersed_volume_rejects_a_waterline_outside_the_body() -> None:
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    with pytest.raises(RightingArmError, match="immersed volume"):
        immersed_volume(mesh, -DRAFT - 1.0)
    # Exactly at the keel the region is a plane with no volume, so there is no centroid
    # to report; refusing beats dividing by zero.
    with pytest.raises(RightingArmError, match="immersed volume"):
        immersed_volume(mesh, -DRAFT)


def test_immersed_volume_rejects_a_non_finite_waterline() -> None:
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    with pytest.raises(RightingArmError, match="finite"):
        immersed_volume(mesh, math.nan)


def test_immersed_volume_counts_the_faces_the_plane_actually_cuts() -> None:
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    # A mid-height waterline cuts the four sides (two triangles each) and leaves the two
    # bottom triangles below it whole: ten contributing faces, none of them the top.
    assert immersed_volume(mesh, 1.0).clipped_faces == 10
    # A waterline at the very top only touches the top face, which contributes nothing
    # to the flux (every integrand carries a factor z - w) and is skipped as degenerate.
    assert immersed_volume(mesh, DEPTH - DRAFT).volume == pytest.approx(
        LENGTH * BEAM * DEPTH, rel=1e-14
    )


# ---------------------------------------------------------------------------
# the rotation convention
# ---------------------------------------------------------------------------
def test_positive_heel_puts_starboard_down() -> None:
    mesh = Mesh(vertices=((0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, 0.0)), faces=((0, 1, 2),))
    heeled = heeled_mesh(mesh, 90.0)
    # The starboard point (y = 1) goes to z = -1: starboard down.
    assert heeled.vertices[0] == pytest.approx((0.0, 0.0, -1.0), abs=1e-15)
    # The masthead point (z = 1) goes to y = +1: the body leans to starboard.
    assert heeled.vertices[1] == pytest.approx((0.0, 1.0, 0.0), abs=1e-15)


def test_heeled_mesh_preserves_the_body() -> None:
    """A rotation moves the body, it does not change how much of it there is.

    (What a fixed *waterline* sees does change, which is the whole point of the draft
    solve: the displaced volume below a plane is not rotation invariant.)
    """
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    heeled = heeled_mesh(mesh, 37.0)
    top = heeled.bounding_box().maximum[2]
    assert immersed_volume(heeled, top).volume == pytest.approx(
        LENGTH * BEAM * DEPTH, rel=1e-13
    )
    with pytest.raises(RightingArmError, match="finite"):
        heeled_mesh(mesh, math.inf)


# ---------------------------------------------------------------------------
# the box against the wall-sided formula
# ---------------------------------------------------------------------------
WALL_SIDED_ANGLES = (2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 44.0)


@pytest.mark.parametrize("kg", (0.6, 1.0, 1.9))
def test_box_reproduces_the_wall_sided_formula_to_round_off(kg: float) -> None:
    """Inside its regime the closed form is exact, so the mesh must be too."""
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    curve = righting_arms(mesh, kg=kg, heel_angles=WALL_SIDED_ANGLES)
    validity = wall_sided_validity(beam=BEAM, draft=DRAFT, depth=DEPTH)
    assert validity.binding_deg == pytest.approx(45.0)  # bilge and deck edge coincide
    for arm in curve.arms:
        reference = box_reference_arm(
            arm.heel_deg, length=LENGTH, beam=BEAM, draft=DRAFT, kg=kg
        )
        assert arm.gz == pytest.approx(reference, abs=1e-12)


def test_wall_sided_formula_leaves_its_regime_and_the_module_does_not() -> None:
    """The reference stops being exact past the binding angle; the module keeps working.

    The check is two-sided on purpose: the departure of the closed form is asserted (so
    the limit is documented rather than discovered later), and the module's own value is
    pinned to the independent 2D section computation instead.
    """
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    kg = 1.9
    validity = wall_sided_validity(beam=BEAM, draft=DRAFT, depth=DEPTH)
    beyond = validity.binding_deg + 15.0
    arm = righting_arm(mesh, kg=kg, heel_deg=beyond)
    closed_form = box_reference_arm(beyond, length=LENGTH, beam=BEAM, draft=DRAFT, kg=kg)
    assert abs(arm.gz - closed_form) > 1e-2  # the formula has no claim here
    reference = section_reference_arm(beyond, beam=BEAM, draft=DRAFT, depth=DEPTH, kg=kg)
    assert arm.gz == pytest.approx(reference, abs=1e-12)


@pytest.mark.parametrize("heel", (1.0, 5.0, 20.0, 44.0, 45.0, 46.0, 55.0, 60.0, 75.0, 90.0))
def test_box_arm_matches_the_independent_section_computation(heel: float) -> None:
    """All heel angles, inside and outside the wall-sided regime."""
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    kg = 1.9
    arm = righting_arm(mesh, kg=kg, heel_deg=heel)
    reference = section_reference_arm(heel, beam=BEAM, draft=DRAFT, depth=DEPTH, kg=kg)
    assert arm.gz == pytest.approx(reference, abs=1e-12)


def test_every_heel_holds_the_displacement() -> None:
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    curve = righting_arms(mesh, kg=0.9, heel_angles=tuple(float(a) for a in range(0, 121, 10)))
    for arm in curve.arms:
        assert arm.displaced_volume == pytest.approx(arm.target_volume, rel=1e-12)
        assert abs(arm.relative_volume_error) < 1e-12


def test_zero_heel_has_no_arm_and_the_curve_is_odd() -> None:
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    angles = (-40.0, -20.0, -10.0, 0.0, 10.0, 20.0, 40.0)
    curve = righting_arms(mesh, kg=0.8, heel_angles=angles)
    assert curve.arm_at(0.0).gz == pytest.approx(0.0, abs=1e-15)
    for heel in (-40.0, -20.0, -10.0):
        assert curve.arm_at(heel).gz == pytest.approx(-curve.arm_at(-heel).gz, abs=1e-13)


def test_raising_g_shifts_the_arm_by_exactly_kg_sine_phi() -> None:
    """``GZ(phi; KG + d) = GZ(phi; KG) - d sin(phi)``: B does not depend on the loading."""
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    angles = (10.0, 20.0, 30.0, 45.0, 60.0)
    base = righting_arms(mesh, kg=0.8, heel_angles=angles)
    raised = righting_arms(mesh, kg=1.3, heel_angles=angles)
    for arm, other in zip(base.arms, raised.arms):
        shift = 0.5 * math.sin(math.radians(arm.heel_deg))
        assert other.gz == pytest.approx(arm.gz - shift, abs=1e-13)


def test_vanishing_angle_of_an_unstable_body_matches_the_closed_form() -> None:
    """With ``GM < 0`` the body lolls: the wall-sided formula predicts the resting heel."""
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    kg = 1.9
    kb, bm = 0.5 * DRAFT, (LENGTH * BEAM**3 / 12.0) / (LENGTH * BEAM * DRAFT)
    closed = wall_sided_vanishing_angle(kb=kb, bm=bm, kg=kg)
    assert closed is not None and closed < wall_sided_validity(
        beam=BEAM, draft=DRAFT, depth=DEPTH
    ).binding_deg
    curve = righting_arms(mesh, kg=kg, heel_angles=tuple(0.5 * i for i in range(0, 181)))
    crossings = zero_crossings(list(curve.heel_deg), list(curve.gz))
    rising = [angle for angle in crossings if angle > 0.0]
    assert rising, crossings
    assert rising[0] == pytest.approx(closed, abs=0.02)


def test_area_under_the_curve_matches_the_closed_form_integral() -> None:
    """``int GZ dphi = GM (1 - cos phi) + BM/2 (cos phi + sec phi - 2)``."""
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    kg, limit = 1.0, 45.0
    kb, bm = 0.5 * DRAFT, (LENGTH * BEAM**3 / 12.0) / (LENGTH * BEAM * DRAFT)
    gm = kb + bm - kg
    phi = math.radians(limit)
    closed = gm * (1.0 - math.cos(phi)) + 0.5 * bm * (math.cos(phi) + 1.0 / math.cos(phi) - 2.0)
    steps = tuple(0.25 * i for i in range(int(limit / 0.25) + 1))
    curve = righting_arms(mesh, kg=kg, heel_angles=steps)
    assert curve.area(limit) == pytest.approx(closed, abs=2e-6)


# ---------------------------------------------------------------------------
# curve bookkeeping
# ---------------------------------------------------------------------------
def test_curve_reports_the_maximum_and_the_range_on_its_own_grid() -> None:
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    curve = righting_arms(mesh, kg=1.0, heel_angles=(0.0, 30.0, 60.0, 90.0, 120.0))
    assert curve.max_gz == max(curve.gz)
    assert curve.angle_of_max_gz in curve.heel_deg
    assert curve.range_of_stability() is None  # a barge with GM > 0 does not capsize here
    assert curve.area() > 0.0


def test_gz_at_interpolates_and_refuses_to_invent_grid_points() -> None:
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    curve = righting_arms(mesh, kg=1.0, heel_angles=(0.0, 10.0, 20.0))
    assert curve.gz_at(5.0) == pytest.approx(0.5 * (curve.gz_at(0.0) + curve.gz_at(10.0)))
    assert curve.gz_at(-5.0) == curve.gz_at(0.0)  # clamped, not extrapolated
    assert curve.gz_at(90.0) == curve.gz_at(20.0)
    with pytest.raises(RightingArmError, match="not on this curve"):
        curve.arm_at(5.0)


def test_arm_carries_the_loading_and_the_moment() -> None:
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    arm = righting_arm(mesh, kg=1.0, heel_deg=30.0)
    assert isinstance(arm, RightingArm)
    assert arm.heel_rad == pytest.approx(math.radians(30.0))
    assert arm.displaced_mass == pytest.approx(arm.displaced_volume * 1025.0)
    assert arm.righting_moment == pytest.approx(arm.displaced_mass * GRAVITY * arm.gz / 1000.0)
    assert arm.iterations >= 1


def test_curve_and_arm_serialise() -> None:
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    curve = righting_arms(mesh, kg=1.0, heel_angles=(0.0, 30.0, 60.0))
    payload = curve.to_dict()
    assert payload["schema"] == "tensorcad.righting_arm.curve/1"
    assert len(payload["arms"]) == 3
    assert payload["arms"][1]["schema"] == "tensorcad.righting_arm.arm/1"
    assert payload["max_relative_volume_error"] < 1e-12
    immersed = immersed_volume(mesh, 0.0)
    assert isinstance(immersed, ImmersedVolume)
    assert immersed.to_dict()["schema"] == "tensorcad.righting_arm.immersed/1"


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------
def test_parameters_are_validated() -> None:
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    with pytest.raises(RightingArmError, match="kg must be finite"):
        righting_arms(mesh, kg=math.nan, heel_angles=(0.0,))
    with pytest.raises(RightingArmError, match="must not be empty"):
        righting_arms(mesh, kg=1.0, heel_angles=())
    with pytest.raises(RightingArmError, match="duplicates"):
        righting_arms(mesh, kg=1.0, heel_angles=(0.0, 10.0, 10.0))
    with pytest.raises(RightingArmError, match="outside"):
        righting_arms(mesh, kg=1.0, heel_angles=(190.0,))
    with pytest.raises(RightingArmError, match="finite"):
        righting_arms(mesh, kg=1.0, heel_angles=(math.nan,))
    with pytest.raises(RightingArmError, match="not closed"):
        holed = Mesh(vertices=mesh.vertices, faces=mesh.faces[:-1])
        righting_arms(holed, kg=1.0, heel_angles=(10.0,))
    with pytest.raises(RightingArmError, match="not positive"):
        righting_arms(mesh, kg=1.0, heel_angles=(10.0,), displacement=-1.0)


def test_a_body_capped_at_its_waterline_has_no_righting_arm() -> None:
    """No freeboard means the displacement is the whole volume, and the arm is zero.

    ``generate_hull`` caps the surface at the design waterline, so this is the shipped
    default geometry -- and the honest answer is that a large-angle curve does not exist
    for it, not a row of zeros.
    """
    mesh = generate_hull().mesh
    with pytest.raises(RightingArmError, match="freeboard"):
        righting_arms(mesh, kg=0.3, heel_angles=(0.0, 10.0))


def test_wall_sided_validity_reports_the_binding_edge() -> None:
    bilge_bound = wall_sided_validity(beam=4.0, draft=1.0, depth=3.0)
    assert bilge_bound.binding == "bilge emergence"
    assert bilge_bound.binding_deg == pytest.approx(math.degrees(math.atan2(2.0, 4.0)))
    deck_bound = wall_sided_validity(beam=4.0, draft=2.0, depth=2.5)
    assert deck_bound.binding == "deck immersion"
    assert deck_bound.binding_deg == pytest.approx(math.degrees(math.atan2(1.0, 4.0)))
    with pytest.raises(RightingArmError, match="exceeds depth"):
        wall_sided_validity(beam=4.0, draft=3.0, depth=2.0)
    with pytest.raises(RightingArmError, match="positive"):
        wall_sided_validity(beam=0.0, draft=1.0, depth=2.0)


def test_wall_sided_reference_helpers() -> None:
    assert wall_sided_righting_arm(0.0, kb=1.0, bm=0.5, kg=0.8) == 0.0
    # Small angles: the slope of the arm at zero heel is GM.
    delta = 1e-6
    slope = wall_sided_righting_arm(math.degrees(delta), kb=1.0, bm=0.5, kg=0.8) / delta
    assert slope == pytest.approx(0.7, rel=1e-9)
    assert wall_sided_vanishing_angle(kb=1.0, bm=0.5, kg=0.8) is None


def test_zero_crossings_finds_both_directions() -> None:
    # down at 5 deg (0.1 -> -0.1), up at 10 + 10/3 deg (-0.1 -> 0.2)
    assert zero_crossings([0.0, 10.0, 20.0], [0.1, -0.1, 0.2]) == pytest.approx((5.0, 40.0 / 3.0))
    assert zero_crossings([0.0, 10.0], [0.0, 0.5]) == pytest.approx((0.0,))
    assert zero_crossings([0.0, 10.0], [0.5, 0.2]) == ()
    with pytest.raises(RightingArmError, match="same length"):
        zero_crossings([0.0, 10.0], [0.5])


# ---------------------------------------------------------------------------
# a hull with topsides
# ---------------------------------------------------------------------------
def test_hull_arm_slope_reproduces_the_waterplane_metacentric_height() -> None:
    """``dGZ/dphi`` at zero heel is ``GM``, reached here by two different integrations.

    The 3D heeled-displacement route of this module and the waterplane-moment route of
    :mod:`tensorcad.stability` (Green's theorem on the plane section) share no code, so
    their agreement is a cross-check and not a tautology. The ratio ``GZ/sin(phi)`` tends
    to ``GM`` from above with an ``O(phi**2)`` error, which is removed by Richardson.
    """
    mesh = flared_hull()
    waterline, kg = -0.9, 0.30
    immersed = immersed_volume(mesh, waterline)
    section = waterplane_moments(mesh, waterline)
    keel = min(vertex[2] for vertex in mesh.vertices)
    kb = immersed.centre_of_buoyancy[2] - keel
    bm = section.transverse_second_moment / immersed.volume
    gm = kb + bm - kg

    ratios = {
        heel: righting_arm(mesh, kg=kg, heel_deg=heel, waterline_z=waterline).gz
        / math.sin(math.radians(heel))
        for heel in (2.0, 4.0)
    }
    slope = (4.0**2 * ratios[2.0] - 2.0**2 * ratios[4.0]) / (4.0**2 - 2.0**2)
    assert slope == pytest.approx(gm, rel=5e-5)


def test_hull_curve_converges_second_order_with_the_mesh() -> None:
    """The arm at a fixed heel is second order in the grid step, measured not asserted."""
    angles = (10.0, 20.0)
    arms = {}
    for nx, nz in ((41, 13), (81, 25), (161, 41)):
        curve = righting_arms(
            flared_hull(nx, nz), kg=0.30, heel_angles=angles, waterline_z=-0.9
        )
        arms[(nx, nz)] = curve.gz

    for index in range(len(angles)):
        coarse = arms[(41, 13)][index]
        medium = arms[(81, 25)][index]
        fine = arms[(161, 41)][index]
        assert coarse != medium != fine
        order = math.log2(abs(coarse - fine) / abs(medium - fine))
        assert order > 1.8, order
    # And the shipped default grid is worth reporting, not hiding: a few per cent.
    assert abs(arms[(41, 13)][1] - arms[(161, 41)][1]) / abs(arms[(161, 41)][1]) > 0.01


def test_hull_curve_is_positive_and_holds_the_displacement() -> None:
    curve = righting_arms(
        flared_hull(), kg=0.30, heel_angles=tuple(float(a) for a in range(0, 91, 10)), waterline_z=-0.9
    )
    assert isinstance(curve, RightingArmCurve)
    assert all(arm.gz > 0.0 for arm in curve.arms[1:])
    assert all(abs(arm.relative_volume_error) < 1e-12 for arm in curve.arms)
    assert curve.area(60.0) < curve.area(90.0)
    # Raising KG tilts the whole curve down by d sin(phi), which is the loading identity
    # again, now on a real hull rather than on a prism.
    raised = righting_arms(
        flared_hull(), kg=0.45, heel_angles=tuple(float(a) for a in range(0, 91, 10)), waterline_z=-0.9
    )
    for base_arm, other_arm in zip(curve.arms, raised.arms):
        shift = 0.15 * math.sin(math.radians(base_arm.heel_deg))
        assert other_arm.gz == pytest.approx(base_arm.gz - shift, abs=1e-12)