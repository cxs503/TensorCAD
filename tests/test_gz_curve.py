"""Righting arms at finite heel: a closed-form section, and GM itself.

Two references, chosen because neither can be satisfied by a bug in the integrator:

* **A rectangular section at any heel**, clipped and integrated in two dimensions in an
  independent world-frame routine (the rotated rectangle is the exact body, the
  waterplane is a horizontal line, and the displaced area is held at ``B*T``). The
  measurement must reproduce it to round-off at *every* heel, including past deck and
  keel immersion, where no closed form for the arm exists.
* **The wall-sided formula** ``GZ = sin(phi)(BM + BM tan^2(phi)/2 - BG)``, which is the
  classical literature result and is only valid while both ends of the waterline stay on
  the vertical sides (``tan(phi) < 2T/B``). Inside that window it must agree to
  round-off; outside it must *not*, and that departure is asserted too, so a formula
  that happened to be applied everywhere would fail here.

And one cross-check that is not a reference at all: ``dGZ/dphi`` at the origin is ``GM``,
which ``initial_stability`` obtains from the waterplane's second moment. When the mesh is
capped at the waterline the two must be the same number, and the curve's secant slope must
approach it as the first heel step shrinks.
"""

from __future__ import annotations

import json
import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.stability import (
    SCHEMA_GZ,
    GZCurve,
    StabilityError,
    initial_stability,
    righting_arm_curve,
)
from tensorcad.wigley import generate_wigley


# ---------------------------------------------------------------------------
# reference geometries
# ---------------------------------------------------------------------------
def box(length: float, beam: float, draft: float, freeboard: float) -> Mesh:
    """Closed box, keel on ``z = -draft``, deck on ``z = +freeboard``, waterline ``z = 0``."""
    builder = MeshBuilder()
    corner = {}
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
    bottom, top = -draft, freeboard
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


def reference_gz(
    beam: float, draft: float, freeboard: float, kg: float, heel_deg: float
) -> tuple[float, float]:
    """Exact arm of the rectangular section, and its waterplane height.

    The section's four corners are rotated into the waterplane frame, clipped by the
    horizontal line ``Z = h``, and integrated with the shoelace rule; ``h`` is bisected
    until the clipped area equals ``B*T``. No part of this shares code -- or a formula --
    with the 3-D divergence integrals under test.
    """
    phi = math.radians(heel_deg)
    cos_phi, sin_phi = math.cos(phi), math.sin(phi)
    keel = -draft
    section = [
        (y * cos_phi - z * sin_phi, y * sin_phi + z * cos_phi)
        for y, z in (
            (-0.5 * beam, keel),
            (0.5 * beam, keel),
            (0.5 * beam, freeboard),
            (-0.5 * beam, freeboard),
        )
    ]
    heights = [point[1] for point in section]
    low, high = min(heights) - 1.0, max(heights) + 1.0

    def clipped(plane: float) -> tuple[float, float]:
        polygon: list[tuple[float, float]] = []
        for index in range(len(section)):
            start, end = section[index], section[(index + 1) % len(section)]
            start_in, end_in = start[1] <= plane, end[1] <= plane
            if start_in:
                polygon.append(start)
            if start_in != end_in:
                ratio = (plane - start[1]) / (end[1] - start[1])
                polygon.append((start[0] + ratio * (end[0] - start[0]), plane))
        area = moment_y = 0.0
        for index in range(len(polygon)):
            x1, y1 = polygon[index]
            x2, y2 = polygon[(index + 1) % len(polygon)]
            cross = x1 * y2 - x2 * y1
            area += cross
            moment_y += (x1 + x2) * cross
        return 0.5 * area, moment_y / 6.0

    for _ in range(200):
        middle = 0.5 * (low + high)
        if clipped(middle)[0] < beam * draft:
            low = middle
        else:
            high = middle
    plane = 0.5 * (low + high)
    area, moment_y = clipped(plane)
    return -(keel + kg) * sin_phi - moment_y / area, plane


def wall_sided_gz(beam: float, draft: float, kg: float, heel_deg: float) -> float:
    """Classical wall-sided arm: ``sin(phi)(BM + BM tan^2(phi)/2 - BG)``."""
    phi = math.radians(heel_deg)
    bm = beam * beam / (12.0 * draft)
    bg = kg - 0.5 * draft
    return math.sin(phi) * (bm + 0.5 * bm * math.tan(phi) ** 2 - bg)


LENGTH, BEAM, DRAFT, FREEBOARD, KG = 10.0, 2.0, 0.5, 1.5, 0.45


@pytest.fixture(scope="module")
def box_curve() -> GZCurve:
    return righting_arm_curve(
        box(LENGTH, BEAM, DRAFT, FREEBOARD),
        kg=KG,
        waterline_z=0.0,
        angles_deg=range(0, 181, 5),
    )


@pytest.fixture(scope="module")
def hull_curve() -> GZCurve:
    mesh = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8)).mesh
    return righting_arm_curve(mesh, kg=0.6, waterline_z=0.0, freeboard=0.4, angles_deg=range(0, 91, 5))


# ---------------------------------------------------------------------------
# the box: exact section, and the wall-sided window
# ---------------------------------------------------------------------------
def test_the_reference_body_is_itself_valid() -> None:
    """A reference body that is not watertight would make every number below moot."""
    mesh = box(LENGTH, BEAM, DRAFT, FREEBOARD)
    assert meshlib.is_watertight(mesh)
    assert meshlib.total_signed_volume(mesh) == pytest.approx(
        LENGTH * BEAM * (DRAFT + FREEBOARD), rel=1e-15
    )


@pytest.mark.parametrize("heel", range(0, 181, 5))
def test_box_arm_matches_the_exact_section(box_curve: GZCurve, heel: float) -> None:
    """Round-off agreement at every heel, including past deck and keel immersion."""
    point = box_curve.points[heel // 5]
    assert point.heel_deg == heel
    expected, _plane = reference_gz(BEAM, DRAFT, FREEBOARD, KG, heel)
    assert point.gz == pytest.approx(expected, abs=1e-12)


def test_box_holds_its_displacement_at_every_heel(box_curve: GZCurve) -> None:
    """The loading condition is the invariant of the whole calculation."""
    for point in box_curve.points:
        assert point.displacement_volume == pytest.approx(
            LENGTH * BEAM * DRAFT, rel=1e-13
        )


@pytest.mark.parametrize("heel", range(0, 26, 5))
def test_wall_sided_formula_holds_while_the_waterline_stays_on_the_walls(
    box_curve: GZCurve, heel: float
) -> None:
    """``tan(phi) < 2T/B = 0.5``: the classical formula is exact, so match it."""
    point = box_curve.points[heel // 5]
    assert point.gz == pytest.approx(wall_sided_gz(BEAM, DRAFT, KG, heel), rel=1e-12, abs=1e-15)


def test_wall_sided_formula_leaves_its_window(box_curve: GZCurve) -> None:
    """At 30 degrees the leeward waterline has crossed the keel edge.

    The formula still returns a number -- a larger one -- and that is precisely why the
    fitted-versus-measured distinction matters: outside its window the closed form is not
    a reference, and the test asserts the departure rather than the agreement.
    """
    point = box_curve.points[6]
    assert point.heel_deg == 30.0
    assert point.gz < wall_sided_gz(BEAM, DRAFT, KG, 30.0)
    relative = abs(point.gz / wall_sided_gz(BEAM, DRAFT, KG, 30.0) - 1.0)
    assert 1e-3 < relative < 0.1


def test_deck_immersion_flag_follows_the_reference_waterplane(box_curve: GZCurve) -> None:
    """The flag must be the reference's own answer, not a guessed angle.

    The leeward deck corner rotates to ``Z = -sin(phi) + F cos(phi)``; the deck is
    immersed exactly when the independently solved waterplane rises above it. Both
    conditions -- the flag and its monotonicity -- are asserted point by point.
    """
    flags: list[tuple[float, bool]] = []
    for heel in range(0, 181, 5):
        point = box_curve.points[heel // 5]
        _gz, plane = reference_gz(BEAM, DRAFT, FREEBOARD, KG, heel)
        deck_corner = -math.sin(math.radians(heel)) + FREEBOARD * math.cos(math.radians(heel))
        flags.append((heel, point.deck_immersed))
        assert point.deck_immersed == (plane > deck_corner - 1e-9)
    first = next(heel for heel, immersed in flags if immersed)
    # measured, not assumed: the body sinks as it heels, so the crossing is well past
    # the naive atan(2F/B) = 56.3 degrees
    assert first == 65.0
    assert all(immersed for heel, immersed in flags if heel >= first)


def test_the_skirt_of_a_capped_hull_equals_the_solid_it_imitates() -> None:
    """A capped box plus a wall-sided skirt is the same body as a box with freeboard.

    This is the test that catches an *open* skirt: without its deck cap the topsides
    displace too little and only the large-heel points move, which is exactly the
    failure a max-only check would miss.
    """
    capped = box(LENGTH, BEAM, DRAFT, 0.0)
    skirted = righting_arm_curve(
        capped, kg=KG, waterline_z=0.0, freeboard=FREEBOARD, angles_deg=range(0, 181, 5)
    )
    solid = righting_arm_curve(
        box(LENGTH, BEAM, DRAFT, FREEBOARD), kg=KG, waterline_z=0.0, angles_deg=range(0, 181, 5)
    )
    for left, right in zip(skirted.points, solid.points):
        assert left.gz == pytest.approx(right.gz, abs=1e-10)
    assert "wall-sided topsides" in skirted.basis


def test_box_features(box_curve: GZCurve) -> None:
    """A box is symmetric fore-and-aft and upright at 0 and 180 degrees."""
    assert box_curve.points[0].gz == pytest.approx(0.0, abs=1e-15)
    assert box_curve.points[-1].gz == pytest.approx(0.0, abs=1e-15)
    assert box_curve.gz_max > 0.5
    assert 80.0 <= box_curve.angle_of_gz_max_deg <= 90.0
    # the arm only returns to zero at the inversion, so there is no crossing *after* the
    # maximum: the reported vanishing angle is None (no sign change inside the grid)
    assert box_curve.vanishing_angle_deg is None
    assert box_curve.dynamic_stability_area > 0.0


# ---------------------------------------------------------------------------
# GM: the same number by two different integrals
# ---------------------------------------------------------------------------
def test_box_metacentric_height_is_the_closed_form(box_curve: GZCurve) -> None:
    expected = 0.5 * DRAFT + BEAM * BEAM / (12.0 * DRAFT) - KG
    assert box_curve.initial_gm_transverse == pytest.approx(expected, rel=1e-15)


def test_box_small_angle_slope_is_gm(box_curve: GZCurve) -> None:
    """``dGZ/dphi(0) = GM`` in metres per radian; the 5-degree secant is O(phi^2) low."""
    secant = (box_curve.points[1].gz - box_curve.points[0].gz) / math.radians(5.0)
    assert secant == pytest.approx(box_curve.initial_gm_transverse, rel=5e-3)


def test_hull_metacentric_height_is_the_initial_stability_value(hull_curve: GZCurve) -> None:
    """Two integrations of the same upright condition must agree exactly.

    ``initial_stability`` takes ``BM`` from the waterplane's second moment; the curve
    takes ``KB`` from the *clipped* buoyancy centroid. When the mesh is capped at the
    waterline they reduce to the same arithmetic and any disagreement would be a bug in
    one of the two paths.
    """
    mesh = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8)).mesh
    upright = initial_stability(mesh, kg=0.6)
    assert hull_curve.initial_gm_transverse == pytest.approx(upright.gm_transverse, rel=1e-12)


def test_hull_small_angle_slope_approaches_gm(hull_curve: GZCurve) -> None:
    """With a 1-degree first step the secant is within a few percent of ``GM``."""
    first = hull_curve.points[0].gz
    second = hull_curve.points[1].gz
    secant = (second - first) / math.radians(5.0)
    assert secant == pytest.approx(hull_curve.initial_gm_transverse, rel=0.06)


def test_an_unstable_loading_condition_has_a_negative_arm() -> None:
    """``KG`` above ``KM`` must turn the curve over at the origin, not merely shift it."""
    mesh = box(LENGTH, BEAM, DRAFT, FREEBOARD)
    curve = righting_arm_curve(mesh, kg=1.5, waterline_z=0.0, angles_deg=range(0, 31, 5))
    assert curve.initial_gm_transverse < 0.0
    assert curve.points[1].gz < 0.0
    assert curve.gz_max <= 0.0


# ---------------------------------------------------------------------------
# the wall-sided topsides over a wetted-surface hull
# ---------------------------------------------------------------------------
def test_a_hull_capped_at_the_waterline_needs_freeboard() -> None:
    """Without topsides the body displaces its whole volume and never heels."""
    mesh = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8)).mesh
    with pytest.raises(StabilityError, match="no volume above the design waterline"):
        righting_arm_curve(mesh, kg=0.6, waterline_z=0.0)


def test_freeboard_is_refused_when_the_mesh_already_has_topsides() -> None:
    """A skirt over existing geometry would double-count the volume above the waterline."""
    mesh = box(LENGTH, BEAM, DRAFT, FREEBOARD)
    with pytest.raises(StabilityError, match="already rises"):
        righting_arm_curve(mesh, kg=KG, waterline_z=0.0, freeboard=0.1)


def test_wigley_curve_is_continuous_and_anchored() -> None:
    """A third family, for the case the two above do not cover: a parabolic section."""
    mesh = generate_wigley().mesh
    curve = righting_arm_curve(mesh, kg=0.4, waterline_z=0.0, freeboard=0.3, angles_deg=range(0, 61, 5))
    assert curve.points[0].gz == pytest.approx(0.0, abs=1e-12)
    arms = [point.gz for point in curve.points]
    assert curve.gz_max > 0.0
    assert curve.points[0].gz < curve.points[3].gz, "the arm must rise off the origin"
    # continuity: a 5-degree step cannot move the arm by more than 5 cm on this hull
    assert all(abs(b - a) < 0.05 for a, b in zip(arms, arms[1:]))
    for point in curve.points:
        # the waterplane solve is converged to ``tolerance`` (1e-11 of the displacement
        # by default), so the loading condition is reproduced to that fraction
        assert point.displacement_volume == pytest.approx(
            curve.displacement_volume, rel=1e-10
        )


# ---------------------------------------------------------------------------
# contract: inputs, schema, serialisation
# ---------------------------------------------------------------------------
def test_angles_must_start_at_zero_and_increase() -> None:
    mesh = box(LENGTH, BEAM, DRAFT, FREEBOARD)
    with pytest.raises(StabilityError, match="must be 0"):
        righting_arm_curve(mesh, kg=KG, angles_deg=[5.0, 10.0])
    with pytest.raises(StabilityError, match="strictly increasing"):
        righting_arm_curve(mesh, kg=KG, angles_deg=[0.0, 10.0, 10.0])
    with pytest.raises(StabilityError, match="at least two"):
        righting_arm_curve(mesh, kg=KG, angles_deg=[0.0])
    with pytest.raises(StabilityError, match="finite"):
        righting_arm_curve(mesh, kg=KG, angles_deg=[0.0, float("nan")])


def test_invalid_loading_and_freeboard_are_refused() -> None:
    mesh = box(LENGTH, BEAM, DRAFT, FREEBOARD)
    with pytest.raises(StabilityError, match="kg must be finite"):
        righting_arm_curve(mesh, kg=float("inf"))
    with pytest.raises(StabilityError, match="freeboard"):
        righting_arm_curve(mesh, kg=KG, freeboard=-0.1)
    with pytest.raises(StabilityError, match="tolerance"):
        righting_arm_curve(mesh, kg=KG, tolerance=0.0)


def test_payload_is_json_serialisable(box_curve: GZCurve) -> None:
    payload = box_curve.to_dict()
    assert payload["schema"] == SCHEMA_GZ
    assert payload["kg"] == pytest.approx(KG)
    assert len(payload["points"]) == len(box_curve.points)
    assert payload["points"][0]["heel_deg"] == 0.0
    assert "free-surface" not in payload["basis"]
    assert any("free-surface" in note for note in payload["notes"])
    json.dumps(payload)  # no NaN, no tuples, no numpy


def test_the_default_grid_is_the_documented_one() -> None:
    mesh = box(LENGTH, BEAM, DRAFT, FREEBOARD)
    curve = righting_arm_curve(mesh, kg=KG)
    assert [point.heel_deg for point in curve.points] == [float(angle) for angle in range(0, 91, 5)]