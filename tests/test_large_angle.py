"""Large-angle stability: every reference here is outside the integrator that is tested.

Four independent yardsticks are used, and none of them is another implementation of
the same integral:

* the **wall-sided closed form** ``sin(phi)(GM + BM tan^2(phi)/2)`` for a box barge --
  exact for a rectangular section, so the agreement is a statement about the
  integrator, not about a tolerance;
* ``GM`` from :func:`tensorcad.stability.initial_stability`, i.e. the **slope of the
  curve at zero heel** measured by this module against a metacentric height measured
  by the upright-waterplane module: two code paths, no shared lines;
* **rotational symmetry**: a circular-section barge ("log") cannot move its centre of
  buoyancy sideways, so its righting arm is exactly ``(R - KG) sin(phi)`` at every
  heel angle, which exercises a curved surface rather than a polyhedron;
* the **closed-form integral** of the wall-sided arm for the area under a curve.

The waterline-neutrality of :func:`tensorcad.large_angle.extend_freeboard` is checked
against :func:`tensorcad.hydrostatics.compute_hydrostatics`, and the failure modes that
would silently produce a number (open mesh, inward mesh, a hull capped at its own
design waterline) are checked to raise.
"""
from __future__ import annotations

import math
from itertools import pairwise

import pytest

import tensorcad
from tensorcad import hydrostatics as hydrostatics_module
from tensorcad import mesh as meshlib
from tensorcad import stability as stability_module
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.large_angle import (
    GzCurve,
    HeeledHydrostatics,
    LargeAngleError,
    RightingArm,
    box_dynamical_stability,
    box_hull_mesh,
    box_righting_arm,
    dynamical_stability,
    extend_freeboard,
    gz_curve,
    heeled_hydrostatics,
    righting_arm,
    wall_sided_limit_deg,
)

LENGTH = 6.0
BREADTH = 2.0
DRAFT = 0.8
FREEBOARD = 0.6
KG = 0.5
HULL_KG = 0.5


@pytest.fixture(scope="module")
def box_mesh() -> Mesh:
    """Box barge with enough freeboard for the wall-sided formula to stay valid to 30 deg."""
    return box_hull_mesh(LENGTH, BREADTH, DRAFT, FREEBOARD)


@pytest.fixture(scope="module")
def capped_hull() -> Mesh:
    """The shipped hull family, capped at its design waterline (no freeboard)."""
    return generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8)).mesh


@pytest.fixture(scope="module")
def hull_mesh(capped_hull: Mesh) -> Mesh:
    """The same hull with 0.6 m of freeboard added, which is what GZ needs."""
    return extend_freeboard(capped_hull, FREEBOARD)


@pytest.fixture(scope="module")
def hull_gz(hull_mesh: Mesh) -> GzCurve:
    return gz_curve(hull_mesh, kg=HULL_KG, up_to_deg=40.0, step_deg=10.0)


def log_prism_mesh(radius: float, length: float, facets: int) -> Mesh:
    """Circular-section barge as a prism: axis along ``x`` at ``z = 0``, keel at ``-R``."""
    builder = MeshBuilder()
    ring: list[int] = []
    for index in range(facets):
        angle = 2.0 * math.pi * index / facets
        ring.append(builder.add_vertex((0.0, radius * math.cos(angle), radius * math.sin(angle))))
    far: list[int] = []
    for index in range(facets):
        angle = 2.0 * math.pi * index / facets
        far.append(
            builder.add_vertex((length, radius * math.cos(angle), radius * math.sin(angle)))
        )
    for index in range(facets):
        following = (index + 1) % facets
        builder.add_quad(ring[index], ring[following], far[following], far[index])
    for index in range(1, facets - 1):
        builder.add_triangle(ring[0], ring[index + 1], ring[index])
        builder.add_triangle(far[0], far[index], far[index + 1])
    return meshlib.orient_outward(builder.build()).mesh


# --------------------------------------------------------------------------- #
# 1. the box barge against the wall-sided closed form
# --------------------------------------------------------------------------- #
def test_box_gz_reproduces_the_wall_sided_formula(box_mesh: Mesh) -> None:
    angles = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
    for angle in angles:
        measured = righting_arm(box_mesh, angle, kg=KG).gz
        reference = box_righting_arm(breadth=BREADTH, draft=DRAFT, kg=KG, heel_deg=angle)
        assert measured == pytest.approx(reference, rel=1e-12), angle


def test_box_gz_is_zero_at_zero_heel(box_mesh: Mesh) -> None:
    point = righting_arm(box_mesh, 0.0, kg=KG)
    assert point.gz == pytest.approx(0.0, abs=1e-15)
    assert point.kn == pytest.approx(0.0, abs=1e-15)


def test_box_holds_its_displacement_at_every_heel(box_mesh: Mesh) -> None:
    for angle in (0.0, 7.5, 15.0, 22.5, 30.0):
        point = righting_arm(box_mesh, angle, kg=KG)
        assert point.volume == pytest.approx(LENGTH * BREADTH * DRAFT, rel=1e-12)
        assert point.hydrostatics.volume_error <= 1e-13


def test_wall_sided_limit_is_the_deck_edge_immersion_angle() -> None:
    limit = wall_sided_limit_deg(breadth=BREADTH, freeboard=FREEBOARD)
    assert limit == pytest.approx(math.degrees(math.atan(2.0 * FREEBOARD / BREADTH)), abs=1e-12)
    assert limit == pytest.approx(30.963756532, abs=1e-9)


def test_beyond_the_wall_sided_limit_the_reference_over_predicts() -> None:
    """Negative control: the formula is an upper bound once the deck edge is wet.

    A shallow-freeboard box is heeled past its limit angle and the measured arm is
    compared with the formula the module is tested against. The formula must be the
    larger of the two -- if the integrator were returning the formula rather than the
    geometry, this test would fail.
    """
    shallow = box_hull_mesh(LENGTH, BREADTH, DRAFT, 0.1)
    limit = wall_sided_limit_deg(breadth=BREADTH, freeboard=0.1)
    assert limit == pytest.approx(5.710593137, abs=1e-9)
    angle = 25.0
    measured = righting_arm(shallow, angle, kg=KG).gz
    reference = box_righting_arm(breadth=BREADTH, draft=DRAFT, kg=KG, heel_deg=angle)
    assert 0.0 < measured < 0.5 * reference


# --------------------------------------------------------------------------- #
# 2. the area under the curve against its closed-form integral
# --------------------------------------------------------------------------- #
def test_dynamical_stability_matches_the_closed_form(box_mesh: Mesh) -> None:
    """Inside the wall-sided validity sector (deck edge clear, < 31 deg) the area matches."""
    curve = gz_curve(box_mesh, kg=KG, up_to_deg=30.0, step_deg=1.25)
    measured = curve.dynamical_stability(30.0)
    reference = box_dynamical_stability(breadth=BREADTH, draft=DRAFT, kg=KG, upto_deg=30.0)
    assert measured == pytest.approx(reference, rel=5e-4)


def test_dynamical_stability_converges_second_order(box_mesh: Mesh) -> None:
    reference = box_dynamical_stability(breadth=BREADTH, draft=DRAFT, kg=KG, upto_deg=30.0)
    coarse = gz_curve(box_mesh, kg=KG, up_to_deg=30.0, step_deg=10.0).dynamical_stability(30.0)
    fine = gz_curve(box_mesh, kg=KG, up_to_deg=30.0, step_deg=5.0).dynamical_stability(30.0)
    coarse_error = abs(coarse - reference)
    fine_error = abs(fine - reference)
    assert fine_error < coarse_error
    assert coarse_error / fine_error > 3.0


def test_dynamical_stability_handles_a_partial_interval(box_mesh: Mesh) -> None:
    curve = gz_curve(box_mesh, kg=KG, up_to_deg=30.0, step_deg=5.0)
    reference = box_dynamical_stability(breadth=BREADTH, draft=DRAFT, kg=KG, upto_deg=17.0)
    assert curve.dynamical_stability(17.0) == pytest.approx(reference, rel=5e-3)


def test_dynamical_stability_needs_two_points(box_mesh: Mesh) -> None:
    curve = gz_curve(box_mesh, angles_deg=[10.0], kg=KG)
    with pytest.raises(LargeAngleError):
        dynamical_stability(curve.points)


# --------------------------------------------------------------------------- #
# 3. the slope at zero heel against the upright-waterplane GM
# --------------------------------------------------------------------------- #
def test_slope_at_zero_heel_equals_the_metacentric_height(hull_mesh: Mesh, capped_hull: Mesh) -> None:
    """``dGZ/dphi(0) = GM``, cross-checked against ``stability.initial_stability``.

    ``GM`` is measured on the *capped* mesh, which is what that module expects (it
    integrates a surface already clipped at the waterline). The slope is measured on
    the freeboarded mesh, by this module. The two share no code, and the agreement is
    also the evidence that adding freeboard above the waterline leaves ``KM`` alone.
    """
    gm = stability_module.initial_stability(capped_hull, kg=HULL_KG).gm_transverse
    step = 1e-4
    plus = righting_arm(hull_mesh, step, kg=HULL_KG).gz
    minus = righting_arm(hull_mesh, -step, kg=HULL_KG).gz
    slope = (plus - minus) / (2.0 * math.radians(step))
    assert slope == pytest.approx(gm, rel=1e-6)


def test_gm_sign_controls_the_small_angle_arm(hull_mesh: Mesh, capped_hull: Mesh) -> None:
    gm = stability_module.initial_stability(capped_hull, kg=HULL_KG).gm_transverse
    assert gm > 0.0
    assert righting_arm(hull_mesh, 2.0, kg=HULL_KG).gz > 0.0


# --------------------------------------------------------------------------- #
# 4. rotational symmetry: the circular-section "log"
# --------------------------------------------------------------------------- #
def test_log_barge_righting_arm_is_exactly_gm_sin_phi() -> None:
    """A circular section cannot move its centre of buoyancy, so ``GZ = (R - KG) sin phi``.

    ``KM = R`` for a half-submerged circular section (``KB + BM = R(1 - 4/3pi) + R4/3pi``),
    so this is also a large-angle statement of the metacentric height. The prism is a
    polygonal approximation, so the residual is a *facet* error and must shrink when
    the facet count rises.
    """
    radius = 1.0
    kg = 0.3 * radius
    errors = {}
    for facets in (32, 128):
        prism = log_prism_mesh(radius=radius, length=4.0, facets=facets)
        worst = 0.0
        for angle in (10.0, 20.0, 30.0, 45.0):
            measured = righting_arm(prism, angle, kg=kg).gz
            reference = (radius - kg) * math.sin(math.radians(angle))
            worst = max(worst, abs(measured - reference) / abs(reference))
        errors[facets] = worst
    assert errors[32] < 5e-3
    assert errors[128] < 1e-3
    assert errors[32] / errors[128] > 4.0


# --------------------------------------------------------------------------- #
# 5. freeboard extension
# --------------------------------------------------------------------------- #
def test_extend_freeboard_reproduces_the_wall_sided_box() -> None:
    """Extruding a waterline-capped box by the freeboard must give the box with freeboard."""
    capped = box_hull_mesh(LENGTH, BREADTH, DRAFT)
    extended = extend_freeboard(capped, FREEBOARD)
    assert meshlib.total_signed_volume(extended) == pytest.approx(
        LENGTH * BREADTH * (DRAFT + FREEBOARD), rel=1e-12
    )
    for angle in (0.0, 10.0, 20.0, 30.0):
        measured = righting_arm(extended, angle, kg=KG).gz
        reference = box_righting_arm(breadth=BREADTH, draft=DRAFT, kg=KG, heel_deg=angle)
        assert measured == pytest.approx(reference, rel=1e-12)


def test_extend_freeboard_is_waterline_neutral(capped_hull: Mesh, hull_mesh: Mesh) -> None:
    reference = hydrostatics_module.compute_hydrostatics(capped_hull, 0.0)
    extended = heeled_hydrostatics(hull_mesh, 0.0, waterline_z=0.0)
    assert extended.volume == pytest.approx(reference.displacement_volume, rel=1e-12)
    assert extended.centre_of_buoyancy[0] == pytest.approx(
        reference.longitudinal_centre_buoyancy, rel=1e-12
    )
    assert extended.centre_of_buoyancy[2] == pytest.approx(
        reference.vertical_centre_buoyancy, rel=1e-12
    )


def test_extended_hull_stays_a_closed_oriented_surface(hull_mesh: Mesh, capped_hull: Mesh) -> None:
    assert meshlib.is_watertight(hull_mesh)
    assert meshlib.is_edge_manifold(hull_mesh)
    euler = len(hull_mesh.vertices) - len({edge for face in hull_mesh.faces
                                           for edge in _face_edges(face)}) + len(hull_mesh.faces)
    assert euler == 2
    assert meshlib.total_signed_volume(hull_mesh) > meshlib.total_signed_volume(capped_hull)
    assert meshlib.zero_area_faces(hull_mesh) == ()


def _face_edges(face: tuple[int, int, int]) -> tuple[tuple[int, int], tuple[int, int], tuple[int, int]]:
    a, b, c = face
    return ((min(a, b), max(a, b)), (min(b, c), max(b, c)), (min(c, a), max(c, a)))


# --------------------------------------------------------------------------- #
# 6. the GZ curve as a curve
# --------------------------------------------------------------------------- #
def test_curve_is_zero_at_zero_and_odd_in_heel(hull_mesh: Mesh) -> None:
    assert righting_arm(hull_mesh, 0.0, kg=HULL_KG).gz == pytest.approx(0.0, abs=1e-14)
    port = righting_arm(hull_mesh, 20.0, kg=HULL_KG).gz
    starboard = righting_arm(hull_mesh, -20.0, kg=HULL_KG).gz
    assert port > 0.0
    assert port == pytest.approx(-starboard, rel=1e-12)


def test_curve_holds_displacement_at_every_point(hull_gz: GzCurve) -> None:
    reference = hull_gz.points[0].volume
    for point in hull_gz.points:
        assert point.volume == pytest.approx(reference, rel=1e-12)
        assert point.hydrostatics.volume_error <= 1e-12


def test_heeled_waterline_sinks_as_the_section_is_lost(hull_gz: GzCurve) -> None:
    """Conserving displacement while losing waterplane area must lower the waterline."""
    waterlines = [point.hydrostatics.waterline_z for point in hull_gz.points]
    assert waterlines[0] == pytest.approx(0.0, abs=1e-12)
    assert all(later <= earlier for earlier, later in pairwise(waterlines))
    assert waterlines[-1] < -0.01


def test_righting_arm_rises_before_it_falls(hull_gz: GzCurve) -> None:
    arms = hull_gz.arms
    assert arms[1] > arms[0]
    assert arms[-1] > arms[1]
    assert hull_gz.peak.gz > 0.0
    assert hull_gz.peak.heel_deg >= 10.0


def test_negative_control_high_kg_makes_the_hull_unstable(hull_mesh: Mesh, capped_hull: Mesh) -> None:
    """Raise ``KG`` above ``KM``: the arm must go negative, not stay positive."""
    km = stability_module.initial_stability(capped_hull, kg=HULL_KG).km_transverse
    assert righting_arm(hull_mesh, 5.0, kg=km + 0.2).gz < 0.0
    assert righting_arm(hull_mesh, 5.0, kg=km - 0.2).gz > 0.0


def test_vanishing_angle_is_interpolated_inside_the_sector(capped_hull: Mesh) -> None:
    """A low-freeboard hull loses its waterplane early, so the arm rises then vanishes.

    ``freeboard = 0.1 m`` on a 2 m beam puts the deck edge in the water at
    ``atan(0.2/2) = 5.7 deg``: the curve peaks around 20 deg and crosses zero before
    60 deg, which is what ``vanishing_angle_deg`` has to bracket by interpolation.
    """
    low_freeboard = extend_freeboard(capped_hull, 0.1)
    curve = gz_curve(low_freeboard, kg=0.6, up_to_deg=60.0, step_deg=5.0)
    vanishing = curve.vanishing_angle_deg
    assert vanishing is not None
    assert 0.0 < vanishing < 60.0
    lower = next(point for point in reversed(curve.points) if point.heel_deg < vanishing)
    upper = next(point for point in curve.points if point.heel_deg > vanishing)
    assert lower.gz > 0.0 > upper.gz
    expected = lower.heel_deg + lower.gz / (lower.gz - upper.gz) * (
        upper.heel_deg - lower.heel_deg
    )
    assert vanishing == pytest.approx(expected, rel=1e-12)


def test_curve_reports_unbracketed_vanishing_angle(hull_gz: GzCurve) -> None:
    assert hull_gz.vanishing_angle_deg is None


def test_at_rejects_an_uncomputed_angle(hull_gz: GzCurve) -> None:
    with pytest.raises(LargeAngleError):
        hull_gz.at(11.5)


def test_curve_to_dict_is_serialisable(hull_gz: GzCurve) -> None:
    payload = hull_gz.to_dict()
    assert payload["schema"] == "tensorcad.stability.gz_curve/1"
    assert len(payload["points"]) == len(hull_gz)
    assert payload["kg_m"] == pytest.approx(HULL_KG)
    assert payload["peak_gz_m"] == pytest.approx(hull_gz.peak.gz, abs=1e-9)
    for key in ("heel_deg", "gz_m", "kn_m", "volume_m3", "waterline_z_m"):
        assert key in payload["points"][0]
    assert isinstance(payload["dynamical_stability_m_rad_upto_40deg"], float)


def test_gz_equals_kn_minus_kg_sin_phi(hull_gz: GzCurve) -> None:
    for point in hull_gz.points:
        assert point.gz == pytest.approx(
            point.kn - HULL_KG * math.sin(math.radians(point.heel_deg)), abs=1e-14
        )


# --------------------------------------------------------------------------- #
# 7. held-waterline convention
# --------------------------------------------------------------------------- #
def test_held_waterline_convention_immerses_freeboard_when_heeled(hull_mesh: Mesh) -> None:
    """With the waterplane held, heeling dips the freeboard in: displacement must rise."""
    upright = heeled_hydrostatics(hull_mesh, 0.0, waterline_z=0.0)
    heeled = heeled_hydrostatics(hull_mesh, 20.0, waterline_z=0.0)
    assert heeled.waterline_z == 0.0
    assert heeled.volume > upright.volume
    assert heeled.target_volume == heeled.volume


# --------------------------------------------------------------------------- #
# 8. failure modes that would otherwise return a plausible number
# --------------------------------------------------------------------------- #
def test_mesh_capped_at_the_waterline_is_refused(capped_hull: Mesh) -> None:
    with pytest.raises(LargeAngleError, match="freeboard"):
        righting_arm(capped_hull, 10.0, kg=HULL_KG)


def test_open_mesh_is_refused(hull_mesh: Mesh) -> None:
    open_mesh = Mesh(vertices=hull_mesh.vertices, faces=hull_mesh.faces[:-1])
    with pytest.raises(LargeAngleError, match="closed"):
        righting_arm(open_mesh, 10.0, kg=HULL_KG)


def test_inward_oriented_mesh_is_refused(hull_mesh: Mesh) -> None:
    flipped = Mesh(
        vertices=hull_mesh.vertices,
        faces=tuple((a, c, b) for (a, b, c) in hull_mesh.faces),  # type: ignore[arg-type]
    )
    with pytest.raises(LargeAngleError):
        righting_arm(flipped, 10.0, kg=HULL_KG)


def test_displacement_and_waterline_are_mutually_exclusive(hull_mesh: Mesh) -> None:
    with pytest.raises(LargeAngleError, match="not both"):
        heeled_hydrostatics(hull_mesh, 10.0, displacement_volume=1.0, waterline_z=0.0)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), 200.0, -181.0])
def test_invalid_heel_is_refused(hull_mesh: Mesh, bad: float) -> None:
    with pytest.raises(LargeAngleError):
        heeled_hydrostatics(hull_mesh, bad)


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_invalid_kg_is_refused(hull_mesh: Mesh, bad: float) -> None:
    with pytest.raises(LargeAngleError):
        righting_arm(hull_mesh, 10.0, kg=bad)


def test_kg_is_a_required_keyword(hull_mesh: Mesh) -> None:
    with pytest.raises(TypeError):
        righting_arm(hull_mesh, 10.0)  # type: ignore[call-arg]


def test_flat_mesh_is_refused() -> None:
    flat = Mesh(
        vertices=((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        faces=((0, 1, 2),),
    )
    with pytest.raises(LargeAngleError):
        heeled_hydrostatics(flat, 0.0)


@pytest.mark.parametrize("bad_angle", [0.0, -5.0])
def test_step_and_limit_must_be_positive(box_mesh: Mesh, bad_angle: float) -> None:
    with pytest.raises(LargeAngleError):
        gz_curve(box_mesh, kg=KG, step_deg=bad_angle)
    with pytest.raises(LargeAngleError):
        gz_curve(box_mesh, kg=KG, up_to_deg=bad_angle)


def test_empty_angle_list_is_refused(box_mesh: Mesh) -> None:
    with pytest.raises(LargeAngleError):
        gz_curve(box_mesh, [], kg=KG)


def test_deck_detection_failure_is_reported() -> None:
    """A mesh with no flat cap cannot have its freeboard extended, and says so."""
    prism = log_prism_mesh(radius=1.0, length=4.0, facets=12)
    with pytest.raises(LargeAngleError):
        extend_freeboard(prism, 0.5)


# --------------------------------------------------------------------------- #
# 9. the public surface
# --------------------------------------------------------------------------- #
def test_module_is_exported_at_package_level() -> None:
    for name in (
        "LargeAngleError",
        "HeeledHydrostatics",
        "RightingArm",
        "GzCurve",
        "box_hull_mesh",
        "box_righting_arm",
        "box_dynamical_stability",
        "dynamical_stability",
        "extend_freeboard",
        "gz_curve",
        "heeled_hydrostatics",
        "righting_arm",
        "wall_sided_limit_deg",
    ):
        assert getattr(tensorcad, name) is tensorcad.large_angle.__dict__[name]


def test_heeled_hydrostatics_to_dict_rounds(hull_mesh: Mesh) -> None:
    payload: HeeledHydrostatics = heeled_hydrostatics(hull_mesh, 15.0)
    data = payload.to_dict()
    assert data["heel_deg"] == pytest.approx(15.0)
    assert data["volume_m3"] == pytest.approx(payload.volume, abs=1e-9)
    assert len(data["centre_of_buoyancy_m"]) == 3
    assert isinstance(RightingArm.__name__, str)
