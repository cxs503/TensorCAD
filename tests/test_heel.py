"""Righting arms at finite heel: the heeled integration against closed forms.

Three references, each one attacking a different way the integration could be wrong:

* a **rectangular barge**, where the wall-sided formula
  ``GZ = sin(phi) (GM + BM tan^2(phi) / 2)`` is exact until the deck edge immerses: the
  mesh integral must reproduce it to round-off, and must *leave* it, on the right side,
  once the deck edge goes under (the formula assumes unbounded walls);
* an **elliptical prism**, same wall-sided reference but with a waterplane whose second
  moment is not the box's ``L B^3 / 12`` -- so ``BM`` has to come out of the actual
  section rather than out of a rectangle shortcut;
* a **sphere**, where rotating the body changes nothing: the centre of buoyancy stays on
  the vertical through the centre, so ``GZ = sin(phi) (z_C - z_G)`` exactly, at *any*
  heel, with no deck to immerse. The residual here is the facet error of the mesh and is
  shown to fall with refinement instead of being asserted small.

Plus the internal consistency checks that no shared bug can fake: ``GZ(0) = 0``,
``GZ(-phi) = -GZ(phi)``, constant displacement at every sampled angle, and
``dGZ/dphi(0) = GM`` where ``GM`` came from the waterplane moments rather than from this
integration.
"""

from __future__ import annotations

import math

import pytest

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.heel import (
    DEFAULT_ANGLES_DEG,
    GZCurve,
    HeeledSection,
    HeelError,
    SubmergedBody,
    gz_curve,
    heeled_displacement,
    righting_arm,
    wall_sided_gz,
    waterplane_normal,
)
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.stability import initial_stability, waterplane_moments


# ---------------------------------------------------------------------------
# reference bodies
# ---------------------------------------------------------------------------
def box(length: float, beam: float, draft: float, depth: float) -> Mesh:
    """Closed box, keel on ``z = -draft``, deck on ``z = depth - draft``, outward-faced."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, depth - draft):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(*points: tuple[float, float, float]) -> None:
        a, b, c, d = (corner[point] for point in points)
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi = -0.5 * beam, 0.5 * beam
    keel, deck = -draft, depth - draft
    quad((0.0, lo, keel), (0.0, hi, keel), (length, hi, keel), (length, lo, keel))
    quad((0.0, lo, deck), (length, lo, deck), (length, hi, deck), (0.0, hi, deck))
    quad((0.0, lo, keel), (length, lo, keel), (length, lo, deck), (0.0, lo, deck))
    quad((0.0, hi, keel), (0.0, hi, deck), (length, hi, deck), (length, hi, keel))
    quad((0.0, lo, keel), (0.0, lo, deck), (0.0, hi, deck), (0.0, hi, keel))
    quad((length, lo, keel), (length, hi, keel), (length, hi, deck), (length, lo, deck))
    return builder.build()


def elliptical_prism(
    length: float, beam: float, draft: float, depth: float, *, facets: int = 64
) -> Mesh:
    """Vertical-sided prism on an elliptical waterplane (wall-sided, but not a box)."""
    builder = MeshBuilder()
    a, b = 0.5 * length, 0.5 * beam
    keel, deck = -draft, depth - draft
    bottom = []
    top = []
    for index in range(facets):
        theta = 2.0 * math.pi * index / facets
        x = a * math.cos(theta)
        y = b * math.sin(theta)
        bottom.append(builder.add_vertex((x, y, keel)))
        top.append(builder.add_vertex((x, y, deck)))
    centre_bottom = builder.add_vertex((0.0, 0.0, keel))
    centre_top = builder.add_vertex((0.0, 0.0, deck))
    for index in range(facets):
        nxt = (index + 1) % facets
        builder.add_triangle(bottom[index], bottom[nxt], top[nxt])
        builder.add_triangle(bottom[index], top[nxt], top[index])
        builder.add_triangle(centre_bottom, bottom[nxt], bottom[index])
        builder.add_triangle(centre_top, top[index], top[nxt])
    return builder.build()


def sphere(radius: float, centre_z: float, *, rings: int = 24, segments: int = 48) -> Mesh:
    """UV sphere, outward-faced, centred on the z axis at ``centre_z``."""
    builder = MeshBuilder()
    north = builder.add_vertex((0.0, 0.0, centre_z + radius))
    bands = []
    for ring in range(1, rings):
        theta = math.pi * ring / rings
        band = []
        for segment in range(segments):
            phi = 2.0 * math.pi * segment / segments
            band.append(
                builder.add_vertex(
                    (
                        radius * math.sin(theta) * math.cos(phi),
                        radius * math.sin(theta) * math.sin(phi),
                        centre_z + radius * math.cos(theta),
                    )
                )
            )
        bands.append(band)
    south = builder.add_vertex((0.0, 0.0, centre_z - radius))
    for segment in range(segments):
        nxt = (segment + 1) % segments
        builder.add_triangle(north, bands[0][segment], bands[0][nxt])
    for ring in range(len(bands) - 1):
        lower, upper = bands[ring], bands[ring + 1]
        for segment in range(segments):
            nxt = (segment + 1) % segments
            builder.add_triangle(lower[segment], upper[nxt], lower[nxt])
            builder.add_triangle(lower[segment], upper[segment], upper[nxt])
    for segment in range(segments):
        nxt = (segment + 1) % segments
        builder.add_triangle(south, bands[-1][nxt], bands[-1][segment])
    return builder.build()


#: Barge used by the wall-sided checks: L=20, B=6, T=2, freeboard 2.
BARGE = (20.0, 6.0, 2.0, 4.0)
BARGE_KB = 1.0
BARGE_BM = 6.0**2 / (12.0 * 2.0)
BARGE_DECK_EDGE_DEG = math.degrees(math.atan(2.0 * 2.0 / 6.0))


# ---------------------------------------------------------------------------
# waterplane normal and the wall-sided reference itself
# ---------------------------------------------------------------------------
def test_waterplane_normal_tilts_towards_the_low_side() -> None:
    up = waterplane_normal(30.0)
    assert up[0] == 0.0
    assert up[1] == pytest.approx(-0.5)
    assert up[2] == pytest.approx(math.cos(math.radians(30.0)))
    # starboard (y > 0) is under water: its height along "up" is negative
    assert up[1] * 1.0 + up[2] * 0.0 < 0.0
    assert waterplane_normal(0.0) == (0.0, -0.0, 1.0)


def test_waterplane_normal_rejects_non_finite() -> None:
    with pytest.raises(HeelError):
        waterplane_normal(float("nan"))


def test_wall_sided_reference_reduces_to_gm_sin_phi_at_small_heel() -> None:
    gm = BARGE_KB + BARGE_BM - 2.0
    # the tan^2 term is O(1e-4) of GM at one degree, so this is a 1e-3-level identity
    assert wall_sided_gz(kb=BARGE_KB, bm=BARGE_BM, kg=2.0, heel_deg=1.0) == pytest.approx(
        math.radians(1.0) * gm, rel=1e-3
    )
    assert wall_sided_gz(kb=BARGE_KB, bm=BARGE_BM, kg=2.0, heel_deg=0.0) == 0.0


def test_wall_sided_reference_rejects_bad_input() -> None:
    with pytest.raises(HeelError):
        wall_sided_gz(kb=1.0, bm=-0.1, kg=1.0, heel_deg=10.0)
    with pytest.raises(HeelError):
        wall_sided_gz(kb=1.0, bm=1.0, kg=float("inf"), heel_deg=10.0)


# ---------------------------------------------------------------------------
# the box: exact to round-off, and honest once the deck edge immerses
# ---------------------------------------------------------------------------
def test_box_reproduces_wall_sided_formula_below_deck_edge() -> None:
    mesh = box(*BARGE)
    kg = 2.0
    angles = (0.0, 5.0, 10.0, 20.0, 30.0)
    curve = gz_curve(mesh, kg=kg, angles_deg=angles)
    assert curve.displacement_volume == pytest.approx(20.0 * 6.0 * 2.0)
    for section in curve.sections:
        reference = wall_sided_gz(kb=BARGE_KB, bm=BARGE_BM, kg=kg, heel_deg=section.heel_deg)
        assert section.gz == pytest.approx(reference, abs=1e-12, rel=1e-11), section.heel_deg


def test_box_heeled_waterline_stays_on_the_upright_waterline() -> None:
    """A wall-sided body keeps its volume with the water surface on ``z = 0``."""
    mesh = box(*BARGE)
    curve = gz_curve(mesh, kg=2.0, angles_deg=(10.0, 20.0, 30.0))
    for section in curve.sections:
        assert section.waterline_offset == pytest.approx(0.0, abs=1e-12)


def test_box_departs_from_wall_sided_after_deck_immersion() -> None:
    mesh = box(*BARGE)
    kg = 2.0
    curve = gz_curve(mesh, kg=kg, angles_deg=(30.0, 40.0, 60.0))
    reference_40 = wall_sided_gz(kb=BARGE_KB, bm=BARGE_BM, kg=kg, heel_deg=40.0)
    gz_40 = curve.at(40.0)
    assert 40.0 > BARGE_DECK_EDGE_DEG
    # the unbounded-wall formula over-predicts once the deck edge is under water
    assert gz_40 < reference_40
    assert (reference_40 - gz_40) / reference_40 > 0.05


def test_box_has_a_finite_range_of_stability() -> None:
    """A shallow-freeboard barge capsizes; the crossing is interpolated, not guessed."""
    mesh = box(20.0, 6.0, 2.0, 3.0)
    curve = gz_curve(mesh, kg=2.2, angles_deg=tuple(float(a) for a in range(0, 95, 5)))
    assert curve.gm_transverse == pytest.approx(BARGE_KB + BARGE_BM - 2.2, rel=1e-12)
    angle, gz_max = curve.max_gz
    assert gz_max > 0.1 and 15.0 <= angle <= 35.0
    vanishing = curve.vanishing_angle_deg
    assert vanishing is not None and 40.0 < vanishing < 50.0
    assert curve.at(vanishing) == pytest.approx(0.0, abs=1e-12)
    assert curve.gz[curve.angles_deg.index(angle)] > 0.0


def test_box_at_ninety_degrees_has_zero_righting_arm() -> None:
    """At 90 deg the waterline is a symmetry plane of the barge, so GZ is exactly zero."""
    mesh = box(*BARGE)
    curve = gz_curve(mesh, kg=2.0, angles_deg=(90.0,))
    assert curve.sections[0].gz == pytest.approx(0.0, abs=1e-12)
    assert curve.sections[0].centre_of_buoyancy[2] == pytest.approx(0.0, abs=1e-12)


def test_negative_metacentric_height_gives_negative_arm_at_small_heel() -> None:
    mesh = box(*BARGE)
    curve = gz_curve(mesh, kg=2.6, angles_deg=(0.0, 5.0, 10.0))
    assert curve.gm_transverse == pytest.approx(BARGE_KB + BARGE_BM - 2.6, rel=1e-12)
    assert curve.gm_transverse < 0.0
    assert curve.sections[0].gz == pytest.approx(0.0, abs=1e-15)
    assert curve.sections[1].gz < 0.0
    assert curve.max_gz == (0.0, pytest.approx(0.0, abs=1e-15))


# ---------------------------------------------------------------------------
# the elliptical prism: the same reference with a different waterplane
# ---------------------------------------------------------------------------
def test_elliptical_prism_reproduces_wall_sided_with_its_own_section() -> None:
    length, beam, draft, depth = 12.0, 4.0, 1.5, 3.0
    facets = 72
    mesh = elliptical_prism(length, beam, draft, depth, facets=facets)
    section = waterplane_moments(mesh, 0.0)
    volume = 0.5 * facets * math.sin(2.0 * math.pi / facets) * 0.5 * length * 0.5 * beam * draft
    kb = 0.5 * draft
    bm = section.transverse_second_moment / volume
    kg = 1.2
    curve = gz_curve(mesh, kg=kg, angles_deg=(0.0, 10.0, 20.0, 30.0))
    assert curve.displacement_volume == pytest.approx(volume, rel=1e-12)
    assert curve.gm_transverse == pytest.approx(kb + bm - kg, rel=1e-9)
    for point in curve.sections:
        reference = wall_sided_gz(kb=kb, bm=bm, kg=kg, heel_deg=point.heel_deg)
        assert point.gz == pytest.approx(reference, abs=1e-11, rel=1e-10), point.heel_deg


def test_elliptical_prism_matches_the_analytic_ellipse_inertia() -> None:
    length, beam, facets = 12.0, 4.0, 512
    mesh = elliptical_prism(length, beam, 1.5, 3.0, facets=facets)
    section = waterplane_moments(mesh, 0.0)
    analytic = math.pi * (0.5 * length) * (0.5 * beam) ** 3 / 4.0
    assert section.transverse_second_moment == pytest.approx(analytic, rel=2e-3)


# ---------------------------------------------------------------------------
# the sphere: an exact large-angle reference with no deck to immerse
# ---------------------------------------------------------------------------
def sphere_case(radius: float = 1.0, draft_fraction: float = 0.6, **kwargs) -> tuple[Mesh, float, float]:
    centre_z = radius - draft_fraction * radius
    return sphere(radius, centre_z, **kwargs), centre_z, centre_z - radius


def test_sphere_gz_matches_the_rotational_symmetry_closed_form() -> None:
    """``GZ = sin(phi) (z_C - z_G)`` -- the body is a sphere, so heel changes nothing."""
    kg = 0.3
    mesh, centre_z, keel_z = sphere_case()
    curve = gz_curve(mesh, kg=kg, angles_deg=(0.0, 10.0, 20.0, 40.0, 60.0, 80.0))
    assert curve.gm_transverse == pytest.approx(centre_z - (keel_z + kg), rel=5e-3)
    for section in curve.sections:
        reference = math.sin(math.radians(section.heel_deg)) * (centre_z - (keel_z + kg))
        assert section.gz == pytest.approx(reference, abs=2e-3), section.heel_deg


def test_sphere_residual_is_the_facet_error_and_falls_with_refinement() -> None:
    kg = 0.3
    heel = 60.0
    errors = []
    for rings in (12, 24, 48):
        mesh, centre_z, keel_z = sphere_case(rings=rings, segments=2 * rings)
        curve = gz_curve(mesh, kg=kg, angles_deg=(heel,))
        reference = math.sin(math.radians(heel)) * (centre_z - (keel_z + kg))
        errors.append(abs(curve.sections[0].gz - reference))
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < 5e-4
    assert errors[0] / errors[2] > 4.0


def test_sphere_waterline_offset_moves_with_heel() -> None:
    """Constant displacement does not mean constant waterline height: the offset is solved."""
    mesh, _, _ = sphere_case()
    curve = gz_curve(mesh, kg=0.3, angles_deg=(0.0, 80.0))
    upright, heeled = curve.sections
    assert heeled.waterline_offset < upright.waterline_offset - 0.2


# ---------------------------------------------------------------------------
# internal consistency: no shared bug can fake these
# ---------------------------------------------------------------------------
def test_odd_symmetry_of_the_righting_arm() -> None:
    mesh = box(*BARGE)
    curve = gz_curve(mesh, kg=2.0, angles_deg=(10.0, 20.0, 30.0, 40.0))
    mirrored = gz_curve(mesh, kg=2.0, angles_deg=(-40.0, -30.0, -20.0, -10.0))
    for positive, negative in zip(curve.sections, reversed(mirrored.sections)):
        assert negative.gz == pytest.approx(-positive.gz, abs=1e-12)


def test_zero_heel_gives_zero_arm_and_the_upright_buoyancy_centre() -> None:
    mesh = box(*BARGE)
    curve = gz_curve(mesh, kg=2.0, angles_deg=(0.0,))
    section = curve.sections[0]
    assert section.gz == pytest.approx(0.0, abs=1e-15)
    assert section.centre_of_buoyancy[1] == pytest.approx(0.0, abs=1e-12)
    assert section.centre_of_buoyancy[2] == pytest.approx(-1.0, abs=1e-12)
    assert section.waterline_offset == pytest.approx(0.0, abs=1e-12)


def test_small_angle_slope_equals_the_metacentric_height() -> None:
    mesh = box(*BARGE)
    kg = 2.0
    step = 0.5
    curve = gz_curve(mesh, kg=kg, angles_deg=(0.0, step, 2.0 * step))
    slope = (curve.at(2.0 * step) - curve.at(0.0)) / math.radians(2.0 * step)
    # a 2 deg chord carries the wall-sided tan^2 term: ~1.4e-3 high, not 0
    assert slope == pytest.approx(curve.gm_transverse, rel=2e-3)
    assert curve.gm_transverse == pytest.approx(BARGE_KB + BARGE_BM - kg, rel=1e-12)


def test_constant_displacement_is_held_at_every_angle() -> None:
    mesh = box(*BARGE)
    target = 0.5 * 20.0 * 6.0 * 2.0
    curve = gz_curve(mesh, kg=2.0, displacement_volume=target, angles_deg=range(0, 80, 10))
    for section in curve.sections:
        assert section.displacement_volume == pytest.approx(target, rel=1e-9)


def test_this_module_computes_its_own_gm_for_bodies_with_freeboard() -> None:
    """``initial_stability`` divides by the whole closed volume; here it must not.

    The barge displaces 240 m^3 but encloses 480 m^3, so a ``BM`` that used the closed
    volume would be half the truth and would report ``GM = -0.25`` m instead of ``0.5`` m.
    """
    mesh = box(*BARGE)
    st = initial_stability(mesh, kg=2.0)
    assert st.displacement_volume == pytest.approx(480.0)
    # KB from the full-body centroid (2.0 instead of 1.0) and BM from the closed volume
    # (0.75 instead of 1.5): two errors that only cancel by accident here.
    assert st.gm_transverse == pytest.approx(2.0 + 0.75 - 2.0, rel=1e-12)
    assert abs(st.gm_transverse - 0.5) > 0.2
    curve = gz_curve(mesh, kg=2.0, angles_deg=(0.0,))
    assert curve.displacement_volume == pytest.approx(240.0)
    assert curve.gm_transverse == pytest.approx(0.5, rel=1e-12)


def test_default_angles_cover_upright_to_sixty_degrees() -> None:
    assert DEFAULT_ANGLES_DEG[0] == 0.0
    assert DEFAULT_ANGLES_DEG[-1] == 60.0
    assert len(DEFAULT_ANGLES_DEG) == 13
    curve = gz_curve(box(*BARGE), kg=2.0)
    assert curve.angles_deg == DEFAULT_ANGLES_DEG


# ---------------------------------------------------------------------------
# the public API surface
# ---------------------------------------------------------------------------
def test_heeled_displacement_reports_the_immersed_body() -> None:
    mesh = box(*BARGE)
    body = heeled_displacement(mesh, heel_deg=30.0, displacement_volume=240.0)
    assert isinstance(body, SubmergedBody)
    assert body.displacement_volume == pytest.approx(240.0, rel=1e-9)
    assert abs(body.residual_volume) < 1e-9
    assert body.centre_of_buoyancy[1] > 0.0  # starboard-down heel moves B to starboard
    data = body.to_dict()
    assert data["schema"] == "tensorcad.heel.submerged/1"
    assert data["heel_deg"] == 30.0


def test_righting_arm_is_callable_on_a_solved_body() -> None:
    mesh = box(*BARGE)
    body = heeled_displacement(mesh, heel_deg=20.0, displacement_volume=240.0)
    assert righting_arm(body, kg=2.0, keel_z=-2.0) == pytest.approx(
        wall_sided_gz(kb=BARGE_KB, bm=BARGE_BM, kg=2.0, heel_deg=20.0), abs=1e-12
    )
    with pytest.raises(HeelError):
        righting_arm(body, kg=float("nan"), keel_z=-2.0)


def test_curve_dict_is_serialisable_and_self_describing() -> None:
    curve = gz_curve(box(20.0, 6.0, 2.0, 3.0), kg=2.2, angles_deg=range(0, 95, 5))
    assert isinstance(curve, GZCurve)
    assert all(isinstance(section, HeeledSection) for section in curve.sections)
    data = curve.to_dict()
    assert data["schema"] == "tensorcad.heel.gz_curve/1"
    assert data["kg_m"] == 2.2
    assert data["displacement_volume_m3"] == pytest.approx(240.0)
    assert data["angle_of_max_gz_deg"] == curve.max_gz[0]
    assert data["vanishing_angle_deg"] == pytest.approx(curve.vanishing_angle_deg)
    assert len(data["sections"]) == len(curve.sections)
    assert "constant displacement" in data["basis"]


def test_truncated_curve_does_not_invent_a_vanishing_angle() -> None:
    curve = gz_curve(box(*BARGE), kg=2.0, angles_deg=(0.0, 10.0, 20.0))
    assert curve.vanishing_angle_deg is None
    assert curve.to_dict()["vanishing_angle_deg"] is None


def test_at_rejects_angles_outside_the_sampled_range() -> None:
    curve = gz_curve(box(*BARGE), kg=2.0, angles_deg=(0.0, 10.0))
    with pytest.raises(HeelError):
        curve.at(11.0)


# ---------------------------------------------------------------------------
# refusals: the interesting failure is a body with no freeboard
# ---------------------------------------------------------------------------
def test_hull_capped_at_the_waterline_has_no_heeled_equilibrium() -> None:
    mesh = generate_hull(HullParameters(nx=21, nz=9)).mesh
    with pytest.raises(HeelError, match="freeboard"):
        gz_curve(mesh, kg=0.4, angles_deg=(10.0,))


def test_displacement_above_the_body_volume_is_refused() -> None:
    mesh = box(*BARGE)
    with pytest.raises(HeelError, match="whole volume"):
        heeled_displacement(mesh, heel_deg=10.0, displacement_volume=1000.0)


def test_bad_inputs_are_refused() -> None:
    mesh = box(*BARGE)
    with pytest.raises(HeelError):
        heeled_displacement(mesh, heel_deg=10.0, displacement_volume=0.0)
    with pytest.raises(HeelError):
        heeled_displacement(mesh, heel_deg=10.0, displacement_volume=240.0, tolerance=-1.0)
    with pytest.raises(HeelError):
        gz_curve(mesh, kg=float("inf"), angles_deg=(0.0,))
    with pytest.raises(HeelError):
        gz_curve(mesh, kg=2.0, angles_deg=())
    with pytest.raises(HeelError):
        gz_curve(mesh, kg=2.0, angles_deg=(0.0,), displacement_volume=-5.0)


def test_inward_faced_mesh_is_refused() -> None:
    mesh = box(*BARGE)
    flipped = Mesh(
        vertices=mesh.vertices,
        faces=tuple((a, c, b) for a, b, c in mesh.faces),
        units=mesh.units,
    )
    with pytest.raises(HeelError, match="outward"):
        heeled_displacement(flipped, heel_deg=10.0, displacement_volume=240.0)