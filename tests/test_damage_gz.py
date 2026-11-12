"""The damaged ship's righting arm: held to the box closed forms and to an exact oracle.

Three independent standards, in increasing strength:

1. **The wall-sided box formula** ``GZ = sin(phi) (GM + BM/2 tan^2 phi)``, which is
   exact for a rectangular hull only while the waterline stays on both vertical sides.
   The tests assert the agreement *inside* that range and assert the disagreement
   outside it, so the comparison cannot pass by being vacuous.
2. **An exact closed-form oracle** (``tests/damage_gz_oracle.py``): the plane is
   independent of ``x`` for a midships compartment, so the net buoyancy of a box is a
   one-dimensional integral of piecewise-linear functions and every moment is
   elementary. No grid, no refinement, no tolerance to negotiate.
Plus the invariants that hold whatever the mesh: displaced volume and longitudinal
centre stay on their targets at every heel, the arm is odd in the heel for a centred
compartment, a zero-permeability compartment is the intact ship, and the whole curve
slides with the hull under a translation.
"""
from __future__ import annotations

import math

import pytest

from damage_gz_oracle import exact_damaged_arm, wall_sided_arm
from tensorcad import HullParameters, generate_hull
from tensorcad.damage import (
    MAX_ARM_HEEL_DEG,
    Compartment,
    DamageError,
    damage_gm,
    flooded_equilibrium,
    flooded_righting_arm,
    flooded_righting_arm_curve,
)
from tensorcad.geometry import Mesh
from tensorcad.gz import with_freeboard

L, B, T = 10.0, 2.0, 0.8
DISPLACEMENT = L * B * T
FREEBOARD = 1.0
DECK_Z = FREEBOARD
KEEL_Z = -T
KG = 0.5
BM0 = B * B / (12.0 * T)
GM0 = 0.5 * T + BM0 - KG
#: ``|y| = B/2`` reaches the keel at this heel: past it the wall-sided form is gone.
WALL_SIDED_LIMIT_DEG = math.degrees(math.atan(2.0 * T / B))


def box_mesh(
    length: float = L,
    beam: float = B,
    draft: float = T,
    freeboard: float = FREEBOARD,
    x_shift: float = 0.0,
) -> Mesh:
    """A watertight rectangular hull, keel at ``-draft``, deck at ``+freeboard``."""
    xs = [x_shift, x_shift + length]
    ys = [-0.5 * beam, 0.5 * beam]
    zs = [-draft, freeboard]
    vertices: list[tuple[float, float, float]] = []
    index: dict[tuple[float, float, float], int] = {}
    faces: list[tuple[int, int, int]] = []

    def vid(point):
        if point not in index:
            index[point] = len(vertices)
            vertices.append(point)
        return index[point]

    def quad(a, b, c, d):
        faces.append((vid(a), vid(b), vid(c)))
        faces.append((vid(a), vid(c), vid(d)))

    quad((xs[0], ys[0], zs[0]), (xs[0], ys[1], zs[0]), (xs[1], ys[1], zs[0]), (xs[1], ys[0], zs[0]))
    quad((xs[0], ys[0], zs[1]), (xs[1], ys[0], zs[1]), (xs[1], ys[1], zs[1]), (xs[0], ys[1], zs[1]))
    quad((xs[0], ys[0], zs[0]), (xs[1], ys[0], zs[0]), (xs[1], ys[0], zs[1]), (xs[0], ys[0], zs[1]))
    quad((xs[1], ys[1], zs[0]), (xs[0], ys[1], zs[0]), (xs[0], ys[1], zs[1]), (xs[1], ys[1], zs[1]))
    quad((xs[0], ys[1], zs[0]), (xs[0], ys[0], zs[0]), (xs[0], ys[0], zs[1]), (xs[0], ys[1], zs[1]))
    quad((xs[1], ys[0], zs[0]), (xs[1], ys[1], zs[0]), (xs[1], ys[1], zs[1]), (xs[1], ys[0], zs[1]))
    return Mesh(vertices=vertices, faces=faces)


@pytest.fixture(scope="module")
def box() -> Mesh:
    return box_mesh()


@pytest.fixture(scope="module")
def midship() -> Compartment:
    """Full-breadth, bottom-to-waterline, 60 % permeable: the classic midship case."""
    return Compartment(
        4.0, 6.0, y0=-1.0, y1=1.0, z0=-0.8, z1=0.0, permeability=0.6, name="midship"
    )


@pytest.fixture(scope="module")
def damaged_curve(box, midship):
    return flooded_righting_arm_curve(
        box, midship, DISPLACEMENT, kg=KG, heels_deg=tuple(range(0, 61, 5))
    )


# --------------------------------------------------------------------------- #
# 1. the intact ship: the closed form, and the exact oracle
# --------------------------------------------------------------------------- #
def test_the_intact_box_reproduces_the_wall_sided_form_where_it_holds(box):
    curve = flooded_righting_arm_curve(
        box, (), DISPLACEMENT, kg=KG, heels_deg=tuple(range(0, 36, 5))
    )
    for point in curve.points:
        reference = wall_sided_arm(point.heel_deg, gm=GM0, bm=BM0)
        assert point.righting_arm == pytest.approx(reference, abs=1e-15)
        assert point.kb == pytest.approx(
            0.5 * T + 0.5 * BM0 * math.tan(math.radians(point.heel_deg)) ** 2, rel=1e-14
        )


def test_the_wall_sided_form_is_not_used_past_its_limit(box):
    """The comparison above must not be vacuous: past keel emergence it is wrong."""
    curve = flooded_righting_arm_curve(box, (), DISPLACEMENT, kg=KG, heels_deg=(60.0,))
    reference = wall_sided_arm(60.0, gm=GM0, bm=BM0)
    assert abs(curve.points[0].righting_arm - reference) > 0.3
    assert WALL_SIDED_LIMIT_DEG == pytest.approx(38.65980825409009, rel=1e-12)


def test_the_intact_box_matches_the_exact_oracle_at_every_heel(box):
    curve = flooded_righting_arm_curve(
        box, (), DISPLACEMENT, kg=KG, heels_deg=tuple(range(0, 61, 5))
    )
    for point in curve.points:
        exact = exact_damaged_arm(
            heel_deg=point.heel_deg,
            displacement_volume=DISPLACEMENT,
            length=L,
            beam=B,
            keel_z=KEEL_Z,
            deck_z=DECK_Z,
            kg=KG,
        )
        assert point.righting_arm == pytest.approx(exact.righting_arm, abs=1e-12)
        assert point.waterline.base_z == pytest.approx(exact.base_z, abs=1e-12)
        assert point.tcb == pytest.approx(exact.tcb, abs=1e-12)
        assert point.kb == pytest.approx(exact.kb, abs=1e-12)


# --------------------------------------------------------------------------- #
# 2. the damaged ship: the exact oracle, because no closed form covers it
# --------------------------------------------------------------------------- #
def test_the_damaged_box_matches_the_exact_oracle_at_every_heel(box, midship, damaged_curve):
    for point in damaged_curve.points:
        exact = exact_damaged_arm(
            heel_deg=point.heel_deg,
            displacement_volume=DISPLACEMENT,
            length=L,
            beam=B,
            keel_z=KEEL_Z,
            deck_z=DECK_Z,
            kg=KG,
            compartment=(4.0, 6.0, -0.8, 0.0, 0.6),
        )
        assert point.righting_arm == pytest.approx(exact.righting_arm, abs=1e-11)
        assert point.waterline.base_z == pytest.approx(exact.base_z, abs=1e-11)
        assert point.lcb == pytest.approx(exact.lcb, abs=1e-11)
        assert point.tcb == pytest.approx(exact.tcb, abs=1e-11)
        assert point.kb == pytest.approx(exact.kb, abs=1e-11)
        assert point.volume == pytest.approx(DISPLACEMENT, rel=1e-12)


def test_the_damaged_arm_is_not_the_intact_arm(box, midship, damaged_curve):
    """Flooding a bottom compartment raises KB, so the arm changes -- visibly."""
    intact = flooded_righting_arm_curve(box, (), DISPLACEMENT, kg=KG, heels_deg=(20.0,))
    assert damaged_curve.points[4].heel_deg == 20.0
    assert abs(damaged_curve.points[4].righting_arm - intact.points[0].righting_arm) > 5e-3


def test_the_small_angle_slope_of_the_curve_is_the_damage_gm(box, midship):
    """``dGZ/dphi`` at upright, from the curve, against ``KB + BM - KG`` of the state."""
    state = flooded_equilibrium(box, midship, DISPLACEMENT, kg=KG)
    gm = damage_gm(state, KG)
    for heel in (1e-2, 1e-3, 1e-4):
        point = flooded_righting_arm(box, midship, DISPLACEMENT, heel_deg=heel, kg=KG)
        slope = point.righting_arm / math.sin(math.radians(heel))
        # the O(phi^2) term is still there at finite heel; it shrinks with the heel
        assert slope == pytest.approx(gm, rel=1e-6)
        assert abs(slope - gm) < 1e-6 * gm * (heel / 1e-2)


def test_the_intact_small_angle_slope_is_the_wall_sided_gm(box):
    point = flooded_righting_arm(box, (), DISPLACEMENT, heel_deg=1e-3, kg=KG)
    slope = point.righting_arm / math.sin(math.radians(1e-3))
    # the wall-sided form adds BM/2 tan^2(phi) = 6e-10 at this heel, so 1e-9 is the
    # tolerance the *series*, not the integrator, sets
    assert slope == pytest.approx(GM0, rel=1e-9)


# --------------------------------------------------------------------------- #
# 3. invariants that hold whatever the mesh
# --------------------------------------------------------------------------- #
def test_volume_and_longitudinal_centre_are_held_at_every_heel(damaged_curve):
    for point in damaged_curve.points:
        assert point.converged
        assert abs(point.volume - DISPLACEMENT) / DISPLACEMENT < 1e-12
        assert point.lcb == pytest.approx(point.lcg, abs=1e-9)
        assert point.kb > 0.0
        assert point.lost_buoyancy_volume > 0.0


def test_a_symmetric_compartment_leaves_no_trim_or_list(box, midship, damaged_curve):
    for point in damaged_curve.points:
        assert point.waterline.trim_tangent == pytest.approx(0.0, abs=1e-10)
        assert point.waterline.heel_tangent == pytest.approx(
            math.tan(math.radians(point.heel_deg)), rel=1e-12
        )


def test_the_arm_is_odd_in_the_heel_for_a_centred_compartment(box, midship):
    port = flooded_righting_arm(box, midship, DISPLACEMENT, heel_deg=-25.0, kg=KG)
    stbd = flooded_righting_arm(box, midship, DISPLACEMENT, heel_deg=25.0, kg=KG)
    assert stbd.righting_arm == pytest.approx(-port.righting_arm, abs=1e-12)
    assert stbd.tcb == pytest.approx(-port.tcb, abs=1e-12)


def test_a_zero_permeability_compartment_is_the_intact_ship(box, midship):
    sealed = Compartment(4.0, 6.0, y0=-1.0, y1=1.0, z0=-0.8, z1=0.0, permeability=0.0)
    for heel in (0.0, 15.0, 40.0):
        flooded = flooded_righting_arm(box, sealed, DISPLACEMENT, heel_deg=heel, kg=KG)
        intact = flooded_righting_arm(box, (), DISPLACEMENT, heel_deg=heel, kg=KG)
        assert flooded.righting_arm == pytest.approx(intact.righting_arm, abs=1e-12)
        assert flooded.volume == pytest.approx(intact.volume, rel=1e-12)


def test_a_side_compartment_lists_the_ship(box):
    """A one-sided compartment makes ``GZ(0)`` non-zero: the ship no longer floats upright.

    At a *prescribed* heel the transverse balance is not solved, so the arm at upright
    is the ship's list, measured -- ``GZ(0) = y_B - y_G``. A centred compartment must
    give exactly zero there; a port-side one must not.
    """
    port = Compartment(
        4.0, 6.0, y0=-1.0, y1=0.0, z0=-0.8, z1=0.0, permeability=0.9, name="port"
    )
    listed = flooded_righting_arm(box, port, DISPLACEMENT, heel_deg=0.0, kg=KG)
    centred = flooded_righting_arm(
        box,
        Compartment(4.0, 6.0, y0=-1.0, y1=1.0, z0=-0.8, z1=0.0, permeability=0.9),
        DISPLACEMENT,
        heel_deg=0.0,
        kg=KG,
    )
    assert centred.righting_arm == pytest.approx(0.0, abs=1e-12)
    assert listed.righting_arm == pytest.approx(listed.tcb, abs=1e-12)
    assert abs(listed.righting_arm) > 1e-2
    # removed buoyancy on the port side pushes B to starboard: the list is to port
    assert listed.tcb > 0.0
    # heeling further to port must therefore *reduce* the arm: the effect has a sign
    deeper = flooded_righting_arm(box, port, DISPLACEMENT, heel_deg=-10.0, kg=KG)
    upright_arm = listed.righting_arm
    assert deeper.righting_arm < upright_arm


def test_an_off_centre_compartment_trims(box):
    """A compartment away from midships drives a trim the two-variable solve must find."""
    forward = Compartment(
        7.0, 9.0, y0=-1.0, y1=1.0, z0=-0.8, z1=0.0, permeability=0.8, name="forward"
    )
    point = flooded_righting_arm(box, forward, DISPLACEMENT, heel_deg=20.0, kg=KG)
    state = flooded_equilibrium(box, forward, DISPLACEMENT, kg=KG)
    assert point.converged
    assert point.lcb == pytest.approx(point.lcg, abs=1e-9)
    assert abs(point.waterline.trim_tangent) > 1e-3
    # the heel-fixed solve at 0 degrees must agree with the three-equation solve
    level = flooded_righting_arm(box, forward, DISPLACEMENT, heel_deg=0.0, kg=KG)
    assert level.waterline.base_z == pytest.approx(state.waterline.base_z, abs=1e-9)
    assert level.waterline.trim_tangent == pytest.approx(
        state.waterline.trim_tangent, abs=1e-9
    )
    assert level.kb == pytest.approx(state.kb, abs=1e-9)


def test_translating_hull_and_compartment_together_changes_nothing(box, midship):
    shifted = box_mesh(x_shift=3.5)
    shifted_compartment = Compartment(
        7.5, 9.5, y0=-1.0, y1=1.0, z0=-0.8, z1=0.0, permeability=0.6, name="midship"
    )
    for heel in (0.0, 20.0, 45.0):
        moved = flooded_righting_arm(shifted, shifted_compartment, DISPLACEMENT, heel_deg=heel, kg=KG)
        stayed = flooded_righting_arm(box, midship, DISPLACEMENT, heel_deg=heel, kg=KG)
        assert moved.righting_arm == pytest.approx(stayed.righting_arm, abs=1e-12)
        assert moved.lcb - 3.5 == pytest.approx(stayed.lcb, abs=1e-9)


# --------------------------------------------------------------------------- #
# 4. the curve as a rule input: peak, range, area -- and the refusal beyond 60 deg
# --------------------------------------------------------------------------- #
def test_the_curve_summary_reads_the_peak_and_the_range(midship):
    """A low-freeboard box loses its arm once the deck edge goes under: it vanishes."""
    shallow = box_mesh(freeboard=0.15)
    steep = flooded_righting_arm_curve(
        shallow, midship, DISPLACEMENT, kg=KG, heels_deg=tuple(range(0, 61, 5))
    )
    assert steep.max_arm > 0.03
    assert 0.0 < steep.heel_at_max_arm_deg <= 60.0
    assert steep.vanishing_angle_deg is not None
    assert steep.heel_at_max_arm_deg < steep.vanishing_angle_deg
    assert steep.positive_range_deg == pytest.approx(steep.vanishing_angle_deg, rel=1e-12)
    assert steep.area_to_max_arm > 0.0
    # the interpolation is between the bracketing ladder points
    arms = {point.heel_deg: point.righting_arm for point in steep.points}
    low = max(h for h in arms if h < steep.vanishing_angle_deg)
    assert arms[low] > 0.0
    assert arms[low + 5.0] <= 0.0


def test_a_curve_that_never_vanishes_says_so(box, midship, damaged_curve):
    assert damaged_curve.vanishing_angle_deg is None
    assert damaged_curve.positive_range_deg == damaged_curve.points[-1].heel_deg
    assert damaged_curve.max_arm == max(point.righting_arm for point in damaged_curve.points)


def test_the_heel_ladder_is_sorted_and_deduplicated(box, midship):
    curve = flooded_righting_arm_curve(
        box, midship, DISPLACEMENT, kg=KG, heels_deg=(30.0, 0.0, 15.0, 15.0)
    )
    assert [point.heel_deg for point in curve.points] == [0.0, 15.0, 30.0]


def test_a_heel_past_the_model_limit_is_refused(box, midship):
    assert MAX_ARM_HEEL_DEG == 60.0
    with pytest.raises(DamageError, match="lost-buoyancy bookkeeping"):
        flooded_righting_arm(box, midship, DISPLACEMENT, heel_deg=61.0, kg=KG)
    with pytest.raises(DamageError, match="displacement_volume must be positive"):
        flooded_righting_arm(box, midship, 0.0, heel_deg=10.0, kg=KG)
    with pytest.raises(DamageError, match="kg must be finite"):
        flooded_righting_arm(box, midship, DISPLACEMENT, heel_deg=10.0, kg=float("nan"))


def test_the_curve_serialises(box, midship, damaged_curve):
    import json

    payload = json.loads(json.dumps(damaged_curve.to_dict()))
    assert payload["schema"] == "tensorcad.damage.flooded_arm/1"
    assert len(payload["righting_arms"]) == len(damaged_curve.points)
    assert payload["vanishing_angle_deg"] is None
    point = payload["righting_arms"][3]
    assert point["heel_deg"] == 15.0
    assert point["converged"] is True


# --------------------------------------------------------------------------- #
# 5. a curved hull: an independent module, and the panel discretisation
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def hull_with_freeboard() -> Mesh:
    model = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13))
    return with_freeboard(model.mesh, 0.4, waterline_z=0.0)


def test_the_panel_discretisation_is_not_the_answer(hull_with_freeboard):
    """Three panel counts over a parabolic hull: the answer must not move."""
    compartment = Compartment(4.0, 6.0, name="midship")
    arms = {}
    for panels in (4, 8, 16):
        curve = flooded_righting_arm_curve(
            hull_with_freeboard,
            compartment,
            5.0,
            kg=0.6,
            panels=panels,
            heels_deg=(10.0, 30.0),
        )
        arms[panels] = [point.righting_arm for point in curve.points]
    for coarse, fine in ((4, 16), (8, 16)):
        for a, b in zip(arms[coarse], arms[fine], strict=True):
            assert a == pytest.approx(b, rel=5e-4)


def test_the_damaged_curve_on_a_curved_hull_conserves_volume(hull_with_freeboard):
    compartment = Compartment(4.0, 6.0, name="midship")
    curve = flooded_righting_arm_curve(
        hull_with_freeboard, compartment, 5.0, kg=0.6, heels_deg=(0.0, 20.0, 40.0)
    )
    for point in curve.points:
        assert point.converged
        assert abs(point.volume - 5.0) / 5.0 < 1e-12
        assert point.lost_buoyancy_volume > 0.0