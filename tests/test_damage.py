"""Damage stability: the box closed forms, the mirror, and a mesh-free hull baseline.

The standard the module is held to is the one the rest of the stability stack uses: a
wall-sided box has no discretisation, so every number has an exact value; a hull with
curvature is checked against the closed-form section area of its own family, integrated
without the mesh.
"""
from __future__ import annotations

import json
import math

import pytest

from tensorcad import HullParameters, analytic_volume, generate_hull
from tensorcad.bonjean import analytic_section_area
from tensorcad.damage import (
    DamageError,
    DamageWaterline,
    Compartment,
    buoyancy_curve,
    damage_gm,
    flooded_equilibrium,
    intact_equilibrium,
)
from tensorcad.geometry import Mesh

L, B, T = 10.0, 2.0, 0.8
DISPLACEMENT = L * B * T  # 16 m^3 of box under the intact waterline


def box_mesh(length: float = L, beam: float = B, draft: float = T, freeboard: float = 1.0) -> Mesh:
    """A watertight rectangular hull, keel at ``-draft``, deck at ``+freeboard``."""
    xs = [0.0, length]
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
    return Compartment(4.0, 6.0, name="midship")


# --------------------------------------------------------------------------- #
# 1. the intact ship, because a damage solver that cannot reproduce it is wrong
# --------------------------------------------------------------------------- #
def test_the_intact_ship_lands_on_its_design_waterline_and_closed_forms(box):
    state = intact_equilibrium(box, DISPLACEMENT)
    assert state.waterline.base_z == pytest.approx(0.0, abs=1e-12)
    assert state.volume_net == pytest.approx(DISPLACEMENT, rel=1e-12)
    assert state.lcb == pytest.approx(0.5 * L, rel=1e-12)
    assert state.tcb == pytest.approx(0.0, abs=1e-12)
    assert state.draft_at_reference == pytest.approx(T, rel=1e-12)
    assert state.kb == pytest.approx(0.5 * T, rel=1e-12)
    assert state.transverse_second_moment == pytest.approx(B**3 * L / 12.0, rel=1e-12)
    assert state.bm_transverse == pytest.approx(B**2 / (12.0 * T), rel=1e-12)
    assert state.trim_angle_deg == 0.0
    assert state.heel_angle_deg == 0.0
    assert state.converged
    assert not state.damaged


# --------------------------------------------------------------------------- #
# 2. a midship compartment: nothing to discretise, so everything is a closed form
# --------------------------------------------------------------------------- #
def test_a_midship_compartment_gives_the_textbook_draft(box, midship):
    """``T' = L T / (L - l)``: the intact waterplane shrinks by the compartment's length."""
    state = flooded_equilibrium(box, midship, DISPLACEMENT)
    expected = L * T / (L - midship.length)
    assert state.draft_at_reference == pytest.approx(expected, rel=1e-12)
    assert state.trim_angle_deg == pytest.approx(0.0, abs=1e-9)
    assert state.heel_angle_deg == pytest.approx(0.0, abs=1e-9)
    assert state.lcb == pytest.approx(0.5 * L, rel=1e-12)
    assert state.sinkage == pytest.approx(expected - T, rel=1e-9)
    assert state.converged
    assert state.volume_residual == pytest.approx(0.0, abs=1e-12)


def test_the_lost_buoyancy_is_the_flooded_volume_and_says_so(box, midship):
    state = flooded_equilibrium(box, midship, DISPLACEMENT)
    draft = L * T / (L - midship.length)
    assert state.lost_buoyancy_volume == pytest.approx(midship.length * B * draft, rel=1e-12)
    # the external volume is the displacement plus exactly what was lost
    assert state.volume_net + state.lost_buoyancy_volume == pytest.approx(
        DISPLACEMENT + midship.length * B * draft, rel=1e-12
    )


def test_the_damaged_waterplane_is_the_intact_one_minus_the_compartment(box, midship):
    state = flooded_equilibrium(box, midship, DISPLACEMENT)
    assert state.waterplane_area == pytest.approx(B * (L - midship.length), rel=1e-12)
    assert state.flooded_waterplane_area == pytest.approx(B * midship.length, rel=1e-12)
    assert state.transverse_second_moment == pytest.approx(
        B**3 * (L - midship.length) / 12.0, rel=1e-12
    )


def test_permeability_scales_the_lost_buoyancy(box):
    """``mu = 0.6`` (a cargo hold): ``T' = L T / (L - mu l)``."""
    mu = 0.6
    state = flooded_equilibrium(box, Compartment(4.0, 6.0, permeability=mu), DISPLACEMENT)
    assert state.draft_at_reference == pytest.approx(L * T / (L - mu * 2.0), rel=1e-12)
    assert state.lost_buoyancy_volume == pytest.approx(
        mu * 2.0 * B * state.draft_at_reference, rel=1e-12
    )


def test_metacentric_height_after_damage_and_a_negative_gm_is_reported(box, midship):
    kg = 1.2
    state = flooded_equilibrium(box, midship, DISPLACEMENT, kg=kg)
    draft = L * T / (L - midship.length)
    expected_kb = 0.5 * draft
    expected_bm = (B**3 * (L - 2.0) / 12.0) / DISPLACEMENT
    assert state.kb == pytest.approx(expected_kb, rel=1e-12)
    assert state.bm_transverse == pytest.approx(expected_bm, rel=1e-12)
    assert state.gm_transverse == pytest.approx(expected_kb + expected_bm - kg, rel=1e-9)
    # a damaged box barge with the weight this high has negative GM, and it is not clamped
    assert state.gm_transverse < 0.0
    assert damage_gm(state, kg) == pytest.approx(state.gm_transverse, rel=0, abs=0)
    # without KG the GM is unknown, not guessed
    assert flooded_equilibrium(box, midship, DISPLACEMENT).gm_transverse is None


# --------------------------------------------------------------------------- #
# 3. asymmetric damage: trim, heel, and the analytic box solution of each
# --------------------------------------------------------------------------- #
def test_an_end_compartment_trims_by_the_analytic_tangent(box):
    """A two-equation linear solve on the box gives ``z0 = 1.1875 m``, ``tan = -0.1875``.

    Volume: ``16 z0 + 16 tan = 16``; moment: ``96 z0 + 181.3333 tan = 80`` (LCG at
    midships, ``x_ref = 5``), both exact for a rectangle.
    """
    state = flooded_equilibrium(box, Compartment(0.0, 2.0, name="aft"), DISPLACEMENT)
    assert state.draft_at_reference == pytest.approx(1.1875, rel=1e-12)
    assert state.waterline.trim_tangent == pytest.approx(-0.1875, rel=1e-12)
    assert state.trim_angle_deg == pytest.approx(math.degrees(math.atan(-0.1875)), rel=1e-12)
    assert state.heel_angle_deg == pytest.approx(0.0, abs=1e-9)
    assert state.lcb == pytest.approx(0.5 * L, rel=1e-12)
    # the flooded end is the deep end: the deepest draft is at x = 0
    assert state.draft_max == pytest.approx(1.1875 + 0.1875 * 5.0, rel=1e-12)
    assert state.draft_min == pytest.approx(1.1875 - 0.1875 * 5.0, rel=1e-12)


def test_a_side_compartment_finds_its_own_heel(box):
    """A port-side box compartment: ``18 z0 + 1 tan = 16`` and ``z0 + 6 tan = 0``."""
    state = flooded_equilibrium(
        box, Compartment(4.0, 6.0, y0=-1.0, y1=0.0, name="port"), DISPLACEMENT
    )
    z0 = 16.0 / (18.0 - 1.0 / 6.0)
    tangent = -z0 / 6.0
    assert state.draft_at_reference == pytest.approx(z0, rel=1e-12)
    assert state.waterline.heel_tangent == pytest.approx(tangent, rel=1e-12)
    assert state.heel_angle_deg == pytest.approx(math.degrees(math.atan(tangent)), rel=1e-12)
    assert state.trim_angle_deg == pytest.approx(0.0, abs=1e-9)
    assert state.tcb == pytest.approx(0.0, abs=1e-12)


def test_the_starboard_mirror_heels_the_other_way_by_the_same_angle(box):
    port = flooded_equilibrium(box, Compartment(4.0, 6.0, y0=-1.0, y1=0.0), DISPLACEMENT)
    starboard = flooded_equilibrium(box, Compartment(4.0, 6.0, y0=0.0, y1=1.0), DISPLACEMENT)
    assert starboard.heel_angle_deg == pytest.approx(-port.heel_angle_deg, rel=1e-12)
    assert starboard.draft_at_reference == pytest.approx(port.draft_at_reference, rel=1e-12)


def test_a_compartment_whose_top_is_under_the_waterline_loses_only_its_height(box):
    """A hold from the keel to ``z = -0.2``: ``20 z0 - 2.4 = 16``, so ``z0 = 0.92``."""
    state = flooded_equilibrium(
        box, Compartment(4.0, 6.0, z0=-T, z1=-0.2, name="lower hold"), DISPLACEMENT
    )
    assert state.draft_at_reference == pytest.approx(0.92, rel=1e-12)
    assert state.lost_buoyancy_volume == pytest.approx(2.0 * B * 0.6, rel=1e-12)


# --------------------------------------------------------------------------- #
# 4. the identities the formulation has to satisfy
# --------------------------------------------------------------------------- #
def test_a_zero_permeability_compartment_reproduces_the_intact_ship_exactly(box):
    """``mu = 0`` floods nothing; the damaged solver must return the intact answer."""
    intact = intact_equilibrium(box, DISPLACEMENT)
    state = flooded_equilibrium(box, Compartment(4.0, 6.0, permeability=0.0), DISPLACEMENT)
    assert state.draft_at_reference == pytest.approx(intact.draft_at_reference, rel=1e-12)
    assert state.lcb == pytest.approx(intact.lcb, rel=1e-12)
    assert state.lost_buoyancy_volume == 0.0 or state.lost_buoyancy_volume < 1e-15
    assert state.transverse_second_moment == pytest.approx(
        intact.transverse_second_moment, rel=1e-12
    )


def test_two_adjacent_compartments_are_one_compartment(box):
    together = flooded_equilibrium(box, Compartment(4.0, 6.0), DISPLACEMENT)
    split = flooded_equilibrium(
        box, (Compartment(4.0, 5.0, name="fwd half"), Compartment(5.0, 6.0, name="aft half")),
        DISPLACEMENT,
    )
    assert split.draft_at_reference == pytest.approx(together.draft_at_reference, rel=1e-12)
    assert split.lost_buoyancy_volume == pytest.approx(together.lost_buoyancy_volume, rel=1e-12)


def test_a_bigger_compartment_sinks_the_ship_further(box):
    drafts = [
        flooded_equilibrium(box, Compartment(4.0, 4.0 + length), DISPLACEMENT).draft_at_reference
        for length in (0.5, 1.0, 2.0, 4.0)
    ]
    assert drafts == sorted(drafts)
    assert drafts[0] < drafts[-1]


def test_the_integration_grid_is_exact_for_a_box(box, midship):
    """Gauss panels make the box answer panel-count independent, to round-off."""
    drafts = [
        flooded_equilibrium(box, midship, DISPLACEMENT, panels=panels).draft_at_reference
        for panels in (1, 2, 4, 12, 32)
    ]
    assert max(drafts) - min(drafts) < 1e-12


# --------------------------------------------------------------------------- #
# 5. refusals: a case that cannot be posed must not produce a number
# --------------------------------------------------------------------------- #
def test_a_compartment_off_the_hull_is_refused(box):
    with pytest.raises(DamageError, match="outside the hull"):
        flooded_equilibrium(box, Compartment(20.0, 22.0), DISPLACEMENT)


def test_overlapping_compartments_are_refused(box):
    with pytest.raises(DamageError, match="overlap"):
        flooded_equilibrium(
            box, (Compartment(2.0, 5.0, name="a"), Compartment(4.0, 6.0, name="b")),
            DISPLACEMENT,
        )


def test_a_bad_compartment_is_refused_at_construction(box):
    with pytest.raises(DamageError, match="must exceed"):
        Compartment(6.0, 4.0)
    with pytest.raises(DamageError, match="permeability"):
        Compartment(4.0, 6.0, permeability=1.5)
    with pytest.raises(DamageError, match="permeability"):
        Compartment(4.0, 6.0, permeability=-0.1)


def test_flooding_the_whole_waterplane_away_is_refused(box):
    with pytest.raises(DamageError, match="vanished|intact waterplane"):
        flooded_equilibrium(box, Compartment(0.0, L), DISPLACEMENT)


def test_a_displacement_larger_than_the_hull_is_refused(box):
    with pytest.raises(DamageError, match="submarine"):
        intact_equilibrium(box, 100.0)


def test_values_must_be_finite(box, midship):
    with pytest.raises(DamageError, match="finite"):
        flooded_equilibrium(box, midship, float("nan"))
    with pytest.raises(DamageError, match="finite"):
        flooded_equilibrium(box, midship, DISPLACEMENT, lcg=float("inf"))
    with pytest.raises(DamageError, match="positive"):
        flooded_equilibrium(box, midship, -1.0)


# --------------------------------------------------------------------------- #
# 6. the buoyancy curve, and the record
# --------------------------------------------------------------------------- #
def test_the_buoyancy_curve_cuts_the_compartment_out(box, midship):
    state = flooded_equilibrium(box, midship, DISPLACEMENT)
    curve = buoyancy_curve(box, midship, waterline=state.waterline, stations=21)
    assert len(curve.stations) == 21
    for station, external, net in zip(
        curve.stations, curve.external_areas, curve.net_areas, strict=True
    ):
        assert external == pytest.approx(B * state.draft_at_reference, rel=1e-12)
        if 4.0 <= station <= 6.0:  # the boundary itself is measure zero, either way
            assert net == pytest.approx(0.0, abs=1e-12)
        else:
            assert net == pytest.approx(external, rel=1e-12)


def test_the_buoyancy_curve_needs_a_waterline(box):
    with pytest.raises(DamageError, match="no default waterline"):
        buoyancy_curve(box)
    curve = buoyancy_curve(box, displacement_volume=DISPLACEMENT)
    assert curve.external_areas[0] == pytest.approx(B * T, rel=1e-12)


def test_the_record_serialises_with_its_schema(box, midship):
    state = flooded_equilibrium(box, midship, DISPLACEMENT, kg=0.9)
    record = state.to_dict()
    assert record["schema"] == "tensorcad.damage.flooded_state/1"
    assert record["converged"] is True
    assert record["compartments"][0]["permeability"] == 1.0
    assert record["waterline"]["trim_angle_deg"] == pytest.approx(0.0, abs=1e-9)
    json.dumps(record)  # must be JSON, not "JSON-ish"


# --------------------------------------------------------------------------- #
# 7. a real curved hull: the same integrals, against a different code path
# --------------------------------------------------------------------------- #
def _simpson(function, low: float, high: float, intervals: int = 2000) -> float:
    if intervals % 2:
        intervals += 1
    step = (high - low) / intervals
    total = function(low) + function(high)
    for index in range(1, intervals):
        total += (4.0 if index % 2 else 2.0) * function(low + index * step)
    return total * step / 3.0


def test_a_curved_hull_reproduces_its_own_hydrostatics():
    """The intact solve on a curved hull against ``tensorcad.hydrostatics``.

    Two independent integrals of the same mesh -- one over plane sections, one over the
    wetted surface -- have to agree at the waterline they agree on. ``hydrostatics``
    fixes the waterline and reports the volume; the damage solver does the reverse, so
    this pins both the section integrator and the bisection that inverts it.
    """
    from tensorcad.hydrostatics import compute_hydrostatics

    mesh = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13)).mesh
    reference = compute_hydrostatics(mesh)
    state = intact_equilibrium(mesh, reference.displacement_volume)
    expected_draft = reference.waterline_z - reference.keel_z
    assert state.draft_at_reference == pytest.approx(expected_draft, rel=1e-2)
    assert state.volume_net == pytest.approx(reference.displacement_volume, rel=1e-12)
    assert state.lcb == pytest.approx(reference.longitudinal_centre_buoyancy, rel=3e-2)
    # hydrostatics reports the centres in the mesh frame; the damage module reports
    # KB above the keel, so the two are offset by exactly keel_z
    assert state.kb == pytest.approx(
        reference.vertical_centre_buoyancy - reference.keel_z, rel=3e-2
    )
    assert state.tcb == pytest.approx(0.0, abs=1e-6)
    assert state.transverse_second_moment == pytest.approx(
        reference.transverse_waterplane_inertia, rel=5e-2
    )
    assert state.converged


def test_the_net_waterplane_subtracts_permeability_times_the_flooded_area(box):
    """The lost-buoyancy method's answer to the free-surface question, at ``mu != 1``.

    Found by the 2026-11-08 gate-power audit (`docs/analysis/2026-11-08-regulatory-gate-power.md`):
    every waterplane/GM assertion in this file used the default ``mu = 1.0``, so dropping the
    permeability factor from the flooded waterplane survived the suite with ``I_T`` 9.09 % high.
    ``mu < 1`` is the case the "no separate free-surface correction" claim is about.
    """
    mu = 0.6
    state = flooded_equilibrium(box, Compartment(4.0, 6.0, permeability=mu), DISPLACEMENT, kg=1.0)
    assert state.waterplane_area == pytest.approx(B * (L - mu * 2.0), rel=1e-12)
    assert state.flooded_waterplane_area == pytest.approx(mu * B * 2.0, rel=1e-12)
    assert state.transverse_second_moment == pytest.approx(
        B**3 * (L - mu * 2.0) / 12.0, rel=1e-12
    )
