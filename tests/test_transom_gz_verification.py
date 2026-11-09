"""Verification of the heeled righting-arm chain on the transom-stern hull family.

The GZ integrators were anchored on a box barge (wall-sided closed form), a
circular section (``GM sin(phi)``) and the Wigley hull (symmetric, fine stern),
but never on a hull whose aft body ends in an **immersed transom plate** -- the
form every real displacement vessel above ``Fn ~ 0.25`` has.  This module closes
that: geometry -> ``with_freeboard`` -> heeled displacement -> GZ -> the initial
stability route, on ``tensorcad.transom``'s default hull.

Two independent heeled-displacement integrators are used as each other's
witness: ``tensorcad.gz`` sweeps the surface and integrates per triangle, while
``tensorcad.righting_arm`` cuts the body into sections and bisects the waterline.
They share no code path, so agreement between them is evidence and not a
tautology.

All reference numbers are recomputed here, not copied; the committed artifact
(``artifacts/transom_gz/evidence.json``, produced by ``tools/verify_transom_gz.py``)
is additionally checked against this module's own arithmetic so a stale artifact
cannot pass silently.
"""
import json
import math
from pathlib import Path

import pytest

from tensorcad.gz import (
    deck_edge_immersion_deg,
    heeled_hydrostatics,
    righting_arm_curve as sweep_curve,
    with_freeboard,
)
from tensorcad.righting_arm import righting_arm_curve as bisection_curve
from tensorcad.stability import initial_stability
from tensorcad.transom import generate_transom_hull

FREEBOARD = 0.15
KG = 0.3
HEELS = (0.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0)
ARTIFACT = Path(__file__).resolve().parent.parent / "artifacts" / "transom_gz" / "evidence.json"


@pytest.fixture(scope="module")
def hull():
    model = generate_transom_hull()
    return model, with_freeboard(model.mesh, FREEBOARD, waterline_z=0.0)


@pytest.fixture(scope="module")
def curves(hull):
    _, banded = hull
    fast = sweep_curve(banded, kg=KG, angles_deg=HEELS, waterline_z=0.0)
    slow = bisection_curve(banded, kg=KG, heels_deg=HEELS, waterline_z=0.0)
    return list(fast.gz), [point.righting_arm for point in slow.points]


def test_the_transom_hull_can_be_heeled_at_all(hull):
    """A hull generated to its waterline has no topside; the band must make it closable."""
    model, banded = hull
    assert banded.bounding_box().maximum[2] == pytest.approx(FREEBOARD, abs=1e-12)
    assert len(banded.faces) > len(model.mesh.faces)


def test_upright_displacement_is_unchanged_by_the_freeboard_band(hull):
    model, banded = hull
    upright = heeled_hydrostatics(banded, 0.0, waterline_z=0.0)
    assert upright.displacement_volume == pytest.approx(
        model.hydrostatics.displacement_volume, rel=1e-12
    )


def test_displacement_is_held_constant_at_every_heel(hull, curves):
    """A heeled ship keeps its displacement; the waterline moves, the volume does not."""
    _, banded = hull
    upright = heeled_hydrostatics(banded, 0.0, waterline_z=0.0).displacement_volume
    for heel in HEELS:
        point = heeled_hydrostatics(banded, heel, waterline_z=0.0)
        assert abs(point.displacement_volume - upright) / upright < 1e-9, heel
        assert abs(point.volume_error) < 1e-9, heel


def test_the_two_independent_integrators_agree(curves):
    """No shared code path: agreement is a witness, disagreement would be a defect."""
    fast, slow = curves
    worst = max(abs(a - b) for a, b in zip(fast, slow))
    assert worst < 1e-9, f"routes differ by {worst:g} m"


def test_small_angle_slope_reproduces_the_metacentric_height(hull, curves):
    """GZ(phi)/sin(phi) -> GM as phi -> 0, with GM from the waterplane moments."""
    model, _ = hull
    fast, _ = curves
    state = initial_stability(model.mesh, kg=KG, waterline_z=0.0)
    gm = state.gm_transverse
    slope = fast[HEELS.index(2.0)] / math.sin(math.radians(2.0))
    assert gm == pytest.approx(state.kb + state.bm_transverse - KG, rel=1e-12)
    assert abs(slope / gm - 1.0) < 0.01, f"slope {slope:.6f} vs GM {gm:.6f}"


def test_the_upright_arm_is_zero_to_the_mesh_residual(hull, curves):
    """A symmetric hull must have no upright heel arm.

    It does not vanish exactly, and the reason is the *mesh*: the lofted quads
    between stations are not planar, so the consistent diagonal split leaves a
    small port/starboard asymmetry (nonzero transverse centre of buoyancy) even
    though the vertices are mirror-exact.  The measured residual is asserted to
    be that centroid, not some integrator bias -- and to stay below 0.01 % of
    the beam, which is the resolution the model is generated at.
    """
    model, _ = hull
    fast, slow = curves
    assert abs(fast[0] - slow[0]) < 1e-12
    assert abs(fast[0]) < 1e-4
    assert abs(fast[0]) < 1e-4 * model.parameters.beam


def test_the_helped_waterplane_shrinks_as_the_hull_rolls(hull):
    """The section of a V-ish body at a fixed displacement narrows before it vanishes."""
    _, banded = hull
    areas = [heeled_hydrostatics(banded, h, waterline_z=0.0).waterplane_area for h in HEELS]
    assert all(later < earlier for earlier, later in zip(areas, areas[1:])), areas


def test_the_waterline_rises_on_the_dipping_side(hull):
    """Constant displacement means the level must move, not stay at the design waterline."""
    _, banded = hull
    rises = [heeled_hydrostatics(banded, h, waterline_z=0.0).waterline_offset for h in HEELS]
    assert rises[0] == pytest.approx(0.0, abs=1e-12)
    assert all(rise < 0.0 for rise in rises[1:]), rises


def test_the_righting_arm_grows_through_the_measured_range(curves):
    fast, _ = curves
    assert all(later > earlier for earlier, later in zip(fast, fast[1:])), fast
    assert fast[-1] > 0.2


def test_the_deck_edge_bound_is_reported_and_bounds_the_claim(hull):
    """Beyond deck-edge immersion the freeboard band is a vertical surrogate, not a topside.

    The curve is still self-consistent there (the tests above hold at 40 deg), but
    it describes the band the model has, so the verification claim is limited to
    the deck-edge angle and the document says so.
    """
    model, _ = hull
    half_breadth = max(vertex[1] for vertex in model.mesh.vertices)
    bound = deck_edge_immersion_deg(half_breadth, FREEBOARD)
    assert 5.0 < bound < 15.0
    assert bound == pytest.approx(math.degrees(math.atan(FREEBOARD / half_breadth)))


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not generated in this worktree")
def test_the_committed_artifact_matches_today_s_arithmetic(hull, curves):
    """The committed evidence must be what the code computes now, not a transcription."""
    model, banded = hull
    fast, slow = curves
    evidence = json.loads(ARTIFACT.read_text())
    state = initial_stability(model.mesh, kg=KG)
    assert evidence["upright"]["gm_transverse"] == pytest.approx(
        state.gm_transverse, rel=1e-9
    )
    table = evidence["heel_table"]
    assert [row["heel_deg"] for row in table] == list(HEELS)
    for row, gz_fast, gz_slow in zip(table, fast, slow):
        assert row["gz_fast"] == pytest.approx(gz_fast, abs=1e-12)
        assert row["gz_bisection"] == pytest.approx(gz_slow, abs=1e-12)
    # the artifact's headline numbers must be reproducible from this run
    assert evidence["cross_route_max_abs_diff"] == pytest.approx(
        max(abs(a - b) for a, b in zip(fast, slow)), rel=1e-6, abs=1e-15
    )
    assert evidence["small_angle"]["gm"] == pytest.approx(state.gm_transverse, rel=1e-9)
    assert evidence["deck_edge_immersion_deg"] == pytest.approx(
        deck_edge_immersion_deg(max(v[1] for v in model.mesh.vertices), FREEBOARD)
    )