"""The PR #52 station integrator, kept as an oracle and held to the canonical route.

``tensorcad.gz`` is the canonical ``GZ`` integration and ``tensorcad.righting_arm`` is its
cross-checked second implementation; both integrate the closed surface. The third code path
lives in ``examples/righting/oracle.py`` -- heeled station areas from Green's theorem with
Simpson along the hull -- because a third public route would be a duplicate, while a third
*independent* route is the only thing that can catch a defect shared by the first two.

These tests pin what the oracle is for: the wall-sided barge to round-off (its clip is exact
for a rectangle), the deck-edge limit, the sign and datum conventions, the input validation,
and agreement with ``tensorcad.gz`` on a real hull. They also guard the published evidence
bundle in ``docs/righting/`` against drifting away from the code that produced it.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest

from tensorcad.geometry import Mesh
from tensorcad.gz import (
    box_wall_sided_limit_deg,
    heeled_hydrostatics,
    righting_arm_curve,
    wall_sided_righting_arm,
    with_freeboard,
)
from tensorcad.hull import HullParameters, generate_hull

ROOT = Path(__file__).resolve().parents[1]
ORACLE_PATH = ROOT / "examples" / "righting" / "oracle.py"
EVIDENCE_DIR = ROOT / "docs" / "righting"


def _load_oracle():
    spec = importlib.util.spec_from_file_location("righting_oracle", ORACLE_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("righting_oracle", module)
    spec.loader.exec_module(module)
    return module


oracle = _load_oracle()

BARGE = {"length": 10.0, "beam": 2.0, "draft": 1.0, "freeboard": 1.0, "kg": 0.8}
HULL = {"length": 10.0, "beam": 2.0, "draft": 1.0, "nx": 41, "nz": 13, "kg": 0.5}
HULL_FREEBOARD = 1.0
STEP = 5.0


def barge_mesh(length: float, beam: float, draft: float, freeboard: float) -> Mesh:
    """Watertight rectangular barge: ``draft`` below ``z = 0``, deck at ``z = freeboard``."""
    half = beam / 2.0
    vertices = (
        (0.0, -half, -draft),
        (length, -half, -draft),
        (length, half, -draft),
        (0.0, half, -draft),
        (0.0, -half, freeboard),
        (length, -half, freeboard),
        (length, half, freeboard),
        (0.0, half, freeboard),
    )
    faces = (
        (0, 2, 1), (0, 3, 2),
        (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4),
        (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6),
        (3, 0, 4), (3, 4, 7),
    )
    return Mesh(vertices=vertices, faces=faces)


@pytest.fixture(scope="module")
def barge() -> Mesh:
    return barge_mesh(BARGE["length"], BARGE["beam"], BARGE["draft"], BARGE["freeboard"])


@pytest.fixture(scope="module")
def banded_hull() -> Mesh:
    hull = generate_hull(
        HullParameters(
            length=HULL["length"],
            beam=HULL["beam"],
            draft=HULL["draft"],
            nx=HULL["nx"],
            nz=HULL["nz"],
        )
    ).mesh
    return with_freeboard(hull, HULL_FREEBOARD)


def _barge_reference(angles):
    bm = BARGE["beam"] ** 2 / (12.0 * BARGE["draft"])
    gm = BARGE["draft"] / 2.0 + bm - BARGE["kg"]
    return [wall_sided_righting_arm(angle, gm=gm, bm=bm) for angle in angles]


# ---------------------------------------------------------------- the closed form it must hit
def test_barge_reproduces_the_wall_sided_closed_form_to_round_off(barge):
    """A polygon barge is the case the station integrator is exact on, by construction."""
    limit = box_wall_sided_limit_deg(BARGE["beam"], BARGE["draft"], BARGE["draft"] + BARGE["freeboard"])
    angles = tuple(STEP * index for index in range(int(limit / STEP) + 1))
    curve = oracle.righting_arm_curve(
        barge, kg=BARGE["kg"], heel_angles_deg=angles, station_intervals=8
    )
    assert curve.displacement_volume == pytest.approx(
        BARGE["length"] * BARGE["beam"] * BARGE["draft"], rel=1e-14
    )
    reference = _barge_reference(angles)
    worst = max(abs(point.gz - expected) for point, expected in zip(curve.points, reference, strict=True))
    assert worst < 1e-12, f"wall-sided barge deviates by {worst:.3e} m"


def test_upright_volume_and_centre_of_buoyancy_of_the_barge(barge):
    """The upright point is the datum: ``B`` sits on the centreline at half the draft."""
    body = oracle.immersed_body(barge, 0.0, station_intervals=8)
    assert body.volume == pytest.approx(BARGE["length"] * BARGE["beam"] * BARGE["draft"], rel=1e-14)
    assert body.y_b == pytest.approx(0.0, abs=1e-12)
    assert body.z_b == pytest.approx(BARGE["draft"] / 2.0, rel=1e-12)


def test_closed_form_is_visibly_wrong_just_above_its_window(barge):
    """The agreement above is not vacuous: one step past the limit the formula is out."""
    limit = box_wall_sided_limit_deg(BARGE["beam"], BARGE["draft"], BARGE["draft"] + BARGE["freeboard"])
    angle = limit + STEP
    curve = oracle.righting_arm_curve(
        barge, kg=BARGE["kg"], heel_angles_deg=(0.0, angle), station_intervals=8
    )
    reference = _barge_reference((angle,))[0]
    assert abs(curve.gz[1] - reference) > 1e-3
    # ... and the integrated arm is the smaller of the two: widening the waterplane past the
    # deck edge does not buy arm, so the closed form over-predicts rather than under.
    assert curve.gz[1] < reference


def test_deck_edge_immersion_onset_matches_the_wall_sided_limit(barge):
    limit = box_wall_sided_limit_deg(BARGE["beam"], BARGE["draft"], BARGE["draft"] + BARGE["freeboard"])
    angles = tuple(STEP * index for index in range(1 + round(90.0 / STEP)))
    curve = oracle.righting_arm_curve(
        barge, kg=BARGE["kg"], heel_angles_deg=angles, station_intervals=8
    )
    onset = next(point.heel_deg for point in curve.points if point.deck_immersion_fraction >= 0.5)
    assert abs(onset - limit) <= STEP


# ------------------------------------------------------------- agreement with the canonical route
def test_heeled_body_moves_buoyancy_to_the_heeled_side(barge):
    """Positive heel puts ``+y`` down, so ``B`` moves to ``+y`` and the arm is positive."""
    point = next(
        point
        for point in oracle.righting_arm_curve(
            barge, kg=BARGE["kg"], heel_angles_deg=(0.0, 20.0), station_intervals=8
        ).points
        if point.heel_deg == 20.0
    )
    assert point.transverse_centre_buoyancy > 0.0
    assert point.vertical_centre_buoyancy > BARGE["draft"] / 2.0
    assert point.gz > 0.0
    assert point.submerged_volume == pytest.approx(
        BARGE["length"] * BARGE["beam"] * BARGE["draft"], rel=1e-10
    )


def test_oracle_agrees_with_the_canonical_route_on_the_hull(banded_hull):
    """Two different integrations of one body, within a fraction of the arm."""
    angles = (0.0, 20.0, 40.0, 60.0)
    canonical = righting_arm_curve(banded_hull, HULL["kg"], list(angles), waterline_z=0.0)
    section = oracle.righting_arm_curve(
        banded_hull, kg=HULL["kg"], heel_angles_deg=angles, station_intervals=64
    )
    worst = max(abs(left - right) for left, right in zip(section.gz, canonical.gz, strict=True))
    assert worst / canonical.max_gz < 5e-4, f"routes differ by {worst:.3e} m"
    assert section.displacement_volume == pytest.approx(
        canonical.displacement_volume, rel=1e-3
    )
    assert section.max_volume_error < 1e-9


def test_upright_displacement_agrees_with_the_closed_surface_integration(banded_hull):
    """The station sum and the divergence-theorem sum of the same body, at zero heel."""
    section = oracle.immersed_body(banded_hull, 0.0, station_intervals=64)
    closed = heeled_hydrostatics(banded_hull, 0.0, waterline_z=0.0)
    assert section.volume == pytest.approx(closed.displacement_volume, rel=1e-3)
    assert section.z_b == pytest.approx(closed.centre_of_buoyancy_z, rel=1e-3)


def test_station_refinement_converges(banded_hull):
    """The longitudinal quadrature is the error source, and it is first order at worst."""
    angles = (0.0, 40.0)
    coarse = oracle.righting_arm_curve(
        banded_hull, kg=HULL["kg"], heel_angles_deg=angles, station_intervals=16
    )
    fine = oracle.righting_arm_curve(
        banded_hull, kg=HULL["kg"], heel_angles_deg=angles, station_intervals=64
    )
    finer = oracle.righting_arm_curve(
        banded_hull, kg=HULL["kg"], heel_angles_deg=angles, station_intervals=128
    )
    error_64 = abs(fine.gz[1] - finer.gz[1])
    error_16 = abs(coarse.gz[1] - finer.gz[1])
    assert error_64 < 5e-5
    assert error_16 > 4.0 * error_64


# --------------------------------------------------------------------------- input validation
def test_station_intervals_below_the_floor_are_rejected(barge):
    with pytest.raises(oracle.OracleError, match="at least"):
        oracle.section_grid(barge, 2)


def test_odd_station_intervals_are_rounded_up(barge):
    assert oracle.section_grid(barge, 9).intervals == 10


def test_heel_angles_must_start_at_zero_and_increase(barge):
    with pytest.raises(oracle.OracleError, match="include 0 degrees"):
        oracle.righting_arm_curve(barge, kg=0.8, heel_angles_deg=(5.0, 10.0))
    with pytest.raises(oracle.OracleError, match="strictly increasing"):
        oracle.righting_arm_curve(barge, kg=0.8, heel_angles_deg=(0.0, 10.0, 5.0))
    with pytest.raises(oracle.OracleError, match="at least one heel angle"):
        oracle.righting_arm_curve(barge, kg=0.8, heel_angles_deg=())


def test_a_waterline_at_the_deck_is_rejected(barge):
    with pytest.raises(oracle.OracleError, match="floods the deck"):
        oracle.righting_arm_curve(
            barge, kg=0.8, heel_angles_deg=(0.0, 10.0), waterline_z=BARGE["freeboard"]
        )


def test_a_displacement_the_body_cannot_hold_is_reported_not_invented():
    """Ask for more displacement than the body has volume: the solve must refuse, not clamp."""
    box = barge_mesh(1.0, 1.0, 1.0, 1.0)  # 1 x 1 x 2 m solid, 2 m3
    with pytest.raises(oracle.OracleError, match="fully immersed"):
        oracle.immersed_body(box, 30.0, station_intervals=8, displacement_volume=2.5)


def test_non_finite_inputs_are_rejected(barge):
    with pytest.raises(oracle.OracleError, match="finite"):
        oracle.righting_arm_curve(barge, kg=float("nan"), heel_angles_deg=(0.0, 10.0))
    with pytest.raises(oracle.OracleError, match="finite"):
        oracle.immersed_body(barge, float("inf"))


# ------------------------------------------------------------------------- the evidence bundle
def test_published_evidence_bundle_matches_a_fresh_evaluation(barge):
    """``docs/righting`` is generated, so it is checked against the code on every run."""
    results = json.loads((EVIDENCE_DIR / "results.json").read_text())
    assert results["checks"]["passed"] is True, results["checks"]["failures"]
    assert all(claim["passed"] for claim in results["claims"])
    for name in ("report.md", "righting_arms.png"):
        assert (EVIDENCE_DIR / name).stat().st_size > 0

    barge_record = results["barge"]
    limit = barge_record["wall_sided_limit_deg"]
    angles = tuple(row["heel_deg"] for row in barge_record["curve"] if row["heel_deg"] <= limit)
    curve = oracle.righting_arm_curve(
        barge, kg=BARGE["kg"], heel_angles_deg=angles, station_intervals=8
    )
    reference = _barge_reference(angles)
    fresh = max(abs(point.gz - expected) for point, expected in zip(curve.points, reference, strict=True))
    recorded = max(
        row["abs_deviation_m"] for row in barge_record["curve"] if row["heel_deg"] <= limit
    )
    assert fresh == pytest.approx(recorded, abs=1e-15)
    assert barge_record["displacement_volume_m3"] == pytest.approx(
        BARGE["length"] * BARGE["beam"] * BARGE["draft"], rel=1e-14
    )
    assert math.isclose(limit, 45.0, abs_tol=1e-9)
    assert results["hull"]["station_intervals"] == 64
