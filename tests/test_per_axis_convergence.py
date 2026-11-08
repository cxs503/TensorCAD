"""Per-axis convergence study of the analytic hull: the claims must stay true.

The committed artifact ``artifacts/per_axis_convergence/per_axis_convergence.json`` is
written by ``tools/per_axis_convergence.py``.  These tests (a) re-derive the study's
conclusions on a *reduced* ladder, so the assertions are not a re-reading of the file,
and (b) pin the committed artifact to its own claims, including the one hypothesis the
study falsifies -- a falsified hypothesis that silently starts passing is also a change
of result and has to be noticed.

Reduced ladder: axial ``nx in (21, 41, 81)`` at ``nz = 41``; vertical
``nz in (9, 13, 21, 41)`` at ``nx = 161``.  The reference grid is 161x41; the study's own
reference is 321x81, which is too slow for a test.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad import HullParameters, analytic_hydrostatics, generate_hull
from tensorcad.hydrostatics import compute_hydrostatics

ARTIFACT = Path(__file__).resolve().parents[1] / "artifacts/per_axis_convergence/per_axis_convergence.json"

REDUCED_AXIAL = (21, 41, 81)
REDUCED_VERTICAL = (9, 13, 21, 41)
N_Z_REF = 41
N_X_REF = 161


def _errors(nx: int, nz: int) -> dict[str, float]:
    parameters = HullParameters(nx=nx, nz=nz)
    reference = analytic_hydrostatics(parameters)
    measured = compute_hydrostatics(generate_hull(parameters).mesh,
                                    water_density=parameters.water_density)
    return {
        "volume": (measured.displacement_volume - reference.displacement_volume)
        / reference.displacement_volume,
        "waterplane_area": (measured.waterplane_area - reference.waterplane_area)
        / reference.waterplane_area,
        "vcb": (measured.vertical_centre_buoyancy - reference.vertical_centre_buoyancy)
        / reference.vertical_centre_buoyancy,
        "section": (measured.maximum_section_area - reference.maximum_section_area)
        / reference.maximum_section_area,
    }


def _slope(hs: list[float], deltas: list[float]) -> float:
    points = [(h, abs(d)) for h, d in zip(hs, deltas) if abs(d) > 1e-13]
    xs = [math.log(h) for h, _ in points]
    ys = [math.log(e) for _, e in points]
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sum((x - mx) ** 2 for x in xs)


@pytest.fixture(scope="module")
def reduced() -> dict[str, object]:
    baseline = _errors(N_X_REF, N_Z_REF)
    axial = {nx: _errors(nx, N_Z_REF) for nx in REDUCED_AXIAL}
    vertical = {nz: _errors(N_X_REF, nz) for nz in REDUCED_VERTICAL}
    return {"baseline": baseline, "axial": axial, "vertical": vertical}


def test_axial_order_is_the_endpoint_exponent_plus_one(reduced):
    hs = [10.0 / (nx - 1) for nx in REDUCED_AXIAL]
    for name in ("volume", "waterplane_area"):
        deltas = [reduced["axial"][nx][name] - reduced["baseline"][name] for nx in REDUCED_AXIAL]
        order = _slope(hs, deltas)
        # predicted 1.7 = BREADTH_EXPONENT (0.7) + 1 for trapezoid error on f ~ x**a
        assert 1.45 < order < 2.05, (name, order)


def test_vertical_order_is_the_keel_exponent_plus_one(reduced):
    """The keel ``sqrt`` taper costs about half an order: p = 0.5 + 1 = 1.5, not 2.

    Coarse levels scatter (the measured 13 -> 21 step is the flattest, 21 -> 41 the
    steepest), so this test asserts the envelope and monotonicity; the exact pinned
    values live in ``TestCommittedArtifact`` on the full 321x81 ladder.
    """
    hs = [0.8 / (nz - 1) for nz in REDUCED_VERTICAL]
    for name in ("volume", "vcb", "section"):
        deltas = [reduced["vertical"][nz][name] - reduced["baseline"][name] for nz in REDUCED_VERTICAL]
        magnitudes = [abs(d) for d in deltas]
        assert all(b <= a * (1.0 + 1e-9) for a, b in zip(magnitudes, magnitudes[1:])), name
        order = _slope(hs, deltas)
        assert 1.2 < order < 2.0, (name, order)


def test_quantities_are_blind_to_the_axis_they_do_not_integrate(reduced):
    for nz in REDUCED_VERTICAL:
        assert reduced["vertical"][nz]["waterplane_area"] == reduced["baseline"]["waterplane_area"]
    for nx in REDUCED_AXIAL:
        assert reduced["axial"][nx]["vcb"] == pytest.approx(reduced["baseline"]["vcb"], abs=1e-12)
        assert reduced["axial"][nx]["section"] == pytest.approx(reduced["baseline"]["section"],
                                                               abs=1e-12)


def test_the_two_axis_errors_add(reduced):
    """Additive model, evaluated at the reduced reference: E(41,13) = E_ref + dx(41) + dz(13)."""
    measured = _errors(41, 13)["volume"]
    predicted = (reduced["baseline"]["volume"]
                 + (reduced["axial"][41]["volume"] - reduced["baseline"]["volume"])
                 + (reduced["vertical"][13]["volume"] - reduced["baseline"]["volume"]))
    assert abs(measured - predicted) <= 0.25 * abs(measured)


def test_vertical_axis_owns_the_default_grid_error(reduced):
    default = _errors(41, 13)["volume"]
    axial_part = reduced["axial"][41]["volume"] - reduced["baseline"]["volume"]
    vertical_part = reduced["vertical"][13]["volume"] - reduced["baseline"]["volume"]
    assert abs(vertical_part) > abs(axial_part)
    assert abs(vertical_part) > 0.5 * abs(default)


def test_default_grid_errors_are_within_the_audit_band(reduced):
    default = _errors(41, 13)
    for name in ("volume", "waterplane_area", "vcb"):
        assert abs(default[name]) < 0.03, (name, default[name])


@pytest.mark.skipif(not ARTIFACT.exists(), reason="artifact not generated in this checkout")
class TestCommittedArtifact:
    @pytest.fixture(scope="class")
    def artifact(self) -> dict[str, object]:
        return json.loads(ARTIFACT.read_text())

    def test_every_invariant_claim_passes(self, artifact):
        failed = [name for name, claim in artifact["claims"].items() if not claim["passed"]]
        assert failed == []

    def test_default_grid_numbers_match_a_fresh_computation(self, artifact):
        default = _errors(41, 13)
        stored = artifact["default_grid"]["relative_error"]
        assert default["volume"] == pytest.approx(stored["displacement_volume"], rel=1e-9)
        assert default["waterplane_area"] == pytest.approx(stored["waterplane_area"], rel=1e-9)
        assert default["vcb"] == pytest.approx(stored["vertical_centre_buoyancy"], rel=1e-9)

    def test_per_axis_orders_are_pinned(self, artifact):
        """Pinned full-ladder orders; the study's headline numbers must not drift."""
        assert artifact["predicted_orders"] == {"axial": 1.7, "vertical": 1.5}
        axial = artifact["axial_orders"]
        vertical = artifact["vertical_orders"]
        assert 1.75 < axial["displacement_volume"]["order"] < 1.90
        assert 1.75 < axial["waterplane_area"]["order"] < 1.90
        assert 1.55 < vertical["displacement_volume"]["order"] < 1.72
        assert 1.50 < vertical["vertical_centre_buoyancy"]["order"] < 1.68
        assert vertical["waterplane_area"]["order"] is None      # z-independent
        assert axial["vertical_centre_buoyancy"]["order"] is None  # x-independent

    def test_falsified_hypothesis_is_recorded_as_falsified(self, artifact):
        """The nz-heavy budget hypothesis did not survive; it must not silently pass."""
        hypothesis = artifact["hypotheses"]["nz_heavy_split_is_budget_optimal"]
        assert hypothesis["passed"] is False
        assert artifact["claims"]["equal_cell_split_beats_default"]["passed"] is True

    def test_reference_volume_is_the_closed_form(self, artifact):
        reference = artifact["reference_grid"]["reference_volume"]
        analytic = analytic_hydrostatics(HullParameters(nx=321, nz=81)).displacement_volume
        assert reference == pytest.approx(analytic, rel=1e-12)