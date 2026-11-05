"""The published GZ bundle is *derived*, not transcribed.

`docs/gz/report.md` quotes deviations of the righting-arm integrator against two closed forms.
These tests rebuild that evidence from the library (`publish_report.build_evidence`) and compare
it with the committed bundle, so the documentation cannot silently drift away from the code: a
change that moves a righting arm of the box barge, the circular section or the hull fails here
before it can ship as a plausible report.
"""

from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path
from types import ModuleType

import pytest

import tensorcad

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs" / "gz"
CLAIM_KEYS = {
    "box_wall_sided",
    "box_negative_control",
    "round_bilge_exact",
    "hull_gm_consistency",
    "displacement_residual",
    "free_surface_flat_shift",
}


def _load_publisher() -> ModuleType:
    spec = importlib.util.spec_from_file_location("tensorcad_publish_report", ROOT / "publish_report.py")
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def publisher() -> ModuleType:
    return _load_publisher()


@pytest.fixture(scope="module")
def evidence(publisher: ModuleType) -> dict[str, object]:
    return publisher.build_evidence()


@pytest.fixture(scope="module")
def claims(evidence: dict[str, object]) -> dict[str, dict[str, object]]:
    return {item["key"]: item for item in evidence["claims"]}  # type: ignore[union-attr]


def test_committed_bundle_matches_the_code(publisher: ModuleType, evidence: dict[str, object]) -> None:
    """``--check`` mode: the committed results.json is what the library produces today."""
    assert publisher.check(evidence) == []


def test_every_claim_holds(evidence: dict[str, object], claims: dict[str, dict[str, object]]) -> None:
    assert set(claims) == CLAIM_KEYS
    assert evidence["verdict"] == {"passed": True, "failed": []}  # type: ignore[index]
    for key, claim in claims.items():
        assert claim["passed"], f"{key}: {claim['measured']} {claim['unit']}"
        if claim["two_sided"]:
            assert abs(claim["measured"]) <= claim["tolerance_abs"]


def test_the_two_exact_benchmarks_are_at_machine_precision(
    claims: dict[str, dict[str, object]],
) -> None:
    box = claims["box_wall_sided"]
    log = claims["round_bilge_exact"]
    assert box["measured"] < 1e-14, box
    assert log["measured"] < 1e-14, log
    # ...and relative to the arm itself, not just in metres
    assert box["relative_percent"] < 1e-12
    assert log["relative_percent"] < 1e-12


def test_wall_sided_formula_has_a_domain_and_the_bundle_says_so(
    evidence: dict[str, object],
) -> None:
    """The negative control: past the deck-edge limit the closed form must stop matching."""
    box_case = evidence["cases"]["box_barge"]  # type: ignore[index]
    assert 30.0 < box_case["wall_sided_limit_deg"] < 45.0
    assert box_case["beyond_limit"][-1]["abs_error_m"] > 1e-2
    assert box_case["beyond_limit"][-1]["heel_deg"] > box_case["wall_sided_limit_deg"]


def test_hull_curve_is_consistent_with_the_waterplane_route(
    evidence: dict[str, object],
) -> None:
    hull = evidence["cases"]["hull"]  # type: ignore[index]
    assert hull["gm_abs_error_m"] < 1e-6
    assert math.isclose(
        hull["gm_from_curve_slope_m"], hull["gm_from_initial_stability_m"], rel_tol=1e-8
    )
    assert hull["max_gz_m"] > 0.0
    assert hull["max_volume_error"] < 1e-10


def test_free_surface_correction_is_the_flat_shift_it_claims(
    evidence: dict[str, object],
) -> None:
    free = evidence["cases"]["hull"]["free_surface"]  # type: ignore[index]
    assert free["max_shift_error_m"] < 1e-12
    assert math.isclose(
        free["max_gz_solid_m"] - free["max_gz_with_free_surface_m"], free["loss_m"], rel_tol=1e-12
    )


def test_results_json_is_the_documented_schema(evidence: dict[str, object]) -> None:
    committed = json.loads((DOCS / "results.json").read_text())
    assert committed["schema"] == "tensorcad/gz-evidence@1"
    assert committed["tensorcad_version"] == tensorcad.__version__
    assert committed["verdict"]["passed"] is True
    for path in (committed["figure"], committed["example"]):
        assert (ROOT / path).is_file(), path
    assert {item["key"] for item in committed["claims"]} == CLAIM_KEYS
    assert evidence["schema"] == committed["schema"]


def test_report_has_the_six_sections_and_only_injected_numbers(
    evidence: dict[str, object], claims: dict[str, dict[str, object]]
) -> None:
    report = (DOCS / "report.md").read_text()
    for section in ("## 1.", "## 2.", "## 3.", "## 4.", "## 5.", "## 6."):
        assert section in report, section
    assert "![GZ curves against the closed forms](gz-curve.png)" in report
    assert f"all {len(CLAIM_KEYS)} claims passed" in report
    for key, claim in claims.items():
        assert f"`{key}`" in report, key
        assert f"{claim['measured']:.2e}" in report, key
    # the hull's two independent GM routes are both quoted, at full precision
    hull = evidence["cases"]["hull"]  # type: ignore[index]
    assert f"{hull['gm_from_initial_stability_m']:.9f}" in report
    assert f"{hull['gm_from_curve_slope_m']:.9f}" in report


def test_report_is_reproducible_from_the_library(publisher: ModuleType, evidence: dict[str, object]) -> None:
    assert (DOCS / "report.md").read_text() == publisher.render_markdown(evidence)


def test_figure_is_a_small_png() -> None:
    payload = (DOCS / "figures" / "gz-curve.png").read_bytes()
    assert payload[:8] == b"\x89PNG\r\n\x1a\n"
    assert 10_000 < len(payload) < 200_000, len(payload)  # embedded in an API report; stays small


def test_the_bundle_agrees_with_the_runnable_example(
    publisher: ModuleType, evidence: dict[str, object]
) -> None:
    """``examples/gz/outputs/results.json`` is committed too: the two derivations must agree."""
    example = json.loads((ROOT / "examples" / "gz" / "outputs" / "results.json").read_text())
    cases = evidence["cases"]
    assert math.isclose(
        example["box_wall_sided"]["worst_abs_error_m"],
        cases["box_barge"]["worst_abs_error_m"],  # type: ignore[index]
        rel_tol=1e-12,
    )
    assert math.isclose(
        example["log_circular_section"]["worst_abs_error_m"],
        cases["round_bilge_log"]["worst_abs_error_m"],  # type: ignore[index]
        rel_tol=1e-12,
    )
    assert math.isclose(
        example["hull"]["gm_from_curve_slope_m"],
        cases["hull"]["gm_from_curve_slope_m"],  # type: ignore[index]
        rel_tol=1e-12,
    )
    assert example["checks"]["passed"] is True


def test_check_detects_drift(publisher: ModuleType, evidence: dict[str, object]) -> None:
    """A tampered committed JSON must be reported, so the test above is not vacuous."""
    committed = json.loads((DOCS / "results.json").read_text())
    tampered = json.loads(json.dumps(committed))
    tampered["cases"]["box_barge"]["worst_abs_error_m"] = 1.0e-3
    mismatches: list[str] = []
    publisher._compare(evidence, tampered, "", mismatches)
    assert mismatches, "a moved righting arm must produce a mismatch"
    assert publisher.check(evidence) == []


def test_cli_check_mode_exits_zero(publisher: ModuleType) -> None:
    assert publisher.main(["--check"]) == 0
