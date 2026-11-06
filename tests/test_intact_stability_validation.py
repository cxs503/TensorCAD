"""Tests for ``tools/intact_stability_validation.py`` and the evidence it commits.

The harness is the evidence behind section 6 of ``docs/intact_stability.md``: it checks the
numbers ``tensorcad.intact_stability``'s weather criterion compares against closed forms.
These tests do three things the script cannot do about itself --

* run the harness end to end (``--no-figure``) in a temporary directory and require every
  claim to hold and the exit code to be zero;
* compare the committed ``docs/intact_stability/intact_stability_validation.json`` against a
  fresh run, so a stale artifact cannot sit next to changed code; the committed figure is
  required to exist as well;
* break one claim on purpose and require a non-zero exit, so a harness that cannot fail is
  not mistaken for a passing one.
"""
from __future__ import annotations

import importlib.util
import json
import pathlib
import sys

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
TOOL = REPO_ROOT / "tools" / "intact_stability_validation.py"
COMMITTED = REPO_ROOT / "docs" / "intact_stability" / "intact_stability_validation.json"
FIGURE = REPO_ROOT / "docs" / "intact_stability" / "intact_stability_validation.png"


def _load_tool():
    spec = importlib.util.spec_from_file_location("intact_stability_validation", TOOL)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def harness():
    return _load_tool()


@pytest.fixture(scope="module")
def fresh_claims(harness) -> list[dict[str, object]]:
    """One evaluation reused by every test in the module (the run is the expensive part)."""
    return harness.evaluate()


def test_every_claim_holds(fresh_claims: list[dict[str, object]]) -> None:
    failed = [claim["name"] for claim in fresh_claims if not claim["holds"]]
    assert failed == [], f"claims that no longer hold: {failed}"


def test_claims_are_distinct_and_measured(fresh_claims: list[dict[str, object]]) -> None:
    names = [claim["name"] for claim in fresh_claims]
    assert len(names) == len(set(names)) >= 10
    for claim in fresh_claims:
        assert claim["claim"], claim["name"]
        assert claim["measured"], claim["name"]
        assert claim["expected"], claim["name"]
        assert float(claim["tolerance"]) >= 0.0


def test_main_writes_artifacts_and_exits_zero(harness, tmp_path: pathlib.Path) -> None:
    assert harness.main(["--output", str(tmp_path), "--no-figure"]) == 0
    payload = json.loads((tmp_path / "intact_stability_validation.json").read_text())
    assert payload["schema"] == harness.SCHEMA
    assert payload["all_claims_hold"] is True
    assert payload["failures"] == []
    assert len(payload["claims"]) == len(harness.CHECKS)
    reference = payload["reference"]
    assert reference["gm0_m"] == pytest.approx(harness.GM0)
    assert reference["bm0_m"] == pytest.approx(harness.BM0)
    assert reference["theta2_cap_deg"] == 50.0


def test_main_returns_one_when_a_claim_breaks(harness, tmp_path: pathlib.Path, monkeypatch) -> None:
    def broken() -> dict[str, object]:
        return {
            "claim": "deliberately broken claim",
            "holds": False,
            "measured": "nothing",
            "expected": "something",
            "error": 1.0,
            "tolerance": 0.0,
        }

    monkeypatch.setattr(harness, "CHECKS", (("broken", broken),))
    assert harness.main(["--output", str(tmp_path), "--no-figure"]) == 1
    payload = json.loads((tmp_path / "intact_stability_validation.json").read_text())
    assert payload["failures"] == ["broken"]
    assert payload["all_claims_hold"] is False


def test_committed_evidence_is_current(fresh_claims: list[dict[str, object]]) -> None:
    assert COMMITTED.is_file(), "docs/intact_stability/intact_stability_validation.json is missing"
    payload = json.loads(COMMITTED.read_text())
    assert payload["all_claims_hold"] is True
    assert payload["failures"] == []
    # runtime_s is wall-clock and is the only field allowed to move between runs.
    assert payload["claims"] == fresh_claims, (
        "the committed evidence no longer matches the harness output: "
        "re-run python tools/intact_stability_validation.py"
    )


def test_committed_figure_exists() -> None:
    assert FIGURE.is_file(), "docs/intact_stability/intact_stability_validation.png is missing"
    assert FIGURE.stat().st_size > 10_000
    assert FIGURE.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def test_harness_uses_no_unverified_constants(fresh_claims: list[dict[str, object]]) -> None:
    """Every claim must carry a number and a tolerance it is compared against."""
    by_name = {claim["name"]: claim for claim in fresh_claims}
    assert by_name["gust_area_a"]["error"] < 1e-12
    assert by_name["second_intercept"]["error"] < 1e-2
    assert by_name["grid_convergence"]["error"] < 0.5
    ratios = by_name["grid_convergence"]["ratios"]
    assert len(ratios) == 3 and all(3.5 <= ratio <= 4.5 for ratio in ratios[-2:])

def test_docs_section_is_generated_not_transcribed(fresh_claims, tmp_path: pathlib.Path) -> None:
    """Section 6 of docs/intact_stability.md must equal what the harness renders now."""
    doc = REPO_ROOT / "docs" / "intact_stability.md"
    text = doc.read_text(encoding="utf-8")
    harness = _load_tool()
    begin = text.find(harness.MARKDOWN_BEGIN)
    end = text.find(harness.MARKDOWN_END)
    assert begin >= 0 and end > begin, "the generated-block markers are missing from the doc"
    committed_block = text[begin : end + len(harness.MARKDOWN_END)]
    payload = json.loads(COMMITTED.read_text())
    rendered = harness.render_section(payload)
    assert committed_block.strip() == rendered.strip(), (
        "docs/intact_stability.md section 6 no longer matches the harness output: "
        "re-run python tools/intact_stability_validation.py --update-docs docs/intact_stability.md"
    )


def test_heading_numbering_of_the_page_is_consistent() -> None:
    """The generated section is numbered 6; the hand-written sections after it must follow."""
    text = (REPO_ROOT / "docs" / "intact_stability.md").read_text(encoding="utf-8")
    headings = [line.split(".", 1)[0][3:] for line in text.splitlines() if line.startswith("## ")]
    assert headings == [str(index) for index in range(1, 10)]
