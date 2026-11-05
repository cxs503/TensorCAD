"""`publish_report.py` must build its payload from the artefacts, not from prose.

The report is what a reviewer reads, so the failure mode worth a test is not "the Markdown
is pretty" but "the numbers in it came from somewhere else". Here the payload is rebuilt
from the committed `docs/gz/results.json` and checked for: the six required sections, every
IMO criterion and every claim present, the mission id, a real PNG inside the figure bundle,
and the endpoint's size ceiling.
"""
from __future__ import annotations

import base64
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import publish_report  # noqa: E402

RESULTS_PATH = ROOT / "docs" / "gz" / "results.json"
FIGURE_PATH = ROOT / "docs" / "gz" / "figures" / "gz-curve.png"

pytestmark = pytest.mark.skipif(
    not (RESULTS_PATH.is_file() and FIGURE_PATH.is_file()),
    reason="the example artefacts (docs/gz/results.json, figures/gz-curve.png) are not built",
)


@pytest.fixture(scope="module")
def payload() -> dict:
    results = json.loads(RESULTS_PATH.read_text())
    return publish_report.build_payload(
        results, publish_report.geometry_checks(), mission_id=4342, figure=FIGURE_PATH
    )


def test_payload_has_the_six_sections(payload: dict) -> None:
    markdown = payload["markdown"]
    for heading in ("## 1. Case", "## 2. Method", "## 3. Results", "## 4. Figure",
                    "## 5. Notes", "## 6. Error judgement"):
        assert heading in markdown, heading


def test_payload_carries_the_measured_numbers(payload: dict) -> None:
    results = json.loads(RESULTS_PATH.read_text())
    markdown = payload["markdown"]
    for check in results["barge"]["imo"]["checks"]:
        assert f"`{check['key']}`" in markdown
        assert f"{check['measured']:.6f}" in markdown
    for row in results["barge"]["table"]:
        assert f"| {row['heel_deg']:g} |" in markdown
    assert results["barge"]["imo_negative_control"]["passed"] is False
    assert "`passed = False`" in markdown


def test_payload_identifies_the_mission_and_fits_the_endpoint(payload: dict) -> None:
    assert payload["mission_id"] == 4342
    assert payload["scenario_key"] == publish_report.SCENARIO_KEY
    size = len(json.dumps(payload))
    assert size <= publish_report.PAYLOAD_LIMIT
    encoded = payload["figures"][0]["png_base64"]
    assert base64.b64decode(encoded)[:8] == b"\x89PNG\r\n\x1a\n"


def test_geometry_checks_recompute_the_barge_volume() -> None:
    barge = publish_report.geometry_checks()["barge"]
    assert barge["closed_box_volume_m3"] == pytest.approx(1800.0)
    assert barge["mesh_volume_m3"] == pytest.approx(1800.0, rel=1e-12)
    assert barge["watertight"] is True
    assert barge["block_coefficient"] == pytest.approx(1.0)
