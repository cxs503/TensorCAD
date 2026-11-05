"""The published ``docs/gz`` evidence bundle must stay a view of a real run.

``docs/gz`` holds the righting-arm evidence: the example's ``results.json``, its figure, and
``tables.md`` rendered from the JSON by ``examples/gz/publish_evidence.py``. Two failure modes
are worth a test of their own:

* the bundle drifts from the code (the numbers were measured once and the module changed);
* a number in the prose is transcribed by hand and stops matching the measurement.

The first is covered by re-deriving stored values from the library and from the closed forms;
the second by re-rendering ``tables.md`` and demanding byte equality.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BUNDLE = ROOT / "docs" / "gz"
sys.path.insert(0, str(ROOT / "examples" / "gz"))

import publish_evidence  # noqa: E402  (path is set up above)
import run_example  # noqa: E402


@pytest.fixture(scope="module")
def results() -> dict:
    return json.loads((BUNDLE / "results.json").read_text())


def test_bundle_is_present() -> None:
    for name in ("results.json", "tables.md", "plot.png"):
        path = BUNDLE / name
        assert path.is_file(), f"missing evidence file {name}; run examples/gz/publish_evidence.py"
        assert path.stat().st_size > 0


def test_plot_is_a_real_png() -> None:
    """A committed figure that is empty or truncated would be worse than none."""
    payload = (BUNDLE / "plot.png").read_bytes()
    assert payload[:8] == b"\x89PNG\r\n\x1a\n"
    assert payload[12:16] == b"IHDR"
    width = int.from_bytes(payload[16:20], "big")
    height = int.from_bytes(payload[20:24], "big")
    assert width >= 1000 and height >= 400, (width, height)


def test_tables_are_rendered_from_the_json(results: dict) -> None:
    """No hand-typed numbers: the committed page must be reproducible byte for byte."""
    assert publish_evidence.render_tables(results) == (BUNDLE / "tables.md").read_text()


def test_claim_ledger_holds(results: dict) -> None:
    checks = results["checks"]
    assert checks["failures"] == []
    assert checks["passed"] is True
    assert results["box_wall_sided"]["worst_abs_error_m"] <= 1e-12
    assert results["log_circular_section"]["worst_abs_error_m"] <= 1e-9
    hull = results["hull"]
    slope, reference = hull["gm_from_curve_slope_m"], hull["gm_from_initial_stability_m"]
    assert abs(slope - reference) / reference < 1e-6
    assert hull["curve"]["max_volume_error"] < 1e-10
    # negative control: the wall-sided formula has to be *wrong* past its limit
    assert results["box_wall_sided"]["beyond_limit"][-1]["abs_error"] > 1e-2


def test_stored_box_barge_row_recomputes(results: dict) -> None:
    """Re-heel the barge from scratch and compare with the stored row (code, not memory)."""
    from tensorcad.gz import righting_arm_curve

    barge = run_example.box(
        (0.0, -run_example.BREADTH / 2, 0.0),
        (run_example.LENGTH, run_example.BREADTH / 2, run_example.DEPTH),
    )
    row = next(row for row in results["box_wall_sided"]["inside"] if row["heel_deg"] == 20.0)
    recomputed = righting_arm_curve(barge, run_example.KG, [20.0], waterline_z=run_example.DRAFT)
    assert recomputed.gz[0] == pytest.approx(row["gz_m"], abs=1e-12)


def test_stored_log_rows_match_the_exact_form(results: dict) -> None:
    """The log reference is analytic, so the stored rows can be checked without the library."""
    section = results["log_circular_section"]
    radius, kg = section["radius_m"], results["hull"]["parameters"]["kg_m"]
    for row in section["rows"]:
        exact = (radius - kg) * math.sin(math.radians(row["heel_deg"]))
        assert row["exact_m"] == pytest.approx(exact, abs=1e-12)
        assert row["gz_m"] == pytest.approx(exact, abs=1e-9)


def test_every_curve_point_floats_at_one_displacement(results: dict) -> None:
    curve = results["hull"]["curve"]
    target = curve["displacement_volume_m3"]
    assert len(curve["points"]) == len(curve["heel_angles_deg"]) == len(curve["gz_m"])
    for point in curve["points"]:
        assert point["target_volume_m3"] == pytest.approx(target, rel=1e-9)
        assert abs(point["volume_error"]) < 1e-10 or point["volume_error"] == 0.0


def test_doc_headlines_match_the_json(results: dict) -> None:
    """The prose in docs/gz.md cites four numbers; they must be the measured ones.

    The page is written by hand (it is documentation), so this is the guard that keeps its
    headline figures from drifting: each token below is *computed* from ``results.json``.
    """
    page = (ROOT / "docs" / "gz.md").read_text()
    box = results["box_wall_sided"]
    hull = results["hull"]
    curve = hull["curve"]
    slope, reference = hull["gm_from_curve_slope_m"], hull["gm_from_initial_stability_m"]
    tokens = [
        f"{box['worst_abs_error_m']:.1e}",
        f"{results['log_circular_section']['worst_abs_error_m']:.1e}",
        f"{abs(slope - reference) / reference:.1e}",
        f"{curve['max_volume_error']:.1e}",
        f"{box['beyond_limit'][-1]['abs_error']:.2f}",
        f"{len(curve['points'])} heel stations",
        f"{box['limit_deg']:.2f} deg",
    ]
    for token in tokens:
        assert token in page, f"docs/gz.md no longer cites {token!r} from results.json"
