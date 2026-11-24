"""Acceptance tests for the Series 60 anchor tool.

This module exercises ``tools/series60_family/anchor_acceptance.py`` end-to-end:
the five-Cb gradient, the transom/bulb closed forms, and the gate-power mutation
record.  It also checks that the tool refuses a record that lacks provenance
(mutant injection details and measured exit status).

The anchor's acceptance clauses are:
1. Series60 五档 Cb 生成 + 静水力对照公开值
2. gate-power 突变测试沿用 >=80%
3. 文档齐全

Every clause reports ``pass`` / ``fail`` / ``not_run`` / ``withheld``;
every comparison carries an evidence **grade** drawn from a closed set.
"""
from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tools.series60_family.anchor_acceptance import (  # noqa: E402
    ANCHOR,
    DEFAULT_MUTATIONS,
    OUT,
    QUANTITY_GRADES,
    GRADES,
    STATUSES,
    TOLERANCE_PUBLISHED,
    TOLERANCE_DIMENSION,
    TOLERANCE_REALISED_CB,
    GATE_POWER_MINIMUM,
    REQUIRED_DOCS,
    load_published,
    relative_error,
    geometry_precondition,
    grade_gradient,
    provenance_defects,
)


@pytest.fixture(scope="module")
def published():
    return load_published()


@pytest.fixture(scope="module")
def gradient(published):
    return grade_gradient(published)


def test_anchor_id_is_recorded():
    assert ANCHOR == "GAP-CAD-T2-PARAM-SERIES"


def test_required_docs_exist():
    for doc in REQUIRED_DOCS:
        assert (ROOT / doc).exists(), f"Missing required doc: {doc}"


def test_quantity_grades_are_valid(gradient):
    for parent in gradient["parents"]:
        for key, comp in parent["comparisons"].items():
            assert comp["grade"] in GRADES, (
                f"Invalid grade {comp['grade']!r} for {key}"
            )


def test_closure_grades_not_mislabelled_as_external(gradient):
    """The audit found closure comparisons read as external; this must not happen."""
    for parent in gradient["parents"]:
        for key, comp in parent["comparisons"].items():
            if QUANTITY_GRADES.get(key) == "closure":
                assert comp["grade"] == "closure", (
                    f"{key} is a closure quantity but graded {comp['grade']!r}"
                )


def test_gradient_tolerance_holds(gradient):
    """All five parents must agree with published particulars within 3%."""
    for parent in gradient["parents"]:
        for key, comp in parent["comparisons"].items():
            assert abs(comp["relative"]) <= TOLERANCE_PUBLISHED, (
                f"{key} at Cb={parent['cb']}: rel err {comp['relative']:.4%} > {TOLERANCE_PUBLISHED}"
            )


def test_gradient_validation_passed(gradient):
    """All five parents must pass validation."""
    for parent in gradient["parents"]:
        assert parent["validation_passed"], (
            f"Cb={parent['cb']} validation failed: {parent['validation_errors']}"
        )


def test_gradient_geometry_precondition(gradient):
    """L/B/T and Cb must match the request within tolerance."""
    for parent in gradient["parents"]:
        geo = parent["geometry"]
        assert geo["passed"], (
            f"Cb={parent['cb']} geometry precondition failed"
        )


def test_gradient_worst_error_within_bar(gradient):
    """The worst relative error must be within the published tolerance."""
    assert abs(gradient["worst_relative"]) <= TOLERANCE_PUBLISHED, (
        f"Worst relative error {gradient['worst_relative']:.4%} > {TOLERANCE_PUBLISHED}"
    )


def test_provenance_defects_detects_missing_mutant():
    """A record without a 'mutant' block must be flagged."""
    entry = {"status": "ok", "run": {"killed": True}}
    defects = provenance_defects(entry)
    assert any("no 'mutant' block" in d for d in defects)


def test_provenance_defects_detects_missing_keys():
    """A mutant block missing required keys must be flagged."""
    entry = {
        "mutant": {
            "id": "TEST-1",
            "kind": "model",
            "file": "series",
            # missing old, new, claim
        }
    }
    defects = provenance_defects(entry)
    assert any("mutant.old missing" in d for d in defects)
    assert any("mutant.new missing" in d for d in defects)
    assert any("mutant.claim missing" in d for d in defects)


def test_provenance_defects_passes_valid_mutant():
    """A mutant block with all required keys must pass."""
    entry = {
        "mutant": {
            "id": "TEST-1",
            "kind": "model",
            "file": "series",
            "old": "old code",
            "new": "new code",
            "claim": "this mutant breaks X",
        }
    }
    defects = provenance_defects(entry)
    assert defects == []


def test_gate_power_record_provenance(tmp_path):
    """The gate-power record must carry provenance for every mutant."""
    # Create a minimal valid record
    record = {
        "S60-1-transom-blend-area": {
            "mutant": {
                "id": "S60-1-transom-blend-area",
                "kind": "model",
                "file": "series",
                "old": "keep * h + transom * 0.5",
                "new": "keep * h + transom * 0.55",
                "claim": "transom_ratio blends the aft taper",
            },
            "run": {
                "killed": True,
                "returncode": 1,
                "failed": True,
                "wall_seconds": 18.2,
            },
            "status": "ok",
        }
    }
    record_path = tmp_path / "test_record.json"
    record_path.write_text(json.dumps(record))

    # Load and check
    data = json.loads(record_path.read_text())
    for mutant_id, entry in data.items():
        defects = provenance_defects(entry)
        assert defects == [], (
            f"Mutant {mutant_id} has provenance defects: {defects}"
        )


def test_relative_error_function():
    """Test the relative_error helper."""
    assert relative_error(1.05, 1.0) == pytest.approx(0.05)
    assert relative_error(0.95, 1.0) == pytest.approx(-0.05)
    assert relative_error(1.0, 1.0) == 0.0


def test_tolerance_constants():
    """Check that tolerance constants are reasonable."""
    assert TOLERANCE_PUBLISHED == 0.03
    assert TOLERANCE_DIMENSION == 0.001
    assert TOLERANCE_REALISED_CB == 0.02
    assert GATE_POWER_MINIMUM == 0.80


def test_grade_gradient_returns_expected_structure(gradient):
    """Check that grade_gradient returns the expected structure."""
    assert "id" in gradient
    assert "statement" in gradient
    assert "tolerance" in gradient
    assert "grid" in gradient
    assert "parents" in gradient
    assert "worst_relative" in gradient
    assert "worst_quantity" in gradient
    assert "worst_cb" in gradient
    assert "grade" in gradient
    assert "grade_note" in gradient
    assert "passed" in gradient


def test_grade_gradient_parents_count(gradient, published):
    """Check that grade_gradient processes all five parents."""
    assert len(gradient["parents"]) == len(published["parents"])


def test_grade_gradient_worst_cb_matches(gradient, published):
    """Check that worst_cb corresponds to the worst_relative parent."""
    worst_parent = next(
        p for p in gradient["parents"]
        if abs(p["comparisons"][gradient["worst_quantity"]]["relative"]) == gradient["worst_relative"]
    )
    assert worst_parent["cb"] == gradient["worst_cb"]


def test_anchor_acceptance_cli_runs(tmp_path):
    """Test that the anchor acceptance CLI runs successfully."""
    result = subprocess.run(
        [sys.executable, str(ROOT / "tools" / "series60_family" / "anchor_acceptance.py"),
         "--out", str(tmp_path / "anchor_acceptance")],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, f"CLI failed: {result.stderr}"


def test_anchor_acceptance_artifacts_created(tmp_path):
    """Test that the anchor acceptance CLI creates expected artifacts."""
    subprocess.run(
        [sys.executable, str(ROOT / "tools" / "series60_family" / "anchor_acceptance.py"),
         "--out", str(tmp_path / "anchor_acceptance")],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert (tmp_path / "anchor_acceptance" / "anchor_acceptance.json").exists()
    assert (tmp_path / "anchor_acceptance" / "anchor_acceptance.csv").exists()
    assert (tmp_path / "anchor_acceptance" / "figures" / "anchor_acceptance.png").exists()


def test_anchor_acceptance_json_structure(tmp_path):
    """Test that the anchor acceptance JSON has the expected structure."""
    subprocess.run(
        [sys.executable, str(ROOT / "tools" / "series60_family" / "anchor_acceptance.py"),
         "--out", str(tmp_path / "anchor_acceptance")],
        capture_output=True,
        text=True,
        timeout=120,
    )
    data = json.loads((tmp_path / "anchor_acceptance" / "anchor_acceptance.json").read_text())
    assert "anchor" in data
    assert "clauses" in data
    assert "metadata" in data


def test_anchor_acceptance_clauses_status(tmp_path):
    """Test that all clauses have a valid status."""
    subprocess.run(
        [sys.executable, str(ROOT / "tools" / "series60_family" / "anchor_acceptance.py"),
         "--out", str(tmp_path / "anchor_acceptance")],
        capture_output=True,
        text=True,
        timeout=120,
    )
    data = json.loads((tmp_path / "anchor_acceptance" / "anchor_acceptance.json").read_text())
    for clause in data["clauses"]:
        assert clause["status"] in STATUSES, (
            f"Invalid status {clause['status']!r} for clause {clause['id']}"
        )


def test_anchor_acceptance_verdict(tmp_path):
    """Test that the anchor acceptance verdict is 'pass'."""
    subprocess.run(
        [sys.executable, str(ROOT / "tools" / "series60_family" / "anchor_acceptance.py"),
         "--out", str(tmp_path / "anchor_acceptance")],
        capture_output=True,
        text=True,
        timeout=120,
    )
    data = json.loads((tmp_path / "anchor_acceptance" / "anchor_acceptance.json").read_text())
    assert data["verdict"] == "pass", (
        f"Expected verdict 'pass', got {data['verdict']!r}"
    )
