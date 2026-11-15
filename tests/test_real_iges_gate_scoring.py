"""The gate's scoring rule: which published numbers a bbox measurement may be judged on.

Regression for a real defect: ``run_gate.py`` compared the mesh's bounding-box x-extent
against the published ``Lpp`` and scored the difference as a precision error.  On the KCS
specimen the shell spans x = -6.000...237.844 m (bulbous bow and stern overhang outside
both perpendiculars) against Lpp 230 m, so a correctly scaled, correctly imported file was
reported as a 6.02 % failure and the whole gate returned "not within tolerance".
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "real_iges_gate", ROOT / "tools" / "real_iges_gate" / "run_gate.py")
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)

KCS_SPEC = {"length_pp": 230.0, "beam": 32.2, "published_volume_m3": 52030.0}
TOLERANCE = 3.0


def test_overhanging_bounding_box_is_not_a_precision_error():
    """bbox longer than Lpp is normal; it must not be scored, and must not block the gate."""
    measured = {"length_bbox": 243.844, "beam_bbox": 32.2234}
    scoring = gate.score_against_published(measured, KCS_SPEC, TOLERANCE)
    assert "length_bbox" not in scoring["relative_errors_pct"]
    assert scoring["within_tolerance"]["beam_bbox"] <= TOLERANCE
    assert scoring["containment"] == {"length_bbox_at_least_lpp": True}
    assert scoring["length_bbox_vs_lpp_pct"] == pytest.approx(6.019, abs=1e-2)
    assert all(scoring["within_tolerance"].values())
    assert all(scoring["containment"].values())


def test_model_shorter_than_its_own_lpp_fails_containment():
    """A millimetre/metre mix-up shrinks the model; the one-sided length test catches it."""
    measured = {"length_bbox": 0.2438, "beam_bbox": 32.2234}
    scoring = gate.score_against_published(measured, KCS_SPEC, TOLERANCE)
    assert scoring["containment"] == {"length_bbox_at_least_lpp": False}


def test_beam_error_is_still_scored():
    """The commensurable quantity keeps its precision test."""
    measured = {"length_bbox": 243.844, "beam_bbox": 40.0}
    scoring = gate.score_against_published(measured, KCS_SPEC, TOLERANCE)
    assert scoring["within_tolerance"]["beam_bbox"] is False
    assert scoring["relative_errors_pct"]["beam_bbox"] == pytest.approx(24.22, abs=0.01)


def test_displacement_and_cb_are_scored_when_present():
    measured = {"length_bbox": 243.844, "beam_bbox": 32.2,
                "volume_relative_error_pct": -2.0, "cb_relative_error_pct": 1.0}
    scoring = gate.score_against_published(measured, KCS_SPEC, TOLERANCE)
    assert scoring["relative_errors_pct"]["displacement"] == pytest.approx(2.0)
    assert all(scoring["within_tolerance"].values())


def test_withheld_displacement_is_absent_not_zero():
    """A specimen with no closable volume must not be credited with a 0 % displacement error."""
    scoring = gate.score_against_published({"length_bbox": 243.844, "beam_bbox": 32.2},
                                           KCS_SPEC, TOLERANCE)
    assert "displacement" not in scoring["relative_errors_pct"]
