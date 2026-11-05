"""Regression tests for the wetted-surface uncertainty evidence (tools/wetted_uncertainty.py)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "tools"))

from wetted_uncertainty import (  # noqa: E402
    DESIGN,
    area_adaptive,
    area_path_a,
    closure_error,
    design_cb,
    holtrop_wetted_area,
    observed_order,
    prolate_spheroid_area,
    superellipse_area,
)

SUMMARY = REPO / "data/wetted-uncertainty-2026-11-04/summary.json"


@pytest.fixture(scope="module")
def summary() -> dict:
    return json.loads(SUMMARY.read_text())


def test_three_integrators_agree():
    """Repo GL, tanh-sinh and adaptive quad must agree far below the closure error."""
    design = dict(DESIGN)
    a = area_path_a(design, 320)
    c = area_adaptive(design)
    assert abs(a / c - 1.0) < 1e-9


def test_quadrature_is_converged():
    design = dict(DESIGN)
    assert abs(area_path_a(design, 160) / area_path_a(design, 320) - 1.0) < 1e-10


def test_closure_reproduces_platform_regression(summary):
    """Local Holtrop implementation vs the platform's own wetted area."""
    worst = summary["cross_checks"]["platform_wetted_vs_local_holtrop_relative"]
    assert worst is not None and worst < 1e-4
    design = dict(DESIGN)
    cb = design_cb(design)
    value = holtrop_wetted_area(design["length_m"], design["beam_m"], design["draft_m"], cb)
    assert abs(value / summary["closure"]["S_platform_m2"] - 1.0) < 1e-6


def test_cross_check_against_report_3871(summary):
    """Two independent closed-form anchors published in tensor_reports#3871."""
    assert abs(summary["cross_checks"]["relative_difference_vs_3871_one_side"]) < 1e-6
    assert abs(summary["cross_checks"]["relative_difference_prolate_vs_3871"]) < 1e-5
    design = dict(DESIGN)
    design.update({"length_m": 100.0, "beam_m": 10.0, "draft_m": 5.0})
    assert abs(superellipse_area(design, 2.0, 2.0, nodes=200) / 2.0 / 1275.2810 - 1.0) < 1e-6
    assert abs(prolate_spheroid_area(100.0, 10.0) / 2478.7758 - 1.0) < 1e-5


def test_observed_order_sign_convention():
    """Refined-by-two sequences with relative ratios give a positive order."""
    f = [1.0 + 4e-3, 1.0 + 1.7e-3, 1.0 + 7.4e-4]
    p = observed_order(f, 2.0)
    assert p == pytest.approx(1.26, abs=0.05)


def test_gci_brackets_true_error(summary):
    mesh = summary["mesh"]
    assert mesh["gci_brackets_all_levels"] is True
    assert 1.0 < mesh["finest_p_observed"] < 1.5
    assert abs(mesh["richardson_vs_analytic_relative"]) < 1e-4


def test_budget_is_far_below_bias(summary):
    budget = summary["budget"]
    bias = abs(summary["closure"]["total_relative"])
    assert budget["expanded_95"] < 0.005
    assert bias / budget["expanded_95"] > 20.0
    assert budget["expanded_95"] <= budget["worst_case_sum"] + 1e-12


def test_propagation_under_three_percent_at_planning_speeds(summary):
    rows = summary["propagation"]
    assert rows, "platform resistance cache missing"
    for row in rows:
        assert abs(row["drt_relative"]) < 0.03
        assert row["rf_identity_relative_residual"] < 1e-3
    # the friction fraction must be the report-3871 one at Fn ~ 0.13
    slow = min(rows, key=lambda r: abs(r["froude"] - 0.1314))
    assert abs(slow["pi"] / 0.3422 - 1.0) < 0.01


def test_design_geometry_precheck(summary):
    """Cb/L/B gate: mesh volume must reproduce LBT*Cb to well under 0.1 %."""
    volume = summary["mesh"]["volume_check"]
    assert abs(volume["mesh_volume_finest_m3"] / volume["LBT_Cb_m3"] - 1.0) < 1e-3
    assert abs(summary["mesh"]["mesh_cb_relative_error_finest"]) < 1e-3


def test_doc_is_current():
    from wetted_uncertainty_doc import DOC, render

    assert DOC.exists()
    assert DOC.read_text() == render(json.loads(SUMMARY.read_text()))
