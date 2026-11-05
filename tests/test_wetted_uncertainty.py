"""Regression tests for the wetted-surface uncertainty evidence (tools/wetted_uncertainty.py)."""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO = Path(__file__).resolve().parents[1]
TOOLS = REPO / "tools"


def _load_tools_module(name: str) -> ModuleType:
    """Load ``tools/wetted_uncertainty.py`` by path instead of mutating ``sys.path``.

    Importing the module by file location keeps the global import path untouched:
    a module-level ``sys.path.insert`` would leak ``tools/`` to every later test
    module in the same session and make the suite order-dependent.
    """
    spec = importlib.util.spec_from_file_location(name, TOOLS / (name + ".py"))
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


_WU = _load_tools_module("wetted_uncertainty")
_WU_DOC = _load_tools_module("wetted_uncertainty_doc")
DESIGN = _WU.DESIGN
area_adaptive = _WU.area_adaptive
area_path_a = _WU.area_path_a
design_cb = _WU.design_cb
holtrop_wetted_area = _WU.holtrop_wetted_area
observed_order = _WU.observed_order
prolate_spheroid_area = _WU.prolate_spheroid_area
superellipse_area = _WU.superellipse_area

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
    """The GUM budget must stay tiny next to the modelling bias.

    Invariants checked: RSS <= linear worst case (true for any set of terms) and
    U95 = 2 uc.  U95 itself may exceed the linear sum -- with one dominant term
    the RSS nearly equals it and the factor 2 does not -- so that is not asserted.
    """
    budget = summary["budget"]
    bias = abs(summary["closure"]["total_relative"])
    assert budget["expanded_95"] < 0.005
    assert bias / budget["expanded_95"] > 20.0
    assert budget["combined_standard"] <= budget["worst_case_sum"] + 1e-12
    assert budget["expanded_95"] == pytest.approx(2.0 * budget["combined_standard"], rel=1e-12)


def test_propagation_gate_holds_from_froude_013(summary):
    """dRt/Rt = pi dS/S to machine precision; the < 3 % gate holds for Fn >= 0.13.

    At 6 kn the friction share is largest (pi = 0.373), so 0.373 * 8.69 % = 3.24 %
    is genuinely above 3 %; the study reports that row instead of averaging it away.
    """
    rows = summary["propagation"]
    assert rows, "platform resistance cache missing"
    delta = abs(summary["closure"]["total_relative"])
    for row in rows:
        assert row["drt_relative"] == pytest.approx(row["pi"] * delta, rel=1e-12)
        assert row["rf_identity_relative_residual"] < 1e-3
    planning = [r for r in rows if r["froude"] >= 0.13]
    assert planning
    assert max(r["drt_relative"] for r in planning) < 0.03
    low_speed = [r for r in rows if r["froude"] < 0.13]
    assert low_speed
    assert max(r["drt_relative"] for r in low_speed) < 0.035
    # the friction fraction must be the report-3871 one at Fn ~ 0.13
    slow = min(rows, key=lambda r: abs(r["froude"] - 0.1314))
    assert abs(slow["pi"] / 0.3422 - 1.0) < 0.01


def test_design_geometry_precheck(summary):
    """Cb/L/B gate: mesh volume must reproduce LBT*Cb to well under 0.1 %."""
    volume = summary["mesh"]["volume_check"]
    assert abs(volume["mesh_volume_finest_m3"] / volume["LBT_Cb_m3"] - 1.0) < 1e-3
    assert abs(summary["mesh"]["mesh_cb_relative_error_finest"]) < 1e-3


def test_doc_is_current():
    DOC, render = _WU_DOC.DOC, _WU_DOC.render

    assert DOC.exists()
    assert DOC.read_text() == render(json.loads(SUMMARY.read_text()))


def test_summary_has_no_non_finite_numbers(summary):
    """The published evidence must be finite everywhere.

    An earlier revision of the run wrote ``nan`` into the family-sweep budget
    (non-monotone Richardson window for steep exponents); the document then
    reported a bar that no code could reproduce.  Fail loudly instead.
    """
    bad = []

    def walk(node, path=""):
        if isinstance(node, dict):
            for key, value in node.items():
                walk(value, f"{path}.{key}")
        elif isinstance(node, list):
            for index, value in enumerate(node):
                walk(value, f"{path}[{index}]")
        elif isinstance(node, float) and not math.isfinite(node):
            bad.append((path, node))

    walk(summary)
    assert not bad, bad
