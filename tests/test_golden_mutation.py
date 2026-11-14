"""Mutation audit of the golden geometry gate (GAP-CAD-D2-GOLDEN-REGRESSION).

``tests/test_golden_regression.py`` already proves the gate *harness* fails closed: a
missing, unreadable or stale baseline, a hand-edited reference, a tree that refuses to bless
itself, the git hook and the CI workflow. Every one of those negative tests injects either a
synthetic probe case or a tampered baseline file, so none of them answers the question the
gate exists for: does it see a regression in the geometry code it is wired to?

These tests answer it, and keep the answer honest:

* :func:`test_library_mutations_are_detected` pushes small mutations into the real
  ``tensorcad`` entry points and requires the untouched gate to report them, localised to
  the mutated case;
* :func:`test_no_mutation_control_passes` requires that a gate which fails everything
  cannot masquerade as a sensitive one;
* the envelope tests require the *declared* tolerance to be the threshold that actually
  decides, and require that a coherent mutation - one that keeps every reported identity
  true - is invisible to the identity layer, so the accuracy column really measures
  accuracy;
* :func:`test_committed_audit_report_matches_a_fresh_run` re-derives the committed report
  and compares, so the numbers in ``docs/golden/mutation_audit.json`` cannot be hand-edited
  or left to rot.

The audit is expensive for a unit test (it measures all four golden hulls once, then sweeps
~900 re-blessed evaluations in pure arithmetic), so the sweep is driven from a session
fixture.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
AUDIT_PATH = ROOT / "tools" / "golden" / "mutation_audit.py"
COMMITTED = ROOT / "docs" / "golden" / "mutation_audit.json"


def _load_audit():
    spec = importlib.util.spec_from_file_location("tensorcad_mutation_audit", AUDIT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    # dataclasses resolve their own module through sys.modules while the class body runs
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


audit_module = _load_audit()


@pytest.fixture(scope="module")
def payload() -> dict:
    return audit_module.audit(quick=True)


@pytest.fixture(scope="module")
def rows(payload: dict) -> dict[str, dict]:
    return {f"{row['case']}/{row['metric']}": row for row in payload["envelope"]}


# ---------------------------------------------------------------------------------------
# the audit itself
# ---------------------------------------------------------------------------------------
def test_the_audit_has_no_policy_violation(payload: dict) -> None:
    assert payload["violations"] == []


def test_every_library_mutation_is_detected(payload: dict) -> None:
    controls = payload["controls"]
    assert len(controls) >= 5
    for control in controls:
        assert control["detected"], control
    mutated = [control for control in controls if control["expected_codes"]]
    assert mutated, "the battery must contain real mutations, not only the sanity control"
    for control in mutated:
        assert set(control["expected_codes"]) & set(control["codes"]), control
        # a mutation that leaks into another golden case means the cases are entangled
        assert control["failed_cases"] == [control["case"]], control


def test_no_mutation_control_passes(payload: dict) -> None:
    sanity = [c for c in payload["controls"] if not c["expected_codes"]]
    assert sanity, "a gate with no sanity control cannot show it is not failing everything"
    for control in sanity:
        assert control["detected"] is True
        assert control["codes"] == []


def test_every_detection_layer_is_exercised(payload: dict) -> None:
    """accuracy, invariant and drift must each be *shown* to work, not assumed."""

    exercised = set(payload["summary"]["layers_exercised"])
    assert {"accuracy", "invariant", "drift"} <= exercised, exercised


def test_the_declared_tolerance_is_the_real_threshold(rows: dict[str, dict]) -> None:
    for label, row in rows.items():
        tolerance = row["tolerance"]
        floor = row["coherent"]["blind_spot_floor"]
        ceiling = row["coherent"]["detection_ceiling"]
        assert ceiling is not None, f"{label} is never caught by accuracy"
        assert ceiling <= 3.0 * tolerance, (label, tolerance, ceiling)
        if floor is not None:
            # The threshold is one-sided where the baseline already sits below its
            # reference: a metric whose measured error is e tolerates a *positive*
            # perturbation up to (1 + tol)/(1 + e) - 1, which is larger than tol for e < 0.
            # The ladder is anchored on that estimate, so a metric whose own tolerance
            # decides brackets it inside the factor; a metric guarded by a tighter member of
            # its dependency cone is caught earlier, never later.
            error = row["reference_error"]
            estimate = max(tolerance, (1.0 + tolerance) / (1.0 + error) - 1.0)
            assert floor <= 1.1 * estimate, (label, tolerance, error, floor)


def test_incoherent_reports_are_caught_by_the_identities(rows: dict[str, dict]) -> None:
    """A number reported inconsistent with its own identity is a plumbing bug, and an
    identity is exact: it must fire orders of magnitude below the accuracy tolerance."""

    eager = [
        (label, row)
        for label, row in rows.items()
        if audit_module.identities_of(row["case"], row["metric"]) > 1
    ]
    assert eager, "no swept metric participates in an identity"
    for label, row in eager:
        floor = row["incoherent"]["blind_spot_floor"]
        assert row["coherent"]["invariants_hold_at_every_delta"], label
        if floor is not None:
            assert floor <= audit_module.INCOHERENT_CEILING, (label, floor)


def test_a_coherent_mutation_is_invisible_to_the_identities(rows: dict[str, dict]) -> None:
    """Otherwise the "accuracy" column would be measuring the invariant layer."""

    for label, row in rows.items():
        assert row["coherent"]["invariants_hold_at_every_delta"], label


def test_blind_spots_are_reported_not_hidden(payload: dict, rows: dict[str, dict]) -> None:
    """Every metric whose tolerance admits a 1 % coherent regression is named, with its
    number. A future tightening of the tolerance shrinks this list; it may never grow
    silently, because the list is compared against the committed report below."""

    named = payload["blind_spots_ge_1pct"]
    expected = [
        label
        for label, row in rows.items()
        if (row["coherent"]["blind_spot_floor"] or 0.0) >= 1e-2
    ]
    assert sorted(named) == sorted(expected)
    for label in named:
        row = rows[label]
        assert row["coherent"]["blind_spot_floor"] >= 1e-2


def test_the_audit_is_not_trivially_green(payload: dict) -> None:
    """A sweep that never perturbs anything would also report zero violations."""

    assert payload["summary"]["metrics_swept"] >= 20
    survived = [
        row
        for row in payload["envelope"]
        if row["coherent"]["blind_spot_floor"] is not None
    ]
    assert len(survived) >= len(payload["envelope"]) - 1, "the sweep caught everything even at 1e-12"
    assert payload["summary"]["worst_blind_spot"] > 1e-4


# ---------------------------------------------------------------------------------------
# the committed report must be the one this code produces
# ---------------------------------------------------------------------------------------
def test_committed_audit_report_matches_a_fresh_run(payload: dict) -> None:
    assert COMMITTED.is_file(), f"{COMMITTED} is missing; regenerate it with the audit tool"
    committed = json.loads(COMMITTED.read_text(encoding="utf-8"))
    assert committed["schema"] == audit_module.SCHEMA
    assert committed["violations"] == payload["violations"] == []

    def key(entry: dict) -> tuple[str, str]:
        return (entry["case"], entry["metric"])

    fresh_rows = {key(row): row for row in payload["envelope"]}
    stored_rows = {key(row): row for row in committed["envelope"]}
    assert set(fresh_rows) == set(stored_rows)
    for label, row in fresh_rows.items():
        stored = stored_rows[label]
        assert row["tolerance"] == pytest.approx(stored["tolerance"], rel=1e-9), label
        assert row["coherent"]["blind_spot_floor"] == pytest.approx(
            stored["coherent"]["blind_spot_floor"], rel=1e-9
        ), label
        fresh_incoherent = row["incoherent"]["blind_spot_floor"]
        stored_incoherent = stored["incoherent"]["blind_spot_floor"]
        assert (fresh_incoherent is None) == (stored_incoherent is None), label
        if fresh_incoherent is not None and stored_incoherent is not None:
            assert fresh_incoherent == pytest.approx(stored_incoherent, rel=1e-9), label

    fresh_controls = {control["name"]: control for control in payload["controls"]}
    stored_controls = {control["name"]: control for control in committed["controls"]}
    assert set(fresh_controls) == set(stored_controls)
    for name, control in fresh_controls.items():
        assert control["detected"] == stored_controls[name]["detected"], name
        assert control["codes"] == stored_controls[name]["codes"], name


def test_committed_markdown_is_rendered_from_the_committed_json() -> None:
    markdown = (ROOT / "docs" / "golden" / "mutation_audit.md").read_text(encoding="utf-8")
    committed = json.loads(COMMITTED.read_text(encoding="utf-8"))
    assert markdown == audit_module.render_markdown(committed)