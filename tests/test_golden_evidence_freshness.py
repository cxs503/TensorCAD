"""The committed golden evidence must describe the tree it ships with.

``docs/golden/results.{json,md}``, ``docs/golden/mutation_audit.{json,md}`` and
``docs/golden/tensorcad_geometry_golden.json`` are generated artifacts committed to the
repository: they are the record a reader consults instead of re-running the suite.  Two
ways such a record can lie, both observed in this repository before this test existed:

**Stale.** ``inclined_box`` was registered as a fifth golden case (code + baseline) while
``docs/golden/results.json`` still reported four cases: the baseline entry and the code
agreed, so ``python -m tensorcad.golden --check`` stayed green, and nothing pointed at the
report that had silently stopped describing the gate.

**Machine-bound.** Both reports recorded the *absolute* path of the worktree that generated
them (``/root/tp-agent/work/<mission>/docs/golden/...``), so re-running the documented
command could never reproduce the committed file and the artifact published the local
layout of the machine that produced it.

The checks below are cheap on purpose: no geometry is re-measured.  The baseline is the
oracle, and ``--check`` (already run by CI and the pre-commit hook) proves the baseline
matches the code.  What is asserted here is the missing link: the reports match the
baseline.  ``test_the_checker_rejects_a_stale_report`` and
``test_the_checker_rejects_a_machine_bound_report`` prove these checks can fail at all.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pytest

from tensorcad import golden

ROOT = Path(__file__).resolve().parents[1]
GOLDEN = ROOT / "docs" / "golden"
RESULTS = GOLDEN / "results.json"
RESULTS_MD = GOLDEN / "results.md"
AUDIT = GOLDEN / "mutation_audit.json"
BASELINE = GOLDEN / "tensorcad_geometry_golden.json"

#: Substrings that only ever appear in a report produced on one particular machine: a
#: private worktree, an NFS share, a home directory.  Their presence means the artifact is
#: not reproducible, whichever machine reads it.
MACHINE_MARKERS = ("/root/", "/home/", "/tmp/", "/nfs/", "tp-agent")


def _problems(report: dict[str, Any], baseline: dict[str, Any], registered: list[str]) -> list[str]:
    """Why the committed report does not describe this tree (empty list: it does)."""
    problems: list[str] = []
    if not report.get("passed"):
        problems.append("the report does not record a passing gate")
    recorded = [case["name"] for case in report.get("cases", [])]
    if sorted(recorded) != sorted(registered):
        problems.append(
            f"the report covers {sorted(recorded)} but the suite registers {sorted(registered)}"
        )
    baseline_path = str(report.get("baseline", ""))
    if os.path.isabs(baseline_path) or any(m in baseline_path for m in MACHINE_MARKERS):
        problems.append(f"the report records a machine-specific baseline path {baseline_path!r}")
    for case in report.get("cases", []):
        entry = (baseline.get("cases") or {}).get(case["name"])
        if entry is None:
            problems.append(f"{case['name']}: no baseline entry to compare against")
            continue
        if case.get("geometry_sha256") != entry.get("geometry_sha256"):
            problems.append(f"{case['name']}: report measures a different geometry")
        for level, measured in (case.get("measurements") or {}).items():
            blessed = (entry.get("measurements") or {}).get(level) or {}
            if set(measured) != set(blessed):
                problems.append(f"{case['name']}/{level}: reported metrics differ from the baseline")
                continue
            for name, value in measured.items():
                if value != blessed[name]:
                    problems.append(f"{case['name']}/{level}/{name}: {value!r} != blessed {blessed[name]!r}")
    return problems


def _registered() -> list[str]:
    return [case.name for case in golden.golden_cases()]


@pytest.fixture(scope="module")
def report() -> dict[str, Any]:
    if not RESULTS.is_file():  # pragma: no cover - the artifact is committed
        pytest.skip("docs/golden/results.json is not present")
    return json.loads(RESULTS.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def baseline() -> dict[str, Any]:
    return json.loads(BASELINE.read_text(encoding="utf-8"))


def test_the_committed_report_describes_the_current_tree(report, baseline) -> None:
    """The one check the CI gate was missing: report vs baseline vs registered cases."""
    assert _problems(report, baseline, _registered()) == []


def test_the_checker_rejects_a_stale_report(report, baseline) -> None:
    """Mutation of the check itself: drop a case and the freshness check must go red."""
    stale = json.loads(json.dumps(report))
    stale["cases"] = [case for case in stale["cases"] if case["name"] != "inclined_box"]
    assert _problems(stale, baseline, _registered())


def test_the_checker_rejects_a_machine_bound_report(report, baseline) -> None:
    machine_bound = json.loads(json.dumps(report))
    machine_bound["baseline"] = "/root/tp-agent/work/m00000/docs/golden/tensorcad_geometry_golden.json"
    assert any("machine-specific" in p for p in _problems(machine_bound, baseline, _registered()))


def test_the_report_recorded_a_repository_relative_baseline(report) -> None:
    expected = "docs/golden/tensorcad_geometry_golden.json"
    assert golden.display_path(golden.BASELINE) == expected
    assert report["baseline"] == expected


def test_display_path_keeps_a_baseline_outside_the_repository_absolute(tmp_path: Path) -> None:
    """Outside the repository the absolute path is the only name that identifies the file."""
    outside = tmp_path / "baseline.json"
    assert golden.display_path(outside) == str(outside)
    assert golden.display_path(golden.BASELINE) == "docs/golden/tensorcad_geometry_golden.json"


def test_no_generated_evidence_names_a_private_machine_path() -> None:
    """Scoped to the artifacts generated by ``tools/golden/*`` in this round.

    ``docs/golden/gate_reach.md`` and ``docs/golden/grid_convergence.md`` are hand-written
    records from earlier missions and still cite ``/tmp`` log paths; they are outside the
    generator this test guards and are listed as residual in
    ``docs/golden/verification_report.md`` rather than silently tolerated here.
    """
    generated = [RESULTS, RESULTS_MD, AUDIT, AUDIT.with_suffix(".md")]
    offenders = {
        path.name: marker
        for path in generated
        if path.is_file()
        for marker in MACHINE_MARKERS
        if marker in path.read_text(encoding="utf-8", errors="replace")
    }
    assert offenders == {}


def test_the_tolerances_in_the_report_are_the_ones_in_force(report, baseline) -> None:
    """A report quoting a wider band than the baseline records would hide a regression."""
    mismatches: list[str] = []
    for case in report["cases"]:
        entry = baseline["cases"][case["name"]]
        for metric in case["metrics"]:
            if metric["kind"] == "drift":
                continue
            # A three-level case records both bands: the two-grid one under ``tolerances``
            # and the GCI one under ``extrapolated``.  The report carries both, so the
            # comparison has to read the arm the metric came from.
            if metric["kind"] == "invariant":
                stored = ((entry.get("invariants") or {}).get(metric["name"]) or {}).get("tolerance")
            elif metric.get("origin") == "richardson":
                stored = ((entry.get("extrapolated") or {}).get(metric["name"]) or {}).get("tolerance")
            else:
                stored = ((entry.get("tolerances") or {}).get(metric["name"]) or {}).get("value")
            if stored is None or stored != metric["tolerance"]:
                mismatches.append(f"{case['name']}/{metric['name']}: {metric['tolerance']!r} != {stored!r}")
            if metric["passed"] != (abs(metric["relative_error"]) <= metric["tolerance"]):
                mismatches.append(f"{case['name']}/{metric['name']}: passed flag contradicts the band")
    assert mismatches == []


def test_the_markdown_table_quotes_every_case_and_assertion(report) -> None:
    text = RESULTS_MD.read_text(encoding="utf-8")
    expected_rows = [f"| {case['name']} | {m['name']} |" for case in report["cases"] for m in case["metrics"]]
    missing = [row for row in expected_rows if row not in text]
    assert missing == []
    assert "gate verdict: **PASS**" in text


def test_the_mutation_audit_report_is_repository_relative() -> None:
    if not AUDIT.is_file():  # pragma: no cover - the artifact is committed
        pytest.skip("docs/golden/mutation_audit.json is not present")
    payload = json.loads(AUDIT.read_text(encoding="utf-8"))
    assert payload["baseline"] == "docs/golden/tensorcad_geometry_golden.json"