"""The PR-#45 cure is checkable, so check it: the union merge, re-audited in the suite.

``tools/verify_merge_m3391.py`` answers the merge question with numbers (file-level union,
line/name-level survival, both routes still computing, tests not below baseline). These
tests run the same audit inside pytest so a later merge cannot quietly reintroduce a
one-sided file set, ``CHANGELOG`` conflict markers, or a public surface that lost a
parent's names. The audit needs a two-parent commit; it is skipped where that does not
hold (e.g. a squashed export).
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools" / "verify_merge_m3391.py"


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout


def _needs_two_parents() -> None:
    try:
        _git("rev-parse", "HEAD^2")
    except subprocess.CalledProcessError:  # pragma: no cover - squashed export
        pytest.skip("HEAD is not a merge commit: nothing to audit")


def test_the_audit_tool_runs_and_reports_a_clean_union() -> None:
    """End to end: the tool exits 0, which is what the merge-governor re-check assumes."""
    _needs_two_parents()
    result = subprocess.run(
        [sys.executable, str(TOOL)], cwd=ROOT, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    data = json.loads((ROOT / "artifacts" / "m3781" / "merge_verification.json").read_text())

    files = data["files"]
    assert files["missing_from_main"] == []
    assert files["missing_from_branch"] == []
    assert files["invented"] == []
    assert data["api"]["main_only_lost"] == []
    assert data["api"]["branch_only_lost"] == []
    for row in data["lines"]:
        assert row["missing"] == 0, row


def test_the_branch_route_and_main_routes_both_compute() -> None:
    """A merge that compiled but neutered one side would fail here, not in review."""
    _needs_two_parents()
    data = json.loads((ROOT / "artifacts" / "m3781" / "merge_verification.json").read_text())

    branch = data["branch_route"]
    assert branch["max_abs_error_m"] < 1e-12  # GZ against the wall-sided closed form
    assert branch["displacement_volume"] > 0.0
    assert branch["max_volume_error"] < 1e-12

    main = data["main_routes"]
    for row in main["hull_vs_closed_form"]:
        assert abs(row["I_T_rel_err"]) < 1e-2
        assert abs(row["I_L_rel_err"]) < 1e-2
    assert main["transom"]["block_coefficient"] > 0.0
    assert main["wigley"]["midship_coefficient"] > 0.0


def test_the_changelog_carries_no_conflict_markers() -> None:
    text = (ROOT / "CHANGELOG.md").read_text()
    for marker in ("<<<<<<<", "=======", ">>>>>>>"):
        assert marker not in text
    assert "0.21.0" in text