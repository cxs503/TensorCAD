"""Tests for the ``GAP-CAD-T2-PARAM-SERIES`` acceptance verdict tool.

These tests exist because the anchor's acceptance used to be checkable only by
reading prose, and because a board audit (m17308) found the machine-readable
``role`` field in the five-Cb artifact let a *closure* check be read as an
external verification.  So the assertions here are about invariants, not about
"the script ran":

* a clause that is ``not_run`` must never be reported as passed;
* a comparison whose reference is an input of the model must be graded
  ``closure``, and the report must contain no clause that upgrades it;
* the external shape anchor must be ``withheld``, not silently dropped;
* negative controls: a perturbed published table and a weak mutation record must
  both flip the verdict.
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
TOOL = ROOT / "tools" / "series60_family" / "anchor_acceptance.py"
AUDIT = ROOT / "tools" / "series60_family" / "audit.py"
DATA = ROOT / "docs" / "data" / "series60_published_parents.json"


def load_tool():
    spec = importlib.util.spec_from_file_location("anchor_acceptance", TOOL)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def tool():
    return load_tool()


@pytest.fixture(scope="module")
def published():
    return json.loads(DATA.read_text())


@pytest.fixture(scope="module")
def gradient(tool, published):
    return tool.grade_gradient(published)


# --------------------------------------------------------------------- clauses
def test_five_parent_forms_generate_and_match_the_published_ladder(tool, gradient):
    assert gradient["passed"], gradient["worst_quantity"]
    assert [parent["cb"] for parent in gradient["parents"]] == [0.60, 0.65, 0.70, 0.75, 0.80]
    for parent in gradient["parents"]:
        assert parent["validation_passed"], parent["validation_errors"]
        assert parent["geometry"]["passed"]
        for key, comparison in parent["comparisons"].items():
            assert comparison["within_bar"], (parent["cb"], key, comparison["relative"])
    # The 3% bar must actually bite: the worst row is well inside it but not zero.
    assert 0.0 < gradient["worst_relative"] < tool.TOLERANCE_PUBLISHED


def test_every_published_comparison_is_graded_closure(gradient):
    for parent in gradient["parents"]:
        for key, comparison in parent["comparisons"].items():
            assert comparison["grade"] == "closure", (parent["cb"], key)
    assert gradient["grade"] == "closure"


def test_grades_are_from_the_closed_vocabulary(tool, gradient):
    for consent in gradient["parents"]:
        for comparison in consent["comparisons"].values():
            assert comparison["grade"] in tool.GRADES


def test_geometry_gate_reports_cb_length_beam_and_draft(tool, published):
    parent = published["parents"][2]
    # A request the loft cannot honour must fail the geometry gate rather than
    # passing silently: a deliberately impossible beam.
    broken = dict(length=tool.L, beam=tool.B * 3.0, draft=tool.T, cb=parent["cb"],
                  nx=41, nz=13, cwp=parent["cwp"], lcb=parent["lcb_over_length"])
    from tensorcad.series60 import Series60Parameters, generate_series60
    model = generate_series60(Series60Parameters(**broken))
    gate = tool.geometry_precondition(model, broken)
    assert not gate["passed"]
    assert {row["quantity"] for row in gate["rows"]} == {
        "length", "beam", "draft", "block_coefficient"}
    assert all(row["grade"] == "closure" for row in gate["rows"])


# ------------------------------------------------------------ negative controls
def test_a_wrong_published_midship_coefficient_fails_the_clause(tool, published):
    # ``cm`` is a published column the shape solve does not consume, so corrupting
    # it is a clean negative control: no feasibility rejection can mask the error.
    corrupted = json.loads(json.dumps(published))
    for parent in corrupted["parents"]:
        parent["cm"] = parent["cm"] * 1.10
    clause = tool.grade_gradient(corrupted)
    assert not clause["passed"], "a 10% Cm error must trip the 3% bar"
    worst = max(abs(p["comparisons"]["midship_coefficient"]["relative"])
                for p in clause["parents"])
    assert worst > tool.TOLERANCE_PUBLISHED


def test_an_input_target_follows_the_table_which_is_why_it_is_closure(tool, published):
    """Moving an *input* target cannot fail the comparison -- it moves the model.

    This is the negative control that pins the grade: ``cwp`` and ``lcb`` are
    solved for, so a corrupted table is reproduced faithfully and the relative
    error stays near zero.  If this ever starts failing, those quantities have
    become genuine predictions and their grade must be raised to ``external``.
    """
    corrupted = json.loads(json.dumps(published))
    for parent in corrupted["parents"]:
        parent["lcb_over_length"] = parent["lcb_over_length"] + 0.02
    clause = tool.grade_gradient(corrupted)
    assert clause["passed"], "an input target is followed, not checked"
    worst = max(abs(p["comparisons"]["lcb_over_length"]["relative"])
                for p in clause["parents"])
    assert worst < tool.TOLERANCE_PUBLISHED


def test_an_infeasible_published_waterplane_is_rejected_loudly(tool, published):
    # The family cannot reach a Cwp 10 % below the ladder: it must refuse, not clamp.
    corrupted = json.loads(json.dumps(published))
    for parent in corrupted["parents"]:
        parent["cwp"] = parent["cwp"] * 0.90
    with pytest.raises(ValueError):
        tool.grade_gradient(corrupted)


# ------------------------------------------------------------------ gate power
def _entry(mutant_id: str, killed: bool) -> dict:
    """One record entry in the shape ``mutate.py`` actually writes."""
    return {
        "status": "ok",
        "mutant": {
            "id": mutant_id, "kind": "model", "file": "series60",
            "old": "    return keep * h + transom * 0.5",
            "new": f"    return keep * h + transom * 0.55  # mutant {mutant_id}",
            "claim": "the transom's share of the aft section area",
            "tests": ["tests/test_series60.py"],
        },
        "run": {
            "returncode": 1 if killed else 0,
            "killed": killed,
            "summary": "3 failed, 53 passed" if killed else "56 passed",
            "failed": ["tests/test_series60.py::test_transom"] if killed else [],
            "wall_seconds": 21.5,
        },
    }


def _record(tmp_path: Path, killed: int, survived: int, name: str = "gate_power.json") -> Path:
    entries = {}
    for index in range(killed):
        entries[f"M-K{index}"] = _entry(f"M-K{index}", killed=True)
    for index in range(survived):
        entries[f"M-S{index}"] = _entry(f"M-S{index}", killed=False)
    path = tmp_path / name
    path.write_text(json.dumps(entries))
    return path


def test_missing_mutation_record_is_not_run_and_not_passed(tool, tmp_path):
    clause = tool.gate_power_clause(tmp_path / "absent.json")
    assert clause["status"] == "not_run"
    assert clause["passed"] is None
    assert clause["kill_rate"] is None


def test_mutation_record_below_the_bar_fails(tool, tmp_path):
    clause = tool.gate_power_clause(_record(tmp_path, killed=11, survived=3))
    assert clause["status"] == "fail"
    assert clause["kill_rate"] == pytest.approx(11 / 14)
    assert clause["survivors"] == ["M-S0", "M-S1", "M-S2"]


def test_mutation_record_at_the_bar_passes(tool, tmp_path):
    clause = tool.gate_power_clause(_record(tmp_path, killed=8, survived=2))
    assert clause["status"] == "pass"
    assert clause["kill_rate"] == pytest.approx(0.80)


def test_a_hand_transcribed_record_is_not_evidence(tool, tmp_path):
    """The defect this guards: a fabricated record used to be read as a kill rate.

    ``tools/hull_gate_power/generate_s60_artifact.py`` wrote
    ``{"status": "ok", "run": {"killed": true}}`` for 24 mutants that it never ran,
    and this clause scored it as an *external* 100 % kill rate with a sha256 stamp.
    A record must carry the injection it made and the exit status it measured.
    """
    fabricated = {f"S60-{i}": {"status": "ok", "run": {"killed": True,
                                                       "wall_seconds": 18.2}}
                  for i in range(1, 15)}
    path = tmp_path / "fabricated.json"
    path.write_text(json.dumps(fabricated))
    clause = tool.gate_power_clause(path)
    assert clause["status"] == "unverified"
    assert clause["passed"] is False
    assert clause["kill_rate"] is None
    assert len(clause["unverifiable"]) == 14
    assert "mutant" in clause["detail"]


def test_a_no_op_injection_cannot_be_a_kill(tool, tmp_path):
    entry = _entry("S60-x", killed=True)
    entry["mutant"]["new"] = entry["mutant"]["old"]
    path = tmp_path / "noop.json"
    path.write_text(json.dumps({"S60-x": entry}))
    clause = tool.gate_power_clause(path)
    assert clause["status"] == "unverified"
    assert "no-op" in clause["unverifiable"]["S60-x"][0]


def test_a_kill_must_agree_with_the_measured_exit_status(tool, tmp_path):
    entry = _entry("S60-y", killed=True)
    entry["run"]["returncode"] = 0       # a clean exit cannot be a kill
    path = tmp_path / "contradiction.json"
    path.write_text(json.dumps({"S60-y": entry}))
    clause = tool.gate_power_clause(path)
    assert clause["status"] == "unverified"
    assert any("contradicts" in d for d in clause["unverifiable"]["S60-y"])


def test_one_unverifiable_entry_taints_the_rate(tool, tmp_path):
    record = {"M-K0": _entry("M-K0", killed=True),
              "M-K1": _entry("M-K1", killed=True),
              "M-S0": _entry("M-S0", killed=False),
              "M-X": {"status": "ok", "run": {"killed": True}}}
    path = tmp_path / "mixed.json"
    path.write_text(json.dumps(record))
    clause = tool.gate_power_clause(path)
    assert clause["status"] == "unverified"
    assert clause["passed"] is False
    assert clause["survivors"] == ["M-S0"]      # measured survivors still reported


def test_no_tool_invents_a_gate_power_record(tool, tmp_path):
    """There must be exactly one producer: the harness that patches and runs."""
    assert not (ROOT / "tools" / "hull_gate_power" / "generate_s60_artifact.py").exists()
    writers = [
        path.name
        for path in (ROOT / "tools").rglob("*.py")
        if "gate_power_series60.json" in path.read_text(errors="ignore")
        and "json.dumps" in path.read_text(errors="ignore")
        and "pytest" not in path.read_text(errors="ignore")
    ]
    assert writers == []


# ------------------------------------------------------------- full assembly
def test_verdict_is_incomplete_when_a_clause_was_not_run(tool, published, tmp_path):
    report = tool.build_report(published, tmp_path / "absent.json")
    assert report["verdict"] == "incomplete"
    assert "gate_power_mutation_kill_rate" in report["clauses_not_run"]
    assert "external_shape_anchor_todd_offsets" in report["clauses_withheld"]


def test_verdict_is_incomplete_when_the_record_is_unverifiable(tool, published, tmp_path):
    fabricated = {f"S60-{i}": {"status": "ok", "run": {"killed": True}} for i in range(3)}
    path = tmp_path / "fabricated.json"
    path.write_text(json.dumps(fabricated))
    report = tool.build_report(published, path)
    assert report["verdict"] == "incomplete"
    assert "gate_power_mutation_kill_rate" in report["clauses_unverified"]
    assert "gate_power_mutation_kill_rate" not in report["failed_clauses"]


def test_no_closure_is_reported_as_an_external_verification(tool, published, tmp_path):
    report = tool.build_report(published, tmp_path / "absent.json")
    assert report["closure_where_external_should_be"] == []
    for clause in report["clauses"]:
        for _quantity, grade in tool._graded_items(clause):
            assert grade in tool.GRADES


def test_external_shape_anchor_is_withheld_not_passed(tool):
    clause = tool.external_anchor_clause()
    assert clause["status"] == "withheld"
    assert clause["passed"] is None
    assert "Todd" in clause["statement"]


def test_documentation_clause_detects_a_missing_page(tool):
    assert tool.documentation_clause()["passed"]
    assert not tool.documentation_clause(("docs/series60.md", "docs/does-not-exist.md"))["passed"]


# ------------------------------------------------- F5 regression: --out outside
def test_audit_tool_accepts_an_output_directory_outside_the_checkout(tmp_path):
    out = tmp_path / "outside_root"
    completed = subprocess.run(
        [sys.executable, str(AUDIT), "--out", str(out)],
        cwd=ROOT, capture_output=True, text=True, timeout=600,
    )
    assert completed.returncode == 0, completed.stderr[-3000:]
    payload = json.loads((out / "five_cb_gradient.json").read_text())
    # Absolute path recorded verbatim instead of raising relative_to().
    assert payload["artifacts"]["json"] == str(out / "five_cb_gradient.json")
    assert (out / "figures" / "five_cb_gradient.png").exists()