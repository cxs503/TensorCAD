"""Guard: the shipped Series-60 gate-power record must still match the harness.

The anchor ``GAP-CAD-T2-PARAM-SERIES`` asks for a gate-power mutation kill rate of
at least 80 %.  The record that carries that number
(``artifacts/series60_gate_power/gate_power_series60.json``) used to be a prose
table in ``docs/series60_gate_power.md``, and the 21 Series-60 mutants it names
were only partly executed: seven of them (``S60-15`` .. ``S60-21``, the entry and
family-comparison claims) had no recorded run at all, so the acceptance tool
reported the clause as ``not_run`` -- *unknown*, which is not a pass.

This module pins the three ways that record can rot:

* the harness grows a mutant and the record is not regenerated (coverage gap);
* the record describes an injection the harness no longer contains, so a "kill"
  no longer refers to anything (drift between claim and code);
* an entry loses the injection or the measured exit status, which turns a kill
  rate into a summary of claims (provenance).

It deliberately does *not* re-run the harness -- a test suite that mutated
production sources on every ``pytest`` invocation would be a hazard.  The
regeneration command is quoted in the failure message instead.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "artifacts" / "series60_gate_power" / "gate_power_series60.json"
HARNESS = ROOT / "tools" / "hull_gate_power" / "mutate.py"
TOOL = ROOT / "tools" / "series60_family" / "anchor_acceptance.py"

REGENERATE = ("python3 tools/hull_gate_power/mutate.py "
              "--out artifacts/series60_gate_power/gate_power_series60.json --only "
              + ",".join(f"S60-{index}-" for index in range(1, 22)))


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def harness():
    return _load(HARNESS, "hull_gate_power_mutate")


@pytest.fixture(scope="module")
def tool():
    return _load(TOOL, "anchor_acceptance_record_guard")


@pytest.fixture(scope="module")
def record():
    assert RECORD.exists(), (
        f"the Series-60 mutation record is missing: {RECORD}. "
        f"Regenerate it with: {REGENERATE}"
    )
    raw = json.loads(RECORD.read_text())
    # the harness also records unmutated baselines; they are controls, not mutants
    return {key: entry for key, entry in raw.items() if not key.startswith("baseline")}


def series60_ids(harness_module) -> set[str]:
    """The Series-60 mutants the harness defines, by id."""
    return {m["id"] for m in harness_module.MUTANTS if m["id"].startswith("S60-")}


def uncovered(record_entries: dict, defined: set[str]) -> list[str]:
    """Mutants the harness defines but the record never ran."""
    return sorted(defined - set(record_entries))


# ------------------------------------------------------------------ coverage
def test_record_covers_every_series60_mutant_the_harness_defines(harness, record):
    defined = series60_ids(harness)
    assert defined, "the harness defines no S60-* mutants: this guard has nothing to pin"
    missing = uncovered(record, defined)
    assert not missing, (
        f"{len(missing)} Series-60 mutant(s) defined by the harness were never run: "
        f"{missing}. Regenerate with: {REGENERATE}"
    )


def test_record_carries_no_mutant_the_harness_does_not_define(harness, record):
    defined = series60_ids(harness)
    stray = sorted(set(record) - defined)
    assert not stray, (
        f"the record names mutants the harness no longer contains: {stray}; a kill "
        "of a deleted mutant is not evidence about this tree"
    )


def test_the_series_is_exactly_the_21_mutants_the_family_claims(record):
    assert len(record) == 21, sorted(record)
    assert {entry["mutant"]["id"] for entry in record.values()} == set(record)


# ---------------------------------------------------------------- provenance
def test_every_entry_carries_the_injection_the_harness_makes(harness, record):
    """A kill must refer to the mutation the harness currently defines."""
    definitions = {m["id"]: m for m in harness.MUTANTS}
    for mutant_id, entry in sorted(record.items()):
        definition = definitions[mutant_id]
        assert entry["mutant"]["old"] == definition["old"], mutant_id
        assert entry["mutant"]["new"] == definition["new"], mutant_id
        assert entry["mutant"]["file"] == definition["file"], mutant_id
        assert entry["mutant"]["claim"] == definition["claim"], mutant_id
        assert entry["mutant"]["tests"] == definition["tests"], mutant_id


def test_every_entry_carries_a_measured_exit_status(tool, record):
    defects = {mutant_id: tool.provenance_defects(entry)
               for mutant_id, entry in sorted(record.items())}
    bad = {k: v for k, v in defects.items() if v}
    assert not bad, bad


def test_every_kill_names_the_test_that_killed_it(record):
    for mutant_id, entry in sorted(record.items()):
        run = entry["run"]
        assert run["killed"] is True, f"{mutant_id} survived: {run['summary']}"
        assert run["returncode"] != 0, mutant_id
        assert run["failed"], f"{mutant_id} was killed but no test is named"


# -------------------------------------------------------------- the anchor bar
def test_kill_rate_meets_the_anchor_bar(tool, record):
    clause = tool.gate_power_clause(RECORD)
    assert clause["status"] == "pass", clause.get("detail") or clause["survivors"]
    assert clause["kill_rate"] >= tool.GATE_POWER_MINIMUM
    assert clause["total"] == len(record)
    assert clause["survivors"] == []


def test_the_anchor_verdict_is_no_longer_blocked_on_this_clause(tool, record):
    published = json.loads((ROOT / "docs" / "data" / "series60_published_parents.json").read_text())
    report = tool.build_report(published, RECORD)
    clause = next(c for c in report["clauses"]
                  if c["id"] == "gate_power_mutation_kill_rate")
    assert clause["status"] == "pass", clause
    assert "gate_power_mutation_kill_rate" not in report["clauses_not_run"]
    # The external shape anchor is still the honest gap; it stays withheld.
    assert "external_shape_anchor_todd_offsets" in report["clauses_withheld"]


# ------------------------------------------------------------ negative control
def test_a_record_missing_one_mutant_is_caught():
    """The coverage check must bite: 21 minus one is not 21."""
    defined = {f"S60-{index}-x" for index in range(1, 22)}
    entries = {mutant_id: {} for mutant_id in sorted(defined)[1:]}
    assert uncovered(entries, defined) == [sorted(defined)[0]]
    assert uncovered({mutant_id: {} for mutant_id in defined}, defined) == []