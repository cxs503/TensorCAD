"""The library's cached blocks must follow from its cells, and the audit must read the dissent.

Anchor ``REPRO-GEOM-FLEET-REFS``.  Two defects motivated this file, both introduced by the
library caching quantities that are *functions of* its cell-level report:

* ``closure`` / ``adjudication`` were left at round 2's answer while round 3 rewrote the cells
  (KVLCC2 ``Bwl`` 580→58 m, DTC/JSR-2017 added, JBC volume refined), so the committed file
  claimed a 900 % source defect that its own report contradicted, and the published table
  printed it as a finding.  ``merge_into`` now refreshes the blocks and
  ``refs_freshness.py check`` fails on any drift.
* the power-of-ten typo detector read only the *winning* family value; once a second lineage
  printed the corrected KVLCC2 beam, the majority cluster closed and the ``580`` typo — still
  printed in two entries — went unaudited.  Two assertions had been parked as ``xfail``.

Assertions here are invariants (the cache equals its inputs; the detector attributes exactly
the slips present and invents none; recomputation does not mutate what it inspects), not
"the script ran".
"""

from __future__ import annotations

import copy
import importlib.util
import json
import pathlib

import pytest

from tensorcad.fleet.crossref import (
    Cell,
    adjudicate,
    closure_check,
    coverage_summary,
    cross_reference,
    refresh_derived,
    stale_derived,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
LIBRARY = ROOT / "benchmarks/verified/fleet_refs/crossref.json"
INGEST = ROOT / "benchmarks/verified/fleet_refs/oa_ingest.json"


def _module(name: str, path: pathlib.Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


oa = _module("crossref_oa_tool", ROOT / "tools/baseline_fleet/crossref_oa.py")
freshness = _module("refs_freshness_tool", ROOT / "tools/baseline_fleet/refs_freshness.py")


@pytest.fixture()
def library() -> dict:
    return json.loads(LIBRARY.read_text())


def _cell(quantity: str, value: float, source: str) -> Cell:
    return Cell(
        quantity=quantity,
        value=value,
        raw=str(value),
        label=f"{quantity} row",
        source=source,
        url="https://example.invalid/table",
        sha256="0" * 64,
        column=0,
    )


def _report(**quantities: tuple[float, float]) -> dict:
    """One hull, two families, from ``quantity=(family_a_value, family_b_value)``."""
    extractions = {"HULL": {}}
    for quantity, (left, right) in quantities.items():
        extractions["HULL"].setdefault("t2015_kcs", {})[quantity] = _cell(quantity, left, "t2015_kcs")
        extractions["HULL"].setdefault("simman2014_kcs", {})[quantity] = _cell(
            quantity, right, "simman2014_kcs"
        )
    return cross_reference(extractions)


# --------------------------------------------------------------- cache freshness (defect 1)


def test_committed_library_has_no_stale_block(library: dict) -> None:
    """The regression guard: a round that rewrites cells must rewrite the cached blocks."""
    view = freshness.derived_view(library)
    assert view["stale"] == [], view["differences"]


def test_check_exits_non_zero_on_a_stale_library(library: dict, tmp_path: pathlib.Path) -> None:
    library["closure"]["KVLCC2"]["SIMMAN-SERIES"]["rel_error"] = 8.99999
    path = tmp_path / "crossref.json"
    path.write_text(json.dumps(library))
    assert freshness.main(["check", "--library", str(path)]) == 1


def test_recomputing_does_not_mutate_the_report_it_inspects(library: dict) -> None:
    """``adjudicate`` stamps ``closure_defect`` into the cells; a staleness *check* must not."""
    payload = copy.deepcopy(library)
    before = json.dumps(payload["report"], sort_keys=True)
    stale_derived(payload)
    assert json.dumps(payload["report"], sort_keys=True) == before


def test_stale_derived_names_the_leaves_that_disagree(library: dict) -> None:
    payload = copy.deepcopy(library)
    payload["closure"]["KCS"]["NMRI-T2015"]["rel_error"] = 0.5
    result = stale_derived(payload)
    assert result["stale"] == ["closure"]
    assert any("/KCS/NMRI-T2015/rel_error" in line for line in result["differences"]["closure"])


def test_refresh_rewrites_every_block_from_the_report(library: dict) -> None:
    payload = copy.deepcopy(library)
    payload["closure"] = {}
    payload["adjudication"] = {}
    payload["coverage"] = {}
    assert refresh_derived(payload) == ["coverage", "closure", "adjudication"]
    assert payload["coverage"] == coverage_summary(payload["report"])
    assert payload["closure"] == closure_check(payload["report"])
    assert stale_derived(payload)["stale"] == []


def test_merge_path_cannot_leave_the_blocks_stale(library: dict) -> None:
    """``crossref_oa.merge_into`` recomputes ``closure``/``adjudication`` itself."""
    payload = copy.deepcopy(library)
    payload["closure"] = {"KVLCC2": {"SIMMAN-SERIES": {"rel_error": 9.0, "consistent": False}}}
    payload["adjudication"] = {"scale_typos": []}
    ingest = json.loads(INGEST.read_text())
    merged = oa.merge_into(payload, ingest, fleet_specs=None)
    assert merged["closure"] == closure_check(merged["report"])
    assert merged["adjudication"]["scale_typos"] == adjudicate(merged["report"])["scale_typos"]
    assert stale_derived(merged)["stale"] == []


# ------------------------------------------------------- dissent-aware slip scan (defect 2)


def test_kvlcc2_beam_typo_is_found_in_the_dissent(library: dict) -> None:
    verdict = adjudicate(copy.deepcopy(library["report"]))
    slips = [row for row in verdict["scale_typos"] if row["ship"] == "KVLCC2"]
    assert len(slips) == 1
    slip = slips[0]
    assert (slip["quantity"], slip["origin"]) == ("beam", "dissent")
    assert (slip["printed"], slip["divide_by"], slip["corrected"]) == (580.0, 10.0, 58.0)
    assert slip["family_value"] == pytest.approx(58.0)
    assert abs(slip["corrected_residual"]) <= 0.01 < abs(slip["raw_residual"])


def test_dissent_slip_is_bound_to_the_cells_that_print_it(library: dict) -> None:
    report = copy.deepcopy(library["report"])
    adjudicate(report)
    entries = report["KVLCC2"]["quantities"]["beam"]["entries"]
    flagged = [entry for entry in entries if "closure_defect" in entry["cell"]]
    assert len(flagged) == 2  # both SIMMAN editions, neither corrected reprint
    assert {f["cell"]["value"] for f in flagged} == {580.0}
    assert {f["cell"]["closure_defect"]["origin"] for f in flagged} == {"dissent"}
    # the votes are untouched by the detector
    assert report["KVLCC2"]["quantities"]["beam"]["family_values"]["SIMMAN-SERIES"] == 58.0


def test_adjudicating_twice_gives_the_same_findings(library: dict) -> None:
    """``closure_defect`` is stamped on the report, so a re-run must be idempotent."""
    report = copy.deepcopy(library["report"])
    first = adjudicate(report)["scale_typos"]
    second = adjudicate(report)["scale_typos"]
    assert first == second


def test_a_disagreement_that_is_not_a_digit_slip_is_not_attributed(report: dict | None = None) -> None:
    """580→58 closes; 63 m instead of 58 m does not, so no typo may be reported."""
    hull = _report(
        displacement_volume=(100.0, 100.0),
        cb=(0.6, 0.6),
        length_pp=(10.0, 10.0),
        beam=(1.6666666666666667, 1.6666666666666667),
        draft=(10.0, 10.0),
    )
    assert adjudicate(hull)["scale_typos"] == []
    hull["HULL"]["quantities"]["beam"]["family_disagreement"] = {
        "SIMMAN-SERIES": {"chosen": 1.6667, "dissent": [6.3], "cluster_size": 1, "values": [1.6667, 6.3]}
    }
    assert adjudicate(hull)["scale_typos"] == []


def test_a_dissent_slip_that_closes_is_reported_once_per_printed_value() -> None:
    """Two entries printing the same typo produce one finding, not one per entry."""
    report = _report(
        displacement_volume=(312622.0, 312622.0),
        cb=(0.8098, 0.8098),
        length_pp=(320.0, 320.0),
        beam=(58.0, 580.0),
        draft=(20.8, 20.8),
    )
    # the tie (one value each) keeps the first cluster; force the majority the way SIMMAN does
    report["HULL"]["quantities"]["beam"]["family_values"]["SIMMAN-SERIES"] = 58.0
    report["HULL"]["quantities"]["beam"]["family_disagreement"] = {
        "SIMMAN-SERIES": {"chosen": 58.0, "dissent": [580.0, 580.0], "cluster_size": 2, "values": [58.0, 58.0, 580.0, 580.0]}
    }
    slips = adjudicate(report)["scale_typos"]
    assert len(slips) == 1
    assert slips[0]["printed"] == 580.0 and slips[0]["origin"] == "dissent"