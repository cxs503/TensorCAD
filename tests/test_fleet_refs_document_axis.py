"""Round 6: the document axis must be computed from the library, not asserted.

The anchor is judged on *lineages* (two documents printing one upstream table are one source).
That rule hid how many documents actually state a numeral: ``JBC`` read ``SINGLE_SOURCE`` while
three publishers print its particulars.  These tests pin the second axis:

* counts come from the voting cells only (a fallback-column reading never counts as a document);
* two documents that collapse to one lineage are ``LINEAGE_LIMITED``, never ``INDEPENDENT``;
* the summary the markdown renders is the summary the JSON carries.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.fleet.crossref import LINEAGE, QUANTITIES  # noqa: E402
from tools.baseline_fleet.refs_round6_documents import (  # noqa: E402
    NAMED,
    document_axis,
    gap_ledger,
    summary,
)

LIBRARY = ROOT / "benchmarks/verified/fleet_refs/crossref.json"
LEDGER = ROOT / "benchmarks/verified/fleet_refs/gap_ledger.json"


def _library() -> dict:
    return json.loads(LIBRARY.read_text())


def _cell(source: str, value: float, votes: bool = True) -> dict:
    return {
        "family": LINEAGE.get(source, source),
        "votes": votes,
        "cell": {"source": source, "value": value, "quantity": "cb", "raw": str(value)},
    }


def _synthetic(entries: list[dict], verdict: str = "SINGLE_SOURCE") -> dict:
    return {
        "report": {
            "SHIP": {
                "quantities": {
                    quantity: {"verdict": verdict, "entries": entries}
                    for quantity in QUANTITIES
                }
            }
        }
    }


def test_axis_is_present_for_every_ship_and_quantity() -> None:
    library = _library()
    axis = document_axis(library)
    assert set(axis) == set(library["report"])
    for ship, per_quantity in axis.items():
        assert set(per_quantity) == set(QUANTITIES), ship
        for quantity, row in per_quantity.items():
            assert row["documents"] == len(row["document_keys"]), (ship, quantity)
            assert row["documents"] <= len(row["entries"]) if "entries" in row else True
            assert row["publishers"] <= row["documents"]
            assert len(row["lineages"]) <= row["publishers"]
            assert row["lineage_limited"] == (row["verdict"] == "LINEAGE_LIMITED")


def test_two_documents_one_lineage_is_lineage_limited_not_independent() -> None:
    # JBC's two ingests: a journal paper and the workshop page, both printing the T2015 table.
    library = _synthetic([_cell("jbc_cfdletters", 0.858), _cell("t2015_jbc", 0.858)])
    axis = document_axis(library)["SHIP"]["cb"]
    assert axis["documents"] == 2
    assert axis["publishers"] == 2
    assert len(axis["lineages"]) == 1
    assert axis["verdict"] == "LINEAGE_LIMITED"
    assert axis["lineage_limited"] is True


def test_two_documents_two_lineages_is_independent() -> None:
    library = _synthetic(
        [_cell("t2015_kcs", 0.6505), _cell("simman2014_kcs", 0.6505)], verdict="CROSS_VERIFIED"
    )
    axis = document_axis(library)["SHIP"]["cb"]
    assert axis["documents"] == 2
    assert len(axis["lineages"]) == 2
    assert axis["verdict"] == "INDEPENDENT"


def test_non_voting_cell_is_not_a_document() -> None:
    library = _synthetic([_cell("t2015_jbc", 0.858), _cell("jbc_scse", 0.858, votes=False)])
    axis = document_axis(library)["SHIP"]["cb"]
    assert axis["documents"] == 1
    assert axis["document_keys"] == ["t2015_jbc"]
    assert axis["verdict"] == "SINGLE_DOCUMENT"


def test_single_or_missing_document_states() -> None:
    assert document_axis(_synthetic([_cell("t2015_jbc", 0.858)]))["SHIP"]["cb"]["verdict"] == (
        "SINGLE_DOCUMENT"
    )
    assert document_axis(_synthetic([]))["SHIP"]["cb"]["verdict"] == "MISSING"


def test_conflict_is_not_laundered_into_independence() -> None:
    library = _synthetic(
        [_cell("t2015_onrt", 8299.5), _cell("simman2020_summary", 8507.0)], verdict="CONFLICT"
    )
    axis = document_axis(library)["SHIP"]["cb"]
    assert axis["documents"] == 2
    assert axis["verdict"] == "CONFLICT"


def test_committed_library_axis_matches_recomputation() -> None:
    library = _library()
    stored = library["documents"]
    fresh = document_axis(library)
    assert stored["cells"] == json.loads(json.dumps(fresh))
    assert stored["summary"] == json.loads(json.dumps(summary(fresh)))


def test_gap_ledger_covers_exactly_the_cells_short_of_two_documents() -> None:
    library = _library()
    axis = document_axis(library)
    ledger = gap_ledger(axis)
    expected = {
        (ship, quantity)
        for ship, per_quantity in axis.items()
        for quantity in NAMED
        if per_quantity[quantity]["documents"] < 2
    }
    assert {(gap["ship"], gap["quantity"]) for gap in ledger["gaps"]} == expected
    # A probe is search provenance: nothing in the ledger may claim a document was pinned.
    for gap in ledger["gaps"]:
        assert gap["documents"] < 2
        if gap["probe"] is not None:
            assert gap["probe"]["outcome"]
            assert gap["probe"]["next_action"]
    assert LEDGER.exists()
    assert json.loads(LEDGER.read_text())["counts"]["gaps"] == len(ledger["gaps"])


def test_summary_never_reports_lineage_limited_as_independent() -> None:
    library = _library()
    axis = document_axis(library)
    rows = summary(axis)
    for ship, row in rows["ships"].items():
        assert row["independent_on_named"] + row["lineage_limited_on_named"] <= row[
            "documents_on_named"
        ], ship
    for ship in rows["ships_independent_on_all_named"]:
        assert rows["ships"][ship]["independent_on_named"] == len(NAMED)
    assert set(rows["ships_with_two_documents_on_all_named"]) >= set(
        rows["ships_independent_on_all_named"]
    )