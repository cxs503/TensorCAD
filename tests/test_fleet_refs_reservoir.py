"""The reservoir sweep reports leads, never admissions -- these tests pin that contract."""

from __future__ import annotations

import json
from pathlib import Path

from tensorcad.fleet.reservoir import (
    VERDICT_ABSENT,
    VERDICT_CANDIDATE,
    VERDICT_SOURCE_GAP,
    blocking_cells,
    classify,
    hulls_mentioned,
    printable_ratio,
    quantity_mentions,
    scan_document,
    summarise,
    text_digest,
)

#: The row shape that actually produced this round's one real lead: MuPDF drops the subscripted
#: label and the numeral lands on the following line.
JBC_TABLE = """
3.1 Test case Japan Bulk Carrier (JBC)
Table 1 Particulars of JBC model
Lpp [m] 4.8276
Lwl [m] 4.9138
B [m] 0.7759
d [m] 0.2845
Displacement [m3] 0.9155
Cb 0.8580
LCB (%Lpp), fwd +
2.5475
"""


def test_hull_names_are_detected():
    assert "JBC" in hulls_mentioned(JBC_TABLE)
    assert "KVLCC2" in hulls_mentioned("the MOERI KVLCC2 tanker")
    assert hulls_mentioned("a paper about propeller cavitation") == {}


def test_lcb_is_found_across_the_line_break():
    mentions = quantity_mentions(JBC_TABLE, ["lcb_percent"])
    assert [m["numeral"] for m in mentions["lcb_percent"]] == ["2.5475"]


def test_bibliography_lines_are_not_mentions():
    text = "Proc. CFD Workshop Tokyo, Vol. 1, pp. 350-358, Series 60 CB=0.60 model"
    assert "cb" not in quantity_mentions(text, ["cb"])


def test_implausible_numerals_are_rejected():
    assert "lcb_percent" not in quantity_mentions("LCB 2021.8 m", ["lcb_percent"])
    assert "cb" not in quantity_mentions("block coefficient 15.0", ["cb"])


def test_binary_and_duplicate_text_are_distinguished():
    assert printable_ratio("normal text\n") > 0.9
    assert printable_ratio("\x00\x01\x02binary") < 0.9
    assert text_digest("a b\n\nc") == text_digest("a  b c")


def test_blocking_cells_reads_the_nested_determinability_block():
    library = {
        "report": {
            "DTC": {"quantities": {"cb": {"verdict": "CROSS_VERIFIED"}, "lcb_percent": {"verdict": "MISSING"}}},
            "JBC": {"quantities": {"lcb_percent": {"verdict": "SINGLE_SOURCE"}}},
        },
        "determinability": {
            "ships": {
                "DTC": {"quantities": {"cb": {"class": "MULTI_LINEAGE_AGREE"}, "lcb_percent": {"class": "NO_DETERMINATION"}}},
                "JBC": {"quantities": {"lcb_percent": {"class": "SINGLE_LINEAGE_SINGLE_DOCUMENT"}}},
            }
        },
    }
    cells = {(c["hull"], c["quantity"]) for c in blocking_cells(library)}
    assert cells == {("DTC", "lcb_percent"), ("JBC", "lcb_percent")}


def test_classify_separates_candidate_from_confirmed_gap():
    document = scan_document(
        JBC_TABLE, name="overset.pdf", sha256="a" * 64, bytes_=10, in_library=False,
        quantities=["lcb_percent", "lwl", "beam"],
    )
    cells = classify(
        [
            {"hull": "JBC", "quantity": "lcb_percent", "library_verdict": "SINGLE_SOURCE", "class": "x"},
            {"hull": "JBC", "quantity": "lwl", "library_verdict": "MISSING", "class": "y"},
            {"hull": "ONRT", "quantity": "lcb_percent", "library_verdict": "MISSING", "class": "y"},
        ],
        [document],
    )
    verdicts = {(c["hull"], c["quantity"]): c["verdict"] for c in cells}
    assert verdicts[("JBC", "lcb_percent")] == VERDICT_CANDIDATE
    # "Lwl [m] 4.9138" is a numeral-bearing statement of the waterline length
    assert verdicts[("JBC", "lwl")] == VERDICT_CANDIDATE
    # a hull the document does not name has nothing to be a candidate from
    assert verdicts[("ONRT", "lcb_percent")] == VERDICT_ABSENT
    assert cells[0]["candidates"][0]["sha256"] == "a" * 64
    assert "JBC" in summarise([document], cells)["by_hull"]


def test_committed_artefact_keeps_its_policy_and_digests():
    path = Path("benchmarks/verified/fleet_refs/reservoir_sweep.json")
    payload = json.loads(path.read_text())
    assert payload["schema"] == "tensorcad.fleet_refs.reservoir_sweep/1.0"
    assert "never an admitted value" in payload["policy"]
    for document in payload["documents"]:
        assert len(document["sha256"]) == 64
    for cell in payload["cells"]:
        # a lead must carry the digest of the document it came from
        for hit in cell["candidates"]:
            assert len(hit["sha256"]) == 64
