"""Round-3 additions to the fleet reference library (anchor REPRO-GEOM-FLEET-REFS).

Two properties are pinned here, both of which a later round could quietly undo:

* the two documents added in this round stay digest-pinned and keep producing the cells they
  were admitted for (``kvlcc2_osaka`` for KVLCC2 ``lwl``/``lcb_percent``, ``dtc_ije`` for DTC);
* the library does **not** upgrade them to independent sources.  Two publishers restating the
  same case-definition table are one lineage; the test asserts the collapse, so a future edit
  that counts them as two cannot pass unnoticed.

The acceptance score is pinned the same way: criteria ②..⑤ pass, criterion ① fails, and the
set of ships that fail it is asserted literally with the reason it must stay that way.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from tensorcad.fleet.crossref import SOURCES, DIGEST_PATTERN, lineage_of

LIBRARY = Path("benchmarks/verified/fleet_refs/crossref.json")
INGEST = Path("benchmarks/verified/fleet_refs/oa_ingest.json")
ACCEPTANCE = Path("benchmarks/verified/fleet_refs/acceptance.json")
DOC = Path("docs/fleet_refs_round3.md")


@pytest.fixture(scope="module")
def library() -> dict:
    return json.loads(LIBRARY.read_text())


@pytest.fixture(scope="module")
def ingest() -> dict:
    return json.loads(INGEST.read_text())["sources"]


@pytest.fixture(scope="module")
def acceptance() -> dict:
    return json.loads(ACCEPTANCE.read_text())


def _spec(key: str):
    return next(s for s in SOURCES if s.key == key)


def test_round3_documents_are_registered_and_pinned(ingest):
    for key, expected in {
        "dtc_ije": {"displacement_volume", "cb", "length_pp", "beam", "draft"},
        "kvlcc2_osaka": {
            "displacement_volume",
            "cb",
            "length_pp",
            "lwl",
            "beam",
            "draft",
            "lcb_percent",
        },
    }.items():
        spec = _spec(key)
        assert DIGEST_PATTERN.match(spec.sha256), f"{key} has no pinned digest"
        assert spec.kind == "pdf" and spec.quotes, f"{key} carries no quote addresses"
        row = ingest[key]
        assert row["status"] == "read", row
        assert row["sha256"] == spec.sha256, f"{key}: ingest digest != pinned digest"
        assert row.get("digest_matches_pin") is not False
        assert expected <= set(row["quantities_read"]), row["quantities_read"]


def test_kvlcc2_lwl_and_lcb_are_transcription_corroborated_not_independent(library):
    for quantity, value in (("lwl", 325.5), ("lcb_percent", 3.48)):
        row = library["report"]["KVLCC2"]["quantities"][quantity]
        documents = {entry["cell"]["source"] for entry in row["entries"]}
        assert documents >= {"simman2014_kvlcc2", "kvlcc2_osaka"} or len(documents) >= 2
        assert any(entry["cell"]["value"] == value for entry in row["entries"])
        # ... and the lineage tally refuses to call that independence.
        lineage = library["lineage"]["KVLCC2"][quantity]
        assert lineage["verdict"] == "SINGLE_SOURCE", lineage
        assert lineage["lineages"] == ["SIMMAN-SERIES"], lineage
    assert lineage_of("kvlcc2_osaka") == lineage_of("simman2014_kvlcc2")


def test_dtc_gains_a_third_publisher_on_the_same_lineage(library):
    row = library["report"]["DTC"]["quantities"]["displacement_volume"]
    documents = {entry["cell"]["source"] for entry in row["entries"]}
    assert {"dtc_strath", "dtc_ije"} <= documents
    assert lineage_of("dtc_ije") == lineage_of("dtc_strath") == "DTC-PRIMARY"
    # DTC still publishes no LCB or Lwl here; the library must keep saying so.
    for quantity in ("lcb_percent", "lwl"):
        assert library["report"]["DTC"]["quantities"][quantity]["verdict"] == "MISSING"


def test_kvlcc2_breadth_typo_is_a_recorded_dissent_not_agreement(acceptance):
    dissent = acceptance["criteria"]["three_disputed_marked_with_both_values"]["dissent_recorded"]
    assert "KVLCC2/beam" in dissent
    values = dissent["KVLCC2/beam"]["SIMMAN-SERIES"]["values"]
    assert 580.0 in values and 58.0 in values


def test_acceptance_scores_and_the_honest_fail_are_pinned(acceptance):
    criteria = acceptance["criteria"]
    assert criteria["two_citation_and_digest"]["pass"] is True
    assert criteria["three_disputed_marked_with_both_values"]["pass"] is True
    assert criteria["four_specs_multi_source_arrays"]["pass"] is True
    assert criteria["five_kvlcc2_reference_wired"]["pass"] is True
    assert criteria["five_kvlcc2_reference_wired"]["kvlcc2_published_volume_m3"] == 312622.0
    # Criterion ① is an honest FAIL: no reachable source gives JBC or Series 60 a second
    # lineage, and KVLCC2's LCB/Lwl are single-lineage.  Do not "fix" this by reclassifying.
    one = criteria["one_two_independent_lineages"]
    assert one["pass"] is False
    assert {ship for ship, row in one["ships"].items() if not row["pass"]} == {
        "DTC",
        "JBC",
        "KVLCC2",
        "ONRT",
        "Series60",
    }
    assert one["ships"]["KCS"]["pass"] is True


def test_document_is_generated_from_the_library():
    text = DOC.read_text()
    library = json.loads(LIBRARY.read_text())
    assert f"cells {library['coverage']['cells']}" in text
    assert "FAIL" in text and "kvlcc2_osaka" in text


def test_ingest_keeps_rows_of_uncached_documents(tmp_path):
    """A document fetched in an earlier round is no longer cached; its row must survive."""
    spec = importlib.util.spec_from_file_location(
        "crossref_oa", Path("tools/baseline_fleet/crossref_oa.py")
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    out = tmp_path / "oa_ingest.json"
    out.write_text(
        json.dumps(
            {
                "schema": "tensorcad.fleet_refs.oa_ingest/1.0",
                "sources": {
                    "dtc_strath": {
                        "status": "read",
                        "sha256": "0" * 64,
                        "quantities_read": ["cb"],
                        "cells": {},
                    }
                },
            }
        )
    )
    empty = tmp_path / "raw"
    empty.mkdir()
    module.main(["ingest", "--raw", str(empty), "--out", str(out), "--keep-existing"])
    kept = json.loads(out.read_text())["sources"]["dtc_strath"]
    assert kept["status"] == "read"
    assert kept["sha256"] == "0" * 64
    assert "reused from the previous ingest" in kept["note"]

    module.main(["ingest", "--raw", str(empty), "--out", str(out)])
    assert json.loads(out.read_text())["sources"]["dtc_strath"]["status"] == "missing_local_file"