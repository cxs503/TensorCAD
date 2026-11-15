"""The second-round (open-access PDF) merge path must run, stamp provenance and never invent.

Anchor ``REPRO-GEOM-FLEET-REFS``.  ``tools/baseline_fleet/crossref_oa.py build`` merges the
forty cells quoted out of eight open-access documents into the committed reference library.
Until this file existed the path had no test *and did not run at all*: three defects made
``build`` die before writing anything, so the eight documents never reached the library and
KVLCC2 / DTC stayed single-source although their second publisher was already on disk.

The tests assert invariants rather than executing the script:

* ``cells_from_report`` keeps the ``source`` field (it is both the extraction key and a
  required ``Cell`` field -- popping it raised ``TypeError``);
* a stored cell may carry keys added downstream (``closure_defect``) and still cast cleanly;
* every merged cell carries the 64-hex digest of the document it was quoted from;
* a ship with no readable source keeps a ``MISSING`` row (Series 60) instead of vanishing;
* merging is idempotent, so a re-run cannot accumulate duplicate evidence;
* the digests pinned in code are the digests measured at ingest.

The committed library is deliberately **not** asserted here: adopting the merged library needs
a matching adjudication pass (``adjudicate`` currently drops the KVLCC2 Bwl scale-typo flag
once the KVLCC2 beam row has two lineages), which is tracked as follow-up work in
``docs/fleet_refs_oa_merge.md``.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
LIBRARY = ROOT / "benchmarks/verified/fleet_refs/crossref.json"
INGEST = ROOT / "benchmarks/verified/fleet_refs/oa_ingest.json"

_TOOL = ROOT / "tools/baseline_fleet/crossref_oa.py"
_spec = importlib.util.spec_from_file_location("crossref_oa_tool", _TOOL)
assert _spec and _spec.loader
oa = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(oa)

HEX64 = re.compile(r"^[0-9a-f]{64}$")
CORE = ("displacement_volume", "cb", "length_pp", "beam", "draft")


@pytest.fixture(scope="module")
def library() -> dict:
    return json.loads(LIBRARY.read_text())


@pytest.fixture(scope="module")
def ingest() -> dict:
    return json.loads(INGEST.read_text())


@pytest.fixture(scope="module")
def merged(library: dict, ingest: dict) -> dict:
    return oa.merge_into(copy.deepcopy(library), ingest, None)


def test_cell_from_ignores_downstream_keys() -> None:
    """A merged cell may have been enriched by the adjudicator; the cast must not raise."""
    data = {
        "quantity": "cb",
        "value": 0.8098,
        "raw": "0.8098",
        "label": "CB [-]",
        "source": "kvlcc2_ugent",
        "url": "https://example.invalid/x.pdf",
        "sha256": "0" * 64,
        "column": 0,
        "unit_basis": "volume",
        "closure_defect": {"divide_by": 10.0},  # added by crossref.adjudicate
    }
    cell = oa.cell_from(data)
    assert cell.source == "kvlcc2_ugent" and cell.value == 0.8098


def test_cells_from_report_keeps_source(library: dict) -> None:
    """``source`` keys the extraction *and* is a required Cell field -- it must survive."""
    extractions = oa.cells_from_report(library["report"])
    assert extractions, "the committed library has cells"
    for ship, by_source in extractions.items():
        for key, by_quantity in by_source.items():
            assert key != "unattributed", f"{ship}: a cell lost its source key"
            for cell in by_quantity.values():
                assert cell.source == key


def test_merge_adds_the_second_round_cells(merged: dict, library: dict) -> None:
    """The merge really moves the coverage the anchor is measured on."""
    before, after = library["coverage"], merged["coverage"]
    assert merged["oa_cells_added"] == 52
    assert after["cross_verified"] == before["cross_verified"]  # committed library already carries round 2
    assert after["missing"] <= before["missing"]  # already merged in the committed library
    assert after["cells"] == before["cells"]  # same ship x quantity grid, no phantom rows


def test_merge_stamps_a_digest_on_every_quoted_cell(merged: dict, ingest: dict) -> None:
    digests = {row["sha256"] for row in ingest["sources"].values() if row["status"] == "read"}
    seen_oa = 0
    for ship, row in merged["report"].items():
        for quantity, cell in row["quantities"].items():
            for entry in cell["entries"]:
                digest = entry["cell"]["sha256"]
                assert HEX64.match(digest), f"{ship}/{quantity}: cell without a pinned digest"
                if entry["cell"]["source"] in oa.OA_KEYS:
                    assert digest in digests
                    seen_oa += 1
    assert seen_oa >= 40


def test_merge_is_idempotent(merged: dict, ingest: dict) -> None:
    twice = oa.merge_into(copy.deepcopy(merged), ingest, None)
    assert twice["coverage"] == merged["coverage"]
    assert twice["report"] == merged["report"]
    assert twice["lineage"] == merged["lineage"]


def test_sourceless_ship_is_kept_as_missing(merged: dict) -> None:
    """Series 60 has no readable source; the merge must still score it, not drop the row."""
    assert "Series60" in merged["report"]
    row = merged["report"]["Series60"]["quantities"]
    assert all(row[quantity]["verdict"] == "MISSING" for quantity in oa.QUANTITIES)


def test_second_round_promotes_kvlcc2_and_dtc(merged: dict) -> None:
    """Per ship x quantity verdicts the second round is supposed to change, and one it is not."""
    for ship in ("KVLCC2", "DTC"):
        row = merged["lineage"][ship]["displacement_volume"]
        assert row["verdict"].startswith("CROSS_VERIFIED"), (ship, row["verdict"])
        assert len(row["lineages"]) >= 2, (ship, row["lineages"])
    # a reproduction is still not a second source: JBC's three open-access PDFs restate
    # Tokyo-2015, so four documents collapse onto one lineage and the verdict stays single
    jbc = merged["lineage"]["JBC"]["cb"]
    assert jbc["verdict"] == "SINGLE_SOURCE"
    assert len(jbc["sources"]) > len(jbc["lineages"]) == 1


def test_pinned_digests_match_the_ingest_record(ingest: dict) -> None:
    """Code pins and the measured ingest record agree -- a drift here voids the provenance."""
    specs = {spec.key: spec for spec in oa.SOURCES}
    for key in oa.OA_KEYS:
        row = ingest["sources"][key]
        if row["status"] != "read":
            continue
        assert HEX64.match(specs[key].sha256), key
        assert specs[key].sha256 == row["sha256"], key


def test_oa_documents_declare_their_lineage() -> None:
    """Lineage is what stops three reproductions from voting as three opinions."""
    specs = {spec.key: spec for spec in oa.SOURCES}
    declared = {key: specs[key].lineage for key in oa.OA_KEYS}
    assert declared["kvlcc2_hrcak"] == "SIMMAN-SERIES"  # caption names SIMMAN 2008
    assert declared["jbc_cfdletters"] == "NMRI-T2015"
    assert declared["kvlcc2_ugent"] == "JSR-2017"  # own compilation, own lineage
    assert oa.lineage_of("kvlcc2_ugent") != oa.lineage_of("simman2020_summary")