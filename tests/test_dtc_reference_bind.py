"""Checks on the DTC reference-binding record (anchor ``GAP-CAD-T1-SOURCE-HUNT``).

The record ``benchmarks/verified/fleet_refs/dtc_reference_bind.json`` replaces an
``offline_recollection`` reference (five particulars typed in from memory of El Moctar et al.
2012) with values bound to digested document bytes.  These tests are the offline arm: they
never touch the network, and they fail if the record stops being a chain.

Two halves:

* the positive case -- the committed record validates, every field is bound, and each bound
  value points at a 64-hex document digest plus the row text it was read out of;
* the negative controls -- a flipped digest, a BOUND field with no source, a deviation outside
  tolerance and a drifted reference value must each be reported.  A validator that cannot fail
  proves nothing.
"""

from __future__ import annotations

import copy
import importlib.util
import json
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
RECORD = REPO / "benchmarks" / "verified" / "fleet_refs" / "dtc_reference_bind.json"
MEASURED = REPO / "benchmarks" / "verified" / "fleet_refs" / "dtc_hull_measured.json"
FIGURE = REPO / "docs" / "figs" / "dtc_reference_bind.png"
TOOL = REPO / "tools" / "baseline_fleet" / "dtc_reference_bind.py"


def load_tool():
    spec = importlib.util.spec_from_file_location("dtc_reference_bind", TOOL)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def tool():
    return load_tool()


@pytest.fixture(scope="module")
def record() -> dict:
    return json.loads(RECORD.read_text())


def test_record_is_committed_with_the_expected_schema(record):
    assert record["schema"] == "tensorcad.fleet_refs.dtc_reference_bind/1.0"
    assert record["anchor"] == "GAP-CAD-T1-SOURCE-HUNT"
    assert record["verdict"] == "BOUND"


def test_committed_record_validates(tool):
    assert tool.verify_record(RECORD) == []


def test_every_field_is_bound_to_a_digested_document(record):
    assert set(record["fields"]) == {
        "length_pp", "beam", "draft", "displacement_volume_m3", "cb"}
    for name, field in record["fields"].items():
        assert field["status"] == "BOUND", name
        assert field["sources"], name
        for source in field["sources"]:
            assert len(source["sha256"]) == 64
            assert source["url"].startswith("https://")
            assert source["page"] >= 1
            assert source["token"] in source["row"], name
            assert abs(source["deviation_percent"]) <= field["tolerance_percent"]


def test_reference_values_are_the_ones_the_measurement_record_uses(record):
    measured = json.loads(MEASURED.read_text())
    reference = measured["published_reference"]
    for name, field in record["fields"].items():
        assert field["reference_value"] == pytest.approx(reference[name], rel=0, abs=0)


def test_record_retires_the_recollection_provenance(record):
    supersedes = record["supersedes"]
    assert supersedes["was"]["retrieval"] == "offline_recollection"
    assert supersedes["was"]["verified"] is False
    assert supersedes["now"]["verified"] is True
    assert supersedes["basis"]


def test_document_digest_is_the_one_in_the_record(record):
    document = record["documents"]["dogrul2018"]
    assert document["http_status"] == 200
    assert document["outcome"] == "particulars_block"
    assert len(document["sha256"]) == 64
    bound = record["fields"]["cb"]["sources"][0]
    assert bound["sha256"] == document["sha256"]
    # the row really contains the Cb pair the paper prints (model column, full-scale column)
    assert bound["row"].count("0.661") == 2


def test_archive_reproduces_the_committed_measurement(record):
    archive = record["archive"]
    measured = json.loads(MEASURED.read_text())
    assert archive["sha256"] == measured["source"]["archive_sha256"]
    assert archive["bytes"] == measured["source"]["archive_bytes"]
    assert archive["reproduces_committed"] is True


def test_registry_chain_matches_for_the_members_the_registry_carries(record):
    chain = record["registry_chain"]
    assert chain["all_registered_match"] is True
    registered = chain["members_registered"]
    assert registered, "no registered members to cross-check"
    for name in registered:
        assert chain["members"][name]["match"] is True
    # the propeller geometry ships in the archive but is not adopted by the fleet registry, so it
    # is digested and listed, not silently dropped
    for name in chain["members_unregistered"]:
        assert chain["members"][name]["in_registry"] is False
        # digested even where the registry has nothing to compare it with
        assert len(record["archive"]["members"][name]["sha256"]) == 64


def test_figure_exists_and_is_a_png():
    assert FIGURE.exists()
    assert FIGURE.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"


def _revalidate(tool, record: dict, tmp_path: Path) -> list[str]:
    path = tmp_path / "record.json"
    path.write_text(json.dumps(record))
    return tool.verify_record(path)


def test_flipped_document_digest_is_reported(tool, record, tmp_path):
    broken = copy.deepcopy(record)
    broken["fields"]["cb"]["sources"][0]["sha256"] = "0" * 64
    problems = _revalidate(tool, broken, tmp_path)
    assert any("does not match the document probe" in problem for problem in problems)


def test_malformed_document_digest_is_reported(tool, record, tmp_path):
    broken = copy.deepcopy(record)
    broken["fields"]["cb"]["sources"][0]["sha256"] = "not-a-digest"
    problems = _revalidate(tool, broken, tmp_path)
    assert any("malformed document digest" in problem for problem in problems)


def test_bound_field_without_source_is_reported(tool, record, tmp_path):
    broken = copy.deepcopy(record)
    broken["fields"]["beam"]["sources"] = []
    problems = _revalidate(tool, broken, tmp_path)
    assert any("BOUND without a source" in problem for problem in problems)


def test_deviation_outside_tolerance_is_reported(tool, record, tmp_path):
    broken = copy.deepcopy(record)
    broken["fields"]["draft"]["sources"][0]["deviation_percent"] = 3.0
    problems = _revalidate(tool, broken, tmp_path)
    assert any("deviates" in problem for problem in problems)


def test_drifted_reference_value_is_reported(tool, record, tmp_path):
    broken = copy.deepcopy(record)
    broken["fields"]["length_pp"]["reference_value"] = 340.0
    problems = _revalidate(tool, broken, tmp_path)
    assert any("drifted" in problem for problem in problems)


def test_archive_digest_must_match_the_measurement_record(tool, record, tmp_path):
    broken = copy.deepcopy(record)
    broken["archive"]["sha256"] = "1" * 64
    problems = _revalidate(tool, broken, tmp_path)
    assert any("archive digest" in problem for problem in problems)