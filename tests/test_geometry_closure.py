"""Offline checks on the geometry-closure record (anchor ``GAP-CAD-T1-SOURCE-HUNT``).

The record ``tensorcad/fleet/data/geometry_closure.json`` is committed; the IGES bytes are not.
These tests assert the properties that make it usable as evidence, without network access:

* every number that came off a publisher page is bound to that page's SHA-256 and keeps its row text;
* every error column is *recomputed* from the published/measured pair, so a claim cannot drift away
  from its own arithmetic, and the ``within_3pct`` flag must follow the standard;
* the hull verdicts are the ones the evidence supports: KCS closes on all four quantities, the JBC
  displacement does not close *and* does not converge (three grids recorded), ONRT is declared
  unmeasurable with the reason, not silently dropped;
* the six-hull reading of the anchor is falsified by recorded 404s with digested bodies;
* ``--verify-only`` runs green on the committed record.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
RECORD = REPO / "tensorcad" / "fleet" / "data" / "geometry_closure.json"
SPECS = REPO / "benchmarks" / "verified" / "baseline_fleet" / "fleet_specs.json"
HEX64 = re.compile(r"[0-9a-f]{64}")

pytestmark = pytest.mark.skipif(not RECORD.is_file(), reason="closure record not generated")


@pytest.fixture(scope="module")
def record() -> dict:
    return json.loads(RECORD.read_text())


def test_every_fetched_part_is_digested(record):
    for name, entry in record["ships"].items():
        for part in ("page", "archive", "iges"):
            blob = entry[part]
            assert HEX64.fullmatch(blob["sha256"]), f"{name}.{part}: no digest"
            assert blob["bytes"] > 0
        assert entry["licence"]


def test_published_values_are_bound_to_page_bytes(record):
    for name, entry in record["ships"].items():
        page_sha = entry["page"]["sha256"]
        for field, value in entry["published"].items():
            assert value["page_sha256"] == page_sha, f"{name}.{field}: wrong page digest"
            assert value["snippet"], f"{name}.{field}: no matched row text"
            assert value["value"] > 0.0


def test_draft_is_the_table_value_not_a_footnote_fragment(record):
    """The KCS draft row also carries `0.34` and `18`; picking the wrong one broke an earlier run."""
    assert record["ships"]["KCS"]["published"]["draft"]["value"] == pytest.approx(10.8)
    assert record["ships"]["JBC"]["published"]["draft"]["value"] == pytest.approx(16.5)
    assert record["ships"]["ONRT"]["published"]["draft"]["value"] == pytest.approx(5.494)


def test_errors_follow_from_the_recorded_numbers(record):
    checked = 0
    for name, entry in record["ships"].items():
        for row in entry["closure"]["rows"]:
            error = 100.0 * (row["measured"] / row["published"] - 1.0)
            assert row["error_percent"] == pytest.approx(error, abs=1e-6), f"{name}.{row['quantity']}"
            assert row["within_3pct"] == (abs(error) <= 3.0), f"{name}.{row['quantity']}"
            checked += 1
    assert checked >= 8


def test_kcs_closes_on_every_quantity(record):
    kcs = record["ships"]["KCS"]["closure"]
    assert kcs["verdict"] == "PASS", kcs.get("failed_quantities")
    assert {row["quantity"] for row in kcs["rows"]} == {"length_wl", "beam_wl", "volume", "cb"}


def test_jbc_displacement_is_declared_non_convergent(record):
    jbc = record["ships"]["JBC"]
    assert jbc["closure"]["verdict"] == "FAIL"
    assert set(jbc["closure"]["failed_quantities"]) == {"volume", "cb"}
    assert jbc["closure"]["reason"]  # the reason travels with the verdict
    errors = [entry["error_percent"] for entry in jbc["volume_sensitivity"]]
    assert len(errors) >= 3 and max(errors) - min(errors) > 20.0, errors


def test_onrt_is_declared_unmeasurable_with_a_reason(record):
    ornt = record["ships"]["ONRT"]
    assert ornt["status"] == "NOT_MEASURABLE"
    assert "disagree" in ornt["closure"]["reason"]
    assert all(not m["measurable"] for m in ornt["measurements"])
    # the digest still matches the committed fleet record: the file was fetched, not skipped
    assert ornt["fleet_specs_consistency"]["digest_matches"] is True
    assert ornt["census"]["surfaces_type128"] > 0


def test_fleet_specs_digest_agreement_is_checked_not_assumed(record):
    specs = json.loads(SPECS.read_text())
    by_name = {ship["name"]: ship for ship in specs["ships"]}
    for name in ("KCS", "JBC", "ONRT"):
        entry = record["ships"][name]
        assert entry["fleet_specs_consistency"]["recorded_iges_sha256"] == by_name[name]["iges_sha256"]
        assert entry["fleet_specs_consistency"]["digest_matches"] is True
    # ONRT's recorded model scale is *not* reproduced by the file and the record says so
    assert record["ships"]["ONRT"]["fleet_specs_consistency"]["consistent"] is False
    assert record["ships"]["KCS"]["fleet_specs_consistency"]["consistent"] is True
    assert record["ships"]["JBC"]["fleet_specs_consistency"]["consistent"] is True


def test_six_hull_reading_is_falsified_by_recorded_404s(record):
    probes = record["probes"]
    for key in ("t2015_dtmb5415_page", "t2015_kvlcc2_page", "t2015_dtc_page"):
        assert probes[key]["http_status"] == 404, key
        assert HEX64.fullmatch(probes[key]["sha256"]), f"{key}: error body not digested"
    assert probes["t2015_index"]["http_status"] == 200


def test_verify_only_mode_passes_on_the_committed_record(monkeypatch):
    monkeypatch.chdir(REPO)
    from tools.baseline_fleet.geometry_closure import main

    assert main(["--verify-only"]) == 0
