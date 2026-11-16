"""The re-binding audit must fail loudly on every way a library digit can stop being true.

The invariant under test is the anchor's acceptance clause read literally: a published numeral
is only "verified" if it comes back out of the bytes whose digest it carries.  Each test below
plants one way that chain breaks and asserts the audit refuses to call it a pass.
"""

from __future__ import annotations

import json

import pytest

from tensorcad.fleet import rebind
from tensorcad.fleet.crossref import Cell


def cell(value: float, raw: str | None = None, basis: str = "volume", **extra) -> dict:
    payload = {
        "quantity": "lcb_percent",
        "value": value,
        "raw": raw if raw is not None else f"{value}",
        "label": "LCB",
        "source": "doc",
        "url": "https://example.invalid/doc",
        "sha256": "",
        "column": 1,
        "unit_basis": basis,
        "note": "",
        "derivation": "",
    }
    payload.update(extra)
    return payload


def test_published_digit_rebinds_from_its_own_bytes():
    stored = {"lcb_percent": cell(-1.68)}
    fresh = {"lcb_percent": Cell(**cell(-1.68))}
    verdicts = rebind.compare_cells("doc", "SHIP", stored, fresh)
    assert [v.status for v in verdicts] == [rebind.MATCH]
    assert verdicts[0].ok


def test_a_transcribed_digit_that_drifted_is_never_auto_repaired():
    stored = {"lcb_percent": cell(-1.70)}
    fresh = {"lcb_percent": Cell(**cell(-1.68))}
    verdicts = rebind.compare_cells("doc", "SHIP", stored, fresh)
    assert [v.status for v in verdicts] == [rebind.TRANSCRIPTION_DRIFT]
    assert rebind.repairable(verdicts) == []
    payload = {"report": {"SHIP": {"quantities": {"lcb_percent": {"entries": [
        {"cell": dict(stored["lcb_percent"])}
    ]}}}}}
    assert rebind.apply_repair(payload, verdicts) == []
    assert payload["report"]["SHIP"]["quantities"]["lcb_percent"]["entries"][0]["cell"]["value"] == -1.70


def test_stale_derived_arithmetic_is_repaired_in_place():
    stored = {"lcb_percent": cell(-1.704545, basis="derived", derivation="ratio 49")}
    fresh = {"lcb_percent": Cell(**cell(-1.6365, basis="derived", derivation="ratio 48.9355"))}
    verdicts = rebind.compare_cells("doc", "ONRT", stored, fresh)
    assert [v.status for v in verdicts] == [rebind.STALE_DERIVED]
    block = {
        "entries": [{"cell": dict(stored["lcb_percent"])}],
        "derived": {"cells": [{"family": "F", "value": -1.704545}], "vs_voters": {"reference": -1.68}},
    }
    payload = {"report": {"ONRT": {"quantities": {"lcb_percent": block}}}}
    done = rebind.apply_repair(payload, verdicts)
    assert [d["after"] for d in done] == [-1.6365]
    repaired = block["entries"][0]["cell"]
    assert repaired["value"] == -1.6365
    assert "48.9355" in repaired["derivation"]
    # the cached arithmetic block followed the cell, and the stale voter comparison is gone
    assert block["derived"]["cells"][0]["value"] == -1.6365
    assert "vs_voters" not in block["derived"]


def test_a_quantity_the_document_no_longer_yields_is_not_a_pass():
    stored = {"lcb_percent": cell(3.46875, basis="derived")}
    verdicts = rebind.compare_cells("doc", "KVLCC2", stored, {})
    assert [v.status for v in verdicts] == [rebind.NOT_IN_BYTES]
    summary = rebind.summarize([{"key": "doc", "status": rebind.DIGEST_MATCH}], verdicts)
    assert summary["all_cells_rebound"] is False
    assert summary["cells_rebound"] == 0


def test_a_numeral_the_library_does_not_carry_is_reported():
    verdicts = rebind.compare_cells("doc", "SHIP", {}, {"cb": Cell(**cell(0.65, basis="volume"))})
    assert [v.status for v in verdicts] == [rebind.NOT_IN_BYTES]
    assert verdicts[0].detail.startswith("the document states this quantity")


def test_uncached_document_is_never_counted_as_verified():
    row = rebind.digest_verdict("a" * 64, None)
    assert row["status"] == rebind.NOT_CACHED
    summary = rebind.summarize([{"key": "doc", "status": rebind.NOT_CACHED}], [])
    assert summary["documents_checked"] == 0
    assert summary["unchecked_sources"] == ["doc"]
    assert summary["all_cells_rebound"] is False


def test_request_counter_moves_the_bytes_but_not_the_document():
    """The SIMMAN-2014 pages emit a per-request counter; pinning raw bytes would always drift."""
    first = b"<p>LCB (%) 3.48</p><span>ML:768-PU:145-PLA:19</span>"
    second = b"<p>LCB (%) 3.48</p><span>ML:768-PU:146-PLA:19</span>"
    pinned = rebind.digest_verdict("0" * 64, first, second)
    assert pinned["cached"] != pinned["live"]  # the bytes really did move
    assert pinned["live_status"] == rebind.DIGEST_VOLATILE_ONLY
    assert "simman_request_counter" in pinned["live_masks"]


def test_real_content_change_is_still_drift():
    first = b"<p>LCB (%) 3.48</p><span>ML:768-PU:145-PLA:19</span>"
    second = b"<p>LCB (%) 3.51</p><span>ML:768-PU:146-PLA:19</span>"
    row = rebind.digest_verdict("0" * 64, first, second)
    assert row["live_status"] == rebind.DIGEST_DRIFT
    assert row["status"] == rebind.DIGEST_DRIFT


def test_pinned_bytes_that_still_hash_to_the_pin_are_a_match():
    payload = b"<table><tr><td>CB</td><td>0.6505</td></tr></table>"
    import hashlib

    pin = hashlib.sha256(payload).hexdigest()
    row = rebind.digest_verdict(pin, payload)
    assert row["status"] == rebind.DIGEST_MATCH
    assert row["cached"] == pin


def test_the_committed_library_is_internally_consistent():
    """Regression guard: the audited library must not go back to a stale derived digit."""
    from pathlib import Path

    path = Path("benchmarks/verified/fleet_refs/rebind_audit.json")
    if not path.exists():
        pytest.skip("audit artifact not built in this checkout")
    payload = json.loads(path.read_text())
    assert payload["schema"] == rebind.SCHEMA
    assert payload["counts"]["documents_by_status"].get(rebind.NOT_CACHED, 0) == 0
    onrt = [
        row
        for row in payload["cells"]
        if row["ship"] == "ONRT" and row["quantity"] == "lcb_percent" and row["status"] != rebind.MATCH
    ]
    assert onrt == [], f"ONRT LCB still not re-derived from its document: {onrt}"