"""Offline checks on the source-hunt record (anchor ``GAP-CAD-T1-SOURCE-HUNT``).

The record under ``artifacts/baseline_fleet/source_hunt.json`` is the only committed
copy of the hunt: the page bytes themselves stay in the git-ignored cache.  These
tests therefore assert the properties that make the record usable as evidence, and
they run **without the network** so that CI (and a reader) can re-check them:

* every 200 response carries a well-formed SHA-256, every failure carries its status;
* every extracted particular is bound to the digest of the page it came from and
  keeps the matched row text, so no digit is hand-typed;
* the ONRT archive served by the T2015 host is byte-identical to the SIMMAN copy;
* the page parser still works on markup it has never seen (a synthetic table), so a
  pass cannot come from the cached page simply being re-read.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
RECORD = REPO / "artifacts" / "baseline_fleet" / "source_hunt.json"
HEX64 = re.compile(r"[0-9a-f]{64}")

pytestmark = pytest.mark.skipif(not RECORD.is_file(), reason="source_hunt.json not generated")


@pytest.fixture(scope="module")
def record() -> dict:
    return json.loads(RECORD.read_text())


def test_every_reachable_source_carries_a_digest(record):
    reachable = {k: v for k, v in record["sources"].items() if v["http_status"] == 200}
    assert reachable, "no source was reachable when the record was written"
    for key, entry in reachable.items():
        assert HEX64.fullmatch(entry["sha256"] or ""), f"{key}: no SHA-256 recorded"
        assert entry["bytes"] > 0 and entry["url"].startswith("http")
        assert entry["licence"], f"{key}: licence state not recorded"


def test_dead_ends_are_recorded_not_dropped(record):
    """A 404 on the URL the DTC paper names is the evidence the anchor rests on."""
    for key in ("dtc_unidue_page", "dtc_unidue_section", "t2015_dtc_page"):
        entry = record["sources"][key]
        assert entry["http_status"] == 404, f"{key} unexpectedly reachable: {entry['http_status']}"
        assert HEX64.fullmatch(entry["sha256"] or ""), f"{key}: error body not digested"
    # The geometry directory is not listable, so a dtc/ subtree cannot be enumerated.
    assert record["sources"]["t2015_geometry_dir"]["http_status"] in (403, 404)


def test_particulars_are_bound_to_page_bytes(record):
    pages = record["sources"]
    assert record["particulars"], "no particulars extracted"
    for ship, fields in record["particulars"].items():
        for field, value in fields.items():
            page = pages[value["page"]]
            assert value["page_sha256"] == page["sha256"], f"{ship}.{field}: page digest mismatch"
            assert value["snippet"].strip(), f"{ship}.{field}: no matched row text"
            assert value["values"], f"{ship}.{field}: no number"


def test_published_principal_particulars(record):
    """Full-scale values first, as the SIMMAN tables list them."""
    p = record["particulars"]
    assert float(p["KCS"]["lpp_m"]["values"][0]) == pytest.approx(230.0)
    assert float(p["KVLCC2"]["lpp_m"]["values"][0]) == pytest.approx(320.0)
    assert float(p["KCS"]["displacement_m3"]["values"][0]) == pytest.approx(52030.0)
    assert float(p["KVLCC2"]["displacement_m3"]["values"][0]) == pytest.approx(312622.0)
    # ONRT full-scale displacement: the T2015 page gives 8507 t, SIMMAN 8507 -> one number.
    assert float(p["ONRT"]["displacement_m3"]["values"][0]) == pytest.approx(8507.0)
    assert float(p["ONRT_t2015"]["displacement"]["values"][1]) == pytest.approx(8507.0)


def test_model_scale_factors(record):
    """The scale row reads ``1 | N``: KCS/KVLCC2/ONRT model scales are published."""
    scales = {s: f["scale"]["values"] for s, f in record["particulars"].items() if "scale" in f}
    assert set(scales) == {"KCS", "KVLCC2", "ONRT"}
    for ship, values in scales.items():
        assert float(values[0]) == 1.0 and float(values[1]) > 1.0, ship
    # cross-checked against Lpp: full-scale 230 m over a 6.0702 m model is 1:37.89
    kcs = record["particulars"]["KCS"]
    assert float(kcs["scale"]["values"][1]) == pytest.approx(
        float(kcs["lpp_m"]["values"][0]) / float(kcs["lpp_m"]["values"][1]), rel=5e-3)


def test_onrt_geometry_identical_across_the_two_hosts(record):
    check = record["geometry_cross_check"]
    assert check, "no cross-check recorded"
    assert check["identical"] is True
    assert check["t2015_sha256"].startswith(check["simman_sha256"][:16])


def test_jbc_scale_factor_still_open(record):
    """Honest bookkeeping: the JBC geometry scale was *not* obtained in this round."""
    items = {q["item"]: q for q in record["open_questions"]}
    jbc = items["JBC model scale factor"]
    assert jbc["status"] == "WITHHELD"
    assert "model scale" in jbc["evidence"]


def test_parser_handles_unseen_markup():
    from tools.baseline_fleet.hunt_sources import parse_simman_tables

    markup = (
        "<html><body><table>"
        "<tr><td>Table 7 Main Particulars of TESTHULL</td></tr>"
        "<tr><td>Scale</td><td>1</td><td>40</td></tr>"
        "<tr><td>Lpp (m)</td><td>200</td><td>5.0</td></tr>"
        "<tr><td>Displacement (m3)</td><td>30000</td><td>0.46875</td></tr>"
        "</table></body></html>"
    )
    parsed = parse_simman_tables(markup)
    assert float(parsed["TESTHULL"]["lpp_m"]["values"][0]) == 200.0
    assert float(parsed["TESTHULL"]["scale"]["values"][1]) == 40.0
    assert parsed["TESTHULL"]["displacement_m3"]["snippet"]


def test_verify_only_is_green():
    from tools.baseline_fleet.hunt_sources import main

    assert main(["--out", "artifacts/baseline_fleet", "--verify-only"]) == 0