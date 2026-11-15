"""Tests for the geometry-source ledger (anchor ``GAP-CAD-T1-SOURCE-HUNT``).

Two layers, deliberately separated:

* **pure-function tests** on synthetic input (a zip made in the test, a two-row HTML table,
  a stub IGES model) -- these prove the parser and the scale arithmetic, and they cannot
  pass by re-reading the cached publisher bytes;
* **record tests** on the committed ledger -- these assert the properties that make the
  record usable as evidence: every live source carries a 64-hex digest and a licence
  statement, dead ends keep their error-body digest, every check states its 3 % verdict
  explicitly, and the withheld items are withheld *with evidence* rather than dropped.

The end-to-end re-measurement against the cached bytes runs only when the cache is
present (``TENSORCAD_FLEET_GEOM_CACHE``), because the publisher archives are not
redistributed in this repository.
"""

from __future__ import annotations

import io
import json
import os
import re
import zipfile
from pathlib import Path

import pytest

from tensorcad.fleet import geometry_provenance as gp

REPO = Path(__file__).resolve().parents[1]
LEDGER = REPO / "tensorcad" / "fleet" / "data" / "geometry_source_ledger.json"
CACHE = Path(
    os.environ.get(
        "TENSORCAD_FLEET_GEOM_CACHE",
        "/root/tp-agent/work/m20473/artifacts/baseline_fleet/raw",
    )
)
HEX64 = re.compile(r"[0-9a-f]{64}")


@pytest.fixture(scope="module")
def ledger() -> dict:
    return json.loads(LEDGER.read_text())


# ------------------------------------------------------------------ pure functions
def test_sha256_is_the_documented_digest_of_the_bytes():
    assert gp.sha256_bytes(b"") == (
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    )
    assert gp.sha256_bytes(b"abc") == (
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    )


def test_kind_of_classifies_by_content_not_only_by_name():
    assert gp.kind_of("x.zip", b"PK\x03\x04rest") == "zip"
    assert gp.kind_of("x.bin", b"%PDF-1.5 ...") == "pdf"
    assert gp.kind_of("x.html", b"<!DOCTYPE HTML PUBLIC>") == "html"
    assert gp.kind_of("hull.igs", b"                                                                        S      1") == "iges"
    assert gp.kind_of("mystery", b"\x00\x01\x02") == "other"


def test_zip_census_and_member_selection_are_digest_bound():
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w") as archive:
        archive.writestr("spacer.igs", b"spacer")
        archive.writestr("hull.igs", b"the hull")
        archive.writestr("notes.pdf", b"%PDF-1.4 notes")
    payload = buffer.getvalue()

    census = gp.zip_members(payload)
    assert [row["name"] for row in census] == ["spacer.igs", "hull.igs", "notes.pdf"]
    by_name = {row["name"]: row for row in census}
    assert by_name["hull.igs"]["sha256"] == gp.sha256_bytes(b"the hull")
    assert by_name["notes.pdf"]["kind"] == "pdf"

    # Named-member selection matters: the first .igs in a publisher archive is not the hull.
    assert gp.zip_member(payload, "hull.igs")[1] == b"the hull"
    with pytest.raises(KeyError):
        gp.zip_member(payload, "absent.igs")


def test_axis_denominators_keep_every_axis_separate():
    extent = [7.542605, 1.130691, 0.630385]
    published = {
        "x": {"value": 280.0, "quantity": "Lpp"},
        "y": {"value": 45.0, "quantity": "B"},
    }
    denominators = gp.axis_denominators(extent, published)
    assert set(denominators) == {"x", "y"}
    assert denominators["y"]["denominator"] == pytest.approx(39.7987, rel=1e-4)
    # The x axis is reported even though it disagrees: a blundered axis must be visible.
    assert abs(denominators["x"]["denominator"] - 40.0) > 2.0


def test_classify_regime_separates_normalised_from_finite_length():
    normalised = gp.classify_regime([1.006323, 0.122081, 0.103128], 154.0)
    assert normalised["regime"] == "unit_length_normalised"
    finite = gp.classify_regime([7.542605, 1.130691, 0.630385], 280.0)
    assert finite["regime"] == "finite_length"
    assert finite["implied_denominator"] == pytest.approx(280.0 / 7.542605, rel=1e-6)


def test_check_row_states_the_tolerance_verdict():
    good = gp.check_row("B", 1.1307, 1.125, unit="m", source="test")
    assert good["within_tolerance"] is True and good["abs_rel_error"] < 0.03
    bad = gp.check_row("B", 1.5, 1.125, unit="m", source="test", note="explained")
    assert bad["within_tolerance"] is False and bad["note"] == "explained"


def test_scrapers_keep_the_row_text_as_evidence():
    markup = (
        "<html><body><table>"
        "<tr><td>Depth \\(D\\) (m)</td><td>25.000</td></tr>"
        "<tr><td>Length of waterline L WL (m)</td><td>3.147</td><td>154.0</td></tr>"
        "</table></body></html>"
    )
    depth = gp.scrape_scalar(markup, r"Depth \\?\(D\\?\) \\?\(m\\?\) ([\d.]+)")
    assert depth["value"] == 25.0 and "25.000" in depth["snippet"]
    pair = gp.scrape_pair(markup, r"Length of waterline L WL \\?\(m\\?\) ([\d.]+) ([\d.]+)")
    assert pair["values"] == [3.147, 154.0]
    with pytest.raises(ValueError):
        gp.scrape_scalar(markup, r"Breadth ([\d.]+)")


def test_quoted_scale_factor_is_recovered_from_the_crossref_cell():
    crossref = {
        "report": {
            "JBC": {
                "quantities": {
                    "lcb_percent": {
                        "entries": [
                            {"cell": {"label": "Block coefficient 0.858 -- Scale factor \u03bb 1 40 --"}}
                        ]
                    }
                }
            }
        }
    }
    hits = gp.find_published_quote_entries(crossref, "Scale factor")
    assert len(hits) == 1 and hits[0]["ship"] == "JBC"
    match = re.search(r"Scale factor [^0-9]{0,8}1\s+(\d+(?:\.\d+)?)", hits[0]["cell_text"])
    assert float(match.group(1)) == 40.0


def test_extent_uses_surface_evaluation_and_converts_units(monkeypatch):
    """A stub model with a known 2 x 2 grid: 1000 mm long, 250 mm wide, 100 mm tall."""

    class StubSurface:
        def sample(self, count_u, count_v):
            assert (count_u, count_v) == (5, 5)
            return [(0.0, 0.0, 0.0), (1000.0, 250.0, 100.0)]

    class StubModel:
        global_fields = (
            ",", ";", "", "f.igs", "system", "preprocessor", "32", "38", "6", "308", "15",
            "", "1.0D0", "2", "MM", "1", "0.254D0", "141031.065553", "0.0001D0",
            "1000.0D0", "", "NMRI", "10", "0", "141031.065553", "",
        )
        model_scale = 1.0
        metres_per_unit = 0.001
        units = "mm"

        def surfaces(self):
            return [StubSurface()]

        def stats(self):
            return {"128": 1}

    class EmptyModel(StubModel):
        def surfaces(self):
            return []

    monkeypatch.setattr(gp, "_iges_file", lambda payload: EmptyModel())
    with pytest.raises(ValueError, match="no B-spline surface"):
        gp.iges_extent_m(b"ignored")

    monkeypatch.setattr(gp, "_iges_file", lambda payload: StubModel())
    fingerprint = gp.iges_fingerprint(b"ignored")
    assert fingerprint["unit_flag"] == 2 and fingerprint["metres_per_unit"] == 0.001
    extent = gp.iges_extent_m(b"ignored", samples=5)
    assert extent["extent_m"] == [1.0, 0.25, 0.1]
    assert extent["header_max_coordinate_m"] == pytest.approx(1.0)
    assert extent["header_matches_extent"] is True


# ------------------------------------------------------------------ the committed record
def test_ledger_is_present_and_shaped(ledger):
    assert ledger["schema"] == "tensorcad.fleet.geometry_source_ledger.v1"
    assert ledger["anchor"] == "GAP-CAD-T1-SOURCE-HUNT"
    assert ledger["cache"]["artefacts"] == len(ledger["sources"]) >= 10
    assert ledger["network"]["state"], "the offline state of the run must be declared"


def test_every_live_source_carries_digest_url_and_licence(ledger):
    live = {k: v for k, v in ledger["sources"].items() if v["http_status"] == 200}
    assert len(live) >= 5
    for key, entry in live.items():
        assert HEX64.fullmatch(entry["sha256"] or ""), key
        assert entry["bytes"] > 0 and entry["url"].startswith("http"), key
        assert entry["licence"], f"{key}: licence statement not recorded"
        assert entry["retrieved_utc"], f"{key}: retrieval time not recorded"


def test_dead_ends_keep_their_error_body_digest(ledger):
    dead = ledger["dead_ends"]
    assert set(dead) >= {"dtc_unidue_page", "t2015_dtc_page"}
    for key, entry in dead.items():
        assert entry["http_status"] >= 400, key
        assert HEX64.fullmatch(entry["sha256"] or ""), f"{key}: error body not digested"
    assert ledger["sources"]["t2015_geometry_dir"]["http_status"] == 403


def test_geometry_archives_are_measured_not_merely_listed(ledger):
    assert set(ledger["geometry"]) == {"JBC", "ONRT"}
    for ship, block in ledger["geometry"].items():
        assert HEX64.fullmatch(block["archive_sha256"]), ship
        assert HEX64.fullmatch(block["hull_member"]["sha256"]), ship
        assert block["hull_member"]["bytes"] > 1_000_000, ship
        assert block["iges_fingerprint"]["unit_flag"] in (2, 6), ship
        assert len(block["extent"]["extent_m"]) == 3, ship
        assert block["extent"]["surfaces_measured"] > 10, ship
        # The header's own maximum-coordinate claim must agree with what we measured,
        # otherwise the extent is being taken from the wrong geometry.
        assert block["extent"]["header_matches_extent"] is True, ship
        assert len(block["members"]) >= 1, ship


def test_jbc_scale_denominator_reproduces_the_published_1_in_40(ledger):
    jbc = ledger["geometry"]["JBC"]
    quoted = ledger["scraped"]["jbc_scale_factor"]
    assert quoted["denominator"] == 40.0
    for axis in ("y", "z"):
        row = jbc["axis_denominators"][axis]
        assert abs(row["denominator"] - 40.0) / 40.0 < 0.03, (axis, row)
    # Lpp is not readable off the bounding box; the record says so instead of hiding it.
    assert jbc["axis_denominators"]["x"]["denominator"] < 38.5
    explained = {q["item"]: q for q in ledger["open_questions"]}
    assert explained["JBC x-axis extent vs Lpp"]["status"] == "EXPLAINED_NOT_MATCHED"


def test_onrt_archive_is_unit_length_normalised(ledger):
    onrt = ledger["geometry"]["ONRT"]
    assert onrt["regime"]["regime"] == "unit_length_normalised"
    assert ledger["scraped"]["onrt_lwl"]["values"] == [3.147, 154.0]
    # Cross-host identity: the same archive from the T2015 host and the SIMMAN host.
    assert ledger["cross_host"]["identical"] is True
    assert len(set(ledger["cross_host"]["hosts"].values())) == 1


def test_every_check_states_its_error_and_verdict(ledger):
    rows = ledger["checks"]
    assert len(rows) >= 6
    for row in rows:
        assert 0.0 <= row["abs_rel_error"] < 10.0, row["quantity"]
        assert isinstance(row["within_tolerance"], bool), row["quantity"]
        assert row["reference_source"], row["quantity"]
    passing = [row for row in rows if row["within_tolerance"]]
    assert len(passing) == ledger["checks_summary"]["within_tolerance"] >= 5
    # Every failing row must explain itself: an unexplained miss is a defect, not a finding.
    for row in rows:
        if not row["within_tolerance"]:
            assert row["note"], row["quantity"]
    assert ledger["acceptance"]["dtc_geometry"]["verdict"] == "WITHHELD"
    assert ledger["acceptance"]["dtc_geometry"]["evidence"], "withheld items need evidence"
    assert ledger["acceptance"]["jbc_scale_factor"]["verdict"] == "MET"


def test_ledger_does_not_embed_publisher_bytes(ledger):
    blob = LEDGER.read_text()
    assert "jbc.igs" in blob  # member names are recorded
    assert "<!DOCTYPE" not in blob and "%PDF" not in blob
    for member in ledger["geometry"]["JBC"]["members"]:
        assert HEX64.fullmatch(member["sha256"])
    assert LEDGER.stat().st_size < 200_000, "the ledger must stay a digest record, not a copy"


# ------------------------------------------------------------------ end to end
@pytest.mark.skipif(not CACHE.is_dir(), reason="publisher cache not present (bytes not in git)")
def test_verify_recomputes_every_digest_from_the_cache():
    from tools.baseline_fleet.geometry_source_ledger import main

    assert main(["--cache", str(CACHE), "--out", str(LEDGER), "--verify"]) == 0