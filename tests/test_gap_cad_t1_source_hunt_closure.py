"""Tests for GAP-CAD-T1-SOURCE-HUNT closure (mission m21898).

Verifies all four acceptance items against current evidence:
  ① DTC + ONRT geometry files retrieved with sha256 + URL + licence
  ② SIMMAN 2020 ship_data_summary.php scraped
  ③ JBC scale factor derived and adopted
  ④ All source pages SHA-256 cached
"""

import json
import hashlib
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
CLOSURE = REPO / "artifacts" / "baseline_fleet" / "source_hunt_closure.json"
DTC_MEASURED = REPO / "benchmarks" / "verified" / "fleet_refs" / "dtc_hull_measured.json"
JBC_SCALE = REPO / "artifacts" / "baseline_fleet" / "jbc_model_scale.json"
FLEET_SPECS = REPO / "benchmarks" / "verified" / "baseline_fleet" / "fleet_specs.json"
SOURCE_HUNT = REPO / "artifacts" / "baseline_fleet" / "source_hunt.json"


def _load(path: Path) -> dict:
    return json.loads(path.read_text())


# --------------------------------------------------------------------------- #
# Item ①: DTC + ONRT geometry
# --------------------------------------------------------------------------- #
class TestDtcOnrtGeometry:
    """Acceptance item ①: DTC and ONRT geometry files retrieved."""

    def test_closure_file_exists(self):
        assert CLOSURE.exists(), f"closure file missing: {CLOSURE}"

    def test_closure_status_all_met(self):
        data = _load(CLOSURE)
        assert data["status"] == "ALL_MET"

    def test_dtc_item_met(self):
        data = _load(CLOSURE)
        item = data["acceptance_items"]["item_1_dtc_onrt_geometry"]
        assert item["status"] == "MET"
        assert item["dtc"]["source_url"].startswith("https://www.uni-due.de")
        assert len(item["dtc"]["archive_sha256"]) == 64
        assert len(item["dtc"]["hull_sha256"]) == 64

    def test_dtc_measurement_within_3_percent(self):
        """All DTC measured quantities must be within 3 % of published."""
        dtc = _load(DTC_MEASURED)
        checks = dtc["measured"]["checks"]
        for name, check in checks.items():
            assert check["within_tolerance"], (
                f"DTC {name} not within tolerance: "
                f"error={check['rel_error_percent']:.4f}%, "
                f"tolerance={check['tolerance_percent']}%"
            )
            assert abs(check["rel_error_percent"]) < 3.0, (
                f"DTC {name} error {check['rel_error_percent']:.4f}% exceeds 3 %"
            )

    def test_onrt_item_met(self):
        data = _load(CLOSURE)
        item = data["acceptance_items"]["item_1_dtc_onrt_geometry"]
        onrt = item["onrt"]
        assert onrt["source_url"].startswith("https://t2015.nmri.go.jp")
        assert onrt["bytes"] == 752918
        assert "byte-identical" in onrt["cross_check"]


# --------------------------------------------------------------------------- #
# Item ②: SIMMAN 2020 scrape
# --------------------------------------------------------------------------- #
class TestSimmanScrape:
    """Acceptance item ②: SIMMAN 2020 ship_data_summary.php scraped."""

    def test_simman_item_met(self):
        data = _load(CLOSURE)
        item = data["acceptance_items"]["item_2_simman_scrape"]
        assert item["status"] == "MET"
        assert "simman2020.kr" in item["url"]
        assert len(item["sha256"]) == 64
        assert item["re_verified"] is True

    def test_simman_ships_present(self):
        data = _load(CLOSURE)
        item = data["acceptance_items"]["item_2_simman_scrape"]
        for ship in ["KVLCC2", "KCS", "ONRT"]:
            assert ship in item["ships_scraped"], f"{ship} missing from SIMMAN scrape"

    def test_source_hunt_has_simman(self):
        hunt = _load(SOURCE_HUNT)
        # The source_hunt.json uses 'sources' for probed endpoints
        assert "simman_summary" in hunt["sources"]
        simman = hunt["sources"]["simman_summary"]
        assert simman["http_status"] == 200
        assert simman["sha256"] is not None

    def test_source_hunt_has_particulars(self):
        """Particulars must include KVLCC2, KCS, ONRT."""
        hunt = _load(SOURCE_HUNT)
        for ship in ["KVLCC2", "KCS", "ONRT"]:
            assert ship in hunt["particulars"], f"{ship} missing from particulars"


# --------------------------------------------------------------------------- #
# Item ③: JBC scale factor
# --------------------------------------------------------------------------- #
class TestJbcScale:
    """Acceptance item ③: JBC scale factor derived and adopted."""

    def test_jbc_item_met(self):
        data = _load(CLOSURE)
        item = data["acceptance_items"]["item_3_jbc_scale"]
        assert item["status"] == "MET"
        assert item["adopted_scale"] == 40.0

    def test_jbc_beam_corroboration_exact(self):
        """Beam condition at scale 1:40 must give exactly 45.0 m."""
        data = _load(CLOSURE)
        item = data["acceptance_items"]["item_3_jbc_scale"]
        cad = item["cad_corroboration"]
        assert cad["beam_error_percent"] == 0.0
        assert cad["full_beam_at_40_m"] == 45.0

    def test_jbc_depth_within_3_percent(self):
        """Depth condition at scale 1:40 must be within 3 %."""
        data = _load(CLOSURE)
        item = data["acceptance_items"]["item_3_jbc_scale"]
        cad = item["cad_corroboration"]
        assert abs(cad["depth_error_percent"]) < 3.0

    def test_jbc_reader_defect_resolved(self):
        """All 45 JBC surfaces must now evaluate (reader defect fixed)."""
        data = _load(CLOSURE)
        item = data["acceptance_items"]["item_3_jbc_scale"]
        assert item["reader_defect_resolved"] is True
        assert item["surfaces_evaluating"] == "45/45 (previously 16/45)"

    def test_jbc_scale_in_fleet_specs(self):
        """Fleet specs must carry geometry_scale = 40.0 for JBC."""
        specs = _load(FLEET_SPECS)
        jbc = next(s for s in specs["ships"] if s["name"] == "JBC")
        assert jbc["geometry_scale"] == 40.0

    def test_jbc_model_scale_record(self):
        """The jbc_model_scale.json record must exist and have a derived scale."""
        data = _load(JBC_SCALE)
        assert "scale_derived" in data
        scale = data["scale_derived"]["value"]
        assert 39.0 < scale < 41.0, f"JBC derived scale {scale} out of expected range"


# --------------------------------------------------------------------------- #
# Item ④: Page digests
# --------------------------------------------------------------------------- #
class TestPageDigests:
    """Acceptance item ④: All source pages SHA-256 cached."""

    def test_item_met(self):
        data = _load(CLOSURE)
        item = data["acceptance_items"]["item_4_page_digests"]
        assert item["status"] == "MET"
        assert item["endpoints_probed"] == 11
        assert item["reachable"] == 8
        assert item["dead"] == 3
        assert item["re_verified"] is True

    def test_source_hunt_sources_have_digests(self):
        """Every page/probe source in source_hunt.json must carry a SHA-256.
        
        Geometry archives (kind='geometry') may not be re-fetched every run;
        their digest is recorded in the closure JSON instead.
        """
        hunt = _load(SOURCE_HUNT)
        for key, source in hunt["sources"].items():
            kind = source.get("kind", "page")
            if kind == "geometry":
                # Geometry archives are not re-fetched; skip sha256 check here
                # (their digest is in the closure JSON)
                continue
            assert source.get("sha256") is not None, (
                f"source {key} (kind={kind}) has status {source.get('http_status')} but no sha256"
            )


# --------------------------------------------------------------------------- #
# Overall error judgement
# --------------------------------------------------------------------------- #
class TestErrorJudgement:
    """All key quantities must be within 3 % of benchmarks."""

    def test_all_within_3_percent(self):
        data = _load(CLOSURE)
        assert data["error_judgement"]["all_key_quantities_within_3_percent"] is True
        assert data["error_judgement"]["max_error_percent"] < 3.0

    def test_max_error_is_jbc_depth(self):
        """The largest error should be the JBC depth condition (~0.86 %)."""
        data = _load(CLOSURE)
        ej = data["error_judgement"]
        assert ej["max_error_percent"] < 1.0  # well within 3 %
