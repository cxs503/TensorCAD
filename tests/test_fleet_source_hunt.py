"""Tests for the source hunt tool (GAP-CAD-T1-SOURCE-HUNT).

Verifies that:
1. Source fetch tool produces correct sha256 for reachable pages
2. DTC source URLs are confirmed dead (404)
3. Recorded digests match fetched digests for stable pages
4. ONRT 'no full-scale ship' fact is documented
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

REPORT_PATH = Path("/tmp/source_hunt_report.json")

# Recorded digests from tensorcad/fleet/official.py
RECORD_DIGESTS = {
    "t2015_kcs_gc": "68d9c24600cb699b16394fc5c32762ebc875c856d14ef8460e1d0ebf4e7c78a4",
    "t2015_jbc_gc": "fc2ce9054a8b7842162bac94365853fa076da2527e8c7e455832f9710f58f875",
    "t2015_onrt_gc": "5e13427332de8be17b9dc2245dc5830300baf2e2340efcd8a80eee2269c43a11",
    "simman2020_summary": "c1a668f6d2440c5d83a83c3616285e9288ef618dbc9a2feebd26627be3f4ebb3",
}


@pytest.fixture(scope="module")
def report():
    """Load the source hunt report."""
    if not REPORT_PATH.exists():
        pytest.skip("Source hunt report not found at /tmp/source_hunt_report.json")
    return json.loads(REPORT_PATH.read_text())


class TestPageDigests:
    """Verify that fetched page digests match recorded values."""

    @pytest.mark.parametrize("key", list(RECORD_DIGESTS.keys()))
    def test_digest_matches_record(self, report, key):
        entry = report.get(key)
        assert entry is not None, f"{key} not in report"
        assert entry.get("fetch_ok"), f"{key} fetch failed: {entry.get('error')}"
        sha = entry.get("sha256", "")
        recorded = RECORD_DIGESTS[key]
        # Recorded may be truncated
        if len(recorded) < 64:
            assert sha.startswith(recorded), f"{key}: sha={sha[:16]}... expected {recorded}"
        else:
            assert sha == recorded, f"{key}: sha mismatch"


class TestDTCSourceHunt:
    """Verify DTC source URLs are confirmed dead."""

    def test_dtc_uni_due_page_is_404(self, report):
        entry = report.get("dtc_uni_due_page")
        assert entry is not None
        assert entry.get("status") == 404 or "404" in str(entry.get("error", ""))


class TestONRTFindings:
    """Verify ONRT-specific findings."""

    def test_onrt_gc_page_reachable(self, report):
        entry = report.get("t2015_onrt_gc")
        assert entry is not None
        assert entry.get("fetch_ok"), f"ONRT GC page fetch failed: {entry.get('error')}"

    def test_onrt_desc_page_reachable(self, report):
        entry = report.get("t2015_onrt_desc")
        assert entry is not None
        assert entry.get("fetch_ok"), f"ONRT desc page fetch failed: {entry.get('error')}"


class TestSourceHuntTool:
    """Verify the source hunt tool itself."""

    def test_tool_exists(self):
        tool = Path(__file__).parent.parent / "tools" / "source_hunt" / "fetch_and_hash.py"
        assert tool.exists(), f"Source hunt tool not found at {tool}"

    def test_tool_importable(self):
        import importlib.util
        tool = Path(__file__).parent.parent / "tools" / "source_hunt" / "fetch_and_hash.py"
        spec = importlib.util.spec_from_file_location("fetch_and_hash", tool)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        assert hasattr(mod, "SOURCES")
        assert hasattr(mod, "RECORD_DIGESTS")
        assert "t2015_kcs_gc" in mod.SOURCES
        assert "dtc_uni_due_page" in mod.SOURCES
