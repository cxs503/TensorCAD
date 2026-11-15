"""Guards on the measured DTC hull geometry (anchor GAP-CAD-T1-SOURCE-HUNT).

Two independent layers:

* **record** -- the committed ``benchmarks/verified/fleet_refs/dtc_hull_measured.json`` must
  describe a hull whose *published reference* is the same one the fleet manifest carries, whose
  archive/member digests are 64 hex characters, and whose verdict is the logical consequence of
  its own check table (no verdict can be written that its numbers do not support);
* **recompute** -- when the source archive is still in this environment's cache, the hull is
  re-sampled at a coarse density and the volume/beam checks are recomputed from scratch, so a
  silent reader regression shows up as a red test instead of a stale green record.

The recompute is skipped (not passed) when the cache is absent: the bytes are publisher material
and are never committed.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:  # pragma: no cover
    sys.path.insert(0, str(REPO))

RECORD = REPO / "benchmarks" / "verified" / "fleet_refs" / "dtc_hull_measured.json"
ARCHIVE = Path("/root/.cache/hull_sources/dtc_igs_ordner.zip")
MANIFEST = REPO / "tensorcad" / "fleet" / "data" / "fleet_manifest.json"

sys.path.insert(0, str(REPO / "tools" / "baseline_fleet"))


def _record() -> dict:
    assert RECORD.is_file(), f"missing record {RECORD}"
    return json.loads(RECORD.read_text())


def test_reference_is_the_fleet_manifest_reference():
    """The published comparison values must be the ones the fleet manifest trusts."""
    record = _record()
    manifest = json.loads(MANIFEST.read_text())
    dtc = next(ship for ship in manifest["ships"] if ship["name"] == "DTC")
    reference = record["published_reference"]
    assert reference["length_pp"] == dtc["length"]
    assert reference["beam"] == dtc["beam"]
    assert reference["draft"] == dtc["draft"]
    assert reference["cb"] == dtc["block_coefficient"]
    assert reference["displacement_volume_m3"] == dtc["displacement_volume"]


def test_provenance_is_pinned():
    """sha256 + source URL + licence: the anchor's manifest requirement, per file."""
    record = _record()
    source = record["source"]
    assert source["archive_url"].startswith("https://")
    assert source["license"].strip(), "a licence statement is mandatory"
    assert len(source["archive_sha256"]) == 64
    assert len(source["hull_sha256"]) == 64
    assert source["hull_member"] == "DTC_hull.igs"
    for name, digest in source["member_sha256"].items():
        assert len(digest) == 64, name
    assert "DTC_hull.igs" in source["member_sha256"]


def test_verdict_follows_from_the_check_table():
    """A VERIFIED stamp is only legal when every check is inside its tolerance."""
    record = _record()
    checks = record["measured"]["checks"]
    assert set(checks) == {"beam", "block_coefficient", "displacement_volume"}
    failing = [name for name, item in checks.items() if not item["within_tolerance"]]
    for name, item in checks.items():
        assert item["within_tolerance"] == (
            abs(item["rel_error_percent"]) <= item["tolerance_percent"]
        ), name
    assert record["verdict"] == ("VERIFIED" if not failing else "WITHHELD")
    assert bool(record["verdict_reason"]) is True


def test_measured_volume_is_inside_the_fleet_volume_tolerance():
    """DTC is now a measured hull, not a withheld one: the volume gate must pass."""
    record = _record()
    volume = record["measured"]["checks"]["displacement_volume"]
    assert volume["within_tolerance"] is True
    assert abs(volume["rel_error_percent"]) < 1.0
    rows = record["measured"]["rows"]
    assert len(rows) >= 1
    # resolution sensitivity: the two densities must agree to well inside the tolerance
    if len(rows) >= 2:
        coarse, fine = rows[0], rows[-1]
        spread = abs(coarse["volume_rel_error_percent"] - fine["volume_rel_error_percent"])
        assert spread < 0.5, f"volume is not converged: {spread:.3f} % between densities"


@pytest.mark.skipif(not ARCHIVE.is_file(), reason="publisher archive not in this cache")
def test_recompute_from_cached_archive(tmp_path: Path):
    """Re-derive the hull from the cached bytes at a coarse density (no record trusted)."""
    import dtc_hull_measure as tool  # noqa: PLC0415  (tools dir on sys.path)

    reference = tool.published_reference(MANIFEST, RECORD.parent / "dtc_geometry.json")
    members = tool.extract_members(ARCHIVE)
    hull = tmp_path / tool.HULL_MEMBER
    hull.write_bytes(members[tool.HULL_MEMBER])
    measured = tool.measure(hull, reference, spacings=(2.0,), dx=1.0, dz=0.2)
    row = measured["rows"][0]
    assert measured["unit"]["unit_flag"] == 6  # metres
    # a point cloud, not the all-zero cloud the earlier round recorded
    assert row["point_cloud"]["distinct"] > 1000
    assert abs(row["point_cloud"]["extent"][1] - reference["beam"]) / reference["beam"] < 0.01
    assert abs(row["volume_rel_error_percent"]) < 1.0