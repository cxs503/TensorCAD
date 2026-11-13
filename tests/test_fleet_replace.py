"""Tests for the placeholder-replacement pipeline (GAP-CAD-T1-BASELINE-FLEET, round 2).

The pipeline is the part of the anchor a gate alone cannot deliver: a gate says "this
file is a placeholder", this module has to put a gated hull in its place and prove the
file it delivered is the one that was gated. The tests therefore assert invariants, not
just "the function returned":

* the source library fixture carries both reported defect classes at their reported
  sizes (6 placeholders, 1 byte-identical pair) — otherwise the fixture proves nothing;
* nothing is silently dropped: every scanned file appears in the ledger exactly once and
  in the delivery directory, either replaced or withheld with a reason;
* every delivered file is *re-read from disk* and independently measured: watertight,
  manifold, consistently oriented, Euler 2, and within the volume tolerance of its
  reference;
* the delivered digest differs from the placeholder digest it replaces;
* re-auditing the delivery directory finds no placeholder left, and the files the
  pipeline could not fix (no manifest record) are still flagged rather than forgotten.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tensorcad.fleet.build import read_stl, write_binary_stl
from tensorcad.fleet.controls import defect_geometry_library, placeholder_box
from tensorcad.fleet.gates import measure, placeholder_reasons
from tensorcad.fleet.manifest import load_manifest, sha256_file
from tensorcad.fleet.replace import (
    DUPLICATE,
    PLACEHOLDER,
    VALID,
    audit_directory,
    plan_replacements,
    replace_library,
)

# Coarse grids keep the suite fast; the delivered production run is gated at the
# mission's 1 % at nx = 481, nz = 97 (docs/gap_cad_t1_baseline_fleet_r2.md).
NX, NZ = 121, 49
TOL = 0.02


@pytest.fixture(scope="module")
def library(tmp_path_factory) -> Path:
    root = tmp_path_factory.mktemp("library")
    defect_geometry_library(root)
    return root


@pytest.fixture(scope="module")
def delivery(library, tmp_path_factory) -> tuple[Path, dict]:
    out = tmp_path_factory.mktemp("delivery")
    ledger = replace_library(library, out, nx=NX, nz=NZ, volume_tolerance=TOL)
    return out, ledger


def test_fixture_reproduces_both_reported_defects(library: Path) -> None:
    reports = audit_directory(library)
    placeholders = [r for r in reports if r.classification == PLACEHOLDER]
    duplicates = [r for r in reports if r.classification == DUPLICATE]
    assert len(reports) == 8
    assert len(placeholders) == 6, [r.name for r in placeholders]
    assert len(duplicates) == 2, [r.name for r in duplicates]
    assert {r.name for r in duplicates} == {"KVLCC1.stl", "DTMB5415.stl"}
    assert len({r.sha256 for r in duplicates}) == 1
    for report in placeholders:
        assert report.triangle_count == 304
        assert any("placeholder_signature" in reason for reason in report.reasons)


def test_audit_flags_identical_bytes_across_different_ships(library: Path) -> None:
    reports = {report.name: report for report in audit_directory(library)}
    assert reports["KVLCC1.stl"].sha256 == reports["DTMB5415.stl"].sha256
    assert any("byte_identical_geometry" in reason for reason in reports["KVLCC1.stl"].reasons)


def test_plan_pairs_defects_with_manifest_records_and_withholds_the_rest(library: Path) -> None:
    reports = audit_directory(library)
    manifest = load_manifest()
    actionable, withheld = plan_replacements(reports, manifest)
    assert {entry["ship"] for entry in actionable} == {
        "KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"
    }
    assert {entry["stem"] for entry in withheld} == {"KVLCC1", "DTMB5415"}
    for entry in withheld:
        assert "no fleet-manifest record" in entry["reason"]
        assert entry["status"] == "WITHHELD"


def test_every_scanned_file_is_accounted_for_once(library: Path, delivery) -> None:
    out, ledger = delivery
    scanned = sorted(report.name for report in audit_directory(library))
    assert sorted(block["file"] for block in ledger["files"]) == scanned
    assert ledger["summary"]["files_scanned"] == 8
    assert ledger["summary"]["withheld"] == 2
    assert len(ledger["summary"]["replaced"]) == 6
    # a drop-in directory: same file set, six of them replaced in place
    assert sorted(path.name for path in out.glob("*.stl")) == scanned


def test_replacements_are_gated_watertight_manifold_and_in_tolerance(delivery) -> None:
    out, ledger = delivery
    replacements = [block for block in ledger["files"] if block["status"] == "REPLACED"]
    assert len(replacements) == 6
    for block in replacements:
        path = out / block["replacement_file"]
        assert path.exists()
        assert sha256_file(path) == block["replacement_sha256"]
        assert block["replacement_sha256"] != block["sha256"], "placeholder bytes survived"
        assert block["sha256_changed"] is True
        mesh = read_stl(path)
        facts = measure(mesh)
        assert facts.boundary_edge_count == 0
        assert facts.nonmanifold_edge_count == 0
        assert facts.inconsistently_oriented_edge_count == 0
        assert facts.euler_characteristic == 2
        assert facts.triangle_count > 1000
        assert placeholder_reasons(facts, mesh=mesh) == []
        assert block["volume_relative_error"] == pytest.approx(
            abs(facts.signed_volume - block["volume_reference"]) / block["volume_reference"]
        )
        assert block["volume_relative_error"] < TOL, block["verification_conclusion"]
        assert block["reingest_ok"] is True
        assert "REPLACED" in block["verification_conclusion"]


def test_delivered_replacements_do_not_share_bytes(delivery) -> None:
    out, ledger = delivery
    digests = [block["replacement_sha256"] for block in ledger["files"] if block["status"] == "REPLACED"]
    assert len(set(digests)) == len(digests), "the same geometry was delivered for two ships"
    # the delivered KCS hull must not equal the bytes of the hull that used to be shared
    shared = next(block["sha256"] for block in ledger["files"] if block["file"] == "KVLCC1.stl")
    assert shared not in digests


def test_delivery_is_a_fixpoint(delivery) -> None:
    """Re-auditing the delivery finds no placeholder; unregistered ships stay flagged."""
    out, _ = delivery
    reports = {report.name: report for report in audit_directory(out)}
    assert len(reports) == 8
    assert [name for name, report in reports.items() if report.classification == PLACEHOLDER] == []
    still_flagged = sorted(name for name, report in reports.items() if report.classification != VALID)
    assert still_flagged == ["DTMB5415.stl", "KVLCC1.stl"], still_flagged


def test_stl_roundtrip_loses_only_float32_precision(tmp_path: Path) -> None:
    mesh = placeholder_box(304, length=100.0, beam=20.0, height=8.0)
    path = tmp_path / "box.stl"
    written_sha = write_binary_stl(mesh, path)
    assert written_sha == sha256_file(path)
    facts_in = measure(mesh)
    facts_out = measure(read_stl(path))
    assert facts_out.triangle_count == facts_in.triangle_count
    # the control box is deliberately not watertight, so the invariant is that writer and
    # reader agree on its topology, not that the topology is good
    assert facts_out.vertex_count == facts_in.vertex_count
    assert facts_out.boundary_edge_count == facts_in.boundary_edge_count
    assert abs(facts_out.signed_volume - facts_in.signed_volume) / abs(facts_in.signed_volume) < 1e-6


def test_ledger_is_json_serialisable_and_records_settings(delivery) -> None:
    out, _ = delivery
    restored = json.loads((out / "ledger.json").read_text())
    assert restored["schema"] == "tensorcad.fleet.replace.v1"
    assert restored["anchor"] == "GAP-CAD-T1-BASELINE-FLEET"
    assert restored["settings"]["manifest_sha256"] == load_manifest().sha256
    assert restored["settings"]["source_registry_sha256"]
    assert restored["summary"]["all_replacements_gated"] is True
    assert restored["summary"]["worst_replacement_volume_error"] < TOL