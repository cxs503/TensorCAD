"""Fleet baseline verification test (GAP-CAD-T1-BASELINE-FLEET, mission m21168).

This test runs the full fleet pipeline and asserts:

1. All 6 ships are admitted (watertight, manifold, Euler=2, volume < 1%).
2. All 5 negative controls are rejected.
3. Cross-route agreement (KCS: series60_form vs parametric_sections) < 1%.
4. Volume integration agreement (mesh vs hydrostatics) < 1e-6.
5. Artifact round-trip (write STL → read STL → re-gate) OK for all ships.
6. Source anchors complete for all ships.

The test also verifies that the per-ship manifest files contain all required
provenance fields and that the evidence.json schema is valid.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import pytest

from tensorcad.fleet import load_manifest, measure
from tensorcad.fleet.gates import VOLUME_TOLERANCE
from tensorcad.fleet.manifest import FleetManifest
from tensorcad.fleet.run import run_fleet

# Coarse resolution for the fixture (fast enough for CI)
COARSE = {"nx": 121, "nz": 25}
# Production resolution for specific tests (run separately)
PRODUCTION = {"nx": 481, "nz": 97}
RICHARDSON_TOLERANCE = 0.003  # 0.3%
ORDER_TOLERANCE = 0.5  # allow p from 0.5 to 1.5


@pytest.fixture(scope="module")
def fleet_evidence(tmp_path_factory):
    """Run the full fleet pipeline once for all tests (coarse resolution)."""
    out_dir = tmp_path_factory.mktemp("m21168_fleet")
    evidence = run_fleet(
        out_dir,
        nx=COARSE["nx"],
        nz=COARSE["nz"],
        volume_tolerance=0.05,  # relaxed tolerance for coarse grid
    )
    return evidence


@pytest.fixture(scope="module")
def manifest():
    return load_manifest()


# --------------------------------------------------------------------------- #
# Core acceptance tests
# --------------------------------------------------------------------------- #


def test_all_six_ships_admitted(fleet_evidence):
    """All 6 ships must be admitted (no refusals, no failures)."""
    summary = fleet_evidence["summary"]
    assert summary["ships_registered"] == 6
    assert summary["ships_admitted"] == 6
    assert summary["ships_refused_or_failed"] == 0


def test_all_negative_controls_rejected(fleet_evidence):
    """All 5 negative controls must be rejected."""
    assert fleet_evidence["summary"]["controls_all_rejected"] is True
    assert fleet_evidence["summary"]["negative_control_count"] == 5
    for control in fleet_evidence["negative_controls"]:
        assert control["rejected"], f"Control {control['control']} was not rejected"


def test_cross_route_agreement(fleet_evidence):
    """KCS built via series60_form and parametric_sections must agree < 1%."""
    assert fleet_evidence["summary"]["cross_route_agrees"] is True
    cross = fleet_evidence["cross_route_check"]
    assert cross["lofted_euler"] == 2
    assert cross["relative_difference"] < 0.01


def test_volume_integration_agreement(fleet_evidence):
    """Mesh volume and hydrostatics volume must agree < 1e-6 for all ships."""
    assert fleet_evidence["summary"]["volume_integration_all_agree"] is True
    for ship in fleet_evidence["ships"]:
        if ship.get("topology") is not None:
            assert ship["volume_integration_ok"], (
                f"{ship['ship']}: volume integration disagreement "
                f"{ship['volume_integration_agreement']}"
            )


def test_artifact_round_trip(fleet_evidence):
    """All 6 STL artifacts must re-ingest OK."""
    assert fleet_evidence["summary"]["artifacts_written"] == 6
    assert fleet_evidence["summary"]["artifacts_reingest_ok"] is True
    for ship in fleet_evidence["ships"]:
        artifact = ship.get("artifact")
        if artifact is not None:
            assert artifact["reingest_ok"], (
                f"{ship['ship']}: artifact re-ingest failed "
                f"(roundtrip error {artifact['roundtrip_relative_error']})"
            )


def test_source_anchors_complete(fleet_evidence):
    """All 6 ships must have complete source anchors."""
    assert fleet_evidence["summary"]["source_anchors_complete"] is True


# --------------------------------------------------------------------------- #
# Per-ship volume tests (coarse grid, relaxed tolerance)
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", ["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"])
def test_coarse_volume_within_five_percent(fleet_evidence, name):
    """Every ship's coarse-resolution mesh must be within 5% of its reference volume."""
    ship_block = next(s for s in fleet_evidence["ships"] if s["ship"] == name)
    assert ship_block["passed"], f"{name} failed gates: {ship_block.get('failed_gates')}"
    vol_err = ship_block["volume_relative_error"]
    assert vol_err < 0.05, (
        f"{name}: volume error {vol_err*100:.4f}% exceeds 5% at coarse resolution"
    )


# --------------------------------------------------------------------------- #
# Topology tests
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", ["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"])
def test_topology_gates(fleet_evidence, name):
    """Every admitted ship must pass all topology gates."""
    ship_block = next(s for s in fleet_evidence["ships"] if s["ship"] == name)
    topo = ship_block.get("topology")
    assert topo is not None, f"{name}: no topology (ship was refused)"
    assert topo["boundary_edge_count"] == 0, f"{name}: not watertight"
    assert topo["nonmanifold_edge_count"] == 0, f"{name}: not edge-manifold"
    assert topo["inconsistently_oriented_edge_count"] == 0, f"{name}: inconsistent orientation"
    assert topo["euler_characteristic"] == 2, f"{name}: Euler characteristic {topo['euler_characteristic']} != 2"
    assert topo["signed_volume"] > 0, f"{name}: signed volume {topo['signed_volume']} <= 0"
    assert topo["zero_area_face_count"] == 0, f"{name}: zero-area faces present"
    assert topo["duplicate_face_count"] == 0, f"{name}: duplicate faces present"


# --------------------------------------------------------------------------- #
# Manifest and provenance tests
# --------------------------------------------------------------------------- #


def test_manifest_has_six_ships_with_provenance(manifest):
    """The fleet manifest must register exactly 6 ships, each with complete provenance."""
    assert manifest.names() == ("KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60")
    for ship in manifest.ships:
        assert ship.provenance.complete, f"{ship.name}: incomplete provenance"
        assert ship.provenance.url.startswith("http"), f"{ship.name}: URL not http(s)"
        assert ship.provenance.license, f"{ship.name}: no license stated"
        assert ship.reference_volume() > 0, f"{ship.name}: reference volume <= 0"


def test_per_ship_manifests_written(fleet_evidence, tmp_path_factory):
    """Per-ship manifest JSON files must exist and contain required fields."""
    out_dir = tmp_path_factory.mktemp("m21168_fleet_check")
    run_fleet(out_dir, nx=61, nz=13, volume_tolerance=0.10)
    for name in ("KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"):
        manifest_path = out_dir / "manifests" / f"{name}.json"
        assert manifest_path.exists(), f"Missing manifest: {manifest_path}"
        data = json.loads(manifest_path.read_text())
        assert "ship" in data
        assert "record" in data
        assert "verdict" in data
        assert "topology" in data
        assert "provenance" in data["record"]


def test_ref_vs_got_csv_format(fleet_evidence, tmp_path_factory):
    """The ref_vs_got.csv must have the correct header and 6 data rows."""
    out_dir = tmp_path_factory.mktemp("m21168_fleet_csv")
    run_fleet(out_dir, nx=61, nz=13, volume_tolerance=0.10)
    csv_path = out_dir / "ref_vs_got.csv"
    assert csv_path.exists()
    rows = csv_path.read_text().strip().splitlines()
    assert len(rows) == 7  # header + 6 ships
    header = rows[0].split(",")
    assert "ship" in header
    assert "reference_volume_m3" in header
    assert "mesh_volume_m3" in header
    assert "volume_rel_error" in header


# --------------------------------------------------------------------------- #
# Evidence schema validation
# --------------------------------------------------------------------------- #


def test_evidence_schema(fleet_evidence):
    """The evidence.json must have the expected top-level keys."""
    required_keys = {
        "schema", "summary", "ships", "negative_controls",
        "cross_route_check", "settings", "provenance",
    }
    for key in required_keys:
        assert key in fleet_evidence, f"Missing key in evidence: {key}"
    assert fleet_evidence["schema"] == "tensorcad.fleet.run.v1"


# --------------------------------------------------------------------------- #
# Production resolution tests (separate, slower)
# --------------------------------------------------------------------------- #


@pytest.mark.slow
@pytest.mark.parametrize("name", ["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"])
def test_production_volume_within_one_percent(manifest, name):
    """Every ship's production-resolution mesh must be within 1% of its reference volume.

    This test is marked 'slow' because it builds each ship at 481×97 resolution.
    Run with: pytest -m slow tests/test_fleet_baseline_verification.py
    """
    from tensorcad.fleet.build import build_ship_geometry

    ship = manifest.by_name(name)
    mesh = build_ship_geometry(ship, nx=PRODUCTION["nx"], nz=PRODUCTION["nz"]).mesh
    facts = measure(mesh)
    error = abs(facts.signed_volume - ship.reference_volume()) / ship.reference_volume()
    assert error < VOLUME_TOLERANCE, (
        f"{name}: volume error {error*100:.4f}% exceeds {VOLUME_TOLERANCE*100}% at production resolution"
    )
    assert facts.boundary_edge_count == 0, f"{name}: not watertight at production resolution"
