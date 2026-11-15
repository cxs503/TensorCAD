"""Acceptance tests for GAP-CAD-T1-BASELINE-FLEET (mission m21232).

These tests directly assert the mission acceptance criteria:
1. All 6 ships watertight + manifold
2. Volume vs published displacement error < 1%
3. Each ship has manifest with provenance
4. Placeholder STL replacement plan verified
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from tensorcad.fleet import (
    build_ship_geometry,
    evaluate_gates,
    load_manifest,
    measure,
    placeholder_box,
    placeholder_reasons,
)
from tensorcad.fleet.sources import load_registry
from tensorcad.fleet.run import run_fleet

# Coarse resolution for fast tests; production resolution tested separately
COARSE = {"nx": 121, "nz": 25}


@pytest.fixture(scope="module")
def manifest():
    return load_manifest()


@pytest.fixture(scope="module")
def registry():
    return load_registry()


@pytest.fixture(scope="module")
def coarse_geometries(manifest):
    return {ship.name: build_ship_geometry(ship, **COARSE) for ship in manifest.ships}


# --------------------------------------------------------------------------- #
# Acceptance Criterion 1: All 6 ships watertight + manifold
# --------------------------------------------------------------------------- #
class TestWatertightManifold:
    """Every registered hull must be watertight, edge-manifold, single-shell (Euler 2)."""

    @pytest.mark.parametrize("name", ["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"])
    def test_watertight(self, manifest, coarse_geometries, name):
        mesh = coarse_geometries[name].mesh
        facts = measure(mesh)
        assert facts.boundary_edge_count == 0, f"{name}: {facts.boundary_edge_count} boundary edges"

    @pytest.mark.parametrize("name", ["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"])
    def test_edge_manifold(self, manifest, coarse_geometries, name):
        mesh = coarse_geometries[name].mesh
        facts = measure(mesh)
        assert facts.nonmanifold_edge_count == 0, f"{name}: {facts.nonmanifold_edge_count} non-manifold edges"

    @pytest.mark.parametrize("name", ["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"])
    def test_consistently_oriented(self, manifest, coarse_geometries, name):
        mesh = coarse_geometries[name].mesh
        facts = measure(mesh)
        assert facts.inconsistently_oriented_edge_count == 0, f"{name}: inconsistent orientation"

    @pytest.mark.parametrize("name", ["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"])
    def test_euler_characteristic(self, manifest, coarse_geometries, name):
        mesh = coarse_geometries[name].mesh
        facts = measure(mesh)
        assert facts.euler_characteristic == 2, f"{name}: Euler = {facts.euler_characteristic}"

    @pytest.mark.parametrize("name", ["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"])
    def test_positive_signed_volume(self, manifest, coarse_geometries, name):
        mesh = coarse_geometries[name].mesh
        facts = measure(mesh)
        assert facts.signed_volume > 0.0, f"{name}: signed volume = {facts.signed_volume}"

    @pytest.mark.parametrize("name", ["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"])
    def test_no_zero_area_faces(self, manifest, coarse_geometries, name):
        mesh = coarse_geometries[name].mesh
        facts = measure(mesh)
        assert facts.zero_area_face_count == 0, f"{name}: {facts.zero_area_face_count} zero-area faces"

    @pytest.mark.parametrize("name", ["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"])
    def test_no_duplicate_faces(self, manifest, coarse_geometries, name):
        mesh = coarse_geometries[name].mesh
        facts = measure(mesh)
        assert facts.duplicate_face_count == 0, f"{name}: {facts.duplicate_face_count} duplicate faces"


# --------------------------------------------------------------------------- #
# Acceptance Criterion 2: Volume error < 1% (at production resolution)
# --------------------------------------------------------------------------- #
class TestVolumeAccuracy:
    """Volume vs published displacement error must be < 1%."""

    VOLUME_TOLERANCE = 0.01

    @pytest.mark.parametrize("name", ["KCS", "KVLCC2", "ONRT", "JBC", "DTC", "Series60"])
    def test_volume_within_one_percent_production(self, manifest, name):
        """Production resolution (481x97) must meet 1% criterion."""
        ship = manifest.by_name(name)
        mesh = build_ship_geometry(ship, nx=481, nz=97).mesh
        facts = measure(mesh)
        ref = ship.reference_volume()
        error = abs(facts.signed_volume - ref) / ref
        assert error < self.VOLUME_TOLERANCE, (
            f"{name}: volume error {error*100:.4f}% exceeds 1% tolerance "
            f"(V_ref={ref:.1f}, V_mesh={facts.signed_volume:.1f})"
        )

    def test_worst_volume_error_below_tolerance(self, manifest):
        """The worst volume error across all ships must be below 1%."""
        worst_error = 0.0
        worst_ship = None
        for ship in manifest.ships:
            mesh = build_ship_geometry(ship, nx=481, nz=97).mesh
            facts = measure(mesh)
            ref = ship.reference_volume()
            error = abs(facts.signed_volume - ref) / ref
            if error > worst_error:
                worst_error = error
                worst_ship = ship.name
        assert worst_error < self.VOLUME_TOLERANCE, (
            f"Worst volume error: {worst_ship} @ {worst_error*100:.4f}%"
        )


# --------------------------------------------------------------------------- #
# Acceptance Criterion 3: Each ship has manifest with provenance
# --------------------------------------------------------------------------- #
class TestProvenance:
    """Every ship must have complete provenance (source URL / license / SHA / verification)."""

    def test_all_ships_have_complete_provenance(self, manifest):
        for ship in manifest.ships:
            assert ship.provenance.complete, f"{ship.name}: incomplete provenance"
            assert ship.provenance.url.startswith("http"), f"{ship.name}: invalid URL"
            assert ship.provenance.license, f"{ship.name}: no license stated"
            assert ship.provenance.source, f"{ship.name}: no source stated"

    def test_source_registry_covers_all_ships(self, manifest, registry):
        missing = registry.missing_for(ship.name for ship in manifest.ships)
        assert not missing, f"Source registry missing anchors for: {missing}"

    def test_manifest_sha256_exists(self, manifest):
        assert manifest.sha256 is not None
        assert len(manifest.sha256) == 64

    def test_run_fleet_writes_per_ship_manifests(self, tmp_path, manifest):
        evidence = run_fleet(tmp_path, nx=61, nz=13, volume_tolerance=0.10, write=True)
        for ship in manifest.ships:
            manifest_path = tmp_path / "manifests" / f"{ship.name}.json"
            assert manifest_path.exists(), f"Missing manifest for {ship.name}"
            data = json.loads(manifest_path.read_text())
            assert "provenance" in data or "record" in data


# --------------------------------------------------------------------------- #
# Acceptance Criterion 4: Placeholder STL replacement verified
# --------------------------------------------------------------------------- #
class TestPlaceholderReplacement:
    """The 304-triangle placeholder must be detected and rejected; replacement meshes pass."""

    def test_placeholder_box_detected(self):
        mesh = placeholder_box(304)
        facts = measure(mesh)
        reasons = placeholder_reasons(facts, mesh=mesh)
        assert any("placeholder_signature" in r for r in reasons), reasons

    def test_placeholder_rejected_by_gate(self, manifest):
        mesh = placeholder_box(304)
        ship = manifest.by_name("KCS")
        verdict = evaluate_gates(ship, mesh, manifest)
        assert "G1_not_placeholder" in verdict["failed_gates"], verdict["failed_gates"]

    def test_replacement_mesh_passes_placeholder_gate(self, manifest, coarse_geometries):
        for ship in manifest.ships:
            mesh = coarse_geometries[ship.name].mesh
            facts = measure(mesh)
            reasons = placeholder_reasons(facts, mesh=mesh)
            assert not reasons, f"{ship.name}: replacement mesh flagged as placeholder: {reasons}"

    def test_replacement_mesh_has_realistic_triangle_count(self, manifest, coarse_geometries):
        for ship in manifest.ships:
            mesh = coarse_geometries[ship.name].mesh
            facts = measure(mesh)
            assert facts.triangle_count > 10_000, (
                f"{ship.name}: only {facts.triangle_count} triangles, expected > 10k at coarse resolution"
            )

    def test_negative_controls_all_rejected(self, tmp_path, manifest):
        evidence = run_fleet(tmp_path, nx=61, nz=13, volume_tolerance=0.10, write=True)
        assert evidence["summary"]["controls_all_rejected"] is True
        assert evidence["summary"]["negative_control_count"] >= 5


# --------------------------------------------------------------------------- #
# Integration: full fleet run
# --------------------------------------------------------------------------- #
class TestFullFleetRun:
    """End-to-end fleet run must produce all required artifacts."""

    def test_fleet_run_all_pass(self, tmp_path, manifest):
        evidence = run_fleet(tmp_path, nx=121, nz=25, volume_tolerance=0.05, write=True)
        assert evidence["summary"]["ships_admitted"] == 6
        assert evidence["summary"]["failed"] == []
        assert evidence["summary"]["controls_all_rejected"] is True

    def test_fleet_run_writes_evidence(self, tmp_path):
        run_fleet(tmp_path, nx=61, nz=13, volume_tolerance=0.10, write=True)
        assert (tmp_path / "evidence.json").exists()
        assert (tmp_path / "ref_vs_got.csv").exists()
        assert (tmp_path / "settings.json").exists()

    def test_fleet_run_writes_stl_artifacts(self, tmp_path, manifest):
        run_fleet(tmp_path, nx=61, nz=13, volume_tolerance=0.10, write=True)
        for ship in manifest.ships:
            stl_path = tmp_path / "manifests" / f"{ship.name}.stl"
            assert stl_path.exists(), f"Missing STL for {ship.name}"
            assert stl_path.stat().st_size > 1000, f"STL for {ship.name} too small"

    def test_fleet_run_artifacts_reingest_ok(self, tmp_path, manifest):
        evidence = run_fleet(tmp_path, nx=61, nz=13, volume_tolerance=0.10, write=True)
        assert evidence["summary"]["artifacts_reingest_ok"] is True
