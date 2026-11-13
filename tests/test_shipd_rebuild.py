"""Tests for tensorcad.shipd_rebuild module."""
import json
import pytest
from pathlib import Path

from tensorcad.shipd_rebuild import (
    rebuild_ship,
    run_shipd_rebuild,
    _select_method,
    _mesh_bbox,
    _euler_characteristic,
    SERIES60_CB_WINDOW,
    PARAMETRIC_CB_WINDOW,
    WIGLEY_CB_WINDOW,
)


class TestSelectMethod:
    """Method selection logic."""

    def test_series60_in_window(self):
        method, reason = _select_method(0.65, 0.98, 0.80)
        assert method == "series60_form"

    def test_series60_at_lower_bound(self):
        method, reason = _select_method(0.55, 0.98, 0.80)
        assert method == "series60_form"

    def test_series60_at_upper_bound(self):
        method, reason = _select_method(0.85, 0.98, 0.80)
        assert method == "series60_form"

    def test_parametric_below_series60(self):
        method, reason = _select_method(0.50, 0.95, 0.75)
        assert method == "parametric_sections"

    def test_parametric_above_series60(self):
        method, reason = _select_method(0.88, 0.99, 0.85)
        assert method == "parametric_sections"

    def test_wigley_fallback(self):
        method, reason = _select_method(0.35, 0.90, 0.70)
        assert method == "wigley_form"

    def test_refused_too_low(self):
        method, reason = _select_method(0.25, 0.90, 0.70)
        assert method == "refused"

    def test_refused_too_high(self):
        method, reason = _select_method(0.96, 0.99, 0.88)
        assert method == "refused"

    def test_cm_near_1_falls_back(self):
        """Cm ≈ 1.0 should fall back to parametric even if Cb is in Series 60 window."""
        method, reason = _select_method(0.65, 1.0, 0.80)
        assert method == "parametric_sections"


class TestMeshHelpers:
    """Mesh utility functions."""

    def test_euler_characteristic_sphere(self):
        """A closed mesh should have Euler characteristic 2."""
        from tensorcad.hull import HullParameters, generate_hull
        model = generate_hull(HullParameters(length=10, beam=2, draft=0.8, nx=21, nz=7))
        euler = _euler_characteristic(model.mesh)
        assert euler == 2

    def test_bbox_symmetric_hull(self):
        """Symmetric hull should have bbox centered at y=0."""
        from tensorcad.hull import HullParameters, generate_hull
        model = generate_hull(HullParameters(length=10, beam=2, draft=0.8, nx=21, nz=7))
        bbox = _mesh_bbox(model.mesh)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox
        assert abs(ymin + ymax) < 1e-10  # symmetric about y=0
        assert abs(xmin) < 1e-10  # starts at x=0
        assert abs(zmax) < 1e-10  # waterline at z=0


class TestRebuildShip:
    """Single ship rebuild."""

    def test_series60_ship_verified(self):
        """A ship with Cb in Series 60 window should rebuild successfully."""
        ship = {
            "ship_id": 1, "ship_type": "container",
            "L": 200.0, "B": 30.0, "T": 10.0,
            "Cb": 0.65, "Cm": 0.98, "Cwp": 0.80,
        }
        result = rebuild_ship(ship, nx=41, nz=13)
        assert result.method == "series60_form"
        assert result.mesh is not None
        assert result.triangle_count is not None
        assert result.triangle_count >= 1000

    def test_parametric_ship(self):
        """A ship with Cb outside Series 60 window uses parametric_sections."""
        ship = {
            "ship_id": 2, "ship_type": "naval",
            "L": 150.0, "B": 18.0, "T": 6.0,
            "Cb": 0.50, "Cm": 0.90, "Cwp": 0.70,
        }
        result = rebuild_ship(ship, nx=41, nz=13)
        assert result.method == "parametric_sections"

    def test_refused_ship(self):
        """A ship with Cb outside all windows should be refused."""
        ship = {
            "ship_id": 3, "ship_type": "special",
            "L": 50.0, "B": 10.0, "T": 3.0,
            "Cb": 0.25, "Cm": 0.85, "Cwp": 0.65,
        }
        result = rebuild_ship(ship, nx=41, nz=13)
        assert result.method == "refused"
        assert result.verdict == "REFUSED"

    def test_gates_present(self):
        """Rebuild result should have 8 gates."""
        ship = {
            "ship_id": 4, "ship_type": "tanker",
            "L": 250.0, "B": 45.0, "T": 15.0,
            "Cb": 0.75, "Cm": 0.99, "Cwp": 0.85,
        }
        result = rebuild_ship(ship, nx=41, nz=13)
        assert len(result.gates) == 8
        gate_names = [g.name for g in result.gates]
        assert "G1_watertight" in gate_names
        assert "G6_volume_consistency" in gate_names


class TestRunShipdRebuild:
    """Full pipeline test."""

    def test_small_fleet(self, tmp_path):
        """Run rebuild on a small fleet and check evidence structure."""
        ships = [
            {"ship_id": i, "ship_type": "container",
             "L": 200.0, "B": 30.0, "T": 10.0,
             "Cb": 0.65, "Cm": 0.98, "Cwp": 0.80}
            for i in range(5)
        ]
        evidence = run_shipd_rebuild(ships, nx=41, nz=13, out_dir=str(tmp_path))
        assert evidence["n_ships"] == 5
        assert "pass_rate" in evidence
        assert "method_distribution" in evidence
        assert "gate_pass_rates" in evidence
        assert "failure_mode_ledger" in evidence
        # Check artefacts written
        assert (tmp_path / "evidence.json").exists()
        assert (tmp_path / "ship_results.json").exists()
        assert (tmp_path / "summary.csv").exists()
        assert (tmp_path / "failure_ledger.json").exists()
