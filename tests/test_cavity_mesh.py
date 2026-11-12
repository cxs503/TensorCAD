"""Tests for tensorcad.cavity_mesh — lid-driven cavity mesh and reference data."""
from __future__ import annotations

import json
import math
import tempfile
from pathlib import Path

import pytest

from tensorcad.cavity_mesh import (
    GHIA_RE5000_U_CENTERLINE,
    GHIA_RE5000_V_CENTERLINE,
    CavityMesh,
    CavityParameters,
    generate_cavity_mesh,
    ghia_u_at_y,
    ghia_v_at_x,
    mesh_summary,
    mesh_to_json,
    write_mesh_json,
    write_reference_data,
)


# ---------------------------------------------------------------------------
# CavityParameters
# ---------------------------------------------------------------------------


class TestCavityParameters:
    def test_default_parameters(self):
        p = CavityParameters()
        assert p.reynolds_number == 5000.0
        assert p.lid_velocity == 1.0
        assert p.domain_size == 1.0

    def test_nu_derived(self):
        p = CavityParameters(reynolds_number=5000.0, lid_velocity=1.0, domain_size=1.0)
        assert abs(p.nu - 2e-4) < 1e-12

    def test_nu_formula(self):
        p = CavityParameters(reynolds_number=1000.0, lid_velocity=2.0, domain_size=0.5)
        expected = 2.0 * 0.5 / 1000.0
        assert abs(p.nu - expected) < 1e-15

    def test_custom_reynolds(self):
        p = CavityParameters(reynolds_number=100.0)
        assert p.reynolds_number == 100.0
        assert abs(p.nu - 0.01) < 1e-12


# ---------------------------------------------------------------------------
# Mesh generation
# ---------------------------------------------------------------------------


class TestGenerateCavityMesh:
    def test_default_mesh(self):
        mesh = generate_cavity_mesh(nx=10, ny=10)
        assert mesh.nx == 10
        assert mesh.ny == 10
        assert mesh.n_nodes == 11 * 11
        assert mesh.n_elements == 10 * 10

    def test_rectangular_mesh(self):
        mesh = generate_cavity_mesh(nx=20, ny=10)
        assert mesh.nx == 20
        assert mesh.ny == 10
        assert mesh.n_nodes == 21 * 11
        assert mesh.n_elements == 20 * 10

    def test_ny_defaults_to_nx(self):
        mesh = generate_cavity_mesh(nx=15)
        assert mesh.ny == 15

    def test_node_coordinates_bounds(self):
        mesh = generate_cavity_mesh(nx=8, ny=8)
        for x, y in mesh.nodes:
            assert 0.0 <= x <= 1.0 + 1e-12
            assert 0.0 <= y <= 1.0 + 1e-12

    def test_corner_nodes(self):
        mesh = generate_cavity_mesh(nx=4, ny=4)
        # Bottom-left
        assert mesh.nodes[0] == (0.0, 0.0)
        # Bottom-right
        assert mesh.nodes[4] == (1.0, 0.0)
        # Top-left
        assert mesh.nodes[-5] == (0.0, 1.0)
        # Top-right
        assert mesh.nodes[-1] == (1.0, 1.0)

    def test_element_connectivity_valid(self):
        mesh = generate_cavity_mesh(nx=5, ny=5)
        for elem in mesh.elements:
            assert len(elem) == 4
            for idx in elem:
                assert 0 <= idx < mesh.n_nodes

    def test_uniform_spacing(self):
        mesh = generate_cavity_mesh(nx=10, ny=10, wall_grading=1.0)
        # Check x-spacing is uniform
        for j in range(mesh.ny + 1):
            for i in range(mesh.nx):
                x0 = mesh.nodes[mesh.node_index(i, j)][0]
                x1 = mesh.nodes[mesh.node_index(i + 1, j)][0]
                assert abs((x1 - x0) - 0.1) < 1e-12

    def test_wall_grading_finer_near_walls(self):
        mesh = generate_cavity_mesh(nx=20, ny=20, wall_grading=3.0)
        # First cell should be smaller than middle cell
        dx_first = mesh.nodes[1][0] - mesh.nodes[0][0]
        dx_mid = mesh.nodes[11][0] - mesh.nodes[10][0]
        assert dx_first < dx_mid

    def test_custom_domain_size(self):
        params = CavityParameters(domain_size=2.0)
        mesh = generate_cavity_mesh(nx=5, ny=5, params=params)
        assert mesh.nodes[-1][0] == pytest.approx(2.0, abs=1e-12)
        assert mesh.nodes[-1][1] == pytest.approx(2.0, abs=1e-12)

    def test_minimum_nx(self):
        with pytest.raises(ValueError):
            generate_cavity_mesh(nx=1)

    def test_minimum_ny(self):
        with pytest.raises(ValueError):
            generate_cavity_mesh(nx=5, ny=1)

    def test_boundary_metadata(self):
        mesh = generate_cavity_mesh(nx=10, ny=10)
        assert "lid" in mesh.boundary
        assert "bottom_wall" in mesh.boundary
        assert "left_wall" in mesh.boundary
        assert "right_wall" in mesh.boundary
        assert mesh.boundary["lid"]["velocity"] == (1.0, 0.0)
        assert mesh.boundary["bottom_wall"]["velocity"] == (0.0, 0.0)

    def test_aspect_ratio_uniform(self):
        mesh = generate_cavity_mesh(nx=10, ny=10, wall_grading=1.0)
        ar_min, ar_max = mesh.aspect_ratio_range()
        assert ar_min == pytest.approx(1.0, abs=1e-10)
        assert ar_max == pytest.approx(1.0, abs=1e-10)


# ---------------------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------------------


class TestGhiaReferenceData:
    def test_u_profile_length(self):
        assert len(GHIA_RE5000_U_CENTERLINE) > 50

    def test_v_profile_length(self):
        assert len(GHIA_RE5000_V_CENTERLINE) > 50

    def test_u_profile_endpoints(self):
        # Bottom wall: u = 0
        assert GHIA_RE5000_U_CENTERLINE[0][0] == 0.0
        assert GHIA_RE5000_U_CENTERLINE[0][1] == 0.0
        # Lid: u = 1
        assert GHIA_RE5000_U_CENTERLINE[-1][0] == 1.0
        assert GHIA_RE5000_U_CENTERLINE[-1][1] == 1.0

    def test_v_profile_endpoints(self):
        # Left wall: v = 0
        assert GHIA_RE5000_V_CENTERLINE[0][0] == 0.0
        assert GHIA_RE5000_V_CENTERLINE[0][1] == 0.0
        # Right wall: v = 0
        assert GHIA_RE5000_V_CENTERLINE[-1][0] == 1.0
        assert GHIA_RE5000_V_CENTERLINE[-1][1] == 0.0

    def test_u_monotonic_y(self):
        for i in range(len(GHIA_RE5000_U_CENTERLINE) - 1):
            y0 = GHIA_RE5000_U_CENTERLINE[i][0]
            y1 = GHIA_RE5000_U_CENTERLINE[i + 1][0]
            assert y1 > y0, f"y not increasing at index {i}"

    def test_v_monotonic_x(self):
        for i in range(len(GHIA_RE5000_V_CENTERLINE) - 1):
            x0 = GHIA_RE5000_V_CENTERLINE[i][0]
            x1 = GHIA_RE5000_V_CENTERLINE[i + 1][0]
            assert x1 > x0, f"x not increasing at index {i}"

    def test_u_interpolation_at_wall(self):
        assert ghia_u_at_y(0.0) == 0.0

    def test_u_interpolation_at_lid(self):
        assert ghia_u_at_y(1.0) == 1.0

    def test_v_interpolation_at_walls(self):
        assert ghia_v_at_x(0.0) == 0.0
        assert ghia_v_at_x(1.0) == 0.0

    def test_u_interpolation_midrange(self):
        # At y=0.5, u should be small and positive (secondary vortex region)
        u_mid = ghia_u_at_y(0.5)
        assert 0.0 < u_mid < 0.05

    def test_v_interpolation_midrange(self):
        # At x=0.5, v should be small
        v_mid = ghia_v_at_x(0.5)
        assert abs(v_mid) < 0.01


# ---------------------------------------------------------------------------
# Export utilities
# ---------------------------------------------------------------------------


class TestExportUtilities:
    def test_mesh_to_json_keys(self):
        mesh = generate_cavity_mesh(nx=4, ny=4)
        d = mesh_to_json(mesh)
        assert "nx" in d
        assert "ny" in d
        assert "n_nodes" in d
        assert "n_elements" in d
        assert "reynolds_number" in d
        assert "nodes" in d
        assert "elements" in d
        assert "boundary" in d

    def test_write_mesh_json(self, tmp_path):
        mesh = generate_cavity_mesh(nx=4, ny=4)
        path = write_mesh_json(mesh, tmp_path / "mesh.json")
        assert path.exists()
        with open(path) as f:
            data = json.load(f)
        assert data["nx"] == 4
        assert data["ny"] == 4
        assert len(data["nodes"]) == 25

    def test_write_reference_data(self, tmp_path):
        path = write_reference_data(tmp_path / "ref.json")
        assert path.exists()
        with open(path) as f:
            data = json.load(f)
        assert data["reynolds_number"] == 5000.0
        assert "u_centerline_x0.5" in data
        assert "v_centerline_y0.5" in data
        assert "Ghia" in data["source"]

    def test_mesh_summary_content(self):
        mesh = generate_cavity_mesh(nx=8, ny=8)
        summary = mesh_summary(mesh)
        assert "Re = 5000" in summary
        assert "8 × 8" in summary
        assert "Lid" in summary
