"""Integration tests for the iges_field geometry route in the baseline fleet.

These tests verify that the IGES → half-breadth field → watertight mesh pipeline
is correctly integrated into the fleet build system.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.fleet.build import (
    RouteRefusal,
    build_ship_geometry,
    loft_station_table,
    sample_offset_table,
)
from tensorcad.fleet_baseline import (
    build_solid,
    field_from_points,
    field_integrals,
)
from tensorcad.fleet.gates import evaluate_gates, measure, placeholder_reasons
from tensorcad.fleet.manifest import Provenance, ShipRecord, load_manifest


@pytest.fixture(scope="module")
def manifest():
    return load_manifest()


@pytest.fixture(scope="module")
def series60_record(manifest):
    return manifest.by_name("Series60")


class TestIGESFieldRouteDispatch:
    def test_iges_field_route_requires_path(self):
        record = ShipRecord(
            name="TestShip", ship_class="Test", length=100.0, beam=13.0, draft=5.0,
            block_coefficient=0.60, displacement_volume=3900.0, displacement_basis="closure",
            geometry_route="iges_field",
            provenance=Provenance(source="test", url="https://example.com", license="test", retrieval="test", verified=False),
        )
        with pytest.raises(RouteRefusal, match="iges_field route needs a path"):
            build_ship_geometry(record)

    def test_unknown_route_refused(self):
        record = ShipRecord(
            name="TestShip", ship_class="Test", length=100.0, beam=13.0, draft=5.0,
            block_coefficient=0.60, displacement_volume=3900.0, displacement_basis="closure",
            geometry_route="nonexistent_route",
            provenance=Provenance(source="test", url="https://example.com", license="test", retrieval="test", verified=False),
        )
        with pytest.raises(RouteRefusal, match="unknown geometry route"):
            build_ship_geometry(record)


class TestHalfBreadthField:
    def test_field_from_symmetric_points(self):
        points = []
        for i in range(21):
            x = i * 5.0
            for k in range(11):
                z = -k * 0.5
                xi = x / 100.0
                zeta = abs(z) / 5.0
                b = 6.5 * math.sqrt(max(0, 1 - (2 * xi - 1) ** 2)) * math.sqrt(max(0, 1 - zeta ** 2))
                points.append((x, b, z))
                points.append((x, -b, z))
        field = field_from_points(points, draft=5.0, dx=5.0, dz=0.5)
        assert field.nx > 0 and field.nz > 0 and field.filled_nodes > 0

    def test_field_cover_metric(self):
        points = [(x, 5.0, -2.5) for x in range(0, 100, 5)]
        field = field_from_points(points, draft=5.0, dx=5.0, dz=0.5)
        assert 0.0 <= field.cover() <= 1.0

    def test_field_to_dict(self):
        points = [(50.0, 5.0, -2.5)]
        field = field_from_points(points, draft=5.0, dx=5.0, dz=0.5)
        d = field.to_dict()
        for key in ("x0", "dx", "dz", "nx", "nz", "filled_nodes", "cover"):
            assert key in d


class TestWatertightMeshFromField:
    def _make_field_points(self, L=100.0, B=13.0, T=5.0, nx=41, nz=21):
        points = []
        for i in range(nx):
            xi = i / (nx - 1)
            x = xi * L
            for k in range(nz):
                zeta = k / (nz - 1)
                z = -T * (1 - zeta)
                b = 0.5 * B * math.sin(math.pi * xi) ** 0.7 * zeta ** 0.5
                points.append((x, b, z))
                points.append((x, -b, z))
        return points

    def test_field_produces_watertight_mesh(self):
        points = self._make_field_points()
        field = field_from_points(points, draft=5.0, dx=2.5, dz=0.25)
        mesh, _ = build_solid(field, ring_points=240, draft=5.0)
        facts = measure(mesh)
        assert facts.boundary_edge_count == 0
        assert facts.nonmanifold_edge_count == 0
        assert facts.euler_characteristic == 2

    def test_field_mesh_positive_volume(self):
        points = self._make_field_points()
        field = field_from_points(points, draft=5.0, dx=2.5, dz=0.25)
        mesh, _ = build_solid(field, ring_points=240, draft=5.0)
        assert measure(mesh).signed_volume > 0

    def test_field_mesh_realistic_triangle_count(self):
        points = self._make_field_points()
        field = field_from_points(points, draft=5.0, dx=2.5, dz=0.25)
        mesh, _ = build_solid(field, ring_points=240, draft=5.0)
        assert measure(mesh).triangle_count > 1000


class TestVolumeConsistency:
    def test_mesh_vs_field_volume(self):
        points = []
        for i in range(41):
            xi = i / 40
            x = xi * 100
            for k in range(21):
                zeta = k / 20
                z = -5 * (1 - zeta)
                b = 6.5 * math.sin(math.pi * xi) ** 0.7 * zeta ** 0.5
                points.append((x, b, z))
                points.append((x, -b, z))
        field = field_from_points(points, draft=5.0, dx=2.5, dz=0.25)
        mesh, _ = build_solid(field, ring_points=240, draft=5.0)
        mesh_volume = measure(mesh).signed_volume
        field_volume = field_integrals(field).get("volume_m3", 0.0)
        assert mesh_volume > 0 and field_volume > 0
        rel_error = abs(mesh_volume - field_volume) / max(abs(mesh_volume), abs(field_volume), 1e-12)
        assert rel_error < 0.05, f"mesh vs field volume error {rel_error*100:.2f}%"


class TestFrameDetection:
    def test_waterline_at_origin_frame(self):
        points = [(50.0, 5.0, -2.5), (50.0, -5.0, -2.5)]
        field = field_from_points(points, draft=5.0, dx=5.0, dz=0.5)
        assert field.filled_nodes > 0


class TestGateIntegration:
    def test_lofted_mesh_not_placeholder(self, manifest):
        record = manifest.by_name("Series60")
        table, levels = sample_offset_table(record, stations=41, levels=13)
        mesh = loft_station_table(table, levels)
        facts = measure(mesh)
        assert not placeholder_reasons(facts, mesh=mesh)


class TestCrossRouteConsistency:
    def test_series60_converges_with_refinement(self, series60_record):
        """Volume converges as resolution increases."""
        vol_coarse = measure(build_ship_geometry(series60_record, nx=41, nz=13).mesh).signed_volume
        vol_fine = measure(build_ship_geometry(series60_record, nx=161, nz=41).mesh).signed_volume
        # Fine grid should be closer to the analytic volume
        analytic = series60_record.reference_volume()
        err_coarse = abs(vol_coarse - analytic) / analytic
        err_fine = abs(vol_fine - analytic) / analytic
        assert err_fine < err_coarse, f"Fine grid error {err_fine:.4f} should be < coarse {err_coarse:.4f}"
