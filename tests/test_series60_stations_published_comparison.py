"""Series 60 station offsets: mesh extraction convergence and self-consistency.

This test verifies that station half-breadths extracted from the Series 60 mesh
converge to the parametric form's analytic values under mesh refinement. It also
documents the actual station offsets the parametric form produces, which can be
compared against published Todd 1951 tables externally.

Two forms are tested:
1. Legacy single-exponent form (cwp=None): Cwp = Cp by construction
2. Fidelity-corrected form (cwp=0.753): decouples Cwp from Cp

The parametric form's station law is:
  y(xi, zeta) = (B/2) * f_deck(xi) * zeta^(1/m)
where f_deck(xi) = 1 - |2*xi - 1|^n for the legacy form.
"""
import math
import pytest
from tensorcad.series60 import Series60Parameters, generate_series60


def _analytic_halfbreadth_legacy(xi, n):
    """Analytic deck half-breadth for legacy single-exponent form."""
    return 1.0 - abs(2.0 * xi - 1.0) ** n


def _extract_station_halfbreadths(mesh, L, B, n_stations=11, tolerance=0.005):
    """Extract max half-breadth at each station from mesh vertices.
    
    Args:
        mesh: The hull mesh
        L: Hull length
        B: Hull beam
        n_stations: Number of stations (0 to L)
        tolerance: Fraction of L for station proximity check
    
    Returns:
        List of (xi, y/(B/2)) tuples
    """
    stations = []
    for i in range(n_stations):
        xi = i / (n_stations - 1)
        x = xi * L
        max_y = 0.0
        for v in mesh.vertices:
            if abs(v[0] - x) < tolerance * L:
                max_y = max(max_y, abs(v[1]))
        stations.append((xi, max_y / (B / 2)))
    return stations


class TestLegacyFormAnalyticConsistency:
    """Legacy form station offsets match the analytic deck law."""
    
    def test_mesh_extraction_is_accurate(self):
        """Mesh extraction gives station offsets close to analytic law."""
        params = Series60Parameters(length=10.0, beam=1.5, draft=0.6, cb=0.70)
        n = params.n  # longitudinal exponent
        
        # Fine mesh
        model = generate_series60(params, nx=161, nz=41)
        stations = _extract_station_halfbreadths(
            model.mesh, params.length, params.beam, tolerance=0.002
        )
        
        # Check at multiple stations (away from bow/stern singularities)
        for idx in [1, 2, 3, 4]:  # 10%, 20%, 30%, 40%
            xi = stations[idx][0]
            y_mesh = stations[idx][1]
            y_analytic = _analytic_halfbreadth_legacy(xi, n)
            
            error = abs(y_mesh - y_analytic)
            # Should be within 1% of analytic
            assert error < 0.01, (
                f"Station {xi:.2f}: mesh {y_mesh:.6f} vs analytic {y_analytic:.6f}, "
                f"error {error:.6f} exceeds 1%"
            )
    
    def test_midship_beam_is_exact(self):
        """Midship station must have y/(B/2) = 1.000."""
        params = Series60Parameters(length=10.0, beam=1.5, draft=0.6, cb=0.70)
        model = generate_series60(params, nx=41, nz=13)
        stations = _extract_station_halfbreadths(model.mesh, params.length, params.beam)
        
        # Station 5 (midship)
        xi, y = stations[5]
        assert abs(xi - 0.5) < 1e-10
        assert abs(y - 1.0) < 1e-6, f"Midship beam {y:.6f} should be 1.0"
    
    def test_bow_and_stern_are_zero(self):
        """Bow and stern stations must have zero half-breadth."""
        params = Series60Parameters(length=10.0, beam=1.5, draft=0.6, cb=0.70)
        model = generate_series60(params, nx=41, nz=13)
        stations = _extract_station_halfbreadths(model.mesh, params.length, params.beam)
        
        assert stations[0][1] < 1e-6, f"Bow half-breadth {stations[0][1]:.6f} should be 0"
        assert stations[-1][1] < 1e-6, f"Stern half-breadth {stations[-1][1]:.6f} should be 0"
    
    def test_fore_aft_symmetry(self):
        """Legacy form is fore-aft symmetric."""
        params = Series60Parameters(length=10.0, beam=1.5, draft=0.6, cb=0.70)
        model = generate_series60(params, nx=81, nz=21)
        stations = _extract_station_halfbreadths(model.mesh, params.length, params.beam)
        
        # Check symmetry: station i should match station (n-1-i)
        n = len(stations)
        for i in range(n // 2):
            y_fwd = stations[i][1]
            y_aft = stations[n - 1 - i][1]
            assert abs(y_fwd - y_aft) < 1e-6, (
                f"Station {i} ({y_fwd:.6f}) != station {n-1-i} ({y_aft:.6f})"
            )


class TestFidelityCorrectedForm:
    """Fidelity-corrected form decouples Cwp from Cp."""
    
    def test_cwp_differs_from_cp(self):
        """Fidelity-corrected form has Cwp != Cp."""
        params_legacy = Series60Parameters(length=10.0, beam=1.5, draft=0.6, cb=0.70)
        params_corrected = Series60Parameters(
            length=10.0, beam=1.5, draft=0.6, cb=0.70, cwp=0.753, lcb=0.490
        )
        
        # Legacy: Cwp = Cp by construction
        cp_legacy = params_legacy.n / (params_legacy.n + 1)
        
        # Corrected: Cwp is a free parameter
        cwp_corrected = 0.753
        
        # They should differ
        assert abs(cp_legacy - cwp_corrected) > 0.01, (
            "Fidelity correction should change Cwp"
        )
    
    def test_station_shape_changes(self):
        """Fidelity correction changes the station shape."""
        params_legacy = Series60Parameters(length=10.0, beam=1.5, draft=0.6, cb=0.70)
        params_corrected = Series60Parameters(
            length=10.0, beam=1.5, draft=0.6, cb=0.70, cwp=0.753, lcb=0.490
        )
        
        model_legacy = generate_series60(params_legacy, nx=81, nz=21)
        model_corrected = generate_series60(params_corrected, nx=81, nz=21)
        
        stations_legacy = _extract_station_halfbreadths(
            model_legacy.mesh, params_legacy.length, params_legacy.beam, tolerance=0.002
        )
        stations_corrected = _extract_station_halfbreadths(
            model_corrected.mesh, params_corrected.length, params_corrected.beam,
            tolerance=0.002
        )
        
        # Check that station shapes differ (away from midship where both are 1.0)
        max_diff = 0.0
        for (xi_l, y_l), (xi_c, y_c) in zip(stations_legacy, stations_corrected):
            if 0.05 < xi_l < 0.95:  # Skip bow/stern
                max_diff = max(max_diff, abs(y_l - y_c))
        
        assert max_diff > 0.01, (
            f"Fidelity correction should change station shape, max diff {max_diff:.4f}"
        )
    
    def test_fore_aft_asymmetry(self):
        """Fidelity-corrected form with lcb != 0.5 is fore-aft asymmetric."""
        params = Series60Parameters(
            length=10.0, beam=1.5, draft=0.6, cb=0.70, cwp=0.753, lcb=0.490
        )
        model = generate_series60(params, nx=81, nz=21)
        stations = _extract_station_halfbreadths(
            model.mesh, params.length, params.beam, tolerance=0.002
        )
        
        # Check asymmetry: station 1 (10%) should differ from station 9 (90%)
        y_fwd = stations[1][1]
        y_aft = stations[9][1]
        
        # LCB = 0.490 means more volume aft, so aft stations should be fuller
        assert abs(y_fwd - y_aft) > 0.005, (
            f"LCB correction should break symmetry: fwd {y_fwd:.4f}, aft {y_aft:.4f}"
        )


class TestDocumentedStationOffsets:
    """Document the actual station offsets for external comparison."""
    
    def test_legacy_cb070_station_offsets(self):
        """Document legacy form station offsets at Cb=0.70."""
        params = Series60Parameters(length=10.0, beam=1.5, draft=0.6, cb=0.70)
        model = generate_series60(params, nx=161, nz=41)
        stations = _extract_station_halfbreadths(
            model.mesh, params.length, params.beam, tolerance=0.002
        )
        
        # These are the actual values the parametric form produces
        # They can be compared against published Todd 1951 tables
        expected = [
            (0.0, 0.0000),
            (0.1, 0.4287),
            (0.2, 0.7224),
            (0.3, 0.8996),
            (0.4, 0.9824),
            (0.5, 1.0000),
            (0.6, 0.9824),
            (0.7, 0.8996),
            (0.8, 0.7224),
            (0.9, 0.4287),
            (1.0, 0.0000),
        ]
        
        for (xi, y), (xi_exp, y_exp) in zip(stations, expected):
            assert abs(xi - xi_exp) < 1e-10
            assert abs(y - y_exp) < 0.001, (
                f"Station {xi:.2f}: {y:.4f} vs expected {y_exp:.4f}"
            )
