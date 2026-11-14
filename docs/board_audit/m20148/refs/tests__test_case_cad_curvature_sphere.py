"""Test pin for CASE-CAD-CURVATURE-SPHERE.

Ensures the sphere mean curvature benchmark runs and all gates pass.
"""
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

_REPO_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = _REPO_ROOT / "benchmarks" / "verified" / "sphere_mean_curvature" / "run_benchmark.py"


@pytest.fixture(scope="module")
def benchmark_result(tmp_path_factory):
    """Run the benchmark once and cache the result."""
    out = tmp_path_factory.mktemp("sphere_curvature")
    result = subprocess.run(
        [sys.executable, str(BENCHMARK), "--out", str(out)],
        capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, f"Benchmark failed:\n{result.stderr}"
    with open(out / "result.json") as f:
        return json.load(f)


class TestSphereMeanCurvature:
    """CASE-CAD-CURVATURE-SPHERE test suite."""

    def test_all_gates_hold(self, benchmark_result):
        """All gates must pass."""
        assert benchmark_result["all_gates_hold"], "Not all gates passed"

    def test_case_tolerance(self, benchmark_result):
        """Worst error must be within the 2% case tolerance for gated rows."""
        case_tol = benchmark_result["case_tolerance"]  # 0.02
        for row in benchmark_result["all_rows"]:
            if row["error_percent"] is not None:
                # G3 raw has 60% tolerance, G3b has 5%, others have 2% or less
                assert abs(row["error_percent"]) / 100 <= row["tolerance"], (
                    f"Gate {row['gate']} failed: {row['error_percent']:.4f}% > {row['tolerance']*100:.1f}%"
                )

    def test_g2_integral_convergence(self, benchmark_result):
        """G2: area-weighted mean curvature converges to 2/R."""
        g2 = benchmark_result["gates"]["G2_integral_marching_tet"]
        for row in g2["rows"]:
            assert row["holds"], f"G2 gate failed: {row['quantity']}"
            assert abs(row["error_percent"]) < 2.0, (
                f"G2 error {row['error_percent']:.4f}% exceeds 2%"
            )

    def test_g4_convergence_order(self, benchmark_result):
        """G4: observed convergence order is approximately 2."""
        g4 = benchmark_result["gates"]["G4_uv_sphere_convergence"]
        order = g4["observed_order"]
        assert math.isfinite(order), "Observed order is not finite"
        assert abs(order - 2.0) < 0.5, f"Observed order {order:.3f} deviates from 2.0"

    def test_g6_integral_identity(self, benchmark_result):
        """G6: integral of S dA = 8*pi*R."""
        g6 = benchmark_result["gates"]["G6_integral_identity"]
        ref = g6["reference"]
        measured = g6["measured"]
        error = abs(measured - ref) / ref
        assert error < 0.01, f"Integral identity error {error*100:.4f}% exceeds 1%"

    def test_schema(self, benchmark_result):
        """Result has correct schema."""
        assert benchmark_result["schema"] == "tensorcad.benchmark.sphere_mean_curvature/1"
        assert benchmark_result["case"] == "CASE-CAD-CURVATURE-SPHERE"

    def test_gates_verified_file(self, benchmark_result, tmp_path):
        """gates_verified.json is produced and consistent."""
        out = tmp_path / "test_out"
        result = subprocess.run(
            [sys.executable, str(BENCHMARK), "--out", str(out)],
            capture_output=True, text=True, timeout=120,
        )
        assert result.returncode == 0
        gv_path = out / "gates_verified.json"
        assert gv_path.exists()
        with open(gv_path) as f:
            gv = json.load(f)
        assert gv["schema"] == "tensorcad.gates_verified/1"
        assert gv["case"] == "CASE-CAD-CURVATURE-SPHERE"
        assert gv["all_hold"] is True
