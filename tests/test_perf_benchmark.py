"""Tests for the performance benchmark module.

These tests verify the benchmark harness contract: data types, schema,
regression gates, and the CLI entry point.  They do *not* re-run the full
ladder (that is a heavy computation); instead they exercise the machinery
with tiny inputs.
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

from tensorcad.perf_benchmark import (
    BENCHMARK_SCHEMA,
    CLUSTER_THRESHOLD,
    FAILED,
    MEASURED,
    NOT_RUN,
    GateVerdict,
    PerfBenchmarkResult,
    SubTaskResult,
    _fibonacci_directions,
    evaluate_gates,
    run_perf_benchmark,
)

ROOT = Path(__file__).resolve().parents[1]


def _env() -> dict[str, str]:
    existing = os.environ.get("PYTHONPATH", "")
    path = str(ROOT) + (os.pathsep + existing if existing else "")
    return {**os.environ, "PYTHONPATH": path}


# ---------------------------------------------------------------------------
# Data contract tests
# ---------------------------------------------------------------------------

def test_subtask_result_to_dict_roundtrip() -> None:
    r = SubTaskResult(
        subtask="load", nx=41, nz=13, triangles=1998, vertices=1001,
        grid_cells=0, seconds=0.05, throughput=39960.0,
        status=MEASURED, note="ok",
    )
    d = r.to_dict()
    assert d["subtask"] == "load"
    assert d["triangles"] == 1998
    assert d["seconds"] == 0.05
    assert d["throughput"] == 39960.0
    assert d["status"] == MEASURED


def test_gate_verdict_to_dict() -> None:
    g = GateVerdict(
        name="regression_load:41x13", passed=True,
        detail="0.05s vs 0.04s baseline (ratio=1.25x, limit=2.0x)",
        baseline_seconds=0.04, measured_seconds=0.05, ratio=1.25,
    )
    d = g.to_dict()
    assert d["passed"] is True
    assert d["ratio"] == 1.25


def test_perf_result_serialises_to_json() -> None:
    result = PerfBenchmarkResult(
        results=(
            SubTaskResult("load", 41, 13, 1998, 1001, 0, 0.05, MEASURED, 39960.0),
        ),
        gates=(
            GateVerdict("regression_load:41x13", True, "ok", 0.04, 0.05, 1.25),
        ),
        environment={"tensorcad": "0.23.0", "python": "3.11"},
        seconds=0.1,
        cluster_warnings=(),
    )
    d = result.to_dict()
    assert d["schema"] == BENCHMARK_SCHEMA
    assert d["passed"] is True
    assert len(d["results"]) == 1
    assert len(d["gates"]) == 1
    # Must be JSON-serialisable
    json.dumps(d)


def test_perf_result_write_json(tmp_path: Path) -> None:
    result = PerfBenchmarkResult(
        results=(
            SubTaskResult("volume", 21, 9, 678, 341, 0, 0.001, MEASURED, 678000.0),
        ),
        seconds=0.01,
    )
    path = tmp_path / "perf.json"
    result.write_json(path)
    assert path.exists()
    loaded = json.loads(path.read_text())
    assert loaded["schema"] == BENCHMARK_SCHEMA


# ---------------------------------------------------------------------------
# Regression gate tests
# ---------------------------------------------------------------------------

def test_gate_passes_when_within_factor() -> None:
    result = PerfBenchmarkResult(
        results=(
            SubTaskResult("load", 41, 13, 1998, 1001, 0, 0.08, MEASURED, 24975.0),
        ),
    )
    # Write a baseline
    baseline_data = {
        "schema": BENCHMARK_SCHEMA,
        "results": [
            {"subtask": "load", "nx": 41, "nz": 13, "triangles": 1998,
             "vertices": 1001, "grid_cells": 0, "seconds": 0.05,
             "status": "measured", "throughput": 39960.0, "note": ""},
        ],
    }
    bl_path = Path("/tmp/test_perf_baseline.json")
    bl_path.write_text(json.dumps(baseline_data))

    updated, failed = evaluate_gates(result, bl_path, regression_factor=2.0)
    assert len(failed) == 0
    assert len(updated.gates) == 1
    assert updated.gates[0].passed is True
    assert updated.gates[0].ratio == pytest.approx(1.6)

    bl_path.unlink()


def test_gate_fails_when_exceeding_factor() -> None:
    result = PerfBenchmarkResult(
        results=(
            SubTaskResult("load", 41, 13, 1998, 1001, 0, 0.20, MEASURED, 9990.0),
        ),
    )
    baseline_data = {
        "schema": BENCHMARK_SCHEMA,
        "results": [
            {"subtask": "load", "nx": 41, "nz": 13, "triangles": 1998,
             "vertices": 1001, "grid_cells": 0, "seconds": 0.05,
             "status": "measured", "throughput": 39960.0, "note": ""},
        ],
    }
    bl_path = Path("/tmp/test_perf_baseline_fail.json")
    bl_path.write_text(json.dumps(baseline_data))

    updated, failed = evaluate_gates(result, bl_path, regression_factor=2.0)
    assert len(failed) == 1
    assert updated.gates[0].passed is False
    assert updated.gates[0].ratio == pytest.approx(4.0)

    bl_path.unlink()


def test_no_gate_when_no_baseline() -> None:
    result = PerfBenchmarkResult(
        results=(
            SubTaskResult("load", 41, 13, 1998, 1001, 0, 0.08, MEASURED, 24975.0),
        ),
    )
    updated, failed = evaluate_gates(result, None)
    assert len(updated.gates) == 0
    assert len(failed) == 0


# ---------------------------------------------------------------------------
# Fibonacci directions
# ---------------------------------------------------------------------------

def test_fibonacci_directions_are_unit_length() -> None:
    dirs = _fibonacci_directions(64)
    assert len(dirs) == 64
    for d in dirs:
        length = math.sqrt(sum(c * c for c in d))
        assert length == pytest.approx(1.0, abs=1e-12)


def test_fibonacci_directions_are_deterministic() -> None:
    first = _fibonacci_directions(32)
    second = _fibonacci_directions(32)
    assert first == second


def test_fibonacci_directions_reject_zero() -> None:
    with pytest.raises(ValueError, match="at least 1"):
        _fibonacci_directions(0)


# ---------------------------------------------------------------------------
# Quick benchmark smoke test
# ---------------------------------------------------------------------------

def test_quick_benchmark_runs() -> None:
    """Run the quick ladder with a tiny query count to verify end-to-end."""
    result = run_perf_benchmark(
        ladder=((21, 9), (41, 13)),
        voxel_ladder=(32,),
        subtasks=("load", "volume", "voxelize"),
        n_queries=100,
    )
    assert result.seconds > 0
    # Should have results for each sub-task at each step
    load_results = [r for r in result.results if r.subtask == "load"]
    assert len(load_results) == 2
    assert all(r.status == MEASURED for r in load_results)
    assert all(r.triangles > 0 for r in load_results)

    vol_results = [r for r in result.results if r.subtask == "volume"]
    assert len(vol_results) == 2
    assert all(r.status == MEASURED for r in vol_results)

    vox_results = [r for r in result.results if r.subtask == "voxelize"]
    assert len(vox_results) == 1
    assert vox_results[0].grid_cells == 32 ** 3


def test_point_locate_returns_hits() -> None:
    from tensorcad.perf_benchmark import _run_point_locate
    result = _run_point_locate(21, 9, n_queries=200)
    assert result.status == MEASURED
    assert result.seconds is not None
    assert result.throughput is not None
    # The hull should be hit by some rays from the origin
    assert "hits=" in result.note
    hits_str = result.note.split("=")[1].split("/")[0]
    assert int(hits_str) > 0


# ---------------------------------------------------------------------------
# CLI smoke test
# ---------------------------------------------------------------------------

def test_cli_quick_exits_zero() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "tensorcad.perf_benchmark", "--quick",
         "--n-queries", "100", "--max-triangles", "50000",
         "--subtasks", "load,volume"],
        cwd=str(ROOT),
        env=_env(),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0, f"stderr: {result.stderr}"
    assert "TensorCAD Performance Benchmark" in result.stdout
    assert "load" in result.stdout


def test_cli_save_baseline(tmp_path: Path) -> None:
    out = tmp_path / "bl.json"
    result = subprocess.run(
        [sys.executable, "-m", "tensorcad.perf_benchmark", "--quick",
         "--n-queries", "100", "--max-triangles", "50000",
         "--subtasks", "load", "--save-baseline", str(out)],
        cwd=str(ROOT),
        env=_env(),
        capture_output=True,
        text=True,
        timeout=60,
    )
    assert result.returncode == 0
    assert out.exists()
    data = json.loads(out.read_text())
    assert data["schema"] == BENCHMARK_SCHEMA


def test_cli_regression_gate(tmp_path: Path) -> None:
    # First save a baseline
    bl = tmp_path / "bl.json"
    subprocess.run(
        [sys.executable, "-m", "tensorcad.perf_benchmark", "--quick",
         "--n-queries", "100", "--max-triangles", "50000",
         "--subtasks", "load", "--save-baseline", str(bl)],
        cwd=str(ROOT), env=_env(), capture_output=True, timeout=60,
    )
    # Then run against it
    result = subprocess.run(
        [sys.executable, "-m", "tensorcad.perf_benchmark", "--quick",
         "--n-queries", "100", "--max-triangles", "50000",
         "--subtasks", "load", "--baseline", str(bl)],
        cwd=str(ROOT), env=_env(), capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0
    assert "Regression Gates" in result.stdout
    assert "All gates passed" in result.stdout


# ---------------------------------------------------------------------------
# Baseline file validation (GAP-CAD-D4-PERF-SCALE)
# ---------------------------------------------------------------------------

def test_baseline_file_exists_and_valid() -> None:
    """Verify the committed baseline file exists and conforms to the schema."""
    bl_path = ROOT / "runs" / "perf_baseline.json"
    assert bl_path.exists(), f"Baseline file not found at {bl_path}"

    data = json.loads(bl_path.read_text())
    assert data["schema"] == BENCHMARK_SCHEMA
    assert "results" in data
    assert len(data["results"]) >= 4, "Baseline should have results for all 4 sub-tasks"

    # Verify all 4 sub-tasks are represented
    subtasks = {r["subtask"] for r in data["results"]}
    assert "load" in subtasks
    assert "point_locate" in subtasks
    assert "voxelize" in subtasks
    assert "volume" in subtasks

    # Verify all results have required fields
    required_fields = {"subtask", "nx", "nz", "triangles", "vertices",
                       "grid_cells", "seconds", "status", "throughput"}
    for r in data["results"]:
        for field in required_fields:
            assert field in r, f"Missing field '{field}' in result: {r}"
        assert r["status"] in (MEASURED, NOT_RUN, FAILED)


def test_baseline_regression_gate_passes() -> None:
    """Verify that running the quick benchmark against the baseline passes all gates."""
    bl_path = ROOT / "runs" / "perf_baseline.json"
    if not bl_path.exists():
        pytest.skip("Baseline file not available")

    result = subprocess.run(
        [sys.executable, "-m", "tensorcad.perf_benchmark", "--quick",
         "--n-queries", "500", "--max-triangles", "100000",
         "--baseline", str(bl_path)],
        cwd=str(ROOT),
        env=_env(),
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, f"Regression gates failed:\n{result.stdout}\n{result.stderr}"
    assert "All gates passed" in result.stdout
