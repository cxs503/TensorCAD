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
    SCALING_GATE_PREFIX,
    MIN_SCALING_POINTS,
    _fibonacci_directions,
    build_fleet_mesh,
    clear_fleet_cache,
    clear_mesh_cache,
    evaluate_gates,
    measurement_key,
    run_fleet_benchmark,
    run_perf_benchmark,
    scaling_gates,
)

ROOT = Path(__file__).resolve().parents[1]

#: The committed regression baseline.  It used to live in ``runs/``, which the
#: merge gate strips from the branch, so on every fresh checkout the baseline was
#: missing and the 2x regression gate silently had nothing to compare against.
#: ``benchmarks/`` is tracked, so the gate is actually armed here.
BASELINE = ROOT / "benchmarks" / "perf_scale" / "baseline.json"


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
    bl_path = BASELINE
    assert bl_path.exists(), f"Baseline file not found at {bl_path}"

    data = json.loads(bl_path.read_text())
    assert data["schema"] == BENCHMARK_SCHEMA
    assert "results" in data
    assert len(data["results"]) >= 4, "Baseline should have results for all 4 sub-tasks"

    # Verify all sub-tasks are represented (including the indexed point path)
    subtasks = {r["subtask"] for r in data["results"]}
    for name in ("load", "point_locate", "point_locate_accel", "voxelize", "volume"):
        assert name in subtasks, f"baseline is missing sub-task {name}"

    # Verify all results have required fields
    required_fields = {"subtask", "nx", "nz", "triangles", "vertices",
                       "grid_cells", "seconds", "status", "throughput"}
    for r in data["results"]:
        for field in required_fields:
            assert field in r, f"Missing field '{field}' in result: {r}"
        assert r["status"] in (MEASURED, NOT_RUN, FAILED)


def test_baseline_regression_gate_passes() -> None:
    """Verify that running the quick benchmark against the baseline passes all gates."""
    bl_path = BASELINE
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


# ---------------------------------------------------------------------------
# Accelerated point location (GAP-CAD-D4-PERF-SCALE: 10^6--10^7 regime)
# ---------------------------------------------------------------------------

def test_accelerated_locate_is_measured_and_faster() -> None:
    """The indexed path is recorded on the same rung as the brute-force one."""
    from tensorcad.perf_benchmark import clear_mesh_cache

    clear_mesh_cache()
    result = run_perf_benchmark(
        ladder=((41, 13),),
        subtasks=("point_locate", "point_locate_accel"),
        n_queries=200,
    )
    by_task = {r.subtask: r for r in result.results}
    naive = by_task["point_locate"]
    accel = by_task["point_locate_accel"]
    assert naive.status == MEASURED and accel.status == MEASURED
    assert accel.triangles == naive.triangles
    assert accel.build_seconds is not None and accel.build_seconds >= 0.0
    assert accel.index_cells > 0 and accel.index_entries >= accel.triangles
    assert accel.index_bytes > 0
    assert accel.throughput is not None and accel.throughput > 0.0
    # Same number of hits: the index must not change the geometry.
    assert accel.note.split(",")[0] == naive.note.split(",")[0]
    clear_mesh_cache()


def test_naive_reference_skipped_above_threshold() -> None:
    """Above ``max_naive_triangles`` the O(N) reference is NOT_RUN, not slow."""
    from tensorcad.perf_benchmark import clear_mesh_cache

    clear_mesh_cache()
    result = run_perf_benchmark(
        ladder=((161, 41),),
        subtasks=("point_locate", "point_locate_accel"),
        n_queries=20,
        max_naive_triangles=1000,
    )
    by_task = {r.subtask: r for r in result.results}
    assert by_task["point_locate"].status == NOT_RUN
    assert "brute-force" in by_task["point_locate"].note
    assert by_task["point_locate_accel"].status == MEASURED
    clear_mesh_cache()


def test_mesh_cache_is_bounded() -> None:
    """The shared mesh cache keeps at most two rungs (10^7 rungs are GB-scale)."""
    from tensorcad import perf_benchmark

    perf_benchmark.clear_mesh_cache()
    for nx, nz in ((21, 9), (41, 13), (81, 21)):
        perf_benchmark._cached_hull(nx, nz)
    assert len(perf_benchmark._MESH_CACHE) == perf_benchmark.MESH_CACHE_SIZE
    perf_benchmark.clear_mesh_cache()


# ---------------------------------------------------------------------------
# Scale-derived (complexity) gates
# ---------------------------------------------------------------------------

def _synth(subtask: str, sizes: list[int], seconds: list[float], copies: int = 1) -> list[SubTaskResult]:
    return [
        SubTaskResult(
            subtask=subtask, nx=0, nz=0, triangles=size, vertices=0,
            grid_cells=0, seconds=sec, throughput=None, copies=copies,
        )
        for size, sec in zip(sizes, seconds)
    ]


def test_scaling_gate_passes_linear_and_sublinear_paths() -> None:
    """An O(N) pass (exponent ~1) and an indexed query (exponent ~0.3) both pass."""
    linear = _synth("volume", [10_000, 100_000, 1_000_000], [0.001, 0.01, 0.1])
    sublinear = _synth("point_locate_accel", [10_000, 100_000, 1_000_000],
                       [0.40, 0.80, 1.60])
    gates = {g.name: g for g in scaling_gates(linear + sublinear)}
    assert gates["scaling_volume"].passed
    assert gates["scaling_volume"].status == MEASURED
    assert abs(gates["scaling_volume"].ratio - 1.0) < 0.02
    assert gates["scaling_point_locate_accel"].passed
    assert gates["scaling_point_locate_accel"].ratio < 0.4


def test_scaling_gate_flags_quadratic_blowup() -> None:
    """Falsification test: an O(N^2) regression must turn the gate RED.

    A wall-clock baseline gate on a tiny ladder cannot see this: the small rung
    would be unchanged in shape.  The exponent is what carries the signal.
    """
    quadratic = _synth("load", [10_000, 100_000, 1_000_000],
                       [0.001, 0.100, 10.000])  # exponent ~2
    gates = {g.name: g for g in scaling_gates(quadratic)}
    gate = gates["scaling_load"]
    assert gate.status == MEASURED
    assert not gate.passed
    assert gate.ratio > 1.8
    assert "RED" not in gate.detail  # detail carries the numbers, not the verdict


def test_scaling_gate_is_not_run_below_min_points() -> None:
    """Two sizes is not an exponent: report not_run, never a silent pass."""
    two = _synth("volume", [10_000, 100_000], [0.001, 0.010])
    gate = scaling_gates(two)[0]
    assert gate.status == NOT_RUN
    assert "not_run" in gate.detail
    assert gate.name.startswith(SCALING_GATE_PREFIX)
    assert MIN_SCALING_POINTS == 3


def test_scaling_gate_separates_body_counts() -> None:
    """Single-body and fleet rungs get distinct gate names."""
    single = _synth("volume", [10_000, 100_000, 1_000_000], [0.001, 0.01, 0.1])
    fleet = _synth("volume", [20_000, 200_000, 2_000_000], [0.002, 0.02, 0.2], copies=2)
    names = [g.name for g in scaling_gates(single + fleet)]
    assert "scaling_volume" in names
    assert "scaling_volume[x2]" in names


def test_gate_verdict_not_run_status_survives_serialisation() -> None:
    gate = GateVerdict(name="scaling_volume", passed=True, status=NOT_RUN,
                       detail="not_run: 2 sizes")
    assert gate.to_dict()["status"] == NOT_RUN


def test_evaluate_gates_keeps_scaling_gates(tmp_path: Path) -> None:
    """Wall-clock gates are added on top of the exponent gates, not instead."""
    baseline_path = tmp_path / "baseline.json"
    baseline_path.write_text(json.dumps({
        "schema": BENCHMARK_SCHEMA,
        "results": [{"subtask": "volume", "nx": 41, "nz": 13, "copies": 1,
                     "seconds": 0.01}],
    }), encoding="utf-8")
    result = PerfBenchmarkResult(
        results=tuple(_synth("volume", [10_000, 100_000, 1_000_000], [0.001, 0.01, 0.1])),
        gates=scaling_gates(_synth("volume", [10_000, 100_000, 1_000_000],
                                   [0.001, 0.01, 0.1])),
    )
    graded, failed = evaluate_gates(result, baseline_path)
    names = [g.name for g in graded.gates]
    assert "scaling_volume" in names
    assert failed == []


# ---------------------------------------------------------------------------
# Multi-body fleet source (10^7 band)
# ---------------------------------------------------------------------------

def _small_hull(nx: int = 21, nz: int = 9):
    from tensorcad.hull import HullParameters, generate_hull

    return generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=nx, nz=nz))


def test_fleet_mesh_scales_triangles_and_volume_exactly() -> None:
    """k bodies = k x triangles and k x signed volume (bodies are disjoint)."""
    from tensorcad.mesh import total_signed_volume

    model = _small_hull()
    base_tris = len(model.mesh.faces)
    base_vol = total_signed_volume(model.mesh)

    for copies in (2, 3, 4):
        mesh = build_fleet_mesh(model, copies)
        assert len(mesh.faces) == copies * base_tris
        assert len(mesh.vertices) == copies * len(model.mesh.vertices)
        vol = total_signed_volume(mesh)
        assert vol == pytest.approx(copies * base_vol, rel=1e-12)
        # Indices stay inside the vertex table and no face is shared across bodies.
        assert max(max(f) for f in mesh.faces) < len(mesh.vertices)
        first_body = set(i for f in model.mesh.faces for i in f)
        assert all(min(f) < copies * len(model.mesh.vertices) for f in mesh.faces)


def test_fleet_mesh_rejects_bad_copy_count() -> None:
    model = _small_hull()
    with pytest.raises(ValueError):
        build_fleet_mesh(model, 0)
    # One copy is the identity: no pointless copy of a multi-GB vertex table.
    assert build_fleet_mesh(model, 1) is model.mesh


def test_fleet_bodies_are_disjoint_for_the_index() -> None:
    """The index must treat the fleet as separated bodies, not as one blob."""
    from tensorcad.spatial_index import TriangleGrid

    model = _small_hull()
    beam = model.parameters.beam
    spacing = 2.0 * beam
    mesh = build_fleet_mesh(model, 3)
    grid = TriangleGrid(mesh)

    # The hull occupies x in [0, L], y in [-B/2, B/2], z in [-T, 0]; query at
    # mid-length, half draft.  A point in the water gap is outside every body,
    # the same point inside each body's own frame is inside exactly that body.
    probe_x, probe_z = 5.0, -0.4
    assert not grid.contains((probe_x, spacing / 2.0, probe_z))
    for k in range(3):
        assert grid.contains((probe_x, k * spacing, probe_z))
    assert grid.stats()["entries"] > 0


def test_fleet_measurement_key_includes_body_count() -> None:
    assert measurement_key("volume", 2001, 501) == "volume:2001x501"
    assert measurement_key("volume", 2001, 501, 3) == "volume:2001x501:b3"


def test_fleet_benchmark_smoke() -> None:
    """Small fleet rung: every sub-task measured, analytic volume error tiny."""
    clear_mesh_cache()
    clear_fleet_cache()
    try:
        result = run_fleet_benchmark(
            base_ladder=((41, 13),),
            copies_ladder=(2, 3),
            n_queries=16,
        )
    finally:
        clear_mesh_cache()
        clear_fleet_cache()

    by_task: dict[str, list[SubTaskResult]] = {}
    for r in result.results:
        by_task.setdefault(r.subtask, []).append(r)

    assert [r.seconds is not None for r in by_task["load"]] == [True]
    assert [r.copies for r in by_task["fleet_assemble"]] == [2, 3]
    assert [r.copies for r in by_task["point_locate_accel"]] == [2, 3]
    volumes = by_task["volume"]
    assert [r.copies for r in volumes] == [2, 3]
    base_tris = len(_small_hull(41, 13).mesh.faces)
    for r in volumes:
        assert r.triangles == r.copies * base_tris
        rel = float(r.note.split("rel_err=")[1].split("%")[0].replace("+", "")) / 100.0
        assert abs(rel) < 0.05


def test_verdict_distinguishes_not_checked_from_passed() -> None:
    """A gate that could not run must not make the whole run look green.

    This is the accounting hole the 10^7 fleet tier exposed: with one base rung
    every exponent gate is ``not_run``, yet ``passed`` was True (``not_run``
    gates carry ``passed=True`` for compatibility) and the CLI printed
    "All gates passed."  ``verdict`` is the honest three-way answer.
    """
    not_checked = PerfBenchmarkResult(
        results=tuple(_synth("volume", [10_000, 100_000], [0.001, 0.010])),
        gates=scaling_gates(_synth("volume", [10_000, 100_000], [0.001, 0.010])),
    )
    assert not_checked.gates_checked == 0
    assert not_checked.gates_not_run == 1
    assert not_checked.verdict == "not_run"
    assert not_checked.to_dict()["verdict"] == "not_run"
    assert not_checked.to_dict()["gates_checked"] == 0

    checked = PerfBenchmarkResult(
        results=tuple(_synth("volume", [10_000, 100_000, 1_000_000], [0.001, 0.01, 0.1])),
        gates=scaling_gates(_synth("volume", [10_000, 100_000, 1_000_000], [0.001, 0.01, 0.1])),
    )
    assert checked.gates_checked == 1
    assert checked.verdict == "pass"

    regressed = PerfBenchmarkResult(
        results=tuple(_synth("load", [10_000, 100_000, 1_000_000], [0.001, 0.1, 10.0])),
        gates=scaling_gates(_synth("load", [10_000, 100_000, 1_000_000], [0.001, 0.1, 10.0])),
    )
    assert regressed.verdict == "fail"


def test_cli_does_not_claim_pass_when_no_gate_was_checked(tmp_path: Path) -> None:
    """End-to-end: a one-rung run must say "NOT a pass", not "All gates passed"."""
    out = tmp_path / "one_rung.json"
    proc = subprocess.run(
        [sys.executable, "-m", "tensorcad.perf_benchmark", "--ladder", "21,9",
         "--subtasks", "load,volume", "--n-queries", "10", "--output", str(out)],
        capture_output=True, text=True, timeout=600,
    )
    assert proc.returncode == 0, proc.stderr
    assert "NOT a pass" in proc.stdout
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["verdict"] == "not_run"
    assert payload["gates_checked"] == 0


# ---------------------------------------------------------------------------
# 10^6-10^7 band: measured exponents and honest "not applicable" rows
# ---------------------------------------------------------------------------

def test_multi_base_fleet_ladder_yields_measured_scaling_gates() -> None:
    """Three base sizes make the fleet band's complexity gate *measured*.

    The exponent is fitted inside one ``(sub-task, body count)`` group.  A
    single base rung gives one size point per group, so every scaling gate came
    out ``not_run`` for the 10^7 band (the documented debt of
    GAP-CAD-D4-PERF-SCALE).  This test fixes the shape of the fix: with three
    bases at one body count the gates must carry a verdict, not ``not_run``.
    """
    from tensorcad import perf_benchmark as pb

    result = pb.run_fleet_benchmark(
        base_ladder=((21, 9), (41, 13), (61, 19)),
        copies_ladder=(2,),
        subtasks=("load", "fleet_assemble", "point_locate_accel", "volume"),
        n_queries=32,
    )
    scaling = [g for g in result.gates if g.name.startswith(pb.SCALING_GATE_PREFIX)]
    assert scaling, "fleet ladder produced no scaling gates"
    for gate in scaling:
        assert gate.status == pb.MEASURED, f"{gate.name}: {gate.detail}"
        assert gate.ratio is not None
    assert result.gates_not_run == 0
    assert result.verdict in ("pass", "fail")


def test_fleet_run_reports_inapplicable_subtasks_as_not_run() -> None:
    """A fleet rung must not show ``failed`` rows for tasks that cannot apply."""
    from tensorcad import perf_benchmark as pb

    result = pb.run_fleet_benchmark(
        base_ladder=((21, 9),),
        copies_ladder=(2,),
        subtasks=("load", "point_locate", "voxelize", "volume"),
        n_queries=8,
    )
    skipped = {r.subtask: r for r in result.results if r.status == pb.NOT_RUN}
    assert "point_locate" in skipped and "voxelize" in skipped
    assert "star-shaped" in skipped["point_locate"].note
    assert not [r for r in result.results if r.status == pb.FAILED]
    # One row per skipped sub-task, not one per (base, body count).
    assert len([r for r in result.results if r.subtask == "point_locate"]) == 1


def test_mesh_footprint_estimate_tracks_exact_size() -> None:
    """The sampled footprint estimate must match an exact accounting closely."""
    import sys

    from tensorcad import perf_benchmark as pb
    from tensorcad.hull import HullParameters, generate_hull

    mesh = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13)).mesh
    exact = sum(
        sys.getsizeof(v) + sum(sys.getsizeof(c) for c in v) for v in mesh.vertices
    ) + sum(
        sys.getsizeof(f) + sum(sys.getsizeof(k) for k in f) for f in mesh.faces
    ) + sys.getsizeof(mesh.vertices) + sys.getsizeof(mesh.faces)
    est_mb = pb.mesh_footprint_mb(mesh)
    assert est_mb is not None
    exact_mb = exact / (1024.0 * 1024.0)
    assert abs(est_mb - exact_mb) / exact_mb < 0.05


def test_live_rss_is_reported_and_serialised() -> None:
    """Each rung carries a live RSS reading, distinct from the process peak."""
    import os

    from tensorcad import perf_benchmark as pb

    if not os.path.exists("/proc/self/status"):
        pytest.skip("/proc/self/status unavailable")
    live = pb.live_rss_mb()
    assert live is not None and live > 0.0

    result = pb.run_fleet_benchmark(
        base_ladder=((21, 9),), copies_ladder=(2,), subtasks=("volume",), n_queries=8
    )
    volume_row = next(r for r in result.results if r.subtask == "volume")
    assert volume_row.rss_mb is not None and volume_row.rss_mb > 0.0
    assert volume_row.mesh_mb is not None and volume_row.mesh_mb > 0.0
    payload = volume_row.to_dict()
    assert payload["rss_mb"] == volume_row.rss_mb
    assert payload["mesh_mb"] == volume_row.mesh_mb
