"""Adversarial validation tests for the perf_benchmark framework.

These tests deliberately inject known defects (quadratic blowup, regression,
truncated data) and verify that the framework's gates catch them.  They are
the "negative test" complement to the existing positive tests in
``test_perf_benchmark.py``.

Coverage gaps addressed:
- Scaling gate with *just over* the limit exponent (boundary condition)
- Scaling gate with *just under* the limit exponent (should pass)
- Payload truncation detection in perf_scale_artifact
- Verdict honesty: not_run is not pass
"""
from __future__ import annotations

import json
import tempfile
from pathlib import Path

import pytest

from tensorcad.perf_benchmark import (
    GateVerdict,
    MEASURED,
    NOT_RUN,
    PerfBenchmarkResult,
    SubTaskResult,
    evaluate_gates,
    scaling_gates,
)
from tensorcad.perf_scale_artifact import (
    PayloadError,
    encode_payload,
    parse_emitted,
)


# ---------------------------------------------------------------------------
# Scaling gate boundary conditions
# ---------------------------------------------------------------------------

def _synth(subtask: str, sizes: list[int], seconds: list[float]) -> list[SubTaskResult]:
    return [
        SubTaskResult(
            subtask=subtask, nx=s, nz=0, triangles=s, vertices=0, grid_cells=0,
            seconds=t, status=MEASURED, throughput=s / t if t else None,
        )
        for s, t in zip(sizes, seconds)
    ]


def test_scaling_gate_just_under_limit_is_pass() -> None:
    """Exponent just under limit should pass."""
    # O(N^1.49): just under the 1.50 limit
    k = 1e-6
    sizes = [10_000, 100_000, 1_000_000]
    seconds = [k * s ** 1.49 for s in sizes]
    results = _synth("volume", sizes, seconds)
    gates = {g.name: g for g in scaling_gates(results)}
    gate = gates["scaling_volume"]
    assert gate.passed, f"exponent={gate.ratio} should pass below limit 1.50"
    assert gate.ratio < 1.50


def test_scaling_gate_exactly_at_limit_is_fail() -> None:
    """Exponent == limit should fail (the gate is strictly less-than)."""
    # O(N^1.5) exactly
    k = 1e-6
    sizes = [10_000, 100_000, 1_000_000]
    seconds = [k * s ** 1.5 for s in sizes]
    results = _synth("volume", sizes, seconds)
    gates = {g.name: g for g in scaling_gates(results)}
    gate = gates["scaling_volume"]
    assert not gate.passed, f"exponent={gate.ratio} should fail at limit 1.50 (strict <)"


def test_scaling_gate_just_over_limit_is_fail() -> None:
    """Exponent > limit should fail."""
    # O(N^1.51): just over the 1.50 limit
    k = 1e-6
    sizes = [10_000, 100_000, 1_000_000]
    seconds = [k * s ** 1.51 for s in sizes]
    results = _synth("volume", sizes, seconds)
    gates = {g.name: g for g in scaling_gates(results)}
    gate = gates["scaling_volume"]
    assert not gate.passed, f"exponent={gate.ratio} should fail above limit 1.50"


def test_scaling_gate_accel_limit_is_strict() -> None:
    """point_locate_accel limit is 0.80, not 1.50."""
    k = 1e-3
    sizes = [10_000, 100_000, 1_000_000]
    # O(N^0.81): just over the 0.80 limit
    seconds = [k * s ** 0.81 for s in sizes]
    results = _synth("point_locate_accel", sizes, seconds)
    gates = {g.name: g for g in scaling_gates(results)}
    gate = gates["scaling_point_locate_accel"]
    assert not gate.passed, f"accel exponent={gate.ratio} should fail above 0.80"


def test_scaling_gate_accel_just_under_limit_is_pass() -> None:
    """point_locate_accel just under 0.80 should pass."""
    k = 1e-3
    sizes = [10_000, 100_000, 1_000_000]
    # O(N^0.79): just under the 0.80 limit
    seconds = [k * s ** 0.79 for s in sizes]
    results = _synth("point_locate_accel", sizes, seconds)
    gates = {g.name: g for g in scaling_gates(results)}
    gate = gates["scaling_point_locate_accel"]
    assert gate.passed, f"accel exponent={gate.ratio} should pass below 0.80"


# ---------------------------------------------------------------------------
# Regression gate edge cases
# ---------------------------------------------------------------------------

def test_evaluate_gates_with_no_baseline_returns_empty_gates() -> None:
    """When no baseline file exists, evaluate_gates returns empty gates (no regression check)."""
    results = _synth("load", [1000, 2000], [0.1, 0.2])
    result = PerfBenchmarkResult(
        results=tuple(results),
        gates=(),
        environment={"tensorcad": "0.23.0", "python": "3.11.6"},
    )
    evaluated, failed = evaluate_gates(result, baseline_path="/nonexistent/baseline.json")
    # No baseline means no regression gates added
    assert len(evaluated.gates) == 0, "expected no gates when baseline is missing"


def test_evaluate_gates_regression_at_2x_passes() -> None:
    """Ratio == 2.0 should pass (the gate is less-than-or-equal)."""
    results = _synth("load", [1000], [0.2])  # exactly 2x the baseline 0.1
    result = PerfBenchmarkResult(
        results=tuple(results),
        gates=(),
        environment={"tensorcad": "0.23.0", "python": "3.11.6"},
    )
    baseline_data = {
        "schema": "tensorcad.perf_benchmark/2",
        "results": [
            {"subtask": "load", "nx": 1000, "nz": 0, "triangles": 1000,
             "vertices": 0, "grid_cells": 0, "seconds": 0.1, "status": "measured"}
        ],
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(baseline_data, f)
        baseline_path = f.name

    try:
        evaluated, failed = evaluate_gates(result, baseline_path=baseline_path)
        load_gates = [g for g in evaluated.gates if "load" in g.name]
        assert len(load_gates) == 1
        assert load_gates[0].passed, f"ratio={load_gates[0].ratio} should pass at 2x limit"
    finally:
        Path(baseline_path).unlink(missing_ok=True)


def test_evaluate_gates_regression_over_2x_fails() -> None:
    """Ratio > 2.0 should fail."""
    results = _synth("load", [1000], [0.25])  # 2.5x the baseline 0.1
    result = PerfBenchmarkResult(
        results=tuple(results),
        gates=(),
        environment={"tensorcad": "0.23.0", "python": "3.11.6"},
    )
    baseline_data = {
        "schema": "tensorcad.perf_benchmark/2",
        "results": [
            {"subtask": "load", "nx": 1000, "nz": 0, "triangles": 1000,
             "vertices": 0, "grid_cells": 0, "seconds": 0.1, "status": "measured"}
        ],
    }
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        json.dump(baseline_data, f)
        baseline_path = f.name

    try:
        evaluated, failed = evaluate_gates(result, baseline_path=baseline_path)
        load_gates = [g for g in evaluated.gates if "load" in g.name]
        assert len(load_gates) == 1
        assert not load_gates[0].passed, f"ratio={load_gates[0].ratio} should fail above 2x"
    finally:
        Path(baseline_path).unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# Payload truncation detection
# ---------------------------------------------------------------------------

def test_parse_emitted_detects_missing_chunk() -> None:
    """A payload with a missing chunk should raise PayloadError, not return partial data."""
    payload = {"results": [{"subtask": "load", "seconds": 1.23}]}
    encoded = encode_payload(payload)
    from tensorcad.perf_scale_artifact import chunk_payload, BEGIN_MARKER, END_MARKER, CHUNK_PREFIX
    chunks = chunk_payload(encoded)
    if len(chunks) > 1:
        # Remove the last chunk
        chunks = chunks[:-1]
        total = len(chunks) + 1  # Original total
        lines = [BEGIN_MARKER, f"# chunks={total} chunk_chars=900"]
        lines += [f"{CHUNK_PREFIX} {i + 1}/{total} {chunk}" for i, chunk in enumerate(chunks)]
        lines.append(END_MARKER)
        truncated_text = "\n".join(lines) + "\n"

        with pytest.raises(PayloadError, match="truncated"):
            parse_emitted(truncated_text)


def test_parse_emitted_detects_missing_end_marker() -> None:
    """Missing end marker should raise PayloadError."""
    payload = {"results": []}
    encoded = encode_payload(payload)
    from tensorcad.perf_scale_artifact import BEGIN_MARKER, CHUNK_PREFIX
    lines = [BEGIN_MARKER, f"{CHUNK_PREFIX} 1/1 {encoded}"]
    # No END_MARKER
    text = "\n".join(lines) + "\n"

    with pytest.raises(PayloadError, match="end marker"):
        parse_emitted(text)


def test_parse_emitted_detects_corrupt_payload() -> None:
    """Corrupt base64 should raise PayloadError."""
    from tensorcad.perf_scale_artifact import BEGIN_MARKER, END_MARKER, CHUNK_PREFIX
    lines = [BEGIN_MARKER, f"# chunks=1 chunk_chars=900", f"{CHUNK_PREFIX} 1/1 NOTVALIDBASE64!!!", END_MARKER]
    text = "\n".join(lines) + "\n"

    with pytest.raises(PayloadError):
        parse_emitted(text)


# ---------------------------------------------------------------------------
# Verdict honesty: not_run is not pass
# ---------------------------------------------------------------------------

def test_verdict_not_run_when_all_gates_not_run() -> None:
    """When all gates are not_run, verdict should be 'not_run', not 'pass'."""
    results = _synth("load", [1000], [0.1])
    result = PerfBenchmarkResult(
        results=tuple(results),
        gates=(GateVerdict(name="scaling_load", passed=True, status=NOT_RUN),),
    )
    assert result.verdict == "not_run", f"expected not_run, got {result.verdict}"
    assert result.gates_checked == 0
    assert result.gates_not_run == 1


def test_verdict_fail_when_any_gate_fails() -> None:
    """When any gate fails, verdict should be 'fail'."""
    results = _synth("load", [1000], [0.1])
    result = PerfBenchmarkResult(
        results=tuple(results),
        gates=(
            GateVerdict(name="scaling_load", passed=True, status=MEASURED),
            GateVerdict(name="regression_load", passed=False, status=MEASURED, detail="3.5x > 2.0x"),
        ),
    )
    assert result.verdict == "fail"
    assert result.gates_checked == 2
