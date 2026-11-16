"""Noise-aware wall-clock gate: estimator, third state, workload comparability.

Every test here is synthetic (injected timings, no sleeping) so it pins the
*judgement* rather than the machine: the gate must (a) still call a sustained
slowdown a RED, (b) refuse to call a straddling measurement a pass, and (c)
refuse to compare timings taken at different workloads.  The measured behaviour
on a loaded box is in ``docs/perf_gate_noise.md``.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tensorcad import perf_benchmark as pb
from tensorcad.perf_benchmark import (
    MEASURED,
    NOISE,
    NOISE_REPEATS,
    NOISE_SPREAD_LIMIT,
    NOT_RUN,
    GateVerdict,
    PerfBenchmarkResult,
    SubTaskResult,
    evaluate_gates,
    measure_with_repeats,
    reduce_seconds,
    scaling_gates,
)

ROOT = Path(__file__).resolve().parents[1]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _result(seconds, samples=(), subtask="load", nx=1000, nz=0, n_queries=0):
    return SubTaskResult(
        subtask=subtask, nx=nx, nz=nz, triangles=1000, vertices=0, grid_cells=0,
        seconds=seconds, n_queries=n_queries, seconds_samples=tuple(samples),
    )


def _baseline_file(tmp_path, entries):
    payload = {
        "schema": "tensorcad.perf_benchmark/2",
        "environment": {"machine": "test", "platform": "test"},
        "results": entries,
    }
    path = tmp_path / "baseline.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def _entry(key_seconds, subtask="load", nx=1000, nz=0, n_queries=0):
    return {
        "subtask": subtask, "nx": nx, "nz": nz, "triangles": 1000, "vertices": 0,
        "grid_cells": 0, "seconds": key_seconds, "status": "measured",
        "n_queries": n_queries,
    }


def _evaluate(results, baseline_path):
    result = PerfBenchmarkResult(
        results=tuple(results),
        gates=(),
        environment={"machine": "test", "platform": "test"},
    )
    return evaluate_gates(result, baseline_path=baseline_path)


# ---------------------------------------------------------------------------
# 1. the estimator
# ---------------------------------------------------------------------------

def test_reduce_seconds_policies():
    samples = [0.4, 0.1, 0.3]
    assert reduce_seconds(samples, "min") == 0.1
    assert reduce_seconds(samples, "max") == 0.4
    assert reduce_seconds(samples, "median") == 0.3


def test_reduce_seconds_rejects_empty_and_unknown():
    with pytest.raises(ValueError):
        reduce_seconds([])
    with pytest.raises(ValueError):
        reduce_seconds([0.1], "average")


def test_seconds_samples_drive_best_worst_spread():
    plain = _result(0.5)
    assert plain.seconds_best == 0.5 and plain.seconds_worst == 0.5
    assert plain.seconds_spread == 1.0, "no repeats must read as zero spread"

    repeated = _result(0.5, samples=(0.5, 0.2, 0.4))
    assert repeated.seconds_best == 0.2
    assert repeated.seconds_worst == 0.5
    assert repeated.seconds_spread == pytest.approx(2.5)
    # seconds keeps its single-shot meaning: the first observation.
    assert repeated.seconds == 0.5


# ---------------------------------------------------------------------------
# 2. repeat policy
# ---------------------------------------------------------------------------

def test_measure_with_repeats_repeats_short_measurements():
    scripted = [0.10, 0.05, 0.07, 0.99]
    calls = []

    def measure():
        value = scripted[len(calls)]
        calls.append(value)
        return _result(value)

    out = measure_with_repeats(measure, repeats=3)
    assert len(calls) == 3, "should stop after repeats=3"
    assert out.seconds == 0.10, "seconds must stay the first observation"
    assert out.seconds_samples == (0.10, 0.05, 0.07)
    assert out.seconds_best == 0.05


def test_measure_with_repeats_skips_long_measurements():
    calls = []

    def measure():
        calls.append(1.0)
        return _result(1.0)

    out = measure_with_repeats(measure, repeats=5, short_seconds=0.5)
    assert calls == [1.0], "a 1.0 s measurement is not a jitter question"
    assert out.seconds_samples == ()


def test_measure_with_repeats_respects_cost_budget():
    calls = []

    def measure():
        calls.append(0.4)
        return _result(0.4)

    # 0.4 s per shot against a 1.0 s budget: shots while the accumulated time is
    # still under budget (0.4, 0.8), then stop - the budget, not repeats=10, ends it.
    out = measure_with_repeats(measure, repeats=10, short_seconds=0.5,
                               repeat_seconds=1.0)
    assert len(calls) == 3, f"budget should stop the loop, got {len(calls)} shots"
    assert out.seconds_samples == (0.4, 0.4, 0.4)


def test_measure_with_repeats_disabled():
    calls = []

    def measure():
        calls.append(0.01)
        return _result(0.01)

    out = measure_with_repeats(measure, repeats=1)
    assert len(calls) == 1
    assert out.seconds_samples == ()


# ---------------------------------------------------------------------------
# 3. gate verdicts: RED survives repetition, ambiguity becomes `noise`
# ---------------------------------------------------------------------------

def test_sustained_regression_is_red_even_with_repeats(tmp_path):
    baseline = _baseline_file(tmp_path, [_entry(0.10)])
    # Every observation is slower: this is a regression, not jitter.
    result = _result(0.33, samples=(0.30, 0.32, 0.33))
    evaluated, failed = _evaluate([result], baseline)
    gate = evaluated.gates[0]
    assert gate.status == MEASURED
    assert not gate.passed
    assert gate.ratio == pytest.approx(3.0)
    assert failed == [gate.name]
    assert evaluated.verdict == "fail"


def test_single_shot_spike_is_still_red(tmp_path):
    baseline = _baseline_file(tmp_path, [_entry(0.10)])
    # Without repeats there is no evidence the spike was noise: RED, as before.
    evaluated, failed = _evaluate([_result(0.26)], baseline)
    gate = evaluated.gates[0]
    assert gate.status == MEASURED and not gate.passed
    assert evaluated.verdict == "fail"


def test_straddling_repeats_are_noise_not_pass(tmp_path):
    baseline = _baseline_file(tmp_path, [_entry(0.10)])
    # Best repeat inside the limit, worst outside, spread above the limit.
    result = _result(0.10, samples=(0.10, 0.11, 0.28))
    evaluated, failed = _evaluate([result], baseline)
    gate = evaluated.gates[0]
    assert gate.status == NOISE, gate.detail
    assert gate.passed, "noise is not a failure"
    assert gate.ratio == pytest.approx(1.0), "ratio reports the best observation"
    assert failed == [], "a noise gate must not be listed as RED"
    assert evaluated.gates_noise == 1
    assert evaluated.verdict == "inconclusive", "must never read as a pass"


def test_consistent_repeats_pass(tmp_path):
    baseline = _baseline_file(tmp_path, [_entry(0.10)])
    result = _result(0.13, samples=(0.13, 0.135, 0.14))
    evaluated, _ = _evaluate([result], baseline)
    gate = evaluated.gates[0]
    assert gate.status == MEASURED and gate.passed
    assert evaluated.verdict == "pass"


def test_spread_below_limit_resolves_as_pass(tmp_path):
    baseline = _baseline_file(tmp_path, [_entry(0.10)])
    # Worst observation is 1.45x the best (< NOISE_SPREAD_LIMIT) and still under
    # the 2x limit: ordinary jitter, decided as a pass.
    result = _result(0.13, samples=(0.13, 0.13, 0.185))
    assert 0.185 / 0.13 < NOISE_SPREAD_LIMIT
    evaluated, _ = _evaluate([result], baseline)
    assert evaluated.gates[0].status == MEASURED
    assert evaluated.gates[0].passed


# ---------------------------------------------------------------------------
# 4. workload comparability (the ray-query count defect)
# ---------------------------------------------------------------------------

def test_query_count_mismatch_is_not_run_not_pass(tmp_path):
    baseline = _baseline_file(
        tmp_path, [_entry(0.05, subtask="point_locate_accel", n_queries=500)]
    )
    result = _result(0.10, subtask="point_locate_accel", n_queries=1000)
    evaluated, failed = _evaluate([result], baseline)
    gate = evaluated.gates[0]
    assert gate.status == NOT_RUN, gate.detail
    assert "workload differs" in gate.detail
    assert failed == []
    assert evaluated.gates_checked == 0
    assert evaluated.verdict == "not_run", "an uncomparable key is not a pass"


def test_query_count_match_is_compared(tmp_path):
    baseline = _baseline_file(
        tmp_path, [_entry(0.05, subtask="point_locate_accel", n_queries=1000)]
    )
    result = _result(0.06, subtask="point_locate_accel", n_queries=1000)
    evaluated, _ = _evaluate([result], baseline)
    gate = evaluated.gates[0]
    assert gate.status == MEASURED and gate.passed
    assert gate.ratio == pytest.approx(1.2)


def test_baseline_without_n_queries_still_compares(tmp_path):
    # v0.3 baselines predate the field: compare, do not silently start skipping.
    entry = _entry(0.05, subtask="point_locate")
    entry.pop("n_queries")
    baseline = _baseline_file(tmp_path, [entry])
    evaluated, _ = _evaluate([_result(0.06, subtask="point_locate", n_queries=1000)], baseline)
    assert evaluated.gates[0].status == MEASURED


# ---------------------------------------------------------------------------
# 5. scaling gate uses the same estimator
# ---------------------------------------------------------------------------

def test_scaling_gate_uses_best_observation():
    sizes = [1000, 10000]
    spike = [
        SubTaskResult(subtask="load", nx=s, nz=0, triangles=s, vertices=0,
                      grid_cells=0, seconds=10.0, seconds_samples=(10.0, 0.1))
        for s in sizes
    ]
    gate = scaling_gates(spike, limits={"load": 5.0}, min_points=2)[0]
    # Both keys measure exactly 0.1 s after the estimator -> exponent 0.
    assert gate.ratio == pytest.approx(0.0, abs=1e-9), gate.detail
    assert gate.passed


# ---------------------------------------------------------------------------
# 6. CLI contract: three exit codes, and the third one says it is not a pass
# ---------------------------------------------------------------------------

def _synthetic_run(results, gates=()):
    return lambda *args, **kwargs: PerfBenchmarkResult(
        results=tuple(results),
        gates=tuple(gates),
        environment={"machine": "test", "platform": "test"},
    )


def test_cli_returns_2_for_noise_only(tmp_path, monkeypatch):
    baseline = _baseline_file(tmp_path, [_entry(0.10)])
    results = [_result(0.10, samples=(0.10, 0.11, 0.28))]
    monkeypatch.setattr(pb, "run_perf_benchmark", _synthetic_run(results))
    code = pb.main(["--quick", "--repeats", "1", "--baseline", str(baseline)])
    monkeypatch.undo()
    assert code == 2, "noise must not exit 0 (false green) or 1 (false RED)"


def test_cli_returns_0_for_pass(tmp_path, monkeypatch):
    baseline = _baseline_file(tmp_path, [_entry(0.10)])
    results = [_result(0.10, samples=(0.10, 0.10, 0.11))]
    monkeypatch.setattr(pb, "run_perf_benchmark", _synthetic_run(results))
    assert pb.main(["--quick", "--repeats", "1", "--baseline", str(baseline)]) == 0


def test_cli_returns_1_for_red(tmp_path, monkeypatch):
    baseline = _baseline_file(tmp_path, [_entry(0.10)])
    monkeypatch.setattr(pb, "run_perf_benchmark", _synthetic_run([_result(0.40)]))
    assert pb.main(["--quick", "--repeats", "1", "--baseline", str(baseline)]) == 1


# ---------------------------------------------------------------------------
# 7. committed artifacts agree with the code
# ---------------------------------------------------------------------------

def test_committed_baseline_uses_the_gate_estimator():
    data = json.loads((ROOT / "benchmarks/perf_scale/baseline.json").read_text(encoding="utf-8"))
    assert data.get("baseline_estimator") == "min", (
        "the stored baseline must use the same estimator as the gate "
        "(see docs/perf_benchmark.md §9)"
    )
    assert "min of" in data["baseline_policy"]


def test_committed_baseline_records_query_counts():
    data = json.loads((ROOT / "benchmarks/perf_scale/baseline.json").read_text(encoding="utf-8"))
    query_rows = [r for r in data["results"] if r["subtask"].startswith("point_locate")]
    assert query_rows, "expected point_locate rows in the baseline"
    for row in query_rows:
        assert row.get("n_queries"), f"{row['subtask']} {row['nx']}x{row['nz']} has no n_queries"
    # And the recorded count matches what the documented CI command runs.
    from tensorcad.perf_benchmark import QUICK_N_QUERIES
    assert {row["n_queries"] for row in query_rows} == {QUICK_N_QUERIES}


def test_gate_noise_document_matches_artifact():
    from scripts.perf_gate_noise_repro import main as repro_main
    assert repro_main(["--check"]) == 0


def test_gate_verdict_status_is_serialised():
    gate = GateVerdict(name="regression_x", passed=True, status=NOISE)
    assert gate.to_dict()["status"] == NOISE
    assert NOISE != MEASURED and NOISE != NOT_RUN