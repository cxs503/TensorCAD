"""Build the committed performance baseline used by the 2x regression gate.

The gate compares a fresh measurement against a stored baseline; a baseline
taken from a single run on a loaded shared box is noisy, so this script runs the
ladder ``--runs`` times and keeps, per (sub-task, resolution) key, the
**slowest** observation.  Taking the max makes the gate tolerant of load noise
(it flags a real 2x regression, not a busy neighbour) and the policy is written
into the JSON so the number is not mistaken for a best case.

Usage::

    python scripts/build_perf_baseline.py                        # quick ladder, 3 runs
    python scripts/build_perf_baseline.py --runs 5 --n-queries 500
    python scripts/build_perf_baseline.py --full \\
        --output benchmarks/perf_scale/baseline_full.json        # 10^6--10^7 rungs
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.perf_benchmark import (  # noqa: E402
    BENCHMARK_SCHEMA,
    DEFAULT_LADDER,
    QUICK_LADDER,
    QUICK_VOXEL_LADDER,
    VOXEL_LADDER,
    _environment,
    clear_mesh_cache,
    run_perf_benchmark,
)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs", type=int, default=3, help="repeat the ladder, keep the slowest")
    ap.add_argument("--n-queries", type=int, default=500)
    ap.add_argument("--full", action="store_true", help="use the full ladder (10^6--10^7)")
    ap.add_argument("--max-triangles", type=int, default=None)
    ap.add_argument("--max-naive-triangles", type=int, default=300_000)
    ap.add_argument("--output", type=str, default="benchmarks/perf_scale/baseline.json")
    args = ap.parse_args(argv)

    ladder = DEFAULT_LADDER if args.full else QUICK_LADDER
    voxel_ladder = VOXEL_LADDER if args.full else QUICK_VOXEL_LADDER

    merged: dict[str, dict] = {}
    for run_index in range(args.runs):
        clear_mesh_cache()
        result = run_perf_benchmark(
            ladder=ladder,
            voxel_ladder=voxel_ladder,
            n_queries=args.n_queries,
            max_triangles=args.max_triangles,
            max_naive_triangles=args.max_naive_triangles,
        )
        for record in result.to_dict()["results"]:
            key = f"{record['subtask']}:{record['nx']}x{record['nz']}"
            if record["seconds"] is None:
                merged.setdefault(key, record)
                continue
            if key in merged and merged[key]["seconds"] is not None:
                if record["seconds"] <= merged[key]["seconds"]:
                    continue
            merged[key] = record
        print(f"run {run_index + 1}/{args.runs}: {len(result.results)} measurements", flush=True)

    payload = {
        "schema": BENCHMARK_SCHEMA,
        "passed": True,
        "seconds": sum(r["seconds"] or 0.0 for r in merged.values()),
        "command": "scripts/build_perf_baseline.py " + " ".join(sys.argv[1:]),
        "environment": _environment(),
        "baseline_policy": (
            f"max of {args.runs} runs per (subtask, nx, nz) key - tolerant of load noise "
            "on a shared machine"
        ),
        "cluster_warnings": [],
        "gates": [],
        "results": [merged[k] for k in sorted(merged)],
    }
    out = Path(args.output)
    if not out.is_absolute():
        out = ROOT / out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"baseline ({len(payload['results'])} records) written to {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())