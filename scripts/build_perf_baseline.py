"""Build the committed performance baseline used by the 2x regression gate.

The gate compares a fresh measurement against a stored baseline; a baseline
taken from a single run on a loaded shared box is noisy, so this script runs the
ladder ``--runs`` times and keeps, per (sub-task, resolution) key, one
observation chosen by ``--policy`` (default ``min``, i.e. the fastest).

**The estimator must match the gate's.**  The gate compares the *best* of its
repeated short measurements against the baseline, because contention only ever
adds seconds.  A baseline recorded as the *max* of N runs is inflated by load
noise by exactly the amount the gate then filters out on the measured side, so a
real 2x regression shows up as ~1.5x and passes - the estimator asymmetry hides
it.  v0.1-v0.3 shipped that asymmetry with a max-of-3 baseline; v0.4 records
min-of-N and says so in ``baseline_policy`` (see ``docs/perf_benchmark.md`` §9).

Usage::

    python scripts/build_perf_baseline.py                        # quick ladder, 3 runs (min)
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
    FLEET_BASE_LADDER,
    FLEET_COPIES_LADDER,
    QUICK_LADDER,
    QUICK_N_QUERIES,
    QUICK_VOXEL_LADDER,
    FULL_N_QUERIES,
    VOXEL_LADDER,
    _environment,
    clear_fleet_cache,
    clear_mesh_cache,
    measurement_key,
    run_fleet_benchmark,
    reduce_seconds,
    run_perf_benchmark,
)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs", type=int, default=3,
                    help="repeat the ladder; --policy picks which observation is kept")
    ap.add_argument("--policy", choices=("min", "max", "median"), default="min",
                    help="estimator for the kept observation per key (default min: the "
                         "fastest, matching the gate's best-of-N; use max only to "
                         "reproduce the pre-0.4 baselines")
    ap.add_argument("--n-queries", type=int, default=None,
                    help="ray queries for point_locate* (default: the same count the "
                         "gate's --quick/--full run uses, so the baseline is comparable "
                         "to the documented CI command)")
    ap.add_argument("--full", action="store_true", help="use the full ladder (10^6--10^7)")
    ap.add_argument("--fleet", action="store_true",
                    help="build the 10^6--10^7 multi-body (fleet) tier baseline "
                         "(cluster-tier project: minutes, several GB)")
    ap.add_argument("--max-triangles", type=int, default=None)
    ap.add_argument("--max-naive-triangles", type=int, default=300_000)
    ap.add_argument("--subtasks", type=str, default=None,
                    help="comma-separated sub-tasks (default: the runner's set; "
                         "use e.g. load,obj_write,obj_read to add the file-import keys)")
    ap.add_argument("--merge", action="store_true",
                    help="keep the seconds already recorded in the existing --output baseline "
                         "and only add keys that are missing; without it the file is replaced")
    ap.add_argument("--output", type=str, default="benchmarks/perf_scale/baseline.json")
    args = ap.parse_args(argv)

    subtasks = None
    if args.subtasks:
        subtasks = tuple(name.strip() for name in args.subtasks.split(",") if name.strip())

    ladder = DEFAULT_LADDER if args.full else QUICK_LADDER
    if args.n_queries is None:
        from tensorcad.perf_benchmark import FULL_N_QUERIES, QUICK_N_QUERIES
        args.n_queries = FULL_N_QUERIES if args.full else QUICK_N_QUERIES
    voxel_ladder = VOXEL_LADDER if args.full else QUICK_VOXEL_LADDER

    collected: dict[str, list[dict]] = {}
    for run_index in range(args.runs):
        clear_mesh_cache()
        clear_fleet_cache()
        if args.fleet:
            result = run_fleet_benchmark(
                base_ladder=FLEET_BASE_LADDER,
                copies_ladder=FLEET_COPIES_LADDER,
                n_queries=args.n_queries,
                max_triangles=args.max_triangles,
            )
        else:
            kwargs = {}
            if subtasks is not None:
                kwargs["subtasks"] = subtasks
            result = run_perf_benchmark(
                ladder=ladder,
                voxel_ladder=voxel_ladder,
                n_queries=args.n_queries,
                max_triangles=args.max_triangles,
                max_naive_triangles=args.max_naive_triangles,
                **kwargs,
            )
        for record in result.to_dict()["results"]:
            key = measurement_key(
                record["subtask"], record["nx"], record["nz"], record.get("copies", 1)
            )
            collected.setdefault(key, []).append(record)
        print(f"run {run_index + 1}/{args.runs}: {len(result.results)} measurements", flush=True)

    # One record per key, chosen by the estimator: a timing is the fastest (or
    # slowest) observation, and the rest of the record travels with it.
    merged: dict[str, dict] = {}
    for key, records in collected.items():
        measured = [r["seconds"] for r in records if r["seconds"] is not None]
        if not measured:
            merged[key] = records[0]
            continue
        target = reduce_seconds(measured, args.policy)
        merged[key] = next(r for r in records if r["seconds"] == target)

    preserved: list[dict] = []
    out_path = Path(args.output)
    if not out_path.is_absolute():
        out_path = ROOT / out_path
    if args.merge and out_path.exists():
        existing = json.loads(out_path.read_text(encoding="utf-8"))
        for record in existing.get("results", []):
            key = measurement_key(
                record["subtask"], record["nx"], record["nz"], record.get("copies", 1)
            )
            if key not in merged:
                preserved.append(record)
        if preserved:
            print(
                f"merge: kept {len(preserved)} existing measurement(s) not re-measured",
                flush=True,
            )

    payload = {
        "schema": BENCHMARK_SCHEMA,
        "passed": True,
        "seconds": sum(r["seconds"] or 0.0 for r in merged.values()),
        "command": "scripts/build_perf_baseline.py " + " ".join(sys.argv[1:]),
        "environment": _environment(),
        "baseline_policy": (
            f"{args.policy} of {args.runs} runs per (subtask, nx, nz) key - the same "
            "estimator the wall-clock gate applies to its own repeats "
            "(docs/perf_benchmark.md §9)"
        ),
        "baseline_estimator": args.policy,
        "cluster_warnings": [],
        "gates": [],
        "results": [merged[k] for k in sorted(merged)] + preserved,
    }
    if args.fleet and args.output == "benchmarks/perf_scale/baseline.json":
        out_path = ROOT / "benchmarks" / "perf_scale" / "baseline_fleet.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if args.merge and preserved:
        payload["baseline_policy"] += (
            f"; --merge kept {len(preserved)} key(s) recorded by an earlier run"
        )
    out_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"baseline ({len(payload['results'])} records) written to {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())