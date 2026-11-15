"""Driver for the 10^7-triangle band of GAP-CAD-D4-PERF-SCALE.

Runs the fleet tier (one generated 4.0e6-triangle hull, assembled into 2 and 3
disjoint copies = 8.0e6 and 1.2e7 triangles) and pushes the result out on three
independent channels, so a packaging quirk cannot lose a multi-minute run:

1. ``--output`` JSON at a version-controlled path (default
   ``benchmarks/perf_scale/ladder_fleet_10m.json``),
2. a gzip+base64 payload between ``PERF10M_JSON_BEGIN/END`` markers on stdout
   (recover with ``python -m tensorcad.perf_scale_artifact <log> -o out.json``),
3. a compact human table on stdout.

Exit code 0 only if every measurement is ``measured`` and every volume check is
within 3% of its analytic reference; otherwise 1, with the failure reason in the
payload.  A failed rung is reported as ``status=failed`` - never dropped.

Usage
-----
::

    python scripts/perf_scale_10m.py                       # 2001x501, copies 2,3
    python scripts/perf_scale_10m.py --base 1281:321 \
        --copies 2,4 --n-queries 2000 \
        --output benchmarks/perf_scale/ladder_fleet_10m.json
"""
from __future__ import annotations

import argparse
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tensorcad import perf_benchmark as pb  # noqa: E402
from tensorcad import perf_scale_artifact as art  # noqa: E402

#: Volume-check tolerance (relative).  Same 3% criterion as the report.
VOLUME_TOL = 0.03


def parse_base(text: str) -> tuple[int, int]:
    nx, _, nz = text.partition(":")
    return int(nx), int(nz)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="10^7-band fleet perf run (GAP-CAD-D4-PERF-SCALE)")
    parser.add_argument("--base", default="2001:501", help="base rung nx:nz (default 2001:501 = 4.0e6 tri)")
    parser.add_argument(
        "--bases", default=None,
        help="comma-separated base rungs nx:nz, e.g. '641:161,1281:321,2001:501'. "
             "Needed for a *measurable* complexity exponent at the 10^6-10^7 band: "
             "the exponent is fitted inside one (sub-task, body-count) group, so a "
             "single base gives one size point per group and every scaling gate "
             "comes out not_run (GAP-CAD-D4-PERF-SCALE debt). Overrides --base.",
    )
    parser.add_argument("--copies", default="2,3", help="comma-separated body counts (default 2,3)")
    parser.add_argument("--n-queries", type=int, default=2000, help="ray queries for point_locate_accel")
    parser.add_argument("--subtasks", default="load,fleet_assemble,point_locate_accel,volume")
    parser.add_argument(
        "--output", default="benchmarks/perf_scale/ladder_fleet_10m.json",
        help="destination JSON; keep it inside a version-controlled directory "
             "(artifacts/ sub-paths are filtered out of merges by the release gate)",
    )
    parser.add_argument("--json-out", action="store_true",
                        help="also emit the raw result JSON on stdout (noisy; payload is canonical)")
    args = parser.parse_args(argv)

    if args.bases:
        base_ladder = tuple(parse_base(item) for item in args.bases.split(",") if item.strip())
    else:
        base_ladder = (parse_base(args.base),)
    nx, nz = base_ladder[0]
    copies = tuple(int(c) for c in args.copies.split(",") if c.strip())
    subtasks = tuple(s.strip() for s in args.subtasks.split(",") if s.strip())

    print(f"[perf10m] bases={[f'{a}x{b}' for a, b in base_ladder]} copies={list(copies)} "
          f"n_queries={args.n_queries} subtasks={list(subtasks)}", flush=True)
    print(f"[perf10m] loadavg={pb._loadavg()} cpus={pb._environment().get('cpu_count')}", flush=True)

    t0 = time.perf_counter()
    try:
        result = pb.run_fleet_benchmark(
            base_ladder=base_ladder,
            copies_ladder=copies,
            subtasks=subtasks,
            n_queries=args.n_queries,
        )
        payload = result.to_dict()
    except Exception as exc:  # noqa: BLE001 - a lost run must still be explained
        payload = {
            "schema": pb.BENCHMARK_SCHEMA,
            "passed": False,
            "crash": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc()[-4000:],
            "command": " ".join(sys.argv),
            "environment": pb._environment(),
            "gates": [],
            "results": [],
            "cluster_warnings": [],
            "seconds": round(time.perf_counter() - t0, 3),
        }
        art.emit(payload)
        print(f"PERF10M_CRASH: {payload['crash']}", file=sys.stderr)
        return 2

    payload["command"] = " ".join(sys.argv)
    payload["seconds"] = round(time.perf_counter() - t0, 3)

    # Volume accuracy check per rung: a fast-but-wrong rung must not pass.
    checks = []
    for r in payload["results"]:
        if r["subtask"] != "volume" or r["status"] != pb.MEASURED:
            continue
        note = r.get("note", "")
        rel = None
        for token in note.split(","):
            token = token.strip()
            if token.startswith("rel_err="):
                rel = float(token.split("=", 1)[1].rstrip("%")) / 100.0
        ok = rel is not None and abs(rel) <= VOLUME_TOL
        checks.append({"triangles": r["triangles"], "copies": r["copies"], "rel_err": rel, "ok": ok})
    payload["volume_checks"] = checks

    failures = [c for c in checks if not c["ok"]]
    not_measured = [r for r in payload["results"] if r["status"] != pb.MEASURED]
    payload["passed"] = payload.get("passed", False) and not failures and not not_measured

    art.write_result(payload, args.output)
    art.emit(payload, chunk_chars=art.CHUNK_CHARS)
    print(f"[perf10m] wrote {args.output} in {payload['seconds']:.1f}s", flush=True)

    if failures:
        print(f"PERF10M_VOLUME_FAIL: {failures}", file=sys.stderr)
    if not_measured:
        print(f"PERF10M_NOT_MEASURED: {[r['subtask'] + '@' + str(r['triangles']) for r in not_measured]}",
              file=sys.stderr)
    return 0 if payload["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())