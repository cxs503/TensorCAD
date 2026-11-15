"""Measure the 10^6-10^7 band on ONE connected body (subdivision source).

The fleet ladder (``scripts/perf_scale_10m.py``) reaches 10^7 triangles by
assembling disjoint copies, which is not the shape an imported solver mesh has.
This driver subdivides a coarser hull instead: one connected closed body per
rung, volume/area preserved by construction, so the divergence-theorem integral
on 10^7 triangles is a *check* rather than a comparison of two discretisations.

Usage::

    # fast mechanism check (seconds)
    python scripts/perf_scale_single10m.py --bases 41:13 --levels 0,1,2 --out /tmp/x.json

    # the real 10^6-10^7 run (minutes; submit through the compute scheduler)
    python scripts/perf_scale_single10m.py \
        --bases 401:129,1601:513 --levels 0,1,2,3 \
        --output benchmarks/perf_scale/artifacts/ladder_single_10m.json

Rungs above ``--max-triangles`` are recorded ``not_run`` with the expected count
instead of being attempted and OOM-killed.
"""
from __future__ import annotations

import argparse
import json
import shlex
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.hull import HullParameters  # noqa: E402
from tensorcad.perf_subdiv import run_single_body_benchmark  # noqa: E402


def _pairs(text: str) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        nx, nz = item.split(":")
        out.append((int(nx), int(nz)))
    if not out:
        raise SystemExit("--bases must list at least one nx:nz pair")
    return out


def _ints(text: str) -> list[int]:
    return [int(piece) for piece in text.split(",") if piece.strip()]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--bases", default="401:129,1601:513",
                        help="comma-separated nx:nz of the parent hulls")
    parser.add_argument("--levels", default="0,1,2,3",
                        help="subdivision levels to measure (comma-separated)")
    parser.add_argument("--length", type=float, default=10.0)
    parser.add_argument("--beam", type=float, default=2.0)
    parser.add_argument("--draft", type=float, default=0.8)
    parser.add_argument("--queries", type=int, default=2000,
                        help="ray count for point_locate_accel")
    parser.add_argument("--max-triangles", type=int, default=16_000_000,
                        help="per-rung cap; larger rungs are recorded not_run")
    parser.add_argument("--output", default=str(
        ROOT / "benchmarks" / "perf_scale" / "artifacts" / "ladder_single_10m.json"))
    args = parser.parse_args(argv)

    artifact = run_single_body_benchmark(
        bases=_pairs(args.bases),
        levels_ladder=_ints(args.levels),
        hull=HullParameters(length=args.length, beam=args.beam, draft=args.draft),
        n_queries=args.queries,
        max_triangles=args.max_triangles,
        command=" ".join(shlex.quote(a) for a in (sys.executable, *sys.argv[:1]))
        + " " + " ".join(shlex.quote(a) for a in (argv if argv is not None else sys.argv[1:])),
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(artifact, indent=2, sort_keys=False) + "\n", encoding="utf-8")

    print(f"{'subtask':<20}{'base':>10}{'lvl':>4}{'triangles':>12}{'seconds':>12}{'status':>10}")
    for row in artifact["results"]:
        seconds = row["seconds"]
        print(
            f"{row['subtask']:<20}{str(row['nx']) + 'x' + str(row['nz']):>10}"
            f"{row['levels']:>4}{row['triangles']:>12,}"
            f"{(f'{seconds:.3f}' if seconds is not None else '-'):>12}"
            f"{row['status']:>10}"
        )
    summary = artifact["summary"]
    print()
    print(f"max triangles          : {summary['max_triangles']:,}")
    print(f"single body >= 1e7     : {summary['single_body_10m_reached']}")
    print(f"invariants ok          : {summary['invariants_ok']}")
    print(f"volume drift vs parent : {summary['max_volume_rel_drift_vs_base']}")
    print(f"area drift vs parent   : {summary['max_surface_area_rel_drift']}")
    print(f"analytic rel err       : {summary['max_analytic_rel_err']}")
    print(f"artifact               : {out}")
    return 0 if summary["invariants_ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())