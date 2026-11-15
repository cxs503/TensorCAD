"""Render the figure embedded in ``docs/perf_benchmark.md`` from committed JSON.

The doc's claim is narrow and quantitative: the naive ray-caster's per-query
cost grows linearly with triangle count, while the uniform-grid path does not,
so the cluster question is algorithmic.  Rather than draw that from memory, this
script reads ``benchmarks/perf_scale/artifacts/ladder_mesh.json`` (the committed
measurement) and plots per-query latency against triangle count, plus the
same-rung speed-up where both paths were measured.

Usage::

    python scripts/gen_perf_benchmark_doc_figure.py           # write the png
    python scripts/gen_perf_benchmark_doc_figure.py --check    # verify it exists

The png lives under ``docs/`` on purpose: the merge gate filters ``artifacts/``,
so a figure referenced from there is a broken link in every clean checkout, and
``docs/`` is where this repository already keeps every committed figure.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "benchmarks" / "perf_scale" / "artifacts" / "ladder_mesh.json"
OUTPUT = ROOT / "docs" / "perf_benchmark" / "figures" / "point_locate_latency.png"


def _latency_us(rows: list[dict], subtask: str) -> tuple[list[int], list[float]]:
    points = [
        (int(r["triangles"]), 1.0e6 / float(r["throughput"]))
        for r in rows
        if r.get("subtask") == subtask
        and r.get("status") == "measured"
        and r.get("throughput")
    ]
    points.sort()
    return [p[0] for p in points], [p[1] for p in points]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--artifact", default=str(ARTIFACT))
    parser.add_argument("--output", default=str(OUTPUT))
    parser.add_argument("--check", action="store_true", help="only verify the png exists")
    args = parser.parse_args(argv)

    output = Path(args.output)
    if args.check:
        if not output.exists() or output.stat().st_size < 4096:
            print(f"missing or truncated figure: {output}")
            return 1
        print(f"ok {output.relative_to(ROOT)} ({output.stat().st_size} bytes)")
        return 0

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = json.loads(Path(args.artifact).read_text(encoding="utf-8"))["results"]
    naive_tri, naive_us = _latency_us(rows, "point_locate")
    accel_tri, accel_us = _latency_us(rows, "point_locate_accel")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ax1.loglog(naive_tri, naive_us, "o-", color="#c0392b", label="brute force O(N) per ray")
    ax1.loglog(accel_tri, accel_us, "s-", color="#1f6fb2", label="uniform grid + DDA (TriangleGrid)")
    ax1.set_xlabel("triangle count")
    ax1.set_ylabel("per-query latency  [us]")
    ax1.set_title("Ray-cast latency vs mesh size (L=10 m hull ladder)")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(loc="upper left", fontsize=9)
    ax1.axhline(60.0e6, ls=":", color="0.4", lw=1.0)
    ax1.annotate("60 s per query", xy=(naive_tri[0], 60.0e6), xytext=(naive_tri[0], 1.2e8),
                 fontsize=8, color="0.35")

    shared = sorted(set(naive_tri) & set(accel_tri))
    lookup = dict(zip(accel_tri, accel_us))
    naive_lookup = dict(zip(naive_tri, naive_us))
    speeds = [naive_lookup[t] / lookup[t] for t in shared]
    ax2.semilogx(shared, speeds, "d-", color="#27ae60")
    ax2.set_xlabel("triangle count")
    ax2.set_ylabel("speed-up (naive / grid)")
    ax2.set_title("Same-query speed-up of the spatial index")
    ax2.grid(True, which="both", alpha=0.3)
    for t, s in zip(shared, speeds):
        if t >= 100_000:
            ax2.annotate(f"{s:.0f}x", xy=(t, s), xytext=(0, 6), textcoords="offset points",
                         fontsize=8, ha="center")

    fig.suptitle(
        "TensorCAD GAP-CAD-D4-PERF-SCALE: point location is an algorithm gap, not a cluster need\n"
        f"(source: {Path(args.artifact).name}, 2000 rays per rung)",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=140)
    print(f"wrote {output.relative_to(ROOT)} ({output.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())