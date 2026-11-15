"""Figures for the single-body 10^6-10^7 benchmark.

Reads ``benchmarks/perf_scale/artifacts/ladder_single_10m.json`` (the same
committed artifact the document tables are rendered from) and writes
``benchmarks/perf_scale/figures/single_body_scaling.png``:

* left  -- seconds vs triangle count, log-log, one connected body per rung;
* right -- point-location query time and index memory vs triangle count.

Nothing here is typed by hand: every plotted point is a measured row.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARTIFACT = ROOT / "benchmarks" / "perf_scale" / "artifacts" / "ladder_single_10m.json"
DEFAULT_OUT = ROOT / "benchmarks" / "perf_scale" / "figures" / "single_body_scaling.png"

COLORS = {
    "load": "#1f77b4",
    "subdivide": "#2ca02c",
    "volume": "#ff7f0e",
    "point_locate_accel": "#d62728",
}


def _series(artifact: dict, subtask: str) -> list[tuple[int, float]]:
    rows = [
        (int(r["triangles"]), float(r["seconds"]))
        for r in artifact["results"]
        if r["subtask"] == subtask
        and r["status"] == "measured"
        and r.get("triangles")
        and r.get("seconds") is not None
    ]
    return sorted(rows)


def _index_series(artifact: dict) -> list[tuple[int, float, float]]:
    rows = [
        (int(r["triangles"]), float(r["seconds"]), (r.get("index_bytes") or 0) / 2**20)
        for r in artifact["results"]
        if r["subtask"] == "point_locate_accel" and r["status"] == "measured"
    ]
    return sorted(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--artifact", default=str(DEFAULT_ARTIFACT))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    artifact = json.loads(Path(args.artifact).read_text(encoding="utf-8"))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.0))

    for subtask in ("load", "subdivide", "volume", "point_locate_accel"):
        pts = _series(artifact, subtask)
        if not pts:
            continue
        xs, ys = zip(*pts)
        ax1.plot(xs, ys, "o-", color=COLORS[subtask], label=subtask, markersize=4)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("triangles (one connected body)")
    ax1.set_ylabel("wall clock [s]")
    ax1.set_title("single-body subdivision ladder")
    ax1.grid(True, which="both", alpha=0.25)
    ax1.legend(fontsize=8)

    idx = _index_series(artifact)
    if idx:
        xs = [r[0] for r in idx]
        query = [r[1] for r in idx]
        index_mb = [r[2] for r in idx]
        ax2.plot(xs, query, "o-", color=COLORS["point_locate_accel"],
                 label="2000 rays [s]", markersize=4)
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_xlabel("triangles")
        ax2.set_ylabel("query batch [s]")
        twin = ax2.twinx()
        twin.plot(xs, index_mb, "s--", color="#555555", label="index [MiB]", markersize=4)
        twin.set_ylabel("index size [MiB]")
        h1, l1 = ax2.get_legend_handles_labels()
        h2, l2 = twin.get_legend_handles_labels()
        ax2.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left")
        ax2.set_title("indexed point location (uniform grid + DDA)")
        ax2.grid(True, which="both", alpha=0.25)

    summary = artifact.get("summary", {})
    fig.suptitle(
        "TensorCAD {}  max {:,.0f} tri, invariants_ok={}".format(
            artifact.get("schema", ""), summary.get("max_triangles") or 0,
            summary.get("invariants_ok"),
        ),
        fontsize=9,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())