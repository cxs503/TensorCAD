#!/usr/bin/env python3
"""Figure for the 10^6--10^7 triangle mesh-derived voxelisation ladder.

Reads ``benchmarks/perf_scale/artifacts/ladder_voxelize_10m.json`` and writes
``benchmarks/perf_scale/figures/voxelize_10m_scaling.png``:

* left  -- wall clock vs triangle count at the fixed 48x24x12 grid (the axis
  that reaches 10^7), with the fitted log-log slope annotated;
* right -- wall clock vs cell count at the fixed 822k-triangle mesh, with the
  fitted slope.

Zero hand-typed numbers: every axis value comes from the artifact.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "benchmarks" / "perf_scale" / "artifacts" / "ladder_voxelize_10m.json"
FIGURE = ROOT / "benchmarks" / "perf_scale" / "figures" / "voxelize_10m_scaling.png"

#: The grid whose triangle axis reaches 10^7 (the finest grid with 4 rungs).
FIXED_GRID = (48, 24, 12)


def _slope(points: list[tuple[float, float]]) -> float:
    xs = [math.log(x) for x, _ in points]
    ys = [math.log(y) for _, y in points]
    n = len(points)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx


def main() -> int:
    if not ARTIFACT.exists():
        print(f"missing artifact {ARTIFACT}", file=sys.stderr)
        return 1
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    measured = [r for r in artifact["results"] if r["status"] == "measured"]
    if not measured:
        print("no measured rungs", file=sys.stderr)
        return 1

    tri_pts = sorted(
        (r["triangles"], r["seconds"])
        for r in measured
        if tuple(r["grid"]) == FIXED_GRID and r["seconds"]
    )
    # fixed mesh = the second-smallest distinct triangle count (the 822k rung)
    distinct = sorted({r["triangles"] for r in measured})
    fixed_mesh = distinct[1] if len(distinct) > 1 else distinct[0]
    cell_pts = sorted(
        (r["grid_cells"], r["seconds"])
        for r in measured
        if r["triangles"] == fixed_mesh and r["seconds"]
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4))
    xs = [p[0] for p in tri_pts]
    ys = [p[1] for p in tri_pts]
    ax1.loglog(xs, ys, "o-", color="tab:blue")
    if len(tri_pts) >= 2:
        ax1.set_title(
            f"wall clock vs triangles (grid {FIXED_GRID[0]}x{FIXED_GRID[1]}x{FIXED_GRID[2]})\n"
            f"slope={_slope(tri_pts):.2f}"
        )
    for x, y in tri_pts:
        ax1.annotate(
            f"{x/1e6:.1f}M", (x, y), textcoords="offset points",
            xytext=(0, 8), ha="center", fontsize=8,
        )
    ax1.set_xlabel("triangles (log)")
    ax1.set_ylabel("seconds (log)")

    ax2.loglog([p[0] for p in cell_pts], [p[1] for p in cell_pts], "s-", color="tab:red")
    if len(cell_pts) >= 2:
        ax2.set_title(
            f"wall clock vs cells (mesh {fixed_mesh:,} tri)\nslope={_slope(cell_pts):.2f}"
        )
    ax2.set_xlabel("voxel cells (log)")
    ax2.set_ylabel("seconds (log)")

    fig.suptitle(
        "Mesh-derived voxelisation, 10^6--10^7 triangle regime (GAP-CAD-D4-PERF-SCALE)"
    )
    fig.tight_layout()
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=140)
    print(f"wrote {FIGURE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
