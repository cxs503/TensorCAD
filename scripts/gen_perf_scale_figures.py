"""Figures for GAP-CAD-D4-PERF-SCALE: 10^6--10^7 triangle performance regime.

Reads the benchmark JSON produced by :mod:`tensorcad.perf_benchmark` and writes
PNGs into ``benchmarks/perf_scale/figures/``.  Every number plotted comes from
the JSON (no hand-typed data).

Usage::

    python scripts/gen_perf_scale_figures.py \
        --mesh artifacts/perf_scale/ladder_mesh.json \
        --voxel artifacts/perf_scale/ladder_voxel.json \
        --out benchmarks/perf_scale/figures
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COLORS = {
    "load": "#1f77b4",
    "point_locate": "#d62728",
    "point_locate_accel": "#2ca02c",
    "volume": "#9467bd",
    "voxelize": "#ff7f0e",
}
LABELS = {
    "load": "load (hull generation)",
    "point_locate": "point_locate (brute force, O(N))",
    "point_locate_accel": "point_locate_accel (grid + DDA)",
    "volume": "volume (divergence theorem)",
    "voxelize": "voxelize (SDF sampling)",
}


def _load(path: str | None) -> dict:
    if not path:
        return {}
    p = Path(path)
    if not p.exists():
        print(f"[skip] {p} not found")
        return {}
    return json.loads(p.read_text())


def _series(data: dict, subtask: str) -> tuple[list[float], list[float], list[float]]:
    xs, secs, tput = [], [], []
    for r in data.get("results", []):
        if r["subtask"] != subtask or r["seconds"] is None:
            continue
        if r["status"] != "measured":
            continue
        x = r["triangles"] if r["triangles"] else r["grid_cells"]
        if not x:
            continue
        xs.append(float(x))
        secs.append(float(r["seconds"]))
        tput.append(float(r["throughput"]) if r["throughput"] else 0.0)
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    return [xs[i] for i in order], [secs[i] for i in order], [tput[i] for i in order]


def fig_scaling(mesh: dict, out: Path) -> Path | None:
    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    plotted = False
    for task in ("load", "point_locate", "point_locate_accel", "volume"):
        xs, secs, _ = _series(mesh, task)
        if not xs:
            continue
        ax.plot(xs, secs, "o-", color=COLORS[task], label=LABELS[task], linewidth=1.6, markersize=4)
        plotted = True
    if not plotted:
        return None
    for scale, name in ((1e6, "10^6"), (1e7, "10^7")):
        ax.axvline(scale, color="grey", linestyle=":", linewidth=1.0)
        ax.text(scale, ax.get_ylim()[1], f" {name}", va="top", fontsize=8, color="grey")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("problem size (triangles for mesh tasks)")
    ax.set_ylabel("wall-clock seconds")
    ax.set_title("TensorCAD geometry sub-task cost, 6.4e2 - 1.0e7 triangles")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = out / "scaling.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def fig_point_locate(mesh: dict, out: Path) -> Path | None:
    """Throughput of brute force vs indexed, plus the speed-up factor."""
    xs_n, _, tp_n = _series(mesh, "point_locate")
    xs_a, _, tp_a = _series(mesh, "point_locate_accel")
    if not xs_a:
        return None
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.6))
    if xs_n:
        ax1.plot(xs_n, tp_n, "o-", color=COLORS["point_locate"], label="brute force (no index)")
    ax1.plot(xs_a, tp_a, "s-", color=COLORS["point_locate_accel"], label="uniform grid + DDA")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("triangles")
    ax1.set_ylabel("ray queries / second")
    ax1.set_title("Point-location throughput")
    ax1.grid(True, which="both", alpha=0.3)
    ax1.legend(fontsize=8)

    # Speed-up where both paths were measured on the same rung.
    naive = {r["triangles"]: r["throughput"] for r in mesh.get("results", [])
             if r["subtask"] == "point_locate" and r["status"] == "measured" and r["throughput"]}
    ax2.plot(xs_a, tp_a, "s-", color=COLORS["point_locate_accel"], label="indexed throughput")
    common_x = [x for x in xs_a if int(x) in naive]
    if common_x:
        speed = [naive[int(x)] for x in common_x]
        ax2.plot(common_x, speed, "o--", color="#7f7f7f", label="brute force (same rungs)")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("triangles")
    ax2.set_ylabel("ray queries / second")
    ax2.set_title("Indexed path reaches the 10^6 band; brute force does not")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=8)
    fig.tight_layout()
    path = out / "point_locate.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def fig_index_cost(mesh: dict, out: Path) -> Path | None:
    """Index build time and memory footprint vs triangles."""
    rows = [r for r in mesh.get("results", [])
            if r["subtask"] == "point_locate_accel" and r["status"] == "measured"
            and r.get("build_seconds") is not None]
    if not rows:
        return None
    rows.sort(key=lambda r: r["triangles"])
    xs = [r["triangles"] for r in rows]
    build = [r["build_seconds"] for r in rows]
    mem = [r["index_bytes"] / 1e6 for r in rows]
    fig, ax1 = plt.subplots(figsize=(8.0, 5.0))
    ax1.plot(xs, build, "o-", color="#2ca02c", label="index build time (s)")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel("triangles")
    ax1.set_ylabel("index build time (s)", color="#2ca02c")
    ax1.tick_params(axis="y", labelcolor="#2ca02c")
    ax1.grid(True, which="both", alpha=0.3)
    ax2 = ax1.twinx()
    ax2.plot(xs, mem, "s--", color="#8c564b", label="index size (MB)")
    ax2.set_ylabel("index size (MB)", color="#8c564b")
    ax2.tick_params(axis="y", labelcolor="#8c564b")
    ax1.set_title("Spatial-index construction cost (one-off) vs mesh size")
    fig.tight_layout()
    path = out / "index_cost.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def fig_voxel(voxel: dict, out: Path) -> Path | None:
    xs, secs, tput = _series(voxel, "voxelize")
    if not xs:
        return None
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    ax.plot(xs, secs, "o-", color=COLORS["voxelize"], label="voxelize wall clock (s)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("grid cells")
    ax.set_ylabel("seconds")
    ax.set_title("SDF voxelization: NumPy path, 3.3e4 - 5.7e7 cells")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = out / "voxelize.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mesh", default="artifacts/perf_scale/ladder_mesh.json")
    ap.add_argument("--voxel", default="artifacts/perf_scale/ladder_voxel.json")
    ap.add_argument("--out", default="benchmarks/perf_scale/figures")
    args = ap.parse_args(argv)

    mesh = _load(args.mesh)
    voxel = _load(args.voxel)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    made = []
    for fn, data in ((fig_scaling, mesh), (fig_point_locate, mesh),
                     (fig_index_cost, mesh), (fig_voxel, voxel)):
        path = fn(data, out)
        if path:
            made.append(path)
            print(f"wrote {path}")
    if not made:
        print("no data; nothing written")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())