"""Figures for the single-body 10^7 triangle subdivision ladder.

Reads the JSON produced by ``scripts/perf_scale_single10m.py`` and writes PNGs
into ``benchmarks/perf_scale/figures/``.  Every plotted number comes from the JSON.

Figure 1 (``single10m_scaling.png``)
    log-log wall-clock seconds vs triangles for load, subdivide, volume,
    point_locate_accel.

Figure 2 (``single10m_memory.png``)
    Mesh footprint, spatial-index size and live RSS vs triangles.

Figure 3 (``single10m_volume_accuracy.png``)
    Volume relative error vs triangles, showing that subdivision preserves
    volume to machine precision.

Usage::

    python scripts/gen_perf_single10m_figures.py \
        --json benchmarks/perf_scale/artifacts/ladder_single_10m.json \
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
    "subdivide": "#8c564b",
    "volume": "#9467bd",
    "point_locate_accel": "#2ca02c",
}
LABELS = {
    "load": "load (hull generation)",
    "subdivide": "subdivide (1-to-4 midpoint)",
    "volume": "volume (divergence theorem)",
    "point_locate_accel": "point_locate_accel (grid + DDA)",
}


def load_data(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def fig_single10m_scaling(data: dict, out: Path) -> Path:
    rows = [r for r in data.get("results", []) if r.get("status") == "measured"]
    fig, ax = plt.subplots(figsize=(7.2, 5.0), dpi=150)
    for subtask in ("load", "subdivide", "volume", "point_locate_accel"):
        entries = sorted(
            [r for r in rows if r["subtask"] == subtask],
            key=lambda r: r["triangles"],
        )
        if not entries:
            continue
        xs = [e["triangles"] for e in entries]
        ys = [e["seconds"] for e in entries]
        ax.plot(xs, ys, "o-", color=COLORS[subtask], label=LABELS[subtask],
                linewidth=1.5, markersize=5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("triangles (-)")
    ax.set_ylabel("wall-clock seconds (-)")
    ax.set_title("Single-body 10^7 triangle ladder: time vs mesh size\n"
                 "(subdivision from a generated parent, not generate_hull)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out.mkdir(parents=True, exist_ok=True)
    path = out / "single10m_scaling.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_single10m_memory(data: dict, out: Path) -> Path:
    rows = [r for r in data.get("results", []) if r.get("status") == "measured"]
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=150)
    seen = set()
    for r in sorted(rows, key=lambda e: e["triangles"]):
        subtask = r["subtask"]
        if r.get("mesh_mb"):
            key = ("mesh", subtask)
            if key not in seen:
                seen.add(key)
                ax.plot(r["triangles"], r["mesh_mb"], "o", ms=6,
                        color=COLORS.get(subtask, "#333"),
                        label=f"mesh footprint ({subtask})")
        if r.get("index_bytes"):
            if ("idx",) not in seen:
                seen.add(("idx",))
                ax.plot(r["triangles"], r["index_bytes"] / (1024.0 * 1024.0), "^", ms=6,
                        color="#ff7f0e", label="spatial index (MB)")
        if r.get("rss_mb"):
            if ("rss",) not in seen:
                seen.add(("rss",))
                ax.plot(r["triangles"], r["rss_mb"], "x", ms=7,
                        color="#7f7f7f", label="live RSS (MiB)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("triangles (-)")
    ax.set_ylabel("MiB (-)")
    ax.set_title("Single-body 10^7: memory vs mesh size")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    path = out / "single10m_memory.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_single10m_volume_accuracy(data: dict, out: Path) -> Path:
    checks = data.get("checks", [])
    if not checks:
        return None
    checks = sorted(checks, key=lambda c: c["triangles"])
    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=150)
    xs = [c["triangles"] for c in checks]
    ys = [abs(c["volume_rel_err"]) * 100 for c in checks]  # as percentage
    ax.plot(xs, ys, "o-", color="#9467bd", linewidth=1.5, markersize=6)
    ax.axhline(0.001, color="red", linestyle="--", linewidth=1.0, label="0.001% threshold")
    ax.axhline(0.01, color="orange", linestyle="--", linewidth=1.0, label="0.01% threshold")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("triangles (-)")
    ax.set_ylabel("|volume rel error| (%)")
    ax.set_title("Single-body 10^7: volume accuracy vs mesh size\n"
                 "(subdivision preserves volume to machine precision)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = out / "single10m_volume_accuracy.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--json", default="benchmarks/perf_scale/artifacts/ladder_single_10m.json")
    ap.add_argument("--out", default="benchmarks/perf_scale/figures")
    args = ap.parse_args(argv)

    data = load_data(Path(args.json))
    out = Path(args.out)
    made = []
    for fn in (fig_single10m_scaling, fig_single10m_memory, fig_single10m_volume_accuracy):
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
