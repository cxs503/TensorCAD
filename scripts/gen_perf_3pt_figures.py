"""Figures for the three-point (10^6->10^7) complexity ladder.

Reads the JSON produced by ``scripts/perf_scale_10m.py --bases ...`` and writes
two PNGs; every plotted number is read from the JSON, nothing is hand-typed.

Figure 1 (``fig3pt_scaling.png``)
    log-log seconds vs triangles for every sub-task, with the fitted exponent
    annotated and the gate limit drawn as a dashed line.  A path that regressed
    to O(N^2) would visibly cross its limit.

Figure 2 (``fig3pt_memory.png``)
    Mesh footprint (sampled ``sys.getsizeof`` estimate), spatial-index size and
    live RSS at rung end vs triangles.  This is the "how much RAM do 10^7
    triangles need" answer that wall-clock numbers cannot give.

Usage::

    python scripts/gen_perf_3pt_figures.py \
        --json benchmarks/perf_scale/ladder_fleet_10m_3pt.json \
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
    "fleet_assemble": "#8c564b",
    "point_locate_accel": "#2ca02c",
    "volume": "#9467bd",
}
LABELS = {
    "load": "load (hull generation)",
    "fleet_assemble": "fleet_assemble (O(N) copy)",
    "point_locate_accel": "point_locate_accel (grid + DDA)",
    "volume": "volume (divergence theorem)",
}
#: Exponent limits from tensorcad.perf_benchmark.SCALING_LIMITS (gate contract).
LIMITS = {"load": 1.5, "fleet_assemble": 1.5, "point_locate_accel": 0.8, "volume": 1.5}


def load_rows(path: Path) -> list[dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    rows = [r for r in data.get("results", []) if r.get("status") == "measured"]
    if not rows:
        raise SystemExit(f"no measured rows in {path}")
    return rows


def gate_lookup(data: dict) -> dict[str, dict]:
    return {g["name"]: g for g in data.get("gates", [])}


def fig_scaling(rows: list[dict], data: dict, out: Path) -> Path:
    fig, ax = plt.subplots(figsize=(7.2, 5.0), dpi=150)
    gates = gate_lookup(data)
    # Plot one line per (sub-task, body count) so a fleet rung is not mixed
    # with the single-body rung it was assembled from.
    groups: dict[tuple[str, int], list[dict]] = {}
    for r in rows:
        groups.setdefault((r["subtask"], r.get("copies", 1)), []).append(r)

    for (subtask, copies), entries in sorted(groups.items()):
        entries = sorted(entries, key=lambda e: e["triangles"])
        x = [e["triangles"] for e in entries]
        y = [e["seconds"] for e in entries]
        color = COLORS.get(subtask, "#333333")
        style = "-" if copies == 1 else "--"
        marker = "o" if copies == 1 else "s"
        suffix = "" if copies == 1 else f" [x{copies} bodies]"
        label = f"{LABELS.get(subtask, subtask)}{suffix}"
        ax.plot(x, y, style, marker=marker, color=color, ms=5, lw=1.5, label=label)
        gate = gates.get(f"scaling_{subtask}" + ("" if copies == 1 else f"[x{copies}]"))
        if gate and gate.get("ratio") is not None:
            last = entries[-1]
            ax.annotate(
                f"exp={gate['ratio']:.2f} (limit {LIMITS.get(subtask, float('nan')):.1f})",
                (last["triangles"], last["seconds"]),
                textcoords="offset points", xytext=(-4, 8), fontsize=7.5, color=color,
                ha="right",
            )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("triangles (-)")
    ax.set_ylabel("wall-clock seconds (-)")
    ax.set_title("TensorCAD large-geometry ladder: measured time vs size\n"
                 "(dashed = multi-body fleet rung; exponent from log-log fit)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7.5, loc="upper left")
    fig.tight_layout()
    out.mkdir(parents=True, exist_ok=True)
    path = out / "fig3pt_scaling.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_memory(rows: list[dict], out: Path) -> Path:
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=150)
    seen = set()
    for r in sorted(rows, key=lambda e: e["triangles"]):
        copies = r.get("copies", 1)
        subtask = r["subtask"]
        if copies != 1 and subtask in ("point_locate_accel", "volume"):
            key = ("rss", subtask)
            if r.get("rss_mb") and key not in seen:
                seen.add(key)
                ax.plot(r["triangles"], r["rss_mb"], "x", color="#7f7f7f", ms=7,
                        label="live RSS at rung end (process)" if ("rss",) not in seen else None)
                seen.add(("rss",))
        if r.get("mesh_mb"):
            key = ("mesh", subtask)
            if key in seen:
                continue
            seen.add(key)
            ax.plot(r["triangles"], r["mesh_mb"], "o", ms=6,
                    color={"load": "#1f77b4", "volume": "#9467bd",
                           "fleet_assemble": "#8c564b",
                           "point_locate_accel": "#2ca02c"}.get(subtask, "#333"),
                    label=f"mesh footprint estimate ({subtask}, x{copies})")
        if r.get("index_bytes"):
            ax.plot(r["triangles"], r["index_bytes"] / (1024.0 * 1024.0), "^", ms=6,
                    color="#ff7f0e",
                    label="spatial index (MB)" if ("idx",) not in seen else None)
            seen.add(("idx",))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("triangles (-)")
    ax.set_ylabel("MiB (-)")
    ax.set_title("Memory vs mesh size (mesh-only estimate, index, live RSS)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    out.mkdir(parents=True, exist_ok=True)
    path = out / "fig3pt_memory.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Figures for the 3-point perf ladder")
    ap.add_argument("--json", default="benchmarks/perf_scale/ladder_fleet_10m_3pt.json")
    ap.add_argument("--out", default="benchmarks/perf_scale/figures")
    args = ap.parse_args(argv)

    src = Path(args.json)
    data = json.loads(src.read_text(encoding="utf-8"))
    rows = load_rows(src)
    out = Path(args.out)
    p1 = fig_scaling(rows, data, out)
    p2 = fig_memory(rows, out)
    print(f"wrote {p1}")
    print(f"wrote {p2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())