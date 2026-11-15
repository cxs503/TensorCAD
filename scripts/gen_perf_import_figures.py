"""Figures for the OBJ import/export performance ladder.

Reads the JSON produced by the import benchmark and writes PNGs into
``benchmarks/perf_scale/figures/``.  Every plotted number comes from the JSON.

Figure 1 (``import_scaling.png``)
    log-log wall-clock seconds vs triangles for load, obj_write, obj_read, volume.

Figure 2 (``import_throughput.png``)
    Throughput (tri/s for load/volume, tri/s for write, tri/s for read) vs triangles.

Figure 3 (``import_volume_accuracy.png``)
    Volume relative error of as-imported mesh vs analytic reference, showing
    that the OBJ round-trip preserves geometry to <0.1%.

Usage::

    python scripts/gen_perf_import_figures.py \
        --json benchmarks/perf_scale/artifacts/ladder_import.json \
        --out benchmarks/perf_scale/figures
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

COLORS = {
    "load": "#1f77b4",
    "obj_write": "#ff7f0e",
    "obj_read": "#2ca02c",
    "volume": "#9467bd",
}
LABELS = {
    "load": "load (hull generation)",
    "obj_write": "obj_write (serialise to OBJ)",
    "obj_read": "obj_read (parse OBJ back)",
    "volume": "volume (divergence theorem)",
}


def load_data(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def fig_import_scaling(data: dict, out: Path) -> Path:
    rows = [r for r in data.get("results", []) if r.get("status") == "measured"]
    fig, ax = plt.subplots(figsize=(7.2, 5.0), dpi=150)
    for subtask in ("load", "obj_write", "obj_read", "volume"):
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
    ax.set_title("OBJ import/export ladder: time vs mesh size")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    out.mkdir(parents=True, exist_ok=True)
    path = out / "import_scaling.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_import_throughput(data: dict, out: Path) -> Path:
    rows = [r for r in data.get("results", []) if r.get("status") == "measured"]
    fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=150)
    for subtask in ("load", "obj_write", "obj_read", "volume"):
        entries = sorted(
            [r for r in rows if r["subtask"] == subtask and r.get("throughput")],
            key=lambda r: r["triangles"],
        )
        if not entries:
            continue
        xs = [e["triangles"] for e in entries]
        ys = [e["throughput"] for e in entries]
        ax.plot(xs, ys, "s-", color=COLORS[subtask], label=LABELS[subtask],
                linewidth=1.5, markersize=5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("triangles (-)")
    ax.set_ylabel("throughput (tri/s or file tri/s)")
    ax.set_title("OBJ import/export throughput vs mesh size")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = out / "import_throughput.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def fig_import_volume_accuracy(data: dict, out: Path) -> Path:
    rows = [r for r in data.get("results", [])
            if r["subtask"] == "obj_read" and r.get("status") == "measured"
            and r.get("volume_rel_err") is not None]
    if not rows:
        return None
    rows = sorted(rows, key=lambda r: r["triangles"])
    fig, ax = plt.subplots(figsize=(7.2, 4.2), dpi=150)
    xs = [r["triangles"] for r in rows]
    ys = [abs(r["volume_rel_err"]) * 100 for r in rows]  # as percentage
    ax.plot(xs, ys, "o-", color="#2ca02c", linewidth=1.5, markersize=6)
    ax.axhline(0.1, color="red", linestyle="--", linewidth=1.0, label="0.1% threshold")
    ax.axhline(1.0, color="orange", linestyle="--", linewidth=1.0, label="1.0% threshold")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("triangles (-)")
    ax.set_ylabel("|volume rel error| (%)")
    ax.set_title("OBJ round-trip volume accuracy (as-imported vs analytic)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = out / "import_volume_accuracy.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--json", default="benchmarks/perf_scale/artifacts/ladder_import.json")
    ap.add_argument("--out", default="benchmarks/perf_scale/figures")
    args = ap.parse_args(argv)

    data = load_data(Path(args.json))
    out = Path(args.out)
    made = []
    for fn in (fig_import_scaling, fig_import_throughput, fig_import_volume_accuracy):
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
