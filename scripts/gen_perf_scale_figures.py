#!/usr/bin/env python3
"""Generate all 12 figures for GAP-CAD-D4-PERF-SCALE from committed artifact JSON.

Reads from benchmarks/perf_scale/artifacts/*.json and
benchmarks/perf_scale/ladder_fleet_10m_3pt.json, writes PNGs to
benchmarks/perf_scale/figures/.

Usage:
    python scripts/gen_perf_scale_figures.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = ROOT / "benchmarks" / "perf_scale" / "artifacts"
#: Figures live under docs/ (not benchmarks/): the push gate admits PNGs only in
#: docs/ (v5.8 F3: 交付图是文档成果非运行记录), and the merge-gate artifact budget
#: (Rendering budget below) keeps the committed total under the per-PR cap.
FIGURES = ROOT / "docs" / "perf_scale" / "figures"

# ---------------------------------------------------------------------------
# Rendering budget
# ---------------------------------------------------------------------------
# The merge gate (owner artifact discipline 2026-11-04, max_artifact_kb=200)
# sums every artifact-class file in a PR -- these 12 PNGs included -- and
# blocks the merge above the cap.  The figures are therefore rendered
# compact (small canvas, modest dpi): small enough to land in git, large
# enough to stay readable and above the 4 KiB validity floor that
# tests/test_perf_scale_evidence.py enforces.
FIGSIZE = (4.4, 2.7)
FIGSIZE_WIDE = (6.2, 2.7)
DPI = 70
#: Per-PR artifact budget of the merge gate, mirrored here so a regeneration
#: that blows the budget fails loudly instead of stranding the next PR.
ARTIFACT_CAP_KB = 200
MIN_FIGURE_BYTES = 4096


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_json(path: Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def fmt(n: float | int) -> str:
    """Format a number for axis labels."""
    if n >= 1e6:
        return f"{n/1e6:.1f}M"
    if n >= 1e3:
        return f"{n/1e3:.0f}k"
    return str(n)


def setup_ax(ax, xlabel: str, ylabel: str, title: str, logx: bool = False, logy: bool = False) -> None:
    import matplotlib
    matplotlib.rcParams.update({
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
    })
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.grid(True, which="major", alpha=0.2, linewidth=0.4)


# ---------------------------------------------------------------------------
# Figure generators
# ---------------------------------------------------------------------------

def gen_scaling(artifact: dict) -> Path:
    """scaling.png: wall-clock vs size log-log for 4 sub-tasks."""
    import matplotlib.pyplot as plt

    results = artifact["results"]
    subtasks = {"load": [], "point_locate": [], "point_locate_accel": [], "volume": []}
    for r in results:
        st = r["subtask"]
        if st in subtasks and r.get("status") == "measured" and r.get("seconds"):
            subtasks[st].append((r["triangles"], r["seconds"]))

    fig, ax = plt.subplots(figsize=FIGSIZE)
    colors = {"load": "#1f77b4", "point_locate": "#ff7f0e",
              "point_locate_accel": "#2ca02c", "volume": "#d62728"}
    markers = {"load": "o", "point_locate": "s", "point_locate_accel": "^", "volume": "D"}
    labels = {"load": "load", "point_locate": "point_locate",
              "point_locate_accel": "point_locate_accel", "volume": "volume"}

    for st, pts in subtasks.items():
        if pts:
            xs, ys = zip(*sorted(pts))
            ax.plot(xs, ys, marker=markers[st], color=colors[st], label=labels[st], linewidth=1.2)

    setup_ax(ax, "Triangle count", "Wall-clock (s)",
             "Scaling: wall-clock vs mesh size (log-log)", logx=True, logy=True)
    ax.legend(loc="upper left")
    ax.set_xticks([1e3, 1e4, 1e5, 1e6, 1e7])
    ax.set_xticklabels(["1k", "10k", "100k", "1M", "10M"])

    out = FIGURES / "scaling.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def gen_point_locate(artifact: dict) -> Path:
    """point_locate.png: brute-force vs accel throughput comparison."""
    import matplotlib.pyplot as plt

    results = artifact["results"]
    naive = []
    accel = []
    for r in results:
        if r.get("status") != "measured" or not r.get("throughput"):
            continue
        if r["subtask"] == "point_locate":
            naive.append((r["triangles"], r["throughput"]))
        elif r["subtask"] == "point_locate_accel":
            accel.append((r["triangles"], r["throughput"]))

    fig, ax = plt.subplots(figsize=FIGSIZE)
    if naive:
        xs, ys = zip(*sorted(naive))
        ax.plot(xs, ys, "s-", color="#ff7f0e", label="brute-force", linewidth=1.2)
    if accel:
        xs, ys = zip(*sorted(accel))
        ax.plot(xs, ys, "^-", color="#2ca02c", label="accel", linewidth=1.2)

    setup_ax(ax, "Triangle count", "Throughput (query/s)",
             "Point locate: brute-force vs accelerated throughput", logx=True, logy=True)
    ax.legend()
    ax.set_xticks([1e3, 1e4, 1e5, 1e6, 1e7])
    ax.set_xticklabels(["1k", "10k", "100k", "1M", "10M"])

    out = FIGURES / "point_locate.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def gen_index_cost(artifact: dict) -> Path:
    """index_cost.png: index build time and size vs triangle count."""
    import matplotlib.pyplot as plt

    results = artifact["results"]
    build_times = []
    index_sizes = []
    for r in results:
        if r["subtask"] != "point_locate_accel" or r.get("status") != "measured":
            continue
        tri = r["triangles"]
        if r.get("build_seconds"):
            build_times.append((tri, r["build_seconds"]))
        if r.get("index_bytes"):
            index_sizes.append((tri, r["index_bytes"] / (1024 * 1024)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIGSIZE_WIDE)

    if build_times:
        xs, ys = zip(*sorted(build_times))
        ax1.plot(xs, ys, "o-", color="#1f77b4", linewidth=1.2)
    setup_ax(ax1, "Triangle count", "Build time (s)",
             "Spatial index build time", logx=True, logy=False)
    ax1.set_xticks([1e3, 1e4, 1e5, 1e6, 1e7])
    ax1.set_xticklabels(["1k", "10k", "100k", "1M", "10M"])

    if index_sizes:
        xs, ys = zip(*sorted(index_sizes))
        ax2.plot(xs, ys, "s-", color="#ff7f0e", linewidth=1.2)
    setup_ax(ax2, "Triangle count", "Index size (MiB)",
             "Spatial index memory footprint", logx=True, logy=False)
    ax2.set_xticks([1e3, 1e4, 1e5, 1e6, 1e7])
    ax2.set_xticklabels(["1k", "10k", "100k", "1M", "10M"])

    fig.suptitle("Spatial index: build cost and memory", fontsize=12, fontweight="bold")
    out = FIGURES / "index_cost.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def gen_voxelize(artifact: dict) -> Path:
    """voxelize.png: voxelisation scaling (32³ → 384³)."""
    import matplotlib.pyplot as plt

    results = artifact["results"]
    pts = [(r["grid_cells"], r["seconds"]) for r in results
           if r["subtask"] == "voxelize" and r.get("status") == "measured" and r.get("seconds")]

    fig, ax = plt.subplots(figsize=FIGSIZE)
    if pts:
        xs, ys = zip(*sorted(pts))
        ax.plot(xs, ys, "o-", color="#9467bd", linewidth=1.2)
        # Add resolution labels
        for tri, sec in pts:
            res = round(tri ** (1/3))
            ax.annotate(f"{res}³", (tri, sec), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=8)

    setup_ax(ax, "Voxel grid cells", "Wall-clock (s)",
             "Voxelisation scaling (SUBOFF SDF)", logx=True, logy=True)

    out = FIGURES / "voxelize.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def gen_fig3pt_scaling(artifact: dict) -> Path:
    """fig3pt_scaling.png: fleet 3-point ladder scaling with fitted exponents."""
    import matplotlib.pyplot as plt

    results = artifact["results"]
    groups: dict[str, list] = {}
    for r in results:
        if r.get("status") != "measured" or not r.get("seconds"):
            continue
        st = r["subtask"]
        groups.setdefault(st, []).append((r["triangles"], r["seconds"]))

    fig, ax = plt.subplots(figsize=FIGSIZE)
    colors = {"load": "#1f77b4", "fleet_assemble": "#ff7f0e",
              "point_locate_accel": "#2ca02c", "volume": "#d62728"}
    labels = {"load": "load", "fleet_assemble": "fleet_assemble",
              "point_locate_accel": "point_locate_accel", "volume": "volume"}

    for st, pts in groups.items():
        if len(pts) >= 2:
            xs, ys = zip(*sorted(pts))
            ax.plot(xs, ys, "o-", color=colors.get(st, "#333"),
                    label=f"{labels.get(st, st)}", linewidth=1.2)

    # Add gate info from artifact
    gates = artifact.get("gates", [])
    gate_text = ""
    for g in gates:
        if g.get("detail"):
            gate_text += f"{g['name']}: {g['detail']}\n"

    setup_ax(ax, "Triangle count", "Wall-clock (s)",
             "10⁶–10⁷ fleet 3-point ladder scaling", logx=True, logy=True)
    ax.legend()
    if gate_text:
        ax.text(0.02, 0.02, gate_text.strip(), transform=ax.transAxes,
                fontsize=7, va="bottom", family="monospace",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

    out = FIGURES / "fig3pt_scaling.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def gen_fig3pt_memory(artifact: dict) -> Path:
    """fig3pt_memory.png: fleet memory (mesh footprint, index, live RSS)."""
    import matplotlib.pyplot as plt

    results = artifact["results"]
    mesh_mb = []
    index_mb = []
    rss_mb = []
    for r in results:
        if r.get("status") != "measured":
            continue
        tri = r["triangles"]
        if r.get("mesh_mb"):
            mesh_mb.append((tri, r["mesh_mb"]))
        if r.get("index_bytes"):
            index_mb.append((tri, r["index_bytes"] / (1024 * 1024)))
        if r.get("rss_mb"):
            rss_mb.append((tri, r["rss_mb"]))

    fig, ax = plt.subplots(figsize=FIGSIZE)
    if mesh_mb:
        xs, ys = zip(*sorted(mesh_mb))
        ax.plot(xs, ys, "o-", color="#1f77b4", label="mesh (MiB)", linewidth=1.2)
    if index_mb:
        xs, ys = zip(*sorted(index_mb))
        ax.plot(xs, ys, "s-", color="#ff7f0e", label="index (MiB)", linewidth=1.2)
    if rss_mb:
        xs, ys = zip(*sorted(rss_mb))
        ax.plot(xs, ys, "^-", color="#2ca02c", label="RSS (MiB)", linewidth=1.2)

    setup_ax(ax, "Triangle count", "Memory (MiB)",
             "10⁶–10⁷ fleet memory breakdown", logx=True, logy=False)
    ax.legend()

    out = FIGURES / "fig3pt_memory.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def gen_import_scaling(artifact: dict) -> Path:
    """import_scaling.png: OBJ import/export wall-clock vs size."""
    import matplotlib.pyplot as plt

    results = artifact["results"]
    groups: dict[str, list] = {}
    for r in results:
        if r.get("status") != "measured" or not r.get("seconds"):
            continue
        groups.setdefault(r["subtask"], []).append((r["triangles"], r["seconds"]))

    fig, ax = plt.subplots(figsize=FIGSIZE)
    colors = {"load": "#1f77b4", "obj_write": "#ff7f0e", "obj_read": "#2ca02c", "volume": "#d62728"}
    for st, pts in groups.items():
        if pts:
            xs, ys = zip(*sorted(pts))
            ax.plot(xs, ys, "o-", color=colors.get(st, "#333"), label=st, linewidth=1.2)

    setup_ax(ax, "Triangle count", "Wall-clock (s)",
             "OBJ import/export scaling", logx=True, logy=True)
    ax.legend()

    out = FIGURES / "import_scaling.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def gen_import_throughput(artifact: dict) -> Path:
    """import_throughput.png: OBJ import/export throughput vs size."""
    import matplotlib.pyplot as plt

    results = artifact["results"]
    groups: dict[str, list] = {}
    for r in results:
        if r.get("status") != "measured" or not r.get("throughput"):
            continue
        groups.setdefault(r["subtask"], []).append((r["triangles"], r["throughput"]))

    fig, ax = plt.subplots(figsize=FIGSIZE)
    colors = {"load": "#1f77b4", "obj_write": "#ff7f0e", "obj_read": "#2ca02c", "volume": "#d62728"}
    for st, pts in groups.items():
        if pts:
            xs, ys = zip(*sorted(pts))
            ax.plot(xs, ys, "o-", color=colors.get(st, "#333"), label=st, linewidth=1.2)

    setup_ax(ax, "Triangle count", "Throughput (tri/s)",
             "OBJ import/export throughput", logx=True, logy=True)
    ax.legend()

    out = FIGURES / "import_throughput.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def gen_import_volume_accuracy(artifact: dict) -> Path:
    """import_volume_accuracy.png: OBJ round-trip volume error vs size."""
    import matplotlib.pyplot as plt

    results = artifact["results"]
    pts = []
    for r in results:
        if r["subtask"] == "volume" and r.get("volume_rel_err") is not None:
            pts.append((r["triangles"], r["volume_rel_err"] * 100))

    fig, ax = plt.subplots(figsize=FIGSIZE)
    if pts:
        xs, ys = zip(*sorted(pts))
        ax.plot(xs, ys, "o-", color="#d62728", linewidth=1.2)
        for tri, err in pts:
            ax.annotate(f"{err:+.3f}%", (tri, err), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=8)

    # 3% threshold line
    ax.axhline(y=-3.0, color="red", linestyle="--", alpha=0.5, label="3% threshold")
    ax.axhline(y=0.0, color="gray", linestyle="-", alpha=0.3)

    setup_ax(ax, "Triangle count", "Volume error (%)",
             "OBJ round-trip volume accuracy (as-imported vs analytic)")
    ax.legend()

    out = FIGURES / "import_volume_accuracy.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def gen_single10m_scaling(artifact: dict) -> Path:
    """single10m_scaling.png: single-body 10⁷ subdivision scaling."""
    import matplotlib.pyplot as plt

    results = artifact["results"]
    groups: dict[str, list] = {}
    for r in results:
        if r.get("status") != "measured" or not r.get("seconds"):
            continue
        groups.setdefault(r["subtask"], []).append((r["triangles"], r["seconds"]))

    fig, ax = plt.subplots(figsize=FIGSIZE)
    colors = {"load": "#1f77b4", "subdivide": "#ff7f0e",
              "point_locate_accel": "#2ca02c", "volume": "#d62728"}
    for st, pts in groups.items():
        if pts:
            xs, ys = zip(*sorted(pts))
            ax.plot(xs, ys, "o-", color=colors.get(st, "#333"), label=st, linewidth=1.2)

    setup_ax(ax, "Triangle count", "Wall-clock (s)",
             "Single-body 10⁷ subdivision scaling", logx=True, logy=True)
    ax.legend()

    out = FIGURES / "single10m_scaling.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def gen_single10m_memory(artifact: dict) -> Path:
    """single10m_memory.png: single-body 10⁷ memory breakdown."""
    import matplotlib.pyplot as plt

    results = artifact["results"]
    mesh_mb = []
    index_mb = []
    rss_mb = []
    for r in results:
        if r.get("status") != "measured":
            continue
        tri = r["triangles"]
        if r.get("mesh_mb"):
            mesh_mb.append((tri, r["mesh_mb"]))
        if r.get("index_bytes"):
            index_mb.append((tri, r["index_bytes"] / (1024 * 1024)))
        if r.get("rss_mb"):
            rss_mb.append((tri, r["rss_mb"]))

    fig, ax = plt.subplots(figsize=FIGSIZE)
    if mesh_mb:
        xs, ys = zip(*sorted(mesh_mb))
        ax.plot(xs, ys, "o-", color="#1f77b4", label="mesh (MiB)", linewidth=1.2)
    if index_mb:
        xs, ys = zip(*sorted(index_mb))
        ax.plot(xs, ys, "s-", color="#ff7f0e", label="index (MiB)", linewidth=1.2)
    if rss_mb:
        xs, ys = zip(*sorted(rss_mb))
        ax.plot(xs, ys, "^-", color="#2ca02c", label="RSS (MiB)", linewidth=1.2)

    setup_ax(ax, "Triangle count", "Memory (MiB)",
             "Single-body 10⁷ memory breakdown", logx=True, logy=False)
    ax.legend()

    out = FIGURES / "single10m_memory.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


def gen_single10m_volume_accuracy(artifact: dict) -> Path:
    """single10m_volume_accuracy.png: subdivision volume conservation."""
    import matplotlib.pyplot as plt

    results = artifact["results"]
    pts = []
    for r in results:
        if r["subtask"] == "volume" and r.get("volume_rel_err") is not None:
            pts.append((r["triangles"], r["volume_rel_err"] * 100))

    fig, ax = plt.subplots(figsize=FIGSIZE)
    if pts:
        xs, ys = zip(*sorted(pts))
        ax.plot(xs, ys, "o-", color="#d62728", linewidth=1.2)
        for tri, err in pts:
            ax.annotate(f"{err:+.4f}%", (tri, err), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=7)

    ax.axhline(y=0.0, color="gray", linestyle="-", alpha=0.3)
    ax.axhline(y=-3.0, color="red", linestyle="--", alpha=0.5, label="3% threshold")

    setup_ax(ax, "Triangle count", "Volume error (%)",
             "Subdivision volume conservation (machine-precision drift)")
    ax.legend()

    out = FIGURES / "single10m_volume_accuracy.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    FIGURES.mkdir(parents=True, exist_ok=True)

    # Load artifacts
    ladder_mesh = load_json(ARTIFACTS / "ladder_mesh.json")
    ladder_voxel = load_json(ARTIFACTS / "ladder_voxel.json")
    ladder_import = load_json(ARTIFACTS / "ladder_import.json")
    ladder_single = load_json(ARTIFACTS / "ladder_single_10m.json")
    ladder_fleet3pt = load_json(ROOT / "benchmarks" / "perf_scale" / "ladder_fleet_10m_3pt.json")

    generators = [
        ("scaling.png", lambda: gen_scaling(ladder_mesh)),
        ("point_locate.png", lambda: gen_point_locate(ladder_mesh)),
        ("index_cost.png", lambda: gen_index_cost(ladder_mesh)),
        ("voxelize.png", lambda: gen_voxelize(ladder_voxel)),
        ("fig3pt_scaling.png", lambda: gen_fig3pt_scaling(ladder_fleet3pt)),
        ("fig3pt_memory.png", lambda: gen_fig3pt_memory(ladder_fleet3pt)),
        ("import_scaling.png", lambda: gen_import_scaling(ladder_import)),
        ("import_throughput.png", lambda: gen_import_throughput(ladder_import)),
        ("import_volume_accuracy.png", lambda: gen_import_volume_accuracy(ladder_import)),
        ("single10m_scaling.png", lambda: gen_single10m_scaling(ladder_single)),
        ("single10m_memory.png", lambda: gen_single10m_memory(ladder_single)),
        ("single10m_volume_accuracy.png", lambda: gen_single10m_volume_accuracy(ladder_single)),
    ]

    for name, gen in generators:
        try:
            path = gen()
            print(f"  ✓ {name} → {path} ({path.stat().st_size:,} bytes)")
        except Exception as e:
            print(f"  ✗ {name}: {e}", file=sys.stderr)
            return 1

    # Merge-gate budget check: the committed figures are artifact-class
    # files, so their total must stay under the cap or the next PR carrying
    # them is held (see Rendering budget above).
    total = 0
    for name, _ in generators:
        size = (FIGURES / name).stat().st_size
        total += size
        if size < MIN_FIGURE_BYTES:
            print(f"  x {name} is {size} B < {MIN_FIGURE_BYTES} B "
                  f"(test floor) -- increase DPI", file=sys.stderr)
            return 1
    if total > ARTIFACT_CAP_KB * 1024:
        print(f"  x figures total {total:,} B > {ARTIFACT_CAP_KB} KiB "
              "merge-gate cap -- reduce FIGSIZE/DPI", file=sys.stderr)
        return 1
    print(f"\nAll 12 figures generated in {FIGURES} "
          f"(total {total:,} B, cap {ARTIFACT_CAP_KB} KiB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
