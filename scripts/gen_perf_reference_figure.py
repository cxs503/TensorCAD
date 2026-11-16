"""Figure for ``docs/perf_scale_reference.md`` (GAP-CAD-D4-PERF-SCALE).

Two panels, both read straight from the committed artifacts -- no number is
retyped, so the picture cannot drift from the tables:

  (a) measured wall-clock vs scale (log-log) for every sub-task that has a
      multi-rung ladder, with the 10^6 / 10^7 acceptance band marked;
  (b) deepest measured scale per sub-task, against the same band lines.

Panel (b) is the visual form of the coverage matrix: it shows at a glance that
``load`` (generator) stops at 8.2e6 while the file path ``obj_read`` crosses
1.0e7 -- i.e. "loading at 10^7" exists only through the import path.

Usage::

    python scripts/gen_perf_reference_figure.py
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tensorcad import perf_evidence as pe  # noqa: E402

PERF = ROOT / "benchmarks" / "perf_scale"
OUT = PERF / "figures" / "reference_scaling_bands.png"

ARTIFACTS = {
    "mesh": PERF / "artifacts" / "ladder_mesh.json",
    "voxel": PERF / "artifacts" / "ladder_voxel.json",
    "fleet3pt": PERF / "ladder_fleet_10m_3pt.json",
    "single10m": PERF / "artifacts" / "ladder_single_10m.json",
    "import_single10m": PERF / "artifacts" / "ladder_import_single_10m.json",
    "import_fleet": PERF / "artifacts" / "ladder_import_fleet_10m.json",
}

#: (artifact key, sub-task, label, linestyle) drawn in panel (a).
SERIES = (
    ("mesh", "load", "load (generator, O(N) Python)", "-"),
    ("single10m", "subdivide", "subdivide (1 connected body)", "--"),
    ("mesh", "point_locate_accel", "point_locate_accel (TriangleGrid)", "-"),
    ("fleet3pt", "point_locate_accel", "point_locate_accel ×3 (fleet)", ":"),
    ("mesh", "volume", "volume (divergence theorem)", "-"),
    ("single10m", "volume", "volume @ 10^7 single body", "o"),
    ("import_single10m", "obj_read", "obj_read (file load, 8.2e6)", "s"),
    ("import_fleet", "obj_read", "obj_read (file load, 1.2e7)", "D"),
)


def main() -> int:
    artifacts = {k: pe.load_artifact(p) for k, p in ARTIFACTS.items()}
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13.5, 5.4))

    for key, subtask, label, style in SERIES:
        rows = [
            r for r in pe._measured(artifacts[key])
            if r.get("subtask") == subtask and pe._row_scale(r) and r.get("seconds")
        ]
        if not rows:
            continue
        rows.sort(key=pe._row_scale)
        xs = [pe._row_scale(r) for r in rows]
        ys = [r["seconds"] for r in rows]
        if style in ("o", "s", "D"):
            ax_a.loglog(xs, ys, style, ms=8, label=label)
        else:
            ax_a.loglog(xs, ys, style, lw=2, marker="o", ms=3.5, label=label)

    for x, name in ((1e6, "10^6"), (1e7, "10^7")):
        ax_a.axvline(x, color="crimson", ls="-.", lw=1.0, alpha=0.7)
        ax_a.text(x, ax_a.get_ylim()[0] * 1.6, name, color="crimson", fontsize=8,
                  ha="right", va="bottom")
    ax_a.set_xlabel("mesh size [triangles]")
    ax_a.set_ylabel("wall clock [s]")
    ax_a.set_title("(a) measured scaling, all committed ladders")
    ax_a.grid(True, which="both", alpha=0.25)
    ax_a.legend(fontsize=7.5, loc="upper left")

    # Panel (b): coverage, straight from the rendered matrix.
    body = pe.render_reference_coverage_table(*[artifacts[k] for k in pe.REFERENCE_ARTIFACT_ORDER])
    rows = [
        [c.strip() for c in line.strip("|").split("|")]
        for line in body.splitlines() if line.startswith("| `")
    ]
    labels = [r[0].strip("`") for r in rows]
    scales = [int(r[2].replace(",", "")) for r in rows]
    colors = ["#2a7" if s >= 1e7 else ("#48a" if s >= 1e6 else "#c44") for s in scales]
    ax_b.barh(range(len(labels)), scales, color=colors)
    ax_b.set_yticks(range(len(labels)))
    ax_b.set_yticklabels(labels, fontsize=8)
    ax_b.set_xscale("log")
    for x, name in ((1e6, "10^6"), (1e7, "10^7")):
        ax_b.axvline(x, color="crimson", ls="-.", lw=1.0, alpha=0.7)
        ax_b.text(x, len(labels) - 0.4, name, color="crimson", fontsize=8, ha="right")
    for i, s in enumerate(scales):
        ax_b.text(s * 1.05, i, f"{s/1e6:.2f}M", va="center", fontsize=7)
    ax_b.set_xlabel("deepest measured scale on its own axis (triangles / voxels)")
    ax_b.set_title("(b) coverage vs the 10^6–10^7 acceptance band")
    ax_b.set_xlim(5e4, 1e8)
    ax_b.grid(True, axis="x", which="both", alpha=0.25)

    fig.suptitle("TensorCAD GAP-CAD-D4-PERF-SCALE: measured scaling and band coverage", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=140)
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())