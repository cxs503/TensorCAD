"""Closure-budget evidence for the IGES/STEP import pipeline (publisher patch sets).

Measures, on a real publisher IGES, how many metres of open boundary are

* a **seam** (a neighbour boundary passes within the budget: sewing closes it),
* a **patch overlap** (the boundary lies on another face's interior: the surface has to
  be trimmed, sewing cannot help),
* a **declared opening** (on a named plane, e.g. a half-hull symmetry cut),
* **unexplained** (nothing within the budget: a real gap or a hole),

both before and after the seam stage, and writes the audit plus one figure.

    python tools/iges_pipe/closure_evidence.py \
        --input /root/tp-agent/assets/geometry/hulls/KCS/KCS_hull_SVA.igs \
        --out artifacts/gap_t1_iges_closure

The numbers in the report are read from the JSON this script writes; nothing is typed
by hand.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.geomquality.trim_mesh import import_patch_set  # noqa: E402

FAMILIES = ("seam", "overlap", "declared_plane", "unexplained")


def _measure(path: Path, *, tolerance: float, seam_tolerance: float | None,
             plane_normal, plane_offset: float, plane_tolerance: float) -> dict:
    planes = [("y0", ((0.0, plane_offset, 0.0), tuple(plane_normal)), plane_tolerance)]
    started = time.time()
    result = import_patch_set(path, closure_tolerance=tolerance, closure_planes=planes,
                              seam_tolerance=seam_tolerance)
    audit = result.report["open_boundary"]
    return {
        "variant": "seam_stage" if seam_tolerance else "as_imported",
        "seam_tolerance": seam_tolerance,
        "closure_tolerance": tolerance,
        "elapsed_seconds": round(time.time() - started, 1),
        "free_edges_after_repair": result.report["after_repair"].get("free_edges", result.report["after_repair"].get("free_edge_count")),
        "watertight": result.report["after_repair"]["watertight"],
        "seam_splits": result.report["seam"].get("splits_total"),
        "open_boundary": audit,
        "family": {name: audit["family"].get(name, {}) for name in FAMILIES},
        "classes": audit["classes"],
        "neighbour_chord_distance": audit["neighbour_chord_distance"],
        "other_surface_distance": audit["other_surface_distance"],
        "unexplained_worst": audit["unexplained_worst"],
    }


def _figure(payload: dict, out_dir: Path) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    runs = payload["runs"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    names = list(FAMILIES)
    width = 0.38
    for index, run in enumerate(runs):
        lengths = [run["family"].get(name, {}).get("length", 0.0) for name in names]
        positions = [i + (index - 0.5) * width for i in range(len(names))]
        axes[0].bar(positions, lengths, width=width, label=run["variant"])
        for position, value in zip(positions, lengths):
            axes[0].text(position, value, f"{value:.0f}", ha="center", va="bottom", fontsize=8)
    axes[0].set_xticks(range(len(names)))
    axes[0].set_xticklabels(names)
    axes[0].set_ylabel("open-boundary length [m]")
    axes[0].set_title("Open boundary by cause")
    axes[0].legend()

    for run in runs:
        audit = run["open_boundary"]
        if not audit.get("unexplained_worst"):
            continue
        lengths = sorted(entry["length"] for entry in audit["unexplained_worst"])
        axes[1].plot(range(len(lengths)), lengths, marker="o", label=run["variant"])
    axes[1].set_xlabel("rank of the longest unexplained chord")
    axes[1].set_ylabel("chord length [m]")
    axes[1].set_title("Longest unexplained chords (witnesses)")
    axes[1].legend()
    fig.tight_layout()
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / "closure_budget.png"
    fig.savefig(target, dpi=130)
    plt.close(fig)
    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", default="artifacts/gap_t1_iges_closure")
    parser.add_argument("--tolerance", type=float, default=5e-3,
                        help="closure coincidence budget in model units")
    parser.add_argument("--seam-tolerance", type=float, default=1.5e-2,
                        help="seam stage budget for the second variant")
    parser.add_argument("--plane-offset", type=float, default=0.0)
    parser.add_argument("--plane-normal", default="0,1,0")
    parser.add_argument("--plane-tolerance", type=float, default=1e-6)
    args = parser.parse_args(argv)

    path = Path(args.input)
    normal = tuple(float(piece) for piece in args.plane_normal.split(","))
    out_dir = Path(args.out)
    runs = [_measure(path, tolerance=args.tolerance, seam_tolerance=None,
                     plane_normal=normal, plane_offset=args.plane_offset,
                     plane_tolerance=args.plane_tolerance),
            _measure(path, tolerance=args.tolerance, seam_tolerance=args.seam_tolerance,
                     plane_normal=normal, plane_offset=args.plane_offset,
                     plane_tolerance=args.plane_tolerance)]
    payload = {
        "tool": "tools/iges_pipe/closure_evidence.py",
        "input": str(path),
        "bytes": path.stat().st_size,
        "runs": runs,
    }
    figure = _figure(payload, out_dir / "figs")
    payload["figure"] = str(figure)
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / f"{path.stem}_closure.json"
    target.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(json.dumps({run["variant"]: run["family"] for run in runs}, indent=2))
    print(f"written: {target}  figure: {figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())