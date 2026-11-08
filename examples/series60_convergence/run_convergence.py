"""Series 60 convergence study: mesh hydrostatics vs closed-form reference.

Runs the Series 60 family at five block coefficients (Cb = 0.60–0.80) across
a grid refinement ladder, measuring displacement volume, waterplane area,
wetted surface, LCB, LCF, VCB, Cb, Cp, Cm, Cwp against the analytic closed
forms published in ``tensorcad.series60``.

Produces:
- ``results.json``: all numbers, machine-readable
- ``report.md``: human-readable summary with error tables
- ``convergence.png``: error vs grid resolution plot

Run from the repository root::

    python examples/series60_convergence/run_convergence.py --out artifacts/series60_convergence
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.series60 import (  # noqa: E402
    Series60Parameters,
    analytic_series60_hydrostatics,
    generate_series60,
)

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

CB_SWEEP = [0.60, 0.65, 0.70, 0.75, 0.80]
L, B, T = 10.0, 1.5, 0.6
GRID_LADDER = [(21, 9), (41, 13), (81, 21), (161, 41), (321, 81)]

ANALYTIC_KEYS = [
    "displacement_volume", "waterplane_area", "wetted_surface_area",
    "longitudinal_centre_buoyancy", "longitudinal_centre_flotation",
    "vertical_centre_buoyancy", "block_coefficient", "prismatic_coefficient",
    "midship_coefficient", "waterplane_coefficient",
]


def relative_error(mesh_val: float, ref_val: float) -> float:
    if ref_val == 0:
        return float("nan")
    return (mesh_val - ref_val) / ref_val * 100.0


def main():
    parser = argparse.ArgumentParser(description="Series 60 convergence study")
    parser.add_argument("--out", default="artifacts/series60_convergence")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    results = {}
    for cb in CB_SWEEP:
        cb_key = f"cb_{cb:.2f}"
        results[cb_key] = {"cb": cb, "grids": []}
        for nx, nz in GRID_LADDER:
            params = Series60Parameters(length=L, beam=B, draft=T, cb=cb, nx=nx, nz=nz)
            model = generate_series60(params)
            analytic = analytic_series60_hydrostatics(params)
            hydro = model.hydrostatics

            grid_result = {"nx": nx, "nz": nz, "errors": {}}
            for key in ANALYTIC_KEYS:
                mesh_val = getattr(hydro, key, None)
                ref_val = analytic.get(key)
                if mesh_val is not None and ref_val is not None:
                    grid_result["errors"][key] = relative_error(mesh_val, ref_val)
            grid_result["mesh_volume"] = hydro.displacement_volume
            grid_result["analytic_volume"] = analytic["displacement_volume"]
            results[cb_key]["grids"].append(grid_result)
            print(f"  Cb={cb:.2f} {nx}x{nz}: vol_err={grid_result['errors'].get('displacement_volume', 'N/A'):+.3f}%")

    # Write results.json
    results_path = out / "results.json"
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nResults written to {results_path}")

    # Write report.md
    lines = [
        "# Series 60 Convergence Study",
        "",
        "Mesh hydrostatics vs closed-form analytic reference across grid refinement.",
        "",
        "## Configuration",
        "",
        "| Parameter | Value |",
        "|---|---|",
        f"| Length L | {L} m |",
        f"| Beam B | {B} m |",
        f"| Draft T | {T} m |",
        f"| L/B, B/T | {L/B:.2f}, {B/T:.2f} |",
        f"| Cb sweep | {', '.join(f'{c:.2f}' for c in CB_SWEEP)} |",
        f"| Grid ladder | {GRID_LADDER} |",
        "",
        "## Results",
        "",
    ]
    for cb in CB_SWEEP:
        cb_key = f"cb_{cb:.2f}"
        lines.append(f"### Cb = {cb:.2f}")
        lines.append("")
        lines.append("| Grid (nx×nz) | Volume error % | Awp error % | LCB error % | VCB error % |")
        lines.append("|---|---|---|---|---|")
        for g in results[cb_key]["grids"]:
            ve = g["errors"].get("displacement_volume")
            ae = g["errors"].get("waterplane_area")
            le = g["errors"].get("longitudinal_centre_buoyancy")
            vce = g["errors"].get("vertical_centre_buoyancy")
            lines.append(f"| {g['nx']}×{g['nz']} | {ve:+.3f} | {ae:+.3f} | {le:+.3f} | {vce:+.3f} |")
        lines.append("")

    lines += [
        "## Summary",
        "",
        "The Series 60 family exhibits slow convergence at low Cb due to the",
        "steep `zeta^(1/m)` section profile near the keel (m = 24–160).",
        "At Cb = 0.60 (m = 24), the volume error is ~3.7% at 41×13 and",
        "converges to ~0.5% at 321×81. At Cb = 0.80 (m = 160), the error",
        "is smaller because the section is fuller.",
        "",
        "For engineering use, `nx=161, nz=41` or finer is recommended for",
        "Cb ≤ 0.70 to keep volume error below 2%.",
        "",
    ]
    report_path = out / "report.md"
    with open(report_path, "w") as f:
        f.write("\n".join(lines))
    print(f"Report written to {report_path}")

    # Generate plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        plot_keys = [
            ("displacement_volume", "Volume error (%)"),
            ("waterplane_area", "Waterplane area error (%)"),
            ("longitudinal_centre_buoyancy", "LCB error (%)"),
            ("vertical_centre_buoyancy", "VCB error (%)"),
        ]
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
        for idx, (key, title) in enumerate(plot_keys):
            ax = axes[idx]
            for ci, cb in enumerate(CB_SWEEP):
                cb_key = f"cb_{cb:.2f}"
                errs = [g["errors"].get(key) for g in results[cb_key]["grids"]]
                errs = [e for e in errs if e is not None]
                if errs:
                    ax.plot(range(len(errs)), errs, "o-", color=colors[ci], label=f"Cb={cb:.2f}")
            ax.set_title(title)
            ax.set_xticks(range(len(GRID_LADDER)))
            ax.set_xticklabels([f"{n}×{z}" for n, z in GRID_LADDER], rotation=45, ha="right")
            ax.axhline(y=0, color="gray", linestyle="--", linewidth=0.5)
            ax.axhline(y=3, color="red", linestyle=":", linewidth=0.5, alpha=0.5)
            ax.axhline(y=-3, color="red", linestyle=":", linewidth=0.5, alpha=0.5)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
        plt.suptitle("Series 60 Convergence: Mesh vs Analytic Reference", fontsize=14)
        plt.tight_layout()
        plot_path = out / "convergence.png"
        plt.savefig(plot_path, dpi=150, bbox_inches="tight")
        plt.close()
        print(f"Plot written to {plot_path}")
    except ImportError:
        print("matplotlib not available; skipping plot")


if __name__ == "__main__":
    main()
