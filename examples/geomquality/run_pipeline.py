"""Run the geometry-quality pipeline end to end and write the evidence bundle.

    python examples/geomquality/run_pipeline.py --out artifacts/geomquality

Produces
--------
``evidence.json``         every number the report quotes (specimens, exchange, library, seams)
``gate.json``             the machine criteria, their values and the verdict
``<hull>/<hull>.step``     the CAD files the chain was actually run on
``figures/*.png``         geometry, parity, grid convergence, defect matrix, seam separation
``report_table.md``       a Markdown table generated from the JSON, never typed by hand
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from tensorcad.geomquality import baselines as baselines_module
from tensorcad.geomquality import continuity
from tensorcad.geomquality.hydro import (
    convergence_study,
    facet_hydrostatics,
)
from tensorcad.geomquality.pipeline import run_evidence


def figure_geometry(specimens: dict, target: Path) -> Path:
    figure, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    for axis, (name, specimen) in zip(axes, specimens.items()):
        hull = [face for face in specimen.brep.faces if face.name.endswith("hull_half")][0]
        for u in np.linspace(0.0, 1.0, 9):
            points = np.array([hull.surface.evaluate(float(u), v / 30) for v in range(31)])
            axis.plot(points[:, 0], points[:, 2], color="C0", lw=0.7)
        for v in (0.02, 0.2, 0.4, 0.6, 0.8, 1.0):
            points = np.array([hull.surface.evaluate(u / 60, float(v)) for u in range(61)])
            axis.plot(points[:, 0], points[:, 2], color="C3", lw=0.8)
        axis.set_title(f"{name}\nL={specimen.length:g} m, B={specimen.beam:g} m, T={specimen.draft:g} m")
        axis.set_xlabel("x [m]")
        axis.set_ylabel("z [m]")
        axis.set_aspect("equal", adjustable="datalim")
        axis.grid(alpha=0.3)
    figure.suptitle("CAD patch set: stations (red) and buttock lines (blue)")
    figure.tight_layout()
    path = target / "fig1_geometry.png"
    figure.savefig(path, dpi=92)
    plt.close(figure)
    return path


def figure_parity(evidence: dict, target: Path) -> Path:
    quantities = ("volume", "lcb_x", "vcb_z", "waterplane_area")
    hulls = list(evidence["specimens"])
    figure, axis = plt.subplots(figsize=(8.4, 4.2))
    width = 0.2
    floor = 1e-16
    for index, hull in enumerate(hulls):
        comparison = evidence["specimens"][hull]["comparison"]
        patch = [max(comparison[q]["patch_relative_error"], floor) for q in quantities]
        facet = [max(comparison[q]["facet_relative_error"], floor) for q in quantities]
        positions = np.arange(len(quantities)) + index * 2.4 * width
        axis.bar(positions - width / 2, patch, width, label=f"{hull} (patch quadrature)")
        axis.bar(positions + width / 2, facet, width, label=f"{hull} (tessellation)")
    axis.set_yscale("log")
    axis.axhline(1e-3, color="k", ls="--", lw=1.0, label="campaign limit 0.1 %")
    axis.set_xticks(np.arange(len(quantities)) + 2.4 * width)
    axis.set_xticklabels(quantities)
    axis.set_ylabel("relative error vs analytic hull")
    axis.set_title("Hydrostatic parity: exact integration vs the mesh a solver consumes")
    axis.legend(fontsize=7, ncol=2)
    axis.grid(alpha=0.3, axis="y")
    figure.tight_layout()
    path = target / "fig2_parity.png"
    figure.savefig(path, dpi=92)
    plt.close(figure)
    return path


def figure_convergence(specimens: dict, target: Path) -> Path:
    figure, axis = plt.subplots(figsize=(6.4, 4.4))
    for name, specimen in specimens.items():
        reference = specimen.reference["volume"]
        h = []
        error = []
        for u_samples, v_samples in ((12, 6), (24, 12), (48, 24), (96, 48)):
            table = facet_hydrostatics(specimen.brep, u_samples=u_samples, v_samples=v_samples)
            h.append(1.0 / u_samples)
            error.append(abs(table.volume - reference) / reference)
        axis.loglog(h, error, "o-", label=name)
        study = convergence_study(
            specimen.brep, levels=((12, 6), (24, 12), (48, 24)), reference=reference
        )
        apex = error[1]
        axis.annotate(
            f"p = {study['observed_order']:.2f}",
            xy=(h[1], apex),
            xytext=(h[1] * 0.55, apex * 2.6),
            fontsize=8,
            arrowprops={"arrowstyle": "->", "lw": 0.7},
        )
        slope = 2.0
        xs = np.array([h[0], h[-1]])
        axis.plot(xs, error[0] * (xs / h[0]) ** slope, "k:", lw=0.9)
    axis.set_xlabel("1 / u samples  (mesh density)")
    axis.set_ylabel("relative displacement error")
    axis.set_title("Grid convergence of the tessellated displacement\n(dotted: O(h^2) reference slope)")
    axis.grid(alpha=0.3, which="both")
    axis.legend(fontsize=8)
    figure.tight_layout()
    path = target / "fig3_convergence.png"
    figure.savefig(path, dpi=92)
    plt.close(figure)
    return path


def figure_defects(evidence: dict, target: Path) -> Path:
    library = evidence["defect_library"]
    rows = library["rows"]
    codes: list[str] = []
    for row in rows:
        for code in [row["expected_code"], *row["detected_codes"]]:
            if code not in codes:
                codes.append(code)
    matrix = np.zeros((len(rows), len(codes)))
    for i, row in enumerate(rows):
        for j, code in enumerate(codes):
            matrix[i, j] = 1.0 if code in row["detected_codes"] else 0.0
    figure, axis = plt.subplots(figsize=(9.2, 5.6))
    axis.imshow(matrix, cmap="Greys", aspect="auto", vmin=0, vmax=1)
    axis.set_xticks(range(len(codes)))
    axis.set_xticklabels(codes, rotation=45, ha="right", fontsize=7)
    axis.set_yticks(range(len(rows)))
    axis.set_yticklabels([f"{row['case']} -> {row['expected_code']}" for row in rows], fontsize=7)
    for i, row in enumerate(rows):
        j = codes.index(row["expected_code"])
        axis.add_patch(
            plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="C3", lw=1.6)
        )
    axis.set_title(
        f"Defect library on {evidence['defect_reference']}: "
        f"{library['detected']}/{library['cases']} detected, "
        f"false alarms on the clean model: {library['false_alarm_codes'] or 'none'}"
    )
    figure.tight_layout()
    path = target / "fig4_defect_matrix.png"
    figure.savefig(path, dpi=92)
    plt.close(figure)
    return path


def figure_continuity(evidence: dict, specimens: dict, target: Path) -> Path:
    rows = continuity.label_seam_dataset(cases=90)
    thresholds = evidence["continuity"]["thresholds"]
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    colours = {"G2": "C0", "G1": "C1", "G0": "C3"}
    for kind, colour in colours.items():
        subset = [row for row in rows if row["kind"] == kind]
        axes[0].scatter(
            [max(row["max_curvature_jump"], 1e-6) for row in subset],
            [max(row["max_normal_angle_deg"], 1e-4) for row in subset],
            s=14,
            color=colour,
            label=kind,
            alpha=0.8,
        )
    axes[0].axvline(thresholds["curvature_limit"], color="k", ls="--", lw=0.9)
    axes[0].axhline(thresholds["angle_limit_deg"], color="k", ls=":", lw=0.9)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("relative curvature jump")
    axes[0].set_ylabel("normal angle [deg]")
    axes[0].set_title("Seam classification with calibrated limits")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3, which="both")

    for name, specimen in specimens.items():
        face = [item for item in specimen.brep.faces if item.name.endswith("hull_half")][0]
        distribution = continuity.curvature_distribution(face.surface, samples=11)
        axes[1].plot(
            [name, name],
            [distribution["gaussian"]["min"], distribution["gaussian"]["max"]],
            "o-",
        )
        axes[1].plot(
            [name + " (H)", name + " (H)"],
            [distribution["mean"]["min"], distribution["mean"]["max"]],
            "s-",
            color="C3",
        )
    axes[1].set_ylabel("curvature [1/m]")
    axes[1].set_title("Gaussian (blue) and mean (red) curvature range per hull patch")
    axes[1].tick_params(axis="x", labelrotation=15, labelsize=8)
    axes[1].grid(alpha=0.3, axis="y")
    figure.tight_layout()
    path = target / "fig5_continuity.png"
    figure.savefig(path, dpi=92)
    plt.close(figure)
    return path


def write_table(evidence: dict, target: Path) -> Path:
    lines = [
        "# Geometry-quality pipeline: measured numbers",
        "",
        "Generated by `examples/geomquality/run_pipeline.py`; every value is read from",
        "`gate.json` / `evidence.json`, so nothing here is transcribed by hand.",
        "",
        "## Machine criteria",
        "",
        "| criterion | value | limit | verdict |",
        "|---|---:|---:|---|",
    ]
    for item in evidence["criteria"]:
        value = item["value"]
        shown = f"{value:.3e}" if isinstance(value, float) else str(value)
        lines.append(
            f"| {item['criterion']} | {shown} | {item['limit']:.3g} | "
            f"{'PASS' if item['passed'] else 'FAIL'} |"
        )
    lines += [
        "",
        "## Hydrostatics against the analytic hull",
        "",
        "| hull | quantity | patch quadrature | tessellation | patch error | facet error |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for hull, item in evidence["specimens"].items():
        for quantity, row in item["comparison"].items():
            lines.append(
                f"| {hull} | {quantity} | {row['patch']:.9g} | {row['facet']:.9g} | "
                f"{row['patch_relative_error']:.2e} | {row['facet_relative_error']:.2e} |"
            )
    lines += [
        "",
        "## CAD exchange",
        "",
        "| hull | format | watertight in | watertight out | max deviation [m] |",
        "|---|---|---|---:|---:|",
    ]
    for hull, trip in evidence["exchange"].items():
        for suffix, result in trip.items():
            lines.append(
                f"| {hull} | {suffix} | {result['watertight_before']} | "
                f"{result['watertight_after']} | {result['max_deviation']:.2e} |"
            )
    lines += [
        "",
        "## Grid convergence",
        "",
        "| hull | observed order p | GCI (fine) | finest error |",
        "|---|---:|---:|---:|",
    ]
    for hull, item in evidence["specimens"].items():
        study = item["convergence"]
        lines.append(
            f"| {hull} | {study['observed_order']:.3f} | {study['gci_percent']:.3f} % | "
            f"{study['error_finest_percent']:.3e} % |"
        )
    path = target / "report_table.md"
    path.write_text("\n".join(lines) + "\n")
    return path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="artifacts/geomquality")
    parser.add_argument("--cases", type=int, default=150)
    args = parser.parse_args(argv)
    target = Path(args.out)
    figures = target / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    evidence = run_evidence(out_dir=target, cases=args.cases)
    specimens = {name: factory() for name, factory in baselines_module.specimens().items()}
    produced = [
        figure_geometry(specimens, figures),
        figure_parity(evidence, figures),
        figure_convergence(specimens, figures),
        figure_defects(evidence, figures),
        figure_continuity(evidence, specimens, figures),
    ]
    table = write_table(evidence, target)
    summary = {
        "gate": evidence["gate"],
        "criteria": evidence["criteria"],
        "figures": [str(path) for path in produced],
        "table": str(table),
        "defect_library": {
            key: value
            for key, value in evidence["defect_library"].items()
            if key != "rows"
        },
        "continuity": {
            "thresholds": evidence["continuity"]["thresholds"],
            "test": evidence["continuity"]["test"],
            "separation": evidence["continuity"]["separation"],
        },
    }
    (target / "summary.json").write_text(json.dumps(summary, indent=1, default=str))
    print(json.dumps(summary["gate"], indent=1))
    for item in evidence["criteria"]:
        print(
            f"{item['criterion']:32s} {'PASS' if item['passed'] else 'FAIL'} "
            f"value={item['value']!r} limit={item['limit']!r}"
        )
    print(f"figures: {[path.name for path in produced]}")
    return 0 if evidence["gate"]["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
