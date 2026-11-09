"""Regenerate the IN383-CAD01 evidence: artifact JSON, figures and result tables.

Run from the repository root::

    python tools/hull_nurbs_report.py

Everything this script prints or plots is measured in the same run from
:mod:`tensorcad.hull_nurbs`; no number is transcribed. It writes

* ``artifacts/in383_cad01/acceptance.json``   -- the machine-readable acceptance record
* ``artifacts/in383_cad01/figures/*.png``     -- the figures embedded in the report
* ``docs/experiments/in383-cad01-hull-nurbs/results.md`` -- the result tables
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tensorcad.hull_nurbs import (  # noqa: E402
    SCALED_JACOBIAN_GOOD,
    SCALED_JACOBIAN_MIN,
    acceptance_report,
    build_body_mesh,
    build_nurbs_hull,
    hex_scaled_jacobian,
    hex_to_tets,
    holtrop_wetted_area,
    mesh_integrals,
    nurbs_hydrostatics,
    schneekluth_lcb_ratio,
)

ARTIFACT_DIR = ROOT / "artifacts" / "in383_cad01"
FIGURE_DIR = ARTIFACT_DIR / "figures"
DOC_DIR = ROOT / "docs" / "experiments" / "in383-cad01-hull-nurbs"


def relative(value: float, reference: float) -> float:
    return abs(value - reference) / reference


# ---------------------------------------------------------------------------
# measurements
# ---------------------------------------------------------------------------
def measure() -> dict:
    """Run the acceptance harness plus the extra measurements the report quotes."""
    report = acceptance_report()

    hull = build_nurbs_hull()
    hydro = nurbs_hydrostatics(hull)
    p = hull.parameters

    ladder = []
    for stations, levels, above, layers in (
        (16, 6, 2, 4),
        (24, 9, 3, 5),
        (32, 12, 4, 6),
        (64, 24, 8, 12),
    ):
        mesh = build_body_mesh(
            hull,
            stations=stations,
            levels_below_waterline=levels,
            levels_above_waterline=above,
            layers=layers,
        )
        integrals = mesh_integrals(mesh)
        ladder.append(
            {
                "stations": stations,
                "levels_below_waterline": levels,
                "levels_above_waterline": above,
                "layers": layers,
                "cells": len(mesh.hexes),
                "nodes": len(mesh.nodes),
                "displaced_volume_m3": integrals.displaced_volume,
                "wetted_surface_area_m2": integrals.wetted_surface_area,
                "vertical_centre_of_buoyancy_m": integrals.vertical_centre_of_buoyancy,
                "volume_relative_error": relative(integrals.displaced_volume, hydro.displaced_volume),
                "wetted_area_relative_error": relative(
                    integrals.wetted_surface_area, hydro.wetted_surface_area
                ),
                "pass_fraction": mesh.quality.pass_fraction,
                "collapsed_cells": mesh.quality.collapsed,
            }
        )

    # Section-wise anatomy: local block coefficient along the hull, NURBS against mesh.
    stations, nz, layers = 32, 12, 6
    mesh = build_body_mesh(
        hull, stations=stations, levels_below_waterline=nz, levels_above_waterline=4, layers=layers
    )
    sections = []
    for i in range(0, stations + 1, 2):
        xi = mesh.node(i, 0, 0)[0] / p.length
        u = min(xi, 0.999999)
        area = 0.0
        samples = 200
        for k in range(samples):
            zeta = (k + 0.5) / samples
            area += 2.0 * hull.patch_half_breadth(u, zeta) * p.draft / samples
        sections.append({"x_over_l": xi, "section_area_m2": area})

    # Quality distributions for the histogram figure.
    mesh_default = build_body_mesh(hull)
    sj = []
    ar = []
    for cell in mesh_default.hexes:
        points = [mesh_default.nodes[index] for index in cell]
        value = hex_scaled_jacobian(points)
        sj.append(value)
        if value > 0.0:
            short, long = 0.0, 0.0
            edges = (
                (0, 1), (1, 2), (2, 3), (3, 0),
                (4, 5), (5, 6), (6, 7), (7, 4),
                (0, 4), (1, 5), (2, 6), (3, 7),
            )
            lengths = [
                math.dist(points[b], points[a]) for a, b in edges
            ]
            short, long = min(lengths), max(lengths)
            ar.append(long / short)

    tets, tet_report = hex_to_tets(mesh_default)

    return {
        "acceptance": report,
        "ladder": ladder,
        "sections": sections,
        "scaled_jacobians": sj,
        "aspect_ratios": ar,
        "tetrahedra": tet_report,
        "tet_count": len(tets),
        "baselines": {
            "holtrop_wetted_area_m2": holtrop_wetted_area(
                p,
                midship_section_coefficient=hydro.midship_section_coefficient,
                waterplane_coefficient=hydro.waterplane_coefficient,
            ),
            "schneekluth_lcb_ratio": schneekluth_lcb_ratio(p),
        },
    }


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def figure_lines(hull, sections, path: Path) -> None:
    p = hull.parameters
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6))

    # body plan: sections
    ax = axes[0]
    for station in range(0, 21):
        xi = station / 20
        xs = [hull.patch_half_breadth(min(xi, 0.999999), k / 40) for k in range(41)]
        zs = [p.draft * k / 40 for k in range(41)]
        ax.plot(xs, zs, color=plt.cm.viridis(xi), lw=1.1)
        ax.plot([-x for x in xs], zs, color=plt.cm.viridis(xi), lw=1.1)
    ax.set_xlabel("y [m]")
    ax.set_ylabel("z [m]")
    ax.set_title("Body plan: 21 NURBS sections (starboard and mirror)")
    ax.set_aspect("equal")
    ax.grid(alpha=0.25)

    # waterlines and buttocks
    ax = axes[1]
    for level in range(0, 6):
        zeta = level / 4
        xs = [i / 200 for i in range(201)]
        ys = [
            (hull.patch_half_breadth(min(xi, 0.999999), min(zeta, 1.0)), p.length * xi) for xi in xs
        ]
        ax.plot([q[0] for q in ys], [q[1] for q in ys], label=f"z/T = {zeta:.2f}")
    ax.set_xlabel("y [m]")
    ax.set_ylabel("x [m]")
    ax.set_title("Waterlines of the NURBS patch")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def figure_quality(values, aspect_ratios, quality, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ax = axes[0]
    ax.hist(values, bins=50, color="#3b6fa0")
    ax.axvline(SCALED_JACOBIAN_MIN, color="#c0392b", ls="--", label=f"admission {SCALED_JACOBIAN_MIN}")
    ax.axvline(SCALED_JACOBIAN_GOOD, color="#e67e22", ls=":", label=f"good {SCALED_JACOBIAN_GOOD}")
    ax.set_xlabel("min corner scaled Jacobian per hex")
    ax.set_ylabel("cells")
    ax.set_title(
        f"Skew: pass {quality.pass_fraction * 100:.2f} % of {quality.cells} cells "
        f"({quality.collapsed} collapsed)"
    )
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.hist([min(value, 1000.0) for value in aspect_ratios], bins=50, color="#4d8354")
    ax.axvline(100.0, color="#c0392b", ls="--", label="cap 100")
    ax.set_xlabel("edge-length aspect ratio (clipped at 1000)")
    ax.set_ylabel("cells")
    ax.set_title(f"Aspect ratio: max {quality.max_aspect_ratio:.0f} at the stem wedge")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def figure_convergence(ladder, hydro, path: Path) -> None:
    cells = [item["cells"] for item in ladder]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ax = axes[0]
    ax.loglog(cells, [item["volume_relative_error"] * 100 for item in ladder], "o-", label="displacement")
    ax.loglog(
        cells,
        [item["wetted_area_relative_error"] * 100 for item in ladder],
        "s-",
        label="wetted area",
    )
    ax.axhline(1.0, color="#c0392b", ls="--", label="1 % consistency band")
    ax.set_xlabel("hex cells")
    ax.set_ylabel("relative error vs NURBS quadrature [%]")
    ax.set_title("Mesh -> patch convergence")
    ax.grid(which="both", alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.semilogx(cells, [item["pass_fraction"] * 100 for item in ladder], "^-", color="#8e44ad")
    ax.axhline(90.0, color="#c0392b", ls="--", label="C4 gate 90 %")
    ax.set_xlabel("hex cells")
    ax.set_ylabel("cells passing the quality gate [%]")
    ax.set_title("Quality against resolution")
    ax.grid(which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def figure_comparison(report, hydro, baselines, path: Path) -> None:
    """Three routes to the same three numbers, with their bands."""
    integrals = report["mesh_integrals"]
    p = report["parameters"]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.0))

    ax = axes[0]
    labels = ["design\nCb·L·B·T", "NURBS\nquadrature", "hex mesh"]
    values = [p["design_displacement_m3"], hydro.displaced_volume, integrals["displaced_volume_m3"]]
    ax.bar(labels, values, color=["#7f8c8d", "#3b6fa0", "#4d8354"])
    for index, value in enumerate(values):
        ax.text(index, value * 1.002, f"{value:.5f}", ha="center", fontsize=8)
    ax.set_ylabel("displacement [m³]")
    ax.set_title(f"C1: {(values[1] / values[0] - 1) * 100:+.3f} % vs design")
    ax.set_ylim(0, max(values) * 1.12)

    ax = axes[1]
    baseline = baselines["holtrop_wetted_area_m2"]
    values = [baseline, hydro.wetted_surface_area, integrals["wetted_surface_area_m2"]]
    ax.bar(["Holtrop-\nMennen", "NURBS", "hex mesh"], values, color=["#7f8c8d", "#3b6fa0", "#4d8354"])
    for index, value in enumerate(values):
        ax.text(index, value * 1.002, f"{value:.4f}", ha="center", fontsize=8)
    ax.set_ylabel("wetted area [m²]")
    ax.set_title(f"C2: {(values[1] / values[0] - 1) * 100:+.2f} % vs empirical")
    ax.set_ylim(0, max(values) * 1.12)

    ax = axes[2]
    baseline = baselines["schneekluth_lcb_ratio"] * 100
    values = [baseline, hydro.lcb_ratio_from_midship * 100, integrals["lcb_ratio_from_midship"] * 100]
    ax.bar(["Schneekluth", "NURBS", "hex mesh"], values, color=["#7f8c8d", "#3b6fa0", "#4d8354"])
    for index, value in enumerate(values):
        ax.text(index, value + (0.15 if value >= 0 else -0.35), f"{value:+.3f}", ha="center", fontsize=8)
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_ylabel("LCB from midship [% L]")
    ax.set_title(f"C3: |Δ| = {abs(values[1] - values[0]):.3f} % L (limit 5 % L)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def figure_sweep(report, hydro, path: Path) -> None:
    sweep = report["criteria"]["C6"]["sweep"]
    cbs = [item["block_coefficient"] for item in sweep]
    volumes = [item["displaced_volume_m3"] for item in sweep]
    p = report["parameters"]
    design_line = [cb * p["length_m"] * p["beam_m"] * p["draft_m"] for cb in cbs]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(cbs, design_line, "k--", label="Cb·L·B·T (analytic)")
    ax.plot(cbs, volumes, "o-", color="#3b6fa0", label="NURBS quadrature")
    ax.set_xlabel("requested block coefficient Cb")
    ax.set_ylabel("displacement [m³]")
    ax.set_title("C6: the form follows the parameter monotonically")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def figure_mesh(hull, mesh, path: Path) -> None:
    stations, nz, layers = mesh.shape
    fig = plt.figure(figsize=(11.0, 4.6))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    for k in (0, layers // 2, layers):
        for j in range(0, nz + 1, 3):
            xs = [mesh.node(i, j, k) for i in range(0, stations + 1, 2)]
            ax.plot([q[0] for q in xs], [q[1] for q in xs], [q[2] for q in xs], lw=0.6)
    for i in range(0, stations + 1, 4):
        for k in (0, layers):
            xs = [mesh.node(i, j, k) for j in range(nz + 1)]
            ax.plot([q[0] for q in xs], [q[1] for q in xs], [q[2] for q in xs], lw=0.6)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_zlabel("z [m]")
    ax.set_title(f"Body-fitted hex mesh {mesh.shape} (starboard half)")
    ax.view_init(elev=22, azim=-60)

    ax = fig.add_subplot(1, 2, 2)
    for k in (0, layers // 2, layers):
        mid = stations // 2
        xs = [mesh.node(mid, j, k) for j in range(nz + 1)]
        ax.plot([q[1] for q in xs], [q[2] for q in xs], "-o", ms=2, lw=1.0, label=f"layer {k}/{layers}")
    ax.axhline(hull.parameters.draft, color="#c0392b", ls="--", lw=1.0, label="design waterline")
    ax.set_xlabel("y [m]")
    ax.set_ylabel("z [m]")
    ax.set_title("Midship column: the waterline is an exact grid level")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# result tables
# ---------------------------------------------------------------------------
def write_results(measured: dict, hydro, path: Path) -> None:
    report = measured["acceptance"]
    p = report["parameters"]
    lines: list[str] = []
    lines.append("# IN383-CAD01 result tables (generated)\n")
    lines.append(
        "Generated by `python tools/hull_nurbs_report.py`; every number below is measured in\n"
        "that run and stored in `artifacts/in383_cad01/acceptance.json`. Do not edit by hand.\n"
    )
    lines.append("## Principal dimensions and form\n")
    lines.append("| quantity | value | unit |\n|---|---|---|\n")
    for key in (
        "length_m",
        "beam_m",
        "draft_m",
        "deck_height_m",
        "block_coefficient",
        "design_displacement_m3",
        "length_to_beam",
        "beam_to_draft",
        "design_froude_number",
    ):
        lines.append(f"| {key} | {p[key]:.6g} | - |\n")
    lines.append(
        f"| calibrated parallel-body ratio | {report['hull']['form']['parallel_body_ratio']:.5f} | - |\n"
    )
    lines.append(f"| aft taper / entrance ratio | {report['hull']['form']['aft_taper_ratio']:.4f} / "
                 f"{report['hull']['form']['entrance_ratio']:.4f} | - |\n")
    lines.append(
        "| realised beam at waterline | "
        f"{report['hull']['realised_dimensions']['beam_at_waterline_m']:.6f} | m |\n"
    )
    lines.append(
        "| realised beam deficit | "
        f"{report['hull']['realised_dimensions']['beam_relative_deficit'] * 100:.3f} | % |\n"
    )

    lines.append("\n## Hydrostatics on the NURBS patch\n")
    lines.append("| quantity | value | unit |\n|---|---|---|\n")
    for key, unit in (
        ("displaced_volume_m3", "m³"),
        ("lcb_from_ap_m", "m"),
        ("lcb_ratio_from_midship", "fraction of L"),
        ("vertical_centre_of_buoyancy_m", "m"),
        ("waterplane_area_m2", "m²"),
        ("waterplane_coefficient", "-"),
        ("maximum_section_area_m2", "m²"),
        ("midship_section_coefficient", "-"),
        ("prismatic_coefficient", "-"),
        ("wetted_surface_area_m2", "m²"),
        ("side_area_m2", "m²"),
        ("bottom_area_m2", "m²"),
        ("transom_area_below_waterline_m2", "m²"),
    ):
        value = report["hydrostatics"][key]
        lines.append(f"| {key} | {value:.8g} | {unit} |\n")

    lines.append("\n## Acceptance criteria\n")
    lines.append("| id | criterion | measured | limit | pass | detail |\n|---|---|---|---|---|---|\n")
    for key in ("C1", "C2", "C3", "C4", "C5", "C6", "C7"):
        item = report["criteria"][key]
        lines.append(
            f"| {key} | {item['criterion']} | {item['value']:.6g} ({item['unit']}) | "
            f"{item['limit']:.6g} | {'PASS' if item['passed'] else 'FAIL'} | {item['detail']} |\n"
        )

    lines.append("\n## Mesh\n")
    quality = report["mesh"]["quality"]
    lines.append("| quantity | value |\n|---|---|\n")
    lines.append(f"| cells (half hull) | {quality['cells']} |\n")
    lines.append(f"| nodes (half hull) | {report['mesh']['nodes']} |\n")
    lines.append(f"| collapsed cells at the stem column | {quality['collapsed_cells']} |\n")
    lines.append(f"| pass fraction | {quality['pass_fraction'] * 100:.3f} % |\n")
    lines.append(f"| good fraction | {quality['good_fraction'] * 100:.3f} % |\n")
    lines.append(f"| mean scaled Jacobian | {quality['mean_scaled_jacobian']:.5f} |\n")
    lines.append(f"| max aspect ratio | {quality['max_aspect_ratio']:.1f} |\n")
    lines.append(f"| tetrahedra after subdivision | {int(measured['tetrahedra']['tetrahedra'])} |\n")
    lines.append(
        f"| tet/hex volume mismatch | {measured['tetrahedra']['relative_volume_mismatch']:.3e} |\n"
    )

    lines.append("\n## Mesh ladder\n")
    lines.append(
        "| stations x levels x layers | cells | V [m³] | S_wet [m²] | V err [%] | S err [%] | pass [%] |\n"
        "|---|---|---|---|---|---|---|\n"
    )
    for item in measured["ladder"]:
        lines.append(
            f"| {item['stations']} x {item['levels_below_waterline']} x {item['layers']} | "
            f"{item['cells']} | {item['displaced_volume_m3']:.6f} | "
            f"{item['wetted_surface_area_m2']:.6f} | "
            f"{item['volume_relative_error'] * 100:.4f} | "
            f"{item['wetted_area_relative_error'] * 100:.4f} | "
            f"{item['pass_fraction'] * 100:.2f} |\n"
        )

    lines.append("\n## Cb sweep (C6)\n")
    lines.append(
        "| Cb requested | Cb realised | V [m³] | V design [m³] | calibration residual | parallel body |\n"
        "|---|---|---|---|---|---|\n"
    )
    lbt = p["length_m"] * p["beam_m"] * p["draft_m"]
    for item in report["criteria"]["C6"]["sweep"]:
        lines.append(
            f"| {item['block_coefficient']:.4f} | {item['displaced_volume_m3'] / lbt:.6f} | "
            f"{item['displaced_volume_m3']:.6f} | "
            f"{item['design_displacement_m3']:.6f} | {item['calibration_residual']:.2e} | "
            f"{item['parallel_body_ratio']:.4f} |\n"
        )
    lines.append(
        f"\nMeasured dV/dCb = {report['criteria']['C6']['measured_dv_dcb_m3']:.6f} m³, "
        f"analytic L·B·T = {report['criteria']['C6']['analytic_dv_dcb_m3']:.6f} m³ "
        f"({report['criteria']['C6']['slope_relative_error'] * 100:.4f} % difference).\n"
    )
    lines.append(
        f"\nReachable Cb band of this form family: "
        f"[{report['criteria']['C6']['reachable_block_coefficient_band'][0]:.4f}, "
        f"{report['criteria']['C6']['reachable_block_coefficient_band'][1]:.4f}].\n"
    )

    lines.append("\n## Empirical baselines\n")
    lines.append("| baseline | value | source |\n|---|---|---|\n")
    lines.append(
        f"| Holtrop-Mennen wetted area | {measured['baselines']['holtrop_wetted_area_m2']:.6f} m² | "
        "Holtrop & Mennen (1984), no bulb (A_BT = 0) |\n"
    )
    lines.append(
        f"| Schneekluth-Bertram LCB | {measured['baselines']['schneekluth_lcb_ratio'] * 100:+.4f} % L | "
        "Schneekluth & Bertram (1998), Fn = "
        f"{p['design_froude_number']:.3f} |\n"
    )
    path.write_text("".join(lines), encoding="utf-8")


def main() -> int:
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    DOC_DIR.mkdir(parents=True, exist_ok=True)

    measured = measure()
    report = measured["acceptance"]
    hull = build_nurbs_hull()
    hydro = nurbs_hydrostatics(hull)
    mesh = build_body_mesh(hull)

    figure_lines(hull, measured["sections"], FIGURE_DIR / "body-plan-waterlines.png")
    figure_quality(
        measured["scaled_jacobians"],
        measured["aspect_ratios"],
        mesh.quality,
        FIGURE_DIR / "mesh-quality.png",
    )
    figure_convergence(measured["ladder"], hydro, FIGURE_DIR / "convergence.png")
    figure_comparison(report, hydro, measured["baselines"], FIGURE_DIR / "routes-comparison.png")
    figure_sweep(report, hydro, FIGURE_DIR / "cb-sweep.png")
    figure_mesh(hull, mesh, FIGURE_DIR / "body-mesh.png")

    artifact = {
        "schema": report["schema"],
        "mission": "IN383-CAD01",
        "repository": "TensorCAD",
        "acceptance": report,
        "mesh_ladder": measured["ladder"],
        "section_areas": measured["sections"],
        "tetrahedra": measured["tetrahedra"],
        "tet_count": measured["tet_count"],
        "baselines": measured["baselines"],
        "figures": sorted(path.name for path in FIGURE_DIR.glob("*.png")),
    }
    out = ARTIFACT_DIR / "acceptance.json"
    out.write_text(json.dumps(artifact, indent=2, sort_keys=False), encoding="utf-8")
    write_results(measured, hydro, DOC_DIR / "results.md")

    print(f"artifact: {out.relative_to(ROOT)}")
    print(f"figures : {len(artifact['figures'])} in {FIGURE_DIR.relative_to(ROOT)}")
    print(f"tables  : {(DOC_DIR / 'results.md').relative_to(ROOT)}")
    print("criteria:")
    for key in ("C1", "C2", "C3", "C4", "C5", "C6", "C7"):
        item = report["criteria"][key]
        print(
            f"  {key} {item['criterion']:<34} {item['value']:.6g} {item['unit']:<18} "
            f"<= {item['limit']:.6g}  {'PASS' if item['passed'] else 'FAIL'}"
        )
    print(f"all_passed = {report['all_passed']}")
    return 0 if report["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
