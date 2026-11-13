"""Regenerate the GAP-CAD-003 evidence bundle: ``artifacts/mesh_boolean_v1/``.

Runs the measured case matrix through ``tensorcad.mesh_boolean`` (no numbers are typed by
hand), writes ``results.json`` and the two-panel figure, and prints the summary that
``docs/mesh_boolean.md`` interpolates.

Run: ``python tools/mesh_boolean_evidence.py``
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

from test_mesh_boolean import box, cylinder, euler_characteristic  # noqa: E402

from tensorcad.mesh import is_edge_manifold, is_watertight, total_signed_volume  # noqa: E402
from tensorcad.mesh_boolean import MeshBooleanError, boolean  # noqa: E402

OUT = ROOT / "artifacts" / "mesh_boolean_v1"
OUT.mkdir(parents=True, exist_ok=True)


def main() -> int:
    prism = cylinder(1.0, 2.0)
    cube = box((-1.5, -1.5, 0.0), (1.5, 1.5, 2.0))
    small = box((1.0, 1.0, 1.0), (3.0, 3.0, 3.0))
    rows: list[dict] = []

    def record(name, op, first, second, reference, chi_expected):
        started = time.time()
        try:
            mesh = boolean(first, second, op)
        except MeshBooleanError as exc:
            rows.append(dict(case=name, operation=op, status="refused",
                             reason=str(exc).split(";")[0],
                             seconds=round(time.time() - started, 3)))
            return None
        volume = total_signed_volume(mesh)
        rows.append(dict(
            case=name, operation=op, status="valid", volume=volume, reference=reference,
            rel_error=abs(volume - reference) / abs(reference),
            chi=euler_characteristic(mesh), chi_expected=chi_expected,
            faces=len(mesh.faces), watertight=is_watertight(mesh),
            edge_manifold=is_edge_manifold(mesh), seconds=round(time.time() - started, 3)))
        return mesh

    unit = box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    record("offset cubes 2x2x2 vs 2x2x2 (+1)", "union", box((0, 0, 0), (2, 2, 2)), small, 15.0, 2)
    record("offset cubes", "intersection", box((0, 0, 0), (2, 2, 2)), small, 1.0, 2)
    record("offset cubes", "difference", box((0, 0, 0), (2, 2, 2)), small, 7.0, 2)
    record("coplanar shared face", "union", unit, box((1, 0, 0), (2, 1, 1)), 2.0, 2)
    record("nested cavity 3^3 - 1^3", "difference", box((0, 0, 0), (3, 3, 3)),
           box((1, 1, 1), (2, 2, 2)), 26.0, 4)
    record("identical solids", "union", unit, unit, 1.0, 2)
    record("identical solids", "intersection", unit, unit, 1.0, 2)
    tessellated = 18.0 - total_signed_volume(prism)
    holed = record("cube 3x3x2 - flush through bolt r=1 (tessellated ref)", "difference",
                   cube, prism, tessellated, 0)
    analytic = 18.0 - 2.0 * math.pi
    rows.append(dict(case="cube 3x3x2 - flush through bolt (analytic 18-2pi)",
                     operation="difference", status="valid",
                     volume=total_signed_volume(holed), reference=analytic,
                     rel_error=abs(total_signed_volume(holed) - analytic) / analytic,
                     chi=euler_characteristic(holed), chi_expected=0,
                     faces=len(holed.faces), watertight=is_watertight(holed),
                     edge_manifold=is_edge_manifold(holed), seconds=0.0))
    corner = box((0.5, 0.5, 0.5), (1.5, 1.5, 1.5))
    record("cylinder - box (unsupported)", "difference", prism,
           box((-1.5, -1.5, 0.5), (1.5, 1.5, 2.0)), None, 0)
    record("cylinder + corner box (unsupported)", "union", prism, corner, None, 0)
    record("cylinder & corner box (unsupported)", "intersection", prism, corner, None, 0)

    valid = [row for row in rows if row["status"] == "valid"]
    summary = dict(
        n_cases=len(rows),
        n_valid=len(valid),
        n_refused=len(rows) - len(valid),
        max_rel_error_on_polyhedral_reference=max(
            row["rel_error"] for row in valid if "analytic 18-2pi" not in row["case"]),
        analytic_flush_cut_rel_error=next(
            row["rel_error"] for row in valid if "analytic 18-2pi" in row["case"]),
        tessellation_error_prism=abs(total_signed_volume(prism) - 2.0 * math.pi) / (2.0 * math.pi),
    )
    payload = dict(generated=time.strftime("%Y-%m-%dT%H:%M:%S"), cases=rows, summary=summary)
    (OUT / "results.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False))

    figure = plt.figure(figsize=(11, 4.6))
    axis = figure.add_subplot(1, 2, 1, projection="3d")
    axis.add_collection3d(Poly3DCollection(
        [[holed.vertices[i] for i in face] for face in holed.faces],
        facecolor="#7fb3d5", edgecolor="#1b4f72", linewidths=0.15, alpha=0.9))
    axis.set_xlim(-1.5, 1.5)
    axis.set_ylim(-1.5, 1.5)
    axis.set_zlim(0, 2)
    axis.set_box_aspect((1, 1, 0.67))
    axis.set_title(f"cube minus flush through bolt\n{len(holed.faces)} tris, "
                   f"chi = {euler_characteristic(holed)} (genus 1)")
    axis.set_xlabel("x [m]")
    axis.set_ylabel("y [m]")
    axis.set_zlabel("z [m]")

    chart = figure.add_subplot(1, 2, 2)
    labels = ["offset union", "offset inter", "offset diff", "shared face", "cavity",
              "identical U", "identical I", "flush cut"]
    values = [max(row["rel_error"], 1e-16) for row in valid[:7]]
    values.append(summary["analytic_flush_cut_rel_error"])
    chart.bar(range(len(values)), values, color=["#1e8449"] * 7 + ["#b9770e"])
    chart.axhline(0.01, color="#c0392b", linestyle="--", linewidth=1.2,
                  label="GAP-CAD-003 gate: 1 %")
    chart.set_yscale("log")
    chart.set_xticks(range(len(labels)))
    chart.set_xticklabels(labels, rotation=35, ha="right", fontsize=7)
    chart.set_ylabel("relative volume error")
    chart.set_title("volume error vs analytic reference")
    chart.legend(fontsize=7)
    figure.tight_layout()
    figure.savefig(OUT / "mesh_boolean_evidence.png", dpi=140)

    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
