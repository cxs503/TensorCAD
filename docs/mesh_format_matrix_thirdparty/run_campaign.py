"""Run the 6-format x 3-geometry third-party matrix and emit table + figure.

Writes (next to this script): matrix.json, matrix.md, matrix_heatmap.png.
Regenerate with:  python docs/mesh_format_matrix_thirdparty/run_campaign.py
"""
from __future__ import annotations

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tensorcad.mesh_interop_independent import ANALYTIC_VOLUME_3, matrix_markdown, run_format_matrix

HERE = pathlib.Path(__file__).resolve().parent
GEOMETRIES = ("cube", "capsule", "wigley_hull")
FORMATS = ("stl", "obj", "vtk", "gmsh", "su2", "cgns")
LABEL = {"cube": "cube (V=1)", "capsule": "capsule", "wigley_hull": "Wigley hull"}


def main() -> int:
    rows = run_format_matrix(workdir=HERE / "files")
    (HERE / "matrix.json").write_text(json.dumps(rows, indent=2, default=str) + "\n")
    (HERE / "matrix.md").write_text(matrix_markdown(rows))

    fig, (ax_v, ax_e) = plt.subplots(1, 2, figsize=(13.5, 4.6))
    width = 0.13
    colors = plt.cm.tab10(np.linspace(0, 1, len(FORMATS)))
    for gi, geom in enumerate(GEOMETRIES):
        vals = [next(r for r in rows if r["geometry"] == geom and r["format"] == f)["volume_independent"] for f in FORMATS]
        for fi, (fmt, v) in enumerate(zip(FORMATS, vals)):
            x = gi + (fi - (len(FORMATS) - 1) / 2) * width
            ax_v.bar(x, v, width=width * 0.95, color=colors[fi], label=fmt if gi == 0 else None)
            ax_e.bar(x, max(abs(v - vals[1]) / abs(vals[1]), 1e-18), width=width * 0.95,
                     color=colors[fi], label=None)
        ax_v.hlines(ANALYTIC_VOLUME_3[geom], gi - 0.42, gi + 0.42, color="k", linestyle="--", linewidth=1.2,
                    label="analytic" if gi == 0 else None)
        ax_v.text(gi, ANALYTIC_VOLUME_3[geom] * 1.01, f"{ANALYTIC_VOLUME_3[geom]:.4f}", ha="center", va="bottom", fontsize=8)
    ax_v.set_xticks(range(len(GEOMETRIES)))
    ax_v.set_xticklabels([LABEL[g] for g in GEOMETRIES])
    ax_v.set_ylabel("enclosed volume V [m$^3$]")
    ax_v.set_title("Volume recovered by an independent parser (6 formats per geometry)")
    ax_v.legend(fontsize=8, ncol=2)

    ax_e.set_yscale("log")
    ax_e.axhline(1e-12, color="r", linestyle=":", linewidth=1.2, label="acceptance 1e-12")
    ax_e.set_xticks(range(len(GEOMETRIES)))
    ax_e.set_xticklabels([LABEL[g] for g in GEOMETRIES])
    ax_e.set_ylabel(r"$|V_{\rm read} - V_{\rm OBJ}| / V_{\rm OBJ}$  (format equivalence)")
    ax_e.set_title("Cross-format volume identity (floor = float64 round-off)")
    ax_e.legend(fontsize=8)
    fig.suptitle("TensorCAD mesh exchange: 6 formats x 3 geometries, verified without meshio", fontsize=11)
    fig.tight_layout()
    fig.savefig(HERE / "matrix_heatmap.png", dpi=140)
    print(matrix_markdown(rows))
    print("rows:", len(rows), "all ok:", all(r["ok"] for r in rows))
    return 0 if all(r["ok"] for r in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
