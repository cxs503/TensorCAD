#!/usr/bin/env python3
"""Figure for the m17634 board audit of m17595 (XVAL six-format matrix).

Panel A: independent recount (this audit) of the volume of every TensorCAD-written
         artefact, as a relative deviation from the geometry's reference volume
         (cube = 8.0 exact; icosphere/torus = median over the six formats).
Panel B: how many of the readers that produced a number are third-party (independent)
         versus subject-authored (self-consistency only) -- the axis the headline
         "OSS cross-read" wording glosses over.

Inputs: docs/board_audit_m17634/audit_evidence.json (produced by audit_independent.py)
Output: docs/board_audit_m17634/figures/audit_volume_deviation.png
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ev = json.loads((HERE / "audit_evidence.json").read_text())["checks"]

FMT_ORDER = ["stl", "obj", "vtk", "msh", "cgns", "su2"]
FMT_LABEL = {"stl": "STL", "obj": "OBJ", "vtk": "VTK", "msh": "Gmsh",
             "cgns": "CGNS", "su2": "SU2"}
OWN_READERS = {"cgns-sids-audit", "su2-audit"}  # written by this audit, *not* by the subject;
# subject-side self-readers are reported separately below from results.json

cells = ev["independent_reread"]
geoms = sorted({k.split("/")[0] for k in cells})
ref = {g: ev["volume_agreement"][g]["reference"] for g in geoms}

mat = np.full((len(geoms), len(FMT_ORDER)), np.nan)
for i, g in enumerate(geoms):
    for j, f in enumerate(FMT_ORDER):
        vals = [v["volume"] for k, v in cells.items()
                if k.startswith(f"{g}/{f}/") and v.get("ok")]
        if vals:
            mat[i, j] = max(abs(v - ref[g]) for v in vals) / ref[g]

fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.2, 4.6), gridspec_kw={"width_ratios": [2.05, 1]})

floor = 1e-18
shown = np.where(np.isnan(mat), floor, np.maximum(mat, floor))
im = ax.imshow(shown, cmap="viridis", aspect="auto",
               norm=matplotlib.colors.LogNorm(vmin=1e-12, vmax=1e-5))
ax.set_xticks(range(len(FMT_ORDER)), [FMT_LABEL[f] for f in FMT_ORDER])
ax.set_yticks(range(len(geoms)), geoms)
ax.set_xlabel("format written by TensorCAD (m17483 writers)")
ax.set_ylabel("geometry")
ax.set_title("A. independent recount of volume: rel. deviation from reference\n"
             "(upper bound over all readers that returned a volume)", fontsize=10)
for i in range(len(geoms)):
    for j in range(len(FMT_ORDER)):
        v = mat[i, j]
        txt = "not read" if np.isnan(v) else ("0" if v < 1e-15 else f"{v:.1e}")
        ax.text(j, i, txt, ha="center", va="center", fontsize=8,
                color="white" if np.isnan(v) or v < 1e-7 else "black")
ax.text(0.0, -0.30, "owner precision line 15 % = 1.5e-1  |  mission's own gate 1e-3  |  colour scale 1e-12..1e-5",
        transform=ax.transAxes, fontsize=8, color="#444")
fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02, label="relative volume error")

# panel B: oracle independence per format
indep = {"STL": ["meshio", "trimesh", "numpy-stl"], "OBJ": ["meshio", "trimesh"],
         "VTK": ["meshio (UNSTRUCTURED_GRID only)"], "Gmsh": ["meshio"],
         "CGNS": ["none (meshio CGNS reader is a Zone1/tetra stub)"],
         "SU2": ["none (meshio has no SU2 reader)"]}
own = {"STL": 0, "OBJ": 0, "VTK": 0, "Gmsh": 0, "CGNS": 1, "SU2": 1}
labels = list(indep)
n3p = [len(indep[l]) if "none" not in indep[l][0] else 0 for l in labels]
bars1 = ax2.bar(range(len(labels)), n3p, color="#2a7ab9", label="third-party readers")
bars2 = ax2.bar(range(len(labels)), [own[l] for l in labels], bottom=n3p,
                color="#d1495b", label="subject-authored reader")
ax2.set_xticks(range(len(labels)), labels, fontsize=9)
ax2.set_ylabel("geometry oracles per format")
ax2.set_title("B. who verified the geometry of each format?\n"
              "(CGNS/SU2 geometry = self-consistency unless an auditor re-derives it)", fontsize=10)
ax2.set_ylim(0, 3.9)
for b, l in zip(bars1, labels):
    if n3p[labels.index(l)] == 0:
        ax2.text(b.get_x() + b.get_width() / 2, 0.12, "0 third-party", ha="center",
                 fontsize=8, rotation=90, color="white")
for b, l in zip(bars2, labels):
    if own[l]:
        ax2.text(b.get_x() + b.get_width() / 2, n3p[labels.index(l)] + 0.12, "+ own",
                 ha="center", fontsize=8, color="#d1495b")
ax2.legend(fontsize=8, loc="upper right")

fig.suptitle("m17595 XVAL six-format matrix — independent recount (audit m17634)", fontsize=12)
fig.tight_layout(rect=(0, 0.02, 1, 0.95))
out = HERE / "figures" / "audit_volume_deviation.png"
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=140)
print("wrote", out, out.stat().st_size, "bytes")
print("matrix:", np.array2string(mat, precision=2))