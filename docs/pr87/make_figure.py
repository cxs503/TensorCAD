"""Figure for the PR #87 CI-parity report: every bar is read from docs/trim/results.json.

Panel (a): relative error of each measured quantity against its own closed form / analytic
reference (box particulars, the hull's two centimetre routes, the trim angles).
Panel (b): the trim angle's convergence with station count, relative to the finest mesh.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
DATA = json.loads((ROOT / "docs/trim/results.json").read_text())

box = DATA["box_particulars"]
hull = DATA["hull_particulars"]
wig = DATA["wigley"]


def rel(measured: float, exact: float) -> float:
    return abs(measured - exact) / abs(exact)


box_items = [(name, abs(row["rel"])) for name, row in sorted(box.items())]
hull_items = [
    ("hull TPC vs dV/dz", rel(hull["tpc"], hull["tpc_from_dV_dz"])),
    ("hull MCT vs I_L", rel(hull["mct"], hull["mct_analytic"])),
    ("hull trim vs small-angle", abs(DATA["hull_trim"][0]["rel"])),
    ("hull volume vs analytic", abs(hull["volume_rel"])),
]
wig_items = [
    ("Wigley TPC vs dV/dz", rel(wig["tpc"], wig["tpc_from_dV_dz"])),
    ("Wigley trim vs small-angle", rel(wig["trim_deg"], wig["small_angle_deg"])),
]
box_trim = abs(DATA["box_trim"][0]["rel"])
extra = [("box trim (8 deg, shift 8 m)", box_trim)]

labels = [k for k, _ in box_items] + [k for k, _ in extra]
values = [v if v > 0 else 1e-17 for _, v in box_items] + [v if v > 0 else 1e-17 for _, v in extra]
for name, value in hull_items + wig_items:
    labels.append(name)
    values.append(value if value > 0 else 1e-17)

fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))
axes[0].barh(range(len(labels)), values, color="#2b6cb0")
axes[0].set_yticks(range(len(labels)))
axes[0].set_yticklabels(labels, fontsize=8)
axes[0].set_xscale("log")
axes[0].axvline(0.03, color="#c53030", linestyle="--", linewidth=1.2, label="3 % criterion")
axes[0].axvline(5e-16, color="#718096", linestyle=":", linewidth=1.0, label="round-off (5e-16)")
axes[0].set_xlabel("relative error vs closed form / analytic route")
axes[0].set_title("(a) Trim module: measured vs reference (log scale)")
axes[0].legend(fontsize=8, loc="lower right")
axes[0].grid(axis="x", alpha=0.3)

conv = sorted(DATA["hull_trim_convergence"], key=lambda row: row["stations"])
stations = [row["stations"] for row in conv]
err = [max(row["rel_to_finest"], 1e-17) for row in conv]
axes[1].loglog(stations, err, "o-", color="#2f855a", label="trim angle vs finest mesh")
axes[1].axhline(0.03, color="#c53030", linestyle="--", linewidth=1.2, label="3 % criterion")
axes[1].set_xlabel("station count (Simpson quadrature in x)")
axes[1].set_ylabel("relative change of trim angle vs finest mesh")
axes[1].set_title("(b) Discretisation: trim angle convergence, parametric hull")
axes[1].legend(fontsize=8)
axes[1].grid(True, which="both", alpha=0.3)

fig.tight_layout()
out = Path(__file__).resolve().parent / "trim_ci_evidence.png"
fig.savefig(out, dpi=140)
print("wrote", out)
print("max box rel:", max(v for _, v in box_items))
print("3% verdict:", all(v < 0.03 for v in values))
