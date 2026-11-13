#!/usr/bin/env python3
"""Figure for the m17405 audit: contract-violation matrix + volume parity plot.

Left  : 8 shipped STL assets x 4 conventions — which convention fired and at
        what severity (blank = clean check).
Right : measured signed volume (recomputed from the shipped STL) vs the
        closed-form analytic volume, with the +/-3% band and per-point relative
        error labels.

All numbers come from artifacts/geom_contract_audit/probe_results.json and the
subject's compliance_report.json — nothing hand-typed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SUBJECT = Path(sys.argv[1] if len(sys.argv) > 1 else "/root/tp-agent/work/m17357")
PROBE = json.loads((ROOT / "artifacts/geom_contract_audit/probe_results.json").read_text())
COMP = json.loads((SUBJECT / "artifacts/geom_contract_assets/compliance_report.json").read_text())

CONVENTIONS = ["coordinate_system", "units", "scale", "orientation"]
SHORT = {"coordinate_system": "coordinate\nsystem", "units": "units",
         "scale": "scale", "orientation": "orientation"}

# ---- matrix: 0 = clean, 1 = warning, 2 = error -----------------------------
assets = [r["source"].rsplit("/", 1)[-1].replace(".stl", "") for r in COMP]
M = np.zeros((len(COMP), len(CONVENTIONS)))
for i, r in enumerate(COMP):
    for v in r["violations"]:
        M[i, CONVENTIONS.index(v["convention"])] = 2 if v["severity"] == "error" else 1

fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2),
                         gridspec_kw={"width_ratios": [1.05, 1.0]})
ax = axes[0]
cmap = ListedColormap(["#eef2f6", "#f6c445", "#cf3b2f"])
ax.imshow(M, cmap=cmap, vmin=0, vmax=2, aspect="auto")
ax.set_xticks(range(len(CONVENTIONS)))
ax.set_xticklabels([SHORT[c] for c in CONVENTIONS], fontsize=9)
ax.set_yticks(range(len(assets)))
ax.set_yticklabels([f"{a}   {'PASS' if COMP[i]['passed'] else 'FAIL'}"
                    for i, a in enumerate(assets)], fontsize=9)
for i, r in enumerate(COMP):
    for j in range(len(CONVENTIONS)):
        if M[i, j]:
            ax.text(j, i, "ERR" if M[i, j] == 2 else "warn", ha="center", va="center",
                    fontsize=8, color="white" if M[i, j] == 2 else "#4a3000", fontweight="bold")
ax.set_title("(a) 8 shipped STL assets × 4 contract conventions\n"
             f"{sum(1 for r in COMP if r['passed'])} PASS / "
             f"{sum(1 for r in COMP if not r['passed'])} FAIL "
             "(red-team assets fail on purpose)", fontsize=10)
ax.legend(handles=[Patch(facecolor="#cf3b2f", label="error"),
                   Patch(facecolor="#f6c445", label="warning"),
                   Patch(facecolor="#eef2f6", edgecolor="#999", label="clean")],
          loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3, fontsize=8, frameon=False)
ax.set_xticks(np.arange(-.5, len(CONVENTIONS), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(assets), 1), minor=True)
ax.grid(which="minor", color="white", linewidth=2)
ax.tick_params(which="minor", length=0)

# ---- parity: got vs analytic ref -------------------------------------------
ax = axes[1]
ref = np.array([a["ref"] for a in PROBE["anchors"]])
got = np.array([a["got"] for a in PROBE["anchors"]])
labels = ["box barge\nL·B·D", "unit sphere\n4/3πr³", "Wigley\n(4/9)L·B·T"]
lo, hi = 1.0, 1400.0
ax.fill_between([lo, hi], [lo * 0.97, hi * 0.97], [lo * 1.03, hi * 1.03],
                color="#cfe3ff", alpha=0.55, label="±3% band")
ax.plot([lo, hi], [lo, hi], "k--", lw=1, label="y = x")
ax.scatter(ref, got, s=70, color="#1f6feb", zorder=3, edgecolor="white")
for x, y, lab, a in zip(ref, got, labels, PROBE["anchors"]):
    ax.annotate(f"{lab}\n{a['rel_err_pct']:+.3f}%", (x, y), textcoords="offset points",
                xytext=(10, -22), fontsize=8)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(lo, hi); ax.set_ylim(lo * 0.9, hi * 0.9)
ax.set_xlabel("closed-form analytic volume [m³]  (log)")
ax.set_ylabel("signed volume recomputed from shipped STL [m³]  (log)")
ax.set_title("(b) validator measures real geometry:\n"
             "divergence-theorem volume vs analytic reference", fontsize=10)
ax.legend(fontsize=8, loc="upper left")
ax.grid(which="both", alpha=0.25, lw=0.5)

fig.suptitle("m17405 audit · GAP-CAD-D3-GEOM-CONTRACT (m17357, TensorCAD) — "
             "contract enforcement is real, but no registered benchmark case ⇒ L2", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.95])
out = ROOT / "artifacts" / "geom_contract_audit" / "contract_audit.png"
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=150)
print("saved", out, out.stat().st_size, "bytes")