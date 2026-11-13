#!/usr/bin/env python3
"""Figures for the m17401 board audit (all data read from verify.json)."""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
AUD = ROOT / "artifacts/board_audit/m17401"
v = json.loads((AUD / "verify.json").read_text())
OUT = AUD / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# ---- figure 1: independent volume cross-check ----------------------------------
cases = ["plate_with_hole", "stepped_shaft", "bushing", "l_bracket"]
subj = [v["C1_rerun"]["per_file"][c + ".step"]["error_percent"] for c in cases]
indep = [v["C2_volume_independent"][c]["error_percent"] for c in cases]

fig, ax = plt.subplots(figsize=(7.2, 4.0))
x = range(len(cases))
w = 0.36
b1 = ax.bar([i - w / 2 for i in x], subj, w, label="subject run (own driver)", color="#4c72b0")
b2 = ax.bar([i + w / 2 for i in x], indep, w, label="auditor independent OBJ re-integration",
            color="#dd8452", hatch="//", edgecolor="k", linewidth=0.4)
ax.axhspan(-1, 1, color="#c7e9c0", alpha=0.35, label="acceptance band |err| < 1 %")
ax.axhline(0, color="k", linewidth=0.8)
for bars in (b1, b2):
    for r in bars:
        h = r.get_height()
        ax.annotate(f"{h:+.3f}%", (r.get_x() + r.get_width() / 2, h),
                    textcoords="offset points", xytext=(0, 6 if h >= 0 else -12),
                    ha="center", fontsize=7.5)
ax.set_xticks(list(x))
ax.set_xticklabels([c.replace("_", "\n") for c in cases], fontsize=8.5)
ax.set_ylabel("volume error vs closed-form nominal [%]")
ax.set_title("m17349 STEP B-rep reader: volume error, subject vs auditor (independent re-integration)")
ax.legend(fontsize=7.5, loc="lower left")
ax.margins(y=0.25)
fig.tight_layout()
fig.savefig(OUT / "audit_volume_crosscheck.png", dpi=150)

# ---- figure 2: convergence order ----------------------------------------------
rec = v["C1_rerun"]["convergence"]
steps = [r["angular_step"] for r in rec]
fig, ax = plt.subplots(figsize=(7.2, 4.2))
colors = {"plate_with_hole.step": "#4c72b0", "stepped_shaft.step": "#dd8452",
          "bushing.step": "#55a868"}
for key, col in colors.items():
    err = [abs(r[key]) for r in rec]
    ax.loglog(steps, err, "o-", label=key.replace(".step", ""), color=col, markersize=4)
ref = [abs(rec[0]["plate_with_hole.step"]) * (s / steps[0]) ** 2 for s in steps]
ax.loglog(steps, ref, "k--", linewidth=1.0, label="slope 2 (second order) reference")
ratios = [abs(rec[i]["plate_with_hole.step"] / rec[i + 1]["plate_with_hole.step"])
          for i in range(len(rec) - 1)]
ax.set_xlabel("angular step [rad]")
ax.set_ylabel("|volume error| [%]")
ax.set_title("Error vs angular step (halving ratio "
             + "/".join(f"{r:.2f}" for r in ratios) + " -> second order)")
ax.invert_xaxis()
ax.grid(True, which="both", alpha=0.3)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "audit_convergence_order.png", dpi=150)
print("wrote", sorted(p.name for p in OUT.glob("*.png")))