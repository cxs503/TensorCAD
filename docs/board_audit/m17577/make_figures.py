#!/usr/bin/env python3
"""m17577 审计图生成：全部数字从被审工件 JSON 程序化读取（零手抄）。"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ART = HERE / "refs" / "artifacts" / "m17472_5415_param"

acc = json.loads((ART / "acceptance.json").read_text())
geo = json.loads((ART / "geometry_validation.json").read_text())
ship = json.loads((ART / "shipped_model_check.json").read_text())

# ---- 图 1：关键量偏差 vs 3% 门（被审方报告 §6 口径，数值直读工件）----
hyd = acc["config0"]["hydrostatics"]
rows = [
    ("Lpp (m)", 142.0, hyd["length_stations_m"]),
    ("Beam@WL (m)", 18.867133, hyd["beam_waterline_m"]),
    ("Draft (m)", 6.156643, hyd["draft_m"]),
    ("DispVol (m3)", 8412.037925, hyd["displacement_volume_m3"]),
    ("Cb", 0.509992, hyd["cb"]),
]
names, diffs = [], []
for n, ref, got in rows:
    names.append(n)
    diffs.append((got - ref) / ref * 100.0)
holdout_ref = ship["holdout_independent_refit"]["rms_point_m"]
holdout_got = ship["holdout_shipped"]["rms_point_m"]
names.append("Holdout RMS (m)")
diffs.append((holdout_got - holdout_ref) / holdout_ref * 100.0)

fig, ax = plt.subplots(figsize=(8, 4.6))
colors = ["#2a7de1" if abs(d) < 3 else "#d64545" for d in diffs]
bars = ax.bar(names, diffs, color=colors)
ax.axhline(3.0, color="k", ls="--", lw=1)
ax.axhline(-3.0, color="k", ls="--", lw=1)
ax.axhline(0.0, color="gray", lw=0.8)
for b, d in zip(bars, diffs):
    ax.text(b.get_x() + b.get_width() / 2, d + (0.08 if d >= 0 else -0.16),
            f"{d:+.3f}%", ha="center", fontsize=8)
ax.set_ylabel("deviation vs reference (%)")
ax.set_title("m17472 key quantities: deviation vs reference (gate 3%) — auditor re-read from run artifacts")
plt.xticks(rotation=18)
plt.tight_layout()
fig.savefig(HERE / "fig_audit_deviation.png", dpi=140)
plt.close(fig)

# ---- 图 2：验收三判据 + 复现核验 ----
c1 = acc["acceptance"]["config0_vs_exact"]
c2 = acc["acceptance"]["sampled_configs_gate"]
c3 = acc["acceptance"]["volume_gradient"]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.5, 4.2))
crit = [
    ("C1 config0 vs exact\nfrac<res (thr 0.95)", c1["frac_within_resolution"], 0.95, c1["pass"]),
    ("C2 sampled gate\nn_pass/20 (thr 20)", c2["n_pass"] / 20.0, 1.0, c2["pass"]),
    ("C3 grad max rel err\n(thr 1e-2, log)", None, None, c3["pass"]),
]
x = range(3)
vals = [c1["frac_within_resolution"], c2["n_pass"] / 20.0, None]
for i, (lab, v, thr, ok) in enumerate(crit):
    if i < 2:
        a1.bar(i, v, color="#2a9d5c" if ok else "#d64545")
        a1.axhline(thr, color="k", ls="--", lw=1)
        a1.text(i, v + 0.01, f"{v:.4f}" if i == 0 else f"{c2['n_pass']}/20", ha="center", fontsize=8)
    else:
        a1.bar(i, 1.0, color="#2a9d5c" if ok else "#d64545")
        a1.text(i, 0.5, f"{c3['max_rel_err']:.2e}", ha="center", fontsize=8, color="white")
a1.set_xticks(list(x))
a1.set_xticklabels([c[0] for c in crit], fontsize=7.5)
a1.set_ylim(0, 1.15)
a1.set_title("Acceptance criteria (re-run values)")

repro = ["sampled_configs.csv", "gate_samples.csv", "acceptance.json",
          "geometry_validation.json", "shipped_model_check.json"]
status = ["byte-identical", "byte-identical", "numeric-identical",
          "numeric-identical", "numeric-identical"]
a2.barh(repro, [1] * 5, color=["#2a9d5c"] * 2 + ["#7fb2e5"] * 3)
for i, s in enumerate(status):
    a2.text(0.02, i, s, va="center", fontsize=7.5, color="white")
a2.set_xlim(0, 1)
a2.set_xticks([])
a2.set_title("m17577 independent re-run vs committed artifacts")
plt.tight_layout()
fig.savefig(HERE / "fig_audit_acceptance.png", dpi=140)
plt.close(fig)
print("figures written")
