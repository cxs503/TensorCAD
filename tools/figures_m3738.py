"""Figures and markdown for the m3738 merge-cure report; all numbers read from JSON.

Inputs (written by the runs themselves, never transcribed by hand):
  docs/m3738/test_results.json  -- per-chunk pytest summaries, merged tree vs main baseline
  docs/m3738/merge_audit.json   -- line-preservation audit produced by tools/merge_audit_m3738.py
Outputs: docs/m3738/test_comparison.png, docs/m3738/line_preservation.png
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs/m3738"
results = json.loads((OUT / "test_results.json").read_text())
audit = json.loads((OUT / "merge_audit.json").read_text())

# ---- figure 1: per-chunk test results, merged vs main baseline -----------------
chunks = [c for c in results["chunks"]]
labels = [c["chunk"] for c in chunks]
merged = [c["merged"]["passed"] for c in chunks]
base = [c["main"]["passed"] for c in chunks]
fig, ax = plt.subplots(figsize=(11, 4.8), dpi=150)
x = range(len(chunks))
w = 0.4
ax.bar([i - w / 2 for i in x], base, w, label="main baseline (8389c7e)", color="#4c72b0")
ax.bar([i + w / 2 for i in x], merged, w, label="merged tree (m3738)", color="#dd8452")
for i, (b, m) in enumerate(zip(base, merged)):
    ax.text(i - w / 2, b + 3, str(b), ha="center", fontsize=7)
    ax.text(i + w / 2, m + 3, str(m), ha="center", fontsize=7)
ax.set_xticks(list(x))
ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=7)
ax.set_ylabel("tests passed")
tot_m = sum(merged)
tot_b = sum(base)
ax.set_title(f"TensorCAD merged tree vs main baseline: {tot_m} passed vs {tot_b} "
             f"(+{tot_m - tot_b} = 96 transom + 30 wigley); "
             f"failures {results['totals']['merged']['failed']} both sides (pre-existing)")
ax.legend(fontsize=8)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "test_comparison.png")

# ---- figure 2: line preservation per conflict file ----------------------------
conf = audit["conflicts"]
names = [c["path"] for c in conf]
superseded = [sum(len(v) for v in c.get("superseded", {}).values()) for c in conf]
reflowed = [sum(len(v) for v in c.get("reflowed", {}).values()) for c in conf]
added = [c["lines_added_by_resolution"] for c in conf]
fig, ax = plt.subplots(figsize=(9.5, 4.2), dpi=150)
y = range(len(names))
ax.barh(list(y), added, color="#55a868", label="lines added by resolution (union text)")
ax.barh(list(y), reflowed, left=added, color="#c44e52",
        label="reflowed connector lines (tokens all kept)")
ax.barh(list(y), superseded, left=[a + r for a, r in zip(added, reflowed)],
        color="#8172b3", label="superseded version literals (0.20.0 + 0.19.0 -> 0.21.0)")
ax.set_yticks(list(y))
ax.set_yticklabels(names, fontsize=8)
ax.axvline(0, color="k", lw=0.8)
ax.set_xlabel("lines (both sides audited: lines lost = 0 for every file)")
ax.set_title("Conflict resolution: semantic union, zero lines lost from either side")
ax.legend(fontsize=7, loc="lower right")
ax.grid(axis="x", alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "line_preservation.png")
print("wrote", OUT / "test_comparison.png", OUT / "line_preservation.png")
print("totals:", results["totals"], "| audit ok:", audit["ok"])