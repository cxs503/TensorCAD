"""Build the m22814 report figure (old vs delivered flip counts) — one PNG, two panels."""
import json
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

NEW = json.load(open("docs/lines_fairing/results.json"))
OLD = json.loads(subprocess.check_output(
    ["git", "show", "061c8a7:docs/lines_fairing/results.json"]))

rows = []
for oc, nc in zip(OLD["cases"], NEW["cases"]):
    for mode in ("off", "auto"):
        for oe, ne in zip(oc["sweeps"][mode], nc["sweeps"][mode]):
            rows.append((oc["form"], mode, ne["tolerance_fraction_of_half_beam"] * 100,
                         oe["flips_after"], ne["flips_after"], ne["flips_before"],
                         ne["rejected_steps"]))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4.8))

changed = [r for r in rows if r[3] != r[4]]
labels = [f"{f.replace('_iii','').replace('_cb060','').replace('_class','')}\n{m} {b:.1f}%"
          for f, m, b, *_ in changed]
x = range(len(changed))
ax1.bar([i - 0.2 for i in x], [r[3] for r in changed], width=0.4,
        label="before the fix (last trial, rejected)", color="tab:red", alpha=0.85)
ax1.bar([i + 0.2 for i in x], [r[4] for r in changed], width=0.4,
        label="delivered table", color="tab:blue")
for i, r in enumerate(changed):
    ax1.text(i, max(r[3], r[4]) + 4, f"{r[3]}→{r[4]}", ha="center", fontsize=7)
ax1.set_xticks(list(x))
ax1.set_xticklabels(labels, fontsize=7, rotation=30, ha="right")
ax1.set_ylabel("curvature sign flips, three families")
ax1.set_title("published flips_after: discarded trial vs delivered table\n"
              f"{len(changed)} of {len(rows)} study runs had a rejected outer step", fontsize=9)
ax1.legend(fontsize=7)
ax1.grid(alpha=0.3, axis="y")

forms = [c["form"] for c in NEW["cases"]]
pairs = []
for c in NEW["cases"]:
    entry = next(e for e in c["sweeps"]["off"]
                 if abs(e["tolerance_fraction_of_half_beam"] - 0.008) < 1e-9)
    fmt = entry["flips_by_direction"]
    before = sum(v["before"] for v in fmt.values())
    after = sum(v["after"] for v in fmt.values())
    pairs.append((c["form"].replace("_iii", "").replace("_cb060", "").replace("_class", ""),
                  before, after))
ax2.bar([i - 0.2 for i in range(len(pairs))], [p[1] for p in pairs], width=0.4,
        label="input table (jitter 0.4 % B/2)", color="tab:gray")
ax2.bar([i + 0.2 for i in range(len(pairs))], [p[2] for p in pairs], width=0.4,
        label="delivered table", color="tab:green")
for i, (name, before, after) in enumerate(pairs):
    ax2.text(i, max(before, after) + 6, f"{before}→{after}", ha="center", fontsize=8)
ax2.set_xticks(range(len(pairs)))
ax2.set_xticklabels([p[0] for p in pairs])
ax2.set_ylabel("curvature sign flips")
ax2.set_title("per-family sum, plain energy at 0.8 % B/2 budget\n"
              "(before = input table, after = delivered table)", fontsize=9)
ax2.legend(fontsize=7)
ax2.grid(alpha=0.3, axis="y")
fig.tight_layout()
fig.savefig("/tmp/m22814_delivered_flips.png", dpi=150)
print("wrote /tmp/m22814_delivered_flips.png")
print(json.dumps([{"form": f, "mode": m, "budget_pct": b, "old": o, "delivered": nd,
                   "rejected": rj} for f, m, b, o, nd, _, rj in changed], indent=0))