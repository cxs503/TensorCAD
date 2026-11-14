#!/usr/bin/env python3
"""Figure for docs/baseline_fleet_type126_falsified.md, drawn from the audit JSON."""
import json, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

report = json.loads(Path(sys.argv[1]).read_text())
out = Path(sys.argv[2] if len(sys.argv) > 2 else "docs/figs/t126_trim_audit.png")
rows = report["files"]
names = [r["file"].split("_")[0].replace(".igs", "") for r in rows]
surf = [r["census"]["surfaces_type128"] for r in rows]
curv = [r["census"]["curves_type126"] for r in rows]
trim = [r["census"]["basis_surfaces_trimmed"] for r in rows]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.0))
ax = axes[0]
x = range(len(rows)); w = 0.27
ax.bar([i - w for i in x], surf, w, label="type 128 surfaces (decoded)", color="#2777b1")
ax.bar(list(x), curv, w, label="type 126 curves (decoded)", color="#e08a1e")
ax.bar([i + w for i in x], trim, w, label="basis surfaces bound to a trim", color="#2c8f57")
ax.set_xticks(list(x)); ax.set_xticklabels(names, fontsize=9)
ax.set_ylabel("entities in file"); ax.set_yscale("symlog")
ax.set_title("Entity identities, publisher archives", fontsize=10)
ax.legend(fontsize=7, frameon=False); ax.grid(axis="y", alpha=.3)

ax = axes[1]
kvl = [r for r in rows if "kvlcc2" in r["file"]][0]["trim_pointer_audit"]["patches"]
patches = [p for p in kvl if p.get("inside_fraction")]
idx = [p["patch"] for p in patches]
bptr = [p["inside_fraction"]["field_bptr"] or 0.0 for p in patches]
crv = [p["inside_fraction"]["field_crv"] or 0.0 for p in patches]
ax.plot(idx, bptr, "o-", ms=4, color="#2c8f57", label="BPTR (field 4, index 3)")
ax.plot(idx, crv, "s--", ms=3, color="#b1342f", label="CRV (field 2, index 1)")
ax.axhline(0.5, color="grey", lw=.8, ls=":")
ax.set_xlabel("KVLCC2 trimmed patch"); ax.set_ylabel("samples inside (u,v) domain")
ax.set_ylim(-0.05, 1.05); ax.set_title("Which pointer is the parameter-space loop", fontsize=10)
ax.legend(fontsize=7, frameon=False); ax.grid(alpha=.3)
fig.tight_layout(); out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=150)
print("wrote", out)
