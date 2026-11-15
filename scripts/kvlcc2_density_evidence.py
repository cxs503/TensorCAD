"""Before/after evidence for the KVLCC2 sample-density defect (anchor GAP-CAD-T1-IGES-PIPE).

Runs ``fleet_baseline.build_ship`` twice on the publisher archive, once with
``density_follows_adopted_factor=False`` (the shipped-before behaviour) and once with the
fix, and writes ``runs/kvlcc2_density/evidence.json`` + ``figure.png``.  Every number in
the report and the figure is read back from these two runs.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tensorcad.fleet_baseline import build_ship, specs_from_json

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs" / "kvlcc2_density"
OUT.mkdir(parents=True, exist_ok=True)
spec = {s.name: s for s in specs_from_json(ROOT / "benchmarks/verified/baseline_fleet/fleet_specs.json")}["KVLCC2"]

def run(flag: bool) -> dict:
    mesh, rep = build_ship(spec, dx=1.0, dz=0.1, density_follows_adopted_factor=flag)
    ts = rep["source"]["trim_sampling"]
    return {
        "density_follows_adopted_factor": flag,
        "verdict": rep["verdict"],
        "volume_m3": rep["derived"]["volume_m3"],
        "volume_error_pct": rep["derived"]["published_volume_rel_error_percent"],
        "cb": rep["derived"].get("cb"),
        "watertight": rep["checks"]["watertight"],
        "edge_manifold": rep["checks"]["edge_manifold"],
        "patches_used": ts["patches_used"],
        "patch_count": ts["patches_total"],
        "samples": ts.get("samples"),
        "patch_grids": [r.get("grid") for r in ts["patches"]],
        "patch_kept": [r.get("kept") for r in ts["patches"]],
        "dz_refinement": ts.get("dz_refinement"),
        "unit_gate": {k: v for k, v in rep["source"]["unit_gate"].items()
                      if k != "candidates"},
        "reasons": rep.get("reasons"),
        "triangles": rep["checks"].get("triangles"),
    }

before, after = run(False), run(True)
published = {"volume_m3": spec.published_volume_m3, "cb": spec.published_cb,
             "Lpp_m": spec.length_pp, "B_m": spec.beam, "T_m": spec.draft}
evidence = {"anchor": "GAP-CAD-T1-IGES-PIPE", "input": spec.iges_path,
            "iges_sha256": spec.iges_sha256, "published": published,
            "before": before, "after": after,
            "repro": "python scripts/kvlcc2_density_evidence.py"}
(OUT / "evidence.json").write_text(json.dumps(evidence, indent=2, default=str))

fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
a = ax[0][0]
labels = ["declared-factor\ndensity (before)", "adopted-factor\ndensity (after)"]
vals = [before["volume_m3"], after["volume_m3"]]
a.bar(labels, vals, color=["#b03030", "#2a7a3a"])
a.axhline(published["volume_m3"], color="k", ls="--", lw=1.2,
          label=f"published {published['volume_m3']:.0f} m$^3$")
for x, v in enumerate(vals):
    a.text(x, v, f"{v:,.0f}\n({100*(v-published['volume_m3'])/published['volume_m3']:+.2f}%)",
           ha="center", va="bottom", fontsize=9)
a.set_ylabel("displacement [m$^3$]"); a.set_title("KVLCC2 volume: unit-factor defect and fix")
a.legend(loc="center right", fontsize=8)

b = ax[0][1]
width = 0.38
idx = range(len(after["patch_kept"]))
b.bar([i - width/2 for i in idx], [g[0]*g[1] if g else 0 for g in before["patch_grids"]],
      width, color="#b03030", label="grid samples (before)")
b.bar([i - width/2 for i in idx], before["patch_kept"], 0.0001, color="#e08080")
b.bar([i + width/2 for i in idx], after["patch_kept"], width, color="#2a7a3a", label="kept samples (after)")
b.set_xlabel("trimmed patch index"); b.set_ylabel("samples inside the trim loops")
b.set_yscale("log"); b.set_title("sampling density per trimmed patch (19/19)")
b.legend(fontsize=8)

c = ax[1][0]
dz = after["dz_refinement"]["dz_m"]; vol = after["dz_refinement"]["volume_m3"]
c.plot(dz, vol, "o-", color="#2a7a3a", label="after fix")
c.axhline(published["volume_m3"], color="k", ls="--", lw=1.2, label="published")
b0 = before["dz_refinement"]
c.plot(b0["dz_m"], b0["volume_m3"], "s--", color="#b03030", label="before fix")
c.set_xscale("log"); c.set_xlabel("vertical cell dz [m]"); c.set_ylabel("displacement [m$^3$]")
c.set_title("dz refinement ladder"); c.legend(fontsize=8)

d = ax[1][1]
keys = ["volume_error_pct", "cb_error_pct"]
cb_before = 100 * (before["cb"] - published["cb"]) / published["cb"] if before["cb"] else float("nan")
cb_after = 100 * (after["cb"] - published["cb"]) / published["cb"]
d.bar([0, 1], [before["volume_error_pct"], cb_before], 0.35, color="#b03030", label="before (∇, Cb)")
d.bar([0.4, 1.4], [after["volume_error_pct"], cb_after], 0.35, color="#2a7a3a", label="after (∇, Cb)")
d.axhline(0, color="k", lw=0.8); d.axhline(-3, color="r", ls=":", lw=1, label="±3 % acceptance")
d.set_xticks([0.2, 1.2]); d.set_xticklabels(["displacement", "$C_B$"])
d.set_ylabel("relative error [%]"); d.set_title("error vs published (SIMMAN 2020)")
d.set_yscale("symlog", linthresh=1.0); d.legend(fontsize=8)
fig.suptitle("GAP-CAD-T1-IGES-PIPE · KVLCC2 publisher IGES → watertight solid (m21995)")
fig.tight_layout()
fig.savefig(OUT / "figure.png", dpi=130)
print(json.dumps({"before": {k: before[k] for k in ("volume_m3", "volume_error_pct", "cb", "verdict")},
                  "after": {k: after[k] for k in ("volume_m3", "volume_error_pct", "cb", "verdict")},
                  "figure": str(OUT / "figure.png")}, indent=1))
