"""Report figure for m24408 (TensorCAD#746 see-through round 2).

Left: the 5415 displacement teeth — measured value, honest published reference and the
10 %-corrupted reference against the declared 3 % band (the gate's reference_mismatch).
Right: GCI — Richardson-extrapolated error vs fine-level error for the four
sphere/ellipsoid quantities (the extrapolated value must be closer to the reference).
All numbers are read from the committed baseline, not typed.
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path("/root/tp-agent/work/m24408")
d = json.loads((ROOT / "docs/golden/tensorcad_geometry_golden.json").read_text())
dt = d["cases"]["dtmb5415"]
se = d["cases"]["sphere_ellipsoid"]

meas = dt["measurements"]["fine"]["displacement_volume_m3"]
ref = dt["reference"]["displacement_volume_m3"]
band = dt["tolerances"]["displacement_volume_m3"]["value"]
corrupted = ref * 1.10

fig, (left, right) = plt.subplots(1, 2, figsize=(11.5, 4.4))

# ---- left: 5415 displacement teeth -------------------------------------------
lo, hi = ref * (1 - band), ref * (1 + band)
left.axhspan(lo, hi, color="#2ca02c", alpha=0.18, label="declared 3 % band")
left.plot([meas, ref], [0, 0], marker="o", ms=8, color="#1f77b4", lw=2,
          label=f"measured {meas:.1f} m$^3$")
left.plot([ref, ref], [-0.012, 0.012], marker="|", ms=14, color="#1f77b4", lw=2)
left.plot([corrupted, corrupted], [-0.012, 0.012], marker="|", ms=14, color="#d62728", lw=2)
left.annotate(f"published ref {ref:.1f}", (ref, 0.016), ha="center", fontsize=9, color="#1f77b4")
left.annotate(f"corrupted (+10 %) {corrupted:.1f}", (corrupted, 0.016), ha="center",
             fontsize=9, color="#d62728")
left.annotate(f"measured {meas:.1f}\n(err vs corrupted {abs(meas - corrupted) / corrupted:+.2%})",
              (meas, -0.016), ha="center", fontsize=9, color="#1f77b4")
left.set_yticks([])
left.set_xlabel("DTMB 5415 displacement volume (full scale, m$^3$)")
left.set_title("corrupted table: the gate's reference_mismatch fires")
left.legend(fontsize=8, loc="lower right")
left.grid(alpha=0.3, axis="x")

# ---- right: GCI extrapolation -------------------------------------------------
keys = ("sphere_volume_m3", "sphere_surface_area_m2",
        "ellipsoid_volume_m3", "ellipsoid_surface_area_m2")
labels = ["sphere\nvolume", "sphere\narea", "ellipsoid\nvolume", "ellipsoid\narea"]
ext_err, fine_err = [], []
for k in keys:
    r = se["reference"][k]
    ext_err.append(abs(se["extrapolated"][k]["value"] - r) / r)
    fine_err.append(abs(se["measurements"]["fine"][k] - r) / r)
x = range(len(keys))
left_b = right.bar([i - 0.19 for i in x], ext_err, width=0.38, color="#1f77b4",
                   label="Richardson extrapolated")
right.bar([i + 0.19 for i in x], fine_err, width=0.38, color="#ff7f0e", label="fine level")
right.set_yscale("log")
right.set_xticks(list(x))
right.set_xticklabels(labels, fontsize=9)
right.set_ylabel("relative error vs analytic reference")
right.set_title("GCI: extrapolated error < fine error (all four)")
right.legend(fontsize=8)
right.grid(alpha=0.3, axis="y")

fig.suptitle("m24408: the m21410 golden deliverable re-certified on current main", fontsize=11)
fig.tight_layout()
out = ROOT / "work" / "m24408_report_figure.png"
fig.savefig(out, dpi=150)
print("saved", out)

# print the numbers for the report table
print("5415: measured=%.4f honest_ref=%.4f corrupted=%.4f band=%.3f" % (meas, ref, corrupted, band))
print("5415 err vs honest ref: %+.4f%%" % (100 * (meas - ref) / ref))
print("5415 err vs corrupted: %+.4f%%" % (100 * (meas - corrupted) / corrupted))
for k, l in zip(keys, labels):
    r = se["reference"][k]
    print(f"{k}: ext={abs(se['extrapolated'][k]['value'] - r) / r:.3e} fine={abs(se['measurements']['fine'][k] - r) / r:.3e}")
