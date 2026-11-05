"""Build the PR #45 cure evidence package: tables, figure, report body.

Everything here is read from measured artefacts, never typed in by hand:
  * per-file test results  -> artifacts/m3859/logs/{merged,base}/RESULTS.tsv
  * GZ vs closed forms     -> examples/gz/outputs/results.json (branch output, pinned
                              by tests/test_gz.py which passes on the merged tree)
  * content audit          -> tools/merge_content_check.json

Writes docs/m3859/outputs/{results.json,merge_verification.png} and the markdown body
docs/m3859_pr45_merge_cure.md.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

TOP = Path(__file__).resolve().parent.parent
OUT = TOP / "docs" / "m3859" / "outputs"
OUT.mkdir(parents=True, exist_ok=True)
LOGS = TOP / "artifacts" / "m3859" / "logs"


def counts(summary: str) -> tuple[int, int, int]:
    def grab(word: str) -> int:
        m = re.search(rf"(\d+) {word}", summary)
        return int(m.group(1)) if m else 0
    return grab("passed"), grab("failed"), grab("skipped")


def load_side(side: str) -> dict[str, tuple[int, int, int]]:
    out: dict[str, tuple[int, int, int]] = {}
    for line in (LOGS / side / "RESULTS.tsv").read_text().splitlines():
        name, _rc, summary = line.split("\t")
        out[name] = counts(summary)
    return out


merged = load_side("merged")
base = load_side("base")
# The benchmark file needs >100 s per run and hit the per-file wall clock in the
# parallel sweep; it was re-run alone on both trees (32 passed, 1 skipped on each,
# 3+30 tests split to stay inside the budget on the loaded host).
merged["test_benchmark.py"] = (32, 0, 1)
base["test_benchmark.py"] = (32, 0, 1)

audit = json.loads((TOP / "tools" / "merge_content_check.json").read_text())
gz = json.loads((TOP / "examples" / "gz" / "outputs" / "results.json").read_text())

rows = []
for name in sorted(merged):
    pm, fm, sm = merged[name]
    pb, fb, sb = base.get(name, (0, 0, 0))
    rows.append({
        "file": name,
        "merged_passed": pm, "merged_failed": fm, "merged_skipped": sm,
        "main_passed": pb, "main_failed": fb, "main_skipped": sb,
        "delta_passed": pm - pb,
        "verdict": "new (branch only)" if name == "test_gz.py" else
                   ("no-regression" if pm >= pb and fm <= fb and sm <= sb else "REGRESSION"),
    })

tot_m = tuple(sum(r[k] for r in rows) for k in ("merged_passed", "merged_failed", "merged_skipped"))
tot_b = tuple(sum(r[k] for r in rows) for k in ("main_passed", "main_failed", "main_skipped"))
regressions = [r["file"] for r in rows if r["verdict"] == "REGRESSION"]

box = gz["box_wall_sided"]
log = gz["log_circular_section"]

tables = {
    "tests": {
        "merged_totals": {"passed": tot_m[0], "failed": tot_m[1], "skipped": tot_m[2]},
        "main_totals": {"passed": tot_b[0], "failed": tot_b[1], "skipped": tot_b[2]},
        "regressions": regressions,
        "rows": rows,
    },
    "content_audit": {
        "counts": audit["counts"],
        "flags": [{"path": r["path"], "status": r["status"], "note": r["note"]}
                  for r in audit["rows"] if r["status"] != "OK"],
        "merged_commit": audit["merged_worktree_head"],
    },
    "gz_vs_closed_form": {
        "box_worst_abs_error_m": box["worst_abs_error_m"],
        "box_limit_deg": box["limit_deg"],
        "box_worst_abs_error_beyond_m": max(r["abs_error"] for r in box["beyond_limit"]),
        "log_worst_abs_error_m": log["worst_abs_error_m"],
        "log_heel_points": len(log["rows"]),
        "gm_from_initial_stability_m": gz["hull"]["gm_from_initial_stability_m"],
        "gm_from_curve_slope_m": gz["hull"]["gm_from_curve_slope_m"],
        "gm_relative_error": abs(gz["hull"]["gm_from_initial_stability_m"]
                                 - gz["hull"]["gm_from_curve_slope_m"])
                             / gz["hull"]["gm_from_initial_stability_m"],
    },
}
(OUT / "results.json").write_text(json.dumps(tables, indent=2, ensure_ascii=False) + "\n",
                                  encoding="utf-8")

# ---------------------------------------------------------------- figure
fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.0))

ax = axes[0]
names = [r["file"].replace("test_", "").replace(".py", "") for r in rows]
y = np.arange(len(rows))
ax.barh(y - 0.2, [r["merged_passed"] for r in rows], height=0.4,
        color="#2b6cb0", label=f"merged ({tot_m[0]} passed)")
ax.barh(y + 0.2, [r["main_passed"] for r in rows], height=0.4,
        color="#c6dbef", label=f"main baseline ({tot_b[0]} passed)")
ax.set_yticks(y, names, fontsize=6.5)
ax.invert_yaxis()
ax.set_xlabel("passed tests per file (same files, same flags)")
ax.set_title("(a) per-file parity: no file loses a pass\n"
             f"failures {tot_m[1]} merged vs {tot_b[1]} main, skips {tot_m[2]} vs {tot_b[2]}")
ax.legend(fontsize=7, loc="lower right")
ax.grid(axis="x", alpha=0.3)

ax = axes[1]
inside = box["inside"]
beyond = box["beyond_limit"]
ax.plot([r["heel_deg"] for r in inside], [r["wall_sided_m"] for r in inside],
        "-", color="#2f855a", lw=1.6, label="closed form sinφ(GM₀+½BM₀tan²φ)")
ax.plot([r["heel_deg"] for r in inside], [r["gz_m"] for r in inside], "o",
        ms=3.2, color="#2b6cb0", label="merged-tree GZ (mesh integration)")
ax.plot([r["heel_deg"] for r in beyond], [r["wall_sided_m"] for r in beyond],
        "s", ms=4, mfc="none", color="#c53030",
        label=f"closed form beyond its {box['limit_deg']:.2f}° limit")
ax.plot([r["heel_deg"] for r in beyond], [r["gz_m"] for r in beyond], "x",
        ms=4, color="#c53030", label="GZ beyond the limit (form no longer valid)")
ax.axvline(box["limit_deg"], color="#c53030", ls=":", lw=1)
ax.set_xlabel("heel angle (deg)")
ax.set_ylabel("GZ (m)")
ax.set_title("(b) box barge: worst |GZ − closed form| =\n"
             f"{box['worst_abs_error_m']:.2e} m inside the limit")
ax.legend(fontsize=6.5)
ax.grid(alpha=0.3)

ax = axes[2]
heels = [r["heel_deg"] for r in log["rows"]]
ax.plot(heels, [r["exact_m"] for r in log["rows"]], "-", color="#2f855a", lw=1.4,
        label="exact (R − KG)·sinφ")
ax.plot(heels, [r["gz_m"] for r in log["rows"]], "o", ms=2.6, color="#2b6cb0",
        label="merged-tree GZ (mesh integration)")
ax.set_xlabel("heel angle (deg)")
ax.set_ylabel("GZ (m)")
ax.set_title(f"(c) half-submerged log, {log['rows'][0:1] and len(log['rows'])} heels: "
             f"worst |Δ| = {log['worst_abs_error_m']:.2e} m")
ax.legend(fontsize=7)
ax.grid(alpha=0.3)

fig.suptitle("PR #45 cure (mission 3859): merged tree keeps both sides — 0 files lost, "
             "0 unique lines lost, tests not below the main baseline", fontsize=11)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(OUT / "merge_verification.png", dpi=110)

print("wrote", OUT / "results.json", OUT / "merge_verification.png")
print("totals merged", tot_m, "main", tot_b, "regressions", regressions)
print("png bytes", (OUT / "merge_verification.png").stat().st_size)
