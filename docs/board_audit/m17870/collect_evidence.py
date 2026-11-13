#!/usr/bin/env python3
"""m17870 board audit of m17833 / GAP-CAD-D2-GOLDEN-REGRESSION: collect evidence, no hand-typed numbers.

Inputs (all read, nothing typed):
  docs/board_audit/m17870/benchmark_rows.json   tensor_benchmark_cases rows id 65-68 (the judge)
  /tmp/audit17870/bench/result.json             fresh audit re-run  (reproduced from the audited tree)
  <work>/artifacts/golden_benchmark/result.json committed artifact of m17833
Outputs:
  docs/board_audit/m17870/audit_reproduction.json
  docs/board_audit/m17870/audit_golden.png
"""
from __future__ import annotations

import json
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
AUDIT = os.environ.get("AUDIT_TREE", "/root/tp-agent/work/m17833")
FRESH = os.environ.get("FRESH_RESULT", "/tmp/audit17870/bench/result.json")

ROWS = json.load(open(os.path.join(HERE, "benchmark_rows.json")))
FRESH_RES = json.load(open(FRESH))
SHIPPED = json.load(open(os.path.join(AUDIT, "artifacts", "golden_benchmark", "result.json")))
GATES = json.load(open(os.path.join(AUDIT, "artifacts", "golden_benchmark", "gates_verified.json")))

METRIC_RE = re.compile(r"([A-Za-z0-9_]+): ref ([-\d.e+]+) vs got ([-\d.e+]+) -> ([-\d.e+]+)")
TOL_RE = re.compile(r"worst-case tolerance ([-\d.e+]+) \(([A-Za-z0-9_]+)\)")

# ---- 1. decode the judge rows (tensor_benchmark_cases) ----------------------
judge = []
for row in ROWS:
    parsed = METRIC_RE.findall(row["metric"])
    worst_tol = float(TOL_RE.search(row["tolerance"]).group(1))
    worst_metric = TOL_RE.search(row["tolerance"]).group(2)
    judge.append(
        {
            "case_row_id": row["id"],
            "case_name": row["case_name"],
            "status": row["status"],
            "metrics_parsed": len(parsed),
            "worst_case_tolerance": worst_tol,
            "worst_case_tolerance_metric": worst_metric,
            "max_abs_rel_error_recorded": max(abs(float(p[3])) for p in parsed),
            "all_within_tolerance": all(abs(float(p[3])) <= worst_tol for p in parsed),
            "case_slug": row["case_name"].replace("CASE-CAD-D2-GOLDEN-", "").lower(),
        }
    )

# ---- 2. fresh run rows vs shipped rows vs judge ----------------------------
fresh_rows = {(r["case"], r["metric"]): r for r in FRESH_RES["rows"]}
shipped_rows = {(r["case"], r["metric"]): r for r in SHIPPED["rows"]}

BIT_IDENTICAL = []
for k, a in shipped_rows.items():
    b = fresh_rows.get(k)
    BIT_IDENTICAL.append(
        {
            "case": k[0],
            "metric": k[1],
            "shipped": a["measured"],
            "fresh": None if b is None else b["measured"],
            "abs_diff": None if b is None else abs(a["measured"] - b["measured"]),
        }
    )
max_diff = max(x["abs_diff"] for x in BIT_IDENTICAL)

# judge got values appear verbatim in the shipped artifact?
judge_vs_artifact = []
for row in ROWS:
    slug = row["case_name"].replace("CASE-CAD-D2-GOLDEN-", "").lower()
    for name, ref, got, rel in METRIC_RE.findall(row["metric"]):
        key = next((k for k in shipped_rows if k[0] == slug and k[1] == name), None)
        if key is None:
            judge_vs_artifact.append({"case": slug, "metric": name, "in_artifact": False})
            continue
        m = shipped_rows[key]["measured"]
        got_f = float(got)
        judge_vs_artifact.append(
            {
                "case": slug,
                "metric": name,
                "in_artifact": True,
                "judge_got": got_f,
                "artifact_measured": m,
                "abs_diff": abs(m - got_f),
                "judge_rel": float(rel),
                "artifact_rel": shipped_rows[key]["relative_error"],
            }
        )
missing = [j for j in judge_vs_artifact if not j["in_artifact"]]
max_judge_diff = max((j["abs_diff"] for j in judge_vs_artifact if j["in_artifact"]), default=None)

# ---- 3. accuracy against the reference (platform line: >15% unacceptable) --
analytic = [r for r in FRESH_RES["rows"] if r["kind"] == "analytic"]
worst_analytic = max(analytic, key=lambda r: abs(r["relative_error"]))
worst_any = max(FRESH_RES["rows"], key=lambda r: abs(r["relative_error"]))
over15 = [r for r in FRESH_RES["rows"] if abs(r["relative_error"]) > 0.15]
over3 = [r for r in analytic if abs(r["relative_error"]) > 0.03]

out = {
    "audit_mission": 17870,
    "audited_mission": 17833,
    "anchor": "GAP-CAD-D2-GOLDEN-REGRESSION",
    "audited_tree": AUDIT,
    "judge": {
        "table": "tensor_benchmark_cases",
        "anchor_slug": "GAP-CAD-D2-GOLDEN-REGRESSION",
        "rows": judge,
        "all_status_passed": all(j["status"] == "passed" for j in judge),
        "n_cases": len(judge),
        "n_metrics": sum(j["metrics_parsed"] for j in judge),
        "all_within_tolerance": all(j["all_within_tolerance"] for j in judge),
    },
    "reproduction": {
        "command": "PYTHONPATH=. python3 benchmarks/verified/golden_geometry/run_benchmark.py --out /tmp/audit17870/bench",
        "rows": len(FRESH_RES["rows"]),
        "gates": [(g["gate"], g["verdict"]) for g in FRESH_RES["gates"]],
        "passed": FRESH_RES["passed"],
        "seconds": FRESH_RES["seconds"],
        "max_abs_diff_vs_shipped_artifact": max_diff,
        "bit_identical_to_shipped": max_diff == 0.0,
        "judge_values_missing_from_artifact": missing,
        "max_abs_diff_judge_vs_artifact": max_judge_diff,
        "geometry_fingerprint": {
            "shipped": SHIPPED["geometry_fingerprint"],
            "fresh": FRESH_RES["geometry_fingerprint"],
            "equal": SHIPPED["geometry_fingerprint"] == FRESH_RES["geometry_fingerprint"],
        },
        "git_commit_shipped_artifact": SHIPPED["git_commit"],
        "git_commit_fresh": FRESH_RES["git_commit"],
    },
    "accuracy": {
        "platform_line": 0.15,
        "analytic_metrics": len(analytic),
        "max_abs_rel_error_analytic": abs(worst_analytic["relative_error"]),
        "worst_analytic": [worst_analytic["case"], worst_analytic["metric"]],
        "max_abs_rel_error_any": abs(worst_any["relative_error"]),
        "metrics_over_3pct_analytic": [[r["case"], r["metric"], r["relative_error"]] for r in over3],
        "metrics_over_15pct": [[r["case"], r["metric"], r["relative_error"]] for r in over15],
        "worst_tolerance_any": worst_any["tolerance"],
    },
    "gates_shipped": GATES["gates"],
    "gates_shipped_passed": GATES["passed"],
}

# ---- figure ---------------------------------------------------------------
fresh_rows_sorted = sorted(FRESH_RES["rows"], key=lambda r: (r["case"], -abs(r["relative_error"])))
labels = [f"{r['case'][:9]}.{r['metric'][:16]}" for r in fresh_rows_sorted]
errs = [abs(r["relative_error"]) for r in fresh_rows_sorted]
tols = [r["tolerance"] for r in fresh_rows_sorted]
colors = ["#2b6cb0" if r["kind"] == "analytic" else "#718096" for r in fresh_rows_sorted]

fig, (ax, bx) = plt.subplots(1, 2, figsize=(13.5, 5.2), gridspec_kw={"width_ratios": [2.1, 1]})
x = range(len(errs))
ax.bar(x, [max(e, 1e-12) for e in errs], color=colors, label="|measured - reference| / |reference|")
ax.plot(list(x), tols, "r_", ms=9, label="derived tolerance (2 x coarse->fine, 1e-4 floor)")
ax.axhline(0.03, color="crimson", ls=":", lw=1.2, label="3% precision line")
ax.axhline(0.15, color="black", ls="--", lw=1.0, label="15% platform reject line")
ax.set_yscale("log")
ax.set_xticks(list(x))
ax.set_xticklabels(labels, rotation=90, fontsize=5.5)
ax.set_ylabel("relative error (log)")
ax.set_title(
    f"m17833 GAP-CAD-D2 golden geometry — {len(errs)} rows re-run by the audit\n"
    f"worst {abs(worst_any['relative_error']):.3e} ({worst_any['case']}/{worst_any['metric']}) vs tolerance {worst_any['tolerance']:.3e}"
)
ax.grid(True, which="both", axis="y", alpha=0.25)
ax.legend(fontsize=6.5, loc="lower left")

names = [j["case_slug"] for j in judge]
bx.bar(names, [j["max_abs_rel_error_recorded"] for j in judge], color="#2f855a", label="max |rel err| in judge row")
bx.plot(names, [j["worst_case_tolerance"] for j in judge], "r_", ms=18, label="judge worst tolerance")
bx.set_yscale("log")
bx.set_ylabel("relative")
bx.set_title("tensor_benchmark_cases rows 65-68 (the judge)")
bx.legend(fontsize=6.5)
bx.grid(True, which="both", axis="y", alpha=0.25)
for i, j in enumerate(judge):
    bx.text(i, j["max_abs_rel_error_recorded"] * 1.6, f"{j['metrics_parsed']}m", ha="center", fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(HERE, "audit_golden.png"), dpi=150)
plt.close(fig)

json.dump(out, open(os.path.join(HERE, "audit_reproduction.json"), "w"), indent=1, sort_keys=True)
print(json.dumps({k: out[k] for k in ("accuracy",)}, indent=1))
print("judge", out["judge"]["n_cases"], "cases,", out["judge"]["n_metrics"], "metrics, all_within_tolerance",
      out["judge"]["all_within_tolerance"])
print("repro bit-identical:", out["reproduction"]["bit_identical_to_shipped"],
      "max_diff", max_diff, "judge-vs-artifact max_diff", max_judge_diff, "missing", len(missing))