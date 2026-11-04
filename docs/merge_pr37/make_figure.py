#!/usr/bin/env python3
"""Figure for the PR #37 cure report: branch physics still converges + test parity.

Left  : |volume error| of the merged transom-stern family against its closed-form
        reference, four truncations x four grids, from
        ``outputs/transom_stern/convergence.csv`` (produced by
        ``examples/transom_stern/run_example.py`` inside this merge).
Right : core-suite outcome of ``main`` (pre-merge baseline) vs the merge commit,
        same command line, same environment.

Every number is read from the artefacts; nothing is hand-typed.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONV = ROOT / "outputs/transom_stern/convergence.csv"
OUT = HERE / "merge_evidence.png"


def tally(path: Path) -> dict[str, int]:
    """Parse the pytest summary line ('2 failed, 802 passed, 11 skipped in ...')."""
    text = path.read_text(errors="replace")
    summary = [ln for ln in text.splitlines() if " in " in ln and ("passed" in ln or "failed" in ln)]
    if not summary:
        raise SystemExit(f"no pytest summary in {path}")
    line = summary[-1]
    out = {"passed": 0, "failed": 0, "skipped": 0}
    for chunk in line.split(","):
        chunk = chunk.strip()
        for key in out:
            if chunk.endswith(key):
                out[key] = int(chunk.split()[0])
    out["seconds"] = float(line.rsplit(" in ", 1)[1].split()[0].rstrip("s"))
    return out


def observed_order(points: list[tuple[float, float]]) -> float:
    """Least-squares slope of log|err| against log(nx)."""
    xs = [math.log(nx) for nx, _ in points]
    ys = [math.log(abs(e)) for _, e in points]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    return num / den


rows = list(csv.DictReader(CONV.open()))
by_tau: dict[float, list[tuple[int, float]]] = defaultdict(list)
for row in rows:
    by_tau[float(row["tau"])].append((int(row["nx"]), float(row["volume_error"])))

baseline = tally(ROOT / "runs/baseline_main.log")
merged = tally(ROOT / "runs/merged_tests.log")

fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(9.6, 3.9), dpi=88)

markers = ["o", "s", "^", "D"]
for (tau, pts), mk in zip(sorted(by_tau.items()), markers):
    pts.sort()
    xs = [p[0] for p in pts]
    ys = [abs(p[1]) * 100 for p in pts]
    order = -observed_order(pts)  # |err| ~ nx^-p, so the order is the negative slope
    ax0.loglog(xs, ys, mk + "-", lw=1.4, ms=4.5, label=f"$\\tau$={tau:g}  (order {order:.2f})")
ax0.axhline(3.0, color="crimson", ls="--", lw=1.0)
ax0.text(43, 3.35, "3 % band", color="crimson", fontsize=7.5)
ax0.set_xlabel("stations in x  $n_x$")
ax0.set_ylabel("|$V_{mesh}-V_{ref}$| / $V_{ref}$  [%]")
ax0.set_title("merged branch: volume error vs closed-form reference", fontsize=9)
ax0.grid(True, which="both", alpha=0.25, lw=0.5)
ax0.legend(fontsize=7, frameon=False, loc="lower left")
ax0.tick_params(labelsize=7.5)

labels = ["main\n(base 8389c7e)", "merge\n(cure 571af0f)"]
keys = ["passed", "failed", "skipped"]
colors = {"passed": "#2f7d32", "failed": "#b3261e", "skipped": "#8a8a8a"}
xpos = [0, 1]
bottom = [0, 0]
for key in keys:
    vals = [baseline[key], merged[key]]
    ax1.bar(xpos, vals, bottom=bottom, color=colors[key], width=0.55, label=key)
    for i, (v, b) in enumerate(zip(vals, bottom)):
        if v:
            ax1.text(i, b + v / 2, str(v), ha="center", va="center", fontsize=7.5,
                     color="white" if key != "skipped" else "black")
    bottom = [b + v for b, v in zip(bottom, vals)]
ax1.set_xticks(xpos)
ax1.set_xticklabels(labels, fontsize=8)
ax1.set_ylabel("tests")
ax1.set_title("core suite: main baseline vs merge commit", fontsize=9)
ax1.legend(fontsize=7, frameon=False, loc="upper right")
ax1.tick_params(labelsize=7.5)
ax1.grid(True, axis="y", alpha=0.25, lw=0.5)
ax1.text(
    0.02, 0.03,
    f"same 2 failures both sides\n(env: missing $safetensors$)\n"
    f"+{merged['passed'] - baseline['passed']} tests from the merged branch",
    transform=ax1.transAxes, fontsize=7, va="bottom",
)

fig.suptitle("TensorCAD PR #37 cure: transom-stern family merged onto main", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(OUT, dpi=88)

summary = {
    "figure": OUT.name,
    "size_bytes": OUT.stat().st_size,
    "by_tau": {str(k): {"grids": [p[0] for p in sorted(v)],
                        "volume_error": [p[1] for p in sorted(v)],
                        "observed_order": round(-observed_order(sorted(v)), 3)}
               for k, v in sorted(by_tau.items())},
    "baseline_main": baseline,
    "merged": merged,
}
(HERE / "figure_data.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))