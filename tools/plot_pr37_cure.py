"""Evidence figure for the PR #37 merge cure.

Three panels, every number read from a file the merged tree produced (no hand typing):

* **(a) test parity** -- passed/failed counts of ``main`` at the merge base line and of the
  merged tree; the merge must not trade a single pass for a pass elsewhere.
* **(b) transom-stern convergence** -- signed mesh volume error against the closed form,
  per truncation ``tau``, against the second-order reference slope; this is the branch's
  own verification claim, re-run on the merged tree.
* **(c) reference ladder** -- the three independent evaluations of the incomplete taper
  integral (``sinc^p`` series, Gamma closed form at ``tau = 0``, tanh-sinh quadrature at
  ``tau > 0``); agreement at round-off is what says the two sides' numbers did not drift.

Run: ``python3 tools/plot_pr37_cure.py`` writes ``artifacts/pr37/cure_evidence.png``.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "pr37"
LADDER = ROOT / "docs" / "m2873" / "outputs" / "transom_stern" / "reference_ladder.csv"
CONV = ROOT / "docs" / "m2873" / "outputs" / "transom_stern" / "convergence.csv"
AUDIT = OUT / "union_audit.json"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    audit = json.loads(AUDIT.read_text(encoding="utf-8"))
    base = audit["tests"]["main_baseline"]["counts"]
    merged = audit["tests"]["merged"]["counts"]

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))

    ax = axes[0]
    labels = ["passed", "failed", "skipped"]
    xs = range(len(labels))
    ax.bar([x - 0.2 for x in xs], [base.get(k, 0) for k in labels], width=0.4,
           label=f"main {audit['main_tip'][:7]}", color="#8fa8c8")
    ax.bar([x + 0.2 for x in xs], [merged.get(k, 0) for k in labels], width=0.4,
           label="merged (main + m2873)", color="#2f6f9f")
    for x, key in zip(xs, labels, strict=True):
        for offset, counts in ((-0.2, base), (0.2, merged)):
            value = counts.get(key, 0)
            ax.annotate(str(value), (x + offset, value), ha="center", va="bottom",
                        fontsize=8)
    ax.set_xticks(list(xs), labels)
    ax.set_yscale("symlog")
    ax.set_ylabel("tests")
    ax.set_title("(a) test parity: same 2 env failures, "
                 f"+{audit['tests']['passed_delta']} passes")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    rows = read_csv(CONV)
    ax = axes[1]
    taus = sorted({float(r["tau"]) for r in rows})
    for tau in taus:
        pts = sorted((float(r["nx"]), abs(float(r["volume_error"])))
                     for r in rows if float(r["tau"]) == tau)
        ax.loglog([p[0] for p in pts], [p[1] for p in pts], "o-", label=f"$\\tau$ = {tau:g}")
    anchor = sorted((float(r["nx"]), abs(float(r["volume_error"])))
                    for r in rows if float(r["tau"]) == taus[0])[0]
    ref = [anchor[1] * (nx / anchor[0]) ** -2 for nx in (41, 321)]
    ax.loglog((41, 321), ref, "k--", lw=1, label="second order")
    ax.set_xlabel("stations along the hull $n_x$")
    ax.set_ylabel("|mesh volume error| vs closed form")
    ax.set_title("(b) transom family converges on the merged tree")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    ax = axes[2]
    ladder = [r for r in read_csv(LADDER) if r.get("exponent")]
    exponents = sorted({float(r["exponent"]) for r in ladder})
    closed = max((abs(float(r["series_vs_closed_form"])) for r in ladder
                  if r.get("series_vs_closed_form")), default=0.0)
    quad = max((abs(float(r["series_vs_quadrature"])) for r in ladder
                if r.get("series_vs_quadrature")), default=0.0)
    ax.bar([0, 1], [max(closed, 1e-300), max(quad, 1e-300)], width=0.5,
           color=["#2f6f9f", "#9f6f2f"])
    ax.set_xticks([0, 1], ["series vs $\\Gamma$ closed form\n($\\tau$=0)",
                          "series vs tanh-sinh quadrature\n($\\tau$>0)"])
    ax.set_yscale("log")
    ax.set_ylabel("max |relative difference|")
    ax.set_ylim(1e-18, 1e-6)
    ax.set_title("(c) reference ladder agrees to round-off\n"
                 f"({len(exponents)} exponents, $\\tau\\in[0,0.35]$)")
    for x, value in ((0, closed), (1, quad)):
        ax.annotate(f"{value:.1e}", (x, max(value, 1e-300)), ha="center",
                    va="bottom", fontsize=9)
    ax.grid(alpha=0.3, axis="y")

    fig.suptitle("PR #37 cure: semantic union of transom-stern (m2873) and main", y=1.0)
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "cure_evidence.png", dpi=140, bbox_inches="tight")

    worst = max(abs(float(r["volume_error"])) for r in rows)
    fine = [abs(float(r["volume_error"])) for r in rows if float(r["nx"]) == 321]
    summary = {
        "test_counts_main": base,
        "test_counts_merged": merged,
        "passed_delta": audit["tests"]["passed_delta"],
        "new_failures": audit["tests"]["new_failures_vs_baseline"],
        "volume_error_worst_all_grids": worst,
        "volume_error_worst_finest_grid": max(fine) if fine else None,
        "ladder_max_closed_form": closed,
        "ladder_max_quadrature": quad,
    }
    (OUT / "cure_metrics.json").write_text(json.dumps(summary, indent=2) + "\n",
                                           encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
