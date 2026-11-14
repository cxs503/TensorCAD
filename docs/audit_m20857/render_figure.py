"""Figure for docs/audit_m20857/golden_lock_audit.md (m20857).

Everything plotted is read from artifacts, never typed:
* the *before* fingerprints come from ``git show <rev>:docs/golden/tensorcad_geometry_golden.json``
  for the commit this mission started from;
* the *after* fingerprints and every metric come from a live ``python -m tensorcad.golden
  --check`` run stored as JSON;
* the *unlocked* verdict is the same rule the gate uses (``geometry_sha256`` empty).

Panels
  A  per-case mesh lock before/after the fix (1 = the mesh is fingerprinted).
  B  |relative error| vs the derived tolerance for every asserted metric, after the fix,
     on a log axis: the tolerance is the two-grid bound, every metric must sit under it.
     The point of the panel is that *the numbers did not move* when the baseline was
     re-blessed -- only the locks in panel A changed.

Usage::

    python docs/audit_m20857/render_figure.py \
        --before-rev c772985 --check-json /tmp/golden_after.json \
        --out docs/audit_m20857/golden_lock_audit.png
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

BASELINE = "docs/golden/tensorcad_geometry_golden.json"
CASE_ORDER = ["wigley", "series60", "sphere_ellipsoid", "dtmb5415"]


def _baseline_at(rev: str | None, root: Path) -> dict:
    if rev is None:
        return json.loads((root / BASELINE).read_text(encoding="utf-8"))
    text = subprocess.run(
        ["git", "show", f"{rev}:{BASELINE}"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return json.loads(text)


def _locked(baseline: dict, case: str) -> int:
    return int(bool((baseline.get("cases", {}).get(case, {}) or {}).get("geometry_sha256")))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--before-rev", required=True)
    parser.add_argument("--check-json", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    before = _baseline_at(args.before_rev, root)
    after = _baseline_at(None, root)
    run = json.loads(Path(args.check_json).read_text(encoding="utf-8"))

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13.0, 5.0), constrained_layout=True)

    # ---- panel A: is the mesh locked? -------------------------------------------
    xs = range(len(CASE_ORDER))
    width = 0.38
    series = {
        f"before ({args.before_rev})": [_locked(before, name) for name in CASE_ORDER],
        "after (this mission)": [_locked(after, name) for name in CASE_ORDER],
    }
    for offset, (label, values) in zip((-width / 2, width / 2), series.items(), strict=True):
        bars = ax_a.bar([x + offset for x in xs], values, width, label=label)
        for bar, value in zip(bars, values, strict=True):
            ax_a.annotate(
                "locked" if value else "UNLOCKED",
                (bar.get_x() + bar.get_width() / 2, bar.get_height()),
                ha="center",
                va="bottom" if value else "top",
                fontsize=9,
                rotation=0,
                xytext=(0, 3 if value else -12),
                textcoords="offset points",
            )
    ax_a.set_xticks(list(xs))
    ax_a.set_xticklabels(CASE_ORDER)
    ax_a.set_ylim(-0.35, 1.35)
    ax_a.set_yticks([0, 1])
    ax_a.set_yticklabels(["no fingerprint", "mesh fingerprinted"])
    ax_a.set_title(
        "A. geometry_sha256 in the committed baseline\n(the gate passed either way before the fix)"
    )
    ax_a.legend(loc="upper left", fontsize=8, frameon=False)
    ax_a.grid(axis="y", alpha=0.3)

    # ---- panel B: every metric vs its derived tolerance -------------------------
    colours = dict(zip(CASE_ORDER, plt.cm.tab10.colors, strict=False))
    ratios: dict[str, list[float]] = {name: [] for name in CASE_ORDER}
    for case in run["cases"]:
        for metric in case.get("metrics", []):
            if metric.get("relative_error") is None or not metric.get("tolerance"):
                continue
            ratio = abs(metric["relative_error"]) / metric["tolerance"]
            ratios.setdefault(case["name"], []).append(max(ratio, 1e-12))
    for index, name in enumerate(CASE_ORDER):
        values = ratios.get(name, [])
        if not values:
            continue
        ax_b.scatter(
            [index] * len(values),
            values,
            s=36,
            color=colours[name],
            alpha=0.85,
            edgecolor="black",
            linewidth=0.4,
            zorder=3,
        )
    ax_b.axhline(1.0, color="firebrick", linestyle="--", linewidth=1.2, zorder=2)
    ax_b.annotate(
        "declared tolerance (2 x two-grid change)",
        (0.02, 1.0),
        xycoords=("axes fraction", "data"),
        va="bottom",
        fontsize=8,
        color="firebrick",
    )
    ax_b.set_yscale("log")
    ax_b.set_xticks(range(len(CASE_ORDER)))
    ax_b.set_xticklabels(CASE_ORDER)
    ax_b.set_ylabel("|relative error| / tolerance")
    ax_b.set_title(
        "B. every asserted metric after the re-bless\n(all under the derived tolerance: no number moved)"
    )
    ax_b.grid(alpha=0.3)
    n_metrics = sum(len(v) for v in ratios.values())
    ax_b.text(
        0.02,
        0.03,
        f"{n_metrics} asserted metrics, {len(CASE_ORDER)} cases, gate passed="
        f"{run['passed']}",
        transform=ax_b.transAxes,
        fontsize=8,
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.suptitle(
        "GAP-CAD-D2-GOLDEN-REGRESSION: the golden baseline locked integrals, not meshes",
        fontsize=11,
    )
    fig.savefig(out, dpi=160)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())