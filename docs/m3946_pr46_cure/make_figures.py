"""Figures for the PR #46 merge cure (mission 3946).

Panel A: wall-sided box GZ -- integration vs closed form, with the deck-corner
         immersion angle marked so the divergence past the formula's validity
         is visible rather than hidden.
Panel B: circular-section prism convergence -- relative error per facet count,
         against the measured second-order slope.
Panel C: merge reconciliation -- per-file side that changed, and the pytest
         counts of main vs the merged tree (parsed from the run logs).

Every number is read from docs/stability/gz_results.json and the pytest logs;
nothing is typed in by hand.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = json.loads((ROOT / "docs" / "stability" / "gz_results.json").read_text())

BASE_LOG = Path("/tmp/pytest_base.log")
MERGED_LOG = Path("/tmp/pytest_merged.log")


def _pytest_counts(path: Path) -> dict[str, int]:
    """Parse the trailing 'N passed, M skipped ...' line of a pytest -q log."""
    if not path.exists():
        return {}
    text = path.read_text(errors="replace")
    tail = text.strip().splitlines()[-1] if text.strip() else ""
    out: dict[str, int] = {}
    for num, word in re.findall(r"(\d+) (passed|failed|error|errors|skipped|xfailed|xpassed)", tail):
        key = "error" if word == "errors" else word
        out[key] = out.get(key, 0) + int(num)
    return out


def panel_a(ax: plt.Axes) -> None:
    ws = DATA["wall_sided"]
    heel = np.array([r["heel_angle_deg"] for r in ws["rows"]])
    integ = np.array([r["righting_arm_m"] for r in ws["rows"]])
    closed = np.array([r["closed_form_m"] for r in ws["rows"]])
    limit = ws["deck_corner_immersion_deg"]
    ax.axvspan(limit, 72, color="0.85", alpha=0.6, lw=0, label="closed form invalid (deck corner immersed)")
    ax.axvline(limit, color="0.35", ls="--", lw=1.0)
    ax.plot(heel, integ, "o", ms=4.5, color="#1f4e79", label="heeled-displacement integration")
    ax.plot(np.r_[heel, [65.0, 70.0]], np.r_[closed, [ws["beyond_the_limit"][k]["closed_form_m"] for k in ("65.0", "70.0")]],
            "-", lw=1.6, color="#c00000", label="wall-sided closed form")
    ax.plot([65.0, 70.0], [ws["beyond_the_limit"][k]["righting_arm_m"] for k in ("65.0", "70.0")],
            "s", ms=4.5, mfc="none", color="#1f4e79")
    ax.annotate(f"max |diff| = {ws['max_abs_difference_m']:.1e} m\nfor heel $\\leq$ {heel.max():.0f}$^\\circ$",
                xy=(30, 1.35), fontsize=8, color="#1f4e79")
    ax.annotate(f"deck corner at {limit:.2f}$^\\circ$", xy=(limit, 0.25), xytext=(40, 0.35), fontsize=8,
                arrowprops=dict(arrowstyle="->", color="0.35", lw=0.8))
    ax.set_xlabel("heel angle [deg]")
    ax.set_ylabel("righting arm GZ [m]")
    ax.set_title("(a) wall-sided box: L=%.0f, B=%.0f, T=%.0f, KG=%.1f m" % (
        ws["geometry"]["length"], ws["geometry"]["beam"], ws["geometry"]["draft"], ws["geometry"]["kg"]), fontsize=9)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="upper left")


def panel_b(ax: plt.Axes) -> None:
    circ = DATA["circular"]
    rows = circ["rows"]
    sides = np.array([r["sides"] for r in rows], dtype=float)
    err = np.array([r["relative_error"] for r in rows])
    ax.loglog(sides, err, "o-", color="#1f4e79", ms=5, label="measured max relative error")
    ref = err[0] * (sides[0] / sides) ** 2
    ax.loglog(sides, ref, "--", color="#c00000", lw=1.2, label="second order reference $\\propto N^{-2}$")
    for n, e in zip(sides, err):
        ax.annotate(f"{e:.2e}", xy=(n, e), xytext=(4, -12), textcoords="offset points", fontsize=7)
    ax.set_xlabel("facets around the circular section")
    ax.set_ylabel("max relative error of $GZ/\\sin\\phi$")
    ax.set_title("(b) inscribed prism, R=%.0f m, exact $GZ/\\sin\\phi$=%.1f m; "
                 "error ratios per doubling: %s" % (
                     circ.get("radius_m", rows[0]["circumradius_m"]), circ["expected_gz_over_sine_m"],
                     ", ".join(f"{r:.2f}" for r in circ["error_ratios_per_doubling"])), fontsize=9)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)


def panel_c(ax: plt.Axes) -> None:
    labels = ["main\n(2ff80b2)", "PR #46 merged"]
    base = _pytest_counts(BASE_LOG)
    merged = _pytest_counts(MERGED_LOG)
    keys = [k for k in ("passed", "skipped", "failed", "error") if k in base or k in merged]
    if not keys:
        ax.text(0.5, 0.5, "pytest logs unavailable", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return
    x = np.arange(len(keys))
    w = 0.36
    bvals = [base.get(k, 0) for k in keys]
    mvals = [merged.get(k, 0) for k in keys]
    ax.bar(x - w / 2, bvals, w, color="#8faadc", label=labels[0])
    ax.bar(x + w / 2, mvals, w, color="#1f4e79", label=labels[1])
    for xi, (b, m) in enumerate(zip(bvals, mvals)):
        ax.annotate(str(b), (xi - w / 2, b), ha="center", va="bottom", fontsize=7)
        ax.annotate(str(m), (xi + w / 2, m), ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x, keys)
    ax.set_ylabel("tests (pytest -q, neural modules excluded)")
    ax.set_title("(c) core suite: merged tree vs main baseline", fontsize=9)
    ax.grid(axis="y", alpha=0.3)
    ax.legend(fontsize=7)


def main() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.4))
    panel_a(axes[0])
    panel_b(axes[1])
    panel_c(axes[2])
    fig.suptitle("TensorCAD PR #46 cure: righting-arm validation carried over the merge + suite parity", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = HERE / "pr46_merge_evidence.png"
    fig.savefig(out, dpi=135)
    print("wrote", out)


if __name__ == "__main__":
    main()