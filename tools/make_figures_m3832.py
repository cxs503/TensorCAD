"""Figures for the m3832 cure report (PR #40): post-merge coexistence of both feature lines.

Every number comes out of a file that the code itself produces, never hand-typed:

* panel (a) - ``docs/righting_arm/results.json`` ``curve_payload`` (this branch's feature,
  ``tensorcad/righting_arm.py``);
* panel (b) - the same file's ``barge`` block: mesh GZ against the wall-sided closed form;
* panel (c) - ``docs/m2873/outputs/transom_stern/convergence.csv`` (main's feature,
  ``tensorcad/transom.py``) volume error against the analytic reference;
* panel (d) - ``artifacts/m3832/pytest_summary.json``, written by the two pytest runs.

Usage: ``python3 tools/make_figures_m3832.py``
"""

from __future__ import annotations

import csv
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "m3832" / "post_merge_verification.png"


def panel_a(ax) -> None:
    payload = json.loads((ROOT / "docs/righting_arm/results.json").read_text())["curve_payload"]
    heels = payload["heels_deg"]
    arms = payload["righting_arms"]
    peak = payload["peak"]
    ax.plot(heels, arms, "o-", color="#1f4e79", ms=4, label="mesh GZ (righting_arm)")
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.plot([peak["heel_deg"]], [peak["righting_arm"]], "r*", ms=11,
            label=f"peak GZ = {peak['righting_arm']:.3f} m @ {peak['heel_deg']:.0f} deg")
    ax.axvline(payload.get("vanishing_crossed") or 0.0, color="0.8", lw=0.0)
    ax.set_xlabel("heel phi [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title(f"(a) this branch: GZ(phi) at constant displacement\n"
                 f"nabla = {payload['displacement_volume']:.3f} m^3, KG = {payload['kg']:.2f} m")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left")


def panel_b(ax) -> None:
    barge = json.loads((ROOT / "docs/righting_arm/results.json").read_text())["barge"]
    inside = barge["inside_window"]
    outside = barge["beyond_window"]
    ax.semilogy([p["heel_deg"] for p in inside if p["heel_deg"] > 0],
                [abs(p["difference_m"]) for p in inside if p["heel_deg"] > 0],
                "o", ms=4, color="#2e7d32", label="inside validity window")
    ax.semilogy([p["heel_deg"] for p in outside if p["heel_deg"] > 0],
                [abs(p["difference_m"]) for p in outside if p["heel_deg"] > 0],
                "s", ms=4, color="#c62828", label="past deck-edge immersion")
    ax.axvline(barge["wall_sided_validity_deg"], color="0.4", ls="--", lw=0.9,
               label=f"deck edge @ {barge['wall_sided_validity_deg']:.2f} deg")
    ax.set_xlabel("heel phi [deg]")
    ax.set_ylabel("|GZ_mesh - GZ_wall-sided| [m]")
    ax.set_title("(b) exact reference, not a fit:\nwall-sided barge, error vs heel")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7, loc="lower right")


def panel_c(ax) -> None:
    rows = list(csv.DictReader((ROOT / "docs/m2873/outputs/transom_stern/convergence.csv").open()))
    by_tau: dict[float, list[tuple[int, float]]] = {}
    for row in rows:
        tau = float(row["tau"])
        by_tau.setdefault(tau, []).append((int(row["nz"]), abs(float(row["volume_error"])) * 100.0))
    for tau, points in sorted(by_tau.items()):
        points.sort()
        ax.loglog([p[0] for p in points], [max(p[1], 1e-4) for p in points], "o-", ms=3.5,
                  label=f"tau = {tau:g}")
    ax.set_xlabel("transverse grid nz")
    ax.set_ylabel("|volume error| vs analytic [%]")
    ax.set_title("(c) main's feature still converges:\ntransom-stern volume error")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)


def panel_d(ax) -> None:
    path = ROOT / "artifacts/m3832/pytest_summary.json"
    data = json.loads(path.read_text())
    labels = ["main\n(2ff80b2)", "merged\n(agent/m3169)"]
    rows = [data["baseline"], data["merged"]]
    passed = [r["passed"] for r in rows]
    x = range(len(rows))
    ax.bar(x, passed, width=0.55, color=["#78909c", "#1f4e79"])
    for xi, row in zip(x, rows):
        ax.text(xi, row["passed"] + 12,
                f"{row['passed']} passed\n{row['skipped']} skipped\n{row['failed']} failed",
                ha="center", fontsize=8)
    ax.set_xticks(list(x), labels)
    ax.set_ylim(0, max(passed) * 1.22)
    ax.set_ylabel("tests")
    ax.set_title("(d) no regression:\ntest suite, baseline vs merge")
    ax.grid(alpha=0.3, axis="y")


def main() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.2), dpi=140)
    panel_a(axes[0][0])
    panel_b(axes[0][1])
    panel_c(axes[1][0])
    panel_d(axes[1][1])
    fig.suptitle("TensorCAD PR #40 cure (m3832): both feature lines alive after the merge",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight")
    print("wrote", OUT, OUT.stat().st_size, "bytes")


if __name__ == "__main__":
    main()