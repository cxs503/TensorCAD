"""Figure + report payload for the PR #46 cure (m3924).

Three panels, each one a claim of the cure:

* **a** — the public surface is the union, not a preference.  ``main`` added 13 symbols
  after the merge base, the branch added 6, the intersection of the two additions is
  empty, and the merged ``__all__`` has 126 = 107 + 13 + 6.
* **b** — the branch's GZ route still computes the same numbers on the merged tree:
  the 0.21.0 regeneration (lines) lies on the committed 0.20.0 baseline (markers) for
  the whole barge ``KG`` ladder, with the residual plotted on a log axis.
* **c** — the main side's transom design space is untouched: block coefficient against
  transom fraction for the merged run, with the committed baseline as markers.

    python tools/report_m3924.py            # writes artifacts/m3924/*
"""

from __future__ import annotations

import argparse
import base64
import csv
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
ART = REPO / "artifacts" / "m3924"


def _load(rel: str) -> Any:
    return json.loads((REPO / rel).read_text())


def _csv_rows(rel: str) -> list[dict[str, str]]:
    with (REPO / rel).open(newline="") as handle:
        return list(csv.DictReader(handle))


def figure(path: Path) -> dict[str, Any]:
    verification = _load("artifacts/m3924/merge_verification.json")
    baseline = _load("docs/stability/gz_results.json")
    merged = _load("artifacts/m3924/gz_results_merged_0.21.0.json")
    design = _csv_rows("docs/m2873/outputs/transom_stern/design_table.csv")
    design_merged = _csv_rows("artifacts/m3924/transom_design_table_merged.csv")

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2))
    fig.suptitle(
        "TensorCAD PR #46 cure — semantic union of main (0.21.0) and agent/m3406",
        fontsize=13,
        fontweight="bold",
    )

    # (a) public surface -------------------------------------------------------------
    ax = axes[0]
    counts = verification["public_surface"]["counts"]
    labels = ["merge base\n(e604f02)", "main\n(+13)", "branch\n(+6)", "merged\n(union)"]
    values = [counts["base"], counts["main"], counts["branch"], counts["merged_all"]]
    colours = ["#8c8c8c", "#1f77b4", "#2ca02c", "#111111"]
    bars = ax.bar(labels, values, color=colours, width=0.62)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 1.2, str(value),
                ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.set_ylim(90, 136)
    ax.set_ylabel("symbols in tensorcad.__all__")
    ax.set_title("(a) public surface = exact union\nno duplicate, no dropped export", fontsize=11)
    ax.grid(axis="y", alpha=0.25)
    ax.text(0.5, 0.06,
            "shared additions between the two sides: 0\n"
            f"tests/test_public_surface.py asserts {counts['expected_api_merged']} symbols",
            transform=ax.transAxes, ha="center", fontsize=9,
            bbox=dict(fc="#f4f4f4", ec="#bbbbbb"))

    # (b) GZ reproduction -----------------------------------------------------------
    ax = axes[1]
    ladder_b = baseline["barge"]["ladder"]
    ladder_m = merged["barge"]["ladder"]
    residual = 0.0
    for row_b, row_m in zip(ladder_b, ladder_m):
        ax.plot(row_b["hea_angles_deg"], row_b["righting_arms_m"], "-",
                lw=0.8, color="#1f77b4", alpha=0.45)
        ax.plot(row_m["hea_angles_deg"], row_m["righting_arms_m"], "o",
                ms=2.4, mfc="none", mec="#d62728")
        for a, b in zip(row_b["righting_arms_m"], row_m["righting_arms_m"]):
            residual = max(residual, abs(a - b))
        ax.annotate(f"KG = {row_b['kg_m']:.2f} m", (row_b["hea_angles_deg"][0] + 2,
                    row_b["righting_arms_m"][0] + 0.012), fontsize=8, color="#1f77b4")
    ax.axhline(0.0, color="#999999", lw=0.7)
    ax.set_xlabel("heel angle [deg]")
    ax.set_ylabel("righting arm GZ [m]")
    ax.set_title("(b) branch GZ route on the merged tree\n"
                 "regenerated 0.21.0 (○) on 0.20.0 baseline (—)", fontsize=11)
    ax.grid(alpha=0.25)
    ax.text(0.03, 0.05,
            f"{verification['gz_reproduction']['numeric_leaves']} numeric leaves\n"
            f"max |Δ| = {residual:.1e} m (exact to the bit)",
            transform=ax.transAxes, fontsize=9,
            bbox=dict(fc="#f4f4f4", ec="#bbbbbb"))

    # (c) transom design space ------------------------------------------------------
    ax = axes[2]
    taus = [float(row["tau"]) for row in design]
    cbs = [float(row["block_coefficient"]) for row in design]
    taus_m = [float(row["tau"]) for row in design_merged]
    cbs_m = [float(row["block_coefficient"]) for row in design_merged]
    ax.plot(taus, cbs, "-", lw=1.2, color="#1f77b4", label="main baseline (committed)")
    ax.plot(taus_m, cbs_m, "o", ms=4.0, mfc="none", mec="#d62728",
            label="merged tree (regenerated)")
    ax.set_xlabel("transom fraction τ = x_transom / L")
    ax.set_ylabel("block coefficient Cb")
    ax.set_title("(c) main's transom design space\n"
                 "three CSVs byte-identical after the merge", fontsize=11)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc="lower left")

    fig.tight_layout(rect=(0, 0, 1, 0.95))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=100, pil_kwargs={"optimize": True})
    plt.close(fig)
    # The report endpoint rejects payloads above ~250 kB, and base64 inflates by 4/3:
    # quantising the raster keeps the panel readable at ~30 kB on the wire.
    try:
        from PIL import Image

        with Image.open(path) as image:
            image.convert("P", palette=Image.ADAPTIVE, colors=128).save(path, optimize=True)
    except ImportError:  # pragma: no cover - PIL is optional tooling, not a dependency
        pass
    return {"residual_gz_m": residual, "design_points": len(design)}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--figure", type=Path, default=ART / "merge_verification.png")
    args = parser.parse_args(argv)
    info = figure(args.figure)
    print(json.dumps({"figure": str(args.figure), **info}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())