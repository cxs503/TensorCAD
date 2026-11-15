"""Figure for the DTMB 5415 published-particulars cross-check (docs/golden_regression.md).

Reads the *committed* baseline and the package's own reference derivation, so the numbers
in the plot cannot be typed by hand: the same call the gate makes produces both bars.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
BASELINE = ROOT / "docs" / "golden" / "tensorcad_geometry_golden.json"

METRICS = {
    "displacement_volume_m3": "displacement\nvolume",
    "block_coefficient": "block\ncoefficient Cb",
    "length_stations_m": "station span\n(vs Lpp)",
    "waterline_beam_m": "beam at WL\n(vs B moulded)",
}


def main() -> None:
    from tensorcad.dtmb5415 import README_PARTICULARS_TOLERANCE

    data = json.loads(BASELINE.read_text(encoding="utf-8"))["cases"]["dtmb5415"]
    fine = data["measurements"]["fine"]
    reference = data["reference"]
    names = [METRICS[key] for key in METRICS]
    errors = [100.0 * (fine[key] - reference[key]) / reference[key] for key in METRICS]
    band = 100.0 * README_PARTICULARS_TOLERANCE

    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2), width_ratios=[1.15, 1])
    axis = axes[0]
    bars = axis.bar(names, errors, color=["#2b6cb0" if abs(e) <= band else "#c53030" for e in errors])
    axis.axhline(0.0, color="#2d3748", linewidth=0.8)
    for sign in (+1, -1):
        axis.axhline(sign * band, color="#c53030", linestyle="--", linewidth=1.0)
    axis.text(len(names) - 0.5, band, f"contract band ±{band:.0f}%", color="#c53030",
              ha="right", va="bottom", fontsize=9)
    for bar, error in zip(bars, errors):
        axis.annotate(f"{error:+.2f}%", (bar.get_x() + bar.get_width() / 2, error),
                      textcoords="offset points", xytext=(0, 4 if error >= 0 else -12),
                      ha="center", fontsize=9)
    axis.set_ylabel("measured vs published [%]")
    axis.set_title("DTMB 5415: asserted against the dataset README table")
    axis.set_ylim(-band * 1.6, band * 1.6)
    axis.grid(axis="y", alpha=0.3)

    axis = axes[1]
    keys = list(METRICS)
    axis.scatter([1.0] * len(keys), [reference[k] for k in keys], marker="o", s=52,
                 facecolors="none", edgecolors="#2b6cb0", label="published (README x λ)")
    axis.scatter([2.0] * len(keys), [fine[k] for k in keys], marker="x", s=46,
                 color="#2d3748", label="shipped reduced model")
    for key in keys:
        axis.plot([1.0, 2.0], [reference[key], fine[key]], color="#a0aec0", linewidth=0.8)
    axis.set_xticks([1.0, 2.0], ["published", "measured"])
    axis.set_yscale("log")
    axis.set_title("same four quantities, log scale")
    axis.legend(fontsize=8, loc="lower right")
    axis.grid(alpha=0.3)

    figure.tight_layout()
    output = ROOT / "docs" / "golden" / "dtmb5415_particulars.png"
    figure.savefig(output, dpi=140)
    print(output, output.stat().st_size)


if __name__ == "__main__":
    main()
