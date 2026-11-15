"""Figures for the mesh-quality door report (``docs/mesh_quality_door.md``).

Reads ``artifacts/mesh_quality_door/results.json`` -- the record the door itself wrote --
so no number on the figure is transcribed by hand. Writes
``artifacts/mesh_quality_door/door.png``.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "artifacts" / "mesh_quality_door" / "results.json"
OUTPUT = ROOT / "artifacts" / "mesh_quality_door" / "door.png"


def main() -> int:
    record = json.loads(RESULTS.read_text(encoding="utf-8"))
    figure, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))

    # (a) closed form vs library ------------------------------------------------------
    rows = record["closed_form"]
    labels = [f"{row['case']}\n{row['quantity']}" for row in rows]
    deviation = [max(row["difference"] / max(abs(row["expected"]), 1e-300), 1e-18) for row in rows]
    axes[0].barh(range(len(rows)), deviation, color="#2b6cb0")
    axes[0].set_yticks(range(len(rows)))
    axes[0].set_yticklabels(labels, fontsize=7)
    axes[0].set_xscale("log")
    axes[0].axvline(record.get("tolerance", 1e-12), color="#c53030", linestyle="--", linewidth=1.2,
                    label=f"door tolerance {record.get('tolerance', 1e-12):g}")
    axes[0].set_xlabel("relative deviation, library vs closed form")
    axes[0].set_title("(a) the volume cases against their closed forms")
    axes[0].legend(fontsize=7)
    axes[0].invert_yaxis()

    # (b) the shell cases against the shipped angle floor ------------------------------
    shells = ["wigley", "series60", "suboff_coarse", "suboff_taubin30", "box_unit"]
    values = [record["cases"][name]["checks"].get("tri.min_angle", math.nan) for name in shells]
    colors = ["#2f855a" if record["cases"][name]["passed"] else "#c53030" for name in shells]
    axes[1].bar(range(len(shells)), values, color=colors)
    axes[1].set_xticks(range(len(shells)))
    axes[1].set_xticklabels(shells, rotation=20, ha="right", fontsize=8)
    axes[1].set_yscale("log")
    axes[1].axhline(10.0, color="#c53030", linestyle="--", linewidth=1.2, label="shipped floor 10 deg")
    axes[1].set_ylabel("min angle (deg, log)")
    axes[1].set_title("(b) shell cases: the angle floor is what decides")
    axes[1].legend(fontsize=7)
    for position, value in enumerate(values):
        axes[1].annotate(f"{value:.4g}", (position, value), textcoords="offset points",
                         xytext=(0, 3), ha="center", fontsize=7)

    # (c) the sheared block: measured vs closed form ----------------------------------
    block = record["cases"]["hex_shear_2"]["checks"]
    quantities = ["volume (each)", "min_sicn", "scaled_jacobian", "non-orthogonality (deg)"]
    measured = [
        record["cases"]["hex_shear_2"]["snapshot"]["measurements"]["summary.hex.summary.volume.mean"],
        block["hex.min_sicn"],
        block["hex.scaled_jacobian"],
        block["volume.max_non_orthogonality"],
    ]
    expected = [
        next(row["expected"] for row in rows if row["quantity"] == "volume_per_cell"),
        next(row["expected"] for row in rows if row["quantity"] == "min_sicn"),
        next(row["expected"] for row in rows if row["quantity"] == "scaled_jacobian"),
        next(row["expected"] for row in rows if row["quantity"] == "max_non_orthogonality"),
    ]
    positions = range(len(quantities))
    axes[2].bar([p - 0.2 for p in positions], expected, width=0.4, label="closed form", color="#4a5568")
    axes[2].bar([p + 0.2 for p in positions], measured, width=0.4, label="library", color="#dd6b20")
    axes[2].set_xticks(list(positions))
    axes[2].set_xticklabels(quantities, rotation=20, ha="right", fontsize=8)
    axes[2].set_yscale("log")
    axes[2].set_ylabel("value (log)")
    axes[2].set_title("(c) sheared block: the four exact quantities")
    axes[2].legend(fontsize=7)

    figure.suptitle(
        "mesh-quality door: closed-form differential and the meshes the verdict is measured on "
        f"({'PASS' if record['passed'] else 'FAIL'})",
        fontsize=10,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=150)
    print(OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())