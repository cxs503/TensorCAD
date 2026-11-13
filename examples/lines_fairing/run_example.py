"""Regenerate docs/lines_fairing/{results.json,report.md,plot.png} for GAP-CAD-C3-LINES-FAIRING.

Every number in the page comes from this run; nothing is transcribed by hand.
"""
from __future__ import annotations

import base64
import json
import math
import random
from pathlib import Path

from tensorcad.lines_fairing import (
    FairingSettings,
    curvature_sign_flips,
    fair_offset_table,
    hand_fairing_reference,
    line_bending_energy,
)
from tensorcad.offset_forms import offset_table
from tensorcad.offset_loft import ShipOffsetTable, table_integrals

OUT = Path("docs/lines_fairing")
GRIDS = {"wigley_iii": (21, 17), "series60_cb060": (25, 21), "kvlcc2_class": (33, 25)}
NOISE = {"wigley_iii": 4.0e-3, "series60_cb060": 4.0e-3, "kvlcc2_class": 2.0e-2}
TOLERANCE = {"wigley_iii": 1.0e-2, "series60_cb060": 1.0e-2, "kvlcc2_class": 5.0e-2}


def scattered(table: ShipOffsetTable, amplitude: float, seed: int = 7) -> ShipOffsetTable:
    generator = random.Random(seed)
    rows = [
        tuple(
            value + generator.uniform(-amplitude, amplitude) if value > 1.0e-9 else value
            for value in row
        )
        for row in table.half_breadths
    ]
    return ShipOffsetTable(
        name=f"{table.name}-scattered", length=table.length, stations=table.stations,
        waterlines=table.waterlines, half_breadths=tuple(rows),
        design_draft=table.design_draft, provenance=table.provenance, units=table.units,
        waterline_exponent=table.waterline_exponent,
    )


def distances(table: ShipOffsetTable, reference: ShipOffsetTable) -> tuple[float, float]:
    errors = [
        table.half_breadths[i][j] - reference.half_breadths[i][j]
        for i in range(len(reference.stations))
        for j in range(len(reference.waterlines))
    ]
    return (
        math.sqrt(math.fsum(e * e for e in errors) / len(errors)),
        max(abs(e) for e in errors),
    )


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    study: dict[str, object] = {"schema": "tensorcad.lines_fairing.study/1", "hulls": {}}
    curves_for_plot: list[tuple[str, list[float], list[float], float]] = []
    claims_ok = True

    for key, (nx, nz) in GRIDS.items():
        clean = offset_table(key, n_stations=nx, n_waterlines=nz)
        rough = scattered(clean, NOISE[key])
        entry: dict[str, object] = {"grid": [nx, nz], "noise": NOISE[key],
                                    "tolerance": TOLERANCE[key], "kernels": {}}
        for kernel in ("local", "energy"):
            result = fair_offset_table(
                rough,
                FairingSettings(offset_tolerance=TOLERANCE[key], max_passes=5,
                                buttock_levels=5, kernel=kernel),
            )
            report = result.to_dict()
            before = distances(rough, clean)
            after = distances(result.table, clean)
            reference, touched = hand_fairing_reference(rough, window=5)
            diffs = [abs(reference[i][j] - result.table.half_breadths[i][j])
                     for i in range(len(clean.stations))
                     for j in range(len(clean.waterlines)) if touched[i][j] > 0]
            history = list(report["total_flips"])
            monotone = all(b <= a for a, b in zip(history, history[1:]))
            within = report["tolerance_respected"]
            claims_ok &= monotone and within
            entry["kernels"][kernel] = {
                "flips": history,
                "families": report["flips"],
                "monotone": monotone,
                "max_deviation": report["max_deviation"],
                "tolerance_respected": within,
                "rms_to_analytic_in": before[0], "rms_to_analytic_out": after[0],
                "max_to_analytic_in": before[1], "max_to_analytic_out": after[1],
                "rms_error_reduction": 1.0 - after[0] / before[0],
                "hand_reference_max": max(diffs) if diffs else None,
                "hand_reference_rms": (
                    math.sqrt(math.fsum(d * d for d in diffs) / len(diffs)) if diffs else None
                ),
                "volume_relative_change": abs(report["volume_after"] - report["volume_before"])
                / report["volume_before"],
                "passes": [{"index": p["index"], "accepted": p["accepted"],
                            "flips": p["total_flips"],
                            "relaxation": p["relaxation"]} for p in report["passes"]],
            }
            if kernel == "local":
                j = len(clean.waterlines) // 2
                curves_for_plot.append(
                    (
                        f"{key} waterline z={clean.waterlines[j]:.2f} m",
                        list(rough.half_breadths[i][j] for i in range(len(clean.stations))),
                        list(result.table.half_breadths[i][j] for i in range(len(clean.stations))),
                        max(row[j] for row in clean.half_breadths),
                    )
                )
        study["hulls"][key] = entry

    (OUT / "results.json").write_text(json.dumps(study, indent=2, sort_keys=True) + "\n")

    figure_reason = ""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as error:  # pragma: no cover - environment dependent
        plt = None
        figure_reason = f"matplotlib unavailable: {error}"

    rows = []
    if plt is not None:
        fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.4))
        for axis, (label, before, after, peak) in zip(axes[:3], curves_for_plot):
            axis.plot(range(len(before)), before, "o", ms=3, alpha=0.55, label="scattered input")
            axis.plot(range(len(after)), after, "-", lw=1.6, label="faired")
            axis.set_title(label, fontsize=9)
            axis.set_xlabel("station index")
            axis.set_ylabel("half-breadth y [m]")
            axis.legend(fontsize=7)
        fig.suptitle("Three-directional fairing: one waterline per hull (local kernel)", fontsize=10)
        fig.tight_layout()
        fig.savefig(OUT / "plot.png", dpi=130)
        plt.close(fig)
        rows.append("![fairing](plot.png)")
    else:
        rows.append(f"(no figure: {figure_reason})")

    lines = [
        "# Three-directional offset-table fairing — measured results",
        "",
        "Regenerated by `python examples/lines_fairing/run_example.py`; every number below",
        "is read from `results.json`, which the runner writes. The runner reports the",
        f"acceptance status (`max_deviation <= tolerance`, monotone flips): **{'PASS' if claims_ok else 'FAIL'}**.",
        "",
        "| hull | grid | noise [m] | tol [m] | kernel | flips (history) | monotone | max dev [m] |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for key, entry in study["hulls"].items():
        for kernel, data in entry["kernels"].items():
            lines.append(
                f"| {key} | {entry['grid'][0]}x{entry['grid'][1]} | {entry['noise']:.3f} |"
                f" {entry['tolerance']:.3f} | {kernel} | {data['flips']} |"
                f" {'yes' if data['monotone'] else 'no'} | {data['max_deviation']:.5f} |"
            )
    lines += [
        "",
        "## Accuracy against the analytic half-breadth surface",
        "",
        "| hull | kernel | RMS in [m] | RMS out [m] | reduction | max in [m] | max out [m] |",
        "|---|---|---|---|---|---|---|",
    ]
    for key, entry in study["hulls"].items():
        for kernel, data in entry["kernels"].items():
            lines.append(
                f"| {key} | {kernel} | {data['rms_to_analytic_in']:.5f} |"
                f" {data['rms_to_analytic_out']:.5f} | {100 * data['rms_error_reduction']:.1f} % |"
                f" {data['max_to_analytic_in']:.5f} | {data['max_to_analytic_out']:.5f} |"
            )
    lines += ["", "## Independent smoother and displacement", "",
              "| hull | kernel | RMS vs hand reference [m] | max vs hand reference [m] | |dV|/V |",
              "|---|---|---|---|---|"]
    for key, entry in study["hulls"].items():
        for kernel, data in entry["kernels"].items():
            hand_rms = data["hand_reference_rms"]
            hand_max = data["hand_reference_max"]
            lines.append(
                f"| {key} | {kernel} |"
                f" {hand_rms:.5f} | {hand_max:.5f} | {data['volume_relative_change']:.2e} |"
            )
    lines += ["", "## Figures", "", *rows, ""]
    (OUT / "report.md").write_text("\n".join(lines) + "\n")
    print(json.dumps({"acceptance": "PASS" if claims_ok else "FAIL"}, indent=2))
    return 0 if claims_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
