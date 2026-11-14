#!/usr/bin/env python3
"""Round-2 evidence for GAP-CAD-T2-PARAM-SERIES: the family through its entry.

Everything here is driven through the *product entry* (``tensorcad-hull series60``)
rather than through an import, because the recurring audit finding on this anchor
was that a kernel nobody can reach is not a capability.  The script

1. runs ``series60 family`` over the five Todd parents on a three-grid ladder
   (81x21 -> 161x41 -> 321x81) and compares the generated coefficients with the
   published particulars, retyping nothing: the published table is read from
   ``docs/data/series60_published_parents.json`` and its sha256 is recorded;
2. runs ``series60 parent`` at two points of a *second* parameter axis (transom
   breadth ratio 0.0 / 0.15 at Cb = 0.70) so the sweep is not a one-parameter
   study;
3. emits ``round2_evidence.json`` + CSVs + figures, and a verdict per acceptance
   criterion with the number that decides it.

The role labels are deliberately blunt.  ``Cwp`` and ``LCB/L`` are handed to the
generator as targets and ``Cp`` is the published ladder the exponents are read
from, so agreement on those is a **closure**, not an external verification.  The
only externally anchored quantity here is ``Cb`` itself plus the published
``Cm``/``Cp`` *values* being reproduced through the geometry; the family's
*shape* against a Todd offset table remains WITHHELD (no such table in-repo, no
network).  Saying so is the point: a closure reported as a verification is how
the previous rounds produced "passing" numbers nobody could trust.

Usage::

    python3 tools/series60_family/round2_evidence.py [--out runs/m19694_series60_family]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "docs" / "data" / "series60_published_parents.json"

GRID_LADDER = ((81, 21), (161, 41), (321, 81))
TRANSOM_POINTS = (0.0, 0.15)
L, B, T = 10.0, 1.5, 0.6

#: Predeclared acceptance bars (unchanged from the earlier rounds: 3 % against the
#: published particulars, 2 % mesh-vs-closed-form).
PUBLISHED_TOLERANCE = 0.03
DISCRETISATION_TOLERANCE = 0.02


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_cli(args: list[str], out: Path) -> dict:
    """Invoke the product entry as a user would; return the record it wrote."""
    out.mkdir(parents=True, exist_ok=True)
    command = [sys.executable, "-m", "tensorcad.cli", "series60", *args, "--out", str(out)]
    completed = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    if completed.returncode not in (0, 1):  # 1 = a claim failed; the record still exists
        raise SystemExit(f"entry failed ({completed.returncode}):\n{completed.stderr}")
    return {"command": " ".join(command[2:]), "returncode": completed.returncode,
            "stderr": completed.stderr.strip()}


def family_record(grid: tuple[int, int], out: Path) -> dict:
    nx, nz = grid
    run = run_cli(["family", "--published", str(DATA), "--nx", str(nx), "--nz", str(nz)], out)
    document = json.loads((out / "five_cb_gradient.json").read_text())
    document["entry_run"] = run
    return document


def parent_record(*extra: str, out: Path) -> dict:
    """One parent: the exported document is the artifact, the record is the audit trail."""
    run = run_cli(["parent", *extra, "--output", str(out / "hull")], out)
    document = json.loads((out / "series60_parent.json").read_text())
    document["entry_run"] = run
    return document


def worst_published(document: dict) -> tuple[float, str]:
    worst, where = 0.0, ""
    for row in document["rows"]:
        for source, entries in row["error_vs_published"].items():
            for name, entry in entries.items():
                value = entry.get("relative")
                if value is None:
                    continue
                if abs(value) > abs(worst):
                    worst, where = value, f"Cb={row['cb_requested']:.2f} {source} {name}"
    return worst, where


def worst_discretisation(document: dict) -> float:
    worst = 0.0
    for row in document["rows"]:
        for value in row["discretisation_relative_errors"].values():
            if value is None:
                continue
            worst = max(worst, abs(value))
    return worst


def figure_gradient(ladder: dict[tuple[int, int], dict], target: Path) -> None:
    grid = GRID_LADDER[-1]
    rows = ladder[grid]["rows"]
    cbs = [row["cb_requested"] for row in rows]
    series = {
        "Cwp": ("waterplane_coefficient", "input target (closure)"),
        "LCB/L": ("lcb_over_length", "input target (closure)"),
        "Cm": ("midship_coefficient", "published ladder (input)"),
        "Cp": ("prismatic_coefficient", "published ladder (input)"),
    }
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    width = 0.2
    for index, (label, (key, _)) in enumerate(series.items()):
        got = [row["measured"][key] for row in rows]
        axes[0].bar([c + (index - 1.5) * width for c in range(len(cbs))], got, width,
                    label=f"{label} (generated)")
        expect = [row["error_vs_published"]["measured"][key]["published"] for row in rows]
        axes[0].plot(range(len(cbs)), expect, "k_", markersize=18,
                     label="published (Todd 1951)" if index == 0 else None)
    axes[0].set_xticks(range(len(cbs)))
    axes[0].set_xticklabels([f"{c:.2f}" for c in cbs])
    axes[0].set_xlabel("Cb (parent form)")
    axes[0].set_ylabel("coefficient")
    axes[0].set_title("Five parents: generated vs published particulars")
    axes[0].legend(fontsize=7, loc="center left", bbox_to_anchor=(0.0, 0.62))

    for label, (key, role) in series.items():
        errors = [100.0 * row["error_vs_published"]["measured"][key]["relative"] for row in rows]
        axes[1].plot(cbs, errors, "o-", label=f"{label}: {role}")
    axes[1].axhline(PUBLISHED_TOLERANCE * 100.0, color="red", ls=":", lw=0.8)
    axes[1].axhline(-PUBLISHED_TOLERANCE * 100.0, color="red", ls=":", lw=0.8,
                    label="±3 % bar")
    axes[1].axhline(0.0, color="gray", lw=0.5)
    axes[1].set_xlabel("Cb (parent form)")
    axes[1].set_ylabel("measured vs published (%)")
    axes[1].set_title("Residuals: two of the four are closures by construction")
    axes[1].legend(fontsize=6)
    axes[1].grid(alpha=0.3)
    fig.suptitle("Series 60 parents through `tensorcad-hull series60 family`", fontsize=13)
    fig.tight_layout()
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def figure_convergence(ladder: dict[tuple[int, int], dict], target: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.4))
    colours = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
    for index, cb in enumerate((0.60, 0.65, 0.70, 0.75, 0.80)):
        volumes, midships = [], []
        for grid in GRID_LADDER:
            row = next(r for r in ladder[grid]["rows"] if abs(r["cb_requested"] - cb) < 1e-9)
            volumes.append(100.0 * abs(row["discretisation_relative_errors"]["displacement_volume"]))
            midships.append(100.0 * abs(row["discretisation_relative_errors"]["midship_coefficient"]))
        cells = [nx * nz for nx, nz in GRID_LADDER]
        axes[0].loglog(cells, volumes, "o-", color=colours[index], label=f"Cb={cb:.2f}")
        axes[1].loglog(cells, midships, "o-", color=colours[index], label=f"Cb={cb:.2f}")
    for axis, title in zip(axes, ("displacement volume", "midship section area")):
        axis.axhline(DISCRETISATION_TOLERANCE * 100.0, color="red", ls=":", lw=0.8,
                     label="2 % bar")
        axis.set_xlabel("surface cells nx*nz")
        axis.set_ylabel(f"|mesh - closed form| on {title} (%)")
        axis.set_title(f"Discretisation: {title}")
        axis.grid(alpha=0.3, which="both")
        axis.legend(fontsize=7)
    fig.suptitle("Series 60 convergence ladder through the entry (3 grids x 5 parents)", fontsize=13)
    fig.tight_layout()
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def figure_transom(points: dict[float, dict], target: Path) -> None:
    ratios = sorted(points)
    cwp = [1000.0 * points[r]["analytic_reference"]["waterplane_coefficient"] for r in ratios]
    cb = [points[r]["analytic_reference"]["block_coefficient"] for r in ratios]
    transom = [points[r]["transom"]["immersed_area_ratio_of_mid_section"] for r in ratios]
    deck = [points[r]["station_law"]["transom_ratio"] for r in ratios]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    axes[0].plot(ratios, cwp, "o-", color="#1f77b4")
    axes[0].set_xlabel("transom breadth ratio tau")
    axes[0].set_ylabel("Cwp x 1000 (closed form)")
    axes[0].set_title("Cwp is held while tau sweeps (the target does not move)")
    axes[0].grid(alpha=0.3)
    axes[1].plot(ratios, transom, "s-", color="#d62728", label="immersed transom area / midship area")
    axes[1].plot(ratios, deck, "o--", color="#2ca02c", label="breadth ratio tau (input)")
    axes[1].set_xlabel("transom breadth ratio tau")
    axes[1].set_ylabel("ratio")
    axes[1].set_title("Second parameter axis: stern fullness, Cb = 0.70")
    axes[1].grid(alpha=0.3)
    axes[1].legend(fontsize=7)
    fig.suptitle("Series 60 transom sweep through `series60 parent`", fontsize=13)
    fig.tight_layout()
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def write_csvs(ladder: dict[tuple[int, int], dict], out: Path) -> list[str]:
    written = []
    path = out / "five_cb_gradient_ladder.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["nx", "nz", "cb", "source", "coefficient", "published", "got",
                         "relative", "role"])
        for grid, document in ladder.items():
            for row in document["rows"]:
                for source, entries in row["error_vs_published"].items():
                    for name, entry in entries.items():
                        writer.writerow([grid[0], grid[1], row["cb_requested"], source, name,
                                         entry["published"], entry["got"],
                                         entry["relative"], entry["role"]])
    written.append(str(path))
    return written


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="runs/m19694_series60_family")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    ladder: dict[tuple[int, int], dict] = {}
    for grid in GRID_LADDER:
        document = family_record(grid, out / f"family_{grid[0]}x{grid[1]}")
        ladder[grid] = document
        worst, where = worst_published(document)
        print(f"  {grid[0]}x{grid[1]}: worst published {100.0 * worst:+.3f} % ({where or 'n/a'}), "
              f"worst discretisation {100.0 * worst_discretisation(document):.3f} %")

    transoms: dict[float, dict] = {}
    for tau in TRANSOM_POINTS:
        transoms[tau] = parent_record("--cb", "0.70", "--cwp", "0.753", "--lcb", "0.490",
                                      "--transom-ratio", f"{tau}", "--nx", "161", "--nz", "41",
                                      "--gate", out=out / f"parent_tau_{tau:.2f}")
        cwp = transoms[tau]["analytic_reference"]["waterplane_coefficient"]
        print(f"  tau={tau:.2f}: Cwp={cwp:.4f} "
              f"Cb={transoms[tau]['analytic_reference']['block_coefficient']:.4f} "
              f"transom/mid={transoms[tau]['transom']['immersed_area_ratio_of_mid_section']:.4f} "
              f"gate={transoms[tau]['gate']['passed']}")

    figure_gradient(ladder, figures / "five_cb_gradient.png")
    figure_convergence(ladder, figures / "convergence_ladder.png")
    figure_transom(transoms, figures / "transom_sweep.png")
    csvs = write_csvs(ladder, out)

    reference_grid = GRID_LADDER[-1]
    best_worst, best_where = worst_published(ladder[reference_grid])
    discretisation_by_grid = {f"{nx}x{nz}": worst_discretisation(document)
                              for (nx, nz), document in ladder.items()}
    refinements = [discretisation_by_grid[f"{nx}x{nz}"] for nx, nz in GRID_LADDER]

    parent_gate_ok = all(point["gate"]["passed"] for point in transoms.values())
    parent_discretisation_ok = all(
        point["discretisation"]["within_tolerance"] for point in transoms.values())

    criteria = {
        "five_cb_generated": {
            "criterion": "five Cb parents generated through the product entry",
            "evidence": [row["cb_requested"] for row in ladder[reference_grid]["rows"]],
            "passed": len(ladder[reference_grid]["rows"]) == 5,
        },
        "hydrostatics_vs_published": {
            "criterion": "generated vs published particulars within 3 %",
            "worst_relative": best_worst,
            "worst_at": best_where,
            "tolerance": PUBLISHED_TOLERANCE,
            "passed": abs(best_worst) <= PUBLISHED_TOLERANCE,
            "caveat": "Cwp and LCB/L are generator targets and Cp is the published ladder: "
                      "agreement there is a closure, not an external verification",
        },
        "mesh_vs_closed_form": {
            "criterion": "mesh integration within 2 % of the module's closed form",
            "worst_by_grid": discretisation_by_grid,
            "tolerance": DISCRETISATION_TOLERANCE,
            "passed": refinements[-1] <= DISCRETISATION_TOLERANCE,
            "monotone_refinement": refinements[0] >= refinements[1] >= refinements[2],
        },
        "second_parameter_axis": {
            "criterion": "transom sweep stays gated and inside the bar",
            "points": {f"{tau:.2f}": {
                "cwp": transoms[tau]["analytic_reference"]["waterplane_coefficient"],
                "block_coefficient": transoms[tau]["analytic_reference"]["block_coefficient"],
                "transom_area_ratio_of_mid_section":
                    transoms[tau]["transom"]["immersed_area_ratio_of_mid_section"],
                "gate_passed": transoms[tau]["gate"]["passed"],
            } for tau in TRANSOM_POINTS},
            "passed": parent_gate_ok and parent_discretisation_ok,
        },
        "product_entry": {
            "criterion": "the family is reachable from the console entry, not only by import",
            "entry": "tensorcad-hull series60 {parent,family}",
            "commands": [document["entry_run"]["command"] for document in ladder.values()]
                        + [point["entry_run"]["command"] for point in transoms.values()],
            "passed": all(document["entry_run"]["returncode"] == 0 for document in ladder.values())
                      and all(point["entry_run"]["returncode"] in (0, 1)
                              for point in transoms.values()),
        },
    }

    evidence = {
        "schema": "tensorcad.series60_family.round2/1.0",
        "mission": "m19694",
        "anchor": "GAP-CAD-T2-PARAM-SERIES",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "case": {"family": "Series 60 parent forms (Todd 1951)",
                 "length_m": L, "beam_m": B, "draft_m": T, "lb": L / B, "bt": B / T,
                 "cb_values": [0.60, 0.65, 0.70, 0.75, 0.80],
                 "grid_ladder": [list(grid) for grid in GRID_LADDER],
                 "transom_points": list(TRANSOM_POINTS)},
        "published_source": {"path": str(DATA.relative_to(ROOT)), "sha256": sha256_of(DATA),
                             "citation": json.loads(DATA.read_text())["provenance"]["citation"]},
        "tolerances": {"published": PUBLISHED_TOLERANCE,
                       "discretisation": DISCRETISATION_TOLERANCE},
        "criteria": criteria,
        "grids": {f"{nx}x{nz}": document for (nx, nz), document in ladder.items()},
        "transom_points": {f"{tau:.2f}": point for tau, point in transoms.items()},
        "artifacts": {"figures": [str(path) for path in sorted(figures.glob("*.png"))],
                      "csv": csvs},
        "repro_cmd": "python3 tools/series60_family/round2_evidence.py --out " + str(out),
        "external_anchor": {
            "available": "published parent particulars Cb/Cm/Cp/Cwp/LCB (Todd 1951, in-repo table)",
            "withheld": "a Todd station-offset table for the five parents: absent in-repo and no "
                        "network in this sandbox, so the family's *shape* is not externally anchored",
            "consequence": "the coefficient closure and the discretisation ladder are verified; "
                           "the sectional-area distribution is NOT",
        },
        "all_passed": all(entry["passed"] for entry in criteria.values()),
    }
    path = out / "round2_evidence.json"
    path.write_text(json.dumps(evidence, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"\nevidence -> {path}")
    for name, entry in criteria.items():
        flag = "PASS" if entry["passed"] else "FAIL"
        print(f"  [{flag}] {name}: {entry['criterion']}")
    return 0 if evidence["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())