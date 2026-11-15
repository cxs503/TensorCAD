#!/usr/bin/env python3
"""Round 4 of the fleet reference library: stop losing numerals the pages do print.

Four publisher tables are pinned in ``tensorcad.fleet.crossref.SOURCES`` as HTML.  Three of
them print a "Model scale | Full scale" pair, and a full-scale cell may say ``NA``.  The
reader used to read the neighbouring model column, flag the cell ``fallback_column`` and then
drop it from the verdict -- correct about the scale, but the numeral was then invisible.

This tool re-reads the pinned bytes (no network) and reports every cell that does not vote,
what the document's own model/full ratio makes of it, and how far the derived reading sits
from the published one.  Outputs: the JSON artefact, the markdown section and the figure.

    python tools/baseline_fleet/refs_round4_derived.py
"""

from __future__ import annotations

import dataclasses
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.fleet.crossref import (  # noqa: E402
    QUANTITIES,
    QUANTITY_LABELS,
    SOURCES,
    TOLERANCE,
    cross_reference,
    parse_source,
    split_ships,
    votes,
)

RAW = Path("artifacts/fleet_refs/raw")
OUT = Path("benchmarks/verified/fleet_refs/derived_cells.json")
DOC = Path("docs/fleet_refs_round4_derived.md")
FIG = Path("docs/figs/fleet_refs_round4_derived.png")


def main() -> int:
    extractions: dict[str, dict[str, dict]] = {}
    read: list[str] = []
    missing: list[str] = []
    for spec in SOURCES:
        if spec.kind != "html":
            continue
        path = RAW / (spec.key + ".html")
        if not path.exists():
            missing.append(spec.key)
            continue
        payload = path.read_bytes()
        read.append(spec.key)
        for ship, block in split_ships(spec, payload):
            cells = parse_source(dataclasses.replace(spec, ship=ship, section=None), block)
            if cells:
                extractions.setdefault(ship, {})[spec.key] = cells
    report = cross_reference(extractions)

    non_voting = []
    for ship in sorted(report):
        for quantity in QUANTITIES:
            block = report[ship]["quantities"][quantity]
            voters = block["family_values"]
            for row in block.get("derived", {}).get("cells", []):
                reference = sorted(voters.values())[len(voters) // 2] if voters else None
                delta = None if reference is None else round(row["value"] - reference, 6)
                non_voting.append(
                    {
                        "ship": ship,
                        "quantity": quantity,
                        "family": row["family"],
                        "value": row["value"],
                        "unit_basis": row["unit_basis"],
                        "derivation": row["derivation"],
                        "note": row["note"],
                        "voters": voters,
                        "delta_vs_voter": delta,
                        "relative_to_voter": None
                        if not reference
                        else round(abs(delta) / abs(reference), 5),
                        "tolerance": TOLERANCE[quantity],
                        "within_tolerance": None
                        if reference is None
                        else abs(delta) <= TOLERANCE[quantity] * abs(reference),
                    }
                )

    payload = {
        "schema": "tensorcad.fleet_refs.derived/1.0",
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "documents_read": read,
        "documents_not_cached": missing,
        "non_voting_cells": non_voting,
        "counts": {
            "non_voting": len(non_voting),
            "agreeing_with_a_voter": sum(1 for r in non_voting if r["within_tolerance"]),
            "disagreeing_with_a_voter": sum(1 for r in non_voting if r["within_tolerance"] is False),
            "no_voter": sum(1 for r in non_voting if r["within_tolerance"] is None),
        },
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n")

    lines = [
        "# Round 4 — the numerals the reference library was dropping",
        "",
        "Anchor **REPRO-GEOM-FLEET-REFS** (TensorCAD). Machine artefact: "
        "`benchmarks/verified/fleet_refs/derived_cells.json` (regenerate: "
        "`python tools/baseline_fleet/refs_round4_derived.py`). Every number below is rendered "
        "from that file.",
        "",
        f"Documents re-read from the pinned bytes (no network): {len(read)} "
        f"({', '.join(read) or 'none'}).",
        "",
    ]
    if missing:
        lines += [f"Not cached in this worktree, so not re-read: {', '.join(missing)}.", ""]
    lines += [
        "## 1. Cells that do not vote, and what they say",
        "",
        "| ship | quantity | basis | derived value | voters | Δ vs voter | within tolerance |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in non_voting:
        voters = ", ".join(f"{k}={v:g}" for k, v in sorted(row["voters"].items())) or "—"
        lines.append(
            f"| {row['ship']} | {QUANTITY_LABELS[row['quantity']]} | `{row['unit_basis']}` | "
            f"{row['value']:g} | {voters} | "
            f"{'—' if row['delta_vs_voter'] is None else f'{row[chr(39)+chr(39)] if False else row[chr(100)+chr(101)+chr(108)+chr(116)+chr(97)+chr(95)+chr(118)+chr(115)+chr(95)+chr(118)+chr(111)+chr(116)+chr(101)+chr(114)]:+g}'} | "
            f"{row['within_tolerance']} |"
        )
    lines += ["", "## 2. Derivations in full", ""]
    for row in non_voting:
        lines.append(f"* {row['ship']} {QUANTITY_LABELS[row['quantity']]}: {row['derivation']}")
        lines.append(f"  * note: {row['note']}")
    lines += [
        "",
        "## 3. Why these cells do not vote",
        "",
        "A published numeral votes only when the publisher prints it in a column whose scale is "
        "known. A model-scale numeral cannot vote for a full-scale quantity (that is the "
        "`fallback_column` rule, enforced by `tensorcad.fleet.crossref.votes`). The derived "
        "readings above use the *document's own* model/full ratio, so they are arithmetic on "
        "one publisher's numbers and one length assumption -- evidence about the document, not "
        "a second opinion about the ship. None of them is counted in criterion ①.",
        "",
    ]
    DOC.write_text("\n".join(lines) + "\n")

    # figure: LCB of ONRT as metres aft of FP and as % of length, from the two publishers
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    onrt = [r for r in non_voting if r["ship"] == "ONRT" and r["quantity"] == "lcb_percent"]
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))
    if onrt:
        row = onrt[0]
        lwl = extractions["ONRT"]["t2015_onrt"]["lwl"].value
        metres_derived = lwl / 2 - row["value"] / 100.0 * lwl
        voters = row["voters"]
        names, pct = [], []
        for family, value in sorted(voters.items()):
            names.append(f"{family}\n(published %L)")
            pct.append(value)
        names.append("NMRI-T2015\n(derived % of Lwl)")
        pct.append(row["value"])
        axes[0].bar(names, pct, color=["#4477aa", "#ee6677"])
        axes[0].axhline(0.0, color="k", lw=0.8)
        axes[0].set_ylabel("LCB (% of length, fwd +)")
        axes[0].set_title("ONRT LCB: published vs derived")
        for i, v in enumerate(pct):
            axes[0].text(i, v, f"{v:.3f}", ha="center", va="bottom" if v >= 0 else "top", fontsize=8)
        metres = [metres_derived]
        for family, value in sorted(voters.items()):
            metres.append(lwl / 2 - value / 100.0 * lwl)
        axes[1].bar(["NMRI-T2015\n(derived)", *[f"{f}\n(inverted)" for f in sorted(voters)]],
                    metres, color=["#ee6677", "#4477aa"])
        axes[1].set_ylabel("LCB (m aft of FP)")
        axes[1].set_title(f"physical value; span {max(metres) - min(metres):.3f} m")
        for i, v in enumerate(metres):
            axes[1].text(i, v, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=160)
    print(f"wrote {OUT}, {DOC}, {FIG}: {len(non_voting)} non-voting cells")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
