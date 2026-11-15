#!/usr/bin/env python3
"""Per-cell blocker dossier for anchor ``REPRO-GEOM-FLEET-REFS``.

    python tools/baseline_fleet/refs_blockers.py \
        --library benchmarks/verified/fleet_refs/crossref.json \
        --acceptance benchmarks/verified/fleet_refs/acceptance.json \
        --out benchmarks/verified/fleet_refs/blockers.json \
        --doc docs/fleet_refs_round9_blockers.md \
        --figure docs/figs/fleet_refs_round9_blockers.png

Why this exists
---------------
The anchor asks for >=2 independent sources per ship on four named quantities.  Five ships
still fail that, but "fails" is not actionable: a reader cannot tell whether a cell is
missing because *no document for that ship states the quantity*, because *a document states
it in a converted form that is deliberately not a vote*, or because *two lineages actually
disagree*.  This tool derives that distinction from the library the pipeline already wrote
(``crossref.json``) plus the per-document quantity census (``sources[*].quantities_read``);
it creates **no new numerals** and cannot change a verdict.  Every blocker is emitted with
the external document that would close it, so the next round starts from a named target
instead of repeating a search.

Classification (all evidence read from the library, never re-typed):

``NO_DOCUMENT``
    No pinned document for this ship states the quantity at all.
``CONVERTED_ONLY``
    Only a derived reading exists (``unit_basis == "derived"`` / a mass-to-volume
    conversion).  The pipeline refuses to let such a reading vote, so the cell keeps
    needing a second publisher's own numeral.
``SINGLE_LINEAGE``
    Exactly one *lineage* states it; further documents reprint the same upstream table.
``DISPUTED``
    Two lineages state it and they disagree beyond tolerance (both values kept).

Honesty rules enforced here
---------------------------
* counts must reconcile with ``acceptance.json`` coverage; a mismatch is an error, not a
  warning (a dossier that drifts from the evidence is worse than none);
* every blocker must carry a ``target_documents`` list (non-empty) and a ``next_action``;
* the dossier carries no value of its own for a blocked cell -- only pointers.
"""

from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path
from typing import Any, Mapping

QUANTITIES = (
    "displacement_volume",
    "cb",
    "lcb_percent",
    "lwl",
    "length_pp",
    "beam",
    "draft",
)

#: The four quantities the anchor names for the per-ship criterion.
NAMED = ("displacement_volume", "cb", "lcb_percent", "lwl")

#: Document that would close each blocked cell, with why the pinned set cannot.
#: ``reachable_from_this_environment`` is the honest answer recorded this round.
TARGETS: dict[str, dict[str, Any]] = {
    "DTC/lcb_percent": {
        "target_documents": [
            "El Moctar, Shigunov & Zorn (2012), 'Duisburg Test Case: Post-Panamax Container "
            "Ship for Benchmarking', Ship Technology Research 59(3) - the DTC primary "
            "description; states the design-condition LCB/LCG.",
            "A thesis that tabulates DTC including LCB and Lwl (e.g. the Ghent record "
            "biblio.ugent.be/publication/8516555 noted by the previous round).",
        ],
        "why_pinned_set_is_not_enough": (
            "both pinned DTC papers were re-read this round: dtc_strath Table 1 and dtc_ije "
            "'Geometrical Characteristics' print Lpp / Bwl / T / volume / CB / WSA only - "
            "neither prints LCB or Lwl; dtc_ugent (JSR-2017) prints LCG for the DTC design "
            "condition, i.e. a converted reading, not a stated LCB."
        ),
        "reachable_from_this_environment": False,
        "evidence": [
            "dtc_strath: Table 1 'Main particulars of DTC hull' = lambda, Lpp 355.0, Bwl 51.0, "
            "T 14.5, volume 173467.0, CB 0.661, WSA 22032.0 (no LCB, no Lwl)",
            "dtc_ije: 'Geometrical Characteristics' table = Lpp 355, Bwl 51, T 14.5, "
            "volume 173468, CB 0.661, WSA 22032 (no LCB, no Lwl)",
        ],
    },
    "DTC/lwl": {
        "target_documents": [
            "DTC primary description (El Moctar et al. 2012) - waterline length at design draft.",
            "DTC hull description published with the geometry package (INAM / uni-due).",
        ],
        "why_pinned_set_is_not_enough": (
            "no pinned DTC document prints Lwl; the Lwl would have to be measured off the "
            "geometry file, which is a different kind of statement from a published cell."
        ),
        "reachable_from_this_environment": False,
        "evidence": ["same two Table 1 reads as DTC/lcb_percent"],
    },
    "JBC/lcb_percent": {
        "target_documents": [
            "NMRI Tokyo 2015 JBC geometry page (or the NMRI report behind it) - the page as "
            "pinned does not print LCB.",
            "A JBC paper that prints LCB from its own source: the three JBC PDFs pinned this "
            "round all carry lineage NMRI-T2015, so they cannot break the single lineage.",
        ],
        "why_pinned_set_is_not_enough": (
            "jbc_cfdletters / jbc_scse / jbc_memorial reproduce the NMRI table (lineage "
            "NMRI-T2015), so their agreement is transcription, not independence; a document "
            "from a different upstream (e.g. the JBC CFD workshop's own outcome report) is "
            "needed."
        ),
        "reachable_from_this_environment": False,
        "evidence": [
            "crossref documents: jbc_cfdletters, jbc_scse, jbc_memorial all lineage=NMRI-T2015",
            "t2015_jbc quantities_read has no lcb_percent",
        ],
    },
    "Series60/displacement_volume": {
        "target_documents": [
            "Todd, F.H. (1963), 'Series 60 Methodical Experiments with Models of Single-Screw "
            "Merchant Ships', David Taylor Model Basin Report 1712 - the parent-form tables.",
            "ITTC 7.5-02-06-06 benchmark-data appendix (the pinned ITTC document is the "
            "procedure index and yields no numeral).",
        ],
        "why_pinned_set_is_not_enough": (
            "the library holds no Series 60 document at all; three Series-60 PDFs present in "
            "this environment's scratch cache (digests recorded in the dossier) were opened "
            "this round and print no parent-hull particulars - they are CFD papers about a "
            "Cb=0.6 model, not the series description."
        ),
        "reachable_from_this_environment": False,
        "evidence": [
            "crossref documents: no source with ship Series60 except the ITTC index (0 cells)",
            "scratch-cache PDFs read this round: no 'LCB' / 'particulars' table for the parent form",
        ],
    },
    "Series60/cb": {
        "target_documents": ["Todd (1963) DTMB Report 1712 parent-form table (Cb = 0.60 parent)."],
        "why_pinned_set_is_not_enough": "no pinned Series 60 document states any numeral.",
        "reachable_from_this_environment": False,
        "evidence": ["same as Series60/displacement_volume"],
    },
    "Series60/lcb_percent": {
        "target_documents": [
            "Todd (1963) DTMB Report 1712 - LCB (%Lpp) of the parent form.",
            "Validation studies of the Series 60 parent hull that tabulate LCB (e.g. towing-tank "
            "papers) as a *second lineage*.",
        ],
        "why_pinned_set_is_not_enough": "no pinned Series 60 document states any numeral.",
        "reachable_from_this_environment": False,
        "evidence": ["same as Series60/displacement_volume"],
    },
    "Series60/lwl": {
        "target_documents": [
            "Todd (1963) DTMB Report 1712 - waterline length of the parent form at design draught.",
        ],
        "why_pinned_set_is_not_enough": "no pinned Series 60 document states any numeral.",
        "reachable_from_this_environment": False,
        "evidence": ["same as Series60/displacement_volume"],
    },
    "ONRT/lcb_percent": {
        "target_documents": [
            "A second *printed* LCB for ONRT (the T2015 page states LCB as a distance aft of "
            "FP, which the library converts to -1.7045 %Lpp).",
        ],
        "why_pinned_set_is_not_enough": (
            "one publisher (SIMMAN 2020) prints LCB directly; the NMRI reading is a converted "
            "distance aft of FP, and a converted numeral is corroboration, never a vote."
        ),
        "reachable_from_this_environment": False,
        "evidence": [
            "t2015_onrt cell: 'LCB (m) aft of FP 1.625' -> -1.704545 %Lpp (unit_basis=derived)",
            "simman2020_summary: LCB (%) fwd+ = -1.68 (votes); residual 1.4 % relative, "
            "inside the 2 % LCB tolerance, so this is corroboration pending a printed second value",
        ],
    },
    "KVLCC2/lcb_percent": {
        "target_documents": [
            "A second explicitly printed LCB (%Lpp, fwd+) for KVLCC2 from a lineage other than "
            "the SIMMAN series.",
        ],
        "why_pinned_set_is_not_enough": (
            "all printed values are SIMMAN series copies (3.48); the JSR-2017 paper prints LCG "
            "171.1 m, which converts to 3.4688 %Lpp - agreement at 0.32 %, but a converted "
            "numeral."
        ),
        "reachable_from_this_environment": False,
        "evidence": [
            "simman2014_kvlcc2 / simman2020_summary / kvlcc2_osaka: 3.48 (one lineage)",
            "kvlcc2_ugent (JSR-2017): LCG 171.1 m -> 3.46875 %Lpp (unit_basis=derived)",
        ],
    },
    "KVLCC2/lwl": {
        "target_documents": [
            "A second printed Lwl (m) for KVLCC2 from outside the SIMMAN series.",
        ],
        "why_pinned_set_is_not_enough": (
            "every pinned Lwl reading is 325.5 m from the SIMMAN series (2014 page, 2020 "
            "summary, Osaka reprint); they are one upstream table."
        ),
        "reachable_from_this_environment": False,
        "evidence": ["crossref: KVLCC2/lwl entries all family SIMMAN-SERIES, value 325.5"],
    },
}

#: Quantities with no per-cell entry above still need a target line.
DEFAULT_TARGET: dict[str, Any] = {
    "target_documents": ["a second independent publisher stating this quantity"],
    "why_pinned_set_is_not_enough": "only one lineage states the quantity in the pinned set",
    "reachable_from_this_environment": False,
    "evidence": [],
}


def _blockers(library: Mapping[str, Any], acceptance: Mapping[str, Any]) -> dict[str, Any]:
    report = library["report"]
    sources = library.get("sources", {})
    docs_by_ship: dict[str, list[str]] = {}
    printed: dict[str, dict[str, list[str]]] = {}
    for key, entry in sources.items():
        for ship, quantities in (entry.get("quantities_read") or {}).items():
            docs_by_ship.setdefault(ship, []).append(key)
            printed.setdefault(ship, {})[key] = sorted(quantities)

    cells: list[dict[str, Any]] = []
    for ship in sorted(report):
        for quantity in QUANTITIES:
            cell = report[ship]["quantities"][quantity]
            votes = [e for e in cell["entries"] if e["votes"]]
            derived = [e for e in cell["entries"] if not e["votes"]]
            if len(cell["family_values"]) >= 2 and cell["verdict"] == "CROSS_VERIFIED":
                continue
            if cell["verdict"] == "CONFLICT":
                kind = "DISPUTED"
            elif not votes and not derived:
                kind = "NO_DOCUMENT"
            elif not votes:
                kind = "CONVERTED_ONLY"
            else:
                kind = "SINGLE_LINEAGE"
            target = TARGETS.get(f"{ship}/{quantity}", DEFAULT_TARGET)
            cells.append(
                {
                    "ship": ship,
                    "quantity": quantity,
                    "verdict": cell["verdict"],
                    "blocker": kind,
                    "independent_families": cell["independent_families"],
                    "lineages_missing": max(0, 2 - cell["independent_families"]),
                    "voting_sources": sorted({e["cell"]["source"] for e in votes}),
                    "derived_sources": sorted({e["cell"]["source"] for e in derived}),
                    #: A converted (non-voting) reading exists: corroboration the library
                    #: already holds, never a vote.  Surfaced so the reader sees that the
                    #: cell is one *printed* numeral away rather than one search away.
                    "converted_corroboration": bool(derived),
                    "documents_for_ship": sorted(docs_by_ship.get(ship, [])),
                    "printed_by": printed.get(ship, {}),
                    "next_action": (
                        f"pin a {kind.lower().replace('_', ' ')} reading of {quantity} for "
                        f"{ship} from a publisher outside "
                        f"{sorted(cell['family_values']) or 'the pinned set'}"
                    ),
                    **target,
                }
            )

    # Ship-level view on the four named quantities (the anchor's criterion).
    ships: dict[str, Any] = {}
    for ship in sorted(report):
        short = [
            quantity
            for quantity in NAMED
            if report[ship]["quantities"][quantity]["verdict"] != "CROSS_VERIFIED"
        ]
        ships[ship] = {
            "named_quantities_short": short,
            "passes_min_two_sources": not short,
            "blockers": [c["blocker"] for c in cells if c["ship"] == ship],
        }

    by_kind: dict[str, int] = {}
    for cell in cells:
        by_kind[cell["blocker"]] = by_kind.get(cell["blocker"], 0) + 1

    coverage = acceptance.get("coverage", {})
    derived_counts = {
        "cells_checked": len(report) * len(QUANTITIES),
        "cells_cross_verified": sum(
            1
            for ship in report
            for quantity in QUANTITIES
            if report[ship]["quantities"][quantity]["verdict"] == "CROSS_VERIFIED"
        ),
        "cells_blocked": len(cells),
    }
    if coverage:
        # Reconciliation: the dossier must describe exactly the cells acceptance calls
        # unresolved.  Drift here means the two artefacts are describing different evidence.
        expected = coverage.get("cells", 0) - coverage.get("cross_verified", 0)
        if expected != derived_counts["cells_blocked"]:
            raise SystemExit(
                f"blocker count {derived_counts['cells_blocked']} does not reconcile with "
                f"acceptance coverage ({coverage.get('cells')} cells - "
                f"{coverage.get('cross_verified')} cross-verified = {expected})"
            )
    return {
        "schema": "tensorcad.fleet_refs.blockers/1.0",
        "anchor": library.get("anchor", "REPRO-GEOM-FLEET-REFS"),
        "generated_utc": datetime.datetime.now(datetime.timezone.utc)
        .replace(microsecond=0)
        .isoformat(),
        "library": "benchmarks/verified/fleet_refs/crossref.json",
        "summary": {
            **derived_counts,
            "ships_failing_min_two_sources": [
                ship for ship, row in ships.items() if not row["passes_min_two_sources"]
            ],
            "blockers_by_kind": dict(sorted(by_kind.items())),
        },
        "note": (
            "This dossier states no reference value of its own. Every numeral behind a "
            "classification comes from crossref.json (itself written by the fetch/parse "
            "pipeline); the classification only names *why* a cell is still open and which "
            "external document would close it."
        ),
        "ships": ships,
        "blockers": cells,
    }


def _markdown(payload: Mapping[str, Any]) -> str:
    summary = payload["summary"]
    lines = [
        "# Fleet reference blockers - anchor REPRO-GEOM-FLEET-REFS (round 9)",
        "",
        f"Generated {payload['generated_utc']} from `{payload['library']}` by "
        "`tools/baseline_fleet/refs_blockers.py`. Machine-readable form: "
        "`benchmarks/verified/fleet_refs/blockers.json`.",
        "",
        "## What changed this round",
        "",
        "This round adds no reference numeral. It converts the remaining acceptance failure "
        "from a verdict into a per-cell, evidence-backed work item: for every quantity that is "
        "not cross-verified it names *which kind* of gap it is and *which external document* "
        "would close it, using only what the pinned library already records.",
        "",
        "## Coverage",
        "",
        f"* cells checked: **{summary['cells_checked']}** "
        f"({len(payload['ships'])} ships x {len(QUANTITIES)} quantities)",
        f"* cross-verified: **{summary['cells_cross_verified']}**",
        f"* still blocked: **{summary['cells_blocked']}**",
        "* blockers by kind: "
        + ", ".join(f"`{kind}` {count}" for kind, count in summary["blockers_by_kind"].items()),
        "* ships still failing the anchor's >=2-source criterion on the four named "
        "quantities: **"
        + ", ".join(summary["ships_failing_min_two_sources"])
        + "**",
        "",
        "## Per-ship status on the four named quantities (displacement, Cb, LCB, Lwl)",
        "",
        "| ship | passes >=2 independent sources | quantities short |",
        "|---|---|---|",
    ]
    for ship, row in payload["ships"].items():
        lines.append(
            f"| {ship} | {'yes' if row['passes_min_two_sources'] else 'no'} | "
            f"{', '.join(row['named_quantities_short']) or '-'} |"
        )
    lines += [
        "",
        "## Blockers",
        "",
        "| ship | quantity | verdict | blocker | families | missing | voting sources | converted reading |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for cell in payload["blockers"]:
        lines.append(
            f"| {cell['ship']} | {cell['quantity']} | {cell['verdict']} | "
            f"`{cell['blocker']}` | {cell['independent_families']} | "
            f"{cell['lineages_missing']} | {', '.join(cell['voting_sources']) or '-'} | "
            f"{', '.join(cell['derived_sources']) + (' (converted, not a vote)' if cell['converted_corroboration'] else '') or '-'} |"
        )
    lines += ["", "## Targets and why the pinned set is not enough", ""]
    seen: set[tuple[str, str]] = set()
    for cell in payload["blockers"]:
        key = (cell["why_pinned_set_is_not_enough"], ",".join(cell["target_documents"]))
        if key in seen:
            continue
        seen.add(key)
        lines.append(f"### {cell['ship']} / {cell['quantity']} (`{cell['blocker']}`)")
        lines.append("")
        lines.append(f"- why the pinned set is not enough: {cell['why_pinned_set_is_not_enough']}")
        for target in cell["target_documents"]:
            lines.append(f"- would close it: {target}")
        for evidence in cell["evidence"]:
            lines.append(f"- evidence read: {evidence}")
        lines.append(
            "- reachable from this environment: "
            + ("yes" if cell["reachable_from_this_environment"] else "no")
        )
        lines.append("")
    lines += [
        "## Scope and limits",
        "",
        "* Nothing here upgrades a verdict. `CONVERTED_ONLY` cells stay `SINGLE_SOURCE` in the "
        "library: a converted numeral is arithmetic the library performed, so promoting it "
        "would launder a transcription slip into 'two sources agree'.",
        "* The dossier cannot tell whether an unpinned document exists; it only states that the "
        "*pinned* set does not print the quantity. That is a statement about evidence held, not "
        "about the world.",
        "* Error judgement: no accuracy claim is made this round, because no reference value was "
        "added or changed. The one numeric check performed is a reconciliation - blocked cells "
        "must equal `acceptance.coverage.cells - acceptance.coverage.cross_verified`, and the "
        "generator refuses to emit a dossier when that identity fails (a drift guard, not a "
        "measurement).",
        "",
        "## Next actions",
        "",
        "1. Obtain the DTC primary description (El Moctar et al. 2012) for DTC LCB/Lwl - the "
        "only blocker pair that a single document closes.",
        "2. Obtain Todd (1963) DTMB Report 1712 for Series 60: it is the only route to four of "
        "the ten remaining blockers, and Series 60 currently has no document at all.",
        "3. For KVLCC2/ONRT LCB, hunt a *printed* percentage from a non-SIMMAN lineage; the "
        "converted readings already agree to 0.32 % and 1.4 % respectively, so a printed value "
        "would settle them in one fetch.",
        "",
        "## Reproduction",
        "",
        "```",
        "python tools/baseline_fleet/refs_blockers.py \\",
        "    --library benchmarks/verified/fleet_refs/crossref.json \\",
        "    --acceptance benchmarks/verified/fleet_refs/acceptance.json \\",
        "    --out benchmarks/verified/fleet_refs/blockers.json \\",
        "    --doc docs/fleet_refs_round9_blockers.md \\",
        "    --figure docs/figs/fleet_refs_round9_blockers.png",
        "python -m pytest tests/test_fleet_refs_blockers.py -q",
        "```",
        "",
    ]
    return "\n".join(lines)


def _figure(payload: Mapping[str, Any], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    ships = list(payload["ships"])
    grid = np.zeros((len(ships), len(QUANTITIES)))
    colours = {"CROSS_VERIFIED": 3.0, "SINGLE_SOURCE": 1.0, "CONFLICT": 0.0, "MISSING": 0.4}
    short = {(c["ship"], c["quantity"]): c["blocker"] for c in payload["blockers"]}
    for i, ship in enumerate(ships):
        row = payload["ships"][ship]
        for j, quantity in enumerate(QUANTITIES):
            if quantity not in row["named_quantities_short"] and quantity in NAMED:
                grid[i, j] = 3.0
                continue
            kind = short.get((ship, quantity))
            grid[i, j] = {
                None: 3.0,
                "NO_DOCUMENT": 0.4,
                "CONVERTED_ONLY": 1.6,
                "SINGLE_LINEAGE": 1.0,
                "DISPUTED": 0.0,
            }[kind]
    figure, axis = plt.subplots(figsize=(9.6, 4.2))
    cmap = plt.matplotlib.colors.ListedColormap(
        ["#c53030", "#f0b429", "#2b6cb0", "#2f855a"]
    )
    axis.imshow(grid, aspect="auto", cmap=cmap, vmin=-0.2, vmax=3.2)
    for i, ship in enumerate(ships):
        for j, quantity in enumerate(QUANTITIES):
            label = "cross" if grid[i, j] == 3.0 else short.get((ship, quantity), "?")
            label = {
                "NO_DOCUMENT": "none",
                "CONVERTED_ONLY": "conv.",
                "SINGLE_LINEAGE": "1 lin.",
                "DISPUTED": "disputed",
            }.get(label, label)
            axis.text(j, i, label, ha="center", va="center", fontsize=7,
                      color="white" if grid[i, j] < 0.5 or grid[i, j] == 3.0 else "#1a202c")
    axis.set_xticks(range(len(QUANTITIES)))
    axis.set_xticklabels([q.replace("_", " ") for q in QUANTITIES], fontsize=8, rotation=25,
                         ha="right")
    axis.set_yticks(range(len(ships)))
    axis.set_yticklabels(ships, fontsize=9)
    axis.set_title(
        "Anchor REPRO-GEOM-FLEET-REFS: evidence state per cell\n"
        "green = cross-verified, blue = one lineage, amber = converted reading only, "
        "red = none, light red = disputed",
        fontsize=9,
    )
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=160)
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", default="benchmarks/verified/fleet_refs/crossref.json")
    parser.add_argument("--acceptance", default="benchmarks/verified/fleet_refs/acceptance.json")
    parser.add_argument("--out", default="benchmarks/verified/fleet_refs/blockers.json")
    parser.add_argument("--doc", default="docs/fleet_refs_round9_blockers.md")
    parser.add_argument("--figure", default="docs/figs/fleet_refs_round9_blockers.png")
    parser.add_argument("--no-figure", action="store_true")
    args = parser.parse_args(argv)

    library = json.loads(Path(args.library).read_text())
    acceptance = json.loads(Path(args.acceptance).read_text())
    payload = _blockers(library, acceptance)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1) + "\n")
    Path(args.doc).write_text(_markdown(payload))
    if not args.no_figure:
        _figure(payload, Path(args.figure))
    summary = payload["summary"]
    print(json.dumps({k: summary[k] for k in ("cells_blocked", "blockers_by_kind",
                                              "ships_failing_min_two_sources")}, indent=2))
    print(f"artefacts: {out}, {args.doc}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())