#!/usr/bin/env python3
"""Round 6 of the fleet reference library: report the *document* axis next to the lineage axis.

The library already collapses sources to *lineages* -- two documents that print one upstream
table are one source -- and the anchor's criterion ① is judged on that strict axis.  The
collapse is right, but it is also lossy in one direction: ``JBC`` shows ``SINGLE_SOURCE (1)``
although three independent publishers (a journal, a conference proceedings and the workshop
page itself) state the same numerals.  A reader cannot tell "only one document mentions this"
from "several documents mention it but they all reprint the workshop table", and the second
situation is very different evidence.

This round adds the missing axis, in place, without re-fetching anything:

* ``documents``   -- per ship/quantity: how many distinct pinned documents vote, which
  publishers they come from, which lineages they collapse to, and whether the agreement is
  ``lineage_limited`` (≥2 documents, 1 lineage: corroboration of transcription, *not* a second
  measurement);
* ``gap_ledger.json`` -- for every cell that still has fewer than two documents, the search
  this mission actually ran and the open-access candidates it returned, so the next round does
  not repeat the hunt blind.

Run it from the repository root:

    python tools/baseline_fleet/refs_round6_documents.py
    python tools/baseline_fleet/refs_acceptance.py check      # re-score the anchor
    python tools/baseline_fleet/backfill_refs.py backfill     # carry the axis into fleet_specs

No numeral is transcribed by hand; every count below is computed from the committed library.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Mapping

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.fleet.crossref import FAMILIES, LINEAGE, QUANTITIES, QUANTITY_LABELS  # noqa: E402

LIBRARY = ROOT / "benchmarks/verified/fleet_refs/crossref.json"
LEDGER = ROOT / "benchmarks/verified/fleet_refs/gap_ledger.json"
DOC = ROOT / "docs/fleet_refs_round6_documents.md"
FIG = ROOT / "docs/figs/fleet_refs_round6_documents.png"

#: Quantities the anchor names for criterion ①.
NAMED = ("displacement_volume", "cb", "lcb_percent", "lwl")

#: Searches this mission ran (OpenAlex, ``is_oa:true``) and the open-access candidates they
#: returned.  Nothing below was fetched, pinned or read in this mission: a candidate is a
#: hunch until its SHA-256 is pinned and its quote read in ``tensorcad.fleet.crossref``.
PROBES: tuple[dict[str, Any], ...] = (
    {
        "ship": "ONRT",
        "quantities": ["lcb_percent", "length_pp"],
        "query": "ONR tumblehome DTMB 5613 main particulars LCB",
        "candidates": [],
        "outcome": "no open-access work returned by the API for this query",
        "next_action": "hunt the SIMMAN 2020 ONRT detail page and DTMB 5613 model reports; the "
        "numerals wanted are Lpp and LCB (% Lpp)",
    },
    {
        "ship": "DTC",
        "quantities": ["lcb_percent", "lwl"],
        "query": "Duisburg Test Case DTC container ship LCB waterline length",
        "candidates": [
            "https://biblio.ugent.be/publication/8516555",
            "https://strathprints.strath.ac.uk/54642/1/Tezdogan_etal_JMST_2015_A_numerical_investigation_of_the_squat_and_resistance_of_ships_advancing.pdf",
        ],
        "outcome": "the Strathclyde paper is already pinned as dtc_strath and does not print LCB/Lwl "
        "in a full-scale column; the Ghent record was not fetched",
        "next_action": "read the DTC primary description (Moctar/Shigunov/Zorn, Ship Technology "
        "Research 2012) or a thesis that tabulates DTC including LCB and Lwl",
    },
    {
        "ship": "Series60",
        "quantities": list(QUANTITIES),
        "query": "Series 60 Cb 0.6 hull main particulars table resistance",
        "candidates": [
            "https://www.jstage.jst.go.jp/article/jjasnaoe1968/1975/138/1975_138_23/_pdf",
            "https://www.jstage.jst.go.jp/article/jjasnaoe1968/1978/144/1978_144_1/_pdf",
            "https://ynu.repo.nii.ac.jp/record/4500/files/11508181-01.pdf",
            "https://www.mdpi.com/2077-1312/10/5/686/pdf?version=1652864259",
        ],
        "outcome": "all four are candidate documents; none fetched in this mission, so Series 60 "
        "still has zero pinned document and the parent-form table in "
        "docs/audit_m17434/data/series60_published_parents.json remains a non-voting source",
        "next_action": "Series 60 is a *methodical series*: a document that reprints Todd's parent "
        "table gives one lineage; a second lineage needs an independent tabulation (e.g. an ITTC "
        "or SNAME procedure reproducing the parent-form coefficients)",
    },
)


def document_axis(library: Mapping[str, Any]) -> dict[str, Any]:
    """Per ship/quantity: documents, publishers, lineages and the ``lineage_limited`` flag."""
    out: dict[str, Any] = {}
    for ship, report in library["report"].items():
        per_quantity: dict[str, Any] = {}
        for quantity in QUANTITIES:
            row = report["quantities"][quantity]
            votes = [entry for entry in row["entries"] if entry.get("votes", True)]
            documents = sorted({entry["cell"]["source"] for entry in votes})
            publishers = sorted({FAMILIES.get(key, LINEAGE.get(key, key)) for key in documents})
            lineages = sorted({LINEAGE.get(key, FAMILIES.get(key, key)) for key in documents})
            if not documents:
                verdict = "MISSING"
            elif len(documents) == 1:
                verdict = "SINGLE_DOCUMENT"
            elif row["verdict"] == "CONFLICT":
                verdict = "CONFLICT"
            elif len(lineages) >= 2:
                verdict = "INDEPENDENT"
            else:
                # Two publishers, one origin: the numerals agree because one table was copied.
                verdict = "LINEAGE_LIMITED"
            per_quantity[quantity] = {
                "documents": len(documents),
                "document_keys": documents,
                "publishers": len(publishers),
                "publisher_families": publishers,
                "lineages": lineages,
                "lineage_limited": verdict == "LINEAGE_LIMITED",
                "verdict": verdict,
                # The strict verdict the library already publishes stays untouched.
                "lineage_verdict": row["verdict"],
            }
        out[ship] = per_quantity
    return out


def summary(axis: Mapping[str, Any]) -> dict[str, Any]:
    ships: dict[str, Any] = {}
    for ship, per_quantity in axis.items():
        named = [per_quantity[q] for q in NAMED]
        ships[ship] = {
            "documents_on_named": sum(1 for row in named if row["documents"] >= 2),
            "independent_on_named": sum(1 for row in named if row["verdict"] == "INDEPENDENT"),
            "lineage_limited_on_named": sum(1 for row in named if row["lineage_limited"]),
            "single_document_on_named": sum(1 for row in named if row["verdict"] == "SINGLE_DOCUMENT"),
            "missing_on_named": sum(1 for row in named if row["verdict"] == "MISSING"),
            "documents_total": len(
                {key for row in per_quantity.values() for key in row["document_keys"]}
            ),
        }
    return {
        "named_quantities": list(NAMED),
        "ships": ships,
        "ships_with_two_documents_on_all_named": sorted(
            ship for ship, row in ships.items() if row["documents_on_named"] == len(NAMED)
        ),
        "ships_independent_on_all_named": sorted(
            ship for ship, row in ships.items() if row["independent_on_named"] == len(NAMED)
        ),
    }


def gap_ledger(axis: Mapping[str, Any]) -> dict[str, Any]:
    gaps: list[dict[str, Any]] = []
    for ship, per_quantity in axis.items():
        for quantity in NAMED:
            row = per_quantity[quantity]
            if row["documents"] >= 2:
                continue
            probe = next(
                (p for p in PROBES if p["ship"] == ship and quantity in p["quantities"]), None
            )
            gaps.append(
                {
                    "ship": ship,
                    "quantity": quantity,
                    "label": QUANTITY_LABELS[quantity],
                    "documents": row["documents"],
                    "verdict": row["verdict"],
                    "probe": probe,
                }
            )
    return {
        "schema": "tensorcad.fleet_refs.gap_ledger/1.0",
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "mission": "m22012",
        "rule": "a probe is search provenance, never evidence: no candidate listed here was "
        "fetched or pinned, so it cannot move a verdict",
        "documents_pinned_total": sum(
            row["documents"] for ship in axis.values() for row in ship.values()
        ),
        "gaps": gaps,
        "counts": {
            "gaps": len(gaps),
            "ships_with_gaps": sorted({gap["ship"] for gap in gaps}),
        },
    }


def markdown(library: Mapping[str, Any], axis: Mapping[str, Any], summary_rows: Mapping[str, Any],
             ledger: Mapping[str, Any]) -> str:
    matrix = ["| ship | " + " | ".join(QUANTITY_LABELS[q] for q in QUANTITIES) + " |",
              "|---" * (len(QUANTITIES) + 1) + "|"]
    for ship in library["report"]:
        cells = []
        for quantity in QUANTITIES:
            row = axis[ship][quantity]
            cells.append(f"{row['verdict']} ({row['documents']}d/{len(row['lineages'])}l)")
        matrix.append(f"| {ship} | " + " | ".join(cells) + " |")
    score = ["| ship | documents on ∇/Cb/LCB/Lwl | independent (≥2 lineages) | lineage-limited | "
             "single document | missing | documents in library |",
             "|---|" + "---|" * 6]
    for ship, row in summary_rows["ships"].items():
        score.append(
            f"| {ship} | {row['documents_on_named']}/4 | {row['independent_on_named']}/4 | "
            f"{row['lineage_limited_on_named']} | {row['single_document_on_named']} | "
            f"{row['missing_on_named']} | {row['documents_total']} |"
        )
    gaps = ["| ship | quantity | documents | verdict | search run this mission | candidates |",
            "|---|---|---|---|---|---|"]
    for gap in ledger["gaps"]:
        probe = gap["probe"] or {}
        gaps.append(
            f"| {gap['ship']} | {gap['label']} | {gap['documents']} | {gap['verdict']} | "
            f"{probe.get('query', '—')} | {len(probe.get('candidates', []))} |"
        )
    return "\n".join(
        [
            "# Fleet reference library — the document axis (round 6)",
            "",
            "Anchor **REPRO-GEOM-FLEET-REFS** (TensorCAD). Machine artefacts:",
            "`benchmarks/verified/fleet_refs/crossref.json` (new `documents` block),",
            "`benchmarks/verified/fleet_refs/gap_ledger.json`.",
            "Regenerate: `python tools/baseline_fleet/refs_round6_documents.py`.",
            "Every digit below is computed from those files.",
            "",
            "## 1. Why a second axis",
            "",
            "Criterion ① is judged on *lineages*: two documents that print one upstream table are",
            "one source, because their agreement is a statement about copying, not about the hull.",
            "That rule is kept. But it hid the documents: `JBC` was reported `SINGLE_SOURCE` while",
            "three publishers state its particulars. The two axes answer different questions —",
            "",
            "* *documents*: is the numeral printed by more than one body?",
            "* *lineages*: is there an origin independent of the first one?",
            "",
            "## 2. Per-cell verdict (`verdict` = documents × lineages, `lineage_verdict` = library)",
            "",
            *matrix,
            "",
            "`Nd/Ml` = N voting documents from M lineages. `INDEPENDENT` needs ≥2 lineages;",
            "`LINEAGE_LIMITED` means ≥2 documents but one origin — corroboration of transcription,",
            "never a second measurement.",
            "",
            "## 3. The anchor's four named quantities",
            "",
            *score,
            "",
            "Ships with ≥2 documents on all four named quantities: "
            + (", ".join(summary_rows["ships_with_two_documents_on_all_named"]) or "none"),
            ".",
            "Ships independent (≥2 lineages) on all four: "
            + (", ".join(summary_rows["ships_independent_on_all_named"]) or "none"),
            ".",
            "The first list is the owner's literal wording (`≥2 独立来源`, journal papers and",
            "workshop pages both named as sources); the second is the strict reading the library",
            "already enforced. Neither is presented as the other.",
            "",
            "## 4. What this round did *not* achieve",
            "",
            f"{ledger['counts']['gaps']} named cells still carry fewer than two documents. The",
            "search actually run is recorded per cell in `gap_ledger.json`; no candidate was",
            "fetched or pinned, so no verdict moved. In particular Series 60 still has no pinned",
            "document, and DTC has no document at all for LCB or Lwl.",
            "",
            *gaps,
            "",
            "## 5. How this was verified",
            "",
            "```",
            "python tools/baseline_fleet/refs_round6_documents.py",
            "python -m pytest tests/test_fleet_refs_document_axis.py tests/test_fleet_refs_crossref.py \\",
            "                 tests/test_fleet_refs_adjudication.py tests/test_fleet_refs_freshness.py -q",
            "python tools/baseline_fleet/refs_acceptance.py check",
            "```",
            "",
            "![documents vs lineages](docs/figs/fleet_refs_round6_documents.png)",
            "",
        ]
    )


def figure(summary_rows: Mapping[str, Any], path: Path) -> str:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ships = list(summary_rows["ships"])
    docs = [summary_rows["ships"][s]["documents_on_named"] for s in ships]
    indep = [summary_rows["ships"][s]["independent_on_named"] for s in ships]
    limited = [summary_rows["ships"][s]["lineage_limited_on_named"] for s in ships]
    x = range(len(ships))
    fig, ax = plt.subplots(figsize=(8, 3.6))
    ax.bar([i - 0.26 for i in x], docs, width=0.26, label="≥2 documents", color="#2c7fb8")
    ax.bar(list(x), indep, width=0.26, label="≥2 lineages (strict)", color="#7fb800")
    ax.bar([i + 0.26 for i in x], limited, width=0.26, label="lineage-limited", color="#c0392b")
    ax.set_xticks(list(x))
    ax.set_xticklabels(ships)
    ax.set_ylabel("named quantities satisfied (of 4)")
    ax.set_ylim(0, 4.8)
    ax.set_title("REPRO-GEOM-FLEET-REFS: document axis vs lineage axis")
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return str(path)


def main() -> int:
    library = json.loads(LIBRARY.read_text())
    axis = document_axis(library)
    rows = summary(axis)
    ledger = gap_ledger(axis)
    library["documents"] = {
        "schema": "tensorcad.fleet_refs.documents/1.0",
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "note": "document-level independence, additive to the lineage verdict in `report`; "
                "LINEAGE_LIMITED = ≥2 documents, 1 origin (transcription corroboration)",
        "summary": rows,
        "cells": axis,
    }
    LIBRARY.write_text(json.dumps(library, indent=1) + "\n")
    LEDGER.write_text(json.dumps(ledger, indent=1) + "\n")
    DOC.write_text(markdown(library, axis, rows, ledger))
    figure(rows, FIG)
    print(f"wrote {LIBRARY.relative_to(ROOT)} (documents block)")
    print(f"wrote {LEDGER.relative_to(ROOT)} ({ledger['counts']['gaps']} open gaps)")
    print(f"wrote {DOC.relative_to(ROOT)}")
    print("≥2 documents on all four named:", ", ".join(rows["ships_with_two_documents_on_all_named"]))
    print("≥2 lineages on all four named:", ", ".join(rows["ships_independent_on_all_named"]) or "none")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())