"""How well do the published reference values of the six hulls hold up, cell by cell?

``crossref`` answers the two questions the library was built for -- *does a second publisher
state this quantity at all*, and *do the publishers that do agree?*  This module answers the
two questions a reader of the benchmark report asks next:

1. **Transcription fidelity.**  A pilot case definition that a later paper reprints should
   come back digit for digit.  Where the corpus prints the same particular more than once,
   how large is the spread?  A cell that only ever had one lineage is not made trustworthy
   by three copies of it; it is made trustworthy by three copies that *agree*, and the size
   of that agreement is a measurable number.
2. **Determinability.**  For each (hull, quantity) cell, classify what the retrievable
   corpus actually determines: nothing, one document, one lineage through many documents
   (transcription), or two or more independent lineages.  A quantity that the corpus
   nowhere prints is a *source* gap; a quantity the corpus prints but the parser drops is a
   *reader* gap.  The two need different work, so :func:`keyword_scan` separates them by
   scanning the fetched bytes for the quantity's own tokens.

Both are pure functions of the library JSON plus the fetched documents, so the numbers in
``docs/fleet_refs_fidelity.md`` are regenerated, never transcribed.
"""

from __future__ import annotations

import re
from typing import Any, Iterable, Mapping

from .crossref import QUANTITIES, TOLERANCE

#: Tokens a publisher uses when it states each quantity, independent of the row-label patterns
#: the parser applies.  Deliberately broad: the scan's job is to find *candidate* mentions and
#: let a human/§parser review decide, so a false positive costs a look while a false negative
#: hides a parser gap.
QUANTITY_TOKENS: dict[str, tuple[str, ...]] = {
    "displacement_volume": (r"displacement", r"volume displacement", r"∇", r"▽"),
    "cb": (r"\bC\s?B\b", r"block coefficient", r"方形係数"),
    "lcb_percent": (r"\bLCB\b", r"longitudinal cent(?:er|re) of buoyancy", r"浮心"),
    "lwl": (r"\bL\s?wl\b", r"length (?:at|on) (?:the )?water\s?line", r"水線長"),
    "length_pp": (r"\bL\s?pp\b", r"length between perpendiculars", r"垂線間長"),
    "beam": (r"\bB\s?wl\b", r"\bbeam\b", r"\bbreadth\b", r"幅"),
    "draft": (r"\bdraught\b", r"\bdraft\b", r"喫水"),
}

#: Cell classes, ordered from least to most determinate.  ``MULTI_LINEAGE_DISPUTED`` is not
#: "worse" than a single lineage -- two publishers disagreeing is information -- but it is not
#: a reference value either, and the anchor counts it as still open.
CLASS_NO_DETERMINATION = "NO_DETERMINATION"
CLASS_SINGLE_DOCUMENT = "SINGLE_LINEAGE_SINGLE_DOCUMENT"
CLASS_TRANSCRIPTION = "SINGLE_LINEAGE_TRANSCRIPTION"
CLASS_MULTI_AGREE = "MULTI_LINEAGE_AGREE"
CLASS_MULTI_DISPUTED = "MULTI_LINEAGE_DISPUTED"

CLASSES: tuple[str, ...] = (
    CLASS_NO_DETERMINATION,
    CLASS_SINGLE_DOCUMENT,
    CLASS_TRANSCRIPTION,
    CLASS_MULTI_AGREE,
    CLASS_MULTI_DISPUTED,
)

#: Cells in these classes still block ">= 2 independent sources" for their hull.
BLOCKING_CLASSES: frozenset[str] = frozenset(
    {CLASS_NO_DETERMINATION, CLASS_SINGLE_DOCUMENT, CLASS_TRANSCRIPTION}
)


def _voting(entries: Iterable[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    """Entries allowed to vote (a fallback-column reading has no known scale).

    ``cross_reference`` keeps such a reading in ``entries`` for the audit trail and denies it a
    vote; the library JSON does not serialise the vote list, so the same rule is re-applied here
    from the recorded ``unit_basis``.  Both spellings are honoured so the function works on
    in-memory reports as well as on the written artefact.
    """
    out = []
    for entry in entries:
        if "votes" in entry:
            if entry["votes"]:
                out.append(entry)
        elif entry.get("cell", {}).get("unit_basis") != "fallback_column":
            out.append(entry)
    return out


def _rel(a: float, b: float) -> float:
    scale = max(abs(a), abs(b))
    return 0.0 if scale == 0 else abs(a - b) / scale


def cluster(values: Iterable[float], tolerance: float) -> list[list[float]]:
    """Group values that agree to ``tolerance``, largest cluster first.

    Single-linkage on purpose: three publishers printing 0.6505 / 0.651 / 0.651 must land in
    one cluster, and a typo of 10x must land in its own.
    """
    ordered = sorted(values)
    groups: list[list[float]] = []
    for value in ordered:
        if groups and _rel(value, groups[-1][-1]) <= tolerance:
            groups[-1].append(value)
        else:
            groups.append([value])
    groups.sort(key=len, reverse=True)
    return groups


def transcription_fidelity(
    library: Mapping[str, Any],
    tolerance: Mapping[str, float] = TOLERANCE,
) -> dict[str, Any]:
    """Per (hull, quantity): how far apart the corpus' statements of that number are.

    ``spread_relative`` is ``(max - min) / min`` over the *voting* published values of one
    cell, i.e. across publishers *and* across lineages.  ``outliers`` are the entries that
    disagree with the largest agreeing cluster by more than the quantity's tolerance; each
    carries the printed ``raw`` text so a typo can be told from a different definition.
    """
    report = library["report"]
    cells: dict[str, Any] = {}
    for ship in sorted(report):
        per_quantity: dict[str, Any] = {}
        for quantity in QUANTITIES:
            cell = report[ship]["quantities"].get(quantity, {})
            entries = _voting(cell.get("entries", []))
            values = [float(entry["cell"]["value"]) for entry in entries]
            tol = tolerance[quantity]
            groups = cluster(values, tol)
            modal = groups[0] if groups else []
            outliers = []
            seen: set[tuple[Any, ...]] = set()
            for group in groups[1:]:
                for value in group:
                    for entry in entries:
                        if float(entry["cell"]["value"]) != value:
                            continue
                        # One row per (source, value, raw): a document that prints the same
                        # outlier on two rows is one finding, not two.
                        key = (entry["cell"].get("source"), value, entry["cell"].get("raw"))
                        if key in seen:
                            continue
                        seen.add(key)
                        outliers.append(
                            {
                                "source": entry["cell"].get("source"),
                                "family": entry.get("family"),
                                "value": entry["cell"]["value"],
                                "raw": entry["cell"].get("raw"),
                                "label": entry["cell"].get("label"),
                                "rel_to_modal": round(_rel(value, modal[0]), 6) if modal else None,
                            }
                        )
            per_quantity[quantity] = {
                "documents": len(entries),
                "lineages": cell.get("independent_families", 0),
                "verdict": cell.get("verdict"),
                "distinct_values": sorted({round(v, 6) for v in values}),
                "spread_relative": round(
                    (max(values) - min(values)) / min(values), 6
                ) if len(values) > 1 and min(values) else 0.0,
                "modal_value": modal[0] if modal else None,
                "modal_cluster_size": len(modal),
                "outliers": outliers,
                "tolerance": tol,
            }
        cells[ship] = per_quantity
    return {
        "schema": "tensorcad.fleet_refs.fidelity/1.0",
        "anchor": library.get("anchor", "REPRO-GEOM-FLEET-REFS"),
        "tolerance_basis": "tensorcad.fleet.crossref.TOLERANCE",
        "cells": cells,
    }


def determinability(
    library: Mapping[str, Any],
    fidelity: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Classify each cell by *what the retrievable corpus determines*, not by who agrees."""
    fidelity = fidelity or transcription_fidelity(library)
    report = library["report"]
    ships: dict[str, Any] = {}
    for ship in sorted(report):
        per_quantity: dict[str, Any] = {}
        for quantity in QUANTITIES:
            cell = report[ship]["quantities"].get(quantity, {})
            stats = fidelity["cells"][ship][quantity]
            lineages = stats["lineages"]
            documents = stats["documents"]
            if lineages >= 2:
                klass = (
                    CLASS_MULTI_AGREE
                    if cell.get("verdict") == "CROSS_VERIFIED"
                    else CLASS_MULTI_DISPUTED
                )
            elif lineages == 1:
                klass = CLASS_TRANSCRIPTION if documents > 1 else CLASS_SINGLE_DOCUMENT
            else:
                klass = CLASS_NO_DETERMINATION
            per_quantity[quantity] = {
                "class": klass,
                "documents": documents,
                "lineages": lineages,
                "spread_relative": stats["spread_relative"],
                "blocking": klass in BLOCKING_CLASSES,
            }
        counts = {klass: 0 for klass in CLASSES}
        for row in per_quantity.values():
            counts[row["class"]] += 1
        ships[ship] = {
            "quantities": per_quantity,
            "counts": counts,
            "blocking": [q for q, row in per_quantity.items() if row["blocking"]],
            #: The anchor's wording is "every hull: >= 2 independent sources" on four named
            #: quantities; the core four are reported separately from the auxiliary three.
            "core_four_resolved": all(
                not per_quantity[q]["blocking"]
                for q in ("displacement_volume", "cb", "lcb_percent", "lwl")
            ),
        }
    return {
        "schema": "tensorcad.fleet_refs.determinability/1.0",
        "anchor": report and library.get("anchor", "REPRO-GEOM-FLEET-REFS"),
        "core_quantities": ["displacement_volume", "cb", "lcb_percent", "lwl"],
        "ships": ships,
        "ships_core_four_resolved": sorted(
            ship for ship, row in ships.items() if row["core_four_resolved"]
        ),
        "classes": list(CLASSES),
        "blocking_classes": sorted(BLOCKING_CLASSES),
    }


def keyword_scan(
    documents: Mapping[str, str],
    extracted: Mapping[str, Mapping[str, Iterable[str]]] | None = None,
) -> dict[str, Any]:
    """Does a fetched document *mention* a quantity, and did the parser extract it?

    ``documents`` maps source key -> full text.  ``extracted`` maps source key -> ship ->
    quantities the parser read out of it.  A mention with no extraction is either a parser
    gap (the document prints the number and the reader missed it) or a mention in prose
    (the number is not stated as a particular) -- both are worth naming, only the first is a
    defect, and the report keeps the evidence so the difference can be settled by reading.
    """
    extracted = extracted or {}
    per_source: dict[str, Any] = {}
    gaps: list[dict[str, Any]] = []
    for key in sorted(documents):
        text = documents[key] or ""
        read: set[str] = set()
        for quantities in (extracted.get(key) or {}).values():
            read.update(quantities)
        rows: dict[str, Any] = {}
        for quantity in QUANTITIES:
            hits = 0
            for token in QUANTITY_TOKENS[quantity]:
                hits += len(re.findall(token, text, flags=re.I))
            rows[quantity] = {"mentions": hits, "extracted": quantity in read}
            if hits and quantity not in read:
                gaps.append(
                    {
                        "source": key,
                        "quantity": quantity,
                        "mentions": hits,
                        "kind": "mention_without_extraction",
                    }
                )
        per_source[key] = rows
    return {
        "schema": "tensorcad.fleet_refs.keyword_scan/1.0",
        "sources": per_source,
        "mention_without_extraction": gaps,
        "documents_scanned": len(per_source),
    }


def summarise(
    fidelity: Mapping[str, Any],
    determinability_report: Mapping[str, Any],
    scan: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """The three numbers a benchmark report should quote about its own reference set."""
    cells = [
        row
        for ship in fidelity["cells"]
        for row in fidelity["cells"][ship].values()
    ]
    reproduced = [row for row in cells if row["documents"] > 1]
    worst = max(reproduced, key=lambda r: r["spread_relative"], default=None)
    outliers = [
        {"ship": ship, "quantity": quantity, **outlier}
        for ship in fidelity["cells"]
        for quantity, row in fidelity["cells"][ship].items()
        for outlier in row["outliers"]
    ]
    return {
        "cells_total": len(cells),
        "cells_reproduced_by_two_or_more_documents": len(reproduced),
        "cells_with_outlier": sum(1 for row in cells if row["outliers"]),
        "outliers": outliers,
        "largest_transcription_spread": None
        if worst is None
        else {
            "spread_relative": worst["spread_relative"],
            "documents": worst["documents"],
            "tolerance": worst["tolerance"],
        },
        "ships_core_four_resolved": determinability_report["ships_core_four_resolved"],
        "ships": sorted(determinability_report["ships"]),
        "mention_without_extraction": 0
        if scan is None
        else len(scan["mention_without_extraction"]),
    }