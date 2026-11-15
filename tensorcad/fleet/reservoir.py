"""What does the *already fetched* corpus state that the reference library does not?

``fidelity`` classifies every (hull, quantity) cell by what the retrievable corpus determines and
lists the cells that still block the anchor's ">= 2 independent sources": five hulls, and for each
of them a handful of quantities that no admitted source states.  The next question is the one the
next worker asks: *is that a retrieval gap or a reader gap?* -- i.e. has anybody even looked at the
bytes that were already downloaded for those cells, or would a new network hunt be a repeat of work
already done?

That question has a cheap, offline answer.  Earlier rounds fetched far more documents than they
admitted: hunts for a second JBC source, the DTC originator paper, Series 60 tables and KVLCC2
editions left PDFs in scratch caches that were never parsed because the hunt ended when the
*admitted* set was still young.  This module scans that **reservoir** -- every fetched document,
admitted or not -- for the exact cells the library is missing, and reports:

1. the roster: each document's SHA-256, size, whether it is already a library source, which hulls
   it mentions at all;
2. per blocking cell, every *numeral-bearing* mention of the quantity in a document that mentions
   that hull, quoted with its line -- a *candidate* for a reader, not an admitted value;
3. a verdict per cell: ``CANDIDATE_FIELD`` when some reservoir document appears to state it (the
   cheap next increment is a parser/citation job, not a network hunt), ``SOURCE_GAP_CONFIRMED``
   when the reservoir mentions the hull but never states the quantity, ``NOT_IN_RESERVOIR`` when no
   reservoir document even mentions the hull.

Nothing here admits a value.  A candidate is a lead with a pinned digest and a quoted line; turning
it into a library cell stays a deliberate act (see ``crossref.SourceSpec``), because a scan that
ingests on its own would make the library's independence claims unfalsifiable.

Both layers are pure functions of (document bytes, library JSON) so ``docs/fleet_refs_reservoir_sweep.md``
is regenerated, never transcribed.
"""

from __future__ import annotations

import hashlib
import re
from typing import Any, Iterable, Mapping, Sequence

#: How a hull is named in the literature.  Deliberately broad -- a false positive costs one line of
#: review, a false negative would hide a document that does state a missing cell.
HULL_PATTERNS: dict[str, str] = {
    "KCS": r"KRISO\s*3600|\(KCS\)|\bKCS\b",
    "KVLCC2": r"\bKVLCC\s*2\b",
    "JBC": r"Japan\s+Bulk\s+Carrier|\bJBC\b|SR\s*266",
    "ONRT": r"ONR\s*Tumblehome|\bONRT\b|DTMB\s*5613|5613",
    "DTC": r"Duisburg\s+Test\s+Case|\(DTC\)|\bDTC\b",
    "Series60": r"Series\s*60\b",
}

#: A reservoir hit must carry a numeral: "LCB" in a sentence about longitudinal centre of buoyancy
#: is prose, "LCB (%) fwd + 3.48" is a particular.  The patterns capture the numeral so the report
#: can quote it, and the number is *never* written into the library from here.
QUANTITY_PATTERNS: dict[str, str] = {
    "displacement_volume": r"(?:displacement|volume)[^\n]{0,40}?([0-9][0-9 ,.]{2,})",
    "cb": r"(?:block\s+coefficient|\bC\s?B\b|\bCb\b)[^\n]{0,20}?([01]?[.,][0-9]{2,4})",
    #: LCB also appears with mangled subscripts in PDF text extraction ("LCB(%)", "L") and as the
    #: right-hand column of a "fwd +" row; both are scanned, neither is trusted without review.
    "lcb_percent": r"(?:\bLCB\b|\bLCB\s*\(|fwd\s*\+)[^\n]{0,40}?(-?[0-9]+[.,][0-9]+)",
    "lwl": r"(?:L\s?wl|waterline\s+length|length\s+(?:at|on)\s+(?:the\s+)?water\s?line)[^\n]{0,20}?([0-9]+[.,]?[0-9]*)",
    "length_pp": r"(?:L\s?pp|LPP|length\s+between\s+perpendiculars)[^\n]{0,20}?([0-9]+[.,]?[0-9]*)",
    "beam": r"(?:B\s?wl|beam|breadth)[^\n]{0,20}?([0-9]+[.,]?[0-9]*)",
    "draft": r"(?:draft|draught)[^\n]{0,20}?([0-9]+[.,]?[0-9]*)",
}

#: Cell verdicts.  ``CANDIDATE_FIELD`` is the only one that implies work here rather than a network.
VERDICT_CANDIDATE = "CANDIDATE_FIELD"
VERDICT_SOURCE_GAP = "SOURCE_GAP_CONFIRMED"
VERDICT_ABSENT = "NOT_IN_RESERVOIR"

#: Lines quoted per (document, quantity) in the artefact -- enough to see the row, not a transcript.
MAX_QUOTES = 4

#: Plausibility window per quantity, in SI / dimensionless units.  A scan that reports "LCB 2021. 8"
#: is worse than no scan: it costs a review to reject and teaches the next reader to ignore the
#: table.  The windows are wide (model scale included) but they are the reason the reported leads
#: are readable.
PLAUSIBLE: dict[str, tuple[float, float]] = {
    "displacement_volume": (0.002, 2.0e6),
    "cb": (0.15, 0.95),
    "lcb_percent": (-20.0, 20.0),
    "lwl": (0.5, 600.0),
    "length_pp": (0.5, 600.0),
    "beam": (0.05, 80.0),
    "draft": (0.01, 40.0),
}

#: A row in a reference list is a citation, not a statement of a particular.  Filtering these is
#: what keeps "Series 60 C_B=0.60 Model, Proceedings ..." out of the Series 60 lead list.
_BIBLIOGRAPHY = re.compile(r"proceedings|vol\.|\bpp\.|doi|journal|conference|isbn|et al", re.I)


def text_digest(text: str) -> str:
    """Digest of the *content*, used to collapse a PDF and its pre-extracted ``.txt`` sidecar."""
    normalised = " ".join(text.split())
    return hashlib.sha256(normalised.encode("utf-8", "replace")).hexdigest()


def printable_ratio(text: str) -> float:
    """Fraction of characters that are printable -- a binary blob is not a document."""
    if not text:
        return 0.0
    good = sum(1 for char in text if char.isprintable() or char in "\n\t")
    return good / len(text)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def hulls_mentioned(text: str) -> dict[str, str]:
    """Which hulls the document names, each with one quoted context line."""
    found: dict[str, str] = {}
    for hull, pattern in HULL_PATTERNS.items():
        match = re.search(pattern, text, flags=re.I)
        if match:
            start = max(0, match.start() - 80)
            found[hull] = " ".join(text[start : match.end() + 100].split())[:220]
    return found


def _plausible(quantity: str, numeral: str) -> bool:
    window = PLAUSIBLE.get(quantity)
    if window is None:
        return True
    try:
        value = float(numeral.replace(" ", "").replace(",", "."))
    except ValueError:
        return False
    return window[0] <= value <= window[1]


def quantity_mentions(text: str, quantities: Sequence[str] | None = None) -> dict[str, list[dict[str, Any]]]:
    """Numeral-bearing mentions per quantity, with the quoted line and the captured numeral."""
    out: dict[str, list[dict[str, Any]]] = {}
    for quantity, pattern in QUANTITY_PATTERNS.items():
        if quantities is not None and quantity not in quantities:
            continue
        hits: list[dict[str, Any]] = []
        lines = text.splitlines()
        for index, line in enumerate(lines):
            if _BIBLIOGRAPHY.search(line):
                continue
            # PDF extraction of a two-column particulars table drops the numeral onto the next line
            # ("LCB (%), fwd +" / "3.48"), so LCB is also matched over a two-line window.  The other
            # quantities are printed inline in every document seen so far.
            candidates = [line]
            if quantity == "lcb_percent" and index + 1 < len(lines):
                candidates.append(line + " " + lines[index + 1])
            for window in candidates:
                if _BIBLIOGRAPHY.search(window):
                    continue
                for match in re.finditer(pattern, window, flags=re.I):
                    if not _plausible(quantity, match.group(1)):
                        continue
                    record = {
                        "line": index + 1,
                        "numeral": match.group(1),
                        "quote": " ".join(window.split())[:200],
                    }
                    if record not in hits:
                        hits.append(record)
                    if len(hits) >= MAX_QUOTES:
                        break
                if len(hits) >= MAX_QUOTES:
                    break
            if len(hits) >= MAX_QUOTES:
                break
        if hits:
            out[quantity] = hits
    return out


def blocking_cells(library: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Every (hull, quantity) the fidelity report leaves blocking, from the artefact itself.

    A cell is blocking when the library holds no value (``MISSING``) or when its determinability
    class is one of :data:`tensorcad.fleet.fidelity.BLOCKING_CLASSES`.  Both sources are read, so a
    cell that gains a source between runs stops being blocking without touching this module.
    """
    from .fidelity import BLOCKING_CLASSES

    cells: list[dict[str, Any]] = []
    report = library.get("report", {})
    determinability = library.get("determinability", {})
    flat_cells = determinability.get("cells") or {}
    ship_cells = determinability.get("ships") or {}
    for hull in sorted(report):
        for quantity in sorted(report[hull].get("quantities", {})):
            core_verdict = report[hull]["quantities"][quantity].get("verdict", "")
            if flat_cells:
                cell_class = flat_cells.get(f"{hull}/{quantity}", {}).get("class")
            else:
                cell_class = (
                    ship_cells.get(hull, {}).get("quantities", {}).get(quantity, {}).get("class")
                )
            if core_verdict != "MISSING" and cell_class not in BLOCKING_CLASSES:
                continue
            cells.append(
                {
                    "hull": hull,
                    "quantity": quantity,
                    "library_verdict": core_verdict,
                    "class": cell_class,
                }
            )
    return cells


def scan_document(
    text: str,
    *,
    name: str,
    sha256: str,
    bytes_: int,
    in_library: bool,
    quantities: Sequence[str] | None = None,
) -> dict[str, Any]:
    """One reservoir document's contribution: hulls named, numerals found, digest pinned."""
    hulls = hulls_mentioned(text)
    return {
        "name": name,
        "sha256": sha256,
        "bytes": bytes_,
        "in_library": in_library,
        "chars": len(text),
        "hulls": hulls,
        "mentions": quantity_mentions(text, quantities) if hulls else {},
    }


def classify(cells: Iterable[Mapping[str, Any]], documents: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Per blocking cell: does any reservoir document state it, and which ones mention the hull?"""
    out: list[dict[str, Any]] = []
    for cell in cells:
        hull, quantity = cell["hull"], cell["quantity"]
        mentioners = [d for d in documents if hull in d["hulls"]]
        hits = [
            {
                "document": d["name"],
                "sha256": d["sha256"],
                "in_library": d["in_library"],
                "mentions": d["mentions"][quantity],
            }
            for d in mentioners
            if quantity in d["mentions"]
        ]
        if hits:
            verdict = VERDICT_CANDIDATE
        elif mentioners:
            verdict = VERDICT_SOURCE_GAP
        else:
            verdict = VERDICT_ABSENT
        out.append(
            {
                **dict(cell),
                "verdict": verdict,
                "documents_mentioning_hull": len(mentioners),
                "documents_stating_quantity": len(hits),
                "candidates": hits,
            }
        )
    return out


def summarise(documents: Sequence[Mapping[str, Any]], cells: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Counts the report quotes; every one of them is derived, none typed in."""
    by_verdict: dict[str, int] = {}
    for cell in cells:
        by_verdict[cell["verdict"]] = by_verdict.get(cell["verdict"], 0) + 1
    by_hull: dict[str, dict[str, int]] = {}
    for cell in cells:
        bucket = by_hull.setdefault(cell["hull"], {})
        bucket[cell["verdict"]] = bucket.get(cell["verdict"], 0) + 1
    return {
        "documents_scanned": len(documents),
        "documents_in_library": sum(1 for d in documents if d["in_library"]),
        "documents_outside_library": sum(1 for d in documents if not d["in_library"]),
        "blocking_cells": len(cells),
        "by_verdict": by_verdict,
        "by_hull": by_hull,
        "candidate_cells": [
            f"{c['hull']}/{c['quantity']}" for c in cells if c["verdict"] == VERDICT_CANDIDATE
        ],
    }