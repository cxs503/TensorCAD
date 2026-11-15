#!/usr/bin/env python3
"""Ingest the second-round open-access documents into the fleet reference library.

    python tools/baseline_fleet/crossref_oa.py ingest --raw /path/to/pdfs
    python tools/baseline_fleet/crossref_oa.py build

Why a second tool
-----------------
``crossref_refs.py`` handles HTML tables: it addresses a value by *cell*.  The documents this
round could actually be obtained are PDFs (journal articles and a thesis on repository hosts,
because the publishers' own pages answer 403/404/RemoteDisconnected from this egress), and a
PDF has no cells -- the address is a *verbatim quote*.  ``tensorcad.fleet.crossref.parse_pdf_text``
applies those quotes; this tool fetches nothing on ``build``, extracts the text, pins each
document's SHA-256, and merges the new cells into the existing library rather than rebuilding
it, because the HTML pages of the first round are no longer cached in the worktree
(``pointer_status: cache-lost``) and a rebuild would silently drop KCS/ONRT coverage.

Independence is reported on two axes, which is the honest reading of "at least two sources":

* ``families``   -- distinct publishers holding the numeral (document-level independence);
* ``lineages``   -- distinct *origins*, after collapsing documents that state in their own text
  that they are reproducing another body's table (e.g. a 2023 journal paper that prints the
  Tokyo-2015 JBC table is one lineage with the workshop page, not a second opinion).

A cell with two publishers but one lineage is reported ``REPRODUCED_SINGLE_LINEAGE``: it is
corroboration of transcription, never evidence of a second measurement.
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Mapping

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.fleet.crossref import (  # noqa: E402
    QUANTITIES,
    QUANTITY_LABELS,
    SOURCES,
    Cell,
    closure_check,
    coverage_summary,
    cross_reference,
    lineage_of,
    parse_pdf_text,
)

#: ``source key -> file name`` of the locally obtained document (kept out of the repository:
#: none of these documents states a redistribution licence, only digest + quote are committed).
RAW_FILES: dict[str, str] = {
    "kvlcc2_ugent": "kvlcc2_ugent.pdf",
    "dtc_ugent": "kvlcc2_ugent.pdf",
    "kvlcc2_hrcak": "kvlcc2_hrcak.pdf",
    "kvlcc2_jstage": "kvlcc2_jstage.pdf",
    "jbc_cfdletters": "jbc_semarak.pdf",
    "jbc_scse": "jbc_scienpg.pdf",
    "jbc_memorial": "jbc_kabonge.pdf",
    "dtc_strath": "dtc_strath.pdf",
    # third round (m21345): documents fetched after the second round was merged
    "dtc_ije": "dtc_ije.pdf",
    "kvlcc2_osaka": "kvlcc2_osaka.pdf",
}

OA_KEYS = tuple(RAW_FILES)

#: Verdict tiers.  ``CROSS_VERIFIED`` from the core module means "two publishers agree";
#: the rewrite below tightens it with the lineage axis.
TIER_LINEAGE = "CROSS_VERIFIED_LINEAGE"
TIER_REPRODUCED = "REPRODUCED_SINGLE_LINEAGE"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def pdf_text(path: Path) -> tuple[str, int]:
    """Extracted text with ``\\f`` page breaks, plus the page count (PyMuPDF, optional)."""
    try:
        import fitz  # type: ignore[import-not-found]
    except ImportError:  # pragma: no cover - documented dependency of this tool only
        raise SystemExit("this tool needs PyMuPDF (pip install pymupdf); the core library does not")
    document = fitz.open(path)
    return "\n\f".join(page.get_text() for page in document), document.page_count


def cmd_ingest(args: argparse.Namespace) -> int:
    raw = Path(args.raw)
    out = Path(args.out)
    previous: dict[str, dict[str, Any]] = {}
    if out.exists():
        previous = json.loads(out.read_text()).get("sources", {})
    seen: dict[str, dict[str, Any]] = {}
    for key in OA_KEYS:
        path = raw / RAW_FILES[key]
        if not path.exists():
            # A document obtained in an earlier round is no longer cached in the worktree.
            # Re-using its last read row keeps coverage from silently shrinking; the row keeps
            # the digest it was read from, so the cell stays auditable.
            cached = previous.get(key)
            if args.keep_existing and cached is not None:
                seen[key] = dict(cached, note="reused from the previous ingest; local file absent")
            else:
                seen[key] = {"status": "missing_local_file", "file": str(path)}
            continue
        spec = next(s for s in SOURCES if s.key == key)
        digest = sha256(path)
        text, pages = pdf_text(path)
        cells = parse_pdf_text(spec, text)
        seen[key] = {
            "status": "read" if cells else "no_quote_matched",
            "file": path.name,
            "bytes": path.stat().st_size,
            "pages": pages,
            "sha256": digest,
            "pinned": spec.sha256 or None,
            "digest_matches_pin": (spec.sha256 == digest) if spec.sha256 else None,
            "quantities_read": sorted(cells),
            "cells": {q: c.as_dict() for q, c in cells.items()},
        }
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "tensorcad.fleet_refs.oa_ingest/1.0",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "sources": seen,
    }
    out.write_text(json.dumps(payload, indent=1) + "\n")
    read = sum(1 for row in seen.values() if row["status"] == "read")
    print(f"wrote {out}: {read}/{len(OA_KEYS)} documents read")
    for key, row in seen.items():
        print(f"  {key:16s} {row['status']:16s} {row.get('quantities_read', '')}")
    return 0


CELL_FIELDS = tuple(f.name for f in dataclasses.fields(Cell))


def cell_from(data: Mapping[str, Any]) -> Cell:
    """``Cell`` from a stored dict, ignoring keys added downstream (adjudication, notes)."""
    return Cell(**{name: data[name] for name in CELL_FIELDS if name in data})


def cells_from_report(report: Mapping[str, Any]) -> dict[str, dict[str, dict[str, Cell]]]:
    """Rebuild the first round's cells from the committed library (no re-fetch required)."""
    extractions: dict[str, dict[str, dict[str, Cell]]] = {}
    for ship, row in report.items():
        for quantity, cell in row["quantities"].items():
            for entry in cell["entries"]:
                if entry["cell"].get("unit_basis") == "fallback_column":
                    continue
                data = dict(entry["cell"])
                # ``source`` both keys the extraction and is a required ``Cell`` field: it is
                # read, never dropped (popping it made ``build`` raise TypeError on the first
                # committed cell, so the merge path had never run).
                key = data.get("source") or "unattributed"
                data["source"] = key
                extractions.setdefault(ship, {}).setdefault(key, {})[quantity] = cell_from(data)
    return extractions


def lineage_tally(report: Mapping[str, Any]) -> dict[str, Any]:
    """Per ship: publisher families vs independent lineages behind the numerals."""
    out: dict[str, Any] = {}
    for ship, row in report.items():
        per_quantity: dict[str, Any] = {}
        for quantity in QUANTITIES:
            cell = row["quantities"][quantity]
            families = sorted({entry["family"] for entry in cell["entries"]})
            lineages = sorted({lineage_of(entry["cell"]["source"]) for entry in cell["entries"]})
            values = sorted({round(entry["cell"]["value"], 4) for entry in cell["entries"] if entry["cell"]["unit_basis"] == "volume"})
            verdict = cell["verdict"]
            if verdict == "CROSS_VERIFIED" and len(lineages) < 2:
                verdict = TIER_REPRODUCED
            elif verdict == "CROSS_VERIFIED":
                verdict = TIER_LINEAGE
            per_quantity[quantity] = {
                "verdict": verdict,
                "core_verdict": cell["verdict"],
                "families": families,
                "lineages": lineages,
                "sources": sorted({entry["cell"]["source"] for entry in cell["entries"]}),
                "published_values": cell["published_values"],
                "distinct_values": values,
                "disputed": cell["disputed"],
            }
        out[ship] = per_quantity
    return out


def merge_into(payload: dict[str, Any], ingest: Mapping[str, Any], fleet_specs: Path | None) -> dict[str, Any]:
    """Add the OA cells to the committed library and recompute report/coverage/lineage."""
    extractions = cells_from_report(payload["report"])
    added = 0
    for key, row in ingest["sources"].items():
        if row["status"] != "read":
            continue
        spec = next(s for s in SOURCES if s.key == key)
        digest = row.get("sha256") or spec.sha256
        for quantity, cell in row["cells"].items():
            data = dict(cell)
            # A cell quoted out of a PDF is only auditable with the digest of the bytes it was
            # read from; the ingest report holds that digest, the cell itself may predate it.
            if not data.get("sha256"):
                data["sha256"] = digest
            extractions.setdefault(spec.ship, {}).setdefault(key, {})[quantity] = Cell(**data)
            added += 1
    report = cross_reference(extractions)
    # A ship whose every quantity is MISSING contributes no extraction, so rebuilding from
    # cells alone would silently drop it from the report (Series 60 has no readable source at
    # all).  Keep such ships with empty rows: "no source found" is a result the coverage table
    # must keep scoring, not a row that disappears.
    for ship in payload.get("report", {}):
        if ship not in report:
            report.update(cross_reference({ship: {}}))
    payload["report"] = report
    payload["coverage"] = coverage_summary(report)
    payload["lineage"] = lineage_tally(report)
    payload["oa_sources"] = {
        key: {k: row[k] for k in ("status", "file", "bytes", "pages", "sha256", "quantities_read") if k in row}
        for key, row in ingest["sources"].items()
    }
    payload["oa_cells_added"] = added
    payload["closure_oa"] = closure_check(report)
    if fleet_specs is not None:
        payload["fleet_backfill"] = backfill(fleet_specs, payload["lineage"], report)
    return payload


def backfill(path: Path, lineage: Mapping[str, Any], report: Mapping[str, Any]) -> dict[str, Any]:
    """Write the multi-source arrays into ``fleet_specs.json`` (additive: old keys untouched)."""
    specs = json.loads(path.read_text())
    summary: dict[str, Any] = {}
    for record in specs["ships"]:
        ship = record["name"]
        if ship not in report:
            summary[ship] = {"status": "no_cross_reference_rows"}
            continue
        refs: dict[str, Any] = {}
        for quantity in QUANTITIES:
            cell = report[ship]["quantities"][quantity]
            if cell["verdict"] == "MISSING":
                continue
            refs[quantity] = [
                {
                    "source": entry["cell"]["source"],
                    "family": entry["family"],
                    "lineage": lineage_of(entry["cell"]["source"]),
                    "value": entry["cell"]["value"],
                    "unit": "m3" if quantity == "displacement_volume" else "-",
                    "url": entry["cell"]["url"],
                    "sha256": entry["cell"]["sha256"],
                    "quote": entry["cell"]["raw"],
                }
                for entry in cell["entries"]
                if entry["cell"]["unit_basis"] != "fallback_column"
            ]
            summary.setdefault(ship, {})[quantity] = {
                "verdict": lineage[ship][quantity]["verdict"],
                "families": lineage[ship][quantity]["families"],
                "lineages": lineage[ship][quantity]["lineages"],
                "disputed": cell["disputed"],
            }
        record["published_cross_refs"] = refs
        record["published_cross_refs_status"] = {
            "quantities": summary.get(ship, {}),
            "library": "benchmarks/verified/fleet_refs/crossref.json",
        }
    path.write_text(json.dumps(specs, indent=2) + "\n")
    return summary


def cmd_build(args: argparse.Namespace) -> int:
    library = Path(args.library)
    ingest_path = Path(args.ingest)
    payload = json.loads(library.read_text())
    ingest = json.loads(ingest_path.read_text())
    fleet = Path(args.fleet_specs) if args.fleet_specs else None
    if fleet is not None and not fleet.exists():
        fleet = None
    merged = merge_into(payload, ingest, fleet)
    merged["generated_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    merged["quantity_labels"] = QUANTITY_LABELS
    library.write_text(json.dumps(merged, indent=1) + "\n")
    cov = merged["coverage"]
    print(f"wrote {library}: cells {cov['cells']}, cross_verified {cov['cross_verified']}, "
          f"conflict {cov['conflict']}, missing {cov['missing']}, oa cells {merged['oa_cells_added']}")
    for ship in cov["ships"]:
        row = merged["lineage"][ship]
        tiers = {q: row[q]["verdict"] for q in QUANTITIES}
        print(f"  {ship:8s} " + " ".join(f"{q[:6]}={v[:12]}" for q, v in tiers.items()))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p_ingest = sub.add_parser("ingest", help="hash + quote-extract the local PDFs")
    p_ingest.add_argument("--raw", default="/tmp/oa")
    p_ingest.add_argument("--out", default="benchmarks/verified/fleet_refs/oa_ingest.json")
    p_ingest.add_argument("--keep-existing", action="store_true",
                          help="keep the rows of documents whose local file is no longer cached")
    p_ingest.set_defaults(func=cmd_ingest)
    p_build = sub.add_parser("build", help="merge the OA cells into the committed library")
    p_build.add_argument("--library", default="benchmarks/verified/fleet_refs/crossref.json")
    p_build.add_argument("--ingest", default="benchmarks/verified/fleet_refs/oa_ingest.json")
    p_build.add_argument("--fleet-specs", default="",
                     help="optional; the canonical specs refresh is tools/baseline_fleet/backfill_refs.py")
    p_build.set_defaults(func=cmd_build)
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())