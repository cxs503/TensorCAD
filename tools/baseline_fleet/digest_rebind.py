#!/usr/bin/env python3
"""Re-read the fleet reference library out of its own documents and audit every digit.

    python tools/baseline_fleet/digest_rebind.py audit --cache /tmp/xref/raw --oa-cache /tmp/oa \
        --out benchmarks/verified/fleet_refs/rebind_audit.json
    python tools/baseline_fleet/digest_rebind.py audit --fetch --out .../rebind_audit.json
    python tools/baseline_fleet/digest_rebind.py repair --out .../rebind_audit.json

``audit`` walks every pinned document of :data:`tensorcad.fleet.crossref.SOURCES`, re-derives
its cells from the cached bytes with the *current* reading rules, and compares them with the
digits the published library carries.  ``--fetch`` additionally re-downloads each document so
a page that moved since the pin is reported as drift instead of being assumed unchanged.

``repair`` rewrites the stale *derived* cells found by the audit (arithmetic on a document's
own numerals) and re-runs the cached derived blocks.  Transcribed numerals are never touched:
a quoting error has to be re-read by a human, not recomputed by a script.

Neither mode fetches anything by default: an audit that needs the network to say "verified"
cannot be re-run by a reviewer offline.
"""

from __future__ import annotations

import argparse
import datetime
import json
import ssl
import sys
import urllib.request
from pathlib import Path
from typing import Any

if __package__ in (None, ""):  # allow ``python tools/baseline_fleet/digest_rebind.py``
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.fleet.crossref import (  # noqa: E402
    SOURCES,
    parse_pdf_text,
    parse_source,
    refresh_derived,
    split_ships,
)
from tensorcad.fleet.rebind import (  # noqa: E402
    SCHEMA,
    apply_repair,
    compare_cells,
    digest_verdict,
    summarize,
)

UA = "Mozilla/5.0 (X11; Linux x86_64) tensorcad-fleet-refs/1.0"
LIBRARY = Path("benchmarks/verified/fleet_refs/crossref.json")
INGEST = Path("benchmarks/verified/fleet_refs/oa_ingest.json")


def pdf_text(path: Path) -> str:
    """Extracted text with ``\\f`` page breaks (PyMuPDF; the extractor the ingest used)."""
    try:
        import fitz  # type: ignore[import-not-found]
    except ImportError:  # pragma: no cover - tool-only dependency
        raise SystemExit("this tool needs PyMuPDF (pip install pymupdf); the core library does not")
    document = fitz.open(str(path))
    return "\n\f".join(page.get_text() for page in document)


def library_cells(payload: dict[str, Any]) -> dict[str, dict[str, dict[str, Any]]]:
    """``source key -> ship -> quantity -> cell`` exactly as the published library has it."""
    out: dict[str, dict[str, dict[str, Any]]] = {}
    for ship, block in payload.get("report", {}).items():
        for quantity, cell_block in (block.get("quantities") or {}).items():
            for row in cell_block.get("entries") or []:
                cell = row.get("cell") or {}
                key = cell.get("source")
                if not key:
                    continue
                out.setdefault(key, {}).setdefault(ship, {})[quantity] = cell
    return out


def ingest_cells(payload: dict[str, Any]) -> dict[str, dict[str, dict[str, Any]]]:
    """``source key -> ship -> quantity -> cell`` for the open-access documents."""
    out: dict[str, dict[str, dict[str, Any]]] = {}
    for key, row in (payload.get("sources") or {}).items():
        for quantity, cell in (row.get("cells") or {}).items():
            out.setdefault(key, {}).setdefault(str(cell.get("ship") or key), {})[quantity] = cell
    return out


def fresh_cells(spec, payload: bytes) -> dict[str, dict[str, Any]]:
    """``ship -> quantity -> Cell`` re-derived from the document bytes."""
    out: dict[str, dict[str, Any]] = {}
    if spec.kind == "html":
        for ship, block in split_ships(spec, payload):
            cells = parse_source(spec, block)
            if cells:
                out.setdefault(ship, {}).update(cells)
        return out
    cells = parse_pdf_text(spec, pdf_text_bytes(payload))
    if cells:
        out.setdefault(spec.ship, {}).update(cells)
    return out


def pdf_text_bytes(payload: bytes) -> str:
    """Text of an in-memory PDF (written to a temporary file for PyMuPDF)."""
    import tempfile

    import fitz  # type: ignore[import-not-found]

    with tempfile.NamedTemporaryFile(suffix=".pdf") as handle:
        handle.write(payload)
        handle.flush()
        document = fitz.open(handle.name)
        return "\n\f".join(page.get_text() for page in document)


def cached_bytes(spec, cache: Path, oa_cache: Path, ingest: dict[str, Any]) -> bytes | None:
    """The copy the library was read from, if this environment still holds it."""
    row = (ingest.get("sources") or {}).get(spec.key) or {}
    names = []
    if row.get("file"):
        names.append(oa_cache / str(row["file"]))
    names.append(cache / f"{spec.key}.{spec.kind}")
    for path in names:
        if path.exists():
            return path.read_bytes()
    return None


def fetch(url: str, timeout: float = 45.0) -> bytes | None:
    """A live copy of a document, or ``None`` when this egress cannot reach it."""
    context = ssl.create_default_context()
    context.check_hostname = False
    context.verify_mode = ssl.CERT_NONE
    try:
        request = urllib.request.Request(url, headers={"User-Agent": UA})
        return urllib.request.urlopen(request, timeout=timeout, context=context).read()
    except Exception:  # noqa: BLE001 - an unreachable document is a recorded outcome
        return None


def run_audit(args: argparse.Namespace) -> dict[str, Any]:
    cache = Path(args.cache)
    oa_cache = Path(args.oa_cache)
    library_path = Path(args.library)
    library = json.loads(library_path.read_text())
    ingest = json.loads(Path(args.ingest).read_text()) if Path(args.ingest).exists() else {}
    stored = library_cells(library)
    spec_by_key = {spec.key: spec for spec in SOURCES}
    for key, per_ship in ingest_cells(ingest).items():
        spec = spec_by_key.get(key)
        # A PDF cell carries no ship field of its own: the document belongs to one hull, and
        # that hull is what the audit must compare against (never the source key).
        ship = spec.ship if spec is not None else next(iter(per_ship))
        for cells in per_ship.values():
            stored.setdefault(key, {}).setdefault(ship, {}).update(cells)

    documents: list[dict[str, Any]] = []
    verdicts = []
    derived_block: list[dict[str, Any]] = []
    for spec in SOURCES:
        payload = cached_bytes(spec, cache, oa_cache, ingest)
        live = fetch(spec.url) if args.fetch else None
        row: dict[str, Any] = {
            "key": spec.key,
            "family": spec.family,
            "kind": spec.kind,
            "url": spec.url,
            "ship": spec.ship,
        }
        row.update(digest_verdict(spec.sha256, payload, live))
        row["cells_stored"] = sum(len(cells) for cells in stored.get(spec.key, {}).values())
        if payload is None:
            documents.append(row)
            continue
        fresh = fresh_cells(spec, payload)
        row["cells_fresh"] = sum(len(cells) for cells in fresh.values())
        row["ships_fresh"] = sorted(fresh)
        for ship in sorted(set(stored.get(spec.key, {})) | set(fresh)):
            verdicts.extend(
                compare_cells(spec.key, ship, stored.get(spec.key, {}).get(ship, {}), fresh.get(ship, {}))
            )
        documents.append(row)

    summary = summarize(documents, verdicts)
    report = {
        "schema": SCHEMA,
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "library": str(library_path),
        "ingest": str(args.ingest),
        "caches": {"documents": str(cache), "open_access": str(oa_cache)},
        "method": (
            "bytes -> sha256 -> text -> numeral (current reading rules) -> library value; "
            "a document that is not cached is NOT_CACHED and never counts as verified"
        ),
        "documents": documents,
        "cells": [verdict.as_dict() for verdict in verdicts],
        "counts": summary,
    }
    if args.repair:
        repaired = apply_repair(library, verdicts)
        if repaired:
            refreshed = refresh_derived(library)
            report["repairs"] = repaired
            report["refreshed_blocks"] = refreshed
            library_path.write_text(json.dumps(library, indent=1) + "\n")
        else:
            report["repairs"] = []
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1) + "\n")
    counts = summary
    print(
        f"{out}: documents {counts['documents_checked']}/{counts['documents']} read, "
        f"cells rebound {counts['cells_rebound']}/{counts['cells']}, "
        f"non-rebound {len(counts['non_rebound'])}, all_rebound={counts['all_cells_rebound']}"
    )
    for item in counts["non_rebound"]:
        print(f"  {item['status']:20s} {item['ship']:9s} {item['quantity']:20s} {item.get('detail','')}")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode",
        nargs="?",
        default="audit",
        choices=("audit", "repair"),
        help="audit = report only (default); repair = also rewrite stale derived cells",
    )
    parser.add_argument("--library", default=str(LIBRARY))
    parser.add_argument("--ingest", default=str(INGEST))
    parser.add_argument("--cache", default="/tmp/xref/raw", help="cached publisher bytes (round one)")
    parser.add_argument("--oa-cache", default="/tmp/oa", help="cached open-access PDFs")
    parser.add_argument(
        "--out", default="benchmarks/verified/fleet_refs/rebind_audit.json"
    )
    parser.add_argument("--fetch", action="store_true", help="also re-download every document")
    args = parser.parse_args(argv)
    args.repair = args.mode == "repair"
    run_audit(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())