#!/usr/bin/env python3
"""Scan the fetched-but-unadmitted corpus for the cells the reference library still lacks.

    python tools/baseline_fleet/refs_reservoir_sweep.py \
        --library benchmarks/verified/fleet_refs/crossref.json \
        --out benchmarks/verified/fleet_refs/reservoir_sweep.json \
        --reservoir artifacts/fleet_refs/raw --reservoir /tmp/oa --reservoir /tmp/hunt --reservoir /tmp/refs

Why this tool exists: the library's blocking cells (``tensorcad.fleet.reservoir.blocking_cells``)
are the cells the anchor REPRO-GEOM-FLEET-REFS still owes.  Before anybody runs a new network hunt
for them, the cheaper question is whether a *previous* hunt already downloaded a document that
states them and simply never parsed it.  The reservoir is every PDF/HTML/text file in the given
roots; documents whose SHA-256 is already a library source are marked as such, the rest are
scanned for hull names and numeral-bearing quantity mentions.

The output admits nothing.  It lists leads -- document digest + quoted line -- and classifies each
blocking cell as ``CANDIDATE_FIELD`` (work here, no network), ``SOURCE_GAP_CONFIRMED`` (the hull is
present in the reservoir, the quantity is nowhere stated) or ``NOT_IN_RESERVOIR``.  A document whose
publisher URL was not recorded at fetch time cannot become a library source from this tool alone:
criterion (2) of the anchor is citation *plus* digest, and half of that is not a citation.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any

if __package__ in (None, ""):  # allow ``python tools/baseline_fleet/refs_reservoir_sweep.py``
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.fleet.reservoir import (  # noqa: E402
    blocking_cells,
    classify,
    printable_ratio,
    scan_document,
    sha256_file,
    summarise,
    text_digest,
)

SUFFIXES = {".pdf", ".html", ".htm", ".txt", ".bin"}


def _text_of(path: Path) -> tuple[str, str | None]:
    """(text, error) for one reservoir file -- PDF via MuPDF, everything else as UTF-8 bytes."""
    if path.suffix.lower() == ".pdf":
        try:
            import fitz  # PyMuPDF

            with fitz.open(path) as document:
                return "\n".join(page.get_text() for page in document), None
        except Exception as error:  # noqa: BLE001 - a broken PDF is a recorded outcome, not a crash
            return "", f"{type(error).__name__}: {error}"
    try:
        return path.read_text(errors="replace"), None
    except OSError as error:
        return "", f"{type(error).__name__}: {error}"


def library_digests(library: dict[str, Any]) -> set[str]:
    """Every SHA-256 the library already cites, collected by walking the JSON."""
    found: set[str] = set()

    def walk(node: Any) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                if key in ("sha256", "page_sha256") and isinstance(value, str):
                    found.add(value.lower())
                walk(value)
        elif isinstance(node, list):
            for item in node:
                walk(item)
        elif isinstance(node, str) and re.fullmatch(r"[0-9a-f]{64}", node):
            found.add(node)

    walk(library)
    return found


def collect(roots: list[Path]) -> dict[str, Path]:
    """Digest -> first path, over every reservoir file (duplicates by content collapse)."""
    documents: dict[str, Path] = {}
    for root in roots:
        if not root.exists():
            continue
        paths = [root] if root.is_file() else sorted(root.rglob("*"))
        for path in paths:
            if not path.is_file() or path.suffix.lower() not in SUFFIXES:
                continue
            if path.stat().st_size == 0:
                continue
            documents.setdefault(sha256_file(str(path)), path)
    return documents


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--library", default="benchmarks/verified/fleet_refs/crossref.json")
    parser.add_argument("--out", default="benchmarks/verified/fleet_refs/reservoir_sweep.json")
    parser.add_argument(
        "--reservoir",
        action="append",
        default=None,
        help="root to scan (repeatable); default: the library's own raw cache",
    )
    parser.add_argument("--fidelity", default="benchmarks/verified/fleet_refs/fidelity.json")
    parser.add_argument("--markdown", default=None, help="also write a generated table here")
    args = parser.parse_args(argv)

    library = json.loads(Path(args.library).read_text())
    fidelity_path = Path(args.fidelity)
    if "determinability" not in library and fidelity_path.exists():
        library["determinability"] = json.loads(fidelity_path.read_text())["determinability"]
    cells = blocking_cells(library)
    quantities = sorted({cell["quantity"] for cell in cells})
    cited = library_digests(library)
    roots = [Path(r) for r in (args.reservoir or ["artifacts/fleet_refs/raw"])]

    documents = []
    errors = []
    seen_text: dict[str, str] = {}
    for digest, path in sorted(collect(roots).items(), key=lambda item: str(item[1])):
        text, error = _text_of(path)
        if error:
            errors.append({"name": path.name, "sha256": digest, "error": error})
            continue
        if printable_ratio(text) < 0.90:
            errors.append({"name": path.name, "sha256": digest, "error": "not text (printable < 90 %)"})
            continue
        content = text_digest(text)
        if content in seen_text:
            # A PDF and its pre-extracted ``.txt`` sidecar are one document, not two statements.
            errors.append({"name": path.name, "sha256": digest, "error": f"duplicate of {seen_text[content]}"})
            continue
        seen_text[content] = path.name
        documents.append(
            scan_document(
                text,
                name=path.name,
                sha256=digest,
                bytes_=path.stat().st_size,
                in_library=digest in cited,
                quantities=quantities,
            )
        )

    cells_report = classify(cells, documents)
    payload = {
        "schema": "tensorcad.fleet_refs.reservoir_sweep/1.0",
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "generated_by": "python tools/baseline_fleet/refs_reservoir_sweep.py",
        "library": args.library,
        "reservoirs": [str(r) for r in roots],
        "policy": "a candidate is a lead with a pinned digest and a quoted line; it is never an "
        "admitted value and never changes a library verdict",
        "quantities_scanned": quantities,
        "summary": summarise(documents, cells_report),
        "errors": errors,
        "documents": documents,
        "cells": cells_report,
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(payload, indent=1, ensure_ascii=False) + "\n")

    if args.markdown:
        Path(args.markdown).write_text(markdown(payload))
    print(json.dumps(payload["summary"], indent=1, ensure_ascii=False))
    return 0


def markdown(payload: dict[str, Any]) -> str:
    """The generated tables of ``docs/fleet_refs_reservoir_sweep.md``."""
    lines = [
        "<!-- generated by tools/baseline_fleet/refs_reservoir_sweep.py; do not edit -->",
        "",
        f"Reservoir: {payload['summary']['documents_scanned']} documents "
        f"({payload['summary']['documents_in_library']} already library sources, "
        f"{payload['summary']['documents_outside_library']} outside it); "
        f"{payload['summary']['blocking_cells']} blocking cells scanned.",
        "",
        "| hull / quantity | library | determinability | documents naming hull | documents stating it | verdict |",
        "|---|---|---|---:|---:|---|",
    ]
    for cell in payload["cells"]:
        lines.append(
            f"| {cell['hull']} / {cell['quantity']} | {cell['library_verdict']} | "
            f"{cell['class'] or '-'} | {cell['documents_mentioning_hull']} | "
            f"{cell['documents_stating_quantity']} | {cell['verdict']} |"
        )
    candidates = [cell for cell in payload["cells"] if cell["candidates"]]
    lines += ["", "## Candidate leads (not admitted values)", ""]
    if not candidates:
        lines.append("None: no reservoir document states any blocking cell.")
    for cell in candidates:
        lines += [
            f"### {cell['hull']} / {cell['quantity']}",
            "",
            "| document | in library | sha256 (first 12) | numeral | quoted line |",
            "|---|---|---|---|---|",
        ]
        for hit in cell["candidates"]:
            for mention in hit["mentions"]:
                lines.append(
                    f"| {hit['document']} | {'yes' if hit['in_library'] else 'no'} | "
                    f"`{hit['sha256'][:12]}` | `{mention['numeral']}` | {mention['quote']} |"
                )
        lines.append("")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())