#!/usr/bin/env python3
"""Hunt open-access documents that may state a *missing* principal particular.

Round 4 of the ``REPRO-GEOM-FLEET-REFS`` reference library.  The library still has cells with
no independent publisher (``Series60`` has no document at all; ``JBC``/``KVLCC2``/``ONRT``/
``DTC`` are short on LCB and Lwl).  This tool does the search part mechanically:

1. ask the OpenAlex API for open-access works matching a query;
2. download each PDF/HTML candidate into a scratch directory with a short timeout;
3. extract text and print every line that mentions the hull *and* the wanted quantity, with
   the candidate URL, so a human reading can decide what to pin.

It never writes into the registry: candidates are hunches until their digest is pinned and
their quote read, which happens in ``crossref_round4.py``.

    python tools/baseline_fleet/hunt_round4.py --ship ONRT --want lcb --query "ONR tumblehome"
"""

from __future__ import annotations

import argparse
import hashlib
import http.client
import json
import re
import ssl
import sys
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

UA = "Mozilla/5.0 (X11; Linux x86_64) tensorcad-fleet-refs/1.0"

WANT_PATTERNS = {
    "lcb": r"lcb|lcg|longitudinal cent(?:er|re) of buoyancy|centre of buoyancy|center of buoyancy",
    "lwl": r"lwl|length (?:on|at) (?:the )?waterline|waterline length",
    "displacement_volume": r"displacement|∇",
    "cb": r"block coefficient|\bcb\b",
}


def context() -> ssl.SSLContext:
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def get(url: str, timeout: int = 12, accept: str = "*/*") -> tuple[int, bytes]:
    request = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": accept})
    with urllib.request.urlopen(request, timeout=timeout, context=context()) as response:
        return response.status, response.read()


def openalex(query: str, limit: int = 8) -> list[dict]:
    url = "https://api.openalex.org/works?" + urllib.parse.urlencode(
        {"search": query, "per-page": limit, "filter": "is_oa:true"}
    )
    _, body = get(url, accept="application/json")
    out = []
    for work in json.loads(body).get("results", []):
        loc = work.get("best_oa_location") or {}
        pdf = loc.get("pdf_url") or (work.get("open_access") or {}).get("oa_url")
        if not pdf:
            continue
        out.append(
            {
                "title": work.get("title"),
                "year": work.get("publication_year"),
                "doi": work.get("doi"),
                "url": pdf,
            }
        )
    return out


def text_of(path: Path) -> str:
    if path.suffix.lower() == ".pdf":
        import fitz

        with fitz.open(path) as document:
            return "\n".join(page.get_text() for page in document)
    return path.read_text("utf-8", "replace")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ship", required=True)
    parser.add_argument("--want", required=True, choices=sorted(WANT_PATTERNS))
    parser.add_argument("--query", action="append", required=True)
    parser.add_argument("--out", default="artifacts/fleet_refs/hunt")
    parser.add_argument("--limit", type=int, default=8)
    parser.add_argument("--context", type=int, default=3)
    args = parser.parse_args(argv)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ship_re = re.compile(args.ship.replace("_", "[ _]?"), re.I)
    want_re = re.compile(WANT_PATTERNS[args.want], re.I)

    seen: set[str] = set()
    for query in args.query:
        try:
            works = openalex(query, args.limit)
        except (urllib.error.URLError, http.client.HTTPException, TimeoutError, ValueError) as exc:
            print(f"SEARCH-ERROR {query}: {type(exc).__name__}: {str(exc)[:100]}")
            continue
        for work in works:
            if work["url"] in seen:
                continue
            seen.add(work["url"])
            suffix = ".pdf" if ".pdf" in work["url"].lower() else ".html"
            name = hashlib.sha256(work["url"].encode()).hexdigest()[:10] + suffix
            target = out / name
            if not target.exists():
                try:
                    _, data = get(work["url"], timeout=20)
                except (urllib.error.URLError, http.client.HTTPException, TimeoutError, OSError) as exc:
                    print(f"FETCH-ERROR {work['url']}: {type(exc).__name__}: {str(exc)[:80]}")
                    continue
                target.write_bytes(data)
            try:
                text = text_of(target)
            except Exception as exc:  # noqa: BLE001 - a broken file is a hunt outcome
                print(f"PARSE-ERROR {target.name}: {type(exc).__name__}")
                continue
            lines = text.splitlines()
            hits = False
            for index, line in enumerate(lines):
                if not want_re.search(line):
                    continue
                window = " ".join(lines[max(0, index - args.context) : index + args.context + 1])
                if not ship_re.search(window):
                    continue
                if not hits:
                    print(f"\n=== {work['year']} {work['url']}\n    {str(work['title'])[:110]}")
                    hits = True
                print(f"  [{target.name}:{index}] {' '.join(line.split())[:200]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())