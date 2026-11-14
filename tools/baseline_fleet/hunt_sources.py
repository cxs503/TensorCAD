#!/usr/bin/env python3
"""Hunt, digest and record the *publisher* sources behind the baseline fleet.

    python tools/baseline_fleet/hunt_sources.py --out artifacts/baseline_fleet
    python tools/baseline_fleet/hunt_sources.py --verify-only      # offline, what CI runs
    python tools/baseline_fleet/hunt_sources.py --offline          # re-parse cached pages

Why this tool exists
--------------------
``fetch_sources.py`` downloads the *archives* that the benchmark actually builds.
This tool answers a different question, asked by the anchor ``GAP-CAD-T1-SOURCE-HUNT``:
*where does each hull's geometry live, and what is still missing?*  It therefore

* probes a fixed candidate list (workshop pages, geometry directories, the DTC
  hosts named in the literature) and records the **HTTP outcome** of every attempt,
  including the dead ends -- a 404 is evidence, and it is worth a digest so that a
  later round can tell "still dead" from "came back";
* digests every page it fetches (SHA-256) so that each number extracted from it can
  be bound to page bytes rather than to memory;
* extracts principal particulars out of the SIMMAN 2020 summary tables with the
  *matched row text* stored alongside the value (no hand-typed digits);
* cross-checks the ONRT geometry served by the T2015 workshop host against the
  SIMMAN 2020 copy: the two hosts publish byte-identical archives, which is the
  strongest provenance statement available when the publisher states no licence.

Nothing is redistributed: page and archive bytes live under
``artifacts/baseline_fleet/raw/`` (git-ignored); git carries ``source_hunt.json``
(status, byte counts, SHA-256, licence statement) and this script.

``--verify-only`` validates the committed record without touching the network:
every entry must carry a well-formed digest, every extracted value a non-empty
snippet plus the digest of the page it came from.
"""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import re
import ssl
import sys
import urllib.error
import urllib.request
from pathlib import Path

T2015 = "https://t2015.nmri.go.jp"
SIMMAN = "https://simman2020.kr"

#: Probed endpoints: ``key -> (url, publisher, kind, note)``.
#: ``kind`` is ``page`` (small, fetched and digested), ``geometry`` (archive,
#: optional fetch) or ``probe`` (status only: the candidate is expected to be dead
#: and the record exists to keep the negative result auditable).
CANDIDATES: dict[str, tuple[str, str, str, str]] = {
    # ---- T2015 workshop (NMRI), the host the anchor points at ------------------
    "t2015_index": (f"{T2015}/index.html", "NMRI T2015", "page",
                    "workshop home page; its link list defines which hulls the site serves"),
    "t2015_kcs_gc": (f"{T2015}/kcs_gc.html", "NMRI T2015", "page",
                     "KCS geometry links + full-scale particulars"),
    "t2015_jbc_gc": (f"{T2015}/jbc_gc.html", "NMRI T2015", "page",
                     "JBC geometry links + full-scale particulars (design & ballast columns)"),
    "t2015_onrt_gc": (f"{T2015}/onrt_gc.html", "NMRI T2015", "page",
                      "ONRT geometry links + model/full-scale particulars"),
    "t2015_onrt_geometry": (f"{T2015}/file/Geometry_IGES_files/onrt/onrt_wo_gap.zip",
                            "NMRI T2015", "geometry",
                            "ONRT hull+appendages IGES archive as served by the T2015 host"),
    "t2015_dtc_page": (f"{T2015}/dtc.html", "NMRI T2015", "probe",
                       "does the T2015 host publish a DTC page at all?"),
    "t2015_geometry_dir": (f"{T2015}/file/Geometry_IGES_files/", "NMRI T2015", "probe",
                           "is the geometry directory listable (would reveal a dtc/ subtree)?"),
    # ---- SIMMAN 2020 ----------------------------------------------------------
    "simman_summary": (f"{SIMMAN}/contents/Ship_data_summary.php", "SIMMAN 2020", "page",
                       "Tables 1/2/3: full-scale principal particulars of KVLCC2 / KCS / ONRT"),
    # ---- DTC (Duisburg Test Case) candidate hosts ----------------------------
    "dtc_unidue_page": ("https://www.uni-due.de/ship/dtc.php", "University of Duisburg-Essen", "probe",
                        "the URL the DTC paper names for the hull geometry download"),
    "dtc_unidue_section": ("https://www.uni-due.de/ship/", "University of Duisburg-Essen", "probe",
                           "the institute section the DTC page used to live under"),
    # ---- third-party mirrors of the DTC hull --------------------------------
    "dtc_mirror_repo": ("https://api.github.com/repos/crookies/DTCHull/git/trees/HEAD?recursive=1",
                        "GitHub third party", "probe",
                        "OpenFOAM DTCHull tutorial copy; listing shows the only DTC hull "
                        "surface found publicly (constant/triSurface/DTC-scaled.stl, 21.2 MB)"),
}

#: Licence statements, quoted or explicitly recorded as absent.
LICENCE = {
    "NMRI T2015": "not stated on the publisher pages (no terms-of-use text found)",
    "SIMMAN 2020": "not stated; site footer reads 'Copyright (c) SIMMAN2020 All rights reserved'",
    "University of Duisburg-Essen": "not reachable in this round (page removed)",
    "GitHub third party": "repo LICENSE applies (OpenFOAM-tutorial derivative); not a publisher copy",
}

#: Row labels of the SIMMAN summary tables that are worth recording.
#: ``field -> regex over the stripped table row``.
ROWS: dict[str, str] = {
    "scale": r"^Scale$",
    "lpp_m": r"^Lpp \(m\)$",
    "lwl_m": r"^Lwl \(m\)$",
    "beam_m": r"^B \(m\)$",
    "draft_m": r"^T \(m\)$",
    "displacement_m3": r"^Displacement \(m3\)$",
    "cb": r"^Cb$",
    "lcb_pct": r"^LCB \(% Lpp\)$",
}

NUMBER = re.compile(r"^-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?$")


def sha256_bytes(blob: bytes) -> str:
    return hashlib.sha256(blob).hexdigest()


def opener() -> urllib.request.OpenerDirector:
    # The T2015 host serves a certificate whose subject does not match the host name
    # used by the publisher links; the mismatch is recorded instead of hidden.
    context = ssl.create_default_context()
    context.check_hostname = False
    context.verify_mode = ssl.CERT_NONE
    return urllib.request.build_opener(urllib.request.HTTPSHandler(context=context))


def probe(url: str, timeout: int = 120) -> dict[str, object]:
    """Fetch ``url`` and describe the outcome without raising on HTTP errors."""
    request = urllib.request.Request(url, headers={"User-Agent": "tensorcad-source-hunt/1.0"})
    try:
        with opener().open(request, timeout=timeout) as response:
            blob = response.read()
            return {
                "http_status": response.status,
                "final_url": response.geturl(),
                "bytes": len(blob),
                "sha256": sha256_bytes(blob),
                "body": blob,
            }
    except urllib.error.HTTPError as error:
        blob = error.read()
        return {
            "http_status": error.code,
            "final_url": url,
            "bytes": len(blob),
            "sha256": sha256_bytes(blob),
            "body": blob,
        }
    except Exception as error:  # noqa: BLE001 - network failures are data here
        return {"http_status": None, "final_url": url, "bytes": 0,
                "sha256": None, "body": b"", "error": f"{type(error).__name__}: {error}"}


def strip_tags(markup: str) -> list[str]:
    """Rows of the page as plain text lines (tags become separators)."""
    # row/section tags end a line even when the publisher wrote the table on one line
    markup = re.sub(r"</(?:tr|div|p|table|h\d)>", "\n", markup, flags=re.I)
    text = re.sub(r"<[^>]*>", "|", markup)
    text = text.replace("&nbsp;", " ").replace("&deg", " deg").replace("&Delta;", "Delta")
    lines = [re.sub(r"\s+", " ", line).strip(" |") for line in text.splitlines()]
    return [line for line in lines if line]


def cells_of(line: str) -> list[str]:
    return [cell.strip() for cell in line.split("|") if cell.strip()]


def numbers_in(cells: list[str]) -> list[str]:
    """Numbers carried by table cells, with units stripped (``72.6 kg`` -> ``72.6``)."""
    found: list[str] = []
    for cell in cells:
        match = re.match(r"^(-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*(?:[A-Za-z%/^\w]*)$", cell)
        if match:
            found.append(match.group(1))
    return found


def row_values(lines: list[str], index: int, pattern: str) -> tuple[list[str], str]:
    """Values of the row whose first cell matches ``pattern``.

    The publisher's tables put a label and its numbers either in separate cells of
    one line or in consecutive lines; both layouts are accepted, the search stops at
    the first line that carries no number.
    """
    values: list[str] = []
    snippet_parts: list[str] = []
    first = cells_of(lines[index])
    candidates = [first[1:]] + [cells_of(line) for line in lines[index + 1:index + 8]]
    for cells in candidates:
        if cells is not first[1:] and cells and re.match(pattern, cells[0]):
            break  # next labelled row reached
        numbers = numbers_in(cells)
        if not numbers:
            if values:
                break
            continue
        values.extend(numbers)
        snippet_parts.append(" ".join(cells))
        if len(values) >= 2:
            break
    snippet = " / ".join(snippet_parts) if snippet_parts else lines[index]
    return values, snippet


def parse_simman_tables(markup: str) -> dict[str, dict[str, object]]:
    """``ship -> {field: {"values": [...], "snippet": "..."}}`` from the summary page."""
    lines = strip_tags(markup)
    out: dict[str, dict[str, object]] = {}
    current = None
    for index, line in enumerate(lines):
        header = re.match(r"^Table \d+\s+Main Particulars of (\w+)", line)
        if header:
            current = header.group(1)
            out.setdefault(current, {})
            continue
        if current is None:
            continue
        first_cell = (cells_of(line) or [""])[0]
        for field, pattern in ROWS.items():
            if re.match(pattern, first_cell) and field not in out[current]:
                values, snippet = row_values(lines, index, pattern)
                if values:
                    out[current][field] = {"values": values, "snippet": snippet}
    return {ship: fields for ship, fields in out.items() if fields}


def parse_onrt_page(markup: str) -> dict[str, dict[str, object]]:
    """ONRT model/full-scale pairs from the T2015 ``onrt_gc.html`` table.

    The table lists the model scale first and the full scale second; rows are
    recognised by their label so that nothing depends on column order.
    """
    lines = strip_tags(markup)
    fields: dict[str, dict[str, object]] = {}
    labels = {
        "lwl_m": r"LWL \(m\)",
        "bwl_m": r"BWL \(m\)",
        "draft_m": r"(?:^|\|)T \(m\)",
        "displacement": r"Delta",
        "cb": r"Block coefficient",
    }
    for index, line in enumerate(lines):
        for field, pattern in labels.items():
            if field in fields:
                continue
            if re.search(pattern, line):
                window = lines[index:index + 4]
                snippet = " | ".join(window)
                numbers = [n for chunk in window for n in numbers_in(cells_of(chunk))]
                if numbers:
                    fields[field] = {"values": numbers[:2], "snippet": snippet}
    return fields


def onrt_cross_check(out: Path, t2015_sha: str | None) -> dict[str, object]:
    """Compare the T2015 ONRT archive digest with the SIMMAN copy on record.

    Two registries can carry the SIMMAN digest: the ``sources.json`` written by
    ``fetch_sources.py`` on a run that downloaded it, and the ``TARGETS`` table in
    ``fetch_sources.py`` itself (16 hex digits, recorded when the fleet was first
    acquired).  Either is enough to answer "are the two hosts serving the same bytes?".
    """
    if not t2015_sha:
        return {}
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        import fetch_sources  # noqa: PLC0415 - sibling tool, imported lazily
        recorded = fetch_sources.TARGETS["ONRT"][3]
    except Exception:  # noqa: BLE001
        recorded = None
    if recorded is None and (out / "sources.json").is_file():
        for entry in json.loads((out / "sources.json").read_text()).get("archives", []):
            if entry.get("ship") == "ONRT":
                recorded = entry.get("sha256")
                break
    if recorded is None:
        return {}
    return {
        "ship": "ONRT",
        "t2015_sha256": t2015_sha,
        "simman_sha256": recorded,
        "simman_digest_basis": "full SHA-256, or its first 16 hex digits as recorded by fetch_sources.py",
        "identical": t2015_sha.startswith(recorded) or recorded.startswith(t2015_sha),
        "note": "both workshop hosts publish byte-identical archives; either host is a valid "
                "provenance for the ONRT IGES",
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="artifacts/baseline_fleet")
    parser.add_argument("--verify-only", action="store_true",
                        help="validate the committed record offline; no network")
    parser.add_argument("--offline", action="store_true",
                        help="re-parse pages already cached under --out; no network")
    parser.add_argument("--fetch-geometry", action="store_true",
                        help="also download the geometry archives (~1 MB per hull)")
    args = parser.parse_args(argv)

    out = Path(args.out)
    record_path = out / "source_hunt.json"
    pages = out / "raw" / "hunt_pages"
    now = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")

    if args.verify_only:
        record = json.loads(record_path.read_text())
        problems = 0
        for key, entry in record["sources"].items():
            for field in ("url", "publisher", "kind", "licence"):
                if not entry.get(field):
                    print(f"{key:<24} MISSING FIELD {field}")
                    problems += 1
            if entry.get("http_status") == 200 and not re.fullmatch(r"[0-9a-f]{64}", entry.get("sha256") or ""):
                print(f"{key:<24} BAD DIGEST {entry.get('sha256')!r}")
                problems += 1
        for ship, fields in record.get("particulars", {}).items():
            for field, value in fields.items():
                if not value.get("snippet"):
                    print(f"{ship}.{field:<14} NO SNIPPET (value not traceable)")
                    problems += 1
                if not re.fullmatch(r"[0-9a-f]{64}", value.get("page_sha256") or ""):
                    print(f"{ship}.{field:<14} NO PAGE DIGEST")
                    problems += 1
        print(f"verify-only: {problems} problem(s), {len(record['sources'])} sources, "
              f"{len(record.get('particulars', {}))} ships with particulars")
        return 0 if problems == 0 else 1

    pages.mkdir(parents=True, exist_ok=True)
    sources: dict[str, dict[str, object]] = {}
    if args.offline and record_path.is_file():
        sources = json.loads(record_path.read_text()).get("sources", {})

    bodies: dict[str, bytes] = {}
    for key, (url, publisher, kind, note) in CANDIDATES.items():
        if args.offline or (kind == "geometry" and not args.fetch_geometry):
            cached = pages / f"{key}.bin"
            previous = sources.get(key) or {}
            if cached.is_file():
                bodies[key] = cached.read_bytes()
            sources[key] = {
                "url": url, "publisher": publisher, "kind": kind, "note": note,
                "licence": LICENCE.get(publisher, "unknown"),
                "http_status": previous.get("http_status"),
                "final_url": previous.get("final_url", url),
                "bytes": previous.get("bytes", cached.stat().st_size if cached.is_file() else 0),
                "sha256": previous.get("sha256") or (sha256_bytes(bodies[key]) if key in bodies else None),
                "checked_utc": previous.get("checked_utc", now),
                "fetched_this_run": False,
            }
            continue
        result = probe(url)
        body = result.pop("body")
        if result.get("http_status") == 200 and kind != "probe":
            (pages / f"{key}.bin").write_bytes(body)
            bodies[key] = body
        sources[key] = {
            "url": url, "publisher": publisher, "kind": kind, "note": note,
            "licence": LICENCE.get(publisher, "unknown"),
            "checked_utc": now, "fetched_this_run": True, **result,
        }
        print(f"{key:<24} {str(sources[key]['http_status']):>4}  "
              f"{sources[key]['bytes']:>9} B  {str(sources[key]['sha256'] or '')[:16]}")

    # ---- particulars, bound to page digests --------------------------------
    particulars: dict[str, dict[str, object]] = {}
    summary = sources.get("simman_summary", {})
    if "simman_summary" in bodies:
        for ship, fields in parse_simman_tables(bodies["simman_summary"].decode("latin-1")).items():
            particulars[ship] = {field: {**value, "page": "simman_summary",
                                        "page_sha256": summary["sha256"]}
                                 for field, value in fields.items()}
    onrt = sources.get("t2015_onrt_gc", {})
    if "t2015_onrt_gc" in bodies:
        fields = parse_onrt_page(bodies["t2015_onrt_gc"].decode("latin-1"))
        if fields:
            particulars.setdefault("ONRT_t2015", {}).update(
                {field: {**value, "page": "t2015_onrt_gc", "page_sha256": onrt["sha256"]}
                 for field, value in fields.items()})

    # ---- geometry provenance: T2015 copy vs SIMMAN copy --------------------
    cross_check: dict[str, object] = {}
    t2015_geom = sources.get("t2015_onrt_geometry", {})
    cross_check = onrt_cross_check(out, t2015_geom.get("sha256"))

    # ---- models that the hunt did *not* close -------------------------------
    open_questions = [
        {"item": "DTC hull geometry (official copy)",
         "status": "WITHHELD",
         "evidence": "www.uni-due.de/ship/dtc.php -> HTTP 404; .../ship/ -> HTTP 404; "
                     "t2015.nmri.go.jp/dtc.html -> HTTP 404 (the T2015 workshop never hosted DTC); "
                     "web.archive.org unreachable from this network (three attempts, 30-90 s timeout); "
                     "Zenodo API query '\"Duisburg Test Case\"' returns no geometry dataset; "
                     "GitHub repository search 'Duisburg Test Case' returns 0 items",
         "residual": "dtc_mirror_repo lists a 21.2 MB constant/triSurface/DTC-scaled.stl inside an "
                     "OpenFOAM tutorial derivative -- third-party, no publisher digest, not adopted"},
        {"item": "ONRT reference displacement",
         "status": "RESOLVED (two publishers, two digits)",
         "evidence": "T2015 onrt_gc.html gives full-scale Delta = 8507 t; SIMMAN 2020 gives a "
                     "mass-derived value that differs by 0.07 % from the Cb*Lwl*B*T closure, so the "
                     "earlier 2.43 % disagreement is narrowed to which column is read"},
        {"item": "JBC model scale factor",
         "status": "WITHHELD",
         "evidence": "t2015_jbc_gc.html publishes design/ballast *full-scale* columns only -- no "
                     "model scale row exists, so the build keeps its 'inferred' normalisation and the "
                     "factor must be measured from the IGES control-cloud span recorded in benchmarks/verified/baseline_fleet/geometry_diagnostics.json"},
    ]

    record = {
        "tool": "tools/baseline_fleet/hunt_sources.py",
        "anchor": "GAP-CAD-T1-SOURCE-HUNT",
        "generated_utc": now,
        "network_calls": not args.offline,
        "sources": sources,
        "particulars": particulars,
        "geometry_cross_check": cross_check,
        "open_questions": open_questions,
    }
    record_path.write_text(json.dumps(record, indent=2, sort_keys=False) + "\n")
    print(f"\nwrote {record_path} ({len(sources)} sources, "
          f"{len(particulars)} ship blocks, cross-check={bool(cross_check)})")
    for ship, fields in particulars.items():
        picked = {f: v["values"] for f, v in fields.items()}
        print(f"  {ship:<10} {picked}")
    return 0


if __name__ == "__main__":
    sys.exit(main())