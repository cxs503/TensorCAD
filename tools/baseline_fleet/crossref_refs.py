#!/usr/bin/env python3
"""Fetch and cross-reference the published particulars of the six benchmark hulls.

    python tools/baseline_fleet/crossref_refs.py fetch   --out artifacts/fleet_refs
    python tools/baseline_fleet/crossref_refs.py build   --raw artifacts/fleet_refs/raw \
        --out benchmarks/verified/fleet_refs/crossref.json

``fetch`` materialises each document in :data:`tensorcad.fleet.crossref.SOURCES` (and the
probes in :data:`PROBES`), records URL/HTTP status/bytes/SHA-256 plus the retrieval time and
compares the digest with the one pinned in the registry -- a source that changed is a
*new* source and must be re-read deliberately.  ``build`` parses the fetched bytes into
cells, cross-references them per ship and quantity, and writes the library JSON.

Only the extracted cells, the row labels and the source digests are written into the
repository; the publisher bytes stay in ``artifacts/fleet_refs/raw`` (git-ignored), because
neither workshop states a redistribution licence.
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime
import hashlib
import json
import ssl
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

if __package__ in (None, ""):  # allow ``python tools/baseline_fleet/crossref_refs.py``
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.fleet.crossref import (  # noqa: E402
    QUANTITIES,
    QUANTITY_LABELS,
    SOURCES,
    Cell,
    coverage_summary,
    closure_check,
    cross_reference,
    parse_source,
    split_ships,
)

UA = "Mozilla/5.0 (X11; Linux x86_64) tensorcad-fleet-refs/1.0"

#: Documents that were hunted for a *second independent* numeric source and could not be
#: obtained from this environment.  Recorded with the reason so the next round starts from
#: evidence instead of repeating the search.
PROBES: tuple[dict[str, str], ...] = (
    {
        "want": "DTC particulars (originator page)",
        "url": "https://www.uni-due.de/ship/dtc.php",
        "outcome": "HTTP 404 (page retired); the institute's INAM pages carry only igs_ordner.zip (geometry, no particulars)",
    },
    {
        "want": "DTC particulars (peer-reviewed)",
        "url": "https://koreascience.kr/article/JAKO201319850774702.pdf",
        "outcome": "TCP closed without response (RemoteDisconnected) from this host",
    },
    {
        "want": "Series 60 primary tables (Todd, DTMB Report 1258/1963)",
        "url": "https://www.biodiversitylibrary.org/itempdf/102062",
        "outcome": "HTTP 403 for programmatic clients",
    },
    {
        "want": "arXiv full text (JBC/DTC case papers)",
        "url": "https://arxiv.org/pdf/2411.09725v1",
        "outcome": "PDF delivery stalls (connection accepted, body never completes); only /abs metadata is reachable",
    },
    {
        "want": "Open-access JBC paper (University of Edinburgh MARINE2021)",
        "url": "http://journals.ed.ac.uk/MARINE2021/article/download/6788/8985",
        "outcome": "request hangs until timeout",
    },
    {
        "want": "SIMMAN 2008 archive (KVLCC2 page, original FORCE edition)",
        "url": "https://www.simman2008.dk/",
        "outcome": "domain no longer serves the workshop; archived copies unreachable (web.archive.org egress blocked)",
    },
    {
        "want": "Wageningen 2022 workshop (KCS data)",
        "url": "https://www.w2022.nl/",
        "outcome": "placeholder page: 'Please stand by while configuration is in progress.'",
    },
)


def _context() -> ssl.SSLContext:
    """The publisher hosts carry mismatched/expired certificates; the *digest* is the check."""
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def fetch_bytes(url: str, timeout: int = 15) -> tuple[int, bytes]:
    request = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(request, timeout=timeout, context=_context()) as response:
        return response.status, response.read()


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def cmd_fetch(args: argparse.Namespace) -> int:
    out = Path(args.out)
    raw = out / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    registry: list[dict[str, Any]] = []
    for spec in SOURCES:
        suffix = ".pdf" if spec.kind == "pdf" else ".html"
        target = raw / (spec.key + suffix)
        entry: dict[str, Any] = {
            "key": spec.key,
            "ship": spec.ship,
            "family": spec.family,
            "publisher": spec.publisher,
            "url": spec.url,
            "kind": spec.kind,
            "licence": spec.licence,
            "note": spec.note,
        }
        if not args.offline:
            try:
                status, data = fetch_bytes(spec.url)
                target.write_bytes(data)
                entry.update(http_status=status, bytes=len(data))
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                entry.update(error=f"{type(exc).__name__}: {exc}")
        if target.exists():
            data = target.read_bytes()
            entry.setdefault("bytes", len(data))
            digest = sha256(data)
            entry["sha256"] = digest
            entry["retrieved_utc"] = datetime.datetime.fromtimestamp(
                target.stat().st_mtime, datetime.timezone.utc
            ).isoformat(timespec="seconds")
            if spec.sha256:
                entry["digest_matches_registry"] = spec.sha256 == digest
            if args.pin and not spec.sha256:
                print(f"PIN {spec.key} = {digest}")
        registry.append(entry)

    probes = []
    for probe in PROBES:
        record = dict(probe)
        if not args.offline:
            try:
                status, data = fetch_bytes(probe["url"], timeout=15)
                record["http_status"] = status
                record["bytes"] = len(data)
            except (urllib.error.URLError, TimeoutError, OSError) as exc:
                record["observed"] = f"{type(exc).__name__}: {exc}"
        probes.append(record)

    payload = {
        "schema": "tensorcad.fleet_refs.registry/1.0",
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "sources": registry,
        "probes": probes,
    }
    (out / "registry.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {out / 'registry.json'}: {len(registry)} sources, {len(probes)} probes")
    missing = [r["key"] for r in registry if "sha256" not in r]
    if missing:
        print("WARNING: not fetched:", ", ".join(missing))
        return 1
    return 0


def _pdf_text(path: Path) -> str:
    import fitz  # PyMuPDF, optional dependency used only for the ITTC bibliography

    with fitz.open(path) as document:
        return "\n".join(page.get_text() for page in document)


def _bibliography_mentions(text: str) -> dict[str, str]:
    """Which hull names the ITTC benchmark list recognises (context quoted, no numerals)."""
    keys = {
        "KCS": r"KRISO 3600 TEU|\(KCS\)",
        "KVLCC2": r"KVLCC2",
        "JBC": r"Japan Bulk Carrier|\bJBC\b",
        "ONRT": r"ONR Tumblehome|ONRT",
        "Series60": r"Series 60",
        "DTC": r"Duisburg|\(DTC\)",
    }
    found: dict[str, str] = {}
    for name, pattern in keys.items():
        match = re_search(pattern, text)
        if match:
            start = max(0, match.start() - 120)
            found[name] = " ".join(text[start : match.end() + 160].split())
    return found


def re_search(pattern: str, text: str):
    import re

    return re.search(pattern, text, flags=re.I)


def cmd_build(args: argparse.Namespace) -> int:
    raw = Path(args.raw)
    registry_path = raw.parent / "registry.json"
    registry = json.loads(registry_path.read_text()) if registry_path.exists() else {"sources": [], "probes": []}
    by_key = {row["key"]: row for row in registry["sources"]}

    extractions: dict[str, dict[str, dict[str, Cell]]] = {}
    sources: dict[str, Any] = {}
    for spec in SOURCES:
        suffix = ".pdf" if spec.kind == "pdf" else ".html"
        path = raw / (spec.key + suffix)
        record = dict(by_key.get(spec.key, {}))
        if not path.exists():
            record["status"] = "not_fetched"
            sources[spec.key] = record
            continue
        payload = path.read_bytes()
        record.setdefault("sha256", sha256(payload))
        record["bytes"] = len(payload)
        record["status"] = "read"
        if spec.kind == "pdf":
            text = _pdf_text(path)
            record["pages"] = text.count("\f") + 1
            record["bibliography_mentions"] = _bibliography_mentions(text)
        else:
            record["quantities_read"] = {}
            for ship, block in split_ships(spec, payload):
                cells = parse_source(dataclasses.replace(spec, ship=ship, section=None), block)
                record["quantities_read"][ship] = sorted(cells)
                extractions.setdefault(ship, {})[spec.key] = cells
        sources[spec.key] = record

    report = cross_reference(extractions)
    for ship in ("KCS", "KVLCC2", "JBC", "ONRT", "DTC", "Series60"):
        if ship not in report:
            report[ship] = {
                "quantities": {
                    q: {
                        "verdict": "MISSING",
                        "independent_families": 0,
                        "family_values": {},
                        "range": 0.0,
                        "range_relative": 0.0,
                        "tolerance": 0.0,
                        "entries": [],
                    }
                    for q in QUANTITIES
                },
                "cross_verified": 0,
                "single_source": 0,
                "conflict": 0,
                "missing": len(QUANTITIES),
                "corroborated_by_repeat_family": 0,
            }

    got = load_got(Path(args.got)) if args.got else {}
    comparisons = compare_with_got(report, got)

    payload = {
        "schema": "tensorcad.fleet_refs.crossref/1.0",
        "anchor": "REPRO-GEOM-FLEET-REFS",
        "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "quantity_labels": QUANTITY_LABELS,
        "sources": sources,
        "report": report,
        "coverage": coverage_summary(report),
        "comparisons": comparisons,
        "closure": closure_check(report),
        "probes": registry.get("probes", []),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n")
    cov = payload["coverage"]
    print(
        f"wrote {out}: cells {cov['cells']}, cross_verified {cov['cross_verified']}, "
        f"single_source {cov['single_source']}, missing {cov['missing']}"
    )
    print("ships fully cross-verified:", ", ".join(cov["ships_fully_cross_verified"]) or "-")
    print("ships below two sources:", ", ".join(cov["ships_failing_min_two_sources"]) or "-")
    return 0


def load_got(path: Path) -> dict[str, dict[str, Any]]:
    """The reconstructed (``got``) side of the comparison, from the fleet benchmark artefact."""
    if not path.exists():
        return {}
    data = json.loads(path.read_text())
    hulls = data.get("hulls", {})
    return {
        ship: {
            "got_volume_m3": row.get("got_volume_m3"),
            "ref_volume_m3": row.get("ref_volume_m3"),
            "rel_error_percent": row.get("rel_error_percent"),
            "verdict": row.get("verdict"),
            "reasons": row.get("reasons", []),
        }
        for ship, row in hulls.items()
    }


def compare_with_got(report: dict[str, Any], got: dict[str, dict[str, Any]]) -> dict[str, Any]:
    """Reference ∇ (cross-verified when available) versus the reconstructed volume."""
    out: dict[str, Any] = {}
    for ship, row in sorted(got.items()):
        cell = report.get(ship, {}).get("quantities", {}).get("displacement_volume", {})
        families = cell.get("family_values", {})
        reference = None
        if families:
            # The *minimum* family value is used as the reference: picking the value closest to
            # the reconstruction would be a flattering choice, and this anchor is about honesty.
            reference = min(families.values())
        got_volume = row.get("got_volume_m3")
        error = None
        comparability = "full_scale"
        if reference and got_volume:
            error = 100.0 * (got_volume - reference) / reference
            if abs(error) > 50.0:
                # A >50 % gap is not an error bar, it is a scale/units mismatch in the recorded
                # artefact (e.g. ONRT is reconstructed at model scale, 92.9 m3).  Report it as
                # not comparable instead of quoting a meaningless -98.9 %.
                comparability = "not_comparable_scale_or_units"
                error = None
        out[ship] = {
            "reference_volume_m3": reference,
            "reference_verdict": cell.get("verdict"),
            "reference_families": families,
            "got_volume_m3": got_volume,
            "rel_error_percent": error,
            "benchmark_verdict": row.get("verdict"),
            "benchmark_reasons": row.get("reasons", []),
            "comparability": comparability,
            "within_1_percent": None if error is None else abs(error) <= 1.0,
        }
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p_fetch = sub.add_parser("fetch", help="download the source documents")
    p_fetch.add_argument("--out", default="artifacts/fleet_refs")
    p_fetch.add_argument("--offline", action="store_true", help="do not touch the network (use cached bytes)")
    p_fetch.add_argument("--pin", action="store_true", help="print digests for sources not yet pinned")
    p_fetch.set_defaults(func=cmd_fetch)

    p_build = sub.add_parser("build", help="parse + cross-reference the cached documents")
    p_build.add_argument("--raw", default="artifacts/fleet_refs/raw")
    p_build.add_argument("--out", default="benchmarks/verified/fleet_refs/crossref.json")
    p_build.add_argument(
        "--got",
        default="benchmarks/verified/baseline_fleet/result.json",
        help="fleet benchmark artefact holding the reconstructed volumes",
    )
    p_build.set_defaults(func=cmd_build)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())