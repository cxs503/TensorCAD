#!/usr/bin/env python3
"""Fetch source pages and geometry archives for the baseline fleet, compute sha256, and report.

This tool is the reproducible entry point for GAP-CAD-T1-SOURCE-HUNT.
It fetches every known source URL, computes sha256, and compares against
the recorded digests in tensorcad/fleet/official.py.

Usage:
    python tools/source_hunt/fetch_and_hash.py [--out DIR] [--skip-geometry]

Output:
    - JSON report with per-URL sha256, byte count, and match status
    - Optional: saved bytes to --out directory
"""

from __future__ import annotations

import argparse
import hashlib
import json
import ssl
import sys
from pathlib import Path
from typing import Any

import urllib.request

# T2015 site has SSL hostname mismatch; disable hostname verification
_CTX = ssl.create_default_context()
_CTX.check_hostname = False
_CTX.verify_mode = ssl.CERT_NONE

#: All known source URLs for the six standard models.
SOURCES: dict[str, dict[str, str]] = {
    # T2015 workshop pages
    "t2015_kcs_gc": {
        "url": "https://www.t2015.nmri.go.jp/kcs_gc.html",
        "type": "page",
        "ship": "KCS",
    },
    "t2015_jbc_gc": {
        "url": "https://www.t2015.nmri.go.jp/jbc_gc.html",
        "type": "page",
        "ship": "JBC",
    },
    "t2015_onrt_gc": {
        "url": "https://www.t2015.nmri.go.jp/onrt_gc.html",
        "type": "page",
        "ship": "ONRT",
    },
    "t2015_onrt_desc": {
        "url": "https://www.t2015.nmri.go.jp/onrt.html",
        "type": "page",
        "ship": "ONRT",
    },
    "t2015_main": {
        "url": "https://www.t2015.nmri.go.jp/",
        "type": "page",
        "ship": None,
    },
    # SIMMAN 2020
    "simman2020_summary": {
        "url": "https://simman2020.kr/contents/Ship_data_summary.php",
        "type": "page",
        "ship": None,
    },
    # Geometry archives (T2015)
    "t2015_kcs_geom": {
        "url": "https://www.t2015.nmri.go.jp/file/Geometry_IGES_files/kcs/KCS_hull_SVA.zip",
        "type": "geometry",
        "ship": "KCS",
    },
    "t2015_jbc_geom": {
        "url": "https://www.t2015.nmri.go.jp/file/Geometry_IGES_files/jbc/JBC_IGES.zip",
        "type": "geometry",
        "ship": "JBC",
    },
    "t2015_onrt_wo_gap": {
        "url": "https://www.t2015.nmri.go.jp/file/Geometry_IGES_files/onrt/onrt_wo_gap.zip",
        "type": "geometry",
        "ship": "ONRT",
    },
    "t2015_onrt_w_gap": {
        "url": "https://www.t2015.nmri.go.jp/file/Geometry_IGES_files/onrt/onrt_w_gap.zip",
        "type": "geometry",
        "ship": "ONRT",
    },
    "t2015_onrt_propeller": {
        "url": "https://www.t2015.nmri.go.jp/file/Geometry_IGES_files/onrt/propeller_IIHR2.zip",
        "type": "geometry",
        "ship": "ONRT",
    },
    # Geometry archives (SIMMAN 2020)
    "simman2020_kvlcc2_geom": {
        "url": "https://simman2020.kr/download/KVLCC2-Hull.zip",
        "type": "geometry",
        "ship": "KVLCC2",
    },
    "simman2020_onrt_geom": {
        "url": "https://simman2020.kr/download/onrt_wo_gap.zip",
        "type": "geometry",
        "ship": "ONRT",
    },
    # DTC sources (known dead)
    "dtc_uni_due_page": {
        "url": "https://www.uni-due.de/ship/dtc.php",
        "type": "page",
        "ship": "DTC",
        "expected_status": 404,
    },
    "dtc_uni_due_geom": {
        "url": "https://www.uni-due.de/ship/dtc_download/DTC_IGES.zip",
        "type": "geometry",
        "ship": "DTC",
        "expected_status": 404,
    },
}

#: Recorded sha256 digests from tensorcad/fleet/official.py (PAGES and ARCHIVES).
RECORD_DIGESTS: dict[str, str] = {
    "t2015_kcs_gc": "68d9c24600cb699b16394fc5c32762ebc875c856d14ef8460e1d0ebf4e7c78a4",
    "t2015_jbc_gc": "fc2ce9054a8b7842162bac94365853fa076da2527e8c7e455832f9710f58f875",
    "t2015_onrt_gc": "5e13427332de8be17b9dc2245dc5830300baf2e2340efcd8a80eee2269c43a11",
    "simman2020_summary": "c1a668f6d2440c5d83a83c3616285e9288ef618dbc9a2feebd26627be3f4ebb3",
    "t2015_kcs_geom": "b9e334853857b819",  # truncated in official.py
    "t2015_jbc_geom": "e19800ee954f4960",
    "t2015_onrt_wo_gap": "bd747eff8f531ba4",
    "simman2020_kvlcc2_geom": "dc032647f00a68d3",
    "simman2020_onrt_geom": "bd747eff8f531ba4",
}


def fetch(url: str, timeout: int = 20) -> tuple[bytes, int, str | None]:
    """Fetch a URL, return (data, status_code, error)."""
    try:
        req = urllib.request.Request(
            url, headers={"User-Agent": "TensorCAD/0.4 source-hunt"}
        )
        resp = urllib.request.urlopen(req, timeout=timeout, context=_CTX)
        data = resp.read()
        return data, resp.status, None
    except Exception as e:
        return b"", 0, str(e)


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def matches_record(key: str, digest: str) -> bool:
    recorded = RECORD_DIGESTS.get(key)
    if recorded is None:
        return False  # no record to compare
    # Recorded digests may be truncated (first 16 hex chars)
    if len(recorded) < 64:
        return digest.startswith(recorded)
    return digest == recorded


def run(out_dir: Path | None = None, skip_geometry: bool = False) -> dict[str, Any]:
    results: dict[str, Any] = {}

    for key, spec in SOURCES.items():
        if skip_geometry and spec["type"] == "geometry":
            results[key] = {"url": spec["url"], "skipped": True, "reason": "skip_geometry"}
            continue

        data, status, error = fetch(spec["url"])

        entry: dict[str, Any] = {
            "url": spec["url"],
            "type": spec["type"],
            "ship": spec.get("ship"),
            "status": status,
            "error": error,
        }

        if error:
            expected = spec.get("expected_status")
            if expected and "404" in error:
                entry["expected_404"] = True
                entry["status"] = 404
            else:
                entry["fetch_ok"] = False
                results[key] = entry
                continue

        entry["fetch_ok"] = True
        entry["bytes"] = len(data)
        entry["sha256"] = sha256_hex(data)
        entry["sha256_short"] = entry["sha256"][:16]
        entry["matches_record"] = matches_record(key, entry["sha256"])

        if out_dir:
            out_dir.mkdir(parents=True, exist_ok=True)
            suffix = ".html" if spec["type"] == "page" else ".zip"
            (out_dir / f"{key}{suffix}").write_bytes(data)
            entry["saved_to"] = str(out_dir / f"{key}{suffix}")

        results[key] = entry

    return results


def main():
    parser = argparse.ArgumentParser(description="Fetch and hash fleet source URLs")
    parser.add_argument("--out", type=Path, default=None, help="Directory to save fetched bytes")
    parser.add_argument("--skip-geometry", action="store_true", help="Skip geometry archives")
    args = parser.parse_args()

    results = run(out_dir=args.out, skip_geometry=args.skip_geometry)

    # Summary
    ok = sum(1 for r in results.values() if r.get("fetch_ok"))
    fail = sum(1 for r in results.values() if not r.get("fetch_ok") and not r.get("expected_404"))
    match = sum(1 for r in results.values() if r.get("matches_record"))
    skipped = sum(1 for r in results.values() if r.get("skipped"))

    print(f"\n=== Source Hunt Summary ===")
    print(f"Total sources: {len(results)}")
    print(f"Fetched OK: {ok}")
    print(f"Failed: {fail}")
    print(f"Matches recorded digest: {match}")
    print(f"Skipped: {skipped}")

    # Print per-source status
    for key, entry in sorted(results.items()):
        status = "OK" if entry.get("fetch_ok") else ("EXPECTED_404" if entry.get("expected_404") else "FAIL")
        sha = entry.get("sha256_short", "N/A")
        match_str = "✓" if entry.get("matches_record") else ("?" if entry.get("sha256") is None else "✗")
        print(f"  {key:30s} {status:12s} sha={sha} match={match_str}")

    # Save JSON report
    report_path = Path("/tmp/source_hunt_report.json")
    report_path.write_text(json.dumps(results, indent=2))
    print(f"\nFull report: {report_path}")

    return 0 if fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
