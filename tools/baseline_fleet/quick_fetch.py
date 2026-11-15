#!/usr/bin/env python3
"""Fetch a *subset* of the cross-reference documents with a short per-request timeout.

The full ``crossref_refs.py fetch`` walks every source plus the probe list with a 15 s
timeout each, which is fine on a workstation but exceeds the 120 s foreground budget of the
agent shell.  This helper fetches only the named sources into the same ``raw/`` layout so
``crossref_refs.py fetch --offline`` can then pin/re-hash them locally.

    python tools/baseline_fleet/quick_fetch.py --out artifacts/fleet_refs t2015_kcs ...
"""

from __future__ import annotations

import argparse
import http.client
import ssl
import sys
import urllib.error
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.fleet.crossref import SOURCES  # noqa: E402
from tools.baseline_fleet.crossref_refs import UA, sha256  # noqa: E402


def context() -> ssl.SSLContext:
    ctx = ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = ssl.CERT_NONE
    return ctx


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("keys", nargs="+")
    parser.add_argument("--out", default="artifacts/fleet_refs")
    parser.add_argument("--timeout", type=int, default=8)
    args = parser.parse_args(argv)

    raw = Path(args.out) / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    by_key = {spec.key: spec for spec in SOURCES}
    for key in args.keys:
        spec = by_key[key]
        suffix = ".pdf" if spec.kind == "pdf" else ".html"
        target = raw / (key + suffix)
        request = urllib.request.Request(spec.url, headers={"User-Agent": UA})
        try:
            with urllib.request.urlopen(request, timeout=args.timeout, context=context()) as r:
                data = r.read()
                status = r.status
        except (urllib.error.URLError, http.client.HTTPException, TimeoutError, OSError) as exc:
            print(f"{key}: ERROR {type(exc).__name__}: {str(exc)[:120]}")
            continue
        target.write_bytes(data)
        digest = sha256(data)
        match = "match" if spec.sha256 == digest else f"MISMATCH registry={spec.sha256[:16]}"
        print(f"{key}: http={status} bytes={len(data)} sha256={digest[:16]} {match}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())