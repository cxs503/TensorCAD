#!/usr/bin/env python3
"""Audit the pinned source digests against the bytes that were actually read.

    python tools/baseline_fleet/refs_audit.py --raw artifacts/fleet_refs/raw \
        --out benchmarks/verified/fleet_refs/digest_audit.json

The first version of the registry pinned four of its seven digests as the empty string and
two with hand-written continuations (they shared the leading eight hex digits with the real
bytes and then diverged).  A digest that is typed rather than fetched is not a fingerprint,
so this audit recomputes every digest from the payload and fails on anything but ``match``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.fleet.crossref import SOURCES, pin_audit  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, default=ROOT / "artifacts/fleet_refs/raw")
    parser.add_argument("--out", type=Path,
                        default=ROOT / "benchmarks/verified/fleet_refs/digest_audit.json")
    args = parser.parse_args(argv)

    fetched = {}
    for spec in SOURCES:
        suffix = ".pdf" if spec.kind == "pdf" else ".html"
        path = args.raw / (spec.key + suffix)
        if path.exists():
            fetched[spec.key] = path.read_bytes()
    audit = pin_audit(fetched)
    payload = {
        "schema": "tensorcad.fleet_refs.digest_audit/1.0",
        "sources": audit,
        "pinned": sum(1 for row in audit.values() if row["verdict"] == "match"),
        "checked": len(fetched),
        "all_pinned_digests_match_fetched_bytes": all(
            row["verdict"] == "match" for row in audit.values()
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {args.out}: {payload['pinned']}/{len(audit)} match, "
          f"all={payload['all_pinned_digests_match_fetched_bytes']}")
    for key, row in audit.items():
        if row["verdict"] != "match":
            print(f"  {row['verdict']}: {key}")
    return 0 if payload["all_pinned_digests_match_fetched_bytes"] else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
