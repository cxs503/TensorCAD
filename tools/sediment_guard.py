#!/usr/bin/env python3
"""Sediment-line guard: fail the build when a verified claim stops being true.

Usage::

    python tools/sediment_guard.py                    # verify, human-readable table
    python tools/sediment_guard.py --json out.json    # also write the machine report
    python tools/sediment_guard.py --entry wigley.cb_rational

Exit code 1 on any DRIFT / MISMATCH / MISSING / BROKEN row, 2 on a malformed
ledger.  A probe that raises is a failure, never a skip: a claim nobody can
re-derive any more is exactly what this guard exists to surface.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.knowledge import KnowledgeError, load_ledger, verify  # noqa: E402

DEFAULT_LEDGER = REPO_ROOT / "knowledge" / "sediment_line.json"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", default=str(DEFAULT_LEDGER), help="sediment ledger JSON")
    parser.add_argument("--root", default=str(REPO_ROOT), help="repository root the paths resolve against")
    parser.add_argument("--entry", action="append", default=None, help="verify only this entry id (repeatable)")
    parser.add_argument("--json", dest="json_out", default=None, help="write the machine report here")
    parser.add_argument("--quiet", action="store_true", help="only print failures and the summary")
    args = parser.parse_args(argv)

    try:
        ledger = load_ledger(args.ledger)
    except (OSError, ValueError, KnowledgeError) as exc:
        print(f"ledger error: {exc}", file=sys.stderr)
        return 2
    if args.entry:
        missing = [entry for entry in args.entry if entry not in ledger.ids()]
        if missing:
            print(f"unknown entry id(s): {missing}", file=sys.stderr)
            return 2
        ledger = type(ledger)(
            entries=tuple(e for e in ledger.entries if e.id in set(args.entry)),
            sha256=ledger.sha256,
            path=ledger.path,
            version=ledger.version,
        )

    report = verify(ledger, args.root)
    for row in report.rows:
        if args.quiet and row["status"] == "PASS":
            continue
        print(f"[{row['status']:8s}] {row['id']} ({row['kind']})")
        for check in row["checks"]:
            marker = "ok " if check["status"] == "PASS" else "!! "
            print(f"    {marker}{check['check']}: {check['detail']}")

    summary = report.as_dict()
    counts = summary["counts"]
    print(
        f"\nledger {ledger.path or args.ledger} sha256={ledger.sha256[:16]}... "
        f"entries={summary['entries']} checks={summary['checks']} "
        f"PASS={counts['PASS']} DRIFT={counts['DRIFT']} MISMATCH={counts['MISMATCH']} "
        f"MISSING={counts['MISSING']} BROKEN={counts['BROKEN']}"
    )
    print("verdict:", "PASS" if report.ok else "FAIL")

    if args.json_out:
        Path(args.json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.json_out).write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
        )
        print(f"report written to {args.json_out}")
    return 0 if report.ok else 1


if __name__ == "__main__":
    raise SystemExit(main())