#!/usr/bin/env python3
"""Baseline fleet gate CLI -- the one command that replaces "the file exists".

Usage
-----
    python tools/baseline_fleet/audit.py audit <geometry_dir> [--published published.json]
    python tools/baseline_fleet/audit.py sources <artifacts/baseline_fleet/sources.json>

``audit`` grades every ``*_Hull.stl`` in a directory (placeholder / watertight /
volume-vs-published).  ``sources`` re-checks the provenance registry: every
archive must carry source URL, publisher, licence, SHA-256 and retrieval time,
and the SHA-256 is recomputed from the bytes on disk -- a manifest that disagrees
with its own file fails here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.fleet_gate import audit_fleet, validate_source, write_json  # noqa: E402


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def cmd_audit(args: argparse.Namespace) -> int:
    published = json.loads(Path(args.published).read_text()) if args.published else {}
    report = audit_fleet(args.directory, published_volumes=published, tolerance=args.tolerance)
    if args.out:
        write_json(report, args.out)
    for ship in report["ships"]:
        mesh = ship["mesh"] or {}
        print(
            f"{ship['ship']:<10} {ship['grade']:<12} tris={mesh.get('triangles', '-'):>7} "
            f"vol={mesh.get('volume_m3', '-')!s:>12} " + "; ".join(ship["reasons"])
        )
    print("counts:", report["counts"], "pass:", report["pass"])
    return 0 if args.soft else (0 if report["pass"] else 1)


def cmd_sources(args: argparse.Namespace) -> int:
    registry = json.loads(Path(args.registry).read_text())
    problems = 0
    for entry in registry.get("archives", []):
        if "error" in entry:
            print(f"{entry.get('ship'):<12} FETCH-FAILED {entry['error']}")
            problems += 1
            continue
        issues = validate_source(entry)
        local = Path(entry.get("local_path", ""))
        if local.is_file():
            actual = _sha256(local)
            if actual != entry.get("sha256"):
                issues.append(f"sha256 on disk {actual[:16]} != manifest {str(entry.get('sha256'))[:16]}")
        else:
            issues.append(f"archive missing at {local}")
        state = "OK" if not issues else "PROBLEM"
        print(f"{entry['ship']:<12} {state:<8} {'; '.join(issues) or entry['source_url']}")
        problems += bool(issues)
    for key, page in registry.get("pages", {}).items():
        if "error" in page:
            print(f"{key:<12} FETCH-FAILED {page['error']}")
            problems += 1
    print(f"provenance problems: {problems}")
    return 0 if problems == 0 else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    audit = sub.add_parser("audit", help="grade a directory of *_Hull.stl")
    audit.add_argument("directory")
    audit.add_argument("--published", help="JSON mapping ship -> published displacement volume (m3)")
    audit.add_argument("--tolerance", type=float, default=0.01)
    audit.add_argument("--out", help="write the JSON report here")
    audit.add_argument("--soft", action="store_true", help="exit 0 even when the fleet fails")
    audit.set_defaults(func=cmd_audit)

    sources = sub.add_parser("sources", help="re-verify the provenance registry")
    sources.add_argument("registry")
    sources.set_defaults(func=cmd_sources)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())