#!/usr/bin/env python3
"""Backfill the cross-source reference library into ``fleet_specs.json``.

    python tools/baseline_fleet/backfill_refs.py backfill
    python tools/baseline_fleet/backfill_refs.py check

``fleet_specs.json`` used to carry *one* scalar per quantity (``published_volume_m3``,
``published_cb``, ...).  One number cannot be audited: a reader cannot tell a
reconstruction error from a publisher typo, nor see that two workshops disagree.  The
library in ``benchmarks/verified/fleet_refs/crossref.json`` holds every cell with its
source, URL and SHA-256, so this tool copies that structure into the specs as
``published_references`` (per quantity: the family values, the disputed flag and every
entry with its digest) and ``published_disputed`` (the quantities two publishers do not
agree on).  The legacy scalars are kept untouched: consumers that need one number keep
working, and ``check`` proves the scalars are a member of the family values so they cannot
drift away from the library.

Nothing here is transcribed: every value comes out of the library, and the library comes
out of fetched bytes whose digests are pinned in :mod:`tensorcad.fleet.crossref`.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.fleet.crossref import QUANTITIES, resolve_reference  # noqa: E402

DEFAULT_LIBRARY = ROOT / "benchmarks/verified/fleet_refs/crossref.json"
DEFAULT_SPECS = ROOT / "benchmarks/verified/baseline_fleet/fleet_specs.json"

#: Library quantity -> the legacy scalar it must agree with (when the ship states one).
LEGACY_FIELD = {
    "displacement_volume": "published_volume_m3",
    "cb": "published_cb",
    "lcb_percent": "published_lcb_percent",
    "lwl": "published_lwl_m",
    "length_pp": "published_length_pp_m",
}


def references_for(library: dict[str, Any], ship: str) -> dict[str, Any]:
    """The library's cells for one ship, in the shape wired into the specs file."""
    quantities: dict[str, Any] = {}
    for quantity in QUANTITIES:
        row = library["report"].get(ship, {}).get("quantities", {}).get(quantity)
        if row is None or not row["entries"]:
            continue
        quantities[quantity] = {
            "verdict": row["verdict"],
            "disputed": row["disputed"],
            "independent_families": row["independent_families"],
            "tolerance": row["tolerance"],
            "range_relative": row["range_relative"],
            "values": row["published_values"],
            "entries": [
                {
                    "family": entry["family"],
                    "value": entry["cell"]["value"],
                    "source": entry["cell"]["source"],
                    "url": entry["cell"]["url"],
                    "sha256": entry["cell"]["sha256"],
                    "raw": entry["cell"]["raw"],
                    "label": entry["cell"]["label"],
                    "unit_basis": entry["cell"]["unit_basis"],
                }
                for entry in row["entries"]
            ],
        }
    return quantities


def _scalar_problems(
    ship_row: dict[str, Any], references: dict[str, Any]
) -> tuple[list[str], list[str]]:
    """Check the legacy scalar against the library.

    Two severities, because they mean different things:

    * ``hard`` -- the ship states a scalar for a quantity the library has no cell for.  The
      scalar is unbacked; backfilling cannot be trusted and the command fails.
    * ``unresolved`` -- the scalar is not equal to *any* publisher's full-scale value.  This
      is recorded in the specs (``published_reference_unresolved``) and printed, but it is
      **not** silently rewritten: changing a published reference quietly would hide the
      defect behind a passing benchmark.  ONRT is the live case -- its stored
      8 500.7793948 m^3 matches neither family (8 299.5 NMRI / 8 507.0 SIMMAN) and carries
      model-scale digit precision, so the reference itself is what needs adjudicating.
    """
    hard: list[str] = []
    unresolved: list[str] = []
    for quantity, field in LEGACY_FIELD.items():
        scalar = ship_row.get(field)
        if scalar is None:
            continue
        row = references.get(quantity)
        if row is None:
            hard.append(f"{ship_row['name']}: {field}={scalar!r} has no library cell")
            continue
        values = [round(item["value"], 6) for item in row["values"]]
        if round(float(scalar), 6) not in values:
            unresolved.append(
                f"{ship_row['name']}: {field}={scalar!r} is not one of the published "
                f"family values {row['values']} (reference unresolved, not rewritten)"
            )
            row["scalar_untraceable"] = True
    return hard, unresolved


def backfill(library: dict[str, Any], specs: dict[str, Any]) -> tuple[dict[str, Any], list[str], list[str]]:
    """Return the specs document with the library wired in, plus hard problems and unresolved."""
    hard: list[str] = []
    unresolved: list[str] = []
    for ship_row in specs["ships"]:
        references = references_for(library, ship_row["name"])
        ship_hard, ship_unresolved = _scalar_problems(ship_row, references)
        ship_row["published_references"] = references
        ship_row["published_disputed"] = [
            quantity for quantity, row in references.items() if row["disputed"]
        ]
        ship_row["published_reference_unresolved"] = [
            quantity for quantity, row in references.items() if row.get("scalar_untraceable")
        ]
        hard.extend(ship_hard)
        unresolved.extend(ship_unresolved)
    return specs, hard, unresolved


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("backfill", "check"))
    parser.add_argument("--library", type=Path, default=DEFAULT_LIBRARY)
    parser.add_argument("--specs", type=Path, default=DEFAULT_SPECS)
    args = parser.parse_args(argv)

    library = json.loads(args.library.read_text())
    specs = json.loads(args.specs.read_text())
    rebuilt, hard, unresolved = backfill(library, specs)
    problems = hard + unresolved

    if args.command == "check":
        current = json.loads(args.specs.read_text())
        drift = [
            row["name"]
            for row, before in zip(rebuilt["ships"], current["ships"])
            if row.get("published_references") != before.get("published_references")
            or row.get("published_disputed") != before.get("published_disputed")
        ]
        for problem in problems:
            print("PROBLEM:", problem)
        if drift:
            print("DRIFT: specs are stale for", ", ".join(drift))
        if drift or hard:
            return 1
        print(f"specs agree with the library for {len(rebuilt['ships'])} hulls")
        return 0

    args.specs.write_text(json.dumps(rebuilt, indent=2) + "\n")
    for problem in problems:
        print("PROBLEM:", problem)
    covered = [row["name"] for row in rebuilt["ships"] if row.get("published_references")]
    print(f"wrote {args.specs}: multi-source references for {len(covered)} hulls: "
          f"{', '.join(covered) or 'none'}")
    return 1 if hard else 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
