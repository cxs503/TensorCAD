#!/usr/bin/env python3
"""Adopt a derived model scale from its record into ``fleet_specs.json``.

    python tools/baseline_fleet/adopt_derived_scale.py --check   # is the manifest in sync?
    python tools/baseline_fleet/adopt_derived_scale.py           # write the derived scale

Why this tool exists
--------------------
``build_fleet.py`` fills ``geometry_scale`` from a *published* scale row when one
exists, and otherwise records ``None``, which makes
:func:`tensorcad.fleet_baseline.build_ship` normalise the hull's **x-extent** to the
published ``L_PP``.  For JBC the x-extent is the deck-level extreme (bow flare plus
stern overhang): 7.5428 m against ``L_PP = 280 m``, so that fit builds the ship
7.0702 / 7.5428 = 93.7 % of its true length -- a systematic **-6.3 %** geometry
bias that no volume gate can see, because the fit also rescales the volume.

The Tokyo-2015 page publishes no model-scale row, so the scale cannot be scraped;
``tools/baseline_fleet/jbc_model_scale.py`` derives it from the publisher's own CAD
instead (1:39.7285, half-spread 0.18 %, two independent published lengths).  This tool
carries that number from the record into the fleet manifest, and nothing else:

* the scale is read from the record -- no value is ever typed into this file;
* ``scale_basis`` keeps the record path plus the digests of the publisher page and
  archive the derivation rests on, so the manifest entry is a pointer into a
  reproducible derivation;
* ``--check`` is what the test suite runs: it fails when the manifest and the record
  disagree, which is how a regeneration with a stale policy (or a hand edit) is caught.

Only ``fleet_specs.json`` is touched; the geometry bytes stay out of the repository
because the publisher states no licence.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

DEFAULT_SPECS = REPO / "benchmarks" / "verified" / "baseline_fleet" / "fleet_specs.json"
DEFAULT_RECORD = REPO / "artifacts" / "baseline_fleet" / "jbc_model_scale.json"
HULL = "JBC"


def derived_scale(record_path: Path) -> tuple[float, str]:
    """Return ``(scale, basis)`` from the record, refusing anything unrecognised."""
    record = json.loads(record_path.read_text())
    if record.get("schema") != "tensorcad.fleet.model_scale.v1":
        raise SystemExit(f"{record_path}: unexpected schema {record.get('schema')!r}")
    scale = float(record["scale_derived"]["value"])
    spread = record["scale_derived"].get("relative_half_spread")
    published = record.get("published_full_scale", {})
    page = record["sources"]["particulars_page"]["sha256"][:16]
    archive = record["sources"]["iges_archive"]["sha256"][:16]
    basis = (
        "model scale derived from the publisher CAD (no scale row is published):"
        f" B={published['beam_m']['value']:g} m / D={published['depth_m']['value']:g} m"
        f" over the measured model spans, half-spread {spread:.2%};"
        f" record {record_path.relative_to(REPO)}"
        f" (page sha256 {page}…, archive sha256 {archive}…)"
    )
    return scale, basis


def sync(specs_path: Path, record_path: Path, *, write: bool) -> int:
    scale, basis = derived_scale(record_path)
    document = json.loads(specs_path.read_text())
    hits = [ship for ship in document.get("ships", []) if ship.get("name") == HULL]
    if len(hits) != 1:
        print(f"FAIL {specs_path}: expected exactly one {HULL} entry, found {len(hits)}")
        return 1
    ship = hits[0]
    current = ship.get("geometry_scale")
    in_sync = current is not None and abs(float(current) - scale) <= 1e-9
    if in_sync:
        print(f"OK {HULL} geometry_scale {scale:.5f} matches the record")
        return 0
    if not write:
        print(
            f"FAIL {HULL} geometry_scale is {current!r}, record says {scale:.5f}:"
            f" the manifest does not use the derived scale"
        )
        return 1
    ship["geometry_scale"] = scale
    ship["scale_basis"] = basis
    specs_path.write_text(json.dumps(document, indent=2) + "\n")
    print(f"wrote {specs_path}: {HULL} geometry_scale {scale:.5f}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--specs", default=str(DEFAULT_SPECS))
    parser.add_argument("--record", default=str(DEFAULT_RECORD))
    parser.add_argument("--check", action="store_true",
                        help="verify the manifest already carries the derived scale")
    args = parser.parse_args(argv)
    return sync(Path(args.specs), Path(args.record), write=not args.check)


if __name__ == "__main__":
    raise SystemExit(main())