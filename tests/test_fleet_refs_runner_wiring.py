"""The runner must see the multi-source reference library (REPRO-GEOM-FLEET-REFS).

The cross-source library (``benchmarks/verified/fleet_refs/crossref.json``) is
backfilled into ``fleet_specs.json`` by ``tools/baseline_fleet/backfill_refs.py``,
but the *runner* reads ``tensorcad/fleet/data/fleet_manifest.json``.  Until the
manifest rows carried the same ``published_references`` cells, the benchmark
compared every hull against a recalled scalar with no provenance -- the exact
failure mode the anchor exists to close.  These tests pin the wiring:

* every manifest row that has a library cell in the specs carries it verbatim;
* the KVLCC2 312 622 m^3 reference resolves through the runner with its
  family values and per-source digests (acceptance item 5 of the anchor);
* a CONFLICT (ONRT) is surfaced as ``disputed`` and never averaged away.
"""

from __future__ import annotations

import json
from pathlib import Path

from tensorcad.fleet.manifest import load_manifest

REPO = Path(__file__).resolve().parent.parent
SPECS = REPO / "benchmarks" / "verified" / "baseline_fleet" / "fleet_specs.json"
MANIFEST = REPO / "tensorcad" / "fleet" / "data" / "fleet_manifest.json"


def _specs_rows() -> dict[str, dict]:
    return {row["name"]: row for row in json.loads(SPECS.read_text())["ships"]}


def test_manifest_rows_carry_the_specs_reference_cells():
    manifest = load_manifest(MANIFEST)
    specs = _specs_rows()
    for record in manifest.ships:
        spec_row = specs.get(record.name)
        if spec_row is None or not spec_row.get("published_references"):
            assert not record.published_references(), record.name
            continue
        assert record.published_references() == spec_row["published_references"], (
            f"{record.name}: manifest published_references drifted from fleet_specs.json"
        )


def test_kvlcc2_reference_resolves_with_provenance():
    """Acceptance 5: the KVLCC2 312 622 m^3 reference is wired into the runner."""
    record = next(s for s in load_manifest(MANIFEST).ships if s.name == "KVLCC2")
    assert record.reference_volume() == 312622.0
    prov = record.published_provenance("displacement_volume")
    assert prov is not None
    assert prov["verdict"] == "CROSS_VERIFIED"
    assert prov["disputed"] is False
    families = {value["family"] for value in prov["values"]}
    # two independent lineages: the SIMMAN series and the Ghent JSR 2017 paper
    assert families == {"SIMMAN-SERIES", "JSR-2017"}
    assert all(len(entry["sha256"]) == 64 for entry in prov["entries"])
    assert {entry["source"] for entry in prov["entries"]} >= {
        "simman2014_kvlcc2",
        "simman2020_summary",
        "kvlcc2_ugent",
    }
    # the runner's reference is a *member* of the published family values, not a mean
    assert record.reference_volume() in {value["value"] for value in prov["values"]}


def test_onrt_conflict_is_surfaced_not_averaged():
    record = next(s for s in load_manifest(MANIFEST).ships if s.name == "ONRT")
    prov = record.published_provenance("displacement_volume")
    assert prov is not None and prov["disputed"] is True
    values = {value["value"] for value in prov["values"]}
    assert values == {8299.5, 8507.0}
    # the stored runner reference is one of the published values or flagged unresolved
    # upstream; it must never be the mean of the two families.
    assert record.reference_volume() != (8299.5 + 8507.0) / 2.0


def test_ships_without_library_cells_stay_empty():
    manifest = load_manifest(MANIFEST)
    specs = _specs_rows()
    for record in manifest.ships:
        spec_row = specs.get(record.name)
        if spec_row is None or not spec_row.get("published_references"):
            assert record.published_provenance("displacement_volume") is None, record.name
