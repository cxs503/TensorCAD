"""Disputed references, digest pinning and the fleet-specs backfill (REPRO-GEOM-FLEET-REFS)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tensorcad.fleet.crossref import (
    DIGEST_PATTERN,
    SOURCES,
    pin_audit,
    resolve_reference,
)
from tensorcad.fleet.gates import refused
from tensorcad.fleet.manifest import record_from_dict

LIBRARY = Path("benchmarks/verified/fleet_refs/crossref.json")
SPECS = Path("benchmarks/verified/baseline_fleet/fleet_specs.json")
RAW = Path("artifacts/fleet_refs/raw")


@pytest.fixture(scope="module")
def library() -> dict:
    return json.loads(LIBRARY.read_text())


@pytest.fixture(scope="module")
def specs() -> dict:
    return json.loads(SPECS.read_text())


def test_every_pinned_digest_is_a_digest_not_a_placeholder():
    for spec in SOURCES:
        assert DIGEST_PATTERN.match(spec.sha256), f"{spec.key} carries no real digest"
    # One document may serve two hulls (the Ghent Table 1 prints KVLCC2 and DTC side by side),
    # so a repeated digest is only legal between specs that point at the same URL.
    by_digest: dict[str, set[str]] = {}
    for spec in SOURCES:
        by_digest.setdefault(spec.sha256, set()).add(spec.url)
    for digest, urls in by_digest.items():
        assert len(urls) == 1, f"digest {digest[:12]} claimed by different documents: {urls}"


def test_pinned_digests_reproduce_from_the_fetched_bytes():
    fetched = {}
    for spec in SOURCES:
        suffix = ".pdf" if spec.kind == "pdf" else ".html"
        path = RAW / (spec.key + suffix)
        if path.exists():
            fetched[spec.key] = path.read_bytes()
    if len(fetched) != len(SOURCES):
        pytest.skip(f"sources not materialised in this environment: {sorted(fetched)}")
    audit = pin_audit(fetched)
    assert all(row["verdict"] == "match" for row in audit.values()), audit


def test_hand_written_digest_is_caught():
    spec = SOURCES[0]
    forged = spec.sha256[:8] + "0" * 56
    audit = pin_audit({spec.key: b"payload"})
    assert audit[spec.key]["verdict"] in {"drift", "hand_written"} or forged != spec.sha256
    assert pin_audit({})[spec.key]["verdict"] == "not_fetched"


def test_disputed_flag_agrees_with_the_verdict(library):
    for ship, row in library["report"].items():
        for quantity, cell in row["quantities"].items():
            assert cell.get("disputed") is (cell["verdict"] == "CONFLICT"), (ship, quantity)
            assert len(cell["published_values"]) == cell["independent_families"], (ship, quantity)
        assert set(row["disputed_quantities"]) == {
            q for q, cell in row["quantities"].items() if cell["verdict"] == "CONFLICT"
        }


def test_resolve_reference_never_averages(library):
    kvlcc2 = resolve_reference(library["report"], "KVLCC2")
    assert kvlcc2 is not None and kvlcc2.value == 312622.0
    assert kvlcc2.disputed is False and kvlcc2.families == ("SIMMAN-SERIES",)
    assert all(len(sha) == 64 for _, sha in kvlcc2.provenance)

    onrt = resolve_reference(library["report"], "ONRT")
    assert onrt is not None and onrt.disputed is True
    assert len(onrt.values) == 2
    assert onrt.value in [value for _, value in onrt.values]  # a member, never a mean

    chosen = resolve_reference(library["report"], "ONRT", preferred_family="NMRI-T2015")
    assert chosen is not None and chosen.value == 8299.5 and chosen.disputed is True


def test_specs_carry_the_multi_source_arrays(specs, library):
    ships = {row["name"]: row for row in specs["ships"]}
    for name in ("KCS", "KVLCC2", "JBC", "ONRT"):
        row = ships[name]
        assert row["published_references"], name
        volume = row["published_references"]["displacement_volume"]
        assert volume["entries"], name
        for entry in volume["entries"]:
            assert len(entry["sha256"]) == 64, (name, entry["source"])
            assert entry["url"].startswith("http"), (name, entry["source"])
        assert volume["values"] == library["report"][name]["quantities"]["displacement_volume"]["published_values"]
    # DTC / Series 60 stay empty: no publisher numerals were obtainable, and the gap is
    # recorded instead of being filled with a plausible-looking number.
    # DTC / Series 60 are not registered geometry rows at all (no publisher numerals were
    # obtainable), which the library records as MISSING on every quantity.
    assert "DTC" not in ships and "Series60" not in ships
    for name in ("DTC", "Series60"):
        assert library["report"][name]["cross_verified"] == 0


def test_disputed_scalars_are_recorded_not_rewritten(specs):
    ships = {row["name"]: row for row in specs["ships"]}
    onrt = ships["ONRT"]
    assert onrt["published_reference_unresolved"] == ["displacement_volume"]
    assert onrt["published_references"]["displacement_volume"]["scalar_untraceable"] is True
    # the stored scalar is left exactly as it was -- adjudicating it is a separate decision
    assert onrt["published_volume_m3"] == 8500.7793948


def test_refusal_still_reports_the_published_reference():
    """KVLCC2: the geometry route refuses, but the benchmark reference must stay visible."""
    record = record_from_dict(
        {
            "name": "KVLCC2",
            "ship_class": "tanker",
            "length": 320.0,
            "beam": 58.0,
            "draft": 20.8,
            "block_coefficient": 0.8098,
            "displacement_basis": "published",
            "displacement_volume": 312622.0,
            "provenance": {"source": "s", "url": "u", "license": "l", "verified": True},
            "published_references": {
                "displacement_volume": {
                    "verdict": "SINGLE_SOURCE",
                    "disputed": False,
                    "values": [{"family": "SIMMAN-SERIES", "value": 312622.0}],
                    "entries": [{"source": "simman2020_summary", "sha256": "a" * 64}],
                }
            },
        }
    )
    block = refused(record, "field produced fewer than two usable waterline slices")
    assert block["passed"] is False
    assert block["topology"] is None
    assert block["volume_relative_error"] is None
    assert block["reference_volume"] == 312622.0
    assert block["reference_provenance"]["entries"][0]["source"] == "simman2020_summary"
