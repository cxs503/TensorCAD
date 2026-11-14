"""Tests for the cross-source reference library (anchor REPRO-GEOM-FLEET-REFS)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tensorcad.fleet.crossref import (
    QUANTITIES,
    SOURCES,
    Cell,
    SourceSpec,
    closure_check,
    coverage_summary,
    cross_reference,
    parse_source,
)

LIBRARY = Path("benchmarks/verified/fleet_refs/crossref.json")
DOC = Path("docs/fleet_refs_crossref.md")


@pytest.fixture(scope="module")
def library() -> dict:
    return json.loads(LIBRARY.read_text())


def _spec(**kwargs) -> SourceSpec:
    base = dict(
        key="synthetic",
        ship="SHIP",
        publisher="p",
        family="F",
        url="u",
        kind="html",
        sha256="0" * 64,
        licence="test",
    )
    base.update(kwargs)
    return SourceSpec(**base)


def _cell(quantity: str, value: float, source: str, basis: str = "volume") -> Cell:
    return Cell(
        quantity=quantity,
        value=value,
        raw=str(value),
        label=quantity,
        source=source,
        url="u",
        sha256="0" * 64,
        column=1,
        unit_basis=basis,
    )


# --------------------------------------------------------------------------- parser behaviour


def test_label_comes_from_the_first_cell_only():
    """A derivative row ("izz/Lpp 0.25") must not be read as Lpp."""
    markup = (
        "<table>"
        "<tr><td>Length between perpendiculars Lpp</td><td>320.0</td></tr>"
        "<tr><td>izz/Lpp</td><td>0.25</td></tr>"
        "</table>"
    )
    cells = parse_source(_spec(), markup)
    assert cells["length_pp"].value == 320.0


def test_fallback_column_is_flagged_and_does_not_vote():
    """When the configured full-scale column is empty, the neighbour is kept but never votes."""
    markup = (
        "<table>"
        "<tr><td>Displacement volume</td><td></td><td>0.0726</td></tr>"
        "<tr><td>Block coefficient</td><td></td><td>0.535</td></tr>"
        "</table>"
    )
    spec = _spec(column=1)  # column 1 is empty on this snippet, column 2 holds the model value
    cells = parse_source(spec, markup)
    flagged = [c for c in cells.values() if c.unit_basis.startswith("fallback_column")]
    assert flagged, "a neighbour column must be flagged, not silently adopted"

    report = cross_reference({"SHIP": {"synthetic": cells}})
    assert report["SHIP"]["quantities"]["displacement_volume"]["verdict"] == "MISSING"


def test_mass_reference_is_converted_with_the_stated_density():
    markup = (
        "<table>"
        "<tr><td>Displacement</td><td>8507 t</td></tr>"
        "<tr><td>Block coefficient</td><td>0.535</td></tr>"
        "</table>"
    )
    cells = parse_source(_spec(), markup)
    assert cells["displacement_volume"].value == pytest.approx(8507 / 1.025, abs=0.1)
    assert "converted" in cells["displacement_volume"].note


# --------------------------------------------------------------------------- classification


def test_two_documents_from_one_lineage_are_not_two_sources():
    """Same lineage twice = SINGLE_SOURCE, which is the failure this anchor exists to expose."""
    cells = {q: _cell(q, 1.0, "a") for q in ("displacement_volume", "cb", "lcb_percent", "lwl")}
    # two documents of the same workshop series (SIMMAN 2014 + 2020) -- one lineage
    report = cross_reference({"SHIP": {"simman2014_kcs": cells, "simman2020_summary": cells}})
    assert report["SHIP"]["quantities"]["cb"]["verdict"] == "SINGLE_SOURCE"
    assert report["SHIP"]["quantities"]["cb"]["independent_families"] == 1


def test_conflict_is_reported_not_averaged():
    left = {q: _cell(q, 100.0, "a") for q in ("displacement_volume", "cb", "lcb_percent", "lwl")}
    right = {q: _cell(q, 110.0, "simman2020_summary") for q in left}
    report = cross_reference({"SHIP": {"a": left, "simman2020_summary": right}})
    cell = report["SHIP"]["quantities"]["displacement_volume"]
    assert cell["verdict"] == "CONFLICT"
    assert cell["family_values"] == {"F": 100.0, "SIMMAN-SERIES": 110.0}


# --------------------------------------------------------------------------- the built library


def test_every_source_is_digest_pinned_the_same_way_after_a_fetch(library):
    unfetched = [key for key, row in library["sources"].items() if row.get("status") != "read"]
    if unfetched:
        pytest.skip(f"sources not materialised in this environment: {unfetched}")
    for key, row in library["sources"].items():
        assert len(row["sha256"]) == 64, key


def test_cross_verified_cells_have_two_lineages_and_agree(library):
    for ship, row in library["report"].items():
        for quantity in QUANTITIES:
            cell = row["quantities"][quantity]
            if cell["verdict"] != "CROSS_VERIFIED":
                continue
            assert cell["independent_families"] >= 2, (ship, quantity)
            assert cell["range_relative"] <= cell["tolerance"], (ship, quantity)


def test_kcs_is_fully_cross_verified(library):
    kcs = library["report"]["KCS"]["quantities"]
    for quantity in QUANTITIES:
        assert kcs[quantity]["verdict"] == "CROSS_VERIFIED", quantity
    assert set(kcs["displacement_volume"]["family_values"]) == {"NMRI-T2015", "SIMMAN-SERIES"}


def test_kvlcc2_collapses_to_one_lineage_and_exposes_the_bwl_typo(library):
    kvlcc2 = library["report"]["KVLCC2"]["quantities"]
    # Two SIMMAN editions state 312622 m3, but they are one lineage: honest single source.
    assert kvlcc2["displacement_volume"]["verdict"] == "SINGLE_SOURCE"
    assert kvlcc2["displacement_volume"]["family_values"] == {"SIMMAN-SERIES": 312622.0}
    closure = library["closure"]["KVLCC2"]["SIMMAN-SERIES"]
    assert not closure["consistent"]
    # Cb*Lpp*Bwl*T with the printed Bwl=580 is 10x the published volume; the implied breadth
    # is the ship's real 58.0 m, i.e. the table carries a factor-of-ten typo.
    assert closure["implied_beam_m"] == pytest.approx(58.0, abs=0.5)
    assert closure["rel_error"] > 8.0


def test_honest_gaps_are_recorded(library):
    """JBC/ONRT/DTC/Series60 must be visibly below the acceptance bar, not quietly rounded up."""
    for ship in ("JBC", "DTC", "Series60"):
        assert ship in library["coverage"]["ships_failing_min_two_sources"]
    assert library["report"]["DTC"]["cross_verified"] == 0
    assert library["report"]["Series60"]["cross_verified"] == 0
    assert library["report"]["JBC"]["quantities"]["displacement_volume"]["verdict"] == "SINGLE_SOURCE"
    # ONRT disagrees between publishers on displacement (mass vs volume) and on LCB.
    assert library["report"]["ONRT"]["quantities"]["displacement_volume"]["verdict"] == "CONFLICT"
    assert library["report"]["ONRT"]["quantities"]["lcb_percent"]["verdict"] == "CONFLICT"


def test_missing_reconstructions_are_not_quoted_as_errors(library):
    """A model-scale or absent reconstruction must not be turned into a percentage error."""
    for ship, row in library["comparisons"].items():
        if row.get("comparability") != "full_scale":
            assert row["rel_error_percent"] is None, ship


def test_closure_uses_each_sources_own_digits(library):
    closure = closure_check(library["report"])
    assert closure["KCS"]["NMRI-T2015"]["consistent"]
    assert closure["KCS"]["SIMMAN-SERIES"]["consistent"]


def test_report_is_rendered_from_the_library(library):
    if not DOC.exists():
        pytest.skip("report not rendered")
    text = DOC.read_text()
    cov = coverage_summary(library["report"])
    assert f"{cov['cross_verified']}/{cov['cells']}" in text
    assert "implied Bwl 58.00 m" in text


def test_registry_covers_all_six_hulls():
    declared = {ship for spec in SOURCES for ship in spec.ship.split("|")}
    assert {"KCS", "KVLCC2", "JBC", "ONRT", "DTC", "Series60"} <= declared