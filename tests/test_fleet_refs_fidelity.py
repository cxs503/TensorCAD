"""Invariants of the reference-library fidelity / determinability measurement."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tensorcad.fleet.crossref import QUANTITIES
from tensorcad.fleet.fidelity import (
    CLASS_MULTI_AGREE,
    CLASS_MULTI_DISPUTED,
    CLASS_NO_DETERMINATION,
    CLASS_SINGLE_DOCUMENT,
    CLASS_TRANSCRIPTION,
    cluster,
    determinability,
    keyword_scan,
    summarise,
    transcription_fidelity,
)

ROOT = Path(__file__).resolve().parents[1]
LIBRARY = ROOT / "benchmarks" / "verified" / "fleet_refs" / "crossref.json"


def _entry(value, source="s1", family="FAM-A", basis="volume", raw=None):
    return {
        "family": family,
        "cell": {
            "quantity": "cb",
            "value": value,
            "raw": raw if raw is not None else str(value),
            "label": "Cb",
            "source": source,
            "unit_basis": basis,
        },
    }


def _library(entries, verdict, lineages):
    return {
        "anchor": "TEST",
        "report": {
            "HULL": {
                "quantities": {
                    q: {
                        "verdict": verdict if q == "cb" else "MISSING",
                        "independent_families": lineages if q == "cb" else 0,
                        "entries": entries if q == "cb" else [],
                    }
                    for q in QUANTITIES
                }
            }
        },
    }


def test_cluster_groups_agreeing_values_and_isolates_a_tenfold_typo():
    groups = cluster([58.0, 58.0, 58.0, 580.0, 580.0], 0.005)
    assert groups[0] == [58.0, 58.0, 58.0]
    assert groups[1] == [580.0, 580.0]


def test_transcription_spread_and_outlier_are_measured_not_asserted():
    entries = [
        _entry(0.6505, "a", "FAM-A"),
        _entry(0.651, "b", "FAM-B"),
        _entry(6.51, "c", "FAM-B", raw="6.51"),
    ]
    report = transcription_fidelity(_library(entries, "CROSS_VERIFIED", 2))
    cell = report["cells"]["HULL"]["cb"]
    assert cell["documents"] == 3
    assert cell["modal_value"] == 0.6505
    assert cell["outliers"][0]["source"] == "c"
    assert cell["outliers"][0]["rel_to_modal"] > 0.8
    assert cell["spread_relative"] == pytest.approx((6.51 - 0.6505) / 0.6505, rel=1e-6)


def test_a_duplicated_outlier_is_one_finding():
    entries = [
        _entry(58.0, "a", "FAM-A"),
        _entry(58.0, "a2", "FAM-A"),
        _entry(580.0, "b", "FAM-B", raw="580"),
        _entry(580.0, "b", "FAM-B", raw="580"),
    ]
    report = transcription_fidelity(_library(entries, "CROSS_VERIFIED", 2))
    assert len(report["cells"]["HULL"]["cb"]["outliers"]) == 1


def test_fallback_column_reading_never_votes():
    entries = [_entry(0.65), _entry(9.99, "s2", basis="fallback_column")]
    report = transcription_fidelity(_library(entries, "SINGLE_SOURCE", 1))
    assert report["cells"]["HULL"]["cb"]["documents"] == 1


@pytest.mark.parametrize(
    "entries,verdict,lineages,expected",
    [
        ([], "MISSING", 0, CLASS_NO_DETERMINATION),
        ([_entry(0.65)], "SINGLE_SOURCE", 1, CLASS_SINGLE_DOCUMENT),
        ([_entry(0.65), _entry(0.65, "s2", "FAM-A")], "SINGLE_SOURCE", 1, CLASS_TRANSCRIPTION),
        (
            [_entry(0.65), _entry(0.65, "s2", "FAM-B")],
            "CROSS_VERIFIED",
            2,
            CLASS_MULTI_AGREE,
        ),
        (
            [_entry(0.65), _entry(0.70, "s2", "FAM-B")],
            "CONFLICT",
            2,
            CLASS_MULTI_DISPUTED,
        ),
    ],
)
def test_determinability_classes(entries, verdict, lineages, expected):
    library = _library(entries, verdict, lineages)
    fidelity = transcription_fidelity(library)
    report = determinability(library, fidelity)
    assert report["ships"]["HULL"]["quantities"]["cb"]["class"] == expected


def test_transcription_and_single_document_still_block():
    library = _library(
        [_entry(0.65), _entry(0.65, "s2", "FAM-A"), _entry(0.65, "s3", "FAM-A")],
        "SINGLE_SOURCE",
        1,
    )
    report = determinability(library)
    hull = report["ships"]["HULL"]
    assert hull["blocking"] == list(QUANTITIES)
    assert hull["core_four_resolved"] is False


def test_keyword_scan_separates_reader_gap_from_source_gap():
    documents = {"doc": "Table 1: LCB 2.55 % Lpp and Cb 0.858 and displacement 178369 m3"}
    extracted = {"doc": {"JBC": ["cb"]}}
    scan = keyword_scan(documents, extracted)
    gaps = {row["quantity"] for row in scan["mention_without_extraction"]}
    assert "displacement_volume" in gaps and "lcb_percent" in gaps
    assert "cb" not in gaps
    rows = scan["sources"]["doc"]
    assert rows["cb"]["extracted"] is True and rows["displacement_volume"]["extracted"] is False


@pytest.mark.skipif(not LIBRARY.exists(), reason="reference library not built in this checkout")
def test_shipped_library_reproduces_the_readme_numbers():
    library = json.loads(LIBRARY.read_text())
    fidelity = transcription_fidelity(library)
    report = determinability(library, fidelity)
    summary = summarise(fidelity, report)
    assert summary["cells_total"] == 42
    #: The adjudicated 10x breadth typo and the ONRT mass/volume disagreement are findings of
    #: this library; if the measurement stops seeing them, the measurement is broken.
    assert {row["quantity"] for row in summary["outliers"]} == {"beam", "displacement_volume"}
    assert report["ships"]["KCS"]["core_four_resolved"] is True
    assert report["ships"]["Series60"]["counts"][CLASS_NO_DETERMINATION] == len(QUANTITIES)
