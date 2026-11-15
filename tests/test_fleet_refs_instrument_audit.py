"""The reference instrument against a designed truth (anchor REPRO-GEOM-FLEET-REFS).

These tests are the regression fence around the audit in
``tools/baseline_fleet/refs_instrument_audit.py``: every fixture there is a synthetic table whose
correct reading is written down, so a wrong row or a wrong unit fails here by name instead of
passing as a plausible numeral.  The pre-guard reading is asserted too -- a fix that cannot be
switched off cannot be shown to have fixed anything.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tensorcad.fleet import crossref as xref
from tensorcad.fleet.crossref import (
    _LCB_AS_LENGTH,
    _MASS_LABEL,
    _RATIO_LABEL,
    Cell,
    SourceSpec,
    cross_reference,
    parse_source,
    votes,
)
from tools.baseline_fleet import refs_instrument_audit as audit

ARTEFACT = Path("artifacts/fleet_refs/instrument_audit.json")
PREFIX = Path("artifacts/fleet_refs/instrument_audit_prefix.json")
DOC = Path("docs/fleet_refs_instrument_audit.md")


def _spec(key: str = "k", family: str = "F", column: int | None = None) -> SourceSpec:
    return SourceSpec(
        key=key,
        ship="SYNTH",
        publisher="synthetic",
        family=family,
        url="fixture://k",
        kind="html",
        sha256="0" * 64,
        licence="synthetic fixture",
        note="",
        column=column,
    )


def _cells(html: str, **kwargs) -> dict[str, Cell]:
    return parse_source(_spec(**kwargs), html)


# ---------------------------------------------------------------------------------------
# the corpus, as shipped
# ---------------------------------------------------------------------------------------
def test_corpus_passes_with_guards() -> None:
    artefact = audit.run_corpus()
    failing = {row["id"]: row["problems"] for row in artefact["fixtures"] if not row["passed"]}
    assert failing == {}
    assert artefact["sensitivity"]["passed"] == artefact["sensitivity"]["fixtures"] >= 5
    assert artefact["specificity"]["passed"] == artefact["specificity"]["fixtures"] >= 3


def test_corpus_catches_the_pre_fix_reading() -> None:
    """With the guards off the same corpus must fail on exactly the defects they address."""
    guards = xref.INSTRUMENT_GUARDS
    xref.INSTRUMENT_GUARDS = False
    try:
        artefact = audit.run_corpus()
    finally:
        xref.INSTRUMENT_GUARDS = guards
    failing = {row["id"] for row in artefact["fixtures"] if not row["passed"]}
    assert {
        "ratio_rows_present",
        "lcb_printed_in_metres",
        "lcb_metres_only",
        "mass_labelled_volume",
    } <= failing
    # controls stay clean either way: the guards must not be doing the work of the tolerance rules
    assert artefact["specificity"]["passed"] == artefact["specificity"]["fixtures"]


def test_artefact_and_prefix_are_present_and_agree_with_the_code() -> None:
    assert ARTEFACT.exists() and PREFIX.exists() and DOC.exists()
    shipped = json.loads(ARTEFACT.read_text())
    assert shipped["guards"] is True
    assert [row["id"] for row in shipped["fixtures"] if not row["passed"]] == []
    prefix = json.loads(PREFIX.read_text())
    assert prefix["guards"] is False
    assert prefix["sensitivity"]["passed"] < prefix["sensitivity"]["fixtures"]


def test_guards_do_not_reclassify_the_shipped_library() -> None:
    """A guard that changes a published verdict must be reported, not shipped quietly."""
    impact = audit.impact_on_library()
    assert impact["checked"] is True
    assert impact["html_cells_scanned"] > 0
    assert impact["reclassified"] == 0
    assert impact["ratio_label_voters"] == []
    assert impact["mass_label_volume_cells"] == []


# ---------------------------------------------------------------------------------------
# one behaviour per test, reading the guards directly
# ---------------------------------------------------------------------------------------
@pytest.mark.parametrize(
    "label",
    ["Lwl/Lpp", "Lwl / Lpp", "Bwl/Lpp", "T/L", "L/B", "B/T", "Cb/Cm"],
)
def test_quotient_labels_are_recognised(label: str) -> None:
    assert _RATIO_LABEL.search(label) is not None


@pytest.mark.parametrize("label", ["kg/m3", "Wetted surface Sw (m2)", "TFP / TAP (m)", "Draft T (m)", "Lpp (m)"])
def test_ordinary_labels_are_not_quotients(label: str) -> None:
    assert _RATIO_LABEL.search(label) is None


def test_ratio_row_is_skipped_and_the_real_row_wins() -> None:
    html = audit.table(
        (("Lwl/Lpp", "1.0109"), ("T/L", "0.0464"), ("Length on waterline Lwl (m)", "232.5"), ("Design draft T (m)", "10.8"))
    )
    cells = _cells(html)
    assert cells["lwl"].value == 232.5
    assert cells["draft"].value == 10.8
    assert cells["lwl"].unit_basis == "volume"


def test_lcb_in_metres_never_votes_and_a_percentage_row_wins() -> None:
    html = audit.table((("LCB (m)", "-3.395"), ("LCB (%), fwd+", "-1.48")))
    cells = _cells(html)
    assert cells["lcb_percent"].value == -1.48
    assert votes(cells["lcb_percent"]) is True
    alone = _cells(audit.table((("LCB (m)", "-3.395"),)))
    assert alone["lcb_percent"].value == -3.395
    assert alone["lcb_percent"].unit_basis == "lcb_as_length"
    assert votes(alone["lcb_percent"]) is False
    assert _LCB_AS_LENGTH.search("LCB (m)") is not None


def test_mass_stated_in_the_label_is_converted_not_read_as_volume() -> None:
    html = audit.table((("Displacement (t)", "52,030"),))
    cell = _cells(html)["displacement_volume"]
    assert cell.unit_basis == "mass"
    assert cell.value == pytest.approx(52030 / 1.025, abs=0.05)  # rounded to 0.1 m3
    assert _MASS_LABEL.search("Displacement (t)") is not None
    assert _MASS_LABEL.search("Displacement volume (m3)") is None


def test_model_only_volume_is_recovered_with_the_cube_of_the_table_ratio() -> None:
    html = (
        "<table><tr><th>Particular</th><th>Model (1:100)</th><th>Full scale</th></tr>"
        "<tr><td>Length between perpendiculars Lpp (m)</td><td>2.300</td><td>230.0</td></tr>"
        "<tr><td>Design draft T (m)</td><td>0.108</td><td>10.8</td></tr>"
        "<tr><td>Displacement volume (m3)</td><td>0.05203</td><td>NA</td></tr></table>"
    )
    cell = _cells(html, column=2)["displacement_volume"]
    assert cell.unit_basis == "derived_model_column_cubed"
    assert cell.value == pytest.approx(52030.0, rel=1e-9)
    assert votes(cell) is False
    assert "(230 m / 2.3 m)^3" in cell.derivation


def test_derived_and_length_cells_never_vote() -> None:
    for basis in ("fallback_column", "derived_model_column_scaled", "lcb_as_length"):
        cell = Cell(
            quantity="lwl",
            value=1.0,
            raw="1",
            label="x",
            source="s",
            url="u",
            sha256="0" * 64,
            column=1,
            unit_basis=basis,
        )
        assert votes(cell) is False