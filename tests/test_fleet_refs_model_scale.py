"""Tests for the published-model-scale cross-check (anchor REPRO-GEOM-FLEET-REFS).

The artefact under test prices a *derived* unit factor against a published statement, so the
tests are written to catch the two ways such a check can lie: reading the wrong hull's table,
and quietly rescaling the benchmark it is supposed to explain.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from tools.baseline_fleet import model_scale_refs as refs

LIBRARY = Path("benchmarks/verified/fleet_refs/model_scale.json")
SPECS = Path("benchmarks/verified/baseline_fleet/fleet_specs.json")
RESULT = Path("benchmarks/verified/baseline_fleet/result.json")
DOC = Path("docs/fleet_refs_jbc_scale.md")

#: The derived factor that ``fleet_specs.json`` admits for JBC; a change here is a change of
#: geometry, not of documentation, and must not happen as a side effect of a reference check.
JBC_DERIVED_SCALE = 39.72849598755336


@pytest.fixture(scope="module")
def library() -> dict:
    return json.loads(LIBRARY.read_text())


@pytest.fixture(scope="module")
def specs() -> dict:
    return {row["name"]: row for row in json.loads(SPECS.read_text())["ships"]}


@pytest.fixture(scope="module")
def result() -> dict:
    return json.loads(RESULT.read_text())


def test_every_claim_carries_a_quote_page_and_digest(library):
    claims = library["comparisons"]["JBC"]["claims"]
    assert claims, "the comparison states no claims at all"
    for name, cell in claims.items():
        assert re.fullmatch(r"[0-9a-f]{64}", cell["sha256"]), name
        assert cell["url"].startswith("https://"), name
        assert cell["doi"], name
        assert cell["page"] >= 1, name
        assert len(cell["quote"]) > 20, name
        # the quotation must actually contain the numeral it is supposed to address
        assert re.search(re.escape(cell["raw"]), cell["quote"]) or float(cell["raw"]) == cell["value"], name


def test_sha_of_the_claim_is_the_sha_of_the_document(library):
    document = library["documents"]["jbc_scale_pr"]
    for cell in library["comparisons"]["JBC"]["claims"].values():
        assert cell["sha256"] == document["sha256"]


def test_published_scale_is_read_from_the_jbc_table_not_the_kcs_table(library):
    """The preprint prints KCS (B=32.2) and JBC (B=45.0) tables back to back."""
    claims = library["comparisons"]["JBC"]["claims"]
    assert claims["beam"]["value"] == pytest.approx(45.0)
    assert claims["draft"]["value"] == pytest.approx(16.5)
    assert claims["displacement_volume"]["value"] == pytest.approx(178369.9, abs=0.05)
    assert claims["cb"]["value"] == pytest.approx(0.858)
    assert claims["model_length_pp"]["value"] == pytest.approx(7.00)
    assert claims["declared_scale_ratio"]["value"] == pytest.approx(40.0)


def test_page_numbers_are_real_pages(library):
    # two statements sit on different pages (table 3 on p.8, the sentence on p.7)
    claims = library["comparisons"]["JBC"]["claims"]
    assert claims["declared_scale_ratio"]["page"] != claims["beam"]["page"]
    assert library["documents"]["jbc_scale_pr"]["pages"] > 10


def test_scale_comparison_arithmetic_is_recomputable(library):
    comparison = library["comparisons"]["JBC"]
    assert comparison["status"] == "cross_checked"
    derived = comparison["derived_scale"]
    published = comparison["published_scale_ratio"]
    assert comparison["volume_factor_published_over_derived"] == pytest.approx((published / derived) ** 3)
    assert comparison["delta_percent"] == pytest.approx((published / derived - 1.0) * 100.0)


def test_benchmark_sensitivity_matches_result_json(library, result):
    comparison = library["comparisons"]["JBC"]["benchmark"]
    hull = result["hulls"]["JBC"]
    assert comparison["ref_volume_m3"] == pytest.approx(hull["ref_volume_m3"])
    assert comparison["got_volume_at_derived_scale_m3"] == pytest.approx(hull["got_volume_m3"])
    assert comparison["error_at_derived_scale_percent"] == pytest.approx(hull["rel_error_percent"], abs=1e-9)
    assert comparison["got_volume_at_published_scale_m3"] == pytest.approx(
        hull["got_volume_m3"] * library["comparisons"]["JBC"]["volume_factor_published_over_derived"]
    )
    assert comparison["points_explained_by_scale"] == pytest.approx(
        comparison["error_at_published_scale_percent"] - comparison["error_at_derived_scale_percent"]
    )


def test_the_scale_does_not_explain_the_withheld_verdict(library, result):
    """The tool must not be usable to excuse a geometry defect with a unit factor."""
    comparison = library["comparisons"]["JBC"]
    assert result["hulls"]["JBC"]["verdict"] == "WITHHELD"
    assert comparison["benchmark"]["verdict_unchanged"] == "WITHHELD"
    assert comparison["benchmark"]["residual_is_geometric"] is True
    assert abs(comparison["benchmark"]["error_at_published_scale_percent"]) > library["tolerance_percent"] * 5


def test_no_silent_rescaling_of_fleet_specs(specs, library):
    jbc = specs["JBC"]
    assert jbc["geometry_scale"] == pytest.approx(JBC_DERIVED_SCALE, rel=0, abs=1e-9)
    cross = jbc["scale_cross_check"]
    assert cross["derived_scale"] == pytest.approx(JBC_DERIVED_SCALE, rel=0, abs=1e-9)
    assert cross["published_scale_ratio"] == pytest.approx(library["comparisons"]["JBC"]["published_scale_ratio"])
    assert cross["evidence"], "fleet_specs must carry the quoted evidence, not just the verdict"


def test_the_document_is_rendered_from_the_artefact(library):
    text = DOC.read_text()
    comparison = library["comparisons"]["JBC"]
    assert f"{comparison['delta_percent']:+.3f} %" in text
    assert f"{comparison['benchmark']['error_at_published_scale_percent']:+.2f} %" in text
    assert comparison["published_scale_basis"] in text


# --------------------------------------------------------------------------------------
# parser behaviour on synthetic text (no network, no PDFs)


def _document(**overrides) -> refs.Document:
    base = dict(
        key="synthetic",
        ship="SHIP",
        publisher="p",
        url="https://example.org/x.pdf",
        doi="10.0000/x",
        licence="test",
        sha256="0" * 64,
        file="x.pdf",
    )
    base.update(overrides)
    return refs.Document(**base)


def test_unanchored_table_row_returns_the_first_hull_and_the_anchor_fixes_it():
    text = (
        "Table 1: KCS hull principal particulars\n"
        "Length, Lpp 230.0 7.2786 m\n"
        "Breadth, B 32.2 1.0190 m\n"
        "Table 3: JBC hull principal particulars\n"
        "Length, Lpp 280.0 7.00 m\n"
        "Breadth, B 45.0 1.125 m\n"
    )
    loose = refs.Claim(
        name="beam", label="B", unit="m", pattern=r"Breadth,\s*B\s*([0-9.]+)\s+[0-9.]+"
    )
    anchored = refs.Claim(
        name="beam",
        label="B",
        unit="m",
        pattern=r"Breadth,\s*B\s*([0-9.]+)\s+[0-9.]+",
        after=r"JBC hull principal particulars",
    )
    assert refs.extract(_document(claims=(loose,)), text)["beam"]["value"] == pytest.approx(32.2)
    assert refs.extract(_document(claims=(anchored,)), text)["beam"]["value"] == pytest.approx(45.0)


def test_quote_records_the_page_it_came_from():
    text = "page one\fpage two\nBreadth, B 45.0 1.125 m\n"
    claim = refs.Claim(name="beam", label="B", unit="m", pattern=r"Breadth,\s*B\s*([0-9.]+)")
    cell = refs.extract(_document(claims=(claim,)), text)["beam"]
    assert cell["page"] == 2
    assert "Breadth" in cell["quote"]


def test_disagreeing_statements_are_not_used():
    claims = {
        "declared_scale_ratio": {"value": 40.0, "source_claim": "declared_scale_ratio"},
        "model_length_pp": {"value": 7.00, "source_claim": "model_length_pp"},
        "full_length_pp": {"value": 284.0, "source_claim": "full_length_pp"},  # 284/7 = 40.57
    }
    spec = {"name": "JBC", "geometry_scale": 39.0, "scale_basis": "derived"}
    out = refs.scale_comparison(claims, spec, {})
    assert out["status"] == "no_published_scale_obtained"
    assert out["published_scale_ratio"] is None


def test_a_claim_the_document_does_not_state_yields_no_cell():
    claim = refs.Claim(name="lcb", label="LCB", unit="%", pattern=r"LCB\s*([0-9.]+)")
    assert refs.extract(_document(claims=(claim,)), "no such row here\f") == {}