"""Closure-based adjudication of disagreeing / defective published cells.

Anchor ``REPRO-GEOM-FLEET-REFS``.  A reference library that only *reports* a disagreement
leaves the reader to guess; one that picks the value closest to the reconstruction flatters
the code under test.  :func:`tensorcad.fleet.crossref.adjudicate` does neither: it asks
whether each publisher's own table closes on its own numbers (``∇ = Cb·L·B·T``) and it never
rewrites a printed cell.

The tests below pin the four behaviours the library claims:

* the ONRT displacement dispute is decided in favour of the publisher whose table closes
  (SIMMAN states 8 507 m³ and closes to 0.07 %; the NMRI cell states 8 507 **t**, which at
  1.025 t/m³ is 8 299.5 m³ and misses by 2.4 %);
* the KVLCC2 ``Bwl 580`` cell is attributed to a factor-of-ten slip, recorded with the factor
  and the corrected residual, and the printed number is left untouched;
* every family that closes is marked consistent and untouched, so the machinery is inert on a
  healthy hull (KCS);
* a hull with no published cells at all (DTC, Series 60) produces no adjudication, no
  exception and no invented value.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tensorcad.fleet.crossref import (
    CLOSURE_TOLERANCE,
    QUANTITIES,
    Cell,
    SourceSpec,
    adjudicate,
    coverage_summary,
    cross_reference,
    parse_source,
)

REPO = Path(__file__).resolve().parents[1]
LIBRARY = REPO / "benchmarks/verified/fleet_refs/crossref.json"


@pytest.fixture(scope="module")
def library() -> dict:
    return json.loads(LIBRARY.read_text())


@pytest.fixture(scope="module")
def verdict(library: dict) -> dict:
    return adjudicate(library["report"])


def _cell(quantity: str, value: float, source: str, family: str) -> dict:
    return {
        "family": family,
        "cell": Cell(
            quantity=quantity,
            value=value,
            raw=str(value),
            label=quantity,
            source=source,
            url=f"https://example.invalid/{source}",
            sha256="0" * 64,
            column=1,
        ).as_dict(),
    }


def _synthetic(rows: dict[str, dict[str, tuple[float, str]]]) -> dict:
    """Build a cross-reference report from ``{ship: {quantity: (value, family)}}``."""
    extractions: dict[str, dict[str, dict[str, Cell]]] = {}
    for ship, quantities in rows.items():
        per_source: dict[str, dict[str, Cell]] = {}
        for quantity, (value, family) in quantities.items():
            key = f"{ship}_{family}"
            per_source.setdefault(key, {})[quantity] = Cell(
                quantity=quantity,
                value=value,
                raw=str(value),
                label=quantity,
                source=key,
                url=f"https://example.invalid/{key}",
                sha256="0" * 64,
                column=1,
            )
        extractions[ship] = per_source
    return cross_reference(extractions)


def test_onrt_displacement_conflict_is_resolved_by_closure(verdict: dict) -> None:
    """The publisher whose own table closes wins; the other keeps its printed value."""
    row = verdict["ships"]["ONRT"]["conflicts"]["displacement_volume"]
    assert row["resolution"] == "resolved_by_closure"
    assert row["preferred_family"] == "SIMMAN-SERIES"
    residuals = row["closure_residual_by_family"]
    assert abs(residuals["SIMMAN-SERIES"]) <= CLOSURE_TOLERANCE
    assert abs(residuals["NMRI-T2015"]) > CLOSURE_TOLERANCE
    assert row["family_values"]["NMRI-T2015"] == pytest.approx(8299.5, abs=1e-6)
    assert row["family_values"]["SIMMAN-SERIES"] == pytest.approx(8507.0, abs=1e-6)


def test_onrt_closure_uses_the_length_the_publisher_states(verdict: dict) -> None:
    """NMRI's ONRT table prints no Lpp; the identity must close on Lwl, not be skipped."""
    closure = verdict["ships"]["ONRT"]["closure"]
    assert set(closure) == {"NMRI-T2015", "SIMMAN-SERIES"}
    assert {row["length_basis"] for row in closure.values()} <= {"lwl", "length_pp"}


def test_kvlcc2_beam_cell_is_attributed_to_a_scale_slip(verdict: dict) -> None:
    """The ``Bwl 580`` typo is found in the *dissent*, which is where round 3 put it.

    One lineage printing 580 (SIMMAN 2014/2020 pages) is outvoted by two printing the
    corrected 58 (Osaka dissertation, UGent reprint), so the family value closes on its own
    ``Cb·Lpp·Bwl·T`` and a detector that only looks at the winning numeral finds nothing --
    which is what had turned this assertion into an ``xfail``.  The typo is still printed in
    two entries of the library, so the scan must audit the dissenting numerals too, and it
    must return exactly the numbers round 1 reported for it.
    """
    slips = [row for row in verdict["scale_typos"] if row["ship"] == "KVLCC2"]
    assert len(slips) == 1
    slip = slips[0]
    assert slip["quantity"] == "beam"
    assert slip["origin"] == "dissent"
    assert slip["family_value"] == pytest.approx(58.0)  # what the family settled on
    assert slip["printed"] == pytest.approx(580.0)
    assert slip["divide_by"] == pytest.approx(10.0)
    assert slip["corrected"] == pytest.approx(58.0)
    assert abs(slip["corrected_residual"]) <= CLOSURE_TOLERANCE
    assert abs(slip["raw_residual"]) > CLOSURE_TOLERANCE


def test_scale_slip_does_not_rewrite_the_printed_cell(library: dict, verdict: dict) -> None:
    """The library keeps the bytes as printed; the finding rides alongside them."""
    entries = library["report"]["KVLCC2"]["quantities"]["beam"]["entries"]
    printed = [entry["cell"]["value"] for entry in entries]
    assert max(printed) == pytest.approx(580.0)  # untouched
    flagged = [entry for entry in entries if "closure_defect" in entry["cell"]]
    assert len(flagged) >= 1
    assert {f["cell"]["closure_defect"]["divide_by"] for f in flagged} == {10.0}
    assert {f["cell"]["closure_defect"]["origin"] for f in flagged} == {"dissent"}
    # Every flagged cell is one that prints the typo, never the corrected numeral.
    assert all(f["cell"]["value"] == pytest.approx(580.0) for f in flagged)
    values = library["report"]["KVLCC2"]["quantities"]["beam"]["family_values"]
    assert values["SIMMAN-SERIES"] == pytest.approx(58.0)  # majority reading, never the mean
    dissent = library["report"]["KVLCC2"]["quantities"]["beam"]["family_disagreement"]
    assert sorted(dissent["SIMMAN-SERIES"]["dissent"]) == [580.0, 580.0]


def test_healthy_hull_is_left_alone(verdict: dict) -> None:
    """KCS closes on every lineage and has no dispute, so nothing is adjudicated."""
    assert verdict["ships"]["KCS"]["conflicts"] == {}
    assert all(row["consistent"] for row in verdict["ships"]["KCS"]["closure"].values())
    assert all(
        row["explained_by_scale_slip"] is None
        for row in verdict["ships"]["KCS"]["closure"].values()
    )
    assert not [row for row in verdict["scale_typos"] if row["ship"] == "KCS"]


def test_hulls_without_published_cells_are_not_invented(verdict: dict) -> None:
    """DTC and Series 60 have no source: no closure, no dispute, no fabricated number."""
    for ship in ("Series60",):
        assert verdict["ships"][ship]["closure"] == {}
        assert verdict["ships"][ship]["conflicts"] == {}


def test_unresolved_conflict_is_reported_not_guessed() -> None:
    """Two families that both close, or neither, must leave the conflict open."""
    # Two families stating the same coefficients but ∇ 1 % apart: both close to within
    # tolerance (10000 vs 10000 exactly would be one family, so perturb the volume).
    report2 = cross_reference(
        {
            "SYNTH": {
                "a": {
                    "displacement_volume": Cell("displacement_volume", 10000.0, "10000", "displacement_volume", "a", "u", "0" * 64, 1),
                    "cb": Cell("cb", 0.5, "0.5", "cb", "a", "u", "0" * 64, 1),
                    "length_pp": Cell("length_pp", 100.0, "100", "length_pp", "a", "u", "0" * 64, 1),
                    "beam": Cell("beam", 20.0, "20", "beam", "a", "u", "0" * 64, 1),
                    "draft": Cell("draft", 10.0, "10", "draft", "a", "u", "0" * 64, 1),
                },
                "b": {
                    "displacement_volume": Cell("displacement_volume", 10080.0, "10080", "displacement_volume", "b", "u", "0" * 64, 1),
                    "cb": Cell("cb", 0.5, "0.5", "cb", "b", "u", "0" * 64, 1),
                    "length_pp": Cell("length_pp", 100.0, "100", "length_pp", "b", "u", "0" * 64, 1),
                    "beam": Cell("beam", 20.0, "20", "beam", "b", "u", "0" * 64, 1),
                    "draft": Cell("draft", 10.0, "10", "draft", "b", "u", "0" * 64, 1),
                },
            }
        }
    )
    assert report2["SYNTH"]["quantities"]["displacement_volume"]["verdict"] == "CONFLICT"
    verdict = adjudicate(report2)
    row = verdict["ships"]["SYNTH"]["conflicts"]["displacement_volume"]
    # 10000/10000 closes exactly; 10080 misses by 0.8 % -- outside the 0.5 % agreement
    # tolerance (so the library reports CONFLICT) but inside the 1 % closure tolerance, so
    # both families close and the library refuses to pick between them.
    assert row["resolution"] == "unresolved"
    assert row["preferred_family"] is None
    assert len(row["family_values"]) == 2


def test_single_scale_slip_only_reports_one_quantity() -> None:
    """A slip in Bwl is attributed to Bwl only -- the other cells stay consistent."""
    verdict = adjudicate(
        _synthetic(
            {
                "SYNTH": {
                    "displacement_volume": (312622.0, "A"),
                    "cb": (0.8098, "A"),
                    "length_pp": (320.0, "A"),
                    "beam": (580.0, "A"),
                    "draft": (20.8, "A"),
                }
            }
        )
    )
    assert [row["quantity"] for row in verdict["scale_typos"]] == ["beam"]
    assert verdict["ships"]["SYNTH"]["closure"]["SYNTH_A"]["consistent"] is False


def test_coverage_summary_counts_stay_consistent(library: dict) -> None:
    cov = coverage_summary(library["report"])
    assert cov["cells"] == len(cov["ships"]) * len(QUANTITIES)
    assert cov["cross_verified"] + cov["single_source"] + cov["conflict"] + cov["missing"] == cov["cells"]