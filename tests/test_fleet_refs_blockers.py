"""Tests for the per-cell blocker dossier (anchor ``REPRO-GEOM-FLEET-REFS``).

The dossier is an *accounting* artefact: it may not introduce a numeral and may not change a
verdict.  These tests pin exactly that: the classification is a function of the stored
library, it reconciles with ``acceptance.json``, and a tampered acceptance is refused.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools" / "baseline_fleet"))

from refs_blockers import QUANTITIES, _blockers, _markdown  # noqa: E402

LIBRARY = json.loads((ROOT / "benchmarks/verified/fleet_refs/crossref.json").read_text())
ACCEPTANCE = json.loads((ROOT / "benchmarks/verified/fleet_refs/acceptance.json").read_text())


@pytest.fixture(scope="module")
def payload() -> dict:
    return _blockers(LIBRARY, ACCEPTANCE)


def test_blocked_cells_reconcile_with_acceptance(payload: dict) -> None:
    coverage = ACCEPTANCE["coverage"]
    assert payload["summary"]["cells_checked"] == coverage["cells"]
    assert payload["summary"]["cells_cross_verified"] == coverage["cross_verified"]
    assert payload["summary"]["cells_blocked"] == coverage["cells"] - coverage["cross_verified"]
    assert payload["summary"]["ships_failing_min_two_sources"] == sorted(
        coverage["ships_failing_min_two_sources"]
    )


def test_every_blocked_cell_is_described_once(payload: dict) -> None:
    seen = [(cell["ship"], cell["quantity"]) for cell in payload["blockers"]]
    assert len(seen) == len(set(seen))
    for cell in payload["blockers"]:
        assert cell["quantity"] in QUANTITIES
        assert cell["verdict"] != "CROSS_VERIFIED"


def test_no_value_is_asserted_by_the_dossier(payload: dict) -> None:
    """The dossier points at evidence; it must never carry a numeral of its own."""
    for cell in payload["blockers"]:
        assert "value" not in cell
        assert "family_values" not in cell


def test_blocker_kinds_are_evidence_derived(payload: dict) -> None:
    by_key = {(cell["ship"], cell["quantity"]): cell for cell in payload["blockers"]}
    # DTC prints LCB only as a converted LCG (JSR-2017 Table 1): converted, never a vote.
    dtc_lcb = by_key[("DTC", "lcb_percent")]
    assert dtc_lcb["blocker"] == "CONVERTED_ONLY"
    assert dtc_lcb["voting_sources"] == []
    assert dtc_lcb["converted_corroboration"] is True
    assert dtc_lcb["lineages_missing"] == 2
    # DTC prints no Lwl anywhere in the pinned set.
    assert by_key[("DTC", "lwl")]["blocker"] == "NO_DOCUMENT"
    # Series 60 has no document at all.
    assert all(
        by_key[("Series60", quantity)]["blocker"] == "NO_DOCUMENT" for quantity in QUANTITIES
    )
    # KVLCC2 / ONRT LCB: one printed lineage plus one converted reading.
    for ship, derived in (("KVLCC2", "kvlcc2_ugent"), ("ONRT", "t2015_onrt")):
        cell = by_key[(ship, "lcb_percent")]
        assert cell["blocker"] == "SINGLE_LINEAGE"
        assert cell["voting_sources"]
        assert cell["derived_sources"] == [derived]
        assert cell["converted_corroboration"] is True
        assert cell["lineages_missing"] == 1
    # The single disputed cell keeps both values instead of being averaged away.
    onrt_disp = by_key[("ONRT", "displacement_volume")]
    assert onrt_disp["blocker"] == "DISPUTED"
    assert onrt_disp["lineages_missing"] == 0


def test_every_blocker_names_a_target_and_an_action(payload: dict) -> None:
    for cell in payload["blockers"]:
        assert cell["target_documents"], cell
        assert cell["next_action"]
        assert cell["why_pinned_set_is_not_enough"]


def test_drift_from_acceptance_is_refused() -> None:
    tampered = json.loads(json.dumps(ACCEPTANCE))
    tampered["coverage"]["cross_verified"] += 1
    with pytest.raises(SystemExit):
        _blockers(LIBRARY, tampered)


def test_markdown_reports_scope_and_no_accuracy_claim(payload: dict) -> None:
    text = _markdown(payload)
    assert "Nothing here upgrades a verdict" in text
    assert "no accuracy claim is made this round" in text
    assert "| ship | passes >=2 independent sources | quantities short |" in text