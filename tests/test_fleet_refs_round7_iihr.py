"""Round 7 artifact checks (offline: the pinned PDF is not committed).

The tests read the committed artifact only, so they hold on a clean checkout: the digit-level
read is reproducible from the pinned bytes, and every derived statement below is recomputed
here from the numbers the tool stored rather than trusted.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.baseline_fleet.refs_round7_iihr import (  # noqa: E402
    IIHR,
    LCB_QUOTE,
    MODEL_SCALE,
    resolve_basis,
)

ART = ROOT / "benchmarks/verified/fleet_refs/round7_onrt_iihr.json"


def art() -> dict:
    return json.loads(ART.read_text())


def test_pin_is_a_sha256_and_matches_observed() -> None:
    d = art()["document"]
    assert len(d["sha256_pinned"]) == 64 and len(d["sha256_observed"]) == 64
    assert d["digest_matches_pin"] is True
    assert d["family"] == "IIHR-2016" and d["lineage"] == "IIHR-2016"


def test_every_read_cell_carries_a_verbatim_quote_and_the_pinned_digest() -> None:
    d = art()["document"]
    cells = art()["cells"]
    assert set(cells) == {"lwl", "beam", "draft", "displacement_volume", "cb"}
    for quantity, cell in cells.items():
        assert cell["raw"], quantity
        assert cell["label"], quantity
        assert cell["sha256"] == d["sha256_observed"], quantity
        assert cell["url"] == d["url"], quantity


def test_main_particulars_are_the_documented_ones() -> None:
    cells = art()["cells"]
    assert cells["lwl"]["value"] == 154.0
    assert cells["beam"]["value"] == 18.78
    assert cells["draft"]["value"] == 5.494
    assert cells["cb"]["value"] == 0.535
    # the '8507 ton' column is stored as a mass, so Cell divided by the seawater density
    assert cells["displacement_volume"]["unit_basis"] == "mass"
    assert abs(cells["displacement_volume"]["value"] - 8507 / 1.025) < 1e-3  # Cell rounds to 4 dp


def test_closure_prefers_the_volume_reading_and_says_why() -> None:
    a = art()
    closure = a["closure"]
    recomputed = resolve_basis(
        a["disputed_cell"]["printed_numeral"],
        a["cells"]["lwl"]["value"],
        a["cells"]["beam"]["value"],
        a["cells"]["draft"]["value"],
        a["cells"]["cb"]["value"],
    )
    assert recomputed == closure
    assert closure["reading"] == "volume"
    assert closure["residual_if_volume_percent"] < 0.5
    assert closure["residual_if_mass_percent"] > 2.0
    # the mass reading converts to the NMRI-T2015 value to within rounding
    assert abs(a["disputed_cell"]["mass_equivalent_of_this_document_m3"] - 8299.5) < 0.2


def test_no_published_value_is_changed_and_the_cell_stays_disputed() -> None:
    a = art()
    stored = {row["family"]: row["value"] for row in a["disputed_cell"]["library_values"]}
    assert stored["NMRI-T2015"] == 8299.5
    assert stored["SIMMAN-SERIES"] == 8507.0
    assert a["disputed_cell"]["verdict"] == "UNIT_BASIS_CONFLICT"
    assert "never" not in a["disputed_cell"]["explanation"].lower()  # the rule is stated in doc, not asserted as an average


def test_lcb_does_not_vote() -> None:
    a = art()["lcb_evidence"]
    assert a["votes"] is False
    assert a["quoted"] == "1.625"
    assert a["full_scale_m_aft_fp"] == round(1.625 * MODEL_SCALE, 3)
    assert a["full_scale_m_aft_fp"] == 79.625
    assert "NA" in a["raw_row"]


def test_lcb_quote_pattern_is_specific_enough() -> None:
    snippet = "Table 2.1: Ship model and full scale dimensions \nLCB \nLCB (m) aft of FP \n1.625 \nNA \n"
    import re

    match = re.search(LCB_QUOTE, snippet)
    assert match is not None and match.group(1) == "1.625"
    assert re.search(LCB_QUOTE, "no LCB here") is None


def test_source_spec_is_registered_shape() -> None:
    assert IIHR.kind == "pdf"
    assert IIHR.ship == "ONRT"
    assert all(len(quote) == 3 for quote in IIHR.quotes)