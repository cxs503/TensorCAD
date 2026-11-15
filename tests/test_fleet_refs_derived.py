"""A published numeral must never be lost, and arithmetic must never vote."""

from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tensorcad.fleet.crossref import (  # noqa: E402
    SOURCES,
    cross_reference,
    parse_source,
    votes,
)

RAW = Path("artifacts/fleet_refs/raw/t2015_onrt.html")


@pytest.mark.skipif(not RAW.exists(), reason="pinned page bytes are not cached in this worktree")
def test_model_scale_lcb_is_derived_not_dropped() -> None:
    spec = next(s for s in SOURCES if s.key == "t2015_onrt")
    cells = parse_source(dataclasses.replace(spec, section=None), RAW.read_bytes())
    lcb = cells["lcb_percent"]
    # the page prints "LCB (m) aft of FP = 1.625" at model scale and "NA" full scale
    assert lcb.unit_basis == "derived_length_aft_fp"
    assert "48.9355" in lcb.derivation  # ratio taken from the page's own Lwl 154.0 / 3.147
    assert lcb.value == pytest.approx(-1.6365, abs=1e-3)
    assert not votes(lcb), "a derived cell must not vote as a published second source"


def test_transcribed_cell_votes() -> None:
    spec = next(s for s in SOURCES if s.key == "t2015_kcs")
    payload = Path("artifacts/fleet_refs/raw/t2015_kcs.html")
    if not payload.exists():
        pytest.skip("pinned page bytes are not cached in this worktree")
    cells = parse_source(spec, payload.read_bytes())
    assert votes(cells["displacement_volume"])
    report = cross_reference({"KCS": {"t2015_kcs": cells}})
    block = report["KCS"]["quantities"]["displacement_volume"]
    assert block["verdict"] == "SINGLE_SOURCE"
    assert "derived" not in block
