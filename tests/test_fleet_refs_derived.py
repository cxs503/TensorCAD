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

#: Committed, licence-compliant reproduction of the NMRI T2015 ONRT table *structure*
#: (model/full-scale columns; LCB printed as a distance from FP with the full-scale cell
#: ``NA``; Lwl printed in both columns so the document's own ratio is pinned).  Unlike the
#: pinned page bytes -- which stay git-ignored because the publisher's bytes are not
#: redistributed -- this fixture is original content present in every worktree, so the
#: round-4 derivation is exercised in CI instead of being skipped.
FIXTURE = Path(__file__).resolve().parent / "fixtures" / "fleet_refs" / "onrt_model_scale.html"


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


def test_derived_lcb_is_reproducible_without_page_bytes() -> None:
    """The round-4 derivation must hold on the committed fixture, not only on cached bytes.

    A fresh worktree has no ``artifacts/fleet_refs/raw`` (the publisher's bytes are
    git-ignored), so the two tests above skip there.  This one runs everywhere: it feeds
    the committed structural fixture through the same reader and demands the same
    derived reading, so a regression in :func:`derive_cells` / :func:`votes` is caught
    even when the pinned bytes are absent.
    """
    spec = next(s for s in SOURCES if s.key == "t2015_onrt")
    cells = parse_source(dataclasses.replace(spec, section=None), FIXTURE.read_bytes())

    # full-scale LCB is NA -> the reader falls back to the model column and must *derive*
    # the full-scale value from the table's own Lwl ratio, never drop the numeral.
    lcb = cells["lcb_percent"]
    assert lcb.unit_basis == "derived_length_aft_fp"
    assert "48.9355" in lcb.derivation  # 154.0 / 3.147, the table's own ratio
    assert lcb.value == pytest.approx(-1.6365, abs=1e-3)
    assert not votes(lcb), "a derived cell must not vote as a published second source"

    # derivation must not leak onto the transcribed cells: they keep their basis and vote.
    lwl = cells["lwl"]
    assert lwl.unit_basis == "volume"
    assert lwl.value == pytest.approx(154.0)
    assert votes(lwl)

    # and in the cross-reference report the derived reading is kept for the audit trail
    # but the quantity still has no published (voting) source.
    report = cross_reference({"ONRT": {"t2015_onrt": cells}})
    block = report["ONRT"]["quantities"]["lcb_percent"]
    assert block["verdict"] == "MISSING"
    assert block["independent_families"] == 0
    derived = block["derived"]["cells"]
    assert len(derived) == 1
    assert derived[0]["value"] == pytest.approx(-1.6365, abs=1e-3)
    assert derived[0]["unit_basis"] == "derived_length_aft_fp"
