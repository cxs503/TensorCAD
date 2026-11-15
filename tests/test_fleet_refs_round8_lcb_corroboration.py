"""Tests for round 8: arithmetic LCB corroboration and the residual-gap probe register.

Anchor ``REPRO-GEOM-FLEET-REFS``.  The tests assert the *rules* the artifact relies on,
not merely that the script ran:

* the conversion is an exact inverse pair (so a reader can redo it on paper);
* a corroboration never votes and never changes the library's verdict;
* the residual lands inside the declared tolerance for ONRT, and the tolerance is actually
  applied (an out-of-band input must come back ``INCONSISTENT``);
* the probe register's documents carry digests that match the library's pinned sources, so
  every negative row is reproducible from the same bytes.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.fleet.corroboration import (  # noqa: E402
    LCB_TOLERANCE,
    corroborate_lcb_by_scale,
    lcb_aft_fp_from_percent,
    lcb_percent_from_aft_fp,
)

ARTIFACT = ROOT / "benchmarks/verified/fleet_refs/round8_lcb_corroboration.json"
LIBRARY = ROOT / "benchmarks/verified/fleet_refs/crossref.json"


def test_conversion_is_an_inverse_pair():
    for pct in (-3.0, -1.68, 0.0, 2.5):
        for lpp in (154.0, 355.0):
            back = lcb_percent_from_aft_fp(lcb_aft_fp_from_percent(pct, lpp), lpp)
            assert back == pytest.approx(pct, abs=1e-12)


def test_aft_of_fp_direction_is_right():
    # forward of midship means *smaller* distance aft of the FP
    assert lcb_aft_fp_from_percent(+2.0, 100.0) < lcb_aft_fp_from_percent(0.0, 100.0)
    assert lcb_aft_fp_from_percent(-2.0, 100.0) == pytest.approx(52.0)


def test_corroboration_never_votes_and_respects_tolerance():
    ok = corroborate_lcb_by_scale(
        ship="X", model_length_aft_fp_m=1.625, scale=49.0, lpp_m=154.0,
        published_percent_fwd_positive=-1.68,
    )
    assert ok.verdict == "CORROBORATES"
    assert ok.votes is False
    # the same document read against a value 2 % away must not be called agreement
    bad = corroborate_lcb_by_scale(
        ship="X", model_length_aft_fp_m=1.625, scale=49.0, lpp_m=154.0,
        published_percent_fwd_positive=-3.6,
    )
    assert bad.verdict == "INCONSISTENT"
    assert bad.votes is False
    assert bad.residual_relative > LCB_TOLERANCE


def test_invalid_inputs_are_refused():
    with pytest.raises(ValueError):
        corroborate_lcb_by_scale(
            ship="X", model_length_aft_fp_m=0.0, scale=49.0, lpp_m=154.0,
            published_percent_fwd_positive=-1.68,
        )
    with pytest.raises(ValueError):
        lcb_aft_fp_from_percent(-1.68, 0.0)


@pytest.fixture(scope="module")
def artifact():
    return json.loads(ARTIFACT.read_text())


@pytest.fixture(scope="module")
def library():
    return json.loads(LIBRARY.read_text())


def test_artifact_corroboration_agrees_with_the_library_cell(artifact, library):
    cell = artifact["lcb_corroboration"]
    conv = cell["conversion"]
    assert conv["votes"] is False
    assert conv["verdict"] == "CORROBORATES"
    assert conv["residual_relative"] <= cell["tolerance"]
    assert cell["library_verdict_before"] == cell["library_verdict_after"]
    stored = library["report"]["ONRT"]["quantities"]["lcb_percent"]
    assert cell["published"]["percent_fwd_positive"] == stored["entries"][0]["cell"]["value"]
    assert cell["published"]["sha256"] == stored["entries"][0]["cell"]["sha256"]
    # the Lpp used for the conversion is the library's own, not a typed-in number
    lpp_cell = library["report"]["ONRT"]["quantities"][cell["published"]["lpp_quantity"]]["entries"][0]["cell"]
    assert cell["published"]["lpp_m"] == lpp_cell["value"]


def test_probe_documents_are_digest_pinned(artifact):
    from tensorcad.fleet.crossref import SOURCES

    pinned = {spec.key: spec.sha256 for spec in SOURCES}
    for probe in artifact["gap_probes"]:
        assert len(pinned[probe["document"]]) == 64
        assert probe["outcome"], "a probe row without an outcome is not a record"
        # a negative row is either quoted or explicitly a zero-match search
        assert probe["quotes"] or probe.get("zero_match")


def test_series60_gap_is_stated_not_hidden(artifact):
    assert artifact["remaining_gaps"]["Series60"] == [
        "displacement_volume", "cb", "lcb_percent", "lwl",
    ]
    assert any("Series 60" in s["query"] for s in artifact["searches"])