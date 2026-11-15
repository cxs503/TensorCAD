"""Tests for the longitudinal reference-frame converter (``tensorcad.fleet.lcb_frames``).

What is asserted here is the rule set, not the arithmetic alone: a frame conversion of a
full-scale numeral may vote, a model-to-full-scale multiplication may not, an LCG-for-LCB
substitution needs even keel, and an unstated datum is refused instead of guessed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tensorcad.fleet import lcb_frames as LF

ROOT = Path(__file__).resolve().parents[1]
REGISTRY = ROOT / "benchmarks/verified/fleet_refs/lcb_frames.json"


def frame(**kw) -> LF.Frame:
    base = {"datum": "ap", "sign": "fwd", "reference": "Lpp", "scale": 1.0}
    base.update(kw)
    return LF.Frame(**base)


# --- geometry of the conversion -------------------------------------------------------


def test_ap_datum_is_the_identity_on_the_length():
    # 171.1 m from the AP on a 320 m hull is 11.1 m forward of the mid-point.
    assert LF.to_percent_lpp(171.1, frame(datum="ap"), 320.0) == pytest.approx(3.46875, abs=1e-6)


def test_fp_datum_with_aft_sign_is_measured_backwards_from_the_forward_perpendicular():
    # 79.625 m aft of the FP on a 154 m hull: 74.375 m from the AP, 2.625 m aft of midship.
    assert LF.to_percent_lpp(79.625, frame(datum="fp", sign="aft"), 154.0) == pytest.approx(-1.7045, abs=1e-3)


def test_lpp_midship_datum_passes_percentages_through_unchanged():
    # a value already referred to the Lpp mid-point is only expressed as metres (1.48 % of 230 m)
    assert LF.to_percent_lpp(-0.0148 * 230.0, frame(datum="lpp_mid", sign="fwd"), 230.0) == pytest.approx(-1.48, abs=1e-9)


def test_lwl_midpoint_datum_shifts_by_half_the_bulb_overhang():
    # KCS: Lwl 232.5 vs Lpp 230 -> the waterline mid-point sits 1.25 m forward of the Lpp one.
    # The same physical point read as 4.654 m aft of the Lwl mid-point is 3.404 m aft of the
    # Lpp mid-point, i.e. the published -1.48 %Lpp.
    got = LF.to_percent_lpp(4.654, frame(datum="lwl_mid", sign="aft"), 230.0, lwl=232.5)
    assert got == pytest.approx(-1.48, abs=5e-3)
    assert LF._lwl_mid_offset(230.0, 232.5, None) == pytest.approx(1.25)


def test_model_scale_is_applied_before_the_datum():
    model = LF.to_x_from_ap(1.625, frame(datum="fp", sign="aft", scale=49.0), 154.0)
    assert model == pytest.approx(154.0 - 1.625 * 49.0)


def test_trimmed_reading_would_move_lcb_so_the_offset_is_fitted_by_the_document():
    # An explicit offset overrides the (Lwl - Lpp)/2 default: the overhang need not be forward.
    assert LF._lwl_mid_offset(230.0, 232.5, -0.4) == pytest.approx(-0.4)


# --- refusals --------------------------------------------------------------------------


def test_uninterpretable_datum_is_refused():
    with pytest.raises(LF.FrameError):
        frame(datum="somewhere near the stern")


def test_lwl_datum_without_a_length_is_refused():
    with pytest.raises(LF.FrameError):
        LF.to_percent_lpp(1.0, frame(datum="lwl_mid", sign="aft"), 230.0)


def test_lcg_may_not_stand_in_for_lcb_without_a_stated_even_keel():
    record = {
        "id": "x", "ship": "KVLCC2", "source": "kvlcc2_ugent", "value": 171.1,
        "frame": {"datum": "ap", "sign": "fwd"}, "lpp_m": 320.0,
        "provenance": LF.IDENTITY_DERIVED, "even_keel_stated": False,
    }
    with pytest.raises(LF.FrameError):
        LF.canonicalise(record)


def test_unstated_provenance_is_refused():
    record = {
        "id": "x", "ship": "KCS", "source": "s", "value": 1.0,
        "frame": {"datum": "ap", "sign": "fwd"}, "lpp_m": 230.0, "provenance": LF.UNSTATED,
    }
    with pytest.raises(LF.FrameError):
        LF.canonicalise(record)


# --- the rule that decides whether a reading may vote ----------------------------------


@pytest.mark.parametrize(
    "provenance,expected",
    [
        (LF.EXPLICIT_FRAME, True),
        (LF.SCALE_DERIVED, False),
        (LF.IDENTITY_DERIVED, False),
        (LF.COMPUTED_NOT_PUBLISHED, False),
    ],
)
def test_only_an_explicit_frame_votes(provenance, expected):
    record = {
        "id": "x", "ship": "KVLCC2", "source": "kvlcc2_ugent", "value": 171.1,
        "frame": {"datum": "ap", "sign": "fwd"}, "lpp_m": 320.0,
        "provenance": provenance, "even_keel_stated": True,
    }
    assert LF.canonicalise(record)["votes"] is expected


# --- the committed registry and its readings -------------------------------------------


def test_registry_records_are_all_convertible_and_label_their_frame():
    registry = json.loads(REGISTRY.read_text())
    assert registry["anchor"] == "REPRO-GEOM-FLEET-REFS"
    for record in registry["records"]:
        assert record["quote"], f"{record['id']} has no quote: a conversion without its source text is unverifiable"
        assert record["sha256"], f"{record['id']} has no digest"
        assert record["frame"]["datum"], f"{record['id']} states no datum"
        converted = LF.canonicalise(record)
        assert converted["provenance"] in LF.PROVENANCE_CLASSES
        assert converted["votes"] == (converted["provenance"] in LF.VOTING_CLASSES)


def test_committed_readings_match_their_registry_arithmetic():
    out = ROOT / "benchmarks/verified/fleet_refs/lcb_frames_out.json"
    if not out.exists():
        pytest.skip("run tools/baseline_fleet/lcb_frames_build.py to produce the artefact")
    payload = json.loads(out.read_text())
    registry = json.loads(REGISTRY.read_text())
    by_id = {row["id"]: row for row in payload["readings"]}
    for record in registry["records"]:
        row = by_id[record["id"]]
        assert row["value_canonical"] == pytest.approx(LF.canonicalise(record)["value_canonical"])
        assert row["votes"] is False, "no committed reading in this round may vote on a published value"


def test_kvlcc2_lcg_corroborates_the_published_lcb_within_the_lcb_tolerance():
    """The one hard number this round buys: JSR-2017's LCG, converted, lands on SIMMAN's LCB."""
    out = ROOT / "benchmarks/verified/fleet_refs/lcb_frames_out.json"
    if not out.exists():
        pytest.skip("artefact absent")
    rows = {row["id"]: row for row in json.loads(out.read_text())["readings"]}
    residual = rows["kvlcc2_jsr2017_lcg"]["residual"]
    assert abs(residual["delta_relative"]) < 0.005, residual
    # and it is corroboration only: the verdict axis must not move
    assert rows["kvlcc2_jsr2017_lcg"]["votes"] is False