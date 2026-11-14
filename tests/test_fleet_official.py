"""Official-source layer of the baseline fleet: page extraction, repair, record provenance.

The fixtures here are *markup shaped like the publishers' tables* with the published digits
in them; they are written for the test rather than copied out of the publisher pages, so the
suite exercises the parser's structure without redistributing the sources. A second group of
tests runs against the real cached page bytes when they are present and skips (never fakes)
when they are not -- the repository does not ship the publisher pages.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
if str(REPO / "tools" / "baseline_fleet") not in sys.path:
    sys.path.insert(0, str(REPO / "tools" / "baseline_fleet"))

import official_fleet  # noqa: E402
from tensorcad.fleet.build import RouteRefusal, build_ship_geometry  # noqa: E402
from tensorcad.fleet.manifest import (  # noqa: E402
    FleetManifest,
    Provenance,
    ShipRecord,
    load_manifest,
)
from tensorcad.fleet.official import (  # noqa: E402
    ARCHIVES,
    PAGES,
    compare_with_record,
    consistency_check,
    parse_particulars,
    published_reference,
)

CACHE = REPO / ".cache" / "baseline_fleet"


def page(rows: str) -> str:
    return f"<html><body><table>{rows}</table></body></html>"


def row(*cells: str) -> str:
    return "<tr>" + "".join(f"<td>{cell}</td>" for cell in cells) + "</tr>"


#: T2015 layout: label | symbol (unit) | full scale | model 1 | model 2
KCS_LIKE = page(
    "".join(
        [
            row("Length between perpendiculars", "LPP (m)", "230.0", "7.2786", "6.0702"),
            row("Length of waterline", "LWL (m)", "232.5", "7.3586", "6.1357"),
            row("Maximum beam of waterline", "BWL (m)", "32.2", "1.0190", "0.8498"),
            row("Draft", "T (m)", "10.8", "0.3418", "0.2850"),
            row("Displacement volume", "&nabla; (m3 )", "52030", "1.6490", "0.9571"),
            row("Block coefficient (CB)", "&nabla; /(LPPBWL T)", "0.6505", "0.6505", "0.651"),
            row("Midship section coefficient (CM)", "&nbsp;", "0.9849", "0.9849", "0.9849"),
            row("LCB (%LPP), fwd+", "&nbsp;", "-1.48", "-1.48", "-2.96"),
        ]
    )
)

#: ONRT layout: label | symbol (unit) | model scale | full scale
ONRT_LIKE = page(
    "".join(
        [
            row("Length of waterline", "LWL (m)", "3.147", "154.0"),
            row("Maximum beam of waterline", "BWL (m)", "0.384", "18.78"),
            row("Draft", "T (m)", "0.112", "5.494"),
            row("Displacement", "&Delta;", "72.6 kg", "8507 ton"),
            row("Block coefficient (CB)", "&nabla; /(LWLBWL T)", "0.535", "0.535"),
            row("LCB", "LCB (m) aft of FP", "1.625", "NA"),
        ]
    )
)

#: SIMMAN summary layout: label (unit) | full scale | model ...
KVLCC2_LIKE = page(
    "".join(
        [
            row("Lpp (m)", "320", "7", "4.2667"),
            row("Lwl (m)", "325.5", "7.1204", "4.34"),
            row("Bwl (m)", "580", "1.1688", "0.7733"),
            row("D (m)", "30", "0.6563", "0.4"),
            row("T (m)", "20.8", "0.455", "0.2773"),
            row("Displacement (m3)", "312622", "3.2724", "0.741"),
            row("CB", "0.8098", "0.8098", "0.8098"),
            row("CM", "0.998", "0.998", "0.998"),
            row("LCB (%), fwd+", "3.48", "3.48", "3.48"),
        ]
    )
)


def test_parse_reads_the_full_scale_column_of_a_symbol_layout():
    extraction = parse_particulars("KCS", KCS_LIKE)
    assert extraction.get("lpp") == pytest.approx(230.0)
    assert extraction.get("bwl") == pytest.approx(32.2)
    assert extraction.get("draft") == pytest.approx(10.8)
    assert extraction.get("displacement_volume") == pytest.approx(52030.0)
    assert extraction.get("cb") == pytest.approx(0.6505)
    assert extraction.get("cm") == pytest.approx(0.9849)
    assert extraction.get("lcb_percent") == pytest.approx(-1.48)
    # The value never leaves its row: the label and the raw cell are recorded with it.
    assert extraction.values["draft"].raw == "10.8"
    assert extraction.values["draft"].scale == "full"
    assert extraction.values["draft"].row[0].startswith("Draft")


def test_parse_reads_the_last_column_when_the_model_column_comes_first():
    extraction = parse_particulars("ONRT", ONRT_LIKE)
    assert extraction.get("lwl") == pytest.approx(154.0)
    assert extraction.get("bwl") == pytest.approx(18.78)
    assert extraction.get("draft") == pytest.approx(5.494)
    assert extraction.get("displacement_mass") == pytest.approx(8507.0)
    # LCB is printed "NA" at full scale: the parser falls back to the model column and says so.
    assert extraction.get("lcb_metres_aft_fp") == pytest.approx(1.625)
    assert extraction.values["lcb_metres_aft_fp"].scale == "reduced_scale"


def test_parse_reads_a_two_column_table():
    extraction = parse_particulars("KVLCC2", KVLCC2_LIKE)
    assert extraction.get("lpp") == pytest.approx(320.0)
    assert extraction.get("bwl") == pytest.approx(580.0)  # the page digit, unrepaired
    assert extraction.get("displacement_volume") == pytest.approx(312622.0)


def test_missing_quantity_is_absent_rather_than_defaulted():
    extraction = parse_particulars("KCS", page(row("Draft", "T (m)", "10.8")))
    assert extraction.get("lpp") is None
    assert extraction.get("displacement_volume") is None


def test_consistency_check_accepts_a_closing_page():
    check = consistency_check(parse_particulars("KCS", KCS_LIKE))
    assert check["closes"] is True
    assert check["relative_difference"] < 1e-5
    assert check["rule"] == "Cb * Lpp * Bwl * T"


def test_consistency_check_exposes_a_decimal_defect_in_the_page():
    check = consistency_check(parse_particulars("KVLCC2", KVLCC2_LIKE))
    assert check["closes"] is False
    assert check["relative_difference"] > 5.0  # 580 m of beam is 900 % off


def test_resolve_digits_repairs_the_beam_and_records_why():
    extraction = parse_particulars("KVLCC2", KVLCC2_LIKE)
    values, corrections = official_fleet.resolve_digits("KVLCC2", extraction)
    assert values["beam"] == pytest.approx(58.0)
    assert len(corrections) == 1
    assert "page prints" in corrections[0] and "58" in corrections[0]


def test_resolve_digits_leaves_a_closing_page_untouched():
    extraction = parse_particulars("KCS", KCS_LIKE)
    values, corrections = official_fleet.resolve_digits("KCS", extraction)
    assert corrections == []
    assert values["beam"] == pytest.approx(32.2)
    assert values["lcb_fraction"] == pytest.approx(-0.0148)


def test_onrt_returns_both_published_reference_currencies():
    extraction = parse_particulars("ONRT", ONRT_LIKE)
    references = published_reference("ONRT", extraction)["candidates"]
    bases = {item["basis"]: item["volume_m3"] for item in references}
    assert bases["published_mass_over_density"] == pytest.approx(8507.0 / 1.025)
    assert bases["closure"] == pytest.approx(0.535 * 154.0 * 18.78 * 5.494)
    # The two published currencies disagree by >2 %: reported, not averaged away.
    assert abs(bases["closure"] - bases["published_mass_over_density"]) / bases["closure"] > 0.02


def test_compare_with_record_flags_a_recalled_digit_that_the_page_contradicts():
    record = ShipRecord(
        name="KCS",
        ship_class="container ship",
        length=230.0,
        beam=32.2,
        draft=10.8,
        block_coefficient=0.651,
        displacement_basis="literature_value_recalled",
        provenance=Provenance("offline", "https://example.invalid", "n/a", "offline_recollection", False),
        displacement_volume=52030.0,
        midship_coefficient=0.971,
    )
    fidelity = compare_with_record(record, parse_particulars("KCS", KCS_LIKE))
    disagreeing = {item["field"] for item in fidelity["disagreements"]}
    assert disagreeing == {"midship_coefficient"}
    assert fidelity["all_agree"] is False


def test_official_record_carries_fetched_provenance_and_external_reference():
    extraction = parse_particulars("KCS", KCS_LIKE, page_sha256="a" * 64)
    record, context = official_fleet.official_record(
        "KCS", extraction, geometry_member="iges/KCS_hull_SVA.igs", member_sha256="b" * 64
    )
    assert record.geometry_route == "iges_field"
    assert record.provenance.retrieval == "fetched"
    assert record.provenance.verified is True
    assert record.reference_is_external() is True
    assert record.displacement_basis == "literature_value_verified"
    assert record.displacement_volume == pytest.approx(52030.0)
    assert record.geometry_sha256 == "b" * 64
    assert context["corrections"] == []


def test_official_manifest_round_trips_through_the_loader(tmp_path):
    record, _ = official_fleet.official_record(
        "KCS",
        parse_particulars("KCS", KCS_LIKE, page_sha256="a" * 64),
        geometry_member="iges/KCS_hull_SVA.igs",
        member_sha256="b" * 64,
    )
    manifest = official_fleet.official_manifest([record], placeholder={"max_placeholder_triangles": 304})
    target = tmp_path / "fleet_manifest_official.json"
    target.write_text(json.dumps(manifest))
    loaded = load_manifest(target)
    assert [ship.name for ship in loaded.ships] == ["KCS"]
    assert loaded.ships[0].provenance.verified is True
    # The shipped manifest is untouched: the official fleet is an additional artefact.
    packaged = load_manifest()
    assert packaged.sha256 != loaded.sha256
    assert all(ship.geometry_route != "iges_field" for ship in packaged.ships)


def test_iges_field_route_refuses_without_a_file():
    record = ShipRecord(
        name="KCS",
        ship_class="container ship",
        length=230.0,
        beam=32.2,
        draft=10.8,
        block_coefficient=0.651,
        displacement_basis="literature_value_verified",
        provenance=Provenance("NMRI T2015", "https://example.invalid", "n/a", "fetched", True),
        displacement_volume=52030.0,
        geometry_route="iges_field",
        geometry_path=None,
    )
    with pytest.raises(RouteRefusal):
        build_ship_geometry(record)


def test_page_registry_digests_are_pinned_sha256():
    for ship, entry in PAGES.items():
        assert len(entry["sha256"]) == 64, ship
        assert entry["url"].startswith("https://"), ship
        assert entry["locator"], ship
    for ship, entry in ARCHIVES.items():
        assert len(entry["sha256"]) == 16, ship
        assert entry["member"].endswith(".igs"), ship


def test_withheld_hulls_state_a_reason():
    assert set(official_fleet.WITHHELD) == {"DTC", "Series60"}
    for reason in official_fleet.WITHHELD.values():
        assert len(reason) > 40


#: The real page bytes are not in the repository (the publishers state no licence); when a
#: fetch has materialised them, the parser is checked against the actual publisher markup.
REAL_PAGES = CACHE / "raw" / "pages"


@pytest.mark.parametrize(
    "ship,file_name,expected",
    [
        ("KCS", "t2015_kcs.html", {"lpp": 230.0, "displacement_volume": 52030.0, "cb": 0.6505}),
        ("JBC", "t2015_jbc.html", {"lpp": 280.0, "displacement_volume": 178369.9, "cb": 0.8580}),
        ("ONRT", "t2015_onrt.html", {"lwl": 154.0, "draft": 5.494, "displacement_mass": 8507.0}),
        ("KVLCC2", "simman_summary.html", {"lpp": 320.0, "displacement_volume": 312622.0, "cb": 0.8098}),
    ],
)
def test_real_publisher_page_matches_the_pinned_digest_and_digits(ship, file_name, expected):
    path = REAL_PAGES / file_name
    if not path.is_file():
        pytest.skip("publisher page not materialised (run tools/baseline_fleet/fetch_sources.py)")
    raw = path.read_bytes()
    extraction = parse_particulars(ship, raw.decode("utf-8", "replace"), page_sha256=hashlib.sha256(raw).hexdigest())
    assert extraction.page_sha256 == PAGES[ship]["sha256"], "source drifted from the pinned digest"
    for quantity, value in expected.items():
        assert extraction.get(quantity) == pytest.approx(value)


def test_real_page_onrt_publishes_two_inconsistent_currencies():
    path = REAL_PAGES / "t2015_onrt.html"
    if not path.is_file():
        pytest.skip("publisher page not materialised")
    extraction = parse_particulars("ONRT", path.read_text(encoding="utf-8", errors="replace"))
    check = consistency_check(extraction)
    assert check["closes"] is False
    assert 0.02 < check["relative_difference"] < 0.03