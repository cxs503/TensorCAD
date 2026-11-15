"""File-frame policy for imported hull geometry (``tensorcad.fleet.frames``).

The published KVLCC2 IGES is a waterline-origin, z-down, stern-negative export; the
published KCS IGES is keel-origin and z-up.  A route that assumes one convention keeps
*nothing* from the other file, which is how KVLCC2 ended up with ``—`` in the benchmark
while its reference (312 622 m³) was already in the library.  These tests pin the
measurement and the transformation, so the policy cannot silently drift into an
assumption.
"""

from __future__ import annotations

import json

import pytest

from tensorcad.fleet import build as fleet_build
from tensorcad.fleet.frames import (
    NORMALISED_FRAME,
    FramePolicy,
    detect,
    normalise,
)

# The KVLCC2 asset measured at coarse patch sampling: 19 patches, z from -28.067 to
# +0.055, and only ~37 % of samples below zero.  A sample-fraction rule calls this "up";
# the depth-ratio rule (below 28.07 vs above 0.06) calls it "down", which is right.
KVLCC2_MEASURED_POINTS = [
    (x, y, z)
    for x in (0.0, -160.0, -328.0)
    for y in (0.0, 29.0)
    for z in (-28.067, -20.0, -10.0, -2.0, 0.0) + (0.055,) * 8
]
from tensorcad.fleet.manifest import ShipRecord, Provenance, load_manifest, record_from_dict


# A keel-origin, z-up half model: z in [0, 20.8], x from 0 to 320, y one-sided.
Z_UP_POINTS = [
    (x, y, z)
    for x in (0.0, 80.0, 160.0, 240.0, 320.0)
    for y in (0.0, 5.0, 15.0, 29.0)
    for z in (0.0, 5.0, 10.0, 20.8)
]

# The same hull as published for KVLCC2: z = 0 at the waterline, z down to -28.067,
# x = 0 near the bow with the stern at x = -328.1.
Z_DOWN_POINTS = [(x - 328.116, y, z - 28.067) for x, y, z in Z_UP_POINTS]


def test_detects_keel_origin_z_up_file():
    evidence = detect(Z_UP_POINTS)
    assert evidence.z_axis == "up"
    assert evidence.x_min == pytest.approx(0.0)
    assert evidence.half_model is True
    assert evidence.full_beam == pytest.approx(58.0)
    assert evidence.x_extent == pytest.approx(320.0)


def test_detects_waterline_origin_z_down_file():
    evidence = detect(Z_DOWN_POINTS)
    assert evidence.z_axis == "down"
    assert evidence.z_below_zero_fraction == pytest.approx(1.0)
    assert evidence.z_min == pytest.approx(-28.067)
    assert evidence.x_min == pytest.approx(-328.116)


def test_auto_policy_restores_the_fleet_frame():
    points, report = normalise(Z_DOWN_POINTS)
    assert report.z_flip is True
    assert report.x_shift == pytest.approx(328.116, abs=1e-6)
    assert report.evidence.z_axis == "down"
    assert min(p[0] for p in points) == pytest.approx(0.0, abs=1e-6)
    assert max(p[0] for p in points) == pytest.approx(320.0, abs=1e-6)
    assert min(p[2] for p in points) == pytest.approx(7.267, abs=1e-6)
    assert max(p[2] for p in points) == pytest.approx(28.067, abs=1e-6)
    assert report.as_dict()["frame"] == NORMALISED_FRAME


def test_auto_policy_is_a_no_op_on_an_already_normalised_file():
    points, report = normalise(Z_UP_POINTS)
    assert report.z_flip is False
    assert report.x_shift == pytest.approx(0.0)
    assert points == [(round(x, 9), y, z) for x, y, z in Z_UP_POINTS]


def test_depth_ratio_beats_the_sample_fraction_on_the_real_export():
    evidence = detect(KVLCC2_MEASURED_POINTS)
    assert evidence.z_below_zero_fraction < 0.5  # a fraction rule would get this wrong
    assert evidence.below_zero == pytest.approx(28.067, abs=1e-6)
    assert evidence.above_zero == pytest.approx(0.055, abs=1e-6)
    assert evidence.z_axis == "down"


def test_forced_policy_overrides_the_measurement():
    _, report = normalise(Z_DOWN_POINTS, {"z_axis": "up", "x_origin": "as_is"})
    assert report.z_flip is False and report.x_shift == 0.0
    assert report.policy == FramePolicy(z_axis="up", x_origin="as_is")


def test_invalid_policy_is_refused():
    with pytest.raises(ValueError):
        FramePolicy(z_axis="sideways")
    with pytest.raises(ValueError):
        FramePolicy(x_origin="guess")


def test_empty_cloud_is_refused():
    with pytest.raises(ValueError):
        detect([])


def test_frame_policy_round_trips_through_the_manifest(tmp_path):
    payload = {
        "name": "KVLCC2",
        "ship_class": "tanker",
        "length": 320.0,
        "beam": 58.0,
        "draft": 20.8,
        "block_coefficient": 0.8098,
        "displacement_basis": "literature_value_verified",
        "displacement_volume": 312622.0,
        "geometry_route": "iges_field",
        "geometry_path": "asset://geometry/hulls/KVLCC2/kvlcc2.igs",
        "geometry_frame": {"z_axis": "auto", "x_origin": "auto", "x_extent": "geometry"},
        "provenance": {"source": "SIMMAN", "url": "https://example.invalid", "license": "n/a",
                       "retrieval": "fetched", "verified": True},
    }
    record = record_from_dict(payload)
    assert record.geometry_frame == payload["geometry_frame"]
    assert record.as_dict()["geometry_frame"] == payload["geometry_frame"]


def test_build_dispatch_hands_the_frame_policy_to_the_iges_route(monkeypatch):
    seen = {}

    def fake_iges_field_mesh(record, path, **kwargs):
        seen["path"] = str(path)
        seen["frame"] = kwargs.get("frame")
        seen["x_extent"] = kwargs.get("x_extent")
        raise RuntimeError("stop after capture")

    monkeypatch.setattr(fleet_build, "iges_field_mesh", fake_iges_field_mesh)
    record = ShipRecord(
        name="KVLCC2",
        ship_class="tanker",
        length=320.0,
        beam=58.0,
        draft=20.8,
        block_coefficient=0.8098,
        displacement_basis="literature_value_verified",
        provenance=Provenance("SIMMAN", "https://example.invalid", "n/a", "fetched", True),
        geometry_route="iges_field",
        geometry_path="asset://geometry/hulls/KVLCC2/kvlcc2.igs",
        geometry_frame={"z_axis": "auto", "x_origin": "auto", "x_extent": "geometry"},
    )
    with pytest.raises(RuntimeError):
        fleet_build.build_ship_geometry(record)
    assert seen["frame"] == record.geometry_frame
    assert seen["x_extent"] == "geometry"
    assert seen["path"].startswith("asset://")


def test_packaged_manifest_still_has_no_iges_route():
    """The offline fleet must keep working without the asset store (documented policy)."""
    packaged = load_manifest()
    assert all(ship.geometry_route != "iges_field" for ship in packaged.ships)
    assert json.dumps([ship.geometry_frame for ship in packaged.ships]) == json.dumps(
        [None] * len(packaged.ships)
    )

# ---------------------------------------------------- x-origin: AP vs forward end
#
# The `auto` x policy has to answer *where the file's own origin is*, not merely *where the
# model starts*.  The published KCS export is the case that exposed the difference: its
# origin is the aft perpendicular, with 6 m of transom overhang behind it and a bulb 7.8 m
# ahead of the FP.  Translating its aft-most point to 0 slides the hull 6 m forward into
# the Lpp window, and the window then cuts the forebody off -- measured length 227 m
# against a declared 230 m, which the principal-dimension gate refuses.

# KCS as published: x = 0 at the AP, transom overhang behind it, bulb ahead of the FP.
KCS_PUBLISHED_POINTS = [
    (x, y, z)
    for x in (-6.0, 0.0, 60.0, 120.0, 180.0, 230.0, 237.84)
    for y in (0.0, 8.0, 16.12)
    for z in (0.0, 5.4, 10.8)
]


def test_aft_perpendicular_origin_is_not_translated():
    points, report = normalise(KCS_PUBLISHED_POINTS)
    assert report.z_flip is False
    assert report.x_shift == 0.0
    assert min(p[0] for p in points) == pytest.approx(-6.0)
    assert "aft perpendicular" in report.x_origin_reason
    assert report.as_dict()["x_origin_reason"] == report.x_origin_reason


def test_model_forward_of_its_own_origin_is_translated_to_zero():
    shifted = [(x + 3.0, y, z) for x, y, z in Z_UP_POINTS]
    _, report = normalise(shifted)
    assert report.x_shift == pytest.approx(-3.0)
    assert "forward of its own x = 0" in report.x_origin_reason


def test_stern_negative_export_is_still_translated():
    _, report = normalise(Z_DOWN_POINTS)
    assert report.x_shift == pytest.approx(328.116, abs=1e-6)
    assert "stern-negative" in report.x_origin_reason


def test_the_kcs_and_kvlcc2_policies_disagree_for_the_right_reason():
    """Same publisher, opposite x decisions -- and each decision is measured, not assumed."""
    _, kcs = normalise(KCS_PUBLISHED_POINTS)
    _, kvlcc2 = normalise(Z_DOWN_POINTS)
    assert kcs.x_shift == 0.0 and kvlcc2.x_shift != 0.0
    assert abs(kcs.evidence.x_min) / kcs.evidence.x_extent < 0.5
    assert abs(kvlcc2.evidence.x_min) / kvlcc2.evidence.x_extent > 0.5
