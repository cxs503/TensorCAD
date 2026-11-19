"""Intake gates for publisher hull archives (anchor GAP-CAD-T1-IGES-PIPE).

Two defects found by running the baseline-fleet intake on the real workshop files:

* the IGES *global section* can contradict the coordinates.  ``KVLCC2/kvlcc2.igs``
  declares unit flag 2 (millimetres) over metre coordinates, so the declared factor
  turned a 320 m ship into a 0.334 m one -- still watertight, still a mesh, of the
  wrong object.  :func:`measured_unit_factor` measures both readings against the
  publisher's own ``L_PP`` and records the contradiction instead of guessing;
* a patch set in which *every* basis surface is trimmed is not described by its
  untrimmed domains: KVLCC2's untrimmed sample grid reaches 7 m below the keel and
  up to 6x outside the trimmed region, which is exactly the data the half-breadth
  field takes ``max |y|`` over.  :func:`sample_trimmed_patches` samples inside the
  trim loops instead.

The pointer decision inside a type-142 record is *measured*, and the measurement has
a trap this file pins down: a degenerate curve (one repeated point) sitting in a
corner of the knot domain scores a perfect inside-fraction, so it must not be allowed
to win.  ``test_degenerate_candidate_cannot_win`` is that case.
"""
from __future__ import annotations

import math
import os
from pathlib import Path

import pytest

from tensorcad.fleet_baseline import (
    UNIT_BAND,
    measured_unit_factor,
    point_in_polygon,
    sample_trimmed_patches,
    trim_loops_uv,
)

ROOT = Path(__file__).resolve().parents[1]
# The library reads $TENSORCAD_ASSET_ROOT (tensorcad.assetlib.store.ENV_VAR);
# TENSORCAD_ASSET_STORE is accepted for backwards compatibility with earlier gates.
STORE = Path(
    os.environ.get("TENSORCAD_ASSET_ROOT")
    or os.environ.get("TENSORCAD_ASSET_STORE")
    or "/srv/tensor/agent/assets"
)
KVLCC2 = STORE / "geometry/hulls/KVLCC2/kvlcc2.igs"
SPECS = ROOT / "benchmarks/verified/baseline_fleet/fleet_specs.json"


# --------------------------------------------------------------------------- #
# unit gate
# --------------------------------------------------------------------------- #
def test_declared_unit_kept_when_it_implies_a_plausible_hull():
    # KCS: metres, 243.9 m of point cloud for a 230 m L_PP -> plausible.
    record = measured_unit_factor(243.854, 1.0, 230.0, metre_factor=1.0)
    assert record["verdict"] == "declared_consistent"
    assert record["factor"] == pytest.approx(1.0)
    assert record["candidates"][0]["ratio_to_published_lpp"] == pytest.approx(1.0602, abs=1e-3)


def test_millimetre_flag_over_metre_coordinates_is_overridden_and_recorded():
    # KVLCC2: flag 2 (mm) over a 333.6 m point cloud for a 320 m L_PP.
    record = measured_unit_factor(333.619, 1e-3, 320.0, metre_factor=1.0)
    assert record["verdict"] == "declared_contradicted_by_geometry"
    assert record["factor"] == pytest.approx(1.0)
    assert record["relative_to_declared"] == pytest.approx(1000.0)
    declared = record["candidates"][0]
    assert declared["implied_length_m"] == pytest.approx(0.3336, abs=1e-4)
    assert declared["plausible"] is False
    assert "contradicts its own coordinates" in record["reason"]


def test_no_plausible_reading_keeps_declared_and_says_so():
    # ONRT as it arrives: neither reading gives a 154 m hull, so no correction is
    # invented -- the downstream scale guard is left to refuse.
    record = measured_unit_factor(1.006, 1.0, 154.0, metre_factor=1.0)
    assert record["verdict"] == "no_plausible_unit"
    assert record["factor"] == pytest.approx(1.0)
    assert record["candidates"][0]["plausible"] is False
    assert len(record["candidates"]) == 1  # the metre reading is the declared one


def test_band_is_documented_and_decade_only():
    low, high = UNIT_BAND
    assert (low, high) == (0.5, 2.0)
    # within the band, the declared reading always wins (no silent rescaling of drafting)
    assert measured_unit_factor(1.6 * 100.0, 1.0, 100.0)["verdict"] == "declared_consistent"
    assert measured_unit_factor(0.6 * 100.0, 1.0, 100.0)["verdict"] == "declared_consistent"


# --------------------------------------------------------------------------- #
# trim loops: point-in-polygon and the pointer decision
# --------------------------------------------------------------------------- #
def test_point_in_polygon_square_and_hole():
    square = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]
    assert point_in_polygon((1.0, 1.0), square) is True
    assert point_in_polygon((3.0, 1.0), square) is False
    # even-odd: two nested loops -> the annulus between them is inside, the hole is not
    hole = [(0.5, 0.5), (1.5, 0.5), (1.5, 1.5), (0.5, 1.5)]
    assert sum(1 for loop in (square, hole) if point_in_polygon((0.25, 0.25), loop)) % 2 == 1
    assert sum(1 for loop in (square, hole) if point_in_polygon((1.0, 1.0), loop)) % 2 == 0


class FlatPatch:
    """A patch whose parameter plane *is* its 3-D plane: 10 m x 1 m."""

    def __init__(self, u=(0.0, 10.0), v=(0.0, 1.0)):
        self._u = u
        self._v = v

    @property
    def domain_u(self):
        return self._u

    @property
    def domain_v(self):
        return self._v

    def evaluate(self, u, v):
        return (u, v, 0.0)


class Curve:
    def __init__(self, curve_uv=None, curve_3d=None):
        self.curve_uv = curve_uv
        self.curve_3d = curve_3d


class StubModel:
    """Duck-typed stand-in for ``IgesFile``: only the trim API the sampler uses."""

    def __init__(self, patches):
        self._patches = patches

    def trimmed_surfaces(self):
        return self._patches

    def trim_polygon_uv(self, pointer, samples=96):
        return list(self._polygons[pointer])

    _polygons: dict = {}


def _rectangle(x0, x1, y0, y1):
    return [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]


def _stub(polygons, patches):
    model = StubModel(patches)
    model._polygons = polygons
    return model


def test_sampling_stays_inside_the_trim_loop():
    half = _rectangle(0.0, 5.0, 0.0, 1.0)
    model = _stub({7: half}, [(object(), FlatPatch(), [Curve(curve_uv=7)])])
    points, record = sample_trimmed_patches(model, spacing=0.5)
    assert record["patches_used"] == 1
    assert points, "a trimmed half-patch must still yield samples"
    assert max(point[0] for point in points) <= 5.0 + 1e-9
    # half of the 10 x 1 domain, at spacing 0.5 -> roughly half of 21 x 3 grid nodes
    assert record["patches"][0]["kept"] == pytest.approx(11 * 3, rel=0.35)


def test_degenerate_candidate_cannot_win():
    """The model-space curve is the real loop; the uv pointer samples to one point.

    A degenerate polygon sits in a corner of the knot domain, so an inside-fraction
    test alone scores it 1.0 -- it used to win and the patch came out with zero
    samples.  The spread check must reject it.
    """
    domain_spanning = _rectangle(0.0, 10.0, 0.0, 1.0)
    degenerate = [(0.0, 0.0)] * 8
    model = _stub(
        {9: degenerate, 11: domain_spanning},
        [(object(), FlatPatch(), [Curve(curve_uv=9, curve_3d=11)])],
    )
    decision = trim_loops_uv(model, FlatPatch(), [Curve(curve_uv=9, curve_3d=11)])
    assert decision["curves"][0]["verdict"] == "decided"
    assert decision["curves"][0]["field"] == "curve_3d"
    assert len(decision["loops"]) == 1
    assert decision["loops"][0] == domain_spanning


def test_undecided_trim_curves_are_reported_not_meshed():
    model = _stub(
        {9: [(0.0, 0.0)] * 4, 11: [(50.0, 50.0), (51.0, 50.0), (51.0, 51.0)]},
        [(object(), FlatPatch(), [Curve(curve_uv=9, curve_3d=11)])],
    )
    points, record = sample_trimmed_patches(model, spacing=0.5)
    assert points == []
    assert record["patches"][0]["verdict"] == "undecided_trim_curves"
    assert record["patches_used"] == 0


def test_inner_loop_is_a_hole():
    outer = _rectangle(0.0, 10.0, 0.0, 1.0)
    hole = _rectangle(4.0, 6.0, 0.25, 0.75)
    model = _stub(
        {1: outer, 2: hole},
        [(object(), FlatPatch(), [Curve(curve_uv=1), Curve(curve_uv=2)])],
    )
    points, record = sample_trimmed_patches(model, spacing=0.25)
    assert record["patches"][0]["loops"] == 2
    assert points
    assert not [p for p in points if 4.0 < p[0] < 6.0 and 0.25 < p[1] < 0.75]


# --------------------------------------------------------------------------- #
# spec loading: one absent hull must not take the whole fleet down
# --------------------------------------------------------------------------- #
def test_specs_from_json_passes_through_a_pointer_that_resolves_nowhere(tmp_path):
    """One absent hull must not take the whole fleet spec down.

    ``specs_from_json`` resolves every geometry pointer through the asset library.  A
    pointer that resolves nowhere on this machine (an ``asset://`` URI whose bytes are
    not in this store, or a dead work-tree path) is passed through unchanged so the
    geometry reader reports the absence with its own message -- instead of the whole
    spec load dying on the first missing hull and taking every present hull down with
    it.  That is what happened when the JBC bytes left the store: the KVLCC2 intake
    could no longer even load its own spec.
    """
    import json as _json

    from tensorcad.fleet_baseline import specs_from_json

    # A pointer that resolves nowhere: an asset:// URI with no bytes in this store.
    absent_uri = "asset://geometry/hulls/NO_SUCH_HULL/absent.igs"
    # A pointer that resolves: a real file on disk.
    present = tmp_path / "present.igs"
    present.write_bytes(b"dummy-iges-bytes")
    spec_file = tmp_path / "fleet_specs.json"
    spec_file.write_text(
        _json.dumps(
            {
                "ships": [
                    {
                        "name": "PRESENT",
                        "draft": 10.0,
                        "length_pp": 100.0,
                        "beam": 20.0,
                        "iges_path": str(present),
                        "iges_sha256": "",
                    },
                    {
                        "name": "ABSENT",
                        "draft": 10.0,
                        "length_pp": 100.0,
                        "beam": 20.0,
                        "iges_path": absent_uri,
                        "iges_sha256": "",
                    },
                ]
            }
        )
    )
    specs = {ship.name: ship for ship in specs_from_json(spec_file)}
    assert set(specs) == {"PRESENT", "ABSENT"}
    # the present pointer resolves to the real path
    assert specs["PRESENT"].iges_path == str(present)
    # the absent pointer is passed through unchanged, not raised away
    assert specs["ABSENT"].iges_path == absent_uri


# --------------------------------------------------------------------------- #
# real-archive regression (skipped when the publisher bytes are absent)
# --------------------------------------------------------------------------- #
@pytest.mark.skipif(not KVLCC2.exists(), reason="KVLCC2 publisher geometry not in the asset store")
def test_kvlcc2_fully_trimmed_set_imports_watertight_within_three_percent():
    from tensorcad.fleet_baseline import build_ship, specs_from_json

    spec = {ship.name: ship for ship in specs_from_json(SPECS)}["KVLCC2"]
    mesh, report = build_ship(spec, dx=1.0, dz=0.1)
    assert report["source"]["unit_gate"]["verdict"] == "declared_contradicted_by_geometry"
    assert report["source"]["route"] == "trimmed_patch_loops"
    assert report["source"]["trim_sampling"]["patches_used"] == 19
    assert report["checks"]["watertight"] is True
    assert report["checks"]["edge_manifold"] is True
    error = report["derived"]["published_volume_rel_error_percent"]
    assert abs(error) < 3.0, f"volume off by {error:+.2f}% against the published displacement"
    assert report["verdict"] == "VERIFIED"
    assert math.isfinite(report["source"]["trim_sampling"]["dz_refinement"]["relative_change_percent"])