"""DTC hull-file probe and the degenerate-sample guard (anchor REPRO-GEOM-FLEET-REFS).

Two things are asserted here:

* the geometry route refuses a sample cloud that has collapsed to a point, instead of
  lofting a zero-volume "hull" that no volume gate can distinguish from a bad integral;
* the committed artefact ``benchmarks/verified/fleet_refs/dtc_geometry.json`` is
  internally consistent -- the DTC cloud is coincident, the KCS control is not, and the
  published reference it would be judged against still has the verdicts the library gives.

The DTC numbers themselves are not re-derived here (that would need the cached archive,
which is outside the repository by design); they are pinned by digest and re-generated with
``tools/baseline_fleet/dtc_archive_probe.py``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from tensorcad.fleet.build import RouteRefusal, _sample_extent, iges_field_mesh
from tensorcad.fleet.manifest import Provenance, ShipRecord

REPO = Path(__file__).resolve().parents[1]
ARTIFACT = REPO / "benchmarks" / "verified" / "fleet_refs" / "dtc_geometry.json"
HEX64 = re.compile(r"^[0-9a-f]{64}$")


def _record() -> ShipRecord:
    return ShipRecord(
        name="DTC",
        ship_class="container",
        length=355.0,
        beam=51.0,
        draft=14.5,
        block_coefficient=0.661,
        displacement_basis="literature_value_verified",
        provenance=Provenance(
            source="uni-due INAM DTC IGES archive",
            url="https://www.uni-due.de/imperia/md/content/inam/igs_ordner.zip",
            license="not stated by publisher",
            retrieval="cached_bytes",
            verified=True,
        ),
    )


def test_sample_extent_of_empty_and_coincident_clouds_is_zero():
    assert _sample_extent([]) == 0.0
    assert _sample_extent([(0.0, 0.0, 0.0)] * 64) == 0.0
    assert _sample_extent([(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)]) == 2.0
    assert _sample_extent([(0.0, -3.0, 0.0), (0.0, 3.0, 1.0)]) == 6.0


class _FakeSurface:
    def sample(self, spacing: float):  # pragma: no cover - only presence is used
        return []


class _FakeIges:
    def surfaces(self):
        return [_FakeSurface(), _FakeSurface()]


def test_route_refuses_a_collapsed_sample_cloud(tmp_path, monkeypatch):
    iges = tmp_path / "degenerate.igs"
    iges.write_text("not really an IGES file: the reader is stubbed")
    import tensorcad.fleet_baseline as baseline
    import tensorcad.iges as iges_module

    monkeypatch.setattr(iges_module, "read_iges", lambda path: _FakeIges())
    monkeypatch.setattr(baseline, "sample_patches", lambda surfaces, spacing: [(0.0, 0.0, 0.0)] * 1000)

    with pytest.raises(RouteRefusal) as excinfo:
        iges_field_mesh(_record(), iges, dx=0.4, dz=0.08, ring_points=40, samples_per_metre=2.5, trim=False)
    message = str(excinfo.value)
    assert "collapsed to a point" in message
    assert "1000 sample" in message
    assert "0.000e+00 m" in message


def test_artifact_pins_the_measured_reason():
    data = json.loads(ARTIFACT.read_text())
    assert data["schema"] == "tensorcad.fleet_refs.dtc_geometry/1.0"
    assert data["ship"] == "DTC"
    assert HEX64.match(data["archive"]["sha256"])
    assert HEX64.match(data["archive"]["hull_sha256"])
    for digest in data["archive"]["members"].values():
        assert HEX64.match(digest)

    census = data["entity_census"]
    assert census["type_128_surfaces"] > 0
    assert census["type_144_trimmed_patches"] == 0  # the zero cloud is not a trim-loop artefact

    dtc = data["dtc_sample_census"]
    control = data["control"]["sample_census"]
    assert dtc["surfaces_read"] == census["type_128_surfaces"]
    assert dtc["sample_points"] > 0 and dtc["unique_points"] == 1
    assert dtc["max_extent_m"] == 0.0
    assert control["sample_points"] > 0 and control["max_extent_m"] > 100.0

    assert data["route"]["outcome"] == "refused"
    assert "collapsed to a point" in data["route"]["reason"]
    assert data["status"] == "WITHHELD_READER"
    assert "refuses" in data["status_reason"]


def test_published_reference_behind_the_measurement_is_two_lineage_where_it_exists():
    data = json.loads(ARTIFACT.read_text())
    reference = data["published_reference"]
    volume = reference["displacement_volume"]
    assert volume["verdict"] == "CROSS_VERIFIED"
    assert volume["independent_families"] >= 2
    assert volume["range_relative"] < volume["tolerance"]
    # The anchor's honest gap for DTC: no second lineage for LCB/Lwl.
    assert reference["lcb_percent"]["verdict"] != "CROSS_VERIFIED"
    assert reference["lwl"]["verdict"] != "CROSS_VERIFIED"
    assert reference["cb"]["verdict"] == "CROSS_VERIFIED"