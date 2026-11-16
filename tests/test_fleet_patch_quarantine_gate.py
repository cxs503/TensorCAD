"""G9: the reader's patch quarantine must reach the acceptance verdict.

The mission that produced the fail-closed type-128 reader ("IGES type-128 reader:
fail closed on undecidable control-net layout") stopped at the reader: a patch whose
control-point layout cannot be decided is reported with ``layout_ambiguous`` and
:func:`tensorcad.fleet_baseline.partition_surfaces` quarantines it -- but the fleet's
own acceptance route (:func:`tensorcad.fleet.build.iges_field_mesh`) sampled every
surface, and nothing in G1-G8 could tell a hull with a dropped patch from a complete
one, because the half-breadth field interpolates the hole and keeps the mesh watertight.

These tests lock down both halves of the repair:

* the route refuses to *mesh* an ambiguous patch (census in ``details``), and
* the verdict refuses to *pass* a hull whose census is non-empty (G9), while a route
  that reports no census is reported as ``geometry_completeness == "unknown"`` rather
  than silently as "complete".

The synthetic hull is a barge: two bilinear patches at half-breadth 5 m spanning
x in [0, 40] m, z in [0, 4] m (keel at z = 0), plus one deliberately undecidable patch
covering the same region, whose removal changes the half-breadth field by nothing at all
(the field is a max-binning of the union of the patches) -- which is the point: the
volume gate cannot see the loss.
"""

from __future__ import annotations

import dataclasses
import math

import pytest

import tensorcad.iges
from tensorcad.fleet.build import iges_field_mesh
from tensorcad.fleet.gates import (
    QUARANTINE_TOLERANCE,
    evaluate_gates,
    gate_patch_quarantine,
    measure,
    quarantine_census,
)
from tensorcad.fleet.manifest import FleetManifest, Provenance, ShipRecord
from tensorcad.iges import BSplineSurface

L, B2, T = 40.0, 5.0, 4.0  # length, half-breadth, draft of the synthetic barge
BARGE_VOLUME = 2.0 * B2 * L * T


def _patch_parameters(x0: float, x1: float, z0: float, z1: float, *, weights: list[float]) -> list[float]:
    """A type-128 parameter list for one 2x2 bilinear patch in weight-block layout."""
    coords = [
        (x0, B2, z0),
        (x1, B2, z0),
        (x0, B2, z1),
        (x1, B2, z1),
    ]
    params: list[float] = [
        128.0, 1.0, 1.0,      # entity type, K, M
        1.0, 1.0,             # degrees
        0.0, 0.0,             # closed u/v
        0.0,                  # PROP3 = rational
        0.0, 0.0,             # PROP4, slot skipped by the reader
    ]
    params += [0.0, 0.0, 1.0, 1.0]  # u knots (count_u + degree + 1)
    params += [0.0, 0.0, 1.0, 1.0]  # v knots
    params += list(weights)
    for point in coords:
        params.extend(point)
    return params


def _clean_patch(x0: float, x1: float) -> BSplineSurface:
    surface = BSplineSurface.from_parameters(_patch_parameters(x0, x1, 0.0, T, weights=[1.0] * 4))
    assert surface.layout == "weight_block" and surface.layout_ambiguous is False
    return surface


def _undecidable_patch(x0: float, x1: float) -> BSplineSurface:
    """A 2x2 net with non-unit weights: the weight test is silent and no interior control
    point exists, so both candidate readings score 0.0 and the reader must fail closed."""
    surface = BSplineSurface.from_parameters(
        _patch_parameters(x0, x1, 0.0, T, weights=[1.0, 1.06, 0.94, 1.0])
    )
    assert surface.layout_ambiguous is True
    return surface


class _FakeIges:
    """Just enough of :class:`tensorcad.iges.IgesFile` for the iges_field route."""

    def __init__(self, surfaces: list[BSplineSurface]) -> None:
        self._surfaces = surfaces

    def surfaces(self) -> list[BSplineSurface]:
        return list(self._surfaces)

    def trimmed_surfaces(self) -> list[tuple[object, object, list[object]]]:
        return []


@pytest.fixture()
def iges_path(tmp_path):
    path = tmp_path / "barge.igs"
    path.write_text("synthetic barge (reader is stubbed by the test)\n")
    return path


def _build(monkeypatch, iges_path, surfaces, *, declared_volume: float = BARGE_VOLUME):
    monkeypatch.setattr(tensorcad.iges, "read_iges", lambda _path: _FakeIges(surfaces))
    record = ShipRecord(
        name="Barge",
        ship_class="synthetic",
        length=L,
        beam=2 * B2,
        draft=T,
        block_coefficient=BARGE_VOLUME / (L * 2 * B2 * T),
        displacement_basis="closure",
        displacement_volume=declared_volume,
        geometry_route="iges_field",
        provenance=Provenance(
            source="synthetic", url="https://example.invalid/barge",
            license="test", retrieval="test", verified=False,
        ),
    )
    built = iges_field_mesh(record, iges_path, dx=2.0, dz=0.5, samples_per_metre=1.0)
    return record, built


# --------------------------------------------------------------------------- #
# the gate itself
# --------------------------------------------------------------------------- #
class TestGateContract:
    def test_complete_census_passes(self):
        gate = gate_patch_quarantine([], patch_total=45)
        assert gate.passed is True
        assert gate.name == "G9_patch_quarantine"
        assert gate.values["quarantined_count"] == 0
        assert gate.values["quarantined_fraction"] == 0.0

    def test_one_dropped_patch_fails_by_default(self):
        gate = gate_patch_quarantine(
            [{"patch_index": 7, "layout_chosen": "interleaved",
              "reason": "ambiguous control-point layout (two candidates within 4x)"}],
            patch_total=45,
        )
        assert gate.passed is False
        assert QUARANTINE_TOLERANCE == 0
        assert gate.values["quarantined_count"] == 1
        assert math.isclose(gate.values["quarantined_fraction"], 1 / 45)
        assert "patch 7" in gate.detail and "45" in gate.detail

    def test_tolerance_must_be_asked_for_explicitly(self):
        census = [{"patch_index": i} for i in range(3)]
        assert gate_patch_quarantine(census, patch_total=45).passed is False
        assert gate_patch_quarantine(census, patch_total=45, tolerance=3).passed is True
        assert gate_patch_quarantine(census, patch_total=45, tolerance=2).passed is False

    def test_unknown_total_is_reported_not_guessed(self):
        gate = gate_patch_quarantine([{"patch_index": 0}], patch_total=None)
        assert gate.passed is False
        assert gate.values["quarantined_fraction"] is None
        assert "?" in gate.detail


class TestCensusDiscovery:
    def test_fleet_route_shape(self):
        details = {"patch_quarantine": {"patches_total": 45, "patches_used": 44,
                                        "patches_quarantined": [{"patch_index": 3}]}}
        census = quarantine_census(details)
        assert census == ([{"patch_index": 3}], 45)

    def test_legacy_build_ship_shape(self):
        details = {"source": {"patches_total": 132, "patches_used": 106,
                              "patches_quarantined": [{"patch_index": 1}, {"patch_index": 2}]}}
        assert quarantine_census(details) == ([{"patch_index": 1}, {"patch_index": 2}], 132)

    def test_absent_census_is_none_not_empty(self):
        assert quarantine_census({"route": "series60_form"}) is None
        assert quarantine_census(None) is None
        assert quarantine_census({"patch_quarantine": {"patches_total": 10}}) is None


# --------------------------------------------------------------------------- #
# the route: does it stop meshing patches the reader could not decide?
# --------------------------------------------------------------------------- #
class TestRouteQuarantines:
    def test_ambiguous_patch_is_not_meshed_and_is_counted(self, monkeypatch, iges_path):
        surfaces = [_clean_patch(0.0, 20.0), _undecidable_patch(0.0, 20.0), _clean_patch(20.0, 40.0)]
        _record, built = _build(monkeypatch, iges_path, surfaces)
        census = built.details["patch_quarantine"]
        assert census["patches_total"] == 3
        assert census["patches_used"] == 2
        assert len(census["patches_quarantined"]) == 1
        assert census["patches_quarantined"][0]["patch_index"] == 1

    def test_dropping_the_patch_leaves_the_volume_unchanged(self, monkeypatch, iges_path):
        """The measured fail-open: a lost patch is invisible to the volume gate.

        The undecidable patch duplicates surface the two clean patches already carry, so
        the half-breadth field, and therefore the mesh volume, is identical with and
        without it.  Nothing in G1-G8 can separate the two hulls; only a census can.
        """
        with_patch = _build(
            monkeypatch, iges_path,
            [_clean_patch(0.0, 20.0), _undecidable_patch(0.0, 20.0), _clean_patch(20.0, 40.0)],
        )[1]
        without_patch = _build(
            monkeypatch, iges_path,
            [_clean_patch(0.0, 20.0), _clean_patch(20.0, 40.0)],
        )[1]
        v_with = measure(with_patch.mesh).signed_volume
        v_without = measure(without_patch.mesh).signed_volume
        assert v_with == pytest.approx(v_without, rel=1e-12)
        # ... and the mesh is a barge of the declared size.  The 6% band is the field/loft
        # discretisation at dx = 2 m (the square ends taper over one cell), not slack for
        # the patch census: the two volumes above agree to 1e-12.
        assert v_with == pytest.approx(BARGE_VOLUME, rel=0.06)

    def test_all_patches_quarantined_refuses_the_file(self, monkeypatch, iges_path):
        from tensorcad.fleet.build import RouteRefusal

        with pytest.raises(RouteRefusal, match="every one of the 1 B-spline surfaces"):
            _build(monkeypatch, iges_path, [_undecidable_patch(0.0, 20.0)])


# --------------------------------------------------------------------------- #
# the verdict: does the census reach the gate chain?
# --------------------------------------------------------------------------- #
@pytest.fixture()
def barge_case(monkeypatch, iges_path):
    """One hull, meshed once; G7 is satisfied by construction so G9 is the only variable.

    The declared volume is set to the volume the field route measures (``closure``
    basis), which isolates the question these tests ask: does the *census* change the
    verdict?  Hydrostatic accuracy is G7's own subject and is covered elsewhere.
    """
    surfaces = [_clean_patch(0.0, 20.0), _undecidable_patch(0.0, 20.0), _clean_patch(20.0, 40.0)]
    record, built = _build(monkeypatch, iges_path, surfaces)
    measured = measure(built.mesh).signed_volume
    # keep the particulars internally consistent (G6/G6b) so that G9 is the only gate
    # whose outcome can differ between the two verdicts below
    record = dataclasses.replace(
        record, displacement_volume=measured, block_coefficient=measured / (L * 2 * B2 * T)
    )
    manifest = FleetManifest(ships=(record,))
    return record, built, manifest


class TestVerdictConsumesCensus:
    def test_census_flips_the_verdict(self, barge_case):
        record, built, manifest = barge_case
        dropped, total = quarantine_census(built.details)
        assert len(dropped) == 1 and total == 3

        without = evaluate_gates(record, built.mesh, manifest, volume_tolerance=0.01)
        assert without["passed"] is True, without["failed_gates"]
        assert "G9_patch_quarantine" not in [g["name"] for g in without["gates"]]
        assert without["geometry_completeness"] == "unknown"

        with_census = evaluate_gates(
            record, built.mesh, manifest, volume_tolerance=0.01,
            quarantine=dropped, patch_total=total,
        )
        assert with_census["passed"] is False
        assert with_census["failed_gates"] == ["G9_patch_quarantine"]
        assert with_census["geometry_completeness"] == "patches_dropped"
        assert with_census["patches_quarantined"]["count"] == 1
        assert with_census["patches_quarantined"]["total"] == 3

    def test_empty_census_is_complete_and_passes(self, barge_case):
        record, built, manifest = barge_case
        verdict = evaluate_gates(
            record, built.mesh, manifest, volume_tolerance=0.01, quarantine=[], patch_total=3,
        )
        names = [g["name"] for g in verdict["gates"]]
        assert "G9_patch_quarantine" in names
        assert verdict["passed"] is True, verdict["failed_gates"]
        assert verdict["geometry_completeness"] == "complete"

class TestCallSiteForwardsCensus:
    """``fleet.run._ship_block`` must hand the census to the gate chain.

    The gate is only as good as the call site: a builder that reports a census nobody
    reads is the same fail-open as no census at all.  ``build_ship_geometry`` is
    replaced here so the test pins the plumbing, not the IGES reader.
    """

    def test_ship_block_reports_dropped_patches(self, barge_case, monkeypatch):
        from tensorcad.fleet import run as fleet_run
        from tensorcad.fleet.build import BuiltGeometry

        record, built, manifest = barge_case
        census = built.details["patch_quarantine"]

        def _fake_build(*_args, **_kwargs):
            return BuiltGeometry(mesh=built.mesh, route="iges_field",
                                 details={"patch_quarantine": census})

        monkeypatch.setattr(fleet_run, "build_ship_geometry", _fake_build)
        block = fleet_run._ship_block(
            record, manifest, nx=41, nz=21, volume_tolerance=0.01, convergence=False,
        )
        assert block["geometry_completeness"] == "patches_dropped"
        assert "G9_patch_quarantine" in block["failed_gates"]
        assert block["patches_quarantined"]["count"] == 1

    def test_ship_block_without_census_says_unknown(self, barge_case, monkeypatch):
        from tensorcad.fleet import run as fleet_run
        from tensorcad.fleet.build import BuiltGeometry

        record, built, manifest = barge_case
        monkeypatch.setattr(
            fleet_run, "build_ship_geometry",
            lambda *_a, **_k: BuiltGeometry(mesh=built.mesh, route="series60_form", details={}),
        )
        block = fleet_run._ship_block(
            record, manifest, nx=41, nz=21, volume_tolerance=0.01, convergence=False,
        )
        assert block["geometry_completeness"] == "unknown"
        assert "G9_patch_quarantine" not in [g["name"] for g in block["gates"]]
