"""The LCB channel through the ShipD rebuild pipeline (round 8).

Round 7 built the LCB control into ``cb_sections`` and verified it on a coefficient
scan; round 8 wires it into the campaign. These tests pin the three properties the
wiring has to hold:

1. the parameter table gains an LCB column *without* disturbing any historical column
   (the 500-ship pass rates of rounds 3/5 must stay comparable);
2. a vector that carries no LCB target rebuilds to a byte-identical surface, so the
   symmetric family is untouched;
3. a vector that *does* carry one is either honoured on the mesh or refused with a
   classified reason -- never silently dropped.
"""

from __future__ import annotations

import hashlib
import json

import pytest

from tensorcad.cb_sections import LCB_OFFSET_LIMIT
from tensorcad.shipd_population import (
    LCB_FRACTION_RANGE,
    build_shipd_database,
    nominal_performance_fields,
)
from tensorcad.shipd_rebuild import LCB_TOLERANCE, rebuild_ship, run_shipd_rebuild

NX, NZ = 41, 13   # small grid; every assertion below is grid-exact, not grid-dependent


def _digest(mesh) -> str:
    payload = repr((tuple(mesh.vertices), tuple(mesh.faces))).encode()
    return hashlib.sha256(payload).hexdigest()


# --------------------------------------------------------------------------- #
# 1. Parameter table
# --------------------------------------------------------------------------- #
def test_population_carries_lcb_inside_the_sampled_band():
    ships = build_shipd_database(40, seed=42)
    assert all("lcb" in ship for ship in ships)
    for ship in ships:
        low, high = LCB_FRACTION_RANGE[ship["ship_type"]]
        assert low <= ship["lcb"] <= high
        # The API guardrail is the hard bound; the type band must sit inside it.
        assert abs(ship["lcb"]) < LCB_OFFSET_LIMIT


def test_historical_columns_are_undisturbed_by_the_lcb_column():
    """Round 3/5 pass rates are only comparable if the old columns did not move."""
    expected = {
        1: dict(ship_type="container", L=280.85, B=43.42, T=17.32,
                Cb=0.579, Cm=0.9626, Cwp=0.7229, design_froude=0.2535),
        2: dict(ship_type="tanker", L=190.23, B=32.11, T=10.85,
                Cb=0.7975, Cm=0.9888, Cwp=0.8555, design_froude=0.128),
        3: dict(ship_type="bulk_carrier", L=235.36, B=37.69, T=13.6,
                Cb=0.7907, Cm=0.9942, Cwp=0.8506, design_froude=0.1522),
    }
    ships = build_shipd_database(3, seed=42)
    for ship in ships:
        assert ship["ship_type"] == expected[ship["ship_id"]]["ship_type"]
        for key, want in expected[ship["ship_id"]].items():
            if key == "ship_type":
                continue
            assert ship[key] == pytest.approx(want, rel=0, abs=1e-9), (ship["ship_id"], key)


def test_population_is_deterministic_with_lcb():
    assert build_shipd_database(12, seed=7) == build_shipd_database(12, seed=7)


# --------------------------------------------------------------------------- #
# 2. The default path is bit-for-bit unchanged
# --------------------------------------------------------------------------- #
def test_symmetric_family_is_byte_identical_when_lcb_absent_or_zero():
    for ship in build_shipd_database(4, seed=42):
        bare = {k: v for k, v in ship.items() if k != "lcb"}
        zero = dict(bare, lcb=0.0)
        a = rebuild_ship(bare, nx=NX, nz=NZ)
        b = rebuild_ship(zero, nx=NX, nz=NZ)
        assert a.mesh is not None and b.mesh is not None
        assert _digest(a.mesh) == _digest(b.mesh)
        # No target -> no LCB gate, and the realised centroid is zero to round-off.
        assert a.lcb_gates == []
        assert abs(a.realised_lcb) < 1e-12


# --------------------------------------------------------------------------- #
# 3. A requested LCB is honoured or classified, never dropped
# --------------------------------------------------------------------------- #
def test_requested_lcb_is_realised_on_the_mesh():
    ships = build_shipd_database(12, seed=42)
    for ship in ships:
        result = rebuild_ship(ship, nx=NX, nz=NZ)
        assert result.mesh is not None
        assert result.requested_lcb == pytest.approx(ship["lcb"], abs=1e-12)
        assert result.lcb_error is not None and result.lcb_error < LCB_TOLERANCE
        assert result.lcb_gates and result.lcb_gates[0].name == "G11_longitudinal_centre"
        assert result.lcb_gates[0].passed
        assert result.strict_verdict == "VERIFIED"
        # Measured, not echoed: the sign must follow the request.
        assert (result.realised_lcb - ship["lcb"]) * (1 if ship["lcb"] >= 0 else -1) == \
            pytest.approx(0.0, abs=LCB_TOLERANCE)


def test_guardrail_rejects_an_out_of_band_request():
    ship = dict(ship_id=1, ship_type="tanker", L=320.0, B=58.0, T=20.8,
                Cb=0.82, Cm=0.9956, Cwp=0.88, lcb=0.06)
    result = rebuild_ship(ship, nx=NX, nz=NZ)
    assert result.verdict == "REFUSED"
    assert result.mesh is None
    assert "lcb_offset must be a finite fraction of L" in result.error_message


def test_unreachable_lcb_is_refused_and_classified():
    """The fullest hulls cannot put their centroid at +4.9 % L (measured, see report)."""
    ship = dict(ship_id=322, ship_type="tanker", L=181.07, B=29.26, T=8.7,
                Cb=0.8591, Cm=0.9926, Cwp=0.8987,
                lcb=LCB_OFFSET_LIMIT - 1.0e-3)
    result = rebuild_ship(ship, nx=NX, nz=NZ)
    assert result.verdict == "REFUSED"
    assert result.mesh is None                      # nothing built -> nothing to measure
    assert result.error_message.startswith("lcb_unreachable:")
    assert result.requested_lcb == pytest.approx(ship["lcb"])


def test_legacy_methods_report_lcb_as_unsupported():
    """``Cp > Cwp`` routes to the parametric chain, which has no LCB degree of freedom."""
    ship = dict(ship_id="ONRT", ship_type="ONRT", L=154.0, B=18.5, T=5.5,
                Cb=0.535, Cm=0.75, Cwp=0.69, lcb=-0.01)
    result = rebuild_ship(ship, nx=NX, nz=NZ)
    assert result.method == "parametric_sections"
    assert result.lcb_gates and result.lcb_gates[0].name == "G11_longitudinal_centre"
    assert not result.lcb_gates[0].passed
    assert "no LCB control" in result.lcb_gates[0].detail
    assert result.strict_verdict == "FAILED"


# --------------------------------------------------------------------------- #
# 4. Campaign-level accounting
# --------------------------------------------------------------------------- #
def test_run_records_the_lcb_contract(tmp_path):
    ships = build_shipd_database(6, seed=42)
    evidence = run_shipd_rebuild(ships, nx=NX, nz=NZ, out_dir=tmp_path)
    contract = evidence["lcb_contract"]
    assert contract["n_requested"] == 6
    assert contract["n_honoured"] + contract["n_refused_unreachable"] == 6
    assert contract["max_abs_error_L"] < LCB_TOLERANCE
    summary = json.loads((tmp_path / "lcb_contract.json").read_text())
    assert len(summary["rows"]) == 6
    assert all(row["lcb_gate_passed"] for row in summary["rows"])


def test_parallel_and_serial_runs_agree(tmp_path):
    ships = build_shipd_database(6, seed=42)
    serial = run_shipd_rebuild(ships, nx=NX, nz=NZ, out_dir=tmp_path / "serial")
    parallel = run_shipd_rebuild(ships, nx=NX, nz=NZ, out_dir=tmp_path / "parallel",
                                 workers=3)
    for key in ("pass_rate", "strict_pass_rate", "lcb_contract"):
        assert serial[key] == parallel[key]


def test_nominal_fields_tolerate_a_missing_froude_number():
    """Benchmark hulls carry ``design_froude: null``; the geometry fields must survive."""
    ship = dict(L=230.0, B=32.2, T=10.8, Cb=0.651, Cm=0.971, Cwp=0.767,
                design_froude=None)
    fields = nominal_performance_fields(ship)
    assert fields["design_froude"] == pytest.approx(0.20)
    assert fields["displacement_volume_m3"] > 0.0
    assert fields["wetted_surface_m2"] > 0.0