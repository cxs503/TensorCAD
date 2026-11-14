"""Tests for :mod:`tensorcad.shipd_population`.

The population model is the input of the whole rebuild campaign, so what matters
is that it is (a) reproducible, (b) inside the published ranges it claims, and
(c) feasible for the reconstruction family that consumes it. The nominal
performance block is checked against its own physics: the form factor must be a
factor (``> 1``), the friction coefficient must sit in the ITTC-57 band, and a
longer ship at the same Froude number must carry a lower ``Cf``.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.cb_sections import CalibratedHullParameters
from tensorcad.shipd_population import (
    SHIP_TYPE_POPULATION,
    build_shipd_database,
    correlation_sanity,
    holtrop_mennen_form_factor,
    holtrop_mennen_wetted_surface,
    ittc57_friction_coefficient,
    nominal_performance_fields,
)


def test_population_is_deterministic():
    first = build_shipd_database(40, seed=7)
    again = build_shipd_database(40, seed=7)
    other = build_shipd_database(40, seed=8)
    assert first == again
    assert first != other


def test_population_covers_every_ship_type():
    ships = build_shipd_database(500)
    assert len(ships) == 500
    assert {ship["ship_type"] for ship in ships} == set(SHIP_TYPE_POPULATION)
    ids = [ship["ship_id"] for ship in ships]
    assert ids == list(range(1, 501))


def test_population_respects_its_ranges():
    for ship in build_shipd_database(500):
        l_min, l_max, lb, bt, cb_min, cb_max, _fn_min, _fn_max = SHIP_TYPE_POPULATION[
            ship["ship_type"]
        ]
        assert l_min <= ship["L"] <= l_max
        assert cb_min <= ship["Cb"] <= cb_max
        assert 0.0 < ship["Cm"] < 1.0
        assert 0.0 < ship["Cwp"] < 1.0
        assert ship["L"] / ship["B"] == pytest.approx(lb, rel=0.15)
        assert ship["B"] / ship["T"] == pytest.approx(bt, rel=0.15)


def test_every_sampled_ship_is_reconstructable():
    """No ship in the population may be refused by the calibrated family."""
    for ship in build_shipd_database(500):
        params = CalibratedHullParameters(
            ship["L"], ship["B"], ship["T"], ship["Cb"], ship["Cm"], ship["Cwp"]
        )
        assert params.problems() == [], ship


def test_invalid_population_size():
    with pytest.raises(ValueError):
        build_shipd_database(0)


def test_wetted_surface_grows_with_size():
    small = {"L": 100.0, "B": 16.0, "T": 5.0, "Cb": 0.65, "Cm": 0.98, "Cwp": 0.77}
    large = {**small, "L": 200.0, "B": 32.0, "T": 10.0}
    ratio = holtrop_mennen_wetted_surface(large) / holtrop_mennen_wetted_surface(small)
    assert ratio == pytest.approx(4.0, rel=0.05)


def test_form_factor_is_a_factor_above_one():
    for cb in (0.45, 0.60, 0.75, 0.85):
        ship = {"L": 200.0, "B": 30.0, "T": 10.0, "Cb": cb, "Cm": 0.98, "Cwp": 0.80}
        factor = holtrop_mennen_form_factor(ship)
        assert factor > 1.0
        # Fuller hulls carry the larger form factor.
    fine = holtrop_mennen_form_factor(
        {"L": 200.0, "B": 30.0, "T": 10.0, "Cb": 0.45, "Cm": 0.90, "Cwp": 0.70}
    )
    full = holtrop_mennen_form_factor(
        {"L": 200.0, "B": 30.0, "T": 10.0, "Cb": 0.85, "Cm": 0.995, "Cwp": 0.90}
    )
    assert full > fine


def test_ittc57_friction_line():
    assert ittc57_friction_coefficient(1e9) == pytest.approx(0.0015, abs=1e-4)
    assert ittc57_friction_coefficient(1e6) > ittc57_friction_coefficient(1e9)
    with pytest.raises(ValueError):
        ittc57_friction_coefficient(0.0)


def test_nominal_fields_are_physical():
    ship = build_shipd_database(1)[0]
    fields = nominal_performance_fields(ship)
    assert fields["displacement_volume_m3"] == pytest.approx(
        ship["Cb"] * ship["L"] * ship["B"] * ship["T"], rel=1e-12
    )
    assert fields["design_speed_mps"] == pytest.approx(
        ship["design_froude"] * math.sqrt(9.80665 * ship["L"]), rel=1e-12
    )
    assert fields["reynolds"] > 1e7
    assert 0.0005 < fields["friction_coefficient"] < 0.004
    assert fields["form_factor_1pk"] > 1.0
    assert fields["viscous_resistance_kN"] > 0.0
    assert fields["prismatic_coefficient"] == pytest.approx(ship["Cb"] / ship["Cm"], rel=1e-12)


def test_correlation_sanity_reports_nothing_for_empty_input():
    assert correlation_sanity([], [])["n_pairs"] == 0