"""Closed-form hydrostatics against the mesh integrator that has to reproduce them."""
from __future__ import annotations

import importlib.util
import json
import math
import pathlib

import pytest

from tensorcad.hull import (
    TAPER_INTEGRAL,
    HullParameters,
    analytic_volume,
    generate_hull,
)
from tensorcad.hydrostatics_analytic import (
    AnalyticHydrostatics,
    analytic_hydrostatics,
    quadrature_evidence,
    tanh_sinh_quadrature,
    taper_integral,
)

#: Resolution ladder used by several tests: (nx, nz) both doubling.
LADDER = ((41, 13), (81, 21), (161, 41))


def _model(nx: int, nz: int):
    return generate_hull(HullParameters(nx=nx, nz=nz))


# ---------------------------------------------------------------- the taper integral


def test_taper_integral_agrees_with_the_volume_derived_constant() -> None:
    """Two independent derivations: gamma function here, volume identity in hull.py."""
    assert taper_integral() == pytest.approx(TAPER_INTEGRAL, rel=0.0, abs=1e-15)


def test_taper_integral_agrees_with_quadrature() -> None:
    numeric = tanh_sinh_quadrature(lambda x: math.sin(math.pi * x) ** 0.7, 0.0, 1.0, levels=5)
    assert taper_integral() == pytest.approx(numeric, rel=1e-12)


@pytest.mark.parametrize(
    ("exponent", "expected"),
    [(0.0, 1.0), (1.0, 2.0 / math.pi), (2.0, 0.5)],
)
def test_taper_integral_reproduces_textbook_values(exponent: float, expected: float) -> None:
    """e=0 is a rectangle, e=1 is 2/pi, e=2 is a half sine (1/2)."""
    assert taper_integral(exponent) == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("exponent", [-1.0, -2.0, float("nan"), float("inf")])
def test_taper_integral_rejects_divergent_exponents(exponent: float) -> None:
    with pytest.raises(ValueError):
        taper_integral(exponent)


# ---------------------------------------------------------------------- quadrature


def test_tanh_sinh_is_exact_on_a_cubic() -> None:
    assert tanh_sinh_quadrature(lambda x: x**3, 0.0, 1.0, levels=4) == pytest.approx(0.25, abs=1e-15)


def test_tanh_sinh_handles_the_keel_style_endpoint_singularity() -> None:
    """``x**-0.5`` is the germ of the keel edge; the rule must not need a special case."""
    assert tanh_sinh_quadrature(lambda x: x**-0.5, 0.0, 1.0, levels=6) == pytest.approx(2.0, abs=1e-7)


def test_tanh_sinh_on_a_non_unit_interval() -> None:
    assert tanh_sinh_quadrature(math.sin, 0.0, math.pi, levels=5) == pytest.approx(2.0, rel=1e-14)


def test_quadrature_evidence_reports_a_vanishing_difference() -> None:
    evidence = quadrature_evidence(math.exp, 0.0, 1.0, levels=5)
    assert evidence["value"] == pytest.approx(math.e - 1.0, rel=1e-14)
    assert evidence["relative_difference"] < 1e-12


@pytest.mark.parametrize("interval", [(1.0, 1.0), (1.0, 0.0)])
def test_tanh_sinh_rejects_an_empty_interval(interval: tuple[float, float]) -> None:
    with pytest.raises(ValueError):
        tanh_sinh_quadrature(math.sin, *interval)


def test_tanh_sinh_rejects_a_zero_level_count() -> None:
    with pytest.raises(ValueError):
        tanh_sinh_quadrature(math.sin, 0.0, 1.0, levels=0)


# ------------------------------------------------------------------- closed forms


@pytest.fixture(scope="module")
def reference() -> AnalyticHydrostatics:
    return analytic_hydrostatics()


def test_the_reference_is_a_closed_form_of_the_surface(reference: AnalyticHydrostatics) -> None:
    p = HullParameters()
    taper = taper_integral()
    assert reference.length == p.length
    assert reference.beam == p.beam
    assert reference.draft == p.draft
    assert reference.taper_integral == pytest.approx(taper, rel=1e-15)
    assert reference.displacement_volume == pytest.approx(
        p.beam * p.length * p.draft * taper * 2.0 / 3.0, rel=1e-15
    )
    assert reference.waterplane_area == pytest.approx(p.beam * p.length * taper, rel=1e-15)
    assert reference.maximum_section_area == pytest.approx(2.0 * p.beam * p.draft / 3.0, rel=1e-15)
    assert reference.station_of_maximum_section == pytest.approx(0.5 * p.length)


def test_the_coefficients_are_consistent_with_each_other(reference: AnalyticHydrostatics) -> None:
    """Cb = Cp * Cm and Cw = Cp: the identities hold exactly, not to a tolerance."""
    assert reference.block_coefficient == pytest.approx(
        reference.prismatic_coefficient * reference.midship_coefficient, rel=1e-15
    )
    assert reference.waterplane_coefficient == pytest.approx(reference.prismatic_coefficient, rel=1e-15)
    assert reference.block_coefficient == pytest.approx(2.0 * taper_integral() / 3.0, rel=1e-15)
    assert reference.midship_coefficient == pytest.approx(2.0 / 3.0, rel=1e-15)


def test_the_centres_follow_from_the_symmetry_of_the_form(reference: AnalyticHydrostatics) -> None:
    p = HullParameters()
    assert reference.longitudinal_centre_buoyancy == pytest.approx(0.5 * p.length, rel=1e-15)
    assert reference.longitudinal_centre_flotation == pytest.approx(0.5 * p.length, rel=1e-15)
    # z-weighted mean of zeta**0.5 on [-T, 0] is -0.4*T, i.e. 0.6*T above the keel.
    assert reference.vertical_centre_buoyancy == pytest.approx(-0.4 * p.draft, rel=1e-15)
    assert reference.vertical_centre_buoyancy + p.draft == pytest.approx(0.6 * p.draft, rel=1e-15)


def test_water_density_scales_the_mass_and_nothing_else() -> None:
    sea = analytic_hydrostatics(water_density=1025.0)
    fresh = analytic_hydrostatics(water_density=1000.0)
    assert fresh.displacement_volume == pytest.approx(sea.displacement_volume, rel=1e-15)
    assert fresh.displacement_mass / sea.displacement_mass == pytest.approx(1000.0 / 1025.0, rel=1e-14)


def test_the_wetted_area_carries_its_own_convergence_evidence(reference: AnalyticHydrostatics) -> None:
    evidence = reference.quadrature
    assert evidence["wetted_area_relative_difference"] < 1e-10
    assert evidence["method"].startswith("tanh-sinh")


def test_the_wetted_area_beats_the_area_estimate_a_flat_section_would_give(
    reference: AnalyticHydrostatics,
) -> None:
    """Sanity bound: the true wetted area sits between the waterplane and 2x the box side."""
    p = HullParameters()
    assert reference.waterplane_area < reference.wetted_surface_area
    assert reference.wetted_surface_area < 3.0 * p.length * p.draft


def test_the_reference_ignores_the_sampling_resolution() -> None:
    coarse = analytic_hydrostatics(HullParameters(nx=11, nz=5))
    fine = analytic_hydrostatics(HullParameters(nx=161, nz=41))
    assert coarse.as_dict() == fine.as_dict()


@pytest.mark.parametrize(
    "parameters",
    [HullParameters(length=-1.0, beam=2.0, draft=0.8), (10.0, 0.0, 0.8), (10.0, 2.0, float("nan"))],
)
def test_invalid_dimensions_are_rejected(parameters) -> None:
    with pytest.raises(ValueError):
        analytic_hydrostatics(parameters)


def test_a_plain_triple_is_accepted() -> None:
    assert analytic_hydrostatics((10.0, 2.0, 0.8)).displacement_volume == pytest.approx(
        analytic_hydrostatics().displacement_volume, rel=1e-15
    )


def test_an_unusable_argument_names_the_expected_types() -> None:
    with pytest.raises(TypeError, match="HullParameters"):
        analytic_hydrostatics("not a hull")


# ------------------------------------------- the mesh integrator against the reference


def test_the_default_grid_is_within_three_percent_on_every_quantity() -> None:
    """The audit criterion, applied to the coefficients the tests used to only bound."""
    reference = analytic_hydrostatics()
    comparison = reference.compare(_model(41, 13).hydrostatics)
    assert comparison["worst_relative_error"] < 0.03
    for name, row in comparison["rows"].items():
        assert abs(row["relative_error"]) < 0.03, f"{name} is {row['relative_error_percent']:+.2f} %"


def test_the_mesh_under_integrates_every_quantity() -> None:
    """An inscribed polygonal surface can only lose volume, area and area-moment."""
    comparison = analytic_hydrostatics().compare(_model(41, 13).hydrostatics)
    for name in (
        "displacement_volume",
        "waterplane_area",
        "wetted_surface_area",
        "maximum_section_area",
        "vertical_centre_buoyancy",
        "block_coefficient",
        "prismatic_coefficient",
    ):
        assert comparison["rows"][name]["relative_error"] < 0.0, name


def test_the_fine_grid_is_within_half_a_percent_on_every_quantity() -> None:
    comparison = analytic_hydrostatics().compare(_model(161, 41).hydrostatics)
    assert comparison["worst_relative_error"] < 0.005
    assert comparison["worst_quantity"]["relative_error"] < 0.0


def test_the_volume_error_is_what_hull_py_already_claimed() -> None:
    """The new module must reproduce the old reference, not replace it with a different one."""
    p = HullParameters()
    measured = _model(p.nx, p.nz).hydrostatics.displacement_volume
    old = analytic_volume(p.length, p.beam, p.draft)
    new = analytic_hydrostatics(p).displacement_volume
    assert old == pytest.approx(new, rel=1e-15)
    assert measured / old - 1.0 == pytest.approx(-0.0088, abs=5e-5)


@pytest.mark.parametrize(
    "name",
    [
        "displacement_volume",
        "waterplane_area",
        "wetted_surface_area",
        "maximum_section_area",
        "vertical_centre_buoyancy",
        "block_coefficient",
    ],
)
def test_every_quantity_converges_monotonically_down_the_ladder(name: str) -> None:
    reference = analytic_hydrostatics()
    errors = []
    for nx, nz in LADDER:
        row = reference.compare(_model(nx, nz).hydrostatics)["rows"][name]
        errors.append(abs(row["relative_error"]))
    assert errors[0] > errors[1] > errors[2], errors
    assert errors[2] < 0.25 * errors[0], errors


def test_the_longitudinal_centres_are_exact_on_every_grid() -> None:
    """Both centres are fixed by symmetry, so discretisation cannot move them."""
    reference = analytic_hydrostatics()
    for nx, nz in LADDER:
        comparison = reference.compare(_model(nx, nz).hydrostatics)
        assert abs(comparison["rows"]["longitudinal_centre_buoyancy"]["relative_error"]) < 1e-12
        assert abs(comparison["rows"]["longitudinal_centre_flotation"]["relative_error"]) < 1e-12


def test_the_wetted_area_reference_is_confirmed_by_richardson_extrapolation() -> None:
    """Independent check: extrapolate the mesh series; it must land on the quadrature value.

    The mesh is an inscribed surface, so its area is a decreasing sequence that
    brackets the truth from below. Estimating the convergence order from the last
    two steps and extrapolating is a reference the quadrature never saw.
    """
    reference = analytic_hydrostatics().wetted_surface_area
    middle = _model(81, 21).hydrostatics.wetted_surface_area
    fine = _model(161, 41).hydrostatics.wetted_surface_area
    finest = _model(321, 81).hydrostatics.wetted_surface_area
    previous_step = fine - middle
    last_step = finest - fine
    assert 0.0 < last_step < previous_step
    order = math.log2(previous_step / last_step)
    assert 1.0 < order < 2.5  # slower than the smooth-surface O(h^2): sharp ends
    extrapolated = finest + last_step / (2.0**order - 1.0)
    assert extrapolated == pytest.approx(reference, rel=1e-3)


def test_compare_reports_both_directions_of_error() -> None:
    class Bogus:
        displacement_volume = 1.0
        displacement_mass = 1.0
        waterplane_area = 1.0
        wetted_surface_area = 1.0
        maximum_section_area = 1.0
        longitudinal_centre_buoyancy = 1.0
        longitudinal_centre_flotation = 1.0
        vertical_centre_buoyancy = 1.0
        block_coefficient = 1.0
        prismatic_coefficient = 1.0

    comparison = analytic_hydrostatics().compare(Bogus(), quantities=["displacement_volume"])
    row = comparison["rows"]["displacement_volume"]
    assert row["reference"] == pytest.approx(analytic_hydrostatics().displacement_volume)
    assert row["relative_error"] == pytest.approx(1.0 / row["reference"] - 1.0)


def test_compare_rejects_an_unknown_quantity() -> None:
    with pytest.raises(ValueError, match="unknown hydrostatic quantity"):
        analytic_hydrostatics().compare(_model(41, 13).hydrostatics, quantities=["metacentric_height"])


def test_within_tolerance_is_the_criterion_in_one_call() -> None:
    reference = analytic_hydrostatics()
    assert reference.within(_model(161, 41).hydrostatics, tolerance=0.005)
    assert not reference.within(_model(41, 13).hydrostatics, tolerance=0.0001)


def test_the_payload_is_json_ready(reference: AnalyticHydrostatics) -> None:
    import json

    payload = json.loads(json.dumps(reference.as_dict()))
    assert payload["schema"] == "tensorcad.hydrostatics.analytic/1"
    assert payload["displacement_volume_m3"] == pytest.approx(reference.displacement_volume)
    # Keys line up with Hydrostatics.to_dict() so a diff between the two is mechanical.
    mesh_payload = _model(41, 13).hydrostatics.to_dict()
    shared = {"displacement_volume_m3", "waterplane_area_m2", "block_coefficient"}
    assert shared <= set(payload)
    assert shared <= set(mesh_payload)

# ------------------------------------------------- the verification example's claims

EXAMPLE = (
    pathlib.Path(__file__).resolve().parents[1]
    / "examples"
    / "hydrostatics_verification"
    / "run_verification.py"
)


def _load_example():
    spec = importlib.util.spec_from_file_location("hydrostatics_verification", EXAMPLE)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_example_fits_the_second_order_it_should() -> None:
    example = _load_example()
    # Halving h quarters the error: order 2.
    assert example._order(0.5, 1.0, 0.25, 1.0) == pytest.approx(2.0)
    assert math.isnan(example._order(1.0, 1.0, 0.1, 1.0))


def test_the_example_extrapolation_recovers_a_known_limit() -> None:
    """Richardson on a synthetic O(h**2) series must land on the limit it is given."""
    example = _load_example()
    limit = 22.912729522831285
    rows = [
        {"rows": {"wetted_surface_area": {"measured": limit - 0.4 * h * h, "reference": limit}}}
        for h in (1.0, 0.5, 0.25)
    ]
    recovered = example.richardson_wetted_area(rows)
    assert recovered["observed_order"] == pytest.approx(2.0, rel=1e-9)
    assert recovered["extrapolated"] == pytest.approx(limit, rel=1e-12)


def test_the_example_runs_on_a_shortened_ladder(tmp_path, monkeypatch) -> None:
    example = _load_example()
    monkeypatch.setattr(example, "LADDER", ((21, 9), (41, 13), (81, 21)))
    # A three-level ladder is too short for the asymptotic claims (the Richardson
    # order estimate has not settled), so the exit code is 0 or 1 by design; the
    # structural contract and the default-grid claim are what this test pins.
    assert example.main(["--output", str(tmp_path)]) in (0, 1)
    payload = json.loads((tmp_path / "hydrostatics_verification.json").read_text())
    assert payload["schema"] == "tensorcad.verification.hydrostatics/1"
    assert set(payload["claims"]) == {
        "default_grid_within_3_percent",
        "fine_grid_within_half_a_percent",
        "richardson_confirms_wetted_area",
        "quadrature_self_convergence",
        "monotone_convergence",
    }
    assert payload["claims"]["default_grid_within_3_percent"]["ok"]
    assert payload["claims"]["quadrature_self_convergence"]["ok"]
    assert len(payload["ladder"]) == 3
