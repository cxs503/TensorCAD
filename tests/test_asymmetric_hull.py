"""The asymmetric family is verified against its own analytic surface, and that
reference is verified against the shipped symmetric one.

Three independent layers, each checked by a different mechanism:

* the closed-form taper integral ``I_p`` against the 0.7 constant in
  :mod:`tensorcad.hull` and against a plain Simpson rule;
* the ``M_p`` moment quadrature against a second substitution power and against
  node refinement (it is the only piece of the reference without a formula);
* the reference as a whole against the mesh path, which integrates the same
  surface with the divergence theorem on triangles -- an independent integrator,
  not a re-run of the quadrature.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import hull
from tensorcad.asymmetric_hull import (
    COMPARABLE_QUANTITIES,
    MAX_EXPONENT,
    MIN_EXPONENT,
    SCHEMA,
    AsymmetricHullParameters,
    analytic_half_taper_integral,
    analytic_hydrostatics,
    analytic_longitudinal_centre,
    analytic_relative_errors,
    analytic_volume,
    analytic_waterplane_area,
    analytic_wetted_surface_area,
    build_mesh,
    generate_hull,
    half_taper_integrals,
    taper,
    wetted_surface_quadrature_drift,
)
from tensorcad.hydrostatics import compute_hydrostatics

DEFAULT = AsymmetricHullParameters()
GRIDS = ((41, 13), (81, 21), (161, 41))


def simpson(function, lower: float, upper: float, intervals: int = 40000) -> float:
    """Plain composite Simpson rule, independent of the shipped quadrature."""
    if intervals % 2:
        intervals += 1
    step = (upper - lower) / intervals
    total = function(lower) + function(upper)
    for k in range(1, intervals):
        weight = 4.0 if k % 2 else 2.0
        total += weight * function(lower + k * step)
    return total * step / 3.0


# --- the reference's own pieces -------------------------------------------------


def test_closed_form_reproduces_the_shipped_taper_constant():
    # hull.TAPER_INTEGRAL is int_0^1 sin(pi*x)^0.7 dx, hard-coded in 0.11.0.
    # The new Beta-function form must hit half of it, not merely be close.
    assert analytic_half_taper_integral(0.7) * 2.0 == pytest.approx(
        hull.TAPER_INTEGRAL, rel=1e-14
    )


@pytest.mark.parametrize("exponent", [0.4, 0.7, 1.0, 1.6, 2.0, 3.0])
def test_closed_form_matches_independent_simpson(exponent):
    closed = analytic_half_taper_integral(exponent)
    numeric = simpson(lambda x: math.sin(math.pi * x) ** exponent, 0.0, 0.5)
    # Simpson is the looser of the two: the integrand has a singular derivative at
    # the tip for p < 1, so its error decays algebraically, not spectrally.
    assert closed == pytest.approx(numeric, rel=1e-7)


@pytest.mark.parametrize("exponent", [0.4, 0.7, 1.6, 2.5, 3.0])
def test_moment_quadrature_is_stable_under_the_substitution(exponent):
    # Two different substitution powers for the same integral: agreement here is
    # evidence that the tau^m mapping is not shaping the answer.
    x_2p, m_2p = half_taper_integrals(exponent, nodes=120)
    x_3p, m_3p = half_taper_integrals(exponent, nodes=120, substitution=3.0 / exponent)
    assert x_2p == pytest.approx(x_3p, rel=1e-9)
    assert m_2p == pytest.approx(m_3p, rel=1e-7)


@pytest.mark.parametrize("exponent", [0.4, 0.7, 1.6, 3.0])
def test_moment_quadrature_converges_with_nodes(exponent):
    coarse = half_taper_integrals(exponent, nodes=100)
    fine = half_taper_integrals(exponent, nodes=400)
    assert coarse[0] == pytest.approx(fine[0], rel=1e-12)
    assert coarse[1] == pytest.approx(fine[1], rel=1e-12)


def test_wetted_surface_quadrature_is_not_the_error_term():
    for exponents in ((0.7, 1.6), (MIN_EXPONENT, MAX_EXPONENT), (1.0, 1.0)):
        parameters = AsymmetricHullParameters(
            fore_exponent=exponents[0], aft_exponent=exponents[1]
        )
        assert wetted_surface_quadrature_drift(parameters) < 1e-9


# --- the reference against the shipped symmetric one ----------------------------


def test_symmetric_limit_matches_the_shipped_analytic_hydrostatics():
    parameters = AsymmetricHullParameters(fore_exponent=0.7, aft_exponent=0.7)
    reference = analytic_hydrostatics(parameters)
    shipped = hull.analytic_hydrostatics(parameters.length, parameters.beam, parameters.draft)
    for name in COMPARABLE_QUANTITIES:
        if name == "wetted_surface_area":
            # Two different quadrature rules (k = 3/p versus the shipped k = 3).
            assert float(getattr(reference, name)) == pytest.approx(
                float(getattr(shipped, name)), rel=1e-9
            )
        else:
            assert float(getattr(reference, name)) == pytest.approx(
                float(getattr(shipped, name)), rel=1e-12
            ), name


def test_symmetric_limit_mesh_is_the_hull_mesh():
    parameters = AsymmetricHullParameters(fore_exponent=0.7, aft_exponent=0.7, nx=41, nz=13)
    mesh, _ = build_mesh(parameters)
    shipped_mesh, _ = hull.build_mesh(hull.HullParameters(10.0, 2.0, 0.8, 41, 13))
    assert len(mesh.vertices) == len(shipped_mesh.vertices)
    assert len(mesh.faces) == len(shipped_mesh.faces)
    worst = max(
        abs(a - b)
        for vertex, shipped_vertex in zip(mesh.vertices, shipped_mesh.vertices, strict=False)
        for a, b in zip(vertex, shipped_vertex, strict=False)
    )
    assert worst < 1e-14
    ours = compute_hydrostatics(mesh, 0.0, closed=True)
    theirs = compute_hydrostatics(shipped_mesh, 0.0, closed=True)
    assert ours.displacement_volume == pytest.approx(theirs.displacement_volume, rel=1e-14)
    assert ours.wetted_surface_area == pytest.approx(theirs.wetted_surface_area, rel=1e-14)


# --- properties only an asymmetric family has ------------------------------------


def test_mirror_identity_puts_the_lcb_back_on_the_other_side_of_midship():
    forward = AsymmetricHullParameters(fore_exponent=0.7, aft_exponent=1.6)
    mirrored = AsymmetricHullParameters(fore_exponent=1.6, aft_exponent=0.7)
    lcb = analytic_longitudinal_centre(forward)
    mirrored_lcb = analytic_longitudinal_centre(mirrored)
    assert lcb + mirrored_lcb == pytest.approx(forward.length, rel=1e-12)
    # Everything that does not depend on where the volume sits must be unchanged.
    assert analytic_volume(forward) == pytest.approx(analytic_volume(mirrored), rel=1e-14)
    assert analytic_waterplane_area(forward) == pytest.approx(
        analytic_waterplane_area(mirrored), rel=1e-14
    )
    assert analytic_wetted_surface_area(forward) == pytest.approx(
        analytic_wetted_surface_area(mirrored), rel=1e-9
    )


def test_lcb_is_no_longer_trivial():
    lcb = analytic_longitudinal_centre(DEFAULT)
    assert abs(lcb - 0.5 * DEFAULT.length) > 0.02 * DEFAULT.length


def test_lcb_moves_towards_the_fuller_half():
    # A smaller exponent is a fuller taper, so the half with the smaller exponent
    # carries more volume and pulls the centre of buoyancy towards itself.
    fuller_fore = AsymmetricHullParameters(fore_exponent=0.7, aft_exponent=1.6)
    fuller_aft = AsymmetricHullParameters(fore_exponent=1.6, aft_exponent=0.7)
    assert analytic_longitudinal_centre(fuller_fore) > 0.5 * fuller_fore.length
    assert analytic_longitudinal_centre(fuller_aft) < 0.5 * fuller_aft.length


def test_lcf_equals_lcb_because_the_family_shares_one_taper():
    # Not an independent check: volume and waterplane are weighted by the same
    # f(x), so the two centroids coincide. Pinned so a future family that breaks
    # that (a waterline-dependent taper) cannot silently keep the claim.
    reference = analytic_hydrostatics(DEFAULT)
    assert reference.longitudinal_centre_flotation == reference.longitudinal_centre_buoyancy


def test_coefficients_follow_from_the_same_taper():
    reference = analytic_hydrostatics(DEFAULT)
    assert reference.maximum_section_area == pytest.approx(
        (2.0 / 3.0) * DEFAULT.beam * DEFAULT.draft, rel=1e-15
    )
    assert reference.vertical_centre_buoyancy == pytest.approx(
        -(2.0 / 5.0) * DEFAULT.draft, rel=1e-15
    )
    assert reference.prismatic_coefficient == pytest.approx(
        1.5 * reference.block_coefficient, rel=1e-14
    )


# --- the mesh path against the reference ----------------------------------------


@pytest.mark.parametrize("nx,nz", GRIDS)
def test_every_quantity_is_within_three_percent(nx, nz):
    model = generate_hull(AsymmetricHullParameters(nx=nx, nz=nz))
    errors = model.errors()
    assert set(errors) == set(COMPARABLE_QUANTITIES)
    for name, error in errors.items():
        assert abs(error) < 0.03, f"{name} at {nx}x{nz}: {100 * error:+.3f} %"


def test_errors_shrink_with_refinement():
    volume = []
    wetted = []
    lcb = []
    for nx, nz in GRIDS:
        model = generate_hull(AsymmetricHullParameters(nx=nx, nz=nz))
        errors = model.errors()
        volume.append(abs(errors["displacement_volume"]))
        wetted.append(abs(errors["wetted_surface_area"]))
        lcb.append(abs(errors["longitudinal_centre_buoyancy"]))
    for series in (volume, wetted, lcb):
        assert series[0] > series[1] > series[2]


def test_lcb_error_keeps_its_sign_across_refinement():
    # The sampled surface loses volume at the tips, which sit symmetrically about
    # midship only in the symmetric family; here the sign of the LCB error is a
    # real prediction, so it must not flip with resolution.
    signs = set()
    for nx, nz in GRIDS:
        model = generate_hull(AsymmetricHullParameters(nx=nx, nz=nz))
        signs.add(model.errors()["longitudinal_centre_buoyancy"] < 0.0)
    assert len(signs) == 1


def test_generated_mesh_is_a_closed_manifold():
    model = generate_hull(AsymmetricHullParameters(nx=81, nz=21))
    assert model.validation.errors == []
    geometry = model.mesh
    assert geometry.units == "m"
    report = {check.name: check for check in hull.inspect_geometry(geometry)}
    assert report["watertight"].passed, report.keys()
    assert report["euler_characteristic_genus_zero"].passed, report.keys()
    assert report["no_degenerate_faces"].passed, report.keys()


def test_even_station_count_warns_about_the_unreachable_beam():
    model = generate_hull(AsymmetricHullParameters(nx=80, nz=21))
    assert model.validation.errors == []
    assert any("beam" in warning for warning in model.validation.warnings)


def test_generation_is_deterministic_and_carries_its_provenance():
    first = generate_hull(AsymmetricHullParameters(nx=41, nz=13))
    second = generate_hull(AsymmetricHullParameters(nx=41, nz=13))
    assert first.mesh.vertices == second.mesh.vertices
    assert first.mesh.faces == second.mesh.faces
    metadata = first.metadata()
    assert metadata["schema"] == SCHEMA
    assert metadata["parameters"]["fore_exponent"] == DEFAULT.fore_exponent
    assert "analytic_reference" in metadata
    assert "relative_errors" in metadata
    assert "fore_exponent" in metadata["generator"]["longitudinal_taper"]


def test_analytic_relative_errors_report_is_self_contained():
    report = analytic_relative_errors(AsymmetricHullParameters(nx=41, nz=13))
    assert report["worst_quantity"] in COMPARABLE_QUANTITIES
    assert abs(report["worst_relative_error"]) < 0.03
    assert report["quadrature"]["wetted_surface_drift"] < 1e-9
    assert report["parameters"]["aft_exponent"] == DEFAULT.aft_exponent


# --- inputs ---------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs",
    [
        {"length": 0.0},
        {"beam": -1.0},
        {"draft": float("nan")},
        {"water_density": 0.0},
        {"fore_exponent": MIN_EXPONENT - 0.01},
        {"aft_exponent": MAX_EXPONENT + 0.01},
        {"fore_exponent": float("inf")},
        {"nx": 2},
        {"nz": 2},
        {"nx": 3.5},
        {"nx": 5000},
        {"nz": 600},
    ],
)
def test_bad_parameters_are_rejected(kwargs):
    parameters = AsymmetricHullParameters(**kwargs)
    with pytest.raises(ValueError):
        parameters.validate()
    with pytest.raises(ValueError):
        generate_hull(parameters)


def test_taper_switch_is_continuous_at_midship():
    parameters = AsymmetricHullParameters(fore_exponent=0.5, aft_exponent=2.8)
    assert taper(0.5, parameters) == pytest.approx(1.0, abs=1e-15)
    assert taper(0.0, parameters) == 0.0
    assert taper(1.0, parameters) == 0.0
    # Just either side of midship the surface must be continuous, not stepped.
    left = taper(0.5 - 1e-9, parameters)
    right = taper(0.5 + 1e-9, parameters)
    assert abs(left - right) < 1e-6
