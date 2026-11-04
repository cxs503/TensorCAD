"""Tests for the transom-stern hull family.

Three independent things are pinned here:

1. the incomplete taper integrals (series) against an independent closed form and
   an independent quadrature rule -- if the reference is wrong every number the
   family reports is wrong, so it is checked first and hardest;
2. the topology of the truncated body: watertight, edge-manifold, Euler 2, and
   the posteriori consistency with the pointed hull at ``tau = 0``;
3. the mesh integration of volume / waterplane / LCB against those references and
   its convergence under refinement, plus the design-side inverse mapping.
"""

import itertools
import math

import pytest

from tensorcad import hull as hull_module
from tensorcad import mesh as meshlib
from tensorcad import transom
from tensorcad.hydrostatics_analytic import tanh_sinh_quadrature, taper_integral

DEFAULT = transom.TransomHullParameters()


def _euler(mesh):
    return len(mesh.vertices) - len(meshlib.edge_face_counts(mesh)) + len(mesh.faces)


# --------------------------------------------------------------------------- #
# 1. The analytic references themselves
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("exponent", [0.25, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0])
def test_series_reproduces_the_closed_form_at_full_length(exponent):
    """tau = 0 must return the gamma-function closed form to round-off.

    This is the check that makes the series usable as a reference: the complete
    integral is known independently, so an error in the coefficient generation
    cannot hide behind agreement with itself.
    """
    series = transom.truncated_taper_integral(0.0, exponent)
    exact = taper_integral(exponent)
    assert series == pytest.approx(exact, rel=1e-14)


@pytest.mark.parametrize("exponent", [0.3, 0.7, 1.0, 2.5])
@pytest.mark.parametrize("tau", [0.02, 0.1, 0.2, 0.35, 0.45])
def test_series_agrees_with_tanh_sinh_quadrature(exponent, tau):
    """A second, method-independent evaluation of the same incomplete integral.

    tanh-sinh is the rule this repository already uses for the integrable
    endpoint singularities of the hull surface, so it is the natural independent
    judge for the truncated range.
    """
    series = transom.truncated_taper_integral(tau, exponent)
    quadrature = tanh_sinh_quadrature(
        lambda x: math.sin(math.pi * x) ** exponent, tau, 1.0, levels=7
    )
    assert series == pytest.approx(quadrature, rel=1e-13)


def test_series_at_the_default_exponent_matches_the_hull_constant():
    assert transom.truncated_taper_integral(0.0, hull_module.BREADTH_EXPONENT) == pytest.approx(
        hull_module.TAPER_INTEGRAL, rel=1e-14
    )


@pytest.mark.parametrize("exponent", [0.5, 0.7, 1.0, 2.0])
def test_first_moment_has_the_symmetry_implied_lcb(exponent):
    """``M(0) / J(0) = 1/2``: the pointed body's LCB is midship, exactly."""
    assert transom.truncated_taper_moment(0.0, exponent) / transom.truncated_taper_integral(
        0.0, exponent
    ) == pytest.approx(0.5, rel=1e-14)


@pytest.mark.parametrize("exponent", [0.7, 1.0, 2.0])
@pytest.mark.parametrize("tau", [0.0, 0.1, 0.3])
def test_moment_agrees_with_tanh_sinh_quadrature(exponent, tau):
    quadrature = tanh_sinh_quadrature(
        lambda x: x * math.sin(math.pi * x) ** exponent, tau, 1.0, levels=7
    )
    assert transom.truncated_taper_moment(tau, exponent) == pytest.approx(quadrature, rel=1e-13)


def test_series_coefficients_are_the_sinc_power_expansion():
    """The coefficients are pinned by their own definition, not by the integral."""
    coefficients = transom.sinc_power_coefficients(1.0)
    assert coefficients[0] == pytest.approx(1.0)
    # (sinc(pi*sqrt(w)))**1 has coefficients (-1)**k * pi**(2k) / (2k+1)!
    for k, value in enumerate(coefficients[:6]):
        expected = (-1.0) ** k * math.pi ** (2 * k) / math.factorial(2 * k + 1)
        assert value == pytest.approx(expected, rel=1e-13)


def test_taper_integral_decreases_with_the_exponent():
    values = [transom.truncated_taper_integral(0.0, e) for e in (0.3, 0.7, 1.3, 2.0)]
    assert all(a > b for a, b in itertools.pairwise(values))


# --------------------------------------------------------------------------- #
# 2. Topology of the truncated body
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("tau", [0.0, 0.01, 0.05, 0.2, 0.35, 0.49])
def test_transom_geometry_is_a_closed_oriented_manifold(tau):
    model = transom.generate_transom_hull(transom.TransomHullParameters(transom_fraction=tau))
    assert meshlib.is_watertight(model.mesh)
    assert meshlib.is_edge_manifold(model.mesh)
    assert meshlib.zero_area_faces(model.mesh) == ()
    assert _euler(model.mesh) == 2
    assert meshlib.total_signed_volume(model.mesh) > 0.0


def test_tau_zero_reproduces_the_pointed_hull():
    """The family is a generalisation: at tau = 0 nothing about the hull changes."""
    reference = hull_module.generate_hull(hull_module.HullParameters())
    mine = transom.generate_transom_hull(transom.TransomHullParameters(transom_fraction=0.0))
    assert mine.mesh.faces == reference.mesh.faces
    assert len(mine.mesh.vertices) == len(reference.mesh.vertices)
    for got, want in zip(mine.mesh.vertices, reference.mesh.vertices, strict=False):
        assert got == pytest.approx(want, abs=1e-12)
    assert mine.hydrostatics.displacement_volume == pytest.approx(
        reference.hydrostatics.displacement_volume, rel=1e-12
    )
    assert mine.hydrostatics.waterplane_area == pytest.approx(
        reference.hydrostatics.waterplane_area, rel=1e-12
    )


def test_transom_adds_one_planar_fan_without_new_vertices():
    """The plate reuses the aft ring, so it adds faces but no geometry.

    At ``tau = 0`` the aft ring sits on the centreline: its starboard and port
    halves are the same points (``nz - 1`` vertices merged by the builder) and the
    first waterline cap collapses from a quad to a triangle. The transom body has
    both back, plus the ``ring - 2`` fan triangles.
    """
    pointed = transom.generate_transom_hull(transom.TransomHullParameters(transom_fraction=0.0))
    truncated = transom.generate_transom_hull(transom.TransomHullParameters(transom_fraction=0.2))
    ring = 2 * DEFAULT.nz - 1
    assert len(truncated.mesh.vertices) - len(pointed.mesh.vertices) == DEFAULT.nz - 1
    assert len(truncated.mesh.faces) - len(pointed.mesh.faces) == (ring - 2) + 1
    # Every transom triangle lies in the plane x = 0.
    plate = [
        face
        for face in truncated.mesh.faces
        if all(abs(truncated.mesh.vertices[index][0]) <= 1e-12 for index in face)
    ]
    assert len(plate) == ring - 2


def test_hull_spans_the_requested_overall_length_with_the_transom_at_zero():
    model = transom.generate_transom_hull(transom.TransomHullParameters(transom_fraction=0.3))
    box = model.mesh.bounding_box()
    assert box.minimum[0] == pytest.approx(0.0, abs=1e-12)
    assert box.maximum[0] == pytest.approx(DEFAULT.length, rel=1e-12)
    assert box.maximum[2] == pytest.approx(0.0, abs=1e-12)
    assert box.minimum[2] == pytest.approx(-DEFAULT.draft, rel=1e-12)


def test_transom_plate_area_converges_to_the_analytic_section():
    """The fan under-reads the parabolic section, and refining closes the gap."""
    errors = []
    for nz in (13, 41, 81):
        parameters = transom.TransomHullParameters(nz=nz, nx=2 * nz - 1, transom_fraction=0.2)
        model = transom.generate_transom_hull(parameters)
        plate = [
            face
            for face in model.mesh.faces
            if all(abs(model.mesh.vertices[index][0]) <= 1e-12 for index in face)
        ]
        area = sum(meshlib.triangle_area(model.mesh, face) for face in plate)
        errors.append(plate and (area / model.transom_immersed_area_analytic - 1.0))
    assert all(error < 0.0 for error in errors), errors
    assert abs(errors[0]) > abs(errors[1]) > abs(errors[2])
    assert abs(errors[-1]) < 1e-3


# --------------------------------------------------------------------------- #
# 3. Hydrostatics against the references
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("tau", [0.0, 0.05, 0.2, 0.35])
def test_default_grid_integrates_volume_waterplane_and_lcb(tau):
    model = transom.generate_transom_hull(transom.TransomHullParameters(transom_fraction=tau))
    errors = model.relative_errors()
    assert abs(errors["displacement_volume"]) < 0.01
    assert abs(errors["waterplane_area"]) < 0.005
    assert abs(errors["longitudinal_centre_buoyancy"]) < 0.005
    assert abs(errors["block_coefficient"]) < 0.01
    assert model.validation.errors == []


@pytest.mark.parametrize("tau", [0.0, 0.2])
def test_refinement_halves_the_volume_error_twice(tau):
    """Uniform refinement must drive every integral towards the reference."""
    factors = []
    errors = []
    for nx, nz in ((41, 13), (81, 25), (161, 41)):
        model = transom.generate_transom_hull(
            transom.TransomHullParameters(nx=nx, nz=nz, transom_fraction=tau)
        )
        error = model.relative_errors()["displacement_volume"]
        errors.append(abs(error))
        factors.append(model)
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < 0.002


def test_analytic_block_coefficient_follows_its_algebra():
    """``Cb = (2/3) * J / (1 - tau)`` over the LOA, recomputed from the integrals.

    The cut removes a thin aft slice (little volume) and shortens the LOA (the
    denominator), so ``Cb`` grows -- but not without bound: past ``tau ~ 0.2`` the
    removed slice starts to matter and the trend flattens, which is why the test
    pins the algebra and the bounds rather than a monotonicity claim.
    """
    values = {}
    for tau in (0.0, 0.05, 0.2, 0.3):
        parameters = transom.TransomHullParameters(transom_fraction=tau)
        analytic = transom.analytic_transom_hydrostatics(parameters)
        expected = (
            (2.0 / 3.0) * transom.truncated_taper_integral(tau, parameters.breadth_exponent) / (1.0 - tau)
        )
        assert analytic["block_coefficient"] == pytest.approx(expected, rel=1e-14)
        assert analytic["block_coefficient"] == pytest.approx(
            analytic["displacement_volume"] / (parameters.length * parameters.beam * parameters.draft),
            rel=1e-14,
        )
        assert 0.0 < analytic["block_coefficient"] < 1.0
        values[tau] = analytic["block_coefficient"]
    assert values[0.0] == pytest.approx((2.0 / 3.0) * taper_integral(0.7), rel=1e-14)
    assert values[0.2] > values[0.0]
    assert values[0.3] > values[0.0]


def test_analytic_volume_matches_an_independent_station_integration():
    """Trapezoidal-free reference: integrate the analytic section areas directly."""
    parameters = transom.TransomHullParameters(nx=41, nz=13, transom_fraction=0.2)
    analytic = transom.analytic_transom_hydrostatics(parameters)
    reference_length = transom.parent_length(parameters)
    stations = 40001
    step = parameters.length / (stations - 1)
    total = 0.0
    for index in range(stations):
        x = index * step
        u = parameters.transom_fraction + x / reference_length
        taper = math.sin(math.pi * u) ** parameters.breadth_exponent if u < 1.0 else 0.0
        area = (2.0 / 3.0) * parameters.beam * parameters.draft * taper
        weight = 0.5 if index in (0, stations - 1) else 1.0
        total += weight * area * step
    assert total == pytest.approx(analytic["displacement_volume"], rel=2e-5)


def test_vertical_centre_of_buoyancy_is_transom_independent():
    for tau in (0.0, 0.2, 0.4):
        analytic = transom.analytic_transom_hydrostatics(
            transom.TransomHullParameters(transom_fraction=tau)
        )
        assert analytic["vertical_centre_buoyancy"] == pytest.approx(-0.4 * DEFAULT.draft)


def test_discretisation_estimate_reports_both_refinement_and_analytic_errors():
    report = transom.estimate_transom_discretisation_error(
        transom.TransomHullParameters(transom_fraction=0.2)
    )
    assert report["status"] == "measured"
    for quantity in ("displacement_volume", "waterplane_area", "longitudinal_centre_buoyancy"):
        entry = report[quantity]
        assert abs(entry["refined_error_vs_analytic"]) < abs(entry["coarse_error_vs_analytic"])


# --------------------------------------------------------------------------- #
# 4. Design-side inverse and parameter validation
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("ratio", [0.01, 0.05, 0.1, 0.2, 0.5])
@pytest.mark.parametrize("exponent", [0.7, 1.5, 3.0])
def test_area_ratio_inverse_round_trips(ratio, exponent):
    tau = transom.transom_fraction_for_area_ratio(ratio, exponent)
    assert transom.transom_area_ratio(tau, exponent) == pytest.approx(ratio, rel=1e-13)


def test_area_ratio_is_the_plate_over_the_midship_section():
    parameters = transom.TransomHullParameters(transom_fraction=0.2, nz=81, nx=161)
    model = transom.generate_transom_hull(parameters)
    assert model.transom_immersed_area_analytic / model.analytic["maximum_section_area"] == (
        pytest.approx(transom.transom_area_ratio(0.2, parameters.breadth_exponent), rel=1e-14)
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"transom_fraction": -0.01},
        {"transom_fraction": 0.5},
        {"transom_fraction": 0.8},
        {"transom_fraction": float("nan")},
        {"breadth_exponent": 0.0},
        {"breadth_exponent": 5.0},
        {"length": 0.0},
        {"beam": -1.0},
        {"draft": float("inf")},
        {"nz": 2},
        {"water_density": 0.0},
    ],
)
def test_invalid_parameters_are_rejected_with_named_errors(kwargs):
    parameters = transom.TransomHullParameters(**kwargs)
    with pytest.raises(ValueError) as error:
        parameters.validate()
    assert str(error.value)


def test_tau_beyond_midship_is_rejected_by_the_generator():
    with pytest.raises(ValueError):
        transom.generate_transom_hull(transom.TransomHullParameters(transom_fraction=0.6))


def test_incomplete_integrals_reject_out_of_range_tau():
    with pytest.raises(ValueError):
        transom.truncated_taper_integral(0.5, 0.7)
    with pytest.raises(ValueError):
        transom.truncated_taper_integral(-0.1, 0.7)
    with pytest.raises(ValueError):
        transom.transom_fraction_for_area_ratio(1.0, 0.7)


# --------------------------------------------------------------------------- #
# 5. Metadata contract
# --------------------------------------------------------------------------- #
def test_metadata_carries_parameters_reference_and_relative_errors():
    model = transom.generate_transom_hull(transom.TransomHullParameters(transom_fraction=0.2))
    metadata = model.metadata()
    assert metadata["parameters"]["transom_fraction"] == 0.2
    assert metadata["transom"]["station_x_m"] == 0.0
    assert metadata["transom"]["immersed_area_analytic_m2"] > 0.0
    assert metadata["analytic_reference"]["displacement_volume"] > 0.0
    assert abs(metadata["relative_error"]["displacement_volume"]) < 0.01
    assert metadata["hydrostatics"]["displacement_volume_m3"] > 0.0


def test_generated_surface_is_deterministic():
    first = transom.generate_transom_hull(transom.TransomHullParameters(transom_fraction=0.2))
    second = transom.generate_transom_hull(transom.TransomHullParameters(transom_fraction=0.2))
    assert first.mesh.vertices == second.mesh.vertices
    assert first.mesh.faces == second.mesh.faces
