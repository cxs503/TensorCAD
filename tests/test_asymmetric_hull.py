"""The asymmetric fore/aft taper: closed forms, and the longitudinal centres it earns.

The symmetric family measures ``LCB = LCF = length/2``, which is an identity, not a
result: it validates the frame and the sign of the ``x`` moments and nothing else.
Giving the two halves their own taper exponents makes the longitudinal centres a
quantity again -- one with an independent reference, an exact involution
(``LCB(p_aft, p_fore) = length - LCB(p_fore, p_aft)``) and a monotonic response to
the asymmetry. That is what these tests measure.
"""
from __future__ import annotations

import itertools
import math

import pytest

from tensorcad import (
    HullParameters,
    analytic_wetted_surface_area,
    build_mesh,
    estimate_discretisation_error,
    generate_hull,
    resolve_resolution,
)
from tensorcad import mesh as meshlib
from tensorcad.hull import (
    BREADTH_EXPONENT,
    TAPER_INTEGRAL,
    _half_first_moment,
    _half_first_moment_exp_sinh,
    analytic_hydrostatics,
    analytic_longitudinal_centre,
    analytic_volume,
    analytic_waterplane_area,
    taper_integral,
    taper_integral_pair,
)

L, B, T = 10.0, 2.0, 0.8
SYMMETRIC = {"fore_exponent": BREADTH_EXPONENT, "aft_exponent": BREADTH_EXPONENT}


def _hydrostatics(fore: float, aft: float, nx: int = 81, nz: int = 25):
    parameters = HullParameters(L, B, T, nx, nz, 1025.0, fore_exponent=fore, aft_exponent=aft)
    model = generate_hull(parameters)
    assert model.hydrostatics is not None
    return parameters, model, model.hydrostatics


def _simpson(f, a: float, b: float, intervals: int) -> float:
    """Composite Simpson on an even number of intervals."""
    assert intervals % 2 == 0
    h = (b - a) / intervals
    total = f(a) + f(b)
    for i in range(1, intervals):
        total += (4.0 if i % 2 else 2.0) * f(a + i * h)
    return total * h / 3.0


# --------------------------------------------------------------------------
# The symmetric form must not move: every stored number was computed with it.
# --------------------------------------------------------------------------
def test_default_parameters_reproduce_the_symmetric_family_bit_for_bit() -> None:
    default = generate_hull()
    explicit = generate_hull(HullParameters(fore_exponent=0.7, aft_exponent=0.7))
    assert default.mesh.vertices == explicit.mesh.vertices
    assert default.mesh.faces == explicit.mesh.faces
    assert default.hydrostatics == explicit.hydrostatics
    assert HullParameters().is_symmetric_fore_and_aft is True
    assert HullParameters(fore_exponent=0.5).is_symmetric_fore_and_aft is False


def test_symmetric_references_are_the_published_constants_exactly() -> None:
    # Exact equality, not a tolerance: the taper pair averages two identical
    # integrals and the centre short-circuits, so no stored digit may move.
    assert analytic_volume(L, B, T) == B * L * T * TAPER_INTEGRAL * (2.0 / 3.0)
    assert analytic_waterplane_area(L, B) == B * L * TAPER_INTEGRAL
    assert analytic_longitudinal_centre(L) == 0.5 * L
    assert taper_integral(BREADTH_EXPONENT) is TAPER_INTEGRAL
    assert taper_integral_pair(BREADTH_EXPONENT, BREADTH_EXPONENT) == TAPER_INTEGRAL
    symmetric = analytic_hydrostatics(L, B, T)
    assert symmetric.basis == (
        "closed forms of the analytic surface (wetted area by Gauss-Legendre quadrature); "
        "longitudinal waterplane inertia by tanh-sinh quadrature"
    )
    assert (
        analytic_wetted_surface_area(L, B, T)
        == analytic_wetted_surface_area(L, B, T, fore_exponent=0.7, aft_exponent=0.7)
    )


def test_taper_integral_matches_an_independent_quadrature() -> None:
    """``taper/2`` re-integrated with Simpson, on a substitution that kills the tip.

    The taper is singular at *both* ends of ``[0, 1]`` (``sin(pi x)^p`` has an
    infinite derivative at the stem and at the post), so the check integrates one
    half where only ``x = 0`` is singular, with ``x = t^3/2``, and doubles it.
    """
    for p in (0.2, 0.4, 0.7, 1.0, 1.4, 2.0, 3.0):
        half = _simpson(
            lambda t, p=p: math.sin(0.5 * math.pi * t ** 3) ** p * 1.5 * t * t,
            0.0,
            1.0,
            200_000,
        )
        assert 2.0 * half == pytest.approx(taper_integral(p), rel=1e-11)


def test_the_half_moment_sits_where_the_two_half_distributions_say_it_does() -> None:
    """The centre of a two-piece taper is the volume-weighted mean of the pieces.

    ``LCB = (h_a*c_a + h_f*c_f) / (h_a + h_f)`` with ``h`` the half integrals and
    ``c`` the centre of each half. This checks the closed-form centre against that
    decomposition, which is the physical statement it encodes.
    """
    for fore, aft in ((0.5, 0.9), (0.4, 1.4), (1.2, 0.6)):
        h_f = 0.5 * taper_integral(fore)
        h_a = 0.5 * taper_integral(aft)
        c_a = L * _half_first_moment(aft) / h_a                   # measured from x = 0
        c_f = L * (1.0 - _half_first_moment(fore) / h_f)           # mirrored: from x = L
        combined = (h_a * c_a + h_f * c_f) / (h_a + h_f)
        assert analytic_longitudinal_centre(
            L, fore_exponent=fore, aft_exponent=aft
        ) == pytest.approx(combined, abs=1e-9)


# --------------------------------------------------------------------------
# The one quantity without a closed form: two schemes, no shared machinery.
# --------------------------------------------------------------------------
@pytest.mark.parametrize("exponent", [0.2, 0.4, 0.7, 1.0, 1.4, 2.0, 3.0])
def test_half_first_moment_agrees_across_two_independent_rules(exponent: float) -> None:
    gauss = _half_first_moment(exponent)
    exp_sinh = _half_first_moment_exp_sinh(exponent)
    assert gauss == pytest.approx(exp_sinh, rel=1e-13)
    # and the Gauss rule is itself converged in its node count
    assert gauss == pytest.approx(_half_first_moment(exponent, nodes=60), rel=1e-13)


def test_half_first_moment_is_bounded_by_its_own_definition() -> None:
    # 0 <= M(p) <= (1/2) * h(p), because 0 <= x <= 1/2 on the interval.
    for p in (0.2, 0.7, 2.0):
        half = 0.5 * taper_integral(p)
        assert 0.0 < _half_first_moment(p) < 0.5 * half


# --------------------------------------------------------------------------
# The reference itself: swap invariance, monotonicity, continuity.
# --------------------------------------------------------------------------
def test_longitudinal_centre_is_antisymmetric_under_the_swap() -> None:
    cases = [(0.5, 0.9), (0.4, 1.4), (0.7, 1.2), (0.35, 0.71)]
    for fore, aft in cases:
        one = analytic_longitudinal_centre(L, fore_exponent=fore, aft_exponent=aft)
        two = analytic_longitudinal_centre(L, fore_exponent=aft, aft_exponent=fore)
        assert one + two == pytest.approx(L, rel=0.0, abs=1e-12)
        assert one != 0.5 * L
    # and the magnitudes that do not depend on the split do not move
    for fore, aft in cases:
        assert analytic_volume(L, B, T, fore_exponent=fore, aft_exponent=aft) == pytest.approx(
            analytic_volume(L, B, T, fore_exponent=aft, aft_exponent=fore), rel=1e-15
        )
        assert analytic_waterplane_area(
            L, B, fore_exponent=fore, aft_exponent=aft
        ) == pytest.approx(
            analytic_waterplane_area(L, B, fore_exponent=aft, aft_exponent=fore), rel=1e-15
        )


def test_fuller_forebody_moves_the_centre_forward_monotonically() -> None:
    values = [
        analytic_longitudinal_centre(L, fore_exponent=p, aft_exponent=0.9)
        for p in (0.3, 0.5, 0.7, 1.0, 1.5)
    ]
    assert all(earlier > later for earlier, later in itertools.pairwise(values))
    assert values[0] > 0.5 * L > values[-1]


def test_longitudinal_centre_is_continuous_at_equal_exponents() -> None:
    epsilon = 1e-9
    near = analytic_longitudinal_centre(
        L, fore_exponent=0.7, aft_exponent=0.7 + epsilon
    )
    assert near == pytest.approx(0.5 * L, abs=1e-6)


# --------------------------------------------------------------------------
# Mesh against reference.
# --------------------------------------------------------------------------
def test_mesh_matches_the_asymmetric_closed_forms() -> None:
    _, _, measured = _hydrostatics(0.5, 0.9, nx=81, nz=25)
    reference = analytic_hydrostatics(L, B, T, 1025.0, fore_exponent=0.5, aft_exponent=0.9)
    assert measured.displacement_volume == pytest.approx(reference.displacement_volume, rel=1e-2)
    assert measured.waterplane_area == pytest.approx(reference.waterplane_area, rel=1e-2)
    assert measured.block_coefficient == pytest.approx(reference.block_coefficient, rel=1e-2)
    assert measured.prismatic_coefficient == pytest.approx(
        reference.prismatic_coefficient, rel=1e-2
    )
    # the point of the whole family: the centre is no longer length/2 by symmetry
    assert abs(reference.longitudinal_centre_buoyancy - 0.5 * L) > 0.2
    assert measured.longitudinal_centre_buoyancy == pytest.approx(
        reference.longitudinal_centre_buoyancy, abs=5e-3
    )


def test_the_mesh_centre_converges_onto_the_asymmetric_reference() -> None:
    reference = analytic_longitudinal_centre(L, fore_exponent=0.5, aft_exponent=0.9)
    errors = []
    for nx, nz in ((41, 13), (81, 25), (161, 49), (321, 97)):
        _, _, measured = _hydrostatics(0.5, 0.9, nx=nx, nz=nz)
        errors.append(measured.longitudinal_centre_buoyancy - reference)
    # same sign at every resolution (a systematic, not a random, bias) and shrinking
    assert all(error < 0.0 for error in errors)
    assert all(abs(a) > abs(b) for a, b in itertools.pairwise(errors))
    assert abs(errors[-1]) < 2e-3
    assert abs(errors[0]) > 5.0 * abs(errors[-1])


def test_the_swap_identity_holds_on_the_mesh_as_well() -> None:
    _, _, one = _hydrostatics(0.5, 0.9, nx=81, nz=25)
    _, _, two = _hydrostatics(0.9, 0.5, nx=81, nz=25)
    total = one.longitudinal_centre_buoyancy + two.longitudinal_centre_buoyancy
    assert total == pytest.approx(L, abs=1e-6)
    # the discretisation error, not the geometry, is what breaks it: it mirrors
    reference = analytic_longitudinal_centre(L, fore_exponent=0.5, aft_exponent=0.9)
    assert (one.longitudinal_centre_buoyancy - reference) == pytest.approx(
        -(two.longitudinal_centre_buoyancy - (L - reference)), abs=1e-6
    )
    assert one.displacement_volume == pytest.approx(two.displacement_volume, rel=1e-15)
    assert one.waterplane_area == pytest.approx(two.waterplane_area, rel=1e-15)


def test_lcb_and_lcf_are_the_same_number_for_an_asymmetric_hull() -> None:
    # The waterplane breadth is proportional to the same taper, so a quantity that
    # was one identity on the symmetric family becomes two independent mesh
    # integrals that have to land on the same number.
    _, _, measured = _hydrostatics(0.4, 1.4, nx=161, nz=49)
    assert measured.longitudinal_centre_buoyancy == pytest.approx(
        measured.longitudinal_centre_flotation, rel=1e-9
    )
    assert measured.longitudinal_centre_buoyancy != pytest.approx(0.5 * L, abs=1e-3)


def test_wetted_area_is_invariant_under_the_swap() -> None:
    one = analytic_wetted_surface_area(L, B, T, fore_exponent=0.5, aft_exponent=0.9)
    two = analytic_wetted_surface_area(L, B, T, fore_exponent=0.9, aft_exponent=0.5)
    assert one == pytest.approx(two, rel=1e-12)


@pytest.mark.parametrize(("nx", "nz"), [(41, 13), (60, 9), (81, 25)])
def test_an_asymmetric_mesh_is_still_a_closed_oriented_solid(nx: int, nz: int) -> None:
    parameters = HullParameters(L, B, T, nx, nz, 1025.0, fore_exponent=0.4, aft_exponent=1.4)
    mesh, offsets = build_mesh(parameters)
    assert meshlib.is_watertight(mesh)
    assert meshlib.is_edge_manifold(mesh)
    assert not meshlib.zero_area_faces(mesh)
    edges = {
        tuple(sorted((face[k], face[(k + 1) % 3])))
        for face in mesh.faces
        for k in range(3)
    }
    assert len(mesh.vertices) - len(edges) + len(mesh.faces) == 2
    assert meshlib.total_signed_volume(mesh) > 0.0
    assert offsets[0][1] == 0.0 and offsets[-1][1] == 0.0


def test_the_maximum_section_stays_at_midship_for_any_split() -> None:
    # sin(pi*x)^p peaks at x = L/2 for every exponent, so the midship section area
    # -- and therefore Cp -- is untouched by the asymmetry.
    reference = analytic_hydrostatics(L, B, T, fore_exponent=0.3, aft_exponent=2.0)
    symmetric = analytic_hydrostatics(L, B, T)
    assert reference.maximum_section_area == symmetric.maximum_section_area
    assert reference.prismatic_coefficient == pytest.approx(
        analytic_volume(L, B, T, fore_exponent=0.3, aft_exponent=2.0)
        / (L * reference.maximum_section_area),
        rel=1e-15,
    )


# --------------------------------------------------------------------------
# Plumbing: validation, metadata, refinement, CLI.
# --------------------------------------------------------------------------
@pytest.mark.parametrize("value", [0.0, -0.5, 6.0, float("nan"), float("inf")])
def test_the_exponents_are_validated(value: float) -> None:
    parameters = HullParameters(fore_exponent=value)
    errors, _ = parameters.problems()
    assert errors and any("fore_exponent" in error for error in errors)
    with pytest.raises(ValueError):
        parameters.validate()


def test_the_extremes_of_the_accepted_range_still_build() -> None:
    for fore, aft in ((0.05, 5.0), (5.0, 0.05)):
        model = generate_hull(
            HullParameters(L, B, T, 41, 13, 1025.0, fore_exponent=fore, aft_exponent=aft)
        )
        assert model.hydrostatics.displacement_volume > 0.0
        assert model.validation.passed


def test_metadata_records_the_fore_aft_split() -> None:
    model = generate_hull(
        HullParameters(L, B, T, 41, 13, 1025.0, fore_exponent=0.5, aft_exponent=0.9)
    )
    metadata = model.metadata()
    assert metadata["parameters"]["fore_exponent"] == 0.5
    assert metadata["parameters"]["aft_exponent"] == 0.9
    assert metadata["generator"]["longitudinal_taper"] == (
        "sin(pi*x/length)^0.9 (aft) / sin(pi*x/length)^0.5 (fore)"
    )
    assert metadata["generator"]["is_symmetric_fore_and_aft"] is False
    # the record's hydrostatics are mesh integrals; the reference says what it is
    reference = analytic_hydrostatics(L, B, T, 1025.0, fore_exponent=0.5, aft_exponent=0.9)
    assert "no elementary form" in reference.basis
    assert metadata["hydrostatics"]["basis"] == "closed surface integrals"


def test_discretisation_estimate_refines_the_same_asymmetric_hull() -> None:
    parameters = HullParameters(L, B, T, 21, 9, 1025.0, fore_exponent=0.5, aft_exponent=0.9)
    report = estimate_discretisation_error(parameters)
    assert report["status"] == "measured"
    assert report["refined_grid"] == {"nx": 41, "nz": 17}
    for value in report["relative_change"].values():
        assert 0.0 < value < 0.1


def test_resolution_search_uses_the_asymmetric_reference() -> None:
    parameters = HullParameters(L, B, T, 21, 9, 1025.0, fore_exponent=0.5, aft_exponent=0.9)
    plan = resolve_resolution(parameters, target_volume_error=0.005)
    assert plan.accepted
    assert abs(plan.relative_error) <= 0.005
    assert plan.parameters.fore_exponent == 0.5
    assert plan.parameters.aft_exponent == 0.9


def test_cli_carries_the_exponents_through(capsys, tmp_path) -> None:
    import json

    from tensorcad.cli import main

    code = main(
        [
            "--length", f"{L}", "--beam", f"{B}", "--draft", f"{T}",
            "--nx", "41", "--nz", "13",
            "--fore-exponent", "0.5", "--aft-exponent", "0.9",
            "--no-export", "--json",
        ]
    )
    assert code == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["parameters"]["fore_exponent"] == 0.5
    assert payload["parameters"]["aft_exponent"] == 0.9
    assert payload["hydrostatics"]["longitudinal_centre_buoyancy_m"] > 0.5 * L


def test_the_mesh_digest_of_an_asymmetric_hull_is_deterministic() -> None:
    parameters = HullParameters(L, B, T, 41, 13, 1025.0, fore_exponent=0.5, aft_exponent=0.9)
    first, _ = build_mesh(parameters)
    second, _ = build_mesh(parameters)
    assert first.vertices == second.vertices
    assert first.faces == second.faces
    # and the swap is a different surface, not the same one relabelled
    swapped, _ = build_mesh(
        HullParameters(L, B, T, 41, 13, 1025.0, fore_exponent=0.9, aft_exponent=0.5)
    )
    assert first.vertices != swapped.vertices
