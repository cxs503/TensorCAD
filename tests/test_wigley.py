"""The Wigley parabolic hull: closed forms, mesh agreement, and end fairness.

The tests are deliberately split by what they can prove:

* the **closed forms** are checked against independent Simpson quadratures of the
  definition, so a transcription slip in the algebra cannot pass (this is how the
  ``VCB`` coefficient was caught -- an earlier draft said ``-9/16 T`` where the
  integral gives ``-3/8 T``, a 33 % error the mesh exposed immediately);
* the **mesh** is checked against the closed forms, with the *sign* of the error
  pinned: an inscribed loft must come out low, so a sign flip means a defect;
* the **fairness** of the ends is checked on the analytic surface, against the
  sine family's divergence, because a coarse mesh cannot resolve it.
"""
from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path

import pytest

from tensorcad import mesh as meshlib
from tensorcad.curvature import assess_fairness
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.wigley import (
    ANALYTIC_COMPARABLE_QUANTITIES,
    WigleyParameters,
    analytic_hydrostatics,
    analytic_longitudinal_centre,
    analytic_maximum_section_area,
    analytic_relative_errors,
    analytic_vertical_centre_buoyancy,
    analytic_volume,
    analytic_waterplane_area,
    analytic_wetted_surface_area,
    build_mesh,
    check_consistency,
    estimate_discretisation_error,
    generate_wigley,
    half_breadth,
    waterline_curvature_at_ends,
    waterline_half_breadth,
)

CASE = WigleyParameters(length=10.0, beam=1.0, draft=0.625, nx=41, nz=13)


def _simpson(f, a: float, b: float, n: int = 4000) -> float:
    """Composite Simpson rule on ``[a, b]`` with ``n`` (even) panels."""
    assert n % 2 == 0
    h = (b - a) / n
    total = f(a) + f(b)
    for i in range(1, n):
        total += (4.0 if i % 2 else 2.0) * f(a + i * h)
    return total * h / 3.0


# --------------------------------------------------------------------------
# The closed forms, against independent quadrature of the definition
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "length,beam,draft",
    [(10.0, 1.0, 0.625), (4.356, 0.4356, 0.27225), (7.5, 2.0, 1.0)],
)
def test_volume_matches_independent_quadrature(length, beam, draft) -> None:
    p = WigleyParameters(length=length, beam=beam, draft=draft)

    def integrand(x: float) -> float:
        return 2.0 * _simpson(lambda z: half_breadth(p, x, z), -draft, 0.0, 200)

    reference = _simpson(integrand, 0.0, length, 400)
    # Simpson on a quartic integrand is exact up to round-off; a 1e-10 tolerance
    # therefore tests the algebra, not the quadrature.
    assert reference == pytest.approx(analytic_volume(length, beam, draft), rel=1e-10)
    assert analytic_volume(length, beam, draft) == pytest.approx(4.0 * length * beam * draft / 9.0)


def test_block_coefficient_is_four_ninths_for_every_hull() -> None:
    for params in (
        WigleyParameters(),
        WigleyParameters(length=4.0, beam=0.9, draft=0.3, nx=11, nz=7),
        WigleyParameters(length=100.0, beam=20.0, draft=8.0, nx=21, nz=9),
    ):
        reference = analytic_hydrostatics(params.length, params.beam, params.draft)
        assert reference.block_coefficient == pytest.approx(4.0 / 9.0, rel=1e-15)
        assert reference.prismatic_coefficient == pytest.approx(2.0 / 3.0, rel=1e-15)
        assert params.block_coefficient == pytest.approx(4.0 / 9.0, rel=0.0)


def test_waterplane_and_section_areas_match_quadrature() -> None:
    L, B, T = CASE.length, CASE.beam, CASE.draft
    waterplane = 2.0 * _simpson(lambda x: half_breadth(CASE, x, 0.0), 0.0, L, 2000)
    assert waterplane == pytest.approx(analytic_waterplane_area(L, B), rel=1e-12)
    section = 2.0 * _simpson(lambda z: half_breadth(CASE, 0.5 * L, z), -T, 0.0, 2000)
    assert section == pytest.approx(analytic_maximum_section_area(B, T), rel=1e-12)


def test_vertical_centre_buoyancy_matches_quadrature() -> None:
    """The quantity an earlier draft got wrong by 33 %: ``-3/8 T``, not ``-9/16 T``."""
    T = CASE.draft
    profile = _simpson(lambda z: 1.0 - (z / T) ** 2, -T, 0.0, 20000)
    moment = _simpson(lambda z: z * (1.0 - (z / T) ** 2), -T, 0.0, 20000)
    expected = moment / profile
    assert analytic_vertical_centre_buoyancy(T) == pytest.approx(expected, rel=1e-9)
    assert analytic_vertical_centre_buoyancy(T) == pytest.approx(-3.0 * T / 8.0)
    # and it is a different number from the sine family's -2/5 T, i.e. the two
    # families are genuinely different sections rather than the same one relabelled
    assert abs(analytic_vertical_centre_buoyancy(T) - (-2.0 * T / 5.0)) > 0.01 * T


def test_longitudinal_centres_are_midship() -> None:
    assert analytic_longitudinal_centre(CASE.length) == pytest.approx(0.5 * CASE.length)


def test_wetted_area_quadrature_is_converged() -> None:
    """The reference must not be the error term: doubling nodes changes <1e-9."""
    coarse = analytic_wetted_surface_area(CASE.length, CASE.beam, CASE.draft, nodes=40)
    fine = analytic_wetted_surface_area(CASE.length, CASE.beam, CASE.draft, nodes=80)
    finer = analytic_wetted_surface_area(CASE.length, CASE.beam, CASE.draft, nodes=160)
    assert abs(fine - coarse) / fine < 1e-8
    assert abs(finer - fine) / finer < 1e-10


def test_wetted_area_matches_independent_simpson_product_rule() -> None:
    L, B, T = CASE.length, CASE.beam, CASE.draft
    hb = 0.5 * B

    def speed(x: float, z: float) -> float:
        u = 2.0 * x / L - 1.0
        v = z / T
        b_x = hb * (-4.0 * u / L) * (1.0 - v * v)
        b_z = hb * (1.0 - u * u) * (-2.0 * z / (T * T))
        return math.sqrt(1.0 + b_x * b_x + b_z * b_z)

    reference = 2.0 * _simpson(
        lambda x: _simpson(lambda z: speed(x, z), -T, 0.0, 200), 0.0, L, 400
    )
    assert analytic_wetted_surface_area(L, B, T) == pytest.approx(reference, rel=1e-7)


def test_hydrostatics_reference_states_its_basis() -> None:
    reference = analytic_hydrostatics(CASE.length, CASE.beam, CASE.draft, 1000.0)
    assert "closed forms" in reference.basis
    assert "quadrature" in reference.basis
    assert reference.displacement_mass == pytest.approx(
        reference.displacement_volume * 1000.0
    )
    assert reference.waterline_z == 0.0


# --------------------------------------------------------------------------
# The mesh, against the closed forms
# --------------------------------------------------------------------------


def test_mesh_is_a_closed_oriented_solid() -> None:
    model = generate_wigley(CASE)
    mesh = model.mesh
    assert meshlib.is_watertight(mesh)
    assert meshlib.is_edge_manifold(mesh)
    assert meshlib.zero_area_faces(mesh) == ()
    assert meshlib.total_signed_volume(mesh) > 0.0
    assert not model.validation.errors
    # V - E + F = 2 for a genus-zero closed surface
    edges = meshlib.edge_face_counts(mesh)
    assert len(mesh.vertices) - len(edges) + len(mesh.faces) == 2


def test_mesh_honours_the_requested_dimensions() -> None:
    model = generate_wigley(CASE)
    box = model.mesh.bounding_box()
    assert box.minimum[0] == pytest.approx(0.0)
    assert box.maximum[0] == pytest.approx(CASE.length)
    assert box.maximum[1] - box.minimum[1] == pytest.approx(CASE.beam, rel=1e-12)
    assert box.maximum[2] == pytest.approx(0.0)
    assert box.minimum[2] == pytest.approx(-CASE.draft)


def test_mesh_agrees_with_the_closed_forms_and_is_low_everywhere() -> None:
    model = generate_wigley(CASE)
    reference = analytic_hydrostatics(CASE.length, CASE.beam, CASE.draft)
    for name in ANALYTIC_COMPARABLE_QUANTITIES:
        got = getattr(model.hydrostatics, name)
        want = getattr(reference, name)
        error = (got - want) / want
        assert abs(error) < 3e-3, f"{name}: {100 * error:+.4f}%"
        if name in ("longitudinal_centre_buoyancy", "longitudinal_centre_flotation"):
            assert error == pytest.approx(0.0, abs=1e-9)
        elif name in ("displacement_volume", "block_coefficient", "waterplane_area",
                      "prismatic_coefficient", "maximum_section_area", "wetted_surface_area"):
            # inscribed loft: the mesh cannot exceed the analytic surface
            assert error <= 0.0, f"{name} came out high ({100 * error:+.4f}%)"


def test_error_falls_with_the_second_order_and_refines_towards_the_reference() -> None:
    quantities = estimate_discretisation_error(CASE)["quantities"]
    assert quantities["coarse_grid" if False else "displacement_volume"]["coarse_error"] < 0.0
    for name in ("displacement_volume", "wetted_surface_area", "waterplane_area"):
        entry = quantities[name]
        assert abs(entry["refined_error"]) < abs(entry["coarse_error"])
        assert entry["relative_change"] > 0.0
    # refinement doubles both counts; a second-order quantity should improve ~4x
    ratio = abs(quantities["waterplane_area"]["coarse_error"]) / abs(
        quantities["waterplane_area"]["refined_error"]
    )
    assert 3.0 < ratio < 5.0, ratio


def test_analytic_relative_errors_are_signed_and_complete() -> None:
    errors = analytic_relative_errors(CASE)
    assert set(errors) == set(ANALYTIC_COMPARABLE_QUANTITIES)
    for name, entry in errors.items():
        assert entry["relative_error"] is not None, name
        assert entry["mesh"] == pytest.approx(
            getattr(generate_wigley(CASE).hydrostatics, name)
        )


def test_offsets_report_the_full_breadth() -> None:
    model = generate_wigley(CASE)
    for x, breadth in model.offsets:
        expected = 2.0 * waterline_half_breadth(CASE, min(max(x, 0.0), CASE.length))
        assert breadth == pytest.approx(expected, abs=1e-12)


def test_generation_is_deterministic() -> None:
    a = generate_wigley(CASE)
    b = generate_wigley(CASE)
    assert tuple(a.mesh.vertices) == tuple(b.mesh.vertices)
    assert tuple(a.mesh.faces) == tuple(b.mesh.faces)


def test_even_station_count_warns_that_beam_is_under_realised() -> None:
    params = WigleyParameters(length=10.0, beam=1.0, draft=0.625, nx=40, nz=13)
    model = generate_wigley(params)
    assert any(not c.passed and c.name == "midship_station" for c in model.validation.checks)
    assert model.validation.warnings


def test_check_consistency_passes_on_the_default_grid() -> None:
    mesh, _ = build_mesh(CASE)
    hydro = compute_hydrostatics(mesh, 0.0, water_density=CASE.water_density, closed=True)
    checks = check_consistency(CASE, mesh, hydro)
    assert all(c.passed for c in checks)
    assert {c.name for c in checks} == {"sampling_volume", "block_coefficient", "midship_station"}


def test_starved_grid_fails_the_sampling_check() -> None:
    params = WigleyParameters(length=10.0, beam=1.0, draft=0.625, nx=5, nz=3)
    model = generate_wigley(params)
    assert any(not c.passed and c.name == "sampling_volume" for c in model.validation.checks)


# --------------------------------------------------------------------------
# Fairness: the reason to carry a second family at all
# --------------------------------------------------------------------------


def test_waterline_has_a_finite_slope_at_both_ends() -> None:
    """``db/dx`` at the stem and stern is exactly ``2B/L`` -- finite, not zero.

    The parabola ``1 - (2x/L - 1)^2`` is steepest at the ends and flat at midship,
    which is the opposite of what a fairing intuition expects but exactly what
    makes the end an ordinary point: a finite slope and a finite curvature.
    """
    L, B = CASE.length, CASE.beam
    for x, sign in ((0.0, +1.0), (L, -1.0)):
        # 1e-12*L is deliberately not probed at the stem: ``1 - u^2`` with
        # ``u = 1 - 2e-12`` is a cancellation, not a geometric question.
        for depth in (1e-6, 1e-9):
            h = depth * L
            step = h if x == 0.0 else -h
            slope = (waterline_half_breadth(CASE, x + step) - waterline_half_breadth(CASE, x)) / step
            # one-sided secant; the analytic limit is 2B/L in magnitude at both ends
            assert slope == pytest.approx(sign * 2.0 * B / L, rel=1e-4), (x, depth, slope)
    assert waterline_half_breadth(CASE, 0.0) == pytest.approx(0.0, abs=1e-15)
    assert waterline_half_breadth(CASE, L) == pytest.approx(0.0, abs=1e-15)
    # midship is where the waterline is flat
    h = 1e-6 * L
    midship = (
        waterline_half_breadth(CASE, 0.5 * L + h) - waterline_half_breadth(CASE, 0.5 * L - h)
    ) / (2.0 * h)
    assert midship == pytest.approx(0.0, abs=1e-9 * B / L)


def test_waterline_curvature_at_the_ends_is_finite_and_equals_the_closed_form() -> None:
    L, B = CASE.length, CASE.beam
    h = 1e-4 * L
    for x in (0.01 * L, 0.5 * L, 0.99 * L):
        second = (
            waterline_half_breadth(CASE, x + h)
            - 2.0 * waterline_half_breadth(CASE, x)
            + waterline_half_breadth(CASE, x - h)
        ) / (h * h)
        # a parabola: the central second difference is exact to round-off
        assert second == pytest.approx(-4.0 * B / (L * L), rel=1e-6)
    assert waterline_curvature_at_ends(CASE) == pytest.approx(-4.0 * B / (L * L))


def test_the_sine_family_diverges_where_the_wigley_family_does_not() -> None:
    """Same probe depth, two families: one bounded slope, one unbounded."""
    L, B = CASE.length, CASE.beam
    x = 1e-9 * L
    wigley_slope = waterline_half_breadth(CASE, x) / x

    def sine(xx: float) -> float:
        return 0.5 * B * math.sin(math.pi * xx / L) ** 0.7

    sine_slope = sine(x) / x
    # the Wigley secant is bounded by its analytic end slope 2B/L
    assert wigley_slope == pytest.approx(2.0 * B / L, rel=1e-4)
    # the sine family's secant grows as x^-0.3: four decades deeper is ~16x steeper
    assert sine_slope > 4.0 * B / L
    # three decades deeper: 10**(0.3*3) = 8x steeper, measured a bit lower
    assert sine(1e-12 * L) / (1e-12 * L) > 6.0 * sine_slope


def test_fairness_report_reads_the_wigley_mesh_as_smooth() -> None:
    model = generate_wigley(WigleyParameters(nx=81, nz=25, length=10.0, beam=1.0, draft=0.625))
    report = assess_fairness(model.mesh)
    assert report.watertight and report.edge_manifold
    assert report.boundary_edges == 0
    assert report.crease_length_fraction < 0.05


# --------------------------------------------------------------------------
# Contract: input rejection, stdlib-only import
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"length": 0.0}, "length"),
        ({"beam": -1.0}, "beam"),
        ({"draft": float("nan")}, "draft"),
        ({"nx": 2}, "nx"),
        ({"nz": 1000}, "nz"),
        ({"water_density": 0.0}, "water_density"),
        ({"nx": True}, "nx"),
    ],
)
def test_parameter_validation_rejects_bad_input(kwargs, match) -> None:
    with pytest.raises((ValueError, TypeError), match=match):
        WigleyParameters(**kwargs).validate()


def test_half_breadth_rejects_points_outside_the_body() -> None:
    with pytest.raises(ValueError, match="outside"):
        half_breadth(CASE, -0.1, -0.1)
    with pytest.raises(ValueError, match="outside"):
        half_breadth(CASE, CASE.length + 0.1, -0.1)
    with pytest.raises(ValueError, match="outside"):
        half_breadth(CASE, 1.0, 0.5)


def test_closed_forms_reject_non_positive_dimensions() -> None:
    with pytest.raises(ValueError):
        analytic_volume(0.0, 1.0, 0.5)
    with pytest.raises(ValueError):
        analytic_wetted_surface_area(10.0, -1.0, 0.5)
    with pytest.raises(ValueError):
        analytic_wetted_surface_area(10.0, 1.0, 0.5, nodes=4)


def test_model_serialises_with_the_family_named() -> None:
    data = generate_wigley(CASE).to_dict()
    assert data["family"] == "wigley_parabolic_hull"
    assert data["generator"]["half_breadth"].startswith("(B/2)*(1-(2x/L-1)^2)")
    assert data["hydrostatics"]["displacement_volume_m3"] == pytest.approx(
        analytic_volume(CASE.length, CASE.beam, CASE.draft), rel=3e-3
    )


def test_import_is_standard_library_only() -> None:
    script = """
import sys
class Blocker:
    def find_module(self, name, path=None):
        if name.split('.')[0] in ('torch', 'numpy', 'scipy'):
            return self
    def load_module(self, name):
        raise ImportError('BLOCKED: ' + name)
sys.meta_path.insert(0, Blocker())
import tensorcad.wigley as w
print(w.analytic_volume(10.0, 1.0, 0.625))
"""
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        [sys.executable, "-c", script], cwd=root, capture_output=True, text=True, check=False
    )
    assert result.returncode == 0, result.stderr
    assert float(result.stdout.strip()) == pytest.approx(analytic_volume(10.0, 1.0, 0.625))
