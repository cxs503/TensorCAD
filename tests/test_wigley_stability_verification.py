"""Wigley initial stability: closed forms, two mesh routes, and the GM amplification.

The Wigley family was already a verification case for volume, areas and coefficients
(``docs/wigley.md``) and for the displacement side of ``GM`` there was only the curved
``tensorcad.hull`` family, whose ``I_T`` needs a Beta integral. This file pins the
*closed forms for the waterplane second moments* and measures both mesh routes of the
merged stability stack against them:

* ``I_T = (B^3 L / 12) K(n)``, ``K(n) = 1 - 3/(n+1) + 3/(2n+1) - 1/(3n+1)``;
* ``I_L = (B L^3 / 4)(1/3 - 1/(n+3))`` about the centre of flotation;
* ``BM = I/nabla``, ``KB = VCB(m)``, ``GM = KB + BM - KG``.

Three independent references are used, none of which is another route through the same
code:

1. **Hand rationals at ``n = 2``**: ``4LB^3/105``, ``BL^3/30``, ``3B^2/(35T)``,
   ``3L^2/(40T)``, written out in the test.
2. **Composite Simpson on the analytic half-breadth** (a second, unrelated quadrature),
   including a fractional exponent where Simpson is not exact.
3. **Two mesh routes that share no code**: the plane cut in ``tensorcad.stability``
   (Green's theorem on closed loops) and the wetted-surface flux identity in
   ``tensorcad.hydrostatics`` (divergence theorem). Both are exact on a polygonal
   waterplane, so they must agree to round-off; what is left is the mesh's chord error.

The claims are the ones ``examples/wigley_stability/run_example.py`` writes into
``artifacts/wigley_stability/results.json``; the thresholds are duplicated here so a
silent loosening of the example is caught.
"""

from __future__ import annotations

import math

import pytest

from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.stability import initial_stability, waterplane_moments
from tensorcad.wigley import (
    WigleyParameters,
    analytic_initial_stability,
    analytic_metacentric_radius_longitudinal,
    analytic_metacentric_radius_transverse,
    analytic_vertical_centre_buoyancy,
    analytic_volume,
    analytic_waterplane_area,
    analytic_waterplane_inertia,
    analytic_waterplane_inertia_longitudinal,
    generate_wigley,
    half_breadth,
)

LENGTH, BEAM, DRAFT = 10.0, 1.0, 0.625  # ITTC reference proportions: L/B = 10, B/T = 1.6
KG = 0.5


def reference(**overrides) -> WigleyParameters:
    values = dict(length=LENGTH, beam=BEAM, draft=DRAFT)
    values.update(overrides)
    return WigleyParameters(**values)


def simpson(values: list[float], step: float) -> float:
    """Composite Simpson on an even number of panels (written here, not imported)."""
    assert (len(values) - 1) % 2 == 0
    total = values[0] + values[-1]
    total += 4.0 * math.fsum(values[1:-1:2])
    total += 2.0 * math.fsum(values[2:-1:2])
    return total * step / 3.0


def quadrature(parameters: WigleyParameters, panels: int = 4000) -> tuple[float, float]:
    """``(I_T, I_L)`` by Simpson on ``b(x) = half_breadth(x, depth=0)``."""
    step = parameters.length / panels
    xs = [i * step for i in range(panels + 1)]
    breadths = [half_breadth(parameters, x, 0.0) for x in xs]
    transverse = simpson([(2.0 / 3.0) * b**3 for b in breadths], step)
    longitudinal = simpson(
        [2.0 * (x - 0.5 * parameters.length) ** 2 * b for x, b in zip(xs, breadths)], step
    )
    return transverse, longitudinal


# --------------------------------------------------------------------------- #
# 1. the closed forms against hand rationals and an unrelated quadrature
# --------------------------------------------------------------------------- #
def test_the_n2_closed_forms_are_the_hand_rationals():
    assert analytic_waterplane_inertia(LENGTH, BEAM, 2.0) == pytest.approx(
        4.0 * LENGTH * BEAM**3 / 105.0, rel=1e-15
    )
    assert analytic_waterplane_inertia_longitudinal(LENGTH, BEAM, 2.0) == pytest.approx(
        BEAM * LENGTH**3 / 30.0, rel=1e-15
    )
    assert analytic_metacentric_radius_transverse(LENGTH, BEAM, DRAFT, 2.0, 2.0) == pytest.approx(
        3.0 * BEAM**2 / (35.0 * DRAFT), rel=1e-15
    )
    assert analytic_metacentric_radius_longitudinal(
        LENGTH, BEAM, DRAFT, 2.0, 2.0
    ) == pytest.approx(3.0 * LENGTH**2 / (40.0 * DRAFT), rel=1e-15)


@pytest.mark.parametrize("exponent_pair", [(2.0, 2.0), (3.0, 2.0), (2.0, 3.0), (3.0, 3.0), (2.5, 2.0)])
def test_the_closed_forms_match_independent_simpson_quadrature(exponent_pair):
    n, m = exponent_pair
    parameters = reference(longitudinal_exponent=n, vertical_exponent=m)
    transverse, longitudinal = quadrature(parameters)
    assert analytic_waterplane_inertia(LENGTH, BEAM, n) == pytest.approx(transverse, rel=1e-8)
    assert analytic_waterplane_inertia_longitudinal(LENGTH, BEAM, n) == pytest.approx(
        longitudinal, rel=1e-8
    )


def test_the_vertical_exponent_moves_kb_and_nabla_but_not_the_inertias():
    area = analytic_waterplane_area(LENGTH, BEAM, 2.0)
    inertia = analytic_waterplane_inertia(LENGTH, BEAM, 2.0)
    volumes = set()
    for m in (1.5, 2.0, 3.0, 4.0):
        volumes.add(round(analytic_volume(LENGTH, BEAM, DRAFT, 2.0, m), 12))
        assert analytic_waterplane_inertia(LENGTH, BEAM, 2.0) == inertia
        # KB rises with m: a fuller section near the keel carries its centroid lower.
        assert analytic_vertical_centre_buoyancy(DRAFT, m + 1.0) < analytic_vertical_centre_buoyancy(
            DRAFT, m
        )
    assert len(volumes) == 4  # m does move the displacement, so BM_T = I/nabla does move
    assert analytic_metacentric_radius_transverse(
        LENGTH, BEAM, DRAFT, 2.0, 3.0
    ) != pytest.approx(analytic_metacentric_radius_transverse(LENGTH, BEAM, DRAFT, 2.0, 2.0))
    assert analytic_waterplane_area(LENGTH, BEAM, 2.0) == area


def test_invalid_exponents_are_refused_by_the_closed_forms():
    for bad in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            analytic_waterplane_inertia(LENGTH, BEAM, bad)
        with pytest.raises(ValueError):
            analytic_waterplane_inertia_longitudinal(LENGTH, BEAM, bad)


def test_kg_is_required_and_never_defaulted_to_zero():
    with pytest.raises(ValueError):
        analytic_initial_stability(reference())
    reference_state = analytic_initial_stability(reference(), kg=KG)
    assert reference_state.kg == KG
    assert reference_state.gm_transverse == pytest.approx(
        reference_state.km_transverse - KG, rel=1e-15
    )
    assert reference_state.gm_longitudinal == pytest.approx(
        reference_state.km_longitudinal - KG, rel=1e-15
    )
    assert set(reference_state.as_dict()) >= {
        "gm_transverse_m",
        "km_longitudinal_m",
        "transverse_waterplane_inertia_m4",
    }


# --------------------------------------------------------------------------- #
# 2. the mesh routes against the closed forms
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("grid", [(41, 13), (81, 25), (161, 41), (321, 81)])
def test_the_waterline_plane_cut_converges_to_the_closed_form_second_order(grid):
    parameters = reference(nx=grid[0], nz=grid[1])
    mesh = generate_wigley(parameters).mesh
    state = analytic_initial_stability(parameters, kg=KG)
    measured = initial_stability(mesh, kg=KG)
    quantities = state.compare(measured)["relative_error"]

    assert quantities["transverse_waterplane_inertia"] < 0.0  # chords cut the convex taper
    assert quantities["longitudinal_waterplane_inertia"] < 0.0
    tolerance = 3e-3 if grid == (41, 13) else 1e-4 if grid == (321, 81) else 1e-2
    assert abs(quantities["transverse_waterplane_inertia"]) < tolerance
    assert abs(quantities["kb"]) < tolerance
    assert abs(quantities["km_transverse"]) < tolerance


def test_the_inertia_error_is_second_order_in_the_grid_step():
    parameters_exponents = (2.0, 2.0)
    errors = []
    steps = []
    for grid in ((41, 13), (81, 25), (161, 41), (321, 81)):
        parameters = reference(nx=grid[0], nz=grid[1])
        state = analytic_initial_stability(parameters, kg=KG)
        measured = initial_stability(generate_wigley(parameters).mesh, kg=KG)
        error = state.compare(measured)["relative_error"]["transverse_waterplane_inertia"]
        errors.append(math.log(abs(error)))
        steps.append(math.log(1.0 / (grid[0] - 1)))
    mean_x = sum(steps) / len(steps)
    mean_y = sum(errors) / len(errors)
    order = sum((x - mean_x) * (y - mean_y) for x, y in zip(steps, errors)) / sum(
        (x - mean_x) ** 2 for x in steps
    )
    assert 1.8 < order < 2.2, f"fitted order {order:.3f} is not second order"
    assert parameters_exponents == (2.0, 2.0)


def test_the_two_mesh_routes_agree_to_round_off_on_a_polygonal_waterplane():
    for grid in ((41, 13), (161, 41), (321, 81)):
        parameters = reference(nx=grid[0], nz=grid[1])
        mesh = generate_wigley(parameters).mesh
        cut = waterplane_moments(mesh, waterline_z=0.0)
        flux = compute_hydrostatics(mesh)
        assert cut.transverse_second_moment == pytest.approx(
            flux.transverse_waterplane_inertia, rel=1e-12
        )
        assert cut.longitudinal_second_moment == pytest.approx(
            flux.longitudinal_waterplane_inertia, rel=1e-12
        )


def test_the_default_wigley_grid_is_inside_the_three_percent_band_on_every_quantity():
    parameters = reference()
    mesh = generate_wigley(parameters).mesh
    state = analytic_initial_stability(parameters, kg=KG)
    comparison = state.compare(initial_stability(mesh, kg=KG))
    worst = max(abs(value) for value in comparison["relative_error"].values())
    assert worst < 0.03, comparison["relative_error"]
    # and the published headline numbers, pinned so a code change cannot move them quietly
    assert comparison["relative_error"]["transverse_waterplane_inertia"] == pytest.approx(
        -1.4574e-03, rel=2e-3
    )
    assert comparison["relative_error"]["gm_transverse"] == pytest.approx(9.3615e-03, rel=2e-3)
    assert state.km_transverse == pytest.approx(0.527767857142857, rel=1e-12)
    assert state.gm_transverse == pytest.approx(0.027767857142857, rel=1e-12)


# --------------------------------------------------------------------------- #
# 3. the amplification: GM is a difference, and the report says so
# --------------------------------------------------------------------------- #
def test_gm_error_is_the_km_error_amplified_by_km_over_gm():
    parameters = reference()
    mesh = generate_wigley(parameters).mesh
    state = analytic_initial_stability(parameters, kg=KG)
    errors = state.compare(initial_stability(mesh, kg=KG))["relative_error"]
    amplification = abs(errors["gm_transverse"]) / abs(errors["km_transverse"])
    assert amplification > 10.0  # the closed-form value is KM/GM = 19.0
    assert amplification == pytest.approx(state.km_transverse / state.gm_transverse, rel=1e-9)


def test_the_amplification_grows_as_the_loading_condition_approaches_the_limit():
    parameters = reference()
    mesh = generate_wigley(parameters).mesh
    previous = 0.0
    for kg in (0.3, 0.45, 0.5, 0.52, 0.527):
        state = analytic_initial_stability(parameters, kg=kg)
        errors = state.compare(initial_stability(mesh, kg=kg))["relative_error"]
        amplification = abs(errors["gm_transverse"]) / abs(errors["km_transverse"])
        assert amplification > previous
        previous = amplification
    assert previous > 100.0  # 687x at KG = 0.527, where GM is 0.77 mm


def test_the_loading_limit_is_km_and_gm_crosses_zero_one_for_one_in_kg():
    parameters = reference()
    limit = analytic_initial_stability(parameters, kg=0.0).km_transverse
    assert limit == pytest.approx(0.527767857142857, rel=1e-12)
    near = analytic_initial_stability(parameters, kg=limit - 1e-6)
    past = analytic_initial_stability(parameters, kg=limit + 1e-6)
    assert near.gm_transverse > 0.0 > past.gm_transverse  # reported, not clamped
    reference_state = analytic_initial_stability(parameters, kg=KG)
    other = analytic_initial_stability(parameters, kg=KG + 0.1)
    assert other.gm_transverse - reference_state.gm_transverse == pytest.approx(-0.1, rel=1e-12)


# --------------------------------------------------------------------------- #
# 4. the artifact the example writes
# --------------------------------------------------------------------------- #
def test_the_committed_artifact_still_makes_its_claims():
    import json
    from pathlib import Path

    path = Path(__file__).resolve().parents[1] / "artifacts" / "wigley_stability" / "results.json"
    if not path.exists():
        pytest.skip("artifact not generated in this checkout")
    payload = json.loads(path.read_text())
    assert payload["claims"], "the artifact carries no claims"
    failed = [claim for claim in payload["claims"] if not claim["holds"]]
    assert not failed, failed
    default = payload["ladders"]["n=2"][0]
    state = analytic_initial_stability(reference(), kg=KG)
    assert default["inertia_transverse"] == pytest.approx(
        state.transverse_waterplane_inertia * (1.0 + default["inertia_transverse_error"]), rel=1e-12
    )
    assert payload["fitted_orders"]["n=2_inertia_transverse"] == pytest.approx(2.0, abs=0.05)