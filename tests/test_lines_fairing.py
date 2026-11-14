"""Acceptance tests for the three-direction offset-table fairing (GAP-CAD-C3)."""
from __future__ import annotations

import json

import pytest

from tensorcad.lines_fairing import (
    FairingError,
    buttock_curves,
    compare_to_reference,
    curvature_sign_flips,
    fair_offset_table,
    fairing_study,
    main,
    perturb_table,
    polyline_sign_flips,
    second_differences,
    _divided_difference_coefficients,
    _flips_per_direction,
    _TangentOperator,
    _grid_of,
    _operators,
    _solve_with_penalty,
    _free_indices,
    _Grid,
)
from tensorcad.offset_forms import offset_table


def _form(key: str = "wigley_iii", stations: int = 15, waterlines: int = 11):
    return offset_table(key, n_stations=stations, n_waterlines=waterlines)


# ---------------------------------------------------------------------------
# the discrete calculus underneath
# ---------------------------------------------------------------------------


def test_second_differences_exact_on_a_quadratic() -> None:
    t = [0.0, 0.1, 0.35, 0.9, 1.0]
    values = [3.0 * x * x - 2.0 * x + 5.0 for x in t]
    d2 = second_differences(values, t)
    assert d2[0] == 0.0 and d2[-1] == 0.0
    for index in range(1, len(t) - 1):
        assert d2[index] == pytest.approx(6.0, rel=1e-12)


def test_divided_difference_coefficients_are_exact_derivatives() -> None:
    xs = [0.0, 1.0, 2.0, 3.0, 4.0]
    cubic = [x**3 for x in xs]
    second = _divided_difference_coefficients(xs, 2, 2)
    third = _divided_difference_coefficients(xs, 2, 3)
    assert second is not None and third is not None
    assert sum(c * cubic[i] for i, c in second) == pytest.approx(12.0, rel=1e-12)
    assert sum(c * cubic[i] for i, c in third) == pytest.approx(6.0, rel=1e-12)
    assert _divided_difference_coefficients(xs, 0, 2) is None  # stencil leaves the curve


def test_sign_flips_count_inflections_not_units() -> None:
    t = [i / 10 for i in range(11)]
    assert curvature_sign_flips([x * x for x in t], t) == 0
    # a zig-zag changes curvature sign between every interior node
    zig = [0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0, 0.1, 0.0]
    assert curvature_sign_flips(zig, t) >= 5
    # same shape, metres instead of tenths -> the same count (unit invariance)
    t_m = [x * 37.0 for x in t]
    assert curvature_sign_flips(zig, t) == curvature_sign_flips(zig, t_m)


def test_buttock_tangent_operator_matches_the_analytic_second_derivative() -> None:
    """``eta_tt`` from the bilinear stencil vs the analytic value on a quadratic surface."""
    xi = [i / 12 for i in range(13)]
    zeta = [j / 10 for j in range(11)]

    def surface(a: float, b: float) -> float:
        return 0.6 + 0.2 * a - 0.1 * b * b + 0.05 * a * b

    eta = [[surface(a, b) for b in zeta] for a in xi]
    operator = _TangentOperator("buttocks", xi, zeta, eta, order=2, boundary_margin=2)
    assert operator.rows, "the tangent operator built no rows"
    values = operator.evaluate(_flat(eta, len(zeta)))
    worst = 0.0
    for node, _coefficients, _h in operator.rows:
        i, j = divmod(node, len(zeta))
        a, b = xi[i], zeta[j]
        gx, gz = 0.2 + 0.05 * b, -0.2 * b + 0.05 * a
        norm = (gx * gx + gz * gz) ** 0.5
        tx, tz = -gz / norm, gx / norm
        analytic = 2.0 * 0.05 * tx * tz + (-0.2) * tz * tz
        value = values[node]
        worst = max(worst, abs(value - analytic))
    assert worst < 5e-3, worst


def _flat(eta, m: int) -> list[float]:
    return [value for row in eta for value in row]


def test_zero_level_gradient_rows_are_reported_not_penalised() -> None:
    xi = [i / 6 for i in range(7)]
    zeta = [j / 6 for j in range(7)]
    eta = [[0.5 for _ in zeta] for _ in xi]  # constant surface: no level-curve tangent
    operator = _TangentOperator("buttocks", xi, zeta, eta, boundary_margin=1)
    assert not operator.rows
    assert operator.nodes_without_row == len(xi) * len(zeta)


# ---------------------------------------------------------------------------
# argument checking
# ---------------------------------------------------------------------------


def test_bad_requests_are_refused() -> None:
    table = _form()
    with pytest.raises(FairingError):
        fair_offset_table(table, tolerance=0.0)
    with pytest.raises(FairingError):
        fair_offset_table(table, tolerance=0.01, iterations=0)
    with pytest.raises(FairingError):
        fair_offset_table(table, tolerance=0.01, order=4)
    with pytest.raises(FairingError):
        fair_offset_table(table, tolerance=0.01, weights=(0.0, 0.0, 0.0))
    with pytest.raises(FairingError):
        perturb_table(table, amplitude=-1.0)


# ---------------------------------------------------------------------------
# the acceptance of GAP-CAD-C3-LINES-FAIRING
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def wigley_case():
    clean = _form("wigley_iii")
    noisy = perturb_table(clean, amplitude=0.004 * 0.5 * clean.beam, seed=3)
    result = fair_offset_table(
        noisy, tolerance=0.008 * 0.5 * clean.beam, iterations=3, reference=clean
    )
    return clean, noisy, result


def test_deviation_never_exceeds_the_tolerance(wigley_case) -> None:
    _, noisy, result = wigley_case
    assert result.within_tolerance
    assert result.deviation_max <= result.tolerance
    assert result.deviation_max > 0.0
    # the budget is actually used, otherwise "fairing" is a no-op dressed up
    assert result.deviation_max > 0.5 * result.tolerance


def test_curvature_sign_flips_fall_monotonically(wigley_case) -> None:
    _, _, result = wigley_case
    assert result.monotone
    assert result.flips_after < result.flips_before
    totals = [result.iterations[0].total_flips_before] + [
        iteration.total_flips_after for iteration in result.iterations if iteration.accepted
    ]
    assert len(totals) >= 2
    for earlier, later in zip(totals, totals[1:]):
        assert later <= earlier


def test_each_family_reports_its_own_evidence(wigley_case) -> None:
    _, _, result = wigley_case
    metrics = result.iterations[-1].metrics
    assert [metric.direction for metric in metrics] == ["waterlines", "stations", "buttocks"]
    for metric in metrics:
        assert metric.curves > 0
        assert metric.energy_before > 0.0
    assert any(metric.flips_reduced for metric in result.iterations[0].metrics)


def test_fixed_edges_and_non_negative_ordinates(wigley_case) -> None:
    _, noisy, result = wigley_case
    before, after = result.table_before, result.table_after
    for station in (0, before.n_stations - 1):
        assert after.half_breadths[station] == before.half_breadths[station]
    for row in after.half_breadths:
        assert row[0] == 0.0
        assert all(value >= 0.0 for value in row)
    # the digest changes, i.e. the table really was edited, and provenance survives
    assert after.digest() != before.digest()
    assert after.provenance.source == before.provenance.source
    assert after.stations == before.stations and after.waterlines == before.waterlines


def test_reference_comparison_moves_towards_the_smooth_table(wigley_case) -> None:
    clean, _, result = wigley_case
    assert result.reference is not None
    assert result.reference.rms_after < result.reference.rms_before
    direct = compare_to_reference(result.table_before, result.table_after, clean, name="analytic")
    assert direct.rms_improvement is not None and direct.rms_improvement > 1.0


def test_hydrostatics_survive_the_fairing(wigley_case) -> None:
    _, _, result = wigley_case
    change = result.hydrostatics()["relative_change"]
    assert abs(change["volume"]) < 0.02
    assert abs(change["waterplane_area"]) < 0.02


def test_deviation_argmax_names_a_cell_of_the_grid(wigley_case) -> None:
    _, _, result = wigley_case
    i, j = result.deviation_argmax
    assert 0 <= i < result.table_before.n_stations
    assert 0 <= j < result.table_before.n_waterlines


def test_two_runs_are_bit_identical(wigley_case) -> None:
    _, noisy, result = wigley_case
    again = fair_offset_table(
        noisy, tolerance=result.tolerance, iterations=3, reference=result.table_before
    )
    assert again.table_after.digest() == result.table_after.digest()


def test_an_already_fair_table_is_left_alone() -> None:
    clean = _form("wigley_iii")
    result = fair_offset_table(clean, tolerance=0.002 * 0.5 * clean.beam, iterations=2)
    assert result.flips_before == 0 and result.flips_after == 0
    assert result.deviation_max <= result.tolerance
    assert not result.flips_reduced  # nothing to reduce, and it does not pretend otherwise


# ---------------------------------------------------------------------------
# buttocks and the study
# ---------------------------------------------------------------------------


def test_buttock_contours_are_circular_level_sets_of_a_dome() -> None:
    xi = [i / 20 for i in range(21)]
    zeta = [j / 20 for j in range(21)]

    def section(a: float, b: float) -> float:
        return (1.0 - 2.0 * (a - 0.5) ** 2) * b

    eta = [[section(a, b) for b in zeta] for a in xi]
    curves = buttock_curves(eta, xi, zeta, levels=(0.5,))
    assert len(curves) == 1
    # inversion of a linear-in-zeta profile is exact: zeta = 0.5 / plan(a)
    for x, z in curves[0]:
        assert z == pytest.approx(0.5 / (1.0 - 2.0 * (x - 0.5) ** 2), rel=1e-9)
    assert polyline_sign_flips(curves[0]) == 0  # a smooth convex buttock


def test_two_tables_meet_the_criterion() -> None:
    accepted = []
    for key in ("wigley_iii", "series60_cb060", "kvlcc2_class"):
        clean = offset_table(key, n_stations=15, n_waterlines=11)
        noisy = perturb_table(clean, amplitude=0.004 * 0.5 * clean.beam, seed=5)
        result = fair_offset_table(
            noisy, tolerance=0.016 * 0.5 * clean.beam, iterations=2, reference=clean
        )
        assert result.within_tolerance, key
        assert result.monotone, key
        assert result.flips_after <= result.flips_before, key
        accepted.append(result.flips_after < result.flips_before)
    assert sum(accepted) >= 2


def test_study_payload_is_json_serialisable_and_complete() -> None:
    report = fairing_study(
        forms=("wigley_iii", "series60_cb060"),
        n_stations=13,
        n_waterlines=9,
        tolerance_fractions=(0.008, 0.016),
        iterations=2,
    )
    text = json.dumps(report)
    assert json.loads(text)["cases_total"] == 2
    for case in report["cases"]:
        assert case["analytic_flips"] == {"waterlines": 0, "stations": 0, "buttocks": 0}
        assert len(case["sweep"]) == 2
        for entry in case["sweep"]:
            assert entry["within_tolerance"]
            assert entry["monotone"]


def test_cli_runs_a_single_case(tmp_path) -> None:
    target = tmp_path / "report.json"
    code = main(
        [
            "--form", "wigley_iii",
            "--stations", "15",
            "--waterlines", "11",
            "--iterations", "2",
            "--noise-fraction", "0.004",
            "--tolerance-fraction", "0.008",
            "--json", str(target),
        ]
    )
    assert code == 0
    payload = json.loads(target.read_text())
    assert payload["accepted"] is True
    assert payload["schema"].startswith("tensorcad.lines_fairing")


def test_operator_rows_are_the_penalty_that_was_solved() -> None:
    """The energy the solver descends must be the energy the report quotes."""
    clean = _form("wigley_iii", 13, 9)
    grid = _grid_of(clean)
    operators = _operators(grid, grid.eta, order=2)
    free = _free_indices(grid.frozen())
    flat = grid.flat()
    size = len(flat)
    before = sum(operator.bend_energy(flat) for operator in operators)
    correction, _, _ = _solve_with_penalty(flat, operators, (1e-3, 1e-3, 1e-3), free, size)
    after = [flat[i] + correction[i] for i in range(size)]
    assert sum(operator.bend_energy(after) for operator in operators) < before


def test_flips_helper_matches_the_operator_directions() -> None:
    clean = _form("wigley_iii", 13, 9)
    noisy = perturb_table(clean, amplitude=0.02 * 0.5 * clean.beam, seed=11)
    grid = _grid_of(noisy)
    counts = _flips_per_direction(grid, grid.flat())
    assert set(counts) == {"waterlines", "stations", "buttocks"}
    assert sum(counts.values()) > 0
    quiet = _grid_of(clean)
    assert _flips_per_direction(quiet, quiet.flat()) == {
        "waterlines": 0,
        "stations": 0,
        "buttocks": 0,
    }