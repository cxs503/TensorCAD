"""Acceptance tests for the three-direction offset-table fairing (GAP-CAD-C3)."""
from __future__ import annotations

import json

import pytest

from tensorcad.lines_fairing import (
    FairingError,
    buttock_coverage,
    buttock_curves,
    buttock_levels,
    compare_to_reference,
    curvature_sign_flips,
    fair_offset_table,
    fairing_study,
    feature_waterlines,
    main,
    perturb_table,
    polyline_sign_flips,
    second_differences,
    _divided_difference_coefficients,
    _flips_per_direction,
    _TangentOperator,
    _grid_of,
    _plot,
    _operators,
    _solve_with_penalty,
    _free_indices,
    _Grid,
    _energy_share_per_waterline,
)
from tensorcad.offset_forms import offset_table
from tensorcad.offset_loft import OffsetProvenance, ShipOffsetTable


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
    with pytest.raises(FairingError):
        fair_offset_table(table, tolerance=0.01, features="sometimes")
    with pytest.raises(FairingError):
        fair_offset_table(table, tolerance=0.01, features=(table.n_waterlines,))
    with pytest.raises(FairingError):
        fair_offset_table(table, tolerance=0.01, features=("1",))  # type: ignore[arg-type]
    with pytest.raises(FairingError):
        feature_waterlines(table, threshold=0.0)
    with pytest.raises(FairingError):
        fairing_study(forms=("wigley_iii",), feature_modes=("sometimes",))


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
        assert len(case["sweeps"]["off"]) == 2
        for entry in case["sweeps"]["off"]:
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


# ---------------------------------------------------------------------------
# the third family: several branches, and levels chosen from the table
# ---------------------------------------------------------------------------


def _tumblehome_table(stations: int = 15, waterlines: int = 11) -> ShipOffsetTable:
    """A hull whose half-breadth falls again above the maximum breadth.

    Tumblehome is ordinary in the aft body of a real ship and it is exactly the case in
    which a level ``y = const`` is crossed twice: once on the way up from the keel and
    once above the maximum breadth.  The earlier contour treated a section as monotone in
    ``z`` and returned nothing for such a level, so the third family quietly lost those
    buttocks.
    """
    length, beam, draft = 10.0, 1.0, 0.65
    xs = [length * i / (stations - 1) for i in range(stations)]
    zs = [draft * j / (waterlines - 1) for j in range(waterlines)]

    def shape(z: float) -> float:
        zeta = z / draft
        return (zeta / 0.65) ** 0.7 if zeta <= 0.65 else 1.0 - 0.55 * (zeta - 0.65) / 0.35

    def plan(x: float) -> float:
        return 0.5 * beam * (1.0 - (2.0 * x / length - 1.0) ** 6)

    rows = tuple(tuple(max(0.0, plan(x) * shape(z)) for z in zs) for x in xs)
    return ShipOffsetTable(
        name="tumblehome test hull",
        length=length,
        stations=tuple(xs),
        waterlines=tuple(zs),
        half_breadths=rows,
        design_draft=draft,
        provenance=OffsetProvenance(
            source="synthetic: half-breadth plan with tumblehome above z/T = 0.65",
            method="y(x, z) = 0.5 B (1 - (2x/L - 1)^6) * shape(z/T) on a uniform grid",
        ),
    )


@pytest.mark.parametrize("key", ["wigley_iii", "series60_cb060", "kvlcc2_class"])
def test_buttock_levels_cover_the_third_family(key: str) -> None:
    """Every level the metric is taken at has a branch worth counting.

    The floor is a property of the table (a third of the stations), not a hard-coded
    fraction of the peak ordinate, so the guarantee survives a change of grid: the
    earlier fixed fractions (0.2 ... 0.8 of the global peak) left the near-keel and
    near-deck part of the buttock plan unmeasured and on a coarse table dropped levels
    entirely.
    """
    table = offset_table(key, n_stations=21, n_waterlines=17)
    grid = _grid_of(table)
    floor = max(5, round(0.3 * grid.n))
    coverage = buttock_coverage(grid.eta, grid.xi, grid.zeta)
    assert coverage["curves"] >= len(coverage["levels"]) == 5
    for detail in coverage["detail"]:
        assert detail["branches"] >= 1, detail
        assert max(detail["samples"]) >= floor, detail
        assert detail["level"] <= 0.95 * max(max(row) for row in grid.eta)
    # the levels span the buttock plan instead of huddling around the midbreadth
    assert coverage["levels"][0] <= 0.1 * max(max(row) for row in grid.eta)
    assert coverage["levels"][-1] >= 0.8 * max(max(row) for row in grid.eta)
    # and they are what the flip count is actually taken on
    curves = buttock_curves(grid.eta, grid.xi, grid.zeta, levels=coverage["levels"])
    assert len(curves) == coverage["curves"]


def test_buttock_contour_keeps_both_branches_of_a_tumblehome_section() -> None:
    xi = [i / 8 for i in range(9)]
    zeta = [j / 10 for j in range(11)]
    # rises to 0.5 at zeta = 0.6, then tumbles home to 0.15 at the deck
    shape = [0.5 * min(1.0, (z / 0.6) ** 6.0) if z > 0.0 else 0.0 for z in zeta]
    shape = [value if z <= 0.6 else 0.5 * (1.0 - 0.7 * (z - 0.6) / 0.4)
             for value, z in zip(shape, zeta)]
    eta = [[value for value in shape] for _ in xi]

    level = 0.25  # crossed below and above the maximum breadth
    curves = buttock_curves(eta, xi, zeta, levels=(level,))
    assert len(curves) == 2
    rising, falling = sorted(curve[0][1] for curve in curves)
    assert 0.0 < rising < 0.6 < falling < 1.0
    for curve in curves:
        assert [point[0] for point in curve] == xi  # one sample per station
        assert len({round(point[1], 12) for point in curve}) == 1  # a buttock is z = const
    # the deck ordinate of this section is 0.15 < level, so the old endpoint test would
    # have reported no buttock at all for this level
    assert max(row[-1] for row in eta) < level


def test_buttock_contour_splits_where_the_level_set_disappears() -> None:
    """A level set that stops and starts again along the length is two buttocks."""
    xi = [i / 10 for i in range(11)]
    zeta = [j / 10 for j in range(11)]
    # a "sponson": wide amidships only, narrow at the ends, so y = 0.4 exists over two
    # separate stretches of x
    widths = [0.6, 0.6, 0.6, 0.6, 0.2, 0.2, 0.2, 0.6, 0.6, 0.6, 0.6]
    eta = [[width * z for z in zeta] for width in widths]
    curves = buttock_curves(eta, xi, zeta, levels=(0.4,))
    assert len(curves) == 2
    assert [len(curve) for curve in curves] == [4, 4]
    assert max(point[0] for point in curves[0]) < min(point[0] for point in curves[1])


def test_fairing_a_tumblehome_table_reports_the_third_family() -> None:
    """End to end on a hull the old contour could not measure.

    The acceptance does not change (deviation inside the tolerance, flip count not
    rising), but the third family now contributes: the level sets are branched, the
    flips are non-zero before the fairing and the faired table keeps them from rising.
    """
    clean = _tumblehome_table()
    half_beam = 0.5 * clean.beam
    unfair = perturb_table(clean, amplitude=0.004 * half_beam, seed=3)
    result = fair_offset_table(
        unfair, tolerance=0.032 * half_beam, iterations=3, reference=clean
    )
    assert result.within_tolerance
    assert result.monotone
    assert result.flips_after <= result.flips_before

    grid = _grid_of(unfair)
    coverage = buttock_coverage(grid.eta, grid.xi, grid.zeta)
    assert coverage["curves"] > len(coverage["levels"])  # tumblehome sections branch
    buttocks = [
        metric for metric in result.iterations[0].metrics if metric.direction == "buttocks"
    ][0]
    assert buttocks.curves > len(coverage["levels"])
    assert buttocks.flips_before > 0  # the old contour reported zero here
    assert buttocks.flips_after <= buttocks.flips_before


# ---------------------------------------------------------------------------
# design-feature rows: the near-keel band of a power-law section
# ---------------------------------------------------------------------------


def test_energy_share_names_the_keel_band_of_a_power_law_section() -> None:
    """One row owning the whole zeta-curvature energy is the detection signal."""
    for key, others in (("series60_cb060", 0.01), ("kvlcc2_class", 0.01)):
        grid = _grid_of(_form(key, stations=15, waterlines=11))
        share = _energy_share_per_waterline(grid)
        assert sum(share.values()) == pytest.approx(1.0, rel=1e-9)
        assert share[1] > 0.5, key
        assert max(value for index, value in share.items() if index != 1) < others, key

    wigley = _energy_share_per_waterline(_grid_of(_form("wigley_iii", stations=15, waterlines=11)))
    assert max(wigley.values()) < 0.3  # a parabola concentrates nothing


def test_feature_waterlines_flags_the_power_law_forms_only() -> None:
    assert feature_waterlines(_form("wigley_iii")) == ()
    assert feature_waterlines(_form("series60_cb060")) == (1,)
    assert feature_waterlines(_form("kvlcc2_class")) == (1,)
    # the frozen keel and deck rows are never "features": they are not unknowns
    assert 0 not in feature_waterlines(_form("series60_cb060"))
    assert 10 not in feature_waterlines(_form("series60_cb060"))


def test_features_off_reproduces_the_plain_energy_bit_for_bit() -> None:
    clean = _form("series60_cb060")
    noisy = perturb_table(clean, amplitude=0.004 * 0.5 * clean.beam, seed=11)
    plain = fair_offset_table(noisy, tolerance=0.016 * 0.5 * clean.beam, iterations=2)
    explicit_off = fair_offset_table(
        noisy, tolerance=0.016 * 0.5 * clean.beam, iterations=2, features=()
    )
    assert plain.table_after.digest() == explicit_off.table_after.digest()
    assert plain.feature_waterlines == ()
    assert all(iteration.penalty_rows_dropped == 0 for iteration in plain.iterations)

    switched_on = fair_offset_table(
        noisy, tolerance=0.016 * 0.5 * clean.beam, iterations=2, features="auto"
    )
    assert switched_on.feature_waterlines == (1,)
    assert all(iteration.penalty_rows_dropped > 0 for iteration in switched_on.iterations)
    assert switched_on.table_after.digest() != plain.table_after.digest()


def _series60_budget(fraction: float = 0.008):
    """The study's Series 60 case: the configuration the docs report."""
    clean = offset_table("series60_cb060", n_stations=21, n_waterlines=17)
    noisy = perturb_table(clean, amplitude=0.004 * 0.5 * clean.beam, seed=20261114)
    return clean, noisy, fraction * 0.5 * clean.beam


def test_declaring_the_keel_band_a_feature_moves_towards_the_smooth_table() -> None:
    """The documented limitation of the plain energy, measured on the same table."""
    clean, noisy, tolerance = _series60_budget(0.008)
    plain = fair_offset_table(noisy, tolerance=tolerance, iterations=3, reference=clean)
    featured = fair_offset_table(
        noisy, tolerance=tolerance, iterations=3, reference=clean, features="auto"
    )
    assert featured.within_tolerance and featured.monotone
    assert featured.flips_after < plain.flips_after          # more sign flips removed
    assert featured.reference.rms_after < plain.reference.rms_after
    # and the improvement is a real one: closer to the smooth table than the input was
    # (the plain energy is further from it than the unfaired table -- the limitation)
    assert featured.reference.rms_after < featured.reference.rms_before
    assert plain.reference.rms_after > plain.reference.rms_before
    assert featured.reference.rms_before == pytest.approx(plain.reference.rms_before)


def test_the_feature_switch_does_not_rescue_a_budget_far_above_the_noise() -> None:
    """Honest limit: at 3.2 % B/2 the max-norm budget is 8x the jitter, and spending it
    moves the table away from the smooth shape in both modes.  The feature rows give the
    largest flip reduction of the whole sweep (180 -> 18) and the worst reference
    agreement, which is why the docs recommend the small budgets for this switch."""
    clean, noisy, tolerance = _series60_budget(0.032)
    plain = fair_offset_table(noisy, tolerance=tolerance, iterations=3, reference=clean)
    featured = fair_offset_table(
        noisy, tolerance=tolerance, iterations=3, reference=clean, features="auto"
    )
    assert featured.flips_after < plain.flips_after
    assert featured.reference.rms_after > featured.reference.rms_before
    assert featured.within_tolerance and featured.monotone


def test_a_wigley_table_gains_nothing_from_the_feature_switch() -> None:
    """No detected row means no change: the switch cannot flatter a smooth form."""
    clean = _form("wigley_iii")
    noisy = perturb_table(clean, amplitude=0.004 * 0.5 * clean.beam, seed=5)
    tolerance = 0.016 * 0.5 * clean.beam
    plain = fair_offset_table(noisy, tolerance=tolerance, iterations=2)
    automatic = fair_offset_table(noisy, tolerance=tolerance, iterations=2, features="auto")
    assert automatic.feature_waterlines == ()
    assert automatic.table_after.digest() == plain.table_after.digest()


def test_study_measures_both_feature_modes_on_the_same_tables() -> None:
    report = fairing_study(
        forms=("series60_cb060", "wigley_iii"),
        n_stations=13,
        n_waterlines=9,
        tolerance_fractions=(0.016,),
        iterations=2,
        feature_modes=("off", "auto"),
    )
    assert report["feature_modes"] == ["off", "auto"]
    for case in report["cases"]:
        assert set(case["sweeps"]) == {"off", "auto"}
        assert case["accepted_at"] == case["accepted_at_by_mode"]["off"]
        for mode, sweep in case["sweeps"].items():
            assert len(sweep) == 1
            assert sweep[0]["feature_mode"] == mode
        assert case["settings"]["feature_modes"] == ["off", "auto"]
        # the analytic table is what both modes are compared against
        assert case["sweeps"]["off"][0]["reference"]["rms_before_m"] == pytest.approx(
            case["sweeps"]["auto"][0]["reference"]["rms_before_m"]
        )
    series = next(case for case in report["cases"] if case["form"] == "series60_cb060")
    assert series["sweeps"]["off"][0]["feature_waterlines"] == []
    assert series["sweeps"]["auto"][0]["feature_waterlines"] == [1]
    wigley = next(case for case in report["cases"] if case["form"] == "wigley_iii")
    assert wigley["sweeps"]["auto"][0]["feature_waterlines"] == []
    assert set(report["cases_meeting_criterion_at_a_tolerance_by_mode"]) == {"off", "auto"}


# ---------------------------------------------------------------------------
# the report describes the delivered table, not a discarded trial
# ---------------------------------------------------------------------------


def _rejected_trial_case():
    """A pinned run whose outer loop rejects its last trial (table comes back unchanged)."""
    clean = offset_table("wigley_iii", n_stations=11, n_waterlines=9)
    noisy = perturb_table(clean, amplitude=0.004 * 0.5 * clean.beam, seed=1)
    result = fair_offset_table(
        noisy, tolerance=3 * 0.004 * 0.5 * clean.beam, iterations=3
    )
    return noisy, result


def _flips_of(table) -> int:
    """Independent re-measurement of the three-family flip count of a table."""
    grid = _grid_of(table)
    return sum(_flips_per_direction(grid, grid.flat()).values())


def test_a_rejected_trial_does_not_leak_into_the_reported_flip_count() -> None:
    noisy, result = _rejected_trial_case()
    assert result.rejected_steps == 1
    assert not result.last_trial_accepted
    # the rejected trial raised the count; the delivered table is the untouched input
    rejected = result.iterations[-1]
    assert rejected.total_flips_after > rejected.total_flips_before
    assert result.table_after.digest() == result.table_before.digest()
    assert result.flips_after == _flips_of(result.table_after) == result.flips_before
    assert result.flips_after != rejected.total_flips_after
    assert not result.flips_reduced  # nothing was delivered, so nothing was reduced


def test_published_flip_counts_match_the_delivered_table() -> None:
    """Every number in the payload is re-derived from the ordinates the caller gets."""
    cases = [
        ("wigley_iii", 1, 0.004, 3),      # rejected trial, table unchanged
        ("kvlcc2_class", 3, 0.030, 3),    # rejected trial, earlier step accepted
        ("series60_cb060", 1, 0.004, 3),  # rejected trial after a measured improvement
    ]
    for form, seed, amplitude, multiplier in cases:
        clean = offset_table(form, n_stations=11, n_waterlines=9)
        noisy = perturb_table(clean, amplitude=amplitude * 0.5 * clean.beam, seed=seed)
        result = fair_offset_table(
            noisy, tolerance=multiplier * amplitude * 0.5 * clean.beam, iterations=3
        )
        payload = result.to_dict()
        assert result.flips_after == _flips_of(result.table_after), form
        assert payload["flips_after"] == result.flips_after
        assert payload["flips_after_by_direction"] == result.flips_after_by_direction
        assert payload["flips_after"] == sum(payload["flips_after_by_direction"].values())
        assert payload["flips_before"] == _flips_of(result.table_before)
        assert payload["last_trial_accepted"] == result.last_trial_accepted
        json.dumps(payload)  # the audit trail stays serialisable


def test_an_unchanged_table_is_never_accepted() -> None:
    """Anti-false-accept invariant: an accepted run must have edited the table."""
    _, result = _rejected_trial_case()
    assert result.table_after.digest() == result.table_before.digest()
    assert not result.accepted()
    for form, seed in (("wigley_iii", 3), ("series60_cb060", 11)):
        clean = offset_table(form, n_stations=11, n_waterlines=9)
        noisy = perturb_table(clean, amplitude=0.004 * 0.5 * clean.beam, seed=seed)
        run = fair_offset_table(
            noisy, tolerance=0.016 * 0.5 * clean.beam, iterations=2
        )
        if run.accepted():
            assert run.table_after.digest() != run.table_before.digest(), form
            assert run.flips_after < run.flips_before, form


def test_study_reports_delivered_flips_and_keeps_the_rejected_trial() -> None:
    report = fairing_study(
        forms=("wigley_iii",),
        n_stations=11,
        n_waterlines=9,
        noise_fraction=0.004,
        tolerance_fractions=(0.012,),
        iterations=3,
        seed=1,
        feature_modes=("off",),
    )
    entry = report["cases"][0]["sweeps"]["off"][0]
    assert entry["rejected_steps"] == 1 and entry["last_trial_accepted"] is False
    assert entry["flips_after"] == sum(entry["flips_by_direction"][key]["after"] or 0
                                      for key in entry["flips_by_direction"])
    for key, pair in entry["flips_by_direction"].items():
        assert pair["after"] == entry["flips_after_by_direction"][key]
    # the discarded trial is still in the audit trail, with its own number
    last = entry["iterations"][-1]
    assert last["accepted"] is False
    assert sum(d["flips_after"] for d in last["directions"]) != entry["flips_after"]


def test_the_figure_survives_a_rejected_trial(tmp_path) -> None:
    """The histogram must plot the delivered value, with the refused trial off the path."""
    pytest.importorskip("matplotlib")
    _, result = _rejected_trial_case()
    assert not result.last_trial_accepted
    target = _plot(result, tmp_path / "rejected.png", form="wigley_iii")
    assert target.exists() and target.stat().st_size > 0
