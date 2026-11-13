"""Offset table -> NURBS loft: integration, provenance and acceptance tests.

The three things that can go wrong silently here are (1) a quadrature rule that is not
what it claims, (2) a table that is no longer the table its provenance says it is, and
(3) a loft that agrees with its own table but not with the hull.  Each has a test:

* the Simpson rules are checked against exact polynomial integrals, including the
  unequal-spacing branch and the 3/8 hybrid;
* the bundle round trip is checked, and a single edited ordinate must be rejected;
* the loft is compared against the closed form of the Wigley hull (``Cb = 4/9``) and,
  for every landed table, against the table Simpson integrals at the ``0.5 %``
  acceptance used by ``examples/offset_loft/run_example.py``.
"""
from __future__ import annotations

import math

import pytest

from tensorcad.offset_forms import builtin_forms, offset_table, reference_for, sample_offset_table
from tensorcad.offset_loft import (
    OffsetProvenance,
    OffsetTableError,
    ShipOffsetTable,
    compare_offset_table,
    loft_offset_table,
    read_offset_bundle,
    read_offset_csv,
    simpson_integrate,
    station_sensitivity,
    table_integrals,
    write_offset_bundle,
)


# ---------------------------------------------------------------------------
# quadrature
# ---------------------------------------------------------------------------
def test_simpson_is_exact_for_cubics_and_handles_odd_counts():
    xs = [0.0 + 0.125 * i for i in range(9)]  # 8 intervals, 1/3 rule
    assert simpson_integrate([x**3 for x in xs], xs) == pytest.approx(0.25, rel=1e-14)
    xs = [0.0 + 0.125 * i for i in range(8)]  # 7 intervals, 3/8 + 1/3 hybrid
    assert simpson_integrate([x**3 for x in xs], xs) == pytest.approx(7**4 / 4 * 0.125**4, rel=1e-14)


def test_unequal_spacing_rule_is_exact_for_quadratics_and_matches_simpson_when_uniform():
    clustered = [0.0, 0.001, 0.004, 0.01, 0.03, 0.07, 0.15, 0.3, 0.5]
    assert simpson_integrate([z**2 for z in clustered], clustered) == pytest.approx(
        0.5**3 / 3, rel=1e-13
    )
    uniform = [0.25 * i for i in range(9)]
    values = [3.0 - 2.0 * x + 5.0 * x**2 for x in uniform]
    classical = simpson_integrate(values, uniform)
    # the same data through the unequal-spacing branch (detected as uniform, so feed a
    # grid that is uniform to within 1e-12 but built by the general path)
    perturbed = [x * (1.0 + 1e-12 * (i % 2)) for i, x in enumerate(uniform)]
    assert simpson_integrate([3.0 - 2.0 * x + 5.0 * x**2 for x in perturbed], perturbed) == pytest.approx(
        classical, rel=1e-9
    )


def test_unequal_spacing_rule_rejects_an_odd_interval_count():
    grid = [0.0, 0.1, 0.25, 0.4, 0.6, 0.85, 1.0]  # 6 intervals -> fine
    simpson_integrate([1.0] * len(grid), grid)
    grid = [0.0, 0.1, 0.25, 0.4, 0.6, 0.85]  # 5 intervals -> rejected
    with pytest.raises(OffsetTableError, match="odd"):
        simpson_integrate([1.0] * len(grid), grid)


# ---------------------------------------------------------------------------
# table validation and provenance
# ---------------------------------------------------------------------------
def test_table_rejects_a_non_uniform_station_grid_and_a_missing_design_waterline():
    form = builtin_forms()["wigley_iii"]
    spec = form.specimen()
    stations = tuple(0.0 + i * 0.5 ** (1 + 0.1 * i) for i in range(6))
    with pytest.raises(OffsetTableError, match="uniform"):
        ShipOffsetTable(
            name="bad",
            length=stations[-1],
            stations=stations,
            waterlines=(0.0, 0.1, 0.2, 0.3, 0.4),
            half_breadths=tuple((0.0, 0.1, 0.2, 0.2, 0.2) for _ in stations),
            design_draft=0.4,
            provenance=form.provenance,
        )
    with pytest.raises(OffsetTableError, match="design draft"):
        ShipOffsetTable(
            name="bad",
            length=2.0,
            stations=(0.0, 1.0, 2.0),
            waterlines=(0.0, 0.1, 0.2, 0.3, 0.4),
            half_breadths=tuple((0.0, 0.1, 0.2, 0.2, 0.2) for _ in range(3)),
            design_draft=0.35,
            provenance=form.provenance,
        )
    assert spec.draft > 0.0


def test_provenance_requires_a_source_and_a_method():
    with pytest.raises(OffsetTableError, match="source"):
        OffsetProvenance(source="  ", method="closed form")
    with pytest.raises(OffsetTableError, match="method"):
        OffsetProvenance(source="Todd 1951", method="")


def test_bundle_round_trip_and_tamper_detection(tmp_path):
    table = offset_table("wigley_iii", n_stations=13, n_waterlines=9)
    csv_path, json_path = write_offset_bundle(table, tmp_path)
    reloaded = read_offset_bundle(csv_path)
    assert reloaded.digest() == table.digest()
    assert reloaded.provenance.source == table.provenance.source
    assert reloaded.half_breadths == table.half_breadths

    # edit one ordinate in the CSV: the recorded digest must reject it
    lines = csv_path.read_text(encoding="utf-8").splitlines()
    header_room = max(i for i, line in enumerate(lines) if line.startswith("#"))
    cells = lines[header_room + 1].split(",")
    cells[3] = repr(float(cells[3]) + 1e-6)
    lines[header_room + 1] = ",".join(cells)
    csv_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    with pytest.raises(OffsetTableError, match="digest mismatch"):
        read_offset_bundle(csv_path)
    assert json_path.exists()


def test_csv_without_provenance_is_rejected(tmp_path):
    path = tmp_path / "naked.csv"
    path.write_text("x\\z,0.0,0.1,0.2\n0.0,0.0,0.0,0.0\n1.0,0.0,0.1,0.2\n2.0,0.0,0.0,0.0\n", encoding="utf-8")
    with pytest.raises(OffsetTableError, match="provenance"):
        read_offset_csv(path)


# ---------------------------------------------------------------------------
# the loft against closed forms
# ---------------------------------------------------------------------------
def test_wigley_loft_reproduces_the_closed_form_volume_and_waterplane():
    table = offset_table("wigley_iii", n_stations=21, n_waterlines=17)
    reference = reference_for("wigley_iii")
    loft = loft_offset_table(table)
    exact_volume = 4.0 / 9.0 * table.length * table.beam * table.design_draft
    exact_waterplane = 2.0 / 3.0 * table.length * table.beam
    assert exact_volume == pytest.approx(reference["volume"], rel=1e-12)
    assert loft.volume() == pytest.approx(exact_volume, rel=1e-9)
    assert loft.waterplane_area() == pytest.approx(exact_waterplane, rel=1e-9)
    assert loft.vcb() == pytest.approx(5.0 / 8.0 * table.design_draft, rel=1e-9)
    assert loft.lcb() == pytest.approx(0.5 * table.length, abs=1e-9)


def test_loft_passes_through_every_tabulated_ordinate():
    table = offset_table("series60_cb060", n_stations=15, n_waterlines=9, waterline_exponent=3.0)
    loft = loft_offset_table(table)
    assert loft.station_deviation() < 1e-9 * table.length


def test_table_integrals_of_wigley_match_the_closed_form():
    table = offset_table("wigley_iii", n_stations=21, n_waterlines=17)
    baseline = table_integrals(table)
    assert baseline.block_coefficient == pytest.approx(4.0 / 9.0, rel=1e-9)
    assert baseline.waterplane_coefficient == pytest.approx(2.0 / 3.0, rel=1e-9)
    assert baseline.vcb == pytest.approx(5.0 / 8.0 * table.design_draft, rel=1e-9)


# ---------------------------------------------------------------------------
# acceptance: the loft against the table's own Simpson integrals
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "key, config",
    [
        ("wigley_iii", dict(n_stations=21, n_waterlines=17, waterline_exponent=1.0)),
        ("series60_cb060", dict(n_stations=21, n_waterlines=17, waterline_exponent=3.0)),
        ("kvlcc2_class", dict(n_stations=41, n_waterlines=17, waterline_exponent=3.0)),
    ],
)
def test_landed_tables_pass_the_half_percent_acceptance(key, config):
    table = offset_table(key, **config)
    comparison = compare_offset_table(
        table, reference=reference_for(key), tolerance=5e-3, mesh=False
    )
    assert comparison.volume_error_nurbs < 5e-3
    assert comparison.waterplane_error_nurbs < 5e-3
    assert comparison.passed


def test_loft_mesh_route_agrees_with_the_table():
    table = offset_table("series60_cb060", n_stations=21, n_waterlines=17, waterline_exponent=3.0)
    comparison = compare_offset_table(
        table, tolerance=5e-3, mesh=True, n_u=97, n_v=33
    )
    assert comparison.volume_error_mesh < 5e-3
    assert comparison.waterplane_error_mesh < 5e-3


def test_fast_basis_matrix_path_matches_the_surface_derivative_path():
    table = offset_table("wigley_iii", n_stations=13, n_waterlines=9)
    loft = loft_offset_table(table)
    assert loft.volume(order=8) == pytest.approx(loft.volume_direct(order=8), rel=1e-10)
    assert loft.waterplane_area(order=8) == pytest.approx(
        loft.waterplane_area_direct(order=8), rel=1e-10
    )
    assert loft.section_half_area(0.5, order=8) == pytest.approx(
        loft.section_half_area_direct(0.5, order=8), rel=1e-10
    )


def test_station_sensitivity_converges_and_reports_both_baselines():
    reference = reference_for("series60_cb060")

    def build(n_stations: int) -> ShipOffsetTable:
        return offset_table(
            "series60_cb060",
            n_stations=n_stations,
            n_waterlines=17,
            waterline_exponent=3.0,
        )

    denominator = table_integrals(build(41))
    rows = station_sensitivity(build, (11, 21, 41), reference=reference, denominator=denominator)
    assert [row["n_stations"] for row in rows] == [11, 21, 41]
    assert rows[-1]["relative_error_volume_vs_denominator"] < 5e-3
    # against the closed form the ladder must be monotone; against another table it is
    # not (two discretisations of a sharp-ended plan form), which is why both are reported
    assert rows[0]["relative_error_volume_vs_reference"] > rows[1]["relative_error_volume_vs_reference"]
    assert rows[1]["relative_error_volume_vs_reference"] > rows[2]["relative_error_volume_vs_reference"]
    # the finest run's loft-vs-table gap is the section-representation difference, not zero
    assert rows[-1]["relative_error_volume_vs_own_table"] < 5e-3
    assert "relative_error_volume_vs_reference" in rows[0]


def test_clustered_waterlines_resolve_the_power_law_section():
    uniform = offset_table("series60_cb060", n_stations=21, n_waterlines=17)
    clustered = offset_table(
        "series60_cb060", n_stations=21, n_waterlines=17, waterline_exponent=3.0
    )
    reference = reference_for("series60_cb060")["volume"]
    uniform_error = abs(loft_offset_table(uniform).volume() - reference) / reference
    clustered_error = abs(loft_offset_table(clustered).volume() - reference) / reference
    assert uniform_error > 1e-3  # a uniform grid on zeta**(1/24) cannot resolve the keel
    assert clustered_error < 1e-4
    assert clustered_error < uniform_error


def test_sample_offset_table_rejects_a_degenerate_grid():
    spec = builtin_forms()["wigley_iii"].specimen()
    with pytest.raises(OffsetTableError, match="n_waterlines"):
        sample_offset_table(spec, provenance=builtin_forms()["wigley_iii"].provenance, n_waterlines=2)
    table = sample_offset_table(
        spec, provenance=builtin_forms()["wigley_iii"].provenance, n_stations=9, n_waterlines=5
    )
    assert table.half_breadths[0] == pytest.approx((0.0,) * 5)
    assert math.isclose(table.half_breadths[-1][-1], 0.0, abs_tol=1e-12)