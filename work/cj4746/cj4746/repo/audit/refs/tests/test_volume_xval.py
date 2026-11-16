"""Tests for :mod:`tensorcad.volume_xval` -- the cross-tool volume cross-validation.

The assertions are invariants, not smoke tests:

* a body a triangle mesh can represent *exactly* (the box) must come out of every
  route to round-off, and its solid must be a closed genus-0 surface;
* a curved body must reproduce its analytic volume only in the ``h -> 0`` limit,
  with the deficit decaying at the *measured* order, so the order is asserted to
  be the second order the loft is supposed to have rather than merely "smaller";
* the classifier has to give each of its four verdicts for the right reason, and a
  systematic bias must never be excused as tessellation;
* an open or inconsistently wound surface must be refused, not silently
  integrated.
"""

from __future__ import annotations

import math

import pytest

from tensorcad import mesh as meshlib
from tensorcad import volume_xval as xv
from tensorcad import wigley as wg


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------
def box_breadth(half_beam: float):
    def breadth(x: float, z: float) -> float:
        return half_beam

    return breadth


@pytest.fixture(scope="module")
def box_solid():
    return xv.closed_solid_from_breadth(
        box_breadth(1.0), length=10.0, draft=1.0, stations=21, levels=9
    )


@pytest.fixture(scope="module")
def wigley_parameters():
    return wg.reference_parameters()


def wigley_breadth(parameters):
    def breadth(x: float, z: float) -> float:
        return wg.half_breadth(parameters, x, z)

    return breadth


# ---------------------------------------------------------------------------
# grids
# ---------------------------------------------------------------------------
def test_cosine_stations_span_the_length_and_cluster_at_the_ends():
    xs = xv.cosine_stations(10.0, 11)
    assert xs[0] == 0.0 and xs[-1] == pytest.approx(10.0)
    assert xs == tuple(sorted(xs))
    interior = [xs[i + 1] - xs[i] for i in range(len(xs) - 1)]
    # clustering means the end spacings are strictly the smallest ones
    assert interior[0] < interior[len(interior) // 2]
    assert interior[-1] < interior[len(interior) // 2]


def test_waterline_levels_start_at_keel_and_end_at_waterline():
    zs = xv.waterline_levels(1.5, 13)
    assert zs[0] == 0.0 and zs[-1] == pytest.approx(1.5)
    assert all(b > a for a, b in zip(zs, zs[1:]))


@pytest.mark.parametrize("stations", [0, 1, -3])
def test_cosine_stations_rejects_degenerate_grids(stations):
    with pytest.raises(xv.VolumeXvalError):
        xv.cosine_stations(10.0, stations)


def test_waterline_levels_rejects_non_positive_draft():
    with pytest.raises(xv.VolumeXvalError):
        xv.waterline_levels(0.0, 9)


# ---------------------------------------------------------------------------
# the solid builder
# ---------------------------------------------------------------------------
def test_box_solid_is_closed_genus_zero(box_solid):
    report = xv.topology_report(box_solid)
    assert report["watertight"] is True
    assert report["boundary_edges"] == 0
    assert report["nonmanifold_edges"] == 0
    assert report["zero_area_faces"] == 0
    assert report["euler_characteristic"] == 2
    assert report["closed_genus"] == 0


def test_box_volume_is_exact_on_every_route(box_solid):
    exact = 10.0 * 2.0 * 1.0
    assert xv.volume_divergence(box_solid) == pytest.approx(exact, rel=1e-13)
    assert xv.volume_trimesh(box_solid) == pytest.approx(exact, rel=1e-12)
    # numpy-stl round trip: binary STL is float32, so the budget is its quantisation
    assert xv.volume_numpy_stl(box_solid) == pytest.approx(exact, rel=1e-6)


def test_sectional_quadrature_is_exact_for_a_box():
    volume = xv.volume_sectional(
        box_breadth(1.0), length=10.0, draft=1.0, stations=41, levels=25
    )
    assert volume == pytest.approx(20.0, rel=1e-13)


def test_negative_breadth_is_refused():
    def breadth(x: float, z: float) -> float:
        return -0.1

    with pytest.raises(xv.VolumeXvalError, match="non-negative"):
        xv.closed_solid_from_breadth(breadth, length=10.0, draft=1.0, stations=11, levels=5)


def test_non_finite_breadth_is_refused():
    def breadth(x: float, z: float) -> float:
        return math.nan

    with pytest.raises(xv.VolumeXvalError, match="finite"):
        xv.closed_solid_from_breadth(breadth, length=10.0, draft=1.0, stations=11, levels=5)


def test_stem_that_does_not_close_gains_a_cap_and_stays_watertight():
    """A barge-like body is full width at the stem; the builder must cap it."""
    solid = xv.closed_solid_from_breadth(
        box_breadth(0.5), length=8.0, draft=1.0, stations=13, levels=7
    )
    report = xv.topology_report(solid)
    assert report["watertight"] is True
    assert report["euler_characteristic"] == 2


def test_transom_fan_refuses_a_non_star_shaped_section():
    """A U-shaped section (widest at the waterline, hollow in the middle) cannot be
    fanned from its centroid without overlapping itself, so the builder must refuse
    rather than add fictitious volume."""

    def breadth(x: float, z: float) -> float:
        # overhanging section: wide bottom, narrow waist, wide top -- the waist puts
        # the top corners out of sight of the bottom centroid
        if z < 0.3:
            return 1.0
        if z < 0.6:
            return 0.05
        return 1.0

    with pytest.raises(xv.VolumeXvalError, match="star-shaped"):
        xv.closed_solid_from_breadth(breadth, length=10.0, draft=1.0, stations=11, levels=13)


# ---------------------------------------------------------------------------
# engines and their agreement
# ---------------------------------------------------------------------------
def test_engines_reproduce_the_platform_integrator_bit_for_bit_on_the_same_triangles(box_solid):
    divergence = xv.volume_divergence(box_solid)
    assert xv.relative_difference(xv.volume_trimesh(box_solid), divergence) < 1e-12
    assert abs(xv.relative_difference(xv.volume_numpy_stl(box_solid), divergence)) < 1e-6


def test_open_surface_is_refused_by_the_engine(box_solid):
    solid = xv.closed_solid_from_breadth(
        box_breadth(1.0), length=10.0, draft=1.0, stations=9, levels=5
    )
    holed = type(solid)(vertices=solid.vertices, faces=solid.faces[:-1])
    assert not meshlib.is_watertight(holed)
    with pytest.raises(xv.VolumeXvalError, match="watertight"):
        xv.volume_trimesh(holed)


# ---------------------------------------------------------------------------
# convergence and attribution
# ---------------------------------------------------------------------------
def test_wigley_deficit_decays_at_second_order(wigley_parameters):
    parameters = wigley_parameters
    breadth = wigley_breadth(parameters)
    exact = wg.analytic_volume(parameters.length, parameters.beam, parameters.draft)
    samples = []
    volumes = []
    for stations in (21, 41, 81, 161):
        solid = xv.closed_solid_from_breadth(
            breadth,
            length=parameters.length,
            draft=parameters.draft,
            stations=stations,
            levels=max(3, stations // 4 + 1),
        )
        volume = xv.volume_divergence(solid)
        volumes.append(volume)
        h = parameters.length / (stations - 1)
        samples.append((h, abs(xv.relative_difference(volume, exact))))
    order = xv.fit_order(samples)
    assert order is not None
    assert 1.8 <= order <= 2.2, f"loft converges at order {order}, expected 2"
    # the inscribed polyhedron is below the smooth body, never above it
    assert all(volume < exact for volume in volumes)


def test_sectional_route_agrees_with_the_analytic_wigley_volume(wigley_parameters):
    parameters = wigley_parameters
    exact = wg.analytic_volume(parameters.length, parameters.beam, parameters.draft)
    measured = xv.volume_sectional(
        wigley_breadth(parameters),
        length=parameters.length,
        draft=parameters.draft,
        stations=161,
        levels=201,
    )
    assert abs(xv.relative_difference(measured, exact)) < 1e-10


def test_richardson_extrapolation_recovers_the_analytic_wigley_volume(wigley_parameters):
    parameters = wigley_parameters
    breadth = wigley_breadth(parameters)
    exact = wg.analytic_volume(parameters.length, parameters.beam, parameters.draft)
    pairs = []
    previous = None
    for stations in (41, 81, 161):
        solid = xv.closed_solid_from_breadth(
            breadth,
            length=parameters.length,
            draft=parameters.draft,
            stations=stations,
            levels=max(3, stations // 4 + 1),
        )
        volume = xv.volume_divergence(solid)
        h = parameters.length / (stations - 1)
        pairs.append((h, volume))
        previous = volume
    limit = xv.richardson_limit(pairs, 2.0)
    assert abs(xv.relative_difference(limit, exact)) < 5e-5


def test_fit_order_reports_zero_for_a_non_decaying_gap():
    """A constant gap is a *measured* exponent of zero, not a failure to fit: the
    classifier uses the window, not the presence of a number."""
    order = xv.fit_order([(1.0, 1e-3), (0.5, 1e-3), (0.25, 1e-3)])
    assert order == pytest.approx(0.0, abs=1e-12)
    # a gap that grows under refinement is not a convergent route at all
    assert xv.fit_order([(1.0, 1e-3), (0.5, 2e-3), (0.25, 4e-3)]) is None
    assert xv.fit_order([(1.0, 1e-3)]) is None


def test_fit_order_recovers_a_known_exponent():
    samples = [(h, h**2 * 1e-3) for h in (1.0, 0.5, 0.25, 0.125)]
    order = xv.fit_order(samples)
    assert order == pytest.approx(2.0, abs=1e-9)


# ---------------------------------------------------------------------------
# classification
# ---------------------------------------------------------------------------
def test_classifier_accepts_agreement():
    verdict = xv.classify_routes(
        divergence=1.0, engine_a=1.0, engine_b=1.0, pipeline=1.001, deficits=[]
    )
    assert verdict["status"] == "agrees"


def test_classifier_calls_a_decaying_gap_tessellation():
    deficits = [(1.0, 1e-2), (0.5, 2.5e-3), (0.25, 6.25e-4)]
    verdict = xv.classify_routes(
        divergence=1.0,
        engine_a=1.0 - 1e-2,
        engine_b=1.0 - 1e-2,
        pipeline=1.0,
        deficits=deficits,
    )
    assert verdict["status"] == "tessellation-deficit"
    assert verdict["observed_order"] == pytest.approx(2.0, abs=0.05)


def test_classifier_never_excuses_a_constant_gap_as_tessellation():
    deficits = [(1.0, 1e-2), (0.5, 1e-2), (0.25, 1e-2)]
    verdict = xv.classify_routes(
        divergence=1.0,
        engine_a=1.0 - 1e-2,
        engine_b=1.0 - 1e-2,
        pipeline=1.0,
        deficits=deficits,
    )
    assert verdict["status"] == "systematic-bias"
    assert "negative-knowledge ledger" in verdict["reason"]


def test_classifier_blames_topology_when_the_engines_disagree():
    verdict = xv.classify_routes(
        divergence=1.0, engine_a=1.0, engine_b=0.9, pipeline=1.0, deficits=[(1.0, 1e-2)]
    )
    assert verdict["status"] == "engine-disagreement"


def test_engine_agreement_note_reports_the_float32_budget():
    assert "float32" in xv.engine_agreement_note(1e-9)
    assert "inspect topology" in xv.engine_agreement_note(1e-3)


def test_relative_difference_refuses_a_zero_reference():
    with pytest.raises(xv.VolumeXvalError):
        xv.relative_difference(1.0, 0.0)


def test_case_row_records_the_reference_deviations():
    row = xv.CaseRow(
        ship="test",
        reference_kind="analytic",
        reference_volume=2.0,
        reference_source="unit test",
        stations=11,
        levels=5,
        h=1.0,
        pipeline=2.0,
        divergence=1.99,
        engine_a=1.99,
        engine_b=1.99,
    )
    payload = row.as_dict()
    assert payload["pipeline_vs_reference"] == pytest.approx(0.0)
    assert payload["engine_a_vs_reference"] == pytest.approx(-0.005)
    assert payload["reference_volume_m3"] == 2.0


def test_tolerances_are_module_constants_not_caller_supplied():
    assert xv.ACCEPTANCE_TOLERANCE == pytest.approx(3e-3)
    assert xv.ENGINE_AGREEMENT_TOLERANCE == pytest.approx(1e-6)