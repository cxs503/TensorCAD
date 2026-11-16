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

import importlib.metadata as _importlib_metadata
import importlib.util as _importlib_util
import math
import sys as _sys
from pathlib import Path

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

# ---------------------------------------------------------------------------
# revolved solids (bodies of revolution: SUBOFF, bulbs, shafts, sonar domes)
# ---------------------------------------------------------------------------
def test_solid_of_revolution_closes_a_sphere_to_genus_zero():
    """A revolved sphere is the primitive where the answer is known in closed form."""
    a = 1.0
    radius = lambda x: math.sqrt(max(0.0, a * a - (x - a) ** 2))  # noqa: E731
    mesh = xv.solid_of_revolution(radius, length=2 * a, stations=201, segments=128)
    topo = xv.topology_report(mesh)
    assert topo["watertight"] and topo["closed_genus"] == 0 and topo["zero_area_faces"] == 0
    volume = xv.volume_divergence(mesh)
    exact = 4.0 / 3.0 * math.pi * a**3
    # Inscribed polygon + inscribed meridian: the deficit must be negative and small.
    assert -1e-3 < xv.relative_difference(volume, exact) < 0.0
    assert abs(xv.relative_difference(xv.volume_trimesh(mesh), volume)) < 1e-9


def test_revolution_sectional_quadrature_is_exact_for_a_cylinder():
    """Simpson is exact for pi r^2 = const, so any error here is a quadrature bug."""
    volume = xv.volume_revolution_sectional(lambda x: 0.5, length=3.0, samples=101)
    assert abs(xv.relative_difference(volume, math.pi * 0.25 * 3.0)) < 1e-12


def test_solid_of_revolution_refuses_unusable_input():
    with pytest.raises(xv.VolumeXvalError):
        xv.solid_of_revolution(lambda x: -0.1, length=1.0)
    with pytest.raises(xv.VolumeXvalError):
        xv.solid_of_revolution(lambda x: float("nan"), length=1.0)
    with pytest.raises(xv.VolumeXvalError):
        xv.solid_of_revolution(lambda x: 0.1, length=1.0, segments=4)
    with pytest.raises(xv.VolumeXvalError):
        xv.volume_revolution_sectional(lambda x: 0.1, length=1.0, samples=100)


def test_half_breadth_loft_cannot_represent_a_body_of_revolution():
    """The applicability limit that forced ``solid_of_revolution`` (negative knowledge).

    A circular section that does not span the global level grid leaves repeated
    centreline vertices in the ring, so the loft is not a valid solid.  This is
    asserted rather than documented so that a future "fix" cannot silently start
    accepting it with a wrong volume.
    """
    from tensorcad import suboff as so

    length, r_max = so.SUBOFF_LENGTH_M, so.SUBOFF_RADIUS_M

    def breadth(x: float, z: float) -> float:
        r = r_max * so.suboff_radius_fraction(x / length)
        return math.sqrt(max(0.0, r * r - (z - r_max) ** 2))

    with pytest.raises(xv.VolumeXvalError):
        xv.closed_solid_from_breadth(breadth, length=length, draft=2 * r_max, stations=21, levels=7)


def test_revolved_suboff_solid_converges_to_the_published_offset_quadrature():
    """SUBOFF has no published displacement volume; its *offsets* are published.

    So the check is: the revolved solid, integrated by two independent engines,
    converges onto the meridian quadrature of the published offset table at the
    measured O(h^2), and the engines agree to the STL float32 budget.
    """
    from tensorcad import suboff as so

    length, r_max = so.SUBOFF_LENGTH_M, so.SUBOFF_RADIUS_M

    def radius(x: float) -> float:
        return r_max * so.suboff_radius_fraction(x / length)

    reference = so.suboff_reference_geometry()["volume"]
    rows = []
    for stations, segments in ((21, 42), (41, 82)):
        mesh = xv.solid_of_revolution(radius, length=length, stations=stations, segments=segments)
        topo = xv.topology_report(mesh)
        assert topo["watertight"] and topo["closed_genus"] == 0
        engine_a = xv.volume_trimesh(mesh)
        engine_b = xv.volume_numpy_stl(mesh)
        assert abs(xv.relative_difference(engine_b, engine_a)) < 1e-6
        rows.append((length / (stations - 1), abs(xv.relative_difference(engine_a, reference))))
    order = xv.fit_order(rows)
    assert order is not None and 1.5 <= order <= 2.6, order
    assert rows[-1][1] < 3e-3


# ---------------------------------------------------------------------------
# the report's claim surface: toolchain versions and the Cb anchor (m23023)
#
# A cross-validation artifact is read by people who were not there.  Two of its
# fields were wrong in a way no volume gate can catch: the second engine's version
# was recorded as "unknown" although it ran, and "calibrated to Cb 0.6505" was
# written as an identity although the calibration and the reported volume are
# integrated at two different quadrature orders.  These tests pin the *reasoning*,
# not just the numbers: the identity question must be answered by a measurement.
# ---------------------------------------------------------------------------
_REPORT_PATH = Path(__file__).resolve().parents[1] / "tools" / "volume_xval_report.py"


def _load_report():
    spec = _importlib_util.spec_from_file_location("volume_xval_report_m23023", _REPORT_PATH)
    assert spec is not None and spec.loader is not None
    module = _importlib_util.module_from_spec(spec)
    _sys.modules[spec.name] = module  # dataclasses look their module up in sys.modules
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def report():
    return _load_report()


@pytest.fixture(scope="module")
def kcs_anchor(report):
    return report.kcs_geometry_anchor()


def test_second_engine_version_comes_from_the_distribution_not_the_module(report):
    """`numpy-stl` imports as `stl` and sets no `__version__`.

    Reading the module attribute therefore recorded "unknown" for an engine that
    did run, and a machine reading the artifact concluded the engine was missing.
    The version must now be a real one and must match the installed distribution.
    """
    pytest.importorskip("stl")
    version = report.distribution_version("numpy-stl", "stl")
    assert version is not None and version != "unknown"
    assert version == _importlib_metadata.version("numpy-stl")
    assert version[0].isdigit()
    engines = report.engine_versions()
    assert engines["numpy-stl"] == version
    assert isinstance(engines["version_disagreements"], dict)


def test_absent_engine_is_none_but_unknown_version_is_not_unknown(report):
    """Three states must stay distinguishable in a machine-read artifact."""
    assert report.distribution_version("no-such-distribution-xyz", "no_such_module_xyz") is None

    class _NoVersion:
        pass

    _sys.modules["engine_without_version_xyz"] = _NoVersion()
    try:
        value = report.distribution_version("no-such-distribution-xyz", "engine_without_version_xyz")
    finally:
        _sys.modules.pop("engine_without_version_xyz", None)
    assert value == report.VERSION_NOT_REPORTED
    assert value != "unknown"


def test_cb_identity_is_measured_at_the_reporting_order_not_assumed(report, kcs_anchor):
    anchor = kcs_anchor
    lbt = anchor["principal_particulars"]["LBT_m3"]
    # the achieved Cb is the reported volume over L B T, not a copy of the request
    assert anchor["achieved_cb"] == pytest.approx(anchor["reported_volume_m3"] / lbt, rel=1e-12)
    assert anchor["requested_volume_m3"] == pytest.approx(anchor["requested_cb"] * lbt, rel=1e-12)
    # the calibration itself does hit the target at the order it bisects at ...
    assert anchor["calibration_hits_target"]
    assert abs(anchor["calibration_residual_rel"]) <= anchor["calibration_tolerance"]
    # ... but that is a different number from the one reported, so the identity
    # claim is false at the reporting order, by a factor that is computed here
    residual = anchor["achieved_cb"] / anchor["requested_cb"] - 1.0
    assert anchor["cb_residual_rel"] == pytest.approx(residual, rel=1e-12)
    assert anchor["identity_within_calibration_tolerance"] is (
        abs(residual) <= anchor["calibration_tolerance"]
    )
    assert anchor["identity_within_calibration_tolerance"] is False
    assert anchor["residual_in_units_of_calibration_tolerance"] > 100.0
    # the gap is still inside the platform's acceptance line: the row may pass, the
    # label may not claim an identity
    assert anchor["within_acceptance"] is True
    assert abs(residual) < xv.ACCEPTANCE_TOLERANCE


def test_anchor_order_ladder_is_converged_and_brackets_the_reported_number(report, kcs_anchor):
    ladder = kcs_anchor["order_ladder"]
    assert len(ladder) >= 3
    for row in ladder:  # the order ladder is a refinement, not a random walk
        assert abs(row["vs_reported"]) < 1e-4, row
        assert row["cb"] == pytest.approx(row["volume_m3"] / kcs_anchor["principal_particulars"]["LBT_m3"])
        assert row["volume_m3"] > 0.0
    finest = max(ladder, key=lambda row: row["order"][0] * row["order"][1])
    assert kcs_anchor["converged_volume_m3"] == pytest.approx(finest["volume_m3"], rel=1e-12)
    # therefore the residual against the published volume is the calibration-order
    # gap, not an unconverged quadrature
    assert abs(finest["vs_reported"]) < 1e-5


def test_full_scale_case_name_reports_the_achieved_cb(report, kcs_anchor):
    case = report.nurbs_hull_case(full_scale=True, anchor=kcs_anchor)
    assert "achieved Cb" in case.name
    assert f"{kcs_anchor['achieved_cb']:.6f}" in case.name
    assert "calibrated to" not in case.name  # the old wording asserted an identity
    assert "identity" in case.caveat and "not an identity" in case.caveat
