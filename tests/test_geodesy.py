"""Verification of :mod:`tensorcad.geodesy` -- two independent routes, plus anchors.

The module exists to put the *Earth model* into the ice channel's error budget, so the
tests here are the ones that can falsify it:

1. **Closed anchors.** The equator arc has the exact length ``a * dlambda`` and the
   meridian quadrant has a published WGS84 value (10 001 965.7293 m). Both are outside
   the code under test.
2. **Two routes, one number.** Vincenty's series and the Clairaut quadrature share only
   the ellipsoid constants, so their agreement is evidence. The bound is 1e-9 relative;
   the measured worst case over the fixed sample is 1.0e-11.
3. **Round trip.** Feeding the inverse solution into the direct problem must land back
   on the second point -- this catches a sign or quadrant slip that a self-consistent
   pair of errors could hide.
4. **Refusals.** The Clairaut reduction must raise where it does not apply instead of
   integrating one branch of a two-branch geodesic.
5. **The claim that motivated the module.** At 70-77 deg N the mean-radius sphere is
   short by 0.3-0.5 %, which is the same size as the 0.46 % projection distortion that
   ``docs/icechannel.md`` budgets. If that inequality ever stops holding, the report
   that cites it has to be withdrawn, so it is a test.
"""
from __future__ import annotations

import math

import pytest

from tensorcad import geodesy
from tensorcad.icechannel import DE_LONG_PEVEK_ROUTE, great_circle_km


def _legs():
    points = [(la, lo) for _, la, lo in DE_LONG_PEVEK_ROUTE]
    return list(zip(points, points[1:]))


def _random_legs(count=40, seed=11):
    import random

    rng = random.Random(seed)
    legs = []
    while len(legs) < count:
        lat1 = rng.uniform(55.0, 83.0)
        lat2 = lat1 + rng.uniform(-8.0, 8.0)
        if not 55.0 <= lat2 <= 84.0:
            continue
        lon1 = rng.uniform(130.0, 180.0)
        lon2 = lon1 + rng.uniform(-6.0, 6.0)
        legs.append((lat1, lon1, lat2, lon2))
    return legs


class TestAnchors:
    def test_equator_arc_is_the_semi_major_axis_times_the_angle(self):
        s = geodesy.vincenty_inverse(0.0, 0.0, 0.0, 1.0).distance_m
        exact = geodesy.WGS84_A * math.pi / 180.0
        assert abs(s / exact - 1.0) < 1e-12

    def test_equator_circumference_matches_the_published_wgs84_value(self):
        # 2*pi*a for a = 6378137 m is 40 075 016.6856 m; the published value is
        # 40 075 016.686 m, so the anchor tests the constant, not the algorithm.
        assert abs(2.0 * math.pi * geodesy.WGS84_A - 40075016.686) < 1e-3

    def test_meridian_quadrant_matches_the_published_wgs84_value(self):
        measured = geodesy.meridian_arc_m(0.0, 90.0)
        assert measured == pytest.approx(10001965.7293, abs=1e-3)
        # the same quadrant through the Vincenty series, i.e. a different algorithm
        assert geodesy.vincenty_inverse(0.0, 0.0, 90.0, 0.0).distance_m == pytest.approx(
            10001965.7293, abs=1e-3
        )

    def test_mean_radius_is_the_iugg_mean_radius(self):
        iugg = (2.0 * geodesy.WGS84_A + geodesy.WGS84_B) / 3.0
        assert abs(iugg - geodesy.MEAN_RADIUS_M) < 0.05

    def test_polar_flattening_matches_the_definition(self):
        assert geodesy.WGS84_B == pytest.approx(6356752.314245, abs=1e-6)


class TestTwoIndependentRoutes:
    @pytest.mark.parametrize("leg", _legs())
    def test_series_and_quadrature_agree_on_every_route_leg(self, leg):
        (lat1, lon1), (lat2, lon2) = leg
        series = geodesy.vincenty_inverse(lat1, lon1, lat2, lon2)
        quadrature = geodesy.clairaut_distance(lat1, lon1, lat2, lon2)
        assert series.converged
        assert abs(quadrature / series.distance_m - 1.0) < 1e-9

    def test_series_and_quadrature_agree_on_a_random_high_latitude_sample(self):
        worst = 0.0
        for lat1, lon1, lat2, lon2 in _random_legs():
            series = geodesy.vincenty_inverse(lat1, lon1, lat2, lon2).distance_m
            quadrature = geodesy.clairaut_distance(lat1, lon1, lat2, lon2)
            worst = max(worst, abs(quadrature / series - 1.0))
        assert worst < 1e-9

    def test_meridian_arc_agrees_across_three_routes(self):
        simpson = geodesy.meridian_arc_m(10.0, 70.0)
        series = geodesy.vincenty_inverse(10.0, 42.0, 70.0, 42.0).distance_m
        quadrature = geodesy.clairaut_distance(10.0, 42.0, 70.0, 42.0)
        assert abs(series / simpson - 1.0) < 1e-10
        assert abs(quadrature / simpson - 1.0) < 1e-10

    def test_distance_is_symmetric(self):
        for lat1, lon1, lat2, lon2 in _random_legs(count=10):
            forward = geodesy.vincenty_inverse(lat1, lon1, lat2, lon2).distance_m
            back = geodesy.vincenty_inverse(lat2, lon2, lat1, lon1).distance_m
            assert abs(forward / back - 1.0) < 1e-12

    def test_coincident_points_are_zero(self):
        solution = geodesy.vincenty_inverse(72.5, 140.0, 72.5, 140.0)
        assert solution.distance_m == 0.0
        assert solution.converged


class TestRoundTrip:
    def test_direct_problem_returns_to_the_inverse_endpoint(self):
        worst = 0.0
        for lat1, lon1, lat2, lon2 in _random_legs(count=25, seed=3):
            inverse = geodesy.vincenty_inverse(lat1, lon1, lat2, lon2)
            lat_b, lon_b, az2, iterations, converged = geodesy.vincenty_direct(
                lat1, lon1, inverse.azimuth1_deg, inverse.distance_m
            )
            assert converged and iterations < 50
            miss = geodesy.vincenty_inverse(lat_b, lon_b, lat2, lon2).distance_m
            worst = max(worst, miss)
            # the end azimuth of the inverse is the back azimuth of the direct problem
            assert abs((az2 - inverse.azimuth2_deg + 180.0) % 360.0 - 180.0) < 1e-6
        assert worst < 1e-3  # metres, over legs up to ~1000 km

    def test_route_legs_round_trip(self):
        worst = 0.0
        for (lat1, lon1), (lat2, lon2) in _legs():
            inverse = geodesy.vincenty_inverse(lat1, lon1, lat2, lon2)
            lat_b, lon_b, _, _, converged = geodesy.vincenty_direct(
                lat1, lon1, inverse.azimuth1_deg, inverse.distance_m
            )
            assert converged
            worst = max(worst, geodesy.vincenty_inverse(lat_b, lon_b, lat2, lon2).distance_m)
        assert worst < 1e-3


#: Vincenty (1975) worked example, Flinders Peak -> Buninyong, as published. The solver
#: never sees these numbers; they are the ellipse's, from the paper's own example.
VINCENTY_EXAMPLE = (-37.9510334167, 144.4248678889, -37.6528211389, 143.9264955278)
VINCENTY_EXAMPLE_DISTANCE_M = 54972.271
VINCENTY_EXAMPLE_AZIMUTH1_DEG = 306.8681583
VINCENTY_EXAMPLE_AZIMUTH2_DEG = 307.1736306

#: The seven route legs in true WGS84 metres from GeographicLib 2.1 (Karney's
#: algorithms, i.e. an implementation that shares no code with this one), plus the total.
ROUTE_REFERENCE = (
    (76.7, 149.0, 75.8, 155.0, 188170.990),
    (75.8, 155.0, 74.6, 161.0, 217168.866),
    (74.6, 161.0, 73.3, 166.0, 211724.487),
    (73.3, 166.0, 72.0, 169.5, 186026.245),
    (72.0, 169.5, 70.8, 170.5, 138545.657),
    (70.8, 170.5, 70.1, 170.6, 78186.647),
    (70.1, 170.6, 69.7, 170.3, 46084.865),
)
ROUTE_REFERENCE_TOTAL_M = 1065907.757


class TestExternalReferences:
    """Numbers this module did not produce, so agreement is evidence rather than a habit."""

    def test_vincenty_published_worked_example(self):
        solution = geodesy.vincenty_inverse(*VINCENTY_EXAMPLE)
        assert solution.converged
        assert solution.distance_m == pytest.approx(VINCENTY_EXAMPLE_DISTANCE_M, abs=1e-3)
        # The module reports azimuths in (-180, 180] deg, clockwise from north; the paper
        # prints them in [0, 360), so compare after wrapping.
        assert solution.azimuth1_deg % 360.0 == pytest.approx(VINCENTY_EXAMPLE_AZIMUTH1_DEG, abs=1e-5)
        assert solution.azimuth2_deg % 360.0 == pytest.approx(VINCENTY_EXAMPLE_AZIMUTH2_DEG, abs=1e-5)

    @pytest.mark.parametrize("lat1,lon1,lat2,lon2,reference", ROUTE_REFERENCE)
    def test_route_legs_against_geographiclib(self, lat1, lon1, lat2, lon2, reference):
        # The bound is 1 mm: GeographicLib prints these to the millimetre, so the test
        # asserts agreement with the printed value, not with the solver's own precision.
        measured = geodesy.vincenty_inverse(lat1, lon1, lat2, lon2).distance_m
        assert measured == pytest.approx(reference, abs=1e-3)

    def test_route_total_against_geographiclib(self):
        total = sum(geodesy.vincenty_inverse(a, b, c, d).distance_m
                    for a, b, c, d, _ in ROUTE_REFERENCE)
        assert total == pytest.approx(ROUTE_REFERENCE_TOTAL_M, abs=1e-3)

    def test_near_antipodal_pairs_stall_is_reported_not_hidden(self):
        # Vincenty's series is known not to converge near antipodes. The contract is that
        # the stall is visible (``converged=False``), so a caller comparing routes knows
        # the length it got is not the geodesic. The pairs that do converge must match.
        reference = {
            (0.0, 0.0, 0.5, 179.5): 19936288.578965,
            (89.0, 0.0, -89.0, 179.9): 20003746.817624,
        }
        for (lat1, lon1, lat2, lon2), distance in reference.items():
            solution = geodesy.vincenty_inverse(lat1, lon1, lat2, lon2)
            assert solution.converged
            assert solution.distance_m == pytest.approx(distance, abs=1e-3)
        for lat1, lon1, lat2, lon2 in ((0.0, 0.0, 0.0, 180.0), (0.0, 0.0, 0.0, 179.99)):
            solution = geodesy.vincenty_inverse(lat1, lon1, lat2, lon2)
            assert not solution.converged
            assert solution.iterations == 200


class TestRefusals:
    def test_quadrature_refuses_latitudes_that_straddle_the_equator(self):
        with pytest.raises(geodesy.GeodesicError):
            geodesy.clairaut_distance(-5.0, 10.0, 5.0, 10.5)

    def test_quadrature_refuses_a_parallel_off_the_equator(self):
        with pytest.raises(geodesy.GeodesicError):
            geodesy.clairaut_distance(75.0, 140.0, 75.0, 141.0)

    def test_quadrature_refuses_a_vertex_inside_the_interval(self):
        # 35 deg of longitude inside 0.1 deg of latitude: the fitted Clairaut parameter
        # puts the vertex either inside the interval (two branches) or close enough
        # outside that the integrand is nearly singular there. Both are refused: before
        # the guards, this case returned a length 1.1 % long instead of refusing, because
        # the singular point fell between Simpson nodes.
        with pytest.raises(geodesy.GeodesicError):
            geodesy.clairaut_distance(60.0, 140.0, 60.1, 175.0)
        with pytest.raises(geodesy.GeodesicError):
            geodesy.clairaut_distance(75.0, 140.0, 75.05, 155.0)

    def test_vertex_latitude_is_consistent_with_the_parallel_radius(self):
        c = geodesy.parallel_radius_m(math.radians(68.0))
        assert math.degrees(geodesy.vertex_latitude(c)) == pytest.approx(68.0, abs=1e-6)
        assert math.degrees(geodesy.vertex_latitude(geodesy.WGS84_A)) == pytest.approx(0.0)

    def test_quadrature_keeps_its_accuracy_just_outside_the_vertex(self):
        # 5 deg of longitude in 1 deg of latitude: vertex outside, branch single.
        series = geodesy.vincenty_inverse(70.0, 140.0, 71.0, 145.0).distance_m
        assert abs(geodesy.clairaut_distance(70.0, 140.0, 71.0, 145.0) / series - 1.0) < 1e-9

    def test_latitude_out_of_range_is_rejected(self):
        with pytest.raises(geodesy.GeodesicError):
            geodesy.vincenty_inverse(95.0, 0.0, 0.0, 0.0)

    def test_quadrature_still_answers_the_monotonic_branch(self):
        # The mirror image of the refusal above: an eastward run is accepted.
        assert geodesy.clairaut_distance(75.0, 140.0, 76.0, 141.0) > 0.0


class TestSphereVersusEllipsoid:
    def test_the_sphere_is_short_at_campaign_latitudes(self):
        for lat1, lon1, lat2, lon2 in _random_legs(count=25):
            error = geodesy.sphere_error(lat1, lon1, lat2, lon2)["relative"]
            assert -0.0055 < error < 0.0

    def test_the_sphere_error_is_the_same_size_as_the_projection_distortion(self):
        # docs/icechannel.md budgets a 0.46 % projection distortion at 1061 km. The
        # Earth-model term over the route is 0.39 %; both are quoted, so neither may be
        # dismissed as negligible relative to the other.
        total_sphere = 0.0
        total_ellipsoid = 0.0
        for (lat1, lon1), (lat2, lon2) in _legs():
            total_sphere += geodesy.great_circle(lat1, lon1, lat2, lon2)
            total_ellipsoid += geodesy.vincenty_inverse(lat1, lon1, lat2, lon2).distance_m
        relative = total_sphere / total_ellipsoid - 1.0
        assert relative == pytest.approx(-0.003862, abs=2e-5)

    def test_icechannel_length_is_reproduced_by_the_great_circle_helper(self):
        module_total = sum(great_circle_km(a[1], a[2], b[1], b[2]) for a, b in
                           zip(DE_LONG_PEVEK_ROUTE, DE_LONG_PEVEK_ROUTE[1:]))
        helper_total = geodesy.polyline_length_m(
            [(la, lo) for _, la, lo in DE_LONG_PEVEK_ROUTE], ellipsoid=False
        ) / 1000.0
        assert helper_total == pytest.approx(module_total, rel=1e-12)

    def test_polyline_length_switches_datum(self):
        points = [(la, lo) for _, la, lo in DE_LONG_PEVEK_ROUTE]
        sphere = geodesy.polyline_length_m(points, ellipsoid=False)
        ellipsoid = geodesy.polyline_length_m(points, ellipsoid=True)
        assert ellipsoid > sphere
        assert (sphere / ellipsoid - 1.0) == pytest.approx(-0.003862, abs=2e-5)


class TestProjection:
    def test_radial_distance_is_exact_by_construction(self):
        x, y = geodesy.azimuthal_equidistant(75.0, 155.0, (76.7, 149.0))
        assert math.hypot(x, y) == pytest.approx(
            geodesy.great_circle(76.7, 149.0, 75.0, 155.0), rel=1e-12
        )

    def test_origin_projects_to_zero(self):
        x, y = geodesy.azimuthal_equidistant(76.7, 149.0, (76.7, 149.0))
        assert abs(x) < 1e-9 and abs(y) < 1e-9

    def test_cross_track_scale_of_a_250_m_corridor(self):
        # The corridor half-width lives in plane metres on a sphere; the true WGS84
        # separation of the two boundary points is larger by the parallel-radius ratio,
        # N cos(phi) / (R cos(phi)) = N / R.
        lat = 75.8
        half = 125.0
        r_parallel = geodesy.MEAN_RADIUS_M * math.cos(math.radians(lat))
        dlon_deg = half / r_parallel * 180.0 / math.pi
        a = (lat, 155.0 - dlon_deg)
        b = (lat, 155.0 + dlon_deg)
        true_width = geodesy.vincenty_inverse(a[0], a[1], b[0], b[1]).distance_m
        sphere_width = geodesy.great_circle(a[0], a[1], b[0], b[1])
        assert sphere_width == pytest.approx(2.0 * half, rel=1e-8)
        ratio = true_width / sphere_width
        assert ratio == pytest.approx(geodesy.prime_vertical_radius(math.radians(lat))
                                      / geodesy.MEAN_RADIUS_M, rel=1e-9)
        assert 1.0042 < ratio < 1.0044