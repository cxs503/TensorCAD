"""Evidence for the ice channel corridor: the projection law, the solid, and the guard.

The corridor is the campaign's ice-going domain (`docs/icechannel.md`). These
tests do not check the module against itself: the projection's closed form is
checked against the *direct geodesic problem*, the solid against an analytic
prism, and the repeated-waypoint guard against a negative control that lifts it
and shows the folded solid the guard exists to prevent.
"""

from __future__ import annotations

import math

import pytest

from tensorcad import icechannel as ice
from tensorcad import mesh as meshlib

ROUTE = ice.DE_LONG_PEVEK_ROUTE
ORIGIN = (ROUTE[0][1], ROUTE[0][2])


def direct(lat: float, lon: float, bearing_deg: float, distance_km: float) -> tuple[float, float]:
    """Destination by the direct geodesic problem on the same sphere."""
    phi1, lam1, theta = math.radians(lat), math.radians(lon), math.radians(bearing_deg)
    delta = distance_km / ice.EARTH_RADIUS_KM
    phi2 = math.asin(
        math.sin(phi1) * math.cos(delta) + math.cos(phi1) * math.sin(delta) * math.cos(theta)
    )
    lam2 = lam1 + math.atan2(
        math.sin(theta) * math.sin(delta) * math.cos(phi1),
        math.cos(delta) - math.sin(phi1) * math.sin(phi2),
    )
    return math.degrees(phi2), math.degrees(lam2)


def plane_azimuth(vector: tuple[float, float]) -> float:
    return math.degrees(math.atan2(vector[0], vector[1])) % 360.0


# --------------------------------------------------------------------------- #
# the scale law
# --------------------------------------------------------------------------- #
def test_factor_is_one_at_the_origin_in_every_direction() -> None:
    for azimuth in range(0, 360, 45):
        assert ice.plane_to_true_factor(0.0, 0.0, float(azimuth)) == 1.0


def test_radial_step_is_exact_and_tangential_step_is_c_over_sin_c() -> None:
    distance_m = 500_000.0
    central = distance_m / (ice.EARTH_RADIUS_KM * 1000.0)
    tangent = central / math.sin(central)
    # a vector pointing straight away from (or at) the origin is radial: no error
    assert ice.plane_to_true_factor(0.0, distance_m, 0.0) == pytest.approx(1.0, rel=1e-15)
    assert ice.plane_to_true_factor(0.0, distance_m, 180.0) == pytest.approx(1.0, rel=1e-15)
    # an east-west vector at that distance is tangential: true = plane * sin(c)/c
    assert ice.plane_to_true_factor(0.0, distance_m, 90.0) == pytest.approx(
        math.sin(central) / central, rel=1e-15
    )


@pytest.mark.parametrize("distance_km", [10.0, 150.0, 500.0, 1061.79, 1500.0])
@pytest.mark.parametrize("bearing_deg", [0.0, 37.0, 90.0, 143.0, 210.0, 305.0])
def test_closed_form_matches_the_direct_geodesic_problem(
    distance_km: float, bearing_deg: float
) -> None:
    """A 1 m true step: predicted true length vs the great circle through it."""
    start = direct(ORIGIN[0], ORIGIN[1], 55.0, distance_km)
    end = direct(start[0], start[1], bearing_deg, 0.001)
    first = ice._project(start[0], start[1], ORIGIN)
    second = ice._project(end[0], end[1], ORIGIN)
    plane = math.dist(first, second)
    truth = ice.great_circle_km(start[0], start[1], end[0], end[1]) * 1000.0
    azimuth = plane_azimuth((second[0] - first[0], second[1] - first[1]))
    predicted = ice.plane_to_true_factor(first[0], first[1], azimuth) * plane
    assert predicted == pytest.approx(truth, rel=1e-9)


def test_factor_rejects_non_finite_input() -> None:
    with pytest.raises(ValueError, match="must be finite"):
        ice.plane_to_true_factor(float("nan"), 0.0, 0.0)
    with pytest.raises(ValueError, match="must be finite"):
        ice.plane_to_true_factor(0.0, 0.0, float("inf"))


def test_tangential_scale_reaches_half_a_percent_at_campaign_scale() -> None:
    """The distortion is 0.46 % at 1061 km -- too large to call it zero, small enough to correct."""
    central = 1061.79 / ice.EARTH_RADIUS_KM
    factor = ice.plane_to_true_factor(0.0, 1061.79 * 1000.0, 90.0)
    assert 1.0 - factor == pytest.approx(1.0 - math.sin(central) / central, rel=1e-12)
    assert 0.004 < 1.0 - factor < 0.005


# --------------------------------------------------------------------------- #
# the route and the solid
# --------------------------------------------------------------------------- #
def test_route_length_error_is_directional_and_de_biases_to_the_geodesic() -> None:
    channel = ice.ice_channel()
    plane_total = 0.0
    corrected = 0.0
    for index in range(len(channel.centre) - 1):
        start, end = channel.centre[index], channel.centre[index + 1]
        plane_km = math.dist(start, end) / 1000.0
        azimuth = plane_azimuth((end[0] - start[0], end[1] - start[1]))
        plane_total += plane_km
        corrected += plane_km * ice.plane_to_true_factor(start[0], start[1], azimuth)
    # the plane over-states the route, because most legs have a tangential component
    assert plane_total > channel.total_km
    assert plane_total / channel.total_km - 1.0 == pytest.approx(2.33e-4, rel=0.05)
    # and the published factor removes it to second order
    assert corrected == pytest.approx(channel.total_km, rel=1e-6)


def test_corridor_volume_matches_the_analytic_prism() -> None:
    channel = ice.ice_channel()
    measurement = ice.channel_measurements(ice.build_channel_mesh(channel), channel)
    assert measurement["watertight"] is True
    assert measurement["edge_manifold"] is True
    assert measurement["euler_characteristic"] == 2
    assert measurement["self_intersection_pairs"] == 0
    assert measurement["volume_relative_error"] == pytest.approx(0.0, abs=1e-9)


def test_true_widths_bracket_the_nominal_width() -> None:
    channel = ice.ice_channel()
    widths = ice.channel_true_widths(channel)
    assert widths["plane_width_m"] == channel.width
    # narrowing from the projection is bounded by the tangential scale, and the
    # mitre is the only thing allowed to go above the nominal width
    assert widths["true_width_m"]["minimum"] > channel.width * 0.995
    assert widths["true_width_m"]["maximum"] < channel.width * 1.02
    assert widths["true_to_plane_ratio"]["maximum"] == pytest.approx(1.0, abs=1e-15)
    assert len(widths["per_station_true_width_m"]) == len(channel.centre)


# --------------------------------------------------------------------------- #
# the repeated-waypoint guard, with its negative control
# --------------------------------------------------------------------------- #
REPEATED = (
    ("start", 70.0, 170.0),
    ("north", 71.0, 171.0),
    ("north again", 71.0, 171.0),
    ("end", 72.0, 172.0),
)


def test_repeated_waypoint_is_refused_by_name() -> None:
    with pytest.raises(ValueError, match="repeated waypoint") as caught:
        ice.ice_channel(REPEATED)
    message = str(caught.value)
    assert "'north'" in message and "'north again'" in message


def test_negative_control_unguarded_duplicate_folds_the_solid(monkeypatch) -> None:
    """Without the guard the same input closes topologically and still self-intersects."""
    monkeypatch.setattr(ice, "MIN_SEGMENT_KM", 0.0)      # no repeated-waypoint floor
    monkeypatch.setattr(ice, "mitre_excursion_m", lambda *_: -1.0)  # no fold guard
    channel = ice.ice_channel(REPEATED)
    measurement = ice.channel_measurements(ice.build_channel_mesh(channel), channel)
    assert measurement["watertight"] is True
    assert measurement["edge_manifold"] is True
    assert measurement["euler_characteristic"] == 2
    assert measurement["self_intersection_pairs"] > 0


def _corner_route(turn_deg: float, short_leg_m: float):
    """200 km north, a turn, a short leg, then the turn repeated and 50 km on."""
    start = ("a", 70.0, 170.0)
    corner = ("b",) + direct(70.0, 170.0, 0.0, 200.0)
    kink = ("c",) + direct(corner[1], corner[2], turn_deg, short_leg_m / 1000.0)
    end = ("d",) + direct(kink[1], kink[2], turn_deg, 50.0)
    return start, corner, kink, end


@pytest.mark.parametrize("turn_deg,short_leg_m", [(20.0, 5.0), (45.0, 40.0), (90.0, 100.0)])
def test_leg_shorter_than_the_mitre_excursion_is_refused(turn_deg, short_leg_m) -> None:
    """A leg shorter than half_width*tan(turn/2) folds the inner wall: refuse it."""
    excursion = ice.mitre_excursion_m(turn_deg, 125.0)
    assert short_leg_m <= excursion  # the case really is inside the fold region
    with pytest.raises(ValueError, match="fold back over the previous strip"):
        ice.ice_channel(_corner_route(turn_deg, short_leg_m))


@pytest.mark.parametrize("turn_deg", [20.0, 45.0, 90.0])
def test_fold_threshold_is_the_mitre_excursion(turn_deg: float) -> None:
    """The guard's criterion is not a guess: just above it the plan is simple."""
    excursion = ice.mitre_excursion_m(turn_deg, 125.0)
    channel = ice.ice_channel(_corner_route(turn_deg, excursion * 1.2))
    measurement = ice.channel_measurements(ice.build_channel_mesh(channel), channel)
    assert measurement["plan_is_simple"] is True
    assert measurement["self_intersection_pairs"] == 0
    assert measurement["watertight"] is True


def test_mitre_excursion_is_zero_on_a_straight_leg_and_capped_at_the_miter_cap() -> None:
    assert ice.mitre_excursion_m(0.0, 125.0) == pytest.approx(0.0, abs=1e-12)
    assert ice.mitre_excursion_m(0.0, 125.0) == 0.0
    capped = ice.mitre_excursion_m(180.0, 125.0)
    assert capped == pytest.approx(125.0 * math.tan(math.radians(60.0)), rel=1e-12)
    with pytest.raises(ValueError):
        ice.mitre_excursion_m(45.0, 0.0)


def test_guard_does_not_disturb_the_public_route() -> None:
    channel = ice.ice_channel()
    assert channel.total_km == pytest.approx(1061.7912, rel=1e-6)
    assert meshlib.is_watertight(ice.build_channel_mesh(channel)) is True