"""Polyline geometry of the de Long -> Pevek ice channel.

The ice-going leg of the campaign needs three things that a hull generator cannot
give: the great-circle length of the route, a corridor (a ribbon of stated width)
around the track, and a closed solid of that corridor so that a DEM ice-field or
a FEM ice-belt model can be built without re-deriving any of it.

The route is a *campaign-scale* navigation polyline in the East Siberian Sea,
given in decimal degrees (WGS84, rounded to 0.1 deg). It is projected onto a local
azimuthal-equidistant plane centred on the first waypoint: along-track and
across-track distances are then true great-circle distances, which a plain
equirectangular projection would not give at 70-77 deg N.

That claim holds for the *radial* direction only. An azimuthal-equidistant
projection keeps every distance from the centre exact but inflates the tangential
scale by ``c/sin c`` in the central angle ``c`` from the origin, so a plane length
is a true length times ``sqrt(cos^2 beta + (sin beta / k)^2)``, ``k = c/sin c``,
``beta`` the angle to the radial direction (:func:`plane_to_true_factor`). At
campaign scale that is 0.46 % at 1061 km, which is why the corridor width is also
reported in true metres (:func:`channel_true_widths`) instead of being read off
the plane number.

Only the Python standard library is required.
"""
from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

from .geometry import Mesh, MeshBuilder

EARTH_RADIUS_KM = 6371.0088
UNITS = "m"
CHANNEL_FRAME = (
    "local azimuthal-equidistant plane centred on the first waypoint: "
    "x=east [m], y=north [m], z=up [m], z=0 at the water/ice surface"
)

#: Coarse campaign polyline, north-west (de Long Islands) to south-east (Pevek).
DE_LONG_PEVEK_ROUTE: tuple[tuple[str, float, float], ...] = (
    ("de Long (Bennett Island)", 76.70, 149.00),
    ("East Siberian Sea north", 75.80, 155.00),
    ("Medvezhyi Islands north", 74.60, 161.00),
    ("Kolyma drift", 73.30, 166.00),
    ("Chaun approach", 72.00, 169.50),
    ("Cape Billings south", 70.80, 170.50),
    ("Pevek roadstead", 70.10, 170.60),
    ("Pevek port", 69.70, 170.30),
)

#: Default corridor: 250 m of broken channel for a 24 m beam icebreaker escort,
#: extending to the hull's design draft so the walls reach the ice keel line.
DEFAULT_WIDTH_M = 250.0
DEFAULT_DEPTH_M = 11.0
MAX_MITER_FACTOR = 2.0

#: A leg shorter than this is a repeated waypoint: the miter normal is undefined
#: (the two legs are parallel) and the zero-length strip folds the wall onto
#: itself, so the solid is refused rather than built.
MIN_SEGMENT_KM = 1.0e-3

#: Half-angle at which the miter cap binds, ``acos(1/MAX_MITER_FACTOR)`` in degrees.
#: Beyond it the wall vertex sits at ``MAX_MITER_FACTOR`` half-widths instead of on
#: the bisector, so the excursion it demands of the legs stops growing.
MITRE_CAP_HALF_ANGLE_DEG = math.degrees(math.acos(1.0 / MAX_MITER_FACTOR))


def great_circle_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Haversine distance on a sphere of radius :data:`EARTH_RADIUS_KM`."""
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dphi = phi2 - phi1
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2.0) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(dlambda / 2.0) ** 2
    return 2.0 * EARTH_RADIUS_KM * math.asin(min(1.0, math.sqrt(a)))


def initial_bearing_deg(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Initial great-circle bearing from point 1 to point 2, degrees from north."""
    phi1, phi2 = math.radians(lat1), math.radians(lat2)
    dlambda = math.radians(lon2 - lon1)
    y = math.sin(dlambda) * math.cos(phi2)
    x = math.cos(phi1) * math.sin(phi2) - math.sin(phi1) * math.cos(phi2) * math.cos(dlambda)
    return math.degrees(math.atan2(y, x)) % 360.0


def _project(lat: float, lon: float, origin: tuple[float, float]) -> tuple[float, float]:
    """Azimuthal-equidistant projection to metres about ``origin`` = (lat, lon)."""
    distance = great_circle_km(origin[0], origin[1], lat, lon) * 1000.0
    bearing = math.radians(initial_bearing_deg(origin[0], origin[1], lat, lon))
    return (distance * math.sin(bearing), distance * math.cos(bearing))


def _unit(vector: tuple[float, float]) -> tuple[float, float]:
    length = math.hypot(vector[0], vector[1])
    if length == 0.0:
        return (1.0, 0.0)
    return (vector[0] / length, vector[1] / length)


def _azimuth_deg(vector: tuple[float, float]) -> float:
    """Compass azimuth of a plane vector: degrees clockwise from +y (north)."""
    return math.degrees(math.atan2(vector[0], vector[1])) % 360.0


def mitre_excursion_m(turn_deg: float, half_width: float) -> float:
    """How far the mitred wall vertex reaches along a leg, in metres.

    A turn of ``turn_deg`` puts the wall vertex ``half_width * tan(turn/2)`` past
    the centreline station, measured along the leg. A leg shorter than that has
    its inner wall laid back over the previous strip, so the corridor folds
    through itself while staying watertight -- the reason :func:`ice_channel`
    refuses such a route. The miter cap bounds the term at
    :data:`MITRE_CAP_HALF_ANGLE_DEG`.
    """
    if not math.isfinite(turn_deg) or not math.isfinite(half_width):
        raise ValueError("turn_deg and half_width must be finite")
    if half_width <= 0.0:
        raise ValueError(f"half_width must be positive, got {half_width!r}")
    half_turn = min(abs(turn_deg) / 2.0, MITRE_CAP_HALF_ANGLE_DEG)
    return half_width * math.tan(math.radians(half_turn))


def plane_to_true_factor(x: float, y: float, plane_azimuth_deg: float) -> float:
    """True length per unit plane length, for a plane vector at ``plane_azimuth_deg``.

    ``x``, ``y`` are metres in the local plane of :data:`CHANNEL_FRAME`. The
    projection is azimuthal-equidistant about the first waypoint, so its metric is
    diagonal: radial scale exactly 1, tangential scale ``k = c/sin c`` in the
    central angle ``c = |(x, y)| / R``. A plane vector of length ``L`` at angle
    ``beta`` to the radial direction therefore has true length

        ``L * sqrt(cos^2(beta) + (sin(beta) / k)^2)``,

    which is the factor returned here. It is first order in ``L/R``: over the
    250 m corridor width the neglected term is second order, ``1.5e-9``. At the
    origin the factor is exactly 1, because every direction there is radial.
    """
    for label, value in (("x", x), ("y", y), ("plane_azimuth_deg", plane_azimuth_deg)):
        if not math.isfinite(value):
            raise ValueError(f"{label} must be finite, got {value!r}")
    radial = math.hypot(x, y)
    central = radial / (EARTH_RADIUS_KM * 1000.0)
    if central <= 0.0:
        return 1.0
    tangent_scale = central / math.sin(central)  # >= 1: plane / true, tangential
    beta = math.radians(plane_azimuth_deg - _azimuth_deg((x, y)))
    return math.hypot(math.cos(beta), math.sin(beta) / tangent_scale)


@dataclass(frozen=True)
class IceChannel:
    """Corridor geometry: waypoints, local-plane centreline and the two offset walls."""

    names: tuple[str, ...]
    latitude: tuple[float, ...]
    longitude: tuple[float, ...]
    centre: tuple[tuple[float, float], ...]
    left: tuple[tuple[float, float], ...]
    right: tuple[tuple[float, float], ...]
    segment_km: tuple[float, ...]
    cumulative_km: tuple[float, ...]
    width: float
    depth: float
    origin: tuple[float, float]

    @property
    def total_km(self) -> float:
        return self.cumulative_km[-1]

    def to_dict(self) -> dict[str, Any]:
        return {
            "route": "de Long (Bennett Island) -> Pevek",
            "frame": CHANNEL_FRAME,
            "units": UNITS,
            "waypoints": [
                {"name": name, "latitude_deg": lat, "longitude_deg": lon}
                for name, lat, lon in zip(self.names, self.latitude, self.longitude, strict=True)
            ],
            "centreline_local_m": [{"x": x, "y": y} for x, y in self.centre],
            "left_wall_local_m": [{"x": x, "y": y} for x, y in self.left],
            "right_wall_local_m": [{"x": x, "y": y} for x, y in self.right],
            "segment_great_circle_km": list(self.segment_km),
            "cumulative_km": list(self.cumulative_km),
            "total_km": self.total_km,
            "total_nautical_miles": self.total_km / 1.852,
            "corridor_width_m": self.width,
            "corridor_depth_m": self.depth,
            "projection_origin": {"latitude_deg": self.origin[0], "longitude_deg": self.origin[1]},
            "segment_heading_deg": [
                initial_bearing_deg(self.latitude[i], self.longitude[i],
                                    self.latitude[i + 1], self.longitude[i + 1])
                for i in range(len(self.latitude) - 1)
            ],
        }


def ice_channel(
    route: Sequence[tuple[str, float, float]] = DE_LONG_PEVEK_ROUTE,
    *,
    width: float = DEFAULT_WIDTH_M,
    depth: float = DEFAULT_DEPTH_M,
) -> IceChannel:
    """Build the corridor of ``route``: centreline, mitered walls, distances.

    Offsets are perpendicular to the local bisector, mitred so the wall stays
    exactly ``width/2`` from the centreline on straight runs; the miter length is
    capped at :data:`MAX_MITER_FACTOR` half-widths so a sharp waypoint cannot
    throw a spike. The cap means a very sharp turn clips the corner rather than
    preserving the width there, which is stated rather than hidden.
    """
    if len(route) < 2:
        raise ValueError("an ice channel needs at least two waypoints")
    if not math.isfinite(width) or width <= 0.0:
        raise ValueError(f"width must be a finite positive number, got {width!r}")
    if not math.isfinite(depth) or depth <= 0.0:
        raise ValueError(f"depth must be a finite positive number, got {depth!r}")
    for name, lat, lon in route:
        if not (-90.0 <= lat <= 90.0) or not math.isfinite(lon):
            raise ValueError(f"waypoint {name!r} has an out-of-range position ({lat}, {lon})")
    for index in range(len(route) - 1):
        leg = great_circle_km(route[index][1], route[index][2],
                              route[index + 1][1], route[index + 1][2])
        if leg < MIN_SEGMENT_KM:
            raise ValueError(
                f"waypoints {route[index][0]!r} and {route[index + 1][0]!r} are "
                f"{leg * 1000.0:.3f} m apart (minimum {MIN_SEGMENT_KM * 1000.0:.0f} m): a "
                "repeated waypoint has no tangent, so the corridor folds through itself "
                "instead of closing"
            )

    origin = (route[0][1], route[0][2])
    centre = tuple(_project(lat, lon, origin) for _, lat, lon in route)
    segment = tuple(
        great_circle_km(route[i][1], route[i][2], route[i + 1][1], route[i + 1][2])
        for i in range(len(route) - 1)
    )
    cumulative = [0.0]
    for length in segment:
        cumulative.append(cumulative[-1] + length)

    half = width / 2.0
    left: list[tuple[float, float]] = []
    right: list[tuple[float, float]] = []
    for index in range(len(centre)):
        if index == 0:
            tangent = _unit((centre[1][0] - centre[0][0], centre[1][1] - centre[0][1]))
            normal = (-tangent[1], tangent[0])
            scale = 1.0
        elif index == len(centre) - 1:
            tangent = _unit(
                (centre[-1][0] - centre[-2][0], centre[-1][1] - centre[-2][1])
            )
            normal = (-tangent[1], tangent[0])
            scale = 1.0
        else:
            incoming = _unit(
                (centre[index][0] - centre[index - 1][0], centre[index][1] - centre[index - 1][1])
            )
            outgoing = _unit(
                (centre[index + 1][0] - centre[index][0], centre[index + 1][1] - centre[index][1])
            )
            into = (-incoming[1], incoming[0])
            out = (-outgoing[1], outgoing[0])
            normal = _unit((into[0] + out[0], into[1] + out[1]))
            cosine = normal[0] * into[0] + normal[1] * into[1]
            scale = 1.0 / cosine if cosine > 1.0 / MAX_MITER_FACTOR else MAX_MITER_FACTOR
        left.append(
            (centre[index][0] + half * scale * normal[0], centre[index][1] + half * scale * normal[1])
        )
        right.append(
            (centre[index][0] - half * scale * normal[0], centre[index][1] - half * scale * normal[1])
        )
    for index in range(1, len(centre) - 1):
        incoming = _unit((centre[index][0] - centre[index - 1][0],
                          centre[index][1] - centre[index - 1][1]))
        outgoing = _unit((centre[index + 1][0] - centre[index][0],
                          centre[index + 1][1] - centre[index][1]))
        cosine = max(-1.0, min(1.0, incoming[0] * outgoing[0] + incoming[1] * outgoing[1]))
        turn = math.degrees(math.acos(cosine))
        excursion = mitre_excursion_m(turn, half)
        for leg_index in (index - 1, index):
            leg_m = segment[leg_index] * 1000.0
            if leg_m <= excursion:
                raise ValueError(
                    f"the {turn:.1f} deg turn at waypoint {route[index][0]!r} mitres "
                    f"{excursion:.1f} m along each leg, but leg {route[leg_index][0]!r} -> "
                    f"{route[leg_index + 1][0]!r} is only {leg_m:.1f} m long: the inner wall "
                    "would fold back over the previous strip (watertight, "
                    "self-intersecting) instead of forming a corridor"
                )

    return IceChannel(
        names=tuple(name for name, _, _ in route),
        latitude=tuple(lat for _, lat, _ in route),
        longitude=tuple(lon for _, _, lon in route),
        centre=centre,
        left=tuple(left),
        right=tuple(right),
        segment_km=segment,
        cumulative_km=tuple(cumulative),
        width=width,
        depth=depth,
        origin=origin,
    )


def build_channel_mesh(channel: IceChannel) -> Mesh:
    """Closed prism of the corridor: top, bottom, two walls, two end caps.

    The corridor is a structured quad strip, so the solid is built directly from
    the four corner columns per station -- no polygon triangulation and no fan,
    which is what keeps it watertight for any (non-self-intersecting) plan shape.
    """
    builder = MeshBuilder(units=UNITS, coordinate_frame=CHANNEL_FRAME)
    top: list[tuple[int, int]] = []
    bottom: list[tuple[int, int]] = []
    for left, right in zip(channel.left, channel.right, strict=True):
        left_top = builder.add_vertex((left[0], left[1], 0.0))
        right_top = builder.add_vertex((right[0], right[1], 0.0))
        left_bottom = builder.add_vertex((left[0], left[1], -channel.depth))
        right_bottom = builder.add_vertex((right[0], right[1], -channel.depth))
        top.append((left_top, right_top))
        bottom.append((left_bottom, right_bottom))

    for index in range(len(top) - 1):
        left_top, right_top = top[index]
        next_left_top, next_right_top = top[index + 1]
        left_bottom, right_bottom = bottom[index]
        next_left_bottom, next_right_bottom = bottom[index + 1]
        builder.add_quad(right_top, next_right_top, next_left_top, left_top)  # +z
        builder.add_quad(right_bottom, left_bottom, next_left_bottom, next_right_bottom)  # -z
        builder.add_quad(left_top, next_left_top, next_left_bottom, left_bottom)  # left wall
        builder.add_quad(right_top, right_bottom, next_right_bottom, next_right_top)  # right wall
    left_top, right_top = top[0]
    left_bottom, right_bottom = bottom[0]
    builder.add_quad(left_top, left_bottom, right_bottom, right_top)  # start cap, -tangent
    last_left_top, last_right_top = top[-1]
    last_left_bottom, last_right_bottom = bottom[-1]
    builder.add_quad(
        last_left_top, last_right_top, last_right_bottom, last_left_bottom
    )  # end cap, +tangent
    return builder.build()


def channel_measurements(mesh: Mesh, channel: IceChannel) -> dict[str, Any]:
    """Topology and plan-shape numbers of the corridor solid, including whether the
    plan is simple (the corridor walls must not cross each other)."""
    from . import mesh as meshlib
    from .intersect import find_self_intersections

    edges = meshlib.edge_face_counts(mesh)
    widths = [
        math.hypot(left[0] - right[0], left[1] - right[1])
        for left, right in zip(channel.left, channel.right, strict=True)
    ]
    length = sum(
        math.hypot(channel.centre[i + 1][0] - channel.centre[i][0],
                   channel.centre[i + 1][1] - channel.centre[i][1])
        for i in range(len(channel.centre) - 1)
    )
    volume = meshlib.total_signed_volume(mesh)
    return {
        "vertices": len(mesh.vertices),
        "faces": len(mesh.faces),
        "edges": len(edges),
        "watertight": meshlib.is_watertight(mesh),
        "edge_manifold": meshlib.is_edge_manifold(mesh),
        "euler_characteristic": len(mesh.vertices) - len(edges) + len(mesh.faces),
        "signed_volume_m3": volume,
        "surface_area_m2": meshlib.surface_area(mesh),
        "corridor_width_m": channel.width,
        "measured_width_m": {
            "minimum": min(widths),
            "maximum": max(widths),
            "mean": sum(widths) / len(widths),
        },
        "plan_length_m": length,
        "plan_length_km": length / 1000.0,
        "great_circle_length_km": channel.total_km,
        "volume_from_plan_m3": length * channel.width * channel.depth,
        "volume_relative_error": (
            volume / (length * channel.width * channel.depth) - 1.0
            if length * channel.width * channel.depth > 0.0
            else float("nan")
        ),
        "self_intersection_pairs": len(find_self_intersections(mesh)),
        "plan_is_simple": not find_self_intersections(mesh),
    }


def channel_true_widths(channel: IceChannel) -> dict[str, Any]:
    """Corridor width re-expressed in true (great-circle) metres, per station.

    The walls are offset ``width/2`` from the centreline *in the projection
    plane*, so ``width`` is a plane length. Each station's width vector is
    converted with :func:`plane_to_true_factor`; that vector is exactly what
    :func:`channel_measurements` reports as ``measured_width_m``, so the two sets
    of numbers are directly comparable. The maximum exceeds ``width`` by the mitre
    factor at a turn -- a widened station chord, not projection error.
    """
    factors: list[float] = []
    true_widths: list[float] = []
    for centre, left in zip(channel.centre, channel.left, strict=True):
        offset = (left[0] - centre[0], left[1] - centre[1])
        factor = plane_to_true_factor(centre[0], centre[1], _azimuth_deg(offset))
        factors.append(factor)
        true_widths.append(2.0 * math.hypot(*offset) * factor)
    plane_chords = [
        math.hypot(left[0] - right[0], left[1] - right[1])
        for left, right in zip(channel.left, channel.right, strict=True)
    ]
    return {
        "plane_width_m": channel.width,
        "plane_chord_m": {
            "minimum": min(plane_chords),
            "maximum": max(plane_chords),
            "mean": sum(plane_chords) / len(plane_chords),
        },
        "true_width_m": {
            "minimum": min(true_widths),
            "maximum": max(true_widths),
            "mean": sum(true_widths) / len(true_widths),
        },
        "true_to_plane_ratio": {
            "minimum": min(factors),
            "maximum": max(factors),
            "mean": sum(factors) / len(factors),
        },
        "per_station_true_width_m": true_widths,
    }


__all__ = [
    "CHANNEL_FRAME",
    "DEFAULT_DEPTH_M",
    "DEFAULT_WIDTH_M",
    "DE_LONG_PEVEK_ROUTE",
    "EARTH_RADIUS_KM",
    "IceChannel",
    "MIN_SEGMENT_KM",
    "build_channel_mesh",
    "channel_measurements",
    "channel_true_widths",
    "great_circle_km",
    "ice_channel",
    "MITRE_CAP_HALF_ANGLE_DEG",
    "initial_bearing_deg",
    "mitre_excursion_m",
    "plane_to_true_factor",
]
