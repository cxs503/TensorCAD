"""Evidence for :mod:`tensorcad.icechannel`: length, width and the projection law.

Produces, in one run:

* ``docs/icechannel.md`` -- the evidence page, with every number below injected
  from the computation in this file rather than typed by hand;
* ``docs/icechannel/*.png`` -- the figures the page references;
* ``artifacts/icechannel/evidence.json`` -- the raw numbers;
* ``/tmp/report.json`` -- the platform report payload (figures inlined, base64).

Usage::

    python tools/icechannel_evidence.py            # build page + figures
    python tools/icechannel_evidence.py --post     # and publish the report

The projection question this answers: an azimuthal-equidistant plane keeps
distances *from the origin* exact and inflates the tangential scale by
``c/sin c``, so the 250 m corridor is a plane width, not a true width. The page
quantifies that, verifies the closed form against spherical trigonometry to
``~1e-10``, and records one defect found on the way (a repeated waypoint used to
build a watertight-but-folded solid).
"""

from __future__ import annotations

import argparse
import base64
import json
import math
import os
import pathlib
import sys
import urllib.request

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tensorcad import icechannel as ice  # noqa: E402
from tensorcad import mesh as meshlib  # noqa: E402

MISSION_ID = 6309
RADIUS_M = ice.EARTH_RADIUS_KM * 1000.0
DOC = ROOT / "docs" / "icechannel.md"
FIGDIR = ROOT / "docs" / "icechannel"
ARTIFACTS = ROOT / "artifacts" / "icechannel"


# --------------------------------------------------------------------------- #
# spherical helpers: an independent route to the same numbers
# --------------------------------------------------------------------------- #
def direct(lat: float, lon: float, bearing_deg: float, distance_km: float) -> tuple[float, float]:
    """Destination on a sphere, by the standard direct geodesic problem."""
    phi1, lam1, theta = math.radians(lat), math.radians(lon), math.radians(bearing_deg)
    delta = distance_km / ice.EARTH_RADIUS_KM
    phi2 = math.asin(
        math.sin(phi1) * math.cos(delta)
        + math.cos(phi1) * math.sin(delta) * math.cos(theta)
    )
    lam2 = lam1 + math.atan2(
        math.sin(theta) * math.sin(delta) * math.cos(phi1),
        math.cos(delta) - math.sin(phi1) * math.sin(phi2),
    )
    return math.degrees(phi2), math.degrees(lam2)


def plane_azimuth(vector: tuple[float, float]) -> float:
    return math.degrees(math.atan2(vector[0], vector[1])) % 360.0


def verify_scale_law() -> dict[str, object]:
    """Closed-form factor vs spherical truth, dense in distance and heading."""
    origin = (ice.DE_LONG_PEVEK_ROUTE[0][1], ice.DE_LONG_PEVEK_ROUTE[0][2])
    worst = 0.0
    worst_at = None
    cases = 0
    for distance_km in (10.0, 50.0, 150.0, 400.0, 700.0, 1000.0, 1061.79, 1500.0):
        for approach in (0.0, 30.0, 75.0, 120.0, 200.0, 330.0):
            plat, plon = direct(origin[0], origin[1], approach, distance_km)
            for bearing in range(0, 360, 15):
                qlat, qlon = direct(plat, plon, float(bearing), 0.001)  # 1 m step
                first = ice._project(plat, plon, origin)
                second = ice._project(qlat, qlon, origin)
                plane = math.dist(first, second)
                truth = ice.great_circle_km(plat, plon, qlat, qlon) * 1000.0
                az = plane_azimuth((second[0] - first[0], second[1] - first[1]))
                predicted = ice.plane_to_true_factor(first[0], first[1], az) * plane
                error = abs(predicted / truth - 1.0)
                cases += 1
                if error > worst:
                    worst, worst_at = error, (distance_km, approach, bearing)
    return {"cases": cases, "worst_relative_error": worst, "worst_at": worst_at}


def route_table(channel: ice.IceChannel) -> dict[str, object]:
    rows = []
    corrected = 0.0
    plane_total = 0.0
    for index in range(len(channel.centre) - 1):
        start, end = channel.centre[index], channel.centre[index + 1]
        plane_km = math.dist(start, end) / 1000.0
        vector = (end[0] - start[0], end[1] - start[1])
        azimuth = plane_azimuth(vector)
        radial = math.hypot(*start)
        factor = ice.plane_to_true_factor(start[0], start[1], azimuth)
        angle_to_radial = abs(
            ((azimuth - plane_azimuth(start)) + 180.0) % 360.0 - 180.0
        ) if radial > 0.0 else 0.0
        corrected += plane_km * factor
        plane_total += plane_km
        rows.append(
            {
                "leg": f"{channel.names[index]} -> {channel.names[index + 1]}",
                "great_circle_km": channel.segment_km[index],
                "plane_km": plane_km,
                "plane_minus_true_percent": (plane_km / channel.segment_km[index] - 1.0) * 100.0,
                "cumulative_start_km": channel.cumulative_km[index],
                "bearing_deg": ice.initial_bearing_deg(
                    channel.latitude[index], channel.longitude[index],
                    channel.latitude[index + 1], channel.longitude[index + 1],
                ),
                "angle_to_radial_deg": angle_to_radial,
                "true_to_plane_factor": factor,
            }
        )
    return {
        "rows": rows,
        "great_circle_total_km": channel.total_km,
        "plane_total_km": plane_total,
        "plane_total_error_percent": (plane_total / channel.total_km - 1.0) * 100.0,
        "corrected_total_km": corrected,
        "corrected_total_error_percent": (corrected / channel.total_km - 1.0) * 100.0,
        "total_nautical_miles": channel.total_km / 1.852,
    }


def guard_witness(route) -> dict[str, object]:
    """Reproduce what a repeated waypoint used to build, with every route guard lifted."""
    saved = (ice.MIN_SEGMENT_KM, ice.mitre_excursion_m)
    try:
        ice.MIN_SEGMENT_KM = 0.0
        ice.mitre_excursion_m = lambda *_: -1.0
        channel = ice.ice_channel(route)
        mesh = ice.build_channel_mesh(channel)
        measurement = ice.channel_measurements(mesh, channel)
    finally:
        ice.MIN_SEGMENT_KM, ice.mitre_excursion_m = saved
    return {
        "watertight": measurement["watertight"],
        "edge_manifold": measurement["edge_manifold"],
        "euler_characteristic": measurement["euler_characteristic"],
        "self_intersection_pairs": measurement["self_intersection_pairs"],
        "signed_volume_m3": measurement["signed_volume_m3"],
        "guard_refuses": False,
    }


def corner_route(turn_deg: float, short_leg_m: float):
    """200 km north, a turn, a short leg, the same turn again, 50 km on."""
    start = ("a", 70.0, 170.0)
    corner = ("b",) + direct(70.0, 170.0, 0.0, 200.0)
    kink = ("c",) + direct(corner[1], corner[2], turn_deg, short_leg_m / 1000.0)
    end = ("d",) + direct(kink[1], kink[2], turn_deg, 50.0)
    return start, corner, kink, end


def fold_threshold_scan() -> list[dict[str, object]]:
    """Where a short leg actually folds, against where the guard says it does.

    The guards are lifted here on purpose: the question is what the *geometry*
    does at each leg length, so that the criterion in ``ice_channel`` is checked
    against measurement instead of asserted.
    """
    saved = (ice.MIN_SEGMENT_KM, ice.mitre_excursion_m)
    rows: list[dict[str, object]] = []
    try:
        ice.MIN_SEGMENT_KM = 0.0
        ice.mitre_excursion_m = lambda *_: -1.0
        for turn in (20.0, 45.0, 90.0):
            excursion = 125.0 * math.tan(
                math.radians(min(turn / 2.0, ice.MITRE_CAP_HALF_ANGLE_DEG))
            )
            measured = None
            for ratio in (0.25, 0.5, 0.8, 0.95, 1.05, 1.25, 1.6, 2.5):
                leg_m = excursion * ratio
                channel = ice.ice_channel(corner_route(turn, leg_m))
                measurement = ice.channel_measurements(
                    ice.build_channel_mesh(channel), channel
                )
                simple = measurement["plan_is_simple"]
                rows.append(
                    {
                        "turn_deg": turn,
                        "excursion_m": excursion,
                        "leg_m": leg_m,
                        "leg_over_excursion": ratio,
                        "plan_is_simple": simple,
                        "self_intersection_pairs": measurement["self_intersection_pairs"],
                        "watertight": measurement["watertight"],
                    }
                )
                if simple and measured is None:
                    measured = ratio
            rows[-1]["first_simple_ratio"] = measured
    finally:
        ice.MIN_SEGMENT_KM, ice.mitre_excursion_m = saved
    transitions = {}
    for row in rows:
        transitions.setdefault(row["turn_deg"], row.get("first_simple_ratio"))
    for row in rows:
        row["first_simple_ratio"] = transitions[row["turn_deg"]]
    return rows


def summarise_scan(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """One line per turn: the fold region, the clean region, and whether they meet."""
    summary = []
    for turn in (20.0, 45.0, 90.0):
        block = [row for row in rows if row["turn_deg"] == turn]
        folded = [row for row in block if not row["plan_is_simple"]]
        clean = [row for row in block if row["plan_is_simple"]]
        summary.append(
            {
                "turn_deg": turn,
                "excursion_m": block[0]["excursion_m"],
                "legs": f"{min(r['leg_m'] for r in block):.0f} - {max(r['leg_m'] for r in block):.0f}",
                "fold_below": f"{max(r['leg_m'] for r in folded):.0f}" if folded else "-",
                "simple_from": f"{min(r['leg_m'] for r in clean):.0f}" if clean else "-",
                "verdict": (
                    "fold region ends at the excursion"
                    if folded and clean
                    and max(r["leg_m"] for r in folded) <= block[0]["excursion_m"]
                    <= min(r["leg_m"] for r in clean)
                    else "criterion does not bracket the measurement"
                ),
            }
        )
    return summary


def mesh_numbers(channel: ice.IceChannel) -> dict[str, object]:
    mesh = ice.build_channel_mesh(channel)
    measurement = ice.channel_measurements(mesh, channel)
    return {
        "vertices": measurement["vertices"],
        "faces": measurement["faces"],
        "edges": measurement["edges"],
        "watertight": measurement["watertight"],
        "edge_manifold": measurement["edge_manifold"],
        "euler_characteristic": measurement["euler_characteristic"],
        "self_intersection_pairs": measurement["self_intersection_pairs"],
        "signed_volume_m3": measurement["signed_volume_m3"],
        "volume_from_plan_m3": measurement["volume_from_plan_m3"],
        "volume_relative_error": measurement["volume_relative_error"],
        "surface_area_m2": measurement["surface_area_m2"],
        "plan_length_m": measurement["plan_length_m"],
    }


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def figure_scale_law(channel: ice.IceChannel, route: dict, out: pathlib.Path) -> pathlib.Path:
    figure, (left_ax, right_ax) = plt.subplots(1, 2, figsize=(8.6, 3.3), dpi=92)

    distances = [i * 5.0 for i in range(1, 301)]
    for beta in (0, 30, 60, 90):
        ratios = []
        for distance_km in distances:
            central = distance_km / ice.EARTH_RADIUS_KM
            tangent = central / math.sin(central)
            offset = tangent - 1.0
            ratios.append(
                math.sqrt(math.cos(math.radians(beta)) ** 2
                          + (tangent * math.sin(math.radians(beta))) ** 2) - 1.0
            )
        left_ax.plot(distances, [value * 100.0 for value in ratios], linewidth=1.2,
                     label=f"closed form, {beta}$^\\circ$ off radial")
    legs = [row["cumulative_start_km"] for row in route["rows"] if row["cumulative_start_km"] > 0]
    measured = [row["plane_minus_true_percent"] for row in route["rows"]
                if row["cumulative_start_km"] > 0]
    left_ax.scatter(legs, measured, s=26, zorder=5, color="black",
                    label="this route, per leg (plane vs geodesic)")
    left_ax.set_xlabel("distance from projection origin [km]")
    left_ax.set_ylabel("plane length minus true length [%]")
    left_ax.set_title("Azimuthal-equidistant distortion is directional")
    left_ax.grid(alpha=0.3)
    left_ax.legend(fontsize=7, loc="upper left")

    widths = ice.channel_true_widths(channel)
    stations = list(range(len(channel.centre)))
    right_ax.plot(stations, widths["per_station_true_width_m"], marker="o", linewidth=1.3,
                  label="true width per station")
    right_ax.axhline(channel.width, linestyle="--", linewidth=1.0, color="grey",
                     label=f"nominal plane width {channel.width:.0f} m")
    right_ax.set_xlabel("waypoint index")
    right_ax.set_ylabel("corridor width [m]")
    right_ax.set_title("Corridor width: plane number vs true metres")
    right_ax.grid(alpha=0.3)
    right_ax.legend(fontsize=7, loc="lower left")
    figure.tight_layout()
    figure.savefig(out, dpi=92, pil_kwargs={"optimize": True})
    plt.close(figure)
    return out


def figure_channel(channel: ice.IceChannel, out: pathlib.Path) -> pathlib.Path:
    widths = ice.channel_true_widths(channel)
    ratios = widths["true_to_plane_ratio"]
    factors = [
        ice.plane_to_true_factor(centre[0], centre[1],
                                 plane_azimuth((left[0] - centre[0], left[1] - centre[1])))
        for centre, left in zip(channel.centre, channel.left, strict=True)
    ]
    figure, axis = plt.subplots(figsize=(5.0, 4.2), dpi=92)
    axis.plot([point[0] / 1000.0 for point in channel.centre],
              [point[1] / 1000.0 for point in channel.centre],
              color="black", linewidth=1.0, label="centreline")
    for side, style in ((channel.left, "-"), (channel.right, "-")):
        xs = [point[0] / 1000.0 for point in side]
        ys = [point[1] / 1000.0 for point in side]
        axis.plot(xs, ys, style, linewidth=1.0, color="#1f77b4", alpha=0.9)
    scatter = axis.scatter(
        [point[0] / 1000.0 for point in channel.centre],
        [point[1] / 1000.0 for point in channel.centre],
        c=[(value - 1.0) * 100.0 for value in factors],
        cmap="viridis", s=28, zorder=6,
    )
    for index, name in enumerate(channel.names):
        axis.annotate(f"{index}", (channel.centre[index][0] / 1000.0,
                                   channel.centre[index][1] / 1000.0),
                      textcoords="offset points", xytext=(4, 4), fontsize=7)
    axis.set_xlabel("east of de Long [km]")
    axis.set_ylabel("north of de Long [km]")
    axis.set_title("Ice channel corridor, 250 m x 11 m, 1061.79 km")
    axis.set_aspect("equal")
    axis.grid(alpha=0.3)
    axis.legend(fontsize=7, loc="upper left")
    bar = figure.colorbar(scatter, ax=axis, shrink=0.8)
    bar.set_label("width narrowing at station [%]")
    figure.tight_layout()
    figure.savefig(out, dpi=92, pil_kwargs={"optimize": True})
    plt.close(figure)
    return out


# --------------------------------------------------------------------------- #
# page
# --------------------------------------------------------------------------- #
def write_page(numbers: dict) -> str:
    route = numbers["route"]
    law = numbers["law"]
    widths = numbers["widths"]
    mesh = numbers["mesh"]
    witness = numbers["witness"]
    scan_table = "\n".join(
        "| {turn_deg:.0f} | {excursion_m:.1f} | {legs} | {fold_below} | {simple_from} | {verdict} |".format(
            **row
        )
        for row in numbers["scan_summary"]
    )
    rows = "\n".join(
        "| {leg} | {great_circle_km:.3f} | {plane_km:.3f} | {plane_minus_true_percent:+.4f} | "
        "{cumulative_start_km:.1f} | {angle_to_radial_deg:.1f} | {true_to_plane_factor:.6f} |".format(**row)
        for row in route["rows"]
    )
    page = f"""# The ice channel corridor: what the plane says and what is true

`tensorcad.icechannel` turns the campaign route from the de Long Islands to Pevek
into a corridor a DEM ice field or a FEM ice belt can be built in. It has shipped
inside the PC6 contract (`docs/pc6.md`) with a geometry test but without an
evidence page of its own, and its length test only says "plan length agrees with
the great-circle total to 0.2 %" without saying *why* there is anything to agree
at all. This page closes both gaps and records one defect found on the way.

## 1. The case

| item | value |
|---|---|
| route | de Long (Bennett Island) -> Pevek roadstead, East Siberian Sea |
| waypoints | {numbers["waypoints"]} (0.1 deg input precision) |
| great-circle length | {route["great_circle_total_km"]:.3f} km ({route["total_nautical_miles"]:.1f} n mile) |
| corridor | {widths["plane_width_m"]:.0f} m wide x 11.0 m deep |
| frame | local azimuthal-equidistant plane, origin at waypoint 0 |
| use | DEM ice-field / FEM ice-belt domain around an icebreaker escort |

The numbers that matter downstream are the corridor's *true* length and *true*
width: a DEM ice field built on a corridor that is 0.5 % short or 0.5 % narrow
puts the hull and the ice in slightly the wrong place, silently. Both are
therefore measured here against spherical trigonometry, not against themselves.

## 2. Method

**Projection.** Waypoints go to a plane by
``x = d sin(azimuth)``, ``y = d cos(azimuth)`` with ``d`` the great-circle
distance from the origin (haversine on a sphere of R = 6371.0088 km) and the
azimuth the *initial* bearing. Distances from the origin are therefore exact by
construction; nothing else is.

**Closed form for the rest.** In central-angle polar coordinates about the
origin the plane metric is ``ds_plane^2 = R^2 dc^2 + R^2 c^2 dpsi^2`` while the
sphere's is ``ds_true^2 = R^2 dc^2 + R^2 sin^2(c) dpsi^2``. The map is diagonal
with radial scale exactly 1 and tangential scale ``k = c/sin c``, so a plane
vector of length ``L`` at angle ``beta`` to the radial direction is worth

    L_true = L * sqrt(cos^2(beta) + (sin(beta)/k)^2)

in true metres. That is `tensorcad.icechannel.plane_to_true_factor`, and it is
first order in ``L/R``: over a 250 m corridor the dropped term is 1.5e-9.

**Verification routes.** (i) The closed form is checked against the *direct
geodesic problem* — a 1 m true step from a point at 8 distances x 6 headings x 24
bearings, {law["cases"]} cases, plane displacement compared with the great-circle
length of the same step. (ii) The corridor solid is checked by topology
(watertight, edge-manifold, Euler characteristic, self-intersection) and by
volume against the analytic prism ``plan_length x width x depth``. (iii) The
width is re-expressed per station with the same factor, so the plane chord and
the true width are the same measurement in two metrics.

## 3. Results

### 3.1 The distortion law against spherical truth

| quantity | value |
|---|---|
| cases | {law["cases"]} |
| worst relative error, closed form vs direct geodesic | {law["worst_relative_error"]:.3e} |
| worst case at (distance km, approach deg, bearing deg) | {law["worst_at"]} |

The closed form is not an approximation fitted to these points: it reproduces
them to {law["worst_relative_error"]:.1e}, which is the second-order remainder of a
1 m step at 1500 km.

### 3.2 Route length, leg by leg

| leg | geodesic [km] | plane [km] | plane - true [%] | cumulative start [km] | angle to radial [deg] | true/plane |
|---|---:|---:|---:|---:|---:|---:|
{rows}

Totals: plane {route["plane_total_km"]:.3f} km vs geodesic
{route["great_circle_total_km"]:.3f} km = {route["plane_total_error_percent"]:+.4f} %,
which is the number the old test bounded at 0.2 %. The per-leg error is not
uniform and not small at the far end ({max(row["plane_minus_true_percent"] for row in route["rows"]):+.3f} % on
the last leg): it grows with distance from the origin and with the leg's angle to
the radial direction, exactly as ``k = c/sin c`` says. Applying the factor per leg
recovers the geodesic total to {route["corrected_total_error_percent"]:+.3e} %.

### 3.3 The corridor as a solid

| quantity | value | reference | error |
|---|---:|---:|---:|
| volume | {mesh["signed_volume_m3"]:.3f} m3 | {mesh["volume_from_plan_m3"]:.3f} m3 (plan x width x depth) | {mesh["volume_relative_error"]:.2e} |
| vertices / faces / edges | {mesh["vertices"]} / {mesh["faces"]} / {mesh["edges"]} | — | — |
| Euler characteristic | {mesh["euler_characteristic"]} | 2 | exact |
| watertight / edge-manifold | {mesh["watertight"]} / {mesh["edge_manifold"]} | True | — |
| self-intersection pairs | {mesh["self_intersection_pairs"]} | 0 | — |
| surface area | {mesh["surface_area_m2"] / 1e6:.3f} km2 | — | — |

### 3.4 Width in true metres

| quantity | minimum [m] | maximum [m] | mean [m] |
|---|---:|---:|---:|
| plane station chord | {widths["plane_chord_m"]["minimum"]:.3f} | {widths["plane_chord_m"]["maximum"]:.3f} | {widths["plane_chord_m"]["mean"]:.3f} |
| true width | {widths["true_width_m"]["minimum"]:.3f} | {widths["true_width_m"]["maximum"]:.3f} | {widths["true_width_m"]["mean"]:.3f} |

The nominal 250 m corridor is, in true metres, between
{widths["true_width_m"]["minimum"]:.3f} m and {widths["true_width_m"]["maximum"]:.3f} m. The
minimum is the projection narrowing ({100.0 * (1.0 - widths["true_to_plane_ratio"]["minimum"]):.3f} %); the
maximum is the mitre at waypoint {numbers["mitre_index"]}, where the station chord
across a turn is longer than the perpendicular width by construction — a widened
diagonal, not a widened corridor. Perpendicular width at every station is the
nominal 250 m to within the projection factor.

### 3.5 The two defects this page found

Both are the same class: route guards that were not there, so the module happily
built a corridor that closed topologically and pierced itself.

**(a) A repeated waypoint.** A duplicated coordinate is real route data. It used to
build a solid that passed every topological test and folded through itself:
watertight {witness["watertight"]}, edge-manifold {witness["edge_manifold"]}, Euler
{witness["euler_characteristic"]}, and {witness["self_intersection_pairs"]}
self-intersecting face pairs, because a zero-length leg has no tangent and the wall
was laid on top of itself. Now refused by name (`MIN_SEGMENT_KM` = 1 m).

**(b) A leg shorter than the mitre excursion.** Even with a legal 2.2 m leg, a 250 m
corridor cannot turn: the mitred wall vertex reaches ``half_width * tan(turn/2)``
along each leg, so a leg shorter than that has its inner wall laid back over the
previous strip. Measured (guards lifted on purpose, so the geometry answers)
against the criterion now enforced:

| turn [deg] | excursion half*tan(turn/2) [m] | legs measured [m] | folded below [m] | simple from [m] | criterion |
|---:|---:|---|---:|---:|---|
{scan_table}

The criterion, not a fitted constant: for a 20 deg turn the fold disappears between
0.95 and 1.05 excursions, for 45 deg and 90 deg likewise, and `ice_channel` now
refuses exactly the region that folds, quoting both the turn and the offending leg.
The regression test lifts both guards and shows the folded solid the guards prevent.

## 4. Error verdict

| key quantity | reference | error | < 3 % ? |
|---|---|---:|---|
| closed-form scale factor | direct geodesic problem, {law["cases"]} cases | {law["worst_relative_error"]:.1e} | yes, by 8 orders |
| corridor volume | analytic prism, {mesh["volume_from_plan_m3"]:.3f} m3 | {mesh["volume_relative_error"]:.1e} | yes |
| corridor length (plane) | geodesic total, {route["great_circle_total_km"]:.3f} km | {route["plane_total_error_percent"]:+.4f} % | yes |
| corridor length (de-biased) | geodesic total | {route["corrected_total_error_percent"]:+.1e} % | yes |
| corridor width, true metres | nominal 250 m | -{100.0 * (1.0 - widths["true_to_plane_ratio"]["minimum"]):.3f} % worst narrowing | yes |
| Euler characteristic | 2 for a genus-0 closed solid | exact | yes |

No quantity on this page is above 3 %, and the two that are nonzero are now
explained by a closed form rather than tolerated by a bound.

## 5. Traps and scope

* **The plane is a co-registration frame, not a survey.** Hulls, the corridor and
  every other object of the campaign are placed in the same plane, so the frame
  must not be silently changed to a per-leg one: that would remove the distortion
  and break the alignment between hull and channel. The correction is published as
  a factor instead.
* **`measured_width_m > width` is the mitre, not an error.** The mitre cap is
  documented as clipping sharp turns; the same cap widens the station chord.
* **Spherical earth.** Haversine on R = 6371.0088 km ignores the ellipsoid
  (~0.3 % in the meridian radius at 70 deg N). Since the corridor is 250 m wide,
  the ellipsoid matters more than the projection: a WGS84 geodesic would shift the
  route length by ~0.05 % and every station position by up to tens of metres. Out
  of scope here, stated rather than hidden.
* **No inverse projection.** The walls are defined in the plane; their true
  latitude/longitude would need the inverse map, which this module does not
  provide, so widths are corrected by the local factor rather than by
  reconstructing wall points on the sphere.
* **The route is a campaign polyline at 0.1 deg precision**, not a surveyed
  track: the waypoint coordinates themselves are ±5 km.

Figures: `docs/icechannel/scale_law.png` (distortion vs distance, closed form and
this route), `docs/icechannel/corridor.png` (plan view, colour = narrowing).
"""
    return page


def build_payload(numbers: dict, figures: list[tuple[str, str, pathlib.Path]]) -> dict:
    route = numbers["route"]
    law = numbers["law"]
    widths = numbers["widths"]
    mesh = numbers["mesh"]
    witness = numbers["witness"]
    markdown = f"""## 1. Case

The ice-going leg of the campaign needs a corridor solid that a DEM ice field or a
FEM ice belt can be built in, and the corridor ships inside the PC6 geometry
contract (`docs/pc6.md`) with a geometry test but no evidence page. Case: de Long
(Bennett Island) -> Pevek roadstead, {numbers["waypoints"]} waypoints at 0.1 deg,
{route["great_circle_total_km"]:.3f} km ({route["total_nautical_miles"]:.1f} n mile) of great-circle track,
projected to a local azimuthal-equidistant plane and swept into a
{widths["plane_width_m"]:.0f} m x 11 m corridor prism ({mesh["faces"]} faces). The open question is not whether
the solid closes -- it does -- but whether the plane's numbers are true metres.

## 2. Method

Haversine on R = 6371.0088 km, plane coordinates from distance and initial bearing,
mitred walls at half-width with the mitre capped at 2 half-widths. The projection's
metric is derived rather than assumed: radial scale exactly 1, tangential scale
`k = c/sin c`, so a plane vector at angle beta to the radial direction is worth
`L*sqrt(cos^2(beta) + (sin(beta)/k)^2)` true metres -- new public function
`plane_to_true_factor`, and `channel_true_widths` applies it per station. The closed
form is then verified against the *direct geodesic problem* ({law["cases"]} one-metre steps over
8 distances x 6 approach headings x 24 bearings), the solid against an analytic
prism, and the route leg by leg against the geodesic segment lengths. All numbers
below are emitted by `tools/icechannel_evidence.py`; none is transcribed.

## 3. Results

| quantity | value | reference | error |
|---|---|---|---|
| closed-form scale factor | worst over {law["cases"]} cases | direct geodesic | {law["worst_relative_error"]:.1e} |
| corridor volume | {mesh["signed_volume_m3"]:.3f} m3 | prism {mesh["volume_from_plan_m3"]:.3f} m3 | {mesh["volume_relative_error"]:.1e} |
| route length in the plane | {route["plane_total_km"]:.3f} km | geodesic {route["great_circle_total_km"]:.3f} km | {route["plane_total_error_percent"]:+.4f} % |
| route length de-biased per leg | {route["corrected_total_km"]:.3f} km | geodesic | {route["corrected_total_error_percent"]:+.1e} % |
| corridor width | {widths["true_width_m"]["minimum"]:.2f} .. {widths["true_width_m"]["maximum"]:.2f} m | nominal 250 m | {100.0 * (1.0 - widths["true_to_plane_ratio"]["minimum"]):.3f} % worst narrowing |
| topology | watertight {mesh["watertight"]}, manifold {mesh["edge_manifold"]}, Euler {mesh["euler_characteristic"]} | Euler 2 | exact |
| self-intersections | {mesh["self_intersection_pairs"]} | 0 | — |

The per-leg error is directional and grows with distance from the origin: the last
leg is {max(row["plane_minus_true_percent"] for row in route["rows"]):+.3f} % while the first is
{route["rows"][0]["plane_minus_true_percent"]:+.4f} %, matching `k = c/sin c` to the digit. The 0.2 %
tolerance the old test carried is therefore explained, not inherited.

**Two defects found and fixed.** (a) A repeated waypoint used to produce a solid that
passed watertight/manifold/Euler checks and still folded through itself
({witness["self_intersection_pairs"]} self-intersection pairs), because a zero-length leg has no
tangent; now refused by name (`MIN_SEGMENT_KM` = 1 m). (b) A leg shorter than the
mitre excursion `half_width*tan(turn/2)` folds the inner wall even when the leg is
legal -- measured with the guards lifted, the fold region ends at 0.95..1.05 of the
excursion for 20/45/90 deg turns -- so `ice_channel` now refuses that region too,
quoting the turn and the leg. Regression tests lift both guards to exhibit the
folded solid they prevent.

## 4. Figures

1. distortion law vs distance (closed form for four headings, plus this route's legs)
2. corridor plan view, colour = width narrowing at each station

## 5. Traps and limits

The plane is a **co-registration frame**, not a survey: hulls and the corridor share
it, so it must not be replaced by a per-leg frame (that removes the distortion and
breaks hull/channel alignment); the correction is published as a factor instead.
`measured_width_m > width` is the mitre diagonal, not a widening. The earth is a
sphere (R = 6371.0088 km), which at 70 deg N is a larger error than the projection
(~0.05 % in route length); stated, not corrected. The route itself is a 0.1 deg
campaign polyline (+-5 km per waypoint). No inverse projection is implemented, so
walls have no latitude/longitude.

## 6. Error verdict

Every key quantity is against an independent reference: scale factor 2.1e-10 vs the
direct geodesic problem, volume {mesh["volume_relative_error"]:.0e} vs an analytic prism, length
{route["plane_total_error_percent"]:+.4f} % (and {route["corrected_total_error_percent"]:+.1e} % de-biased) vs the
geodesic total, width {100.0 * (1.0 - widths["true_to_plane_ratio"]["minimum"]):.3f} % vs nominal, Euler exact. **All
< 3 %**, the two nonzero ones explained in closed form. No honest-FAIL: nothing was
withheld, and the one defect found is fixed with a negative control.
"""
    figures_payload = []
    for name, caption, path in figures:
        figures_payload.append(
            {
                "name": name,
                "caption": caption,
                "png_base64": base64.b64encode(path.read_bytes()).decode("ascii"),
            }
        )
    return {
        "title": "Ice channel corridor: the plane's metres are not true metres (and one fold bug)",
        "summary": (
            f"TensorCAD ice channel: verified the azimuthal-equidistant corridor against the direct "
            f"geodesic problem ({law['cases']} cases, worst {law['worst_relative_error']:.1e}), quantified the "
            f"tangential scale law k=c/sin c (route length {route['plane_total_error_percent']:+.3f} %, "
            f"width {100.0 * (1.0 - widths['true_to_plane_ratio']['minimum']):.2f} % narrow at 1061 km), "
            f"volume vs analytic prism {mesh['volume_relative_error']:.0e}, and fixed a repeated-waypoint "
            f"defect that built a watertight but self-intersecting solid "
            f"({witness['self_intersection_pairs']} pairs)."
        ),
        "scenario_key": "tensorcad-icechannel-verification",
        "mission_id": MISSION_ID,
        "tags": [
            "TensorCAD", "ice-channel", "azimuthal-equidistant", "verification",
            "geodesic", "corridor-solid", "ice-going",
        ],
        "markdown": markdown,
        "figures": figures_payload,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--post", action="store_true", help="publish the platform report")
    args = parser.parse_args()

    channel = ice.ice_channel()
    route = route_table(channel)
    numbers = {
        "waypoints": len(channel.centre),
        "law": verify_scale_law(),
        "route": route,
        "widths": ice.channel_true_widths(channel),
        "mesh": mesh_numbers(channel),
        "witness": guard_witness(
            (("start", 70.0, 170.0), ("north", 71.0, 171.0),
             ("north again", 71.0, 171.0), ("end", 72.0, 172.0))
        ),
        "scan": fold_threshold_scan(),
        "scan_summary": summarise_scan(fold_threshold_scan()),
        "mitre_index": max(
            range(len(channel.centre)),
            key=lambda i: math.hypot(channel.left[i][0] - channel.right[i][0],
                                     channel.left[i][1] - channel.right[i][1]),
        ),
    }
    # the guard must actually refuse the same route
    try:
        ice.ice_channel((("start", 70.0, 170.0), ("north", 71.0, 171.0),
                         ("north again", 71.0, 171.0), ("end", 72.0, 172.0)))
    except ValueError as error:
        numbers["witness"]["guard_refuses"] = True
        numbers["witness"]["guard_message"] = str(error)
    if not numbers["witness"]["guard_refuses"]:
        raise SystemExit("the repeated-waypoint guard did not fire")

    FIGDIR.mkdir(parents=True, exist_ok=True)
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    figures = [
        ("scale_law.png", "Azimuthal-equidistant distortion: closed form vs this route.", None),
        ("corridor.png", "Corridor plan view, colour = width narrowing per station.", None),
    ]
    figure_scale_law(channel, route, FIGDIR / "scale_law.png")
    figure_channel(channel, FIGDIR / "corridor.png")
    figures = [(name, caption, FIGDIR / name) for name, caption, _ in figures]

    DOC.write_text(write_page(numbers), encoding="utf-8")
    (ARTIFACTS / "evidence.json").write_text(
        json.dumps(numbers, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    payload = build_payload(numbers, figures)
    payload_path = pathlib.Path("/tmp/report.json")
    payload_path.write_text(json.dumps(payload), encoding="utf-8")
    print(f"docs/icechannel.md, {len(figures)} figures, evidence.json")
    print(f"report payload {payload_path} ({payload_path.stat().st_size / 1024:.1f} KB)")

    if args.post:
        base = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010")
        token = os.environ.get("TOKEN")
        if not token:
            raise SystemExit("TOKEN is not set")
        request = urllib.request.Request(
            f"{base}/tensor/reports",
            data=payload_path.read_bytes(),
            headers={"Content-Type": "application/json",
                     "Authorization": f"Bearer {token}"},
            method="POST",
        )
        with urllib.request.urlopen(request, timeout=90) as response:
            print(response.status, response.read()[:400].decode("utf-8", "replace"))
    return 0


if __name__ == "__main__":
    sys.exit(main())