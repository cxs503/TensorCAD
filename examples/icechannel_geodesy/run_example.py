"""Earth-model budget for the ice-channel route: sphere (as used) vs WGS84 (as true).

The campaign corridor is built on a mean-radius sphere. This example measures what that
choice costs, and checks the ellipsoidal geodesic solver it uses to measure it, by two
independent routes (Vincenty's series and a Clairaut quadrature).

Writes ``artifacts/icechannel_geodesy/results.json`` (every number in the report, injected
programmatically) and two PNGs when matplotlib is importable.

Run: ``python examples/icechannel_geodesy/run_example.py``
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_REPO_ROOT = _Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in _sys.path:  # allow running from a checkout without installing
    _sys.path.insert(0, str(_REPO_ROOT))

import json
import math
import pathlib

from tensorcad import geodesy
from tensorcad.icechannel import DE_LONG_PEVEK_ROUTE, ice_channel

OUT = _Path(__file__).resolve().parents[2] / "artifacts" / "icechannel_geodesy"


def route_points():
    return [(name, lat, lon) for name, lat, lon in DE_LONG_PEVEK_ROUTE]


def legs():
    points = route_points()
    return list(zip(points, points[1:]))


def measure():
    rows = []
    total_sphere = total_ellipsoid = 0.0
    worst_route_agreement = 0.0
    for (n1, la1, lo1), (n2, la2, lo2) in legs():
        sphere = geodesy.great_circle(la1, lo1, la2, lo2)
        inverse = geodesy.vincenty_inverse(la1, lo1, la2, lo2)
        quadrature = geodesy.clairaut_distance(la1, lo1, la2, lo2)
        total_sphere += sphere
        total_ellipsoid += inverse.distance_m
        agreement = abs(quadrature / inverse.distance_m - 1.0)
        worst_route_agreement = max(worst_route_agreement, agreement)
        rows.append(
            {
                "leg": f"{n1} -> {n2}",
                "sphere_km": sphere / 1000.0,
                "wgs84_km": inverse.distance_m / 1000.0,
                "sphere_minus_wgs84_m": sphere - inverse.distance_m,
                "relative": sphere / inverse.distance_m - 1.0,
                "clairaut_km": quadrature / 1000.0,
                "route_agreement": agreement,
                "vincenty_iterations": inverse.iterations,
            }
        )

    channel = ice_channel()
    widths = []
    for name, lat, lon in route_points()[1:4]:
        radius_parallel = geodesy.MEAN_RADIUS_M * math.cos(math.radians(lat))
        half_deg = 125.0 / radius_parallel * 180.0 / math.pi
        west = (lat, lon - half_deg)
        east = (lat, lon + half_deg)
        true_width = geodesy.vincenty_inverse(west[0], west[1], east[0], east[1]).distance_m
        sphere_width = geodesy.great_circle(west[0], west[1], east[0], east[1])
        widths.append(
            {
                "station": name,
                "latitude_deg": lat,
                "plane_width_m": 250.0,
                "sphere_width_m": sphere_width,
                "wgs84_width_m": true_width,
                "relative": true_width / sphere_width - 1.0,
            }
        )

    latitude_sweep = []
    for lat in (0.0, 15.0, 30.0, 45.0, 60.0, 65.0, 70.0, 75.0, 77.0, 80.0, 85.0):
        meridional = geodesy.meridional_radius(math.radians(lat)) / geodesy.MEAN_RADIUS_M - 1.0
        transverse = geodesy.prime_vertical_radius(math.radians(lat)) / geodesy.MEAN_RADIUS_M - 1.0
        latitude_sweep.append(
            {"latitude_deg": lat, "meridional_relative": meridional, "transverse_relative": transverse}
        )

    return {
        "route": "de Long (Bennett Island) -> Pevek, East Siberian Sea",
        "waypoints": len(route_points()),
        "ellipsoid": "WGS84 a=6378137.0 m f=1/298.257223563",
        "sphere_radius_m": geodesy.MEAN_RADIUS_M,
        "legs": rows,
        "totals": {
            "sphere_km": total_sphere / 1000.0,
            "wgs84_km": total_ellipsoid / 1000.0,
            "difference_km": (total_sphere - total_ellipsoid) / 1000.0,
            "relative": total_sphere / total_ellipsoid - 1.0,
            "icechannel_total_km": channel.total_km,
            "icechannel_vs_sphere_relative": channel.total_km / (total_sphere / 1000.0) - 1.0,
            "worst_route_agreement": worst_route_agreement,
        },
        "corridor_width": widths,
        "latitude_sweep": latitude_sweep,
        "anchors": {
            "equator_arc_m": geodesy.vincenty_inverse(0.0, 0.0, 0.0, 1.0).distance_m,
            "equator_arc_exact_m": geodesy.WGS84_A * math.pi / 180.0,
            "meridian_quadrant_m": geodesy.meridian_arc_m(0.0, 90.0),
            "meridian_quadrant_published_m": 10001965.7293,
            "equator_circumference_m": 2.0 * math.pi * geodesy.WGS84_A,
        },
    }


def figures(results):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as error:  # pragma: no cover - optional dependency
        return {"matplotlib": f"unavailable: {error}"}

    OUT.mkdir(parents=True, exist_ok=True)
    legs_ = results["legs"]
    labels = [row["leg"].split(" -> ")[0][:18] for row in legs_]

    figure, axis = plt.subplots(figsize=(9.5, 4.6))
    index = range(len(legs_))
    axis.bar([i - 0.2 for i in index], [row["sphere_km"] for row in legs_], width=0.4,
             label="mean-radius sphere (as used)", color="#9ecae1")
    axis.bar([i + 0.2 for i in index], [row["wgs84_km"] for row in legs_], width=0.4,
             label="WGS84 geodesic (true)", color="#08519c")
    axis.set_ylabel("leg length [km]")
    axis.set_xticks(list(index))
    axis.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    twin = axis.twinx()
    twin.plot(list(index), [100.0 * row["relative"] for row in legs_], "o-", color="#cb181d",
              markersize=4, label="sphere / WGS84 - 1 [%]")
    twin.set_ylabel("relative error [%]", color="#cb181d")
    twin.tick_params(axis="y", colors="#cb181d")
    axis.set_title("Ice-channel legs: sphere datum vs WGS84 (%d legs, %.3f km)"
                   % (len(legs_), results["totals"]["wgs84_km"]))
    axis.legend(loc="upper right", fontsize=8)
    twin.legend(loc="lower left", fontsize=8)
    axis.grid(axis="y", alpha=0.3)
    figure.tight_layout()
    figure.savefig(OUT / "route_legs_sphere_vs_wgs84.png", dpi=140)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(7.5, 4.4))
    lats = [row["latitude_deg"] for row in results["latitude_sweep"]]
    axis.plot(lats, [100.0 * row["meridional_relative"] for row in results["latitude_sweep"]],
              "o-", label="meridional  M/R - 1")
    axis.plot(lats, [100.0 * row["transverse_relative"] for row in results["latitude_sweep"]],
              "s-", label="transverse  N/R - 1")
    axis.axhspan(-0.5, 0.0, color="#fee0d2", alpha=0.4)
    axis.axhline(100.0 * results["totals"]["relative"], color="#cb181d", linestyle="--",
                 label="this route (%.3f %%)" % (100.0 * results["totals"]["relative"]))
    axis.axhline(0.46, color="#31a354", linestyle=":", label="projection distortion 0.46 %")
    axis.set_xlabel("latitude [deg]")
    axis.set_ylabel("sphere shortfall against WGS84 [%]")
    axis.set_title("Two error terms of the same size: Earth model vs projection")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(OUT / "earth_model_vs_projection.png", dpi=140)
    plt.close(figure)
    return {"figures": [str(OUT / "route_legs_sphere_vs_wgs84.png"),
                        str(OUT / "earth_model_vs_projection.png")]}


def main() -> int:
    results = measure()
    OUT.mkdir(parents=True, exist_ok=True)
    results["figures"] = figures(results)
    (OUT / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    totals = results["totals"]
    print(f"legs                 : {len(results['legs'])}")
    print(f"sphere total         : {totals['sphere_km']:.3f} km")
    print(f"WGS84 total          : {totals['wgs84_km']:.3f} km")
    print(f"shortfall            : {totals['difference_km']:+.3f} km "
          f"({100.0 * totals['relative']:+.4f} %)")
    print(f"icechannel total_km  : {totals['icechannel_total_km']:.3f} km "
          f"(vs sphere {100.0 * totals['icechannel_vs_sphere_relative']:+.2e} %)")
    print(f"worst route agreement: {totals['worst_route_agreement']:.2e} (Vincenty vs Clairaut)")
    for row in results["corridor_width"]:
        print(f"corridor @ {row['station'][:24]:<24} plane 250.000 m -> WGS84 "
              f"{row['wgs84_width_m']:.3f} m ({100.0 * row['relative']:+.4f} %)")
    print(f"artifacts            : {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())