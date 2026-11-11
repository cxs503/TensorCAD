"""Surface degree elevation: exactness, growth and the rational trap.

Evidence behind ``docs/nurbs_surface_degree_elevation.md``. For the exact rational
cylinder and sphere, a non-rational bicubic patch, and a deliberately rational-blind
control surface it reports:

* degree, control-net shape and knot-multiplicity growth per direction;
* the maximum point deviation ``max |S_elev(u, v) - S(u, v)|`` on a parameter grid;
* the analytic-quantity errors: cylinder radius, sphere radius, ``2 pi R H``,
  ``4 pi R^2``;
* that the two nominal directions commute and that a non-rational patch stays
  non-rational (no silent switch onto the rational evaluation path);
* what happens if the weights are dropped (the trap the homogeneous solve avoids).

Run from the repository root::

    python examples/nurbs_surface_degree_elevation/run_example.py --out examples/nurbs_surface_degree_elevation/outputs
"""
from __future__ import annotations

import argparse
import json
import math
import random
from pathlib import Path

from tensorcad.nurbs import NURBSSurface, cylinder_surface, sphere_surface

CYLINDER_RADIUS = 1.37
CYLINDER_HEIGHT = 5.2
SPHERE_RADIUS = 2.2
U_SAMPLES = 25
V_SAMPLES = 13
GRID = [(i / (U_SAMPLES - 1), j / (V_SAMPLES - 1)) for i in range(U_SAMPLES) for j in range(V_SAMPLES)]


def max_deviation(left: NURBSSurface, right: NURBSSurface) -> float:
    return max(math.dist(left.evaluate(u, v), right.evaluate(u, v)) for u, v in GRID)


def net_shape(surface: NURBSSurface) -> tuple[int, int]:
    return (len(surface.control_net), len(surface.control_net[0]))


def bicubic_patch(seed: int = 20261111) -> NURBSSurface:
    rng = random.Random(seed)
    knots = ((0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0),) * 2
    net = tuple(
        tuple(
            (
                rng.uniform(-1.0, 1.0),
                rng.uniform(-1.0, 1.0),
                0.15 * (i + j) + rng.uniform(-0.2, 0.2),
            )
            for j in range(4)
        )
        for i in range(4)
    )
    return NURBSSurface(degree=(3, 3), knots=knots, control_net=net)


def row(base: str, times_u: int, times_v: int, surface: NURBSSurface, elevated: NURBSSurface) -> dict:
    return {
        "case": base,
        "times_u": times_u,
        "times_v": times_v,
        "degree_before": list(surface.degree),
        "degree_after": list(elevated.degree),
        "net_before": list(net_shape(surface)),
        "net_after": list(net_shape(elevated)),
        "knots_after": [len(elevated.knots[0]), len(elevated.knots[1])],
        "rational": elevated.rational,
        "max_point_deviation_m": max_deviation(elevated, surface),
    }


def cylinder_cases() -> list[dict]:
    analytic_area = 2.0 * math.pi * CYLINDER_RADIUS * CYLINDER_HEIGHT
    base = cylinder_surface(CYLINDER_RADIUS, CYLINDER_HEIGHT)
    rows = []
    for times_u, times_v in ((1, 0), (2, 0), (3, 0), (0, 1), (0, 2), (1, 1), (2, 1), (3, 2)):
        elevated = base.elevate_degree(times_u, times_v)
        entry = row("cylinder", times_u, times_v, base, elevated)
        entry.update(
            {
                "radius_error_m": max(
                    abs(math.hypot(*elevated.evaluate(u, v)[:2]) - CYLINDER_RADIUS)
                    for u, v in GRID
                ),
                "area_analytic_m2": analytic_area,
                "area_elevated_m2": elevated.area(),
                "area_relative_error": abs(elevated.area() / analytic_area - 1.0),
            }
        )
        rows.append(entry)
    return rows


def sphere_cases() -> list[dict]:
    analytic_area = 4.0 * math.pi * SPHERE_RADIUS**2
    base = sphere_surface(SPHERE_RADIUS)
    rows = []
    for times_u, times_v in ((1, 0), (0, 1), (1, 1), (2, 1), (1, 2)):
        elevated = base.elevate_degree(times_u, times_v)
        entry = row("sphere", times_u, times_v, base, elevated)
        entry.update(
            {
                "radius_error_m": max(
                    abs(math.sqrt(sum(c * c for c in elevated.evaluate(u, v))) - SPHERE_RADIUS)
                    for u, v in GRID
                ),
                "area_analytic_m2": analytic_area,
                "area_elevated_m2": elevated.area(),
                "area_relative_error": abs(elevated.area() / analytic_area - 1.0),
            }
        )
        rows.append(entry)
    return rows


def patch_cases() -> dict:
    base = bicubic_patch()
    rows = []
    for times_u, times_v in ((1, 1), (2, 3), (0, 2)):
        elevated = base.elevate_degree(times_u, times_v)
        entry = row("bicubic_patch", times_u, times_v, base, elevated)
        entry["deviation_map_max_m"] = entry["max_point_deviation_m"]
        rows.append(entry)

    u_first = base.elevate_degree(2, 1)
    v_first = base.elevate_degree(0, 1).elevate_degree(2, 0)
    twice = base.elevate_degree(1, 0).elevate_degree(1, 0)
    direct = base.elevate_degree(2, 0)
    return {
        "rows": rows,
        "commutation_max_deviation_m": max_deviation(u_first, v_first),
        "additivity_max_deviation_m": max_deviation(twice, direct),
        "non_rational_preserved": all(not entry["rational"] for entry in rows),
    }


def trap_case() -> dict:
    """Drop the weights and elevate the 3-D net: the surface stops being the cylinder."""
    truth = cylinder_surface(CYLINDER_RADIUS, CYLINDER_HEIGHT)
    blind = NURBSSurface(
        degree=truth.degree, knots=truth.knots, control_net=truth.control_net
    )
    elevated = blind.elevate_degree(2, 0)
    analytic_area = 2.0 * math.pi * CYLINDER_RADIUS * CYLINDER_HEIGHT
    return {
        "radius_error_m": max(
            abs(math.hypot(*elevated.evaluate(u, v)[:2]) - CYLINDER_RADIUS) for u, v in GRID
        ),
        "area_relative_error": abs(elevated.area() / analytic_area - 1.0),
        "honest_route_radius_error_m": max(
            abs(math.hypot(*truth.elevate_degree(2, 0).evaluate(u, v)[:2]) - CYLINDER_RADIUS)
            for u, v in GRID
        ),
    }


def deviation_map(surface: NURBSSurface, elevated: NURBSSurface) -> list[float]:
    return [math.dist(elevated.evaluate(u, v), surface.evaluate(u, v)) for u, v in GRID]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="examples/nurbs_surface_degree_elevation/outputs")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    cylinder = cylinder_cases()
    sphere = sphere_cases()
    patch = patch_cases()
    trap = trap_case()

    sphere_base = sphere_surface(SPHERE_RADIUS)
    field = deviation_map(sphere_base, sphere_base.elevate_degree(2, 1))

    results = {
        "inputs": {
            "cylinder": {
                "radius_m": CYLINDER_RADIUS,
                "height_m": CYLINDER_HEIGHT,
                "degree": [2, 1],
                "area_analytic_m2": 2.0 * math.pi * CYLINDER_RADIUS * CYLINDER_HEIGHT,
            },
            "sphere": {
                "radius_m": SPHERE_RADIUS,
                "degree": [2, 2],
                "area_analytic_m2": 4.0 * math.pi * SPHERE_RADIUS**2,
            },
            "sample_grid": [U_SAMPLES, V_SAMPLES],
        },
        "cylinder": cylinder,
        "sphere": sphere,
        "non_rational_patch": patch,
        "rational_blind_trap": trap,
        "sphere_deviation_field": {
            "times": [2, 1],
            "grid": [U_SAMPLES, V_SAMPLES],
            "max_m": max(field),
            "mean_m": sum(field) / len(field),
        },
    }
    (out / "results.json").write_text(json.dumps(results, indent=2) + "\n")

    worst_deviation = max(
        entry["max_point_deviation_m"] for entry in cylinder + sphere + patch["rows"]
    )
    worst_area = max(entry["area_relative_error"] for entry in cylinder + sphere)
    print(f"worst point deviation      : {worst_deviation:.3e} m")
    print(f"worst analytic area error  : {worst_area:.3e} (relative)")
    print(f"trap: radius lost by       : {trap['radius_error_m']:.3e} m")
    print(f"op -> {out / 'results.json'}")


if __name__ == "__main__":
    main()
