"""Fillet/chamfer evidence run: exact values, convergence orders, and a bilge-radius table.

Usage: ``python examples/fillet_corners/run_example.py --out artifacts/fillet_corners``

Writes ``fillet_report.json`` (every number in ``docs/fillet_chamfer.md`` comes from it) and
``fillet_corners.png`` (the rounded sections, and the two measured convergence orders).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.fillet import (  # noqa: E402
    chamfer_section,
    chamfer_solid,
    fillet_section,
    fillet_solid,
    max_fillet_radius,
    polygon_moments,
    rectangle_section,
    section_turns,
)

SQUARE = rectangle_section(2.0, 2.0, centred=True)
L_SHAPE = ((0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0))
HEXAGON = tuple(
    (1.0 * math.cos(2.0 * math.pi * index / 6), 1.0 * math.sin(2.0 * math.pi * index / 6))
    for index in range(6)
)


def _order(errors: list[float]) -> list[float]:
    return [
        math.log(errors[index] / errors[index + 1]) / math.log(2.0)
        for index in range(len(errors) - 1)
    ]


def build_report() -> dict[str, object]:
    fixtures = {}
    radius = 0.25
    for name, ring, size in (
        ("square", SQUARE, 0.5),
        ("hexagon", HEXAGON, 0.2),
        ("l_shape", L_SHAPE, 0.2),
    ):
        shape = fillet_section(ring, size, arc_points=32)
        fixtures[name] = {
            "points": len(ring),
            "turn_deg": [round(math.degrees(turn), 9) for turn in section_turns(ring)],
            "radius": shape.achieved,
            "base_area": shape.base_area,
            "area": shape.moments.area,
            "area_closed_form": shape.base_area + shape.area_coefficient * shape.achieved**2,
            "area_crosscheck_error": shape.area_crosscheck_error,
            "perimeter": shape.moments.perimeter,
            "g1_error_deg": shape.g1_error_deg,
            "max_feasible_radius": max_fillet_radius(ring),
            "outline": [[u, v] for u, v in shape.points],
        }

    # The circle limit: a square of side 2r filleted with r.
    circle = fillet_section(SQUARE, 1.0, arc_points=64)
    circle_limit = {
        "radius": circle.achieved,
        "area": circle.moments.area,
        "area_exact": math.pi,
        "area_rel_error": abs(circle.moments.area - math.pi) / math.pi,
        "perimeter": circle.moments.perimeter,
        "perimeter_exact": 2.0 * math.pi,
        "perimeter_rel_error": abs(circle.moments.perimeter - 2.0 * math.pi) / (2.0 * math.pi),
        "second_u": circle.moments.second_u,
        "second_u_exact": math.pi / 4.0,
        "second_u_rel_error": abs(circle.moments.second_u - math.pi / 4.0) / (math.pi / 4.0),
    }

    # Chamfer against the same closed forms.
    chamfer = chamfer_section(SQUARE, 0.25)
    chamfer_block = {
        "setback": chamfer.achieved,
        "area": chamfer.moments.area,
        "area_closed_form": 4.0 - 4.0 * 0.5 * 0.25**2,
        "area_crosscheck_error": chamfer.area_crosscheck_error,
        "perimeter": chamfer.moments.perimeter,
        "perimeter_closed_form": 8.0 - 8.0 * 0.25 + 4.0 * math.sqrt(2.0) * 0.25,
        "crease_deg": chamfer.g1_error_deg,
    }

    # Discretisation: sampled area and mesh volume against the exact values.
    exact_area = fillet_section(SQUARE, 0.5, arc_points=512).moments.area
    area_errors = []
    volume_errors = []
    for points in (4, 8, 16, 32):
        sampled = fillet_section(SQUARE, 0.5, arc_points=points)
        area_errors.append(abs(sampled.sampled_moments().area - exact_area) / exact_area)
        solid = fillet_solid(SQUARE, 3.0, 0.5, arc_points=points)
        volume_errors.append(abs(solid.volume_mesh - solid.volume_exact) / solid.volume_exact)
    convergence = {
        "arc_points": [4, 8, 16, 32],
        "sampled_area_rel_error": area_errors,
        "sampled_area_order": _order(area_errors),
        "mesh_volume_rel_error": volume_errors,
        "mesh_volume_order": _order(volume_errors),
    }

    # Tapered fillet: tangency survives a radius change, and the exact volume is quadratic.
    taper = fillet_solid(SQUARE, 2.0, 0.1, radius_end=0.6, stations=4, arc_points=16)
    taper_block = {
        "size_start": taper.size_start,
        "size_end": taper.size_end,
        "volume_exact": taper.volume_exact,
        "volume_mesh": taper.volume_mesh,
        "mesh_volume_rel_error": taper.mesh_volume_error,
        "g1_error_deg": taper.g1_error_deg,
        "watertight": taper.watertight,
        "euler_characteristic": taper.euler_characteristic,
        "zero_area_faces": taper.zero_area_faces,
    }

    # Short-edge fallback: a 2 x 0.2 strip cannot host r = 0.5, the clamp names the limit.
    strip = rectangle_section(2.0, 0.2)
    clamped = fillet_section(strip, 0.5, allow_clamp=True)
    fallback = {
        "requested": 0.5,
        "max_feasible": max_fillet_radius(strip),
        "achieved": clamped.achieved,
        "clamped": clamped.clamped,
        "area": clamped.moments.area,
    }

    # Ship application: bilge radius on a box barge waterplane (B = 10 m, T = 4 m).
    beam, draft = 10.0, 4.0
    bilge = []
    for ratio in (0.0, 0.25, 0.5, 1.0):
        radius_m = ratio * draft / 2.0
        shape = fillet_section(
            rectangle_section(beam, draft, centred=True), radius_m, arc_points=32
        )
        bilge.append(
            {
                "bilge_radius": radius_m,
                "waterplane_area": shape.moments.area,
                "i_transverse": shape.moments.second_u,
                "i_over_area": shape.moments.second_u / shape.moments.area,
                "bm_m": shape.moments.second_u / shape.moments.area / draft,
                "bm_box_m": beam**2 / (12.0 * draft),
            }
        )

    return {
        "schema": "tensorcad.fillet.example/1",
        "fixtures": fixtures,
        "circle_limit": circle_limit,
        "chamfer": chamfer_block,
        "convergence": convergence,
        "taper": taper_block,
        "fallback": fallback,
        "bilge_radius": bilge,
        "chamfer_solid": {
            "volume_exact": chamfer_solid(SQUARE, 3.0, 0.25).volume_exact,
            "volume_mesh": chamfer_solid(SQUARE, 3.0, 0.25).volume_mesh,
        },
        "polygon_check": polygon_moments(SQUARE).to_dict(),
    }


def plot(report: dict[str, object], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.6))
    for name, colour in (("square", "tab:blue"), ("hexagon", "tab:green"), ("l_shape", "tab:red")):
        outline = report["fixtures"][name]["outline"]
        us = [point[0] for point in outline] + [outline[0][0]]
        vs = [point[1] for point in outline] + [outline[0][1]]
        left.plot(us, vs, "-", color=colour, lw=1.6, label=f"{name} (r={report['fixtures'][name]['radius']})")
    left.set_aspect("equal")
    left.set_title("Filleted sections (exact arcs, chorded outline)")
    left.set_xlabel("u")
    left.set_ylabel("v")
    left.grid(True, alpha=0.3)
    left.legend(fontsize=8)

    convergence = report["convergence"]
    right.loglog(convergence["arc_points"], convergence["sampled_area_rel_error"], "o-", label="sampled area")
    right.loglog(convergence["arc_points"], convergence["mesh_volume_rel_error"], "s-", label="mesh volume")
    reference = [4e-2 / (points**2) for points in convergence["arc_points"]]
    right.loglog(convergence["arc_points"], reference, "k--", lw=0.8, label="slope -2")
    right.set_xlabel("chords per arc")
    right.set_ylabel("relative error vs exact")
    right.set_title("Chord discretisation converges at order 2")
    right.grid(True, which="both", alpha=0.3)
    right.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(path, dpi=140)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("artifacts/fillet_corners"))
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    report = build_report()
    payload = json.dumps(report, indent=2, sort_keys=True)
    (args.out / "fillet_report.json").write_text(payload + "\n", encoding="utf-8")
    plot(report, args.out / "fillet_corners.png")
    print(payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())