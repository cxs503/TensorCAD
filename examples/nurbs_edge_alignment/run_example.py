"""Evidence behind ``docs/nurbs_edge_alignment.md``.

Three measurements, all of them produced here (nothing in the document is typed by
hand):

1. **Knot insertion** on the exact rational circle: the curve must not move, the
   representation must grow by exactly one control point per insertion, and the
   inserted curve must still have radius ``R`` and length ``2 pi R``.
2. **Alignment of a real patch seam**: a sphere is split at ``u = 0.5`` so the two
   halves share the interface exactly, then the neighbour is re-represented (one
   extra knot plus one degree elevation in ``v``). Before alignment the edges are
   geometrically coincident but have no shared data; after alignment they have equal
   degree, equal knots, equal control points and equal weights.
3. **Two controls that must fail**: a seam whose two edges are *different* curves
   (radii 1.0 and 1.002) must keep its gap, and two representations of the same point
   set with different parametrisations (three 120-degree arcs against six 60-degree
   arcs) must come out with a weight mismatch instead of a watertight claim.

Run::

    PYTHONPATH=$PWD python examples/nurbs_edge_alignment/run_example.py
    # -> examples/nurbs_edge_alignment/outputs/edge_alignment.json
    # -> examples/nurbs_edge_alignment/outputs/edge_alignment.png
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_REPO_ROOT = _Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in _sys.path:  # allow running from a checkout without installing
    _sys.path.insert(0, str(_REPO_ROOT))

import argparse
import json
import math
from pathlib import Path

from tensorcad.edge_alignment import (
    align_edges,
    boundary_curve,
    edge_report,
    insert_knot,
    reverse_curve,
)
from tensorcad.highorder.nurbs_geom import insert_knot as insert_surface_knot
from tensorcad.highorder.nurbs_geom import split_patch
from tensorcad.nurbs import circle_curve, sphere_surface

RADIUS = 1.0
CIRCLE_RADIUS = 2.0


def _parameters(curve, samples: int = 257) -> list[float]:
    low, high = curve.domain
    return [low + (high - low) * index / (samples - 1) for index in range(samples)]


def _max_deviation(first, second) -> float:
    return max(math.dist(first.evaluate(u), second.evaluate(u)) for u in _parameters(first))


def _radius_error(curve, radius: float) -> float:
    return max(abs(math.hypot(*curve.evaluate(u)[:2]) - radius) for u in _parameters(curve))


def _gap_curve(first, second, samples: int = 257) -> tuple[list[float], list[float]]:
    parameters = _parameters(first, samples)
    return parameters, [math.dist(first.evaluate(u), second.evaluate(u)) for u in parameters]


def insertion_evidence() -> dict[str, object]:
    original = circle_curve(CIRCLE_RADIUS, arcs=3)
    inserted = insert_knot(original, 0.25, 2)
    inserted = insert_knot(inserted, 0.7, 1)
    analytic_length = 2.0 * math.pi * CIRCLE_RADIUS
    return {
        "radius": CIRCLE_RADIUS,
        "degree": original.degree,
        "control_points": [len(original.control_points), len(inserted.control_points)],
        "knots": [len(original.knots), len(inserted.knots)],
        "insertions": [{"value": 0.25, "times": 2}, {"value": 0.7, "times": 1}],
        "max_point_deviation": _max_deviation(original, inserted),
        "radius_error_before": _radius_error(original, CIRCLE_RADIUS),
        "radius_error_after": _radius_error(inserted, CIRCLE_RADIUS),
        "arc_length": [original.arc_length(), inserted.arc_length()],
        "analytic_arc_length": analytic_length,
        "arc_length_relative_error": [
            abs(original.arc_length() - analytic_length) / analytic_length,
            abs(inserted.arc_length() - analytic_length) / analytic_length,
        ],
        "rational_after": inserted.rational,
    }


def seam_evidence() -> dict[str, object]:
    surface = sphere_surface(RADIUS)
    left, right = split_patch(surface, 0, 0.5)
    neighbour = insert_surface_knot(right, 1, 0.37, 1).elevate_degree(0, 1)
    first = boundary_curve(left, "u1")
    second = boundary_curve(neighbour, "u0")

    report = {
        "before": edge_report(first, second),
        "after": edge_report(*align_edges(first, second)),
        "deviation_first_edge": _max_deviation(first, align_edges(first, second)[0]),
        "deviation_second_edge": _max_deviation(second, align_edges(first, second)[1]),
        "neighbour_elevation_deviation": max(
            math.dist(
                surface.evaluate(u, v),
                surface.elevate_degree(0, 1).evaluate(u, v),
            )
            for u in (0.0, 0.25, 0.5, 0.75, 1.0)
            for v in (0.0, 0.25, 0.5, 0.75, 1.0)
        ),
    }
    # A left/right orientation check: the neighbour's edge reversed is the same seam.
    reversed_second = reverse_curve(second)
    report["reversed_seam_gap"] = edge_report(first, reversed_second)["max_gap"]
    return report


def controls_evidence() -> dict[str, object]:
    small = circle_curve(1.0, arcs=3)
    coarse = circle_curve(1.002, arcs=3)
    fine = circle_curve(1.0, arcs=6)
    return {
        "different_radius": edge_report(*align_edges(small, coarse)),
        "same_point_set_other_parametrisation": edge_report(*align_edges(small, fine)),
    }


def figure(payload: dict[str, object], out: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    surface = sphere_surface(RADIUS)
    left, right = split_patch(surface, 0, 0.5)
    neighbour = insert_surface_knot(right, 1, 0.37, 1).elevate_degree(0, 1)
    first = boundary_curve(left, "u1")
    second = boundary_curve(neighbour, "u0")
    aligned_first, _aligned_second = align_edges(first, second)
    small = circle_curve(1.0, arcs=3)
    coarse = circle_curve(1.002, arcs=3)
    fine = circle_curve(1.0, arcs=6)

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.8))

    axis = axes[0]
    dense = _parameters(aligned_first, 400)
    axis.plot(
        [aligned_first.evaluate(u)[0] for u in dense],
        [aligned_first.evaluate(u)[2] for u in dense],
        "-",
        color="black",
        linewidth=1.4,
        label="aligned seam edge",
    )
    for curve, style, colour, name in (
        (first, "o--", "tab:blue", f"left patch edge, degree {first.degree}, {len(first.control_points)} pts"),
        (
            second,
            "s--",
            "tab:orange",
            f"neighbour edge, degree {second.degree}, {len(second.control_points)} pts",
        ),
    ):
        axis.plot(
            [point[0] for point in curve.control_points],
            [point[2] for point in curve.control_points],
            style,
            color=colour,
            markersize=4,
            label=name,
        )
    axis.plot(
        [point[0] for point in aligned_first.control_points],
        [point[2] for point in aligned_first.control_points],
        "k+",
        markersize=7,
        label=f"after alignment, both edges: {len(aligned_first.control_points)} pts",
    )
    axis.set_title("(a) sphere seam at u = 0.5\none curve, two representations")
    axis.set_xlabel("x [m]")
    axis.set_ylabel("z [m]")
    axis.legend(fontsize=7)
    axis.set_aspect("equal")

    axis = axes[1]
    for pair, colour, name in (
        ((first, second), "tab:green", "coincident seam (before/after alignment)"),
        ((small, coarse), "tab:red", "radii 1.000 vs 1.002 m"),
        ((small, fine), "tab:purple", "3 arcs vs 6 arcs (same circle)"),
    ):
        parameters, values = _gap_curve(pair[0], pair[1])
        axis.semilogy(
            parameters,
            [max(value, 1e-18) for value in values],
            color=colour,
            label=name,
        )
    axis.axhline(1e-12, color="grey", linestyle=":", label="watertight tolerance 1e-12 m")
    axis.set_title("(b) gap at the same parameter")
    axis.set_xlabel("u")
    axis.set_ylabel("|A(u) - B(u)| [m]")
    axis.legend(fontsize=7)

    axis = axes[2]
    labels = ["left edge", "neighbour edge"]
    before = payload["seam"]["before"]["control_points"]
    after = payload["seam"]["after"]["control_points"]
    positions = range(len(labels))
    axis.bar([p - 0.18 for p in positions], before, width=0.36, label="before", color="tab:blue")
    axis.bar([p + 0.18 for p in positions], after, width=0.36, label="after", color="tab:green")
    axis.set_xticks(list(positions), labels)
    axis.set_ylabel("control points on the shared edge")
    axis.set_title("(c) shared data on the seam")
    for index, (left_value, right_value) in enumerate(zip(before, after, strict=True)):
        axis.text(index - 0.18, left_value + 0.1, str(left_value), ha="center", fontsize=8)
        axis.text(index + 0.18, right_value + 0.1, str(right_value), ha="center", fontsize=8)
    axis.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(out / "edge_alignment.png", dpi=160)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="examples/nurbs_edge_alignment/outputs")
    arguments = parser.parse_args()
    out = Path(arguments.out)
    out.mkdir(parents=True, exist_ok=True)

    payload: dict[str, object] = {
        "insertion": insertion_evidence(),
        "seam": seam_evidence(),
        "controls": controls_evidence(),
        "notes": {
            "radius_m": RADIUS,
            "sphere": "sphere_surface(1.0), degree (2, 2), split at u = 0.5",
            "neighbour": "insert_knot(v, 0.37, 1) then elevate_degree(0, 1)",
            "tolerance_m": 1e-12,
        },
    }
    (out / "edge_alignment.json").write_text(json.dumps(payload, indent=2) + "\n")
    figure(payload, out)

    insertion = payload["insertion"]
    print("insertion:", insertion["control_points"], "max dev", insertion["max_point_deviation"])
    print("seam before:", payload["seam"]["before"]["degree"], payload["seam"]["before"]["control_points"])
    print(
        "seam after:",
        payload["seam"]["after"]["degree"],
        payload["seam"]["after"]["control_points"],
        "watertight",
        payload["seam"]["after"]["watertight"],
        "cp diff",
        payload["seam"]["after"]["max_control_difference"],
    )
    print(
        "radius control gap:",
        payload["controls"]["different_radius"]["max_gap"],
        "watertight",
        payload["controls"]["different_radius"]["watertight"],
    )
    print(
        "same point set, other parametrisation: weight diff",
        payload["controls"]["same_point_set_other_parametrisation"]["max_weight_difference"],
        "gap",
        payload["controls"]["same_point_set_other_parametrisation"]["max_gap"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
