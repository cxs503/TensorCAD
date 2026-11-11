"""Degree elevation benchmark: geometry, continuity and representation growth.

The script is the evidence behind ``docs/nurbs_degree_elevation.md``. It reports, for
an exact rational circle and a cubic C2 B-spline:

* the degree / control-point / knot-multiplicity growth;
* the maximum point deviation between the original and the elevated curve, sampled at
  every knot and between them;
* the deviation of the elevated circle from the analytic radius;
* the relative arc-length error against ``2 * pi * R``;
* the measured derivative jump at the interior knot (declared continuity preserved).

Run from the repository root::

    python examples/nurbs_degree_elevation/run_example.py --out artifacts/m14481
"""
from __future__ import annotations

import argparse
import json
import math
from itertools import pairwise
from pathlib import Path

from tensorcad.nurbs import NURBSCurve, circle_curve, curve_derivative_continuity

RADIUS = 2.0


def sample_parameters(curve: NURBSCurve, per_span: int = 33) -> list[float]:
    low, high = curve.domain
    distinct = sorted({value for value in curve.knots if low <= value <= high})
    parameters = list(distinct)
    for left, right in pairwise(distinct):
        for index in range(1, per_span):
            parameters.append(left + (right - left) * index / per_span)
    return parameters


def max_deviation(first: NURBSCurve, second: NURBSCurve, per_span: int = 33) -> float:
    return max(
        math.dist(first.evaluate(u), second.evaluate(u))
        for u in sample_parameters(first, per_span)
    )


def circle_row(times: int) -> dict[str, float | int]:
    original = circle_curve(RADIUS, arcs=3)
    elevated = original.elevate_degree(times)
    deviation = max(
        abs(math.hypot(*elevated.evaluate(u)[:2]) - RADIUS)
        for u in sample_parameters(elevated)
    )
    analytic_length = 2.0 * math.pi * RADIUS
    return {
        "times": times,
        "degree": elevated.degree,
        "control_points": len(elevated.control_points),
        "knots": len(elevated.knots),
        "interior_multiplicity": elevated.knots.count(1.0 / 3.0),
        "max_point_deviation_m": max_deviation(original, elevated),
        "max_radius_deviation_m": deviation,
        "arc_length_m": elevated.arc_length(),
        "arc_length_relative_error": abs(elevated.arc_length() - analytic_length)
        / analytic_length,
    }


def cubic_continuity() -> dict[str, float | int]:
    original = NURBSCurve(
        degree=3,
        knots=(0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0),
        control_points=((0, 0, 0), (1, 2, 0), (2, 0, 0), (3, 3, 0), (4, 1, 0)),
    )
    elevated = original.elevate_degree(1)
    records = curve_derivative_continuity(elevated, order=2)
    record = records[0]
    return {
        "degree_before": original.degree,
        "degree_after": elevated.degree,
        "control_points_before": len(original.control_points),
        "control_points_after": len(elevated.control_points),
        "interior_multiplicity_after": elevated.knots.count(0.4),
        "declared_order": record.declared_order,
        "measured_max_jump": max(record.jumps),
        "max_point_deviation_m": max_deviation(original, elevated),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="artifacts/m14481")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    rows = [circle_row(times) for times in (0, 1, 2, 3)]
    cubic = cubic_continuity()

    # The operation is only useful if it preserves the geometry, so fail loudly.
    for row in rows:
        assert row["max_point_deviation_m"] < 1e-12, row
        assert row["max_radius_deviation_m"] < 1e-12, row
        assert row["arc_length_relative_error"] < 1e-12, row
    assert cubic["measured_max_jump"] < 1e-12, cubic
    assert cubic["declared_order"] == 2, cubic
    assert cubic["control_points_after"] == 7, cubic

    payload = {
        "case": {
            "geometry": "exact rational circle, three 120-degree arcs",
            "radius_m": RADIUS,
            "analytic_arc_length_m": 2.0 * math.pi * RADIUS,
            "cubic": "single-span-pair cubic, interior knot 0.4 of multiplicity 1 (C2)",
        },
        "circle": rows,
        "cubic": cubic,
    }
    (out / "degree_elevation.json").write_text(json.dumps(payload, indent=2) + "\n")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - plotting is evidence, not a runtime dep
        print(f"matplotlib unavailable, JSON only: {exc}")
        return 0

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    original = circle_curve(RADIUS, arcs=3)
    elevated = original.elevate_degree(2)
    for curve, label, colour in (
        (original, f"original degree {original.degree}", "#1f77b4"),
        (elevated, f"elevated degree {elevated.degree}", "#d62728"),
    ):
        xs = [point[0] for point in curve.control_points]
        ys = [point[1] for point in curve.control_points]
        left.plot(xs, ys, "o-", ms=3, lw=1.0, color=colour, label=label)
    theta = [2.0 * math.pi * index / 720 for index in range(721)]
    left.plot(
        [RADIUS * math.cos(t) for t in theta],
        [RADIUS * math.sin(t) for t in theta],
        "k--",
        lw=0.8,
        label="analytic circle",
    )
    left.set_aspect("equal")
    left.set_xlabel("x [m]")
    left.set_ylabel("y [m]")
    left.set_title("control polygons vs the exact circle")
    left.legend(fontsize=7)

    times = [row["times"] for row in rows]
    deviations = [max(row["max_point_deviation_m"], 1e-18) for row in rows]
    right.semilogy(times, deviations, "o-", color="#2ca02c")
    right.set_xlabel("degree elevations applied")
    right.set_ylabel("max |C(u) - C_elevated(u)| [m]")
    right.set_title("geometry is preserved to round-off")
    right.grid(True, which="both", alpha=0.3)
    figure.tight_layout()
    figure.savefig(out / "degree_elevation.png", dpi=160)

    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
