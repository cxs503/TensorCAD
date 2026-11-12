"""Piecewise (break-aware) degree elevation: evidence generator for the docs.

The object under test is a degree-2 curve with a genuine break: an exact rational
quarter circle of radius 2 followed by a straight segment, joined with multiplicity
``degree + 1 = 3``. The script measures, for ``times = 0..3``:

* how many pieces the split produces, their degrees and control-point counts;
* the maximum point deviation between each elevated piece and the original curve
  (sampled between the knots, never on a break, where evaluation is one-sided);
* the deviation of the raised arc from the analytic radius (``|P(u)| - R``);
* the arc length of the pieces against the analytic ``pi R / 2 + 2 sqrt 2``;
* the kink angle between the two one-sided tangents at the break;
* for a three-piece curve whose joins do not even share a point, the size of each
  jump before and after elevation (a break stays a break).

Outputs (``--out``): ``metrics.json``, ``table.md`` (the table quoted in
``docs/nurbs_piecewise_degree_elevation.md``) and ``figure.png``.

Run from the repository root::

    python examples/nurbs_piecewise_elevation/run_example.py --out artifacts/m15428
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
from itertools import pairwise
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # backend must be set before this import

from tensorcad.nurbs import NURBSCurve, NURBSError

RADIUS = 2.0
ARC_LENGTH = math.pi * RADIUS / 2.0
LINE_LENGTH = math.hypot(2.0, 2.0)
TOTAL_LENGTH = ARC_LENGTH + LINE_LENGTH


def composite_curve() -> NURBSCurve:
    """Exact rational quarter circle (R = 2) + straight segment, joined at a break."""
    return NURBSCurve(
        degree=2,
        knots=(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0),
        control_points=(
            (2.0, 0.0, 0.0),
            (2.0, 2.0, 0.0),
            (0.0, 2.0, 0.0),
            (0.0, 2.0, 0.0),
            (-1.0, 1.0, 0.0),
            (-2.0, 0.0, 0.0),
        ),
        weights=(1.0, math.cos(math.pi / 4.0), 1.0, 1.0, 1.0, 1.0),
    )


def jump_curve() -> NURBSCurve:
    """Degree 3, three pieces, two breaks that do not share a point."""
    knots = (0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.0, 2.0, 2.0)
    points = tuple((float(index), float((index * index) % 5), 0.0) for index in range(len(knots) - 4))
    return NURBSCurve(degree=3, knots=knots, control_points=points)


def breaks_of(curve: NURBSCurve) -> list[float]:
    low, high = curve.domain
    values = sorted({value for value in curve.knots if low < value < high})
    return [value for value in values if sum(1 for knot in curve.knots if knot == value) > curve.degree]


def sample(segment: NURBSCurve, per_span: int = 61) -> list[tuple[float, tuple[float, float, float]]]:
    low, high = segment.domain
    distinct = sorted({value for value in segment.knots if low <= value <= high})
    parameters = [
        left + (right - left) * index / per_span
        for left, right in pairwise(distinct)
        for index in range(per_span + 1)
    ]
    return [(u, segment.evaluate(u)) for u in parameters]


def max_deviation(curve: NURBSCurve, piece: NURBSCurve) -> float:
    worst = 0.0
    for u, point in sample(piece):
        reference = curve.evaluate(u)
        worst = max(worst, math.dist(point, reference))
    return worst


def elevation_report(times: int) -> dict[str, object]:
    curve = composite_curve()
    pieces = curve.split_at_discontinuities()
    raised = curve.elevate_degree_piecewise(times) if times else pieces
    deviations = [max_deviation(curve, piece) for piece in raised]
    arc = raised[0]
    radius_error = max(
        abs(math.hypot(*arc.evaluate(index / 200.0)[:2]) - RADIUS) for index in range(201)
    )
    lengths = [piece.arc_length() for piece in raised]
    left, right = raised[0].tangent(1.0), raised[1].tangent(1.0)
    cosine = max(-1.0, min(1.0, sum(a * b for a, b in zip(left, right, strict=True))))
    return {
        "times": times,
        "pieces": len(raised),
        "degrees": [piece.degree for piece in raised],
        "control_points": [len(piece.control_points) for piece in raised],
        "max_deviation": max(deviations),
        "radius_error": radius_error,
        "arc_length": lengths[0],
        "line_length": lengths[1],
        "total_length": sum(lengths),
        "total_length_error": sum(lengths) - TOTAL_LENGTH,
        "kink_angle_deg": math.degrees(math.acos(cosine)),
    }


def jump_report() -> list[dict[str, float]]:
    curve = jump_curve()
    pieces = curve.split_at_discontinuities()
    raised = curve.elevate_degree_piecewise(2)
    rows = []
    for index in range(len(pieces) - 1):
        before = math.dist(
            pieces[index].evaluate(pieces[index].domain[1]),
            pieces[index + 1].evaluate(pieces[index + 1].domain[0]),
        )
        after = math.dist(
            raised[index].evaluate(raised[index].domain[1]),
            raised[index + 1].evaluate(raised[index + 1].domain[0]),
        )
        rows.append({"break": pieces[index].domain[1], "jump_before": before, "jump_after": after})
    return rows


def figure(out: Path, reports: list[dict[str, object]]) -> None:
    curve = composite_curve()
    pieces = curve.split_at_discontinuities()
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6))

    ax = axes[0]
    designs = [
        ("original (break at u=1)", pieces, "-"),
        ("elevated x1 (degree 3)", curve.elevate_degree_piecewise(1), "--"),
        ("elevated x3 (degree 5)", curve.elevate_degree_piecewise(3), ":"),
    ]
    for label, segments, style in designs:
        for segment in segments:
            points = sample(segment, per_span=24)
            ax.plot([p[0] for _, p in points], [p[1] for _, p in points], style, lw=1.6, label=label)
            label = None  # one legend entry per design
    ax.plot([-2, 2], [0, 0], "k.", ms=3, alpha=0.25)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("quarter circle R=2 + straight segment, split at the break")
    ax.legend(loc="lower left", fontsize=7)
    ax.grid(alpha=0.3)

    ax = axes[1]
    parameters = [index / 200.0 for index in range(201)]
    for report in reports:
        raised = curve.elevate_degree_piecewise(int(report["times"])) if report["times"] else pieces
        deviations = []
        for u in parameters:
            segment = raised[0] if u <= 1.0 else raised[1]
            deviations.append(max(math.dist(segment.evaluate(u), curve.evaluate(u)), 1e-17))
        ax.semilogy(parameters, deviations, lw=1.4, label=f"degree {report['degrees'][0]} (x{report['times']})")
    ax.axvline(1.0, color="k", lw=0.8, ls="--")
    ax.set_xlabel("parameter u")
    ax.set_ylabel("|elevated - original| [m]")
    ax.set_title("point deviation (machine precision), break at u=1")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    fig.savefig(out / "figure.png", dpi=160)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("artifacts/m15428"))
    arguments = parser.parse_args()
    arguments.out.mkdir(parents=True, exist_ok=True)

    reports = [elevation_report(times) for times in range(4)]
    rows = jump_report()

    try:
        composite_curve().elevate_degree(1)
        refusal = "NOT RAISED"
    except NURBSError as error:
        refusal = str(error)

    payload = {
        "case": {
            "degree": 2,
            "pieces": 2,
            "break_multiplicity": 3,
            "radius": RADIUS,
            "analytic_arc_length": ARC_LENGTH,
            "analytic_line_length": LINE_LENGTH,
            "analytic_total_length": TOTAL_LENGTH,
        },
        "elevation": reports,
        "jumps": rows,
        "single_curve_refusal": refusal,
    }
    (arguments.out / "metrics.json").write_text(json.dumps(payload, indent=2) + "\n")

    header = (
        "| times | degree | control points | max |elev - original| [m] | radius error [m] "
        "| arc length [m] | line length [m] | total - analytic [m] | kink angle [deg] |\n"
        "|---|---|---|---|---|---|---|---|---|\n"
    )
    body = "".join(
        f"| {report['times']} | {report['degrees']} | {report['control_points']} "
        f"| {report['max_deviation']:.3e} | {report['radius_error']:.3e} "
        f"| {report['arc_length']:.15f} | {report['line_length']:.15f} "
        f"| {report['total_length_error']:+.3e} | {report['kink_angle_deg']:.12f} |\n"
        for report in reports
    )
    jump_table = (
        "\n| break | jump before [m] | jump after x2 [m] |\n|---|---|---|\n"
        + "".join(
            f"| {row['break']} | {row['jump_before']:.15f} | {row['jump_after']:.15f} |\n"
            for row in rows
        )
    )
    (arguments.out / "table.md").write_text(header + body + jump_table)

    figure(arguments.out, reports)

    print(header + body + jump_table)
    print(f"analytic total length = {TOTAL_LENGTH:.15f} m")
    print(f"single-curve elevation refused with: {refusal}")
    print(f"artifacts written to {arguments.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
