"""Righting arms at finite heel: a closed-form section, a literature formula, then a hull.

Everything printed here is measured, and the script exits non-zero if one of its own
claims stops holding -- so the tables in ``docs/gz_curve.md`` cannot drift away from the
code the way a hand-copied table does.

What it does, in order:

1. **A box with freeboard**, where the section is a rectangle and an independent 2-D
   clipping gives the exact arm at *any* heel. The 3-D integral must reproduce it to
   round-off, past deck and keel immersion where no closed form exists.
2. **The wall-sided formula** ``GZ = sin(phi)(BM + BM tan^2(phi)/2 - BG)``, which is only
   valid while both ends of the waterline stay on the vertical sides
   (``tan(phi) < 2T/B``). Inside the window it must agree to round-off; outside it must
   not, and the angle at which it leaves is reported rather than assumed.
3. **The parametric hull** of :mod:`tensorcad.hull`, with wall-sided topsides added over
   its design waterline. ``GM`` from the waterplane moment and the curve's slope at the
   origin are two different integrals; they must land on the same number.
4. **The Wigley** parabolic hull, as a second family with a different section shape.

Run it with ``python examples/gz_curve/run_example.py --output docs/gz_curve``.
matplotlib is optional: the JSON is written either way and the PNGs are skipped.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from tensorcad import __version__
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.stability import initial_stability, righting_arm_curve
from tensorcad.wigley import generate_wigley


# ---------------------------------------------------------------------------
# reference bodies and references
# ---------------------------------------------------------------------------
def box(length: float, beam: float, draft: float, freeboard: float) -> Mesh:
    """Closed box: keel on ``z = -draft``, deck on ``z = +freeboard``, waterline ``z = 0``."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, freeboard):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    left, right = 0.0, length
    lo, hi = -0.5 * beam, 0.5 * beam
    bottom, top = -draft, freeboard
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


def reference_arm(
    beam: float, draft: float, freeboard: float, kg: float, heel_deg: float
) -> tuple[float, float]:
    """Exact arm and waterplane height of the rectangular section (independent 2-D code)."""
    phi = math.radians(heel_deg)
    cos_phi, sin_phi = math.cos(phi), math.sin(phi)
    keel = -draft
    section = [
        (y * cos_phi - z * sin_phi, y * sin_phi + z * cos_phi)
        for y, z in (
            (-0.5 * beam, keel),
            (0.5 * beam, keel),
            (0.5 * beam, freeboard),
            (-0.5 * beam, freeboard),
        )
    ]
    heights = [point[1] for point in section]
    low, high = min(heights) - 1.0, max(heights) + 1.0

    def clipped(plane: float) -> tuple[float, float]:
        polygon: list[tuple[float, float]] = []
        for index in range(len(section)):
            start, end = section[index], section[(index + 1) % len(section)]
            start_in, end_in = start[1] <= plane, end[1] <= plane
            if start_in:
                polygon.append(start)
            if start_in != end_in:
                ratio = (plane - start[1]) / (end[1] - start[1])
                polygon.append((start[0] + ratio * (end[0] - start[0]), plane))
        area = moment = 0.0
        for index in range(len(polygon)):
            x1, y1 = polygon[index]
            x2, y2 = polygon[(index + 1) % len(polygon)]
            cross = x1 * y2 - x2 * y1
            area += cross
            moment += (x1 + x2) * cross
        return 0.5 * area, moment / 6.0

    for _ in range(200):
        middle = 0.5 * (low + high)
        if clipped(middle)[0] < beam * draft:
            low = middle
        else:
            high = middle
    plane = 0.5 * (low + high)
    area, moment = clipped(plane)
    return -(keel + kg) * sin_phi - moment / area, plane


def wall_sided(beam: float, draft: float, kg: float, heel_deg: float) -> float:
    phi = math.radians(heel_deg)
    bm = beam * beam / (12.0 * draft)
    bg = kg - 0.5 * draft
    return math.sin(phi) * (bm + 0.5 * bm * math.tan(phi) ** 2 - bg)


def secant_slope(curve) -> float:
    """``GZ`` derivative at the origin from the first two points, in m/rad.

    A secant over the first step, not a fitted derivative: the grid is the caller's, and
    a one-degree first step keeps the ``O(phi^2)`` truncation small enough to compare
    against ``GM`` without pretending the limit was taken.
    """
    first, second = curve.points[0], curve.points[1]
    return (second.gz - first.gz) / math.radians(second.heel_deg - first.heel_deg)


# ---------------------------------------------------------------------------
# the three cases
# ---------------------------------------------------------------------------
LENGTH, BEAM, DRAFT, FREEBOARD, KG = 10.0, 2.0, 0.5, 1.5, 0.45


def box_case() -> dict:
    mesh = box(LENGTH, BEAM, DRAFT, FREEBOARD)
    curve = righting_arm_curve(
        mesh, kg=KG, waterline_z=0.0, angles_deg=range(0, 181, 5)
    )
    wall_sided_limit = math.degrees(math.atan(2.0 * DRAFT / BEAM))
    rows = []
    for point in curve.points:
        expected, plane = reference_arm(BEAM, DRAFT, FREEBOARD, KG, point.heel_deg)
        closed = wall_sided(BEAM, DRAFT, KG, point.heel_deg)
        rows.append(
            {
                "heel_deg": point.heel_deg,
                "gz": point.gz,
                "gz_reference": expected,
                "absolute_error": abs(point.gz - expected),
                "gz_wall_sided": closed,
                "wall_sided_relative_error": (
                    abs(point.gz / closed - 1.0) if closed else 0.0
                ),
                "waterplane_height": plane,
                "deck_immersed": point.deck_immersed,
                "displacement": point.displacement_volume,
            }
        )
    inside = [row for row in rows if row["heel_deg"] < wall_sided_limit]
    outside = [row for row in rows if row["heel_deg"] > wall_sided_limit + 3.0]
    worst_reference = max(row["absolute_error"] for row in rows)
    worst_wall_sided = max(row["wall_sided_relative_error"] for row in inside)
    first_outside = min(row["wall_sided_relative_error"] for row in outside)
    assert worst_reference < 1e-11, worst_reference
    assert worst_wall_sided < 1e-11, worst_wall_sided
    assert first_outside > 1e-3, first_outside
    expected_gm = 0.5 * DRAFT + BEAM * BEAM / (12.0 * DRAFT) - KG
    assert abs(curve.initial_gm_transverse - expected_gm) < 1e-12
    return {
        "length": LENGTH,
        "beam": BEAM,
        "draft": DRAFT,
        "freeboard": FREEBOARD,
        "kg": KG,
        "displacement_volume": curve.displacement_volume,
        "enclosed_volume": curve.enclosed_volume,
        "gm_transverse": curve.initial_gm_transverse,
        "gm_closed_form": expected_gm,
        "small_angle_slope": secant_slope(curve),
        "gz_max": curve.gz_max,
        "angle_of_gz_max_deg": curve.angle_of_gz_max_deg,
        "vanishing_angle_deg": curve.vanishing_angle_deg,
        "dynamic_stability_area": curve.dynamic_stability_area,
        "wall_sided_limit_deg": wall_sided_limit,
        "worst_reference_error": worst_reference,
        "worst_wall_sided_relative_error_inside": worst_wall_sided,
        "first_wall_sided_relative_error_outside": first_outside,
        "deck_immersion_deg": next(
            (row["heel_deg"] for row in rows if row["deck_immersed"]), None
        ),
        "rows": rows,
    }


def hull_case() -> dict:
    mesh = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8)).mesh
    kg, freeboard = 0.6, 0.4
    upright = initial_stability(mesh, kg=kg)
    curve = righting_arm_curve(
        mesh,
        kg=kg,
        waterline_z=0.0,
        freeboard=freeboard,
        angles_deg=(0.0, 1.0, *range(5, 91, 5)),
    )
    slope = secant_slope(curve)
    rows = [
        {
            "heel_deg": point.heel_deg,
            "gz": point.gz,
            "immersion_depth": point.immersion_depth,
            "centre_of_buoyancy_y": point.centre_of_buoyancy_y,
            "deck_immersed": point.deck_immersed,
        }
        for point in curve.points
    ]
    assert curve.initial_gm_transverse == pytest_approx(upright.gm_transverse, 1e-12)
    assert abs(slope / upright.gm_transverse - 1.0) < 0.01, slope / upright.gm_transverse
    return {
        "length": 10.0,
        "beam": 2.0,
        "draft": 0.8,
        "kg": kg,
        "freeboard": freeboard,
        "displacement_volume": curve.displacement_volume,
        "enclosed_volume": curve.enclosed_volume,
        "gm_transverse": curve.initial_gm_transverse,
        "gm_initial_stability": upright.gm_transverse,
        "small_angle_slope": slope,
        "slope_relative_error": abs(slope / upright.gm_transverse - 1.0),
        "gz_max": curve.gz_max,
        "angle_of_gz_max_deg": curve.angle_of_gz_max_deg,
        "vanishing_angle_deg": curve.vanishing_angle_deg,
        "dynamic_stability_area": curve.dynamic_stability_area,
        "deck_immersion_deg": next(
            (row["heel_deg"] for row in rows if row["deck_immersed"]), None
        ),
        "rows": rows,
    }


def wigley_case() -> dict:
    mesh = generate_wigley().mesh
    kg, freeboard = 0.4, 0.3
    curve = righting_arm_curve(
        mesh, kg=kg, waterline_z=0.0, freeboard=freeboard, angles_deg=range(0, 61, 5)
    )
    assert curve.points[0].gz == pytest_approx(0.0, 1e-12)
    assert curve.gz_max > 0.0
    return {
        "freeboard": freeboard,
        "kg": kg,
        "displacement_volume": curve.displacement_volume,
        "gm_transverse": curve.initial_gm_transverse,
        "gz_max": curve.gz_max,
        "angle_of_gz_max_deg": curve.angle_of_gz_max_deg,
        "dynamic_stability_area": curve.dynamic_stability_area,
        "rows": [
            {"heel_deg": point.heel_deg, "gz": point.gz} for point in curve.points
        ],
    }


def pytest_approx(value: float, relative: float) -> float:
    """A tolerance object without importing pytest: exact JSON numbers are what matter."""
    return _Approx(value, relative)


class _Approx:
    def __init__(self, value: float, relative: float) -> None:
        self.value = value
        self.relative = relative

    def __eq__(self, other: object) -> bool:
        return abs(float(other) - self.value) <= self.relative * max(abs(self.value), 1.0)


# ---------------------------------------------------------------------------
# reporting
# ---------------------------------------------------------------------------
def plot(box_data: dict, hull_data: dict, wigley_data: dict, output: Path) -> str | None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # pragma: no cover - the plot is optional
        return None
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    heels = [row["heel_deg"] for row in box_data["rows"]]
    left.plot(heels, [row["gz"] for row in box_data["rows"]], "k-", lw=2, label="measured (3-D clip)")
    left.plot(
        heels,
        [row["gz_reference"] for row in box_data["rows"]],
        "r--",
        lw=1,
        label="exact section (2-D)",
    )
    left.plot(
        heels,
        [row["gz_wall_sided"] for row in box_data["rows"]],
        "b:",
        lw=1,
        label="wall-sided formula",
    )
    left.axvline(box_data["wall_sided_limit_deg"], color="b", lw=0.6, alpha=0.6)
    left.set_xlabel("heel angle [deg]")
    left.set_ylabel("GZ [m]")
    left.set_title(f"box {LENGTH:g}x{BEAM:g}x{DRAFT:g}+{FREEBOARD:g} m, KG={KG:g} m")
    left.legend(fontsize=8)
    left.grid(alpha=0.3)
    for name, data, style in (
        ("parametric hull + 0.4 m topsides", hull_data, "o-"),
        ("Wigley + 0.3 m topsides", wigley_data, "s--"),
    ):
        right.plot(
            [row["heel_deg"] for row in data["rows"]],
            [row["gz"] for row in data["rows"]],
            style,
            ms=3,
            lw=1.4,
            label=name,
        )
    right.set_xlabel("heel angle [deg]")
    right.set_ylabel("GZ [m]")
    right.set_title("hull families (wall-sided topsides)")
    right.legend(fontsize=8)
    right.grid(alpha=0.3)
    figure.tight_layout()
    path = output / "gz_curves.png"
    figure.savefig(path, dpi=130)
    plt.close(figure)
    return path.name


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/gz_curve"))
    arguments = parser.parse_args(argv)
    arguments.output.mkdir(parents=True, exist_ok=True)

    box_data = box_case()
    hull_data = hull_case()
    wigley_data = wigley_case()
    figure = plot(box_data, hull_data, wigley_data, arguments.output)

    print(f"tensorcad {__version__} - righting arms at finite heel")
    print(
        f"\nbox {LENGTH:g} x {BEAM:g} x {DRAFT:g} m + {FREEBOARD:g} m freeboard, "
        f"KG = {KG:g} m, displacement {box_data['displacement_volume']:.6f} m^3"
    )
    print(f"{'heel':>6} {'GZ [m]':>12} {'exact':>12} {'|err|':>9} {'wall-sided':>11} {'rel':>9}")
    for row in box_data["rows"]:
        if int(row["heel_deg"]) % 15 and row["heel_deg"] not in (0.0, 60.0):
            continue
        print(
            f"{row['heel_deg']:6.0f} {row['gz']:12.9f} {row['gz_reference']:12.9f} "
            f"{row['absolute_error']:9.2e} {row['gz_wall_sided']:11.6f} "
            f"{row['wall_sided_relative_error']:9.2e}"
        )
    print(
        f"worst |measured - exact section| = {box_data['worst_reference_error']:.3e} m "
        "over 0..180 deg"
    )
    print(
        "wall-sided agreement inside tan(phi) < 2T/B "
        f"(phi < {box_data['wall_sided_limit_deg']:.2f} deg): "
        f"{box_data['worst_wall_sided_relative_error_inside']:.3e} relative"
    )
    print(
        "wall-sided departure outside it: "
        f"{box_data['first_wall_sided_relative_error_outside']:.3e} relative "
        f"(deck edge immerses at {box_data['deck_immersion_deg']} deg)"
    )
    print(
        f"GZ_max {box_data['gz_max']:.6f} m at {box_data['angle_of_gz_max_deg']:.0f} deg, "
        f"dynamic area {box_data['dynamic_stability_area']:.6f} m*rad"
    )
    print(
        f"GM {box_data['gm_transverse']:.9f} m (closed form {box_data['gm_closed_form']:.9f}), "
        f"secant slope at 0 {box_data['small_angle_slope']:.6f} m/rad"
    )

    print(
        f"\nparametric hull 10x2x0.8 m + {hull_data['freeboard']:g} m topsides, "
        f"KG = {hull_data['kg']:g} m, displacement {hull_data['displacement_volume']:.6f} m^3"
    )
    print(f"{'heel':>6} {'GZ [m]':>12} {'immersion [m]':>14} {'y_B [m]':>10} {'deck':>6}")
    for row in hull_data["rows"]:
        print(
            f"{row['heel_deg']:6.0f} {row['gz']:12.9f} {row['immersion_depth']:14.9f} "
            f"{row['centre_of_buoyancy_y']:10.6f} {row['deck_immersed']!s:>6}"
        )
    print(
        f"GM {hull_data['gm_transverse']:.9f} m == initial_stability "
        f"{hull_data['gm_initial_stability']:.9f} m; slope {hull_data['small_angle_slope']:.6f} "
        f"m/rad ({hull_data['slope_relative_error'] * 100:.2f} % off, 1 deg secant)"
    )
    print(
        f"GZ_max {hull_data['gz_max']:.6f} m at {hull_data['angle_of_gz_max_deg']:.0f} deg, "
        f"vanishing angle {hull_data['vanishing_angle_deg']} deg, "
        f"area {hull_data['dynamic_stability_area']:.6f} m*rad"
    )
    print(
        f"\nWigley + {wigley_data['freeboard']:g} m topsides: GM "
        f"{wigley_data['gm_transverse']:.6f} m, GZ_max {wigley_data['gz_max']:.6f} m at "
        f"{wigley_data['angle_of_gz_max_deg']:.0f} deg"
    )

    payload = {
        "tensorcad_version": __version__,
        "box": box_data,
        "hull": hull_data,
        "wigley": wigley_data,
        "figure": figure,
        "claims": {
            "box_matches_exact_section": box_data["worst_reference_error"] < 1e-11,
            "wall_sided_inside_window": box_data["worst_wall_sided_relative_error_inside"] < 1e-11,
            "wall_sided_departs_outside_window": box_data[
                "first_wall_sided_relative_error_outside"
            ]
            > 1e-3,
            "hull_gm_matches_initial_stability": abs(
                hull_data["gm_transverse"] - hull_data["gm_initial_stability"]
            )
            < 1e-12,
            "hull_slope_matches_gm": hull_data["slope_relative_error"] < 0.01,
        },
    }
    (arguments.output / "results.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"\nwrote {arguments.output / 'results.json'}" + (f" and {figure}" if figure else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
