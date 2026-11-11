#!/usr/bin/env python3
"""Hull-girder section properties: exact geometry against closed forms and a real hull.

Five studies, each with a reference the answer can be wrong against:

1. **Closed forms.** A rectangle, a double-walled box, an I-beam and a plus-shaped
   section are computed both by the polygon sums in :mod:`tensorcad.girder` and by hand
   formulas (``b h^3 / 12``, ``b h^2 / 6``, ``b h^2 / 4``, Steiner sums of rectangles).
   The measured residual is round-off, not discretisation -- there is nothing to converge.
2. **The plastic neutral axis.** For a trapezoid the equal-area level has a closed form,
   and the plastic modulus is compared against a Simpson quadrature of the width law that
   never bisects. This is the only quantity in the module that is searched for.
3. **An inscribed polygon's circle.** The second moment of a regular n-gon approaches
   ``pi R^4 / 4`` with an error coefficient that theory fixes at ``4 pi^2 / 3``; the
   sequence is reported against that prediction rather than against a tolerance.
4. **A real hull.** :func:`tensorcad.generate_hull` at ten stations: the section area is
   cross-checked against the independent clipped integration in :mod:`tensorcad.bonjean`,
   the loft's symmetry shows up as a vanishing product of inertia, and the hull-girder
   second moment and section modulus along the length are reported so the governing
   section can be named.
5. **Torsion.** Bredt's closed-cell constant against the classical box formula, and the
   open-section constant for the same steel -- the number behind "a ship needs a closed
   deck".

Writes ``results.json``, ``report.md`` and ``plot.png`` into ``--output`` and exits
non-zero if any claim stops holding.

Usage::

    python examples/girder_section/run_example.py [--output artifacts/girder_section]
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

from tensorcad import HullParameters, generate_hull
from tensorcad.bonjean import section_loops, section_properties
from tensorcad.girder import (
    box_girder_section,
    closed_cell_torsion_constant,
    girder_section,
    open_section_torsion_constant,
)

LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8


def relative(measured: float, reference: float) -> float:
    if reference == 0.0:
        return abs(measured)
    return abs(measured - reference) / abs(reference)


def simpson(function, low: float, high: float, intervals: int = 200_000) -> float:
    if intervals % 2:
        intervals += 1
    step = (high - low) / intervals
    total = function(low) + function(high)
    for index in range(1, intervals):
        total += (4.0 if index % 2 else 2.0) * function(low + index * step)
    return total * step / 3.0


def i_beam_polygon(
    flange_width: float, flange_thickness: float, web_thickness: float, depth: float
) -> list[tuple[float, float]]:
    half_f, half_w = 0.5 * flange_width, 0.5 * web_thickness
    return [
        (-half_f, 0.0),
        (half_f, 0.0),
        (half_f, flange_thickness),
        (half_w, flange_thickness),
        (half_w, depth - flange_thickness),
        (half_f, depth - flange_thickness),
        (half_f, depth),
        (-half_f, depth),
        (-half_f, depth - flange_thickness),
        (-half_w, depth - flange_thickness),
        (-half_w, flange_thickness),
        (-half_f, flange_thickness),
    ]


def study_closed_forms() -> dict[str, object]:
    beam, depth = 6.0, 4.0
    rectangle = girder_section(box_girder_section(beam, depth))
    thickness = 1.0
    inner_beam, inner_depth = beam - 2.0 * thickness, depth - 2.0 * thickness
    hollow = girder_section(box_girder_section(beam, depth, thickness=thickness))
    i_beam = girder_section(i_beam_polygon(4.0, 0.5, 0.4, 6.0))
    plus = girder_section(
        [
            (-1.0, -3.0),
            (1.0, -3.0),
            (1.0, -1.0),
            (3.0, -1.0),
            (3.0, 1.0),
            (1.0, 1.0),
            (1.0, 3.0),
            (-1.0, 3.0),
            (-1.0, 1.0),
            (-3.0, 1.0),
            (-3.0, -1.0),
            (-1.0, -1.0),
        ]
    )

    references = {
        "rectangle_area": (rectangle.area, beam * depth),
        "rectangle_i_horizontal": (rectangle.second_moment_horizontal, beam * depth**3 / 12.0),
        "rectangle_section_modulus": (rectangle.section_modulus_deck, beam * depth**2 / 6.0),
        "rectangle_plastic_modulus": (rectangle.plastic_modulus, beam * depth**2 / 4.0),
        "rectangle_shape_factor": (rectangle.plastic_to_elastic_ratio, 1.5),
        "hollow_area": (hollow.area, beam * depth - inner_beam * inner_depth),
        "hollow_i_horizontal": (
            hollow.second_moment_horizontal,
            (beam * depth**3 - inner_beam * inner_depth**3) / 12.0,
        ),
        "hollow_plastic_modulus": (
            hollow.plastic_modulus,
            (beam * depth**2 - inner_beam * inner_depth**2) / 4.0,
        ),
        "i_beam_area": (i_beam.area, 4.0 * 0.5 + 0.4 * 5.0 + 4.0 * 0.5),
        "plus_area": (plus.area, 2.0 * 6.0 * 2.0 - 4.0),
        "plus_plastic_modulus": (plus.plastic_modulus, 2.0 * 3.0 * 1.0 + 2.0 * 1.0 * 9.0 - 2.0),
    }
    measured = {name: value for name, (value, _in) in references.items()}
    errors = {name: relative(value, reference) for name, (value, reference) in references.items()}
    return {"measured": measured, "reference": {n: r for n, (_v, r) in references.items()}, "relative_error": errors}


def study_plastic_level() -> dict[str, object]:
    """Trapezoid: equal-area level in closed form, plastic modulus by quadrature."""
    rows = []
    for bottom, top, depth in ((2.0, 6.0, 4.0), (6.0, 2.0, 3.0), (1.0, 8.0, 5.0)):
        polygon = [
            (-0.5 * bottom, 0.0),
            (0.5 * bottom, 0.0),
            (0.5 * top, depth),
            (-0.5 * top, depth),
        ]
        section = girder_section(polygon)
        area = 0.5 * (bottom + top) * depth
        slope = (top - bottom) / (2.0 * depth)
        level = (
            0.5 * depth
            if slope == 0.0
            else (-bottom + math.sqrt(bottom * bottom + 4.0 * slope * area / 2.0)) / (2.0 * slope)
        )
        def width(z: float, b: float = bottom, t: float = top, d: float = depth) -> float:
            return b + (t - b) * z / d

        below = simpson(lambda z, lv=level, w=width: (lv - z) * w(z), 0.0, level, 100_000)
        above = simpson(lambda z, lv=level, w=width: (z - lv) * w(z), level, depth, 100_000)
        plastic = below + above
        rows.append(
            {
                "bottom_width": bottom,
                "top_width": top,
                "depth": depth,
                "level_measured": section.plastic_neutral_axis_z,
                "level_reference": level,
                "plastic_modulus_measured": section.plastic_modulus,
                "plastic_modulus_reference": plastic,
                "level_relative_error": relative(section.plastic_neutral_axis_z, level),
                "plastic_modulus_relative_error": relative(section.plastic_modulus, plastic),
            }
        )
    return {"rows": rows}


def study_circle() -> dict[str, object]:
    radius = 2.0
    exact = math.pi * radius**4 / 4.0
    rows = []
    for sides in (16, 64, 256, 1024):
        polygon = [
            (
                radius * math.cos(2.0 * math.pi * k / sides),
                radius * math.sin(2.0 * math.pi * k / sides),
            )
            for k in range(sides)
        ]
        section = girder_section(polygon)
        error = relative(section.second_moment_horizontal, exact)
        rows.append(
            {
                "sides": sides,
                "second_moment": section.second_moment_horizontal,
                "relative_error": error,
                "error_times_n_squared": error * sides * sides,
            }
        )
    return {"exact": exact, "rows": rows, "predicted_coefficient": 4.0 * math.pi**2 / 3.0}


def study_hull() -> dict[str, object]:
    model = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=41, nz=13))
    stations = [LENGTH * fraction for fraction in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)]
    rows = []
    for station in stations:
        loops = section_loops(model.mesh, station)
        section = girder_section(loops[0])
        bonjean_area = section_properties(model.mesh, station).area
        rows.append(
            {
                "station": station,
                "area": section.area,
                "area_bonjean": bonjean_area,
                "area_relative_gap": relative(section.area, bonjean_area),
                "centroid_z": section.centroid_z,
                "second_moment_horizontal": section.second_moment_horizontal,
                "section_modulus_deck": section.section_modulus_deck,
                "section_modulus_keel": section.section_modulus_keel,
                "product_of_inertia": section.product_of_inertia,
                "plastic_modulus": section.plastic_modulus,
                "plastic_to_elastic_ratio": section.plastic_to_elastic_ratio,
                "beam": section.beam,
                "depth": section.depth,
            }
        )
    midship = max(rows, key=lambda row: row["area"])
    levels = [midship["centroid_z"] + offset * midship["depth"] for offset in (-0.5, -0.25, 0.0, 0.25, 0.5)]
    loops = section_loops(model.mesh, midship["station"])
    section = girder_section(loops[0])
    static_moments = [
        {"z": level, "Q": section.first_moment(level)} for level in levels
    ]
    return {
        "length": LENGTH,
        "beam": BEAM,
        "draft": DRAFT,
        "rows": rows,
        "midship": midship,
        "midship_points": [list(point) for point in section.points],
        "midship_static_moments": static_moments,
        "max_area_gap": max(row["area_relative_gap"] for row in rows),
        "max_abs_product": max(abs(row["product_of_inertia"]) for row in rows),
    }


def study_torsion() -> dict[str, object]:
    beam, depth, thickness = 6.0, 4.0, 0.02
    enclosed = (beam - thickness) * (depth - thickness)
    perimeter = 2.0 * ((beam - thickness) + (depth - thickness))
    bredt = closed_cell_torsion_constant(enclosed, thickness, perimeter)
    classical = (
        2.0
        * thickness
        * (beam - thickness) ** 2
        * (depth - thickness) ** 2
        / ((beam - thickness) + (depth - thickness))
    )
    open_strips = [(beam, thickness), (depth, thickness), (beam, thickness)]
    open_constant = open_section_torsion_constant(open_strips)
    return {
        "closed_cell": bredt,
        "closed_cell_reference": classical,
        "closed_cell_relative_error": relative(bredt, classical),
        "open_section": open_constant,
        "closed_to_open_ratio": bredt / open_constant,
    }


def plot(results: dict[str, object], output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    hull = results["hull"]
    figure, axes = plt.subplots(2, 2, figsize=(11.0, 9.0))

    # (a) the midship section
    points = hull["midship_points"]
    ys = [p[0] for p in points]
    zs = [p[1] for p in points]
    axes[0][0].fill([*ys, ys[0]], [*zs, zs[0]], color="#9fc6e0", edgecolor="#1f4e79")
    axes[0][0].axhline(hull["midship"]["centroid_z"], color="#c0392b", linestyle="--", linewidth=1.0)
    axes[0][0].axhline(
        hull["midship"]["centroid_z"], color="#27864a", linestyle=":", linewidth=1.0
    )
    axes[0][0].set_title("Midship section (hull girder)")
    axes[0][0].set_xlabel("y [m]")
    axes[0][0].set_ylabel("z [m]")
    axes[0][0].set_aspect("equal")
    axes[0][0].grid(alpha=0.3)

    # (b) second moment and section modulus along the hull
    stations = [row["station"] for row in hull["rows"]]
    axes[0][1].plot(stations, [row["second_moment_horizontal"] for row in hull["rows"]], "o-", label="$I_{hh}$")
    axes[0][1].plot(stations, [row["section_modulus_deck"] for row in hull["rows"]], "s-", label="$Z_{deck}$")
    axes[0][1].plot(stations, [row["section_modulus_keel"] for row in hull["rows"]], "^-", label="$Z_{keel}$")
    axes[0][1].set_title("Hull-girder stiffness along the length")
    axes[0][1].set_xlabel("station x [m]")
    axes[0][1].set_ylabel("m$^4$ / m$^3$")
    axes[0][1].legend()
    axes[0][1].grid(alpha=0.3)

    # (c) static moment profile: computed vs the rectangle's closed form scaled out
    profile_z = [row["z"] for row in hull["midship_static_moments"]]
    profile_q = [row["Q"] for row in hull["midship_static_moments"]]
    axes[1][0].plot(profile_q, profile_z, "o-", color="#8e44ad")
    axes[1][0].set_title("Static moment $Q(z)$ at midship")
    axes[1][0].set_xlabel("Q [m$^3$]")
    axes[1][0].set_ylabel("z [m]")
    axes[1][0].grid(alpha=0.3)

    # (d) circle convergence against the 1/n^2 prediction
    circle = results["circle"]
    sides = [row["sides"] for row in circle["rows"]]
    errors = [row["relative_error"] for row in circle["rows"]]
    axes[1][1].loglog(sides, errors, "o-", label="measured")
    axes[1][1].loglog(
        sides,
        [circle["predicted_coefficient"] / n**2 for n in sides],
        "--",
        label=r"$4\pi^2/3 \, n^{-2}$",
    )
    axes[1][1].set_title("Inscribed n-gon vs $\\pi R^4/4$")
    axes[1][1].set_xlabel("sides n")
    axes[1][1].set_ylabel("relative error")
    axes[1][1].legend()
    axes[1][1].grid(alpha=0.3, which="both")

    figure.suptitle("TensorCAD hull-girder section properties")
    figure.tight_layout()
    figure.savefig(output / "plot.png", dpi=150)
    plt.close(figure)


def write_report(results: dict[str, object], output: Path) -> None:
    closed = results["closed_forms"]
    circle = results["circle"]
    hull = results["hull"]
    torsion = results["torsion"]
    plastic = results["plastic"]
    worst = max(closed["relative_error"].items(), key=lambda item: item[1])
    lines = [
        "# Hull-girder section properties -- measured",
        "",
        "## Closed forms",
        "",
        "| quantity | measured | reference | relative error |",
        "|---|---|---|---|",
    ]
    for name, value in closed["measured"].items():
        reference = closed["reference"][name]
        error = closed["relative_error"][name]
        lines.append(f"| {name} | {value:.12g} | {reference:.12g} | {error:.3e} |")
    lines += [
        "",
        f"Worst closed-form residual: **{worst[0]}** at {worst[1]:.3e} relative.",
        "",
        "## Plastic neutral axis (trapezoids)",
        "",
        "| bottom | top | depth | level (measured) | level (closed form) | Zp (measured) | Zp (quadrature) |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in plastic["rows"]:
        lines.append(
            f"| {row['bottom_width']} | {row['top_width']} | {row['depth']} | "
            f"{row['level_measured']:.9f} | {row['level_reference']:.9f} | "
            f"{row['plastic_modulus_measured']:.9f} | {row['plastic_modulus_reference']:.9f} |"
        )
    lines += [
        "",
        "## Inscribed n-gon",
        "",
        "| sides | I | relative error | error * n^2 |",
        "|---|---|---|---|",
    ]
    for row in circle["rows"]:
        lines.append(
            f"| {row['sides']} | {row['second_moment']:.9f} | {row['relative_error']:.3e} | "
            f"{row['error_times_n_squared']:.4f} |"
        )
    lines += [
        "",
        f"Predicted coefficient `4 pi^2 / 3` = {circle['predicted_coefficient']:.4f}.",
        "",
        "## Hull (L=10 m, B=2 m, T=0.8 m)",
        "",
        "| station | A | A (bonjean) | I_hh | Z_deck | Z_keel | Zp | Zp/Z |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in hull["rows"]:
        lines.append(
            f"| {row['station']:.2f} | {row['area']:.6f} | {row['area_bonjean']:.6f} | "
            f"{row['second_moment_horizontal']:.6f} | {row['section_modulus_deck']:.6f} | "
            f"{row['section_modulus_keel']:.6f} | {row['plastic_modulus']:.6f} | "
            f"{row['plastic_to_elastic_ratio']:.4f} |"
        )
    lines += [
        "",
        f"Largest area gap against bonjean: {hull['max_area_gap']:.3e} relative. "
        f"Largest |product of inertia|: {hull['max_abs_product']:.3e} m^4.",
        f"Governing (largest) section at x = {hull['midship']['station']:.2f} m.",
        "",
        "## Torsion",
        "",
        f"Bredt closed cell {torsion['closed_cell']:.6e} m^4 vs classical box "
        f"{torsion['closed_cell_reference']:.6e} m^4 "
        f"({torsion['closed_cell_relative_error']:.3e} relative); "
        f"open-steel version {torsion['open_section']:.6e} m^4, "
        f"ratio {torsion['closed_to_open_ratio']:.1f}.",
        "",
    ]
    (output / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(REPO_ROOT / "artifacts" / "girder_section"))
    arguments = parser.parse_args()
    output = Path(arguments.output)
    output.mkdir(parents=True, exist_ok=True)

    results: dict[str, object] = {
        "closed_forms": study_closed_forms(),
        "plastic": study_plastic_level(),
        "circle": study_circle(),
        "hull": study_hull(),
        "torsion": study_torsion(),
    }

    failures: list[str] = []
    worst_closed = max(results["closed_forms"]["relative_error"].values())
    if worst_closed > 1e-12:
        failures.append(f"closed-form residual {worst_closed:.3e} above 1e-12")
    for row in results["plastic"]["rows"]:
        if row["level_relative_error"] > 1e-12:
            failures.append(f"equal-area level off by {row['level_relative_error']:.3e}")
        if row["plastic_modulus_relative_error"] > 1e-6:
            failures.append(f"plastic modulus off by {row['plastic_modulus_relative_error']:.3e}")
    if results["hull"]["max_area_gap"] > 1e-9:
        failures.append(f"hull section area differs from bonjean by {results['hull']['max_area_gap']:.3e}")
    if results["hull"]["max_abs_product"] > 1e-9:
        failures.append(f"hull section product of inertia {results['hull']['max_abs_product']:.3e} not zero")
    if results["torsion"]["closed_cell_relative_error"] > 1e-12:
        failures.append("Bredt constant disagrees with the box formula")
    last = results["circle"]["rows"][-1]
    if abs(last["error_times_n_squared"] - results["circle"]["predicted_coefficient"]) > 0.01:
        failures.append("circle convergence coefficient drifted from 4 pi^2 / 3")

    (output / "results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True), encoding="utf-8"
    )
    write_report(results, output)
    plot(results, output)

    print(f"wrote {output / 'results.json'}, report.md and plot.png")
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
