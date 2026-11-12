"""Independent re-check of :mod:`tensorcad.section_gz` (merged 2026-11-12, mission m16358).

The module ships its own evidence: exact vs mesh agreement, wall-sided agreement, chord
convergence.  This run asks the *other* question a verifier seat has to ask -- does an
independently written route, which shares no algorithm with the module, reproduce the same
arms?  :mod:`tests.section_gz_oracle` supplies that route: horizontal chord slices integrated
by Simpson quadrature for the ellipse (the conic is recovered numerically, the disk formula is
never used) and exact slice integration for polygons (no clipping, no Green's theorem).

Every number below is computed in this run and written to ``outputs/metrics.json``; the
markdown table in ``outputs/report.md`` is generated from that JSON.

Usage::

    python examples/section_gz_independent/run_example.py [--out DIR]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:  # a clean checkout has no installed package
    sys.path.insert(0, str(REPO_ROOT))
if str(REPO_ROOT / "tests") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "tests"))

from section_gz_oracle import (  # noqa: E402
    EllipseSliceOracle,
    ellipse_prism_arms_independent,
    polygon_prism_arms_independent,
    polygon_properties_below,
    wall_sided_arm,
    wall_sided_limit_deg,
    waterplane_gm,
)
from tensorcad import section_gz  # noqa: E402
from tensorcad.gz import GzError, righting_arm_curve  # noqa: E402

ANGLES = [float(angle) for angle in range(0, 91, 5)]
TRIANGLE = ((-1.0, 0.0), (1.0, 0.0), (0.0, 2.0))
ELLIPSE = {"semi_beam": 2.0, "semi_height": 1.0, "draft": 1.2, "kg": 0.8, "length": 4.0}
BOX = {"breadth": 4.0, "draft": 1.0, "depth": 2.0, "kg": 0.6, "length": 5.0}
RECTANGLE = ((-2.0, 0.0), (2.0, 0.0), (2.0, 2.0), (-2.0, 2.0))


def _max_abs(left, right) -> float:
    return max(abs(a - b) for a, b in zip(left, right, strict=True))


def _max_rel(left, right) -> float:
    return max(abs(a - b) / max(abs(b), 1e-12) for a, b in zip(left, right, strict=True))


def case_ellipse() -> dict[str, object]:
    """2:1 ellipse: closed form vs slice quadrature, both heel directions, GM and landmarks."""
    a, b = ELLIPSE["semi_beam"], ELLIPSE["semi_height"]
    draft, kg, length = ELLIPSE["draft"], ELLIPSE["kg"], ELLIPSE["length"]
    library = section_gz.ellipse_prism_righting_arms(
        a, b, ANGLES, kg=kg, waterline_z=draft, length=length
    )
    independent = ellipse_prism_arms_independent(a, b, ANGLES, kg=kg, waterline_z=draft)
    both = ANGLES + [-angle for angle in ANGLES[1:]]
    published = _published_ellipse()
    library_both = section_gz.ellipse_prism_righting_arms(
        a, b, both, kg=kg, waterline_z=draft, length=length
    ).gz
    independent_both = ellipse_prism_arms_independent(a, b, both, kg=kg, waterline_z=draft)

    # GM by three routes: waterplane moments integrated here, the library's closed-form
    # slope, and a finite difference of the independent arms.  The GM the merged example
    # *published* is read back out of its committed artifact, not retyped.
    oracle = EllipseSliceOracle(a, b, 0.0)
    level = draft
    area, _, kb = oracle.properties_below(level)
    inertia = 2.0 * (0.5 * oracle.width(level)) ** 3 / 3.0
    gm_moments = kb + inertia / area - kg
    kb_library = section_gz.ellipse_prism_righting_arms(
        a, b, [0.0], kg=kg, waterline_z=draft, length=length
    ).arms[0].centre_z
    tiny = 1e-6
    slope_library = (
        section_gz.ellipse_prism_righting_arms(
            a, b, [tiny], kg=kg, waterline_z=draft, length=length
        ).gz[0]
        / math.radians(tiny)
    )
    delta = 0.5
    plus = ellipse_prism_arms_independent(a, b, [delta], kg=kg, waterline_z=draft)[0]
    minus = ellipse_prism_arms_independent(a, b, [-delta], kg=kg, waterline_z=draft)[0]
    slope_independent = (plus - minus) / (2.0 * math.radians(delta))

    # landmarks: GZ_max on a 1 deg grid, and the vanishing angle at a raised KG
    fine = [float(angle) for angle in range(0, 91)]
    fine_arms = ellipse_prism_arms_independent(a, b, fine, kg=kg, waterline_z=draft)
    best = max(zip(fine_arms, fine))
    high_kg = 1.6
    grid = [float(angle) for angle in range(0, 91, 2)]
    high_arms = ellipse_prism_arms_independent(a, b, grid, kg=high_kg, waterline_z=draft)
    library_high = section_gz.ellipse_prism_righting_arms(
        a, b, grid, kg=high_kg, waterline_z=draft, length=length
    )
    library_vanishing = library_high.vanishing_stability_angle_deg()
    vanishing = None
    for left, right, arm_left, arm_right in zip(
        grid, grid[1:], high_arms, high_arms[1:], strict=False
    ):
        if arm_left > 0.0 >= arm_right:
            vanishing = left + arm_left / (arm_left - arm_right) * (right - left)
            break

    # mesh route convergence measured against the *independent* oracle, not the closed form
    convergence = []
    for sides in (16, 32, 64, 128, 256):
        mesh = section_gz.prism_mesh(section_gz.ellipse_section(a, b, sides), length)
        mesh_arms = righting_arm_curve(mesh, kg, ANGLES, waterline_z=draft).gz
        convergence.append(
            {
                "sides": sides,
                "max_abs_error_m": _max_abs(mesh_arms, independent),
                "max_rel_error": _max_rel(mesh_arms, independent),
            }
        )
    for index in range(1, len(convergence)):
        ratio = convergence[index - 1]["max_abs_error_m"] / convergence[index]["max_abs_error_m"]
        convergence[index]["observed_order"] = math.log2(ratio)

    return {
        "name": "2:1 ellipse a=2 m b=1 m, draft 1.2 m, KG 0.8 m, L=4 m",
        "max_abs_error_m": _max_abs(library.gz, independent),
        "max_rel_error": _max_rel(library.gz, independent),
        "max_abs_error_both_heels_m": _max_abs(library_both, independent_both),
        "symmetry_residual_m": _max_abs(
            independent_both[len(ANGLES) :], [-value for value in independent[1:]]
        ),
        "gm_published_m": published.get("gm_from_waterplane_m"),
        "gm_waterplane_independent_m": waterplane_gm(
            0.5 * oracle.width(level), area, kb, kg
        ),
        "gm_integrated_independent_m": gm_moments,
        "gm_slope_library_m": slope_library,
        "gm_slope_independent_m": slope_independent,
        "kb_library_m": kb_library,
        "kb_independent_m": kb,
        "gz_max_published_m": published.get("gz_max_m"),
        "gz_max_published_heel_deg": published.get("gz_max_heel_deg"),
        "gz_max_independent_m": max(independent),
        "gz_max_independent_heel_deg": ANGLES[independent.index(max(independent))],
        "gz_max_refined_m": best[0],
        "gz_max_refined_heel_deg": best[1],
        "vanishing_angle_library_deg": library_vanishing,
        "vanishing_angle_published_deg": published.get("vanishing_angle_at_kg_1p6_deg"),
        "vanishing_angle_independent_deg": vanishing,
        "table": [
            {"heel_deg": angle, "gz_library": l, "gz_independent": i, "diff": l - i}
            for angle, l, i in zip(ANGLES, library.gz, independent, strict=True)
        ],
        "convergence": convergence,
        "wall_sided_relative_error_at_30": abs(
            wall_sided_arm(a, draft, kg, 30.0) - library.angle_of(30.0).righting_arm
        )
        / abs(library.angle_of(30.0).righting_arm),
    }


def _published_ellipse() -> dict[str, float]:
    """Read the GM / GZ_max / vanishing angle that the merged example committed."""
    path = Path(__file__).resolve().parents[1] / "section_gz_exact" / "outputs" / "results.json"
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError):
        return {}
    return dict(data.get("results", {}).get("ellipse", {}))


def case_box() -> dict[str, object]:
    """Box barge: independent slice oracle vs the polygon oracle, and vs the wall-sided form."""
    breadth, draft, kg = BOX["breadth"], BOX["draft"], BOX["kg"]
    library = section_gz.polygon_prism_righting_arms(
        RECTANGLE, ANGLES, kg=kg, waterline_z=draft, length=BOX["length"]
    )
    independent = polygon_prism_arms_independent(
        RECTANGLE, ANGLES, kg=kg, waterline_z=draft
    )
    limit = wall_sided_limit_deg(breadth, draft)
    rows = []
    for angle, library_arm, independent_arm in zip(
        ANGLES, library.gz, independent, strict=True
    ):
        wall = wall_sided_arm(breadth, draft, kg, angle)
        rows.append(
            {
                "heel_deg": angle,
                "gz_library": library_arm,
                "gz_independent": independent_arm,
                "gz_wall_sided": wall,
                "library_minus_wall_sided": library_arm - wall,
                "independent_minus_wall_sided": independent_arm - wall,
                "wall_sided_region": "inside" if angle <= limit else "outside",
            }
        )
    inside = [row for row in rows if row["wall_sided_region"] == "inside"]
    outside = [row for row in rows if row["wall_sided_region"] == "outside"]
    return {
        "name": "box barge B=4 m, T=1 m, D=2 m, KG=0.6 m, L=5 m",
        "wall_sided_limit_deg": limit,
        "max_abs_error_vs_independent_m": _max_abs(library.gz, independent),
        "max_abs_error_wall_sided_inside_m": max(
            abs(row["library_minus_wall_sided"]) for row in inside
        ),
        "min_abs_error_wall_sided_outside_m": min(
            abs(row["library_minus_wall_sided"]) for row in outside
        ),
        "transom_heel_deg": 90.0,
        "gz_at_90_independent_m": independent[-1],
        "table": rows,
    }


def case_polygons() -> dict[str, object]:
    """Two polygonal sections the library uses in its own evidence, re-derived independently."""
    cases = {
        "v_section": {
            "section": TRIANGLE,
            "kg": 0.6,
            "waterline": 1.0,
            "angles": ANGLES,
            "length": 4.0,
            "name": "isosceles V section B=2 m, height 2 m, draft 1 m, KG 0.6 m",
        },
        "rotated_hexagon": {
            "section": section_gz.regular_polygon_section(6, 1.0, rotation_deg=17.0),
            "kg": 0.7,
            "waterline": 0.8,
            "angles": [float(angle) for angle in range(0, 46, 5)],
            "length": 2.5,
            "name": "regular hexagon R=1 m rotated 17 deg, draft 0.8 m, KG 0.7 m",
        },
    }
    out: dict[str, object] = {}
    for key, case in cases.items():
        library = section_gz.polygon_prism_righting_arms(
            case["section"],
            case["angles"],
            kg=case["kg"],
            waterline_z=case["waterline"],
            length=case["length"],
        )
        independent = polygon_prism_arms_independent(
            case["section"], case["angles"], kg=case["kg"], waterline_z=case["waterline"]
        )
        out[key] = {
            "name": case["name"],
            "max_abs_error_m": _max_abs(library.gz, independent),
            "max_rel_error": _max_rel(library.gz, independent),
            "gz_library": list(library.gz),
            "gz_independent": independent,
            "heels_deg": case["angles"],
        }
    return out


def case_inscribed_polygon() -> dict[str, object]:
    """Does ``ellipse_section`` really return a polygon inside the ellipse, for every `sides`?

    The defect this run pins: the closed-form doubling of the sample set once left one sample
    family without a vertex at the keel, and the ``- min(z)`` shift then pushed the boundary
    outside the exact section (up to 3.0 % for ``sides = 18``).
    """
    a, b = ELLIPSE["semi_beam"], ELLIPSE["semi_height"]
    rows = []
    for sides in (5, 16, 17, 18, 20, 32):
        points = section_gz.ellipse_section(a, b, sides)
        worst = 0.0
        count = len(points)
        for index in range(count):
            y0, z0 = points[index]
            y1, z1 = points[(index + 1) % count]
            for step in range(21):
                fraction = step / 20.0
                y = y0 + fraction * (y1 - y0)
                z = z0 + fraction * (z1 - z0)
                worst = max(worst, (y / a) ** 2 + ((z - b) / b) ** 2)
        area = polygon_properties_below(points, max(z for (_, z) in points))[0]
        rows.append(
            {
                "sides": sides,
                "keel_offset_m": min(z for (_, z) in points),
                "outside_ellipse_quadric_max": worst,
                "area_ratio_to_ellipse": area / (math.pi * a * b),
            }
        )
    return {
        "rows": rows,
        "worst_outside": max(row["outside_ellipse_quadric_max"] for row in rows),
        "worst_overfill": max(row["area_ratio_to_ellipse"] for row in rows),
    }


def case_rejections() -> dict[str, object]:
    """The input contract: every rejection the docstring promises, exercised."""
    concave = ((0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0))
    probes = {
        "concave_section": lambda: section_gz.polygon_prism_righting_arms(
            concave, [10.0], kg=0.5, waterline_z=0.5
        ),
        "zero_area_section": lambda: section_gz.polygon_prism_righting_arms(
            ((0.0, 0.0), (1.0, 0.0), (2.0, 0.0)), [10.0], kg=0.5, waterline_z=0.5
        ),
        "non_finite_vertex": lambda: section_gz.polygon_prism_righting_arms(
            ((0.0, 0.0), (1.0, float("nan")), (0.0, 1.0)), [10.0], kg=0.5, waterline_z=0.5
        ),
        "heel_at_180": lambda: section_gz.polygon_prism_righting_arms(
            TRIANGLE, [180.0], kg=0.5, waterline_z=1.0
        ),
        "waterline_at_keel": lambda: section_gz.polygon_prism_righting_arms(
            TRIANGLE, [10.0], kg=0.5, waterline_z=0.0
        ),
        "displaced_area_too_large": lambda: section_gz.polygon_prism_righting_arms(
            TRIANGLE, [10.0], kg=0.5, waterline_z=1.0, displaced_area=5.0
        ),
        "negative_length": lambda: section_gz.prism_mesh(TRIANGLE, -1.0),
        "ellipse_draft_at_top": lambda: section_gz.ellipse_prism_righting_arms(
            2.0, 1.0, [10.0], kg=0.8, waterline_z=2.0
        ),
        "negative_semi_axis": lambda: section_gz.ellipse_prism_righting_arms(
            -2.0, 1.0, [10.0], kg=0.8, waterline_z=1.2
        ),
    }
    results = {}
    for name, probe in probes.items():
        try:
            probe()
            results[name] = {"rejected": False, "error": None}
        except GzError as error:
            results[name] = {"rejected": True, "error": str(error)}
    # clockwise winding must be accepted and give the same arms
    clockwise = tuple(reversed(TRIANGLE))
    forward = section_gz.polygon_prism_righting_arms(
        TRIANGLE, [20.0], kg=0.6, waterline_z=1.0
    ).gz[0]
    reversed_arm = section_gz.polygon_prism_righting_arms(
        clockwise, [20.0], kg=0.6, waterline_z=1.0
    ).gz[0]
    return {"probes": results, "clockwise_winding_difference_m": abs(forward - reversed_arm)}


def run(out_dir: Path) -> dict[str, object]:
    ellipse = case_ellipse()
    box = case_box()
    polygons = case_polygons()
    polygon = case_inscribed_polygon()
    rejections = case_rejections()

    checks = [
        {
            "check": "ellipse: closed form vs independent slice quadrature, 0-90 deg",
            "value": ellipse["max_abs_error_m"],
            "unit": "m",
            "criterion": "< 1e-9 m",
            "passed": ellipse["max_abs_error_m"] < 1e-9,
        },
        {
            "check": "ellipse: both heel directions vs independent",
            "value": ellipse["max_abs_error_both_heels_m"],
            "unit": "m",
            "criterion": "< 1e-9 m",
            "passed": ellipse["max_abs_error_both_heels_m"] < 1e-9,
        },
        {
            "check": "box / V / rotated hexagon: polygon oracle vs independent",
            "value": max(
                box["max_abs_error_vs_independent_m"],
                polygons["v_section"]["max_abs_error_m"],
                polygons["rotated_hexagon"]["max_abs_error_m"],
            ),
            "unit": "m",
            "criterion": "< 1e-12 m",
            "passed": max(
                box["max_abs_error_vs_independent_m"],
                polygons["v_section"]["max_abs_error_m"],
                polygons["rotated_hexagon"]["max_abs_error_m"],
            )
            < 1e-12,
        },
        {
            "check": "box: wall-sided formula reproduced inside its limit (independent formula)",
            "value": box["max_abs_error_wall_sided_inside_m"],
            "unit": "m",
            "criterion": "< 1e-12 m",
            "passed": box["max_abs_error_wall_sided_inside_m"] < 1e-12,
        },
        {
            "check": "box: wall-sided formula genuinely fails outside its limit",
            "value": box["min_abs_error_wall_sided_outside_m"],
            "unit": "m",
            "criterion": "> 0.01 m (non-vacuous)",
            "passed": box["min_abs_error_wall_sided_outside_m"] > 0.01,
        },
        {
            "check": "ellipse GM: published artifact vs independent waterplane moments",
            "value": abs(
                (ellipse["gm_published_m"] or 0.0) - ellipse["gm_waterplane_independent_m"]
            ),
            "unit": "m",
            "criterion": "< 1e-6 m",
            "passed": abs(
                (ellipse["gm_published_m"] or 0.0) - ellipse["gm_waterplane_independent_m"]
            )
            < 1e-6,
        },
        {
            "check": "ellipse GM: library closed-form slope vs independent waterplane moments",
            "value": abs(ellipse["gm_slope_library_m"] - ellipse["gm_waterplane_independent_m"]),
            "unit": "m",
            "criterion": "< 1e-6 m",
            "passed": abs(
                ellipse["gm_slope_library_m"] - ellipse["gm_waterplane_independent_m"]
            )
            < 1e-6,
        },
        {
            "check": "ellipse GM: independent finite difference vs independent waterplane",
            "value": abs(ellipse["gm_slope_independent_m"] - ellipse["gm_waterplane_independent_m"]),
            "unit": "m",
            "criterion": "< 1e-3 m (finite-difference truncation)",
            "passed": abs(
                ellipse["gm_slope_independent_m"] - ellipse["gm_waterplane_independent_m"]
            )
            < 1e-3,
        },
        {
            "check": "ellipse KB: library centroid vs independent centroid",
            "value": abs(ellipse["kb_library_m"] - ellipse["kb_independent_m"]),
            "unit": "m",
            "criterion": "< 1e-9 m",
            "passed": abs(ellipse["kb_library_m"] - ellipse["kb_independent_m"]) < 1e-9,
        },
        {
            "check": "ellipse GZ_max: published artifact vs independent 5 deg grid",
            "value": abs(
                (ellipse["gz_max_published_m"] or 0.0) - ellipse["gz_max_independent_m"]
            ),
            "unit": "m",
            "criterion": "< 1e-6 m",
            "passed": abs(
                (ellipse["gz_max_published_m"] or 0.0) - ellipse["gz_max_independent_m"]
            )
            < 1e-6,
        },
        {
            "check": "ellipse vanishing angle at KG=1.6 m: library vs independent",
            "value": abs(
                (ellipse["vanishing_angle_library_deg"] or 0.0)
                - (ellipse["vanishing_angle_independent_deg"] or 0.0)
            ),
            "unit": "deg",
            "criterion": "< 0.05 deg (2 deg sampling + linear interpolation)",
            "passed": abs(
                (ellipse["vanishing_angle_library_deg"] or 0.0)
                - (ellipse["vanishing_angle_independent_deg"] or 0.0)
            )
            < 0.05,
        },
        {
            "check": "mesh route converges to the independent oracle at second order",
            "value": ellipse["convergence"][-1]["observed_order"],
            "unit": "-",
            "criterion": "1.9 < order < 2.1",
            "passed": 1.9 < ellipse["convergence"][-1]["observed_order"] < 2.1,
        },
        {
            "check": "ellipse_section is inscribed in the exact ellipse for every side count",
            "value": polygon["worst_outside"],
            "unit": "-",
            "criterion": "<= 1 + 1e-12 (quadric)",
            "passed": polygon["worst_outside"] <= 1.0 + 1e-12,
        },
        {
            "check": "input contract: every documented rejection fires",
            "value": sum(1 for row in rejections["probes"].values() if row["rejected"]),
            "unit": "of",
            "criterion": f"{len(rejections['probes'])} probes",
            "passed": all(row["rejected"] for row in rejections["probes"].values()),
        },
        {
            "check": "clockwise winding accepted and identical",
            "value": rejections["clockwise_winding_difference_m"],
            "unit": "m",
            "criterion": "< 1e-15 m",
            "passed": rejections["clockwise_winding_difference_m"] < 1e-15,
        },
    ]

    metrics = {
        "mission": "m16430",
        "scope": "independent verification of tensorcad.section_gz (merged m16358)",
        "ellipse": ellipse,
        "box": box,
        "polygons": polygons,
        "inscribed_polygon": polygon,
        "rejections": rejections,
        "checks": checks,
        "all_passed": all(row["passed"] for row in checks),
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _write_figure(metrics, out_dir / "plot.png")
    (out_dir / "report.md").write_text(_report(metrics) + "\n")
    return metrics


def _write_figure(metrics: dict[str, object], path: Path) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # pragma: no cover - the numbers are the evidence, the figure is a view
        return
    ellipse = metrics["ellipse"]
    box = metrics["box"]
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    heels = [row["heel_deg"] for row in ellipse["table"]]
    left.plot(heels, [row["gz_library"] for row in ellipse["table"]], "k-", lw=2, label="ellipse: closed form (library)")
    left.plot(
        heels,
        [row["gz_independent"] for row in ellipse["table"]],
        "o",
        ms=4,
        mfc="none",
        label="ellipse: independent slice quadrature",
    )
    for sides, style in ((32, ":"), (256, "--")):
        row = next(entry for entry in ellipse["convergence"] if entry["sides"] == sides)
        left.plot([], [], style, label=f"mesh route, {sides} chords (err {row['max_abs_error_m']:.1e} m)")
    box_heels = [row["heel_deg"] for row in box["table"]]
    left.plot(box_heels, [row["gz_library"] for row in box["table"]], "C3-", lw=1.5, label="box: exact polygon oracle")
    left.plot(
        box_heels,
        [row["gz_wall_sided"] for row in box["table"]],
        "C3--",
        lw=1.0,
        label="box: wall-sided formula",
    )
    left.axvline(box["wall_sided_limit_deg"], color="0.6", lw=0.8, ls=":")
    left.set_xlabel("heel angle [deg]")
    left.set_ylabel("GZ [m]")
    left.set_title("exact arms and an independent route")
    left.grid(alpha=0.3)
    left.legend(fontsize=7)

    sides = [row["sides"] for row in ellipse["convergence"]]
    errors = [row["max_abs_error_m"] for row in ellipse["convergence"]]
    right.loglog(sides, errors, "ko-", label="mesh vs independent oracle")
    reference = [errors[0] * (sides[0] / value) ** 2 for value in sides]
    right.loglog(sides, reference, "--", color="0.6", label="slope -2 reference")
    for row in ellipse["convergence"]:
        right.annotate(
            f"{row['observed_order']:.2f}" if "observed_order" in row else "",
            (row["sides"], row["max_abs_error_m"]),
            textcoords="offset points",
            xytext=(4, 4),
            fontsize=7,
        )
    right.set_xlabel("chords on the inscribed ellipse")
    right.set_ylabel("max |GZ_mesh - GZ_independent| [m]")
    right.set_title("chord convergence against an independent reference")
    right.grid(alpha=0.3, which="both")
    right.legend(fontsize=7)
    figure.tight_layout()
    figure.savefig(path, dpi=130)
    plt.close(figure)


def _report(metrics: dict[str, object]) -> str:
    ellipse = metrics["ellipse"]
    box = metrics["box"]
    polygon = metrics["inscribed_polygon"]
    lines = [
        "# Independent re-check of `tensorcad.section_gz` (mission m16430)",
        "",
        "Generated by `examples/section_gz_independent/run_example.py`; numbers are read from",
        "`outputs/metrics.json`, nothing is typed by hand.",
        "",
        "## Checks",
        "",
        "| check | value | unit | criterion | result |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in metrics["checks"]:
        value = row["value"]
        shown = f"{value:.3e}" if isinstance(value, float) else str(value)
        lines.append(
            f"| {row['check']} | {shown} | {row['unit']} | {row['criterion']} | "
            f"{'PASS' if row['passed'] else 'FAIL'} |"
        )
    lines += [
        "",
        f"All checks passed: **{metrics['all_passed']}**",
        "",
        "## Ellipse 2:1 (a = 2 m, b = 1 m, draft 1.2 m, KG = 0.8 m, L = 4 m)",
        "",
        "| heel [deg] | GZ closed form [m] | GZ independent [m] | difference [m] |",
        "| --- | --- | --- | --- |",
    ]
    for row in ellipse["table"]:
        lines.append(
            f"| {row['heel_deg']:.0f} | {row['gz_library']:.9f} | {row['gz_independent']:.9f} | "
            f"{row['diff']:.3e} |"
        )
    lines += [
        "",
        "| quantity | merged artifact | library today | independent route |",
        "| --- | --- | --- | --- |",
        f"| GM from waterplane moments [m] | {ellipse['gm_published_m']:.9f} | "
        f"(closed-form slope {ellipse['gm_slope_library_m']:.9f}) | "
        f"{ellipse['gm_waterplane_independent_m']:.9f} |",
        f"| GM from dGZ/dphi(0) [m] | - | {ellipse['gm_slope_library_m']:.9f} | "
        f"{ellipse['gm_slope_independent_m']:.9f} |",
        f"| KB [m] | - | {ellipse['kb_library_m']:.9f} | {ellipse['kb_independent_m']:.9f} |",
        f"| GZ_max [m] | {ellipse['gz_max_published_m']:.9f} at "
        f"{ellipse['gz_max_published_heel_deg']:.0f} deg (5 deg grid) | - | "
        f"{ellipse['gz_max_independent_m']:.9f} at "
        f"{ellipse['gz_max_independent_heel_deg']:.0f} deg (same grid); "
        f"{ellipse['gz_max_refined_m']:.9f} at {ellipse['gz_max_refined_heel_deg']:.0f} deg "
        f"on a 1 deg grid |",
        f"| vanishing angle at KG = 1.6 m [deg] | {ellipse['vanishing_angle_published_deg']:.3f} | "
        f"{ellipse['vanishing_angle_library_deg']:.3f} | {ellipse['vanishing_angle_independent_deg']:.3f} |",
        "",
        "## Chord convergence measured against the independent oracle",
        "",
        "| chords | max abs error [m] | max rel error | observed order |",
        "| --- | --- | --- | --- |",
    ]
    for row in ellipse["convergence"]:
        order = f"{row['observed_order']:.2f}" if "observed_order" in row else "-"
        lines.append(
            f"| {row['sides']} | {row['max_abs_error_m']:.3e} | {row['max_rel_error']:.3e} | {order} |"
        )
    lines += [
        "",
        "## Box B = 4 m, T = 1 m, D = 2 m, KG = 0.6 m",
        "",
        f"Wall-sided validity limit: {box['wall_sided_limit_deg']:.3f} deg. Inside it the",
        "independently transcribed wall-sided formula reproduces the oracle; outside it the",
        "independent slice oracle still furnishes the true curve, which the library had no",
        "reference for.",
        "",
        "| heel [deg] | exact [m] | independent [m] | wall-sided [m] | exact - wall-sided [m] |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in box["table"]:
        lines.append(
            f"| {row['heel_deg']:.0f} | {row['gz_library']:.9f} | {row['gz_independent']:.9f} | "
            f"{row['gz_wall_sided']:.9f} | {row['library_minus_wall_sided']:.3e} |"
        )
    lines += [
        "",
        "## Inscribed-polygon defect (fixed in this mission)",
        "",
        "| sides | keel offset [m] | max quadric on boundary | area / ellipse area |",
        "| --- | --- | --- | --- |",
    ]
    for row in polygon["rows"]:
        lines.append(
            f"| {row['sides']} | {row['keel_offset_m']:.3e} | "
            f"{row['outside_ellipse_quadric_max']:.9f} | {row['area_ratio_to_ellipse']:.9f} |"
        )
    lines += [
        "",
        "A quadric value above 1 means the polygon boundary left the ellipse. Before the fix",
        "`sides = 17` reached 1.00851 and `sides = 18` reached 1.03015, i.e. up to 3.0 % outside:",
        "the sampling had no vertex at the keel and the `- min(z)` shift translated the polygon",
        "through the exact section, so the stated 'inscribed, known-sign chord error' did not hold.",
        "",
        "## Rejection contract",
        "",
        "| probe | rejected |",
        "| --- | --- |",
    ]
    for name, row in metrics["rejections"]["probes"].items():
        lines.append(f"| {name} | {'yes' if row['rejected'] else 'NO'} |")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", type=Path, default=Path(__file__).resolve().parent / "outputs"
    )
    arguments = parser.parse_args(argv)
    metrics = run(arguments.out)
    for row in metrics["checks"]:
        status = "PASS" if row["passed"] else "FAIL"
        print(f"{status}  {row['check']}: {row['value']}")
    print(f"artifacts written to {arguments.out}")
    return 0 if metrics["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())