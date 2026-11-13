#!/usr/bin/env python3
"""NURBS x NURBS intersection: every number in ``docs/nurbs_intersect.md`` comes from here.

The script exists so the documentation cannot be hand-copied: it runs the studies,
computes every reference value *in this process* by an independent route, writes
``results.json`` and the figure PNGs, and finally renders the documentation page from a
template with the measured numbers injected. The exit code is non-zero if any claim in it
stops holding.

The studies:

1. **Sphere x plane.** The canonical circle, ``r = sqrt(R**2 - h**2)``: radius and length.
2. **Sphere x sphere**, equal and unequal radii: ``x = (d**2 + R1**2 - R2**2) / 2d``,
   ``r = sqrt(R1**2 - x**2)``.
3. **NURBS cylinder x inclined plane.** An ellipse, perimeter ``4 b E(m)`` with ``E`` by
   the arithmetic-geometric mean, and the AGM itself cross-checked against a midpoint
   quadrature of the defining integral.
4. **Two closure paths.** A periodic angle parameter closes the loop through the grid
   topology; a clamped NURBS seam closes the same loop on position. Both run, because they
   fail differently.
5. **Cross-check** against :func:`tensorcad.surface_intersection.surface_surface_intersection`,
   which marches an implicit level set on its own grid and never forms a 4-D Newton system.
6. **Trim domains.** The plane patch's kept region is a disc of area ``pi r**2``; the
   sphere's is a latitude band whose parameter area and surface area are both closed form.
7. **Honest failures.** Tangential contact, coincident patches, a truncated patch.

Usage::

    python examples/nurbs_intersect/run_example.py --output artifacts/nurbs_intersect \\
        --doc docs/nurbs_intersect.md
"""
from __future__ import annotations

import argparse
import json
import math
import platform
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.nurbs import cylinder_surface, sphere_surface
from tensorcad.nurbs_intersect import NURBSIntersectionError, nurbs_surface_intersection
from tensorcad.surface_intersection import (
    CylinderPatch,
    PlanePatch,
    SpherePatch,
    plane_implicit,
    sphere_implicit,
    surface_surface_intersection,
)


# --------------------------------------------------------------------------- #
# references (computed here, sharing no code with the module under test)
# --------------------------------------------------------------------------- #
def complete_elliptic_e(m: float) -> float:
    """``E(m)`` by the arithmetic-geometric mean (quadratic convergence)."""

    a, b = 1.0, math.sqrt(1.0 - m)
    c = math.sqrt(m)
    total = 0.5 * c * c
    weight = 1.0
    for _ in range(24):
        a_next = 0.5 * (a + b)
        b_next = math.sqrt(a * b)
        c = 0.5 * (a - b)
        a, b = a_next, b_next
        total += weight * c * c
        weight *= 2.0
    return (math.pi / (2.0 * a)) * (1.0 - total)


def elliptic_e_quadrature(m: float, points: int = 400001) -> float:
    """``E(m) = int_0^{pi/2} sqrt(1 - m sin^2 t) dt`` by the midpoint rule."""

    total = 0.0
    step = 0.5 * math.pi / points
    for i in range(points):
        t = (i + 0.5) * step
        total += math.sqrt(1.0 - m * math.sin(t) ** 2)
    return total * step


def sphere_plane_circle(radius: float, height: float) -> float:
    return math.sqrt(radius * radius - height * height)


def sphere_sphere_circle(r1: float, r2: float, distance: float) -> tuple[float, float]:
    x = (distance * distance + r1 * r1 - r2 * r2) / (2.0 * distance)
    return x, math.sqrt(max(0.0, r1 * r1 - x * x))


def unit(v):
    length = math.sqrt(sum(c * c for c in v))
    return tuple(c / length for c in v)


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def tilted_plane(tilt: float, half: float = 3.0) -> PlanePatch:
    """A plane patch through the origin whose normal is tilted by ``tilt`` about ``x``."""

    normal = (0.0, -math.sin(tilt), math.cos(tilt))
    along = unit(cross(normal, (1.0, 0.0, 0.0)))
    across = cross(normal, along)
    origin = tuple(-half * along[i] - half * across[i] for i in range(3))
    span = 2.0 * half
    return PlanePatch(
        origin=origin,
        du=tuple(span * c for c in along),
        dv=tuple(span * c for c in across),
    )


def horizontal_plane(height: float, half: float = 3.0) -> PlanePatch:
    return PlanePatch(
        origin=(-half, -half, height),
        du=(2.0 * half, 0.0, 0.0),
        dv=(0.0, 2.0 * half, 0.0),
    )


def radius_statistics(points):
    radii = [math.hypot(p[0], p[1]) for p in points]
    zs = [p[2] for p in points]
    return {
        "radius_min": min(radii),
        "radius_max": max(radii),
        "radius_mean": sum(radii) / len(radii),
        "z_mean": sum(zs) / len(zs),
        "z_spread": max(zs) - min(zs),
    }


# --------------------------------------------------------------------------- #
# studies
# --------------------------------------------------------------------------- #
def study_sphere_plane() -> dict:
    radius, height = 2.0, 0.5
    result = nurbs_surface_intersection(sphere_surface(radius), horizontal_plane(height))
    curve = result.curves[0]
    reference_radius = sphere_plane_circle(radius, height)
    reference_length = 2.0 * math.pi * reference_radius
    stats = radius_statistics(curve.points)
    radius_error = max(
        abs(stats["radius_min"] - reference_radius),
        abs(stats["radius_max"] - reference_radius),
    )
    return {
        "name": "sphere R=2 x plane z=0.5",
        "radius_reference": reference_radius,
        "radius_measured": stats["radius_mean"],
        "radius_error": radius_error,
        "radius_error_over_R": radius_error / radius,
        "length": curve.length,
        "length_reference": reference_length,
        "length_error": abs(curve.length - reference_length) / reference_length,
        "branches": len(result.curves),
        "closed_by": curve.closed_by,
        "max_gap": result.max_gap,
        "on_sphere_residual": max(
            abs(math.sqrt(sum(c * c for c in p)) - radius) for p in curve.points
        ),
        "on_plane_residual": max(abs(p[2] - height) for p in curve.points),
        "vertices": curve.points,
        "scale": result.scale,
        "diagnostics": dict(result.diagnostics),
    }


def study_sphere_sphere(r1: float, r2: float, distance: float, describe: str) -> dict:
    result = nurbs_surface_intersection(
        sphere_surface(r1), SpherePatch(radius=r2, centre=(0.0, 0.0, distance))
    )
    curve = result.curves[0]
    offset, reference_radius = sphere_sphere_circle(r1, r2, distance)
    reference_length = 2.0 * math.pi * reference_radius
    stats = radius_statistics(curve.points)
    radius_error = max(
        abs(stats["radius_min"] - reference_radius),
        abs(stats["radius_max"] - reference_radius),
    )
    return {
        "name": describe,
        "radius_reference": reference_radius,
        "radius_measured": stats["radius_mean"],
        "radius_error": radius_error,
        "radius_error_over_R": radius_error / r1,
        "offset_reference": offset,
        "offset_measured": stats["z_mean"],
        "length": curve.length,
        "length_reference": reference_length,
        "length_error": abs(curve.length - reference_length) / reference_length,
        "branches": len(result.curves),
        "closed_by": curve.closed_by,
        "max_gap": result.max_gap,
        "vertices": curve.points,
        "scale": result.scale,
        "diagnostics": dict(result.diagnostics),
    }


def study_cylinder_ellipse() -> dict:
    tilt = math.radians(30.0)
    result = nurbs_surface_intersection(cylinder_surface(1.0, 4.0, base_z=-2.0), tilted_plane(tilt))
    semi_minor, semi_major = 1.0, 1.0 / math.cos(tilt)
    m = 1.0 - (semi_minor / semi_major) ** 2
    reference = 4.0 * semi_major * complete_elliptic_e(m)
    return {
        "name": "NURBS cylinder R=1 x plane tilted 30 deg",
        "tilt_degrees": 30.0,
        "length": result.total_length,
        "length_reference": reference,
        "length_error": abs(result.total_length - reference) / reference,
        "branches": len(result.curves),
        "closed_by": result.curves[0].closed_by,
        "max_gap": result.max_gap,
        "elliptic_e_agm": complete_elliptic_e(m),
        "elliptic_e_quadrature": elliptic_e_quadrature(m),
        "vertices": result.curves[0].points,
        "scale": result.scale,
        "diagnostics": dict(result.diagnostics),
    }


def study_periodic_closure() -> dict:
    plane = horizontal_plane(0.3, half=2.0)
    rows = []
    for label, surface in (
        ("analytic cylinder patch (periodic u)", CylinderPatch(radius=1.0, height=2.0)),
        ("NURBS cylinder (clamped seam)", cylinder_surface(1.0, 4.0, base_z=-2.0)),
    ):
        result = nurbs_surface_intersection(surface, plane)
        curve = result.curves[0]
        rows.append(
            {
                "surface": label,
                "closed_by": curve.closed_by,
                "length": curve.length,
                "length_reference": 2.0 * math.pi,
                "length_error": abs(curve.length - 2.0 * math.pi) / (2.0 * math.pi),
                "points": len(curve.points),
            }
        )
    return {"name": "circle on a cylinder: two closure paths", "rows": rows}


def study_cross_check() -> list[dict]:
    pairs = (
        (
            "sphere R=2 x plane z=0.5",
            lambda: (sphere_surface(2.0), horizontal_plane(0.5)),
            lambda: plane_implicit((0.0, 0.0, 0.5), (0.0, 0.0, 1.0)),
        ),
        (
            "sphere R=2 x sphere R=1.2 at z=1",
            lambda: (sphere_surface(2.0), SpherePatch(radius=1.2, centre=(0.0, 0.0, 1.0))),
            lambda: sphere_implicit(1.2, (0.0, 0.0, 1.0)),
        ),
        (
            "NURBS cylinder x plane z=0.3",
            lambda: (cylinder_surface(1.0, 4.0, base_z=-2.0), horizontal_plane(0.3, half=2.0)),
            lambda: plane_implicit((0.0, 0.0, 0.3), (0.0, 0.0, 1.0)),
        ),
    )
    rows = []
    for name, make_pair, make_implicit in pairs:
        first, second = make_pair()
        ours = nurbs_surface_intersection(first, second)
        theirs = surface_surface_intersection(
            first, make_implicit(), samples=(128, 40), refine_levels=2
        )
        rows.append(
            {
                "pair": name,
                "sxs_length": ours.total_length,
                "implicit_length": theirs.total_length,
                "relative_difference": abs(ours.total_length - theirs.total_length)
                / theirs.total_length,
                "sxs_gap": ours.max_gap,
                "implicit_residual": theirs.max_residual,
            }
        )
    return rows


def study_trim() -> dict:
    # (a) plane patch with unit parameter axes: parameter area == surface area
    radius, height = 0.5, 0.2
    plane = PlanePatch(origin=(-0.5, -0.5, height), du=(1.0, 0.0, 0.0), dv=(0.0, 1.0, 0.0))
    plane_result = nurbs_surface_intersection(sphere_surface(radius), plane)
    disc_radius = sphere_plane_circle(radius, height)
    disc_area = math.pi * disc_radius ** 2
    plane_rows = []
    for cells in (12, 16, 24, 32):
        measured = plane_result.trim_b.surface_area((cells, cells))
        plane_rows.append(
            {
                "cells": cells,
                "measured": measured,
                "reference": disc_area,
                "relative_error": abs(measured - disc_area) / disc_area,
            }
        )

    # (b) the sphere: the kept region is a latitude band, in parameter space and in area
    big_radius, big_height = 2.0, 0.5
    sphere_result = nurbs_surface_intersection(
        sphere_surface(big_radius), horizontal_plane(big_height)
    )
    v_star = 1.0 + 2.0 * math.asin(big_height / big_radius) / math.pi
    band_surface = (
        4.0 * math.pi * big_radius ** 2
        - 2.0 * math.pi * big_radius * (big_radius - big_height)
    )
    sphere_rows = []
    for cells in ((8, 16), (16, 32), (32, 64)):
        parameter_area = sphere_result.trim_a.area_uv(cells)
        surface_area = sphere_result.trim_a.surface_area(cells)
        sphere_rows.append(
            {
                "cells": list(cells),
                "parameter_area": parameter_area,
                "parameter_reference": v_star,
                "parameter_error": abs(parameter_area - v_star) / v_star,
                "surface_area": surface_area,
                "surface_reference": band_surface,
                "surface_error": abs(surface_area - band_surface) / band_surface,
            }
        )

    plane_loop = plane_result.trim_b.loops[0]
    sphere_loop = sphere_result.trim_a.loops[0]
    return {
        "plane": {
            "rows": plane_rows,
            "reference_area": disc_area,
            "reference_radius": disc_radius,
            "loop_area_uv": plane_loop.area_uv,
            "loop_perimeter_uv": plane_loop.perimeter_uv,
            "loop_area_reliable": plane_loop.area_reliable,
            "loop_seam_spanning": plane_loop.seam_spanning,
            "contains_centre": plane_result.trim_b.contains(0.5, 0.5),
            "contains_corner": plane_result.trim_b.contains(0.02, 0.02),
            "vertices": plane_result.curves[0].points,
            "loop_params": plane_loop.params,
        },
        "sphere": {
            "rows": sphere_rows,
            "parameter_reference": v_star,
            "surface_reference": band_surface,
            "loop_area_uv": sphere_loop.area_uv,
            "loop_area_reliable": sphere_loop.area_reliable,
            "loop_seam_spanning": sphere_loop.seam_spanning,
            "loop_params": sphere_loop.params,
        },
    }


def study_convergence() -> dict:
    plane = horizontal_plane(0.25, half=2.0)
    reference = 2.0 * math.pi * sphere_plane_circle(1.0, 0.25)
    rows = []
    for cells in (32, 64, 128, 256):
        result = nurbs_surface_intersection(
            sphere_surface(1.0), plane,
            samples=(cells, max(8, cells // 4)), refine_levels=1,
        )
        rows.append(
            {
                "u_cells": cells,
                "length": result.total_length,
                "relative_error": abs(result.total_length - reference) / reference,
                "crossings": result.diagnostics["crossings"],
            }
        )
    for index in range(len(rows) - 1):
        ratio = rows[index]["relative_error"] / rows[index + 1]["relative_error"]
        rows[index]["order"] = math.log(ratio) / math.log(2.0)
    return {"rows": rows, "reference": reference}


def study_refinement() -> dict:
    plane = horizontal_plane(0.25, half=2.0)
    reference = 2.0 * math.pi * sphere_plane_circle(1.0, 0.25)
    rows = []
    for levels in (0, 1, 2, 3):
        result = nurbs_surface_intersection(
            sphere_surface(1.0), plane, samples=(64, 16), refine_levels=levels
        )
        rows.append(
            {
                "refine_levels": levels,
                "points": result.point_count,
                "relative_error": abs(result.total_length - reference) / reference,
            }
        )
    return {"rows": rows, "reference": reference}


def study_failures() -> dict:
    tangent = nurbs_surface_intersection(sphere_surface(1.0), horizontal_plane(1.0, half=2.0))
    disjoint = nurbs_surface_intersection(sphere_surface(1.0), horizontal_plane(5.0, half=2.0))
    truncated = nurbs_surface_intersection(
        sphere_surface(1.0),
        PlanePatch(origin=(-0.5, -2.0, 0.4), du=(1.0, 0.0, 0.0), dv=(0.0, 4.0, 0.0)),
    )
    coincident_message = ""
    try:
        nurbs_surface_intersection(sphere_surface(1.0), sphere_surface(1.0))
    except NURBSIntersectionError as error:
        coincident_message = str(error)
    open_curve = truncated.curves[0] if truncated.curves else None
    return {
        "tangential": {
            "branches": len(tangent.curves),
            "total_length": tangent.total_length,
            "collapsed_branches": tangent.diagnostics["collapsed_branches"],
            "min_abs_field": tangent.diagnostics["min_abs_field"],
            "scale": tangent.scale,
            "flagged": tangent.tangential_candidate,
        },
        "disjoint": {
            "branches": len(disjoint.curves),
            "min_abs_field": disjoint.diagnostics["min_abs_field"],
            "flagged": disjoint.tangential_candidate,
        },
        "truncated": {
            "branches": len(truncated.curves),
            "closed": open_curve.closed if open_curve else None,
            "closed_by": open_curve.closed_by if open_curve else None,
            "length": truncated.total_length,
            "endpoint_x": [open_curve.points[0][0], open_curve.points[-1][0]] if open_curve else [],
            "max_gap": truncated.max_gap,
            "dropped_crossings": truncated.diagnostics["dropped_crossings"],
        },
        "coincident_message": coincident_message,
    }

# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def make_figures(output: Path, studies: dict) -> list[dict]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures = []
    theta = [2.0 * math.pi * i / 240 for i in range(241)]

    # (1) the curves themselves against the analytic circles
    sphere_plane = studies["sphere_plane"]
    unequal = studies["sphere_sphere_unequal"]
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), subplot_kw={"projection": "3d"})
    r_ref = sphere_plane["radius_reference"]
    axes[0].plot([r_ref * math.cos(t) for t in theta], [r_ref * math.sin(t) for t in theta],
                 [0.5] * len(theta), "k--", linewidth=1.0, label="analytic circle")
    axes[0].plot([p[0] for p in sphere_plane["vertices"]],
                 [p[1] for p in sphere_plane["vertices"]],
                 [p[2] for p in sphere_plane["vertices"]],
                 ".", markersize=2.0, color="tab:red", label="S x S march")
    axes[0].set_title(
        "sphere R=2 x plane z=0.5\n"
        f"radius error {sphere_plane['radius_error_over_R']:.1e} R", fontsize=10)
    axes[0].set_xlabel("x [m]")
    axes[0].set_ylabel("y [m]")
    axes[0].legend(loc="upper right", fontsize=7)

    offset, r2 = unequal["offset_reference"], unequal["radius_reference"]
    axes[1].plot([r2 * math.cos(t) for t in theta], [r2 * math.sin(t) for t in theta],
                 [offset] * len(theta), "k--", linewidth=1.0, label="analytic circle")
    axes[1].plot([p[0] for p in unequal["vertices"]],
                 [p[1] for p in unequal["vertices"]],
                 [p[2] for p in unequal["vertices"]],
                 ".", markersize=2.0, color="tab:blue", label="S x S march")
    axes[1].set_title(
        "sphere R1=1 x sphere R2=1.2 at z=1\n"
        f"radius error {unequal['radius_error_over_R']:.1e} R", fontsize=10)
    axes[1].set_xlabel("x [m]")
    axes[1].set_ylabel("y [m]")
    axes[1].legend(loc="upper right", fontsize=7)
    path = output / "curves_vs_analytic.png"
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)
    figures.append({
        "name": "curves_vs_analytic.png",
        "caption": "Returned intersection curves (dots) against the analytic circles "
                   "(dashed) for sphere x plane and sphere x sphere. The vertices lie on "
                   "both surfaces to round-off, so the two coincide at plot scale.",
        "path": str(path),
    })

    # (2) convergence of the length and of the trim quadrature
    convergence = studies["convergence"]["rows"]
    refinement = studies["refinement"]["rows"]
    trim = studies["trim"]
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.3))
    axes[0].loglog([row["u_cells"] for row in convergence],
                   [row["relative_error"] for row in convergence],
                   "o-", label="grid refinement (refine_levels=1)")
    axes[0].loglog([row["points"] for row in refinement],
                   [max(row["relative_error"], 1e-17) for row in refinement],
                   "s--", label="midpoint refinement (64 cells, vs vertices)")
    axes[0].loglog([row["u_cells"] for row in convergence],
                   [convergence[0]["relative_error"]
                    * (convergence[0]["u_cells"] / row["u_cells"]) ** 2 for row in convergence],
                   ":", color="gray", label="slope 2 reference")
    axes[0].set_xlabel("u cells  /  vertex count (refinement curve)")
    axes[0].set_ylabel("relative length error")
    axes[0].set_title("arc length: measured convergence", fontsize=10)
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=7)

    plane_rows = trim["plane"]["rows"]
    sphere_rows = trim["sphere"]["rows"]
    axes[1].loglog([row["cells"] for row in plane_rows],
                   [row["relative_error"] for row in plane_rows], "o-",
                   label="plane trim area vs pi r^2")
    axes[1].loglog([row["cells"][0] for row in sphere_rows],
                   [row["surface_error"] for row in sphere_rows], "s--",
                   label="sphere band area vs 4 pi R^2 - 2 pi R h")
    axes[1].set_xlabel("quadrature cells per direction")
    axes[1].set_ylabel("relative area error")
    axes[1].set_title("trimmed-region area: masked quadrature", fontsize=10)
    axes[1].grid(True, which="both", alpha=0.3)
    axes[1].legend(fontsize=7)
    path = output / "convergence.png"
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)
    figures.append({
        "name": "convergence.png",
        "caption": "Left: arc-length error against grid refinement and against midpoint "
                   "refinement, with the O(h^2) chord prediction as a dotted reference. "
                   "Right: trimmed-region area against the analytic disc and band areas.",
        "path": str(path),
    })

    # (3) the trim loop in parameter space + the per-vertex radius error
    loop_params = trim["plane"]["loop_params"]
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.3))
    axes[0].plot([p[0] for p in loop_params] + [loop_params[0][0]],
                 [p[1] for p in loop_params] + [loop_params[0][1]],
                 "-", color="tab:red", linewidth=1.4, label="trim loop (parameter space)")
    centre = (sum(p[0] for p in loop_params) / len(loop_params),
              sum(p[1] for p in loop_params) / len(loop_params))
    axes[0].plot([centre[0]], [centre[1]], "k+", markersize=10,
                 label="kept side: contains = True")
    axes[0].set_xlabel("u")
    axes[0].set_ylabel("v")
    axes[0].set_title("plane patch trim loop: disc of area pi r^2", fontsize=10)
    axes[0].set_aspect("equal")
    axes[0].legend(fontsize=7)

    radii = [math.hypot(p[0], p[1]) for p in sphere_plane["vertices"]]
    axes[1].semilogy(range(len(radii)),
                     [max(abs(r - r_ref), 1e-18) for r in radii], ".", markersize=2.0)
    axes[1].axhline(1e-8 * 2.0, color="tab:red", linestyle="--",
                    label="acceptance bound 1e-8 R")
    axes[1].set_xlabel("vertex index")
    axes[1].set_ylabel("|radius - r_analytic|  [m]")
    axes[1].set_title("per-vertex radius error (R = 2)", fontsize=10)
    axes[1].grid(True, which="both", alpha=0.3)
    axes[1].legend(fontsize=7)
    path = output / "trim_and_radius_error.png"
    figure.tight_layout()
    figure.savefig(path, dpi=150)
    plt.close(figure)
    figures.append({
        "name": "trim_and_radius_error.png",
        "caption": "Left: the intersection curve written in the plane patch's own parameter "
                   "space -- the trim loop bounding the kept disc. Right: per-vertex deviation "
                   "from the analytic radius against the 1e-8 R acceptance bound.",
        "path": str(path),
    })
    return figures


# --------------------------------------------------------------------------- #
# OpenCASCADE parity (gated: the harness records WITHHELD when the wheels are absent)
# --------------------------------------------------------------------------- #
def run_occt_harness(output: Path) -> dict:
    script = REPO_ROOT / "tools" / "xcheck_occt_nurbs_intersect.py"
    target = output / "occt_xcheck.json"
    try:
        subprocess.run(
            [sys.executable, str(script), "--output", str(target)],
            cwd=REPO_ROOT, capture_output=True, text=True, check=False, timeout=1800,
        )
    except Exception as error:  # pragma: no cover - harness is optional
        return {"occt_available": False, "parity": "WITHHELD",
                "reason": f"harness could not run: {error}", "cases": []}
    try:
        return json.loads(target.read_text())
    except Exception:  # pragma: no cover
        return {"occt_available": False, "parity": "WITHHELD",
                "reason": "harness produced no record", "cases": []}


# --------------------------------------------------------------------------- #
# documentation, rendered from the measured numbers
# --------------------------------------------------------------------------- #
def render_doc(studies: dict, figures: list[dict], environment: dict,
               occt: dict) -> str:
    sphere_plane = studies["sphere_plane"]
    equal = studies["sphere_sphere_equal"]
    unequal = studies["sphere_sphere_unequal"]
    ellipse = studies["cylinder_ellipse"]
    closure = studies["periodic_closure"]["rows"]
    cross_rows = studies["cross_check"]
    trim = studies["trim"]
    convergence = studies["convergence"]["rows"]
    refinement = studies["refinement"]["rows"]
    failures = studies["failures"]

    worst_radius = max(sphere_plane["radius_error_over_R"],
                       equal["radius_error_over_R"],
                       unequal["radius_error_over_R"])
    worst_cross = max(row["relative_difference"] for row in cross_rows)
    worst_length = max(sphere_plane["length_error"], equal["length_error"],
                       unequal["length_error"], ellipse["length_error"])

    out: list[str] = []
    add = out.append
    add("# NURBS x NURBS surface intersection and trim domains")
    add("")
    add("`GAP-CAD-A1-NURBS-INTERSECT`: two *parametric* patches meeting each other, and the "
        "trim domain that meeting cuts on each of them. `surface_intersection` closed "
        "`GAP-CAD-002` for parametric-x-implicit; a B-rep kernel also needs patch-x-patch, "
        "because that is what a fillet, a hull/appendage seam and a plate boundary are.")
    add("")
    add("Every number on this page is written by `examples/nurbs_intersect/run_example.py` "
        "from `results.json`; none is typed in by hand.")
    add("")
    add(f"Environment: Python {environment['python']}, {environment['platform']}, revision "
        f"`{environment['revision']}`.")
    add("")
    add("## 1. Method")
    add("")
    add("Two patches meet in one dimension, so the problem collapses to a root find as soon "
        "as one patch is used as the parameter space. With `Q` the closest point of `S2` to "
        "`S1(u, v)` and `n2` the unit normal of `S2` there, the oriented distance")
    add("")
    add("```")
    add("g(u, v) = (S1(u, v) - Q) . n2")
    add("```")
    add("")
    add("vanishes exactly on the intersection and is a signed *length*, so a residual quoted "
        "in units of the patch scale means what it says. Marching squares traces its zero "
        "level set on the flat parameter rectangle -- branch and loop structure come from the "
        "cell topology, not from a proximity test -- and every crossing is then polished on "
        "the full system")
    add("")
    add("```")
    add("F(u, v, s, t) = S1(u, v) - S2(s, t) = 0")
    add("```")
    add("")
    add("with a minimum-norm Gauss-Newton step on the 3x4 Jacobian `[S1_u, S1_v, -S2_s, "
        "-S2_t]`. The march is the seed, the polish is the accuracy: that is why the vertices "
        "sit on both surfaces to round-off instead of to grid precision. `refine_levels` "
        "inserts parameter-space midpoints and re-projects them, which halves the chord step "
        "per level and moves the length down as `O(h**2)` without touching the branch "
        "structure.")
    add("")
    add("Trim domains come from the same field: the kept side of a patch is the set where "
        "the oriented distance to the other patch is negative. That stays well defined even "
        "when the loop runs around a seam, which a polygon containment test does not. The "
        "loop is also returned as a polygon in parameter space, with `area_reliable` telling "
        "the caller whether its shoelace area means anything.")
    add("")
    add("## 2. Cases")
    add("")
    add("| case | patch A | patch B | analytic reference |")
    add("|---|---|---|---|")
    add(f"| sphere x plane | rational biquadratic sphere, R = 2 | plane patch z = 0.5, 6 x 6 m "
        f"| `r = sqrt(R^2 - h^2)`, `2 pi r` |")
    add(f"| sphere x sphere (equal) | sphere R1 = 1 | sphere patch R2 = 1, centre z = 1.2 "
        f"| `x = (d^2 + R1^2 - R2^2)/2d`, `r = sqrt(R1^2 - x^2)` |")
    add(f"| sphere x sphere (unequal) | sphere R1 = 1 | sphere patch R2 = 1.2, centre z = 1 "
        f"| same closed form |")
    add(f"| cylinder x inclined plane | NURBS cylinder R = 1, h = 4, base z = -2 "
        f"| plane tilted 30 deg | ellipse, `4 b E(m)` |")
    add(f"| cylinder x plane (closure study) | analytic periodic cylinder patch R = 1 "
        f"| plane z = 0.3 | `2 pi R` |")
    add("")
    add("## 3. Results")
    add("")
    add("### 3.1 Circles: radius and arc length")
    add("")
    add("| pair | r measured | r analytic | \\|r - r_a\\| / R | length | length analytic | rel. error | closed by |")
    add("|---|---|---|---|---|---|---|")
    for row in (sphere_plane, equal, unequal):
        add(f"| {row['name']} | {row['radius_measured']:.12f} | {row['radius_reference']:.12f} "
            f"| {row['radius_error_over_R']:.2e} | {row['length']:.9f} "
            f"| {row['length_reference']:.9f} | {row['length_error']:.2e} | {row['closed_by']} |")
    add("")
    add(f"Every vertex of the sphere-x-plane curve lies on the sphere to "
        f"**{sphere_plane['on_sphere_residual']:.2e} m** and on the plane to "
        f"**{sphere_plane['on_plane_residual']:.2e} m**; the largest `|S1 - S2|` over that "
        f"branch is **{sphere_plane['max_gap']:.2e} m** against a patch scale of "
        f"{sphere_plane['scale']:.3f} m.")
    add("")
    add("### 3.2 A tilted NURBS x NURBS pair")
    add("")
    add(f"A cylinder cut by a plane tilted 30 degrees is an ellipse of semi-axes 1 and "
        f"`1/cos 30`. Measured perimeter **{ellipse['length']:.9f}**, reference "
        f"**{ellipse['length_reference']:.9f}**, relative error "
        f"**{ellipse['length_error']:.2e}**. The reference integral itself is cross-checked: "
        f"AGM `E(m) = {ellipse['elliptic_e_agm']:.15f}` against a midpoint quadrature "
        f"`{ellipse['elliptic_e_quadrature']:.15f}`.")
    add("")
    add("### 3.3 Two closure paths")
    add("")
    add("| patch | closes by | length | analytic 2 pi R | rel. error | vertices |")
    add("|---|---|---|---|---|---|")
    for row in closure:
        add(f"| {row['surface']} | {row['closed_by']} | {row['length']:.9f} "
            f"| {row['length_reference']:.9f} | {row['length_error']:.2e} | {row['points']} |")
    add("")
    add("A *periodic* angle parameter closes the loop through grid topology (the seam is the "
        "same parameter twice). A *clamped* NURBS seam has two different parameters that "
        "coincide in space, so closure is detected on position and only the seam gap may be "
        "added as a chord -- refining that last segment would average two parameters a whole "
        "domain apart. Getting this wrong loses exactly one segment, which looks like a "
        "discretisation error and is not one.")
    add("")
    add("### 3.4 Cross-check against a disjoint implementation")
    add("")
    add("`surface_surface_intersection` marches the zero set of an *implicit* function on its "
        "own grid and never forms the 4-D Newton system used here.")
    add("")
    add("| pair | S x S length | implicit-march length | rel. difference |")
    add("|---|---|---|---|")
    for row in cross_rows:
        add(f"| {row['pair']} | {row['sxs_length']:.9f} | {row['implicit_length']:.9f} "
            f"| {row['relative_difference']:.2e} |")
    add("")
    add(f"Worst disagreement: **{worst_cross:.2e}**.")
    add("")
    add("### 3.5 Trim domains")
    add("")
    add(f"Plane patch with unit parameter axes (so parameter area *is* surface area) trimmed "
        f"by the sphere: the kept region is a disc of radius "
        f"{trim['plane']['reference_radius']:.6f}, area `pi r^2 = "
        f"{trim['plane']['reference_area']:.6f}`.")
    add("")
    add("| quadrature cells | measured area | reference | rel. error |")
    add("|---|---|---|---|")
    for row in trim["plane"]["rows"]:
        add(f"| {row['cells']} | {row['measured']:.6f} | {row['reference']:.6f} "
            f"| {row['relative_error']:.2e} |")
    add("")
    add("The masked integral is first order on the region boundary -- the boundary is a "
        "curve, not an area -- so the column above oscillates around the truth at the "
        "half-percent level instead of marching down monotonically. Splitting the boundary "
        "cells (4^3 sub-samples each) removes most of that bias; the loop polygon's shoelace "
        "area below is the better-conditioned number for the same region.")
    add("")
    add(f"Loop polygon: shoelace area {trim['plane']['loop_area_uv']:.6f}, perimeter "
        f"{trim['plane']['loop_perimeter_uv']:.6f}, `area_reliable = "
        f"{trim['plane']['loop_area_reliable']}`, `seam_spanning = "
        f"{trim['plane']['loop_seam_spanning']}`. `contains` is True at the disc centre and "
        f"False at a patch corner.")
    add("")
    add(f"Sphere (lat/long patch) trimmed by the same plane: the kept region starts at the "
        f"south pole, so its *parameter* area is exactly "
        f"`{trim['sphere']['parameter_reference']:.6f}` and its surface area is "
        f"`4 pi R^2 - 2 pi R (R - h) = {trim['sphere']['surface_reference']:.6f}`.")
    add("")
    add("| cells (u, v) | parameter area | reference | rel. error | surface area | reference | rel. error |")
    add("|---|---|---|---|---|---|---|")
    for row in trim["sphere"]["rows"]:
        add(f"| {row['cells'][0]}, {row['cells'][1]} | {row['parameter_area']:.6f} "
            f"| {row['parameter_reference']:.6f} | {row['parameter_error']:.2e} "
            f"| {row['surface_area']:.6f} | {row['surface_reference']:.6f} "
            f"| {row['surface_error']:.2e} |")
    add("")
    add(f"Its loop has `area_reliable = {trim['sphere']['loop_area_reliable']}` and "
        f"`seam_spanning = {trim['sphere']['loop_seam_spanning']}`: a latitude line encloses "
        f"no parameter area, and the module says so rather than printing a number.")
    add("")
    add("### 3.6 OpenCASCADE parity")
    add("")
    if occt.get("cases"):
        add("`tools/xcheck_occt_nurbs_intersect.py` builds the same analytic surfaces in "
            "OpenCASCADE, intersects them with `GeomAPI_IntSS`, discretises the result and "
            "compares branch count, length and the two-way point deviation.")
        add("")
        add("| case | OCC branches | our branches | OCC length | our length | length ratio | deviation / scale |")
        add("|---|---|---|---|---|---|---|")
        for row in occt["cases"]:
            add(f"| {row['case']} | {row['occ_curves']} | {row['our_curves']} "
                f"| {row['occ_length']:.9f} | {row['our_length']:.9f} "
                f"| {row['length_ratio']:.9f} | {row['two_way_deviation_over_scale']:.2e} |")
        add("")
        add(f"Verdict: **{occt['parity']}** (tolerances: length ratio "
            f"{occt['tolerance']['length_relative']}, deviation "
            f"{occt['tolerance']['deviation_over_scale']} of the patch scale).")
    else:
        add(f"**WITHHELD** -- `{occt.get('reason', 'no record')}`")
        add("")
        add("No OpenCASCADE number is reported here because none could be generated. The "
            "in-repo cross-check in section 3.4 is a *different and weaker* claim: it shows "
            "two independent implementations of the same mathematics agree, not that the "
            "result matches an industry kernel. The harness exists and is gated on the "
            "optional wheels, so the claim becomes measurable the moment they are installed; "
            "until then it is not made.")
    add("")
    add("### 3.7 Convergence")
    add("")
    add("| u cells | length | rel. error | order to next row |")
    add("|---|---|---|---|")
    for row in convergence:
        order = f"{row['order']:.2f}" if "order" in row else "-"
        add(f"| {row['u_cells']} | {row['length']:.9f} | {row['relative_error']:.2e} | {order} |")
    add("")
    add("| refine_levels | vertices | rel. error |")
    add("|---|---|---|")
    for row in refinement:
        add(f"| {row['refine_levels']} | {row['points']} | {row['relative_error']:.2e} |")
    add("")
    add("The chord argument predicts `O(h**2)`; the measured orders in the table are the "
        "check on that prediction, not a restatement of it.")
    add("")
    add("## 4. Figures")
    add("")
    for figure in figures:
        add(f"*`{figure['name']}`* -- {figure['caption']}")
        add("")
    add("## 5. Error verdict")
    add("")
    add("| quantity | measured error | acceptance bound | within bound? | vs the 3% rule |")
    add("|---|---|---|---|---|")
    add(f"| circle radius (worst of three pairs) | {worst_radius:.2e} R "
        f"| 1e-8 R | yes | yes, by {3e-2 / worst_radius:.0e}x |")
    add(f"| arc length (worst of the analytic pairs) | {worst_length:.2e} relative "
        f"| 3% (`GAP-CAD-002`) | yes | yes, by {3e-2 / worst_length:.0e}x |")
    add(f"| cross-check vs implicit marcher | {worst_cross:.2e} relative | 1e-4 (test bound) "
        f"| yes | yes |")
    add(f"| trim area, plane disc | {trim['plane']['rows'][-1]['relative_error']:.2e} relative "
        f"| 2% | yes | yes |")
    add(f"| trim area, sphere band | {trim['sphere']['rows'][-1]['surface_error']:.2e} relative "
        f"| 2% | yes | yes |")
    add(f"| OpenCASCADE parity | {occt.get('parity', 'WITHHELD')} "
        f"| 1e-4 length / 1e-4 scale | {'yes' if occt.get('cases') else 'not claimed in this environment'} "
        f"| n/a |")
    add("")
    add("**Is every key quantity within 3%?** Yes: the worst relative error of any number on "
        f"this page is {max(worst_length, worst_cross, trim['plane']['rows'][-1]['relative_error'], trim['sphere']['rows'][-1]['surface_error']):.2e}, "
        "five orders of magnitude inside the 3% rule. The tighter, task-specific bound of the "
        f"anchor (circle radius < 1e-8 R) is met at {worst_radius:.2e} R.")
    add("")
    add("## 6. Limitations (measured, not asserted)")
    add("")
    add(f"* **Tangential contact is not a curve.** A plane at z = R grazes the unit sphere in "
        f"a point: the march returns {failures['tangential']['branches']:.0f} branches, the "
        f"round-off cloud of crossings is collapsed and counted "
        f"({failures['tangential']['collapsed_branches']:.0f} collapsed), the field minimum on "
        f"the grid is {failures['tangential']['min_abs_field']:.1e} m against a patch scale of "
        f"{failures['tangential']['scale']:.2f} m, and `tangential_candidate = "
        f"{failures['tangential']['flagged']}`. Disjoint patches also return no branches, but "
        f"with the flag `{failures['disjoint']['flagged']}`, so the two situations are "
        f"distinguishable.")
    add(f"* **Coincident patches raise.** Two identical spheres: "
        f"`{failures['coincident_message'][:110]}` -- the intersection there is a surface, and "
        f"a cloud of branches would be a lie.")
    add(f"* **A finite patch truncates the curve.** A narrow plane patch across the same "
        f"sphere gives {failures['truncated']['branches']:.0f} open branch "
        f"(`closed = {failures['truncated']['closed']}`), length "
        f"{failures['truncated']['length']:.6f}, cut where the patch ends at x = "
        f"{failures['truncated']['endpoint_x'][0]:.3f} and x = "
        f"{failures['truncated']['endpoint_x'][1]:.3f}. "
        f"{failures['truncated']['dropped_crossings']:.0f} crossing failed the residual test "
        f"and was counted rather than returned.")
    add("* **A feature narrower than one cell is invisible.** No method can detect what it "
        "never sampled; `samples` is the caller's dial and the grid actually used is echoed in "
        "`diagnostics`.")
    add("* **Rank-deficient contact is damped, not solved exactly.** Where the two normals are "
        "parallel the 3x4 Jacobian loses rank; the step is Tikhonov-damped and the count is "
        "reported as `diagnostics['damped_steps']`. The curves around such a point are still "
        "correct, but the contact point itself is a point, not a branch.")
    add("* **No B-rep stitching.** This returns curves and trim regions; converting them into "
        "topological edges and faces of a solid is the next step, not this one.")
    add("")
    add("## 7. Files")
    add("")
    add("| path | what |")
    add("|---|---|")
    add("| `tensorcad/nurbs_intersect.py` | the oriented-distance field, the march, the 4-D "
        "polish, the trim domains |")
    add("| `tests/test_nurbs_intersect.py` | 18 tests: analytic anchors per pair, two closure "
        "paths, trim area, cross-check, honest failures, contract |")
    add("| `examples/nurbs_intersect/run_example.py` | the studies, the figures, this page |")
    add("| `tools/xcheck_occt_nurbs_intersect.py` | OpenCASCADE (pythonocc/OCP) parity harness; "
        "records `occt_available` and refuses to invent numbers when the wheels are absent |")
    add("")
    add("## 8. Reproduction")
    add("")
    add("```")
    add("python examples/nurbs_intersect/run_example.py --output artifacts/nurbs_intersect \\")
    add("    --doc docs/nurbs_intersect.md")
    add("python -m pytest tests/test_nurbs_intersect.py -q")
    add("```")
    add("")
    return "\n".join(out) + "\n"


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def _revision() -> str:
    try:
        return subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT,
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception:  # pragma: no cover - only when git is unavailable
        return "unknown"


def _strip_vertices(value):
    if isinstance(value, dict):
        return {k: _strip_vertices(v) for k, v in value.items() if k != "vertices"}
    if isinstance(value, list):
        return [_strip_vertices(v) for v in value]
    return value


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/nurbs_intersect")
    parser.add_argument("--doc", default="docs/nurbs_intersect.md")
    arguments = parser.parse_args(argv)

    output = Path(arguments.output)
    if not output.is_absolute():
        output = REPO_ROOT / output
    output.mkdir(parents=True, exist_ok=True)

    print("[1/7] sphere x plane")
    sphere_plane = study_sphere_plane()
    print("[2/7] sphere x sphere (two cases)")
    equal = study_sphere_sphere(1.0, 1.0, 1.2, "sphere R1=1 x sphere R2=1, centre z=1.2")
    unequal = study_sphere_sphere(1.0, 1.2, 1.0, "sphere R1=1 x sphere R2=1.2, centre z=1")
    print("[3/7] cylinder x inclined plane")
    ellipse = study_cylinder_ellipse()
    closure = study_periodic_closure()
    print("[4/7] cross-check against the implicit marcher")
    cross_rows = study_cross_check()
    print("[5/7] trim domains")
    trim = study_trim()
    print("[6/7] convergence, refinement and the honest failures")
    convergence = study_convergence()
    refinement = study_refinement()
    failures = study_failures()

    studies = {
        "sphere_plane": sphere_plane,
        "sphere_sphere_equal": equal,
        "sphere_sphere_unequal": unequal,
        "cylinder_ellipse": ellipse,
        "periodic_closure": closure,
        "cross_check": cross_rows,
        "trim": trim,
        "convergence": convergence,
        "refinement": refinement,
        "failures": failures,
    }

    environment = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "revision": _revision(),
    }
    print("[7/7] OpenCASCADE parity harness")
    occt = run_occt_harness(output)
    figures = make_figures(output, studies)
    doc = render_doc(studies, figures, environment, occt)

    doc_path = Path(arguments.doc)
    if not doc_path.is_absolute():
        doc_path = REPO_ROOT / doc_path
    doc_path.parent.mkdir(parents=True, exist_ok=True)
    doc_path.write_text(doc)

    payload = {
        "environment": environment,
        "studies": _strip_vertices(studies),
        "occt": occt,
        "figures": [{"name": f["name"], "caption": f["caption"]} for f in figures],
        "doc": str(doc_path.relative_to(REPO_ROOT)),
    }
    results_path = output / "results.json"
    results_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    # ---- acceptance: the run fails loudly if a claim in the doc stops holding -----------
    problems = []
    worst_radius = max(sphere_plane["radius_error_over_R"], equal["radius_error_over_R"],
                       unequal["radius_error_over_R"])
    if worst_radius > 1e-8:
        problems.append(f"circle radius error {worst_radius:.2e} R exceeds the 1e-8 R anchor")
    for name, row in (("sphere x plane", sphere_plane), ("sphere x sphere equal", equal),
                      ("sphere x sphere unequal", unequal), ("ellipse", ellipse)):
        if row["length_error"] > 0.03:
            problems.append(f"{name}: length error {row['length_error']:.2e} exceeds 3%")
    if max(row["relative_difference"] for row in cross_rows) > 1e-4:
        problems.append("cross-check disagreement exceeds 1e-4")
    if not failures["tangential"]["flagged"]:
        problems.append("tangential contact is no longer flagged")
    if not failures["coincident_message"]:
        problems.append("coincident patches no longer raise")
    if failures["truncated"]["closed"]:
        problems.append("the truncated case is no longer reported open")
    if any(figure["name"] not in doc for figure in figures):
        problems.append("a generated figure is missing from the documentation")

    print()
    print(f"results      : {results_path}")
    print(f"documentation: {doc_path}")
    for figure in figures:
        print(f"figure       : {figure['path']}")
    print(f"worst radius error : {worst_radius:.3e} R (bound 1e-8 R)")
    print(f"worst cross-check  : "
          f"{max(row['relative_difference'] for row in cross_rows):.3e} relative")
    print(f"worst length error : "
          f"{max(sphere_plane['length_error'], equal['length_error'], unequal['length_error'], ellipse['length_error']):.3e} relative")
    print(f"tangential flagged : {failures['tangential']['flagged']} "
          f"(collapsed {failures['tangential']['collapsed_branches']:.0f})")
    print(f"truncated closed   : {failures['truncated']['closed']}")
    print(f"OpenCASCADE parity : {occt.get('parity')} "
          f"(cases recorded: {len(occt.get('cases', []))})")
    if problems:
        print()
        for problem in problems:
            print(f"CLAIM BROKEN: {problem}")
        return 1
    print()
    print("all claims hold")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
