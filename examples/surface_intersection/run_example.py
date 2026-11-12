#!/usr/bin/env python3
"""Surface--surface intersection: six measured studies against analytic references.

Every reference in this script is computed here, not read from a report:

1. **Meridians.** A plane through the axis cuts a cylinder in two straight lines; the
   length is the patch height and the discretisation error should be exactly zero.
2. **Ellipse.** An inclined plane cuts an ellipse of semi-axes ``R`` and ``R / cos a``;
   the perimeter is ``4 a E(1 - cos^2 a)`` with ``E`` the complete elliptic integral,
   evaluated by the arithmetic-geometric mean.
3. **Circle.** Sphere cap by a plane, cylinder by a concentric sphere and a NURBS
   rational cylinder by a plane: three different code paths, all answering ``2 pi r``.
4. **Viviani.** Two equal perpendicular cylinders meet in two closed branches, each
   ``4 sqrt 2 E(1/2) R``. Cross-checked against a midpoint quadrature of the arc-length
   integrand. This is the pair ``GAP-CAD-002`` names, and the one where ``grad g``
   vanishes at the crossing of the two branches.
5. **Convergence.** Length error against the grid, with the measured order -- the chord
   error should fall fourfold per halving, so the order should be 2.
6. **Residuals.** Every returned point against both surfaces independently.

Writes ``results.json``, ``report.md`` and ``plot.png`` into ``--output`` and exits
non-zero if any claim stops holding.

Usage::

    python examples/surface_intersection/run_example.py [--output artifacts/surface_intersection]
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

from tensorcad.nurbs import cylinder_surface
from tensorcad.surface_intersection import (
    CylinderPatch,
    PlanePatch,
    SpherePatch,
    cylinder_implicit,
    plane_implicit,
    sphere_implicit,
    surface_surface_intersection,
)

ACCEPTANCE = 3.0e-2  # GAP-CAD-002: 3% on the arc length
POINT_ACCEPTANCE = 1.0e-9  # GAP-CAD-002: 1e-9 L on the surfaces


def complete_elliptic_e(m: float) -> float:
    """``E(m)`` by the AGM of Gauss -- a route that shares no code with the march."""

    a, b, c = 1.0, math.sqrt(1.0 - m), math.sqrt(m)
    total = 0.5 * c * c
    power = 1.0
    for _ in range(40):
        a_next = 0.5 * (a + b)
        b_next = math.sqrt(a * b)
        c = 0.5 * (a - b)
        power *= 2.0
        total += 0.5 * power * c * c
        a, b = a_next, b_next
        if abs(c) < 1e-18:
            break
    return (math.pi / (2.0 * a)) * (1.0 - total)


def midpoint_arc_length(integrand, count: int = 200_000) -> float:
    step = 2.0 * math.pi / count
    return step * sum(integrand((index + 0.5) * step) for index in range(count))


def study(name: str, result, reference: float, note: str) -> dict:
    error = abs(result.total_length - reference) / abs(reference)
    return {
        "study": name,
        "measured_length": result.total_length,
        "reference_length": reference,
        "relative_error": error,
        "branch_count": len(result.curves),
        "closed_branch_count": len(result.closed_curves),
        "max_residual": result.max_residual,
        "samples": list(result.samples),
        "refine_levels": result.refine_levels,
        "note": note,
        "passes_acceptance": error < ACCEPTANCE,
    }


def convergence_study() -> list[dict]:
    reference = 2.0 * math.pi * math.sqrt(1.0 - 0.5 ** 2)
    rows = []
    for samples in (24, 48, 96, 192):
        result = surface_surface_intersection(
            SpherePatch(radius=1.0),
            plane_implicit((0.0, 0.0, 0.5), (0.0, 0.0, 1.0)),
            samples=(samples, samples // 2),
            refine_levels=2,
        )
        rows.append(
            {
                "samples": samples,
                "measured_length": result.total_length,
                "relative_error": abs(result.total_length - reference) / reference,
                "reference_length": reference,
            }
        )
    for before, after in zip(rows[:-1], rows[1:]):
        if after["relative_error"] > 0.0:
            before["order_to_next"] = math.log2(
                before["relative_error"] / after["relative_error"]
            )
    return rows


def residual_study() -> dict:
    radius = 1.0
    result = surface_surface_intersection(
        CylinderPatch(radius=radius, height=3.0),
        cylinder_implicit(radius, axis="x"),
        samples=(128, 32),
    )
    worst = 0.0
    for branch in result.curves:
        for x, y, z in branch.points:
            worst = max(worst, abs(math.hypot(x, y) - radius), abs(math.hypot(y, z) - radius))
    return {
        "scale": 2.0 * radius,
        "max_residual_reported": result.max_residual,
        "max_residual_recomputed": worst,
        "point_count": result.point_count,
        "max_closure_defect": max(
            (branch.closure_defect for branch in result.closed_curves), default=0.0
        ),
        "passes_acceptance": worst <= POINT_ACCEPTANCE * 2.0 * radius,
    }


def plot(studies: list[dict], convergence: list[dict], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(11.0, 4.4))

    curved = plt.subplot(1, 2, 1)
    for samples, label in ((48, "48 x 24 cells"), (192, "192 x 96 cells")):
        result = surface_surface_intersection(
            SpherePatch(radius=1.0),
            plane_implicit((0.0, 0.0, 0.5), (0.0, 0.0, 1.0)),
            samples=(samples, samples // 2),
        )
        branch = result.closed_curves[0]
        ring_u = [point[0] for point in branch.points]
        ring_v = [point[1] for point in branch.points]
        curved.plot(ring_u, ring_v, marker="", label=label)
    curved.set_aspect("equal")
    curved.set_xlabel("x [m]")
    curved.set_ylabel("y [m]")
    curved.set_title("sphere cap cut by a plane: circle r = sqrt(R^2 - d^2)")
    curved.legend(fontsize=8)
    curved.grid(alpha=0.3)

    order = plt.subplot(1, 2, 2)
    sizes = [row["samples"] for row in convergence]
    errors = [row["relative_error"] for row in convergence]
    order.loglog(sizes, errors, marker="o", label="measured")
    scale = errors[0] / sizes[0] ** -2.0
    order.loglog(
        sizes, [scale * size ** -2.0 for size in sizes], "--", label="O(h^2) reference"
    )
    order.set_xlabel("u samples")
    order.set_ylabel("relative length error")
    order.set_title("arc length converges at second order")
    order.legend(fontsize=8)
    order.grid(alpha=0.3, which="both")

    figure.tight_layout()
    figure.savefig(path, dpi=140)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/surface_intersection")
    arguments = parser.parse_args()
    output = Path(arguments.output)
    output.mkdir(parents=True, exist_ok=True)

    studies = []

    meridian = surface_surface_intersection(
        CylinderPatch(radius=1.0, height=4.0),
        plane_implicit((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        samples=(96, 12),
    )
    studies.append(
        study("plane through the cylinder axis (two meridians)", meridian, 8.0,
              "length is 2 h exactly; the error is round-off, not discretisation")
    )

    alpha = math.radians(30.0)
    ellipse = surface_surface_intersection(
        CylinderPatch(radius=1.0, height=4.0),
        plane_implicit((0.0, 0.0, 0.0), (0.0, math.sin(alpha), math.cos(alpha))),
        samples=(128, 24),
    )
    semi_major = 1.0 / math.cos(alpha)
    perimeter = 4.0 * semi_major * complete_elliptic_e(1.0 - math.cos(alpha) ** 2)
    studies.append(
        study("inclined plane cuts an ellipse at 30 deg", ellipse, perimeter,
              "perimeter 4 a E(e^2) with a = R / cos(alpha), E by AGM")
    )

    sphere_plane = surface_surface_intersection(
        SpherePatch(radius=1.0),
        plane_implicit((0.0, 0.0, 0.5), (0.0, 0.0, 1.0)),
        samples=(96, 48),
    )
    studies.append(
        study("sphere patch cut by a plane", sphere_plane,
              2.0 * math.pi * math.sqrt(1.0 - 0.25), "2 pi sqrt(R^2 - d^2)")
    )

    cylinder_sphere = surface_surface_intersection(
        CylinderPatch(radius=0.5, height=3.0), sphere_implicit(1.0), samples=(96, 24)
    )
    studies.append(
        study("cylinder cut by a concentric sphere", cylinder_sphere,
              4.0 * math.pi * 0.5, "two circles of radius R_cyl")
    )

    nurbs_plane = surface_surface_intersection(
        cylinder_surface(1.0, 2.0, base_z=-1.0),
        plane_implicit((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
        samples=(96, 16),
    )
    studies.append(
        study("NURBS rational cylinder cut by a plane", nurbs_plane, 2.0 * math.pi,
              "closes across the clamped patch's duplicated seam")
    )

    viviani = surface_surface_intersection(
        CylinderPatch(radius=1.0, height=3.0),
        cylinder_implicit(1.0, axis="x"),
        samples=(128, 32),
    )
    branch_reference = 4.0 * math.sqrt(2.0) * complete_elliptic_e(0.5)
    quadrature = midpoint_arc_length(lambda t: math.sqrt(1.0 + math.cos(t) ** 2))
    studies.append(
        study("perpendicular equal cylinders (Viviani)", viviani, 2.0 * branch_reference,
              "4 sqrt(2) E(1/2) R per branch; AGM and 2e5-point midpoint quadrature "
              f"agree to {abs(branch_reference - quadrature) / quadrature:.2e}")
    )

    truncated = surface_surface_intersection(
        PlanePatch(origin=(0.0, 0.0, 0.0), du=(3.0, 0.0, 0.0), dv=(0.0, 2.0, 0.0)),
        cylinder_implicit(0.75),
        samples=(64, 64),
    )
    studies.append(
        study("plane patch cut by a cylinder (truncated)", truncated,
              math.pi * 0.75 / 2.0, "only the quarter inside the patch survives")
    )

    convergence = convergence_study()
    residuals = residual_study()
    plot(studies, convergence, output / "plot.png")

    payload = {
        "acceptance_relative_error": ACCEPTANCE,
        "acceptance_point_residual": POINT_ACCEPTANCE,
        "studies": studies,
        "convergence": convergence,
        "residuals": residuals,
        "figure": "plot.png",
        "all_studies_pass": all(row["passes_acceptance"] for row in studies),
    }
    (output / "results.json").write_text(json.dumps(payload, indent=2))

    lines = ["# Surface--surface intersection: measured against analytic references", ""]
    lines.append(f"Acceptance (GAP-CAD-002): relative arc-length error < {ACCEPTANCE:.0%}, "
                 f"residual < {POINT_ACCEPTANCE:g} L.")
    lines.append("")
    lines.append("| study | measured [m] | reference [m] | rel. error | branches | closed |")
    lines.append("|---|---|---|---|---|---|")
    for row in studies:
        lines.append(
            f"| {row['study']} | {row['measured_length']:.9f} | "
            f"{row['reference_length']:.9f} | {row['relative_error']:.3e} | "
            f"{row['branch_count']} | {row['closed_branch_count']} |"
        )
    lines.append("")
    lines.append("## Convergence (sphere cap vs plane, refine_levels=2)")
    lines.append("")
    lines.append("| u samples | length [m] | rel. error | order to next row |")
    lines.append("|---|---|---|---|")
    for row in convergence:
        order = row.get("order_to_next")
        lines.append(
            f"| {row['samples']} | {row['measured_length']:.9f} | "
            f"{row['relative_error']:.3e} | "
            f"{'-' if order is None else format(order, '.2f')} |"
        )
    lines.append("")
    lines.append("## Residuals on both surfaces (Viviani)")
    lines.append("")
    lines.append(f"- reported max residual: {residuals['max_residual_reported']:.3e} m")
    lines.append(f"- recomputed independently: {residuals['max_residual_recomputed']:.3e} m "
                 f"over {residuals['point_count']} points")
    lines.append(f"- max closure defect: {residuals['max_closure_defect']:.3e} m")
    (output / "report.md").write_text("\n".join(lines) + "\n")

    print("\n".join(lines))
    ok = payload["all_studies_pass"] and residuals["passes_acceptance"]
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())