"""Shelling a solid: the box closed form, a pentagon erosion, a sphere and the refusals.

Run:  python examples/shelling/run_example.py --output examples/shelling/outputs

Writes ``results.json``, ``report.md`` and ``shelling.png`` and exits non-zero if a claim stops
holding, so CI can run it as a check rather than as a demo. The claims:

1. ``GAP-CAD-011``'s acceptance criterion, run literally: a box shelled over ``t/L`` from 0.01
   to 0.10 matches the analytic ``(L - 2t)^3`` inner volume, far inside the 2 % threshold;
2. a convex pentagon prism -- facets not mutually orthogonal -- matches the independent 2-D
   half-plane clipping oracle in ``docs/shelling.md``;
3. a faceted sphere converges to the analytic spherical-shell volume at second order, so the
   residual really is faceting and not a modelling error;
4. the wall measurement is not a tautology: an under-offset surface is measured as the
   distance it actually is;
5. the failures fail: a plate thinner than twice the thickness and a collapsing taper raise
   with a reason instead of returning a plausible-looking mesh.
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

from tensorcad.section_gz import prism_mesh
from tensorcad.shell import ShellError, hollow_body, min_wall_thickness, offset_inward, shell
import itertools

RECTANGLE = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
PENTAGON = ((0.0, 0.0), (1.0, 0.0), (1.3, 0.6), (0.6, 1.2), (-0.2, 0.7))
SWEEP = [0.01 * step for step in range(1, 11)]
SPHERE = ((12, 6), (24, 12), (48, 24))


def _expect(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"FAIL: {message}")


def _box_sweep() -> dict[str, object]:
    rows = []
    worst = 0.0
    for thickness in SWEEP:
        result = shell(prism_mesh(RECTANGLE, 1.0), thickness)
        expected_inner = (1.0 - 2.0 * thickness) ** 3
        error = abs(result.inner_volume - expected_inner) / expected_inner
        worst = max(worst, error)
        rows.append(
            {
                "t_over_L": thickness,
                "inner_volume": result.inner_volume,
                "inner_volume_analytic": expected_inner,
                "relative_error": error,
                "shell_volume": result.shell_volume,
                "shell_volume_analytic": 1.0 - expected_inner,
                "min_wall_thickness": result.min_wall_thickness,
                "max_plane_residual": result.max_plane_residual,
                "min_normal_dot": result.min_normal_dot,
            }
        )
    return {"rows": rows, "worst_relative_error": worst, "limit": 0.02}


def _clipped_section(section: tuple[tuple[float, float], ...], thickness: float):
    polygon = list(section)
    total = len(section)
    for index in range(total):
        start, end = section[index], section[(index + 1) % total]
        edge_y, edge_z = end[0] - start[0], end[1] - start[1]
        length = math.hypot(edge_y, edge_z)
        normal = (edge_z / length, -edge_y / length)
        bound = normal[0] * start[0] + normal[1] * start[1] - thickness
        result = []
        for position in range(len(polygon)):
            current = polygon[position]
            following = polygon[(position + 1) % len(polygon)]
            dc = normal[0] * current[0] + normal[1] * current[1] - bound
            df = normal[0] * following[0] + normal[1] * following[1] - bound
            if dc <= 0.0:
                result.append(current)
            if (dc < 0.0 < df) or (df < 0.0 < dc):
                ratio = dc / (dc - df)
                result.append(
                    (
                        current[0] + ratio * (following[0] - current[0]),
                        current[1] + ratio * (following[1] - current[1]),
                    )
                )
        polygon = result
    area = 0.0
    for position in range(len(polygon)):
        y0, z0 = polygon[position]
        y1, z1 = polygon[(position + 1) % len(polygon)]
        area += 0.5 * (y0 * z1 - y1 * z0)
    return abs(area), len(polygon)


def _pentagon(thicknesses: tuple[float, ...] = (0.02, 0.05, 0.10, 0.15)) -> dict[str, object]:
    rows = []
    worst = 0.0
    for thickness in thicknesses:
        result = shell(prism_mesh(PENTAGON, 1.0), thickness)
        area, corners = _clipped_section(PENTAGON, thickness)
        reference = area * (1.0 - 2.0 * thickness)
        error = abs(result.inner_volume - reference) / reference
        worst = max(worst, error)
        rows.append(
            {
                "thickness": thickness,
                "inner_volume": result.inner_volume,
                "clipped_reference": reference,
                "relative_error": error,
                "eroded_section_corners": corners,
                "min_wall_thickness": result.min_wall_thickness,
            }
        )
    return {"rows": rows, "worst_relative_error": worst}


def _sphere(radius: float = 1.0, thickness: float = 0.2) -> dict[str, object]:
    analytic = 4.0 / 3.0 * math.pi * (radius**3 - (radius - thickness) ** 3)

    def build(meridians: int, parallels: int) -> object:
        from tensorcad.geometry import MeshBuilder

        builder = MeshBuilder()
        north = builder.add_vertex((0.0, 0.0, radius))
        south = builder.add_vertex((0.0, 0.0, -radius))
        rings = []
        for ring in range(1, parallels):
            phi = math.pi * ring / parallels
            z, rho = radius * math.cos(phi), radius * math.sin(phi)
            rings.append(
                [
                    builder.add_vertex(
                        (
                            rho * math.cos(2.0 * math.pi * column / meridians),
                            rho * math.sin(2.0 * math.pi * column / meridians),
                            z,
                        )
                    )
                    for column in range(meridians)
                ]
            )
        for column in range(meridians):
            builder.add_triangle(north, rings[0][column], rings[0][(column + 1) % meridians])
        for upper, lower in itertools.pairwise(rings):
            for column in range(meridians):
                nxt = (column + 1) % meridians
                builder.add_quad(upper[column], lower[column], lower[nxt], upper[nxt])
        for column in range(meridians):
            builder.add_triangle(south, rings[-1][(column + 1) % meridians], rings[-1][column])
        return builder.build()

    rows = []
    for meridians, parallels in SPHERE:
        mesh = build(meridians, parallels)
        result = shell(mesh, thickness)
        rows.append(
            {
                "meridians": meridians,
                "parallels": parallels,
                "faces": len(mesh.faces),
                "shell_volume": result.shell_volume,
                "shell_volume_analytic": analytic,
                "relative_error": abs(result.shell_volume - analytic) / analytic,
                "outer_volume": result.outer_volume,
                "outer_sphere_error": abs(result.outer_volume - 4.0 / 3.0 * math.pi * radius**3)
                / (4.0 / 3.0 * math.pi * radius**3),
                "min_wall_thickness": result.min_wall_thickness,
            }
        )
    orders = [
        math.log(coarse["relative_error"] / fine["relative_error"]) / math.log(2.0)
        for coarse, fine in itertools.pairwise(rows)
    ]
    return {"rows": rows, "analytic": analytic, "observed_orders": orders, "thickness": thickness}


def _wall_measurement() -> dict[str, object]:
    outer = prism_mesh(RECTANGLE, 1.0)
    target = 0.1
    rows = []
    for factor in (0.2, 0.5, 0.7, 1.0):
        under = offset_inward(outer, factor * target)
        measured = min_wall_thickness(outer, under)
        rows.append(
            {
                "offset_factor": factor,
                "requested_offset": factor * target,
                "measured_wall": measured,
                "relative_error": abs(measured - factor * target) / (factor * target),
            }
        )
    return {"rows": rows}


def _refusals() -> dict[str, object]:
    cases = {
        "plate_0.1_shelled_at_0.08": (
            prism_mesh(((0.0, 0.0), (1.0, 0.0), (1.0, 0.1), (0.0, 0.1)), 1.0),
            0.08,
        ),
        "taper_shelled_at_0.10": (
            prism_mesh(((0.0, 0.0), (1.0, 0.0), (1.0, 0.1), (0.0, 0.5)), 1.0),
            0.10,
        ),
    }
    rows = []
    for name, (mesh, thickness) in cases.items():
        try:
            shell(mesh, thickness)
        except ShellError as error:
            rows.append({"case": name, "refused": True, "reason": str(error).split(";")[0]})
        else:  # pragma: no cover - a silent success is the failure this example exists for
            rows.append({"case": name, "refused": False, "reason": "no error raised"})
    return {"rows": rows}


def _plot(report: dict[str, object], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    box_rows = report["box"]["rows"]
    sphere_rows = report["sphere"]["rows"]
    wall_rows = report["wall"]["rows"]
    figure, axes = plt.subplots(1, 3, figsize=(15.5, 4.4))

    axes[0].semilogy(
        [row["t_over_L"] for row in box_rows],
        [max(row["relative_error"], 1e-18) for row in box_rows],
        "o-",
    )
    axes[0].axhline(0.02, color="crimson", linestyle="--", label="GAP-CAD-011 limit 2 %")
    axes[0].set_xlabel("t / L")
    axes[0].set_ylabel("|inner volume - (L-2t)^3| / (L-2t)^3")
    axes[0].set_title("box: exact erosion over the sweep")
    axes[0].legend(loc="best")

    axes[1].loglog(
        [1.0 / row["meridians"] for row in sphere_rows],
        [row["relative_error"] for row in sphere_rows],
        "s-",
        label="measured",
    )
    guide = sphere_rows[0]["relative_error"]
    axes[1].loglog(
        [1.0 / row["meridians"] for row in sphere_rows],
        [guide * (row["meridians"] / sphere_rows[0]["meridians"]) ** -2 for row in sphere_rows],
        "k--",
        label="slope -2 guide",
    )
    axes[1].set_xlabel("1 / meridians")
    axes[1].set_ylabel("relative shell-volume error")
    axes[1].set_title(f"sphere: faceting error, t = {report['sphere']['thickness']} m")
    axes[1].legend(loc="best")

    axes[2].plot(
        [row["requested_offset"] for row in wall_rows],
        [row["measured_wall"] for row in wall_rows],
        "o-",
        label="measured wall",
    )
    axes[2].plot([0.0, 0.1], [0.0, 0.1], "k--", label="identity")
    axes[2].set_xlabel("offset applied to the inner surface (m)")
    axes[2].set_ylabel("min wall reported by min_wall_thickness (m)")
    axes[2].set_title("the measurement is not a tautology")
    axes[2].legend(loc="best")

    figure.suptitle("TensorCAD shelling: exact box erosion, second-order sphere faceting, falsifiable wall")
    figure.tight_layout()
    figure.savefig(path, dpi=130)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("outputs"))
    arguments = parser.parse_args()
    output = arguments.output
    output.mkdir(parents=True, exist_ok=True)

    report = {
        "box": _box_sweep(),
        "pentagon": _pentagon(),
        "sphere": _sphere(),
        "wall": _wall_measurement(),
        "refusals": _refusals(),
    }
    body = hollow_body(shell(prism_mesh(RECTANGLE, 1.0), 0.1))
    report["hollow_body"] = body.as_dict()

    _expect(report["box"]["worst_relative_error"] < 1e-12, "box erosion is not exact")
    _expect(report["pentagon"]["worst_relative_error"] < 1e-12, "pentagon erosion disagrees")
    _expect(report["sphere"]["rows"][-1]["relative_error"] < 0.01, "sphere shell misses by >1 %")
    _expect(
        all(order > 1.6 for order in report["sphere"]["observed_orders"]),
        "sphere convergence is below second order",
    )
    _expect(all(row["relative_error"] < 1e-12 for row in report["wall"]["rows"]), "wall measurement")
    _expect(all(row["refused"] for row in report["refusals"]["rows"]), "a bad shell was accepted")
    _expect(body.as_dict()["watertight"], "the hollow export is not watertight")

    (output / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    plot_path = output / "shelling.png"
    _plot(report, plot_path)

    lines = [
        "# Shelling: exact erosion, measured wall, loud refusals",
        "",
        "## 1. Box (GAP-CAD-011 acceptance criterion)",
        "",
        "| t/L | inner volume [m^3] | analytic (1-2t)^3 [m^3] | relative error | wall [m] |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in report["box"]["rows"]:
        lines.append(
            f"| {row['t_over_L']:.2f} | {row['inner_volume']:.15f} | "
            f"{row['inner_volume_analytic']:.15f} | {row['relative_error']:.2e} | "
            f"{row['min_wall_thickness']:.12f} |"
        )
    lines += [
        "",
        f"Worst relative error over the sweep: **{report['box']['worst_relative_error']:.2e}** "
        f"(threshold {report['box']['limit']:.0%}).",
        "",
        "## 2. Pentagon prism vs independent 2-D clipping",
        "",
        "| t [m] | inner volume [m^3] | clipped reference [m^3] | relative error | eroded corners |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in report["pentagon"]["rows"]:
        lines.append(
            f"| {row['thickness']:.2f} | {row['inner_volume']:.15f} | "
            f"{row['clipped_reference']:.15f} | {row['relative_error']:.2e} | "
            f"{row['eroded_section_corners']} |"
        )
    lines += [
        "",
        f"Worst relative error: **{report['pentagon']['worst_relative_error']:.2e}**.",
        "",
        f"## 3. Sphere, t = {report['sphere']['thickness']} m",
        "",
        "| meridians x parallels | faces | shell volume [m^3] | analytic [m^3] | rel. error | observed order |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for index, row in enumerate(report["sphere"]["rows"]):
        order = "" if index == 0 else f"{report['sphere']['observed_orders'][index - 1]:.2f}"
        lines.append(
            f"| {row['meridians']} x {row['parallels']} | {row['faces']} | "
            f"{row['shell_volume']:.6f} | {row['shell_volume_analytic']:.6f} | "
            f"{row['relative_error']:.3%} | {order} |"
        )
    lines += [
        "",
        "## 4. The wall measurement, falsified on purpose",
        "",
        "| offset applied [m] | measured wall [m] | relative error |",
        "|---:|---:|---:|",
    ]
    for row in report["wall"]["rows"]:
        lines.append(
            f"| {row['requested_offset']:.3f} | {row['measured_wall']:.12f} | {row['relative_error']:.2e} |"
        )
    lines += [
        "",
        "## 5. Refusals",
        "",
        "| case | refused | reason |",
        "|---|---|---|",
    ]
    for row in report["refusals"]["rows"]:
        lines.append(f"| {row['case']} | {row['refused']} | {row['reason']} |")
    lines += [
        "",
        f"Hollow export: watertight {report['hollow_body']['watertight']}, "
        f"material {report['hollow_body']['material_volume']:.12f} m^3, "
        f"cavity {report['hollow_body']['cavity_volume']:.12f} m^3, "
        f"signed volume of the exported surface {report['hollow_body']['signed_volume']:.12f} m^3.",
        "",
        f"![shelling]({plot_path.name})",
        "",
        "## 6. Error judgement",
        "",
        "| key quantity | reference | error | below 3 %? |",
        "|---|---|---|---|",
        f"| box inner volume, t/L 0.01-0.10 | (L-2t)^3 | {report['box']['worst_relative_error']:.2e} | yes |",
        f"| pentagon prism inner volume | 2-D half-plane clipping | "
        f"{report['pentagon']['worst_relative_error']:.2e} | yes |",
        f"| sphere shell volume, 48 x 24 | 4/3 pi (R^3-(R-t)^3) | "
        f"{report['sphere']['rows'][-1]['relative_error']:.2%} | yes |",
        f"| sphere shell volume, 12 x 6 | 4/3 pi (R^3-(R-t)^3) | "
        f"{report['sphere']['rows'][0]['relative_error']:.2%} | no (faceting, converges at order ~2) |",
        "",
        "No quantisation of the mesh, no smoothing of the inner surface and no tolerance "
        "tuning is applied in any of the numbers above.",
    ]
    (output / "report.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {output / 'results.json'}, {output / 'report.md'} and {plot_path}")
    print(
        f"box worst {report['box']['worst_relative_error']:.2e}; "
        f"pentagon worst {report['pentagon']['worst_relative_error']:.2e}; "
        f"sphere last {report['sphere']['rows'][-1]['relative_error']:.2%}; "
        f"orders {[round(order, 2) for order in report['sphere']['observed_orders']]}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
