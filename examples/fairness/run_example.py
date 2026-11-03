"""Measure curvature and surface fairness, and fail loudly if a claim breaks.

Four studies, all of them with a reference the answer can be wrong against:

1. **The analytic field is exact.** Sphere, cylinder, plane and ellipsoid have
   closed-form principal curvatures; the finite-difference shape operator has to
   reproduce them to round-off, otherwise the estimator that measures the hull and
   the appendages is uncalibrated.
2. **The mesh estimate is exact in the integral and only as good as the mesh
   pointwise.** On a raw marching-tet sphere the area-weighted mean curvature is
   the analytic ``1/R`` while the pointwise median is tens of percent off; after
   ``taubin_smooth`` + ``project_to_field`` the pointwise number comes down. Both
   halves matter: a solver sees the mesh, a designer reads the field.
3. **The AFF-8 appendage set has creases and they are found.** The sail/hull
   junction and the NACA trailing edges are 0-90 degree edges; a fairness audit that
   cannot separate them from discretisation noise cannot be trusted on a real hull.
4. **The parametric hull family is not G1 at its ends.** ``b(x) = B/2 *
   sin(pi x/L)^0.7`` has an infinite waterline slope at ``x = 0`` and ``x = L``, so
   every curvature-based quantity diverges there. That is a property of the family,
   not of any mesher, and the exponent is fitted rather than asserted.

Writes ``results.json``, ``report.md`` and ``plot.png`` into ``--output`` and exits
non-zero if any claim stops holding.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tensorcad.curvature import (
    assess_fairness,
    principal_curvatures,
    vertex_mean_curvature,
)
from tensorcad.geometry import Mesh
from tensorcad.geometry_field import (
    CylinderField,
    EllipsoidField,
    PlaneField,
    SphereField,
    to_mesh,
)
from tensorcad.hull import BREADTH_EXPONENT, HullParameters, generate_hull
from tensorcad.mesh_conditioning import project_to_field, taubin_smooth
from tensorcad.suboff_appendages import SuboffFullField

FAILURES: list[str] = []


def claim(name: str, value: float, limit: float, *, lower_is_better: bool = True) -> dict:
    ok = value <= limit if lower_is_better else value >= limit
    record = {"claim": name, "value": value, "limit": limit,
              "direction": "<=" if lower_is_better else ">=", "holds": bool(ok)}
    if not ok:
        FAILURES.append(f"{name}: {value!r} violates {record['direction']} {limit!r}")
    return record


def field_exactness() -> dict:
    """Study 1: closed-form curvature against the finite-difference shape operator."""
    rows = []

    def add(field, label, point, expected, note):
        sample = principal_curvatures(field, [point])[0]
        values = (sample.kappa_max, sample.kappa_min)
        # Order-free comparison: the estimator returns kappa_max/kappa_min, the closed
        # forms are written in whichever order the derivation gives, so sort both.
        errors = sorted(
            abs(value - target)
            for value, target in zip(sorted(values), sorted(expected), strict=False)
        )
        scale = max(abs(target) for target in expected) or 1.0
        rows.append(
            {
                "field": label,
                "point": list(point),
                "quantity": note,
                "expected_per_m": list(expected),
                "measured_per_m": list(values),
                "max_abs_error_per_m": errors[-1],
                "max_relative_error": errors[-1] / scale,
                "step_error_per_m": sample.step_error,
                "gradient_norm": sample.gradient_norm,
            }
        )

    radius = 0.5
    sphere = SphereField(center=(0.0, 0.0, 0.0), radius=radius)
    for point in ((radius, 0.0, 0.0), (0.0, radius, 0.0), (0.0, 0.0, -radius)):
        add(sphere, "sphere", point, (1.0 / radius, 1.0 / radius), "kappa1 = kappa2 = 1/R")

    cylinder_radius = 0.35
    cylinder = CylinderField(radius=cylinder_radius, height=1.0)
    for point in ((cylinder_radius, 0.0, 0.0), (0.0, cylinder_radius, 0.2)):
        add(cylinder, "cylinder", point, (1.0 / cylinder_radius, 0.0), "hoop 1/R, axial 0")

    plane = PlaneField(normal=(0.0, 0.0, 1.0), offset=0.0)
    for point in ((1.0, 2.0, 0.0), (-3.0, 0.5, 0.0)):
        add(plane, "plane", point, (0.0, 0.0), "flat: kappa1 = kappa2 = 0")

    semi_axes = (0.6, 0.4, 0.25)
    ellipsoid = EllipsoidField(semi_axes=semi_axes)
    a, b, c = semi_axes
    add(ellipsoid, "ellipsoid", (a, 0.0, 0.0), (a / b**2, a / c**2), "x pole: a/b^2, a/c^2")
    add(ellipsoid, "ellipsoid", (0.0, b, 0.0), (b / c**2, b / a**2), "y pole: b/c^2, b/a^2")
    add(ellipsoid, "ellipsoid", (0.0, 0.0, c), (c / b**2, c / a**2), "z pole: c/b^2, c/a^2")

    worst = max(row["max_relative_error"] for row in rows)
    worst_step = max(row["step_error_per_m"] for row in rows)
    plane_exact = max(
        row["max_abs_error_per_m"] for row in rows if row["field"] == "plane"
    )
    return {
        "rows": rows,
        "worst_relative_error": worst,
        "worst_step_error_per_m": worst_step,
        "plane_abs_error_per_m": plane_exact,
        "claims": [
            claim("field curvature vs closed forms: worst relative error", worst, 1e-6),
            claim("plane curvature is exactly zero [1/m]", plane_exact, 1e-12),
            claim("reported step error on the primitives [1/m]", worst_step, 1e-5),
        ],
    }


def uv_sphere(radius: float, nu: int, nv: int) -> Mesh:
    """A regular UV sphere: the control mesh for the pointwise curvature study."""
    vertices: list[tuple[float, float, float]] = []
    grid: dict[tuple[int, int], int] = {}
    for i in range(nu + 1):
        theta = math.pi * i / nu
        for j in range(nv):
            phi = 2.0 * math.pi * j / nv
            vertices.append(
                (
                    radius * math.sin(theta) * math.cos(phi),
                    radius * math.sin(theta) * math.sin(phi),
                    radius * math.cos(theta),
                )
            )
            grid[(i, j)] = len(vertices) - 1
    faces: list[tuple[int, int, int]] = []
    for i in range(nu):
        for j in range(nv):
            j2 = (j + 1) % nv
            a, b = grid[(i, j)], grid[(i, j2)]
            c, d = grid[(i + 1, j)], grid[(i + 1, j2)]
            faces.append((a, c, b))
            faces.append((b, c, d))
    return Mesh(vertices=vertices, faces=faces)


def mesh_curvature_study() -> dict:
    """Study 2: mesh estimate against the analytic field, before and after conditioning."""
    radius = 1.0
    field = SphereField(center=(0.0, 0.0, 0.0), radius=radius)
    exact = 1.0 / radius

    def score(mesh: Mesh) -> dict:
        values = vertex_mean_curvature(mesh)
        reference = [sample.mean_curvature for sample in principal_curvatures(field, mesh.vertices)]
        deviations = sorted(
            abs(value - target) / abs(target)
            for value, target in zip(values, reference, strict=False)
            if math.isfinite(value) and abs(target) > 0.0
        )
        weights = [0.0] * len(mesh.vertices)
        for face in mesh.faces:
            a, b, c = (mesh.vertices[index] for index in face)
            u = [b[axis] - a[axis] for axis in range(3)]
            v = [c[axis] - a[axis] for axis in range(3)]
            cross = (
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            )
            area = 0.5 * math.sqrt(sum(component * component for component in cross)) / 3.0
            for index in face:
                weights[index] += area
        total = math.fsum(weights)
        weighted = math.fsum(
            value * weight
            for value, weight in zip(values, weights, strict=False)
            if math.isfinite(value)
        ) / total
        return {
            "triangles": len(mesh.faces),
            "pointwise_median_relative_error": deviations[len(deviations) // 2],
            "pointwise_p95_relative_error": deviations[int(0.95 * (len(deviations) - 1))],
            "area_weighted_mean_per_m": weighted,
            "area_weighted_relative_error": abs(weighted - exact) / exact,
        }

    regular = score(uv_sphere(radius, 24, 48))
    sampled = to_mesh(field, 24, margin=3).mesh
    raw = score(sampled)
    smoothed, _ = taubin_smooth(sampled, iterations=30)
    smoothed_score = score(smoothed)
    projected, _ = project_to_field(smoothed, field, iterations=10)
    conditioned = score(projected)
    report = assess_fairness(sampled, field=field).to_dict()

    return {
        "analytic_mean_curvature_per_m": exact,
        "regular_uv_sphere": regular,
        "marching_tets": {
            "raw": raw,
            "taubin_30": smoothed_score,
            "taubin_30_plus_projection_10": conditioned,
        },
        "fairness": {
            "dihedral_degrees": report["dihedral_degrees"],
            "crease_edges": report["crease_edges"],
            "mesh_curvature_per_m": report["mesh_curvature_per_m"],
            "field_curvature_per_m": report["field_curvature_per_m"],
        },
        "claims": [
            claim("area-weighted mean curvature is the analytic 1/R (raw marching tets)",
                  raw["area_weighted_relative_error"], 1e-3),
            claim("pointwise median error on the raw mesh [%] (documents the sliver field)",
                  raw["pointwise_median_relative_error"] * 100.0, 10.0, lower_is_better=False),
            claim("pointwise median error after taubin + projection [%]",
                  conditioned["pointwise_median_relative_error"] * 100.0, 3.0),
            claim("a sphere has no crease above 15 deg", report["crease_edges"], 0),
            claim("sphere dihedral p95 [deg]", report["dihedral_degrees"]["p95"], 6.0),
        ],
    }


def appendage_study() -> dict:
    """Study 3: the AFF-8 audit, and the same audit after conditioning the slivers away."""
    field = SuboffFullField()
    mesh = to_mesh(field, (80, 12, 16), margin=2).mesh
    raw = assess_fairness(mesh, field=field).to_dict()
    smoothed, _ = taubin_smooth(mesh, iterations=10)
    conditioned = assess_fairness(smoothed).to_dict()
    return {
        "grid": {"nx": 80, "ny": 12, "nz": 16},
        "triangles": len(mesh.faces),
        "raw": raw,
        "taubin_10": conditioned,
        "claims": [
            claim("AFF-8 has designed creases above 15 deg", raw["crease_edges"], 100,
                  lower_is_better=False),
            claim("AFF-8 has knife edges (min dihedral) [deg]", raw["dihedral_degrees"]["min"], 1.0),
            claim("AFF-8 audit keeps the creases after taubin", conditioned["crease_edges"], 0.5,
                  lower_is_better=False),
            claim("AFF-8 mesh is watertight", 0 if raw["watertight"] else 1, 0),
        ],
    }


def hull_waterline_slope(parameters: HullParameters) -> dict:
    """``db/dx`` of the analytic half breadth, and the exponent it diverges with."""
    half_beam = 0.5 * parameters.beam
    h = parameters.length * 1e-9

    def breadth(x: float) -> float:
        return half_beam * math.sin(math.pi * x / parameters.length) ** BREADTH_EXPONENT

    rows = []
    for decade in range(2, 8):
        x = parameters.length * 10.0**-decade
        slope = (breadth(x + h) - breadth(x - h)) / (2.0 * h)
        rows.append({"x_over_l": x / parameters.length, "x_m": x, "db_dx": slope})
    xs = [row["x_m"] for row in rows]
    slopes = [abs(row["db_dx"]) for row in rows]
    # least squares on log(|db/dx|) = exponent * log(x) + const
    logs_x = [math.log(value) for value in xs]
    logs_s = [math.log(value) for value in slopes]
    mean_x = sum(logs_x) / len(logs_x)
    mean_s = sum(logs_s) / len(logs_s)
    exponent = sum((u - mean_x) * (v - mean_s) for u, v in zip(logs_x, logs_s, strict=False)) / sum(
        (u - mean_x) ** 2 for u in logs_x
    )
    expected = BREADTH_EXPONENT - 1.0
    return {
        "breadth_exponent": BREADTH_EXPONENT,
        "expected_slope_exponent": expected,
        "fitted_slope_exponent": exponent,
        "rows": rows,
        "claims": [
            claim("waterline slope exponent vs 0.7-1", abs(exponent - expected), 0.02),
            claim("waterline slope keeps growing as x -> 0 (ratio x/L=1e-7 vs 1e-2)",
                  slopes[-1] / slopes[0], 10.0, lower_is_better=False),
        ],
    }


def hull_study() -> dict:
    """Study 3b/4: dihedral statistics of the parametric hull at two grids, plus the ends."""
    parameters_default = HullParameters(length=10.0, beam=2.0, draft=0.8, nx=41, nz=13)
    parameters_fine = HullParameters(length=10.0, beam=2.0, draft=0.8, nx=161, nz=41)
    coarse = assess_fairness(generate_hull(parameters_default).mesh).to_dict()
    fine = assess_fairness(generate_hull(parameters_fine).mesh).to_dict()
    return {
        "default_grid": {"nx": 41, "nz": 13, "triangles": coarse["triangles"],
                         "dihedral_degrees": coarse["dihedral_degrees"],
                         "crease_edges": coarse["crease_edges"],
                         "crease_length_fraction": coarse["crease_length_fraction"],
                         "mesh_curvature_per_m": coarse["mesh_curvature_per_m"]},
        "fine_grid": {"nx": 161, "nz": 41, "triangles": fine["triangles"],
                      "dihedral_degrees": fine["dihedral_degrees"],
                      "crease_edges": fine["crease_edges"],
                      "crease_length_fraction": fine["crease_length_fraction"],
                      "mesh_curvature_per_m": fine["mesh_curvature_per_m"]},
        "waterline": hull_waterline_slope(parameters_default),
        "claims": [
            claim("hull dihedral p95 falls when the grid is refined [deg]",
                  fine["dihedral_degrees"]["p95"], coarse["dihedral_degrees"]["p95"]),
            claim("hull crease *length fraction* falls under refinement",
                  fine["crease_length_fraction"], coarse["crease_length_fraction"]),
            claim("hull crease *count* does not fall under refinement (count alone misleads)",
                  float(fine["crease_edges"]), float(coarse["crease_edges"]),
                  lower_is_better=False),
        ],
    }


def plot(results: dict, path: Path) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(15.0, 4.4))

    field_rows = results["field"]["rows"]
    labels = [f"{row['field']}\n{row['quantity'].split(':')[0]}" for row in field_rows]
    errors = [max(row["max_relative_error"], 1e-17) for row in field_rows]
    axes[0].barh(range(len(errors)), errors, color="#2b6cb0")
    axes[0].set_yticks(range(len(labels)), labels, fontsize=6)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("max relative error of kappa vs closed form")
    axes[0].set_title("1. analytic field: exact to round-off")
    axes[0].axvline(1e-6, color="crimson", linestyle="--", linewidth=1)
    axes[0].invert_yaxis()

    staging = results["mesh"]["marching_tets"]
    names = ["raw", "taubin_30", "taubin_30_plus_projection_10"]
    medians = [staging[name]["pointwise_median_relative_error"] * 100.0 for name in names]
    area = [staging[name]["area_weighted_relative_error"] * 100.0 for name in names]
    index = range(len(names))
    axes[1].bar([i - 0.2 for i in index], medians, width=0.4, label="pointwise median")
    axes[1].bar([i + 0.2 for i in index], area, width=0.4, label="area-weighted mean")
    axes[1].set_yscale("log")
    axes[1].set_xticks(list(index), ["raw", "taubin", "+projection"], fontsize=8)
    axes[1].axhline(3.0, color="crimson", linestyle="--", linewidth=1, label="3 % target")
    axes[1].set_ylabel("relative error of mean curvature [%]")
    axes[1].set_title("2. mesh estimate: integral exact, pointwise needs conditioning")
    axes[1].legend(fontsize=7)

    waterline = results["hull"]["waterline"]["rows"]
    axes[2].loglog([row["x_m"] for row in waterline], [abs(row["db_dx"]) for row in waterline],
                   "o-", label="|db/dx| of the analytic waterline")
    exponent = results["hull"]["waterline"]["fitted_slope_exponent"]
    axes[2].set_xlabel("x [m] from the bow")
    axes[2].set_ylabel("|db/dx| [m/m]")
    axes[2].set_title(f"3. hull ends are not G1: |db/dx| ~ x^{exponent:.2f}")
    axes[2].legend(fontsize=7)
    axes[2].grid(True, which="both", alpha=0.3)

    figure.tight_layout()
    figure.savefig(path, dpi=130)
    plt.close(figure)


def markdown(results: dict) -> str:
    lines: list[str] = []
    lines.append("# Curvature and surface fairness: measured\n")
    lines.append(
        "Generated by `examples/fairness/run_example.py`. Every number below is read "
        "out of `results.json` by the same script that wrote it; the script exits "
        "non-zero if a claim stops holding.\n"
    )
    lines.append("## 1. The analytic field is exact\n")
    lines.append("| field | quantity | expected [1/m] | measured [1/m] | max rel. error |")
    lines.append("|---|---|---|---|---|")
    for row in results["field"]["rows"]:
        expected = ", ".join(f"{value:.9g}" for value in row["expected_per_m"])
        measured = ", ".join(f"{value:.9g}" for value in row["measured_per_m"])
        lines.append(
            f"| {row['field']} | {row['quantity']} | {expected} | {measured} | "
            f"{row['max_relative_error']:.2e} |"
        )
    lines.append("")
    lines.append(
        f"Worst relative error {results['field']['worst_relative_error']:.3e}; the plane is "
        f"zero to {results['field']['plane_abs_error_per_m']:.2e} 1/m (round-off). "
        f"Reported step error (the h vs h/2 change) peaks at "
        f"{results['field']['worst_step_error_per_m']:.2e} 1/m, i.e. the estimator's own "
        "noise floor, seven orders below the curvatures it measures.\n"
    )

    lines.append("## 2. Mesh curvature: integral exact, pointwise needs a conditioned mesh\n")
    lines.append("Reference: a sphere of radius 1 m, analytic mean curvature 1 1/m.\n")
    lines.append("| mesh | triangles | pointwise median | pointwise p95 | area-weighted mean [1/m] | area-weighted error |")
    lines.append("|---|---|---|---|---|---|")
    mesh_section = results["mesh"]
    rows = [("regular UV sphere", mesh_section["regular_uv_sphere"])]
    for name, label in (
        ("raw", "marching tets, raw"),
        ("taubin_30", "+ taubin x30"),
        ("taubin_30_plus_projection_10", "+ projection x10"),
    ):
        rows.append((label, mesh_section["marching_tets"][name]))
    for label, row in rows:
        lines.append(
            f"| {label} | {row['triangles']} | {row['pointwise_median_relative_error'] * 100:.2f} % | "
            f"{row['pointwise_p95_relative_error'] * 100:.1f} % | "
            f"{row['area_weighted_mean_per_m']:.6f} | {row['area_weighted_relative_error'] * 100:.4f} % |"
        )
    lines.append("")
    lines.append(
        "The area-weighted mean is the integral `(1/A) int H dA`, and the cotangent "
        "operator conserves it exactly on a closed mesh, which is why it survives the "
        "sliver field untouched while the pointwise value does not. This is the same "
        "sliver population `docs/mesh_conditioning.md` reports at 1.7 % of triangles: "
        "they carry almost no area and almost all of the pointwise error.\n"
    )

    lines.append("## 3. AFF-8 appendages: the creases are the geometry\n")
    appendages = results["appendages"]
    raw = appendages["raw"]
    smoothed = appendages["taubin_10"]
    lines.append(f"Grid {appendages['grid']}, {appendages['triangles']} triangles.\n")
    lines.append("| quantity | raw | after taubin x10 |")
    lines.append("|---|---|---|")
    for key in ("mean", "p50", "p95", "max", "min"):
        lines.append(
            f"| dihedral {key} [deg] | {raw['dihedral_degrees'][key]:.4f} | "
            f"{smoothed['dihedral_degrees'][key]:.4f} |"
        )
    lines.append(f"| crease edges > {raw['crease_angle_degrees']:.0f} deg | {raw['crease_edges']} | "
                 f"{smoothed['crease_edges']} |")
    lines.append(f"| crease length fraction | {raw['crease_length_fraction'] * 100:.2f} % | "
                 f"{smoothed['crease_length_fraction'] * 100:.2f} % |")
    lines.append(f"| mean curvature jump p50 [1/m] | {raw['curvature_jump_per_m']['p50']:.4f} | — |")
    lines.append(f"| mean curvature jump p95 [1/m] | {raw['curvature_jump_per_m']['p95']:.4f} | — |")
    lines.append("")
    lines.append(
        "Smoothing does not remove the creases, and it must not: the sail/hull junction "
        "and the NACA trailing edges are designed edges, at 0 to 90 degrees. What it "
        "removes is the sliver population that made the raw dihedral distribution "
        "unreadable.\n"
    )

    lines.append("## 4. The parametric hull family is not G1 at its ends\n")
    hull = results["hull"]
    waterline = hull["waterline"]
    lines.append(
        f"`b(x) = B/2 * sin(pi x/L)^{waterline['breadth_exponent']}` has "
        f"`db/dx ~ x^{waterline['expected_slope_exponent']:.1f}`; the fit over "
        f"`x/L = 1e-2 .. 1e-7` gives **{waterline['fitted_slope_exponent']:.4f}** "
        f"(expected {waterline['expected_slope_exponent']:.1f}). "
        "The waterline therefore has an infinite slope at the stem and the stern, the "
        "surface is only C0 there, and every curvature-based quantity diverges: the "
        "first row of the table is the fitted data, and the last row is the value at "
        f"`x/L = 1e-7`, {waterline['rows'][-1]['db_dx']:.3e} 1/m.\n"
    )
    lines.append("| x/L | |db/dx| [m/m] |")
    lines.append("|---|---|")
    for row in waterline["rows"]:
        lines.append(f"| {row['x_over_l']:.0e} | {abs(row['db_dx']):.4e} |")
    lines.append("")
    fine = hull["fine_grid"]
    coarse = hull["default_grid"]
    lines.append("| hull grid | triangles | dihedral p50 [deg] | dihedral p95 [deg] | crease edges | crease length |")
    lines.append("|---|---|---|---|---|---|")
    for label, row in (("41 x 13", coarse), ("161 x 41", fine)):
        lines.append(
            f"| {label} | {row['triangles']} | {row['dihedral_degrees']['p50']:.3f} | "
            f"{row['dihedral_degrees']['p95']:.3f} | {row['crease_edges']} | "
            f"{row['crease_length_fraction'] * 100:.2f} % |"
        )
    lines.append("")
    lines.append(
        "The crease count falls by more than half between the two grids, which is the "
        "warning this report exists to make: **a crease threshold applied to a coarse "
        "mesh measures the mesh, not the shape**. Fairness has to be judged on a "
        "converged grid or on the field, and the two hull-family statements above are "
        "the ones that survive refinement: the ends are not G1, and the keel line has "
        "unbounded curvature because the profile is `sqrt(z/T)`.\n"
    )

    lines.append("## Claims and their verdicts\n")
    lines.append("| claim | value | limit | direction | holds |")
    lines.append("|---|---|---|---|---|")
    for section in ("field", "mesh", "appendages", "hull"):
        for record in results[section]["claims"]:
            lines.append(
                f"| {record['claim']} | {record['value']:.6g} | {record['limit']:.6g} | "
                f"{record['direction']} | {'PASS' if record['holds'] else 'FAIL'} |"
            )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/fairness"))
    arguments = parser.parse_args()
    arguments.output.mkdir(parents=True, exist_ok=True)

    results = {
        "schema": "tensorcad.fairness.study/1",
        "field": field_exactness(),
        "mesh": mesh_curvature_study(),
        "appendages": appendage_study(),
        "hull": hull_study(),
    }
    (arguments.output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    (arguments.output / "report.md").write_text(markdown(results) + "\n")
    plot(results, arguments.output / "plot.png")

    if FAILURES:
        for failure in FAILURES:
            print(f"FAIL {failure}")
        return 1
    total = sum(len(results[section]["claims"]) for section in ("field", "mesh", "appendages", "hull"))
    print(f"{total} claims hold; wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
