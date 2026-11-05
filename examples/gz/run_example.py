"""Large-angle stability measured against a closed form and against itself.

Everything on this page is *measured*: the script exits non-zero if one of its own
claims stops holding, so the numbers in ``docs/gz.md`` cannot drift away from the code
the way a hand-copied table does.

What it does, in order:

1. **A wall-sided box.** Its righting arm has the textbook closed form
   ``sin(phi) * (GM + BM tan^2(phi) / 2)``, and because a box's faces are exactly
   representable the clipped-volume integration must reproduce it to round-off.
2. **The parametric hull, two independent routes.** The 3-D heeled-displacement
   integration and the classical cross curves from station cuts share nothing but the
   mesh; their difference is the script's honest error bar on ``GZ``.
3. **The metacentric-height cross-check.** The curve's slope at the origin is a
   *different* calculation from ``GM = KB + BM - KG`` (one differentiates a heeled
   volume, the other divides a waterplane's second moment by the volume), so their
   agreement ties this module to :mod:`tensorcad.stability` without sharing code with it.
4. **A stability ladder over ``KG``**, down to the loading where the upright position
   stops being stable and the curve has no range to report.

Run it with ``python examples/gz/run_example.py --output docs/gz``. matplotlib is
optional: the JSON is written either way and the PNGs are skipped.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from tensorcad import __version__
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.gz import (
    cross_curve_righting_arm_curve,
    keel_height,
    righting_arm,
    righting_arm_curve,
    submerged_properties,
    wall_sided_righting_arm,
)
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.stability import waterplane_moments

HEEL_GRID = (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 75.0, 90.0)

#: Above this heel the station-cut route stops resolving the same geometry as the 3-D
#: clip (the deck edge is immersed and the sections it integrates are nearly degenerate),
#: so the two-route agreement is claimed below it and the difference is still reported.
ROUTE_COMPARISON_LIMIT = 60.0

BOX_LENGTH, BOX_BEAM, BOX_DEPTH, BOX_DRAFT, BOX_KG = 40.0, 10.0, 8.0, 4.0, 3.0
BOX_VOLUME = BOX_LENGTH * BOX_BEAM * BOX_DRAFT
BOX_BM = BOX_BEAM**2 / (12.0 * BOX_DRAFT)
BOX_KB = 0.5 * BOX_DRAFT
BOX_GM = BOX_KB + BOX_BM - BOX_KG

HULL_LENGTH, HULL_BEAM, HULL_DRAFT = 12.0, 2.4, 0.9
HULL_WATERLINE = -0.31
HULL_KG = 0.55


def box(length: float, beam: float, depth: float) -> Mesh:
    """Axis-aligned closed box with the keel on ``z = 0``; outward winding."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (0.0, depth):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    left, right = 0.0, length
    lo, hi = -0.5 * beam, 0.5 * beam
    bottom, top = 0.0, depth
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


def hull(resolution: tuple[int, int]) -> Mesh:
    nx, nz = resolution
    return generate_hull(
        HullParameters(length=HULL_LENGTH, beam=HULL_BEAM, draft=HULL_DRAFT, nx=nx, nz=nz)
    ).mesh


def metacentric_height(mesh: Mesh, kg: float) -> dict[str, float]:
    submerged = submerged_properties(mesh, (0.0, 0.0, 1.0), HULL_WATERLINE)
    keel = keel_height(mesh)
    kb = submerged.centroid[2] - keel
    section = waterplane_moments(mesh, HULL_WATERLINE)
    bm = section.transverse_second_moment / submerged.volume
    return {
        "waterline_z_m": HULL_WATERLINE,
        "displacement_volume_m3": submerged.volume,
        "kb_m": kb,
        "bm_m": bm,
        "km_m": kb + bm,
        "kg_m": kg,
        "gm_m": kb + bm - kg,
    }


def slope_at_zero(mesh: Mesh, volume: float, kg: float, delta: float = 0.25) -> float:
    starboard = righting_arm(
        mesh, delta, kg=kg, displacement_volume=volume, waterline_z=HULL_WATERLINE
    ).righting_arm
    port = righting_arm(
        mesh, -delta, kg=kg, displacement_volume=volume, waterline_z=HULL_WATERLINE
    ).righting_arm
    return (starboard - port) / (2.0 * math.radians(delta))


def relative_error(measured: float, reference: float) -> float:
    if reference == 0.0:
        return 0.0 if measured == 0.0 else float("inf")
    return (measured - reference) / abs(reference)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="docs/gz", help="directory for results.json and figures")
    parser.add_argument("--no-figures", action="store_true", help="skip matplotlib output")
    arguments = parser.parse_args(argv)
    output = Path(arguments.output)
    output.mkdir(parents=True, exist_ok=True)

    claims: list[dict[str, object]] = []

    def claim(name: str, holds: bool, detail: str) -> None:
        claims.append({"claim": name, "holds": bool(holds), "detail": detail})

    # -- 1. wall-sided box against the closed form --------------------------
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    box_curve = righting_arm_curve(
        mesh, HEEL_GRID, kg=BOX_KG, displacement_volume=BOX_VOLUME
    )
    box_rows = []
    worst = 0.0
    # The wall-sided closed form assumes the waterline still meets the two vertical
    # sides: no deck edge immersed, no bilge emerged. For this box that is
    # tan(phi) < T / (B/2), i.e. phi < 38.66 deg, and the claim is only made inside it.
    box_valid_limit = math.degrees(math.atan(BOX_DRAFT / (0.5 * BOX_BEAM)))
    for point in box_curve.points:
        reference = wall_sided_righting_arm(point.heel_degrees, BOX_GM, BOX_BM)
        error = relative_error(point.righting_arm, reference)
        if reference != 0.0 and point.heel_degrees < box_valid_limit:
            worst = max(worst, abs(error))
        box_rows.append(
            {
                "heel_degrees": point.heel_degrees,
                "wall_sided_closed_form_m": reference,
                "measured_righting_arm_m": point.righting_arm,
                "relative_error": error,
                "submerged_volume_m3": point.submerged_volume,
                "volume_relative_error": point.volume_error,
                "transverse_centre_buoyancy_m": point.transverse_centre_buoyancy,
                "centre_buoyancy_above_keel_m": point.centre_buoyancy_above_keel,
                "waterline_iterations": point.iterations,
            }
        )
    claim(
        "box GZ reproduces the wall-sided closed form to <1e-8 relative",
        worst < 1e-8,
        f"worst relative error {worst:.3e} over the heel angles below the wall-sided "
        f"limit {box_valid_limit:.2f} deg (beyond it the closed form does not apply)",
    )
    volume_balance = max(abs(row["volume_relative_error"]) for row in box_rows)
    claim(
        "box heeled displacement stays on the target volume",
        volume_balance < 1e-9,
        f"worst relative volume error {volume_balance:.3e}",
    )

    # -- 2. the hull, both routes -------------------------------------------
    hull_curves: dict[str, object] = {}
    hull_rows = []
    worst_route = 0.0
    for resolution in ((41, 13), (81, 25), (161, 41)):
        hull_mesh = hull(resolution)
        volume = submerged_properties(hull_mesh, (0.0, 0.0, 1.0), HULL_WATERLINE).volume
        volume_route = righting_arm_curve(
            hull_mesh, HEEL_GRID, kg=HULL_KG, displacement_volume=volume, waterline_z=HULL_WATERLINE
        )
        cross_route = cross_curve_righting_arm_curve(
            hull_mesh, HEEL_GRID, kg=HULL_KG, displacement_volume=volume, waterline_z=HULL_WATERLINE
        )
        key = f"{resolution[0]}x{resolution[1]}"
        hull_curves[key] = {
            "faces": len(hull_mesh.faces),
            "volume_route": volume_route.to_dict(),
            "cross_curve_route": cross_route.to_dict(),
        }
        for first, second in zip(volume_route.points, cross_route.points, strict=True):
            if first.heel_degrees > ROUTE_COMPARISON_LIMIT:
                continue
            reference = first.righting_arm if abs(first.righting_arm) > 1e-9 else 0.0
            difference = second.righting_arm - first.righting_arm
            if reference != 0.0:
                worst_route = max(worst_route, abs(difference / reference))
            else:
                worst_route = max(worst_route, abs(difference))
            hull_rows.append(
                {
                    "resolution": key,
                    "heel_degrees": first.heel_degrees,
                    "gz_volume_route_m": first.righting_arm,
                    "gz_cross_curve_route_m": second.righting_arm,
                    "difference_m": difference,
                    "difference_relative": None if reference == 0.0 else difference / reference,
                    "submerged_volume_m3": first.submerged_volume,
                    "target_volume_m3": first.target_volume,
                    "waterline_iterations": first.iterations,
                }
            )
    claim(
        "3-D heeled-displacement and cross-curve routes agree within 0.5%",
        worst_route < 5e-3,
        f"worst difference {worst_route:.3e} (relative, or absolute metres near the origin) "
        f"for heels up to {ROUTE_COMPARISON_LIMIT:.0f} deg",
    )

    # -- 3. slope at zero against the metacentric height --------------------
    slope_rows = []
    worst_slope = 0.0
    for resolution in ((41, 13), (81, 25), (161, 41)):
        hull_mesh = hull(resolution)
        measurements = metacentric_height(hull_mesh, HULL_KG)
        volume = measurements["displacement_volume_m3"]
        slope = slope_at_zero(hull_mesh, volume, HULL_KG)
        error = relative_error(slope, measurements["gm_m"])
        worst_slope = max(worst_slope, abs(error))
        row = dict(measurements)
        row["resolution"] = f"{resolution[0]}x{resolution[1]}"
        row["curvature_slope_dgz_dphi"] = slope
        row["relative_error"] = error
        slope_rows.append(row)
    claim(
        "GZ slope at the origin equals GM = KB + BM - KG within 0.5%",
        worst_slope < 5e-3,
        f"worst relative error {worst_slope:.3e} on three resolutions",
    )

    # -- 4. symmetry ---------------------------------------------------------
    symmetry_mesh = hull((81, 25))
    symmetry_volume = submerged_properties(symmetry_mesh, (0.0, 0.0, 1.0), HULL_WATERLINE).volume
    symmetry_error = 0.0
    for heel in (5.0, 10.0, 20.0, 30.0):
        starboard = righting_arm(
            symmetry_mesh, heel, kg=HULL_KG, displacement_volume=symmetry_volume, waterline_z=HULL_WATERLINE
        )
        port = righting_arm(
            symmetry_mesh, -heel, kg=HULL_KG, displacement_volume=symmetry_volume, waterline_z=HULL_WATERLINE
        )
        symmetry_error = max(
            symmetry_error,
            abs(starboard.righting_arm + port.righting_arm),
            abs(starboard.transverse_centre_buoyancy + port.transverse_centre_buoyancy),
        )
    claim(
        "GZ is an odd function of the heel angle to <1e-8 m",
        symmetry_error < 1e-8,
        f"worst |GZ(phi) + GZ(-phi)| = {symmetry_error:.3e} m",
    )

    # -- 5. a KG ladder ------------------------------------------------------
    base = metacentric_height(hull((81, 25)), HULL_KG)
    ladder_rows = []
    for offset in (0.0, 0.1, 0.2, 0.3, 0.5, 1.0):
        kg = HULL_KG + offset
        curve = righting_arm_curve(
            hull((41, 13)),
            HEEL_GRID,
            kg=kg,
            displacement_volume=submerged_properties(
                hull((41, 13)), (0.0, 0.0, 1.0), HULL_WATERLINE
            ).volume,
            waterline_z=HULL_WATERLINE,
        )
        ladder_rows.append(
            {
                "kg_m": kg,
                "gm_m": base["km_m"] - kg,
                "maximum_righting_arm_m": curve.maximum_righting_arm,
                "angle_of_maximum_degrees": curve.angle_of_maximum,
                "range_of_stability_degrees": curve.range_of_stability,
                "area_under_curve_0_30_m_rad": curve.area_under_curve(30.0),
                "area_under_curve_0_40_m_rad": curve.area_under_curve(40.0),
                "gz_at_5_degrees_m": curve.righting_arms[1],
            }
        )
    upright_stable = ladder_rows[0]["gm_m"] > 0.0 and ladder_rows[0]["gz_at_5_degrees_m"] > 0.0
    claim(
        "a stable loading has positive arms and a vanishing one does not",
        upright_stable
        and ladder_rows[-1]["gm_m"] < 0.0
        and ladder_rows[-1]["gz_at_5_degrees_m"] < 0.0,
        f"GM {ladder_rows[0]['gm_m']:+.3f} m -> GZ(5deg) {ladder_rows[0]['gz_at_5_degrees_m']:+.4f} m; "
        f"GM {ladder_rows[-1]['gm_m']:+.3f} m -> GZ(5deg) {ladder_rows[-1]['gz_at_5_degrees_m']:+.4f} m",
    )
    stable_rows = [row for row in ladder_rows if row["gm_m"] > 0.05]
    monotone = all(
        stable_rows[index + 1]["maximum_righting_arm_m"] < stable_rows[index]["maximum_righting_arm_m"]
        for index in range(len(stable_rows) - 1)
    )
    claim(
        "raising KG lowers the maximum righting arm for the loadings still stable",
        monotone,
        "; ".join(f"KG={row['kg_m']:.2f} -> GZmax={row['maximum_righting_arm_m']:.4f}" for row in ladder_rows),
    )

    # -- write the evidence --------------------------------------------------
    payload = {
        "schema": "tensorcad.gz_example/1",
        "version": __version__,
        "heel_grid_degrees": list(HEEL_GRID),
        "box": {
            "length_m": BOX_LENGTH,
            "beam_m": BOX_BEAM,
            "depth_m": BOX_DEPTH,
            "draft_m": BOX_DRAFT,
            "kg_m": BOX_KG,
            "displacement_volume_m3": BOX_VOLUME,
            "kb_m": BOX_KB,
            "bm_m": BOX_BM,
            "gm_m": BOX_GM,
            "curve": box_curve.to_dict(),
            "rows": box_rows,
            "worst_relative_error": worst,
            "closed_form_valid_below_degrees": box_valid_limit,
        },
        "hull": {
            "length_m": HULL_LENGTH,
            "beam_m": HULL_BEAM,
            "draft_m": HULL_DRAFT,
            "waterline_z_m": HULL_WATERLINE,
            "kg_m": HULL_KG,
            "metacentric_height": metacentric_height(hull((81, 25)), HULL_KG),
            "curves": hull_curves,
            "route_comparison": hull_rows,
            "worst_route_difference": worst_route,
        },
        "slope_check": {"rows": slope_rows, "worst_relative_error": worst_slope},
        "symmetry": {"worst_absolute_error_m": symmetry_error},
        "kg_ladder": ladder_rows,
        "claims": claims,
    }
    (output / "results.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {output / 'results.json'}")
    write_report(output, payload)

    if not arguments.no_figures:
        try:
            write_figures(output, payload)
        except ImportError as error:  # pragma: no cover - optional dependency
            print(f"matplotlib unavailable ({error}); figures skipped", file=sys.stderr)

    failed = [entry for entry in claims if not entry["holds"]]
    for entry in claims:
        status = "PASS" if entry["holds"] else "FAIL"
        print(f"[{status}] {entry['claim']}: {entry['detail']}")
    if failed:
        print(f"{len(failed)} claim(s) failed", file=sys.stderr)
        return 1
    return 0


def write_report(output: Path, payload: dict) -> None:
    """A markdown summary whose numbers are read from ``payload``, never re-typed."""
    box = payload["box"]
    lines = [
        f"# Righting-arm evidence (TensorCAD {payload['version']})",
        "",
        "Generated by `examples/gz/run_example.py`; every number below is read from",
        "`results.json` in this directory, and the script exits non-zero if a claim stops",
        "holding.",
        "",
        "## Claims",
        "",
        "| claim | verdict | detail |",
        "|---|---|---|",
    ]
    for entry in payload["claims"]:
        lines.append(f"| {entry['claim']} | {'PASS' if entry['holds'] else 'FAIL'} | {entry['detail']} |")
    lines += [
        "",
        "## Wall-sided box against the closed form",
        "",
        f"L B D T = {box['length_m']:.1f} x {box['beam_m']:.1f} x {box['depth_m']:.1f} x "
        f"{box['draft_m']:.1f} m, KG = {box['kg_m']:.2f} m, "
        f"KB = {box['kb_m']:.4f} m, BM = {box['bm_m']:.4f} m, GM = {box['gm_m']:.4f} m.",
        "",
        "| heel [deg] | closed form GZ [m] | measured GZ [m] | rel. error | volume rel. error |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in box["rows"]:
        lines.append(
            f"| {row['heel_degrees']:.1f} | {row['wall_sided_closed_form_m']:+.6f} | "
            f"{row['measured_righting_arm_m']:+.6f} | {row['relative_error']:+.2e} | "
            f"{row['volume_relative_error']:+.2e} |"
        )
    lines += [
        "",
        f"Closed form valid below {box['closed_form_valid_below_degrees']:.2f} deg for this box.",
        "",
        "## Parametric hull: two independent routes",
        "",
        "| resolution | heel [deg] | GZ, 3-D clip [m] | GZ, cross curves [m] | difference [m] |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in payload["hull"]["route_comparison"]:
        lines.append(
            f"| {row['resolution']} | {row['heel_degrees']:.1f} | {row['gz_volume_route_m']:+.6f} | "
            f"{row['gz_cross_curve_route_m']:+.6f} | {row['difference_m']:+.2e} |"
        )
    lines += [
        "",
        "## Slope at the origin against the metacentric height",
        "",
        "| resolution | volume [m^3] | KB [m] | BM [m] | GM [m] | dGZ/dphi [m/rad] | rel. error |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in payload["slope_check"]["rows"]:
        lines.append(
            f"| {row['resolution']} | {row['displacement_volume_m3']:.6f} | {row['kb_m']:.6f} | "
            f"{row['bm_m']:.6f} | {row['gm_m']:.6f} | {row['curvature_slope_dgz_dphi']:.6f} | "
            f"{row['relative_error']:+.2e} |"
        )
    lines += [
        "",
        "## KG ladder",
        "",
        "| KG [m] | GM [m] | max GZ [m] | heel of max [deg] | range [deg] | A(30) [m rad] | A(40) [m rad] |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in payload["kg_ladder"]:
        rng = "none" if row["range_of_stability_degrees"] is None else f"{row['range_of_stability_degrees']:.2f}"
        lines.append(
            f"| {row['kg_m']:.2f} | {row['gm_m']:+.4f} | {row['maximum_righting_arm_m']:.4f} | "
            f"{row['angle_of_maximum_degrees']:.0f} | {rng} | {row['area_under_curve_0_30_m_rad']:.4f} | "
            f"{row['area_under_curve_0_40_m_rad']:.4f} |"
        )
    lines.append("")
    (output / "report.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {output / 'report.md'}")


def write_figures(output: Path, payload: dict) -> None:  # pragma: no cover - plotting
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # 1. box: measured against the closed form, and the heeled geometry behind it
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    rows = payload["box"]["rows"]
    heels = [row["heel_degrees"] for row in rows]
    axes[0].plot(
        heels, [row["wall_sided_closed_form_m"] for row in rows], "-", lw=2.0,
        label="wall-sided closed form",
    )
    axes[0].plot(
        heels, [row["measured_righting_arm_m"] for row in rows], "o", ms=5.0,
        label="TensorCAD clipped-volume GZ",
    )
    axes[0].set_xlabel("heel angle [deg]")
    axes[0].set_ylabel("GZ [m]")
    axes[0].set_title(
        f"box {BOX_LENGTH:.0f}x{BOX_BEAM:.0f}x{BOX_DEPTH:.0f} m, T={BOX_DRAFT:.0f} m, "
        f"GM={BOX_GM:.3f} m"
    )
    axes[0].legend(loc="upper left")
    axes[0].grid(alpha=0.3)

    residuals = [
        abs(row["relative_error"]) for row in rows if abs(row["wall_sided_closed_form_m"]) > 0.0
    ]
    axes[1].semilogy(
        [row["heel_degrees"] for row in rows if abs(row["wall_sided_closed_form_m"]) > 0.0],
        residuals,
        "s-",
    )
    axes[1].axhline(3e-2, color="tab:red", ls="--", lw=1.0, label="3% acceptance line")
    axes[1].set_xlabel("heel angle [deg]")
    axes[1].set_ylabel("|measured - closed form| / |closed form|")
    axes[1].set_title("the clip is exact for flat faces: residual is round-off")
    axes[1].legend(loc="lower right")
    axes[1].grid(alpha=0.3, which="both")
    figure.tight_layout()
    figure.savefig(output / "box_vs_closed_form.png", dpi=140)
    plt.close(figure)

    # 2. hull: both routes, plus the KG ladder
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    for key, entry in payload["hull"]["curves"].items():
        points = entry["volume_route"]["points"]
        axes[0].plot(
            [point["heel_degrees"] for point in points],
            [point["righting_arm_m"] for point in points],
            "-", lw=1.6, label=f"3-D clip, {key}",
        )
        cross_points = entry["cross_curve_route"]["points"]
        axes[0].plot(
            [point["heel_degrees"] for point in cross_points],
            [point["righting_arm_m"] for point in cross_points],
            "x", ms=5.0, alpha=0.8,
        )
    axes[0].plot([], [], "kx", label="cross curves (station cuts)")
    axes[0].set_xlabel("heel angle [deg]")
    axes[0].set_ylabel("GZ [m]")
    axes[0].set_title(
        f"parametric hull {HULL_LENGTH:.0f} m, T_w={-HULL_WATERLINE:.2f} m, KG={HULL_KG:.2f} m"
    )
    axes[0].legend(loc="upper left", fontsize=8)
    axes[0].grid(alpha=0.3)

    ladder = payload["kg_ladder"]
    axes[1].plot(
        [row["gm_m"] for row in ladder],
        [row["maximum_righting_arm_m"] for row in ladder],
        "o-",
    )
    for row in ladder:
        axes[1].annotate(
            f"KG={row['kg_m']:.2f}",
            (row["gm_m"], row["maximum_righting_arm_m"]),
            textcoords="offset points",
            xytext=(4, 6),
            fontsize=8,
        )
    axes[1].axvline(0.0, color="tab:red", ls="--", lw=1.0, label="GM = 0")
    axes[1].set_xlabel("GM [m]")
    axes[1].set_ylabel("max GZ [m]")
    axes[1].set_title("stability ladder over KG")
    axes[1].legend(loc="upper left")
    axes[1].grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(output / "hull_curves_and_ladder.png", dpi=140)
    plt.close(figure)

    # 3. convergence of the two routes and of the GM slope
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.0))
    resolutions = list(payload["hull"]["curves"])
    faces = [payload["hull"]["curves"][key]["faces"] for key in resolutions]
    axes[0].loglog(
        faces,
        [payload["hull"]["curves"][key]["volume_route"]["maximum_righting_arm_m"] for key in resolutions],
        "o-",
        label="max GZ, 3-D clip",
    )
    axes[0].loglog(
        faces,
        [
            payload["hull"]["curves"][key]["cross_curve_route"]["maximum_righting_arm_m"]
            for key in resolutions
        ],
        "s--",
        label="max GZ, cross curves",
    )
    for face, key in zip(faces, resolutions, strict=True):
        axes[0].annotate(key, (face, payload["hull"]["curves"][key]["volume_route"]["maximum_righting_arm_m"]),
                         textcoords="offset points", xytext=(4, -10), fontsize=8)
    axes[0].set_xlabel("mesh faces")
    axes[0].set_ylabel("max GZ [m]")
    axes[0].set_title("both routes converge with the mesh")
    axes[0].legend()
    axes[0].grid(alpha=0.3, which="both")

    slope_rows = payload["slope_check"]["rows"]
    axes[1].semilogy(
        [row["gm_m"] for row in slope_rows],
        [abs(row["relative_error"]) for row in slope_rows],
        "o-",
    )
    axes[1].axhline(3e-2, color="tab:red", ls="--", lw=1.0, label="3% acceptance line")
    for row in slope_rows:
        axes[1].annotate(row["resolution"], (row["gm_m"], abs(row["relative_error"])),
                         textcoords="offset points", xytext=(4, 6), fontsize=8)
    axes[1].set_xlabel("GM = KB + BM - KG [m]")
    axes[1].set_ylabel("|slope(GZ at 0) - GM| / GM")
    axes[1].set_title("the curve's slope at the origin is the metacentric height")
    axes[1].legend(loc="lower right")
    axes[1].grid(alpha=0.3, which="both")
    figure.tight_layout()
    figure.savefig(output / "convergence.png", dpi=140)
    plt.close(figure)
    print(f"wrote 3 figures to {output}")


if __name__ == "__main__":
    raise SystemExit(main())
