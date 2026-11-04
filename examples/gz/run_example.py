"""Large-angle stability: ``GZ(phi)`` measured against closed forms, and plotted.

Run from the repository root::

    python examples/gz/run_example.py

Writes ``docs/gz/results.json`` (every number below, injected by code), ``docs/gz/tables.md``
(the same numbers as markdown tables, generated rather than typed) and ``docs/gz/plot.png``.
Exits non-zero if any claim stops holding, so ``docs/gz.md`` cannot rot silently.
"""
from __future__ import annotations

import json
import math
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.gz import GZError, gz_curve, immersed_moments, righting_arm
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.stability import waterplane_moments

LENGTH, BEAM, DRAFT, FREEBOARD = 10.0, 2.0, 0.8, 0.4
BARGE_KG = 0.6
CYLINDER_RADIUS, CYLINDER_LENGTH, CYLINDER_G = 1.0, 4.0, -0.3
HULL_WATERLINE, HULL_KG = -0.11, 0.7
ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "gz"


def barge(length: float, beam: float, draft: float, freeboard: float) -> Mesh:
    builder = MeshBuilder()
    corner: dict[tuple[float, float, float], int] = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, freeboard):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    keel, deck = -draft, freeboard
    low, high = -0.5 * beam, 0.5 * beam
    quad((0.0, low, keel), (0.0, high, keel), (length, high, keel), (length, low, keel))
    quad((0.0, low, deck), (length, low, deck), (length, high, deck), (0.0, high, deck))
    quad((0.0, low, keel), (length, low, keel), (length, low, deck), (0.0, low, deck))
    quad((0.0, high, keel), (0.0, high, deck), (length, high, deck), (length, high, keel))
    quad((0.0, low, keel), (0.0, low, deck), (0.0, high, deck), (0.0, high, keel))
    quad((length, low, keel), (length, high, keel), (length, high, deck), (length, low, deck))
    return builder.build()


def cylinder(radius: float, length: float, segments: int) -> Mesh:
    builder = MeshBuilder()
    phase = math.pi / segments  # keep every vertex off the waterline
    rings: list[list[int]] = []
    for x in (0.0, length):
        rings.append(
            [
                builder.add_vertex(
                    (
                        x,
                        radius * math.cos(2.0 * math.pi * index / segments + phase),
                        radius * math.sin(2.0 * math.pi * index / segments + phase),
                    )
                )
                for index in range(segments)
            ]
        )
    for index in range(segments):
        following = (index + 1) % segments
        a, b = rings[0][index], rings[0][following]
        c, d = rings[1][following], rings[1][index]
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)
    for index in range(1, segments - 1):
        builder.add_triangle(rings[0][0], rings[0][index + 1], rings[0][index])
    for index in range(1, segments - 1):
        builder.add_triangle(rings[1][0], rings[1][index], rings[1][index + 1])
    return builder.build()


def wall_sided(phi_deg: float, gm: float, bm: float) -> float:
    phi = math.radians(phi_deg)
    return math.sin(phi) * (gm + 0.5 * bm * math.tan(phi) ** 2)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results: dict[str, object] = {"example": "gz", "version": "0.21.0"}
    claims: list[tuple[str, bool, str]] = []
    decks: dict[str, list[str]] = {}

    # ------------------------------------------------------- 1. wall-sided barge --
    body = barge(LENGTH, BEAM, DRAFT, FREEBOARD)
    bm = BEAM * BEAM / (12.0 * DRAFT)
    gm = bm - (BARGE_KG - 0.5 * DRAFT)
    deck_limit = math.degrees(math.atan(2.0 * FREEBOARD / BEAM))
    rows = []
    for phi in (1.0, 5.0, 10.0, 15.0, 20.0):
        measured = righting_arm(body, phi, kg=BARGE_KG)
        reference = wall_sided(phi, gm, bm)
        rows.append(
            {
                "heel_deg": phi,
                "gz_measured_m": measured,
                "gz_wall_sided_m": reference,
                "relative_error": measured / reference - 1.0,
            }
        )
    worst = max(abs(row["relative_error"]) for row in rows)
    claims.append(("wall-sided barge matches the closed form to 1e-9", worst < 1e-9, f"worst {worst:.2e}"))
    beyond = []
    for phi in (25.0, 30.0, 40.0):
        measured = righting_arm(body, phi, kg=BARGE_KG)
        reference = wall_sided(phi, gm, bm)
        beyond.append(
            {
                "heel_deg": phi,
                "gz_measured_m": measured,
                "gz_wall_sided_m": reference,
                "relative_error": measured / reference - 1.0,
                "formula_valid": False,
            }
        )
    claims.append(
        (
            "past deck-edge immersion the wall-sided formula over-predicts",
            all(row["relative_error"] < -0.01 for row in beyond),
            f"deck edge at {deck_limit:.2f} deg",
        )
    )
    results["barge"] = {
        "length_m": LENGTH,
        "beam_m": BEAM,
        "draft_m": DRAFT,
        "freeboard_m": FREEBOARD,
        "kg_m": BARGE_KG,
        "bm_m": bm,
        "kb_m": 0.5 * DRAFT,
        "gm_m": gm,
        "deck_edge_immersion_deg": deck_limit,
        "within_wall_sided_validity": rows,
        "beyond_wall_sided_validity": beyond,
    }
    decks["barge"] = [
        "| heel [deg] | GZ measured [m] | GZ wall-sided [m] | relative error |",
        "|---:|---:|---:|---:|",
    ] + [
        f"| {row['heel_deg']:.0f} | {row['gz_measured_m']:.9f} | {row['gz_wall_sided_m']:.9f} | {row['relative_error']:+.2e} |"
        for row in rows
    ] + [
        "",
        f"Deck-edge immersion is at `tan(phi) = 2 f / B = {deck_limit:.2f} deg`; past it the formula has no",
        "standing, and the measured arm falls below it:",
        "",
        "| heel [deg] | GZ measured [m] | GZ wall-sided [m] | relative error |",
        "|---:|---:|---:|---:|",
    ] + [
        f"| {row['heel_deg']:.0f} | {row['gz_measured_m']:.9f} | {row['gz_wall_sided_m']:.9f} | {row['relative_error']:+.2e} |"
        for row in beyond
    ]

    # ------------------------------------------------- 2. the barge curve itself --
    curve = gz_curve(body, (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0), kg=BARGE_KG)
    peak_angle, peak_arm = curve.max_arm
    vanishing = curve.vanishing_angle()
    area30, area40 = curve.area_under_curve(30.0), curve.area_under_curve(40.0)
    results["barge"]["curve"] = {
        "max_arm_m": peak_arm,
        "angle_of_max_arm_deg": peak_angle,
        "vanishing_angle_deg": vanishing,
        "area_to_30_deg_m_rad": area30,
        "area_to_40_deg_m_rad": area40,
        "gm_from_slope_m": curve.gm_from_slope(),
        "max_volume_error": curve.max_volume_error,
    }
    claims.append(("barge displacement is held to 1e-12 at every station", curve.max_volume_error < 1e-12, f"{curve.max_volume_error:.1e}"))
    claims.append(("barge curve peaks before it vanishes", 25.0 < peak_angle < 55.0 and peak_arm > 0.1, f"{peak_arm:.4f} m at {peak_angle:.0f} deg"))
    claims.append(("barge vanishing angle is 90 deg", vanishing is not None and abs(vanishing - 90.0) < 2.0, f"{vanishing}"))
    decks["barge_curve"] = [
        "| quantity | value |",
        "|---|---:|",
        f"| maximum arm | {peak_arm:.6f} m at {peak_angle:.0f} deg |",
        f"| vanishing angle | {vanishing:.3f} deg |",
        f"| area under the curve to 30 deg | {area30:.6f} m*rad |",
        f"| area under the curve to 40 deg | {area40:.6f} m*rad |",
        f"| GM from the slope at 5 deg | {curve.gm_from_slope():.6f} m |",
        f"| worst displacement residual | {curve.max_volume_error:.1e} (relative) |",
    ]

    # ------------------------------------------------------- 3. circular cylinder --
    cylinder_rows = []
    errors = []
    for segments in (24, 48, 96):
        tank = cylinder(CYLINDER_RADIUS, CYLINDER_LENGTH, segments)
        worst_error = 0.0
        for phi in (5.0, 20.0, 45.0, 60.0, 90.0, 120.0):
            exact = -CYLINDER_G * math.sin(math.radians(phi))
            measured = righting_arm(tank, phi, kg=CYLINDER_RADIUS + CYLINDER_G)
            worst_error = max(worst_error, abs(measured / exact - 1.0))
        errors.append(worst_error)
        cylinder_rows.append({"segments": segments, "worst_relative_error": worst_error})
    claims.append(("cylinder reproduces the exact arm at every heel", errors[-1] < 2e-2, f"{errors[-1]:.2e} at 96 segments"))
    claims.append(("cylinder error is second order in the segment size", errors[0] / errors[-1] > 8.0, f"{errors[0]:.2e} -> {errors[-1]:.2e}"))
    exact_arms = []
    tank = cylinder(CYLINDER_RADIUS, CYLINDER_LENGTH, 48)
    for phi in (0.0, 10.0, 20.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0):
        exact_arms.append(
            {
                "heel_deg": phi,
                "gz_exact_m": -CYLINDER_G * math.sin(math.radians(phi)),
                "gz_measured_m": righting_arm(tank, phi, kg=CYLINDER_RADIUS + CYLINDER_G),
            }
        )
    results["cylinder"] = {
        "radius_m": CYLINDER_RADIUS,
        "length_m": CYLINDER_LENGTH,
        "centre_of_gravity_above_centre_m": CYLINDER_G,
        "exact_arm_form": "GZ(phi) = -(z_G - z_centre) sin(phi) for every phi",
        "arms": exact_arms,
        "segments": cylinder_rows,
    }
    decks["cylinder"] = [
        "| heel [deg] | GZ exact [m] | GZ 48-gon [m] |",
        "|---:|---:|---:|",
    ] + [
        f"| {row['heel_deg']:.0f} | {row['gz_exact_m']:.9f} | {row['gz_measured_m']:.9f} |" for row in exact_arms
    ] + [
        "",
        "| segments | worst relative error over 5..120 deg |",
        "|---:|---:|",
    ] + [
        f"| {row['segments']} | {row['worst_relative_error']:.2e} |" for row in cylinder_rows
    ]

    # --------------------------------------------- 4. the hull, with freeboard --
    hull_rows = []
    hull_curves = {}
    for nx, nz in ((41, 13), (81, 25)):
        mesh = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=nx, nz=nz)).mesh
        volume, centroid = immersed_moments(mesh, HULL_WATERLINE)
        section = waterplane_moments(mesh, HULL_WATERLINE)
        bm_hull = section.transverse_second_moment / volume
        z_g = min(vertex[2] for vertex in mesh.vertices) + HULL_KG
        gm_waterplane = bm_hull - (z_g - centroid[1])
        hull_curve = gz_curve(
            mesh, (0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0), kg=HULL_KG, waterline_z=HULL_WATERLINE
        )
        hull_curves[(nx, nz)] = hull_curve
        hull_rows.append(
            {
                "nx": nx,
                "nz": nz,
                "faces": len(mesh.faces),
                "displacement_m3": volume,
                "bm_m": bm_hull,
                "gm_waterplane_m": gm_waterplane,
                "gm_from_slope_m": hull_curve.gm_from_slope(),
                "gm_relative_error": hull_curve.gm_from_slope() / gm_waterplane - 1.0,
                "gz_10deg_m": hull_curve.arm_at(10.0),
                "max_arm_m": hull_curve.max_arm[1],
                "angle_of_max_arm_deg": hull_curve.max_arm[0],
                "vanishing_angle_deg": hull_curve.vanishing_angle(),
                "max_volume_error": hull_curve.max_volume_error,
            }
        )
    claims.append(
        (
            "hull GM from the curve slope matches the waterplane integral",
            all(abs(row["gm_relative_error"]) < 1e-2 for row in hull_rows),
            f"worst {max(abs(row['gm_relative_error']) for row in hull_rows):.2e}",
        )
    )
    claims.append(
        (
            "hull arm is grid-converged at 10 deg",
            abs(hull_rows[1]["gz_10deg_m"] / hull_rows[0]["gz_10deg_m"] - 1.0) < 0.1,
            f"{hull_rows[0]['gz_10deg_m']:.5f} -> {hull_rows[1]['gz_10deg_m']:.5f} m",
        )
    )
    results["hull"] = {"waterline_z_m": HULL_WATERLINE, "kg_m": HULL_KG, "grids": hull_rows}
    decks["hull"] = [
        "| grid | faces | displacement [m3] | BM [m] | GM (waterplane) [m] | GM (curve slope) [m] | relative error | GZ(10 deg) [m] | vanishing angle [deg] |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ] + [
        f"| {row['nx']}x{row['nz']} | {row['faces']} | {row['displacement_m3']:.6f} | {row['bm_m']:.6f} | "
        f"{row['gm_waterplane_m']:.6f} | {row['gm_from_slope_m']:.6f} | {row['gm_relative_error']:+.2e} | "
        f"{row['gz_10deg_m']:.6f} | {row['vanishing_angle_deg']:.2f} |"
        for row in hull_rows
    ]

    # ------------------------------------- 5. the degenerate body must be refused --
    degenerate = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=21, nz=9)).mesh
    try:
        gz_curve(degenerate, (0.0, 10.0), kg=HULL_KG)
    except GZError as error:
        refusal = str(error)
        claims.append(("a body at exactly full submergence is refused", True, refusal[:80]))
    else:  # pragma: no cover
        refusal = ""
        claims.append(("a body at exactly full submergence is refused", False, "no error raised"))
    results["degenerate"] = {"refused_with": refusal}

    # ---------------------------------------------------------------- plot --
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    barge_angles = [row["heel_deg"] for row in rows] + [row["heel_deg"] for row in beyond]
    barge_measured = [row["gz_measured_m"] for row in rows] + [row["gz_measured_m"] for row in beyond]
    barge_reference = [row["gz_wall_sided_m"] for row in rows] + [row["gz_wall_sided_m"] for row in beyond]
    left.plot(barge_angles, barge_measured, "o-", color="tab:blue", label="barge, clipped (measured)")
    left.plot(barge_angles, barge_reference, "--", color="tab:blue", alpha=0.6, label="barge, wall-sided formula")
    left.axvline(deck_limit, color="tab:grey", ls=":", lw=1.0)
    left.text(deck_limit + 1.0, 0.02, "deck edge immerses", fontsize=7, color="tab:grey")
    cylinder_angles = [row["heel_deg"] for row in exact_arms]
    left.plot(cylinder_angles, [row["gz_exact_m"] for row in exact_arms], "--", color="tab:orange", alpha=0.6, label="cylinder, exact")
    left.plot(cylinder_angles, [row["gz_measured_m"] for row in exact_arms], "o", ms=3.5, color="tab:orange", label="cylinder, clipped (measured)")
    hull_curve = hull_curves[(81, 25)]
    left.plot(hull_curve.angles_deg, hull_curve.arms, "s-", ms=3.5, color="tab:green", label="hull 81x25, 0.11 m freeboard")
    left.axhline(0.0, color="black", lw=0.8)
    left.set_xlabel("heel angle [deg]")
    left.set_ylabel("righting arm GZ [m]")
    left.set_title("Righting-arm curves against closed forms")
    left.grid(alpha=0.3)
    left.legend(fontsize=7)

    right.semilogy(
        [row["heel_deg"] for row in rows],
        [max(abs(row["relative_error"]), 1e-16) for row in rows],
        "o-",
        color="tab:blue",
        label="barge vs wall-sided",
    )
    right.semilogy(
        [row["heel_deg"] for row in exact_arms if row["heel_deg"] > 0.0],
        [
            abs(row["gz_measured_m"] / row["gz_exact_m"] - 1.0) if row["gz_exact_m"] else float("nan")
            for row in exact_arms
            if row["heel_deg"] > 0.0
        ],
        "s-",
        color="tab:orange",
        label="cylinder vs exact",
    )
    right.semilogy(
        [row["segments"] for row in cylinder_rows],
        [row["worst_relative_error"] for row in cylinder_rows],
        "^-",
        color="tab:red",
        label="cylinder worst error vs segments",
    )
    right.axhline(1e-3, color="tab:grey", ls=":", lw=1.0)
    right.text(30.0, 1.2e-3, "0.1 %", fontsize=7, color="tab:grey")
    right.set_xlabel("heel angle [deg]  /  segments (red curve)")
    right.set_ylabel("relative error vs reference")
    right.set_title("Error against the analytic reference")
    right.grid(alpha=0.3, which="both")
    right.legend(fontsize=7)
    figure.tight_layout()
    figure.savefig(OUT / "plot.png", dpi=140)
    plt.close(figure)

    # --------------------------------------------------------------- output --
    results["claims"] = [{"claim": name, "holds": bool(ok), "detail": detail} for name, ok, detail in claims]
    (OUT / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True))
    sections = [
        "# Generated tables (do not edit by hand)", "",
        f"Written by `examples/gz/run_example.py`. Geometry: barge `{LENGTH} x {BEAM} x {DRAFT} m` with "
        f"`{FREEBOARD} m` freeboard and `KG = {BARGE_KG} m`; cylinder `R = {CYLINDER_RADIUS} m`, "
        f"`L = {CYLINDER_LENGTH} m`, `G` `{-CYLINDER_G} R` below the axis; hull `{LENGTH} x {BEAM} x {DRAFT} m` "
        f"at waterline `z = {HULL_WATERLINE} m` ({-HULL_WATERLINE:.2f} m freeboard) and `KG = {HULL_KG} m`.", "",
    ]
    for title, deck in (
        ("Barge against the wall-sided formula", decks["barge"]),
        ("Barge righting-arm curve", decks["barge_curve"]),
        ("Circular cylinder against the exact arm", decks["cylinder"]),
        ("Parametric hull with freeboard", decks["hull"]),
    ):
        sections += [f"## {title}", ""] + deck + [""]
    sections += ["## Claims", "", "| claim | holds | detail |", "|---|---|---|"] + [
        f"| {item['claim']} | {'yes' if item['holds'] else 'no'} | {item['detail']} |" for item in results["claims"]
    ]
    (OUT / "tables.md").write_text("\n".join(sections) + "\n")

    print(json.dumps({item["claim"]: item["holds"] for item in results["claims"]}, indent=2))
    failures = [item["claim"] for item in results["claims"] if not item["holds"]]
    if failures:
        print("CLAIMS FAILED:", failures, file=sys.stderr)
        return 1
    print(f"all {len(claims)} claims hold; wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())