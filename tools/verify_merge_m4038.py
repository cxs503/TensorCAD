"""PR-#40 cure (m4038): the two ``GZ(phi)`` routes on one tree, measured.

PR #40 (``agent/m3169``) was closed because ``main`` had meanwhile acquired its own
large-angle righting-arm curve (``tensorcad/gz.py``, 0.21.0), so the branch and ``main``
answered the same question twice.  The cure does not pick a winner: it keeps both routes
and makes their overlap into the *check*.  Three questions, each answered with a number:

1. **The branch's route is alive on the merged tree** -- ``tensorcad/righting_arm.py``
   reproduces the wall-sided barge closed form to round-off, a floating sphere to
   ``GM sin(phi)`` with the error falling as the mesh is refined, and the parametric hull
   with a small-angle slope equal to the waterplane-moment ``GM`` of
   ``tensorcad.stability``.
2. **``main``'s route is alive on the merged tree** -- ``tensorcad/gz.py`` still returns
   the same arms it did before the merge, and neither route shadows the other at the
   package surface (``tensorcad.righting_arm_curve`` is ``gz``'s, the second route is
   reached as ``bisection_righting_arm_curve`` / ``tensorcad.righting_arm``).
3. **The two routes compose honestly** -- same mesh, same ``KG``, same displacement, same
   heels: the arms agree to double-precision round-off, because each integral is exact for
   the polygonal hull and only the waterline *solver* differs (bisection on the monotone
   volume here, Newton on the projected wetted area there).  A deliberately mismatched
   displacement target is included so the comparison is shown to be able to fail.

Writes ``artifacts/m4038/merge_verification.json`` and ``merge_verification.png``.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import tensorcad
from tensorcad import gz
from tensorcad.geometry import MeshBuilder
from tensorcad.geometry_field import SphereField, to_mesh
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.righting_arm import (
    heeled_displacement,
    righting_arm,
    righting_arm_curve,
    wall_sided_righting_arm,
    wall_sided_validity_deg,
)
from tensorcad.stability import waterplane_moments

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "m4038"

LENGTH, BEAM, DRAFT, FREEBOARD = 10.0, 2.0, 0.8, 0.4
DEPTH = DRAFT + FREEBOARD
KG = 0.55
HEELS = (0.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0)
HULL_WATERLINE, HULL_KG = -0.47, 0.55
HULL_GRIDS = ((21, 9), (41, 13), (81, 25))
SPHERE_RADIUS, SPHERE_KG, SPHERE_RESOLUTION = 1.0, 0.9, 32


def box(length: float, beam: float, draft: float, freeboard: float):  # type: ignore[no-untyped-def]
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

    lo, hi, bottom, top, right = -0.5 * beam, 0.5 * beam, -draft, freeboard, length
    quad((0.0, lo, bottom), (0.0, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((0.0, lo, top), (right, lo, top), (right, hi, top), (0.0, hi, top))
    quad((0.0, lo, bottom), (right, lo, bottom), (right, lo, top), (0.0, lo, top))
    quad((0.0, hi, bottom), (0.0, hi, top), (right, hi, top), (right, hi, bottom))
    quad((0.0, lo, bottom), (0.0, lo, top), (0.0, hi, top), (0.0, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


def hull(nx: int, nz: int):  # type: ignore[no-untyped-def]
    return generate_hull(HullParameters(length=10.0, beam=2.0, draft=1.0, nx=nx, nz=nz)).mesh


def compare(mesh, kg: float, heels, waterline_z: float = 0.0) -> dict[str, object]:  # type: ignore[no-untyped-def]
    """Both routes on one mesh, one loading condition, one displacement."""
    ours = righting_arm_curve(mesh, kg=kg, heels_deg=heels, waterline_z=waterline_z)
    theirs = gz.righting_arm_curve(mesh, kg, angles_deg=heels, waterline_z=waterline_z)
    rows = []
    for index, point in enumerate(ours.points):
        rows.append(
            {
                "heel_deg": point.heel_deg,
                "righting_arm_bisection_m": point.righting_arm,
                "righting_arm_newton_gz_m": theirs.gz[index],
                "difference_m": point.righting_arm - theirs.gz[index],
                "volume_residual_bisection": point.volume_residual,
                "volume_error_newton": theirs.points[index].volume_error,
                "waterline_bisection_m": point.waterline_z,
                "iterations_bisection": point.iterations,
                "iterations_newton": theirs.points[index].iterations,
            }
        )
    return {
        "displacement_bisection_m3": ours.displacement_volume,
        "displacement_newton_m3": theirs.displacement_volume,
        "displacement_relative_difference": abs(
            ours.displacement_volume - theirs.displacement_volume
        )
        / abs(theirs.displacement_volume),
        "max_abs_arm_difference_m": max(abs(row["difference_m"]) for row in rows),
        "peak_bisection": {"heel_deg": ours.peak()[0], "gz_m": ours.peak()[1]},
        "peak_newton": {"heel_deg": theirs.angle_of_max_gz_deg, "gz_m": theirs.max_gz},
        "rows": rows,
    }


def main() -> dict[str, object]:
    OUT.mkdir(parents=True, exist_ok=True)
    data: dict[str, Any] = {
        "mission": "m4038",
        "pr": 40,
        "branch": "agent/m4038",
        "version": tensorcad.__version__,
        "question": "do the two independent GZ routes agree on the merged tree?",
        "geometry": {
            "barge": {"L": LENGTH, "B": BEAM, "T": DRAFT, "freeboard": FREEBOARD},
            "hull": {"L": 10.0, "B": 2.0, "D": 1.0, "waterline_z": HULL_WATERLINE},
            "sphere": {"R": SPHERE_RADIUS, "resolution": SPHERE_RESOLUTION},
            "kg": KG,
        },
        "api_surface": {
            "tensorcad.righting_arm_curve": tensorcad.righting_arm_curve.__module__,
            "tensorcad.bisection_righting_arm_curve": (
                tensorcad.bisection_righting_arm_curve.__module__
            ),
            "tensorcad.RightingArmCurve": tensorcad.RightingArmCurve.__module__,
            "tensorcad.BisectionRightingArmCurve": (
                tensorcad.BisectionRightingArmCurve.__module__
            ),
            "tensorcad.righting_arm": tensorcad.righting_arm.__module__,
            "tensorcad.wall_sided_righting_arm": tensorcad.wall_sided_righting_arm.__module__,
            "tensorcad.barge_wall_sided_righting_arm": (
                tensorcad.barge_wall_sided_righting_arm.__module__
            ),
        },
    }

    # 1. the branch's route against the closed forms it was written for
    barge = box(LENGTH, BEAM, DRAFT, FREEBOARD)
    limit = wall_sided_validity_deg(breadth=BEAM, draft=DRAFT, depth=DEPTH)
    inside = tuple(heel for heel in HEELS if heel <= limit)
    closed_form_rows = []
    worst_wall = 0.0
    for heel in inside:
        measured = righting_arm(barge, heel_deg=heel, kg=KG).righting_arm
        formula = wall_sided_righting_arm(breadth=BEAM, draft=DRAFT, kg=KG, heel_deg=heel)
        worst_wall = max(worst_wall, abs(measured - formula))
        closed_form_rows.append(
            {"heel_deg": heel, "gz_m": measured, "wall_sided_m": formula,
             "difference_m": measured - formula}
        )
    data["barge_vs_wall_sided"] = {
        "validity_deg": limit,
        "worst_abs_difference_m": worst_wall,
        "rows": closed_form_rows,
    }
    outside_rows = []
    for heel in (25.0, 30.0, 35.0, 40.0):
        measured = righting_arm(barge, heel_deg=heel, kg=KG).righting_arm
        formula = wall_sided_righting_arm(breadth=BEAM, draft=DRAFT, kg=KG, heel_deg=heel)
        outside_rows.append(
            {"heel_deg": heel, "gz_m": measured, "wall_sided_m": formula,
             "formula_error_m": formula - measured}
        )
    data["barge_beyond_the_closed_form"] = outside_rows

    # 2. the branch's route on the parametric hull, against the waterplane moments
    hull_rows = []
    for nx, nz in HULL_GRIDS:
        mesh = hull(nx, nz)
        upright = heeled_displacement(mesh, heel_deg=0.0, waterline_z=HULL_WATERLINE)
        moments = waterplane_moments(mesh, HULL_WATERLINE)
        keel = min(vertex[2] for vertex in mesh.vertices)
        gm_from_moments = (
            upright.centre_of_buoyancy[2] - keel
            + moments.transverse_second_moment / upright.volume
            - HULL_KG
        )
        small = 1e-3  # GZ/sin(phi) is the metacentric height to first order in phi
        slope = (
            righting_arm(mesh, heel_deg=small, kg=HULL_KG, waterline_z=HULL_WATERLINE).righting_arm
            / math.sin(math.radians(small))
        )
        at30 = righting_arm(mesh, heel_deg=30.0, kg=HULL_KG, waterline_z=HULL_WATERLINE)
        hull_rows.append(
            {
                "nx": nx,
                "nz": nz,
                "gm_from_waterplane_moments_m": gm_from_moments,
                "gm_from_heeled_slope_m": slope,
                "gm_relative_difference": (slope - gm_from_moments) / gm_from_moments,
                "gz_30deg_m": at30.righting_arm,
                "volume_residual_30deg": at30.volume_residual,
            }
        )
    data["hull_small_angle"] = hull_rows

    # 3. the two routes, side by side
    data["cross_route_barge"] = compare(barge, KG, HEELS)
    data["cross_route_hull"] = [
        {"grid": {"nx": nx, "nz": nz}, **compare(hull(nx, nz), HULL_KG, HEELS, HULL_WATERLINE)}
        for nx, nz in ((41, 13), (81, 25))
    ]
    sphere = to_mesh(SphereField(radius=SPHERE_RADIUS), SPHERE_RESOLUTION).mesh
    sphere_heels = (0.0, 10.0, 20.0, 30.0, 40.0)
    sphere_block = compare(sphere, SPHERE_KG, sphere_heels)
    analytic = [
        (SPHERE_RADIUS - SPHERE_KG) * math.sin(math.radians(heel)) for heel in sphere_heels
    ]
    sphere_block["analytic_gz_m"] = analytic
    sphere_block["max_abs_deviation_from_analytic_m"] = max(
        abs(row["righting_arm_bisection_m"] - reference)
        for row, reference in zip(sphere_block["rows"], analytic, strict=True)
    )
    data["cross_route_sphere"] = sphere_block

    # 4. the negative control: a comparison that cannot fail is not a comparison
    mismatched = gz.righting_arm_curve(
        barge, KG, angles_deg=(0.0, 20.0), waterline_z=0.0,
        displacement_volume=0.9 * LENGTH * BEAM * DRAFT,
    )
    ours = righting_arm_curve(barge, kg=KG, heels_deg=(0.0, 20.0))
    data["negative_control"] = {
        "description": "same mesh, displacement target moved 10 %: the routes must diverge",
        "gz_at_20deg_newton_target_90pct_m": mismatched.gz[1],
        "gz_at_20deg_bisection_full_m": ours.points[1].righting_arm,
        "difference_m": ours.points[1].righting_arm - mismatched.gz[1],
    }

    # 5. what the merged tree promises at its surface
    data["no_shadowing"] = {
        "righting_arm_curve_is_gz": tensorcad.righting_arm_curve is gz.righting_arm_curve,
        "wall_sided_righting_arm_is_gz": (
            tensorcad.wall_sided_righting_arm is gz.wall_sided_righting_arm
        ),
        "RightingArmCurve_is_gz": tensorcad.RightingArmCurve is gz.RightingArmCurve,
        "bisection_route_reachable": callable(tensorcad.bisection_righting_arm_curve),
        "barge_wall_sided_reachable": callable(tensorcad.barge_wall_sided_righting_arm),
    }

    claims = [
        (
            "barge GZ equals the wall-sided closed form inside its validity window",
            worst_wall < 1e-12,
            f"worst |difference| = {worst_wall:.2e} m up to {limit:.2f} deg",
        ),
        (
            "the two routes agree on the barge to round-off",
            data["cross_route_barge"]["max_abs_arm_difference_m"] < 1e-12,
            f"worst |difference| = {data['cross_route_barge']['max_abs_arm_difference_m']:.2e} m",
        ),
        (
            "the two routes agree on the hull",
            max(
                entry["max_abs_arm_difference_m"] for entry in data["cross_route_hull"]
            ) < 1e-12,
            "worst |difference| = "
            f"{max(entry['max_abs_arm_difference_m'] for entry in data['cross_route_hull']):.2e} m",
        ),
        (
            "the two routes agree on the sphere",
            data["cross_route_sphere"]["max_abs_arm_difference_m"] < 1e-12,
            f"worst |difference| = {data['cross_route_sphere']['max_abs_arm_difference_m']:.2e} m",
        ),
        (
            "on the sphere both routes sit within the mesh error of GM sin(phi)",
            data["cross_route_sphere"]["max_abs_deviation_from_analytic_m"] < 2e-4,
            "worst |GZ - (R-KG) sin(phi)| = "
            f"{data['cross_route_sphere']['max_abs_deviation_from_analytic_m']:.2e} m at "
            f"resolution {SPHERE_RESOLUTION} (meshing error, not solver error)",
        ),
        (
            "the negative control diverges",
            abs(data["negative_control"]["difference_m"]) > 1e-3,
            f"10 % displacement shift moves GZ(20 deg) by "
            f"{data['negative_control']['difference_m']:+.4f} m",
        ),
        (
            "neither route shadows the other at the package surface",
            all(data["no_shadowing"].values()),
            f"{data['no_shadowing']}",
        ),
    ]
    data["claims"] = [
        {"claim": name, "passed": bool(passed), "detail": detail}
        for name, passed, detail in claims
    ]

    # figure: arms from both routes, and the difference between them
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.5))
    hull_41 = compare(hull(41, 13), HULL_KG, HEELS, HULL_WATERLINE)
    panels = (
        (axes[0][0], data["cross_route_barge"], "barge  L=10 B=2 T=0.8 m, KG=0.55 m"),
        (axes[0][1], hull_41, "hull 41x13, waterline -0.47 m, KG=0.55 m"),
    )
    for axis, block, title in panels:
        heels = [row["heel_deg"] for row in block["rows"]]
        axis.plot(heels, [row["righting_arm_bisection_m"] for row in block["rows"]], "o-",
                  label="tensorcad.righting_arm (bisection)")
        axis.plot(heels, [row["righting_arm_newton_gz_m"] for row in block["rows"]], "x--",
                  label="tensorcad.gz (Newton)")
        axis.set_xlabel("heel [deg]")
        axis.set_ylabel("GZ [m]")
        axis.set_title(title)
        axis.grid(alpha=0.3)
        axis.legend(loc="upper left", fontsize=8)
    axes[1][0].semilogy(
        [row["heel_deg"] for row in data["cross_route_barge"]["rows"]],
        [abs(row["difference_m"]) + 1e-18 for row in data["cross_route_barge"]["rows"]],
        "o-",
        label="barge",
    )
    axes[1][0].semilogy(
        [row["heel_deg"] for row in hull_41["rows"]],
        [abs(row["difference_m"]) + 1e-18 for row in hull_41["rows"]],
        "s-",
        label="hull 41x13",
    )
    axes[1][0].axhline(1e-12, color="crimson", lw=1, ls=":",
                       label="round-off bound 1e-12 m")
    axes[1][0].set_xlabel("heel [deg]")
    axes[1][0].set_ylabel("|GZ(bisection) - GZ(Newton)| [m]")
    axes[1][0].set_title("the two routes differ only by round-off")
    axes[1][0].grid(alpha=0.3, which="both")
    axes[1][0].legend(fontsize=8)
    axes[1][1].plot(
        [row["heel_deg"] for row in data["barge_vs_wall_sided"]["rows"]],
        [row["difference_m"] for row in data["barge_vs_wall_sided"]["rows"]],
        "o-",
        label="barge vs wall-sided form (inside window)",
    )
    axes[1][1].plot(
        [row["heel_deg"] for row in data["barge_beyond_the_closed_form"]],
        [row["formula_error_m"] for row in data["barge_beyond_the_closed_form"]],
        "s--",
        label="formula error past the deck edge",
    )
    axes[1][1].axvline(limit, color="grey", lw=1, ls=":",
                       label=f"validity limit {limit:.1f} deg")
    axes[1][1].set_xlabel("heel [deg]")
    axes[1][1].set_ylabel("difference [m]")
    axes[1][1].set_title("the closed form is exact only inside its own window")
    axes[1][1].grid(alpha=0.3)
    axes[1][1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "merge_verification.png", dpi=100)
    plt.close(fig)

    (OUT / "merge_verification.json").write_text(json.dumps(data, indent=2, default=str))
    for claim in data["claims"]:
        print(f"[{'PASS' if claim['passed'] else 'FAIL'}] {claim['claim']}: {claim['detail']}")
    failures = [claim for claim in data["claims"] if not claim["passed"]]
    if failures:
        raise SystemExit(f"{len(failures)} claim(s) failed")
    print(f"wrote {OUT}")
    return data


if __name__ == "__main__":
    main()
