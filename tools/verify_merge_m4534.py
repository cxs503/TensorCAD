"""PR-#52 re-open (m4534): the third ``GZ(phi)`` route on the *current* ``main``.

PR #52 (``agent/m3581``) was closed by the owner instead of merged.  It was not a
textual conflict: the branch was cut from 0.20.0 and, between the merge base and
``main``, two other righting-arm modules landed -- ``tensorcad/gz.py`` (0.21.0,
divergence theorem on the rotated closed surface, Newton on the wetted area) and
``tensorcad/righting_arm.py`` (0.22.0, clipped faces in the body frame,
bracket-and-bisect on the monotone ``V(h)``).  Both publish ``righting_arm_curve``,
``RightingArmCurve`` and ``wall_sided_righting_arm``, which is exactly what the
branch's ``tensorcad/righting.py`` wanted to publish too, so a straight merge
would have left one of the three routes unreachable from the package surface.

The cure keeps all three, gives the branch's route the ``section_``/``Section``
prefix at the package surface, and turns the overlap into a check.  This script
measures that check, on the merged tree:

1. **One barge, three routes.**  A rectangular barge is a polyhedron, so each
   route is exact and only the waterline solver and the summation differ; the
   three arms must agree to round-off.
2. **One floating hull, two routes.**  ``sin^0.7 sqrt`` hull at ``nx=41, nz=13``
   with 0.3 m of freeboard: the station-Simpson route against ``gz``'s closed
   surface integral.  Here the two discretisations differ (station grid vs.
   triangles), so the gap is ~1e-7 m rather than round-off -- still far below
   the 15 % band the platform calls unacceptable, and measured, not asserted.
3. **The closed form.**  The wall-sided barge formula, exact up to deck-edge
   immersion, pins the whole chain (half-plane sign, Green moments, Simpson,
   bisection, lever) to round-off.
4. **The integral underneath.**  Station-area volume equals the closed-surface
   tetrahedron volume of the same body to round-off.
5. **Self-convergence.**  40 vs. 80 station intervals on the hull.
6. **A negative control**, so the agreement above is shown to be able to fail: a
   5 mm shift of ``KG`` moves ``GZ(20 deg)`` by ~1.6 mm, four orders of magnitude
   above the inter-route agreement.

Writes ``artifacts/m4534/merge_verification.json`` and ``merge_verification.png``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:  # runnable from a checkout without installing
    sys.path.insert(0, str(REPO_ROOT))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import tensorcad  # noqa: E402
from tensorcad import gz  # noqa: E402
from tensorcad.hull import HullParameters, generate_hull  # noqa: E402
from tensorcad.righting import righting_arm_curve as section_curve  # noqa: E402
from tensorcad.righting import submerged_body, wall_sided_righting_arm  # noqa: E402

OUT = REPO_ROOT / "artifacts" / "m4534"

LENGTH, BEAM, DRAFT, FREEBOARD = 10.0, 2.0, 1.0, 1.0
KG = 0.8
BARGE_ANGLES = tuple(float(angle) for angle in range(0, 46, 5))
HULL_GRID = (41, 13)
HULL_WATERLINE, HULL_KG = -0.3, 0.75
HULL_ANGLES = (0.0, 10.0, 20.0, 30.0, 40.0)
KG_NUDGE = 0.005


def barge_mesh():  # type: ignore[no-untyped-def]
    """Watertight rectangular barge: ``DRAFT`` below ``z = 0``, deck at ``z = FREEBOARD``."""
    half = BEAM / 2.0
    vertices = (
        (0.0, -half, -DRAFT),
        (LENGTH, -half, -DRAFT),
        (LENGTH, half, -DRAFT),
        (0.0, half, -DRAFT),
        (0.0, -half, FREEBOARD),
        (LENGTH, -half, FREEBOARD),
        (LENGTH, half, FREEBOARD),
        (0.0, half, FREEBOARD),
    )
    faces = (
        (0, 2, 1), (0, 3, 2),
        (4, 5, 6), (4, 6, 7),
        (0, 1, 5), (0, 5, 4),
        (1, 2, 6), (1, 6, 5),
        (2, 3, 7), (2, 7, 6),
        (3, 0, 4), (3, 4, 7),
    )
    return tensorcad.Mesh(vertices=vertices, faces=faces)


def hull_mesh(nx: int, nz: int):  # type: ignore[no-untyped-def]
    return generate_hull(HullParameters(length=10.0, beam=2.0, draft=1.0, nx=nx, nz=nz)).mesh


def main() -> dict[str, object]:
    OUT.mkdir(parents=True, exist_ok=True)
    data: dict[str, Any] = {
        "mission": "m4534",
        "pr": 52,
        "branch": "agent/m4534",
        "version": tensorcad.__version__,
        "question": "do the three GZ routes agree on the tree that now carries all of them?",
        "geometry": {
            "barge": {"L": LENGTH, "B": BEAM, "T": DRAFT, "freeboard": FREEBOARD, "kg": KG},
            "hull": {
                "L": 10.0, "B": 2.0, "D": 1.0, "nx": HULL_GRID[0], "nz": HULL_GRID[1],
                "waterline_z": HULL_WATERLINE, "kg": HULL_KG,
            },
        },
        "api_surface": {
            "tensorcad.righting_arm_curve": tensorcad.righting_arm_curve.__module__,
            "tensorcad.RightingArmCurve": tensorcad.RightingArmCurve.__module__,
            "tensorcad.bisection_righting_arm_curve": (
                tensorcad.bisection_righting_arm_curve.__module__
            ),
            "tensorcad.section_righting_arm_curve": (
                tensorcad.section_righting_arm_curve.__module__
            ),
            "tensorcad.SectionRightingArmCurve": tensorcad.SectionRightingArmCurve.__module__,
            "tensorcad.submerged_body": tensorcad.submerged_body.__module__,
            "tensorcad.wall_sided_righting_arm": tensorcad.wall_sided_righting_arm.__module__,
        },
    }

    # 1. one barge, three routes
    barge = barge_mesh()
    ours = section_curve(barge, kg=KG, waterline_z=0.0, heel_angles_deg=BARGE_ANGLES)
    newton = gz.righting_arm_curve(barge, KG, angles_deg=BARGE_ANGLES, waterline_z=0.0)
    bisected = tensorcad.bisection_righting_arm_curve(barge, kg=KG, heels_deg=BARGE_ANGLES)
    barge_rows = []
    worst_three = 0.0
    for index, heel in enumerate(BARGE_ANGLES):
        arm = ours.gz_at(heel)
        gap = max(abs(arm - newton.gz[index]), abs(arm - bisected.gz_at(heel)))
        worst_three = max(worst_three, gap)
        barge_rows.append(
            {
                "heel_deg": heel,
                "section_m": arm,
                "newton_m": newton.gz[index],
                "bisection_m": bisected.gz_at(heel),
                "worst_difference_m": gap,
            }
        )
    data["barge_three_routes"] = {
        "max_difference_m": worst_three,
        "displacement_section_m3": ours.displacement_volume,
        "displacement_newton_m3": newton.displacement_volume,
        "displacement_bisection_m3": bisected.displacement_volume,
        "rows": barge_rows,
    }

    # 2. the closed form, inside its validity window
    gm = DRAFT / 2.0 + BEAM**2 / (12.0 * DRAFT) - KG
    bm = BEAM**2 / (12.0 * DRAFT)
    closed_rows = []
    worst_closed = 0.0
    for heel in BARGE_ANGLES:
        formula = wall_sided_righting_arm(heel, gm=gm, bm=bm)
        gap = abs(ours.gz_at(heel) - formula)
        if formula != 0.0:
            worst_closed = max(worst_closed, gap / abs(formula))
        closed_rows.append({"heel_deg": heel, "gz_m": ours.gz_at(heel),
                            "wall_sided_m": formula, "difference_m": ours.gz_at(heel) - formula})
    data["barge_vs_wall_sided"] = {"worst_relative": worst_closed, "rows": closed_rows}

    # 3. one floating hull, two routes
    hull = hull_mesh(*HULL_GRID)
    hull_section = section_curve(hull, kg=HULL_KG, waterline_z=HULL_WATERLINE,
                                 heel_angles_deg=HULL_ANGLES)
    hull_newton = gz.righting_arm_curve(hull, HULL_KG, angles_deg=HULL_ANGLES,
                                        waterline_z=HULL_WATERLINE)
    hull_rows = []
    worst_hull = 0.0
    for index, heel in enumerate(HULL_ANGLES):
        arm = hull_section.gz_at(heel)
        gap = abs(arm - hull_newton.gz[index])
        worst_hull = max(worst_hull, gap)
        hull_rows.append(
            {
                "heel_deg": heel,
                "section_m": arm,
                "newton_m": hull_newton.gz[index],
                "difference_m": arm - hull_newton.gz[index],
                "volume_residual_section": (
                    hull_section.points[index].submerged_volume
                    / hull_section.displacement_volume
                    - 1.0
                ),
                "volume_error_newton": hull_newton.points[index].volume_error,
            }
        )
    data["hull_two_routes"] = {
        "max_difference_m": worst_hull,
        "displacement_section_m3": hull_section.displacement_volume,
        "displacement_newton_m3": hull_newton.displacement_volume,
        "rows": hull_rows,
    }

    # 4. the integral underneath: station areas vs closed-surface tetrahedra
    # The waterline for this identity is the mesh's own cap (z = 0): below it the two
    # integrations describe the same solid, above it ``compute_hydrostatics`` would still
    # return the whole capped volume while the section integral stops at the waterline.
    body = submerged_body(hull, heel_deg=0.0, waterline_z=0.0)
    hydrostatics = tensorcad.compute_hydrostatics(hull, 0.0)
    data["submerged_body_vs_closed_surface"] = {
        "volume_section_m3": body.volume,
        "volume_closed_surface_m3": hydrostatics.displacement_volume,
        "relative_difference": abs(body.volume - hydrostatics.displacement_volume)
        / hydrostatics.displacement_volume,
    }

    # 5. self-convergence in the station grid
    coarse = section_curve(hull, kg=HULL_KG, waterline_z=HULL_WATERLINE,
                           heel_angles_deg=HULL_ANGLES, station_intervals=40)
    convergence_rows = []
    worst_arm_gap = 0.0
    for heel in HULL_ANGLES:
        fine, rough = hull_section.gz_at(heel), coarse.gz_at(heel)
        gap = abs(fine - rough)
        worst_arm_gap = max(worst_arm_gap, gap)
        convergence_rows.append({"heel_deg": heel, "intervals_80_m": fine,
                                 "intervals_40_m": rough, "difference_m": gap})
    # Relative measures are taken on the curve's peak and on the displacement; ``GZ``
    # passes through zero near 40 deg here, so a per-point relative gap is not a number.
    data["station_convergence"] = {
        "max_abs_arm_difference_m": worst_arm_gap,
        "gz_max_relative": abs(hull_section.gz_max - coarse.gz_max) / hull_section.gz_max,
        "displacement_relative": abs(
            hull_section.displacement_volume - coarse.displacement_volume
        )
        / hull_section.displacement_volume,
        "gz_max_m": {"intervals_80": hull_section.gz_max, "intervals_40": coarse.gz_max},
        "rows": convergence_rows,
    }

    # 6. negative control: the comparison must be able to fail
    nudged = section_curve(hull, kg=HULL_KG + KG_NUDGE, waterline_z=HULL_WATERLINE,
                           heel_angles_deg=(0.0, 20.0))
    shift = abs(nudged.gz_at(20.0) - hull_section.gz_at(20.0))
    data["negative_control"] = {
        "kg_shift_m": KG_NUDGE,
        "gz_20deg_shift_m": shift,
        "ratio_to_inter_route_gap": shift / max(worst_hull, 1e-18),
    }

    # ------------------------------------------------------------------ figure
    figure, axes = plt.subplots(1, 3, figsize=(15.0, 4.4))
    ax = axes[0]
    ax.plot(BARGE_ANGLES, [ours.gz_at(heel) for heel in BARGE_ANGLES], "o", ms=4,
            label="section Simpson (this PR)")
    ax.plot(BARGE_ANGLES, newton.gz, "-", lw=1.2, label="gz: closed surface + Newton")
    ax.plot(BARGE_ANGLES, [bisected.gz_at(heel) for heel in BARGE_ANGLES], "--", lw=1.0,
            label="righting_arm: bisection")
    ax.plot(BARGE_ANGLES, [wall_sided_righting_arm(heel, gm=gm, bm=bm) for heel in BARGE_ANGLES],
            ":", lw=1.4, label="wall-sided closed form")
    ax.set_title(f"barge L={LENGTH:g} B={BEAM:g} T={DRAFT:g} fb={FREEBOARD:g} m")
    ax.set_xlabel("heel [deg]")
    ax.set_ylabel("GZ [m]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)

    ax = axes[1]
    ax.plot(HULL_ANGLES, [hull_section.gz_at(heel) for heel in HULL_ANGLES], "o", ms=5,
            label="section Simpson")
    ax.plot(HULL_ANGLES, hull_newton.gz, "-", lw=1.2, label="gz: closed surface + Newton")
    ax.set_title(f"hull {HULL_GRID[0]}x{HULL_GRID[1]}, waterline {HULL_WATERLINE:g} m, KG {HULL_KG:g} m")
    ax.set_xlabel("heel [deg]")
    ax.set_ylabel("GZ [m]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[2]
    ax.semilogy([row["heel_deg"] for row in barge_rows],
                [max(row["worst_difference_m"], 1e-18) for row in barge_rows], "o-",
                label="barge: 3 routes")
    ax.semilogy(HULL_ANGLES, [max(row["difference_m"], 1e-18) for row in hull_rows], "s-",
                label="hull: section vs Newton")
    ax.axhline(shift, color="crimson", ls=":", lw=1.2,
               label=f"negative control: KG+{KG_NUDGE * 1000:g} mm at 20 deg")
    ax.set_title("|GZ difference| between routes")
    ax.set_xlabel("heel [deg]")
    ax.set_ylabel("[m]")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)

    figure.tight_layout()
    figure.savefig(OUT / "merge_verification.png", dpi=110)
    plt.close(figure)

    (OUT / "merge_verification.json").write_text(json.dumps(data, indent=2) + "\n")
    return data


if __name__ == "__main__":
    payload = main()
    print(json.dumps({k: v for k, v in payload.items() if not isinstance(v, dict) or "rows" not in v},
                     indent=2))
    print("max three-route difference (barge):", payload["barge_three_routes"]["max_difference_m"])
    print("max section-vs-newton difference (hull):", payload["hull_two_routes"]["max_difference_m"])
    print("worst closed-form relative:", payload["barge_vs_wall_sided"]["worst_relative"])
    print("station convergence (gz_max relative):", payload["station_convergence"]["gz_max_relative"])
    print("negative control shift:", payload["negative_control"]["gz_20deg_shift_m"])
    print("wrote", OUT)
