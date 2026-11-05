"""Cure-PR-#44 verification: the union of ``main`` and ``agent/m3347``, measured.

Three questions, each answered with a number instead of an assertion of intent:

1. The branch's route is alive on the merged tree -- ``righting_arm_curve`` produces the
   finite-heel righting arms of a closed box, and those arms agree with an independent
   two-dimensional section integration and, inside its own window, with the wall-sided
   formula.
2. ``main``'s route is alive on the merged tree -- the waterplane second moments, the
   asymmetric fore/aft exponents, and the Wigley and transom families still build and
   still agree with their own references.
3. The two routes compose -- ``righting_arm_curve``'s upright ``GM`` equals the one
   ``initial_stability`` reports from the same mesh, and every symbol either side added
   is exported once from ``tensorcad.__all__``.

Writes ``artifacts/m3916/merge_verification.json`` and the report figure.

Everything here is deterministic: no sampling, no random seeds, no tolerance loosening.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import tensorcad as tc
from tensorcad import transom, wigley
from tensorcad.geometry import MeshBuilder
from tensorcad.hull import HullParameters, analytic_hydrostatics, generate_hull
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.stability import initial_stability, righting_arm_curve

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "m3916"

L, B, T, FREEBOARD, KG = 10.0, 2.0, 0.5, 1.5, 0.45

# Symbols each side of the merge introduced; every one of them must survive the union.
BRANCH_SYMBOLS = (
    "GZCurve",
    "GZPoint",
    "SCHEMA_GZ",
    "DEFAULT_HEEL_ANGLES",
    "righting_arm_curve",
)
MAIN_SYMBOLS = (
    "MAX_TRANSOM_FRACTION",
    "TransomHullParameters",
    "TransomModel",
    "analytic_transom_hydrostatics",
    "build_transom_mesh",
    "estimate_transom_discretisation_error",
    "generate_transom_hull",
    "parent_length",
    "sinc_power_coefficients",
    "transom_area_ratio",
    "transom_fraction_for_area_ratio",
    "truncated_taper_integral",
    "truncated_taper_moment",
    "waterplane_inertia_longitudinal",
    "waterplane_inertia_transverse",
)


def box_mesh(length: float, beam: float, draft: float, freeboard: float):
    """Closed box: keel on ``z = -draft``, deck on ``z = +freeboard``, waterline ``z = 0``."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, freeboard):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[p] for p in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi = -0.5 * beam, 0.5 * beam
    bottom, top = -draft, freeboard
    quad((0.0, lo, bottom), (0.0, hi, bottom), (length, hi, bottom), (length, lo, bottom))
    quad((0.0, lo, top), (length, lo, top), (length, hi, top), (0.0, hi, top))
    quad((0.0, lo, bottom), (length, lo, bottom), (length, lo, top), (0.0, lo, top))
    quad((0.0, hi, bottom), (0.0, hi, top), (length, hi, top), (length, hi, bottom))
    quad((0.0, lo, bottom), (0.0, lo, top), (0.0, hi, top), (0.0, hi, bottom))
    quad((length, lo, bottom), (length, hi, bottom), (length, hi, top), (length, lo, top))
    return builder.build()


def section_arm(heel_deg: float, beam: float, draft: float, freeboard: float, kg: float):
    """Exact arm of the rectangular section and its waterplane height.

    The section's four corners are rotated into the heeled frame, clipped by the
    horizontal line ``Z = plane`` and integrated with the shoelace rule; ``plane`` is
    bisected until the clipped area is ``B*T``. Nothing here shares a formula or a code
    path with the three-dimensional divergence integrals under test: this is a 2-D
    polygon cut, the library cuts a triangulated solid.
    """
    phi = math.radians(heel_deg)
    cos_phi, sin_phi = math.cos(phi), math.sin(phi)
    keel = -draft
    corners = [(-0.5 * beam, keel), (0.5 * beam, keel), (0.5 * beam, freeboard), (-0.5 * beam, freeboard)]
    polygon = [(y * cos_phi - z * sin_phi, y * sin_phi + z * cos_phi) for y, z in corners]
    heights = [p[1] for p in polygon]
    low, high = min(heights) - 1.0, max(heights) + 1.0

    def clipped(plane: float):
        kept = []
        for index in range(len(polygon)):
            start, end = polygon[index], polygon[(index + 1) % len(polygon)]
            start_in, end_in = start[1] <= plane, end[1] <= plane
            if start_in:
                kept.append(start)
            if start_in != end_in:
                ratio = (plane - start[1]) / (end[1] - start[1])
                kept.append((start[0] + ratio * (end[0] - start[0]), plane))
        area = moment = 0.0
        for index in range(len(kept)):
            y1, z1 = kept[index]
            y2, z2 = kept[(index + 1) % len(kept)]
            cross = y1 * z2 - y2 * z1
            area += cross
            moment += (y1 + y2) * cross
        return 0.5 * area, moment / 6.0

    for _ in range(200):
        middle = 0.5 * (low + high)
        if clipped(middle)[0] < beam * draft:
            low = middle
        else:
            high = middle
    plane = 0.5 * (low + high)
    area, moment = clipped(plane)
    centre_y = moment / area
    return -(keel + kg) * sin_phi - centre_y, plane


def wall_sided_arm(heel_deg: float, beam: float, draft: float, kg: float) -> float:
    """Classical wall-sided arm ``sin(phi)(BM + BM tan^2(phi)/2 - BG)``."""
    phi = math.radians(heel_deg)
    bm = beam * beam / (12.0 * draft)
    bg = kg - 0.5 * draft
    return math.sin(phi) * (bm + 0.5 * bm * math.tan(phi) ** 2 - bg)


def export_contract() -> dict:
    exported = list(tc.__all__)
    present = set(exported)
    missing_branch = [name for name in BRANCH_SYMBOLS if name not in present]
    missing_main = [name for name in MAIN_SYMBOLS if name not in present]
    unimportable = [name for name in BRANCH_SYMBOLS + MAIN_SYMBOLS if not hasattr(tc, name)]
    duplicates = sorted({name for name in exported if exported.count(name) > 1})
    return {
        "version": tc.__version__,
        "n_exported": len(exported),
        "missing_from_branch_side": missing_branch,
        "missing_from_main_side": missing_main,
        "declared_but_not_importable": unimportable,
        "duplicates": duplicates,
    }


def branch_route() -> dict:
    """m3347's route (finite-heel GZ) on the merged tree, against two references."""
    mesh = box_mesh(L, B, T, FREEBOARD)
    grid = range(0, 181, 1)
    curve = righting_arm_curve(mesh, kg=KG, waterline_z=0.0, angles_deg=grid)

    ref_gz = {point.heel_deg: section_arm(point.heel_deg, B, T, FREEBOARD, KG) for point in curve.points}
    rows = []
    for point in curve.points:
        reference, plane = ref_gz[point.heel_deg]
        rows.append(
            {
                "heel_deg": point.heel_deg,
                "gz": point.gz,
                "gz_section": reference,
                "abs_err": abs(point.gz - reference),
                "gz_wall_sided": wall_sided_arm(point.heel_deg, B, T, KG),
                "waterplane_z": plane,
                "deck_immersed": point.deck_immersed,
                "displacement_volume": point.displacement_volume,
            }
        )
    worst = max(rows, key=lambda row: row["abs_err"])
    # Wall-sided formula is only a reference while the waterline stays on both vertical
    # sides: tan(phi) < 2T/B. One degree past that is already outside the window.
    window = math.degrees(math.atan(2.0 * T / B))
    inside = [row for row in rows if row["heel_deg"] < window]
    worst_inside = max(inside, key=lambda row: abs(row["gz"] - row["gz_wall_sided"]))
    # Past the wall-sided window the formula is expected to fail, so it is measured only
    # up to 60 deg where it is still finite (``tan`` blows up at 90 deg: 1e31 in the doc).
    outside = [row for row in rows if window + 1.0 <= row["heel_deg"] <= 60.0]

    ref_curve = [ref_gz[a][0] for a in range(0, 181)]
    ref_max = max(ref_curve)
    ref_angle = ref_curve.index(ref_max)
    ref_area = math.radians(1.0) * (0.5 * ref_curve[0] + sum(ref_curve[1:-1]) + 0.5 * ref_curve[-1])

    hull_mesh = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8)).mesh
    upright = initial_stability(hull_mesh, kg=0.6)
    hull_curve = righting_arm_curve(
        hull_mesh, kg=0.6, waterline_z=0.0, freeboard=0.4, angles_deg=range(0, 91, 5)
    )
    return {
        "box": {"L": L, "B": B, "T": T, "freeboard": FREEBOARD, "kg": KG,
                "displacement_volume": curve.displacement_volume,
                "declared_displacement_volume": L * B * T},
        "table": [{k: rows[i][k] for k in ("heel_deg", "gz", "gz_section", "abs_err", "gz_wall_sided")}
                  for i in (0, 15, 30, 45, 60, 90, 120, 180)],
        "worst_abs_err": worst["abs_err"],
        "worst_abs_err_heel": worst["heel_deg"],
        "wall_sided_window_deg": window,
        "wall_sided_worst_rel_inside": abs(worst_inside["gz"] - worst_inside["gz_wall_sided"])
        / abs(worst_inside["gz"]),
        "wall_sided_worst_rel_inside_heel": worst_inside["heel_deg"],
        "wall_sided_worst_rel_outside": max(
            abs(row["gz"] - row["gz_wall_sided"]) / abs(row["gz"]) for row in outside
        ),
        "wall_sided_rel_outside_60deg": next(
            abs(row["gz"] - row["gz_wall_sided"]) / abs(row["gz"])
            for row in rows
            if row["heel_deg"] == 60
        ),
        "gz_max": curve.gz_max,
        "gz_max_ref": ref_max,
        "gz_max_rel_err": (curve.gz_max - ref_max) / ref_max,
        "angle_of_gz_max_deg": curve.angle_of_gz_max_deg,
        "angle_of_gz_max_ref_deg": ref_angle,
        "dynamic_stability_area": curve.dynamic_stability_area,
        "dynamic_stability_area_ref": ref_area,
        "dynamic_stability_area_rel_err": (curve.dynamic_stability_area - ref_area) / ref_area,
        "vanishing_angle_deg": curve.vanishing_angle_deg,
        "deck_immersed_first_heel_deg": next(
            (row["heel_deg"] for row in rows if row["deck_immersed"]), None
        ),
        "max_displacement_drift": max(
            abs(row["displacement_volume"] - L * B * T) for row in rows
        ),
        "hull_gm_curve": hull_curve.initial_gm_transverse,
        "hull_gm_initial_stability": upright.gm_transverse,
        "hull_gm_rel_err": (hull_curve.initial_gm_transverse - upright.gm_transverse)
        / upright.gm_transverse,
        "basis": curve.basis,
        "_rows": rows,
        "_hull_curve": hull_curve,
    }


def main_route() -> dict:
    """main's route (hydrostatics second moments, asymmetric taper, Wigley, transom)."""
    barge = compute_hydrostatics(box_mesh(L, B, 0.8, 0.0))
    ref_it, ref_il = L * B**3 / 12.0, B * L**3 / 12.0
    data = {
        "barge": {
            "I_T": barge.transverse_waterplane_inertia,
            "I_T_ref": ref_it,
            "I_T_rel_err": (barge.transverse_waterplane_inertia - ref_it) / ref_it,
            "I_L": barge.longitudinal_waterplane_inertia,
            "I_L_ref": ref_il,
            "I_L_rel_err": (barge.longitudinal_waterplane_inertia - ref_il) / ref_il,
            "keel_z": barge.keel_z,
            "keel_z_ref": -0.8,
        }
    }

    rows = []
    for nx, nz in ((21, 9), (41, 13), (81, 25)):
        model = generate_hull(HullParameters(10.0, 2.0, 0.8, nx, nz))
        mesh = compute_hydrostatics(model.mesh)
        ref = analytic_hydrostatics(10.0, 2.0, 0.8)
        rows.append(
            {
                "nx": nx,
                "nz": nz,
                "I_T_mesh": mesh.transverse_waterplane_inertia,
                "I_T_ref": ref.transverse_waterplane_inertia,
                "I_T_rel_err": (mesh.transverse_waterplane_inertia - ref.transverse_waterplane_inertia)
                / ref.transverse_waterplane_inertia,
                "I_L_mesh": mesh.longitudinal_waterplane_inertia,
                "I_L_ref": ref.longitudinal_waterplane_inertia,
                "I_L_rel_err": (mesh.longitudinal_waterplane_inertia - ref.longitudinal_waterplane_inertia)
                / ref.longitudinal_waterplane_inertia,
                "keel_z_mesh": mesh.keel_z,
            }
        )
    data["hull_vs_closed_form"] = rows

    asym = analytic_hydrostatics(10.0, 2.0, 0.8, fore_exponent=0.5, aft_exponent=0.9)
    sym = analytic_hydrostatics(10.0, 2.0, 0.8)
    data["asymmetric"] = {
        "LCB": asym.longitudinal_centre_buoyancy,
        "LCF": asym.longitudinal_centre_flotation,
        "LCB_symmetric": sym.longitudinal_centre_buoyancy,
        "mid": 5.0,
        "Cb": asym.block_coefficient,
        "I_T_is_nan": math.isnan(asym.transverse_waterplane_inertia),
        "I_L_is_nan": math.isnan(asym.longitudinal_waterplane_inertia),
        "keel_z": asym.keel_z,
        "basis": asym.basis,
    }

    wparams = wigley.WigleyParameters()
    wmodel = wigley.generate_wigley(wparams)
    wmesh = wigley.compute_hydrostatics(wmodel.mesh)
    wref = wigley.wigley_analytic_hydrostatics(wparams)
    data["wigley"] = {
        "Cb_mesh": wmesh.block_coefficient,
        "Cb_ref": wref.block_coefficient,
        "Cb_rel_err": (wmesh.block_coefficient - wref.block_coefficient) / wref.block_coefficient,
        "volume_rel_err": (wmesh.displacement_volume - wref.displacement_volume)
        / wref.displacement_volume,
    }

    tparams = transom.TransomHullParameters()
    tmodel = transom.generate_transom_hull(tparams)
    tmesh = transom.compute_hydrostatics(tmodel.mesh)
    tref = transom.analytic_transom_hydrostatics(tparams)
    data["transom"] = {
        "tau": getattr(tparams, "transom_fraction", None),
        "Cb_mesh": tmesh.block_coefficient,
        "Cb_ref": tref["block_coefficient"],
        "Cb_rel_err": (tmesh.block_coefficient - tref["block_coefficient"])
        / tref["block_coefficient"],
        "volume_rel_err": (tmesh.displacement_volume - tref["displacement_volume"])
        / tref["displacement_volume"],
    }
    return data


def figure(branch: dict, main_data: dict) -> None:
    rows = branch["_rows"]
    heels = [row["heel_deg"] for row in rows]
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.0))

    ax.plot(heels, [row["gz"] for row in rows], "-", lw=2.0, label="merged tree: `righting_arm_curve`")
    ax.plot(heels, [row["gz_section"] for row in rows], "--", lw=1.2, label="2-D section reference")
    ax.plot(heels, [row["gz_wall_sided"] for row in rows], ":", lw=1.2,
            label="wall-sided formula (valid < 26.6° only)")
    ax.axhline(0.0, color="k", lw=0.8)
    ax.set_xlabel("heel [deg]")
    ax.set_ylabel("GZ [m]")
    ax.set_title("m3347 route after the merge: 10×2×(0.5+1.5) m box, KG = 0.45 m")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    inset = ax.inset_axes([0.56, 0.12, 0.40, 0.34])
    inset.semilogy(heels, [max(row["abs_err"], 1e-18) for row in rows], "-", color="#c44e52", lw=1.0)
    inset.set_title("|GZ − section| [m]", fontsize=7)
    inset.tick_params(labelsize=6)
    inset.grid(alpha=0.3)

    grids = main_data["hull_vs_closed_form"]
    nodes = [row["nx"] * row["nz"] for row in grids]
    ax2.plot(nodes, [100.0 * row["I_T_rel_err"] for row in grids], "o-", label="I_T (∫y² dA)")
    ax2.plot(nodes, [100.0 * row["I_L_rel_err"] for row in grids], "s-", label="I_L (∫(x−LCF)² dA)")
    ax2.axhline(0.0, color="k", lw=0.8, ls=":")
    ax2.set_xscale("log")
    ax2.set_xlabel("grid resolution  nx·nz")
    ax2.set_ylabel("mesh vs closed form  [%]")
    ax2.set_title("main route after the merge: waterplane second moments")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(OUT / "merge_verification.png", dpi=110)
    plt.close(fig)


def main() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    contract = export_contract()
    branch = branch_route()
    main_data = main_route()

    payload = {
        "export_contract": contract,
        "branch_route": {k: v for k, v in branch.items() if not k.startswith("_")},
        "main_route": main_data,
    }
    figure(branch, main_data)

    (OUT / "merge_verification.json").write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8"
    )

    print(json.dumps(contract, indent=2))
    print(
        f"box GZ worst |err| vs 2-D section: {branch['worst_abs_err']:.3e} m "
        f"(at {branch['worst_abs_err_heel']} deg)"
    )
    print(
        f"wall-sided rel. inside window: {branch['wall_sided_worst_rel_inside']:.3e} at "
        f"{branch['wall_sided_worst_rel_inside_heel']} deg; outside: "
        f"{branch['wall_sided_worst_rel_outside']:.3f}"
    )
    print(
        f"GZ_max {branch['gz_max']:.9f} vs {branch['gz_max_ref']:.9f} "
        f"({100 * branch['gz_max_rel_err']:+.2e} %), area "
        f"{branch['dynamic_stability_area']:.9f} vs {branch['dynamic_stability_area_ref']:.9f} "
        f"({100 * branch['dynamic_stability_area_rel_err']:+.2e} %)"
    )
    print(
        f"hull GM: curve {branch['hull_gm_curve']:.12f} vs initial_stability "
        f"{branch['hull_gm_initial_stability']:.12f} ({100 * branch['hull_gm_rel_err']:+.2e} %)"
    )
    print(
        f"max displacement drift over the curve: {branch['max_displacement_drift']:.3e} m^3"
    )
    barge = main_data["barge"]
    print(
        f"barge I_T {barge['I_T_rel_err']:+.2e} %, I_L {barge['I_L_rel_err']:+.2e} %, "
        f"keel_z {barge['keel_z']:.6f}"
    )
    for row in main_data["hull_vs_closed_form"]:
        print(
            f"nx={row['nx']:3d} nz={row['nz']:2d}  I_T {100 * row['I_T_rel_err']:+8.4f}%  "
            f"I_L {100 * row['I_L_rel_err']:+8.4f}%  keel_z={row['keel_z_mesh']:.6f}"
        )
    print(
        f"asymmetric LCB/LCF {main_data['asymmetric']['LCB']:.6f} / "
        f"{main_data['asymmetric']['LCF']:.6f} (midship 5.0), inertias nan: "
        f"{main_data['asymmetric']['I_T_is_nan']}, {main_data['asymmetric']['I_L_is_nan']}"
    )
    for family in ("wigley", "transom"):
        print(f"{family} Cb rel err {100 * main_data[family]['Cb_rel_err']:+.4f} %")
    return payload


if __name__ == "__main__":
    main()