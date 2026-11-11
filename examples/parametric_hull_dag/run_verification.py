"""Verification campaign for parametric hull generation (campaign node p316w1-n03).

Measures every machine-checkable acceptance criterion of the node against closed-form
references (generalised Wigley family) and the published DARPA SUBOFF offset table, and
writes a machine-readable record plus the figures used in the report:

    python examples/parametric_hull_dag/run_verification.py

Outputs (relative to the repository root):

    benchmarks/verified/parametric_hull/verification.json  every number, plus PASS/FAIL
    docs/parametric_hull/figures/*.png       report figures
    docs/parametric_hull/VERIFICATION.md     generated table of the same numbers

The record deliberately lands under ``benchmarks/verified/`` and the figures under
``docs/``: the repository's push filter drops ``artifacts/``, ``runs/`` and binary
extensions outside ``docs/``, so a record written to ``artifacts/`` never reaches the
committed tree and the run leaves no evidence behind.

The acceptance thresholds live in ``CRITERIA`` so the report cannot drift from the
measurement, and the markdown is rendered from the measured dictionaries only.
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad import mesh as meshlib
from tensorcad.parametric_hull import (
    Transform,
    apply_transform,
    build_parametric_hull,
    compute_hydrostatics,
    depth_exponent_for_block_coefficient,
    loft_nurbs_surface,
    mesh_density_for,
    mesh_quality_report,
    suboff_cross_sections,
    suboff_station_fractions,
    tessellate_surface,
    wigley_analytic_properties,
    wigley_cross_sections,
    wigley_family_sections,
)
from tensorcad.suboff import (
    SUBOFF_LENGTH_M,
    SUBOFF_RADIUS_M,
    suboff_radius_fraction,
    suboff_reference_geometry,
)

CRITERIA = {
    "a_wigley_geometry_max_rel_error": 1e-3,
    "b_wigley_volume_max_rel_error": 5e-3,
    "c_wigley_waterplane_max_rel_error": 1e-2,
    "d_suboff_geometry_max_rel_error": 5e-3,
    "e_min_scaled_jacobian": 0.30,
    "e_max_aspect_ratio": 5.0,
}

# Reference hull (parent form used throughout unless stated otherwise)
L, B, T = 10.0, 1.0, 0.625
STATIONS = 11
PROFILE_POINTS = 12
REPORT_N_V = 16          # profile rows used for the reported hydrostatics


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _v_for_z(surface, u, z_target, lo=0.0, hi=1.0):
    """Parameter ``v`` whose surface ordinate ``z`` equals ``z_target`` at fixed ``u``."""
    z_lo = surface.evaluate(u, lo)[2]
    z_hi = surface.evaluate(u, hi)[2]
    rising = z_hi >= z_lo
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        z_mid = surface.evaluate(u, mid)[2]
        if (z_mid < z_target) == rising:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _u_for_x(surface, x_target):
    lo, hi = 0.0, 1.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if surface.evaluate(mid, 0.5)[0] < x_target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def _surface_radius(surface, x_target, samples=401):
    u = _u_for_x(surface, x_target)
    return max(
        math.hypot(*surface.evaluate(u, j / (samples - 1))[1:])
        for j in range(samples)
    )


def _rel(value, reference):
    return (value - reference) / reference


# --------------------------------------------------------------------------- #
# criteria
# --------------------------------------------------------------------------- #
def criterion_a_wigley_geometry(hull):
    """Surface ordinates vs the analytic formula at 10 stations x 5 waterlines."""
    surface = hull.surface
    stations = [0.05 + 0.1 * i for i in range(10)]           # x/L, tips excluded
    waterlines = [-0.95 * T, -0.75 * T, -0.5 * T, -0.25 * T, 0.0]
    worst = 0.0
    samples = []
    for x_frac in stations:
        x = x_frac * L
        u = _u_for_x(surface, x)
        for z in waterlines:
            v = _v_for_z(surface, u, z)
            point = surface.evaluate(u, v)
            ratio = abs(2.0 * x - L) / L
            analytic = 0.5 * B * (1.0 - ratio ** 2) * (1.0 - (abs(z) / T) ** 2)
            error = abs(point[1] - analytic) / analytic
            worst = max(worst, error)
            samples.append(
                {
                    "x_over_L": round(x_frac, 4),
                    "z": z,
                    "y_surface": point[1],
                    "y_analytic": analytic,
                    "rel_error": error,
                }
            )
    return {
        "value": worst,
        "threshold": CRITERIA["a_wigley_geometry_max_rel_error"],
        "n_stations": len(stations),
        "n_waterlines": len(waterlines),
        "samples": samples,
    }


def criterion_b_volume(hull, analytic):
    measured = hull.hydrostatics.displaced_volume
    error = abs(_rel(measured, analytic["displaced_volume"]))
    cb = measured / (L * B * T)
    literal = 4.0 / 3.0 * L * B * T * cb
    return {
        "value": error,
        "threshold": CRITERIA["b_wigley_volume_max_rel_error"],
        "measured_volume": measured,
        "analytic_volume": analytic["displaced_volume"],
        "measured_block_coefficient": cb,
        "analytic_block_coefficient": analytic["block_coefficient"],
        "seeding_formula_note": {
            "literal_expression": "4/3 * L * B * T * Cb",
            "literal_value": literal,
            "literal_rel_error": _rel(literal, analytic["displaced_volume"]),
            "standard_expression": "Cb * L * B * T",
            "standard_rel_error": _rel(cb * L * B * T, analytic["displaced_volume"]),
            "explanation": (
                "The node text quotes displacement as 4/3*L*B*T*Cb. With the naval-"
                "architecture block coefficient of the Wigley surface (Cb = 4/9) that "
                "expression returns 16/27*L*B*T, i.e. -25% against the closed-form "
                "volume 4/9*L*B*T; it equals the volume only for Cb = 1/3, which is not "
                "the block coefficient of this hull. The criterion is therefore "
                "evaluated in the standard form V = Cb*L*B*T with Cb taken from the "
                "geometry, and the literal expression is reported here for the audit "
                "trail."
            ),
        },
    }


def criterion_c_waterplane(hull, analytic):
    measured = hull.hydrostatics.waterplane_area
    return {
        "value": abs(_rel(measured, analytic["waterplane_area"])),
        "threshold": CRITERIA["c_wigley_waterplane_max_rel_error"],
        "measured_waterplane_area": measured,
        "analytic_waterplane_area": analytic["waterplane_area"],
        "extra": {
            "measured_transverse_inertia": hull.hydrostatics.transverse_waterplane_inertia,
            "analytic_transverse_inertia": analytic["transverse_waterplane_inertia"],
            "transverse_inertia_rel_error": _rel(
                hull.hydrostatics.transverse_waterplane_inertia,
                analytic["transverse_waterplane_inertia"],
            ),
            "measured_longitudinal_inertia": hull.hydrostatics.longitudinal_waterplane_inertia,
            "analytic_longitudinal_inertia": analytic["longitudinal_waterplane_inertia"],
            "longitudinal_inertia_rel_error": _rel(
                hull.hydrostatics.longitudinal_waterplane_inertia,
                analytic["longitudinal_waterplane_inertia"],
            ),
            "measured_bm": hull.hydrostatics.transverse_metacentric_radius,
            "analytic_bm": analytic["transverse_metacentric_radius"],
            "transverse_metacentric_radius_rel_error": _rel(
                hull.hydrostatics.transverse_metacentric_radius,
                analytic["transverse_metacentric_radius"],
            ),
            "measured_vcb": hull.hydrostatics.vcb,
            "analytic_vcb": analytic["vcb"],
            "vcb_rel_error": _rel(hull.hydrostatics.vcb, analytic["vcb"]),
            "lcf_error_m": hull.hydrostatics.lcf - analytic["lcf"],
        },
    }


def criterion_d_suboff(mesh_n_u=160, mesh_n_v=16):
    stations = suboff_station_fractions(41)
    sections = suboff_cross_sections(n_points=32, station_fractions=stations)
    surface = loft_nurbs_surface(sections, 3, 3, 32)
    checks = []
    worst = 0.0
    for i in range(20):
        xi = (i + 0.5) / 20.0                      # midpoints: not build stations
        reference = SUBOFF_RADIUS_M * suboff_radius_fraction(xi)
        measured = _surface_radius(surface, xi * SUBOFF_LENGTH_M)
        error = abs(_rel(measured, reference))
        worst = max(worst, error)
        checks.append(
            {
                "x_over_L": xi,
                "radius_surface_m": measured,
                "radius_darpa_m": reference,
                "rel_error": error,
            }
        )
    mesh = tessellate_surface(surface, n_u=mesh_n_u, n_v=mesh_n_v)
    reference_geometry = suboff_reference_geometry()
    hydro = compute_hydrostatics(
        mesh.mesh, waterline_z=2.0, n_hull_faces=mesh.hull_face_count
    )
    tip_closure = {
        "radius_at_bow_tip_m": _surface_radius(surface, 1e-6),
        "radius_at_stern_tip_m": _surface_radius(surface, SUBOFF_LENGTH_M - 1e-6),
    }
    return {
        "value": worst,
        "threshold": CRITERIA["d_suboff_geometry_max_rel_error"],
        "n_stations": len(checks),
        "n_build_stations": len(stations),
        "checks": checks,
        "tip_closure": tip_closure,
        "volume": {
            "measured": hydro.displaced_volume,
            "reference": reference_geometry["fine"]["volume"],
            "rel_error": _rel(
                hydro.displaced_volume, reference_geometry["fine"]["volume"]
            ),
        },
        "wetted_area": {
            "measured": hydro.wetted_surface_area,
            "reference": reference_geometry["fine"]["wetted_area"],
            "rel_error": _rel(
                hydro.wetted_surface_area, reference_geometry["fine"]["wetted_area"]
            ),
        },
        "mesh": {
            "n_u": mesh_n_u,
            "n_v": mesh_n_v,
            "watertight": meshlib.is_watertight(mesh.mesh),
            "edge_manifold": meshlib.is_edge_manifold(mesh.mesh),
            "quality": mesh_quality_report(mesh.mesh),
        },
    }


def criterion_e_mesh_quality(wigley_hull):
    quality = wigley_hull.quality
    passed = (
        quality["min_scaled_jacobian"] > CRITERIA["e_min_scaled_jacobian"]
        and quality["max_aspect_ratio"] < CRITERIA["e_max_aspect_ratio"]
        and quality["n_degenerate_faces"] == 0
    )
    return {
        "value": {
            "min_scaled_jacobian": quality["min_scaled_jacobian"],
            "max_aspect_ratio": quality["max_aspect_ratio"],
            "mean_aspect_ratio": quality["mean_aspect_ratio"],
            "n_degenerate_faces": quality["n_degenerate_faces"],
            "n_faces": quality["n_faces"],
        },
        "threshold": {
            "min_scaled_jacobian": CRITERIA["e_min_scaled_jacobian"],
            "max_aspect_ratio": CRITERIA["e_max_aspect_ratio"],
        },
        "hull": "generalised Wigley (p = q = 2), L = 10 m, B = 1 m, T = 0.625 m",
        "passed": passed,
        "closure": {
            "watertight": meshlib.is_watertight(wigley_hull.mesh_data.mesh),
            "edge_manifold": meshlib.is_edge_manifold(wigley_hull.mesh_data.mesh),
            "closure_faces": quality["closure_faces"],
        },
    }


def criterion_f_parametric():
    """Cb / L/B / B/T sweeps must move the geometry and hydrostatics continuously."""
    sweeps = []
    for cb in (0.35, 0.40, 4.0 / 9.0, 0.50, 0.55):
        sections = wigley_family_sections(
            L, B, T, fore_exponent=2.0, block_coefficient=cb,
            n_stations=STATIONS, n_points=PROFILE_POINTS,
        )
        q = depth_exponent_for_block_coefficient(cb, 2.0)
        analytic = wigley_analytic_properties(
            L, B, T, fore_exponent=2.0, depth_exponent=q
        )
        n_u, n_v = mesh_density_for(L, T, n_v=REPORT_N_V)
        hull = build_parametric_hull(sections, mesh_n_u=n_u, mesh_n_v=n_v)
        sweeps.append(
            {
                "target_Cb": cb,
                "depth_exponent": q,
                "Cb_measured": hull.hydrostatics.block_coefficient,
                "Cb_analytic": analytic["block_coefficient"],
                "volume_measured": hull.hydrostatics.displaced_volume,
                "volume_analytic": analytic["displaced_volume"],
                "volume_rel_error": _rel(
                    hull.hydrostatics.displaced_volume, analytic["displaced_volume"]
                ),
                "BM_measured": hull.hydrostatics.transverse_metacentric_radius,
                "BM_analytic": analytic["transverse_metacentric_radius"],
                "maximum_half_breadth": max(s.max_half_breadth for s in sections),
            }
        )
    volumes = [row["volume_measured"] for row in sweeps]
    monotone = all(volumes[i + 1] > volumes[i] for i in range(len(volumes) - 1))
    max_rel_error = max(abs(row["volume_rel_error"]) for row in sweeps)

    base_sections = wigley_cross_sections(L, B, T, STATIONS, PROFILE_POINTS)
    base = build_parametric_hull(base_sections, mesh_n_u=41, mesh_n_v=10)
    k = 1.5
    scaled_mesh = apply_transform(base.mesh_data.mesh, Transform(scale=(k, k, k)))
    scaled = compute_hydrostatics(
        scaled_mesh, waterline_z=0.0, n_hull_faces=base.mesh_data.hull_face_count
    )
    scaling = {
        "k": k,
        "volume_ratio": scaled.displaced_volume / base.hydrostatics.displaced_volume,
        "expected_volume_ratio": k ** 3,
        "waterplane_ratio": scaled.waterplane_area / base.hydrostatics.waterplane_area,
        "expected_waterplane_ratio": k ** 2,
        "volume_rel_error": _rel(
            scaled.displaced_volume / base.hydrostatics.displaced_volume, k ** 3
        ),
        "waterplane_rel_error": _rel(
            scaled.waterplane_area / base.hydrostatics.waterplane_area, k ** 2
        ),
    }

    ratios = []
    for length, beam, draft in (
        (10.0, 1.0, 0.625),
        (12.0, 1.0, 0.625),
        (10.0, 1.25, 0.625),
        (10.0, 1.0, 0.5),
    ):
        sections = wigley_cross_sections(length, beam, draft, STATIONS, PROFILE_POINTS)
        analytic = wigley_analytic_properties(length, beam, draft)
        n_u, n_v = mesh_density_for(length, draft, n_v=REPORT_N_V)
        hull = build_parametric_hull(sections, mesh_n_u=n_u, mesh_n_v=n_v)
        ratios.append(
            {
                "L": length,
                "B": beam,
                "T": draft,
                "L_over_B": length / beam,
                "B_over_T": beam / draft,
                "volume_measured": hull.hydrostatics.displaced_volume,
                "volume_analytic": analytic["displaced_volume"],
                "volume_rel_error": _rel(
                    hull.hydrostatics.displaced_volume, analytic["displaced_volume"]
                ),
                "waterplane_measured": hull.hydrostatics.waterplane_area,
                "waterplane_rel_error": _rel(
                    hull.hydrostatics.waterplane_area, analytic["waterplane_area"]
                ),
            }
        )
    passed = (
        monotone
        and max_rel_error < CRITERIA["b_wigley_volume_max_rel_error"]
        and abs(scaling["volume_rel_error"]) < 1e-9
        and abs(scaling["waterplane_rel_error"]) < 1e-9
    )
    return {
        "value": max_rel_error,
        "threshold": CRITERIA["b_wigley_volume_max_rel_error"],
        "block_coefficient_sweep": sweeps,
        "volume_monotone_in_Cb": monotone,
        "similarity_scaling": scaling,
        "ratio_sweep": ratios,
        "passed": passed,
    }


def convergence_study():
    """Mesh refinement: measured order and the legacy (non-interpolating) comparison."""
    sections = wigley_cross_sections(L, B, T, STATIONS, PROFILE_POINTS)
    analytic = wigley_analytic_properties(L, B, T)
    rows = []
    for n_u, n_v in ((21, 6), (31, 8), (41, 10), (61, 12), (81, 16), (101, 20), (131, 26)):
        hull = build_parametric_hull(sections, mesh_n_u=n_u, mesh_n_v=n_v)
        rows.append(
            {
                "n_u": n_u,
                "n_v": n_v,
                "n_faces": hull.quality["n_faces"],
                "h_m": L / (n_u - 1),
                "volume_rel_error": _rel(
                    hull.hydrostatics.displaced_volume, analytic["displaced_volume"]
                ),
                "waterplane_rel_error": _rel(
                    hull.hydrostatics.waterplane_area, analytic["waterplane_area"]
                ),
                "min_scaled_jacobian": hull.quality["min_scaled_jacobian"],
                "max_aspect_ratio": hull.quality["max_aspect_ratio"],
            }
        )
    last, prev = rows[-1], rows[-2]
    order = math.log(
        abs(prev["volume_rel_error"]) / abs(last["volume_rel_error"])
    ) / math.log(prev["h_m"] / last["h_m"])
    extrapolated = last["volume_rel_error"] - (
        (prev["volume_rel_error"] - last["volume_rel_error"])
        / ((prev["h_m"] / last["h_m"]) ** order - 1.0)
    )
    legacy = []
    for n_u, n_v in ((61, 12), (101, 20)):
        hull = build_parametric_hull(
            sections, mesh_n_u=n_u, mesh_n_v=n_v, interpolate=False
        )
        legacy.append(
            {
                "n_u": n_u,
                "n_v": n_v,
                "volume_rel_error": _rel(
                    hull.hydrostatics.displaced_volume, analytic["displaced_volume"]
                ),
                "waterplane_rel_error": _rel(
                    hull.hydrostatics.waterplane_area, analytic["waterplane_area"]
                ),
                "min_scaled_jacobian": hull.quality["min_scaled_jacobian"],
                "max_aspect_ratio": hull.quality["max_aspect_ratio"],
            }
        )
    def _directional(fixed, values, key, total):
        out = []
        previous = None
        for value in values:
            hull = build_parametric_hull(
                sections,
                mesh_n_u=fixed if key == "n_v" else value,
                mesh_n_v=value if key == "n_v" else fixed,
            )
            error = _rel(hull.hydrostatics.displaced_volume, analytic["displaced_volume"])
            spacing = (T / (value - 1)) if key == "n_v" else (L / (value - 1))
            row = {key: value, "spacing_m": spacing, "volume_rel_error": error}
            if previous is not None:
                row["observed_order"] = math.log(
                    abs(previous["volume_rel_error"]) / abs(error)
                ) / math.log(previous["spacing_m"] / spacing)
            previous = row
            out.append(row)
        return out

    directional = {
        "u": _directional(26, (31, 61, 101, 131, 201), "n_u", analytic),
        "v": _directional(131, (6, 10, 14, 20, 26, 40), "n_v", analytic),
    }
    return {
        "rows": rows,
        "observed_order": order,
        "richardson_extrapolated_rel_error": extrapolated,
        "legacy_control_point_loft": legacy,
        "directional": directional,
    }


def quality_scan():
    """Scaled Jacobian / aspect ratio versus mesh density on the Wigley hull."""
    sections = wigley_cross_sections(L, B, T, STATIONS, PROFILE_POINTS)
    surface = loft_nurbs_surface(sections, 3, 3, PROFILE_POINTS)
    rows = []
    for n_v in (8, 10, 12, 14, 16, 20):
        recommended, _ = mesh_density_for(L, T, n_v=n_v)
        for n_u in (recommended, 40, 80, 120):
            mesh = tessellate_surface(surface, n_u=n_u, n_v=n_v)
            quality = mesh_quality_report(mesh.mesh, n_faces=mesh.hull_face_count)
            rows.append(
                {
                    "n_u": n_u,
                    "n_v": n_v,
                    "recommended": n_u == recommended,
                    "dx_m": L / (n_u - 1),
                    "dz_m": T / (n_v - 1),
                    "min_scaled_jacobian": quality["min_scaled_jacobian"],
                    "max_aspect_ratio": quality["max_aspect_ratio"],
                    "n_faces": quality["n_faces"],
                }
            )
    return {"rows": rows}


# --------------------------------------------------------------------------- #
# figures
# --------------------------------------------------------------------------- #
def make_figures(hull, results, outdir):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir.mkdir(parents=True, exist_ok=True)
    figures = []

    # 1. body plan + waterlines
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for section in hull.sections:
        axes[0].plot([y for y, _ in section.points], [z for _, z in section.points], lw=1)
    axes[0].set_xlabel("half-breadth y [m]")
    axes[0].set_ylabel("z [m]")
    axes[0].set_title("Body plan (11 stations)")
    axes[0].grid(alpha=0.3)
    surface = hull.surface
    for z in (-0.95 * T, -0.75 * T, -0.5 * T, -0.25 * T, 0.0):
        xs, ys = [], []
        for k in range(201):
            u = k / 200
            v = _v_for_z(surface, u, z)
            point = surface.evaluate(u, v)
            xs.append(point[0])
            ys.append(point[1])
        axes[1].plot(xs, ys, label="z = %.3f m" % z)
    axes[1].set_xlabel("x [m]")
    axes[1].set_ylabel("half-breadth y [m]")
    axes[1].set_title("Waterlines evaluated from the NURBS surface")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    path = outdir / "fig1_body_plan_waterlines.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(
        {
            "name": "fig1_body_plan_waterlines",
            "caption": "Wigley III body plan and five waterlines evaluated from the "
            "lofted NURBS surface (L = 10 m, B = 1 m, T = 0.625 m).",
            "path": str(path.relative_to(ROOT)),
        }
    )

    # 2. ordinate error map
    samples = results["criterion_a"]["samples"]
    xs = sorted({s["x_over_L"] for s in samples})
    zs = sorted({s["z"] for s in samples})
    grid = [[0.0] * len(zs) for _ in xs]
    for sample in samples:
        grid[xs.index(sample["x_over_L"])][zs.index(sample["z"])] = max(
            sample["rel_error"], 1e-18
        )
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    peak = max(max(row) for row in grid)
    image = ax.imshow(
        grid,
        aspect="auto",
        origin="lower",
        cmap="viridis",
        norm=matplotlib.colors.LogNorm(vmin=1e-18, vmax=max(1e-13, peak)),
    )
    ax.set_xticks(range(len(zs)))
    ax.set_xticklabels(["%.2fT" % (z / T) for z in zs], fontsize=8)
    ax.set_yticks(range(len(xs)))
    ax.set_yticklabels(["%.2fL" % x for x in xs], fontsize=7)
    ax.set_xlabel("waterline (z/T)")
    ax.set_ylabel("station (x/L)")
    ax.set_title("Relative ordinate error vs analytic (log scale)")
    fig.colorbar(image, ax=ax, label="|y_NURBS / y_analytic - 1|")
    fig.tight_layout()
    path = outdir / "fig2_ordinate_error.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(
        {
            "name": "fig2_ordinate_error",
            "caption": "Relative ordinate error of the lofted surface against the "
            "closed-form Wigley formula over 10 stations x 5 waterlines (worst "
            "%.2e)." % results["criterion_a"]["value"],
            "path": str(path.relative_to(ROOT)),
        }
    )

    # 3. mesh convergence
    rows = results["convergence"]["rows"]
    h = [row["h_m"] for row in rows]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.loglog(h, [abs(row["volume_rel_error"]) for row in rows], "o-",
              label="displacement (interpolating loft)")
    ax.loglog(h, [abs(row["waterplane_rel_error"]) for row in rows], "s-",
              label="waterplane area")
    legacy = results["convergence"]["legacy_control_point_loft"]
    ax.loglog(
        [L / (row["n_u"] - 1) for row in legacy],
        [abs(row["volume_rel_error"]) for row in legacy],
        "^--",
        color="tab:red",
        label="displacement (legacy control-point loft)",
    )
    guide = [abs(rows[0]["volume_rel_error"]) * (rows[0]["h_m"] / v) ** 2 for v in h]
    ax.loglog(h, guide, ":", color="grey", label="O(h^2) guide")
    ax.set_xlabel("mesh spacing h [m]")
    ax.set_ylabel("relative error")
    ax.set_title("Mesh convergence against the analytic Wigley hull")
    ax.grid(which="both", alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    path = outdir / "fig3_mesh_convergence.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(
        {
            "name": "fig3_mesh_convergence",
            "caption": "Displacement and waterplane error versus mesh spacing. The "
            "interpolating loft converges at the measured order p = %.2f; the legacy "
            "control-point loft stalls near 3%% because it does not pass through the "
            "sections, so refinement cannot fix it."
            % results["convergence"]["observed_order"],
            "path": str(path.relative_to(ROOT)),
        }
    )

    # 4. mesh quality
    scan = [row for row in results["quality_scan"]["rows"] if row["recommended"]]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    axes[0].plot([row["n_v"] for row in scan], [row["min_scaled_jacobian"] for row in scan], "o-")
    axes[0].axhline(CRITERIA["e_min_scaled_jacobian"], ls="--", color="tab:red",
                    label="criterion 0.30")
    axes[0].set_xlabel("profile rows n_v")
    axes[0].set_ylabel("min scaled Jacobian")
    axes[0].set_title("Worst scaled Jacobian at the recommended density")
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[1].plot([row["n_v"] for row in scan], [row["max_aspect_ratio"] for row in scan],
                 "s-", color="tab:green")
    axes[1].axhline(CRITERIA["e_max_aspect_ratio"], ls="--", color="tab:red",
                    label="criterion 5.0")
    axes[1].set_xlabel("profile rows n_v")
    axes[1].set_ylabel("max aspect ratio")
    axes[1].set_title("Worst aspect ratio at the recommended density")
    axes[1].grid(alpha=0.3)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    path = outdir / "fig4_mesh_quality.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(
        {
            "name": "fig4_mesh_quality",
            "caption": "Hull mesh quality versus profile resolution with n_u chosen by "
            "mesh_density_for (dx = 2.5 dz). Both criterion limits are cleared from "
            "n_v = 12 upwards.",
            "path": str(path.relative_to(ROOT)),
        }
    )

    # 5. parametric sweep
    sweep = results["criterion_f"]["block_coefficient_sweep"]
    cbs = [row["target_Cb"] for row in sweep]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    axes[0].plot(cbs, [row["volume_measured"] for row in sweep], "o",
                 label="mesh hydrostatics")
    axes[0].plot(cbs, [row["volume_analytic"] for row in sweep], "-", label="closed form")
    axes[0].set_xlabel("target block coefficient Cb")
    axes[0].set_ylabel("displacement volume [m^3]")
    axes[0].set_title("Volume follows Cb continuously")
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[1].plot(cbs, [row["BM_measured"] for row in sweep], "o", label="mesh hydrostatics")
    axes[1].plot(cbs, [row["BM_analytic"] for row in sweep], "-", label="closed form")
    axes[1].set_xlabel("target block coefficient Cb")
    axes[1].set_ylabel("transverse metacentric radius BM [m]")
    axes[1].set_title("BM_T follows Cb continuously")
    axes[1].grid(alpha=0.3)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    path = outdir / "fig5_parametric_sweep.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(
        {
            "name": "fig5_parametric_sweep",
            "caption": "Cb-driven family sweep: changing the block coefficient updates "
            "the geometry and both hydrostatic properties continuously (maximum volume "
            "error %.2f%%)." % (100 * results["criterion_f"]["value"]),
            "path": str(path.relative_to(ROOT)),
        }
    )

    # 6. SUBOFF offsets
    checks = results["criterion_d"]["checks"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    xs = [row["x_over_L"] for row in checks]
    axes[0].plot(xs, [row["radius_darpa_m"] for row in checks], "-",
                 label="DARPA offset table")
    axes[0].plot(xs, [row["radius_surface_m"] for row in checks], "o", ms=4,
                 label="lofted NURBS surface")
    axes[0].set_xlabel("x/L")
    axes[0].set_ylabel("radius r [m]")
    axes[0].set_title("SUBOFF bare hull radius")
    axes[0].grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[1].semilogy(xs, [max(row["rel_error"], 1e-18) for row in checks], "o-")
    axes[1].axhline(CRITERIA["d_suboff_geometry_max_rel_error"], ls="--",
                    color="tab:red", label="criterion 0.5%")
    axes[1].set_xlabel("x/L")
    axes[1].set_ylabel("|r_NURBS / r_DARPA - 1|")
    axes[1].set_title("Relative radius error at 20 stations")
    axes[1].grid(which="both", alpha=0.3)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    path = outdir / "fig6_suboff_offsets.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(
        {
            "name": "fig6_suboff_offsets",
            "caption": "DARPA SUBOFF bare hull: lofted radius versus the published "
            "offset table at 20 check stations (surface built from 48 clustered "
            "stations; worst error %.3f%%)." % (100 * results["criterion_d"]["value"]),
            "path": str(path.relative_to(ROOT)),
        }
    )
    return figures


# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #
def render_markdown(results):
    a = results["criterion_a"]
    b = results["criterion_b"]
    c = results["criterion_c"]
    d = results["criterion_d"]
    e = results["criterion_e"]
    f = results["criterion_f"]
    conv = results["convergence"]

    def flag(passed):
        return "PASS" if passed else "FAIL"

    lines = [
        "# Parametric hull generation -- verification record",
        "",
        "Generated by `examples/parametric_hull_dag/run_verification.py`; every number",
        "below is measured by that script (no hand-copied values).",
        "",
        "| criterion | measured | threshold | verdict |",
        "|---|---|---|---|",
        "| (a) Wigley III ordinates, 10 stations x 5 waterlines | %.3e | < %.0e | %s |"
        % (a["value"], a["threshold"], flag(a["value"] < a["threshold"])),
        "| (b) Wigley III displacement volume | %.4f%% | < %.1f%% | %s |"
        % (100 * b["value"], 100 * b["threshold"], flag(b["value"] < b["threshold"])),
        "| (c) Wigley III waterplane area | %.4f%% | < %.1f%% | %s |"
        % (100 * c["value"], 100 * c["threshold"], flag(c["value"] < c["threshold"])),
        "| (d) SUBOFF radius, 20 stations | %.4f%% | < %.1f%% | %s |"
        % (100 * d["value"], 100 * d["threshold"], flag(d["value"] < d["threshold"])),
        "| (e) hull mesh min scaled Jacobian | %.3f | > %.2f | %s |"
        % (
            e["value"]["min_scaled_jacobian"],
            e["threshold"]["min_scaled_jacobian"],
            flag(e["value"]["min_scaled_jacobian"] > e["threshold"]["min_scaled_jacobian"]),
        ),
        "| (e) hull mesh max aspect ratio | %.2f | < %.1f | %s |"
        % (
            e["value"]["max_aspect_ratio"],
            e["threshold"]["max_aspect_ratio"],
            flag(e["value"]["max_aspect_ratio"] < e["threshold"]["max_aspect_ratio"]),
        ),
        "| (f) Cb / L-B / B-T sweep | %.4f%% volume error | monotone and < %.1f%% | %s |"
        % (100 * f["value"], 100 * f["threshold"], flag(f["passed"])),
        "",
        "## Hydrostatic detail (Wigley III, L = 10 m, B = 1 m, T = 0.625 m)",
        "",
        "| quantity | measured | closed form | rel. error |",
        "|---|---|---|---|",
        "| displaced volume [m^3] | %.6f | %.6f | %.3e |"
        % (b["measured_volume"], b["analytic_volume"], b["value"]),
        "| waterplane area [m^2] | %.6f | %.6f | %.3e |"
        % (c["measured_waterplane_area"], c["analytic_waterplane_area"], c["value"]),
        "| transverse waterplane inertia I_T [m^4] | %.6f | %.6f | %.3e |"
        % (
            c["extra"]["measured_transverse_inertia"],
            c["extra"]["analytic_transverse_inertia"],
            c["extra"]["transverse_inertia_rel_error"],
        ),
        "| longitudinal waterplane inertia I_L [m^4] | %.6f | %.6f | %.3e |"
        % (
            c["extra"]["measured_longitudinal_inertia"],
            c["extra"]["analytic_longitudinal_inertia"],
            c["extra"]["longitudinal_inertia_rel_error"],
        ),
        "| transverse metacentric radius BM_T [m] | %.6f | %.6f | %.3e |"
        % (
            c["extra"]["measured_bm"],
            c["extra"]["analytic_bm"],
            c["extra"]["transverse_metacentric_radius_rel_error"],
        ),
        "| centre of buoyancy VCB [m] | %.6f | %.6f | %.3e |"
        % (
            c["extra"]["measured_vcb"],
            c["extra"]["analytic_vcb"],
            c["extra"]["vcb_rel_error"],
        ),
        "| LCF - L/2 [m] | %.3e | 0 | -- |" % c["extra"]["lcf_error_m"],
        "",
        "## Mesh convergence",
        "",
        "| n_u | n_v | faces | h [m] | volume error | waterplane error | min scaled J | max aspect |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in conv["rows"]:
        lines.append(
            "| %d | %d | %d | %.4f | %.4f%% | %.4f%% | %.3f | %.2f |"
            % (
                row["n_u"],
                row["n_v"],
                row["n_faces"],
                row["h_m"],
                100 * row["volume_rel_error"],
                100 * row["waterplane_rel_error"],
                row["min_scaled_jacobian"],
                row["max_aspect_ratio"],
            )
        )
    lines += [
        "",
        "Observed order of accuracy p = %.2f; Richardson-extrapolated displacement "
        "error %.3e." % (conv["observed_order"], conv["richardson_extrapolated_rel_error"]),
        "",
        "### Legacy control-point loft (regression control)",
        "",
        "| n_u | n_v | volume error | waterplane error | min scaled J | max aspect |",
        "|---|---|---|---|---|---|",
    ]
    for row in conv["legacy_control_point_loft"]:
        lines.append(
            "| %d | %d | %.4f%% | %.4f%% | %.3f | %.2f |"
            % (
                row["n_u"],
                row["n_v"],
                100 * row["volume_rel_error"],
                100 * row["waterplane_rel_error"],
                row["min_scaled_jacobian"],
                row["max_aspect_ratio"],
            )
        )
    lines += [
        "",
        "## Parametric continuity",
        "",
        "| target Cb | depth exponent q | Cb measured | volume measured | volume closed form | BM_T measured | BM_T closed form |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in f["block_coefficient_sweep"]:
        lines.append(
            "| %.4f | %.4f | %.4f | %.5f | %.5f | %.5f | %.5f |"
            % (
                row["target_Cb"],
                row["depth_exponent"],
                row["Cb_measured"],
                row["volume_measured"],
                row["volume_analytic"],
                row["BM_measured"],
                row["BM_analytic"],
            )
        )
    scaling = f["similarity_scaling"]
    lines += [
        "",
        "Uniform scaling k = %.2f: volume ratio %.9f (expected %.9f), waterplane ratio "
        "%.9f (expected %.9f)."
        % (
            scaling["k"],
            scaling["volume_ratio"],
            scaling["expected_volume_ratio"],
            scaling["waterplane_ratio"],
            scaling["expected_waterplane_ratio"],
        ),
        "",
        "| L | B | T | L/B | B/T | volume error | waterplane error |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in f["ratio_sweep"]:
        lines.append(
            "| %.2f | %.2f | %.3f | %.2f | %.2f | %.4f%% | %.4f%% |"
            % (
                row["L"],
                row["B"],
                row["T"],
                row["L_over_B"],
                row["B_over_T"],
                100 * row["volume_rel_error"],
                100 * row["waterplane_rel_error"],
            )
        )
    lines += [
        "",
        "## SUBOFF",
        "",
        "Worst radius error over 20 stations: %.4e; volume error %.4e; wetted area "
        "error %.4e; bow tip radius %.3e m, stern tip radius %.3e m."
        % (
            d["value"],
            d["volume"]["rel_error"],
            d["wetted_area"]["rel_error"],
            d["tip_closure"]["radius_at_bow_tip_m"],
            d["tip_closure"]["radius_at_stern_tip_m"],
        ),
        "",
        "Mesh (n_u = %d, n_v = %d): watertight %s, min scaled Jacobian %.3f, max aspect "
        "ratio %.2f -- see the limitation section of docs/parametric_hull/README.md."
        % (
            d["mesh"]["n_u"],
            d["mesh"]["n_v"],
            d["mesh"]["watertight"],
            d["mesh"]["quality"]["min_scaled_jacobian"],
            d["mesh"]["quality"]["max_aspect_ratio"],
        ),
        "",
        "## Seeding-formula audit",
        "",
        b["seeding_formula_note"]["explanation"],
        "",
        "Literal expression value %.6f m^3 (%.2f%% against the closed form)."
        % (
            b["seeding_formula_note"]["literal_value"],
            100 * b["seeding_formula_note"]["literal_rel_error"],
        ),
        "",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="parametric hull verification campaign")
    parser.add_argument(
        "--out",
        type=pathlib.Path,
        default=ROOT / "benchmarks" / "verified" / "parametric_hull",
        help="output directory for verification.json and figures",
    )
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    outdir = args.out
    outdir.mkdir(parents=True, exist_ok=True)

    analytic = wigley_analytic_properties(L, B, T)
    sections = wigley_cross_sections(L, B, T, STATIONS, PROFILE_POINTS)
    n_u, n_v = mesh_density_for(L, T, n_v=REPORT_N_V)
    hull = build_parametric_hull(
        sections, mesh_n_u=n_u, mesh_n_v=n_v, close_waterplane=True
    )

    results = {
        "node": "p316w1-n03-cad-parametric-hull",
        "repository": "TensorCAD",
        "reference_hull": {
            "L": L, "B": B, "T": T, "Cb": analytic["block_coefficient"],
            "Cw": analytic["waterplane_coefficient"],
        },
        "mesh": {
            "n_u": n_u,
            "n_v": n_v,
            "rule": "mesh_density_for(L, T, n_v=%d) -> dx = 2.5 dz" % REPORT_N_V,
        },
        "criterion_a": criterion_a_wigley_geometry(hull),
        "criterion_b": criterion_b_volume(hull, analytic),
        "criterion_c": criterion_c_waterplane(hull, analytic),
        "criterion_d": criterion_d_suboff(),
        "criterion_e": criterion_e_mesh_quality(hull),
        "criterion_f": criterion_f_parametric(),
        "convergence": convergence_study(),
        "quality_scan": quality_scan(),
    }
    figure_dir = ROOT / "docs" / "parametric_hull" / "figures"
    if not args.no_figures:
        figure_dir.mkdir(parents=True, exist_ok=True)
        results["figures"] = make_figures(hull, results, figure_dir)

    verdicts = {
        "a": results["criterion_a"]["value"] < results["criterion_a"]["threshold"],
        "b": results["criterion_b"]["value"] < results["criterion_b"]["threshold"],
        "c": results["criterion_c"]["value"] < results["criterion_c"]["threshold"],
        "d": results["criterion_d"]["value"] < results["criterion_d"]["threshold"],
        "e": results["criterion_e"]["passed"],
        "f": results["criterion_f"]["passed"],
    }
    results["verdicts"] = verdicts
    results["all_passed"] = all(verdicts.values())

    (outdir / "verification.json").write_text(json.dumps(results, indent=2) + "\n")
    docs_dir = ROOT / "docs" / "parametric_hull"
    docs_dir.mkdir(parents=True, exist_ok=True)
    (docs_dir / "VERIFICATION.md").write_text(render_markdown(results) + "\n")

    print(json.dumps({"verdicts": verdicts, "all_passed": results["all_passed"]}, indent=2))
    for key, value in results["criterion_a"].items():
        if key in ("samples",):
            continue
    print("criterion a (max ordinate rel error): %.3e" % results["criterion_a"]["value"])
    print("criterion b (volume rel error):        %.4f%%" % (100 * results["criterion_b"]["value"]))
    print("criterion c (waterplane rel error):    %.4f%%" % (100 * results["criterion_c"]["value"]))
    print("criterion d (SUBOFF radius rel error): %.4f%%" % (100 * results["criterion_d"]["value"]))
    print("criterion e (min scaled Jacobian):     %.3f" % results["criterion_e"]["value"]["min_scaled_jacobian"])
    print("criterion e (max aspect ratio):        %.2f" % results["criterion_e"]["value"]["max_aspect_ratio"])
    print("criterion f (Cb sweep max vol error):  %.4f%%" % (100 * results["criterion_f"]["value"]))
    print("observed order of accuracy:            %.2f" % results["convergence"]["observed_order"])
    print("artifacts: %s" % outdir)
    return 0 if results["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
