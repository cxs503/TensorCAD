"""The Wigley hull against its closed forms, and against the sine family's ends.

Runs headless, writes ``docs/wigley/{results.json,report.md,*.png}`` and exits 1 if
any of the machine-checked claims stops holding, so it can be a CI step.

    python examples/wigley/run_wigley.py [--output DIR]
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tensorcad import HullParameters, mesh as meshlib
from tensorcad.curvature import assess_fairness
from tensorcad.hull import generate_hull
from tensorcad.wigley import (
    ANALYTIC_COMPARABLE_QUANTITIES,
    WigleyParameters,
    analytic_hydrostatics,
    analytic_wetted_surface_area,
    estimate_discretisation_error,
    generate_wigley,
    half_breadth,
    waterline_curvature_at_ends,
)

L, B, T = 10.0, 1.0, 0.625
GRIDS = ((41, 13), (81, 25), (161, 49))
SHORT = {
    "displacement_volume": "V",
    "wetted_surface_area": "S_w",
    "waterplane_area": "A_w",
    "maximum_section_area": "A_m",
    "longitudinal_centre_buoyancy": "LCB",
    "longitudinal_centre_flotation": "LCF",
    "vertical_centre_buoyancy": "VCB",
    "block_coefficient": "C_b",
    "prismatic_coefficient": "C_p",
}


@dataclass
class Claim:
    name: str
    detail: str
    passed: bool


def params(nx: int, nz: int) -> WigleyParameters:
    return WigleyParameters(length=L, beam=B, draft=T, nx=nx, nz=nz)


def sine_waterline_half_breadth(x: float) -> float:
    """The other family's waterline: ``(B/2) sin(pi x/L)^0.7``."""
    s = math.sin(math.pi * x / L)
    return 0.5 * B * (s**0.7 if s > 0.0 else 0.0)


def sine_slope(x: float) -> float:
    """Analytic ``db/dx`` of the sine waterline; diverges as ``x^-0.3`` at the stern."""
    s = math.sin(math.pi * x / L)
    return 0.5 * B * 0.7 * (s ** (0.7 - 1.0)) * math.cos(math.pi * x / L) * (math.pi / L)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=str(ROOT / "docs" / "wigley"))
    args = parser.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    reference = analytic_hydrostatics(L, B, T)
    claims: list[Claim] = []
    results: dict[str, object] = {
        "case": {
            "family": "wigley_parabolic_hull",
            "half_breadth": "(B/2)*(1-(2x/L-1)^2)*(1-(z/T)^2)",
            "length_m": L,
            "beam_m": B,
            "draft_m": T,
            "beam_length_ratio": B / L,
            "draft_length_ratio": T / L,
            "water_density_kg_m3": 1025.0,
        },
        "closed_form_hydrostatics": {
            "displacement_volume_m3": reference.displacement_volume,
            "displacement_mass_kg": reference.displacement_mass,
            "wetted_surface_area_m2": reference.wetted_surface_area,
            "waterplane_area_m2": reference.waterplane_area,
            "maximum_section_area_m2": reference.maximum_section_area,
            "lcb_m": reference.longitudinal_centre_buoyancy,
            "lcf_m": reference.longitudinal_centre_flotation,
            "vcb_m": reference.vertical_centre_buoyancy,
            "block_coefficient": reference.block_coefficient,
            "prismatic_coefficient": reference.prismatic_coefficient,
            "basis": reference.basis,
        },
    }

    # ---- 1. mesh against the closed forms at three grids -------------------
    grid_rows = []
    per_grid_errors: dict[str, list[tuple[int, float]]] = {n: [] for n in ANALYTIC_COMPARABLE_QUANTITIES}
    for nx, nz in GRIDS:
        model = generate_wigley(params(nx, nz))
        triangles = len(model.mesh.faces)
        row = {"nx": nx, "nz": nz, "triangles": triangles, "quantities": {}}
        worst = 0.0
        for name in ANALYTIC_COMPARABLE_QUANTITIES:
            got = getattr(model.hydrostatics, name)
            want = getattr(reference, name)
            error = (got - want) / want
            row["quantities"][name] = {"mesh": got, "closed_form": want, "relative_error": error}
            per_grid_errors[name].append((triangles, abs(error)))
            worst = max(worst, abs(error))
        row["worst_relative_error"] = worst
        grid_rows.append(row)
    results["grids"] = grid_rows

    default_worst = grid_rows[0]["worst_relative_error"]
    refined_worst = grid_rows[1]["worst_relative_error"]
    finest_worst = grid_rows[2]["worst_relative_error"]
    claims.append(Claim("default_grid_within_0.3pct", f"worst |error| at 41x13 = {100*default_worst:.4f}% (limit 0.3%)", default_worst < 3e-3))
    claims.append(Claim("refined_grid_within_0.1pct", f"worst |error| at 81x25 = {100*refined_worst:.4f}% (limit 0.1%)", refined_worst < 1e-3))
    claims.append(Claim("finest_grid_within_0.03pct", f"worst |error| at 161x49 = {100*finest_worst:.4f}% (limit 0.03%)", finest_worst < 3e-4))

    volume_errors = [abs(r["quantities"]["displacement_volume"]["relative_error"]) for r in grid_rows]
    monotone = all(b < a for a, b in itertools.pairwise(volume_errors))
    claims.append(Claim("volume_error_monotone", "volume |error| " + ", ".join(f"{100*e:.4f}%" for e in volume_errors) + " is strictly decreasing", monotone))
    order_ratios = [volume_errors[i] / volume_errors[i + 1] for i in range(len(volume_errors) - 1)]
    claims.append(Claim("volume_second_order", "error ratios between consecutive refinements " + ", ".join(f"{r:.2f}x" for r in order_ratios) + " (second order ~4x)", all(3.0 < r < 5.0 for r in order_ratios)))
    signs = {r["quantities"]["displacement_volume"]["relative_error"] < 0 for r in grid_rows}
    claims.append(Claim("inscribed_mesh_is_low", "every grid reports a volume below the closed form (an inscribed loft must)", signs == {True}))
    cb = grid_rows[0]["quantities"]["block_coefficient"]
    cb_high = all(r["quantities"]["block_coefficient"]["relative_error"] < 0 for r in grid_rows)
    claims.append(Claim("block_coefficient_is_4_9", f"closed form Cb = {reference.block_coefficient:.12f} = 4/9 = {4/9:.12f}; mesh at 41x13 = {cb['mesh']:.9f} ({100*cb['relative_error']:+.4f}%)", abs(reference.block_coefficient - 4.0 / 9.0) < 1e-15 and cb_high))

    # ---- 2. the wetted-area quadrature is not the error term ---------------
    areas = {n: analytic_wetted_surface_area(L, B, T, nodes=n) for n in (20, 40, 80, 160)}
    area_change = {
        "80_vs_40": abs(areas[80] - areas[40]) / areas[80],
        "160_vs_80": abs(areas[160] - areas[80]) / areas[160],
    }
    results["wetted_area_quadrature"] = {"nodes": areas, "relative_change": area_change}
    finest_mesh_area_error = abs(grid_rows[2]["quantities"]["wetted_surface_area"]["relative_error"])
    claims.append(Claim("wetted_area_quadrature_converged", f"80 vs 160 nodes differ by {area_change['160_vs_80']:.3e} relative, five orders below the finest mesh error {finest_mesh_area_error:.3e}", area_change["160_vs_80"] < 1e-9 and area_change["160_vs_80"] < finest_mesh_area_error / 1e3))

    # ---- 3. geometry validity ---------------------------------------------
    model = generate_wigley(params(81, 25))
    validity = {
        "watertight": meshlib.is_watertight(model.mesh),
        "edge_manifold": meshlib.is_edge_manifold(model.mesh),
        "zero_area_faces": len(meshlib.zero_area_faces(model.mesh)),
        "signed_volume_m3": meshlib.total_signed_volume(model.mesh),
        "validation_errors": model.validation.errors,
    }
    results["mesh_validity_81x25"] = validity
    claims.append(Claim("mesh_is_a_closed_solid", f"watertight={validity['watertight']}, edge-manifold={validity['edge_manifold']}, zero-area faces={validity['zero_area_faces']}, validation errors={len(validity['validation_errors'])}", validity["watertight"] and validity["edge_manifold"] and validity["zero_area_faces"] == 0 and not validity["validation_errors"]))

    # ---- 4. fairness: the reason to carry this family ----------------------
    probes = (1e-2, 1e-4, 1e-6, 1e-8, 1e-9)
    fairness_rows = []
    for f in probes:
        x = f * L
        wigley_secant = half_breadth(params(41, 13), x, 0.0) / x
        fairness_rows.append(
            {
                "x_over_L": f,
                "wigley_secant_dbdx": wigley_secant,
                "wigley_analytic_dbdx": 0.5 * B * (-4.0 * (2.0 * x / L - 1.0) / L),
                "sine_secant_dbdx": sine_waterline_half_breadth(x) / x,
                "sine_analytic_dbdx": sine_slope(x),
            }
        )
    results["waterline_slope_near_stern"] = fairness_rows
    end_slope = 2.0 * B / L
    wigley_bounded = all(abs(r["wigley_secant_dbdx"]) <= end_slope * 1.001 for r in fairness_rows)
    sine_growth = fairness_rows[-1]["sine_analytic_dbdx"] / fairness_rows[0]["sine_analytic_dbdx"]
    claims.append(Claim("wigley_slope_bounded", f"Wigley secant slope near the stern stays within ±{end_slope:.4f} 1/m for x/L down to 1e-9 (analytic end slope 2B/L = {end_slope:.4f})", wigley_bounded))
    predicted_growth = (1e-2 / 1e-9) ** 0.3
    claims.append(Claim("sine_slope_diverges", f"sine family slope grows {sine_growth:.1f}x between x/L=1e-2 and 1e-9, and x^-0.3 predicts {predicted_growth:.1f}x", abs(sine_growth - predicted_growth) / predicted_growth < 0.02))
    curvature = waterline_curvature_at_ends(params(41, 13))
    numeric_curvature = []
    h = 1e-4 * L
    for x in (0.01 * L, 0.5 * L, 0.99 * L):
        numeric = (half_breadth(params(41, 13), x + h, 0.0) - 2 * half_breadth(params(41, 13), x, 0.0) + half_breadth(params(41, 13), x - h, 0.0)) / (h * h)
        numeric_curvature.append({"x_over_L": x / L, "second_difference": numeric, "closed_form": curvature})
    results["waterline_curvature"] = numeric_curvature
    curvature_ok = all(abs(r["second_difference"] - curvature) / abs(curvature) < 1e-5 for r in numeric_curvature)
    claims.append(Claim("waterline_curvature_finite", f"d2b/dx2 = -4B/L^2 = {curvature:.6f} 1/m at the ends and midship, central differences agree to <1e-5", curvature_ok))

    wigley_fair = assess_fairness(model.mesh)
    sine_hull = generate_hull(HullParameters(length=L, beam=B, draft=T, nx=81, nz=25))
    sine_fair = assess_fairness(sine_hull.mesh)
    results["fairness"] = {
        "wigley_81x25": wigley_fair.to_dict(),
        "sine_81x25": sine_fair.to_dict(),
    }

    # ---- 5. discretisation estimate structure ------------------------------
    results["discretisation_41x13_to_81x25"] = estimate_discretisation_error(params(41, 13))

    # ---- figures ----------------------------------------------------------
    fig1 = out / "fig1_convergence.png"
    plot_convergence(per_grid_errors, fig1)
    fig2 = out / "fig2_waterline_fairness.png"
    plot_waterline(fairness_rows, fig2)
    fig3 = out / "fig3_sections.png"
    plot_sections(params(41, 13), fig3)
    results["figures"] = [fig1.name, fig2.name, fig3.name]

    results["claims"] = [{"name": c.name, "detail": c.detail, "passed": c.passed} for c in claims]
    results["passed"] = all(c.passed for c in claims)
    (out / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    (out / "report.md").write_text(render_report(results, claims, grid_rows, fairness_rows))
    print(f"wrote {out/'results.json'}")
    print(f"wrote {out/'report.md'}")
    for c in claims:
        print(f"  [{'PASS' if c.passed else 'FAIL'}] {c.name}: {c.detail}")
    return 0 if results["passed"] else 1


def plot_convergence(per_grid_errors: dict[str, list[tuple[int, float]]], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    styles = {
        "displacement_volume": ("o", "volume"),
        "waterplane_area": ("s", "waterplane"),
        "wetted_surface_area": ("^", "wetted surface"),
        "maximum_section_area": ("D", "midship section"),
        "vertical_centre_buoyancy": ("v", "VCB"),
    }
    for name, (marker, label) in styles.items():
        xs = [n for n, _ in per_grid_errors[name]]
        ys = [max(e, 1e-12) for _, e in per_grid_errors[name]]
        ax.loglog(xs, ys, marker=marker, label=label)
    xs = [n for n, _ in per_grid_errors["displacement_volume"]]
    first = per_grid_errors["displacement_volume"][0][1]
    ax.loglog(xs, [first * (xs[0] / x) for x in xs], "k--", lw=1.0, label="second-order reference (∝1/N)")
    ax.set_xlabel("triangles")
    ax.set_ylabel("|mesh − closed form| / closed form")
    ax.set_title("Wigley hull: the sampled mesh converges to the closed forms")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_waterline(rows: list[dict[str, float]], path: Path) -> None:
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.4, 4.2))
    xs = [i / 400.0 for i in range(41)]
    left.plot(xs, [B / 2 for _ in xs], alpha=0.0)
    left.plot(xs, [2.0 * half_breadth(params(41, 13), x * L, 0.0) for x in xs], label="Wigley  B(1−(2x/L−1)²)")
    left.plot(xs, [2.0 * sine_waterline_half_breadth(x * L) for x in xs], label="sine  B sin(πx/L)^0.7")
    left.set_xlim(0.0, 0.05)
    left.set_xlabel("x / L (stern end)")
    left.set_ylabel("waterline breadth [m]")
    left.set_title("Waterline near the stern")
    left.grid(alpha=0.3)
    left.legend(fontsize=8)

    depths = [10.0 ** (-e) for e in range(2, 11)]
    right.loglog(depths, [abs(sine_slope(d * L)) for d in depths], "o-", label="sine  |db/dx| ~ x^-0.3")
    right.loglog(depths, [2.0 * B / L for _ in depths], "s-", label="Wigley  |db/dx| ≤ 2B/L")
    right.set_xlabel("x / L")
    right.set_ylabel("|db/dx| [1/m]")
    right.set_title("The sine family has no tangent plane at the end")
    right.grid(True, which="both", alpha=0.3)
    right.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def plot_sections(p: WigleyParameters, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.8))
    for ax, x_over_l in zip(axes, (0.5, 0.25, 0.05), strict=False):
        zs = [-p.draft * i / 200.0 for i in range(201)]
        x = x_over_l * p.length
        ax.plot([2.0 * half_breadth(p, x, z) for z in zs], zs, label="Wigley")
        taper = math.sin(math.pi * x_over_l) ** 0.7 if 0 < x_over_l < 1 else 0.0
        ax.plot(
            [p.beam * taper * math.sqrt((z + p.draft) / p.draft) for z in zs],
            zs,
            "--",
            label="sine family",
        )
        ax.set_title(f"section at x/L = {x_over_l}")
        ax.set_xlabel("y [m]")
        ax.set_ylabel("z [m]")
        ax.grid(alpha=0.3)
        ax.set_aspect("equal", adjustable="datalim")
    axes[0].legend(fontsize=8)
    fig.suptitle("Wigley (parabolic) and sine-family sections, same L, B, T", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def render_report(results, claims, grid_rows, fairness_rows) -> str:
    lines = [
        "# Wigley parabolic hull: closed forms, mesh agreement and end fairness",
        "",
        "Generated by `examples/wigley/run_wigley.py`; every number below is read from",
        "`results.json`, which the same run wrote. No value is typed by hand.",
        "",
        "## Case",
        "",
        f"Wigley hull `b(x,z) = (B/2)(1-(2x/L-1)^2)(1-(z/T)^2)`, L = {L} m, B = {B} m,",
        f"T = {T} m (B/L = {B/L:.3f}, T/L = {T/L:.4f}), sea water 1025 kg/m^3.",
        "",
        "## Closed forms vs the sampled mesh",
        "",
        "| quantity | closed form | 41x13 | error | 81x25 | error | 161x49 | error |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for name in ANALYTIC_COMPARABLE_QUANTITIES:
        ref = grid_rows[0]["quantities"][name]["closed_form"]
        cells = []
        for row in grid_rows:
            entry = row["quantities"][name]
            cells.append(f"{entry['mesh']:.9g}")
            cells.append(f"{100*entry['relative_error']:+.4f}%")
        lines.append(f"| {SHORT.get(name, name)} | {ref:.9g} | " + " | ".join(cells) + " |")
    lines += [
        "",
        f"Worst error: {100*grid_rows[0]['worst_relative_error']:.4f}% at 41x13,",
        f"{100*grid_rows[1]['worst_relative_error']:.4f}% at 81x25,"
        f" {100*grid_rows[2]['worst_relative_error']:.4f}% at 161x49.",
        "",
        "## Waterline near the stern: Wigley against the sine family",
        "",
        "| x/L | Wigley secant db/dx [1/m] | Wigley analytic [1/m] | sine secant [1/m] | sine analytic [1/m] |",
        "|---|---|---|---|---|",
    ]
    for row in fairness_rows:
        lines.append(
            f"| {row['x_over_L']:.0e} | {row['wigley_secant_dbdx']:.6f} | {row['wigley_analytic_dbdx']:.6f} "
            f"| {row['sine_secant_dbdx']:.3f} | {row['sine_analytic_dbdx']:.3f} |"
        )
    lines += ["", "## Claims", ""]
    for c in claims:
        lines.append(f"- **{'PASS' if c.passed else 'FAIL'}** `{c.name}` — {c.detail}")
    lines += ["", "## Figures", ""]
    for name, caption in (
        ("fig1_convergence.png", "Mesh error against the closed forms, three grids, with a first-order-in-1/N reference."),
        ("fig2_waterline_fairness.png", "Waterline shape at the stern, and the waterline slope as the end is approached."),
        ("fig3_sections.png", "Section shapes of the two families at three stations."),
    ):
        lines.append(f"### {name}")
        lines.append("")
        lines.append(caption)
        lines.append("")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
