#!/usr/bin/env python3
"""Self-check of the sampled hull against the closed forms of the analytic surface.

The generator samples a surface whose half-breadth factorises,

    y(x, z) = (B/2) * sin(pi*x/L)^0.7 * sqrt((z + T)/T),

so the volume, waterplane area, midship section area, LCB, LCF, VCB and form
coefficients of that surface all have closed forms, and the wetted area has a
converged quadrature. ``tensorcad.hull.analytic_hydrostatics`` returns them as the
same ``Hydrostatics`` record the mesh path returns, which is what makes a
one-to-one comparison possible.

What this script measures, in order:

1. **The closed forms themselves** against an independent Simpson quadrature of
   the same integrands, so a formula is never checked only against itself.
2. **The mesh integrals** against those closed forms on five grids from 11x5 to
   161x41: signed relative error of every quantity, and whether it converges.
3. **The discretisation estimator**: how the coarse/refined relative change from
   ``estimate_discretisation_error`` compares with the error the closed form
   reveals. It is the number people quote, so how optimistic it is deserves a row.
4. **The wetted area's quadrature**, which is the only quantity without a formula.

Writes ``results.json``, ``report.md`` (tables generated from the numbers, never
typed by hand) and ``convergence.png``. Exits 1 if a claim stops holding.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from tensorcad import HullParameters
from tensorcad.hull import (
    BREADTH_EXPONENT,
    TAPER_INTEGRAL,
    analytic_hydrostatics,
    analytic_relative_errors,
    analytic_volume,
    analytic_wetted_surface_area,
    estimate_discretisation_error,
)

LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8
GRIDS = ((11, 5), (21, 9), (41, 13), (81, 21), (161, 41))
DEFAULT_GRID = (41, 13)
#: Claim: at the shipped default grid every quantity is inside 3 % of its closed
#: form, and at 161x41 every one is inside 0.2 %.
DEFAULT_TOLERANCE = 0.03
FINE_TOLERANCE = 0.002
QUANTITIES = (
    "displacement_volume",
    "wetted_surface_area",
    "waterplane_area",
    "maximum_section_area",
    "longitudinal_centre_buoyancy",
    "longitudinal_centre_flotation",
    "vertical_centre_buoyancy",
    "block_coefficient",
    "prismatic_coefficient",
)


def simpson(function, lower: float, upper: float, panels: int) -> float:
    if panels % 2:
        panels += 1
    step = (upper - lower) / panels
    total = function(lower) + function(upper)
    for index in range(1, panels):
        total += (4.0 if index % 2 else 2.0) * function(lower + index * step)
    return total * step / 3.0


def taper_integral_by_simpson(panels: int = 4000) -> float:
    """Independent quadrature of int_0^1 sin(pi x)^0.7 dx (theta substitution)."""

    def integrand(theta: float) -> float:
        return 0.5 * math.sin(theta) * math.sin(0.5 * math.pi * (1.0 - math.cos(theta))) ** BREADTH_EXPONENT

    return simpson(integrand, 0.0, math.pi, panels)


def wetted_area_by_simpson(length: float, beam: float, draft: float, panels: int = 400) -> float:
    """The same surface element as the module, integrated by Simpson instead of Gauss."""

    def integrand(tau: float, u: float) -> float:
        psi = 0.5 * math.pi * tau ** 3
        dpsi_dtau = 1.5 * math.pi * tau * tau
        sine = math.sin(psi)
        if sine <= 0.0:
            return 0.0
        p = BREADTH_EXPONENT
        s = sine ** p
        ds_dx = p * sine ** (p - 1.0) * math.cos(psi) * (math.pi / length)
        element = math.sqrt(
            4.0 * u * u * draft * draft
            + (draft * beam * ds_dx * u * u) ** 2
            + (0.5 * beam * s) ** 2
        )
        return dpsi_dtau * element

    nodes = [index / panels for index in range(panels + 1)]
    weights = [
        (1.0 / panels) / 3.0 * (1.0 if index in (0, panels) else (4.0 if index % 2 else 2.0))
        for index in range(panels + 1)
    ]
    total = 0.0
    for i, tau in enumerate(nodes):
        for j, u in enumerate(nodes):
            total += weights[i] * weights[j] * integrand(tau, u)
    return (4.0 * length / math.pi) * total


def closed_form_audit() -> dict[str, object]:
    """The closed forms against an independent quadrature (a negative control)."""
    taper_simpson = taper_integral_by_simpson()
    gauss_area = analytic_wetted_surface_area(LENGTH, BEAM, DRAFT)
    simpson_area = wetted_area_by_simpson(LENGTH, BEAM, DRAFT)
    reference = analytic_hydrostatics(LENGTH, BEAM, DRAFT)
    return {
        "taper_integral": TAPER_INTEGRAL,
        "taper_integral_simpson": taper_simpson,
        "taper_integral_relative_difference": abs(TAPER_INTEGRAL / taper_simpson - 1.0),
        "volume_closed_form": analytic_volume(LENGTH, BEAM, DRAFT),
        "volume_rederived": LENGTH * BEAM * DRAFT * taper_simpson * (2.0 / 3.0),
        "wetted_area_gauss": gauss_area,
        "wetted_area_simpson": simpson_area,
        "wetted_area_relative_difference": abs(gauss_area / simpson_area - 1.0),
        "wetted_area_gauss_nodes_60": analytic_wetted_surface_area(LENGTH, BEAM, DRAFT, nodes=60),
        "wetted_area_gauss_nodes_120": analytic_wetted_surface_area(LENGTH, BEAM, DRAFT, nodes=120),
        "reference": reference.to_dict(),
        "checks": {
            "taper_integral_matches_simpson": abs(TAPER_INTEGRAL / taper_simpson - 1.0) < 1e-10,
            "wetted_area_quadrature_families_agree": abs(gauss_area / simpson_area - 1.0) < 1e-6,
            "wetted_area_gauss_converged": abs(
                analytic_wetted_surface_area(LENGTH, BEAM, DRAFT, nodes=120) / gauss_area - 1.0
            ) < 1e-8,
        },
    }


def sweep() -> list[dict[str, object]]:
    rows = []
    for nx, nz in GRIDS:
        report = analytic_relative_errors(HullParameters(nx=nx, nz=nz))
        rows.append(
            {
                "grid": f"{nx}x{nz}",
                "nx": nx,
                "nz": nz,
                "spacing": LENGTH / (nx - 1),
                "relative_error": report["relative_error"],
                "max_abs_relative_error": report["max_abs_relative_error"],
                "measured": report["measured"],
                "reference": report["reference"],
            }
        )
    return rows


def estimator_comparison(rows: list[dict[str, object]], grid: str) -> dict[str, object]:
    nx, nz = next((row["nx"], row["nz"]) for row in rows if row["grid"] == grid)
    parameters = HullParameters(nx=nx, nz=nz)
    estimate = estimate_discretisation_error(parameters)
    errors = next(row["relative_error"] for row in rows if row["grid"] == grid)
    rows_out = {}
    for name, reported in estimate.get("relative_change", {}).items():
        measured = abs(errors[name])
        rows_out[name] = {
            "reported_relative_change": reported,
            "measured_relative_error": measured,
            "ratio_measured_over_reported": (measured / reported) if reported else float("nan"),
        }
    return {"grid": grid, "quantities": rows_out}


def write_figure(rows: list[dict[str, object]], estimate: dict[str, object], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, (left, right) = plt.subplots(1, 2, figsize=(9.4, 3.8), dpi=110)
    spacings = [row["spacing"] for row in rows]
    plotted = {
        "displacement_volume": ("volume", "#1f77b4", "o"),
        "wetted_surface_area": ("wetted area", "#ff7f0e", "s"),
        "waterplane_area": ("waterplane area", "#2ca02c", "^"),
        "vertical_centre_buoyancy": ("VCB", "#d62728", "D"),
    }
    for name, (label, colour, marker) in plotted.items():
        values = [abs(row["relative_error"][name]) for row in rows]
        left.loglog(spacings, values, marker=marker, color=colour, label=label, lw=1.2)
    guide = [abs(rows[0]["relative_error"]["displacement_volume"]) * (s / spacings[0]) ** 2 for s in spacings]
    left.loglog(spacings, guide, "k--", lw=1.0, label="slope 2 guide")
    left.set_xlabel("station spacing [m]")
    left.set_ylabel("|error| vs closed form")
    left.set_title("sampled hull vs closed-form hydrostatics")
    left.grid(True, which="both", alpha=0.3)
    left.legend(fontsize=7)

    names = list(estimate["quantities"])
    reported = [estimate["quantities"][name]["reported_relative_change"] for name in names]
    measured = [estimate["quantities"][name]["measured_relative_error"] for name in names]
    positions = range(len(names))
    right.barh([p - 0.2 for p in positions], measured, height=0.4, label="error vs closed form", color="#1f77b4")
    right.barh([p + 0.2 for p in positions], reported, height=0.4, label="coarse/refined estimate", color="#ff7f0e")
    right.set_yticks(list(positions))
    right.set_yticklabels([name.replace("_", " ") for name in names], fontsize=7)
    right.set_xlabel("relative change")
    right.set_title(f"estimator vs truth at {estimate['grid']}")
    right.grid(True, axis="x", alpha=0.3)
    right.legend(fontsize=7)
    figure.tight_layout()
    figure.savefig(path)
    plt.close(figure)


def markdown_report(audit, rows, estimate) -> str:
    lines = []
    lines.append("# Sampled hull vs closed-form hydrostatics")
    lines.append("")
    lines.append(
        "Generated by `examples/analytic_hydrostatics/run_analytic_check.py`; every number "
        "below is written straight from the computation, none is typed by hand."
    )
    lines.append("")
    lines.append(
        f"Case: `L = {LENGTH} m`, `B = {BEAM} m`, `T = {DRAFT} m`, "
        f"half-breadth `(B/2)*sin(pi*x/L)^{BREADTH_EXPONENT}*sqrt((z+T)/T)`, seawater 1025 kg/m^3."
    )
    reference = audit["reference"]
    lines.append("")
    lines.append("## The reference itself")
    lines.append("")
    lines.append("| Quantity | Closed form / quadrature |")
    lines.append("|---|---|")
    for key, label in (
        ("displacement_volume_m3", "displacement volume [m³]"),
        ("waterplane_area_m2", "waterplane area [m²]"),
        ("maximum_section_area_m2", "midship section area [m²]"),
        ("longitudinal_centre_buoyancy_m", "LCB [m from aft end]"),
        ("longitudinal_centre_flotation_m", "LCF [m from aft end]"),
        ("vertical_centre_buoyancy_m", "VCB [m, +z up from the waterline]"),
        ("block_coefficient", "block coefficient"),
        ("prismatic_coefficient", "prismatic coefficient"),
        ("wetted_surface_area_m2", "wetted surface area [m²]"),
    ):
        lines.append(f"| {label} | {reference[key]:.9g} |")
    lines.append("")
    lines.append(
        f"Independent checks: `TAPER_INTEGRAL` vs Simpson "
        f"`{audit['taper_integral_simpson']:.12f}` → `{audit['taper_integral_relative_difference']:.2e}`; "
        f"wetted area Gauss vs Simpson `{audit['wetted_area_simpson']:.9f}` → "
        f"`{audit['wetted_area_relative_difference']:.2e}`."
    )
    lines.append("")
    lines.append("## Signed relative error of the mesh integrals")
    lines.append("")
    header = "| Quantity | " + " | ".join(row["grid"] for row in rows) + " |"
    lines.append(header)
    lines.append("|---" * (len(rows) + 1) + "|")
    for name in QUANTITIES:
        cells = " | ".join(f"{100 * row['relative_error'][name]:+.4f} %" for row in rows)
        lines.append(f"| {name.replace('_', ' ')} | {cells} |")
    lines.append("")
    lines.append("## Discretisation estimate vs the true error")
    lines.append("")
    lines.append("| Quantity | coarse/refined relative change | error vs closed form | ratio |")
    lines.append("|---|---|---|---|")
    for name, row in estimate["quantities"].items():
        lines.append(
            f"| {name.replace('_', ' ')} | {100 * row['reported_relative_change']:.4f} % | "
            f"{100 * row['measured_relative_error']:.4f} % | {row['ratio_measured_over_reported']:.2f} |"
        )
    lines.append("")
    lines.append("![convergence](convergence.png)")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "docs" / "analytic_hydrostatics",
        help="directory for results.json, report.md and convergence.png",
    )
    arguments = parser.parse_args()
    output = arguments.output
    output.mkdir(parents=True, exist_ok=True)

    audit = closed_form_audit()
    rows = sweep()
    estimate = estimator_comparison(rows, "41x13")

    default_row = next(row for row in rows if row["grid"] == "41x13")
    fine_row = next(row for row in rows if row["grid"] == "161x41")

    checks = dict(audit["checks"])
    checks["default_grid_inside_three_percent"] = default_row["max_abs_relative_error"] < DEFAULT_TOLERANCE
    checks["fine_grid_inside_point_two_percent"] = fine_row["max_abs_relative_error"] < FINE_TOLERANCE
    volume_series = [abs(row["relative_error"]["displacement_volume"]) for row in rows]
    checks["volume_error_is_monotone_in_h"] = volume_series == sorted(volume_series, reverse=True)
    checks["volume_is_under_reported"] = all(
        row["relative_error"]["displacement_volume"] < 0.0 for row in rows
    )
    checks["longitudinal_centre_is_the_midpoint"] = all(
        abs(row["measured"]["longitudinal_centre_buoyancy"] - LENGTH / 2.0) < 1e-9 for row in rows
    )

    results = {
        "case": {
            "length_m": LENGTH,
            "beam_m": BEAM,
            "draft_m": DRAFT,
            "breadth_exponent": BREADTH_EXPONENT,
            "water_density_kg_per_m3": 1025.0,
            "grids": ["x".join(map(str, grid)) for grid in GRIDS],
        },
        "closed_form_audit": audit,
        "sweep": rows,
        "estimator": estimate,
        "checks": checks,
        "passed": all(checks.values()),
    }
    (output / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True))
    (output / "report.md").write_text(markdown_report(audit, rows, estimate))
    write_figure(rows, estimate, output / "convergence.png")

    print(f"{'quantity':30s} " + " ".join(f"{row['grid']:>10s}" for row in rows))
    for name in QUANTITIES:
        cells = " ".join(f"{100 * row['relative_error'][name]:>+9.4f}%" for row in rows)
        print(f"{name:30s} {cells}")
    print()
    for name, passed in checks.items():
        print(f"{'PASS' if passed else 'FAIL'}  {name}")
    print()
    print("results.json, report.md and convergence.png written to", output)
    return 0 if results["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
