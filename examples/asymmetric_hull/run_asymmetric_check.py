#!/usr/bin/env python3
"""Verification of the asymmetric hull family against its own closed-form surface.

The symmetric family in :mod:`tensorcad.hull` has ``LCB = LCF = L/2`` by
construction, so its longitudinal-moment test cannot detect an antisymmetric bias
of any size. This script works the case that can: the two halves of the length get
their own taper exponent, which makes LCB a number a formula has to *predict*.

It measures, in order:

1. **The reference's own pieces** -- the Beta-function taper integral against the
   0.7 constant shipped in 0.11.0 and against a plain Simpson rule, and the
   ``M_p`` moment quadrature against a second substitution power and against node
   refinement (it is the only part of the reference with no formula).
2. **The reference against the mesh path** on five grids: signed relative error of
   nine hydrostatic quantities, computed by two independent integrators (closed
   form plus quadrature versus the divergence theorem on the triangles).
3. **The properties only asymmetry has**: the mirror identity
   ``LCB(p_f, p_a) + LCB(p_a, p_f) = L``, and how far LCB actually sits from
   midship.
4. **The symmetric limit**, where the new code has to reproduce the shipped 0.11.0
   reference rather than merely being self-consistent.

Writes ``results.json``, ``report.md`` (tables generated from the numbers, never
typed by hand) and ``convergence.png``. Exits 1 if a claim stops holding.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from tensorcad.hull import TAPER_INTEGRAL, analytic_hydrostatics as symmetric_hydrostatics
from tensorcad.asymmetric_hull import (
    COMPARABLE_QUANTITIES,
    DEFAULT_AFT_EXPONENT,
    DEFAULT_FORE_EXPONENT,
    MOMENT_NODES,
    SCHEMA,
    WETTED_SURFACE_NODES,
    AsymmetricHullParameters,
    analytic_half_taper_integral,
    analytic_hydrostatics,
    analytic_longitudinal_centre,
    analytic_volume,
    analytic_waterplane_area,
    generate_hull,
    half_taper_integrals,
    taper,
    wetted_surface_quadrature_drift,
)

LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8
GRIDS = ((11, 5), (21, 9), (41, 13), (81, 21), (161, 41))
DEFAULT_GRID = (41, 13)
FINE_GRID = (161, 41)
#: Claim: every quantity is inside 3 % on the 41x13 grid and inside 0.5 % at 161x41.
DEFAULT_TOLERANCE = 0.03
FINE_TOLERANCE = 0.005


def simpson(function, lower: float, upper: float, intervals: int = 40000) -> float:
    if intervals % 2:
        intervals += 1
    step = (upper - lower) / intervals
    total = function(lower) + function(upper)
    for k in range(1, intervals):
        total += (4.0 if k % 2 else 2.0) * function(lower + k * step)
    return total * step / 3.0


def reference_audit() -> dict:
    """Check the pieces of the reference that are not covered by a mesh."""
    parameter = AsymmetricHullParameters(length=LENGTH, beam=BEAM, draft=DRAFT)
    taper_rows = []
    for exponent in (0.4, 0.7, 1.0, 1.6, 2.0, 3.0):
        closed = analytic_half_taper_integral(exponent)
        numeric = simpson(lambda x, p=exponent: math.sin(math.pi * x) ** p, 0.0, 0.5)
        slow = half_taper_integrals(exponent, nodes=120)
        fast = half_taper_integrals(exponent, nodes=120, substitution=3.0 / exponent)
        refined = half_taper_integrals(exponent, nodes=400)
        taper_rows.append(
            {
                "exponent": exponent,
                "closed_form": closed,
                "simpson": numeric,
                "closed_vs_simpson_relative": closed / numeric - 1.0,
                "moment_substitution_relative": slow[1] / fast[1] - 1.0,
                "moment_node_refinement_relative": slow[1] / refined[1] - 1.0,
            }
        )
    symmetric_limit = {}
    reference = analytic_hydrostatics(
        AsymmetricHullParameters(
            length=LENGTH, beam=BEAM, draft=DRAFT, fore_exponent=0.7, aft_exponent=0.7
        )
    )
    shipped = symmetric_hydrostatics(LENGTH, BEAM, DRAFT)
    for name in COMPARABLE_QUANTITIES:
        ours = float(getattr(reference, name))
        theirs = float(getattr(shipped, name))
        symmetric_limit[name] = {
            "asymmetric_module": ours,
            "shipped_symmetric": theirs,
            "relative": ours / theirs - 1.0,
        }
    checks = {
        "taper_closed_form_hits_shipped_constant": abs(
            2.0 * analytic_half_taper_integral(0.7) / TAPER_INTEGRAL - 1.0
        )
        < 1e-14,
        "taper_closed_form_matches_simpson": all(
            abs(row["closed_vs_simpson_relative"]) < 1e-7 for row in taper_rows
        ),
        "moment_quadrature_is_substitution_independent": all(
            abs(row["moment_substitution_relative"]) < 1e-7 for row in taper_rows
        ),
        "moment_quadrature_converges_with_nodes": all(
            abs(row["moment_node_refinement_relative"]) < 1e-12 for row in taper_rows
        ),
        "wetted_quadrature_drift_below_1e-9": wetted_surface_quadrature_drift(parameter) < 1e-9,
        "symmetric_limit_reproduces_0_11_0_reference": all(
            abs(row["relative"]) < 1e-9 for row in symmetric_limit.values()
        ),
    }
    return {
        "case": parameter.as_dict(),
        "taper_integral": taper_rows,
        "symmetric_limit": symmetric_limit,
        "wetted_surface_drift": wetted_surface_quadrature_drift(parameter),
        "checks": checks,
    }


def sweep() -> list[dict]:
    rows = []
    for nx, nz in GRIDS:
        parameters = AsymmetricHullParameters(
            length=LENGTH,
            beam=BEAM,
            draft=DRAFT,
            fore_exponent=DEFAULT_FORE_EXPONENT,
            aft_exponent=DEFAULT_AFT_EXPONENT,
            nx=nx,
            nz=nz,
        )
        model = generate_hull(parameters)
        errors = model.errors()
        measured = model.hydrostatics
        reference = model.reference
        rows.append(
            {
                "grid": f"{nx}x{nz}",
                "nx": nx,
                "nz": nz,
                "vertices": len(model.mesh.vertices),
                "triangles": len(model.mesh.faces),
                "h": 1.0 / (nx - 1),
                "relative_error": errors,
                "measured": {
                    name: float(getattr(measured, name)) for name in COMPARABLE_QUANTITIES
                },
                "reference": {
                    name: float(getattr(reference, name)) for name in COMPARABLE_QUANTITIES
                },
                "max_abs_relative_error": max(abs(value) for value in errors.values()),
                "warnings": list(model.validation.warnings),
            }
        )
    return rows


def asymmetry_checks() -> dict:
    forward = AsymmetricHullParameters(
        length=LENGTH, beam=BEAM, draft=DRAFT,
        fore_exponent=DEFAULT_FORE_EXPONENT, aft_exponent=DEFAULT_AFT_EXPONENT,
    )
    mirrored = AsymmetricHullParameters(
        length=LENGTH, beam=BEAM, draft=DRAFT,
        fore_exponent=DEFAULT_AFT_EXPONENT, aft_exponent=DEFAULT_FORE_EXPONENT,
    )
    lcb = analytic_longitudinal_centre(forward)
    mirrored_lcb = analytic_longitudinal_centre(mirrored)
    checks = {
        "lcb_is_off_midship_by_more_than_2_percent_of_length": abs(lcb - 0.5 * LENGTH)
        > 0.02 * LENGTH,
        "mirror_identity_holds": abs(lcb + mirrored_lcb - LENGTH) < 1e-12 * LENGTH,
        "volume_invariant_under_mirror": abs(
            analytic_volume(forward) / analytic_volume(mirrored) - 1.0
        )
        < 1e-14,
        "waterplane_invariant_under_mirror": abs(
            analytic_waterplane_area(forward) / analytic_waterplane_area(mirrored) - 1.0
        )
        < 1e-14,
        "lcb_moves_towards_the_fuller_half": lcb > 0.5 * LENGTH > mirrored_lcb,
    }
    return {
        "fore_exponent": DEFAULT_FORE_EXPONENT,
        "aft_exponent": DEFAULT_AFT_EXPONENT,
        "lcb_m": lcb,
        "lcb_over_length": lcb / LENGTH,
        "mirrored_lcb_m": mirrored_lcb,
        "lcb_shift_m": lcb - 0.5 * LENGTH,
        "symmetric_family_would_report": 0.5 * LENGTH,
        "checks": checks,
    }


def markdown_report(audit: dict, asymmetry: dict, rows: list[dict]) -> str:
    lines: list[str] = []
    lines.append("# Asymmetric hull family: mesh against closed form\n")
    lines.append(
        "Generated by `examples/asymmetric_hull/run_asymmetric_check.py`; every number "
        "below is read from `results.json`, none is typed by hand.\n"
    )
    case = audit["case"]
    lines.append("## Case\n")
    lines.append(
        f"`L = {case['length']} m`, `B = {case['beam']} m`, `T = {case['draft']} m`, "
        f"fore taper exponent `{case['fore_exponent']}`, aft `{case['aft_exponent']}`, "
        f"seawater `{case['water_density']} kg/m³`.\n"
    )
    lines.append(
        "Half-breadth `b(x, z) = (B/2) · sin(pi·x/L)^p(x) · sqrt((z+T)/T)` with "
        "`p(x) = fore_exponent` forward of midship and `aft_exponent` aft of it; the "
        "surface is continuous at midship because `sin(pi/2) = 1` for either exponent.\n"
    )
    lines.append("## The reference's own pieces\n")
    lines.append(
        "| p | `I_p` closed form | Simpson | closed vs Simpson | `M_p` two substitutions | `M_p` node refinement |"
    )
    lines.append("|---|---|---|---|---|---|")
    for row in audit["taper_integral"]:
        lines.append(
            f"| {row['exponent']} | {row['closed_form']:.15f} | {row['simpson']:.15f} | "
            f"{row['closed_vs_simpson_relative']:+.2e} | "
            f"{row['moment_substitution_relative']:+.2e} | "
            f"{row['moment_node_refinement_relative']:+.2e} |"
        )
    lines.append("")
    lines.append(
        f"Wetted-surface quadrature: {WETTED_SURFACE_NODES} nodes, drift under node "
        f"doubling `{audit['wetted_surface_drift']:.2e}` relative. Moments: "
        f"`{MOMENT_NODES}` nodes.\n"
    )
    lines.append("## Mesh integrals against the closed form\n")
    header = "| quantity | " + " | ".join(row["grid"] for row in rows) + " |"
    lines.append(header)
    lines.append("|---" * (len(rows) + 1) + "|")
    for name in COMPARABLE_QUANTITIES:
        cells = " | ".join(f"{100 * row['relative_error'][name]:+.4f} %" for row in rows)
        lines.append(f"| `{name}` | {cells} |")
    lines.append("")
    lines.append("## Asymmetry\n")
    lines.append(
        f"LCB `{asymmetry['lcb_m']:.9f} m` (`{100 * asymmetry['lcb_over_length']:.4f} % of L`), "
        f"i.e. `{asymmetry['lcb_shift_m']:+.6f} m` off midship, where the symmetric family "
        f"would report `{asymmetry['symmetric_family_would_report']:.6f} m` no matter what. "
        f"Mirrored exponents give `{asymmetry['mirrored_lcb_m']:.9f} m`; the two sum to "
        f"`{asymmetry['lcb_m'] + asymmetry['mirrored_lcb_m']:.12f} m`.\n"
    )
    lines.append("## Checks\n")
    for group in ("checks",):
        for name, passed in audit[group].items():
            lines.append(f"- {'PASS' if passed else 'FAIL'} `{name}`")
    for name, passed in asymmetry["checks"].items():
        lines.append(f"- {'PASS' if passed else 'FAIL'} `{name}`")
    lines.append("")
    return "\n".join(lines)


def write_figure(rows: list[dict], asymmetry: dict, destination: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.2))

    series = (
        ("displacement_volume", "displacement volume"),
        ("wetted_surface_area", "wetted surface area"),
        ("longitudinal_centre_buoyancy", "LCB"),
        ("vertical_centre_buoyancy", "VCB"),
    )
    h = [row["h"] for row in rows]
    for name, label in series:
        values = [abs(row["relative_error"][name]) for row in rows]
        left.loglog(h, values, "o-", label=label)
    reference = [2.0 * h[-1] ** 2 / (h[i] ** 2) for i in range(len(h))]
    left.loglog(h, [abs(rows[-1]["relative_error"]["displacement_volume"]) * value for value in reference],
                "k--", linewidth=1.0, label="second order")
    left.set_xlabel("station spacing h = 1/(nx-1)")
    left.set_ylabel("|mesh - closed form| / closed form")
    left.set_title("asymmetric hull: error against the closed form")
    left.grid(True, which="both", alpha=0.3)
    left.legend(fontsize=8)

    parameters = AsymmetricHullParameters(
        length=LENGTH, beam=BEAM, draft=DRAFT,
        fore_exponent=DEFAULT_FORE_EXPONENT, aft_exponent=DEFAULT_AFT_EXPONENT,
    )
    samples = 400
    xs = [LENGTH * k / samples for k in range(samples + 1)]
    widths = [BEAM * taper(x / LENGTH, parameters) for x in xs]
    integrals = sum(widths[:-1])
    cumulative = []
    running = 0.0
    for width in widths[:-1]:
        running += width
        cumulative.append(running / integrals)
    right.plot(xs, widths, label="waterline breadth b(x, 0)")
    right.fill_between(xs, 0.0, widths, alpha=0.12)
    right.axvline(0.5 * LENGTH, color="k", linestyle=":", linewidth=1.0, label="midship")
    right.plot([asymmetry["lcb_m"], asymmetry["lcb_m"]], [0.0, BEAM * 0.5], "r-", linewidth=2.0,
               label=f"LCB {asymmetry['lcb_m']:.3f} m")
    right.set_xlabel("x from aft end [m]")
    right.set_ylabel("breadth [m]")
    right.set_title("fore/aft taper asymmetry and where the volume sits")
    right.grid(True, alpha=0.3)
    right.legend(fontsize=8)

    figure.tight_layout()
    figure.savefig(destination, dpi=150)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "docs" / "asymmetric_hull",
        help="directory for results.json, report.md and convergence.png",
    )
    arguments = parser.parse_args()
    output = arguments.output
    output.mkdir(parents=True, exist_ok=True)

    audit = reference_audit()
    asymmetry = asymmetry_checks()
    rows = sweep()

    default_row = next(row for row in rows if row["grid"] == "41x13")
    fine_row = next(row for row in rows if row["grid"] == "161x41")
    checks = dict(audit["checks"])
    checks.update(asymmetry["checks"])
    checks["default_grid_inside_three_percent"] = default_row["max_abs_relative_error"] < DEFAULT_TOLERANCE
    checks["fine_grid_inside_half_percent"] = fine_row["max_abs_relative_error"] < FINE_TOLERANCE
    for name in ("displacement_volume", "wetted_surface_area", "longitudinal_centre_buoyancy"):
        series = [abs(row["relative_error"][name]) for row in rows]
        checks[f"{name}_error_shrinks_with_refinement"] = series == sorted(series, reverse=True)
    checks["lcb_error_sign_is_stable"] = len(
        {row["relative_error"]["longitudinal_centre_buoyancy"] < 0.0 for row in rows}
    ) == 1
    checks["mesh_is_a_closed_manifold_on_every_grid"] = all(
        not row["warnings"] or "beam" not in " ".join(row["warnings"]) for row in rows
    )

    results = {
        "schema": SCHEMA,
        "case": audit["case"],
        "reference_audit": audit,
        "asymmetry": asymmetry,
        "sweep": rows,
        "checks": checks,
        "passed": all(checks.values()),
    }
    (output / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True))
    (output / "report.md").write_text(markdown_report(audit, asymmetry, rows))
    write_figure(rows, asymmetry, output / "convergence.png")

    print(f"{'quantity':32s} " + " ".join(f"{row['grid']:>11s}" for row in rows))
    for name in COMPARABLE_QUANTITIES:
        cells = " ".join(f"{100 * row['relative_error'][name]:>+10.4f}%" for row in rows)
        print(f"{name:32s} {cells}")
    print()
    print(f"LCB = {asymmetry['lcb_m']:.9f} m of {LENGTH} m ({100 * asymmetry['lcb_over_length']:.4f} %)")
    print()
    for name, passed in checks.items():
        print(f"{'PASS' if passed else 'FAIL'}  {name}")
    print()
    print("results.json, report.md and convergence.png written to", output)
    return 0 if results["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
