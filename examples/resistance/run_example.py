#!/usr/bin/env python3
"""Resistance from geometry: friction lines, a Prohaska form factor, and scaling to a ship.

Four studies, each with a reference the answer can be wrong against:

1. **Three flat-plate lines, three independent estimates of one number.** The 1957 ITTC
   line is checked against Schoenherr's implicit line (solved by bisection here, not
   tabulated) and against Prandtl-Schlichting's ``0.455/(log10 Re)^2.58``. The measured
   distance between them as a function of Reynolds number *is* the result; the two
   divergences at low ``Re`` (16 % at ``1e5``) are a property of the lines, not a defect,
   and are reported rather than smoothed away.
2. **The wetted area that feeds the friction comes off the mesh.** Wigley hull, three
   grids: shell area minus the waterline deck patch against the analytic flank integral
   (:func:`tensorcad.wigley.analytic_wetted_surface_area`), which is a quadrature of the
   analytic surface and never touches the mesh.
3. **A form factor recovered from a resistance curve.** Synthetic "model test" totals are
   generated from a known ``(1 + k)`` and wave coefficient; the Prohaska fit has to give
   them back. This is the method's own self-test: if the transform or the regression is
   wrong, the recovered intercept is wrong by a number a designer would act on.
4. **Model to ship.** A 10 m Wigley at ``Fn = 0.30`` carried to 100 m by Froude scaling
   and the ITTC-78 two-dimensional split, with the residual coefficient taken from the
   model total. The reported scale effect is the reason a ship is not ``lambda^3`` times
   its model.

Writes ``results.json``, ``report.md`` and ``plot.png`` into ``--output`` and exits
non-zero if any claim stops holding.

Usage::

    python examples/resistance/run_example.py [--output artifacts/resistance]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:  # allow running from a checkout without installing
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.mesh import surface_area, triangle_areas
from tensorcad.resistance import (
    FRESH_WATER_DENSITY,
    FRESH_WATER_KINEMATIC_VISCOSITY,
    SEA_WATER_DENSITY,
    SEA_WATER_KINEMATIC_VISCOSITY,
    bare_hull_resistance,
    froude_number,
    ittc57_friction_coefficient,
    prandtl_schlichting_friction_coefficient,
    prohaska_fit,
    reynolds_number,
    schoenherr_friction_coefficient,
    ship_from_model_test,
)
from tensorcad.wigley import (
    analytic_wetted_surface_area,
    generate_wigley,
    reference_parameters,
)

REYNOLDS_GRID = (1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8, 1.0e9)
CAREY_GRID = (41, 81, 161)  # nx; nz chosen as nx // 4 + ... below
MODEL_LENGTH = 10.0
SHIP_LENGTH = 100.0
INJECTED_FORM_FACTOR = 1.15
INJECTED_WAVE_COEFFICIENT = 3.0e-3
PROHASKA_FROUDES = (0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24)
SHIP_FROUDE = 0.30


def _deck_area(mesh) -> float:
    areas = triangle_areas(mesh)
    return math.fsum(
        area
        for face, area in zip(mesh.faces, areas)
        if max(abs(mesh.vertices[index][2]) for index in face) <= 1e-12
    )


def friction_study() -> dict:
    rows = []
    for re in REYNOLDS_GRID:
        ittc = ittc57_friction_coefficient(re)
        schoenherr = schoenherr_friction_coefficient(re)
        prandtl = prandtl_schlichting_friction_coefficient(re)
        rows.append(
            {
                "reynolds": re,
                "ittc57": ittc,
                "schoenherr": schoenherr,
                "prandtl_schlichting": prandtl,
                "relative_difference_schoenherr": (ittc - schoenherr) / schoenherr,
                "relative_difference_prandtl_schlichting": (ittc - prandtl) / prandtl,
                "schoenherr_equation_residual": 0.242 / math.sqrt(schoenherr)
                - math.log10(re * schoenherr),
            }
        )
    ship_range = [row for row in rows if row["reynolds"] >= 1.0e7]
    return {
        "rows": rows,
        "worst_relative_difference_schoenherr_ship_range": max(
            abs(row["relative_difference_schoenherr"]) for row in ship_range
        ),
        "worst_relative_difference_schoenherr_full_range": max(
            abs(row["relative_difference_schoenherr"]) for row in rows
        ),
        "worst_relative_difference_prandtl_schlichting_ship_range": max(
            abs(row["relative_difference_prandtl_schlichting"]) for row in ship_range
        ),
        "worst_absolute_schoenherr_residual": max(
            abs(row["schoenherr_equation_residual"]) for row in rows
        ),
    }


def wetted_area_study() -> dict:
    parameters = reference_parameters()
    reference = analytic_wetted_surface_area(
        parameters.length,
        parameters.beam,
        parameters.draft,
        parameters.longitudinal_exponent,
        parameters.vertical_exponent,
    )
    rows = []
    for nx in CAREY_GRID:
        nz = nx // 4 + 3
        model = generate_wigley(parameters, nx=nx, nz=nz)
        shell = surface_area(model.mesh)
        deck = _deck_area(model.mesh)
        wetted = shell - deck
        rows.append(
            {
                "nx": nx,
                "nz": nz,
                "shell_area_m2": shell,
                "deck_area_m2": deck,
                "deck_share_of_shell": deck / shell,
                "wetted_area_m2": wetted,
                "analytic_wetted_area_m2": reference,
                "relative_error": (wetted - reference) / reference,
                "friction_inflation_if_shell_used": shell / wetted - 1.0,
            }
        )
    return {
        "parameters": {
            "length_m": parameters.length,
            "beam_m": parameters.beam,
            "draft_m": parameters.draft,
            "longitudinal_exponent": parameters.longitudinal_exponent,
            "vertical_exponent": parameters.vertical_exponent,
        },
        "rows": rows,
        "analytic_wetted_area_m2": reference,
        "finest_relative_error": rows[-1]["relative_error"],
        "worst_deck_share": max(row["deck_share_of_shell"] for row in rows),
    }


def prohaska_study() -> dict:
    froudes = PROHASKA_FROUDES
    friction = tuple(
        ittc57_friction_coefficient(
            reynolds_number(
                MODEL_LENGTH,
                fn * math.sqrt(9.80665 * MODEL_LENGTH),
                SEA_WATER_KINEMATIC_VISCOSITY,
            )
        )
        for fn in froudes
    )
    totals = tuple(
        INJECTED_FORM_FACTOR * cf + INJECTED_WAVE_COEFFICIENT * fn**4
        for fn, cf in zip(froudes, friction)
    )
    fit = prohaska_fit(froudes, totals, length=MODEL_LENGTH)
    return {
        "model_length_m": MODEL_LENGTH,
        "froude_numbers": list(froudes),
        "measured_total_coefficients": list(totals),
        "friction_coefficients": list(fit.friction_coefficients),
        "injected_form_factor": INJECTED_FORM_FACTOR,
        "injected_wave_coefficient": INJECTED_WAVE_COEFFICIENT,
        "recovered_form_factor": fit.form_factor,
        "recovered_wave_coefficient": fit.wave_coefficient,
        "relative_error_form_factor": abs(fit.form_factor - INJECTED_FORM_FACTOR)
        / INJECTED_FORM_FACTOR,
        "max_relative_residual": fit.max_relative_residual,
        "prohaska_points": [list(point) for point in fit.points],
        "fit": fit.to_dict(),
    }


def scaling_study(wetted_area: float) -> dict:
    model_velocity = SHIP_FROUDE * math.sqrt(9.80665 * MODEL_LENGTH)
    model_reynolds = reynolds_number(
        MODEL_LENGTH, model_velocity, SEA_WATER_KINEMATIC_VISCOSITY
    )
    model_total = INJECTED_FORM_FACTOR * ittc57_friction_coefficient(
        model_reynolds
    ) + INJECTED_WAVE_COEFFICIENT * SHIP_FROUDE**4
    scale = SHIP_LENGTH / MODEL_LENGTH
    prediction = ship_from_model_test(
        model_length=MODEL_LENGTH,
        ship_length=SHIP_LENGTH,
        model_velocity=model_velocity,
        model_wetted_surface=wetted_area,
        ship_wetted_surface=wetted_area * scale**2,
        form_factor=INJECTED_FORM_FACTOR,
        model_total_coefficient=model_total,
        correlation_allowance=0.0,
    )
    check = bare_hull_resistance(
        length=SHIP_LENGTH,
        velocity=model_velocity * math.sqrt(scale),
        wetted_surface=wetted_area * scale**2,
        form_factor=INJECTED_FORM_FACTOR,
        residual_coefficient=INJECTED_WAVE_COEFFICIENT * SHIP_FROUDE**4,
    )
    speeds = [SHIP_FROUDE * math.sqrt(9.80665 * MODEL_LENGTH) * f for f in (0.5, 0.75, 1.0, 1.25, 1.5)]
    curve = [
        bare_hull_resistance(
            length=MODEL_LENGTH,
            velocity=velocity,
            wetted_surface=wetted_area,
            form_factor=INJECTED_FORM_FACTOR,
            residual_coefficient=INJECTED_WAVE_COEFFICIENT
            * froude_number(MODEL_LENGTH, velocity) ** 4,
        ).to_dict()
        for velocity in speeds
    ]
    return {
        "model_length_m": MODEL_LENGTH,
        "ship_length_m": SHIP_LENGTH,
        "froude": SHIP_FROUDE,
        "model_velocity_m_s": model_velocity,
        "model_reynolds": model_reynolds,
        "model_total_coefficient": model_total,
        "prediction": prediction.to_dict(),
        "direct_ship_route_total_resistance_n": check.total_resistance,
        "routes_agree_relative": abs(check.total_resistance - prediction.ship.total_resistance)
        / prediction.ship.total_resistance,
        "naive_lambda_cubed_ratio": scale**3,
        "model_speed_sweep": curve,
    }


def claims(results: dict) -> list[dict]:
    friction = results["friction"]
    area = results["wetted_area"]
    prohaska = results["prohaska"]
    scaling = results["scaling"]
    total = [row["total_resistance_n"] for row in scaling["model_speed_sweep"]]
    return [
        {
            "claim": "ITTC-57 within 2.3 % of Schoenherr for Re >= 1e7",
            "measured": friction["worst_relative_difference_schoenherr_ship_range"],
            "passed": friction["worst_relative_difference_schoenherr_ship_range"] < 0.023,
        },
        {
            "claim": "the Schoenherr bisection satisfies its own implicit equation to 1e-12",
            "measured": friction["worst_absolute_schoenherr_residual"],
            "passed": friction["worst_absolute_schoenherr_residual"] < 1e-12,
        },
        {
            "claim": "the low-Re divergence of the two lines is present, not hidden (16 % at 1e5)",
            "measured": abs(friction["rows"][1]["relative_difference_schoenherr"]),
            "passed": 0.15 < abs(friction["rows"][1]["relative_difference_schoenherr"]) < 0.17,
        },
        {
            "claim": "mesh wetted area converges to the analytic flank integral, finest grid < 0.01 %",
            "measured": area["finest_relative_error"],
            "passed": abs(area["finest_relative_error"]) < 1.0e-4,
        },
        {
            "claim": "the deck patch is a third of the shell area, so shell != wetted",
            "measured": area["worst_deck_share"],
            "passed": 0.30 < area["worst_deck_share"] < 0.32,
        },
        {
            "claim": "Prohaska recovers the injected form factor to 1e-12 relative",
            "measured": prohaska["relative_error_form_factor"],
            "passed": prohaska["relative_error_form_factor"] < 1e-12,
        },
        {
            "claim": "model and ship Froude numbers match under sqrt(lambda) scaling",
            "measured": abs(
                scaling["prediction"]["ship"]["froude"] - scaling["prediction"]["model"]["froude"]
            ),
            "passed": abs(
                scaling["prediction"]["ship"]["froude"] - scaling["prediction"]["model"]["froude"]
            )
            < 1e-13,
        },
        {
            "claim": "the scaled ship is not lambda^3 times the model",
            "measured": scaling["prediction"]["ship_resistance_ratio"],
            "passed": scaling["prediction"]["ship_resistance_ratio"] < scaling["naive_lambda_cubed_ratio"],
        },
        {
            "claim": "resistance rises monotonically with speed",
            "measured": total,
            "passed": all(later > earlier for earlier, later in zip(total, total[1:])),
        },
    ]


def plot(results: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    friction = results["friction"]["rows"]
    re = [row["reynolds"] for row in friction]
    figure, (left, right) = plt.subplots(1, 2, figsize=(12.0, 4.8))

    left.loglog(re, [row["ittc57"] for row in friction], "o-", label="ITTC-57 (used)")
    left.loglog(re, [row["schoenherr"] for row in friction], "s--", label="Schoenherr 1932 (implicit)")
    left.loglog(re, [row["prandtl_schlichting"] for row in friction], "^:", label="Prandtl-Schlichting")
    left.set_xlabel("Reynolds number Re")
    left.set_ylabel("friction coefficient $C_F$")
    left.set_title("Three flat-plate lines")
    left.grid(True, which="both", alpha=0.3)
    left.legend(loc="lower left", fontsize=8)

    points = results["prohaska"]["prohaska_points"]
    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    fit = results["prohaska"]["fit"]
    right.plot(xs, ys, "o", label="synthetic model test")
    span = [min(xs), max(xs)]
    right.plot(
        span,
        [fit["form_factor"] + fit["wave_coefficient"] * x for x in span],
        "-",
        label=f"fit: 1+k = {fit['form_factor']:.9f}",
    )
    right.axhline(fit["form_factor"], color="grey", ls=":", lw=1)
    right.set_xlabel("$Fn^4 / C_F$")
    right.set_ylabel("$C_T / C_F$")
    right.set_title("Prohaska plot, model length 10 m")
    right.grid(True, alpha=0.3)
    right.legend(fontsize=8)

    figure.tight_layout()
    figure.savefig(path, dpi=140)
    plt.close(figure)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(REPO_ROOT / "artifacts" / "resistance"))
    arguments = parser.parse_args(argv)
    output = Path(arguments.output)
    output.mkdir(parents=True, exist_ok=True)

    results = {
        "friction": friction_study(),
        "wetted_area": wetted_area_study(),
        "prohaska": prohaska_study(),
        "fluids": {
            "sea_water_density": SEA_WATER_DENSITY,
            "sea_water_kinematic_viscosity": SEA_WATER_KINEMATIC_VISCOSITY,
            "fresh_water_density": FRESH_WATER_DENSITY,
            "fresh_water_kinematic_viscosity": FRESH_WATER_KINEMATIC_VISCOSITY,
        },
    }
    results["scaling"] = scaling_study(results["wetted_area"]["analytic_wetted_area_m2"])
    results["claims"] = claims(results)

    (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    plot(results, output / "plot.png")

    failed = [claim for claim in results["claims"] if not claim["passed"]]
    lines = ["# Resistance example: measured claims", ""]
    lines += ["| claim | measured | pass |", "|---|---|---|"]
    for claim in results["claims"]:
        measured = claim["measured"]
        shown = f"{measured:.6g}" if isinstance(measured, float) else str(measured)
        lines.append(f"| {claim['claim']} | {shown} | {'PASS' if claim['passed'] else 'FAIL'} |")
    (output / "report.md").write_text("\n".join(lines) + "\n")

    for claim in results["claims"]:
        print(f"[{'PASS' if claim['passed'] else 'FAIL'}] {claim['claim']}: {claim['measured']}")
    if failed:
        print(f"\n{len(failed)} claim(s) failed", file=sys.stderr)
        return 1
    print(f"\nwrote {output / 'results.json'}, {output / 'plot.png'}, {output / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())