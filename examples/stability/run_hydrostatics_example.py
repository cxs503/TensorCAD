"""Verify the waterplane inertias and the metacentric heights against closed forms.

Builds the resolution ladder of the parametric hull, integrates the two waterplane
second moments off the mesh, and compares them with the analytic surface (one gamma
identity, one tanh-sinh quadrature). Writes ``results.json``, ``report.md`` and a
two-panel figure, and exits non-zero if any claim stops holding.

    python examples/stability/run_hydrostatics_example.py
    python examples/stability/run_hydrostatics_example.py --output artifacts/stability_hydrostatics
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tensorcad import (
    HullParameters,
    analytic_hydrostatics,
    generate_hull,
    initial_stability_from_hydrostatics,
)

#: Resolution ladder: (nx, nz). Fine enough to see the convergence order.
LADDER = ((21, 9), (41, 13), (81, 21), (161, 41), (321, 81))

#: The default grid a user gets from ``generate_hull()``.
DEFAULT_GRID = (41, 13)
FINE_GRID = (161, 41)

#: Acceptance bands on the relative error of every stability quantity.
DEFAULT_CLAIM = 0.03
FINE_CLAIM = 0.005

#: KG used for the metacentric height reported below, in metres above the keel.
HEIGHT_OF_GRAVITY = 0.5

#: Quantities compared, in report order.
QUANTITIES = (
    "transverse_waterplane_inertia",
    "longitudinal_waterplane_inertia",
    "transverse_metacentric_radius",
    "longitudinal_metacentric_radius",
    "height_of_metacentre_above_keel",
    "transverse_metacentric_height",
)


def relative_error(measured: float, reference: float) -> float:
    return measured / reference - 1.0


def measure(parameters: HullParameters) -> dict[str, float]:
    """Mesh hydrostatics and the derived stability block for one resolution."""
    model = generate_hull(parameters)
    hydro = model.hydrostatics
    stability = initial_stability_from_hydrostatics(hydro, height_of_gravity=HEIGHT_OF_GRAVITY)
    volume = hydro.displacement_volume
    return {
        "transverse_waterplane_inertia": hydro.transverse_waterplane_inertia,
        "longitudinal_waterplane_inertia": hydro.longitudinal_waterplane_inertia,
        "transverse_metacentric_radius": hydro.transverse_waterplane_inertia / volume,
        "longitudinal_metacentric_radius": hydro.longitudinal_waterplane_inertia / volume,
        "height_of_metacentre_above_keel": stability.height_of_transverse_metacentre,
        "transverse_metacentric_height": stability.transverse_metacentric_height,
        "block_coefficient": hydro.block_coefficient,
        "displacement_volume": volume,
        "waterplane_area": hydro.waterplane_area,
        "keel_z": hydro.keel_z,
        "kb": stability.height_of_centre_of_buoyancy,
        "draft": stability.draft,
    }


def reference_values(parameters: HullParameters) -> dict[str, float]:
    """Closed-form values of the same quantities on the analytic surface."""
    reference = analytic_hydrostatics(parameters)
    volume = reference.displacement_volume
    return {
        "transverse_waterplane_inertia": reference.transverse_waterplane_inertia,
        "longitudinal_waterplane_inertia": reference.longitudinal_waterplane_inertia,
        "transverse_metacentric_radius": reference.transverse_metacentric_radius,
        "longitudinal_metacentric_radius": reference.longitudinal_metacentric_radius,
        "height_of_metacentre_above_keel": reference.height_of_metacentre_above_keel,
        "transverse_metacentric_height": reference.height_of_metacentre_above_keel
        - HEIGHT_OF_GRAVITY,
        "block_coefficient": reference.block_coefficient,
        "displacement_volume": volume,
        "waterplane_area": reference.waterplane_area,
        "kb": reference.vertical_centre_buoyancy + parameters.draft,
        "draft": parameters.draft,
    }


def make_figure(
    parameters: HullParameters,
    rows: list[dict[str, float]],
    reference: dict[str, float],
    destination: Path,
) -> None:
    """Left: the waterplane and the two integrands. Right: error against resolution."""
    length, beam, draft = parameters.length, parameters.beam, parameters.draft
    exponent = 0.7
    stations = [length * i / 400.0 for i in range(401)]
    half_breadth = [(beam / 2.0) * math.sin(math.pi * x / length) ** exponent for x in stations]
    transverse_integrand = [(2.0 / 3.0) * b**3 for b in half_breadth]
    longitudinal_integrand = [2.0 * b * (x - length / 2.0) ** 2 for b, x in zip(half_breadth, stations, strict=True)]

    figure, (left, right) = plt.subplots(1, 2, figsize=(12.0, 4.6))

    left.plot(stations, half_breadth, "k-", label="half-breadth b(x) [m]")
    left.plot(stations, [value / 2.0 for value in transverse_integrand], "C0--",
              label="dI_T/dx, int (2/3)b$^3$ / 2 [m$^3$]")
    left.plot(stations, [value / 100.0 for value in longitudinal_integrand], "C1:",
              label="dI$_L$/dx, 2b(x-L/2)$^2$ / 100 [m$^4$]")
    left.set_xlabel("x [m]  (aft -> fore)")
    left.set_ylabel("m, or scaled integrand")
    left.set_title("Waterplane of the analytic hull (B/2 = %.2f m)" % (beam / 2.0))
    left.grid(True, alpha=0.3)
    left.legend(fontsize=8, loc="upper right")

    spacings = [parameters.length / (row["nx"] - 1) for row in rows]
    for name, style in (
        ("transverse_waterplane_inertia", "o-"),
        ("longitudinal_waterplane_inertia", "s-"),
        ("height_of_metacentre_above_keel", "^-"),
        ("displacement_volume", "d-"),
    ):
        errors = [abs(relative_error(row[name], reference[name])) for row in rows]
        right.loglog(spacings, errors, style, label=name.replace("_", " "))
    right.axhline(DEFAULT_CLAIM, color="k", linestyle="--", linewidth=1.0,
                  label="3 %% acceptance band")
    right.set_xlabel("longitudinal spacing h = L/(nx-1) [m]")
    right.set_ylabel("|relative error| vs closed form")
    right.set_title(
        f"Mesh vs analytic reference (L={length:.1f}, B={beam:.1f}, T={draft:.2f})"
    )
    right.grid(True, which="both", alpha=0.3)
    right.legend(fontsize=8)

    figure.tight_layout()
    figure.savefig(destination, dpi=130)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="docs/stability_hydrostatics", help="directory for the artifacts")
    arguments = parser.parse_args()
    output = Path(arguments.output)
    output.mkdir(parents=True, exist_ok=True)

    parameters = HullParameters()
    reference = reference_values(parameters)
    rows: list[dict[str, float]] = []
    failures: list[str] = []
    for nx, nz in LADDER:
        measured = measure(HullParameters(nx=nx, nz=nz))
        row: dict[str, float] = {"nx": nx, "nz": nz}
        for name in QUANTITIES:
            row[name] = measured[name]
            row[f"{name}_relative_error"] = relative_error(measured[name], reference[name])
        row["displacement_volume"] = measured["displacement_volume"]
        row["displacement_volume_relative_error"] = relative_error(
            measured["displacement_volume"], reference["displacement_volume"]
        )
        rows.append(row)
        status = "  ".join(
            f"{name.split('_')[0][:4]}:{100.0 * row[f'{name}_relative_error']:+.3f}%"
            for name in QUANTITIES
        )
        print(f"nx={nx:4d} nz={nz:3d}  {status}")

    default_row = next(row for row in rows if (row["nx"], row["nz"]) == DEFAULT_GRID)
    fine_row = next(row for row in rows if (row["nx"], row["nz"]) == FINE_GRID)
    for name in QUANTITIES:
        default_error = abs(default_row[f"{name}_relative_error"])
        fine_error = abs(fine_row[f"{name}_relative_error"])
        print(f"{name:34s} default {100 * default_error:7.3f} %   fine {100 * fine_error:7.3f} %")
        if default_error > DEFAULT_CLAIM:
            failures.append(f"{name}: default-grid error {100 * default_error:.3f} % > 3 %")
        if fine_error > FINE_CLAIM:
            failures.append(f"{name}: fine-grid error {100 * fine_error:.3f} % > 0.5 %")
        if fine_error > default_error + 1e-12:
            failures.append(f"{name}: refinement did not reduce the error")

    stability = initial_stability_from_hydrostatics(generate_hull().hydrostatics, height_of_gravity=HEIGHT_OF_GRAVITY)
    if abs(stability.transverse_metacentric_height - reference["transverse_metacentric_height"]) > (
        DEFAULT_CLAIM * reference["transverse_metacentric_height"]
    ):
        failures.append("GM_T is outside the 3 % band")
    if not stability.longitudinal_metacentric_radius > stability.transverse_metacentric_radius:
        failures.append("BM_L must exceed BM_T on a slender hull")

    print(f"KB = {stability.height_of_centre_of_buoyancy:.6f} m, "
          f"BM_T = {stability.transverse_metacentric_radius:.6f} m, "
          f"KM_T = {stability.height_of_transverse_metacentre:.6f} m, "
          f"GM_T(KG={HEIGHT_OF_GRAVITY}) = {stability.transverse_metacentric_height:.6f} m")

    payload = {
        "schema": "tensorcad.stability.verification/1",
        "parameters": parameters.as_dict(),
        "height_of_gravity_m": HEIGHT_OF_GRAVITY,
        "reference": reference,
        "reference_notes": {
            "transverse_waterplane_inertia": "(2/3)*(B/2)**3*L*int sin(pi*xi)**2.1 dxi, gamma identity",
            "longitudinal_waterplane_inertia": "B*L**3*int sin(pi*xi)**0.7*(xi-1/2)**2 dxi, tanh-sinh",
            "height_of_metacentre_above_keel": "0.6*T + I_T/V",
        },
        "ladder": rows,
        "stability": stability.to_dict(),
        "claims": {
            "default_grid": [DEFAULT_GRID[0], DEFAULT_GRID[1]],
            "default_grid_max_relative_error": DEFAULT_CLAIM,
            "fine_grid": [FINE_GRID[0], FINE_GRID[1]],
            "fine_grid_max_relative_error": FINE_CLAIM,
            "failures": failures,
        },
    }
    (output / "results.json").write_text(json.dumps(payload, indent=2) + "\n")
    make_figure(parameters, rows, reference, output / "plot.png")

    lines = [
        "# Waterplane inertia and initial stability — verification",
        "",
        "Generated by `python examples/stability/run_hydrostatics_example.py`; every number below is",
        "read from `results.json` in the same run, none is typed by hand.",
        "",
        "## Reference",
        "| Quantity | Closed form |",
        "|---|---|",
        "| $I_T$ | $\\frac{2}{3}\\left(\\frac{B}{2}\\right)^3 L \\int_0^1 \\sin^{2.1}(\\pi\\xi)\\,d\\xi$ |",
        "| $I_L$ | $B L^3 \\int_0^1 \\sin^{0.7}(\\pi\\xi)\\,(\\xi-\\tfrac12)^2\\,d\\xi$ (tanh-sinh) |",
        "| $KB$ | $0.6T$ |",
        "| $BM_T$, $BM_L$ | $I_T/\\nabla$, $I_L/\\nabla$ |",
        "| $KM_T$ | $KB + BM_T$ |",
        "",
        "## Measured vs reference",
        "| nx | nz | $I_T$ | $I_L$ | $KM_T$ | $GM_T$ | vol |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            "| %d | %d | %+.3f %% | %+.3f %% | %+.3f %% | %+.3f %% | %+.3f %% |"  # noqa: UP031
            % (
                row["nx"],
                row["nz"],
                100 * row["transverse_waterplane_inertia_relative_error"],
                100 * row["longitudinal_waterplane_inertia_relative_error"],
                100 * row["height_of_metacentre_above_keel_relative_error"],
                100 * row["transverse_metacentric_height_relative_error"],
                100 * row["displacement_volume_relative_error"],
            )
        )
    lines += [
        "",
        f"KB = {stability.height_of_centre_of_buoyancy:.6f} m, "
        f"BM_T = {stability.transverse_metacentric_radius:.6f} m, "
        f"BM_L = {stability.longitudinal_metacentric_radius:.6f} m, "
        f"KM_T = {stability.height_of_transverse_metacentre:.6f} m, "
        f"GM_T (KG = {HEIGHT_OF_GRAVITY} m) = {stability.transverse_metacentric_height:.6f} m.",
        "",
        "## Verdict",
        "Every quantity is inside the 3 % band at the default 41x13 grid and inside",
        "0.5 % at 161x41: " + ("**PASS**" if not failures else "**FAIL** — " + "; ".join(failures)),
        "",
        "![error vs resolution](plot.png)",
        "",
    ]
    (output / "report.md").write_text("\n".join(lines))

    if failures:
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("all claims hold; wrote", output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
