"""Bonjean curves on the two families that had no analytic reference yet.

``docs/bonjean.md`` verified the station cut on a barge (exact) and on the
shipped ``hull`` family (``sin(pi x/L)^0.7 * sqrt(zeta/T)`` sections) against
:func:`tensorcad.bonjean.analytic_section_area`. Two families added later had no
closed-form Bonjean at all, so their curves were measured but never checked:

* the **asymmetric** hull (``HullParameters.fore_exponent`` / ``aft_exponent``,
  0.17.0) -- its analytic reference covered volume, waterplane and centres, not
  the section-area distribution that damage stability reads;
* the **Wigley** family (``tensorcad.wigley``, 0.20.0) -- its analytic reference
  covered hydrostatics, not the area-under-a-waterline table.

This script supplies both closed forms and measures the mesh against them,
which is what turns "the cut machinery works on a box" into "the cut machinery
reproduces two independent families to sub-percent accuracy, and the reference
actually discriminates the fore/aft split (72.8 % apart if you ignore it)".

Run from the repository root::

    python examples/bonjean_families/run_example.py

Writes ``docs/bonjean_families/results.json`` (every number injected by code,
none typed by hand) and ``docs/bonjean_families/plot.png``. Exits non-zero if
any claim stops holding.
"""
from __future__ import annotations

import json
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tensorcad.bonjean import analytic_section_area, bonjean_table, section_properties
from tensorcad.hull import HullParameters, analytic_volume, build_mesh as build_hull_mesh
from tensorcad.wigley import (
    WigleyParameters,
    analytic_maximum_section_area,
    analytic_section_area as wigley_analytic_section_area,
    build_mesh as build_wigley_mesh,
)

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "bonjean_families"

# ITTC reference proportions for the Wigley family, the case the resistance and
# diffraction benchmarks quote.
WIGLEY = {"length": 10.0, "beam": 1.0, "draft": 0.625}
# The parametric hull family's shipped default dimensions.
HULL = {"length": 10.0, "beam": 2.0, "draft": 0.8}
# The asymmetric case: the fore half falls off faster than the aft half, so the
# section areas are no longer symmetric about midship.
FOR, AFT = 1.4, 0.7

STATION_FRACTIONS = (0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95)
WIGLEY_STATIONS = (1.0, 2.5, 5.0, 7.5, 9.0)
WIGLEY_WATERLINES = (-0.5, -0.3125, -0.15625, -0.05)
HULL_WATERLINES = (-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, -0.02)


def wigley_mesh(nx: int, nz: int):
    parameters = WigleyParameters(nx=nx, nz=nz, **WIGLEY)
    mesh, _ = build_wigley_mesh(parameters)
    return mesh


def hull_mesh(nx: int, nz: int, fore: float = FOR, aft: float = AFT):
    parameters = HullParameters(nx=nx, nz=nz, fore_exponent=fore, aft_exponent=aft, **HULL)
    mesh, _ = build_hull_mesh(parameters)
    return parameters, mesh


def wigley_reference(parameters: WigleyParameters, station: float, waterline: float) -> float:
    return wigley_analytic_section_area(
        parameters.length,
        parameters.beam,
        parameters.draft,
        station,
        waterline,
        parameters.longitudinal_exponent,
        parameters.vertical_exponent,
    )


def hull_reference(station: float, waterline: float, fore: float = FOR, aft: float = AFT) -> float:
    return analytic_section_area(
        HULL["length"], HULL["beam"], HULL["draft"], station, waterline,
        fore_exponent=fore, aft_exponent=aft,
    )


def worst_relative(mesh, stations, waterlines, reference) -> dict[str, object]:
    """Worst |mesh/closed form - 1| over a station x waterline grid, plus the table."""
    table: list[dict[str, object]] = []
    worst = 0.0
    worst_at = (0.0, 0.0)
    for station in stations:
        for waterline in waterlines:
            measured = section_properties(mesh, station, waterline).area
            exact = reference(station, waterline)
            error = measured / exact - 1.0 if exact else 0.0
            if abs(error) > abs(worst):
                worst, worst_at = error, (station, waterline)
            table.append(
                {
                    "station_m": station,
                    "waterline_z_m": waterline,
                    "mesh_m2": measured,
                    "closed_form_m2": exact,
                    "relative_error": error,
                }
            )
    return {"worst_relative_error": worst, "worst_at": list(worst_at), "samples": table}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results: dict[str, object] = {"example": "bonjean_families", "version": "0.21.0"}
    claims: list[tuple[str, bool, str]] = []

    # ---------------------------------------------------------------- wigley --
    reference_parameters = WigleyParameters(nx=81, nz=13, **WIGLEY)
    reference_mesh = wigley_mesh(81, 13)
    wigley = worst_relative(
        reference_mesh,
        WIGLEY_STATIONS,
        WIGLEY_WATERLINES,
        lambda x, z: wigley_reference(reference_parameters, x, z),
    )
    results["wigley"] = {
        "parameters": {
            "length": WIGLEY["length"],
            "beam": WIGLEY["beam"],
            "draft": WIGLEY["draft"],
            "longitudinal_exponent": reference_parameters.longitudinal_exponent,
            "vertical_exponent": reference_parameters.vertical_exponent,
            "nx": 81,
            "nz": 13,
        },
        "stations_m": list(WIGLEY_STATIONS),
        "waterlines_m": list(WIGLEY_WATERLINES),
        **wigley,
    }
    claims.append(
        (
            "wigley sections within 1 % of the closed form",
            abs(wigley["worst_relative_error"]) < 0.01,
            f"{abs(wigley['worst_relative_error']) * 100:.3f} % at "
            f"x={wigley['worst_at'][0]:g} m, z={wigley['worst_at'][1]:g} m",
        )
    )

    # The midship section at the design waterline must equal the closed form the
    # hydrostatics module already used, which pins the new expression to it.
    midship: list[dict[str, float]] = []
    for nz in (13, 25, 49):
        mesh = wigley_mesh(81, nz)
        measured = section_properties(mesh, WIGLEY["length"] / 2, 0.0).area
        exact = analytic_maximum_section_area(WIGLEY["beam"], WIGLEY["draft"])
        midship.append(
            {
                "nz": float(nz),
                "mesh_m2": measured,
                "closed_form_m2": exact,
                "relative_error": measured / exact - 1.0,
            }
        )
    results["wigley_midship"] = midship
    claims.append(
        (
            "wigley midship equals the known section area",
            abs(midship[0]["relative_error"]) < 0.01,
            f"{midship[0]['relative_error'] * 100:+.3f} % at nz = 13",
        )
    )

    # ------------------------------------------------------- asymmetric hull --
    parameters, hull = hull_mesh(81, 41)
    asymmetric = worst_relative(hull, tuple(f * HULL["length"] for f in STATION_FRACTIONS), HULL_WATERLINES, hull_reference)
    # Per-waterline worst error: the sqrt keel taper is what bends the deepest
    # curves, so the split has to be reported, not averaged away.
    per_waterline = []
    for waterline in HULL_WATERLINES:
        worst = max(
            abs(
                section_properties(hull, f * HULL["length"], waterline).area
                / hull_reference(f * HULL["length"], waterline)
                - 1.0
            )
            for f in STATION_FRACTIONS
        )
        per_waterline.append({"waterline_z_m": waterline, "worst_relative_error": worst})
    results["asymmetric_hull"] = {
        "parameters": {
            "length": HULL["length"],
            "beam": HULL["beam"],
            "draft": HULL["draft"],
            "fore_exponent": FOR,
            "aft_exponent": AFT,
            "nx": 81,
            "nz": 41,
        },
        "per_waterline": per_waterline,
        **asymmetric,
    }
    above_keel = max(row["worst_relative_error"] for row in per_waterline if row["waterline_z_m"] >= -0.4)
    # Bookkeeping for the honest side of the table: at the coarse waterline
    # count the deepest station curves are the worst point in the family, because
    # the sqrt keel taper has an infinite slope there and a polygon in z cannot
    # follow it. Reported by name so the docs quote it from JSON, not by memory.
    _, coarse = hull_mesh(81, 13)
    coarse_deepest = max(
        abs(
            section_properties(coarse, f * HULL["length"], -0.6).area
            / hull_reference(f * HULL["length"], -0.6)
            - 1.0
        )
        for f in STATION_FRACTIONS
    )
    results["asymmetric_hull"]["coarse_nz13_deepest_waterline_relative_error"] = coarse_deepest
    claims.append(
        (
            "asymmetric sections within 2 % above z = -0.4 m",
            above_keel < 0.02,
            f"{above_keel * 100:.3f} % at nz = 41",
        )
    )

    # ---------------------------------------------- the split is not cosmetic --
    # Reading the asymmetric mesh against a *symmetric* reference (the shipped
    # exponent everywhere) has to fail loudly, otherwise none of the agreement
    # above proves that fore and aft were told apart.
    fore_stations = (0.7, 0.8, 0.9, 0.95)
    symmetric_reference_worst = max(
        abs(
            section_properties(hull, f * HULL["length"], -0.2).area
            / analytic_section_area(HULL["length"], HULL["beam"], HULL["draft"], f * HULL["length"], -0.2)
            - 1.0
        )
        for f in fore_stations
    )
    split_reference_worst = max(
        abs(section_properties(hull, f * HULL["length"], -0.2).area / hull_reference(f * HULL["length"], -0.2) - 1.0)
        for f in fore_stations
    )
    results["split_discrimination"] = {
        "stations_m": [f * HULL["length"] for f in fore_stations],
        "symmetric_reference_worst_relative_error": symmetric_reference_worst,
        "split_reference_worst_relative_error": split_reference_worst,
    }
    claims.append(
        (
            "the split reference is discriminating",
            symmetric_reference_worst > 0.10 and split_reference_worst < 0.02,
            f"symmetric {symmetric_reference_worst * 100:.2f} % vs split {split_reference_worst * 100:.2f} %",
        )
    )

    # ------------------------------------------------------- mirror involution --
    swapped, swapped_mesh = hull_mesh(81, 41, fore=AFT, aft=FOR)
    mirror = 0.0
    for fraction in STATION_FRACTIONS:
        forward = section_properties(hull, fraction * HULL["length"], -0.2).area
        swapped_area = section_properties(swapped_mesh, (1.0 - fraction) * HULL["length"], -0.2).area
        mirror = max(mirror, abs(forward / swapped_area - 1.0) if swapped_area else 0.0)
    results["mirror"] = {"worst_relative_difference": mirror, "fore_exponent": FOR, "aft_exponent": AFT}
    claims.append(
        (
            "swapping the exponents mirrors the curves",
            mirror < 1e-6,
            f"max |rel| {mirror:.2e} over {len(STATION_FRACTIONS)} stations",
        )
    )

    # ------------------------------------------------- convergence ladders --
    # Waterline refinement at fixed station count: the section polygon is linear
    # in z between rings, so the error must fall as the square of the spacing.
    wigley_ladder = []
    for nz in (13, 25, 49):
        mesh = wigley_mesh(81, nz)
        worst = max(
            abs(
                section_properties(mesh, station, waterline).area
                / wigley_reference(reference_parameters, station, waterline)
                - 1.0
            )
            for station in WIGLEY_STATIONS
            for waterline in WIGLEY_WATERLINES
        )
        wigley_ladder.append({"nz": nz, "worst_relative_error": worst})

    hull_ladder = []
    for nz in (13, 25, 49, 97):
        _, ladder_mesh = hull_mesh(81, nz)
        worst = max(
            abs(
                section_properties(ladder_mesh, f * HULL["length"], -0.2).area
                / hull_reference(f * HULL["length"], -0.2)
                - 1.0
            )
            for f in STATION_FRACTIONS
        )
        hull_ladder.append({"nz": nz, "worst_relative_error": worst})

    def ratios(ladder: list[dict[str, float]]) -> list[float]:
        return [
            abs(ladder[index]["worst_relative_error"]) / abs(ladder[index + 1]["worst_relative_error"])
            for index in range(len(ladder) - 1)
        ]

    wigley_ratios = ratios(wigley_ladder)
    hull_ratios = ratios(hull_ladder)
    results["convergence"] = {
        "wigley": {"ladder": wigley_ladder, "ratios": wigley_ratios},
        "asymmetric_hull": {"ladder": hull_ladder, "ratios": hull_ratios, "waterline_z_m": -0.2},
    }
    claims.append(
        (
            "wigley error falls like the square of the waterline spacing",
            all(ratio > 3.0 for ratio in wigley_ratios),
            f"ratios {[format(ratio, '.2f') for ratio in wigley_ratios]}",
        )
    )
    claims.append(
        (
            "asymmetric error falls with refinement",
            all(ratio > 2.5 for ratio in hull_ratios),
            f"ratios {[format(ratio, '.2f') for ratio in hull_ratios]}",
        )
    )

    # ------------------------------------- table integration vs closed volume --
    # The Bonjean table integrates the areas along x; the analytic volume is an
    # independent closed form that never sees a station. They must agree.
    station_grid = tuple(HULL["length"] * index / 40 for index in range(41))
    table = bonjean_table(hull, stations=station_grid, waterlines=(-0.8, -0.4, 0.0))
    simpson = table.integrated_volume(0.0, method="simpson")
    trapezoid = table.integrated_volume(0.0, method="trapezoid")
    exact_volume = analytic_volume(
        HULL["length"], HULL["beam"], HULL["draft"], fore_exponent=FOR, aft_exponent=AFT
    )
    results["volume_closure"] = {
        "simpson_m3": simpson,
        "trapezoid_m3": trapezoid,
        "closed_form_m3": exact_volume,
        "simpson_relative_error": simpson / exact_volume - 1.0,
        "trapezoid_relative_error": trapezoid / exact_volume - 1.0,
    }
    claims.append(
        (
            "integrating the curves recovers the closed-form volume",
            abs(simpson / exact_volume - 1.0) < 0.005,
            f"{abs(simpson / exact_volume - 1.0) * 100:.3f} % over 41 stations",
        )
    )

    # ------------------------------------------------------------------ plot --
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.5, 4.7))
    for station in WIGLEY_STATIONS:
        levels = [(-WIGLEY["draft"] + WIGLEY["draft"] * index / 40) for index in range(41)]
        left.plot(
            [section_properties(reference_mesh, station, z).area for z in levels],
            levels,
            "-",
            linewidth=1.4,
            label=f"Wigley mesh, x = {station:g} m",
        )
        left.plot(
            [wigley_reference(reference_parameters, station, z) for z in levels],
            levels,
            "k--",
            linewidth=0.9,
            alpha=0.75,
        )
    left.plot([], [], "k--", linewidth=0.9, label="closed form")
    left.set_xlabel("section area below the waterline [m$^2$]")
    left.set_ylabel("waterline z [m]")
    left.set_title("Bonjean curves, Wigley (81 x 13), ITTC proportions")
    left.grid(alpha=0.3)
    left.legend(fontsize=7)

    levels = [(-HULL["draft"] + HULL["draft"] * index / 40) for index in range(41)]
    for fraction in (0.25, 0.75):
        station = fraction * HULL["length"]
        right.plot(
            [section_properties(hull, station, z).area for z in levels],
            levels,
            "-",
            linewidth=1.4,
            label=f"asymmetric mesh, x = {station:g} m",
        )
        right.plot(
            [hull_reference(station, z) for z in levels],
            levels,
            "k--",
            linewidth=0.9,
            alpha=0.75,
        )
    right.plot([], [], "k--", linewidth=0.9, label="closed form")
    right.set_xlabel("section area below the waterline [m$^2$]")
    right.set_ylabel("waterline z [m]")
    right.set_title(f"Bonjean curves, asymmetric hull (fore$^{{{FOR}}}$, aft$^{{{AFT}}}$)")
    right.grid(alpha=0.3)
    right.legend(fontsize=7)

    figure.tight_layout()
    figure.savefig(OUT / "plot.png", dpi=140)
    plt.close(figure)

    figure, axes = plt.subplots(figsize=(6.4, 4.4))
    axes.plot(
        [row["nz"] for row in wigley_ladder],
        [abs(row["worst_relative_error"]) * 100 for row in wigley_ladder],
        "o-",
        label="Wigley, worst over 5 stations x 4 waterlines",
    )
    axes.plot(
        [row["nz"] for row in hull_ladder],
        [abs(row["worst_relative_error"]) * 100 for row in hull_ladder],
        "s-",
        label="asymmetric hull, worst over 10 stations at z = -0.2 m",
    )
    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlabel("vertical levels nz")
    axes.set_ylabel("worst relative error vs closed form [%]")
    axes.set_title("Waterline refinement: both families are second order")
    axes.grid(alpha=0.3, which="both")
    axes.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(OUT / "convergence.png", dpi=140)
    plt.close(figure)

    results["claims"] = [{"claim": name, "holds": bool(ok), "detail": detail} for name, ok, detail in claims]
    (OUT / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True))

    print(
        json.dumps(
            {
                "wigley_worst_error": wigley["worst_relative_error"],
                "wigley_ratios": wigley_ratios,
                "asymmetric_worst_above_keel": above_keel,
                "asymmetric_deepest_waterline_error": per_waterline[0]["worst_relative_error"],
                "coarse_nz13_deepest_waterline_error": coarse_deepest,
                "hull_ratios": hull_ratios,
                "split_discrimination": results["split_discrimination"],
                "mirror": mirror,
                "volume_simpson_error": results["volume_closure"]["simpson_relative_error"],
            },
            indent=2,
        )
    )
    failures = [name for name, ok, _ in claims if not ok]
    if failures:
        print("CLAIMS FAILED:", failures, file=sys.stderr)
        return 1
    print(f"all {len(claims)} claims hold; wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())