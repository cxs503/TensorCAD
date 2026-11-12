"""Bonjean curves: measured, compared with the closed form, and plotted.

Run from the repository root::

    python examples/bonjean/run_example.py

Writes ``docs/bonjean/results.json`` (every number below, injected by code, not
typed by hand) and ``docs/bonjean/plot.png``. Exits non-zero if any claim stops
holding, so the numbers on ``docs/bonjean.md`` cannot rot silently.
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_REPO_ROOT = _Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in _sys.path:  # allow running from a checkout without installing
    _sys.path.insert(0, str(_REPO_ROOT))

import json
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tensorcad import (
    HullParameters,
    analytic_volume,
    bonjean_table,
    compute_hydrostatics,
    generate_hull,
    section_properties,
)
from tensorcad import mesh as meshlib
from tensorcad.bonjean import analytic_section_area
from tensorcad.geometry import Mesh, MeshBuilder

LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8
ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "bonjean"


def box(minimum, maximum) -> Mesh:
    x0, y0, z0 = minimum
    x1, y1, z1 = maximum
    builder = MeshBuilder()
    corners = [
        (x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0),
        (x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1),
    ]
    ids = [builder.add_vertex(point) for point in corners]
    for face in (
        (0, 3, 2), (0, 2, 1), (4, 5, 6), (4, 6, 7), (0, 1, 5), (0, 5, 4),
        (2, 3, 7), (2, 7, 6), (1, 2, 6), (1, 6, 5), (0, 4, 7), (0, 7, 3),
    ):
        builder.add_triangle(*[ids[index] for index in face])
    return meshlib.orient_outward(builder.build()).mesh


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results: dict[str, object] = {"example": "bonjean", "version": "0.13.0"}
    claims: list[tuple[str, bool, str]] = []

    # ---------------------------------------------------------------- barge --
    barge = box((0.0, -BEAM / 2, -DRAFT), (LENGTH, BEAM / 2, 0.0))
    barge_table = bonjean_table(barge, station_count=11, waterline_count=5)
    barge_volume = barge_table.integrated_volume(0.0)
    exact_volume = LENGTH * BEAM * DRAFT
    results["barge"] = {
        "section_area_m2": section_properties(barge, LENGTH / 2).area,
        "section_centroid_z_m": section_properties(barge, LENGTH / 2).centroid_z,
        "integrated_volume_m3": barge_volume,
        "exact_volume_m3": exact_volume,
        "relative_error": barge_volume / exact_volume - 1.0,
    }
    claims.append(("barge area exact", abs(results["barge"]["section_area_m2"] - BEAM * DRAFT) < 1e-12, ""))
    claims.append(("barge volume exact", abs(barge_volume / exact_volume - 1.0) < 1e-12, ""))

    # ----------------------------------------------------------------- hull --
    hull = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=41, nz=13))
    stations = (1.0, 2.5, 5.0, 7.5, 9.0)
    waterlines = (-0.6, -0.4, -0.2, 0.0)
    table: dict[str, object] = {"stations_m": list(stations), "waterlines_m": list(waterlines), "rows": []}
    worst = 0.0
    for station in stations:
        row = {"station_m": station, "areas": []}
        for waterline in waterlines:
            measured = section_properties(hull.mesh, station, waterline).area
            reference = analytic_section_area(LENGTH, BEAM, DRAFT, station, waterline)
            error = measured / reference - 1.0
            worst = max(worst, abs(error))
            row["areas"].append(
                {
                    "waterline_z_m": waterline,
                    "mesh_m2": measured,
                    "closed_form_m2": reference,
                    "relative_error": error,
                }
            )
        table["rows"].append(row)
    table["worst_relative_error"] = worst
    results["section_table"] = table
    claims.append(("hull section within 6 %", worst < 0.06, f"{worst * 100:.2f} %"))

    # ------------------------------------------------- vertical convergence --
    convergence = []
    for nz in (13, 21, 41, 81):
        mesh = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=21, nz=nz)).mesh
        measured = section_properties(mesh, LENGTH / 2, -0.6).area
        reference = analytic_section_area(LENGTH, BEAM, DRAFT, LENGTH / 2, -0.6)
        convergence.append({"nz": nz, "relative_error": measured / reference - 1.0})
    results["vertical_convergence"] = convergence
    ratios = [
        abs(convergence[index]["relative_error"]) / abs(convergence[index + 1]["relative_error"])
        for index in range(len(convergence) - 1)
    ]
    results["vertical_convergence_ratio"] = ratios
    claims.append(("vertical convergence monotone", all(r > 1.5 for r in ratios), f"ratios {[format(r, '.2f') for r in ratios]}"))

    # ------------------------------------------- two integrations, one body --
    volume_closure = []
    analytic = analytic_volume(LENGTH, BEAM, DRAFT)
    for nx in (11, 21, 41, 81):
        mesh = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=nx, nz=41)).mesh
        station_grid = tuple(LENGTH * index / (nx - 1) for index in range(nx))
        grid_table = bonjean_table(mesh, stations=station_grid, waterline_count=3)
        simpson = grid_table.integrated_volume(0.0, method="simpson")
        trapezoid = grid_table.integrated_volume(0.0, method="trapezoid")
        flux = compute_hydrostatics(mesh, 0.0).displacement_volume
        volume_closure.append(
            {
                "nx": nx,
                "nz": 41,
                "simpson_m3": simpson,
                "trapezoid_m3": trapezoid,
                "surface_flux_m3": flux,
                "analytic_m3": analytic,
                "simpson_error": simpson / analytic - 1.0,
                "trapezoid_error": trapezoid / analytic - 1.0,
                "flux_error": flux / analytic - 1.0,
                "sections_minus_flux": trapezoid / flux - 1.0,
            }
        )
    results["volume_closure"] = volume_closure
    identity = max(abs(row["sections_minus_flux"]) for row in volume_closure)
    claims.append(("sections == surface flux", identity < 1e-12, f"max |rel| {identity:.2e}"))

    # ----------------------------------------------------------------- plot --
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.6))
    for station in stations:
        curve = [
            (section_properties(hull.mesh, station, waterline).area, waterline)
            for waterline in [(-DRAFT + DRAFT * index / 40) for index in range(41)]
        ]
        areas = [area for area, _ in curve]
        levels = [waterline for _, waterline in curve]
        left.plot(areas, levels, "-", linewidth=1.4, label=f"mesh x = {station:g} m")
        reference = [analytic_section_area(LENGTH, BEAM, DRAFT, station, z) for z in levels]
        left.plot(reference, levels, "k--", linewidth=0.9, alpha=0.75)
    left.plot([], [], "k--", linewidth=0.9, label="closed form")
    left.set_xlabel("section area below the waterline [m$^2$]")
    left.set_ylabel("waterline z [m]")
    left.set_title("Bonjean curves, parametric hull (41 x 13)")
    left.grid(alpha=0.3)
    left.legend(fontsize=7)

    right.plot([row["nz"] for row in convergence], [abs(row["relative_error"]) * 100 for row in convergence],
               "o-", label="section area at z = -0.6 m")
    right.plot([row["nx"] for row in volume_closure], [abs(row["simpson_error"]) * 100 for row in volume_closure],
               "s-", label="displacement (Simpson over stations)")
    right.set_xscale("log")
    right.set_yscale("log")
    right.set_xlabel("grid resolution (nz for areas, nx for volume)")
    right.set_ylabel("relative error vs closed form [%]")
    right.set_title("Discretisation error against the analytic reference")
    right.grid(alpha=0.3, which="both")
    right.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(OUT / "plot.png", dpi=140)
    plt.close(figure)

    results["claims"] = [{"claim": name, "holds": bool(ok), "detail": detail} for name, ok, detail in claims]
    (OUT / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True))

    print(json.dumps(
        {
            "barge_volume_error": results["barge"]["relative_error"],
            "worst_section_error": worst,
            "vertical_convergence_ratio": ratios,
            "identity_residual": identity,
            "volume_error_at_nx41": volume_closure[2]["simpson_error"],
        },
        indent=2,
    ))
    failures = [name for name, ok, _ in claims if not ok]
    if failures:
        print("CLAIMS FAILED:", failures, file=sys.stderr)
        return 1
    print(f"all {len(claims)} claims hold; wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
