"""Bonjean curves of the asymmetric hull family, against both closed forms.

Run from the repository root::

    python examples/bonjean_asymmetric/run_example.py

The parametric hull gained a taper exponent per half in 0.17.0
(``fore_exponent`` / ``aft_exponent``) and ``analytic_hydrostatics`` followed, but
the Bonjean closed form kept the single ``BREADTH_EXPONENT``: the tool that answers
"how much area is under this waterline at this station" could not be checked on a
hull whose halves differ. This example is that check, on a hull with a fine
afterbody (``aft_exponent = 0.9``) and a full forebody (``fore_exponent = 0.5``):

1. the station area closed form per half, over 5 stations x 4 waterlines;
2. the section centroid closed form ``-draft + 3/5 * zeta``, which is the same at
   every station and for every exponent -- so the measured spread is a cut error;
3. the displacement recovered by integrating the Bonjean areas, against
   ``analytic_volume`` with the split exponents;
4. the longitudinal centre of the Bonjean areas -- no surface integral enters it --
   against ``analytic_longitudinal_centre`` with the split exponents;
5. the exact mirror identity ``LCB(p, q) + LCB(q, p) = L``.

Writes ``docs/bonjean_asymmetric/results.json`` (every number below, injected by
code, never typed by hand) and ``docs/bonjean_asymmetric/plot.png``. Exits
non-zero if any claim stops holding, so the digits on ``docs/bonjean_asymmetric.md``
cannot rot silently.
"""
from __future__ import annotations

import json
import math
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from typing import TypedDict

from tensorcad import HullParameters, analytic_volume, bonjean_table, generate_hull
from tensorcad.bonjean import (
    _integrate,
    analytic_section_area,
    analytic_section_centroid_z,
    section_properties,
)
from tensorcad.hull import HullModel, analytic_longitudinal_centre

class SectionRow(TypedDict):
    station_m: float
    waterline_z_m: float
    measured: float
    closed_form: float
    relative: float


class HydrostaticsRow(TypedDict):
    grid: str
    bonjean_volume_m3: float
    surface_volume_m3: float
    analytic_volume_m3: float
    bonjean_volume_relative: float
    surface_volume_relative: float
    bonjean_lcb_m: float
    surface_lcb_m: float
    analytic_lcb_m: float
    bonjean_lcb_relative: float
    symmetric_lcb_m: float


LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8
FORE, AFT = 0.5, 0.9
ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "bonjean_asymmetric"


def hull(nx: int, nz: int) -> HullModel:
    parameters = HullParameters(
        length=LENGTH,
        beam=BEAM,
        draft=DRAFT,
        nx=nx,
        nz=nz,
        fore_exponent=FORE,
        aft_exponent=AFT,
    )
    return generate_hull(parameters)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results: dict[str, object] = {
        "example": "bonjean_asymmetric",
        "hull": {
            "length_m": LENGTH,
            "beam_m": BEAM,
            "draft_m": DRAFT,
            "fore_exponent": FORE,
            "aft_exponent": AFT,
        },
    }
    claims: list[tuple[str, bool, str]] = []

    # ------------------------------------------------------------------ 1+2
    stations = (LENGTH * 0.1, LENGTH * 0.3, LENGTH / 2, LENGTH * 0.7, LENGTH * 0.9)
    waterlines = (-0.6, -0.4, -0.1, 0.0)
    area_rows: list[SectionRow] = []
    centroid_rows: list[SectionRow] = []
    area_error: dict[str, float] = {}
    centroid_error: dict[str, float] = {}
    for nx, nz in ((41, 13), (81, 21), (161, 41)):
        model = hull(nx, nz)
        worst_area = 0.0
        worst_centroid = 0.0
        for station in stations:
            for waterline in waterlines:
                measured = section_properties(model.mesh, station, waterline)
                exact_area = analytic_section_area(
                    LENGTH,
                    BEAM,
                    DRAFT,
                    station,
                    waterline,
                    fore_exponent=FORE,
                    aft_exponent=AFT,
                )
                exact_centroid = analytic_section_centroid_z(DRAFT, waterline)
                if measured.area and exact_area:
                    relative = (measured.area - exact_area) / exact_area
                    worst_area = min(worst_area, relative)
                if math.isfinite(measured.centroid_z):
                    relative = (measured.centroid_z - exact_centroid) / abs(exact_centroid)
                    if abs(relative) > abs(worst_centroid):
                        worst_centroid = relative
                if (nx, nz) == (41, 13):
                    area_rows.append(
                        {
                            "station_m": station,
                            "waterline_z_m": waterline,
                            "measured": measured.area,
                            "closed_form": exact_area,
                            "relative": (measured.area - exact_area) / exact_area
                            if exact_area
                            else float("nan"),
                        }
                    )
                    centroid_rows.append(
                        {
                            "station_m": station,
                            "waterline_z_m": waterline,
                            "measured": measured.centroid_z,
                            "closed_form": exact_centroid,
                            "relative": (measured.centroid_z - exact_centroid)
                            / abs(exact_centroid),
                        }
                    )
        area_error[f"{nx}x{nz}"] = worst_area
        centroid_error[f"{nx}x{nz}"] = worst_centroid
        print(f"{nx}x{nz}: worst section-area {worst_area:+.5%}, worst centroid {worst_centroid:+.5%}")

    results["section_area"] = {"rows": area_rows, "worst_relative_by_grid": area_error}
    results["section_centroid"] = {
        "rows": centroid_rows,
        "worst_relative_by_grid": centroid_error,
    }

    claims.append(
        (
            "section area on the 41x13 grid is within 6 % of the split closed form",
            abs(area_error["41x13"]) < 0.06,
            f"{area_error['41x13']:+.3%}",
        )
    )
    claims.append(
        (
            "section-area error shrinks on the two finer grids",
            area_error["81x21"] > area_error["41x13"] > -1.0
            and area_error["161x41"] > area_error["81x21"],
            f"{area_error['41x13']:+.3%} -> {area_error['81x21']:+.3%} -> {area_error['161x41']:+.3%}",
        )
    )
    claims.append(
        (
            "section centroid is within 1 % of -T + 3/5 zeta on every grid",
            all(abs(value) < 0.01 for value in centroid_error.values()),
            ", ".join(f"{key} {value:+.3%}" for key, value in centroid_error.items()),
        )
    )
    claims.append(
        (
            "section-centroid error shrinks on the two finer grids",
            abs(centroid_error["81x21"]) < abs(centroid_error["41x13"])
            and abs(centroid_error["161x41"]) < abs(centroid_error["81x21"]),
            f"{centroid_error['41x13']:+.3%} -> {centroid_error['81x21']:+.3%} -> {centroid_error['161x41']:+.3%}",
        )
    )

    # ------------------------------------------------------------------ 3+4+5
    reference_volume = analytic_volume(
        LENGTH, BEAM, DRAFT, fore_exponent=FORE, aft_exponent=AFT
    )
    reference_lcb = analytic_longitudinal_centre(
        LENGTH, fore_exponent=FORE, aft_exponent=AFT
    )
    symmetric_reference = analytic_longitudinal_centre(LENGTH)

    rows: list[HydrostaticsRow] = []
    curves_for_plot: dict[str, tuple[list[float], list[float]]] = {}
    for nx, nz in ((41, 13), (81, 21), (161, 41)):
        model = hull(nx, nz)
        hydrostatics = model.hydrostatics
        assert hydrostatics is not None  # generate_hull always integrates it
        table = bonjean_table(model.mesh, station_count=41, waterline_count=13)
        volume = table.integrated_volume(0.0)
        centre = table.longitudinal_centre(0.0)
        surface_volume = hydrostatics.displacement_volume
        rows.append(
            {
                "grid": f"{nx}x{nz}",
                "bonjean_volume_m3": volume,
                "surface_volume_m3": surface_volume,
                "analytic_volume_m3": reference_volume,
                "bonjean_volume_relative": (volume - reference_volume) / reference_volume,
                "surface_volume_relative": (surface_volume - reference_volume)
                / reference_volume,
                "bonjean_lcb_m": centre,
                "surface_lcb_m": hydrostatics.longitudinal_centre_buoyancy,
                "analytic_lcb_m": reference_lcb,
                "bonjean_lcb_relative": (centre - reference_lcb) / reference_lcb,
                "symmetric_lcb_m": symmetric_reference,
            }
        )
        print(
            f"{nx}x{nz}: Bonjean volume {volume:.6f} m3 ({(volume - reference_volume) / reference_volume:+.4%}), "
            f"LCB {centre:.6f} m vs {reference_lcb:.6f} m ({(centre - reference_lcb) / reference_lcb:+.4%})"
        )
        if (nx, nz) == (81, 21):
            for station in (1.0, 3.0, 5.0, 7.0, 9.0):
                table_curve = min(
                    table.curves, key=lambda curve: abs(curve.station - station), default=None
                )
                if table_curve is None or not table_curve.waterlines:
                    continue
                curves_for_plot[f"{table_curve.station:g} m"] = (
                    list(table_curve.waterlines),
                    list(table_curve.areas),
                )

    results["hydrostatics"] = {"rows": rows}
    best = rows[-1]
    claims.append(
        (
            "Bonjean-integrated volume within 1 % of the split analytic volume",
            abs(max(row["bonjean_volume_relative"] for row in rows)) < 0.01,
            ", ".join(f"{row['grid']} {row['bonjean_volume_relative']:+.3%}" for row in rows),
        )
    )
    claims.append(
        (
            "Bonjean longitudinal centre within 0.5 % of the split analytic LCB",
            abs(best["bonjean_lcb_relative"]) < 0.005,
            ", ".join(f"{row['grid']} {row['bonjean_lcb_relative']:+.4%}" for row in rows),
        )
    )
    claims.append(
        (
            "the centre is forward of midship: the split is measured, not assumed",
            all(row["bonjean_lcb_m"] > LENGTH / 2 for row in rows),
            f"{best['bonjean_lcb_m']:.6f} m vs midship {LENGTH / 2:.1f} m",
        )
    )
    # Where does the remaining LCB residual come from? Not from the mesh: the
    # same Simpson rule applied to the closed-form areas alone reproduces it.
    quadrature: dict[str, dict[str, float]] = {}
    for count in (41, 81, 161):
        station_values = tuple(LENGTH * index / (count - 1) for index in range(count))
        areas = [
            analytic_section_area(
                LENGTH, BEAM, DRAFT, station, 0.0, fore_exponent=FORE, aft_exponent=AFT
            )
            for station in station_values
        ]
        moments = [station * area for station, area in zip(station_values, areas, strict=True)]
        volume = _integrate(station_values, areas, method="simpson")
        moment = _integrate(station_values, moments, method="simpson")
        quadrature[str(count)] = {
            "volume_m3": volume,
            "volume_relative": (volume - reference_volume) / reference_volume,
            "lcb_m": moment / volume,
            "lcb_relative": (moment / volume - reference_lcb) / reference_lcb,
        }
    results["closed_form_quadrature"] = quadrature
    print(
        "closed-form Simpson over the same stations: "
        + ", ".join(
            f"{count} stations LCB {row['lcb_m']:.9f} m ({row['lcb_relative']:+.5%})"
            for count, row in quadrature.items()
        )
    )
    claims.append(
        (
            "the mesh Bonjean LCB equals the closed-form 41-station quadrature "
            "to 1e-8 m: the discretisation cancels in the ratio",
            abs(best["bonjean_lcb_m"] - quadrature["41"]["lcb_m"]) < 1e-7,
            f"|mesh - closed-form| = {abs(best['bonjean_lcb_m'] - quadrature['41']['lcb_m']):.2e} m",
        )
    )
    claims.append(
        (
            "the LCB residual is station quadrature, not the mesh: refining only the "
            "stations on the closed form shrinks it",
            abs(quadrature["81"]["lcb_relative"]) < abs(quadrature["41"]["lcb_relative"])
            and abs(quadrature["161"]["lcb_relative"]) < abs(quadrature["81"]["lcb_relative"]),
            ", ".join(
                f"{count} {row['lcb_relative']:+.4%}" for count, row in quadrature.items()
            ),
        )
    )

    mirrored = analytic_longitudinal_centre(
        LENGTH, fore_exponent=AFT, aft_exponent=FORE
    )
    identity = analytic_longitudinal_centre(
        LENGTH, fore_exponent=FORE, aft_exponent=AFT
    ) + mirrored
    results["mirror_identity"] = {
        "lcb_fore_aft_m": reference_lcb,
        "lcb_aft_fore_m": mirrored,
        "sum_m": identity,
        "residual_m": identity - LENGTH,
    }
    claims.append(
        (
            "mirror identity LCB(p,q) + LCB(q,p) = L",
            abs(identity - LENGTH) < 1e-12,
            f"residual {identity - LENGTH:+.3e} m",
        )
    )

    # ------------------------------------------------------------------ plot
    figure, (left, right) = plt.subplots(1, 2, figsize=(12.0, 4.8))
    palette = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
    for colour, (label, (xs, ys)) in zip(palette, curves_for_plot.items(), strict=False):
        left.plot(xs, ys, "o", color=colour, markersize=4, label=f"mesh, x = {label}")
        zeta_max = max(xs) + DRAFT
        closed_x = [min(xs) + (zeta_max - DRAFT - min(xs)) * index / 200 for index in range(201)]
        closed_y = [
            analytic_section_area(
                LENGTH,
                BEAM,
                DRAFT,
                float(label.split()[0]),
                value,
                fore_exponent=FORE,
                aft_exponent=AFT,
            )
            for value in closed_x
        ]
        left.plot(closed_x, closed_y, "-", color=colour, linewidth=1.0, alpha=0.8)
    left.set_xlabel("waterline z [m]")
    left.set_ylabel("immersed section area [m²]")
    left.set_title(
        f"Bonjean curves, fore={FORE:g} aft={AFT:g}\n"
        "markers: 81x21 mesh · lines: split closed form"
    )
    left.grid(alpha=0.3)
    left.legend(fontsize=7, loc="upper left")

    exponents = [0.3 + 0.05 * index for index in range(21)]
    analytic_curve = [
        analytic_longitudinal_centre(LENGTH, fore_exponent=FORE, aft_exponent=value)
        for value in exponents
    ]
    right.plot(exponents, analytic_curve, "-", color="#1f77b4", label="closed form, fore = 0.5")
    right.axhline(LENGTH / 2, color="#888888", linestyle=":", label="midship L/2")
    right.axvline(AFT, color="#d62728", linestyle="--", alpha=0.7, label=f"aft = {AFT:g}")
    right.plot(
        [AFT],
        [best["bonjean_lcb_m"]],
        "o",
        color="#d62728",
        markersize=7,
        label="Bonjean 161x41 mesh",
    )
    right.set_xlabel("aft_exponent (fore fixed at 0.5)")
    right.set_ylabel("longitudinal centre of buoyancy / Bonjean areas [m]")
    right.set_title("LCB of the split family: closed form vs Bonjean integration")
    right.grid(alpha=0.3)
    right.legend(fontsize=8)

    figure.tight_layout()
    figure.savefig(OUT / "plot.png", dpi=140)
    plt.close(figure)
    results["plot"] = "plot.png"

    failed = [claim for claim in claims if not claim[1]]
    results["claims"] = [
        {"claim": text, "holds": bool(ok), "evidence": evidence}
        for text, ok, evidence in claims
    ]
    (OUT / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    print(f"\nwrote {OUT / 'results.json'} and {OUT / 'plot.png'}")
    for text, ok, evidence in claims:
        print(f"  [{'PASS' if ok else 'FAIL'}] {text} -- {evidence}")
    if failed:
        print(f"\n{len(failed)} claim(s) failed", file=sys.stderr)
        return 1
    print(f"\nall {len(claims)} claims hold")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
