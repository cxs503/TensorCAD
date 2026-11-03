"""Verify the sampled hull's hydrostatics against the closed-form reference.

Runs a resolution ladder, compares every reported hydrostatic quantity with
``tensorcad.analytic_hydrostatics``, fits the observed convergence order, and
cross-checks the wetted-area reference with a Richardson extrapolation of the mesh
series (a reference the quadrature never saw). Writes JSON plus figures and exits
non-zero if any claim stops holding.

    python examples/hydrostatics_verification/run_verification.py
    python examples/hydrostatics_verification/run_verification.py --output artifacts/hydrostatics_verification
"""
from __future__ import annotations

import argparse
from itertools import pairwise
import json
import math
import sys
from pathlib import Path

from tensorcad import HullParameters, analytic_hydrostatics, generate_hull

#: Resolution ladder: (nx, nz). Both axes refine together, so ``h`` below is read as
#: the longitudinal spacing ``length/(nx-1)`` and the fitted slopes are order
#: estimates, not per-axis convergence rates.
LADDER = ((21, 9), (41, 13), (81, 21), (161, 41), (321, 81))

#: Quantities compared, in report order.
QUANTITIES = (
    "displacement_volume",
    "displacement_mass",
    "waterplane_area",
    "wetted_surface_area",
    "maximum_section_area",
    "longitudinal_centre_buoyancy",
    "longitudinal_centre_flotation",
    "vertical_centre_buoyancy",
    "block_coefficient",
    "prismatic_coefficient",
)

#: Acceptance bands, applied to the relative error of every quantity.
DEFAULT_GRID_CLAIM = 0.03
FINE_GRID_CLAIM = 0.005
RICHARDSON_CLAIM = 2e-3
QUADRATURE_CLAIM = 1e-9


def _grid_index(ladder: tuple[tuple[int, int], ...], target: tuple[int, int], fallback: int) -> int:
    """Position of ``target`` in the ladder, or ``fallback`` when it is absent."""
    return ladder.index(target) if target in ladder else min(fallback, len(ladder) - 1)


def _h(parameters: HullParameters) -> float:
    """Nominal cell size of a grid: longitudinal spacing, the coarser of the two axes."""
    return parameters.length / (parameters.nx - 1)


def _order(h_fine: float, h_coarse: float, error_fine: float, error_coarse: float) -> float:
    """Observed convergence order from two refinement levels on a log-log line."""
    if error_fine <= 0.0 or error_coarse <= 0.0 or h_fine >= h_coarse:
        return float("nan")
    return math.log(error_coarse / error_fine) / math.log(h_coarse / h_fine)


def run_ladder(reference) -> list[dict]:
    rows = []
    for nx, nz in LADDER:
        parameters = HullParameters(nx=nx, nz=nz)
        model = generate_hull(parameters)
        comparison = reference.compare(model.hydrostatics, quantities=list(QUANTITIES))
        rows.append(
            {
                "nx": nx,
                "nz": nz,
                "vertices": len(model.mesh.vertices),
                "triangles": len(model.mesh.faces),
                "h_m": _h(parameters),
                "worst_relative_error": comparison["worst_relative_error"],
                "rows": comparison["rows"],
            }
        )
    return rows


def convergence_summary(rows: list[dict]) -> dict:
    """Per-quantity error decay: observed order, final error, monotone or not."""
    summary = {}
    for name in QUANTITIES:
        errors = [abs(row["rows"][name]["relative_error"]) for row in rows]
        hs = [row["h_m"] for row in rows]
        orders = [
            _order(hs[i + 1], hs[i], errors[i + 1], errors[i]) for i in range(len(rows) - 1)
        ]
        finite_orders = [order for order in orders if math.isfinite(order)]
        summary[name] = {
            "relative_errors": errors,
            "observed_orders": orders,
            "median_observed_order": sorted(finite_orders)[len(finite_orders) // 2]
            if finite_orders
            else float("nan"),
            "monotone": all(
                later < earlier for earlier, later in pairwise(errors)
            ),
            # The two longitudinal centres are fixed at L/2 by symmetry, so the mesh
            # reproduces them to round-off and their error cannot decay monotonically:
            # it is already at the float64 floor.
            "at_round_off_floor": max(errors) < 1e-12,
            "error_default_grid": errors[_grid_index(LADDER, (41, 13), 1)],
            "error_fine_grid": errors[_grid_index(LADDER, (161, 41), len(errors) - 1)],
        }
    return summary


def richardson_wetted_area(rows: list[dict]) -> dict:
    """Extrapolate the mesh wetted-area series and compare with the quadrature value.

    An inscribed surface never over-estimates area, so the mesh series is a
    decreasing sequence bounded by the truth. The last two steps give the order;
    the extrapolation is the third estimate of the same number.
    """
    areas = [row["rows"]["wetted_surface_area"]["measured"] for row in rows]
    reference = rows[0]["rows"]["wetted_surface_area"]["reference"]
    previous_step = areas[-2] - areas[-3]
    last_step = areas[-1] - areas[-2]
    order = math.log2(previous_step / last_step) if last_step else float("nan")
    extrapolated = areas[-1] + last_step / (2.0**order - 1.0) if order > 1.0 else float("nan")
    return {
        "mesh_areas": areas,
        "quadrature_reference": reference,
        "observed_order": order,
        "extrapolated": extrapolated,
        "extrapolation_relative_difference": extrapolated / reference - 1.0,
    }


def write_figures(rows: list[dict], summary: dict, output: Path) -> list[str]:
    """PNG figures, or an empty list when matplotlib is not installed.

    The library itself is standard-library only and its benchmark module writes SVG
    by hand for that reason; the figures here are a report artifact, so they are
    optional rather than a new dependency of the example.
    """
    try:
        import matplotlib
    except ImportError:
        return []

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    written = []

    # Figure 1: error decay against the nominal cell size, one line per quantity.
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    axes = list(axes)
    hs = [row["h_m"] for row in rows]
    for name in QUANTITIES:
        errors = [100.0 * value for value in summary[name]["relative_errors"]]
        axes[0].loglog(hs, errors, marker="o", linewidth=1.2, label=name.replace("_", " "))
    axes[0].set_ylim(1e-16, 5.0)
    axes[0].axhline(1.0, color="0.5", linestyle="--", linewidth=0.8)
    axes[0].axhline(3.0, color="0.5", linestyle=":", linewidth=0.8)
    axes[0].set_xlabel("nominal cell size h = L/(nx-1)  [m]")
    axes[0].set_ylabel("|mesh - closed form| / closed form  [%]")
    axes[0].set_title("Hydrostatic error against the closed-form reference")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend(fontsize=5.5, ncol=2, loc="lower left")

    reference = rows[0]["rows"]["wetted_surface_area"]["reference"]
    areas = [row["rows"]["wetted_surface_area"]["measured"] for row in rows]
    axes[1].semilogx(hs, [100.0 * (area / reference - 1.0) for area in areas], marker="s",
                     color="tab:blue", label="mesh wetted area")
    for name, colour in (("displacement_volume", "tab:orange"), ("vertical_centre_buoyancy", "tab:green")):
        points = [100.0 * row["rows"][name]["relative_error"] for row in rows]
        axes[1].semilogx(hs, points, marker="^", linewidth=1.0, color=colour,
                         label=name.replace("_", " "))
    axes[1].axhline(0.0, color="0.3", linewidth=0.8)
    axes[1].set_xlabel("nominal cell size h  [m]")
    axes[1].set_ylabel("signed error  [%]")
    axes[1].set_title("Every quantity is under-integrated, and converges")
    axes[1].grid(True, which="both", alpha=0.25)
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    path = output / "convergence.png"
    fig.savefig(path, dpi=110)
    plt.close(fig)
    written.append(path.name)

    # Figure 2: the offset field the reference describes, with the sampled stations on it.
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.0))
    axes = list(axes)
    parameters = HullParameters(nx=161, nz=41)
    length, beam, draft = parameters.length, parameters.beam, parameters.draft
    stations = 121
    levels = 61
    grid_x = [length * index / (stations - 1) for index in range(stations)]
    grid_z = [-draft * (1.0 - level / (levels - 1)) for level in range(levels)]
    field = [
        [0.5 * beam * math.sin(math.pi * x / length) ** 0.7 * math.sqrt((z + draft) / draft)
         for x in grid_x]
        for z in grid_z
    ]
    contour = axes[0].contourf(grid_x, grid_z, field, levels=24, cmap="viridis")
    fig.colorbar(contour, ax=axes[0], label="half-breadth b(x, z)  [m]")
    model = generate_hull(HullParameters(nx=21, nz=9))
    sampled_x = [point[0] for point in model.mesh.vertices]
    sampled_z = [point[2] for point in model.mesh.vertices]
    axes[0].plot(sampled_x, sampled_z, ".", color="white", markersize=1.5, alpha=0.65,
                 label="sampled grid vertices (21x9)")
    axes[0].set_xlabel("x  [m]  (bow at x = L)")
    axes[0].set_ylabel("z  [m]  (waterline z = 0)")
    axes[0].set_title("Analytic offset field b(x,z) the reference integrates")
    axes[0].legend(fontsize=7, loc="upper left")

    default = rows[_grid_index(LADDER, (41, 13), 1)]
    fine = rows[_grid_index(LADDER, (161, 41), len(rows) - 1)]
    axes[1].bar(
        [name.replace("_", "\n") for name in QUANTITIES],
        [100.0 * default["rows"][name]["relative_error"] for name in QUANTITIES],
        color="tab:blue", label=f"{default['nx']}x{default['nz']} (default)",
    )
    axes[1].bar(
        [name.replace("_", "\n") for name in QUANTITIES],
        [100.0 * fine["rows"][name]["relative_error"] for name in QUANTITIES],
        color="tab:green", label=f"{fine['nx']}x{fine['nz']}",
    )
    axes[1].axhline(-1.0, color="0.4", linestyle="--", linewidth=0.8)
    axes[1].axhline(-3.0, color="tab:red", linestyle=":", linewidth=0.9, label="-3 % acceptance")
    axes[1].set_ylabel("signed relative error  [%]")
    axes[1].set_title("Per-quantity error against the closed form")
    axes[1].tick_params(axis="x", labelsize=5.5)
    axes[1].grid(True, axis="y", alpha=0.25)
    axes[1].legend(fontsize=7)
    fig.tight_layout()
    path = output / "offset_field_and_error.png"
    fig.savefig(path, dpi=110)
    plt.close(fig)
    written.append(path.name)
    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("artifacts/hydrostatics_verification"))
    args = parser.parse_args(argv)
    args.output.mkdir(parents=True, exist_ok=True)

    reference = analytic_hydrostatics()
    rows = run_ladder(reference)
    summary = convergence_summary(rows)
    richardson = richardson_wetted_area(rows)

    default_index = _grid_index(LADDER, (41, 13), 1)
    fine_index = _grid_index(LADDER, (161, 41), 3)
    default_worst = rows[default_index]["worst_relative_error"]
    fine_worst = rows[fine_index]["worst_relative_error"]
    quadrature_error = reference.quadrature["wetted_area_relative_difference"]
    round_off = [name for name in QUANTITIES if summary[name]["at_round_off_floor"]]
    monotone = [
        name
        for name in QUANTITIES
        if not summary[name]["monotone"] and name not in round_off
    ]

    claims = {
        "default_grid_within_3_percent": {
            "value": default_worst,
            "claim": DEFAULT_GRID_CLAIM,
            "ok": default_worst < DEFAULT_GRID_CLAIM,
            "detail": f"{rows[default_index]['nx']}x{rows[default_index]['nz']} grid, worst quantity",
        },
        "fine_grid_within_half_a_percent": {
            "value": fine_worst,
            "claim": FINE_GRID_CLAIM,
            "ok": fine_worst < FINE_GRID_CLAIM,
            "detail": f"{rows[fine_index]['nx']}x{rows[fine_index]['nz']} grid, worst quantity",
        },
        "richardson_confirms_wetted_area": {
            "value": richardson["extrapolation_relative_difference"],
            "claim": RICHARDSON_CLAIM,
            "ok": abs(richardson["extrapolation_relative_difference"]) < RICHARDSON_CLAIM,
            "detail": "321x81 mesh series extrapolated onto the quadrature value",
        },
        "quadrature_self_convergence": {
            "value": quadrature_error,
            "claim": QUADRATURE_CLAIM,
            "ok": quadrature_error < QUADRATURE_CLAIM,
            "detail": "difference between the last two tanh-sinh refinements",
        },
        "monotone_convergence": {
            "value": sorted(monotone),
            "claim": [],
            "ok": not monotone,
            "detail": (
                "quantities whose error did not fall at every refinement; "
                f"excluded as already at the round-off floor: {sorted(round_off)}"
            ),
        },
    }

    payload = {
        "schema": "tensorcad.verification.hydrostatics/1",
        "generated_by": "examples/hydrostatics_verification/run_verification.py",
        "reference": reference.as_dict(),
        "ladder": rows,
        "convergence": summary,
        "richardson_wetted_area": richardson,
        "claims": claims,
    }
    figures = write_figures(rows, summary, args.output)
    payload["figures"] = figures
    (args.output / "hydrostatics_verification.json").write_text(json.dumps(payload, indent=2))

    print(f"closed form  V={reference.displacement_volume:.6f} m3  Cb={reference.block_coefficient:.6f} "
          f"Cp={reference.prismatic_coefficient:.6f}  VCB={reference.vertical_centre_buoyancy:.6f} m")
    print(f"{'nx':>4} {'nz':>3} {'V err %':>9} {'S err %':>9} {'worst %':>9} {'order(vol)':>10}")
    for index, row in enumerate(rows):
        order = summary["displacement_volume"]["observed_orders"][index - 1] if index else float("nan")
        print(f"{row['nx']:>4} {row['nz']:>3} "
              f"{100.0 * row['rows']['displacement_volume']['relative_error']:>9.4f} "
              f"{100.0 * row['rows']['wetted_surface_area']['relative_error']:>9.4f} "
              f"{100.0 * row['worst_relative_error']:>9.4f} {order:>10.3f}")
    print(f"richardson wetted area: {richardson['extrapolated']:.6f} m2 vs quadrature "
          f"{richardson['quadrature_reference']:.6f} m2 "
          f"({100.0 * richardson['extrapolation_relative_difference']:+.2e} %)")
    for name, claim in claims.items():
        print(f"{'PASS' if claim['ok'] else 'FAIL'}  {name}: {claim['value']}")
    return 0 if all(claim["ok"] for claim in claims.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
