"""Transom-stern family: verification ladder, design table, and figures.

Run with the standard library plus matplotlib::

    python3 examples/transom_stern/run_example.py --output outputs/transom_stern

Everything here is driven by the closed-form references in
``tensorcad.transom``; no number is typed by hand.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:  # allow running from a checkout without installing
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad import transom
from tensorcad.hydrostatics_analytic import taper_integral, tanh_sinh_quadrature

GRIDS = ((41, 13), (81, 25), (161, 41), (321, 81))
TAUS = (0.0, 0.05, 0.2, 0.35)
EXPONENTS = (0.7, 1.5, 3.0)


def reference_ladder() -> list[dict[str, float]]:
    """The incomplete taper integrals, three ways: series, closed form, quadrature."""
    rows: list[dict[str, float]] = []
    for exponent in (0.3, 0.7, 1.0, 1.5, 2.0, 3.0):
        closed = taper_integral(exponent)
        series = transom.truncated_taper_integral(0.0, exponent)
        rows.append(
            {
                "exponent": exponent,
                "series_at_tau0": series,
                "gamma_closed_form": closed,
                "series_vs_closed_form": series / closed - 1.0,
            }
        )
    for exponent in (0.3, 0.7, 1.5, 3.0):
        for tau in (0.02, 0.1, 0.2, 0.35, 0.45):
            series = transom.truncated_taper_integral(tau, exponent)
            quadrature = tanh_sinh_quadrature(
                lambda x, e=exponent: math.sin(math.pi * x) ** e, tau, 1.0, levels=7
            )
            rows.append(
                {
                    "exponent": exponent,
                    "tau": tau,
                    "series_incomplete": series,
                    "tanh_sinh_incomplete": quadrature,
                    "series_vs_quadrature": series / quadrature - 1.0,
                }
            )
    return rows


def convergence_rows() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for tau in TAUS:
        for nx, nz in GRIDS:
            model = transom.generate_transom_hull(
                transom.TransomHullParameters(nx=nx, nz=nz, transom_fraction=tau)
            )
            errors = model.relative_errors()
            rows.append(
                {
                    "tau": tau,
                    "nx": nx,
                    "nz": nz,
                    "faces": len(model.mesh.faces),
                    "volume_m3": model.hydrostatics.displacement_volume,
                    "volume_analytic_m3": model.analytic["displacement_volume"],
                    "volume_error": errors["displacement_volume"],
                    "waterplane_error": errors["waterplane_area"],
                    "lcb_error": errors["longitudinal_centre_buoyancy"],
                    "lcb_m": model.hydrostatics.longitudinal_centre_buoyancy,
                    "lcb_analytic_m": model.analytic["longitudinal_centre_buoyancy"],
                }
            )
    return rows


def design_table() -> list[dict[str, float]]:
    """The family as a design space: immersed transom area against Cb, Cp and LCB."""
    rows: list[dict[str, float]] = []
    for exponent in EXPONENTS:
        for ratio in (0.02, 0.05, 0.10, 0.15, 0.20):
            tau = transom.transom_fraction_for_area_ratio(ratio, exponent)
            parameters = transom.TransomHullParameters(
                transom_fraction=tau, breadth_exponent=exponent
            )
            analytic = transom.analytic_transom_hydrostatics(parameters)
            rows.append(
                {
                    "exponent": exponent,
                    "transom_area_ratio": ratio,
                    "tau": tau,
                    "transom_breadth_m": analytic["transom_breadth_at_waterline"],
                    "volume_m3": analytic["displacement_volume"],
                    "block_coefficient": analytic["block_coefficient"],
                    "prismatic_coefficient": analytic["block_coefficient"] / (2.0 / 3.0),
                    "waterplane_coefficient": analytic["waterplane_coefficient"],
                    "lcb_over_length": analytic["longitudinal_centre_buoyancy"] / 10.0,
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    if not rows:
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def make_figures(output: Path, convergence: list[dict[str, float]], design: list[dict[str, float]]):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures = []

    # 1. convergence of the mesh integration against the analytic reference
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))
    for tau in TAUS:
        subset = [row for row in convergence if row["tau"] == tau]
        spacing = [1.0 / (row["nx"] - 1) for row in subset]
        axes[0].loglog(spacing, [abs(row["volume_error"]) for row in subset], "o-",
                       label=f"tau = {tau:g}")
        axes[1].loglog(spacing, [abs(row["lcb_error"]) or 1e-16 for row in subset], "o-",
                       label=f"tau = {tau:g}")
    axes[0].set_xlabel("longitudinal spacing / LOA")
    axes[0].set_ylabel("|volume error| vs analytic")
    axes[0].set_title("Displacement converges to the closed form")
    axes[1].set_xlabel("longitudinal spacing / LOA")
    axes[1].set_ylabel("|LCB error| vs analytic")
    axes[1].set_title("LCB converges for every truncation")
    for axis in axes:
        axis.grid(True, which="both", alpha=0.3)
        axis.legend()
    fig.tight_layout()
    path = output / "convergence.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(("convergence.png", path))

    # 2. the geometry itself: profile and body plan of a transom body
    parameters = transom.TransomHullParameters(transom_fraction=0.2, nx=161, nz=41)
    transom.generate_transom_hull(parameters)  # validates the case being plotted
    ring_size = 2 * parameters.nz - 1

    def ring(station: int):
        x = transom._station_x(parameters, station)
        taper = transom._taper_at(parameters, x)
        return [transom._ring_point(parameters, x, taper, j) for j in range(ring_size)]

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.0))
    for station in range(parameters.nx):
        xs = ring(station)
        axes[0].plot([point[0] for point in xs], [point[2] for point in xs], color="0.75", lw=0.3)
    axes[0].plot([0.0, 0.0], [-parameters.draft, 0.0], color="tab:red", lw=2.0,
                 label="transom plate")
    axes[0].axhline(0.0, color="tab:blue", lw=1.0, ls="--", label="design waterline")
    axes[0].set_xlabel("x [m] (transom at 0)")
    axes[0].set_ylabel("z [m]")
    axes[0].set_title("Profile: stations and the transom cut")
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    rear = (parameters.nx - 1) // 4
    for label, station in (("transom", 0), ("1/4 LOA", rear), ("midship", (parameters.nx - 1) // 2)):
        xs = ring(station)
        axes[1].plot([point[1] for point in xs], [point[2] for point in xs], label=label)
    axes[1].set_xlabel("y [m]")
    axes[1].set_ylabel("z [m]")
    axes[1].set_title("Body plan: sections, including the transom")
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].set_aspect("equal")
    fig.tight_layout()
    path = output / "geometry.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(("geometry.png", path))

    # 3. block coefficient vs transom station for several exponents
    fig, axis = plt.subplots(figsize=(6.4, 4.2))
    taus = [0.005 + 0.005 * index for index in range(90)]
    for exponent in EXPONENTS:
        values = [
            transom.analytic_transom_hydrostatics(
                transom.TransomHullParameters(transom_fraction=tau, breadth_exponent=exponent)
            )["block_coefficient"]
            for tau in taus
        ]
        axis.plot(taus, values, label=f"p = {exponent:g}")
    axis.set_xlabel("transom fraction tau (LOA trim, aft cut)")
    axis.set_ylabel("block coefficient Cb over LOA")
    axis.set_title("The two parameters reach different block coefficients")
    axis.grid(True, alpha=0.3)
    axis.legend()
    fig.tight_layout()
    path = output / "design_space.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    figures.append(("design_space.png", path))
    return figures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("outputs/transom_stern"))
    arguments = parser.parse_args()
    output: Path = arguments.output
    output.mkdir(parents=True, exist_ok=True)

    ladder = reference_ladder()
    convergence = convergence_rows()
    design = design_table()
    write_csv(output / "reference_ladder.csv", ladder)
    write_csv(output / "convergence.csv", convergence)
    write_csv(output / "design_table.csv", design)
    figures = make_figures(output, convergence, design)

    print("== reference ladder: incomplete taper integral, three ways ==")
    print(f"{'p':>5} {'tau':>6} {'series':>20} {'independent':>20} {'rel':>10}")
    for row in ladder:
        if "tau" in row:
            print(f"{row['exponent']:5.2f} {row['tau']:6.2f} {row['series_incomplete']:20.12f} "
                  f"{row['tanh_sinh_incomplete']:20.12f} {row['series_vs_quadrature']:+10.2e}")
        else:
            print(f"{row['exponent']:5.2f} {0.0:6.2f} {row['series_at_tau0']:20.12f} "
                  f"{row['gamma_closed_form']:20.12f} {row['series_vs_closed_form']:+10.2e}")
    worst = max(abs(row.get("series_vs_quadrature", row.get("series_vs_closed_form", 0.0)))
                for row in ladder)
    print(f"worst reference disagreement: {worst:.2e}")

    print("\n== convergence of the mesh integration against the analytic reference ==")
    print(f"{'tau':>5} {'grid':>9} {'V mesh':>14} {'V analytic':>14} {'V err %':>9} "
          f"{'Awp err %':>10} {'LCB err %':>10}")
    for row in convergence:
        print(f"{row['tau']:5.2f} {row['nx']:4}x{row['nz']:<4} {row['volume_m3']:14.6f} "
              f"{row['volume_analytic_m3']:14.6f} {100 * row['volume_error']:+9.3f} "
              f"{100 * row['waterplane_error']:+10.3f} {100 * row['lcb_error']:+10.3f}")

    print("\n== design space: immersed transom area -> coefficients ==")
    print(f"{'p':>5} {'A_t/A_m':>8} {'tau':>8} {'V [m3]':>10} {'Cb':>8} {'Cp':>8} "
          f"{'Cwp':>8} {'LCB/LOA':>8}")
    for row in design:
        print(f"{row['exponent']:5.2f} {row['transom_area_ratio']:8.2f} {row['tau']:8.5f} "
              f"{row['volume_m3']:10.4f} {row['block_coefficient']:8.4f} "
              f"{row['prismatic_coefficient']:8.4f} {row['waterplane_coefficient']:8.4f} "
              f"{row['lcb_over_length']:8.4f}")

    print("\n== artefacts ==")
    for name, path in figures:
        print(f"  {name}: {path} ({path.stat().st_size} bytes)")
    for name in ("reference_ladder.csv", "convergence.csv", "design_table.csv"):
        print(f"  {name}: {output / name}")

    # One last invariant, printed rather than asserted: the topology of each case.
    print("\n== topology ==")
    from tensorcad import mesh as meshlib

    for tau in TAUS:
        model = transom.generate_transom_hull(transom.TransomHullParameters(transom_fraction=tau))
        euler = len(model.mesh.vertices) - len(meshlib.edge_face_counts(model.mesh)) + len(model.mesh.faces)
        print(f"  tau={tau:4}: watertight={meshlib.is_watertight(model.mesh)} "
              f"manifold={meshlib.is_edge_manifold(model.mesh)} euler={euler} "
              f"faces={len(model.mesh.faces)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
