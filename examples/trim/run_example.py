"""Trim and the centimetre constants: every claim in docs/trim.md, re-derived here.

Run from the repository root::

    python examples/trim/run_example.py

The script exits non-zero if any claim fails, writes ``results.json`` next to the
figure, and prints the tables the document quotes. It is deliberately cheap: a
box, the parametric hull and the Wigley hull, each a few hundred station cuts.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad import (  # noqa: E402
    HullParameters,
    analytic_volume,
    generate_hull,
    generate_wigley,
)
from tensorcad.bonjean import section_properties  # noqa: E402
from tensorcad.geometry import Mesh  # noqa: E402
from tensorcad.hydrostatics_analytic import waterplane_inertia_longitudinal  # noqa: E402
from tensorcad.stability import waterplane_moments  # noqa: E402
from tensorcad.trim import (  # noqa: E402
    SEAWATER_DENSITY,
    box_analytic_particulars,
    box_trim_tangent,
    equilibrium_trim,
    hydrostatic_particulars,
    hydrostatic_table,
    small_angle_trim,
)
from tensorcad.wigley import half_breadth  # noqa: E402

OUT = Path(__file__).resolve().parent
FAILURES: list[str] = []
ROWS: dict[str, object] = {}


def claim(name: str, value: float, limit: float, *, relative: bool = True) -> float:
    bad = (abs(value) > limit) if relative else (value > limit)
    status = "FAIL" if bad else "ok"
    if bad:
        FAILURES.append(f"{name}: {value:+.3e} (limit {limit:.1e})")
    print(f"  [{status}] {name}: {value:+.3e}  (limit {limit:.1e})")
    return value


def box_mesh(length: float, beam: float, draft: float, freeboard: float) -> Mesh:
    """A watertight rectangular hull: 8 vertices, 12 triangles, outward wound."""
    xs = [0.0, length]
    ys = [-0.5 * beam, 0.5 * beam]
    zs = [-draft, freeboard]
    vertices: list[tuple[float, float, float]] = []
    index: dict[tuple[float, float, float], int] = {}
    faces: list[tuple[int, int, int]] = []

    def vid(point: tuple[float, float, float]) -> int:
        if point not in index:
            index[point] = len(vertices)
            vertices.append(point)
        return index[point]

    def quad(a, b, c, d) -> None:
        faces.append((vid(a), vid(b), vid(c)))
        faces.append((vid(a), vid(c), vid(d)))

    quad((xs[0], ys[0], zs[0]), (xs[0], ys[1], zs[0]), (xs[1], ys[1], zs[0]), (xs[1], ys[0], zs[0]))
    quad((xs[0], ys[0], zs[1]), (xs[1], ys[0], zs[1]), (xs[1], ys[1], zs[1]), (xs[0], ys[1], zs[1]))
    quad((xs[0], ys[0], zs[0]), (xs[1], ys[0], zs[0]), (xs[1], ys[0], zs[1]), (xs[0], ys[0], zs[1]))
    quad((xs[1], ys[1], zs[0]), (xs[0], ys[1], zs[0]), (xs[0], ys[1], zs[1]), (xs[1], ys[1], zs[1]))
    quad((xs[0], ys[1], zs[0]), (xs[0], ys[0], zs[0]), (xs[0], ys[0], zs[1]), (xs[0], ys[1], zs[1]))
    quad((xs[1], ys[0], zs[0]), (xs[1], ys[1], zs[0]), (xs[1], ys[1], zs[1]), (xs[1], ys[0], zs[1]))
    return Mesh(vertices=tuple(vertices), faces=tuple(faces))


def relative(measured: float, reference: float) -> float:
    return (measured - reference) / reference


def wigley_reference_i_l(parameters, samples: int = 2000) -> float:
    """I_L about midship by Simpson on the analytic half-breadth, not on the mesh."""
    length = parameters.length
    values = []
    for index in range(samples + 1):
        x = length * index / samples
        y = half_breadth(parameters, x / length, 1.0)  # depth fraction 1 = design waterline
        values.append((2.0 / 3.0) * y**3 * (x - 0.5 * length) ** 2)
    total = values[0] + values[-1]
    for index in range(1, samples):
        total += (4.0 if index % 2 else 2.0) * values[index]
    return total * (length / samples) / 3.0


def main() -> int:
    print("== 1. box hull: the mesh route against the closed forms ==")
    length, beam, draft = 20.0, 4.0, 5.0
    box = box_mesh(length, beam, draft, freeboard=2.0)
    measured = hydrostatic_particulars(box, 0.0, station_count=41)
    exact = box_analytic_particulars(length, beam, draft)
    box_row = {}
    for field in (
        "volume",
        "lcb",
        "lcf",
        "kb",
        "waterplane_area",
        "transverse_second_moment",
        "longitudinal_second_moment",
        "tpc",
        "mct",
    ):
        error = relative(getattr(measured, field), getattr(exact, field))
        box_row[field] = {"measured": getattr(measured, field), "exact": getattr(exact, field), "rel": error}
        claim(f"box {field}", error, 1e-12)
    ROWS["box_particulars"] = box_row

    print("== 2. box hull: trim at finite angle against the wallsided closed form ==")
    print("   (only angles whose waterline stays on the vertical sides: beyond that the")
    print("    closed form itself stops applying, and the mesh route is the honest one)")
    trim_rows = []
    for weight, shift in ((24000.0, 4.0), (48000.0, 8.0)):
        solution = equilibrium_trim(box, weight=weight, shift=shift, lcb_tolerance=1e-10)
        closed = math.degrees(math.atan(box_trim_tangent(length, beam, weight, shift)))
        error = relative(solution.trim_angle_deg, closed)
        trim_rows.append(
            {
                "weight_kg": weight,
                "shift_m": shift,
                "mesh_deg": solution.trim_angle_deg,
                "closed_form_deg": closed,
                "rel": error,
                "iterations": solution.iterations,
                "lcb_residual_m": solution.residual_m,
            }
        )
        claim(f"box trim {weight:.0f} kg x {shift:.0f} m", error, 1e-6)
    ROWS["box_trim"] = trim_rows

    _, _, change_cm = small_angle_trim(measured, 48000.0, 8.0)
    loading_route_cm = 48000.0 / 1000.0 * 8.0 / measured.mct
    claim("box two centimetre routes", relative(change_cm, loading_route_cm), 1e-12)
    ROWS["box_centimetre_routes"] = {"from_I_L_cm": change_cm, "from_MCT_cm": loading_route_cm}

    print("== 3. parametric hull: particulars at a part load waterline ==")
    hull_parameters = HullParameters(length=10.0, beam=2.0, draft=0.8)
    hull = generate_hull(hull_parameters)
    z_part = -0.2 * hull_parameters.draft
    particulars = hydrostatic_particulars(hull.mesh, z_part, station_count=41)
    analytic_v = (
        analytic_volume(hull_parameters.length, hull_parameters.beam, hull_parameters.draft)
        * (particulars.draft / hull_parameters.draft) ** 1.5
    )
    claim("hull part-load volume vs analytic", relative(particulars.volume, analytic_v), 1.5e-2)
    step = 1e-3
    up = hydrostatic_particulars(hull.mesh, z_part + step, station_count=41).volume
    down = hydrostatic_particulars(hull.mesh, z_part - step, station_count=41).volume
    tpc_from_derivative = 0.5 * (up - down) / step * SEAWATER_DENSITY / 1e5
    claim("hull TPC vs dV/dz", relative(particulars.tpc, tpc_from_derivative), 5e-3)
    i_l, quadrature_evidence = waterplane_inertia_longitudinal(
        hull_parameters.length, hull_parameters.beam
    )
    mct_analytic = SEAWATER_DENSITY / 1000.0 * i_l / (100.0 * hull_parameters.length)
    at_design = hydrostatic_particulars(hull.mesh, 0.0, station_count=41)
    claim("hull MCT at the design waterline vs analytic I_L", relative(at_design.mct, mct_analytic), 1.5e-2)
    ROWS["hull_particulars"] = {
        "waterline_z": z_part,
        "draft": particulars.draft,
        "volume": particulars.volume,
        "analytic_volume": analytic_v,
        "volume_rel": relative(particulars.volume, analytic_v),
        "tpc": particulars.tpc,
        "tpc_from_dV_dz": tpc_from_derivative,
        "mct": at_design.mct,
        "mct_analytic": mct_analytic,
        "analytic_quadrature_evidence": quadrature_evidence,
        "lcb": particulars.lcb,
        "lcf": particulars.lcf,
        "kb": particulars.kb,
        "km_transverse": particulars.km_transverse,
        "lcb_minus_lcf": particulars.lcb_relative_to_lcf,
    }

    print("== 4. parametric hull: the trim solve against the small-angle answer ==")
    hull_solutions = []
    for weight, shift in ((300.0, 1.0), (150.0, -1.0)):
        solution = equilibrium_trim(
            hull.mesh, weight=weight, shift=shift, waterline_z=z_part, station_count=41
        )
        error = relative(solution.trim_angle_deg, solution.small_angle_deg)
        hull_solutions.append(
            {
                "weight_kg": weight,
                "shift_m": shift,
                "mesh_deg": solution.trim_angle_deg,
                "small_angle_deg": solution.small_angle_deg,
                "rel": error,
                "change_of_trim_cm": solution.change_of_trim_cm,
                "mct_t_m_per_cm": solution.mct,
                "waterline_at_reference_m": solution.waterline_at_reference,
                "iterations": solution.iterations,
                "lcb_residual_m": solution.residual_m,
            }
        )
        claim(f"hull trim {weight:.0f} kg x {shift:+.0f} m", error, 5e-3)
    ROWS["hull_trim"] = hull_solutions

    print("== 5. hull: change of trim against stations (Simpson convergence) ==")
    convergence = []
    reference_angle = None
    for station_count in (21, 41, 81, 161):
        solution = equilibrium_trim(
            hull.mesh, weight=300.0, shift=1.0, waterline_z=z_part, station_count=station_count
        )
        if station_count == 161:
            reference_angle = solution.trim_angle_deg
        convergence.append(
            {"stations": station_count, "trim_deg": solution.trim_angle_deg, "iterations": solution.iterations}
        )
    for row in convergence:
        row["rel_to_finest"] = relative(row["trim_deg"], reference_angle)
    ROWS["hull_trim_convergence"] = convergence
    print("   ", [(row["stations"], round(row["trim_deg"], 9)) for row in convergence])

    print("== 6. Wigley: immersion constant and trim against the analytic waterline ==")
    wigley = generate_wigley()
    wigley_z = -0.2 * wigley.parameters.draft
    wigley_particulars = hydrostatic_particulars(wigley.mesh, wigley_z, station_count=41)
    wigley_step = 1e-3
    w_up = hydrostatic_particulars(wigley.mesh, wigley_z + wigley_step, station_count=41).volume
    w_down = hydrostatic_particulars(wigley.mesh, wigley_z - wigley_step, station_count=41).volume
    wigley_dvdz = 0.5 * (w_up - w_down) / wigley_step * SEAWATER_DENSITY / 1e5
    claim("wigley TPC vs dV/dz", relative(wigley_particulars.tpc, wigley_dvdz), 5e-3)
    wigley_solution = equilibrium_trim(wigley.mesh, weight=300.0, shift=0.5, waterline_z=wigley_z, station_count=41)
    claim(
        "wigley trim mesh vs small-angle",
        relative(wigley_solution.trim_angle_deg, wigley_solution.small_angle_deg),
        5e-3,
    )
    ROWS["wigley"] = {
        "draft": wigley_particulars.draft,
        "volume": wigley_particulars.volume,
        "tpc": wigley_particulars.tpc,
        "tpc_from_dV_dz": wigley_dvdz,
        "i_l_mesh": waterplane_moments(wigley.mesh, wigley_z).longitudinal_second_moment,
        "trim_deg": wigley_solution.trim_angle_deg,
        "small_angle_deg": wigley_solution.small_angle_deg,
        "mct": hydrostatic_particulars(wigley.mesh, wigley_z, station_count=41).mct,
    }

    print("== 7. hydrostatic table: TPC and MCT against draft, and their ratios ==")
    drafts = [-0.08 * hull_parameters.draft * factor for factor in (9, 8, 7, 6, 5, 4, 3, 2, 1)]
    table = hydrostatic_table(hull.mesh, drafts, station_count=41)
    table_rows = [
        {
            "draft": row.draft,
            "volume": row.volume,
            "tpc": row.tpc,
            "mct": row.mct,
            "lcb_minus_lcf": row.lcb_relative_to_lcf,
            "km_transverse": row.km_transverse,
            "nudged": row.waterline_nudged,
        }
        for row in table
    ]
    ROWS["hydrostatic_table"] = table_rows
    monotone_tpc = all(
        table_rows[index + 1]["tpc"] > table_rows[index]["tpc"] for index in range(len(table_rows) - 1)
    )
    monotone_mct = all(
        table_rows[index + 1]["mct"] > table_rows[index]["mct"] for index in range(len(table_rows) - 1)
    )
    print(f"    TPC increasing with draft: {monotone_tpc}; MCT increasing: {monotone_mct}")
    if not (monotone_tpc and monotone_mct):
        FAILURES.append("TPC/MCT are not monotone in draft")

    plot(convergence, table_rows, hull, z_part, hull_solutions[0], box, trim_rows)
    ROWS["box_trim_reference"] = {
        "weight_kg": 48000.0,
        "shift_m": 8.0,
        "closed_form_deg": math.degrees(math.atan(box_trim_tangent(length, beam, 48000.0, 8.0))),
    }
    (OUT / "results.json").write_text(json.dumps(ROWS, indent=2, sort_keys=True))
    print(f"\nwrote {OUT / 'results.json'} and {OUT / 'trim.png'}")

    if FAILURES:
        print("\nFAILED claims:")
        for failure in FAILURES:
            print(" -", failure)
        return 1
    print("\nAll claims held.")
    return 0


def plot(convergence, table_rows, hull, z_part, hull_solution, box, trim_rows) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 2, figsize=(11.0, 8.0))

    drafts = [row["draft"] for row in table_rows]
    axes[0][0].plot(drafts, [row["tpc"] for row in table_rows], "o-", label="hull (mesh)")
    box_row = box_analytic_particulars(20.0, 4.0, 5.0)
    axes[0][0].plot(
        [box_row.draft], [box_row.tpc], "s", color="tab:red", label=f"box closed form {box_row.tpc:.4f}"
    )
    axes[0][0].set_xlabel("draft from keel [m]")
    axes[0][0].set_ylabel("TPC [t/cm]")
    axes[0][0].set_title("immersion constant against draft")
    axes[0][0].legend()
    axes[0][0].grid(alpha=0.3)

    axes[0][1].plot(drafts, [row["mct"] for row in table_rows], "o-", color="tab:green", label="hull (mesh)")
    axes[0][1].set_xlabel("draft from keel [m]")
    axes[0][1].set_ylabel("MCT [t m/cm]")
    axes[0][1].set_title("trim constant against draft")
    axes[0][1].legend()
    axes[0][1].grid(alpha=0.3)

    box_length = 20.0
    mesh_values = [row["mesh_deg"] for row in trim_rows]
    closed_values = [row["closed_form_deg"] for row in trim_rows]
    axes[1][0].plot(closed_values, mesh_values, "o", label="box trim: mesh solve")
    limit = max(closed_values) * 1.1
    axes[1][0].plot([0.0, limit], [0.0, limit], "k--", linewidth=1, label="mesh = wallsided closed form")
    for row in trim_rows:
        axes[1][0].annotate(
            f"{row['rel']:+.1e}", (row["closed_form_deg"], row["mesh_deg"]), textcoords="offset points", xytext=(6, -10), fontsize=8
        )
    axes[1][0].set_xlabel("wallsided closed form trim [deg]")
    axes[1][0].set_ylabel("mesh solve trim [deg]")
    axes[1][0].set_title(f"box {box_length:.0f} m: exact at finite trim (rel. error annotated)")
    axes[1][0].legend()
    axes[1][0].grid(alpha=0.3)

    station_axis = [row["stations"] for row in convergence]
    axes[1][1].semilogy(
        station_axis, [abs(row["rel_to_finest"]) for row in convergence], "o-", color="tab:purple"
    )
    axes[1][1].set_xlabel("stations used for the section integral")
    axes[1][1].set_ylabel("|change of trim - finest| / finest")
    axes[1][1].set_title("trim solve converges as the section integral refines")
    axes[1][1].grid(alpha=0.3, which="both")

    figure.suptitle(
        f"TensorCAD trim: hull L/B/T = {hull.parameters.length:.1f}/"
        f"{hull.parameters.beam:.1f}/{hull.parameters.draft:.1f} m, "
        f"trim {hull_solution['mesh_deg']:+.4f} deg under {hull_solution['weight_kg']:.0f} kg x "
        f"{hull_solution['shift_m']:.0f} m"
    )
    figure.tight_layout()
    figure.savefig(Path(__file__).resolve().parent / "trim.png", dpi=130)


if __name__ == "__main__":
    raise SystemExit(main())