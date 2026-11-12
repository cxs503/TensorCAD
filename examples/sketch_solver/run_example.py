"""The GAP-CAD-006 acceptance run: fixtures, restarts, unit scaling, and the figure.

Everything in ``docs/sketch_constraints.md`` comes from this script. It does four things,
in the order a sceptical reader would want them:

1. **Solves every fixture** and measures the answer against the analytic value that the
   fixture itself states (law of cosines, ``h/2``, ``L sin 30``, ...).
2. **Re-checks the residuals independently**: the returned parameter vector is handed back
   to :meth:`Sketch.residuals`, so the "residual is zero" claim does not come from the
   solver's own bookkeeping.
3. **Restarts from perturbed coordinates** -- 25 deterministic perturbations per solvable
   fixture -- because a solver that only converges from the coordinates its own author
   typed has not been tested. The convergence rate is reported, not assumed to be 1.
4. **Repeats the rectangle in millimetres, metres and kilometres** to show the rank/DOF
   answer does not depend on the drawing's units (the Jacobian is row-normalised before
   the rank is taken).

Outputs: ``outputs/results.json`` (every number, program-generated) and ``outputs/plot.png``
(four panels: three solved sketches with the analytic answer marked, and the residual
history of four cases on a log scale). Run with ``--write-doc`` to splice the generated
table into ``docs/sketch_constraints.md`` between its markers.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import statistics
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad import (
    Angle,
    Distance,
    Fix,
    FixX,
    Horizontal,
    InconsistentConstraintsError,
    LineCircleTangency,
    Sketch,
    Vertical,
)
from tensorcad.sketch_cases import CASES, run_all

HERE = Path(__file__).resolve().parent
OUTPUTS = HERE / "outputs"
DOC = HERE.parents[1] / "docs" / "sketch_constraints.md"
RESTARTS = 25
SEED = 20261112
SCALES = (1e-3, 1.0, 1e3)


# -- studies ----------------------------------------------------------------------------
def study_cases() -> list[dict]:
    """Solve every fixture; locate the outcome against the analytic answer."""
    records = []
    for outcome in run_all():
        case = next(c for c in CASES if c.name == outcome.case)
        # Independent re-check: rebuild, solve, and evaluate the residuals at the returned
        # parameter vector through the public residual entry point.
        sketch = case.build()
        try:
            result = sketch.solve()
        except InconsistentConstraintsError as error:
            assert error.result is not None
            result = error.result
        recheck = max((abs(v) for v in sketch.residuals(result.values)), default=0.0)
        records.append(
            {
                "case": outcome.case,
                "description": outcome.descriptions,
                "status": outcome.status,
                "dof": outcome.dof,
                "rank": outcome.rank,
                "parameters": outcome.parameters,
                "residual_rows": outcome.residual_rows,
                "redundant_rows": outcome.redundant_rows,
                "null_space_dimension": outcome.null_space_dimension,
                "iterations": outcome.iterations,
                "max_residual": outcome.max_residual,
                "recheck_residual": recheck,
                "max_position_error": outcome.max_position_error,
                "max_radius_error": outcome.max_radius_error,
                "blocking": list(outcome.blocking),
                "raised": outcome.raised,
                "analytic_positions": {k: list(v) for k, v in outcome.analytic_positions.items()},
                "solved_positions": {k: list(v) for k, v in outcome.solved_positions.items()},
                "analytic_radii": {k: v for k, v in case.radii.items()},
                "solved_radii": {k: v for k, v in result.radii.items()},
                "note": case.note,
            }
        )
    return records


def study_restarts() -> dict:
    """Perturbed starts: does the solve depend on the fixture author's coordinates?"""
    rng = random.Random(SEED)
    records = []
    total = 0
    converged = 0
    iterations: list[int] = []
    worst_position_error = 0.0
    for case in CASES:
        if case.status == "inconsistent":
            continue
        succeeded = 0
        case_iterations: list[int] = []
        for _ in range(RESTARTS):
            sketch = case.build()
            sketch.perturb(0.25, rng)
            result = sketch.solve(strict=False)
            total += 1
            if result.max_residual <= 1e-9:
                succeeded += 1
                converged += 1
                case_iterations.append(result.iterations)
                iterations.append(result.iterations)
            for name, expected in case.positions.items():
                solved = result.positions[name]
                worst_position_error = max(
                    worst_position_error,
                    math.hypot(solved[0] - expected[0], solved[1] - expected[1]),
                )
        records.append(
            {
                "case": case.name,
                "restarts": RESTARTS,
                "converged": succeeded,
                "median_iterations": statistics.median(case_iterations) if case_iterations else 0,
            }
        )
    return {
        "restarts_per_case": RESTARTS,
        "perturbation": 0.25,
        "seed": SEED,
        "attempts": total,
        "converged": converged,
        "rate": converged / total if total else 0.0,
        "median_iterations": statistics.median(iterations) if iterations else 0,
        "worst_position_error": worst_position_error,
        "cases": records,
    }


def _rectangle(scale: float) -> tuple[Sketch, dict[str, tuple[float, float]]]:
    sketch = Sketch()
    a = sketch.point(0.0, 0.0, "A")
    b = sketch.point(4.0 * scale, 0.3 * scale, "B")
    c = sketch.point(4.2 * scale, 3.0 * scale, "C")
    d = sketch.point(0.1 * scale, 3.0 * scale, "D")
    ab, bc, cd, da = (sketch.line(a, b), sketch.line(b, c), sketch.line(c, d), sketch.line(d, a))
    sketch.add(
        Fix(a, 0.0, 0.0),
        Horizontal(ab),
        Vertical(bc),
        Horizontal(cd),
        Vertical(da),
        Distance(a, b, 4.0 * scale),
        Distance(b, c, 3.0 * scale),
    )
    expected = {
        "A": (0.0, 0.0),
        "B": (4.0 * scale, 0.0),
        "C": (4.0 * scale, 3.0 * scale),
        "D": (0.0, 3.0 * scale),
    }
    return sketch, expected


def study_scales() -> list[dict]:
    """The same rectangle drawn in millimetres, metres and kilometres."""
    records = []
    for scale in SCALES:
        sketch, expected = _rectangle(scale)
        result = sketch.solve(tol=1e-12 * scale)
        error = 0.0
        for name, want in expected.items():
            got = result.positions[name]
            error = max(error, math.hypot(got[0] - want[0], got[1] - want[1]))
        records.append(
            {
                "scale": scale,
                "status": result.status,
                "dof": result.dof,
                "rank": result.rank,
                "residual_rows": result.residual_rows,
                "max_residual": result.max_residual,
                "max_position_error": error,
                "relative_error": error / (4.0 * scale),
            }
        )
    return records


def study_history() -> dict[str, list[float]]:
    """Residual history of four cases, for the log-scale panel."""
    wanted = ("sss_triangle", "angle_arm", "circle_between_lines", "point_on_circle")
    history: dict[str, list[float]] = {}
    for case in CASES:
        if case.name not in wanted:
            continue
        sketch = case.build()
        result = sketch.solve()
        # Clip at the double-precision floor so a log axis can draw the last step.
        history[case.name] = [max(value, 1e-18) for value in result.history]
    return history


# -- figure -----------------------------------------------------------------------------
def draw_figure(results: list[dict], history: dict[str, list[float]]) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    solved = {record["case"]: record for record in results}
    figure, axes = plt.subplots(2, 2, figsize=(11.0, 9.0))
    figure.suptitle(
        "Sketch constraint solving (GAP-CAD-006): solved geometry and Newton convergence",
        fontsize=13,
    )

    def draw_case(ax, name: str, title: str) -> None:
        case = next(c for c in CASES if c.name == name)
        sketch = case.build()
        sketch.solve()
        for line in sketch.lines:
            (x0, y0), (x1, y1) = sketch.position(line.start), sketch.position(line.end)
            ax.plot([x0, x1], [y0, y1], "-", color="tab:blue", lw=1.8)
        for circle in sketch.circles:
            cx, cy = sketch.position(circle.center)
            theta = [i * 2 * math.pi / 200 for i in range(201)]
            radius = sketch.radius(circle)
            ax.plot([cx + radius * math.cos(t) for t in theta],
                    [cy + radius * math.sin(t) for t in theta], "--", color="tab:green", lw=1.2)
        for point in sketch.points:
            x, y = sketch.position(point)
            ax.plot([x], [y], "o", color="tab:blue", ms=5)
        for name_of_point, (x, y) in solved[name]["analytic_positions"].items():
            ax.plot([x], [y], "+", color="tab:red", ms=11, mew=1.6)
        ax.plot([], [], "o", color="tab:blue", label="solved")
        ax.plot([], [], "+", color="tab:red", label="analytic")
        ax.set_title(title)
        ax.set_aspect("equal", adjustable="datalim")
        ax.grid(alpha=0.3)
        ax.legend(loc="best", fontsize=8)

    draw_case(axes[0][0], "rectangle",
              "4 x 3 rectangle: DOF 0, max |position error| = "
              f"{solved['rectangle']['max_position_error']:.1e} m")
    draw_case(axes[0][1], "sss_triangle",
              "SSS triangle (4, 3, 2 m): apex from the law of cosines, error "
              f"{solved['sss_triangle']['max_position_error']:.1e} m")
    circle_radius = solved["circle_between_lines"]["solved_radii"]["C"]
    draw_case(axes[1][0], "circle_between_lines",
              "circle tangent to two lines 2 m apart: r = h/2 = "
              f"{circle_radius:.6f} m (analytic 1 m)")

    axis = axes[1][1]
    for name, values in history.items():
        axis.semilogy(range(len(values)), values, marker="o", ms=3, label=name)
    axis.set_xlabel("Levenberg-Marquardt iteration")
    axis.set_ylabel("max |constraint residual|  [m]")
    axis.set_title("Residual history (exact Jacobian, so the drop is fast)")
    axis.grid(alpha=0.3, which="both")
    axis.legend(fontsize=8)

    figure.tight_layout(rect=(0, 0, 1, 0.96))
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    path = OUTPUTS / "plot.png"
    figure.savefig(path, dpi=130)
    plt.close(figure)
    return path


# -- output -----------------------------------------------------------------------------
def summary(results: list[dict], restarts: dict, scales: list[dict]) -> dict:
    solvable = [r for r in results if r["status"] != "inconsistent"]
    return {
        "fixtures": len(results),
        "acceptance_bound": 1e-10,
        "max_position_error": max(r["max_position_error"] for r in results),
        "max_radius_error": max(r["max_radius_error"] for r in results),
        "max_residual": max(r["max_residual"] for r in solvable),
        "max_recheck_residual": max(r["recheck_residual"] for r in solvable),
        "all_converged": all(r["max_residual"] <= 1e-9 for r in solvable),
        "inconsistent_fixtures": [r["case"] for r in results if r["status"] == "inconsistent"],
        "redundant_rows_reported": sum(r["redundant_rows"] for r in results),
        "restart_rate": restarts["rate"],
        "restart_attempts": restarts["attempts"],
        "max_position_error_over_scales": max(s["relative_error"] for s in scales),
        "scale_dof": [s["dof"] for s in scales],
    }


def doc_table(results: list[dict], restarts: dict, scales: list[dict]) -> str:
    lines = [
        "| fixture | what it checks | status | DOF | rank / rows | iterations | max residual [m] "
        "| max position error [m] |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for record in results:
        lines.append(
            "| `{case}` | {description} | {status} | {dof} | {rank} / {rows} | {it} | "
            "{res:.1e} | {err:.1e} |".format(
                case=record["case"],
                description=record["description"],
                status=record["status"],
                dof=record["dof"],
                rank=record["rank"],
                rows=record["residual_rows"],
                it=record["iterations"],
                res=record["max_residual"],
                err=max(record["max_position_error"], record["max_radius_error"]),
            )
        )
    lines.append("")
    lines.append(
        "Perturbed-restart study: **{conv} / {att}** solves converged "
        "({per_case} deterministic restarts per solvable fixture, perturbation ±0.25 m), "
        "median {med:g} iterations, worst position drift {worst:.1e} m.".format(
            conv=restarts["converged"],
            att=restarts["attempts"],
            per_case=restarts["restarts_per_case"],
            med=restarts["median_iterations"],
            worst=restarts["worst_position_error"],
        )
    )
    lines.append("")
    lines.append("| drawing scale | DOF | rank / rows | max residual [m] | relative position error |")
    lines.append("|---:|---:|---:|---:|---:|")
    for record in scales:
        lines.append(
            "| {scale:g} | {dof} | {rank} / {rows} | {res:.1e} | {rel:.1e} |".format(
                scale=record["scale"],
                dof=record["dof"],
                rank=record["rank"],
                rows=record["residual_rows"],
                res=record["max_residual"],
                rel=record["relative_error"],
            )
        )
    return "\n".join(lines)


def write_doc(table: str) -> None:
    if not DOC.exists():
        raise SystemExit(f"cannot splice: {DOC} is missing")
    text = DOC.read_text()
    start, end = "<!-- TABLE:START -->", "<!-- TABLE:END -->"
    if start not in text or end not in text:
        raise SystemExit("the document lost its table markers")
    head, _, rest = text.partition(start)
    _, _, tail = rest.partition(end)
    DOC.write_text(f"{head}{start}\n{table}\n{end}{tail}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-doc", action="store_true",
                        help="splice the generated table into docs/sketch_constraints.md")
    args = parser.parse_args()

    results = study_cases()
    restarts = study_restarts()
    scales = study_scales()
    history = study_history()
    totals = summary(results, restarts, scales)

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    payload = {
        "generated_by": "examples/sketch_solver/run_example.py",
        "acceptance": totals,
        "cases": results,
        "restarts": restarts,
        "scales": scales,
        "residual_history": history,
    }
    (OUTPUTS / "results.json").write_text(json.dumps(payload, indent=2) + "\n")
    figure = draw_figure(results, history)
    table = doc_table(results, restarts, scales)
    (OUTPUTS / "table.md").write_text(table + "\n")
    if args.write_doc:
        write_doc(table)

    under = sum(1 for r in results if r["status"] == "under_constrained")
    inconsistent = sum(1 for r in results if r["status"] == "inconsistent")
    print(f"fixtures            : {totals['fixtures']} "
          f"({under} under-constrained, {inconsistent} inconsistent, rest fully constrained)")
    print(f"max position error  : {totals['max_position_error']:.2e} m "
          f"(bound {totals['acceptance_bound']:.0e})")
    print(f"max residual        : {totals['max_residual']:.2e} m")
    print(f"independent recheck : {totals['max_recheck_residual']:.2e} m")
    print(f"restart convergence : {totals['restart_rate'] * 100:.1f}% "
          f"({restarts['converged']}/{restarts['attempts']})")
    print(f"worst relative error over scales: {totals['max_position_error_over_scales']:.2e}")
    print(f"figure              : {figure.relative_to(HERE)}")
    assert totals["all_converged"], "a consistent fixture failed to converge"
    assert totals["max_position_error"] <= totals["acceptance_bound"]
    return 0


if __name__ == "__main__":
    raise SystemExit(main())