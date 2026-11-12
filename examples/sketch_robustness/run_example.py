"""Evidence run for the two sketch-solver defects of m16386 (degenerate start, status).

Writes ``artifacts/sketch_robustness/results.json`` (every number below, program
generated) and ``artifacts/sketch_robustness/plot.png``. Run from the repository root:

    python examples/sketch_robustness/run_example.py

Three studies, each of them the probe that found the defect rather than a restatement of
the fix:

1. **Coincident start.** The same sketches the acceptance fixtures use -- a single
   distance, an SSS triangle, a radius -- but started from coordinates that put a length
   row's argument at zero. Reported: the analytic answer, the solved answer, the error,
   whether the deterministic escape fired, and the iteration count.
2. **Restart robustness from a degenerate start.** 25 deterministic perturbations of the
   *coincident* start (the previous mission's restart study perturbed a well-posed start
   by +-0.25 m, which never reaches the singular set). Reported: convergence rate and
   worst error.
3. **The status contract.** For a system that cannot be satisfied from its start, the
   structural verdict (``dof``, rank) is the same whether the residual is 1e-16 or 0.67;
   only the status may differ. Reported for a sweep of ``tol``/``residual_tol`` bands:
   the residual reached, the status now returned, and what the previous ladder (a
   ``dof == 0`` verdict with no convergence test) would have claimed.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.sketch import (
    SINGULAR_START_NUDGE,
    Distance,
    Fix,
    Radius,
    Sketch,
)


def _cases() -> dict[str, dict]:
    """The three degenerate-start cases, each with its closed-form answer."""
    cases: dict[str, dict] = {}

    def single_distance() -> tuple[Sketch, list[str], float]:
        sketch = Sketch()
        origin = sketch.point(0.0, 0.0)
        free = sketch.point(0.0, 0.0)
        sketch.add(Fix(origin, 0.0, 0.0), Distance(free, origin, 5.0))
        return sketch, ["P0", "P1"], 5.0

    def sss_triangle() -> tuple[Sketch, list[str], float]:
        sketch = Sketch()
        a = sketch.point(0.0, 0.0)
        b = sketch.point(0.0, 0.0)
        c = sketch.point(0.0, 0.0)
        sketch.add(Fix(a, 0.0, 0.0), Distance(b, a, 4.0), Distance(c, a, 3.0), Distance(c, b, 2.0))
        return sketch, ["P0", "P1", "P2"], 4.0

    for name, builder, targets in (
        ("single_distance", single_distance, (("P0", "P1", 5.0),)),
        ("sss_triangle", sss_triangle, (("P0", "P1", 4.0), ("P0", "P2", 3.0), ("P1", "P2", 2.0))),
    ):
        sketch, _, _ = builder()
        result = sketch.solve()
        rows = []
        for first, second, target in targets:
            solved = math.dist(result.positions[first], result.positions[second])
            rows.append(
                {
                    "pair": f"{first}-{second}",
                    "target": target,
                    "solved": solved,
                    "error": solved - target,
                }
            )
        cases[name] = {
            "kind": "coincident_start",
            "status": result.status,
            "dof": result.dof,
            "converged": result.converged,
            "nudged_from_singular_start": result.nudged_from_singular_start,
            "iterations": result.iterations,
            "max_residual": result.max_residual,
            "rows": rows,
            "worst_abs_error": max(abs(row["error"]) for row in rows),
        }

    sketch = Sketch()
    centre = sketch.point(0.0, 0.0)
    circle = sketch.circle(centre, 0.0)
    sketch.add(Fix(centre, 0.0, 0.0), Radius(circle, 2.5))
    result = sketch.solve()
    radius = result.radii["C0"]
    cases["zero_radius"] = {
        "kind": "zero_radius_start",
        "status": result.status,
        "dof": result.dof,
        "converged": result.converged,
        "nudged_from_singular_start": result.nudged_from_singular_start,
        "iterations": result.iterations,
        "max_residual": result.max_residual,
        "solved": radius,
        "target": 2.5,
        "worst_abs_error": abs(radius - 2.5),
    }
    return cases


def _restarts(trials: int = 25, scale: float = 1e-3, seed: int = 20261112) -> dict:
    """Restart the SSS triangle from a *coincident* start plus a tiny jitter."""
    rng = random.Random(seed)
    converged = 0
    worst = 0.0
    failures: list[str] = []
    for _ in range(trials):
        jitter = [(rng.uniform(-scale, scale), rng.uniform(-scale, scale)) for _ in range(3)]
        sketch = Sketch()
        a = sketch.point(0.0, 0.0)
        b = sketch.point(*jitter[1])
        c = sketch.point(*jitter[2])
        sketch.add(Fix(a, 0.0, 0.0), Distance(b, a, 4.0), Distance(c, a, 3.0), Distance(c, b, 2.0))
        result = sketch.solve(strict=False)
        error = max(
            abs(math.dist(result.positions["P0"], result.positions["P1"]) - 4.0),
            abs(math.dist(result.positions["P0"], result.positions["P2"]) - 3.0),
            abs(math.dist(result.positions["P1"], result.positions["P2"]) - 2.0),
        )
        worst = max(worst, error)
        if result.converged:
            converged += 1
        else:
            failures.append(result.status)
    return {
        "trials": trials,
        "jitter_scale": scale,
        "seed": seed,
        "converged": converged,
        "rate": converged / trials,
        "worst_abs_error": worst,
        "failures": failures,
    }


def _status_contract() -> list[dict]:
    """An unsatisfiable set swept over the ``tol`` / ``residual_tol`` band."""
    bands = [
        (1e-12, 1e-9, 100, "defaults"),
        (1e-3, 1.0, 100, "wide residual_tol"),
        (1e-3, 10.0, 100, "wider residual_tol"),
        (1e-30, 1e3, 1, "tol unreachable, one iteration"),
    ]
    rows = []
    for tol, residual_tol, max_iterations, label in bands:
        sketch = Sketch()
        origin = sketch.point(0.0, 0.0)
        point = sketch.point(3.0, 0.0)
        sketch.add(Fix(origin, 0.0, 0.0), Fix(point, 3.0, 0.0), Distance(point, origin, 5.0))
        result = sketch.solve(strict=False, tol=tol, residual_tol=residual_tol,
                              max_iterations=max_iterations)
        rows.append(
            {
                "band": label,
                "tol": tol,
                "residual_tol": residual_tol,
                "max_residual": result.max_residual,
                "dof": result.dof,
                "converged": result.converged,
                "status": result.status,
                "fully_constrained": result.fully_constrained,
                "blocking": list(result.blocking),
                # the ladder this mission replaced: dof == 0 was reported as
                # "fully_constrained" with no convergence test at all.
                "previous_ladder_would_say": "fully_constrained" if result.dof == 0 else "under_constrained",
            }
        )
    return rows


def _plot(payload: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    names = [key for key, value in payload["cases"].items()]
    errors = [payload["cases"][key]["worst_abs_error"] for key in names]
    axes[0].bar(names, errors, color="#2b6cb0")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("|solved - analytic| [m]")
    axes[0].set_title(f"degenerate start, after a {SINGULAR_START_NUDGE:g} m escape")
    for index, value in enumerate(errors):
        axes[0].text(index, value, f"{value:.1e}", ha="center", va="bottom", fontsize=8)

    bands = [row["band"] for row in payload["status_contract"]]
    residuals = [row["max_residual"] for row in payload["status_contract"]]
    colours = ["#2f855a" if row["status"] == "not_converged" else "#c53030"
               for row in payload["status_contract"]]
    axes[1].bar(bands, residuals, color=colours)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("max |residual| [m]")
    axes[1].set_title("status returned for an unsatisfiable set")
    axes[1].tick_params(axis="x", labelrotation=20, labelsize=8)
    for index, row in enumerate(payload["status_contract"]):
        axes[1].text(index, row["max_residual"], row["status"], ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default=str(REPO_ROOT / "artifacts" / "sketch_robustness"))
    args = parser.parse_args()
    out = Path(args.out)
    payload = {
        "singular_start_nudge": SINGULAR_START_NUDGE,
        "cases": _cases(),
        "restarts": _restarts(),
        "status_contract": _status_contract(),
    }
    _plot(payload, out / "plot.png")
    out.mkdir(parents=True, exist_ok=True)
    (out / "results.json").write_text(json.dumps(payload, indent=2) + "\n")

    for name, case in payload["cases"].items():
        print(f"{name:16s} status={case['status']:18s} nudged={case['nudged_from_singular_start']!s:5s} "
              f"worst error={case['worst_abs_error']:.2e} m")
    restarts = payload["restarts"]
    print(f"restarts         {restarts['converged']}/{restarts['trials']} converged, "
          f"worst error {restarts['worst_abs_error']:.2e} m")
    for row in payload["status_contract"]:
        print(f"{row['band']:36s} residual={row['max_residual']:.2e} status={row['status']:16s} "
              f"(old ladder: {row['previous_ladder_would_say']})")
    print(f"figure           {out / 'plot.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
