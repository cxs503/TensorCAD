#!/usr/bin/env python3
"""Three-grid convergence study of the golden cases' analytic metrics.

Question this tool answers: is the golden acceptance band (currently ``2 x |fine -
coarse|``) the tightest band the measurements can justify, or does a third grid level
buy a tighter *derived* band?

For each analytic golden case (``wigley``, ``series60``, ``sphere_ellipsoid``) the same
metrics the gate measures are evaluated on a ladder of grid levels and compared with the
case's own reference (a closed form / converged quadrature, not a stored number).  From
the observed errors ``e(h) ~ C h^p`` the tool reports, per metric:

* the relative error at every level;
* the observed order ``p`` fitted from three consecutive levels (non-uniform refinement
  ratios handled by solving the ratio equation numerically);
* the Richardson estimate of the *remaining* error at the acceptance level;
* the band a Richardson estimate would give (``2 x |R|``, floored) next to the band the
  two-grid rule gives today, and the baseline's own error (a band tighter than that is
  a gate that is red on the honest tree).

Usage::

    python tools/golden/convergence_study.py --json docs/golden/convergence_study.json \
        --md docs/golden/convergence_study.md
"""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
from typing import Any, Callable

import tensorcad.golden as G

SCHEMA = "tensorcad.golden.convergence/1"


def _rel(value: float, reference: float) -> float:
    return (value - reference) / reference if reference else value - reference


# --------------------------------------------------------------------------------------
# cases: (level name, grid) ladders.  ``accept`` marks the level the gate asserts on and
# ``coarse`` the level the two-grid rule compares it with, so the ladder contains the
# gate's own pair plus one finer level and (for the order fit) one coarser level.
# --------------------------------------------------------------------------------------
def _wigley(level: tuple[int, int]) -> tuple[dict[str, float], dict[str, float]]:
    from tensorcad import wigley

    nx, nz = level
    hydro = wigley.generate_wigley(wigley.reference_parameters(nx=nx, nz=nz)).hydrostatics
    return (
        G._hydrostatics_values(hydro, G.WigleyGolden.METRICS),
        G.WigleyGolden().reference(),
    )


def _series60(level: tuple[int, int]) -> tuple[dict[str, float], dict[str, float]]:
    from tensorcad import series60

    nx, nz = level
    hydro = series60.generate_series60(series60.Series60Parameters(nx=nx, nz=nz)).hydrostatics
    return (
        G._hydrostatics_values(hydro, G.Series60Golden.METRICS),
        G.Series60Golden().reference(),
    )


def _sphere_ellipsoid(resolution: int) -> tuple[dict[str, float], dict[str, float]]:
    """The SDF case's own measurement, at an arbitrary iso-surface resolution."""
    from tensorcad.geometry_field import EllipsoidField, SphereField, to_mesh
    from tensorcad.hydrostatics import compute_hydrostatics

    case = G.SphereEllipsoidGolden()
    values: dict[str, float] = {}
    for name, geometry in (
        ("sphere", SphereField(radius=case.RADIUS, center=(0.0, 0.0, 0.0))),
        ("ellipsoid", EllipsoidField(semi_axes=case.AXES, center=(0.0, 0.0, 0.0))),
    ):
        mesh = to_mesh(geometry, resolution=resolution).mesh
        hydro = compute_hydrostatics(mesh, waterline_z=-1e3, closed=True)
        values[f"{name}_volume_m3"] = abs(hydro.displacement_volume)
        values[f"{name}_surface_area_m2"] = G._closed_surface_area(mesh)
        for axis, extent in zip("xyz", G._extents(mesh)):
            values[f"{name}_extent_{axis}_m"] = extent
    return values, case.reference()


CASES: dict[str, dict[str, Any]] = {
    # the gate's own coarse/fine pair, plus one finer level ("finer") for the order fit
    "wigley": {
        "fn": _wigley,
        "levels": [("coarse", (41, 13)), ("fine", (161, 41)), ("finer", (321, 81))],
        "coarse": "coarse",
        "accept": "fine",
        "grid_of": lambda level: 1.0 / math.sqrt(level[0] * level[1]),
    },
    "series60": {
        "fn": _series60,
        "levels": [("coarse", (81, 25)), ("fine", (161, 49)), ("finer", (321, 97))],
        "coarse": "coarse",
        "accept": "fine",
        "grid_of": lambda level: 1.0 / math.sqrt(level[0] * level[1]),
    },
    "sphere_ellipsoid": {
        "fn": _sphere_ellipsoid,
        "levels": [("coarse", 24), ("fine", 48), ("finer", 96)],
        "coarse": "coarse",
        "accept": "fine",
        "grid_of": lambda level: 1.0 / level,
    },
}


def _solve_order(q: float, r1: float, r2: float) -> float | None:
    """Solve ``(r1^p - 1) / (1 - r2^-p) = q`` for the observed order ``p``.

    ``q`` is the ratio of the two consecutive error changes, ``r1`` the refinement ratio
    between the coarser pair and ``r2`` the one between the finer pair.  ``q`` is
    monotone in ``p``, so a bisection is enough; a pair that does not move at all (or moves
    the wrong way) has no positive order and the caller falls back to the two-grid rule.
    """
    if not (q > 0.0):
        return None

    def residual(p: float) -> float:
        return (r1**p - 1.0) / (1.0 - r2 ** (-p)) - q

    low, high = 1e-3, 12.0
    if residual(low) > 0.0 or residual(high) < 0.0:
        return None
    for _ in range(200):
        mid = 0.5 * (low + high)
        if residual(mid) > 0.0:
            high = mid
        else:
            low = mid
    return 0.5 * (low + high)


def study(max_metrics: int | None = None) -> dict[str, Any]:
    report: dict[str, Any] = {"schema": SCHEMA, "cases": {}, "verdict": {}}
    for name, spec in CASES.items():
        rows: list[dict[str, Any]] = []
        table: dict[str, dict[str, float]] = {}
        reference: dict[str, float] = {}
        started = time.time()
        for label, level in spec["levels"]:
            t0 = time.time()
            values, ref = spec["fn"](level)
            reference = ref
            table[label] = values
            rows.append({"label": label, "level": list(level) if isinstance(level, tuple) else level,
                         "seconds": round(time.time() - t0, 2)})
            print(f"{name} {label} {level} {rows[-1]['seconds']}s", flush=True)
        elapsed = round(time.time() - started, 2)
        metrics: dict[str, Any] = {}
        for key, ref_value in reference.items():
            errors = {label: _rel(table[label][key], ref_value) for label, _ in spec["levels"] if key in table[label]}
            entry: dict[str, Any] = {"errors": errors}
            labels = [label for label, _ in spec["levels"] if label in errors]
            coarse, accept, finer = spec["coarse"], spec["accept"], "finer"
            if all(label in errors for label in (coarse, accept, finer)):
                h = {label: spec["grid_of"](dict(spec["levels"])[label]) for label in labels}
                e_c, e_a, e_f = abs(errors[coarse]), abs(errors[accept]), abs(errors[finer])
                q = abs(e_c - e_a) / abs(e_a - e_f) if abs(e_a - e_f) > 0 else float("inf")
                r1 = h[coarse] / h[accept]
                r2 = h[accept] / h[finer]
                p = _solve_order(q, r1, r2) if math.isfinite(q) else None
                entry["change_two_grid"] = abs(errors[accept] - errors[coarse])
                entry["change_three_grid"] = abs(errors[finer] - errors[accept])
                entry["order"] = None if p is None else round(p, 3)
                if p is not None:
                    # remaining error at the acceptance level: |e_accept - e_finer| / (1 - r2^-p)
                    remaining = abs(errors[accept] - errors[finer]) / (1.0 - r2 ** (-p))
                    entry["richardson_remaining"] = remaining
                    entry["richardson_band"] = max(2.0 * remaining, G.TOLERANCE_FLOOR)
                two_grid = max(
                    2.0 * abs(errors[accept] - errors[coarse]),
                    G.TOLERANCE_FLOOR if abs(errors[accept] - errors[coarse]) <= G.TOLERANCE_FLOOR / 2 else 0.0,
                )
                entry["two_grid_band"] = two_grid
                entry["observed_error_at_accept"] = abs(errors[accept])
                if "richardson_band" in entry:
                    entry["band_ratio"] = entry["richardson_band"] / two_grid if two_grid else None
                    entry["sound"] = entry["richardson_band"] > abs(errors[accept])
            metrics[key] = entry
            if max_metrics is not None and len(metrics) >= max_metrics:
                break
        report["cases"][name] = {"rows": rows, "seconds": elapsed, "metrics": metrics}
    return report


def _markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Golden acceptance band: does a third grid level tighten it?",
        "",
        f"- schema: `{SCHEMA}`",
        "- reference = each case's own closed form / converged quadrature (no stored number)",
        "- `two-grid band` = the rule in force: 2 x |fine - coarse|; `Richardson band` = 2 x |remaining error at fine|",
        "- `sound` = the Richardson band still admits the honest fine measurement (else it would turn the gate red on a correct tree)",
        "",
    ]
    for name, case in report["cases"].items():
        lines += [f"## {name} ({case['seconds']} s)", "",
                  "| metric | err coarse | err fine | order p | two-grid band | Richardson band | honest err | sound |",
                  "| --- | --- | --- | --- | --- | --- | --- | --- |"]
        for key, m in case["metrics"].items():
            errs = m["errors"]
            lines.append(
                "| `{}` | {:+.2e} | {:+.2e} | {} | {} | {} | {} | {} |".format(
                    key,
                    errs.get("coarse", float("nan")),
                    errs.get("fine", float("nan")),
                    m.get("order", "n/a"),
                    f"{m['two_grid_band']:.3e}" if "two_grid_band" in m else "n/a",
                    f"{m['richardson_band']:.3e}" if "richardson_band" in m else "n/a",
                    f"{m['observed_error_at_accept']:.3e}" if "observed_error_at_accept" in m else "n/a",
                    m.get("sound", "n/a"),
                )
            )
        lines.append("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--json", type=Path, default=None)
    parser.add_argument("--md", type=Path, default=None)
    parser.add_argument("--max-metrics", type=int, default=None)
    args = parser.parse_args(argv)

    report = study(args.max_metrics)
    text = json.dumps(report, indent=1, sort_keys=True)
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text + "\n", encoding="utf-8")
    else:
        print(text)
    if args.md:
        args.md.write_text(_markdown(report) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())