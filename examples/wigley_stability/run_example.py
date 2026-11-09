"""Wigley initial stability measured against its own closed forms.

``tensorcad.wigley`` has always been a *verification case*, not just a hull form: its
volume, waterplane area, block/prismatic/midship coefficients and ``VCB`` are exact
rationals, so the mesh integrator is measured, not trusted. What it never had was a
closed form for the **waterplane second moments** -- and therefore none for ``BM``,
``KM`` or ``GM``. That is the one hydrostatic quantity a naval architect reads first,
so this script derives the references and measures the two independent mesh routes
against them:

* ``I_T = (B^3 L / 12) K(n)``, ``K(n) = 1 - 3/(n+1) + 3/(2n+1) - 1/(3n+1)``,
* ``I_L = (B L^3 / 4) (1/3 - 1/(n+3))``  (about the centre of flotation -- the
  parallel-axis term removed, matching what the mesh routes report),
* ``BM = I / nabla`` with ``nabla = L B T n m / ((n+1)(m+1))``, ``KB`` from the
  existing ``VCB`` closed form, ``GM = KB + BM - KG``.

Two mesh routes are compared, neither of which shares code with the other:

1. :func:`tensorcad.stability.waterplane_moments` / ``initial_stability`` -- the
   surface is cut with the plane ``z = 0`` and each closed loop integrated by Green's
   theorem;
2. :func:`tensorcad.hydrostatics.compute_hydrostatics` -- the same moment obtained by
   the divergence theorem over the *wetted surface* (no waterplane cap at all).

On a polygonal waterplane both are exact, so they must agree to round-off; the
remaining error is the mesh's own chord approximation of the analytic taper, which is
what the refinement ladder measures.

The script exits non-zero if one of its claims stops holding, so the tables in
``docs/wigley_stability_verification.md`` cannot drift away from the code.
matplotlib is optional: the JSON is written either way and the PNGs are skipped.

Run it with ``python examples/wigley_stability/run_example.py`` (writes
``artifacts/wigley_stability/``) or ``--output docs/wigley_stability``.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from tensorcad import __version__
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.stability import initial_stability, waterplane_moments
from tensorcad.wigley import (
    WigleyParameters,
    analytic_initial_stability,
    analytic_metacentric_radius_longitudinal,
    analytic_metacentric_radius_transverse,
    analytic_volume,
    analytic_waterplane_inertia,
    analytic_waterplane_inertia_longitudinal,
    generate_wigley,
)

#: ITTC reference proportions (L/B = 10, B/T = 1.6) -- the case every resistance and
#: diffraction benchmark quotes, so a user can reproduce it by name.
REFERENCE = WigleyParameters(
    length=10.0, beam=1.0, draft=0.625, nx=41, nz=13, longitudinal_exponent=2.0, vertical_exponent=2.0
)

#: KG above the keel used for the headline loading condition. Chosen so GM stays
#: positive but small: at ``B/T = 1.6`` the Wigley hull is stiff only near the keel.
KG = 0.5

LADDER = ((41, 13), (81, 25), (161, 41), (321, 81))

#: Closed-form cases: the classic Wigley, each exponent moved on its own, and a
#: fractional exponent -- which is the only one where Simpson is not integrating a
#: polynomial, so it is the case that actually exercises the quadrature.
EXPONENT_CASES = ((2.0, 2.0), (3.0, 2.0), (2.0, 3.0), (3.0, 3.0), (2.5, 2.0))

KG_LADDER = (0.3, 0.45, 0.5, 0.52, 0.527)


# ---------------------------------------------------------------------------
# independent quadrature route (stdlib only, shares nothing with the closed forms)
# ---------------------------------------------------------------------------
def simpson(values: list[float], step: float) -> float:
    """Composite Simpson on an even number of panels."""
    n = len(values) - 1
    if n < 2 or n % 2:
        raise ValueError("Simpson needs an even number of panels")
    total = values[0] + values[-1]
    total += 4.0 * math.fsum(values[1:-1:2])
    total += 2.0 * math.fsum(values[2:-1:2])
    return total * step / 3.0


def quadrature_inertias(parameters: WigleyParameters, panels: int = 20000) -> dict[str, float]:
    """``I_T`` and ``I_L`` by composite Simpson on the analytic half-breadth.

    The integrands are ``(2/3) b(x)^3`` and ``2 (x - L/2)^2 b(x)`` with
    ``b(x) = half_breadth(x, depth=0)`` -- the waterline of the analytic surface,
    evaluated independently of both mesh routes and of the closed forms. Simpson is
    fourth order, so with 20000 panels the quadrature error is far below the
    tolerance the test asks for.
    """
    from tensorcad.wigley import half_breadth

    length = parameters.length
    step = length / panels
    xs = [i * step for i in range(panels + 1)]
    breadths = [half_breadth(parameters, x, 0.0) for x in xs]
    transverse = simpson([(2.0 / 3.0) * b**3 for b in breadths], step)
    longitudinal = simpson(
        [2.0 * (x - 0.5 * length) ** 2 * b for x, b in zip(xs, breadths)], step
    )
    return {"transverse": transverse, "longitudinal": longitudinal}


# ---------------------------------------------------------------------------
# measurements
# ---------------------------------------------------------------------------
def measure(exponent_pair: tuple[float, float], grid: tuple[int, int]) -> dict[str, object]:
    """One mesh, two routes, one loading condition."""
    n, m = exponent_pair
    nx, nz = grid
    parameters = WigleyParameters(
        length=REFERENCE.length,
        beam=REFERENCE.beam,
        draft=REFERENCE.draft,
        nx=nx,
        nz=nz,
        longitudinal_exponent=n,
        vertical_exponent=m,
    )
    mesh = generate_wigley(parameters).mesh
    reference = analytic_initial_stability(parameters, kg=KG)
    route_cut = initial_stability(mesh, kg=KG)
    route_flux = compute_hydrostatics(mesh, water_density=parameters.water_density)
    exactly = waterplane_moments(mesh, waterline_z=0.0)

    def rel(measured: float, exact: float) -> float:
        return (measured - exact) / abs(exact)

    return {
        "grid": f"{nx}x{nz}",
        "nx": nx,
        "nz": nz,
        "longitudinal_exponent": n,
        "vertical_exponent": m,
        "volume": route_flux.displacement_volume,
        "reference_volume": reference.displacement_volume,
        "volume_error": rel(route_flux.displacement_volume, reference.displacement_volume),
        "waterplane_area": exactly.area,
        "reference_waterplane_area": None,  # filled per exponent below
        "inertia_transverse": route_cut.waterplane.transverse_second_moment,
        "inertia_longitudinal": route_cut.waterplane.longitudinal_second_moment,
        "inertia_transverse_flux": route_flux.transverse_waterplane_inertia,
        "inertia_longitudinal_flux": route_flux.longitudinal_waterplane_inertia,
        "reference_inertia_transverse": reference.transverse_waterplane_inertia,
        "reference_inertia_longitudinal": reference.longitudinal_waterplane_inertia,
        "inertia_transverse_error": rel(
            route_cut.waterplane.transverse_second_moment, reference.transverse_waterplane_inertia
        ),
        "inertia_longitudinal_error": rel(
            route_cut.waterplane.longitudinal_second_moment,
            reference.longitudinal_waterplane_inertia,
        ),
        "route_agreement_transverse": abs(
            route_cut.waterplane.transverse_second_moment - route_flux.transverse_waterplane_inertia
        )
        / reference.transverse_waterplane_inertia,
        "route_agreement_longitudinal": abs(
            route_cut.waterplane.longitudinal_second_moment
            - route_flux.longitudinal_waterplane_inertia
        )
        / reference.longitudinal_waterplane_inertia,
        "kb": route_cut.kb,
        "bm_transverse": route_cut.bm_transverse,
        "km_transverse": route_cut.km_transverse,
        "gm_transverse": route_cut.gm_transverse,
        "bm_longitudinal": route_cut.bm_longitudinal,
        "km_longitudinal": route_cut.km_longitudinal,
        "gm_longitudinal": route_cut.gm_longitudinal,
        "reference": reference.as_dict(),
        "error": reference.compare(route_cut)["relative_error"],
    }


def fit_order(rungs: list[dict[str, object]], key: str) -> float:
    """Least-squares slope of ``log|error|`` against ``log(1/(nx-1))``."""
    points = [
        (math.log(1.0 / (rung["nx"] - 1)), math.log(abs(rung[key])))
        for rung in rungs
        if abs(rung[key]) > 0.0
    ]
    count = len(points)
    mean_x = math.fsum(p[0] for p in points) / count
    mean_y = math.fsum(p[1] for p in points) / count
    numerator = math.fsum((x - mean_x) * (y - mean_y) for x, y in points)
    denominator = math.fsum((x - mean_x) ** 2 for x, _ in points)
    return numerator / denominator


def kg_ladder() -> list[dict[str, float]]:
    """GM error against KG on the default grid: the amplification, measured.

    ``GM = KM - KG`` is a difference of nearly equal numbers when ``KG`` approaches
    ``KM``, so the same ``KM`` error becomes a much larger ''GM'' error. Reported as a
    table rather than a warning, because how close a loading condition may sit to the
    limit is a design decision, not a numerical one.
    """
    parameters = WigleyParameters(
        length=REFERENCE.length,
        beam=REFERENCE.beam,
        draft=REFERENCE.draft,
        nx=LADDER[0][0],
        nz=LADDER[0][1],
        longitudinal_exponent=2.0,
        vertical_exponent=2.0,
    )
    mesh = generate_wigley(parameters).mesh
    rows: list[dict[str, float]] = []
    for kg in KG_LADDER:
        reference = analytic_initial_stability(parameters, kg=kg)
        measured = initial_stability(mesh, kg=kg)
        rows.append(
            {
                "kg": kg,
                "km": reference.km_transverse,
                "gm_reference": reference.gm_transverse,
                "gm_measured": measured.gm_transverse,
                "gm_error": (measured.gm_transverse - reference.gm_transverse)
                / abs(reference.gm_transverse),
                "km_error": (measured.km_transverse - reference.km_transverse)
                / reference.km_transverse,
                "amplification_km_over_gm": reference.km_transverse / reference.gm_transverse,
            }
        )
    return rows


def closed_form_table() -> list[dict[str, object]]:
    """Closed forms for four exponent pairs, plus their hand-derived ``n = 2`` values."""
    rows: list[dict[str, object]] = []
    for n, m in EXPONENT_CASES:
        parameters = WigleyParameters(
            length=REFERENCE.length,
            beam=REFERENCE.beam,
            draft=REFERENCE.draft,
            longitudinal_exponent=n,
            vertical_exponent=m,
        )
        reference = analytic_initial_stability(parameters, kg=KG)
        quadrature = quadrature_inertias(parameters)
        rows.append(
            {
                "longitudinal_exponent": n,
                "vertical_exponent": m,
                "volume": reference.displacement_volume,
                "inertia_transverse": reference.transverse_waterplane_inertia,
                "inertia_longitudinal": reference.longitudinal_waterplane_inertia,
                "inertia_transverse_quadrature": quadrature["transverse"],
                "inertia_longitudinal_quadrature": quadrature["longitudinal"],
                "quadrature_error_transverse": (
                    quadrature["transverse"] - reference.transverse_waterplane_inertia
                )
                / reference.transverse_waterplane_inertia,
                "quadrature_error_longitudinal": (
                    quadrature["longitudinal"] - reference.longitudinal_waterplane_inertia
                )
                / reference.longitudinal_waterplane_inertia,
                "kb": reference.kb,
                "bm_transverse": reference.bm_transverse,
                "bm_longitudinal": reference.bm_longitudinal,
                "km_transverse": reference.km_transverse,
                "km_longitudinal": reference.km_longitudinal,
                "gm_transverse": reference.gm_transverse,
                "gm_longitudinal": reference.gm_longitudinal,
            }
        )
    return rows


# ---------------------------------------------------------------------------
# claims
# ---------------------------------------------------------------------------
def evaluate(results: dict[str, object]) -> list[dict[str, object]]:
    claims: list[dict[str, object]] = []

    def add(claim: str, value: float, threshold: float) -> None:
        """``|value| <= threshold``."""
        claims.append(
            {
                "claim": claim,
                "value": float(value),
                "threshold": float(threshold),
                "holds": abs(value) <= threshold,
            }
        )

    def add_at_least(claim: str, value: float, threshold: float) -> None:
        """``value >= threshold`` (used where the sign is the point)."""
        claims.append(
            {
                "claim": claim,
                "value": float(value),
                "threshold": float(threshold),
                "holds": float(value) >= threshold,
            }
        )

    def add_negative(claim: str, value: float) -> None:
        """``value < 0`` on every sample: the sign of a mesh error is a prediction."""
        claims.append(
            {"claim": claim, "value": float(value), "threshold": 0.0, "holds": float(value) < 0.0}
        )

    closed = results["closed_forms"]
    worst_quadrature = max(
        max(abs(row["quadrature_error_transverse"]), abs(row["quadrature_error_longitudinal"]))
        for row in closed
    )
    add(
        "closed forms agree with independent Simpson quadrature of the analytic waterline "
        "(5 exponent pairs, 20000 panels)",
        worst_quadrature,
        1e-8,
    )

    n2 = next(row for row in closed if row["longitudinal_exponent"] == 2.0 and row["vertical_exponent"] == 2.0)
    L, B, T = REFERENCE.length, REFERENCE.beam, REFERENCE.draft
    hand = {
        "inertia_transverse": 4.0 * L * B**3 / 105.0,
        "inertia_longitudinal": B * L**3 / 30.0,
        "bm_transverse": 3.0 * B**2 / (35.0 * T),
        "bm_longitudinal": 3.0 * L**2 / (40.0 * T),
    }
    worst_hand = max(
        abs(n2[key] - value) / value for key, value in hand.items()
    )
    add(
        "n = 2 closed forms reproduce the hand rationals 4LB^3/105, BL^3/30, 3B^2/(35T), 3L^2/(40T)",
        worst_hand,
        1e-15,
    )

    for label, exponent_pair in (("n=2", (2.0, 2.0)), ("n=3", (3.0, 2.0))):
        rungs = results["ladders"][label]
        default = rungs[0]
        add(
            f"{label}: default grid (41x13) I_T within 3 % of the closed form",
            default["inertia_transverse_error"],
            0.03,
        )
        add(
            f"{label}: default grid (41x13) GM_T within 3 % of the closed form",
            default["error"]["gm_transverse"],
            0.03,
        )
        fine = rungs[-1]
        add(
            f"{label}: fine grid (321x81) I_T within 0.01 % of the closed form",
            fine["inertia_transverse_error"],
            1e-4,
        )
        add(
            f"{label}: fine grid (321x81) GM_T within 0.2 % of the closed form",
            fine["error"]["gm_transverse"],
            2e-3,
        )
        add_negative(
            f"{label}: I_T is under-predicted on every rung (waterline chords cut the convex taper)",
            max(rung["inertia_transverse_error"] for rung in rungs),
        )
        order = fit_order(rungs, "inertia_transverse_error")
        results.setdefault("fitted_orders", {})[f"{label}_inertia_transverse"] = order
        add(f"{label}: I_T error is second order in the grid step (fitted order 1.8-2.2)", order - 2.0, 0.2)
        worst_route = max(
            max(rung["route_agreement_transverse"], rung["route_agreement_longitudinal"])
            for rung in rungs
        )
        add(
            f"{label}: the plane-cut and wetted-flux routes agree to round-off on every rung",
            worst_route,
            1e-12,
        )

    ladder = results["kg_ladder"]
    middle = ladder[2]  # KG = KG (0.5)
    amplifications = [
        abs(row["gm_error"]) / abs(row["km_error"]) for row in ladder
    ]
    add_at_least(
        "at KG = 0.5 the measured GM error exceeds 10x the KM error (KM/GM amplification)",
        amplifications[2],
        10.0,
    )
    add_at_least(
        "at KG = 0.527 (GM = 0.77 mm) the amplification exceeds 100x",
        amplifications[-1],
        100.0,
    )
    increasing = all(
        later > earlier for earlier, later in zip(amplifications, amplifications[1:])
    )
    claims.append(
        {
            "claim": "the amplification grows monotonically as GM -> 0 (KG ladder)",
            "value": amplifications[-1] - amplifications[0],
            "threshold": 0.0,
            "holds": increasing,
        }
    )
    add(
        "at 321x81 the whole KM chain (I_T, I_L, KB, BM, KM) is within 1 %",
        max(
            max(abs(value) for value in rung["error"].values())
            for rung in results["ladders"]["n=2"]
            if rung["grid"] == "321x81"
        ),
        0.01,
    )
    return claims


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def figure_convergence(results: dict[str, object], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))
    series = (
        ("inertia_transverse_error", "I_T", "o"),
        ("inertia_longitudinal_error", "I_L", "s"),
    )
    for axis, label in zip(axes, ("n = 2", "n = 3")):
        key = "n=2" if label == "n = 2" else "n=3"
        rungs = results["ladders"][key]
        steps = [1.0 / (rung["nx"] - 1) for rung in rungs]
        for field, name, marker in series:
            axis.loglog(
                steps,
                [abs(rung[field]) for rung in rungs],
                marker + "-",
                label=f"{name} (plane cut)",
            )
        axis.loglog(
            steps,
            [abs(rung["error"]["gm_transverse"]) for rung in rungs],
            "^-",
            label="GM_T (KB+BM-KG)",
        )
        guide = [abs(rungs[0]["inertia_transverse_error"]) * (s / steps[0]) ** 2 for s in steps]
        axis.loglog(steps, guide, "k--", linewidth=1.0, label="slope 2 guide")
        axis.set_title(f"Wigley {label} (L/B=10, B/T=1.6)")
        axis.set_xlabel("grid step 1/(nx-1)")
        axis.set_ylabel("relative error vs closed form")
        axis.grid(True, which="both", alpha=0.3)
        axis.legend(fontsize=8)
    figure.suptitle("Mesh routes converge to the closed forms at second order")
    figure.tight_layout()
    figure.savefig(path, dpi=140)
    plt.close(figure)


def figure_amplification(results: dict[str, object], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ladder = results["kg_ladder"]
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))
    axes[0].plot([row["kg"] for row in ladder], [row["gm_reference"] for row in ladder], "o-")
    axes[0].axhline(0.0, color="k", linewidth=0.8)
    axes[0].set_xlabel("KG above keel [m]")
    axes[0].set_ylabel("GM_T [m]")
    axes[0].set_title("GM falls one-for-one with KG (KM = 0.527768 m)")
    axes[0].grid(alpha=0.3)

    axes[1].plot(
        [row["kg"] for row in ladder],
        [row["amplification_km_over_gm"] for row in ladder],
        "s-",
        label="KM/GM (predicted amplification)",
    )
    axes[1].plot(
        [row["kg"] for row in ladder],
        [abs(row["gm_error"]) / abs(row["km_error"]) for row in ladder],
        "o--",
        label="measured rel(GM)/rel(KM)",
    )
    axes[1].set_xlabel("KG above keel [m]")
    axes[1].set_ylabel("amplification")
    axes[1].set_yscale("log")
    axes[1].set_title("default 41x13 grid: GM error = KM error x KM/GM")
    axes[1].grid(alpha=0.3, which="both")
    axes[1].legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(path, dpi=140)
    plt.close(figure)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def render_report(results: dict[str, object], claims: list[dict[str, object]]) -> str:
    lines = [
        "# Wigley initial stability: closed forms and the two mesh routes",
        "",
        f"TensorCAD {results['version']} · reference case L = {REFERENCE.length} m, "
        f"B = {REFERENCE.beam} m, T = {REFERENCE.draft} m · KG = {KG} m above the keel",
        "",
        "## Closed forms",
        "",
        "| n | m | I_T [m^4] | I_L [m^4] | KB [m] | BM_T [m] | BM_L [m] | KM_T [m] | GM_T [m] |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for row in results["closed_forms"]:
        lines.append(
            "| {:.0f} | {:.0f} | {:.9f} | {:.6f} | {:.6f} | {:.9f} | {:.6f} | {:.9f} | {:.6f} |".format(
                row["longitudinal_exponent"],
                row["vertical_exponent"],
                row["inertia_transverse"],
                row["inertia_longitudinal"],
                row["kb"],
                row["bm_transverse"],
                row["bm_longitudinal"],
                row["km_transverse"],
                row["gm_transverse"],
            )
        )
    for label in ("n=2", "n=3"):
        lines += [
            "",
            f"## Refinement ladder, {label.replace('n=', 'n = ')}",
            "",
            "| grid | I_T [m^4] | err I_T | err I_L | err KB | err BM_T | err KM_T | err GM_T | routes |",
            "|---|---|---|---|---|---|---|---|---|",
        ]
        for rung in results["ladders"][label]:
            lines.append(
                "| {} | {:.9f} | {:+.4f} % | {:+.4f} % | {:+.4f} % | {:+.4f} % | {:+.4f} % | {:+.4f} % | {:.1e} |".format(
                    rung["grid"],
                    rung["inertia_transverse"],
                    100.0 * rung["inertia_transverse_error"],
                    100.0 * rung["inertia_longitudinal_error"],
                    100.0 * rung["error"]["kb"],
                    100.0 * rung["error"]["bm_transverse"],
                    100.0 * rung["error"]["km_transverse"],
                    100.0 * rung["error"]["gm_transverse"],
                    max(rung["route_agreement_transverse"], rung["route_agreement_longitudinal"]),
                )
            )
        lines.append(
            f"\nFitted order of the ``I_T`` error: **{results['fitted_orders'][label + '_inertia_transverse']:.2f}**."
        )
    lines += [
        "",
        "## GM error against the loading condition (default 41x13 grid)",
        "",
        "| KG [m] | GM ref [m] | GM mesh [m] | rel GM | rel KM | KM/GM |",
        "|---|---|---|---|---|---|",
    ]
    for row in results["kg_ladder"]:
        lines.append(
            "| {:.3f} | {:+.6f} | {:+.6f} | {:+.3f} % | {:+.4f} % | {:.1f} |".format(
                row["kg"],
                row["gm_reference"],
                row["gm_measured"],
                100.0 * row["gm_error"],
                100.0 * row["km_error"],
                row["amplification_km_over_gm"],
            )
        )
    lines += ["", "## Claims", "", "| claim | value | threshold | holds |", "|---|---|---|---|"]
    for claim in claims:
        lines.append(
            "| {} | {:.3e} | {:.0e} | {} |".format(
                claim["claim"], claim["value"], claim["threshold"], "PASS" if claim["holds"] else "FAIL"
            )
        )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/wigley_stability"),
        help="directory for results.json / report.md / figures",
    )
    parser.add_argument("--no-figures", action="store_true", help="skip the PNGs")
    arguments = parser.parse_args(argv)
    output: Path = arguments.output
    output.mkdir(parents=True, exist_ok=True)

    results: dict[str, object] = {
        "version": __version__,
        "reference_case": REFERENCE.as_dict(),
        "kg": KG,
        "closed_forms": closed_form_table(),
        "ladders": {
            "n=2": [measure((2.0, 2.0), grid) for grid in LADDER],
            "n=3": [measure((3.0, 2.0), grid) for grid in LADDER],
        },
        "kg_ladder": kg_ladder(),
    }
    claims = evaluate(results)
    results["claims"] = claims
    (output / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True))
    (output / "report.md").write_text(render_report(results, claims))
    if not arguments.no_figures:
        try:
            figure_convergence(results, output / "fig_convergence.png")
            figure_amplification(results, output / "fig_amplification.png")
        except ImportError:  # pragma: no cover - matplotlib is optional
            print("[wigley-stability] matplotlib missing: figures skipped", file=sys.stderr)

    failed = [claim for claim in claims if not claim["holds"]]
    for claim in claims:
        status = "PASS" if claim["holds"] else "FAIL"
        print(f"[wigley-stability] {status}  {claim['claim']}  ({claim['value']:.3e} vs {claim['threshold']:.0e})")
    if failed:
        print(f"[wigley-stability] {len(failed)} claim(s) do not hold", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())