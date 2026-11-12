#!/usr/bin/env python3
"""Adversarial audit of `tensorcad.savitsky` -- evidence generator (m16954).

Run:  python examples/savitsky_audit/run_example.py

Writes ``artifacts/savitsky_audit/results.json``, ``report.md`` and four figures.
Every claim is asserted; the script exits non-zero if a claim stops holding, so the
numbers in ``docs/savitsky_audit_m16954.md`` cannot drift away from the code.

What is measured
----------------
1. Does the trim solve produce a wetted length ratio at all?  (finding F1)
2. Is the lift equation satisfied exactly, and is an out-of-window case reported as
   such instead of as "non-convergence"?                       (finding F2/F5a)
3. Does the hull-length cap preserve lift equilibrium?         (finding F3)
4. How far apart are the module's lift correlation and the classical Cv-dependent
   Savitsky correlation?                                       (finding F6)
5. Does the new two-equation closure actually close?           (new capability)
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (backend must be set first)

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.savitsky import (  # noqa: E402
    DEFAULT_WETTED_LENGTH_RATIO,
    PlaningPredictionError,
    SEA_WATER_DENSITY,
    center_of_pressure_ratio,
    classic_lift_coefficient,
    equilibrium_trim,
    lift_coefficient,
    planing_prediction,
    running_trim,
    speed_coefficient,
    speed_sweep,
    trim_for_wetted_length,
)

OUT = Path(__file__).resolve().parents[2] / "artifacts" / "savitsky_audit"

# Reference craft: the same family the module's own example uses.
WEIGHT = 15000.0      # N
BEAM = 2.5            # m
LENGTH = 12.0         # m
DEADRISE = 20.0       # deg
LCG = 3.0             # m forward of the transom

claims: list[dict] = []


def claim(name: str, expected: str, actual, error, passed: bool):
    claims.append(
        {
            "name": name,
            "expected": expected,
            "actual": actual,
            "error": error,
            "pass": bool(passed),
        }
    )


# ---------------------------------------------------------------- 1. finding F1
def audit_wetted_length() -> list[dict]:
    rows = []
    for weight, speed, beam, deadrise in [
        (8000.0, 8.0, 2.5, 20.0),
        (12000.0, 10.0, 2.0, 10.0),
        (15000.0, 12.0, 2.5, 20.0),
        (20000.0, 15.0, 3.0, 25.0),
        (26000.0, 9.0, 2.2, 30.0),
        (40000.0, 18.0, 3.0, 15.0),
    ]:
        trim, lam = running_trim(weight, speed, beam, deadrise)
        rows.append(
            {
                "weight_n": weight,
                "speed_ms": speed,
                "beam_m": beam,
                "deadrise_deg": deadrise,
                "trim_deg": trim,
                "wetted_length_ratio_returned": lam,
            }
        )
    lams = {row["wetted_length_ratio_returned"] for row in rows}
    claim(
        "F1 the legacy trim solve returns one single lambda for every input",
        f"one value: {DEFAULT_WETTED_LENGTH_RATIO}",
        f"{len(lams)} distinct value(s)",
        None,
        lams == {DEFAULT_WETTED_LENGTH_RATIO},
    )
    return rows


# ---------------------------------------------------------- 2. closed-form trim
def audit_closed_form() -> dict:
    lam = 3.0
    trim = trim_for_wetted_length(WEIGHT, 12.0, BEAM, DEADRISE, lam)
    lift = lift_coefficient(trim, lam, DEADRISE) * (
        0.5 * SEA_WATER_DENSITY * 12.0 ** 2 * BEAM ** 2
    )
    relative = abs(lift - WEIGHT) / WEIGHT
    claim(
        "F2 closed-form trim satisfies the lift equation",
        "|CL q B^2 - W| / W < 1e-12",
        f"{relative:.3e}",
        relative,
        relative < 1e-12,
    )

    window_message = ""
    try:
        trim_for_wetted_length(40000.0, 6.0, 2.0, DEADRISE, lam)
    except PlaningPredictionError as exc:
        window_message = str(exc)
    honest = "required trim" in window_message and "failed to converge" not in window_message
    claim(
        "F2/F5a out-of-window trim is reported as a window violation",
        "message names the required trim, not a convergence failure",
        window_message[:90],
        None,
        honest,
    )
    return {
        "trim_deg": trim,
        "lift_n": lift,
        "relative_lift_error": relative,
        "out_of_window_message": window_message,
    }


# --------------------------------------------------------- 3. clamp honesty F3
def audit_clamp() -> dict:
    pred = planing_prediction(WEIGHT, 12.0, 5.0, BEAM, DEADRISE)  # L/B = 2
    lift = pred.lift_coefficient * (
        0.5 * SEA_WATER_DENSITY * 12.0 ** 2 * BEAM ** 2
    )
    relative = abs(lift - pred.weight_n) / pred.weight_n
    claim(
        "F3 after the hull-length cap the lift balance still holds",
        "|CL q B^2 - W| / W < 1e-9",
        f"{relative:.3e}",
        relative,
        relative < 1e-9,
    )
    claim(
        "F3 the cap is reported in warnings",
        "at least one warning mentions 'clamped'",
        str([w for w in pred.warnings if "clamped" in w])[:70],
        None,
        any("clamped" in w for w in pred.warnings),
    )
    return {
        "hull_length_m": 5.0,
        "wetted_length_ratio": pred.wetted_length_ratio,
        "lift_n": lift,
        "weight_n": pred.weight_n,
        "relative_lift_error": relative,
        "warnings": list(pred.warnings),
    }


# ------------------------------------------------ 4. correlation gap (F6)
def audit_lift_models() -> dict:
    grid = [
        (tau, lam, beta)
        for tau in (2.0, 3.0, 4.0, 6.0, 8.0, 10.0)
        for lam in (1.0, 2.0, 3.0, 4.0)
        for beta in (10.0, 20.0, 30.0)
    ]

    def worst(cv: float) -> float:
        return max(
            abs(lift_coefficient(t, lam, b) - classic_lift_coefficient(t, lam, cv, b))
            / classic_lift_coefficient(t, lam, cv, b)
            for t, lam, b in grid
        )

    scan = [(0.5 + 0.05 * i, worst(0.5 + 0.05 * i)) for i in range(91)]
    best_cv, best_worst = min(scan, key=lambda item: item[1])
    claim(
        "F6 no single speed coefficient reconciles the two lift correlations",
        "best-fit worst-case disagreement > 50%",
        f"{100 * best_worst:.1f}% at Cv={best_cv:.2f}",
        best_worst,
        best_worst > 0.5,
    )

    rows = []
    for tau in (2.0, 4.0, 6.0, 8.0, 10.0):
        for lam in (2.0, 3.0, 4.0):
            module = lift_coefficient(tau, lam, DEADRISE)
            classic = classic_lift_coefficient(tau, lam, 2.0, DEADRISE)
            rows.append(
                {
                    "trim_deg": tau,
                    "wetted_length_ratio": lam,
                    "module_cl": module,
                    "classic_cl_cv2": classic,
                    "relative_difference": (module - classic) / classic,
                }
            )
    return {
        "grid_size": len(grid),
        "best_fit_cv": best_cv,
        "best_fit_worst_case_relative_difference": best_worst,
        "table_beta20": rows,
    }


# ------------------------------------------- 5. equilibrium closure (new path)
def audit_equilibrium() -> tuple[dict, list[dict]]:
    speeds = [8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0]
    rows = []
    for speed in speeds:
        solution = equilibrium_trim(
            WEIGHT, speed, BEAM, DEADRISE, LCG, max_wetted_length_m=LENGTH
        )
        legacy = planing_prediction(WEIGHT, speed, LENGTH, BEAM, DEADRISE)
        solved = planing_prediction(
            WEIGHT, speed, LENGTH, BEAM, DEADRISE, lcg_m=LCG
        )
        rows.append(
            {
                "speed_ms": speed,
                "froude_number": solved.froude_number,
                "speed_coefficient": solution.speed_coefficient,
                "trim_deg": solution.trim_deg,
                "wetted_length_ratio": solution.wetted_length_ratio,
                "lift_residual": solution.lift_residual,
                "moment_residual_m": solution.moment_residual_m,
                "iterations": solution.iterations,
                "resistance_n": solved.total_drag_n,
                "legacy_trim_deg": legacy.trim_deg,
                "legacy_resistance_n": legacy.total_drag_n,
                "trim_ratio": solved.trim_deg / legacy.trim_deg,
                "resistance_ratio": solved.total_drag_n / legacy.total_drag_n,
            }
        )

    max_lift_residual = max(r["lift_residual"] for r in rows)
    max_moment_residual = max(r["moment_residual_m"] for r in rows)
    claim(
        "F4-equilibrium lift residual",
        "< 1e-9 (relative)",
        f"{max_lift_residual:.3e}",
        max_lift_residual,
        max_lift_residual < 1e-9,
    )
    claim(
        "F4-equilibrium centre-of-pressure residual",
        "< 1e-8 m",
        f"{max_moment_residual:.3e} m",
        max_moment_residual,
        max_moment_residual < 1e-8,
    )

    lams = [r["wetted_length_ratio"] for r in rows]
    monotone = all(b < a for a, b in zip(lams, lams[1:]))
    claim(
        "equilibrium wetted length ratio falls with speed",
        "strictly decreasing over V = 8..20 m/s",
        " -> ".join(f"{lam:.3f}" for lam in lams),
        None,
        monotone,
    )

    lcg_sweep = []
    for lcg in (2.5, 2.8, 3.0, 3.2, 3.5, 4.0, 4.5):
        try:
            solution = equilibrium_trim(
                WEIGHT, 12.0, BEAM, DEADRISE, lcg, max_wetted_length_m=LENGTH
            )
            lcg_sweep.append(
                {
                    "lcg_m": lcg,
                    "lcg_over_beam": lcg / BEAM,
                    "trim_deg": solution.trim_deg,
                    "wetted_length_ratio": solution.wetted_length_ratio,
                }
            )
        except PlaningPredictionError as exc:
            lcg_sweep.append({"lcg_m": lcg, "error": str(exc)})

    trims = [row["trim_deg"] for row in lcg_sweep if "trim_deg" in row]
    claim(
        "moving the LCG aft raises the running trim",
        "trim strictly decreasing with lcg",
        " -> ".join(f"{t:.2f}" for t in trims),
        None,
        all(b < a for a, b in zip(trims, trims[1:])),
    )

    return {
        "speed_sweep": rows,
        "lcg_sweep": lcg_sweep,
        "max_lift_residual": max_lift_residual,
        "max_moment_residual_m": max_moment_residual,
    }, rows


# ---------------------------------------------------------------- figures
def make_figures(lift_models: dict, equilibrium: dict, rows: list[dict]) -> list[dict]:
    OUT.mkdir(parents=True, exist_ok=True)
    figures = []

    # Fig 1: the two lift correlations.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    taus = [1.0 + 0.1 * i for i in range(101)]
    for lam in (2.0, 3.0, 4.0):
        axes[0].plot(
            taus,
            [lift_coefficient(t, lam, DEADRISE) for t in taus],
            label=f"module CL, λ={lam:.0f}",
            linestyle="-",
        )
        axes[0].plot(
            taus,
            [classic_lift_coefficient(t, lam, 2.0, DEADRISE) for t in taus],
            label=f"classic CL (Cv=2), λ={lam:.0f}",
            linestyle="--",
        )
    axes[0].set_xlabel("trim τ [deg]")
    axes[0].set_ylabel("CL [-]")
    axes[0].set_title(f"Lift correlations, β={DEADRISE:.0f} deg")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.3)

    cvs = [0.8 + 0.05 * i for i in range(85)]
    for lam in (2.0, 3.0, 4.0):
        axes[1].plot(
            cvs,
            [
                abs(
                    lift_coefficient(4.0, lam, DEADRISE)
                    - classic_lift_coefficient(4.0, lam, cv, DEADRISE)
                )
                / classic_lift_coefficient(4.0, lam, cv, DEADRISE)
                for cv in cvs
            ],
            label=f"λ={lam:.0f}",
        )
    axes[1].set_xlabel("speed coefficient Cv [-]")
    axes[1].set_ylabel("|CL_module - CL_classic| / CL_classic")
    axes[1].set_title("Disagreement at τ=4 deg")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    path = OUT / "fig1_lift_models.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    figures.append(
        {
            "name": "fig1_lift_models.png",
            "caption": (
                "Left: the module's lift correlation (solid) against the classical "
                "Cv-dependent Savitsky statement (dashed). Right: their relative "
                "disagreement versus speed coefficient; no single Cv reconciles them."
            ),
        }
    )

    # Fig 2: solved versus assumed wetted length and trim.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    speeds = [row["speed_ms"] for row in rows]
    axes[0].plot(
        speeds,
        [row["wetted_length_ratio"] for row in rows],
        "o-",
        label="solved λ (two-equation closure)",
    )
    axes[0].axhline(
        DEFAULT_WETTED_LENGTH_RATIO,
        color="C1",
        linestyle="--",
        label=f"legacy λ = {DEFAULT_WETTED_LENGTH_RATIO} (assumed)",
    )
    axes[0].set_xlabel("speed V [m/s]")
    axes[0].set_ylabel("wetted length ratio λ [-]")
    axes[0].set_title("Wetted length ratio: solved vs assumed")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)

    axes[1].plot(
        speeds, [row["trim_deg"] for row in rows], "o-", label="solved τ (closure)"
    )
    axes[1].plot(
        speeds,
        [row["legacy_trim_deg"] for row in rows],
        "s--",
        label="legacy τ (λ = 3 assumed)",
    )
    axes[1].set_xlabel("speed V [m/s]")
    axes[1].set_ylabel("running trim τ [deg]")
    axes[1].set_title(f"Trim, Δ={WEIGHT / 1000:.0f} kN, B={BEAM} m, β={DEADRISE:.0f} deg")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    fig.tight_layout()
    path = OUT / "fig2_equilibrium_vs_assumed.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    figures.append(
        {
            "name": "fig2_equilibrium_vs_assumed.png",
            "caption": (
                "Left: λ solved from lift plus centre-of-pressure moment equilibrium "
                "(new path) against the constant λ=3 the legacy trim solve echoed back. "
                "Right: the resulting running trim at the same displacement and beam."
            ),
        }
    )

    # Fig 3: resistance consequence.
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    ax.plot(
        speeds, [row["resistance_n"] for row in rows], "o-", label="solved closure"
    )
    ax.plot(
        speeds,
        [row["legacy_resistance_n"] for row in rows],
        "s--",
        label="legacy (λ = 3, module CL)",
    )
    ax.set_xlabel("speed V [m/s]")
    ax.set_ylabel("total resistance R_T [N]")
    ax.set_title("Resistance: closure vs legacy path")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    path = OUT / "fig3_resistance_paths.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    figures.append(
        {
            "name": "fig3_resistance_paths.png",
            "caption": (
                "Total resistance from the two paths. The gap is the combined effect of "
                "solving λ and of the lift correlation, and it is far outside the ±10% "
                "band the module claims for the method."
            ),
        }
    )

    # Fig 4: the moment condition.
    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    lams = [0.5 + 0.02 * i for i in range(201)]
    for cv in (1.5, 2.0, 3.0, 4.0):
        ax.plot(
            lams,
            [center_of_pressure_ratio(lam, cv) * lam * BEAM for lam in lams],
            label=f"Cv={cv}",
        )
    ax.axhline(LCG, color="k", linestyle=":", label=f"LCG = {LCG} m")
    ax.set_xlabel("wetted length ratio λ [-]")
    ax.set_ylabel("centre of pressure from transom l_cp [m]")
    ax.set_title("Moment closure: l_cp(λ, Cv) = LCG")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    path = OUT / "fig4_moment_closure.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    figures.append(
        {
            "name": "fig4_moment_closure.png",
            "caption": (
                "Centre of pressure against wetted length ratio for several speed "
                "coefficients; each crossing with the LCG line is a moment-equilibrium "
                "solution, and its λ is what the legacy path could not produce."
            ),
        }
    )
    return figures


# ---------------------------------------------------------------- report
def write_outputs(payload: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "results.json", "w") as handle:
        json.dump(payload, handle, indent=2)

    lines = ["# Savitsky module audit -- measured evidence (m16954)", "", "## Claims", "",
             "| # | Claim | Expected | Actual | Error | Pass |", "|---:|---|---|---|---|---|"]
    for index, item in enumerate(payload["claims"], start=1):
        error = "—" if item["error"] is None else f"{item['error']:.2e}"
        lines.append(
            f"| {index} | {item['name']} | {item['expected']} | {item['actual']} | "
            f"{error} | {'PASS' if item['pass'] else 'FAIL'} |"
        )
    lines += ["", "## Solved-versus-assumed sweep", "",
              "| V [m/s] | Fr | Cv | τ solved [deg] | λ solved | R_T solved [N] | τ legacy [deg] | R_T legacy [N] |",
              "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in payload["equilibrium"]["speed_sweep"]:
        lines.append(
            f"| {row['speed_ms']:.0f} | {row['froude_number']:.2f} | {row['speed_coefficient']:.2f} | "
            f"{row['trim_deg']:.3f} | {row['wetted_length_ratio']:.3f} | {row['resistance_n']:.0f} | "
            f"{row['legacy_trim_deg']:.3f} | {row['legacy_resistance_n']:.0f} |"
        )
    lines += ["", "## Lift-model disagreement (β=20 deg)", "",
              "| τ [deg] | λ | module CL | classic CL (Cv=2) | rel. diff |",
              "|---:|---:|---:|---:|---:|"]
    for row in payload["lift_models"]["table_beta20"]:
        lines.append(
            f"| {row['trim_deg']:.0f} | {row['wetted_length_ratio']:.0f} | {row['module_cl']:.4f} | "
            f"{row['classic_cl_cv2']:.4f} | {100 * row['relative_difference']:+.1f}% |"
        )
    verdict = "ALL CLAIMS PASS" if all(c["pass"] for c in payload["claims"]) else "SOME CLAIMS FAIL"
    lines += ["", f"**Verdict**: {verdict}", ""]
    (OUT / "report.md").write_text("\n".join(lines))


def main() -> int:
    wetted_rows = audit_wetted_length()
    closed_form = audit_closed_form()
    clamp = audit_clamp()
    lift_models = audit_lift_models()
    equilibrium, sweep_rows = audit_equilibrium()
    figures = make_figures(lift_models, equilibrium, sweep_rows)

    payload = {
        "reference_craft": {
            "weight_n": WEIGHT,
            "beam_m": BEAM,
            "length_m": LENGTH,
            "deadrise_deg": DEADRISE,
            "lcg_m": LCG,
            "water_density_kg_m3": SEA_WATER_DENSITY,
        },
        "wetted_length_legacy": wetted_rows,
        "closed_form_trim": closed_form,
        "clamp_case": clamp,
        "lift_models": lift_models,
        "equilibrium": equilibrium,
        "claims": claims,
        "figures": figures,
        "provenance": {
            "mission": "m16954",
            "module": "tensorcad.savitsky",
            "note": (
                "The Cv-dependent coefficients of classic_lift_coefficient and "
                "center_of_pressure_ratio are not verified against the primary "
                "Savitsky (1964) paper in this environment (evidence level WITHHELD)."
            ),
        },
    }
    write_outputs(payload)

    print(f"wrote {OUT}/results.json and {len(figures)} figures")
    for item in claims:
        marker = "PASS" if item["pass"] else "FAIL"
        print(f"  [{marker}] {item['name']}: {item['actual']}")
    return 0 if all(item["pass"] for item in claims) else 1


if __name__ == "__main__":
    raise SystemExit(main())