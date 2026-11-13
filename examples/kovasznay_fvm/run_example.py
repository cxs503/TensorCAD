#!/usr/bin/env python3
"""Kovasznay 1948 exact solution — reference bundle and order-of-accuracy evidence.

The case-pool entry ``FVM-KOVASZNAY-ORDER`` asks for a measured log-log
convergence order ``>= 1.8`` (nominal second order) from 32/64/128 meshes.  The
solver lives in TensorFVM; TensorCAD supplies the geometry/reference side of the
benchmark, exactly as it does for the lid-driven cavity
(``examples/cavity_re5000``).  This script therefore produces, end to end and
without any hand-typed number:

1. verification of the exact solution itself (analytic Navier-Stokes residual,
   continuity, and two independent vorticity routes),
2. the 32/64/128 structured mesh family with quality metrics,
3. the *consistency* orders a second-order discretisation can deliver
   (midpoint-rule cell-average error, central-difference NS operator residual),
4. a self-test of the order-measurement harness
   (:func:`tensorcad.kovasznay.convergence_study`) against injected solvers of
   known order 1, 2 and 3 — so the harness is shown to *discriminate* order,
5. a solver-consumable JSON bundle (meshes, exact fields, probe table, spec).

Outputs (``--outdir``, default ``artifacts/kovasznay_fvm``)::

    checks.json                  every measured number + pass/fail
    case_spec.json               machine-readable benchmark specification
    mesh_{32,64,128}.json        nodes / quads / boundary tags / exact fields
    probe_table.json             exact solution on a probe lattice
    kovasznay_fields.png         u, v, p of the exact solution
    kovasznay_order.png          log-log error vs h with fitted orders
    kovasznay_profiles.png       centreline profiles at two Reynolds numbers

Usage::

    python examples/kovasznay_fvm/run_example.py [--outdir DIR] [--no-figures]

Exit code is 0 when every check passes, 1 otherwise.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.kovasznay import (
    KOVASZNAY_1948,
    KOVASZNAY_ORDER_ACCEPTANCE,
    KovasznayFlow,
    case_spec,
    convergence_study,
    generate_structured_mesh,
    midpoint_consistency_study,
    probe_table,
    truncation_error_study,
    write_case_bundle,
)

LEVELS = (32, 64, 128)
MIN_ORDER = KOVASZNAY_ORDER_ACCEPTANCE["min_acceptable_order"]


# ---------------------------------------------------------------------------
# 1. The exact solution verifies itself
# ---------------------------------------------------------------------------
def exact_solution_checks(flow: KovasznayFlow) -> dict:
    probe = probe_table(flow)
    max_ru = max(abs(r["ns_residual_u"]) for r in probe)
    max_rv = max(abs(r["ns_residual_v"]) for r in probe)
    max_div = max(abs(r["continuity"]) for r in probe)
    max_vort = max(r["vorticity_mismatch"] for r in probe)
    lam = flow.lam
    quadratic = lam * lam - flow.reynolds_number * lam - 4.0 * math.pi ** 2
    return {
        "reynolds_number": flow.reynolds_number,
        "lam": lam,
        "lam_quadratic_residual": quadratic,
        "n_probe_points": len(probe),
        "max_abs_momentum_residual_u": max_ru,
        "max_abs_momentum_residual_v": max_rv,
        "max_abs_continuity_residual": max_div,
        "max_abs_vorticity_route_mismatch": max_vort,
        "passes": max(max_ru, max_rv, max_div) < 1e-10 and max_vort < 1e-12,
    }


# ---------------------------------------------------------------------------
# 4. Harness self-test: inject solvers of known order
# ---------------------------------------------------------------------------
def injected_solver(flow: KovasznayFlow, order: int, amplitude: float = 0.37):
    """Cell-centre velocity = exact + a smooth defect of size ``amplitude*h**order``."""

    def solver(mesh):
        h = mesh.hx
        out = []
        for x, y in mesh.centroids():
            u, v = flow.velocity(x, y)
            du = amplitude * h ** order * math.sin(3.0 * x) * math.cos(2.0 * y)
            dv = 0.57 * amplitude * h ** order * math.cos(x) * math.sin(4.0 * y)
            out.append((u - du, v + dv))
        return out

    return solver


def harness_self_test(flow: KovasznayFlow, levels=LEVELS) -> dict:
    rows = []
    for order in (1, 2, 3):
        study = convergence_study(flow, levels, injected_solver(flow, order))
        rows.append({
            "injected_order": order,
            "measured_order": study["order"],
            "r2": study["r2"],
            "absolute_error": abs(study["order"] - order),
        })
    exact = convergence_study(
        flow, levels, lambda mesh: [flow.velocity(x, y) for x, y in mesh.centroids()])
    return {
        "injected": rows,
        "max_abs_order_error": max(r["absolute_error"] for r in rows),
        "exact_solver_order": exact["order"],
        "exact_solver_note": exact.get("note"),
        # The harness must recover the injected order to well inside 0.1, and must
        # report the first-order solver as failing the 1.8 acceptance.
        "passes": (max(r["absolute_error"] for r in rows) < 0.1
                   and exact["order"] is None),
    }


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def make_figures(flow: KovasznayFlow, outdir: Path,
                 midpoint: dict, truncation: dict,
                 self_test: dict, levels=LEVELS) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception:  # pragma: no cover - optional dependency
        print("  [SKIP] matplotlib/numpy unavailable; no figures written")
        return []

    written: list[str] = []

    # -- field contours -----------------------------------------------------
    n = 129
    xs = np.linspace(flow.x_min, flow.x_max, n)
    ys = np.linspace(flow.y_min, flow.y_max, n)
    X, Y = np.meshgrid(xs, ys)
    U = np.empty_like(X)
    V = np.empty_like(X)
    P = np.empty_like(X)
    for j in range(n):
        for i in range(n):
            u, v = flow.velocity(float(X[j, i]), float(Y[j, i]))
            U[j, i] = u
            V[j, i] = v
            P[j, i] = flow.pressure(float(X[j, i]), float(Y[j, i]))
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    for ax, field, label, cmap in (
        (axes[0], U, "u", "RdBu_r"),
        (axes[1], V, "v", "RdBu_r"),
        (axes[2], P, "p", "viridis"),
    ):
        cs = ax.contourf(X, Y, field, levels=24, cmap=cmap)
        ax.set_title(f"Kovasznay {label} (Re={flow.reynolds_number:g})")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_aspect("equal")
        fig.colorbar(cs, ax=ax, shrink=0.85)
    fig.tight_layout()
    path = outdir / "kovasznay_fields.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    written.append(path.name)

    # -- order plot ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    series = [
        ("midpoint rule (cell avg vs centroid)", midpoint["h"],
         [r["l2_velocity"] for r in midpoint["levels"]], midpoint["order_velocity"]),
        ("NS operator truncation (u-momentum)", truncation["h"],
         [r["l2_momentum_x"] for r in truncation["levels"]],
         truncation["order_momentum_x"]),
    ]
    for label, hs, errs, order in series:
        ax.loglog(hs, errs, "o-", label=f"{label}: p={order:.3f}")
    # injected-solver lines from the harness self-test
    for row in self_test["injected"]:
        if row["injected_order"] == 2:
            study = convergence_study(flow, levels, injected_solver(flow, 2))
            ax.loglog(study["h"], study["l2_error"], "s--",
                      label=f"harness self-test (injected p=2): p={row['measured_order']:.3f}")
            break
    ref_h = np.array(midpoint["h"], dtype=float)
    ref_e = np.array([r["l2_velocity"] for r in midpoint["levels"]], dtype=float)
    ax.loglog(ref_h, ref_e[0] * (ref_h / ref_h[0]) ** 2, "k:", lw=1.0,
              label="slope 2 reference")
    ax.set_xlabel("h")
    ax.set_ylabel("L2 error")
    ax.set_title("Kovasznay order evidence vs h (32/64/128)")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7.5, loc="lower right")
    fig.tight_layout()
    path = outdir / "kovasznay_order.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    written.append(path.name)

    # -- centreline profiles at two Reynolds numbers ------------------------
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
    for re_ in (40.0, 100.0):
        fl = KovasznayFlow(reynolds_number=re_)
        yy = [flow.y_min + (flow.y_max - flow.y_min) * k / 100.0 for k in range(101)]
        uu = [fl.velocity(-0.25, y)[0] for y in yy]
        axes[0].plot(uu, yy, label=f"Re={re_:g}")
        xx = [flow.x_min + (flow.x_max - flow.x_min) * k / 100.0 for k in range(101)]
        axes[1].plot(xx, [fl.velocity(x, 0.25)[1] for x in xx], label=f"Re={re_:g}")
    axes[0].set_xlabel("u at x=-0.25")
    axes[0].set_ylabel("y")
    axes[0].set_title("Streamwise velocity profile")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("v at y=0.25")
    axes[1].set_title("Transverse velocity profile")
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    path = outdir / "kovasznay_profiles.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    written.append(path.name)
    return written


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", default=str(REPO_ROOT / "artifacts" / "kovasznay_fvm"))
    parser.add_argument("--re", type=float, default=40.0, dest="reynolds")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    flow = KovasznayFlow(reynolds_number=args.reynolds)

    print("=" * 78)
    print("Kovasznay 1948 exact solution — reference bundle / order evidence")
    print("=" * 78)
    print(f"reference: {KOVASZNAY_1948['citation']}")
    print(f"           doi:{KOVASZNAY_1948['doi']}")
    print(f"Re = {flow.reynolds_number:g}   lam = {flow.lam:.12f}   nu = {flow.nu:.6e}")
    print(f"domain = [{flow.x_min}, {flow.x_max}] x [{flow.y_min}, {flow.y_max}]")
    print()

    # 1. exact solution -----------------------------------------------------
    exact = exact_solution_checks(flow)
    print("[1] exact-solution verification (analytic derivatives)")
    print(f"    max |R_u| = {exact['max_abs_momentum_residual_u']:.3e}")
    print(f"    max |R_v| = {exact['max_abs_momentum_residual_v']:.3e}")
    print(f"    max |div| = {exact['max_abs_continuity_residual']:.3e}")
    print(f"    max vorticity route mismatch = {exact['max_abs_vorticity_route_mismatch']:.3e}")
    print(f"    lam quadratic residual = {exact['lam_quadratic_residual']:.3e}")
    print()

    # 2. mesh family --------------------------------------------------------
    print("[2] mesh family")
    mesh_rows = []
    for n in LEVELS:
        summary = generate_structured_mesh(flow, n).summary()
        mesh_rows.append(summary)
        print(f"    n={summary['nx']:>3}  h={summary['hx']:.6f}  cells={summary['n_cells']:>5} "
              f"nodes={summary['n_nodes']:>5}  area_err={summary['area_error']:.2e}  "
              f"neg_jac={summary['negative_jacobians']}")
    mesh_ok = all(r["negative_jacobians"] == 0 and abs(r["area_error"]) < 1e-12
                  for r in mesh_rows)
    print()

    # 3. consistency orders -------------------------------------------------
    midpoint = midpoint_consistency_study(flow, LEVELS)
    truncation = truncation_error_study(flow, LEVELS)
    print("[3] consistency orders (what a 2nd-order discretisation can deliver)")
    print("    grid        cell-avg vs centroid L2      central-diff NS residual L2")
    for row_m, row_t in zip(midpoint["levels"], truncation["levels"], strict=False):
        print(f"    n={row_m['n']:>3}   u {row_m['l2_u']:.6e}  v {row_m['l2_v']:.6e}"
              f"      Rx {row_t['l2_momentum_x']:.6e}  div {row_t['l2_continuity']:.6e}")
    print(f"    fitted order: midpoint(u,v) = {midpoint['order_u']:.4f}, "
          f"{midpoint['order_v']:.4f} (R2={midpoint['r2_velocity']:.6f})")
    print(f"    fitted order: truncation(Rx,Ry,div) = {truncation['order_momentum_x']:.4f}, "
          f"{truncation['order_momentum_y']:.4f}, {truncation['order_continuity']:.4f}")
    print()

    # 4. harness self-test --------------------------------------------------
    self_test = harness_self_test(flow, LEVELS)
    print("[4] order-measurement harness self-test (injected known-order solvers)")
    for row in self_test["injected"]:
        print(f"    injected p={row['injected_order']}  ->  measured p="
              f"{row['measured_order']:.4f}  |Δ|={row['absolute_error']:.2e}")
    print(f"    exact solver -> order={self_test['exact_solver_order']} "
          f"({self_test['exact_solver_note']})")
    print()

    # 5. artefacts ----------------------------------------------------------
    manifest = write_case_bundle(flow, outdir, LEVELS)
    print("[5] case bundle written")
    for name, digest in manifest["files"].items():
        print(f"    {name:<18} sha256={digest[:16]}…")
    print()

    checks = {
        "case_code": "FVM-KOVASZNAY-ORDER",
        "reference": KOVASZNAY_1948,
        "flow": flow.as_dict(),
        "acceptance": KOVASZNAY_ORDER_ACCEPTANCE,
        "exact_solution": exact,
        "meshes": mesh_rows,
        "midpoint_consistency": midpoint,
        "truncation_error": truncation,
        "harness_self_test": self_test,
        "bundle": manifest,
        "conclusions": {
            "exact_solution_real": exact["passes"],
            "mesh_family_valid": bool(mesh_ok),
            "midpoint_order_meets_1p8": bool(midpoint["order_velocity"] >= MIN_ORDER),
            "truncation_order_meets_1p8": bool(
                min(truncation["order_momentum_x"], truncation["order_momentum_y"],
                    truncation["order_continuity"]) >= MIN_ORDER),
            "harness_discriminates_order": self_test["passes"],
            "solver_error_order": "WITHHELD (requires the TensorFVM solve; not run here)",
        },
    }
    if not args.no_figures:
        figures = make_figures(flow, outdir, midpoint, truncation, self_test, LEVELS)
        checks["figures"] = figures
        if figures:
            print(f"[6] figures: {', '.join(figures)}")
            print()

    (outdir / "checks.json").write_text(
        json.dumps(checks, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (outdir / "case_spec.json").write_text(
        json.dumps(case_spec(flow, LEVELS), indent=2, sort_keys=True) + "\n",
        encoding="utf-8")

    verdict = all([
        checks["conclusions"]["exact_solution_real"],
        checks["conclusions"]["mesh_family_valid"],
        checks["conclusions"]["midpoint_order_meets_1p8"],
        checks["conclusions"]["truncation_order_meets_1p8"],
        checks["conclusions"]["harness_discriminates_order"],
    ])
    print("-" * 78)
    for key, value in checks["conclusions"].items():
        print(f"    {key:<32} {value}")
    print("-" * 78)
    print(f"RESULT: {'PASS' if verdict else 'FAIL'}  (checks.json -> {outdir / 'checks.json'})")
    return 0 if verdict else 1


if __name__ == "__main__":
    raise SystemExit(main())
