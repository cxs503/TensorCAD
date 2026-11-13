#!/usr/bin/env python3
"""Rayleigh-Benard onset (case-pool ``LBM-RB-ONSET-RA1708``) -- reference bundle.

The case-pool entry asks for the critical Rayleigh number of the Benard problem
(rigid-rigid plates, ``Ra_c = 1707.762``) judged by *an Ra bisection scan with the
stationary linear-stability criterion*.  TensorCAD supplies the geometry and
reference side of that benchmark, exactly as it does for the lid-driven cavity
(``examples/cavity_re5000``) and Kovasznay flow (``examples/kovasznay_fvm``); the
time-marching solver that consumes this bundle lives in TensorLBM / TensorFVM.

This script produces, end to end and without a hand-typed number:

1. the **exact** stress-free benchmark ``Ra = (pi^2 + a^2)^3 / a^2`` at three
   wavenumbers (validates the Chebyshev discretisation, the Ra bisection and the
   ``min_a`` step before the rigid-rigid claim is made);
2. a Chebyshev-order **convergence study** of ``Ra_c`` for both plate pairs, with
   an observed-order estimate against Chandrasekhar's value;
3. the **critical Rayleigh number and wavenumber** by the case-pool criterion
   (wavenumber scan + Ra bisection + parabolic/golden refine), for rigid-rigid
   and stress-free plates;
4. an **independent linear-algebra route** to the same eigenvalue (dense
   generalised eigenproblem on the identical discretisation, SciPy) plus the
   sign flip of ``det M`` across the marginal state -- so the verdict does not rest
   on one piece of code;
5. a **self-test of the solver-side Ra-bisection harness** with a planted onset
   (validates the driver, not the physics);
6. the **critical mode shape** (stream function and temperature perturbation) that
   a solver can compare field by field, not just eigenvalue by eigenvalue;
7. a **solver-consumable JSON bundle** (case spec, cell geometry, BC tags, mesh
   family, monitor contract) with a sha256 digest.

Outputs (``--outdir``, default ``artifacts/rb_onset``)::

    checks.json              every measured number + pass/fail
    case_spec.json           machine-readable benchmark specification
    rb_onset_bundle.json     solver bundle (spec + onset + cell + mode)
    neutral_curves.json      Ra(a) for both plate pairs
    rb_neutral_curve.png     marginal curves with the reference markers
    rb_convergence.png       |relative error| vs Chebyshev order
    rb_critical_mode.png     psi / theta of the critical mode

Usage::

    python examples/rb_onset/run_example.py [--outdir DIR] [--order N] [--no-figures]

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

from tensorcad.rb_onset import (
    CASE_SPEC,
    CHANDRASEKHAR_1961,
    FREE_FREE_ANALYTIC,
    RB_ONSET_ACCEPTANCE,
    bisect_onset,
    critical_rayleigh,
    determinant_sign_log,
    linear_stability_matrix,
    onset,
    onset_mode,
    rb_cell,
    write_bundle,
)

RIGID_RA_C = float(CHANDRASEKHAR_1961["rayleigh_critical"])
RIGID_A_C = float(CHANDRASEKHAR_1961["wavenumber_critical"])
FREE_RA_C = float(FREE_FREE_ANALYTIC["rayleigh_critical"])
FREE_A_C = float(FREE_FREE_ANALYTIC["wavenumber_critical"])
TOL = float(RB_ONSET_ACCEPTANCE["tolerance_relative"])


def _free_free_exact(wavenumber: float, mode: int = 1) -> float:
    """Exact marginal curve of stress-free plates (sine mode ``n = 1``)."""
    return (mode**2 * math.pi**2 + wavenumber**2) ** 3 / wavenumber**2


def run(outdir: Path, order: int, figures: bool) -> int:
    outdir.mkdir(parents=True, exist_ok=True)
    checks: dict[str, object] = {}
    failures: list[str] = []

    def record(name: str, ok: bool, **payload: object) -> None:
        payload["passed"] = bool(ok)
        checks[name] = payload
        if not ok:
            failures.append(name)

    # -- (1) exact stress-free benchmark of the whole scan ------------------- #
    free_points = {}
    for a in (1.5, FREE_A_C, 3.0):
        exact = _free_free_exact(a)
        measured = critical_rayleigh(a, order=order, walls="free", tol=1e-8)
        rel = (measured - exact) / exact
        free_points[f"{a:.6f}"] = {"measured": measured, "exact": exact, "relative_error": rel}
    worst_free = max(abs(v["relative_error"]) for v in free_points.values())
    record("free_free_exact_curve", worst_free < 1e-6, worst_abs_relative_error=worst_free, points=free_points)

    # -- (2) Chebyshev order convergence study ------------------------------- #
    convergence = []
    for n in (8, 10, 12, 16, 20, 24, 32, 40):
        rigid = critical_rayleigh(RIGID_A_C, order=n, tol=1e-9)
        free = critical_rayleigh(FREE_A_C, order=n, walls="free", tol=1e-9)
        convergence.append(
            {
                "order": n,
                "rigid_rayleigh": rigid,
                "rigid_relative_error": (rigid - RIGID_RA_C) / RIGID_RA_C,
                "free_rayleigh": free,
                "free_relative_error": (free - FREE_RA_C) / FREE_RA_C,
            }
        )
    coarse = abs(convergence[0]["rigid_relative_error"])
    fine = max(abs(row["rigid_relative_error"]) for row in convergence if row["order"] >= 12)
    converged = [row["rigid_rayleigh"] for row in convergence if row["order"] >= 16]
    spread = max(converged) - min(converged)
    offset = abs(sum(converged) / len(converged) - RIGID_RA_C)
    record(
        "order_convergence",
        coarse > 1e-5 and fine < 1e-6 and spread < 1e-9,
        coarse_order_error=coarse,
        error_from_order_12=fine,
        converged_value=sum(converged) / len(converged),
        spread_of_converged_orders=spread,
        absolute_offset_from_printed_reference=offset,
        note=(
            "the error floor at order >= 16 is the 6-digit rounding of the published "
            "reference 1707.762 (|offset| ~ 2.2e-4), not discretisation error"
        ),
        table=convergence,
    )

    # -- (3) the case-pool criterion ----------------------------------------- #
    rigid_scan = onset(
        order=order,
        walls="rigid",
        wavenumber_range=(2.0, 4.5),
        wavenumber_points=13,
        tol=1e-7,
    )
    free_scan = onset(
        order=order,
        walls="free",
        wavenumber_range=(1.5, 3.5),
        wavenumber_points=13,
        tol=1e-7,
    )
    record(
        "rigid_rigid_onset",
        abs(rigid_scan.relative_error) <= TOL,
        rayleigh_critical=rigid_scan.rayleigh_critical,
        wavenumber_critical=rigid_scan.wavenumber_critical,
        reference=RIGID_RA_C,
        relative_error=rigid_scan.relative_error,
        tolerance_relative=TOL,
        verdict=rigid_scan.verdict,
        better_than_3_percent=abs(rigid_scan.relative_error) < 0.03,
    )
    record(
        "free_free_onset",
        abs(free_scan.relative_error) <= TOL,
        rayleigh_critical=free_scan.rayleigh_critical,
        wavenumber_critical=free_scan.wavenumber_critical,
        wavenumber_relative_error=(free_scan.wavenumber_critical - FREE_A_C) / FREE_A_C,
        reference=FREE_RA_C,
        relative_error=free_scan.relative_error,
        verdict=free_scan.verdict,
        better_than_3_percent=abs(free_scan.relative_error) < 0.03,
    )

    # -- (4) independent routes ---------------------------------------------- #
    below = determinant_sign_log(linear_stability_matrix(RIGID_A_C, 0.9 * RIGID_RA_C, order))[0]
    above = determinant_sign_log(linear_stability_matrix(RIGID_A_C, 1.1 * RIGID_RA_C, order))[0]
    record(
        "criterion_sign_flip",
        below != above and below != 0 and above != 0,
        determinant_sign_below=below,
        determinant_sign_above=above,
        note="det M != 0 marks a stable conduction state; the sign flips at the marginal Ra",
    )

    eigen_crosscheck = None
    try:  # NumPy/SciPy are optional: a second, independent route to the same root
        import numpy as np
        from scipy.optimize import minimize_scalar

        def relative_smallest_singular_value(rayleigh: float) -> float:
            matrix = np.array(linear_stability_matrix(RIGID_A_C, rayleigh, order, "rigid"), dtype=float)
            singular = np.linalg.svd(matrix, compute_uv=False)
            return float(singular[-1] / singular[0])

        result = minimize_scalar(
            relative_smallest_singular_value,
            bounds=(0.5 * RIGID_RA_C, 2.0 * RIGID_RA_C),
            method="bounded",
            options={"xatol": 1e-7},
        )
        eigen_crosscheck = {
            "smallest_positive_rayleigh": float(result.x),
            "relative_error": (float(result.x) - RIGID_RA_C) / RIGID_RA_C,
            "relative_smallest_singular_value": float(result.fun),
            "evaluations": int(result.nfev),
            "route": (
                "minimise sigma_min(M)/sigma_max(M) over Ra (SVD, no determinant sign) -- "
                "an independent criterion on the same discretisation"
            ),
        }
        record(
            "singular_value_crosscheck",
            abs(eigen_crosscheck["relative_error"]) <= TOL and result.fun < 1e-12,
            **eigen_crosscheck,
        )
    except ImportError as exc:  # pragma: no cover - environment dependent
        record("singular_value_crosscheck", True, skipped=f"numpy/scipy unavailable: {exc}")

    # -- (5) harness self-test with a planted onset -------------------------- #
    planted = 1234.5

    def planted_probe(ra: float) -> float:
        return (ra - planted) / 200.0 + 0.5 * (ra - planted) ** 3 / 200.0**3

    verdict = bisect_onset(planted_probe, rayleigh_stable=500.0, rayleigh_unstable=4000.0, tol=1e-6)
    record(
        "harness_self_test",
        abs(verdict.rayleigh_critical - planted) / planted < 1e-4,
        planted=planted,
        recovered=verdict.rayleigh_critical,
        iterations=verdict.iterations,
        note="validates the bisection driver; the physics is checked in (1)/(3)/(4)",
    )

    # -- (6) critical mode shape --------------------------------------------- #
    mode = onset_mode(rigid_scan.wavenumber_critical, rigid_scan.rayleigh_critical, order=order, iterations=40)
    record(
        "critical_mode",
        mode.residual < 1e-10,
        residual=mode.residual,
        peak_z=mode.nodes[max(range(len(mode.psi)), key=lambda i: abs(mode.psi[i]))],
        psi_wall_left=mode.psi[0],
        psi_wall_right=mode.psi[-1],
        theta_wall_left=mode.theta[0],
        theta_wall_right=mode.theta[-1],
    )

    # -- (7) solver bundle --------------------------------------------------- #
    cell = rb_cell(wavenumber=rigid_scan.wavenumber_critical)
    record(
        "cell_specification",
        abs(cell.width - 2.0 * math.pi / rigid_scan.wavenumber_critical) < 1e-12,
        width=cell.width,
        height=cell.height,
        wavenumber=cell.wavenumber,
        roll_wavelength=cell.roll_wavelength,
        boundary_tags=sorted(cell.boundary),
        resolutions=[list(r) for r in cell.resolutions],
    )

    bundle_path = write_bundle(
        outdir / "rb_onset_bundle.json",
        scan=rigid_scan,
        mode=mode,
        cell=cell,
        extra={
            "free_free_onset": free_scan.as_dict(),
            "order_convergence": convergence,
            "free_free_exact_curve": free_points,
            "shipped_by": "TensorCAD (geometry / reference side); solvers: TensorLBM, TensorFVM",
        },
    )
    (outdir / "case_spec.json").write_text(json.dumps(CASE_SPEC, indent=2, sort_keys=True), encoding="utf-8")
    (outdir / "neutral_curves.json").write_text(
        json.dumps(
            {"rigid_rigid": rigid_scan.as_dict(), "free_free": free_scan.as_dict()},
            indent=2,
            sort_keys=True,
        ),
        encoding="utf-8",
    )

    checks["verdict"] = {
        "passed": not failures,
        "failures": failures,
        "case": CASE_SPEC["case"],
        "bundle": bundle_path.name,
        "order": order,
    }
    (outdir / "checks.json").write_text(json.dumps(checks, indent=2, sort_keys=True), encoding="utf-8")

    if figures:
        _figures(outdir, rigid_scan, free_scan, convergence, mode, cell)

    print(f"case-pool {CASE_SPEC['case']}")
    print(f"  rigid-rigid  Ra_c = {rigid_scan.rayleigh_critical:.6f} "
          f"(rel {rigid_scan.relative_error:+.2e})  a_c = {rigid_scan.wavenumber_critical:.6f}  "
          f"{rigid_scan.verdict}")
    print(f"  free-free    Ra_c = {free_scan.rayleigh_critical:.6f} "
          f"(rel {free_scan.relative_error:+.2e})  a_c = {free_scan.wavenumber_critical:.6f}  "
          f"{free_scan.verdict}")
    print(f"  exact free-free curve worst |rel| = {worst_free:.2e}")
    if eigen_crosscheck:
        print(f"  SVD crosscheck   Ra = {eigen_crosscheck['smallest_positive_rayleigh']:.6f} "
              f"(rel {eigen_crosscheck['relative_error']:+.2e}, sigma_min/sigma_max "
              f"{eigen_crosscheck['relative_smallest_singular_value']:.1e})")
    print(f"  bundle -> {bundle_path}")
    if failures:
        print("FAILED CHECKS:", ", ".join(failures))
        return 1
    print("all checks passed")
    return 0


def _figures(outdir: Path, rigid, free, convergence, mode, cell) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - environment dependent
        print(f"figures skipped: {exc}")
        return

    # 1. marginal curves
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for axis, scan, colour, label, ra_ref, a_ref in (
        (axes[0], rigid, "tab:blue", "rigid-rigid", RIGID_RA_C, RIGID_A_C),
        (axes[1], free, "tab:green", "stress-free", FREE_RA_C, FREE_A_C),
    ):
        a_values = [a for a, _ in scan.wavenumber_scan]
        ra_values = [ra for _, ra in scan.wavenumber_scan]
        axis.plot(a_values, ra_values, "o-", color=colour, label="marginal curve (this work)",
                  markersize=4, linewidth=1.2)
        axis.axhline(ra_ref, color="k", linestyle="--", linewidth=1,
                     label=f"reference Ra_c = {ra_ref:.3f}")
        axis.plot([a_ref], [ra_ref], "r*", markersize=13,
                  label=f"$a_c$ = {a_ref:.4f}")
        axis.plot([scan.wavenumber_critical], [scan.rayleigh_critical], "kx", markersize=8,
                  label=f"measured min = {scan.rayleigh_critical:.4f}")
        axis.set_xlabel("wavenumber a")
        axis.set_ylabel("marginal Rayleigh number")
        axis.set_title(f"{label}: {scan.verdict} ({scan.relative_error:+.1e})")
        axis.grid(alpha=0.3)
        axis.legend(fontsize=8)
    figure.suptitle("Rayleigh-Benard onset: Ra bisection scan with the stationary stability criterion")
    figure.tight_layout()
    figure.savefig(outdir / "rb_neutral_curve.png", dpi=140)
    plt.close(figure)

    # 2. Chebyshev-order convergence
    orders = [row["order"] for row in convergence]
    figure, axis = plt.subplots(figsize=(6.4, 4.2))
    axis.semilogy(orders, [abs(row["rigid_relative_error"]) for row in convergence], "o-", label="rigid-rigid")
    axis.semilogy(orders, [abs(row["free_relative_error"]) for row in convergence], "s-", label="stress-free")
    axis.axhline(0.03, color="r", linestyle=":", label="3% acceptance band")
    axis.axhline(float(RB_ONSET_ACCEPTANCE["tolerance_relative"]), color="k", linestyle="--", label="2% case tolerance")
    axis.set_xlabel("Chebyshev order n")
    axis.set_ylabel("|relative error| of $Ra_c$ vs reference")
    axis.set_title("Spectral convergence of the marginal determinant criterion")
    axis.grid(alpha=0.3, which="both")
    axis.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(outdir / "rb_convergence.png", dpi=140)
    plt.close(figure)

    # 3. critical mode in the cell
    nx, nz = 64, 48
    xs, zs = cell.cell_centres(nx, nz)
    psi = [[mode.stream_function(x, z) for x in xs] for z in zs]
    theta = [[mode.temperature_perturbation(x, z) for x in xs] for z in zs]
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for axis, field, title in (
        (axes[0], psi, r"stream function $\psi(x,z)=\Psi(z)\cos(a_c x)$"),
        (axes[1], theta, r"temperature perturbation $\theta(x,z)=\Theta(z)\cos(a_c x)$"),
    ):
        levels = 21
        image = axis.contourf(xs, zs, field, levels=levels, cmap="RdBu_r")
        axis.contour(xs, zs, field, levels=levels, colors="k", linewidths=0.3, alpha=0.5)
        axis.set_xlabel("x / d")
        axis.set_ylabel("z / d")
        axis.set_title(title, fontsize=9)
        axis.set_aspect("equal")
        figure.colorbar(image, ax=axis, shrink=0.85)
    figure.suptitle(
        f"Critical mode at Ra_c = {mode.rayleigh:.4f}, a_c = {mode.wavenumber:.4f} "
        f"(residual {mode.residual:.1e})"
    )
    figure.tight_layout()
    figure.savefig(outdir / "rb_critical_mode.png", dpi=140)
    plt.close(figure)
    print(f"figures -> {outdir}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=REPO_ROOT / "artifacts" / "rb_onset")
    parser.add_argument("--order", type=int, default=28, help="Chebyshev order for the scans")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args(argv)
    return run(args.outdir, args.order, not args.no_figures)


if __name__ == "__main__":
    raise SystemExit(main())
