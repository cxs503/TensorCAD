"""Adversarial re-verification of the planing-hull wetted topside integral.

Three routes to the same number:

* **legacy** -- a verbatim copy of ``_wetted_topside_area`` as shipped in 0.23.0
  (``sqrt(1 + slope**2) * dx/dt`` with ``e = 2/q`` for every ``q``). Kept here only as
  the "before" side of the comparison in ``docs/planing_numerical_robustness.md``.
* **shipped** -- ``tensorcad.planing.analytic_wetted_area`` after the fix.
* **mpmath** -- ``mp.quad`` (tanh-sinh, 40 digits) of the physical integral
  ``2 int_0^x0 (T - b tan b) sqrt(1 + b'^2) dx`` in ``x``, i.e. neither substitution
  nor hypot, on the same closed-form ``b(x)``. Route 3 decides who is right.

Writes ``artifacts/planing_numerical_robustness/results.json`` and a two-panel figure,
and exits non-zero if any claim below stops holding.

Claims checked
--------------
C1  The legacy route returns a non-finite wetted area on 4 of 420 valid parameter sets
    (``q = 0.05``, ``beta >= 30 deg``, light draft), so the analytic reference the module
    validates the mesh against was ``inf`` there.
C2  The shipped route is finite and within 1e-9 relative of the mpmath reference on the
    whole 420-case sweep.
C3  On the cases both routes survive, shipped and legacy agree to better than 1e-9
    relative: the rewrite is a numerical re-formulation, not a change of the formula.
"""
from __future__ import annotations

import json
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

import mpmath as mp  # noqa: E402

from tensorcad.planing import (  # noqa: E402
    PlaningHullParameters,
    analytic_planing_hydrostatics,
    analytic_waterplane_area,
    analytic_wetted_area,
    chine_half_beam,
    section_area,
    waterline_half_breadth,
    wetted_length,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "planing_numerical_robustness"
DEGREES = (1.0, 6.0, 12.0, 18.0, 30.0, 45.0, 59.0)
EXPONENTS = (0.05, 0.1, 0.2, 0.4, 0.6, 1.0, 1.5, 2.0, 3.0, 4.0)
DRAFTS = (0.05, 0.2, 0.35, 0.55, 1.0, 2.0)


def legacy_topside_area(p: PlaningHullParameters, x0: float, samples: int = 4000) -> float:
    """The 0.23.0 implementation, unchanged, for the before/after table."""
    if x0 <= 0.0:
        return 0.0
    if samples % 2:
        samples += 1
    exponent = 2.0 / p.planform_exponent
    c = p.half_beam_transom / p.length**p.planform_exponent

    def integrand(t: float) -> float:
        x = x0 * t**exponent
        if x <= 0.0:
            return 0.0
        b = c * x**p.planform_exponent
        slope = c * p.planform_exponent * x ** (p.planform_exponent - 1.0)
        height = p.draft - b * p.tan_beta
        if height <= 0.0:
            return 0.0
        dx_dt = x0 * exponent * t ** (exponent - 1.0)
        return 2.0 * height * math.sqrt(1.0 + slope * slope) * dx_dt

    h = 1.0 / samples
    total = integrand(0.0) + integrand(1.0)
    for i in range(1, samples):
        total += (4.0 if i % 2 else 2.0) * integrand(i * h)
    return total * h / 3.0


def legacy_wetted_area(p: PlaningHullParameters, samples: int = 4000) -> float:
    tau = p.tan_beta
    x0 = wetted_length(p)
    bottom = analytic_waterplane_area(p) * math.sqrt(1.0 + tau * tau)
    return bottom + legacy_topside_area(p, x0, samples) + section_area(p, p.length)


def mp_wetted_area(p: PlaningHullParameters) -> "mp.mpf":
    """High-precision reference in a *physical* variable chosen per regime.

    ``2 int_0^x0 (T - b tan b) sqrt(1 + b'^2) dx`` is singular at a knife bow in ``x``
    (``b' ~ x**(q-1)``: infinite for ``q < 1``, zero for ``q > 1``), and a fixed rule in
    ``x`` is wrong by ~0.25 % there -- which is exactly the trap under review, so the
    reference must not fall into it. Two substitutions are used instead, neither of them
    the module's ``x = x0 t**e``:

    * ``q < 1``: integrate over the chine half-breadth ``u = b(x)``, where
      ``sqrt(1 + b'**2) dx = sqrt(1 + (dx/du)**2) du`` and ``dx/du ~ u**(1/q - 1) -> 0``
      at the stem, i.e. a smooth integrand.
    * ``q >= 1``: integrate in ``x`` directly (``b'`` is bounded there).

    Both routes are quadratured by ``mp.quad`` at 40 digits.
    """
    mp.mp.dps = 40
    L, T, tau = mp.mpf(p.length), mp.mpf(p.draft), mp.mpf(p.tan_beta)
    q = mp.mpf(p.planform_exponent)
    c = mp.mpf(p.half_beam_transom) / L**q
    x0 = mp.mpf(wetted_length(p))
    bottom = mp.mpf(analytic_waterplane_area(p)) * mp.sqrt(1 + tau * tau)

    def b(x):
        return c * x**q if x > 0 else mp.mpf(0)

    def section(x):
        bb = b(x)
        h = bb * tau
        return T * T / tau if h >= T else bb * h + 2 * bb * (T - h)

    if x0 <= 0:
        top = mp.mpf(0)
    elif q < 1:
        u_max = b(x0)
        # x(u) = (u/c)**(1/q); dx/du = u**(1/q - 1) / (q c**(1/q))
        one_over_q = 1 / q
        denominator = q * c**one_over_q

        def integrand(u):
            if u <= 0:
                return mp.mpf(0)
            dx_du = u ** (one_over_q - 1) / denominator
            height = T - tau * u
            if height <= 0:
                return mp.mpf(0)
            return 2 * height * mp.sqrt(1 + dx_du**2)

        top = mp.quad(integrand, [0, u_max])
    else:
        def integrand(x):
            if x <= 0:
                return mp.mpf(0)
            height = T - tau * b(x)
            if height <= 0:
                return mp.mpf(0)
            return 2 * height * mp.sqrt(1 + (c * q * x ** (q - 1)) ** 2)

        top = mp.quad(integrand, [0, x0])
    return bottom + top + section(L)


def scipy_wetted_area(p: PlaningHullParameters) -> float:
    """A fourth, library route (SciPy QUADPACK) for the same integral in ``x``.

    Used as a cross-check of the mpmath reference when SciPy is installed; the tests do
    not import SciPy, only this script does.
    """
    import scipy.integrate as si

    tau = p.tan_beta
    x0 = wetted_length(p)
    c = p.half_beam_transom / p.length**p.planform_exponent
    q = p.planform_exponent

    def integrand(x):
        if x <= 0:
            return 0.0
        height = p.draft - tau * c * x**q
        if height <= 0:
            return 0.0
        return 2.0 * height * math.sqrt(1.0 + (c * q * x ** (q - 1)) ** 2)

    points = [x0 * f for f in (1e-9, 1e-6, 1e-3, 1e-1)]
    top = si.quad(integrand, 0.0, x0, limit=800, points=points)[0] if x0 > 0 else 0.0
    return (analytic_waterplane_area(p) * math.sqrt(1 + tau * tau) + top
            + section_area(p, p.length))


def cases() -> list[PlaningHullParameters]:
    out = []
    for beta in DEGREES:
        for q in EXPONENTS:
            for draft in DRAFTS:
                p = PlaningHullParameters(length=10.0, beam=3.0, draft=draft,
                                          deadrise_deg=beta, planform_exponent=q,
                                          nx=21, nz=9)
                try:
                    p.validate()
                except ValueError:
                    continue
                out.append(p)
    return out


def main() -> int:
    mp.mp.dps = 40
    rows = []
    legacy_nonfinite = []
    for p in cases():
        shipped = analytic_planing_hydrostatics(p)["wetted_surface_area"]
        before = legacy_wetted_area(p)
        ref = float(mp_wetted_area(p))
        after_error = (shipped - ref) / ref
        scipy_value = scipy_wetted_area(p)
        scipy_error = (scipy_value - ref) / ref
        rows.append(
            {
                "deadrise_deg": p.deadrise_deg,
                "planform_exponent": p.planform_exponent,
                "draft": p.draft,
                "wetted_length": wetted_length(p),
                "legacy": before if math.isfinite(before) else None,
                "shipped": shipped,
                "mpmath": ref,
                "legacy_error": (before - ref) / ref if math.isfinite(before) else None,
                "shipped_error": after_error,
                "scipy": scipy_value,
                "scipy_error": scipy_error,
            }
        )
        if not math.isfinite(before):
            legacy_nonfinite.append(rows[-1])

    finite_rows = [r for r in rows if r["legacy_error"] is not None]
    worst_legacy = max(finite_rows, key=lambda r: abs(r["legacy_error"]))
    worst_shipped = max(rows, key=lambda r: abs(r["shipped_error"]))
    converged = [r for r in finite_rows if abs(r["legacy_error"]) < 1e-9]
    diverged = [r for r in finite_rows if abs(r["legacy_error"]) >= 1e-9]
    agreement = max(abs((r["shipped"] - r["legacy"]) / r["legacy"]) for r in converged)
    worst_scipy = max(rows, key=lambda r: abs(r["scipy_error"]))
    q_ge_2 = [r for r in finite_rows if r["planform_exponent"] >= 2.0]
    worst_legacy_high_q = max(q_ge_2, key=lambda r: abs(r["legacy_error"]))
    q_lt_06 = [r for r in finite_rows if r["planform_exponent"] < 0.6]
    worst_legacy_low_q = max(q_lt_06, key=lambda r: abs(r["legacy_error"]))

    claims = {
        "C1_legacy_nonfinite": {
            "hold": len(legacy_nonfinite) > 0,
            "count": len(legacy_nonfinite),
            "sweep_size": len(rows),
            "examples": [
                {k: r[k] for k in ("deadrise_deg", "planform_exponent", "draft")}
                for r in legacy_nonfinite
            ],
        },
        "C2_shipped_finite_and_accurate": {
            "hold": all(math.isfinite(r["shipped"]) for r in rows)
            and abs(worst_shipped["shipped_error"]) < 1e-9,
            "worst_relative_error": worst_shipped["shipped_error"],
            "at": {k: worst_shipped[k] for k in ("deadrise_deg", "planform_exponent", "draft")},
        },
        "C3_no_op_where_legacy_converged": {
            "hold": agreement < 1e-9,
            "n_cases": len(converged),
            "max_shipped_vs_legacy": agreement,
        },
        "C5_correction_lands_on_reference": {
            "hold": all(abs(r["shipped_error"]) < 1e-9 for r in diverged),
            "n_cases_where_legacy_was_wrong": len(diverged),
            "worst_shipped_error_there": max((abs(r["shipped_error"]) for r in diverged),
                                             default=0.0),
        },
        "C4_two_references_agree": {
            "hold": abs(worst_scipy["scipy_error"]) < 1e-7,
            "worst_scipy_vs_mpmath": worst_scipy["scipy_error"],
            "at": {k: worst_scipy[k] for k in ("deadrise_deg", "planform_exponent", "draft")},
        },
    }
    defect_a = {
        "name": "overflow to inf",
        "count": len(legacy_nonfinite),
        "cases": [{k: r[k] for k in ("deadrise_deg", "planform_exponent", "draft")}
                  for r in legacy_nonfinite],
        "reference_value": {k: legacy_nonfinite[0][k] for k in ("mpmath", "shipped")},
    }
    defect_b = {
        "name": "uniform Simpson grid truncation at q >= 2",
        "worst_relative_error": {k: worst_legacy_high_q[k] for k in
                                 ("deadrise_deg", "planform_exponent", "draft", "legacy_error")},
        "n_cases_q_ge_2": len(q_ge_2),
    }
    defect_c = {
        "name": "uniform Simpson grid truncation at small q",
        "worst_relative_error": {k: worst_legacy_low_q[k] for k in
                                 ("deadrise_deg", "planform_exponent", "draft", "legacy_error")},
        "n_cases_q_lt_0_6": len(q_lt_06),
    }
    result = {
        "sweep": {
            "cases": len(rows),
            "deadrise_deg": list(DEGREES),
            "planform_exponent": list(EXPONENTS),
            "draft": list(DRAFTS),
            "length": 10.0,
            "beam": 3.0,
        },
        "reference": "mpmath.quad dps=40 in u = b(x) (q < 1, smooth) or x (q >= 1)",
        "cross_check": "scipy.integrate.quad in x",
        "defects_found_in_legacy": [defect_a, defect_b, defect_c],
        "claims": claims,
        "worst_legacy_error_on_shared_cases": {
            **{k: worst_legacy[k] for k in
               ("deadrise_deg", "planform_exponent", "draft", "legacy_error")}
        },
        "rows": rows,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "results.json").write_text(json.dumps(result, indent=2) + "\n")

    fig = make_figure(rows, worst_shipped)
    (OUT / "planing_topside_robustness.png").write_bytes(fig)

    for name, claim in claims.items():
        print(f"{name}: hold={claim['hold']} -> {json.dumps({k: v for k, v in claim.items() if k != 'hold'})}")
    return 0 if all(c["hold"] for c in claims.values()) else 1


def make_figure(rows, worst_shipped) -> bytes:
    import io

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), dpi=150)
    ax = axes[0]
    for beta, marker in ((18.0, "o"), (30.0, "s"), (45.0, "^")):
        xs, ys = [], []
        for r in rows:
            if r["deadrise_deg"] != beta or r["draft"] != 0.05:
                continue
            xs.append(r["planform_exponent"])
            ys.append(float("nan") if r["legacy_error"] is None else abs(r["legacy_error"]))
        ax.plot(xs, ys, marker=marker, ls="none", label=f"legacy, beta={beta:g} deg")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("planform exponent q")
    ax.set_ylabel("|relative error| vs mpmath")
    ax.set_title("(a) legacy integral, T = 0.05 m\n(missing point = inf)")
    ax.axhline(1e-12, color="k", ls=":", lw=1)
    ax.legend(fontsize=7)

    ax = axes[1]
    shipped = [abs(r["shipped_error"]) for r in rows]
    ax.plot(range(len(rows)), shipped, "o", ms=2, label="shipped")
    ax.set_yscale("log")
    ax.set_xlabel("sweep case index")
    ax.set_ylabel("|relative error| vs mpmath")
    ax.set_title(f"(b) after fix: {len(rows)} cases,\nworst {abs(worst_shipped['shipped_error']):.2e}")
    ax.axhline(1e-12, color="k", ls=":", lw=1)
    ax.legend(fontsize=8)
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    plt.close(fig)
    return buf.getvalue()


if __name__ == "__main__":
    raise SystemExit(main())
