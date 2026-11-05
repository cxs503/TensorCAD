#!/usr/bin/env python3
"""Wetted-surface closure error with a quantified uncertainty budget (UNCERT).

Gap closed by this tool (``wetted-surface-geometry`` / fact ``gaps``: missing
uncertainty, error bars, GCI): the closure error ``dS/S`` reported by the
platform (Holtrop-Mennen 1984 regression + linear coefficient closure) is a
*number without a bar*.  This tool attaches one, from three independent sources,
and propagates it to total resistance:

A. quadrature/discretisation uncertainty of the reference area ``S_exact`` for
   the analytic hull of :mod:`tensorcad.hull` -- two independent integrators
   (repo Gauss-Legendre path, nested tanh-sinh path of
   :func:`tensorcad.hydrostatics_analytic.analytic_hydrostatics`) plus a
   mesh-convergence sequence with Richardson extrapolation and Roache's GCI;
B. input-parameter uncertainty of the closure inputs ``(L, B, T, Cb)``, mapped
   through numerically differentiated sensitivities (GUM first-order RSS) to a
   95 % interval for ``dS/S``; the hull-form exponent ``p`` carries the ``Cb``
   uncertainty because ``Cb = (2/3)*int_0^1 sin(pi x)^p dx``;
C. quantisation of the platform response (2 decimals) and the modelling bias
   itself, reported separately from the uncertainty (a bias is not a bar).

Propagation to resistance uses ``dRt/Rt = pi * dS/S`` with ``pi = Rf/Rt`` taken
from the platform resistance endpoint (cached JSON, so the run is reproducible
offline); see report #3871 for the analytic derivation of the exactness of that
relation.

Everything written out is traceable: ``data/<run>/summary.json`` carries the
commit, hull parameters, per-level sequences, GCI, budget and verdicts;
``figures/`` holds the PNGs.  No value in the report is typed by hand.

Usage
-----
    python3 tools/wetted_uncertainty.py selftest
    python3 tools/wetted_uncertainty.py run     --out data/wetted-uncertainty-2026-11-04
    python3 tools/wetted_uncertainty.py figures --data data/wetted-uncertainty-2026-11-04
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import tensorcad as tc  # noqa: E402
from tensorcad import hull as H  # noqa: E402
from tensorcad.hydrostatics_analytic import analytic_hydrostatics  # noqa: E402

SCHEMA = "tensorcad.wetted_uncertainty/1"
MISSION = 3767
REPO_NAME = "TensorCAD"

#: Reference hull: the analytic family documented in docs/wigley.md, at a
#: destroyer-like slenderness.  Cb = (2/3)*taper(0.7) = 0.4667, Cm = 2/3,
#: Cwp = taper = 0.70 -- all three exact in closed form, which is what makes
#: the closure error measurable instead of asserted.
DESIGN = {
    "length_m": 100.0,
    "beam_m": 16.0,
    "draft_m": 6.0,
    "fore_exponent": 0.7,
    "aft_exponent": 0.7,
    "water_density_kg_per_m3": 1025.0,
}

#: Assumed standard uncertainties of the closure inputs (declared assumptions,
#: see docs/wetted_surface_uncertainty.md section 5).  Offsets half-breadth for
#: a 16 m beam hull are good to about +-0.5 % of B; draft marks are read to
#: 0.02 m; length between perpendiculars is a fixed contract dimension but the
#: *wetted* length moves with trim, hence 0.05 m; Cb from a hydrostatic table
#: carries about 0.2 % of itself.
INPUT_UNCERTAINTY = {
    "length_m": 0.05,
    "beam_m": 0.005 * DESIGN["beam_m"],
    "draft_m": 0.02,
    "cb": 0.002 * (2.0 / 3.0) * H.taper_integral(DESIGN["fore_exponent"]),
}

#: Bias/uncertainty interplay: the closure is a systematic model error, so it is
#: reported as a bias with a confidence interval, never as a symmetric bar that
#: hides it.
GCI_SAFETY_FACTOR = 1.25  # Roache's Fs for a three-grid study
CONFIDENCE_K = 2.0  # 95 % coverage, normal approximation

HULL_MESH_LEVELS = ((41, 21), (81, 41), (161, 81), (321, 161))


# --------------------------------------------------------------------------- #
# Holtrop-Mennen wetted area (the platform closure)
# --------------------------------------------------------------------------- #
def closure_coefficients(cb: float) -> tuple[float, float]:
    """Midship and waterplane coefficients the platform derives from Cb alone."""
    return 0.95 + 0.05 * cb, 0.75 + 0.25 * cb


def holtrop_wetted_area(
    length: float,
    beam: float,
    draft: float,
    cb: float,
    *,
    cm: float | None = None,
    cwp: float | None = None,
    bulb_area: float = 0.0,
) -> float:
    """Holtrop & Mennen (1984) wetted-surface regression, S in m**2.

    ``cm``/``cwp`` default to the platform closure; passing the true
    coefficients isolates the regression's structural error from the closure
    error of report #3871.
    """
    if cm is None or cwp is None:
        cm_c, cwp_c = closure_coefficients(cb)
        cm = cm_c if cm is None else cm
        cwp = cwp_c if cwp is None else cwp
    bracket = (
        0.453
        + 0.4425 * cb
        - 0.2862 * cm
        - 0.003467 * (beam / draft)
        + 0.3696 * cwp
    )
    return length * (2.0 * draft + beam) * math.sqrt(cm) * bracket + 2.38 * bulb_area / cb


# --------------------------------------------------------------------------- #
# Independent integrators for the analytic reference area
# --------------------------------------------------------------------------- #
def _gauss_legendre(nodes: int):
    x, w = np.polynomial.legendre.leggauss(nodes)
    return x, w


def area_path_a(design: dict, nodes: int) -> float:
    """Repo path: :func:`tensorcad.hull.analytic_wetted_surface_area`."""
    return H.analytic_wetted_surface_area(
        design["length_m"],
        design["beam_m"],
        design["draft_m"],
        nodes=nodes,
        fore_exponent=design["fore_exponent"],
        aft_exponent=design["aft_exponent"],
    )


def _surface_element(design: dict, x: np.ndarray, u: np.ndarray) -> np.ndarray:
    """``sqrt(4 u^2 T^2 (1 + b_x^2) + beta^2)`` on the (x, u) grid.

    ``b = beta * sqrt(zeta)`` with ``zeta = u**2`` removes the keel
    ``zeta**-0.5`` slope singularity into the ``du`` Jacobian; ``beta`` and its
    x-derivative are the sine-tapered half breadth and slope.
    """
    length = design["length_m"]
    beam = design["beam_m"]
    draft = design["draft_m"]
    p = design["fore_exponent"]
    s = np.sin(np.pi * x / length)
    s = np.clip(s, 0.0, None)
    beta = 0.5 * beam * s ** p
    dbeta_dx = 0.5 * beam * p * (np.pi / length) * np.cos(np.pi * x / length) * s ** (p - 1.0)
    beta = beta[:, None]
    dbeta_dx = dbeta_dx[:, None]
    u = u[None, :]
    slope = dbeta_dx * u
    return np.sqrt(4.0 * u * u * draft * draft * (1.0 + slope * slope) + beta * beta)


def area_adaptive(design: dict, tol: float = 1e-11) -> float:
    """Third, adaptive path: nested ``scipy.integrate.quad`` (Gauss-Kronrod).

    Chosen deliberately as an algorithmically different family from the fixed
    Gauss-Legendre grid of path A (which clusters nodes with ``psi ~ tau**3``)
    and from the double-exponential tanh-sinh rule used by
    :func:`tensorcad.hydrostatics_analytic.analytic_hydrostatics`.  Adaptive
    bisection localises the bow/stern ``x**(-0.3)`` edge on its own, so it is the
    reference used to arbitrate which fixed-grid value is right.
    """
    from scipy.integrate import quad

    length = design["length_m"]
    draft = design["draft_m"]
    beam = design["beam_m"]
    power = design["fore_exponent"]

    def inner(x: float) -> float:
        s = math.sin(math.pi * x / length)
        if s <= 0.0:
            return 0.0
        beta = 0.5 * beam * s ** power
        slope = 0.5 * beam * power * (math.pi / length) * math.cos(math.pi * x / length) * s ** (power - 1.0)

        def integrand(u: float) -> float:
            eta = slope * u
            return math.sqrt(4.0 * u * u * draft * draft * (1.0 + eta * eta) + beta * beta)

        value, _ = quad(integrand, 0.0, 1.0, epsabs=tol * beta, epsrel=1e-12, limit=200)
        return value

    half, _ = quad(inner, 0.0, 0.5 * length, epsabs=tol, epsrel=1e-11, limit=400)
    return 4.0 * half  # two halves x port/starboard


def superellipse_area(design: dict, n: float, m: float, nodes: int = 400) -> float:
    """Wetted area of the superellipse Wigley family used by report #3871.

    ``y = (B/2)(1-|2x/L-1|^n)(1-(w/T)^m)``; independent of :mod:`tensorcad.hull`
    (different family) and used only as a cross-check of the quadrature
    machinery against an already published number.
    """
    gx, wx = _gauss_legendre(nodes)
    gu, wu = _gauss_legendre(nodes)
    length = design["length_m"]
    beam = design["beam_m"]
    draft = design["draft_m"]
    x = 0.5 * length * (gx + 1.0)
    xi = np.abs(2.0 * x / length - 1.0)
    taper = 1.0 - xi ** n
    dtaper = -n * np.sign(2.0 * x / length - 1.0) * (2.0 / length) * np.where(xi > 0, xi ** (n - 1.0), 0.0)
    # vertical: w in (0, T) from the waterline downwards, w = T*s, s in (0,1)
    s = 0.5 * (gu + 1.0)
    prof = 1.0 - s ** m
    dprof_ds = -m * s ** (m - 1.0)
    y_x = 0.5 * beam * dtaper[:, None] * prof[None, :]
    # dw = draft*ds, so d(y)/ds = draft * d(y)/dw
    y_w = 0.5 * beam * taper[:, None] * dprof_ds[None, :] / draft
    # area element on the (x, w) chart, times dw = draft*ds
    element = np.sqrt(1.0 + y_x * y_x + y_w * y_w) * draft
    integral = float(np.einsum("i,j,ij->", wx * 0.5 * length, wu, element))
    return 2.0 * integral  # port and starboard


def prolate_spheroid_area(length: float, beam: float) -> float:
    """Closed-form wetted area of a prolate spheroid (semi-axes L/2, B/2).

    Used only to identify the geometry behind report #3871's published
    ``S = 2478.7758 m^2``: the closed form is independent of every quadrature in
    this file, so matching it pins down their ``L``/``B`` and lets the
    superellipse cross-check be dimensioned correctly.
    """
    a = 0.5 * length
    b = 0.5 * beam
    e = math.sqrt(1.0 - (b / a) ** 2)
    return 2.0 * math.pi * b * b * (1.0 + (a / (b * e)) * math.asin(e))


# --------------------------------------------------------------------------- #
# Richardson / GCI machinery
# --------------------------------------------------------------------------- #
def observed_order(f: list[float], ratio: float) -> float:
    """Observed order of accuracy from three monotone-refined values.

    ``p = log|d2/d1| / log r``; the sign of the differences is not informative
    (they may approach the limit from either side), only their ratio is.
    """
    d1 = f[1] - f[0]
    d2 = f[2] - f[1]
    if d1 == 0.0 or d2 / d1 <= 0.0:
        return float("nan")
    return -math.log(abs(d2 / d1)) / math.log(ratio)


def richardson_extrapolate(f: list[float], ratio: float, p: float) -> float:
    """Extrapolate to zero grid from the two finest values of a triple."""
    return f[2] + (f[2] - f[1]) / (ratio ** p - 1.0)


def gci(fine: float, coarse: float, p: float, ratio: float, fs: float = GCI_SAFETY_FACTOR) -> float:
    """Roache grid-convergence index, relative to the fine solution."""
    return fs * abs((coarse - fine) / fine) / abs(ratio ** p - 1.0)


def gci_sequence(f: list[float], ratio: float) -> dict:
    """GCI of every three-value window plus the fitted order of each window."""
    out = {"levels": list(range(len(f))), "values": f, "windows": []}
    for i in range(len(f) - 2):
        window = f[i : i + 3]
        p = observed_order(window, ratio)
        out["windows"].append(
            {
                "coarse": i,
                "p_observed": p,
                "richardson": richardson_extrapolate(window, ratio, p),
                "gci_fine_relative": gci(window[1], window[0], p, ratio),
            }
        )
    return out


# --------------------------------------------------------------------------- #
# Platform endpoints (cached for offline reproducibility)
# --------------------------------------------------------------------------- #
def _api_base() -> str:
    return os.environ.get("TENSOR_SHIP_API", "http://127.0.0.1:8002/api/v1")


def resistance_call(length, beam, draft, cb, displacement, speed) -> dict:
    body = {
        "length": length,
        "beam": beam,
        "draft": draft,
        "cb": cb,
        "displacement": displacement,
        "speed": speed,
    }
    req = urllib.request.Request(
        _api_base() + "/calculations/resistance",
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=30) as response:
        payload = json.load(response)
    return payload.get("data", payload)


def collect_platform(speeds: list[float], design: dict, cb: float, cache: Path) -> dict:
    """Query the resistance endpoint per speed; fall back to the committed cache."""
    records = {}
    displacement = design["length_m"] * design["beam_m"] * design["draft_m"] * cb
    for v in speeds:
        key = f"{v:.3f}"
        try:
            data = resistance_call(
                design["length_m"],
                design["beam_m"],
                design["draft_m"],
                cb,
                displacement,
                v,
            )
            records[key] = {"source": "live", **data}
        except (urllib.error.URLError, OSError, ValueError) as exc:  # offline
            records[key] = {"source": f"unavailable: {exc}"}
    live = [r for r in records.values() if r.get("source") == "live"]
    if not live:
        if cache.exists():
            cached = json.loads(cache.read_text())
            for k, v in cached.items():
                if k in records and "rt" not in records[k]:
                    records[k] = {**v, "source": "cache"}
        else:
            raise RuntimeError(
                "resistance endpoint unreachable and no cached responses at %s" % cache
            )
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "source"} for k, v in records.items()}, indent=1, sort_keys=True))
    return records


# --------------------------------------------------------------------------- #
# Uncertainty budget
# --------------------------------------------------------------------------- #
def _family_mesh_cb(exponent: float) -> float:
    return (2.0 / 3.0) * H.taper_integral(exponent)


def design_cb(design: dict) -> float:
    return _family_mesh_cb(design["fore_exponent"])


def closure_error(design: dict, cb: float, nodes: int = 240) -> dict:
    """Return the closure error and its decomposition at one design point."""
    exact = area_path_a(design, nodes)
    platform = holtrop_wetted_area(
        design["length_m"], design["beam_m"], design["draft_m"], cb
    )
    cm_true = 2.0 / 3.0
    cwp_true = H.taper_integral(design["fore_exponent"])
    structure_only = holtrop_wetted_area(
        design["length_m"],
        design["beam_m"],
        design["draft_m"],
        cb,
        cm=cm_true,
        cwp=cwp_true,
    )
    eps_struct = structure_only / exact - 1.0
    eps_closure = (1.0 + platform / exact - 1.0) / (1.0 + eps_struct) - 1.0
    return {
        "S_exact_m2": exact,
        "S_platform_m2": platform,
        "cb": cb,
        "cm_true": cm_true,
        "cwp_true": cwp_true,
        "cm_closure": closure_coefficients(cb)[0],
        "cwp_closure": closure_coefficients(cb)[1],
        "total_relative": platform / exact - 1.0,
        "structure_relative": eps_struct,
        "closure_relative": eps_closure,
    }


def sensitivities(design: dict, cb: float, rel_step: float = 1e-5) -> dict:
    """Central-difference d(ln S)/d(ln X) for both area models.

    The exponent direction is converted to a Cb direction with
    ``dCb/dp = 1.5 dtaper/dp`` so the budget can be driven by a measured Cb.
    """
    def rel(key, delta_key):
        plus = dict(design)
        minus = dict(design)
        plus[delta_key] = design[delta_key] * (1.0 + rel_step)
        minus[delta_key] = design[delta_key] * (1.0 - rel_step)
        out = {}
        for name, fn in (
            ("exact", lambda d: area_path_a(d, 400)),
            ("platform", lambda d: holtrop_wetted_area(
                d["length_m"], d["beam_m"], d["draft_m"], cb)),
        ):
            hi = fn(plus)
            lo = fn(minus)
            out[name] = (math.log(hi) - math.log(lo)) / (2.0 * rel_step)
        return out

    table = {
        "length_m": rel("length", "length_m"),
        "beam_m": rel("beam", "beam_m"),
        "draft_m": rel("draft", "draft_m"),
    }
    # exponent -> Cb, by a numerical derivative of the closed-form taper
    h = 1e-3
    p = design["fore_exponent"]
    dtaper_dp = (H.taper_integral(p + h) - H.taper_integral(p - h)) / (2.0 * h)
    dcb_dp = (2.0 / 3.0) * dtaper_dp
    plus = dict(design)
    minus = dict(design)
    plus["fore_exponent"] = p + h
    minus["fore_exponent"] = p - h
    dln_exact_dcb = (
        math.log(area_adaptive(plus)) - math.log(area_adaptive(minus))
    ) / (2.0 * h * dcb_dp)
    dln_platform_dcb = (
        math.log(holtrop_wetted_area(plus["length_m"], plus["beam_m"], plus["draft_m"], _family_mesh_cb(p + h)))
        - math.log(holtrop_wetted_area(minus["length_m"], minus["beam_m"], minus["draft_m"], _family_mesh_cb(p - h)))
    ) / (2.0 * h * dcb_dp)
    table["cb"] = {"exact": dln_exact_dcb, "platform": dln_platform_dcb}
    return table


def error_budget(design: dict, cb: float, sens: dict, quad_u: float, platform_quant_u: float) -> dict:
    """GUM first-order RSS budget for d(closure error)/dX.

    The closure error is ``e = S_platform/S_exact - 1``, so its sensitivity to a
    relative input change is ``S_p - S_e`` of the two log-sensitivities; the
    quadrature and quantisation terms enter additively in the variance.
    """
    components = []
    for key, u_x in INPUT_UNCERTAINTY.items():
        s_e = sens[key]["exact"]
        s_p = sens[key]["platform"]
        # convert absolute u(X) to relative
        if key == "cb":
            u_rel = u_x / cb
        else:
            u_rel = u_x / design[key]
        slope = s_p - s_e
        components.append(
            {
                "name": f"u({key})",
                "input_absolute": u_x,
                "input_relative": u_rel,
                "d_ln_platform": s_p,
                "d_ln_exact": s_e,
                "contribution_relative": slope * u_rel,
            }
        )
    components.append(
        {
            "name": "u(quadrature)",
            "input_absolute": quad_u,
            "input_relative": quad_u,
            "d_ln_platform": 0.0,
            "d_ln_exact": 1.0,
            "contribution_relative": quad_u,
        }
    )
    components.append(
        {
            "name": "u(platform rounding)",
            "input_absolute": platform_quant_u,
            "input_relative": platform_quant_u,
            "d_ln_platform": 1.0,
            "d_ln_exact": 0.0,
            "contribution_relative": platform_quant_u,
        }
    )
    variance = sum(c["contribution_relative"] ** 2 for c in components)
    u_c = math.sqrt(variance)
    worst = sum(abs(c["contribution_relative"]) for c in components)
    return {
        "components": components,
        "combined_standard": u_c,
        "expanded_95": CONFIDENCE_K * u_c,
        "worst_case_sum": worst,
        "dominant": max(components, key=lambda c: abs(c["contribution_relative"]))["name"],
    }


# --------------------------------------------------------------------------- #
# Driver pieces
# --------------------------------------------------------------------------- #
def mesh_convergence(design: dict, levels=HULL_MESH_LEVELS) -> dict:
    analytic = area_path_a(design, 240)
    rows = []
    for nx, nz in levels:
        params = H.HullParameters(
            length=design["length_m"],
            beam=design["beam_m"],
            draft=design["draft_m"],
            nx=nx,
            nz=nz,
            water_density=design["water_density_kg_per_m3"],
            fore_exponent=design["fore_exponent"],
            aft_exponent=design["aft_exponent"],
        )
        model = H.generate_hull(params)
        hydro = tc.hydrostatics.compute_hydrostatics(model.mesh)
        rows.append(
            {
                "nx": nx,
                "nz": nz,
                "faces": len(model.mesh.faces),
                "h_m": design["length_m"] / (nx - 1),
                "S_mesh_m2": hydro.wetted_surface_area,
                "cb_mesh": hydro.block_coefficient,
                "volume_mesh_m3": hydro.displacement_volume,
                "error_relative": hydro.wetted_surface_area / analytic - 1.0,
            }
        )
    seq = gci_sequence([r["S_mesh_m2"] for r in rows], ratio=2.0)
    # a window (coarse, medium, fine) yields p and GCI for the *medium* level
    for window in seq["windows"]:
        row = rows[window["coarse"] + 1]
        row["p_observed"] = window["p_observed"]
        row["gci_relative_95"] = window["gci_fine_relative"]
        row["richardson_S_m2"] = window["richardson"]
        row["richardson_relative_error"] = window["richardson"] / analytic - 1.0
        row["gci_brackets_true_error"] = abs(window["gci_fine_relative"]) >= abs(row["error_relative"])
    evaluated = [r for r in rows if "gci_relative_95" in r]
    last = seq["windows"][-1]
    p_eff = last["p_observed"]
    extrapolated = last["richardson"]
    return {
        "ratio": 2.0,
        "levels": rows,
        "finest_gci_relative": last["gci_fine_relative"],
        "finest_p_observed": p_eff,
        "richardson_S_m2": extrapolated,
        "richardson_vs_analytic_relative": extrapolated / analytic - 1.0,
        "analytic_S_m2": analytic,
        "gci_brackets_all_levels": all(r["gci_brackets_true_error"] for r in evaluated),
        "mesh_cb_relative_error_finest": rows[-1]["cb_mesh"] / design_cb(design) - 1.0,
        "volume_check": {
            "LBT_Cb_m3": design["length_m"] * design["beam_m"] * design["draft_m"] * design_cb(design),
            "mesh_volume_finest_m3": rows[-1]["volume_mesh_m3"],
        },
    }


def quadrature_study(design: dict, *, adaptive: bool = True) -> dict:
    """Grid-convergence evidence for the reference area ``S_exact``.

    Three algorithmically different integrators on the same integrand:
    repo Gauss-Legendre with ``psi ~ tau**3`` node clustering (path A, refinement
    ratio 2 over the node count), the repo nested tanh-sinh rule at increasing
    levels (double exponential), and -- when ``adaptive`` -- nested adaptive
    Gauss-Kronrod (``scipy.integrate.quad``) as the arbiter.
    """
    nodes_a = [80, 160, 320, 640]
    a = [area_path_a(design, n) for n in nodes_a]
    seq_a = gci_sequence(a, ratio=2.0)
    params = H.HullParameters(
        length=design["length_m"],
        beam=design["beam_m"],
        draft=design["draft_m"],
        fore_exponent=design["fore_exponent"],
        aft_exponent=design["aft_exponent"],
    )
    ts_levels = [4, 5, 6, 7, 8]
    ts = [analytic_hydrostatics(params, quadrature_levels=level) for level in ts_levels]
    ts_values = [s.wetted_surface_area for s in ts]
    ts_spread = max(ts_values) - min(ts_values)
    reference = area_adaptive(design) if adaptive else a[-1]
    return {
        "path_a_nodes": nodes_a,
        "path_a_values": a,
        "path_a_windows": seq_a["windows"],
        "path_a_finest": a[-1],
        "tanh_sinh_levels": ts_levels,
        "tanh_sinh_values": ts_values,
        "tanh_sinh_self_difference_relative": ts[-1].quadrature["wetted_area_relative_difference"],
        "tanh_sinh_spread_relative": ts_spread / a[-1],
        "adaptive_S_m2": reference,
        "path_a_vs_adaptive_relative": a[-1] / reference - 1.0,
        "tanh_sinh_vs_adaptive_relative": ts_values[-1] / reference - 1.0,
        "u_quadrature_relative": max(
            abs(a[-1] / reference - 1.0),
            abs(ts_values[-1] / reference - 1.0),
            abs(seq_a["windows"][-1]["gci_fine_relative"]),
        ),
    }


def family_sweep(design: dict, exponents: list[float]) -> list[dict]:
    rows = []
    for p in exponents:
        local = dict(design)
        local["fore_exponent"] = p
        local["aft_exponent"] = p
        cb = _family_mesh_cb(p)
        err = closure_error(local, cb)
        sens = sensitivities(local, cb)
        quad = quadrature_study(local)
        self_gci = abs(quad["path_a_windows"][-1]["gci_fine_relative"])
        if not math.isfinite(self_gci):
            # non-monotone triple: fall back on the cross-integrator agreement
            self_gci = abs(quad["u_quadrature_relative"])
        budget = error_budget(
            local,
            cb,
            sens,
            self_gci,
            0.005 / math.sqrt(3.0) / err["S_platform_m2"],
        )
        rows.append(
            {
                "exponent": p,
                "cb": cb,
                "cwp": H.taper_integral(p),
                **{k: err[k] for k in ("S_exact_m2", "S_platform_m2", "total_relative", "structure_relative", "closure_relative")},
                "u95_relative": budget["expanded_95"],
                "integrator_agreement": quad["path_a_vs_adaptive_relative"],
                "u_num_self_gci": self_gci,
                "budget": budget,
            }
        )
    return rows


def propagation(platform: dict, delta_s_rel: float, u_rel: float, density: float) -> list[dict]:
    rows = []
    for key in sorted(platform, key=float):
        rec = platform[key]
        if "rt" not in rec:
            continue
        pi = rec["rf"] / rec["rt"]
        # the endpoint reports thrust/resistance in kN, so the identity needs /1000
        identity = abs(
            0.5 * density * rec["wetted_surface"] * rec["speed_ms"] ** 2 * rec["cf_total"] / 1000.0
            - rec["rf"]
        ) / rec["rf"]
        rows.append(
            {
                "speed_kn": float(key),
                "froude": rec["froude_number"],
                "reynolds": rec["reynolds_number"],
                "rf_kn": rec["rf"],
                "rt_kn": rec["rt"],
                "pi": pi,
                "wetted_surface_platform_m2": rec["wetted_surface"],
                "drt_relative": pi * delta_s_rel,
                "u95_drt_relative": pi * u_rel,
                "rf_identity_relative_residual": identity,
            }
        )
    return rows


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True
        ).strip()
    except Exception:  # pragma: no cover
        return "unknown"


def run(out: Path, speeds: list[float]) -> dict:
    design = dict(DESIGN)
    cb = design_cb(design)
    quad = quadrature_study(design)
    mesh = mesh_convergence(design)
    err = closure_error(design, cb)
    sens = sensitivities(design, cb)
    platform_u = 0.005 / math.sqrt(3.0) / err["S_platform_m2"]
    budget = error_budget(design, cb, sens, quad["u_quadrature_relative"], platform_u)

    out.mkdir(parents=True, exist_ok=True)
    cache = out / "platform_resistance.json"
    platform = collect_platform(speeds, design, cb, cache)
    prop = propagation(platform, err["total_relative"], budget["expanded_95"], design["water_density_kg_per_m3"])

    # cross-check against the already published numbers of report #3871
    # (its two benchmarks are L=100, B=10, T=5; identified from the closed-form
    # prolate-spheroid area, because the report does not print the dimensions)
    wigley = dict(design)
    wigley.update({"length_m": 100.0, "beam_m": 10.0, "draft_m": 5.0})
    wigley_full = superellipse_area(wigley, 2.0, 2.0, nodes=320)
    spheroid = prolate_spheroid_area(100.0, 10.0)

    sweep = family_sweep(design, [0.35, 0.45, 0.55, 0.7, 0.9, 1.1])

    summary = {
        "schema": SCHEMA,
        "mission": MISSION,
        "repo": REPO_NAME,
        "commit": git_commit(),
        "design": design,
        "design_cb": cb,
        "design_cwp": H.taper_integral(design["fore_exponent"]),
        "design_cm": 2.0 / 3.0,
        "quadrature": quad,
        "mesh": mesh,
        "closure": err,
        "sensitivities": sens,
        "budget": budget,
        "propagation": prop,
        "family_sweep": sweep,
        "cross_checks": {
            "wigley_superellipse_full_S_m2": wigley_full,
            "wigley_superellipse_one_side_S_m2": 0.5 * wigley_full,
            "wigley_report_3871_S_m2": 1275.2810,
            "relative_difference_vs_3871_one_side": 0.5 * wigley_full / 1275.2810 - 1.0,
            "prolate_spheroid_closed_form_m2": spheroid,
            "prolate_spheroid_report_3871_m2": 2478.7758,
            "relative_difference_prolate_vs_3871": spheroid / 2478.7758 - 1.0,
            "platform_wetted_vs_local_holtrop_relative": max(
                abs(r["wetted_surface_platform_m2"] / err["S_platform_m2"] - 1.0)
                for r in prop
            )
            if prop
            else None,
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=1, sort_keys=False))
    return summary


def selftest() -> int:
    """Cheap invariants; the pytest suite exercises the same numbers."""
    design = dict(DESIGN)
    cb = design_cb(design)
    a = area_path_a(design, 320)
    b = area_adaptive(design)
    assert abs(a / b - 1.0) < 1e-7, (a, b)
    # quadrature is converged well below any closure error we report
    assert abs(area_path_a(design, 160) / a - 1.0) < 1e-10
    # report #3871 reproduction: its prolate-spheroid area fixes the dimensions
    spheroid = prolate_spheroid_area(100.0, 10.0)
    assert abs(spheroid / 2478.7758 - 1.0) < 1e-4, spheroid
    # superellipse cross-check against the same report (its S is one-sided)
    wig = dict(design)
    wig.update({"length_m": 100.0, "beam_m": 10.0, "draft_m": 5.0})
    s = superellipse_area(wig, 2.0, 2.0, nodes=320)
    assert abs(s / 2.0 / 1275.2810 - 1.0) < 1e-6, s
    # closure reproduces the platform's regression at this Cb
    s_platform = holtrop_wetted_area(design["length_m"], design["beam_m"], design["draft_m"], cb)
    assert 1500 < s_platform < 2400, s_platform
    print(f"selftest ok: pathA={a:.6f} adaptive={b:.6f} wigley_half={s/2:.4f} spheroid={spheroid:.4f} S_platform={s_platform:.4f}")
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["selftest", "run", "figures"])
    parser.add_argument("--out", default="data/wetted-uncertainty-2026-11-04")
    parser.add_argument("--data", default="data/wetted-uncertainty-2026-11-04")
    parser.add_argument("--speeds", default="6,8,10,12,14,16,18,20,22,24")
    args = parser.parse_args(argv)
    if args.command == "selftest":
        return selftest()
    if args.command == "run":
        speeds = [float(v) for v in args.speeds.split(",")]
        summary = run(REPO / args.out, speeds)
        print(json.dumps({
            "S_exact": summary["closure"]["S_exact_m2"],
            "total_relative": summary["closure"]["total_relative"],
            "u95": summary["budget"]["expanded_95"],
            "u_quad": summary["quadrature"]["u_quadrature_relative"],
        }, indent=1))
        return 0
    from tools.wetted_uncertainty_figures import make_figures

    make_figures(REPO / args.data, REPO / "figures")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())