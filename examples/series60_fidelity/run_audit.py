"""Audit the Series 60 family against its published particulars, and close the gap.

Two questions, both answered by numbers this script computes rather than quotes:

1. **Does the legacy family describe a Series 60 hull?**  The legacy loft uses a single
   longitudinal exponent, which forces ``Cwp = Cp`` and ``LCB = L/2`` identically.  The
   script measures what that costs against the tabulated parent forms and against an
   independent regression, and states the verdict per block coefficient.
2. **Does the ``cwp``/``lcb`` extension fix it?**  It re-measures the same quantities on
   hulls whose waterplane coefficient and longitudinal centre of buoyancy are *inputs*,
   and refines the grid to show the loft still converges to its own closed forms.

Run from the repository root::

    python examples/series60_fidelity/run_audit.py --out artifacts/series60_fidelity
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.series60 import (  # noqa: E402
    Series60Parameters,
    analytic_series60_hydrostatics,
    gauss_legendre_unit,
    generate_series60,
    series60_section_law,
)

L, B, T = 10.0, 1.5, 0.6

#: Tabulated Series 60 parent-form particulars (Todd 1951, as reproduced in the standard
#: resistance references): midship coefficient, prismatic coefficient and waterplane
#: coefficient.  ``Cb = Cp * Cm`` holds in the table to the digit shown, which is a
#: self-consistency check on the transcription (see ``check_table``).
PUBLISHED_PARENTS: dict[float, dict[str, float]] = {
    0.60: {"cm": 0.960, "cp": 0.625, "cwp": 0.700},
    0.65: {"cm": 0.970, "cp": 0.670, "cwp": 0.727},
    0.70: {"cm": 0.980, "cp": 0.715, "cwp": 0.753},
    0.75: {"cm": 0.985, "cp": 0.760, "cwp": 0.780},
    0.80: {"cm": 0.990, "cp": 0.805, "cwp": 0.807},
}

#: Published parents carry their LCB aft of midship.  These are the order-of-magnitude
#: values quoted with the same tables (per cent of L, negative = aft of midship); they
#: are used as demonstration targets, not as verification references -- the point the
#: audit makes about LCB is structural (the legacy family can only produce 0.0), so it
#: does not depend on the second digit of these numbers.
PUBLISHED_LCB_PERCENT_L: dict[float, float] = {
    0.60: -2.5,
    0.65: -1.7,
    0.70: -1.0,
    0.75: -0.4,
    0.80: +0.2,
}


def schneekluth_cwp(cb: float) -> float:
    """``Cwp = (1 + 2 Cb) / 3`` -- Schneekluth & Bertram's regression."""
    return (1.0 + 2.0 * cb) / 3.0


def check_table() -> dict[str, float]:
    """How far the published table itself is from ``Cb = Cp * Cm``."""
    return {f"{cb:.2f}": cb - p["cp"] * p["cm"] for cb, p in PUBLISHED_PARENTS.items()}


def legacy_audit(verbose: bool = True) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for cb, published in PUBLISHED_PARENTS.items():
        p = Series60Parameters(length=L, beam=B, draft=T, cb=cb)
        ref = analytic_series60_hydrostatics(p)
        cwp = ref["waterplane_coefficient"]
        rows.append({
            "cb": cb,
            "cb_closed_form": ref["block_coefficient"],
            "cp": ref["prismatic_coefficient"],
            "cm": ref["midship_coefficient"],
            "cwp_model": cwp,
            "cwp_published": published["cwp"],
            "cwp_error_vs_published": (cwp - published["cwp"]) / published["cwp"],
            "cwp_schneekluth": schneekluth_cwp(cb),
            "cwp_error_vs_schneekluth": (cwp - schneekluth_cwp(cb)) / schneekluth_cwp(cb),
            "lcb_over_length": ref["longitudinal_centre_buoyancy"] / L,
            "lcf_over_length": ref["longitudinal_centre_flotation"] / L,
            "vcb_over_draft": ref["vertical_centre_buoyancy"] / T,
        })
        if verbose:
            print(f"  Cb={cb:.2f}  Cwp(model)={cwp:.4f}  published={published['cwp']:.4f} "
                  f"({rows[-1]['cwp_error_vs_published']:+.2%})  "
                  f"Schneekluth={rows[-1]['cwp_schneekluth']:.4f} "
                  f"({rows[-1]['cwp_error_vs_schneekluth']:+.2%})  LCB/L={rows[-1]['lcb_over_length']:.3f}")
    return rows


def corrected_audit(verbose: bool = True) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for cb, published in PUBLISHED_PARENTS.items():
        lcb_target = 0.5 + PUBLISHED_LCB_PERCENT_L[cb] / 100.0
        p = Series60Parameters(length=L, beam=B, draft=T, cb=cb,
                               cwp=published["cwp"], lcb=lcb_target)
        law = series60_section_law(p)
        ref = analytic_series60_hydrostatics(p)
        rows.append({
            "cb": cb,
            "n_aft": law.n_aft,
            "n_fore": law.n_fore,
            "section_taper_gain": law.gain,
            "cm_mid": law.cm_mid,
            "cwp_model": ref["waterplane_coefficient"],
            "cwp_requested": published["cwp"],
            "cwp_error_vs_requested": (ref["waterplane_coefficient"] - published["cwp"]) / published["cwp"],
            "lcb_over_length": law.lcb_ratio,
            "lcb_target": lcb_target,
            "lcb_error": law.lcb_ratio - lcb_target,
            "lcf_over_length": law.lcf_ratio,
            "vcb_over_draft": law.vcb_ratio,
        })
        if verbose:
            print(f"  Cb={cb:.2f}  n_aft={law.n_aft:6.3f} n_fore={law.n_fore:6.3f} "
                  f"gain={law.gain:+.4f}  Cwp={ref['waterplane_coefficient']:.6f} "
                  f"({rows[-1]['cwp_error_vs_requested']:+.2e})  LCB/L={law.lcb_ratio:.6f} "
                  f"({rows[-1]['lcb_error']:+.2e})  VCB/T={law.vcb_ratio:.5f}")
    return rows


# ---------------------------------------------------------------------------- independent check
def _tanh_sinh(fn, a: float = 0.0, b: float = 1.0, levels: int = 7) -> float:
    """Double-exponential (tanh-sinh) quadrature.

    Chosen for the independent check because it is *not* the module's Gauss-Legendre
    rule and because it resolves the two things that break a uniform rule here: the
    ``zeta**(1/m)`` keel singularity (infinite slope at ``zeta = 0``) and the fact that
    for a full hull the section area is concentrated in a layer of thickness ``1/m``
    under the waterline.
    """
    half = 0.5 * (b - a)
    mid = 0.5 * (a + b)
    h = 1.0 / (2 ** levels)
    total = 0.0
    k = 0
    while True:
        t = k * h
        u = 0.5 * math.pi * math.sinh(t)
        if u > 40.0:
            break
        ch = math.cosh(u)
        w = 0.5 * math.pi * math.cosh(t) / (ch * ch)
        # tanh(u) is odd, so both signs of t have to be summed: doubling the positive
        # half only is right for an even integrand, and the hull is not symmetric.
        total += w * fn(mid + half * math.tanh(u))
        if k:
            total += w * fn(mid - half * math.tanh(u))
        k += 1
    return total * half * h


def independent_hydrostatics(p: Series60Parameters) -> dict[str, float]:
    """Integrate the analytic surface from its definition, by a different rule.

    No closed form from the module is used: the half-breadth ``b(xi, zeta) =
    (B/2) f(xi) zeta**(1/m(xi))`` is integrated directly, vertically by tanh-sinh
    (keel singularity) and longitudinally by tanh-sinh (|2xi-1| kink at midship).
    """
    law = series60_section_law(p)

    def vertical(xi: float) -> tuple[float, float]:
        """``(int zeta**(1/m) dzeta, int (1-zeta) zeta**(1/m) dzeta)`` at one station."""
        q = law.fullness(xi)
        m = q / (1.0 - q)
        area = _tanh_sinh(lambda z: math.pow(z, 1.0 / m))
        moment = _tanh_sinh(lambda z: (1.0 - z) * math.pow(z, 1.0 / m))
        return area, moment

    def volume_density(xi: float) -> float:
        return law.deck_taper(xi) * vertical(xi)[0]

    def vertical_moment_density(xi: float) -> float:
        return law.deck_taper(xi) * vertical(xi)[1]

    def wpa_density(xi: float) -> float:
        return law.deck_taper(xi)

    def wpa_moment_density(xi: float) -> float:
        return xi * law.deck_taper(xi)

    vol_hatch = _tanh_sinh(volume_density)
    vert_hatch = _tanh_sinh(vertical_moment_density)
    awp_hatch = _tanh_sinh(wpa_density)
    awp_moment = _tanh_sinh(wpa_moment_density)
    return {
        "displacement_volume": L * B * T * vol_hatch,
        "waterplane_area": L * B * awp_hatch,
        "vertical_centre_buoyancy": T * (1.0 - vert_hatch / vol_hatch),
        "longitudinal_centre_flotation": L * awp_moment / awp_hatch,
        "block_coefficient": vol_hatch,
        "waterplane_coefficient": awp_hatch,
    }


def cross_check(verbose: bool = True) -> list[dict[str, object]]:
    cases = [
        ("legacy Cb=0.60", Series60Parameters(length=L, beam=B, draft=T, cb=0.60)),
        ("legacy Cb=0.80", Series60Parameters(length=L, beam=B, draft=T, cb=0.80)),
        ("cwp=0.700 Cb=0.60", Series60Parameters(length=L, beam=B, draft=T, cb=0.60, cwp=0.700)),
        ("cwp=0.700 lcb=0.475", Series60Parameters(length=L, beam=B, draft=T, cb=0.60,
                                                  cwp=0.700, lcb=0.475)),
        ("cwp=0.780 lcb=0.496", Series60Parameters(length=L, beam=B, draft=T, cb=0.75,
                                                  cwp=0.780, lcb=0.496)),
    ]
    rows: list[dict[str, object]] = []
    for name, p in cases:
        ref = analytic_series60_hydrostatics(p)
        ind = independent_hydrostatics(p)
        row = {"case": name}
        for key in ("displacement_volume", "waterplane_area", "vertical_centre_buoyancy",
                    "longitudinal_centre_flotation"):
            row[f"{key}_closed"] = ref[key]
            row[f"{key}_independent"] = ind[key]
            row[f"{key}_error"] = (ref[key] - ind[key]) / ind[key]
        rows.append(row)
        if verbose:
            worst = max(abs(row[f"{k}_error"]) for k in
                        ("displacement_volume", "waterplane_area",
                         "vertical_centre_buoyancy", "longitudinal_centre_flotation"))
            print(f"  {name:22s} worst closed-form vs independent = {worst:.2e}")
    return rows


def convergence(verbose: bool = True) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for label, kwargs in (("legacy", {}), ("cwp=0.700/lcb=0.475", {"cwp": 0.700, "lcb": 0.475})):
        for nx, nz in ((41, 13), (81, 21), (161, 41), (321, 81)):
            p = Series60Parameters(length=L, beam=B, draft=T, cb=0.60, nx=nx, nz=nz, **kwargs)
            model = generate_series60(p)
            ref = analytic_series60_hydrostatics(p)
            hydro = model.hydrostatics
            rows.append({
                "hull": label, "nx": nx, "nz": nz,
                "triangles": len(model.mesh.faces),
                "volume_error": hydro.displacement_volume / ref["displacement_volume"] - 1.0,
                "waterplane_error": hydro.waterplane_area / ref["waterplane_area"] - 1.0,
                "lcb_error_over_length": ((hydro.longitudinal_centre_buoyancy
                                           - ref["longitudinal_centre_buoyancy"]) / L),
                "vcb_error": ((hydro.vertical_centre_buoyancy + T) / ref["vertical_centre_buoyancy"] - 1.0),
                "closed_manifold": not model.validation.errors,
            })
            if verbose:
                print(f"  {label:22s} {nx:4d}x{nz:3d}: dV={rows[-1]['volume_error']:+.3%} "
                      f"dLCB={rows[-1]['lcb_error_over_length']:+.2e} L "
                      f"dVCB={rows[-1]['vcb_error']:+.3%}")
    return rows


# ---------------------------------------------------------------------------- figures
def make_figures(out_dir: Path, legacy, corrected, conv) -> list[dict[str, str]]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures: list[dict[str, str]] = []
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1: Cwp against the two baselines
    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=160)
    cbs = [r["cb"] for r in legacy]
    ax.plot(cbs, [r["cwp_model"] for r in legacy], "o-", label="TensorCAD legacy (Cwp = Cp)")
    ax.plot(cbs, [r["cwp_published"] for r in legacy], "s--", label="tabulated Series 60 parents")
    ax.plot(cbs, [r["cwp_schneekluth"] for r in legacy], "^:", label="Schneekluth (1+2Cb)/3")
    ax.plot([r["cb"] for r in corrected], [r["cwp_model"] for r in corrected], "x",
            markersize=9, label="with cwp = published")
    ax.set_xlabel("block coefficient Cb [-]")
    ax.set_ylabel("waterplane coefficient Cwp [-]")
    ax.set_title("Series 60 waterplane coefficient: one exponent cannot carry it")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    path = out_dir / "cwp_vs_baselines.png"
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    figures.append({"name": path.name, "path": str(path),
                    "caption": "Cwp implied by the single-exponent loft (Cwp = Cp) against the "
                               "tabulated parents and the Schneekluth regression; crosses are the "
                               "hull with Cwp as an input."})

    # 2: station shapes, legacy vs corrected
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.4), dpi=160, sharey=True)
    p_legacy = Series60Parameters(length=L, beam=B, draft=T, cb=0.60)
    p_fixed = Series60Parameters(length=L, beam=B, draft=T, cb=0.60, cwp=0.700, lcb=0.475)
    for ax, xi in zip(axes, (0.25, 0.45, 0.65)):
        for p, style, label in ((p_legacy, "-", "legacy"), (p_fixed, "--", "cwp/lcb")):
            law = series60_section_law(p)
            f = law.deck_taper(xi)
            q = law.fullness(xi)
            m = q / (1.0 - q)
            zetas = [i / 200.0 for i in range(201)]
            ys = [0.5 * B * f * z ** (1.0 / m) for z in zetas]
            ax.plot(ys, [-T * (1.0 - z) for z in zetas], style, label=label)
        ax.set_title(f"x/L = {xi:.2f}", fontsize=9)
        ax.set_xlabel("half-breadth [m]")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("z [m]")
    axes[0].legend(fontsize=8)
    fig.suptitle("Cb=0.60 sections: the corrected hull keeps the midship section and "
                 "thins the ends to buy waterplane fullness", fontsize=9)
    path = out_dir / "sections_legacy_vs_cwp.png"
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    figures.append({"name": path.name, "path": str(path),
                    "caption": "Station shapes at three stations, legacy hull against the same "
                               "hull with Cwp and LCB as inputs."})

    # 3: convergence
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), dpi=160)
    for label, marker in (("legacy", "o"), ("cwp=0.700/lcb=0.475", "s")):
        rows = [r for r in conv if r["hull"] == label]
        nx = [r["nx"] for r in rows]
        axes[0].loglog(nx, [abs(r["volume_error"]) for r in rows], marker + "-", label=label)
        axes[1].loglog(nx, [max(abs(r["lcb_error_over_length"]), 1e-13) for r in rows],
                       marker + "-", label=label)
    axes[0].set_title("|volume error| vs analytic", fontsize=9)
    axes[1].set_title("|LCB error| / L vs analytic", fontsize=9)
    for ax in axes:
        ax.set_xlabel("nx"); ax.grid(alpha=0.3, which="both"); ax.legend(fontsize=8)
    path = out_dir / "convergence.png"
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    figures.append({"name": path.name, "path": str(path),
                    "caption": "First-order convergence of the lofted mesh to its own closed "
                               "forms, for the legacy hull and for the Cwp/LCB-controlled hull."})

    # 4: plan view of the waterplane for both hulls
    fig, ax = plt.subplots(figsize=(7.6, 3.2), dpi=160)
    for p, style, label in ((p_legacy, "-", "legacy (Cwp = 0.625)"),
                            (p_fixed, "--", "cwp = 0.700, lcb = 0.475")):
        law = series60_section_law(p)
        xs = [i / 400.0 for i in range(401)]
        ax.plot([L * x for x in xs], [0.5 * B * law.deck_taper(x) for x in xs], style, label=label)
        ax.plot([L * x for x in xs], [-0.5 * B * law.deck_taper(x) for x in xs], style)
    ax.axvline(0.5 * L, color="k", linewidth=0.6, alpha=0.5)
    ax.set_xlabel("x [m] from the aft end")
    ax.set_ylabel("half-breadth at DWL [m]")
    ax.set_title("Design waterplane: fuller plan, aft-shifted centroid", fontsize=9)
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    path = out_dir / "waterplane_plan.png"
    fig.tight_layout(); fig.savefig(path); plt.close(fig)
    figures.append({"name": path.name, "path": str(path),
                    "caption": "Design waterplane of both hulls; the corrected one has the full "
                               "plan of a real parent form and its centroid moved aft."})
    return figures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="artifacts/series60_fidelity")
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args()
    out_dir = Path(args.out)
    if not out_dir.is_absolute():
        out_dir = ROOT / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    print("published table self-consistency (Cb - Cp*Cm):",
          {k: f"{v:+.4f}" for k, v in check_table().items()})
    print("\n1. legacy family (single longitudinal exponent)")
    legacy = legacy_audit()
    print("\n2. family with cwp / lcb as inputs")
    corrected = corrected_audit()
    print("\n3. closed forms vs an independent tanh-sinh integration of the same surface")
    cross = cross_check()
    print("\n4. loft convergence to the closed forms")
    conv = convergence()

    payload = {
        "case": {"length_m": L, "beam_m": B, "draft_m": T,
                 "published_table_source": "Todd 1951 Series 60 parent forms, as reproduced in "
                                           "the standard resistance references",
                 "regression": "Schneekluth & Bertram, Cwp = (1 + 2 Cb)/3"},
        "published_table_consistency": check_table(),
        "legacy": legacy,
        "corrected": corrected,
        "cross_check": cross,
        "convergence": conv,
    }
    (out_dir / "audit.json").write_text(json.dumps(payload, indent=2))

    figures = [] if args.no_figures else make_figures(out_dir, legacy, corrected, conv)
    (out_dir / "figures.json").write_text(json.dumps(figures, indent=2))
    print(f"\nwrote {out_dir/'audit.json'} and {len(figures)} figures")

    worst_legacy = max(abs(r["cwp_error_vs_published"]) for r in legacy)
    worst_legacy_reg = max(abs(r["cwp_error_vs_schneekluth"]) for r in legacy)
    print(f"legacy worst Cwp error: {worst_legacy:.2%} (published table), "
          f"{worst_legacy_reg:.2%} (regression)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())