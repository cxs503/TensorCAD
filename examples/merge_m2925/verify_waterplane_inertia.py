"""Cure of PR #39: evidence for the merged waterplane-inertia references.

The conflict in ``tensorcad/hull.py`` was semantic: ``agent/m2925`` reports the two
waterplane second moments (``I_T``, ``I_L``) on the ``Hydrostatics`` record and in
``analytic_hydrostatics``, while ``main`` grew a split (asymmetric) sine-taper
exponent.  A textual union that kept either side's call would have been wrong --
m2925's call pinned ``breadth_exponent`` to the symmetric constant, so the
asymmetric family would have reported a symmetric-exponent inertia with no error
signal at all.  This script measures the union that was actually committed:

* the analytic ``I_T`` / ``I_L`` against **independent quadrature** of the very
  surface the generator samples (not against the closed form itself), on the
  symmetric family and on an asymmetric one;
* the **mesh** record of :func:`tensorcad.hydrostatics.compute_hydrostatics`
  against those references on a four-grid ladder, so the convergence is visible;
* bit-for-bit back-compatibility of the symmetric family: the merged reference has
  to return the same numbers the pre-merge tree did.

Writes ``docs/merge_m2925/results.json`` and ``docs/merge_m2925/inertia_split.png``
and exits non-zero if any claim below breaks.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from tensorcad import HullParameters, generate_hull  # noqa: E402
from tensorcad.hull import (  # noqa: E402
    BREADTH_EXPONENT,
    analytic_hydrostatics,
    taper_integral,
)
from tensorcad.hydrostatics_analytic import (  # noqa: E402
    tanh_sinh_quadrature,
    waterplane_inertia_longitudinal,
    waterplane_inertia_transverse,
)

OUT = Path(__file__).resolve().parents[2] / "docs" / "merge_m2925"
LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8
ASYMMETRIC = {"fore_exponent": 0.5, "aft_exponent": 0.9}
GRIDS = ((41, 13), (81, 25), (161, 49), (321, 97))

failures: list[str] = []


def claim(ok: bool, text: str) -> None:
    print(("PASS " if ok else "FAIL ") + text)
    if not ok:
        failures.append(text)


def half_breadth(xi: float, fore: float, aft: float) -> float:
    """Half breadth of the sampled surface at ``xi = x/L``, per half of the hull."""
    p = aft if xi <= 0.5 else fore
    return 0.5 * BEAM * math.sin(math.pi * xi) ** p


def quad_transverse(fore: float, aft: float) -> float:
    """``I_T = int b**3/12 dx``, integrated with tanh-sinh over the whole length."""

    def integrand(xi: float) -> float:
        return (2.0 * half_breadth(xi, fore, aft)) ** 3 / 12.0

    return LENGTH * tanh_sinh_quadrature(integrand, 0.0, 1.0, levels=8)


def quad_longitudinal(fore: float, aft: float, centre: float) -> float:
    """``I_L = int (x - LCF)**2 * 2 * half_breadth dx`` about the given station."""

    def integrand(xi: float) -> float:
        x = xi * LENGTH
        return (x - centre) ** 2 * 2.0 * half_breadth(xi, fore, aft)

    return LENGTH * tanh_sinh_quadrature(integrand, 0.0, 1.0, levels=9)


def reference(fore: float, aft: float):
    hydro = analytic_hydrostatics(
        LENGTH, BEAM, DRAFT, fore_exponent=fore, aft_exponent=aft
    )
    return hydro


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    results: dict[str, object] = {
        "case": {
            "length_m": LENGTH,
            "beam_m": BEAM,
            "draft_m": DRAFT,
            "breadth_exponent_symmetric": BREADTH_EXPONENT,
            "asymmetric": ASYMMETRIC,
            "grids": [list(g) for g in GRIDS],
        }
    }

    # ---------------------------------------------------------------- symmetric
    symmetric = reference(BREADTH_EXPONENT, BREADTH_EXPONENT)
    closed_it = waterplane_inertia_transverse(LENGTH, BEAM, breadth_exponent=BREADTH_EXPONENT)
    claim(
        symmetric.transverse_waterplane_inertia == closed_it,
        "symmetric I_T bit-for-bit equal to the pre-merge reference call",
    )
    # ``hull.taper_integral`` is the gamma identity too, but the two paths differ by
    # one ulp on ``3*0.7`` (a published-digit short circuit on the shipped exponent);
    # the record is bit-for-bit the merged reference, which is what back-compatibility
    # means, so the written-out form is held to the last few bits rather than to ``==``.
    closed_form = BEAM**3 * LENGTH / 12.0 * taper_integral(3.0 * BREADTH_EXPONENT)
    claim(
        abs(closed_it / closed_form - 1.0) < 1e-15,
        "pre-merge reference agrees with the written-down closed form to 1e-15",
    )
    it_quad = quad_transverse(BREADTH_EXPONENT, BREADTH_EXPONENT)
    il_quad = quad_longitudinal(
        BREADTH_EXPONENT, BREADTH_EXPONENT, symmetric.longitudinal_centre_flotation
    )
    claim(
        abs(symmetric.transverse_waterplane_inertia / it_quad - 1.0) < 1e-12,
        "symmetric I_T vs independent quadrature within 1e-12 relative",
    )
    claim(
        abs(symmetric.longitudinal_waterplane_inertia / il_quad - 1.0) < 1e-10,
        "symmetric I_L vs independent quadrature within 1e-10 relative",
    )
    results["symmetric"] = {
        "longitudinal_centre_flotation_m": symmetric.longitudinal_centre_flotation,
        "transverse_waterplane_inertia_m4": symmetric.transverse_waterplane_inertia,
        "longitudinal_waterplane_inertia_m4": symmetric.longitudinal_waterplane_inertia,
        "transverse_quadrature_m4": it_quad,
        "longitudinal_quadrature_m4": il_quad,
        "transverse_relative_difference": symmetric.transverse_waterplane_inertia / it_quad - 1.0,
        "longitudinal_relative_difference": (
            symmetric.longitudinal_waterplane_inertia / il_quad - 1.0
        ),
        "basis": symmetric.basis,
    }

    # --------------------------------------------------------------- asymmetric
    asymmetric = reference(ASYMMETRIC["fore_exponent"], ASYMMETRIC["aft_exponent"])
    lcf = asymmetric.longitudinal_centre_flotation
    it_quad = quad_transverse(ASYMMETRIC["fore_exponent"], ASYMMETRIC["aft_exponent"])
    il_quad = quad_longitudinal(ASYMMETRIC["fore_exponent"], ASYMMETRIC["aft_exponent"], lcf)
    claim(
        abs(asymmetric.transverse_waterplane_inertia / it_quad - 1.0) < 1e-12,
        "asymmetric I_T vs independent quadrature within 1e-12 relative",
    )
    claim(
        abs(asymmetric.longitudinal_waterplane_inertia / il_quad - 1.0) < 1e-10,
        "asymmetric I_L (about the LCF) vs independent quadrature within 1e-10 relative",
    )
    claim(
        lcf > 0.5 * LENGTH,
        "asymmetric LCF forward of midship for fore exponent < aft exponent",
    )
    # A symmetric exponent would sit 2 % low here -- the error the union avoids.
    symmetric_exponent_it = waterplane_inertia_transverse(LENGTH, BEAM)
    claim(
        abs(symmetric_exponent_it / it_quad - 1.0) > 1e-3,
        "keeping the symmetric exponent on an asymmetric hull would be off by > 1e-3",
    )
    results["asymmetric"] = {
        "fore_exponent": ASYMMETRIC["fore_exponent"],
        "aft_exponent": ASYMMETRIC["aft_exponent"],
        "longitudinal_centre_flotation_m": lcf,
        "transverse_waterplane_inertia_m4": asymmetric.transverse_waterplane_inertia,
        "longitudinal_waterplane_inertia_m4": asymmetric.longitudinal_waterplane_inertia,
        "transverse_quadrature_m4": it_quad,
        "longitudinal_quadrature_m4": il_quad,
        "transverse_relative_difference": asymmetric.transverse_waterplane_inertia / it_quad - 1.0,
        "longitudinal_relative_difference": (
            asymmetric.longitudinal_waterplane_inertia / il_quad - 1.0
        ),
        "symmetric_exponent_transverse_relative_error": symmetric_exponent_it / it_quad - 1.0,
        "longitudinal_inertia_quadrature_evidence": (
            waterplane_inertia_longitudinal(
                LENGTH,
                BEAM,
                fore_exponent=ASYMMETRIC["fore_exponent"],
                aft_exponent=ASYMMETRIC["aft_exponent"],
                longitudinal_centre=lcf,
            )[1]
        ),
        "basis": asymmetric.basis,
    }

    # ------------------------------------------------------- mesh ladder (both)
    ladder: dict[str, list[dict[str, float]]] = {}
    for label, exponents, ref in (
        ("symmetric", (BREADTH_EXPONENT, BREADTH_EXPONENT), symmetric),
        (
            "asymmetric",
            (ASYMMETRIC["fore_exponent"], ASYMMETRIC["aft_exponent"]),
            asymmetric,
        ),
    ):
        rows = []
        for nx, nz in GRIDS:
            hydro = generate_hull(
                HullParameters(
                    length=LENGTH,
                    beam=BEAM,
                    draft=DRAFT,
                    nx=nx,
                    nz=nz,
                    fore_exponent=exponents[0],
                    aft_exponent=exponents[1],
                )
            ).hydrostatics
            rows.append(
                {
                    "nx": nx,
                    "nz": nz,
                    "transverse_waterplane_inertia_m4": hydro.transverse_waterplane_inertia,
                    "longitudinal_waterplane_inertia_m4": hydro.longitudinal_waterplane_inertia,
                    "transverse_relative_error": (
                        hydro.transverse_waterplane_inertia
                        / ref.transverse_waterplane_inertia
                        - 1.0
                    ),
                    "longitudinal_relative_error": (
                        hydro.longitudinal_waterplane_inertia
                        / ref.longitudinal_waterplane_inertia
                        - 1.0
                    ),
                    "block_coefficient": hydro.block_coefficient,
                }
            )
        ladder[label] = rows

    for label, rows in ladder.items():
        magnitudes = [abs(r["transverse_relative_error"]) for r in rows]
        claim(
            magnitudes[-1] < magnitudes[0],
            f"{label}: transverse mesh error shrinks with the grid",
        )
        claim(
            magnitudes[-1] < 5e-3,
            f"{label}: transverse mesh error below 0.5 % on the finest grid",
        )
        longitudinal = [abs(r["longitudinal_relative_error"]) for r in rows]
        claim(
            longitudinal[-1] < longitudinal[0],
            f"{label}: longitudinal mesh error shrinks with the grid",
        )
    results["mesh_ladder"] = ladder

    # ------------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), dpi=160)
    for axis, quantity, title in (
        (axes[0], "transverse", r"$I_T=\int y^2\,dA$"),
        (axes[1], "longitudinal", r"$I_L=\int (x-LCF)^2\,dA$"),
    ):
        for label, style in (("symmetric", "o-"), ("asymmetric", "s-.")):
            rows = ladder[label]
            h = [LENGTH / r["nx"] for r in rows]
            err = [abs(r[f"{quantity}_relative_error"]) for r in rows]
            axis.loglog(h, err, style, label=label)
        axis.set_xlabel(r"station spacing $L/n_x$  [m]")
        axis.set_ylabel("relative error vs closed form")
        axis.set_title(title)
        axis.grid(True, which="both", alpha=0.3)
        axis.legend()
    fig.suptitle(
        "Merged waterplane-inertia references vs the hull mesh (length 10 m, beam 2 m)"
    )
    fig.tight_layout()
    fig.savefig(OUT / "inertia_split.png")
    plt.close(fig)

    results["claims"] = {"failed": failures, "passed": not failures}
    (OUT / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    print(f"wrote {OUT / 'results.json'} and {OUT / 'inertia_split.png'}")
    if failures:
        print(f"{len(failures)} claim(s) failed")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())