#!/usr/bin/env python3
"""Make the longitudinal centres measurable instead of an identity.

The shipped hull family is mirror-symmetric, so ``LCB = LCF = length/2`` is a
consequence of the symmetry and cannot detect an antisymmetric bias in the moment
integration. This example gives each half its own sine-taper exponent, computes
the analytic reference of the resulting asymmetric surface, and compares the
sampled mesh against it: volume and waterplane area against closed forms, the
longitudinal centres against a quadrature half-moment that is cross-checked by a
second, unrelated scheme.

Writes ``results.json``, ``report.md`` and ``plot.png`` and exits 1 if a claim
stops holding, so it is usable as a CI step.

    python3 examples/asymmetric_hull/run_example.py --output docs/asymmetric_hull
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_REPO_ROOT = _Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in _sys.path:  # allow running from a checkout without installing
    _sys.path.insert(0, str(_REPO_ROOT))

import argparse
import itertools
import json
from pathlib import Path

from tensorcad import HullParameters, analytic_hydrostatics, generate_hull
from tensorcad.hull import (
    TAPER_INTEGRAL,
    _half_first_moment,
    _half_first_moment_exp_sinh,
    analytic_longitudinal_centre,
    analytic_volume,
    analytic_waterplane_area,
)

L, BEAM, DRAFT = 10.0, 2.0, 0.8
DENSITY = 1025.0
GRIDS = ((41, 13), (81, 25), (161, 49), (321, 97))
ASYM = ((0.5, 0.9), (0.4, 1.4))
EXPLORED = (0.3, 0.4, 0.5, 0.7, 1.0, 1.4, 2.0)


def measure(fore: float, aft: float, nx: int, nz: int) -> dict[str, float]:
    parameters = HullParameters(
        L, BEAM, DRAFT, nx, nz, DENSITY, fore_exponent=fore, aft_exponent=aft
    )
    model = generate_hull(parameters)
    measured = model.hydrostatics
    reference = analytic_hydrostatics(
        L, BEAM, DRAFT, DENSITY, fore_exponent=fore, aft_exponent=aft
    )
    return {
        "nx": nx,
        "nz": nz,
        "volume_error_percent": 100.0 * (measured.displacement_volume / reference.displacement_volume - 1.0),
        "waterplane_error_percent": 100.0 * (measured.waterplane_area / reference.waterplane_area - 1.0),
        "wetted_area_error_percent": 100.0 * (measured.wetted_surface_area / reference.wetted_surface_area - 1.0),
        "lcb_mesh_m": measured.longitudinal_centre_buoyancy,
        "lcb_reference_m": reference.longitudinal_centre_buoyancy,
        "lcb_error_mm": 1000.0 * (measured.longitudinal_centre_buoyancy - reference.longitudinal_centre_buoyancy),
        "lcf_error_mm": 1000.0 * (measured.longitudinal_centre_flotation - reference.longitudinal_centre_flotation),
        "block_coefficient_error": measured.block_coefficient - reference.block_coefficient,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(Path(__file__).resolve().parents[2] / "docs" / "asymmetric_hull"))
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    # --- the symmetric family must not have moved -------------------------
    symmetric_reference = {
        "volume_m3": analytic_volume(L, BEAM, DRAFT),
        "volume_from_constant_m3": BEAM * L * DRAFT * TAPER_INTEGRAL * (2.0 / 3.0),
        "waterplane_m2": analytic_waterplane_area(L, BEAM),
        "longitudinal_centre_m": analytic_longitudinal_centre(L),
        "taper_integral": TAPER_INTEGRAL,
    }

    # --- the one quantity without a closed form ---------------------------
    moment_check = []
    for p in (0.2, 0.4, 0.7, 1.0, 1.4, 2.0, 3.0):
        gauss = _half_first_moment(p)
        exp_sinh = _half_first_moment_exp_sinh(p)
        moment_check.append(
            {
                "exponent": p,
                "gauss_legendre": gauss,
                "exp_sinh": exp_sinh,
                "relative_difference": abs(gauss / exp_sinh - 1.0),
            }
        )

    # --- reference behaviour in the asymmetry parameter --------------------
    sweep = []
    for fore in EXPLORED:
        reference = analytic_longitudinal_centre(L, fore_exponent=fore, aft_exponent=0.9)
        mesh = measure(fore, 0.9, 81, 25)
        sweep.append(
            {
                "fore_exponent": fore,
                "aft_exponent": 0.9,
                "lcb_reference_m": reference,
                "lcb_mesh_m": mesh["lcb_mesh_m"],
                "lcb_error_mm": mesh["lcb_error_mm"],
                "volume_m3": analytic_volume(L, BEAM, DRAFT, fore_exponent=fore, aft_exponent=0.9),
            }
        )

    # --- mesh against reference on a refinement ladder --------------------
    convergence = []
    for fore, aft in ASYM:
        for nx, nz in GRIDS:
            row = measure(fore, aft, nx, nz)
            row["fore_exponent"] = fore
            row["aft_exponent"] = aft
            row["h_over_length"] = 1.0 / (nx - 1)
            convergence.append(row)

    # --- involution: swapping the halves mirrors the centre about midship --
    involution = []
    for fore, aft in ((0.5, 0.9), (0.4, 1.4), (0.7, 1.2)):
        one = analytic_longitudinal_centre(L, fore_exponent=fore, aft_exponent=aft)
        two = analytic_longitudinal_centre(L, fore_exponent=aft, aft_exponent=fore)
        mesh_one = measure(fore, aft, 81, 25)["lcb_mesh_m"]
        mesh_two = measure(aft, fore, 81, 25)["lcb_mesh_m"]
        involution.append(
            {
                "fore": fore,
                "aft": aft,
                "reference_pair_sum_m": one + two,
                "reference_residual_m": one + two - L,
                "mesh_pair_sum_m": mesh_one + mesh_two,
                "mesh_residual_m": mesh_one + mesh_two - L,
            }
        )

    claims = {
        "symmetric_volume_uses_the_published_constant": (
            symmetric_reference["volume_m3"] == symmetric_reference["volume_from_constant_m3"]
        ),
        "symmetric_centre_is_exactly_midship": symmetric_reference["longitudinal_centre_m"] == 0.5 * L,
        "half_moment_two_schemes_agree_below_1e-13": all(
            row["relative_difference"] < 1e-13 for row in moment_check
        ),
        "involution_is_exact_in_the_reference": all(
            row["reference_residual_m"] < 1e-12 for row in involution
        ),
        "involution_holds_on_the_mesh_below_1um": all(
            row["mesh_residual_m"] < 1e-6 for row in involution
        ),
        "the_centre_leaves_midship_as_soon_as_the_halves_differ": all(
            abs(row["lcb_reference_m"] - 0.5 * L) > 0.05 for row in sweep
        ),
        "a_two_tenths_asymmetry_moves_the_centre_over_0.1_m": all(
            abs(row["lcb_reference_m"] - 0.5 * L) > 0.1
            for row in sweep
            if abs(row["fore_exponent"] - row["aft_exponent"]) >= 0.2
        ),
        "a_four_tenths_asymmetry_moves_the_centre_over_0.2_m": all(
            abs(row["lcb_reference_m"] - 0.5 * L) > 0.2
            for row in sweep
            if abs(row["fore_exponent"] - row["aft_exponent"]) >= 0.4
        ),
        "the_reference_centre_is_monotone_in_the_fore_exponent": all(
            earlier["lcb_reference_m"] > later["lcb_reference_m"]
            for earlier, later in itertools.pairwise(sweep)
        ),
        "volume_error_below_1_percent_at_81x25": all(
            abs(row["volume_error_percent"]) < 1.0
            for row in convergence
            if row["nx"] == 81
        ),
        "centre_error_below_1_percent_of_length_at_161x49": all(
            abs(row["lcb_error_mm"]) < 10.0
            for row in convergence
            if row["nx"] == 161
        ),
        "centre_error_is_monotone_in_the_grid": all(
            abs(row["lcb_error_mm"]) > abs(next_row["lcb_error_mm"])
            for fore in (a for a, _ in ASYM)
            for row, next_row in zip(
                [r for r in convergence if r["fore_exponent"] == fore],
                [r for r in convergence if r["fore_exponent"] == fore][1:],
                strict=False,
            )
        ),
        "lcb_equals_lcf_on_the_mesh": all(abs(row["lcb_error_mm"] - row["lcf_error_mm"]) < 1e-6 for row in convergence),
    }
    failed = [name for name, ok in claims.items() if not ok]

    results = {
        "case": {
            "length_m": L,
            "beam_m": BEAM,
            "draft_m": DRAFT,
            "water_density_kg_per_m3": DENSITY,
            "family": "half-breadth (B/2)*sin(pi*x/L)^p*sqrt((z+T)/T), p = aft/fore per half",
            "grids": [{"nx": nx, "nz": nz} for nx, nz in GRIDS],
        },
        "symmetric_reference": symmetric_reference,
        "half_moment": moment_check,
        "asymmetry_sweep": sweep,
        "convergence": convergence,
        "involution": involution,
        "claims": claims,
        "failed_claims": failed,
    }
    (output / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True), encoding="utf-8")

    plot(output, sweep, convergence)
    (output / "report.md").write_text(render(results), encoding="utf-8")

    print(json.dumps({"claims": len(claims), "failed": failed, "output": str(output)}, indent=2))
    return 1 if failed else 0


def plot(output: Path, sweep: list[dict], convergence: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, (left, right) = plt.subplots(1, 2, figsize=(9.0, 3.6))
    fore = [row["fore_exponent"] for row in sweep]
    left.plot(fore, [row["lcb_reference_m"] for row in sweep], "k-", label="analytic reference")
    left.plot(fore, [row["lcb_mesh_m"] for row in sweep], "o", mfc="none", label="mesh 81x25")
    left.axhline(0.5 * L, color="0.6", ls=":", label="midship L/2")
    left.set_xlabel("fore exponent  $p_{fore}$   (aft exponent fixed at 0.9)")
    left.set_ylabel("LCB = LCF  [m]")
    left.set_title("asymmetry moves the longitudinal centre")
    left.legend()
    left.grid(alpha=0.3)

    for fore_value, aft_value in ASYM:
        rows = [row for row in convergence if row["fore_exponent"] == fore_value]
        right.loglog(
            [row["h_over_length"] for row in rows],
            [abs(row["lcb_error_mm"]) for row in rows],
            "o-",
            label=f"$p_f$={fore_value}, $p_a$={aft_value}",
        )
    right.set_xlabel("station spacing  $1/(n_x-1)$")
    right.set_ylabel("|mesh - reference| LCB  [mm]")
    right.set_title("the centre converges onto the reference")
    right.grid(alpha=0.3, which="both")
    right.legend()
    figure.tight_layout()
    figure.savefig(output / "plot.png", dpi=95)
    plt.close(figure)


def render(results: dict) -> str:
    symmetric = results["symmetric_reference"]
    lines = [
        "# Asymmetric fore/aft taper: the longitudinal centres become measurable",
        "",
        "Generated by `examples/asymmetric_hull/run_example.py`; every number below is",
        "read from `results.json`, which the same run wrote.",
        "",
        "Family: half-breadth `(B/2)*sin(pi*x/L)^p*sqrt((z+T)/T)` with `p = aft_exponent`",
        f"abaft midship and `p = fore_exponent` forward of it; L={L} m, B={BEAM} m, T={DRAFT} m,",
        f"seawater {DENSITY} kg/m^3.",
        "",
        "## 1. The symmetric form did not move",
        "",
        "| quantity | value |",
        "|---|---|",
        f"| analytic volume [m^3] | {symmetric['volume_m3']!r} |",
        f"| `B*L*T*TAPER_INTEGRAL*2/3` [m^3] | {symmetric['volume_from_constant_m3']!r} |",
        f"| analytic waterplane [m^2] | {symmetric['waterplane_m2']!r} |",
        f"| analytic LCB/LCF [m] | {symmetric['longitudinal_centre_m']!r} |",
        "",
        "## 2. The half moment, by two unrelated schemes",
        "",
        "| p | Gauss-Legendre (x = t^3/2) | exp-sinh (x = e^-v/2) | relative difference |",
        "|---:|---:|---:|---:|",
    ]
    for row in results["half_moment"]:
        lines.append(
            f"| {row['exponent']} | {row['gauss_legendre']:.17g} | {row['exp_sinh']:.17g} | "
            f"{row['relative_difference']:.2e} |"
        )
    lines += [
        "",
        "## 3. Reference and mesh against the fore exponent (aft fixed at 0.9)",
        "",
        "| p_fore | reference LCB [m] | mesh 81x25 LCB [m] | error [mm] | volume [m^3] |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in results["asymmetry_sweep"]:
        lines.append(
            f"| {row['fore_exponent']:.2f} | {row['lcb_reference_m']:.6f} | {row['lcb_mesh_m']:.6f} | "
            f"{row['lcb_error_mm']:+.4f} | {row['volume_m3']:.6f} |"
        )
    lines += [
        "",
        "## 4. Refinement ladder",
        "",
        "| p_fore | p_aft | grid | volume err [%] | waterplane err [%] | wetted err [%] | LCB err [mm] |",
        "|---:|---:|---|---:|---:|---:|---:|",
    ]
    for row in results["convergence"]:
        lines.append(
            f"| {row['fore_exponent']:.2f} | {row['aft_exponent']:.2f} | {row['nx']}x{row['nz']} | "
            f"{row['volume_error_percent']:+.4f} | {row['waterplane_error_percent']:+.4f} | "
            f"{row['wetted_area_error_percent']:+.4f} | {row['lcb_error_mm']:+.4f} |"
        )
    lines += [
        "",
        "## 5. Swapping the halves mirrors the centre about midship",
        "",
        "| p_fore | p_aft | reference sum [m] | reference residual [m] | mesh residual [m] |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in results["involution"]:
        lines.append(
            f"| {row['fore']} | {row['aft']} | {row['reference_pair_sum_m']:.15g} | "
            f"{row['reference_residual_m']:+.2e} | {row['mesh_residual_m']:+.2e} |"
        )
    lines += ["", "## Claims", ""]
    for name, ok in results["claims"].items():
        lines.append(f"- {'PASS' if ok else 'FAIL'} — `{name}`")
    lines += [
        "",
        "## Reading the numbers",
        "",
        "* The asymmetric reference is closed form for volume, waterplane area, maximum",
        "  section area, VCB and both coefficients. Only the x-weighted half-moment has no",
        "  elementary form, and it is computed two independent ways that agree to machine",
        "  precision, so the reference is not the weak link in any comparison here.",
        "* The mesh is inscribed between stations, so volume, waterplane area, wetted area",
        "  and the centre all come out systematically on one side of the reference and",
        "  converge to it; a sign flip would be a defect rather than a refinement level.",
        "* The mesh LCB error is antisymmetric under swapping the exponents while the mesh",
        "  volume is unchanged: the discretisation error, not the geometry, is what breaks",
        "  the involution on the mesh.",
        "",
        "`plot.png` (left) shows the centre against the fore exponent, reference and mesh;",
        "the dotted line is midship, which the symmetric family never leaves.",
        "",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    raise SystemExit(main())
