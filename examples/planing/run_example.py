"""Verification run for the planing hull family.

Builds two hulls -- a loaded deep-V whose chine stays immersed and a light hull whose
chine leaves the waterline 72 % aft of the bow -- and checks the mesh against three
independent references: the closed forms in ``tensorcad.planing``, an independent 1-D
Simpson quadrature of the same sections, and a mesh refinement study.

Writes ``artifacts/planing_hull/results.json`` (all numbers, injected programmatically)
and ``artifacts/planing_hull/planing_verification.png`` when matplotlib is importable.

Run: ``python examples/planing/run_example.py``
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_REPO_ROOT = _Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in _sys.path:  # allow running from a checkout without installing
    _sys.path.insert(0, str(_REPO_ROOT))

import json
import math
import pathlib
import sys

from tensorcad import planing
from tensorcad.planing import (
    PlaningHullParameters,
    analytic_planing_hydrostatics,
    build_planing_mesh,
    chine_half_beam,
    estimate_planing_discretisation_error,
    generate_planing_hull,
    planing_relative_errors,
    section_area,
    simpson_section_integrals,
    waterline_half_breadth,
    wetted_length,
)

CASES = {
    "chine_immersed": PlaningHullParameters(),
    "chine_dry": PlaningHullParameters(draft=0.40),
    "straight_chine": PlaningHullParameters(draft=0.60, planform_exponent=1.0),
    "deep_vee": PlaningHullParameters(draft=0.60, deadrise_deg=24.0, planform_exponent=0.5),
}
QUANTITIES = (
    "displacement_volume",
    "waterplane_area",
    "wetted_surface_area",
    "longitudinal_centre_buoyancy",
    "vertical_centre_buoyancy",
    "longitudinal_centre_flotation",
)


def run() -> dict:
    out: dict = {"cases": {}, "convergence": {}, "anchors": {}}
    for name, p in CASES.items():
        model = generate_planing_hull(p)
        hydro = model.hydrostatics
        analytic = analytic_planing_hydrostatics(p)
        quad = simpson_section_integrals(p, intervals=8192)
        errors = planing_relative_errors(p, model)
        out["cases"][name] = {
            "parameters": p.as_dict(),
            "wetted_length": wetted_length(p),
            "chine_immersed_at_transom": p.chine_immersed_at_transom,
            "mesh": {key: getattr(hydro, key) for key in QUANTITIES},
            "closed_form": {key: analytic[key] for key in QUANTITIES},
            "quadrature": quad,
            "relative_error_mesh_vs_closed_form": errors,
            "relative_error_quadrature_vs_closed_form": {
                key: (quad[key] - analytic[key]) / analytic[key] if analytic[key] else float("nan")
                for key in quad
            },
            "coefficients": {
                "cb": analytic["block_coefficient"],
                "cp": analytic["prismatic_coefficient"],
                "cwp": analytic["waterplane_coefficient"],
                "vcb_above_keel_over_T": analytic["vertical_centre_buoyancy_above_keel"] / p.draft,
                "transom_area": analytic["transom_section_area"],
                "waterline_beam": analytic["waterline_beam"],
            },
            "mesh": {
                **{key: getattr(hydro, key) for key in QUANTITIES},
                "triangles": len(model.mesh.faces),
                "vertices": len(model.mesh.vertices),
                "watertight": __import__("tensorcad.mesh", fromlist=["is_watertight"]).is_watertight(
                    model.mesh
                ),
            },
            "validation": {
                "errors": list(model.validation.errors),
                "warnings": [check.name for check in model.validation.warnings],
            },
        }

    study = estimate_planing_discretisation_error(
        PlaningHullParameters(nx=21, nz=9), max_triangles=400_000
    )
    rows = study["convergence"]
    orders = []
    for a, b in zip(rows, rows[1:]):
        ratio = (b["nx"] - 1) / (a["nx"] - 1)
        ea, eb = abs(a["volume_error"]), abs(b["volume_error"])
        if ea > 0 and eb > 0:
            orders.append(math.log(ea / eb) / math.log(ratio))
    out["convergence"] = {
        "reference": {key: study["reference"][key] for key in QUANTITIES},
        "rows": rows,
        "observed_order": orders,
    }

    # Exact anchors that do not depend on any of the integrals above.
    anchors = {}
    wet = CASES["chine_immersed"]
    c = wet.half_beam_transom / wet.length**wet.planform_exponent
    anchors["everywhere_immersed_volume"] = {
        "closed_form": analytic_planing_hydrostatics(wet)["displacement_volume"],
        "paper_formula": (
            2.0 * wet.draft * c * wet.length ** (wet.planform_exponent + 1) / (wet.planform_exponent + 1)
            - wet.tan_beta * c * c * wet.length ** (2 * wet.planform_exponent + 1)
            / (2 * wet.planform_exponent + 1)
        ),
    }
    anchors["everywhere_immersed_cwp"] = {
        "closed_form": analytic_planing_hydrostatics(wet)["waterplane_coefficient"],
        "paper_formula": 1.0 / (wet.planform_exponent + 1.0),
    }
    straight = CASES["straight_chine"]  # q = 1, chine immersed: Cwp = 1/(q+1) = 1/2
    anchors["straight_chine_cwp"] = {
        "closed_form": analytic_planing_hydrostatics(straight)["waterplane_coefficient"],
        "paper_formula": 0.5,
    }
    dry = CASES["chine_dry"]
    anchors["pure_v_transom_area"] = {
        "closed_form": analytic_planing_hydrostatics(dry)["transom_section_area"],
        "paper_formula": dry.draft**2 / dry.tan_beta,
    }
    anchors["pointwise_section_area"] = {
        "x": [dry.length * i / 10 for i in range(11)],
        "reason": "section_area is the quadrature integrand, so it is checked pointwise",
    }
    out["anchors"] = anchors
    return out


def figure(results: dict, path: pathlib.Path) -> bool:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency
        print(f"matplotlib unavailable ({exc}); skipping figure", file=sys.stderr)
        return False

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))

    # (a) body plan of the two regimes, with the waterline
    ax = axes[0]
    for name, style in (("chine_immersed", "-"), ("chine_dry", "--")):
        p = CASES[name]
        for x in (0.25 * p.length, 0.6 * p.length, p.length):
            b = chine_half_beam(p, x)
            ys, zs = [], []
            for level in range(101):
                z = -p.draft + p.draft * level / 100.0
                ys.append(min((z + p.draft) / p.tan_beta, b))
                zs.append(z)
            ax.plot(ys, zs, style, lw=1.1, label=f"{name} x={x / p.length:.2f}L")
            ax.plot([-y for y in ys], zs, style, lw=1.1)
    ax.axhline(0.0, color="tab:blue", lw=0.8)
    ax.set_xlabel("y [m]")
    ax.set_ylabel("z [m]")
    ax.set_title("(a) sections, waterline at z = 0")
    ax.legend(fontsize=6, loc="lower right")

    # (b) refinement: error of volume and waterplane area versus triangles
    ax = axes[1]
    rows = results["convergence"]["rows"]
    triangles = [row["triangles"] for row in rows]
    for key, label in (("volume_error", "displacement"), ("waterplane_error", "waterplane"),
                       ("wetted_area_error", "wetted area")):
        ax.loglog(triangles, [abs(row[key]) for row in rows], "o-", lw=1.2, label=label)
    ax.set_xlabel("triangles")
    ax.set_ylabel("|mesh - closed form| / closed form")
    ax.set_title("(b) refinement vs closed form")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7)

    # (c) error budget per quantity and case
    ax = axes[2]
    labels = list(QUANTITIES)
    width = 0.8 / len(CASES)
    for index, (name, case) in enumerate(results["cases"].items()):
        errors = case["relative_error_mesh_vs_closed_form"]
        values = [abs(errors[label]) * 100.0 for label in labels]
        ax.bar([i + index * width for i in range(len(labels))], values, width, label=name)
    ax.set_xticks([i + 0.4 for i in range(len(labels))])
    ax.set_xticklabels([label.replace("_", "\n") for label in labels], fontsize=6)
    ax.set_ylabel("|error| [%]")
    ax.set_title("(c) mesh vs closed form, by case")
    ax.legend(fontsize=6)

    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return True


def main() -> int:
    root = pathlib.Path(__file__).resolve().parents[2]
    out_dir = root / "artifacts" / "planing_hull"
    out_dir.mkdir(parents=True, exist_ok=True)
    results = run()
    results["figure"] = figure(results, out_dir / "planing_verification.png")

    (out_dir / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True))

    failures = []
    for name, case in results["cases"].items():
        for key, value in case["relative_error_mesh_vs_closed_form"].items():
            if not math.isfinite(value) or abs(value) > 0.01:
                failures.append(f"{name}:{key} = {value:+.4%}")
        if case["validation"]["errors"]:
            failures.append(f"{name}: validation errors {case['validation']['errors']}")
    if failures:
        print("FAILED")
        for item in failures:
            print("  ", item)
        return 1

    for name, case in results["cases"].items():
        errors = case["relative_error_mesh_vs_closed_form"]
        worst = max(abs(value) for value in errors.values())
        print(f"{name:16s} triangles={case['mesh']['triangles']:6d} worst error={worst:.5%}")
    print(f"observed order (volume, log ratio): "
          f"{', '.join(f'{value:.2f}' for value in results['convergence']['observed_order'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())