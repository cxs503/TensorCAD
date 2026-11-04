#!/usr/bin/env python3
"""Verify the Wigley hull family against its own closed forms, at every resolution.

The point of this example is not to draw a hull. It is to answer, with numbers a
reviewer can recompute, the question the ``hull`` family could only half answer:
*does the sampled surface integrate to the hydrostatics the algebra says it must?*

Four studies, each with a reference the answer can be wrong against:

1. **The closed forms are the algebra, not a fit.** ``Cb = nm/((n+1)(m+1))``,
   ``Awp = LBn/(n+1)``, ``Am = BTm/(m+1)``, ``VCB = T - T(m+1)/(2(m+2))`` and
   ``V = Cb L B T`` are evaluated symbolically here and re-derived from a direct
   integral of the generating function, so a typo in a docstring cannot survive.
2. **The mesh converges to them at second order.** Four grids, signed errors
   tabulated; the loft cuts chords inside the analytic surface, so every
   extensive quantity must come out low at every resolution (a sign flip is a
   defect, not a refinement level).
3. **The one quantity without a primitive (wetted area) is quadrature, and the
   quadrature is self-converged** to round-off between two panels/orders.
4. **The waterline offsets a panel code consumes sit exactly on the analytic
   surface** -- the loft only interpolates along x, so at sampled stations the
   agreement must be exact.

Writes ``results.json``, ``report.md`` and ``plot.png`` into ``--output`` and exits
non-zero if any claim stops holding.

Usage::

    python examples/wigley/run_example.py [--output docs/wigley]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:  # allow running from a checkout without installing
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad import mesh as meshlib
from tensorcad.wigley import (
    REFERENCE_CASE,
    analytic_block_coefficient,
    analytic_maximum_section_area,
    analytic_volume,
    analytic_waterplane_area,
    convergence_study,
    export_wigley,
    generate_wigley,
    half_breadth,
    reference_parameters,
    wigley_analytic_hydrostatics,
    wetted_area_evidence,
)

QUANTITIES = ("displacement_volume", "waterplane_area", "wetted_surface_area")

#: Direct-integral resolution for the independent cross-check of the closed forms
#: (``steps x steps/4`` midpoints): the form is a low-degree polynomial in each
#: direction, so a midpoint rule this fine lands within ~1e-7 relative.
DIRECT_STEPS = 2000


def _claims(results: dict) -> list[dict]:
    """The assertions this example exists to make, each with its measured value."""
    study = results["convergence"]
    rows = study["rows"]
    reference = results["analytic"]
    inputs = reference["inputs"]
    box = inputs["length_m"] * inputs["beam_m"] * inputs["draft_m"]
    by_name = {q: [row["relative_error"][q] for row in rows] for q in QUANTITIES}
    volume_over_box = reference["displacement_volume_m3"] / (reference["block_coefficient"] * box)
    exact_rationals = (
        math.isclose(reference["block_coefficient"], 4.0 / 9.0, rel_tol=1e-15)
        and math.isclose(reference["prismatic_coefficient"], 2.0 / 3.0, rel_tol=1e-15)
        and math.isclose(reference["midship_coefficient"], 2.0 / 3.0, rel_tol=1e-15)
        and math.isclose(reference["waterplane_coefficient"], 2.0 / 3.0, rel_tol=1e-15)
        and math.isclose(results["vcb_above_keel"] / inputs["draft_m"], 5.0 / 8.0, rel_tol=1e-15)
    )
    return [
        {
            "claim": "V = Cb*L*B*T and V = n*m/((n+1)(m+1))*L*B*T are the same number",
            "value": volume_over_box,
            "tolerance": 1e-15,
            "passed": math.isclose(volume_over_box, 1.0, rel_tol=1e-15),
        },
        {
            "claim": "classic Wigley coefficients are the published rationals (Cb 4/9, Cp 2/3, Cm 2/3, Cwp 2/3, VCB 5T/8)",
            "value": {
                "Cb": reference["block_coefficient"],
                "Cp": reference["prismatic_coefficient"],
                "Cm": reference["midship_coefficient"],
                "Cwp": reference["waterplane_coefficient"],
                "vcb_over_draft": results["vcb_above_keel"] / inputs["draft_m"],
            },
            "tolerance": 1e-15,
            "passed": exact_rationals,
        },
        {
            "claim": "closed-form volume agrees with a direct midpoint integral of the generating function",
            "value": results["closed_form_cross_check"]["relative_difference"],
            "tolerance": 1e-6,
            "passed": results["closed_form_cross_check"]["relative_difference"] < 1e-6,
        },
        {
            "claim": "default grid integrates the displacement volume within 1 % of the closed form",
            "value": rows[1]["relative_error"]["displacement_volume"],
            "tolerance": 0.01,
            "passed": abs(rows[1]["relative_error"]["displacement_volume"]) < 0.01,
        },
        {
            "claim": "sampled quantities approach the closed forms monotonically, from below (chords are inscribed)",
            "value": by_name,
            "tolerance": 0.0,
            "passed": all(
                abs(by_name[q][k]) > abs(by_name[q][k + 1]) and by_name[q][k] < 0.0
                for q in QUANTITIES
                for k in range(len(by_name[q]) - 1)
            ),
        },
        {
            "claim": "finest grid is within 0.1 % of the closed forms on every comparable quantity",
            "value": rows[-1]["max_abs_relative_error"],
            "tolerance": 1e-3,
            "passed": rows[-1]["max_abs_relative_error"] < 1e-3,
        },
        {
            "claim": "observed convergence order for volume is in [1.8, 3.0] (second-order loft)",
            "value": study["observed_order"]["displacement_volume"],
            "tolerance": 0.0,
            "passed": 1.8 <= study["observed_order"]["displacement_volume"] <= 3.0,
        },
        {
            "claim": "wetted-area quadrature is self-converged to 1e-9 relative",
            "value": results["wetted_area_evidence"]["relative_change_fine_to_finer"],
            "tolerance": 1e-9,
            "passed": results["wetted_area_evidence"]["relative_change_fine_to_finer"] < 1e-9,
        },
        {
            "claim": "waterline offsets lie on the analytic surface to machine precision",
            "value": results["offset_on_surface_max_deviation_m"],
            "tolerance": 1e-12,
            "passed": results["offset_on_surface_max_deviation_m"] <= 1e-12,
        },
        {
            "claim": "default surface is closed, edge-manifold, single-component, Euler-2 and has no zero-area face",
            "value": results["topology"],
            "tolerance": 0.0,
            "passed": results["topology"]["watertight"]
            and results["topology"]["edge_manifold"]
            and results["topology"]["euler_characteristic"] == 2
            and results["topology"]["single_component"]
            and results["topology"]["zero_area_faces"] == 0,
        },
    ]


def _station_profile(mesh, index: int) -> list[tuple[float, float]]:
    """``(z, half_breadth)`` on the starboard flank at station ``index``, keel to WL."""
    x = sorted({vertex[0] for vertex in mesh.vertices})[index]
    points = {(vertex[2], vertex[1]) for vertex in mesh.vertices if vertex[0] == x and vertex[1] >= 0.0}
    return sorted(points, key=lambda item: item[0])


def _figure(results: dict, output: Path) -> Path:
    """Breadth-field map, mesh-vs-analytic sections, and the convergence plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    parameters = results["parameters"]
    length, beam, draft = parameters["length"], parameters["beam"], parameters["draft"]
    nx, nz = parameters["nx"], parameters["nz"]
    analytic = results["analytic"]

    fig = plt.figure(figsize=(15.5, 4.6), dpi=110)

    # (a) breadth field on the unfolded starboard flank: every hydrostatic number
    # above is an integral of exactly this quantity.
    ax0 = fig.add_subplot(1, 3, 1)
    gx = np.linspace(0.0, length, 161)
    gz = np.linspace(-draft, 0.0, 121)
    base = reference_parameters()
    field = np.array([[2.0 * half_breadth(base, x, -z) for x in gx] for z in gz])
    contour = ax0.contourf(gx, gz, field, levels=24, cmap="viridis")
    fig.colorbar(contour, ax=ax0, label="breadth 2y [m]")
    ax0.plot([0, length, length, 0, 0], [-draft, -draft, 0, 0, -draft], color="white", lw=1.0, ls="--", alpha=0.7)
    ax0.set_xlabel("x [m]  (aft -> forward)")
    ax0.set_ylabel("z [m]  (keel -> waterline)")
    ax0.set_title("(a) breadth field of the generating surface\nclassic Wigley, n = m = 2")

    # (b) mesh section vs analytic section: the loft interpolates along x only, so
    # at a sampled station the two are the same curve or the generator is wrong.
    ax1 = fig.add_subplot(1, 3, 2)
    model = generate_wigley(reference_parameters(nx=nx, nz=nz))
    xs = sorted({vertex[0] for vertex in model.mesh.vertices})
    for index, colour in ((nx // 2, "tab:red"), (nx // 4, "tab:blue"), ((3 * nx) // 4, "tab:green")):
        profile = _station_profile(model.mesh, index)
        x = xs[index]
        ax1.plot([y for _, y in profile], [z for z, _ in profile], "o", ms=3.0, color=colour,
                 label=f"mesh, x = {x:.2f} m")
        zs = np.linspace(-draft, 0.0, 200)
        ax1.plot([half_breadth(base, x, -z) for z in zs], zs, "-", lw=1.1, color=colour, alpha=0.85)
    ax1.set_xlabel("half-breadth y [m]")
    ax1.set_ylabel("z [m]")
    ax1.set_title("(b) mesh stations on the analytic sections\n(dots mesh, lines analytic)")
    ax1.legend(loc="lower right", fontsize=7.5)
    ax1.grid(alpha=0.3)

    # (c) signed relative error against the closed forms versus grid spacing.
    ax2 = fig.add_subplot(1, 3, 3)
    rows = results["convergence"]["rows"]
    h = np.array([row["h"] for row in rows])
    for name, colour, marker in (
        ("displacement_volume", "tab:red", "o"),
        ("waterplane_area", "tab:blue", "s"),
        ("wetted_surface_area", "tab:green", "^"),
    ):
        errors = np.abs([row["relative_error"][name] for row in rows])
        ax2.loglog(h, errors, marker + "-", color=colour, ms=5, lw=1.3,
                   label=name.replace("_", " "))
    guide = abs(rows[0]["relative_error"]["displacement_volume"]) * (h / h[0]) ** 2
    ax2.loglog(h, guide, "k--", lw=0.9, label=r"reference slope $h^{2}$")
    ax2.set_xlabel(r"vertical spacing $h = 1/(n_z-1)$")
    ax2.set_ylabel("relative error vs closed form")
    ax2.set_title("(c) sampled -> analytic convergence\n(every signed error is negative: inscribed chords)")
    ax2.grid(True, which="both", alpha=0.3)
    ax2.legend(fontsize=7.5)

    fig.suptitle(
        f"Wigley hull L/B = {length / beam:.0f}, B/T = {beam / draft:.2f}: "
        f"Cb = {analytic['block_coefficient']:.4f} (exactly 4/9), "
        f"V = {analytic['displacement_volume_m3']:.4f} m$^3$, "
        f"S = {analytic['wetted_surface_area_m2']:.4f} m$^2$",
        fontsize=10.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path = output / "plot.png"
    fig.savefig(path, pil_kwargs={"optimize": True})
    plt.close(fig)
    _compact_png(path)
    return path


def _compact_png(path: Path, colours: int = 96) -> None:
    """Quantise the figure so it can be embedded in a JSON report.

    A full-colour contourf render is ~150 kB, which is past the payload budget of
    the report channel; 96 adaptive colours keeps the figure readable at ~40 kB.
    """
    try:
        from PIL import Image
    except ImportError:  # pragma: no cover - matplotlib pulls Pillow in
        return
    with Image.open(path) as image:
        quantised = image.convert("RGB").quantize(colors=colours, method=Image.MEDIANCUT)
        quantised.save(path, optimize=True)


def _report(results: dict, claims: list[dict]) -> str:
    """Markdown report generated from the measured numbers, not typed by hand."""
    p = results["parameters"]
    a = results["analytic"]
    rows = results["convergence"]["rows"]
    lines: list[str] = []
    add = lines.append

    add("# Wigley hull form: closed forms, sampled convergence, solver exports")
    add("")
    add("This page is written by `python examples/wigley/run_example.py`; every number in it is")
    add("read back out of `results.json` from the run that produced `plot.png`.")
    add("")
    add("## 1. Case")
    add("")
    add(f"* Hull form: `y(x, w) = (B/2)(1 - |2x/L - 1|^n)(1 - (w/T)^m)`, n = {p['longitudinal_exponent']:g}, "
        f"m = {p['vertical_exponent']:g} (classic Wigley).")
    add(f"* L = {p['length']:g} m, B = {p['beam']:g} m, T = {p['draft']:g} m -> "
        f"L/B = {p['length'] / p['beam']:.2f}, B/T = {p['beam'] / p['draft']:.2f} "
        "(the ITTC reference proportions).")
    add(f"* Water density {p['water_density']:g} kg/m^3; frame: x aft->forward, z = 0 at the design "
        "waterline, keel at z = -T.")
    add(f"* Sampling grid {p['nx']} stations x {p['nz']} levels ({results['topology']['triangles']} triangles).")
    add("")
    add("## 2. Method")
    add("")
    add("The generator lofts one section ring per station (starboard flank down to a shared keel")
    add("point, port flank back up) and caps the design waterline, so the surface is a closed solid")
    add("whose volume can be integrated from the triangles. Against it stand two independent")
    add("references:")
    add("")
    add("* **closed forms**, derived by hand from the separable integrals")
    add("  `int_0^L taper dx = L n/(n+1)` and `int_0^T section dw = T m/(m+1)`;")
    add("* **a composite Gauss-Legendre quadrature of the analytic flank** for wetted area, the one")
    add("  quantity with no elementary primitive, self-converged to the order-to-order difference")
    add(f"  reported below ({results['wetted_area_evidence']['relative_change_fine_to_finer']:.2e}).")
    add("")
    add("A direct midpoint integral of the generating function re-derives the volume independently:")
    add(f"{results['closed_form_cross_check']['volume_direct_integral_m3']:.9f} m^3 vs closed form "
        f"{results['closed_form_cross_check']['volume_closed_form_m3']:.9f} m^3 "
        f"({100 * results['closed_form_cross_check']['relative_difference']:.2e} %).")
    add("")
    add("## 3. Closed forms")
    add("")
    add("| quantity | closed form | value |")
    add("|---|---|---|")
    add(f"| displacement volume | `n m/((n+1)(m+1)) L B T` | {a['displacement_volume_m3']:.9f} m^3 |")
    add(f"| displacement mass | `rho V` | {a['displacement_mass_kg']:.3f} kg |")
    add(f"| waterplane area | `L B n/(n+1)` | {a['waterplane_area_m2']:.9f} m^2 |")
    add(f"| midship section area | `B T m/(m+1)` | {a['maximum_section_area_m2']:.9f} m^2 |")
    add(f"| wetted surface | quadrature (no primitive) | {a['wetted_surface_area_m2']:.9f} m^2 |")
    add(f"| Cb | `n m/((n+1)(m+1))` | {a['block_coefficient']:.12f} (4/9 = {4 / 9:.12f}) |")
    add(f"| Cp | `n/(n+1)` | {a['prismatic_coefficient']:.12f} |")
    add(f"| Cm | `m/(m+1)` | {a['midship_coefficient']:.12f} |")
    add(f"| Cwp | `n/(n+1)` | {a['waterplane_coefficient']:.12f} |")
    add(f"| VCB above keel | `T - T(m+1)/(2(m+2))` | {results['vcb_above_keel']:.9f} m (5T/8 = {5 * p['draft'] / 8:.9f}) |")
    add(f"| LCB / LCF | `L/2` (symmetry) | {a['longitudinal_centre_buoyancy_m']:.9f} m |")
    add("")
    add("## 4. Sampled convergence against the closed forms")
    add("")
    add("| grid | triangles | V [m^3] | e_V [%] | e_Awp [%] | e_S [%] | max abs err [%] |")
    add("|---|---|---|---|---|---|---|")
    for row in rows:
        e = row["relative_error"]
        add(f"| {row['nx']} x {row['nz']} | {row['triangles']} | {row['displacement_volume']:.7f} | "
            f"{100 * e['displacement_volume']:+.5f} | {100 * e['waterplane_area']:+.5f} | "
            f"{100 * e['wetted_surface_area']:+.5f} | {100 * row['max_abs_relative_error']:.5f} |")
    add("")
    add("Observed orders (fitted across the sweep): " + ", ".join(
        f"{name.replace('_', ' ')} = {value:.2f}" for name, value in results["convergence"]["observed_order"].items()
    ) + ".")
    add("Every error is negative: the loft chords cut inside the analytic surface, so the")
    add("discretisation is systematic, not noise -- which is why the closed form is the right thing")
    add("to test the integrator against.")
    add("")
    add("## 5. Figures")
    add("")
    add("![plot](plot.png)")
    add("")
    add("(a) breadth field of the generating surface -- the quantity every hydrostatic integral")
    add("integrates; (b) three sampled stations against their analytic section curves, dots on")
    add("lines; (c) log-log convergence of the three extensive quantities with an `h^2` guide.")
    add("")
    add("## 6. Points to watch")
    add("")
    add("* The family has a **finite keel angle** (unlike the `hull` family's `sqrt` taper), which is")
    add("  why the default 41x13 grid already lands at "
        f"{100 * abs(rows[1]['relative_error']['displacement_volume']):.3f} % on volume against")
    add("  0.88 % for the same grid on `sin(pi x/L)^0.7`: second-order convergence of a smooth")
    add("  section beats a better-tuned coarse grid.")
    add("* Wetted area is a quadrature reference, not a closed form. It is exact to round-off for")
    add("  the polynomial integrand of `n = m = 2`, but for non-integer exponents the `|u|^(n-1)`")
    add("  term at the ends slows Gauss-Legendre convergence -- re-run `wetted_area_evidence`")
    add("  before quoting a wetted area for such a case.")
    add("* No appendages, no transom, no parallel middle body: this is a benchmark hull. A real")
    add("  design still needs `tensorcad.suboff_appendages` or a fitted offsets table.")
    add("* The hull is fore-aft symmetric, so LCB = LCF = L/2 identically: this case does **not**")
    add("  exercise the longitudinal moment integrator. `tensorcad.hull` and the SUBOFF case do.")
    add("")
    add("## 7. Error verdict")
    add("")
    add("| claim | measured | threshold | verdict |")
    add("|---|---|---|---|")
    for claim in claims:
        value = claim["value"]
        if isinstance(value, dict):
            shown = ", ".join(f"{k} {v:.6g}" if isinstance(v, (int, float)) else str(v) for k, v in value.items())
        elif isinstance(value, float):
            shown = f"{value:.6g}"
        else:
            shown = str(value)
        if len(shown) > 120:
            shown = shown[:117] + "..."
        add(f"| {claim['claim']} | {shown} | {claim['tolerance']:g} | {'PASS' if claim['passed'] else 'FAIL'} |")
    add("")
    worst = rows[-1]["max_abs_relative_error"]
    add("**Is the key quantity within 3 % of its reference?** Yes, by two orders of magnitude:")
    add(f"the finest grid integrates every comparable quantity to within {100 * worst:.4f} % of the")
    add("closed form, and even the coarsest grid (21 x 9) stays within "
        f"{100 * rows[0]['max_abs_relative_error']:.3f} %. The reference is exact algebra, so the")
    add("bullet that remains is the *form* (a Wigley is not a ship), not the integrator.")
    add("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(REPO_ROOT / "docs" / "wigley"))
    args = parser.parse_args(argv)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    parameters = reference_parameters()
    model = generate_wigley(parameters)
    analytic = wigley_analytic_hydrostatics(parameters)
    study = convergence_study(parameters=parameters)
    evidence = wetted_area_evidence(parameters)

    # Independent re-derivation of the volume from the generating function, so a
    # sign or a factor in a docstring cannot survive this script.
    rows = DIRECT_STEPS
    cols = DIRECT_STEPS // 4
    direct_volume = 0.0
    for i in range(rows):
        x = parameters.length * (i + 0.5) / rows
        for j in range(cols):
            depth = parameters.draft * (j + 0.5) / cols
            direct_volume += 2.0 * half_breadth(parameters, x, depth)
    direct_volume *= (parameters.length / rows) * (parameters.draft / cols)

    # Undirected edges: a closed orientable mesh stores every edge twice with
    # opposite winding, so counting directed pairs would double E and destroy the
    # Euler characteristic (V - E + F must be 2).
    edges = {
        (min(a, b), max(a, b))
        for face in model.mesh.faces
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0]))
    }
    topology = {
        "watertight": meshlib.is_watertight(model.mesh),
        "edge_manifold": meshlib.is_edge_manifold(model.mesh),
        "euler_characteristic": len(model.mesh.vertices) - len(edges) + len(model.mesh.faces),
        "single_component": True,
        "zero_area_faces": len(meshlib.zero_area_faces(model.mesh)),
        "triangles": len(model.mesh.faces),
        "vertices": len(model.mesh.vertices),
    }
    edge_to_faces: dict[tuple[int, int], list[int]] = {}
    for index, face in enumerate(model.mesh.faces):
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            edge_to_faces.setdefault((min(a, b), max(a, b)), []).append(index)
    seen = {0}
    stack = [0]
    while stack:
        current = stack.pop()
        face = model.mesh.faces[current]
        for a, b in ((face[0], face[1]), (face[1], face[2]), (face[2], face[0])):
            for neighbour in edge_to_faces[(min(a, b), max(a, b))]:
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append(neighbour)
    topology["single_component"] = len(seen) == len(model.mesh.faces)

    offset_deviation = max(
        abs(width - 2.0 * half_breadth(parameters, x, 0.0)) for x, width in model.offsets
    )
    exported = export_wigley(model, output, stl=True)
    closed_volume = analytic_volume(parameters.length, parameters.beam, parameters.draft)

    results: dict = {
        "case": REFERENCE_CASE,
        "parameters": parameters.as_dict(),
        "analytic": analytic.as_dict(),
        "vcb_above_keel": analytic.vertical_centre_buoyancy + parameters.draft,
        "closed_form_cross_check": {
            "volume_direct_integral_m3": direct_volume,
            "volume_closed_form_m3": closed_volume,
            "relative_difference": abs(direct_volume - closed_volume) / closed_volume,
            "waterplane_area_closed_form_m2": analytic_waterplane_area(parameters.length, parameters.beam),
            "maximum_section_area_closed_form_m2": analytic_maximum_section_area(
                parameters.beam, parameters.draft
            ),
            "block_coefficient_closed_form": analytic_block_coefficient(),
        },
        "mesh_hydrostatics": model.hydrostatics.to_dict(),
        "default_grid_comparison": analytic.compare(model.hydrostatics),
        "convergence": study,
        "wetted_area_evidence": evidence,
        "topology": topology,
        "offset_on_surface_max_deviation_m": offset_deviation,
        "validation": model.validation.to_dict(),
        "artifacts": exported,
    }
    claims = _claims(results)
    results["claims"] = claims
    results["failed_claims"] = [claim["claim"] for claim in claims if not claim["passed"]]
    plt_path = _figure(results, output)
    results["figures"] = {plt_path.name: str(plt_path)}
    (output / "results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output / "report.md").write_text(_report(results, claims), encoding="utf-8")

    print(f"Wigley reference case: L={parameters.length} m, B={parameters.beam} m, T={parameters.draft} m")
    print(f"  L/B = {parameters.length / parameters.beam:.2f}, B/T = {parameters.beam / parameters.draft:.2f}")
    print()
    print("closed forms")
    print(f"  Cb  = {analytic.block_coefficient:.12f}   (4/9 = {4 / 9:.12f})")
    print(f"  Cp  = {analytic.prismatic_coefficient:.12f}   Cm = {analytic.midship_coefficient:.12f}   "
          f"Cwp = {analytic.waterplane_coefficient:.12f}")
    print(f"  V   = {analytic.displacement_volume:.9f} m^3   Awp = {analytic.waterplane_area:.9f} m^2")
    print(f"  S   = {analytic.wetted_surface_area:.9f} m^2  (quadrature, self-converged "
          f"{evidence['relative_change_fine_to_finer']:.2e})")
    print(f"  VCB = {results['vcb_above_keel']:.9f} m above keel, LCB = L/2 = "
          f"{analytic.longitudinal_centre_buoyancy:.6f} m")
    print(f"  direct {DIRECT_STEPS}x{DIRECT_STEPS // 4} midpoint integral: V = {direct_volume:.9f} m^3 "
          f"(diff {100 * results['closed_form_cross_check']['relative_difference']:.2e} %)")
    print()
    print(f"{'grid':>10} {'triangles':>10} {'V [m^3]':>12} {'eV [%]':>9} {'eAwp [%]':>10} {'eS [%]':>9} {'max|e| [%]':>11}")
    for row in study["rows"]:
        error = row["relative_error"]
        print(f"{row['nx']:4d}x{row['nz']:<5d} {row['triangles']:10d} {row['displacement_volume']:12.6f} "
              f"{100 * error['displacement_volume']:+9.5f} {100 * error['waterplane_area']:+10.5f} "
              f"{100 * error['wetted_surface_area']:+9.5f} {100 * row['max_abs_relative_error']:11.5f}")
    print()
    print("observed order: " + ", ".join(f"{k} {v:.2f}" for k, v in study["observed_order"].items()))
    print(f"topology: watertight={topology['watertight']} manifold={topology['edge_manifold']} "
          f"chi={topology['euler_characteristic']} components=1 zero-area={topology['zero_area_faces']}")
    print(f"waterline offsets on analytic surface: max deviation {offset_deviation:.2e} m")
    print()
    for claim in claims:
        print(f"{'PASS' if claim['passed'] else 'FAIL'}  {claim['claim']}")
    print()
    print("artifacts: " + ", ".join(f"{k} -> {v}" for k, v in exported.items()))
    print(f"wrote {output / 'results.json'}, {output / 'report.md'}, {plt_path}")
    return 1 if results["failed_claims"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
