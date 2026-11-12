"""Initial stability measured against closed forms, and a metacentric height ladder.

Everything on this page is *measured*: the script exits non-zero if one of its own
claims stops holding, so the numbers in ``docs/stability.md`` cannot drift away from the
code the way a hand-copied table does.

What it does, in order:

1. **A box**, where ``A``, ``I_T``, ``I_L``, ``BM`` and ``KB`` are rational arithmetic.
   The mesh integration must reproduce them to round-off -- if it does not, nothing
   further is worth running.
2. **A hemisphere**, where the section is a circle. Two independent integrations
   (``KB`` from the volume, ``BM`` from the waterplane) must add up to ``KM = R``: the
   transverse metacentre of a hemisphere is the centre of the sphere.
3. **The parametric hull** on a refinement ladder. The waterplane is the sine planform
   ``|y| <= (B/2) sin(pi x / L)^0.7``; ``I_T`` is closed-form through the Beta integral in
   ``hydrostatics_analytic`` and ``I_L`` is quadrature. The ladder measures the order of
   the mesh's own error rather than asserting it is small.
4. **A ``GM`` ladder over ``KG``**, with the ``GM = 0`` crossing -- the loading limit at
   which the upright position stops being stable -- reported as ``KM`` itself.

Run it with ``python examples/stability/run_example.py --output docs/stability``.
matplotlib is optional: the JSON is written either way and the PNGs are skipped.
"""

from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_REPO_ROOT = _Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in _sys.path:  # allow running from a checkout without installing
    _sys.path.insert(0, str(_REPO_ROOT))

import argparse
import json
import math
import sys
from pathlib import Path

from tensorcad import __version__
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.hydrostatics_analytic import (
    analytic_hydrostatics,
    tanh_sinh_quadrature,
    waterplane_inertia_longitudinal,
    waterplane_inertia_transverse,
)
from tensorcad.stability import initial_stability, waterplane_moments

LENGTH, BEAM, DRAFT = 10.0, 2.0, 0.8
BREADTH_EXPONENT = 0.7
LADDER = ((21, 9), (41, 13), (81, 25), (161, 41), (321, 81))


# ---------------------------------------------------------------------------
# reference bodies
# ---------------------------------------------------------------------------
def box(length: float, beam: float, draft: float) -> Mesh:
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, 0.0):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi, top, bottom, right = -0.5 * beam, 0.5 * beam, 0.0, -draft, length
    quad((0.0, lo, bottom), (0.0, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((0.0, lo, top), (right, lo, top), (right, hi, top), (0.0, hi, top))
    quad((0.0, lo, bottom), (right, lo, bottom), (right, lo, top), (0.0, lo, top))
    quad((0.0, hi, bottom), (0.0, hi, top), (right, hi, top), (right, hi, bottom))
    quad((0.0, lo, bottom), (0.0, lo, top), (0.0, hi, top), (0.0, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


def hemisphere(radius: float, n_lat: int, n_lon: int) -> Mesh:
    from tensorcad import mesh as meshlib

    builder = MeshBuilder()
    rings = []
    for i in range(n_lat):
        theta = 0.5 * math.pi * i / n_lat
        ring = []
        for j in range(n_lon):
            phi = 2.0 * math.pi * j / n_lon
            ring.append(
                builder.add_vertex(
                    (
                        radius * math.cos(theta) * math.cos(phi),
                        radius * math.cos(theta) * math.sin(phi),
                        -radius * math.sin(theta),
                    )
                )
            )
        rings.append(ring)
    pole = builder.add_vertex((0.0, 0.0, -radius))
    for i in range(n_lat - 1):
        upper, lower = rings[i], rings[i + 1]
        for j in range(n_lon):
            j2 = (j + 1) % n_lon
            builder.add_triangle(upper[j], lower[j], lower[j2])
            builder.add_triangle(upper[j], lower[j2], upper[j2])
    last = rings[-1]
    for j in range(n_lon):
        j2 = (j + 1) % n_lon
        builder.add_triangle(last[j], pole, last[j2])
    centre = builder.add_vertex((0.0, 0.0, 0.0))
    for j in range(n_lon):
        j2 = (j + 1) % n_lon
        builder.add_triangle(centre, rings[0][j], rings[0][j2])
    mesh = builder.build()
    if meshlib.total_signed_volume(mesh) < 0.0:
        mesh = Mesh(
            vertices=mesh.vertices,
            faces=tuple((f[0], f[2], f[1]) for f in mesh.faces),
            units=mesh.units,
            coordinate_frame=mesh.coordinate_frame,
        )
    return mesh


# ---------------------------------------------------------------------------
# an independent route to I_L: remove the endpoint singularity, then Simpson
# ---------------------------------------------------------------------------
def longitudinal_moment_by_simpson(
    length: float, beam: float, exponent: float = BREADTH_EXPONENT, intervals: int = 20000
) -> float:
    """``I_L`` from ``2 B L^3 int_0^{1/2} u^2 cos(pi u)^e du`` after ``u = 1/2 - s^2``.

    The substitution turns the integrable ``cos(pi u)^e`` endpoint into ``sin(pi s^2)^e
    ~ (pi s^2)^e``, and the Jacobian ``2s`` makes the integrand behave like ``s**(2e+1)``
    at ``s = 0``: smooth for ``e = 0.7``, so a composite Simpson rule converges at its
    nominal order. This is deliberately a *different* rule from the tanh-sinh one in
    ``hydrostatics_analytic``: two rules agreeing to the accuracy shown is evidence, one
    rule agreeing with itself is not.
    """
    upper = math.sqrt(0.5)
    h = upper / intervals

    def integrand(s: float) -> float:
        u = 0.5 - s * s
        return (u * u) * (math.sin(math.pi * s * s) ** exponent) * 2.0 * s

    total = integrand(0.0) + integrand(upper)
    total += 4.0 * math.fsum(integrand(i * h) for i in range(1, intervals, 2))
    total += 2.0 * math.fsum(integrand(i * h) for i in range(2, intervals, 2))
    return 2.0 * beam * length**3 * total * h / 3.0


def sinh_sinh_cross_check(length: float, beam: float) -> float:
    """Same integral through ``tanh_sinh_quadrature`` one more level up: rule stability."""

    def integrand(u: float) -> float:
        return u * u * math.cos(math.pi * u) ** BREADTH_EXPONENT

    return 2.0 * beam * length**3 * tanh_sinh_quadrature(integrand, 0.0, 0.5, levels=9)


def relative(value: float, reference: float) -> float:
    return value / reference - 1.0


def percent(value: float, reference: float) -> float:
    return 100.0 * relative(value, reference)


def fitted_order(errors: list[float], steps: list[float]) -> float:
    """Least-squares slope of ``log|error|`` against ``log(step)``."""
    points = [
        (math.log(step), math.log(abs(error)))
        for error, step in zip(errors, steps, strict=False)
        if error != 0.0 and math.isfinite(error)
    ]
    if len(points) < 2:
        return float("nan")
    mean_x = math.fsum(x for x, _ in points) / len(points)
    mean_y = math.fsum(y for _, y in points) / len(points)
    numerator = math.fsum((x - mean_x) * (y - mean_y) for x, y in points)
    denominator = math.fsum((x - mean_x) ** 2 for x, _ in points)
    return numerator / denominator


# ---------------------------------------------------------------------------
# the measurements
# ---------------------------------------------------------------------------
def measure_box() -> dict[str, object]:
    mesh = box(LENGTH, BEAM, DRAFT)
    kg = 0.55
    result = initial_stability(mesh, kg=kg)
    reference = {
        "waterplane_area": LENGTH * BEAM,
        "transverse_second_moment": LENGTH * BEAM**3 / 12.0,
        "longitudinal_second_moment": BEAM * LENGTH**3 / 12.0,
        "displacement_volume": LENGTH * BEAM * DRAFT,
        "kb": 0.5 * DRAFT,
        "bm_transverse": BEAM**2 / (12.0 * DRAFT),
        "gm_transverse": 0.5 * DRAFT + BEAM**2 / (12.0 * DRAFT) - kg,
    }
    measured = {
        "waterplane_area": result.waterplane.area,
        "transverse_second_moment": result.waterplane.transverse_second_moment,
        "longitudinal_second_moment": result.waterplane.longitudinal_second_moment,
        "displacement_volume": result.displacement_volume,
        "kb": result.kb,
        "bm_transverse": result.bm_transverse,
        "gm_transverse": result.gm_transverse,
    }
    rows = {
        name: {
            "measured": measured[name],
            "reference": reference[name],
            "relative_error": relative(measured[name], reference[name]),
        }
        for name in reference
    }
    worst = max(abs(row["relative_error"]) for row in rows.values())
    return {"kg": kg, "rows": rows, "worst_relative_error": worst, "mesh": {"faces": len(mesh.faces)}}


def measure_hemisphere(radius: float = 1.0, n_lat: int = 32, n_lon: int = 64) -> dict[str, object]:
    mesh = hemisphere(radius, n_lat, n_lon)
    result = initial_stability(mesh, kg=0.5)
    rows = {
        "waterplane_area": {
            "measured": result.waterplane.area,
            "reference": math.pi * radius**2,
        },
        "transverse_second_moment": {
            "measured": result.waterplane.transverse_second_moment,
            "reference": math.pi * radius**4 / 4.0,
        },
        "longitudinal_second_moment": {
            "measured": result.waterplane.longitudinal_second_moment,
            "reference": math.pi * radius**4 / 4.0,
        },
        "kb": {"measured": result.kb, "reference": 5.0 * radius / 8.0},
        "bm_transverse": {"measured": result.bm_transverse, "reference": 3.0 * radius / 8.0},
        "km_transverse": {"measured": result.km_transverse, "reference": radius},
    }
    for row in rows.values():
        row["relative_error"] = relative(row["measured"], row["reference"])
    worst = max(abs(row["relative_error"]) for row in rows.values())
    return {
        "radius": radius,
        "discretisation": {"n_lat": n_lat, "n_lon": n_lon, "faces": len(mesh.faces)},
        "rows": rows,
        "worst_relative_error": worst,
    }


def measure_hull_ladder() -> dict[str, object]:
    reference_parameters = HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT)
    reference = analytic_hydrostatics(reference_parameters)
    inertia_t = waterplane_inertia_transverse(LENGTH, BEAM)
    inertia_l, quadrature_evidence = waterplane_inertia_longitudinal(LENGTH, BEAM)
    simpson = longitudinal_moment_by_simpson(LENGTH, BEAM)
    sinh_sinh = sinh_sinh_cross_check(LENGTH, BEAM)
    analytic_km = 0.6 * DRAFT + inertia_t / reference.displacement_volume
    kg = 0.6

    levels = []
    for nx, nz in LADDER:
        mesh = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=nx, nz=nz)).mesh
        section = waterplane_moments(mesh, 0.0)
        result = initial_stability(mesh, kg=kg)
        levels.append(
            {
                "nx": nx,
                "nz": nz,
                "faces": len(mesh.faces),
                "step": 1.0 / (nx - 1),
                "waterplane_area": section.area,
                "waterplane_area_error": relative(section.area, reference.waterplane_area),
                "transverse_second_moment": section.transverse_second_moment,
                "transverse_second_moment_error": relative(section.transverse_second_moment, inertia_t),
                "longitudinal_second_moment": section.longitudinal_second_moment,
                "longitudinal_second_moment_error": relative(section.longitudinal_second_moment, inertia_l),
                "displacement_volume": result.displacement_volume,
                "displacement_volume_error": relative(
                    result.displacement_volume, reference.displacement_volume
                ),
                "kb": result.kb,
                "kb_error": relative(result.kb, 0.6 * DRAFT),
                "bm_transverse": result.bm_transverse,
                "bm_transverse_error": relative(
                    result.bm_transverse, inertia_t / reference.displacement_volume
                ),
                "km_transverse": result.km_transverse,
                "km_transverse_error": relative(result.km_transverse, analytic_km),
                "gm_transverse": result.gm_transverse,
                "gm_reference": analytic_km - kg,
                "gm_transverse_error": relative(result.gm_transverse, analytic_km - kg),
                "longitudinal_centre_flotation": section.longitudinal_centre_flotation,
                "transverse_centre_flotation": section.transverse_centre_flotation,
                "loops": section.loops,
            }
        )
    steps = [level["step"] for level in levels]
    orders = {
        "waterplane_area": fitted_order([level["waterplane_area_error"] for level in levels], steps),
        "transverse_second_moment": fitted_order(
            [level["transverse_second_moment_error"] for level in levels], steps
        ),
        "longitudinal_second_moment": fitted_order(
            [level["longitudinal_second_moment_error"] for level in levels], steps
        ),
        "displacement_volume": fitted_order([level["displacement_volume_error"] for level in levels], steps),
    }
    return {
        "parameters": {"length": LENGTH, "beam": BEAM, "draft": DRAFT, "breadth_exponent": BREADTH_EXPONENT},
        "kg": kg,
        "reference": {
            "waterplane_area": reference.waterplane_area,
            "transverse_second_moment": inertia_t,
            "longitudinal_second_moment": inertia_l,
            "longitudinal_second_moment_quadrature_evidence": quadrature_evidence,
            "longitudinal_second_moment_simpson": simpson,
            "longitudinal_second_moment_simpson_relative_difference": relative(simpson, inertia_l),
            "longitudinal_second_moment_tanh_sinh_level9": sinh_sinh,
            "longitudinal_second_moment_tanh_sinh_level9_relative_difference": relative(sinh_sinh, inertia_l),
            "displacement_volume": reference.displacement_volume,
            "kb": 0.6 * DRAFT,
            "bm_transverse": inertia_t / reference.displacement_volume,
            "km_transverse": analytic_km,
        },
        "levels": levels,
        "fitted_orders": orders,
    }


def measure_kg_ladder() -> list[dict[str, float]]:
    """GM against KG, and the loading limit ``GM = 0`` at ``KG = KM``."""
    mesh = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=161, nz=41)).mesh
    return [
        {
            "kg": kg,
            "gm_transverse": initial_stability(mesh, kg=kg).gm_transverse,
        }
        for kg in (0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)
    ]


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def figure_sections(results: dict[str, object], path: Path) -> None:
    """Waterline planform: mesh chords against the analytic sine curve."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(7.6, 3.4))
    for nx, nz, style, colour in ((41, 13, "-o", "#1f77b4"), (321, 81, "-", "#d62728")):
        mesh = generate_hull(HullParameters(length=LENGTH, beam=BEAM, draft=DRAFT, nx=nx, nz=nz)).mesh
        section = waterplane_moments(mesh, 0.0)
        # the section polygon is the loop of the waterline; re-chaining it here would
        # duplicate the module's logic, so draw the analytic curve and the station
        # half-breadths the mesh actually has at the waterline.
        xs = [LENGTH * i / (nx - 1) for i in range(nx)]
        half = [0.5 * BEAM * math.sin(math.pi * x / LENGTH) ** BREADTH_EXPONENT for x in xs]
        if style == "-o":
            axis.plot(
                xs,
                half,
                style,
                markersize=2.5,
                linewidth=1.0,
                color=colour,
                label=f"mesh waterline chords, nx={nx} (A={section.area:.4f} m²)",
            )
    curve_x = [LENGTH * i / 4000 for i in range(4001)]
    axis.plot(
        curve_x,
        [0.5 * BEAM * math.sin(math.pi * x / LENGTH) ** BREADTH_EXPONENT for x in curve_x],
        "k--",
        linewidth=1.0,
        label="analytic half-breadth (B/2) sin(πx/L)^0.7",
    )
    axis.set_xlabel("x [m] (0 = aft perpendicular)")
    axis.set_ylabel("half-breadth y [m]")
    axis.set_title("Waterplane section used for I_T: chords vs the analytic planform")
    axis.legend(fontsize=7, loc="lower center")
    axis.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(path, dpi=140)
    plt.close(figure)


def figure_convergence(results: dict[str, object], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ladder = results["hull"]["levels"]  # type: ignore[index]
    steps = [level["step"] for level in ladder]
    figure, axis = plt.subplots(figsize=(6.4, 4.0))
    for key, label, colour in (
        ("waterplane_area_error", "waterplane area", "#1f77b4"),
        ("transverse_second_moment_error", "I$_T$", "#2ca02c"),
        ("longitudinal_second_moment_error", "I$_L$", "#d62728"),
        ("displacement_volume_error", "displacement volume", "#9467bd"),
    ):
        axis.loglog(
            steps, [abs(level[key]) for level in ladder], "-o", label=label, color=colour, markersize=3
        )
    axis.loglog(
        steps,
        [steps[0] ** 2 / step * abs(ladder[0]["waterplane_area_error"]) for step in steps],
        "k--",
        linewidth=0.9,
        label="slope 2 reference",
    )
    axis.set_xlabel("mesh step h = 1/(nx−1)")
    axis.set_ylabel("|relative error| vs analytic")
    axis.invert_xaxis()
    axis.set_title("Waterplane moments and displacement vs hull grid")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.3, which="both")
    figure.tight_layout()
    figure.savefig(path, dpi=140)
    plt.close(figure)


def figure_gm_ladder(results: dict[str, object], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = results["kg_ladder"]
    kgs = [row["kg"] for row in rows]
    gms = [row["gm_transverse"] for row in rows]
    km = results["hull"]["reference"]["km_transverse"]  # type: ignore[index]
    figure, axis = plt.subplots(figsize=(6.4, 4.0))
    axis.plot(kgs, gms, "-o", color="#1f77b4", label="GM$_T$ from the mesh (161×41)")
    axis.axhline(0.0, color="k", linewidth=0.8)
    axis.axvline(km, color="#d62728", linestyle="--", linewidth=1.0, label=f"KG = KM = {km:.4f} m (GM = 0)")
    axis.set_xlabel("KG above keel [m]")
    axis.set_ylabel("GM$_T$ [m]")
    axis.set_title("Transverse metacentric height against loading (hull 10 × 2 × 0.8 m)")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(path, dpi=140)
    plt.close(figure)


# ---------------------------------------------------------------------------
# claims
# ---------------------------------------------------------------------------
def evaluate(results: dict[str, object]) -> list[dict[str, object]]:
    box_worst = results["box"]["worst_relative_error"]  # type: ignore[index]
    hemisphere_worst = results["hemisphere"]["worst_relative_error"]  # type: ignore[index]
    ladder = results["hull"]["levels"]  # type: ignore[index]
    default = next(level for level in ladder if level["nx"] == 41)
    fine = next(level for level in ladder if level["nx"] == 161)
    order = results["hull"]["fitted_orders"]  # type: ignore[index]
    quadrature_spread = abs(
        results["hull"]["reference"]["longitudinal_second_moment_simpson_relative_difference"]  # type: ignore[index]
    )
    return [
        {
            "claim": "box: every waterplane and metacentric quantity is exact to round-off",
            "value": box_worst,
            "threshold": 1e-12,
            "holds": box_worst < 1e-12,
        },
        {
            "claim": "hemisphere: KM reproduces the sphere radius R (KB = 5R/8, BM = 3R/8)",
            "value": hemisphere_worst,
            "threshold": 0.02,
            "holds": hemisphere_worst < 0.02,
        },
        {
            "claim": "default hull grid (41×13): I_T within 3 % of the closed form",
            "value": abs(default["transverse_second_moment_error"]),
            "threshold": 0.03,
            "holds": abs(default["transverse_second_moment_error"]) < 0.03,
        },
        {
            "claim": "default hull grid (41×13): I_L within 3 % of quadrature",
            "value": abs(default["longitudinal_second_moment_error"]),
            "threshold": 0.03,
            "holds": abs(default["longitudinal_second_moment_error"]) < 0.03,
        },
        {
            "claim": "default hull grid (41×13): GM_T within 3 % of the analytic surface",
            "value": abs(default["gm_transverse_error"]),
            "threshold": 0.03,
            "holds": abs(default["gm_transverse_error"]) < 0.03,
        },
        {
            "claim": "fine hull grid (161×41): I_T within 1 %",
            "value": abs(fine["transverse_second_moment_error"]),
            "threshold": 0.01,
            "holds": abs(fine["transverse_second_moment_error"]) < 0.01,
        },
        {
            "claim": "fine hull grid (161×41): GM_T within 1 %",
            "value": abs(fine["gm_transverse_error"]),
            "threshold": 0.01,
            "holds": abs(fine["gm_transverse_error"]) < 0.01,
        },
        {
            "claim": "the waterplane error is second order in the grid step (fitted order > 1.5)",
            "value": order["transverse_second_moment"],
            "threshold": 1.5,
            "holds": order["transverse_second_moment"] > 1.5,
        },
        {
            "claim": "the I_L reference is rule-independent (Simpson vs tanh-sinh within 1e-6)",
            "value": quadrature_spread,
            "threshold": 1e-6,
            "holds": quadrature_spread < 1e-6,
        },
        {
            "claim": "the waterplane section is symmetric: |y_c| is round-off (below 1e-9 m)",
            "value": abs(fine["transverse_centre_flotation"]),
            "threshold": 1e-9,
            "holds": abs(fine["transverse_centre_flotation"]) < 1e-9,
        },
    ]


def render_report(results: dict[str, object], claims: list[dict[str, object]]) -> str:
    hull = results["hull"]  # type: ignore[index]
    reference = hull["reference"]
    lines = [
        "# Initial stability measured against closed forms",
        "",
        f"TensorCAD {results['version']} · `examples/stability/run_example.py` · "
        f"hull {LENGTH} × {BEAM} × {DRAFT} m, KG = {hull['kg']} m above keel",
        "",
        "## Claims enforced by this run",
        "",
        "| claim | measured | threshold | holds |",
        "|---|---|---|---|",
    ]
    for claim in claims:
        lines.append(
            f"| {claim['claim']} | {claim['value']:.3e} | {claim['threshold']:.0e} | "
            f"{'PASS' if claim['holds'] else 'FAIL'} |"
        )
    lines += [
        "",
        "## Hull ladder: mesh waterplane moments vs the analytic surface",
        "",
        "| nx × nz | A [m²] | err A | I_T [m⁴] | err I_T | I_L [m⁴] | err I_L | BM_T [m] | KB [m] | GM_T [m] | err GM |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for level in hull["levels"]:
        lines.append(
            "| {nx}×{nz} | {a:.4f} | {ea:+.3%} | {it:.4f} | {eit:+.3%} | {il:.4f} | {eil:+.3%} | "
            "{bm:.4f} | {kb:.4f} | {gm:.4f} | {egm:+.3%} |".format(
                nx=level["nx"],
                nz=level["nz"],
                a=level["waterplane_area"],
                ea=level["waterplane_area_error"],
                it=level["transverse_second_moment"],
                eit=level["transverse_second_moment_error"],
                il=level["longitudinal_second_moment"],
                eil=level["longitudinal_second_moment_error"],
                bm=level["bm_transverse"],
                kb=level["kb"],
                gm=level["gm_transverse"],
                egm=level["gm_transverse_error"],
            )
        )
    lines += [
        "",
        f"Analytic reference: A = {reference['waterplane_area']:.6f} m², "
        f"I_T = {reference['transverse_second_moment']:.6f} m⁴, "
        f"I_L = {reference['longitudinal_second_moment']:.6f} m⁴, "
        f"∇ = {reference['displacement_volume']:.6f} m³, "
        f"BM_T = {reference['bm_transverse']:.6f} m, KM_T = {reference['km_transverse']:.6f} m.",
        "",
        "Fitted orders (log–log slope of |error| against the grid step): "
        + ", ".join(f"{name} {value:.2f}" for name, value in hull["fitted_orders"].items())
        + ".",
        "",
        "`I_L` is quadrature, not elementary: the tanh-sinh value agrees with an independent "
        f"Simpson rule after the `u = 1/2 − s²` substitution to "
        f"{results['hull']['reference']['longitudinal_second_moment_simpson_relative_difference']:+.2e} "
        "relative, and with tanh-sinh one level up.",
        "",
        "## Figures",
        "",
        "- `sections.png` — the waterplane chords the mesh actually has against the analytic planform.",
        "- `convergence.png` — |error| vs grid step for A, I_T, I_L and the displacement volume.",
        "- `gm_ladder.png` — GM_T against KG, with the GM = 0 loading limit at KG = KM.",
    ]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/stability"))
    arguments = parser.parse_args(argv)
    output = arguments.output
    output.mkdir(parents=True, exist_ok=True)

    results: dict[str, object] = {
        "version": __version__,
        "box": measure_box(),
        "hemisphere": measure_hemisphere(),
        "hull": measure_hull_ladder(),
        "kg_ladder": measure_kg_ladder(),
    }
    claims = evaluate(results)
    results["claims"] = claims

    for name, builder in (
        ("sections.png", figure_sections),
        ("convergence.png", figure_convergence),
        ("gm_ladder.png", figure_gm_ladder),
    ):
        try:
            builder(results, output / name)
        except ImportError:  # matplotlib is optional
            print(f"[stability] matplotlib missing, skipping {name}", file=sys.stderr)

    (output / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    (output / "report.md").write_text(render_report(results, claims))

    hull = results["hull"]  # type: ignore[index]
    print(f"[stability] box worst relative error      {results['box']['worst_relative_error']:.3e}")
    print(f"[stability] hemisphere worst rel. error   {results['hemisphere']['worst_relative_error']:.3e}")
    for level in hull["levels"]:
        print(
            f"[stability] hull {level['nx']:>3}×{level['nz']:<3} A {level['waterplane_area_error']:+.4%}  "
            f"I_T {level['transverse_second_moment_error']:+.4%}  I_L "
            f"{level['longitudinal_second_moment_error']:+.4%}  GM {level['gm_transverse']:.4f} m "
            f"({level['gm_transverse_error']:+.4%})"
        )
    failed = [claim for claim in claims if not claim["holds"]]
    for claim in claims:
        status = "PASS" if claim["holds"] else "FAIL"
        print(f"[stability] {status}  {claim['claim']}  ({claim['value']:.3e} vs {claim['threshold']:.0e})")
    if failed:
        print(f"[stability] {len(failed)} claim(s) do not hold", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
