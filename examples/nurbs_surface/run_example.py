"""NURBS surface service: analytic anchors, continuity, tessellation, STL export.

Runs the whole acceptance chain of campaign node IN16-K5 and writes every number and
figure it reports to ``artifacts/nurbs_surface/`` so the report body never has to
copy a value by hand.

    python3 examples/nurbs_surface/run_example.py [output-directory]

Sequence:

1. exact circular arc and full circle -> length against ``r * theta``;
2. exact rational cylinder -> parametric area against ``2 pi r h`` and the
   Gauss-Legendre convergence curve;
3. tessellation -> watertightness, normal consistency, area convergence;
4. closed solid -> volume against ``pi r^2 h``, STL round-trip;
5. parameter-derivative continuity across the interior knots (declared vs measured);
6. the campaign's negative controls, each asserted to raise;
7. the additive PC6 contract block, with the L/B/Cb geometry gate.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    # The ambient PYTHONPATH may carry another checkout; this example must measure the
    # code it ships with.
    sys.path.insert(0, str(ROOT))

from tensorcad.export import read_stl, write_stl
from tensorcad import mesh as meshlib
from tensorcad.nurbs import (
    NURBSCurve,
    NURBSError,
    NURBS_SCHEMA,
    circle_curve,
    circular_arc,
    cylinder_solid_mesh,
    cylinder_surface,
    pc6_compat_block,
    tessellation_report,
    mesh_solid_report,
)
from tensorcad.pc6 import PC6Parameters, analytic_form
from tensorcad.pc6_contract import CONTRACT_SCHEMA

RADIUS = 12.0
HEIGHT = 11.0


def _relative(measured: float, reference: float) -> float:
    return (measured - reference) / reference


def _negative_controls() -> list[dict[str, object]]:
    """The campaign's negative controls: each must be rejected, with a reason."""
    cases: list[dict[str, object]] = []

    def record(name: str, expectation: str, action) -> None:
        try:
            action()
        except NURBSError as error:
            cases.append({"case": name, "expected": expectation, "rejected": True,
                          "message": str(error)})
        else:
            cases.append({"case": name, "expected": expectation, "rejected": False,
                          "message": "NOT REJECTED - the guard is missing"})

    good = circle_curve(1.0)
    record(
        "non_monotone_knots",
        "knot vector must be non-decreasing",
        lambda: NURBSCurve(
            degree=2,
            knots=(0.0, 0.0, 0.0, 0.6, 0.4, 1.0, 1.0),
            control_points=((0, 0, 0), (1, 1, 0), (2, 0, 0), (3, 1, 0)),
            weights=(1.0, 1.0, 1.0, 1.0),
        ),
    )
    record(
        "zero_weight",
        "weights must be strictly positive",
        lambda: NURBSCurve(
            degree=2,
            knots=good.knots,
            control_points=good.control_points,
            weights=(1.0,) + (0.0,) + (1.0,) * (len(good.weights) - 2),
        ),
    )
    record(
        "negative_weight",
        "weights must be strictly positive",
        lambda: NURBSCurve(
            degree=2,
            knots=good.knots,
            control_points=good.control_points,
            weights=(-0.5,) + (0.5,) * (len(good.weights) - 1),
        ),
    )
    record(
        "parameter_outside_domain",
        "parameter outside [U[p], U[n+1]] is an error, not a silent clamp",
        lambda: good.evaluate(1.5),
    )
    record(
        "unclamped_end",
        "a non-clamped end would miss its end control point",
        lambda: NURBSCurve(
            degree=2,
            knots=(-0.5, 0.0, 0.0, 1.0, 1.0, 1.0),
            control_points=((0, 0, 0), (1, 1, 0), (2, 0, 0)),
            weights=(1.0, 0.5, 1.0),
        ),
    )
    record(
        "knot_multiplicity_too_high",
        "multiplicity above degree + 1",
        lambda: NURBSCurve(
            degree=2,
            knots=(0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0),
            control_points=((0, 0, 0), (1, 1, 0), (2, 0, 0)),
            weights=(1.0, 0.5, 1.0),
        ),
    )
    return cases


def _continuity_section() -> dict[str, object]:
    circle = circle_curve(RADIUS)
    circle_report = circle.continuity(2)
    # A cubic polynomial B-spline with one interior knot of multiplicity 1: declared
    # C2, so the third derivative is the first one that must jump. This is the positive
    # control for the jump measurement itself.
    cubic = NURBSCurve(
        degree=3,
        knots=(0.0, 0.0, 0.0, 0.0, 0.4, 1.0, 1.0, 1.0, 1.0),
        control_points=((0, 0, 0), (1, 2, 0), (2, 0, 0), (3, 3, 0), (4, 1, 0)),
        weights=(1.0, 1.0, 1.0, 1.0, 1.0),
    )
    cubic_report = cubic.continuity(3)
    surface = cylinder_surface(RADIUS, HEIGHT)
    surface_report = surface.continuity(1)
    return {
        "circle": {
            "declared_order": circle_report["declared_minimum_order"],
            "knots": circle_report["knots"],
            "measured_minimum_order": min(
                record["measured_continuous_order"] for record in circle_report["knots"]
            ),
            "max_first_derivative_jump": max(
                record["jumps"][1] for record in circle_report["knots"]
            ),
        },
        "cubic_positive_control": {
            "declared_order": cubic_report["declared_minimum_order"],
            "knots": cubic_report["knots"],
        },
        "cylinder_surface": {
            "u": surface_report["u"],
            "v": surface_report["v"],
            "passed": surface_report["passed"],
        },
    }


def _state(mesh) -> dict[str, object]:
    return {
        "triangles": len(mesh.faces),
        "vertices": len(mesh.vertices),
        "watertight": meshlib.is_watertight(mesh),
    }


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    outdir = Path(argv[0]) if argv else Path("artifacts/nurbs_surface")
    outdir.mkdir(parents=True, exist_ok=True)

    # -- 1. exact arcs -----------------------------------------------------
    arc = circular_arc(2.0, 1.0)
    arc_samples = arc.polyline(4001)
    arc_area_deviation = max(
        abs(math.hypot(point[0], point[1]) - 2.0) for point in arc_samples
    )
    circle = circle_curve(RADIUS)
    circle_deviation = max(
        abs(math.hypot(point[0], point[1]) - RADIUS) for point in circle.polyline(4001)
    )
    arcs = {
        "radius": 2.0,
        "sweep_rad": 1.0,
        "arc_length": arc.arc_length(),
        "arc_length_reference": 2.0,
        "arc_length_relative_error": _relative(arc.arc_length(), 2.0),
        "max_radius_deviation": arc_area_deviation,
        "circle_radius": RADIUS,
        "circle_length": circle.arc_length(),
        "circle_length_reference": 2.0 * math.pi * RADIUS,
        "circle_length_relative_error": _relative(circle.arc_length(), 2.0 * math.pi * RADIUS),
        "circle_max_radius_deviation": circle_deviation,
        "circle_control_points": len(circle.control_points),
        "circle_weights": list(circle.weights),
    }

    # -- 2. cylinder parametric area --------------------------------------
    surface = cylinder_surface(RADIUS, HEIGHT)
    area = surface.area()
    area_reference = 2.0 * math.pi * RADIUS * HEIGHT
    unit_surface = cylinder_surface(1.0, 1.0)
    area_rows = []
    for radius, height in ((RADIUS, HEIGHT), (1.0, 1.0), (0.5, 40.0)):
        candidate = cylinder_surface(radius, height)
        measured = candidate.area()
        reference = 2.0 * math.pi * radius * height
        area_rows.append(
            {
                "radius_m": radius,
                "height_m": height,
                "parametric_area_m2": measured,
                "analytic_area_m2": reference,
                "relative_error": _relative(measured, reference),
            }
        )
    quadrature = surface.area_convergence()

    # -- 3. tessellation: topology, normals, area --------------------------
    tessellation_rows = []
    for u_samples in (12, 24, 48, 96, 192):
        tessellation = surface.tessellate(u_samples, 8)
        report = tessellation_report(surface, tessellation)
        tessellation_rows.append(
            {
                "u_samples": u_samples,
                "v_samples": 8,
                "triangles": report["triangles"],
                "watertight": report["watertight"],
                "boundary_edges": report["boundary_edges"],
                "area_relative_error": report["area_relative_error"],
                "max_normal_angle_deg": report["max_normal_angle_deg"],
                "mean_normal_angle_deg": report["mean_normal_angle_deg"],
                "inward_faces": report["inward_faces"],
                "orientation_flips": report["orientation_flips"],
            }
        )

    # -- 4. closed solid, STL round-trip ----------------------------------
    solid_rows = []
    volume_reference = math.pi * RADIUS * RADIUS * HEIGHT
    lateral_reference = area_reference
    cap_reference = 2.0 * math.pi * RADIUS * RADIUS
    stl_path = outdir / "cylinder_solid.stl"
    for u_samples in (96, 192, 384):
        solid = cylinder_solid_mesh(RADIUS, HEIGHT, u_samples=u_samples, v_samples=6)
        report = mesh_solid_report(
            solid,
            volume_reference=volume_reference,
            area_reference=lateral_reference + cap_reference,
        )
        solid_rows.append(
            {
                "u_samples": u_samples,
                "triangles": report["triangles"],
                "watertight": report["watertight"],
                "components": report["components"],
                "volume_m3": report["signed_volume_m3"],
                "volume_relative_error": report["volume_relative_error"],
                "area_relative_error": report["area_relative_error"],
                "orientation_consistent": report["orientation_consistent"],
            }
        )
        if u_samples == 384:
            write_stl(solid, stl_path)
    recovered = read_stl(stl_path)
    stl_report = mesh_solid_report(recovered, volume_reference=volume_reference)
    stl_round_trip = {
        "path": str(stl_path),
        "watertight": stl_report["watertight"],
        "components": stl_report["components"],
        "volume_m3": stl_report["signed_volume_m3"],
        "volume_relative_error": stl_report["volume_relative_error"],
        "note": (
            "binary STL stores float32, so the recovered volume cannot be bit-identical; "
            "the shell is re-welded by exact position and stays closed"
        ),
    }

    # -- 5. continuity -----------------------------------------------------
    continuity = _continuity_section()

    # -- 6. negative controls ---------------------------------------------
    negative = _negative_controls()

    # -- 7. PC6 contract compatibility ------------------------------------
    parameters = PC6Parameters()
    form = analytic_form(parameters)
    block = pc6_compat_block(surface, parameters=parameters, tessellation=surface.tessellate(48, 6))
    block_path = outdir / "pc6_nurbs_block.json"
    block_path.write_text(json.dumps(block, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    tamper = None
    try:
        pc6_compat_block(
            surface,
            parameters=parameters,
            declared={"block_coefficient": form.block_coefficient * 0.9},
        )
    except NURBSError as error:
        tamper = str(error)
    compatibility = {
        "contract_schema": CONTRACT_SCHEMA,
        "block_schema": block["schema"],
        "block_representation": block["representation"],
        "additive": block["additive"],
        "geometry_gate": block["geometry_gate"],
        "declared_cb_tamper_rejected": tamper is not None,
        "declared_cb_tamper_message": tamper,
        "block_path": str(block_path),
    }

    metrics = {
        "schema": "tensorcad.nurbs.report/1.0",
        "representation_schema": NURBS_SCHEMA,
        "arcs": arcs,
        "cylinder": {
            "radius_m": RADIUS,
            "height_m": HEIGHT,
            "parametric_area_m2": area,
            "analytic_area_m2": area_reference,
            "area_relative_error": _relative(area, area_reference),
            "unit_cylinder_area_m2": unit_surface.area(),
            "area_rows": area_rows,
            "quadrature_convergence": quadrature,
            "reference_gl_order": max(row["order"] for row in quadrature["rows"]),
        },
        "tessellation": tessellation_rows,
        "solid": {
            "volume_reference_m3": volume_reference,
            "rows": solid_rows,
            "stl_round_trip": stl_round_trip,
        },
        "continuity": continuity,
        "negative_controls": negative,
        "pc6_compatibility": compatibility,
        "criteria": {
            "cylinder_area_error_below_1e-6": abs(_relative(area, area_reference)) < 1e-6,
            "tessellation_max_normal_angle_deg_below_1": all(
                row["max_normal_angle_deg"] < 1.0 for row in tessellation_rows
            ),
            "tessellation_no_inward_faces": all(
                row["inward_faces"] == 0 for row in tessellation_rows
            ),
            "solid_watertight": all(row["watertight"] for row in solid_rows),
            "solid_volume_error_below_1e-3": all(
                abs(row["volume_relative_error"]) < 1e-3 for row in solid_rows
            ),
            "stl_round_trip_watertight": stl_round_trip["watertight"],
            "continuity_declared_holds": continuity["cylinder_surface"]["passed"],
            "all_negative_controls_rejected": all(case["rejected"] for case in negative),
            "pc6_gate_passed": compatibility["geometry_gate"]["passed"],
            "pc6_tamper_rejected": compatibility["declared_cb_tamper_rejected"],
        },
    }
    metrics_path = outdir / "metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    _figures(surface, metrics, outdir)

    print(json.dumps(metrics["criteria"], indent=2, sort_keys=True))
    print(f"metrics   -> {metrics_path}")
    print(f"pc6 block -> {block_path}")
    print(f"stl       -> {stl_path}")
    return 0 if all(metrics["criteria"].values()) else 1


def _figures(surface, metrics: dict, outdir: Path) -> None:
    """Four figures: geometry, convergence, continuity, exactness."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as error:  # pragma: no cover - plotting is optional evidence
        print(f"matplotlib unavailable ({error}); figures skipped")
        return

    # figure 1: the tessellated cylinder with its outward normals
    tessellation = surface.tessellate(36, 6)
    mesh = tessellation.mesh
    fig = plt.figure(figsize=(7.2, 5.4))
    axes = fig.add_subplot(111, projection="3d")
    for face in mesh.faces:
        ring = [mesh.vertices[index] for index in face] + [mesh.vertices[face[0]]]
        axes.plot(
            [point[0] for point in ring],
            [point[1] for point in ring],
            [point[2] for point in ring],
            color="#4477aa",
            linewidth=0.25,
        )
    step = max(len(mesh.faces) // 60, 1)
    for face_index in range(0, len(mesh.faces), step):
        face = mesh.faces[face_index]
        centroid = [sum(mesh.vertices[index][axis] for index in face) / 3.0 for axis in range(3)]
        normal = meshlib.face_normal(mesh, face)
        length = 2.2
        axes.quiver(
            centroid[0], centroid[1], centroid[2],
            normal[0], normal[1], normal[2],
            length=length, color="#cc3311", linewidth=0.8,
        )
    axes.set_xlabel("x [m]")
    axes.set_ylabel("y [m]")
    axes.set_zlabel("z [m]")
    axes.set_title(
        f"NURBS cylinder R={metrics['cylinder']['radius_m']:g} m, "
        f"H={metrics['cylinder']['height_m']:g} m\n"
        f"triangles={len(mesh.faces)}, inward faces="
        f"{metrics['tessellation'][2]['inward_faces']}"
    )
    fig.tight_layout()
    fig.savefig(outdir / "fig1_cylinder_normals.png", dpi=150)
    plt.close(fig)

    # figure 2: quadrature and tessellation convergence
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.4, 4.2))
    rows = metrics["cylinder"]["quadrature_convergence"]["rows"]
    left.semilogy([row["order"] for row in rows],
                  [max(row["relative_error"], 1e-17) for row in rows],
                  "o-", color="#4477aa")
    left.set_xlabel("Gauss-Legendre order per knot span [-]")
    left.set_ylabel("|A_quad - A_ref| / A_ref [-]")
    left.set_title(f"Parametric area: exact to round-off\n"
                   f"A = {metrics['cylinder']['parametric_area_m2']:.6f} m$^2$ vs "
                   f"$2\\pi r h$ = {metrics['cylinder']['analytic_area_m2']:.6f} m$^2$")
    left.grid(True, which="both", alpha=0.3)

    tessellation_rows = metrics["tessellation"]
    inverse = [1.0 / row["u_samples"] for row in tessellation_rows]
    error = [abs(row["area_relative_error"]) for row in tessellation_rows]
    right.loglog(inverse, error, "s-", label="tessellated vs parametric area",
                 color="#cc3311")
    reference = [error[0] * (value / inverse[0]) ** 2 for value in inverse]
    right.loglog(inverse, reference, "k--", label="slope -2 (inscribed facets)")
    right.set_xlabel("1 / u samples [-]")
    right.set_ylabel("relative area error [-]")
    right.set_title("Tessellation converges from below as O(h$^2$)")
    right.legend()
    right.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "fig2_convergence.png", dpi=150)
    plt.close(fig)

    # figure 3: derivative continuity, declared vs measured
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.4, 4.2))
    circle_rows = metrics["continuity"]["circle"]["knots"]
    orders = list(range(len(circle_rows[0]["jumps"])))
    width = 0.35
    for offset, record in enumerate(circle_rows):
        left.bar([order + offset * width for order in orders],
                 [max(jump, 1e-18) for jump in record["jumps"]],
                 width=width,
                 label=f"circle knot {record['knot']:.3f}")
    cubic = metrics["continuity"]["cubic_positive_control"]["knots"][0]
    base = len(orders) + 0.6
    left.bar([base + order for order in range(len(cubic["jumps"]))],
             [max(jump, 1e-18) for jump in cubic["jumps"]],
             color="#228833", label="cubic knot 0.4 (C2 control)")
    left.set_yscale("log")
    left.set_xticks(list(range(len(orders))) + [base + order for order in range(len(cubic["jumps"]))])
    left.set_xticklabels([f"arc\nk={order}" for order in orders]
                         + [f"cubic\nk={order}" for order in range(len(cubic["jumps"]))])
    left.set_ylabel("|C$^{(k)}$(t+0) - C$^{(k)}$(t-0)| [m$^{1-k}$]")
    left.set_title("One-sided derivative jumps: measured, not assumed")
    left.legend(fontsize=7)
    left.grid(True, axis="y", which="both", alpha=0.3)

    curve = circle_curve(metrics["arcs"]["circle_radius"])
    samples = [index / 400 for index in range(401)]
    speeds = [math.sqrt(sum(value * value for value in curve.derivatives(u, 1)[1]))
              for u in samples]
    right.plot(samples, speeds, color="#4477aa")
    for knot in curve.interior_knots:
        right.axvline(knot, color="#cc3311", linestyle=":", linewidth=1.0)
    right.set_xlabel("parameter u [-]")
    right.set_ylabel("|C'(u)| [m]")
    right.set_title("Three-arc circle: |C'| is symmetric at the joins,\nhence C1 despite multiplicity 2")
    right.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "fig3_continuity.png", dpi=150)
    plt.close(fig)

    # figure 4: exactness of the arc as a map (radius deviation over the domain)
    fig = plt.figure(figsize=(7.2, 4.6))
    axes = fig.add_subplot(111)
    u_values = [index / 60 for index in range(61)]
    v_values = [index / 8 for index in range(9)]
    deviations = []
    for v in v_values:
        row = []
        for u in u_values:
            x, y, _ = surface.evaluate(u, v)
            row.append(math.hypot(x, y) - metrics["cylinder"]["radius_m"])
        deviations.append(row)
    image = axes.pcolormesh(u_values, v_values, deviations, shading="auto", cmap="viridis")
    bar = fig.colorbar(image, ax=axes)
    bar.set_label("sqrt(x²+y²) - r [m]")
    axes.set_xlabel("u (exact circle direction)")
    axes.set_ylabel("v (linear height)")
    axes.set_title(
        "Deviation from the analytic cylinder: max "
        f"{max(abs(value) for row in deviations for value in row):.2e} m "
        "(round-off, not discretisation)"
    )
    fig.tight_layout()
    fig.savefig(outdir / "fig4_exactness.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
