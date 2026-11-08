"""Closed NURBS patch with singular pole rings: the sphere anchor, measured.

Runs the whole acceptance chain of the pole-surface verification and writes every
number and figure it reports to ``artifacts/nurbs_sphere/``, so the report and
``docs/nurbs_pole_surface.md`` never have to copy a value by hand.

    python3 examples/nurbs_sphere/run_example.py [output-directory]

Sequence:

1. exactness -- ``|S(u, v)| == R`` over the whole parameter domain;
2. parametric area -- ``4 pi R^2`` by Gauss-Legendre quadrature per knot span, plus the
   order sweep that shows where the round-off floor is;
3. pole welding -- how far the per-sample evaluation at ``v = v_min`` drifts from the
   pole (in ulps) and what that does to vertex welding if the row is not keyed;
4. tessellation -- watertightness, orientation, per-density volume and area
   convergence against ``4/3 pi R^3`` and ``4 pi R^2`` (``O(h^2)`` reference slope);
5. STL round-trip of the closed shell;
6. two negative controls -- the reversed-pole meridian (watertight but uniformly
   inward, negative volume) and a ``u``-pinched control net.
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

from tensorcad import mesh as meshlib
from tensorcad.export import read_stl, write_stl
from tensorcad.geometry import MeshBuilder
from tensorcad.nurbs import (
    NURBSError,
    NURBSSurface,
    SQRT_HALF,
    circle_curve,
    mesh_solid_report,
    sphere_surface,
    tessellation_report,
)

#: Scale of the demonstration: the PC6 bilge radius B/2 = 12 m, as in the tube anchor.
RADIUS = 12.0
DENSITIES = ((12, 4), (24, 8), (48, 16), (96, 32), (192, 64))
VOLUME_THRESHOLD = 1e-3
AREA_THRESHOLD = 1e-2


def volume_exact(radius: float) -> float:
    return 4.0 / 3.0 * math.pi * radius**3


def area_exact(radius: float) -> float:
    return 4.0 * math.pi * radius**2


def reversed_meridian_sphere(radius: float) -> NURBSSurface:
    """The same sphere with the meridian running north -> south (negative control)."""
    circle = circle_curve(1.0, arcs=3)
    meridian = ((0.0, 1.0), (1.0, 1.0), (1.0, 0.0), (1.0, -1.0), (0.0, -1.0))
    weights = (1.0, SQRT_HALF, 1.0, SQRT_HALF, 1.0)
    net = tuple(
        tuple((radius * r * p[0], radius * r * p[1], radius * z) for r, z in meridian)
        for p in circle.control_points
    )
    combined = tuple(
        tuple(circle.weights[i] * weights[j] for j in range(len(meridian)))
        for i in range(len(circle.control_points))
    )
    return NURBSSurface(
        degree=(2, 2),
        knots=(circle.knots, (0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 2.0)),
        control_net=net,
        weights=combined,
    )


def formation_ulps(values: list[float], reference: float) -> int:
    """How many representable steps separate ``values`` from ``reference``."""
    step = math.ulp(reference) if hasattr(math, "ulp") else abs(reference) * 2.0**-52
    return max(int(round(abs(value - reference) / step)) for value in values)


def main(output_directory: str | Path | None = None) -> dict[str, object]:
    out = Path(output_directory) if output_directory else ROOT / "artifacts" / "nurbs_sphere"
    out.mkdir(parents=True, exist_ok=True)
    surface = sphere_surface(RADIUS)
    v_low, v_high = surface.domain[1]
    exact_volume = volume_exact(RADIUS)
    exact_area = area_exact(RADIUS)

    # 1. exactness -----------------------------------------------------------
    worst_radial = 0.0
    for i in range(17):
        for j in range(17):
            point = surface.evaluate(i / 16.0, v_low + (v_high - v_low) * j / 16.0)
            worst_radial = max(worst_radial, abs(math.sqrt(sum(c * c for c in point)) - RADIUS))

    # 2. parametric area -----------------------------------------------------
    area_orders = {order: surface.area(order) for order in (4, 6, 8, 12, 16, 24, 32)}
    area_errors = {order: value / exact_area - 1.0 for order, value in area_orders.items()}

    # 3. pole welding --------------------------------------------------------
    south = surface.singular_boundary_point("v")
    north = surface.singular_boundary_point("v", upper=True)
    assert south is not None and north is not None
    u_low, u_high = surface.domain[0]
    raw_south = [
        surface.evaluate(u_low + (u_high - u_low) * index / 48, v_low) for index in range(49)
    ]
    raw_north = [
        surface.evaluate(u_low + (u_high - u_low) * index / 48, v_high) for index in range(49)
    ]
    pole_drift_ulps = max(
        formation_ulps([p[2] for p in raw_south], south[2]),
        formation_ulps([p[2] for p in raw_north], north[2]),
    )
    distinct_raw = len({p for p in raw_south})
    # What a naive (non-keyed) tessellator would produce: a ragged pole ring.
    naive = MeshBuilder()
    for index in range(49):
        naive.add_vertex(surface.evaluate(u_low + (u_high - u_low) * index / 48, v_low))
    naive_vertices = naive.vertex_count
    keyed = surface.tessellate(48, 16)
    keyed_poles = [
        point
        for point in keyed.mesh.vertices
        if abs(abs(point[2]) - RADIUS) < 1e-12 and math.hypot(point[0], point[1]) < 1e-12
    ]

    # 4. tessellation convergence -------------------------------------------
    rows: list[dict[str, object]] = []
    for u_samples, v_samples in DENSITIES:
        tess = surface.tessellate(u_samples, v_samples)
        solid = mesh_solid_report(
            tess.mesh, volume_reference=exact_volume, area_reference=exact_area
        )
        audit = tessellation_report(surface, tess)
        rows.append(
            {
                "u_samples": u_samples,
                "v_samples": v_samples,
                "triangles": len(tess.mesh.faces),
                "vertices": len(tess.mesh.vertices),
                "expected_vertices": u_samples * (v_samples - 1) + 2,
                "degenerate_quads_skipped": tess.degenerate_skipped,
                "watertight": solid["watertight"],
                "edge_manifold": solid["edge_manifold"],
                "boundary_edges": solid["boundary_edges"],
                "non_manifold_edges": solid["non_manifold_edges"],
                "components": solid["components"],
                "zero_area_faces": len(meshlib.zero_area_faces(tess.mesh)),
                "signed_volume_m3": solid["signed_volume_m3"],
                "volume_relative_error": solid["volume_relative_error"],
                "tessellated_area_m2": audit["tessellated_area_m2"],
                "area_relative_error": audit["area_relative_error"],
                "max_normal_angle_deg": audit["max_normal_angle_deg"],
                "inward_faces": audit["inward_faces"],
                "orientation_flips_walk_order": audit["orientation_flips"],
                "orientation_consistent": solid["orientation_consistent"],
                "faces_per_quad_aligned": sum(
                    len(quad) for quad in tess.faces_per_quad
                )
                == len(tess.mesh.faces),
            }
        )

    # 5. STL round-trip ------------------------------------------------------
    stl_path = out / "sphere_solid.stl"
    tess_mid = surface.tessellate(48, 16)
    write_stl(tess_mid.mesh, stl_path)
    recovered = read_stl(stl_path)
    stl = {
        "triangles": len(recovered.faces),
        "watertight": meshlib.is_watertight(recovered),
        "volume_relative_error": meshlib.total_signed_volume(recovered) / exact_volume - 1.0,
        "source_volume_relative_error": meshlib.total_signed_volume(tess_mid.mesh)
        / exact_volume
        - 1.0,
    }

    # 6. negative controls ---------------------------------------------------
    reversed_sphere = reversed_meridian_sphere(RADIUS)
    reversed_tess = reversed_sphere.tessellate(24, 8)
    reversed_report = mesh_solid_report(reversed_tess.mesh, volume_reference=exact_volume)
    negative_controls = {
        "reversed_poles_watertight": reversed_report["watertight"],
        "reversed_poles_signed_volume_m3": reversed_report["signed_volume_m3"],
        "reversed_poles_orientation_consistent": reversed_report["orientation_consistent"],
        "reversed_poles_flipped_faces": reversed_report["orientation_flipped_faces"],
        "reversed_poles_triangles": len(reversed_tess.mesh.faces),
        "reversed_poles_radius_error": abs(
            math.sqrt(sum(c * c for c in reversed_sphere.evaluate(0.3, 1.1))) - RADIUS
        ),
        "zero_radius_rejected": False,
        "bad_axis_rejected": False,
    }
    for label, call in (
        ("zero_radius_rejected", lambda: sphere_surface(0.0)),
        ("bad_axis_rejected", lambda: surface.singular_boundary_point("w")),
    ):
        try:
            call()
        except NURBSError:
            negative_controls[label] = True

    criteria = {
        "sphere_exact_within_1e-12": worst_radial < 1e-12,
        "parametric_area_within_1e-12": abs(area_errors[24]) < 1e-12,
        "poles_weld_to_two_vertices": len(keyed_poles) == 2,
        "vertex_count_closed_form_holds": all(
            row["vertices"] == row["expected_vertices"] for row in rows
        ),
        "faces_per_quad_aligned_at_every_density": all(
            row["faces_per_quad_aligned"] for row in rows
        ),
        "no_zero_area_faces_at_every_density": all(
            row["zero_area_faces"] == 0 for row in rows
        ),
        "watertight_at_every_density": all(row["watertight"] for row in rows),
        "outward_at_every_density": all(
            row["signed_volume_m3"] > 0.0 and row["inward_faces"] == 0 for row in rows
        ),
        "volume_converges_o_h2": _converges(rows, "volume_relative_error"),
        "area_converges_o_h2": _converges(rows, "area_relative_error"),
        "volume_error_below_threshold_at_192x64": abs(
            rows[-1]["volume_relative_error"]
        )
        < VOLUME_THRESHOLD,
        "stl_round_trip_watertight": stl["watertight"],
        "reversed_poles_rejected_as_inward": (
            negative_controls["reversed_poles_orientation_consistent"] is False
            and negative_controls["reversed_poles_signed_volume_m3"] < 0.0
        ),
        "bad_inputs_rejected": negative_controls["zero_radius_rejected"]
        and negative_controls["bad_axis_rejected"],
    }

    payload: dict[str, object] = {
        "case": {
            "object": "rational biquadratic sphere, one patch, two singular pole rings",
            "radius_m": RADIUS,
            "degree": list(surface.degree),
            "control_net": [len(surface.control_net), len(surface.control_net[0])],
            "weights": "circle 1, 1/2, 1, ... x meridian 1, sqrt(1/2), 1, sqrt(1/2), 1",
            "poles": {"south": south, "north": north},
            "why": "singular pole rows had no committed anchor (m6410 follow-up)",
        },
        "exactness": {
            "max_radial_error_m": worst_radial,
            "area_m2": surface.area(),
            "area_exact_m2": exact_area,
            "area_relative_error": area_errors[24],
            "volume_exact_m3": exact_volume,
        },
        "area_quadrature": {
            "orders": {str(order): value for order, value in area_orders.items()},
            "relative_errors": {str(order): value for order, value in area_errors.items()},
        },
        "pole_welding": {
            "drift_ulps_before_keying": pole_drift_ulps,
            "distinct_raw_pole_vertices_49_samples": distinct_raw,
            "naive_pole_vertices_49_samples": naive_vertices,
            "keyed_tessellation_pole_vertices": len(keyed_poles),
            "degenerate_quads_skipped_48x16": keyed.degenerate_skipped,
        },
        "tessellation": rows,
        "stl_round_trip": stl,
        "negative_controls": negative_controls,
        "criteria": criteria,
        "all_criteria_pass": all(criteria.values()),
    }
    (out / "metrics.json").write_text(json.dumps(payload, indent=2) + "\n")

    _write_figures(out, payload)
    print(json.dumps(payload["criteria"], indent=2))
    print(f"all criteria pass: {payload['all_criteria_pass']}")
    print(f"artifacts -> {out}")
    return payload


def _converges(rows: list[dict[str, object]], key: str) -> bool:
    """Each doubling of the grid must cut the error by roughly four (O(h^2))."""
    errors = [abs(float(row[key])) for row in rows]
    for coarse, fine in zip(errors, errors[1:]):
        if coarse <= 1e-12:
            return False
        ratio = coarse / fine
        if not 3.0 < ratio < 5.0:
            return False
    return True


def _write_figures(out: Path, payload: dict[str, object]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    except Exception as error:  # pragma: no cover - figures are illustrative only
        print(f"matplotlib unavailable, figures skipped: {error}")
        return

    surface = sphere_surface(RADIUS)
    tess = surface.tessellate(48, 16)

    # fig1: the closed shell with its outward normals around the pole ring.
    figure = plt.figure(figsize=(11.0, 5.2))
    axis = figure.add_subplot(1, 2, 1, projection="3d")
    triangles = [[tess.mesh.vertices[i] for i in face] for face in tess.mesh.faces]
    axis.add_collection3d(
        Poly3DCollection(
            triangles, facecolors="#cfe3f7", edgecolors="#5b7fa6", linewidths=0.12, alpha=0.95
        )
    )
    axis.set_xlim(-RADIUS, RADIUS)
    axis.set_ylim(-RADIUS, RADIUS)
    axis.set_zlim(-RADIUS, RADIUS)
    axis.set_xlabel("x [m]")
    axis.set_ylabel("y [m]")
    axis.set_zlabel("z [m]")
    axis.set_title("one closed NURBS patch, 48 x 16, watertight")
    axis.view_init(elev=18.0, azim=32.0)

    normals = figure.add_subplot(1, 2, 2)
    sample = []
    for face in tess.mesh.faces:
        centre = meshlib.face_centroid(tess.mesh, face)
        if abs(centre[2]) > 0.93 * RADIUS:
            sample.append((centre, meshlib.face_normal(tess.mesh, face)))
    for centre, raw in sample[:: max(1, len(sample) // 60)]:
        length = math.sqrt(sum(component * component for component in raw)) or 1.0
        unit = [component / length for component in raw]
        normals.quiver(
            centre[0],
            centre[1],
            unit[0],
            unit[1],
            color="#b3541e",
            angles="xy",
            scale_units="xy",
            scale=8.0,
            width=0.004,
        )
    ring = [tess.mesh.vertices[i] for i in range(len(tess.mesh.vertices))]
    normals.scatter(
        [p[0] for p in ring if abs(p[2]) > 0.93 * RADIUS],
        [p[1] for p in ring if abs(p[2]) > 0.93 * RADIUS],
        s=9.0,
        color="#1f3f5b",
    )
    normals.set_aspect("equal")
    normals.set_title("outward normals on both pole rings")
    normals.set_xlabel("x [m]")
    normals.set_ylabel("y [m]")
    normals.grid(alpha=0.3)
    figure.tight_layout()
    figure.savefig(out / "fig1_pole_shell.png", dpi=150)
    plt.close(figure)

    # fig2: convergence in volume and area, plus the quadrature round-off floor.
    rows = payload["tessellation"]
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.4))
    h = [1.0 / row["u_samples"] for row in rows]
    left.loglog(h, [abs(row["volume_relative_error"]) for row in rows], "o-", label="volume")
    left.loglog(h, [abs(row["area_relative_error"]) for row in rows], "s-", label="area")
    reference = [abs(rows[0]["volume_relative_error"]) * (value / h[0]) ** 2 for value in h]
    left.loglog(h, reference, "k--", linewidth=1.0, label="O(h^2) reference")
    left.set_xlabel("1 / u_samples")
    left.set_ylabel("relative error")
    left.set_title("tessellation converges as O(h^2)")
    left.grid(alpha=0.3, which="both")
    left.legend()

    orders = payload["area_quadrature"]["relative_errors"]
    axis_orders = [int(order) for order in orders]
    right.semilogy(
        axis_orders,
        [max(abs(value), 1e-17) for value in orders.values()],
        "o-",
        label="|area error| vs 4 pi R^2",
    )
    right.axhline(2.0**-52, color="k", linestyle=":", label="double round-off 2^-52")
    right.set_xlabel("Gauss points per knot span")
    right.set_ylabel("relative error")
    right.set_title("parametric area: quadrature floor")
    right.grid(alpha=0.3, which="both")
    right.legend()
    figure.tight_layout()
    figure.savefig(out / "fig2_convergence.png", dpi=150)
    plt.close(figure)

    # fig3: why the pole row must be keyed, and what the pole order does to the sign.
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    v_low = surface.domain[1][0]
    u_low, u_high = surface.domain[0]
    step = math.ulp(RADIUS)
    raw = [
        surface.evaluate(u_low + (u_high - u_low) * index / 96, v_low)[2] for index in range(97)
    ]
    left.plot(
        range(97),
        [(value - RADIUS * -1.0) / step for value in raw],
        "o-",
        markersize=3.0,
        label="raw S(u, v_min) z-offset [ulps]",
    )
    left.set_xlabel("u sample index")
    left.set_ylabel("offset from the south pole [ulps]")
    left.set_title("the pole ring does not evaluate bit-for-bit")
    left.grid(alpha=0.3)
    left.legend()

    reversed_report = payload["negative_controls"]
    labels = ["poles south$\\to$north", "poles north$\\to$south"]
    values = [float(payload["tessellation"][1]["signed_volume_m3"]), float(
        reversed_report["reversed_poles_signed_volume_m3"]
    )]
    right.bar(labels, values, color=["#2b6cb0", "#b3541e"])
    right.axhline(0.0, color="k", linewidth=0.8)
    right.set_ylabel("signed volume [m$^3$]")
    right.set_title("watertight in both cases; only the sign changes")
    right.grid(alpha=0.3, axis="y")
    figure.tight_layout()
    figure.savefig(out / "fig3_pole_physics.png", dpi=150)
    plt.close(figure)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)