"""Evidence behind ``docs/patch_network.md``: the multi-patch stitcher, measured.

Nine measurements, every number in the document produced here:

1. **Box** (six bilinear patches, 2 x 3 x 4 m): the tessellation is *exact*, so the
   closed shell's ``V = abc`` and ``A = 2(ab + bc + ca)`` to round-off -- any error
   would be a seam error, not a discretisation error.
2. **Torus** (16 rational biquadratic patches, R = 2, r = 0.5): ``V = 2 pi^2 R r^2``,
   ``A = 4 pi^2 R r``, refinement 8/16/32 and the measured convergence order.
3. **Double cone** (8 patches, two apex poles collapsed): ``V = 2 pi R^2 H / 3``,
   ``A = 2 pi R sqrt(R^2 + H^2)``, same refinement study.
4. **Seams**: per-seam positional gap and dihedral angle for the torus (all G1) and
   the box (all G0 at 90 degrees).
5. **Negative controls**: a patch moved by 1e-3 m (kept topology, lost geometry), a
   missing patch (a real hole), an unclaimed but bit-identical seam (closed by
   position whether or not it was declared), a wrong ``reversed`` flag (gap jumps to
   the edge length), and a sample-count mismatch (refused).

Run::

    PYTHONPATH=$PWD python examples/patch_network/run_example.py
    # -> examples/patch_network/outputs/results.json
    # -> examples/patch_network/outputs/report.md
    # -> examples/patch_network/outputs/figures/*.png
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path

_REPO_ROOT = _Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in _sys.path:  # run from a checkout without installing
    _sys.path.insert(0, str(_REPO_ROOT))

import json
import math
from pathlib import Path

from tensorcad import mesh as meshlib
from tensorcad.nurbs import NURBSCurve, circular_arc
from tensorcad.patches import (
    PatchNetwork,
    PatchNetworkError,
    Seam,
    bilinear_patch,
    derive_seams,
    revolve_patch,
    seam_network,
    stitch,
)

OUTPUTS = Path(__file__).resolve().parent / "outputs"
FIGURES = OUTPUTS / "figures"
BOX = (2.0, 3.0, 4.0)
TORUS = (2.0, 0.5)
CONE = (1.0, 1.0)


# ---------------------------------------------------------------------------
# networks
# ---------------------------------------------------------------------------
def box_network(a: float = BOX[0], b: float = BOX[1], c: float = BOX[2]) -> PatchNetwork:
    faces = (
        ((0.0, 0.0, 0.0), (0.0, b, 0.0), (a, 0.0, 0.0)),
        ((0.0, 0.0, c), (a, 0.0, 0.0), (0.0, b, 0.0)),
        ((0.0, 0.0, 0.0), (a, 0.0, 0.0), (0.0, 0.0, c)),
        ((0.0, b, 0.0), (0.0, 0.0, c), (a, 0.0, 0.0)),
        ((0.0, 0.0, 0.0), (0.0, 0.0, c), (0.0, b, 0.0)),
        ((a, 0.0, 0.0), (0.0, b, 0.0), (0.0, 0.0, c)),
    )
    surfaces = tuple(bilinear_patch(origin, du, dv) for origin, du, dv in faces)
    return seam_network(surfaces, derive_seams(surfaces))


def meridian_arc(major: float, minor: float, start: float, sweep: float) -> NURBSCurve:
    arc = circular_arc(minor, sweep, start_angle=start)
    return NURBSCurve(
        degree=arc.degree,
        knots=arc.knots,
        control_points=tuple((major + p[0], p[1], 0.0) for p in arc.control_points),
        weights=arc.weights,
    )


def torus_network(major: float = TORUS[0], minor: float = TORUS[1]) -> PatchNetwork:
    profiles = [
        meridian_arc(major, minor, -0.5 * math.pi + k * 0.5 * math.pi, 0.5 * math.pi)
        for k in range(4)
    ]
    surfaces = tuple(
        revolve_patch(profiles[tube], quadrant * 0.5 * math.pi, 0.5 * math.pi)
        for quadrant in range(4)
        for tube in range(4)
    )
    return seam_network(surfaces, derive_seams(surfaces))


def double_cone_network(radius: float = CONE[0], half_height: float = CONE[1]) -> PatchNetwork:
    up = NURBSCurve(1, (0.0, 0.0, 1.0, 1.0), ((radius, 0.0, 0.0), (0.0, half_height, 0.0)))
    down = NURBSCurve(1, (0.0, 0.0, 1.0, 1.0), ((0.0, -half_height, 0.0), (radius, 0.0, 0.0)))
    surfaces = tuple(
        revolve_patch(profile, quadrant * 0.5 * math.pi, 0.5 * math.pi)
        for profile in (up, down)
        for quadrant in range(4)
    )
    return seam_network(surfaces, derive_seams(surfaces))


def torus_volume(major: float = TORUS[0], minor: float = TORUS[1]) -> float:
    return 2.0 * math.pi**2 * major * minor**2


def torus_area(major: float = TORUS[0], minor: float = TORUS[1]) -> float:
    return 4.0 * math.pi**2 * major * minor


def cone_volume(radius: float = CONE[0], half: float = CONE[1]) -> float:
    return 2.0 * math.pi * radius**2 * half / 3.0


def cone_area(radius: float = CONE[0], half: float = CONE[1]) -> float:
    return 2.0 * math.pi * radius * math.hypot(radius, half)


# ---------------------------------------------------------------------------
# measurements
# ---------------------------------------------------------------------------
def refine(network: PatchNetwork, volume: float, area: float,
           samples: tuple[int, ...]) -> list[dict]:
    rows = []
    for n in samples:
        stitched = stitch(network, nu=n, nv=n)
        report = stitched.solid_report(reference_volume=volume, reference_area=area)
        rows.append(
            {
                "samples_per_patch": n,
                "triangles": report["triangles"],
                "watertight": report["watertight"],
                "edge_manifold": report["edge_manifold"],
                "components": report["components"],
                "volume_m3": report["signed_volume_m3"],
                "volume_error": report["volume_relative_error"],
                "area_m2": report["surface_area_m2"],
                "area_error": report["area_relative_error"],
                "seam_gap_max_m": max(r.gap_max for r in stitched.seams),
                "seam_angle_max_deg": max(r.angle_max_deg for r in stitched.seams),
            }
        )
    return rows


def observed_order(rows: list[dict], key: str) -> list[float]:
    orders = []
    for coarse, fine in zip(rows, rows[1:]):
        ratio = abs(fine[key]) / abs(coarse[key]) if coarse[key] else float("inf")
        orders.append(math.log2(1.0 / ratio))
    return orders


def negative_controls() -> dict:
    controls: dict[str, object] = {}

    # (a) one face moved by 1e-3 m: the derived topology refuses it, the declared
    #     topology welds it and keeps the residual gap in the report
    delta = 1e-3
    base = box_network()
    surfaces = list(base.surfaces)
    surfaces[5] = bilinear_patch((BOX[0] + delta, 0.0, 0.0), (0.0, BOX[1], 0.0), (0.0, 0.0, BOX[2]))
    derived = PatchNetwork(tuple(surfaces), derive_seams(tuple(surfaces)))
    declared = stitch(PatchNetwork(tuple(surfaces), base.seams), nu=4, nv=4)
    controls["moved_patch"] = {
        "offset_m": delta,
        "sides_left_unclaimed_by_derive": len(derived.unclaimed_sides()),
        "watertight_with_declared_seams": meshlib.is_watertight(declared.mesh),
        "measured_gap_m": max(r.gap_max for r in declared.seams),
        "volume_error": declared.solid_report(reference_volume=BOX[0] * BOX[1] * BOX[2])[
            "volume_relative_error"
        ],
    }

    # (b) a missing patch is a real hole
    network = box_network()
    holed = PatchNetwork(
        network.surfaces[:5],
        tuple(s for s in network.seams if 5 not in (s.first[0], s.second[0])),
    )
    stitched = stitch(holed, nu=4, nv=4)
    controls["missing_patch"] = {
        "open_sides": len(stitched.open_sides),
        "watertight": meshlib.is_watertight(stitched.mesh),
        "boundary_edges": stitched.solid_report()["boundary_edges"],
    }

    # (c) an unclaimed seam whose sides coincide bitwise is welded by position anyway
    unclaimed = stitch(PatchNetwork(network.surfaces, network.seams[:-1]), nu=4, nv=4)
    controls["unclaimed_but_coincident_seam"] = {
        "open_sides": len(unclaimed.open_sides),
        "watertight": meshlib.is_watertight(unclaimed.mesh),
    }

    # (d) a wrong direction flag takes the ends away from each other
    first = bilinear_patch((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    second = bilinear_patch((1.0, 1.0, 0.0), (0.0, -1.0, 0.0), (1.0, 0.0, 0.0))
    right = stitch(PatchNetwork((first, second), derive_seams((first, second))), nu=4, nv=4)
    wrong = stitch(PatchNetwork((first, second), (Seam((0, "u1"), (1, "v0")),)), nu=4, nv=4)
    controls["reversed_flag"] = {
        "gap_declared_correctly_m": right.seams[0].gap_max,
        "gap_flag_omitted_m": wrong.seams[0].gap_max,
        "edge_length_m": 1.0,
    }

    # (e) a u-run seam against a v-run seam needs matching grids
    try:
        stitch(PatchNetwork((first, second), (Seam((0, "u1"), (1, "v0")),)), nu=4, nv=8)
        controls["sample_mismatch"] = {"refused": False, "message": ""}
    except PatchNetworkError as error:
        controls["sample_mismatch"] = {"refused": True, "message": str(error)}
    return controls


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def plot_convergence(torus_rows: list[dict], cone_rows: list[dict]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGURES.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))
    for rows, label, marker in ((torus_rows, "torus (16 patches)", "o"), (cone_rows, "double cone (8 patches)", "s")):
        widths = [row["samples_per_patch"] for row in rows]
        axes[0].loglog(widths, [abs(row["volume_error"]) for row in rows], marker + "-", label=label)
        axes[1].loglog(widths, [abs(row["area_error"]) for row in rows], marker + "-", label=label)
    # reference slopes: h^2 with h = 1/samples
    widths = [row["samples_per_patch"] for row in torus_rows]
    anchor = abs(torus_rows[0]["volume_error"])
    axes[0].loglog(widths, [anchor * (widths[0] / w) ** 2 for w in widths], "k--", lw=1.0, label="$h^2$ reference")
    axes[0].set_xlabel("samples per patch direction")
    axes[0].set_ylabel("|volume relative error|")
    axes[0].set_title("enclosed volume vs analytic")
    axes[0].grid(True, which="both", alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[1].set_xlabel("samples per patch direction")
    axes[1].set_ylabel("|area relative error|")
    axes[1].set_title("surface area vs analytic")
    axes[1].grid(True, which="both", alpha=0.3)
    axes[1].legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(FIGURES / "convergence.png", dpi=140)
    plt.close(figure)


def plot_seams(torus_records, box_records) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGURES.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.0))
    axes[0].bar(range(len(torus_records)), [r.angle_max_deg for r in torus_records], color="#1f77b4")
    axes[0].set_xlabel("torus seam index")
    axes[0].set_ylabel("dihedral angle across the seam [deg]")
    axes[0].set_title("16-patch torus: every seam is G1")
    axes[0].set_ylim(0.0, max(1e-5, 1.2 * max(r.angle_max_deg for r in torus_records)))
    axes[1].bar(range(len(box_records)), [r.angle_max_deg for r in box_records], color="#ff7f0e")
    axes[1].set_xlabel("box seam index")
    axes[1].set_ylabel("dihedral angle across the seam [deg]")
    axes[1].set_title("6-patch box: every seam is G0 at 90 deg")
    axes[1].set_ylim(0.0, 100.0)
    figure.tight_layout()
    figure.savefig(FIGURES / "seams.png", dpi=140)
    plt.close(figure)


def plot_shell(stitched, name: str, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    FIGURES.mkdir(parents=True, exist_ok=True)
    mesh = stitched.mesh
    triangles = [
        [mesh.vertices[index] for index in face] for face in mesh.faces
    ]
    figure = plt.figure(figsize=(6.4, 5.2))
    axes = figure.add_subplot(111, projection="3d")
    collection = Poly3DCollection(triangles, linewidths=0.05, edgecolors="#333333",
                                  facecolors="#7fb3d5", alpha=0.95)
    axes.add_collection3d(collection)
    points = mesh.vertices
    axes.auto_set_scale = True
    axes.set_xlim(min(p[0] for p in points), max(p[0] for p in points))
    axes.set_ylim(min(p[1] for p in points), max(p[1] for p in points))
    axes.set_zlim(min(p[2] for p in points), max(p[2] for p in points))
    axes.set_xlabel("x [m]")
    axes.set_ylabel("y [m]")
    axes.set_zlabel("z [m]")
    axes.set_title(title)
    axes.view_init(elev=22.0, azim=35.0)
    figure.tight_layout()
    figure.savefig(FIGURES / f"shell_{name}.png", dpi=140)
    plt.close(figure)


# ---------------------------------------------------------------------------
def main() -> int:
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    a, b, c = BOX
    box_stitched = stitch(box_network(), nu=4, nv=4)
    box_report = box_stitched.solid_report(
        reference_volume=a * b * c, reference_area=2 * (a * b + b * c + c * a)
    )

    torus_rows = refine(torus_network(), torus_volume(), torus_area(), (8, 16, 32, 64))
    cone_rows = refine(double_cone_network(), cone_volume(), cone_area(), (8, 16, 32, 64))
    torus_stitched = stitch(torus_network(), nu=16, nv=16)
    cone_stitched = stitch(double_cone_network(), nu=16, nv=16)

    results = {
        "box": {
            "dimensions_m": list(BOX),
            "analytic_volume_m3": a * b * c,
            "analytic_area_m2": 2 * (a * b + b * c + c * a),
            "patches": len(box_network().surfaces),
            "seams": len(box_network().seams),
            "report": box_report,
            "seam_angle_deg": sorted({round(r.angle_max_deg, 12) for r in box_stitched.seams}),
            "seam_gap_max_m": max(r.gap_max for r in box_stitched.seams),
        },
        "torus": {
            "major_radius_m": TORUS[0],
            "minor_radius_m": TORUS[1],
            "analytic_volume_m3": torus_volume(),
            "analytic_area_m2": torus_area(),
            "patches": len(torus_network().surfaces),
            "seams": len(torus_network().seams),
            "refinement": torus_rows,
            "observed_volume_order": observed_order(torus_rows, "volume_error"),
            "observed_area_order": observed_order(torus_rows, "area_error"),
            "seam_angle_max_deg": max(r.angle_max_deg for r in torus_stitched.seams),
            "seam_gap_max_m": max(r.gap_max for r in torus_stitched.seams),
            "inset_normals": sum(r.inset_normals for r in torus_stitched.seams),
        },
        "double_cone": {
            "radius_m": CONE[0],
            "half_height_m": CONE[1],
            "analytic_volume_m3": cone_volume(),
            "analytic_area_m2": cone_area(),
            "patches": len(double_cone_network().surfaces),
            "seams": len(double_cone_network().seams),
            "pole_sides": len(cone_stitched.pole_sides),
            "open_sides": len(cone_stitched.open_sides),
            "refinement": cone_rows,
            "observed_volume_order": observed_order(cone_rows, "volume_error"),
            "observed_area_order": observed_order(cone_rows, "area_error"),
            "inset_normals": sum(r.inset_normals for r in cone_stitched.seams),
            "pole_side_samples": [list(side) for side in cone_stitched.pole_sides],
        },
        "seams": {
            "torus": [record.to_dict() for record in torus_stitched.seams],
            "box": [record.to_dict() for record in box_stitched.seams],
        },
        "negative_controls": negative_controls(),
        "environment": {
            "python": _sys.version.split()[0],
            "generated_by": "examples/patch_network/run_example.py",
            "grid": "common (nu, nv) grid per patch, seams welded by index",
        },
    }

    plot_convergence(torus_rows, cone_rows)
    plot_seams(torus_stitched.seams, box_stitched.seams)
    plot_shell(cone_stitched, "double_cone", "8-patch double cone, 2 apex poles collapsed")
    plot_shell(torus_stitched, "torus", "16-patch torus, 32 G1 seams")

    (OUTPUTS / "results.json").write_text(json.dumps(results, indent=2) + "\n")

    def row(cells) -> str:
        return "| " + " | ".join(str(cell) for cell in cells) + " |"

    lines = [
        "# patch_network: numbers injected from `outputs/results.json`",
        "",
        "## Box 2 x 3 x 4 m (6 bilinear patches, exact tessellation)",
        "",
        row(["quantity", "measured", "analytic", "relative error"]),
        row(["---"] * 4),
        row(["volume [m3]", f"{box_report['signed_volume_m3']:.17g}", f"{a*b*c:.17g}",
             f"{box_report['volume_relative_error']:+.3e}"]),
        row(["area [m2]", f"{box_report['surface_area_m2']:.17g}", f"{2*(a*b+b*c+c*a):.17g}",
             f"{box_report['area_relative_error']:+.3e}"]),
        row(["watertight", box_report["watertight"], "-", "-"]),
        row(["seam gap max [m]", f"{max(r.gap_max for r in box_stitched.seams):.3e}", "0", "-"]),
        "",
        "## Torus R = 2 m, r = 0.5 m (16 patches, 32 seams)",
        "",
        row(["samples", "triangles", "volume [m3]", "volume err", "area [m2]", "area err", "order (vol)"]),
        row(["---"] * 7),
    ]
    for index, record in enumerate(torus_rows):
        order = results["torus"]["observed_volume_order"]
        lines.append(row([
            record["samples_per_patch"], record["triangles"],
            f"{record['volume_m3']:.9f}", f"{record['volume_error']:+.3e}",
            f"{record['area_m2']:.9f}", f"{record['area_error']:+.3e}",
            f"{order[index-1]:.2f}" if index else "-",
        ]))
    lines += [
        "",
        f"analytic: V = 2 pi^2 R r^2 = {torus_volume():.9f} m3, A = 4 pi^2 R r = {torus_area():.9f} m2",
        f"seam gap max = {results['torus']['seam_gap_max_m']:.3e} m, "
        f"seam dihedral max = {results['torus']['seam_angle_max_deg']:.3e} deg",
        "",
        "## Double cone R = 1 m, H = 2 x 1 m (8 patches, 2 apex poles)",
        "",
        row(["samples", "triangles", "volume [m3]", "volume err", "order (vol)"]),
        row(["---"] * 5),
    ]
    for index, record in enumerate(cone_rows):
        order = results["double_cone"]["observed_volume_order"]
        lines.append(row([
            record["samples_per_patch"], record["triangles"],
            f"{record['volume_m3']:.9f}", f"{record['volume_error']:+.3e}",
            f"{order[index-1]:.2f}" if index else "-",
        ]))
    lines += [
        "",
        f"analytic: V = 2 pi R^2 H / 3 = {cone_volume():.9f} m3, "
        f"A = 2 pi R sqrt(R^2 + H^2) = {cone_area():.9f} m2",
        f"pole sides collapsed = {len(cone_stitched.pole_sides)}, open sides = {len(cone_stitched.open_sides)}",
        "",
        "## Negative controls",
        "",
        "```json",
        json.dumps(results["negative_controls"], indent=2),
        "```",
        "",
    ]
    (OUTPUTS / "report.md").write_text("\n".join(lines))

    print(json.dumps({
        "box_volume_error": box_report["volume_relative_error"],
        "box_area_error": box_report["area_relative_error"],
        "box_watertight": box_report["watertight"],
        "torus_volume_errors": [r["volume_error"] for r in torus_rows],
        "torus_area_errors": [r["area_error"] for r in torus_rows],
        "torus_orders_volume": results["torus"]["observed_volume_order"],
        "torus_orders_area": results["torus"]["observed_area_order"],
        "torus_seam_gap_max_m": results["torus"]["seam_gap_max_m"],
        "torus_seam_angle_max_deg": results["torus"]["seam_angle_max_deg"],
        "cone_volume_errors": [r["volume_error"] for r in cone_rows],
        "cone_area_errors": [r["area_error"] for r in cone_rows],
        "cone_orders_volume": results["double_cone"]["observed_volume_order"],
        "cone_watertight": [r["watertight"] for r in cone_rows],
        "negative_controls": results["negative_controls"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())