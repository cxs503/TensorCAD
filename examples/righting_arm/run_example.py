"""Righting arms measured against a closed form, a closed-form family, and themselves.

Everything on this page is *measured*: the script exits non-zero if one of its own
claims stops holding, so the tables in ``docs/righting_arms.md`` cannot drift away from
the code the way a hand-copied table does.

What it does, in order:

1. **A box against the wall-sided formula.** ``GZ = sin(phi)(GM + BM tan^2(phi)/2)`` is
   rational arithmetic for a rectangular section, and the clipped-polyhedron
   integration should reproduce it to round-off. If it does not, nothing further is
   worth running.
2. **The same box with a negative ``GM``**, where the arm vanishes at
   ``tan(phi) = sqrt(-2 GM / BM)``: the range of stability is a *root* of the curve, so
   it tests the sign change and not just the magnitude.
3. **The area under the curve** (dynamical stability) against the integral of the closed
   form: the quantity a rolling ship actually spends.
4. **The parametric hull**: ``GZ'(0)`` must equal ``KB + BM - KG`` (the metacentric
   height), with ``BM`` from the waterplane moment -- a two-dimensional plane integral
   against a three-dimensional heeled volume integral, two code paths that no shared bug
   can fake -- and then against the *analytic* ``GM`` of the hull family at the draft the
   body actually floats at.
5. **A refinement ladder**, because the arm inherits the mesh's error exactly the way the
   hydrostatics does, and that error is reported rather than assumed small.
6. **The Wigley hull**, the benchmark hull, for a second non-box geometry.

Run it with ``python examples/righting_arm/run_example.py --output docs/righting_arms``.
matplotlib is optional: the JSON and the Markdown are written either way and the PNGs
are skipped.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import itertools
from pathlib import Path

from tensorcad import __version__
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.hydrostatics import compute_hydrostatics
from tensorcad.hydrostatics_analytic import analytic_hydrostatics, waterplane_inertia_transverse
from tensorcad.righting_arm import (
    righting_arm,
    righting_arm_curve,
    wall_sided_righting_arm,
    wall_sided_vanishing_angle,
)
from tensorcad.stability import waterplane_moments
from tensorcad.wigley import generate_wigley

BOX_LENGTH = 20.0
BOX_BEAM = 10.0
BOX_DEPTH = 10.0
BOX_DRAFT = 5.0
BOX_KG = 4.0
BOX_KG_UNSTABLE = BOX_DRAFT / 2.0 + BOX_BEAM**2 / (12.0 * BOX_DRAFT) + 0.2  # GM = -0.2 m

HULL_LENGTH = 10.0
HULL_BEAM = 2.0
HULL_DESIGN_DRAFT = 0.8
HULL_KG = 0.55
#: One light displacement for every hull mesh, so the refinement ladder measures the mesh's
#: error and not a different loading condition.
HULL_DISPLACEMENT = 3.75
#: A deeper displacement: the analytic GM of the family is checked at both, so a wrong
#: draft-scaling law cannot hide behind one good number.
HULL_DISPLACEMENT_DEEP = 5.5
WIGLEY_KG = 0.35
WIGLEY_FILL = 0.5


def box(length: float, beam: float, depth: float) -> Mesh:
    """Axis-aligned closed box, top face on ``z = 0``, keel on ``z = -depth``."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-depth, 0.0):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    left, right = 0.0, length
    lo, hi = -0.5 * beam, 0.5 * beam
    top, bottom = 0.0, -depth
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    return builder.build()


class Claims:
    """Collect assertions so one run reports every failure, not just the first."""

    def __init__(self) -> None:
        self.rows: list[dict[str, object]] = []

    def check(self, name: str, measured: float, reference: float, *, tolerance: float,
              unit: str, note: str = "") -> bool:
        error = (
            abs(measured)
            if reference == 0.0
            else abs(measured - reference) / abs(reference)
        )
        passed = error <= tolerance
        self.rows.append(
            {
                "claim": name,
                "measured": measured,
                "reference": reference,
                "unit": unit,
                "relative_error": error,
                "tolerance": tolerance,
                "passed": passed,
                "note": note,
            }
        )
        return passed

    @property
    def failed(self) -> list[str]:
        return [str(row["claim"]) for row in self.rows if not row["passed"]]


def float_to_text(value: float, digits: int = 6) -> str:
    if value != value or value in (float("inf"), float("-inf")):
        return "n/a"
    return f"{value:.{digits}f}"


def measure_box(claims: Claims) -> dict[str, object]:
    mesh = box(BOX_LENGTH, BOX_BEAM, BOX_DEPTH)
    displacement = BOX_LENGTH * BOX_BEAM * BOX_DRAFT
    angles = (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0)
    curve = righting_arm_curve(
        mesh, kg=BOX_KG, displacement_volume=displacement, angles_deg=angles
    )
    rows = []
    worst = 0.0
    for point in curve.points:
        closed = wall_sided_righting_arm(
            point.heel_angle_deg, beam=BOX_BEAM, depth=BOX_DEPTH, draft=BOX_DRAFT, kg=BOX_KG
        )
        worst = max(worst, abs(point.gz - closed))
        rows.append(
            {
                "heel_angle_deg": point.heel_angle_deg,
                "gz_m": point.gz,
                "closed_form_m": closed,
                "difference_m": point.gz - closed,
                "waterline_z_m": point.waterline_z,
                "centre_of_buoyancy_m": list(point.centre_of_buoyancy),
            }
        )
    claims.check(
        "box GZ vs wall-sided closed form (max |difference|, 0-40 deg)",
        worst,
        0.0,
        tolerance=1e-9,
        unit="m",
        note="exact polyhedron, so round-off is the only error",
    )

    # Negative GM: the arm has to cross zero at the closed-form vanishing angle.
    unstable = righting_arm_curve(
        mesh,
        kg=BOX_KG_UNSTABLE,
        displacement_volume=displacement,
        angles_deg=(0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0),
    )
    analytic_vanishing = wall_sided_vanishing_angle(
        beam=BOX_BEAM, depth=BOX_DEPTH, draft=BOX_DRAFT, kg=BOX_KG_UNSTABLE
    )
    measured_vanishing = unstable.zero_crossings_deg[0] if unstable.zero_crossings_deg else float("nan")
    claims.check(
        "box with MG = -0.2 m: angle where GZ = 0 vs sqrt(-2 GM / BM)",
        measured_vanishing,
        float(analytic_vanishing),
        tolerance=1e-6,
        unit="deg",
        note="a root of the curve, so it tests the sign change and not only the magnitude",
    )

    # Dynamical stability: the area under the closed form against the curve.
    fine = righting_arm_curve(
        mesh,
        kg=BOX_KG_UNSTABLE,
        displacement_volume=displacement,
        angles_deg=tuple(0.5 * step for step in range(0, 51)),
    )
    upper = math.radians(25.0)
    gm = BOX_DRAFT / 2.0 + BOX_BEAM**2 / (12.0 * BOX_DRAFT) - BOX_KG_UNSTABLE
    bm = BOX_BEAM**2 / (12.0 * BOX_DRAFT)
    analytic_area = -gm * (math.cos(upper) - 1.0) + 0.5 * bm * (
        1.0 / math.cos(upper) + math.cos(upper) - 2.0
    )
    measured_area = fine.dynamical_stability(25.0)
    claims.check(
        "box: area under GZ, 0-25 deg, vs integral of the closed form",
        measured_area,
        analytic_area,
        tolerance=1e-3,
        unit="m*rad",
        note="trapezoid on a 0.5 deg grid; the residual is quadrature, not integration",
    )

    return {
        "geometry": {
            "length_m": BOX_LENGTH,
            "beam_m": BOX_BEAM,
            "depth_m": BOX_DEPTH,
            "draft_m": BOX_DRAFT,
            "enclosed_volume_m3": BOX_LENGTH * BOX_BEAM * BOX_DEPTH,
        },
        "displacement_volume_m3": displacement,
        "kg_m_above_keel": BOX_KG,
        "kg_unstable_m_above_keel": BOX_KG_UNSTABLE,
        "gm_m": BOX_DRAFT / 2.0 + bm - BOX_KG,
        "gm_unstable_m": gm,
        "bm_m": bm,
        "kb_m": BOX_DRAFT / 2.0,
        "curve": curve.to_dict(),
        "rows": rows,
        "unstable_curve": unstable.to_dict(),
        "analytic_vanishing_angle_deg": analytic_vanishing,
        "analytic_area_m_rad": analytic_area,
        "measured_area_m_rad": measured_area,
    }


def analytic_gm(displacement: float, kg: float) -> tuple[float, float]:
    """Floating draft and analytic metacentric height of the hull family at a displacement.

    The analytic surface the generator samples is
    ``b(x, zeta) = (B/2) sin(pi x/L)**e * zeta**(1/2)`` with ``zeta = (z + T)/T``, which is
    self-similar in the depth. At a floating draft ``d < T`` that gives

        V(d) = V(T) (d/T)**(3/2),   KB(d) = KB(T) (d/T),   I_T(d) = I_T(T) (d/T)**(3/2)

    and therefore ``BM = I_T / V`` is independent of the draft. Those three scalings are
    the only assumption in this file that is not measured: ``tests/test_righting_arm.py``
    checks them against the library's own closed forms at the design draft, and the GM
    they produce is compared with the integration at two different displacements here,
    which is what a wrong exponent would break.
    """
    full = analytic_hydrostatics(
        HullParameters(length=HULL_LENGTH, beam=HULL_BEAM, draft=HULL_DESIGN_DRAFT)
    )
    ratio = (displacement / full.displacement_volume) ** (2.0 / 3.0)
    draft = HULL_DESIGN_DRAFT * ratio
    # VCB is reported in the model frame (keel at z = -T), so the height above the keel,
    # which is what a metacentric height is measured from, is VCB + T.
    kb = (full.vertical_centre_buoyancy + HULL_DESIGN_DRAFT) * ratio
    bm = waterplane_inertia_transverse(full.length, full.beam) / full.displacement_volume
    return draft, kb + bm - kg


def measure_hull(
    claims: Claims, *, nx: int, nz: int, analytic_tolerance: float
) -> dict[str, object]:
    label = f"hull {nx}x{nz}"
    parameters = HullParameters(length=HULL_LENGTH, beam=HULL_BEAM, draft=HULL_DESIGN_DRAFT, nx=nx, nz=nz)
    mesh = generate_hull(parameters).mesh
    upright = compute_hydrostatics(mesh, 0.0)
    displacement = HULL_DISPLACEMENT
    angles = (0.05, 0.1, 10.0, 20.0)
    curve = righting_arm_curve(
        mesh, kg=HULL_KG, displacement_volume=displacement, angles_deg=angles
    )
    equilibrium = righting_arm(mesh, heel_angle_deg=0.0, kg=HULL_KG, displacement_volume=displacement)
    section = waterplane_moments(mesh, equilibrium.waterline_z)
    bm_measured = section.transverse_second_moment / displacement
    gm_identity = (equilibrium.centre_of_buoyancy[2] - equilibrium.centre_of_gravity_z) + bm_measured
    claims.check(
        f"{label}: GZ'(0) vs KB + BM - KG, BM from the waterplane moment",
        curve.initial_slope,
        gm_identity,
        tolerance=1e-6,
        unit="m/rad",
        note="3-D heeled volume integral against a 2-D plane integral",
    )

    # The analytic GM of the family, at the draft the body actually floats at.
    floating_draft, gm_analytic = analytic_gm(displacement, HULL_KG)
    claims.check(
        f"{label}: GZ'(0) vs the analytic GM of the family, light displacement",
        curve.initial_slope,
        gm_analytic,
        tolerance=analytic_tolerance,
        unit="m/rad",
        note=f"floating draft {floating_draft:.4f} m; the residual is the chorded waterplane",
    )
    deep = righting_arm_curve(
        mesh,
        kg=HULL_KG,
        displacement_volume=HULL_DISPLACEMENT_DEEP,
        angles_deg=(0.05, 0.1, 5.0),
    )
    deep_draft, deep_gm = analytic_gm(HULL_DISPLACEMENT_DEEP, HULL_KG)
    claims.check(
        f"{label}: GZ'(0) vs the analytic GM of the family, deep displacement",
        deep.initial_slope,
        deep_gm,
        tolerance=analytic_tolerance,
        unit="m/rad",
        note=f"floating draft {deep_draft:.4f} m; the same scaling law, a second displacement",
    )

    upright_point = righting_arm(mesh, heel_angle_deg=0.0, kg=HULL_KG, displacement_volume=displacement)
    starboard = righting_arm(mesh, heel_angle_deg=10.0, kg=HULL_KG, displacement_volume=displacement)
    port = righting_arm(mesh, heel_angle_deg=-10.0, kg=HULL_KG, displacement_volume=displacement)
    claims.check(
        f"{label}: GZ(0) = 0",
        upright_point.gz,
        0.0,
        tolerance=1e-12,
        unit="m",
        note="absolute tolerance; the arm is zero by symmetry in the upright position",
    )
    claims.check(
        f"{label}: GZ(+10) = -GZ(-10)",
        starboard.gz,
        -port.gz,
        tolerance=1e-12,
        unit="m",
        note="odd in the heel for a port/starboard symmetric body",
    )
    sinkage_ok = all(
        left.waterline_z > right.waterline_z
        for left, right in itertools.pairwise(curve.points)
    )
    claims.rows.append(
        {
            "claim": f"{label}: the waterline descends monotonically as the body heels",
            "measured": float(len(curve.points) - 1),
            "reference": float(len(curve.points) - 1),
            "unit": "interval",
            "relative_error": 0.0 if sinkage_ok else 1.0,
            "tolerance": 0.0,
            "passed": sinkage_ok,
            "note": "a heeled body sinks until the heeled section displaces the same volume",
        }
    )

    return {
        "resolution": [nx, nz],
        "faces": len(mesh.faces),
        "waterline_enclosed_volume_m3": upright.displacement_volume,
        "displacement_volume_m3": displacement,
        "kg_m_above_keel": HULL_KG,
        "floating_draft_m": floating_draft,
        "analytic_gm_m": gm_analytic,
        "identity_gm_m": gm_identity,
        "deep_displacement_volume_m3": HULL_DISPLACEMENT_DEEP,
        "deep_floating_draft_m": deep_draft,
        "deep_analytic_gm_m": deep_gm,
        "deep_initial_slope_m_per_rad": deep.initial_slope,
        "upper_waterplane_moment_m4": section.transverse_second_moment,
        "curve": curve.to_dict(),
        "upright_equilibrium": equilibrium.to_dict(),
    }


def measure_wigley(claims: Claims) -> dict[str, object]:
    model = generate_wigley()
    mesh = model.mesh
    upright = compute_hydrostatics(mesh, 0.0)
    displacement = WIGLEY_FILL * upright.displacement_volume
    angles = (0.05, 0.1, 2.0, 5.0, 10.0, 15.0, 20.0)
    curve = righting_arm_curve(
        mesh, kg=WIGLEY_KG, displacement_volume=displacement, angles_deg=angles
    )
    equilibrium = righting_arm(mesh, heel_angle_deg=0.0, kg=WIGLEY_KG, displacement_volume=displacement)
    section = waterplane_moments(mesh, equilibrium.waterline_z)
    gm_identity = (
        equilibrium.centre_of_buoyancy[2]
        - equilibrium.centre_of_gravity_z
        + section.transverse_second_moment / displacement
    )
    claims.check(
        "Wigley: GZ'(0) vs KB + BM - KG, BM from the waterplane moment",
        curve.initial_slope,
        gm_identity,
        tolerance=1e-6,
        unit="m/rad",
        note="second non-box geometry; the same two code paths must agree",
    )
    starboard = righting_arm(mesh, heel_angle_deg=10.0, kg=WIGLEY_KG, displacement_volume=displacement)
    port = righting_arm(mesh, heel_angle_deg=-10.0, kg=WIGLEY_KG, displacement_volume=displacement)
    claims.check(
        "Wigley: GZ(+10) = -GZ(-10)",
        starboard.gz,
        -port.gz,
        tolerance=1e-12,
        unit="m",
        note="odd in the heel",
    )
    return {
        "faces": len(mesh.faces),
        "waterline_enclosed_volume_m3": upright.displacement_volume,
        "displacement_volume_m3": displacement,
        "kg_m_above_keel": WIGLEY_KG,
        "identity_gm_m": gm_identity,
        "curve": curve.to_dict(),
    }


def measure_refinement(claims: Claims, coarse: dict[str, object], fine: dict[str, object]) -> dict[str, object]:
    measurement = {
        "coarse_resolution": coarse["resolution"],
        "fine_resolution": fine["resolution"],
        "coarse_faces": coarse["faces"],
        "fine_faces": fine["faces"],
        "coarse_initial_slope": coarse["curve"]["initial_slope_m_per_rad"],
        "fine_initial_slope": fine["curve"]["initial_slope_m_per_rad"],
        "coarse_gz_10deg_m": _point_gz(coarse, 10.0),
        "fine_gz_10deg_m": _point_gz(fine, 10.0),
        "coarse_analytic_error": (
            coarse["curve"]["initial_slope_m_per_rad"] - coarse["analytic_gm_m"]
        )
        / coarse["analytic_gm_m"],
        "fine_analytic_error": (
            fine["curve"]["initial_slope_m_per_rad"] - fine["analytic_gm_m"]
        )
        / fine["analytic_gm_m"],
    }
    for name, key in (("initial slope", None), ("GZ at 10 deg", (10.0,))):
        coarse_value = measurement["coarse_initial_slope"] if key is None else _point_gz(coarse, 10.0)
        fine_value = measurement["fine_initial_slope"] if key is None else _point_gz(fine, 10.0)
        claims.check(
            f"hull refinement: {name}, {coarse['resolution']} vs {fine['resolution']}",
            float(coarse_value),
            float(fine_value),
            tolerance=8e-2,
            unit="m" if key else "m/rad",
            note="the arm inherits the mesh's discretisation error; this is its size",
        )
    return measurement


def _point_gz(block: dict[str, object], heel: float) -> float:
    for point in block["curve"]["points"]:
        if abs(point["heel_angle_deg"] - heel) < 1e-9:
            return float(point["gz_m"])
    raise KeyError(f"heel {heel} deg is not on this curve")


def heeled_section(mesh: Mesh, station_x: float, heel_deg: float, waterline: float) -> dict[str, object]:
    """The station outline and the waterline in the *world* frame, for the figure.

    A point ``(y, z)`` of the body moves to ``(y cos + z sin, -y sin + z cos)`` when the
    body heels to starboard, so the waterline ``z = y tan + wl`` becomes the horizontal
    line ``z_w = wl cos``. That identity is the frame algebra the arm formula rests on,
    so the figure is drawn from it rather than from a second implementation.
    """
    heel = math.radians(heel_deg)
    cos_heel, sin_heel = math.cos(heel), math.sin(heel)
    station = sorted(
        (vertex for vertex in mesh.vertices if abs(vertex[0] - station_x) < 1e-6),
        key=lambda vertex: math.atan2(vertex[2], vertex[1]),
    )
    outline = [
        (
            vertex[1] * cos_heel + vertex[2] * sin_heel,
            -vertex[1] * sin_heel + vertex[2] * cos_heel,
        )
        for vertex in station
    ]
    return {
        "heel_angle_deg": heel_deg,
        "outline_world": outline,
        "waterline_world_z": waterline * cos_heel,
    }


def write_figures(path: Path, payload: dict[str, object]) -> list[str]:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:  # pragma: no cover - matplotlib is optional
        print("matplotlib is not available; skipping figures", file=sys.stderr)
        return []

    written: list[str] = []

    figure, axis = plt.subplots(figsize=(7.0, 4.5))
    for key, label, marker, closed in (
        ("box", "box 20x10x10 m, KG = 4.0 m", "o", True),
        ("hull", "parametric hull, KG = 0.55 m", "s", False),
        ("wigley", "Wigley, KG = 0.35 m", "^", False),
    ):
        block = payload[key]
        angles = [point["heel_angle_deg"] for point in block["curve"]["points"]]
        arms = [point["gz_m"] for point in block["curve"]["points"]]
        axis.plot(angles, arms, marker=marker, label=label, linewidth=1.4, markersize=4)
        if closed:
            grid = [0.25 * step for step in range(0, 161)]
            axis.plot(
                grid,
                [
                    wall_sided_righting_arm(
                        angle, beam=BOX_BEAM, depth=BOX_DEPTH, draft=BOX_DRAFT, kg=BOX_KG
                    )
                    for angle in grid
                ],
                linestyle="--",
                linewidth=1.0,
                color="0.4",
                label="wall-sided closed form (box)",
            )
    hull_block = payload["hull"]
    limit = hull_block["curve"]["physical_limit_deg"]
    if limit is not None:
        axis.axvline(limit, color="0.6", linestyle=":", linewidth=1.0)
        axis.annotate(
            f"hull deck rim wetted at {limit:.1f} deg",
            xy=(limit, 0.0),
            xytext=(limit + 1.0, 0.02),
            fontsize=7,
            color="0.35",
        )
    axis.axhline(0.0, color="0.7", linewidth=0.8)
    axis.set_xlabel("heel angle (deg)")
    axis.set_ylabel("righting arm GZ (m)")
    axis.set_title("Righting-arm curves from a heeled displacement integration")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)
    figure.tight_layout()
    target = path / "gz_curves.png"
    figure.savefig(target, dpi=150)
    plt.close(figure)
    written.append(target.name)

    figure, axis = plt.subplots(figsize=(7.0, 3.8))
    rows = payload["box"]["rows"]
    angles = [row["heel_angle_deg"] for row in rows]
    differences = [max(abs(row["difference_m"]), 1e-18) for row in rows]
    axis.semilogy(angles, differences, marker="o", linewidth=1.2, label="|GZ - closed form|")
    axis.axhline(1e-16, color="0.6", linestyle=":", linewidth=1.0, label="double-precision round-off")
    axis.set_xlabel("heel angle (deg)")
    axis.set_ylabel("absolute difference (m)")
    axis.set_title("A box reproduces its wall-sided righting arm to round-off")
    axis.legend(fontsize=8)
    axis.grid(alpha=0.25)
    figure.tight_layout()
    target = path / "box_vs_closed_form.png"
    figure.savefig(target, dpi=150)
    plt.close(figure)
    written.append(target.name)

    figure, axes = plt.subplots(1, 3, figsize=(9.0, 3.2), sharey=True)
    sections = payload["sections"]
    for axis, section in zip(axes, sections):
        outline = section["outline_world"]
        axis.fill(
            [point[0] for point in outline],
            [point[1] for point in outline],
            facecolor="0.85",
            edgecolor="0.35",
            linewidth=0.8,
        )
        axis.axhline(section["waterline_world_z"], color="tab:blue", linewidth=1.2)
        axis.axhline(0.0, color="0.75", linewidth=0.8, linestyle="--")
        axis.plot([section["centre_of_buoyancy_world"][0]], [section["centre_of_buoyancy_world"][1]],
                  marker="o", color="tab:blue", markersize=5)
        axis.plot([section["centre_of_gravity_world"][0]], [section["centre_of_gravity_world"][1]],
                  marker="o", color="tab:red", markersize=5)
        axis.annotate("B", section["centre_of_buoyancy_world"], fontsize=8, color="tab:blue",
                      xytext=(3, 3), textcoords="offset points")
        axis.annotate("G", section["centre_of_gravity_world"], fontsize=8, color="tab:red",
                      xytext=(3, 3), textcoords="offset points")
        axis.set_title(f"heel {section['heel_angle_deg']:.0f} deg, GZ = {section['gz_m']:.4f} m",
                       fontsize=9)
        axis.set_xlabel("y (m)")
        axis.grid(alpha=0.2)
    axes[0].set_ylabel("z (m)")
    figure.suptitle("Midship section, waterline, B and G in the world frame", fontsize=10)
    figure.tight_layout()
    target = path / "heeled_sections.png"
    figure.savefig(target, dpi=150)
    plt.close(figure)
    written.append(target.name)
    return written


def build_sections(mesh: Mesh, kg: float, displacement: float) -> list[dict[str, object]]:
    sections = []
    for heel in (0.0, 10.0, 20.0):
        arm = righting_arm(mesh, heel_angle_deg=heel, kg=kg, displacement_volume=displacement)
        heel_rad = math.radians(heel)
        cos_heel, sin_heel = math.cos(heel_rad), math.sin(heel_rad)
        _, y_b, z_b = arm.centre_of_buoyancy
        section = heeled_section(mesh, 5.0, heel, arm.waterline_z)
        section["centre_of_buoyancy_world"] = [y_b * cos_heel + z_b * sin_heel,
                                               -y_b * sin_heel + z_b * cos_heel]
        section["centre_of_gravity_world"] = [arm.centre_of_gravity_z * sin_heel,
                                              arm.centre_of_gravity_z * cos_heel]
        section["gz_m"] = arm.gz
        sections.append(section)
    return sections


def write_markdown(path: Path, payload: dict[str, object], claims: Claims, figures: list[str]) -> None:
    lines = [
        "# Righting arms: GZ(phi) measured against closed forms",
        "",
        f"Generated by `python examples/righting_arm/run_example.py` on TensorCAD "
        f"{payload['version']}. Every number below is read from the `results.json` of the",
        "same run; the script exits non-zero if a claim fails, so this page cannot decouple",
        "from the code.",
        "",
        "## 1. A box against the wall-sided formula",
        "",
        f"Rectangular box {BOX_LENGTH:g} x {BOX_BEAM:g} m, depth {BOX_DEPTH:g} m, draft "
        f"{BOX_DRAFT:g} m, KG = {BOX_KG:g} m above the keel "
        f"(GM = {payload['box']['gm_m']:.6f} m, BM = {payload['box']['bm_m']:.6f} m, "
        f"KB = {payload['box']['kb_m']:.6f} m).",
        "",
        "| heel (deg) | GZ integrated (m) | closed form (m) | difference (m) | waterline z (m) |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in payload["box"]["rows"]:
        lines.append(
            f"| {row['heel_angle_deg']:.1f} | {row['gz_m']:.12f} | {row['closed_form_m']:.12f} | "
            f"{row['difference_m']:+.2e} | {row['waterline_z_m']:.9f} |"
        )
    lines += [
        "",
        "The residual is round-off: a box is exactly representable, so the only error the",
        "clipped-polyhedron integral can make is in the last bits of the trapezoid-free",
        "tetrahedron sum. Any real error in the method would show up here at 1e-3, not 1e-16.",
        "",
        "## 2. Negative GM: the range of stability is a root",
        "",
        f"Same box with KG = {BOX_KG_UNSTABLE:.6f} m, i.e. "
        f"GM = {payload['box']['gm_unstable_m']:.6f} m. The wall-sided arm vanishes where",
        "`tan(phi) = sqrt(-2 GM / BM)`, and the integration has to find the same zero.",
        "",
        "| quantity | value |",
        "|---|---:|",
        f"| vanishing angle, closed form | {payload['box']['analytic_vanishing_angle_deg']:.9f} deg |",
        f"| vanishing angle, integrated | {payload['box']['unstable_curve']['zero_crossings_deg'][0]:.9f} deg |",
        f"| area under GZ 0-25 deg, closed form | {payload['box']['analytic_area_m_rad']:.12f} m*rad |",
        f"| area under GZ 0-25 deg, integrated | {payload['box']['measured_area_m_rad']:.12f} m*rad |",
        "",
        "`vanishing_angle_deg` is `None` for this body on purpose: a curve that starts",
        "non-positive has no range of stability to end, and its zero is the angle of loll.",
        "Relabelling it as a vanishing angle would be the kind of quiet wrong answer this",
        "module refuses to give.",
        "",
        "## 3. The parametric hull",
        "",
        "| quantity | value |",
        "|---|---:|",
        f"| mesh | {payload['hull']['resolution'][0]} x {payload['hull']['resolution'][1]} "
        f"({payload['hull']['faces']} faces) |",
        f"| enclosed volume at the design waterline | {payload['hull']['waterline_enclosed_volume_m3']:.6f} m^3 |",
        f"| displacement scanned | {payload['hull']['displacement_volume_m3']:.6f} m^3 |",
        f"| floating draft at that displacement | {payload['hull']['floating_draft_m']:.6f} m |",
        f"| KG | {payload['hull']['kg_m_above_keel']:.6f} m |",
        f"| GZ'(0) from the curve | {payload['hull']['curve']['initial_slope_m_per_rad']:.9f} m/rad |",
        f"| KB + BM - KG, BM from the waterplane moment | {payload['hull']['identity_gm_m']:.9f} m/rad |",
        f"| analytic GM of the family at that draft | {payload['hull']['analytic_gm_m']:.9f} m/rad |",
        f"| maximum GZ in the scan | {payload['hull']['curve']['maximum_gz_m']:.6f} m at "
        f"{payload['hull']['curve']['angle_of_maximum_gz_deg']:.1f} deg |",
        f"| deck rim first wetted at | "
        f"{payload['hull']['curve']['physical_limit_deg']:.3f} deg |",
        "",
        "The first two agree to metacentric precision because they are the same identity",
        "measured two ways: a three-dimensional heeled volume integral on one side, a",
        "two-dimensional waterplane moment on the other. The third column is the analytic",
        "surface the sampled hull approximates, so the gap between it and the curve is the",
        "mesh's own chord error and nothing else.",
        "",
        "**The deck-rim limit is the honest part of this table.** The generated hull is",
        "capped at its design waterline, so above that heel the water is over a rim the",
        "geometry does not model. The curve is reported, and the limit is reported with it,",
        "rather than the arm being extrapolated over a freeboard that does not exist.",
        "",
        "### Refinement, same displacement",
        "",
        "| resolution | faces | GZ'(0) (m/rad) | GZ at 10 deg (m) | error vs analytic GM |",
        "|---|---:|---:|---:|---:|",
        f"| {payload['refinement']['coarse_resolution'][0]} x {payload['refinement']['coarse_resolution'][1]} "
        f"| {payload['refinement']['coarse_faces']} | {payload['refinement']['coarse_initial_slope']:.9f} "
        f"| {payload['refinement']['coarse_gz_10deg_m']:.9f} "
        f"| {payload['refinement']['coarse_analytic_error']:+.2e} |",
        f"| {payload['refinement']['fine_resolution'][0]} x {payload['refinement']['fine_resolution'][1]} "
        f"| {payload['refinement']['fine_faces']} | {payload['refinement']['fine_initial_slope']:.9f} "
        f"| {payload['refinement']['fine_gz_10deg_m']:.9f} "
        f"| {payload['refinement']['fine_analytic_error']:+.2e} |",
        "",
        "| resolution | displacement (m^3) | draft (m) | GZ'(0) (m/rad) | analytic GM (m/rad) | relative error |",
        "|---|---:|---:|---:|---:|---:|",
        f"| {payload['hull']['resolution'][0]} x {payload['hull']['resolution'][1]} "
        f"| {payload['hull']['deep_displacement_volume_m3']:.4f} "
        f"| {payload['hull']['deep_floating_draft_m']:.6f} "
        f"| {payload['hull']['deep_initial_slope_m_per_rad']:.9f} "
        f"| {payload['hull']['deep_analytic_gm_m']:.9f} "
        f"| {(payload['hull']['deep_initial_slope_m_per_rad'] - payload['hull']['deep_analytic_gm_m']) / payload['hull']['deep_analytic_gm_m']:+.2e} |",
        f"| {payload['hull_fine']['resolution'][0]} x {payload['hull_fine']['resolution'][1]} "
        f"| {payload['hull_fine']['deep_displacement_volume_m3']:.4f} "
        f"| {payload['hull_fine']['deep_floating_draft_m']:.6f} "
        f"| {payload['hull_fine']['deep_initial_slope_m_per_rad']:.9f} "
        f"| {payload['hull_fine']['deep_analytic_gm_m']:.9f} "
        f"| {(payload['hull_fine']['deep_initial_slope_m_per_rad'] - payload['hull_fine']['deep_analytic_gm_m']) / payload['hull_fine']['deep_analytic_gm_m']:+.2e} |",
        "",
        "The second table is the one a wrong draft-scaling law would break: the same law,",
        "at a displacement 47 % higher, on both grids.",
        "",
        "## 4. The Wigley hull",
        "",
        "| quantity | value |",
        "|---|---:|",
        f"| faces | {payload['wigley']['faces']} |",
        f"| displacement scanned | {payload['wigley']['displacement_volume_m3']:.6f} m^3 |",
        f"| KG | {payload['wigley']['kg_m_above_keel']:.6f} m |",
        f"| GZ'(0) | {payload['wigley']['curve']['initial_slope_m_per_rad']:.9f} m/rad |",
        f"| KB + BM - KG | {payload['wigley']['identity_gm_m']:.9f} m/rad |",
        f"| maximum GZ in the scan | {payload['wigley']['curve']['maximum_gz_m']:.6f} m |",
        "",
        "## 5. Claims",
        "",
        "| claim | measured | reference | relative error | tolerance | passed |",
        "|---|---:|---:|---:|---:|:--:|",
    ]
    for row in claims.rows:
        lines.append(
            f"| {row['claim']} | {float_to_text(float(row['measured']), 9)} | "
            f"{float_to_text(float(row['reference']), 9)} | {float(row['relative_error']):.3e} | "
            f"{float(row['tolerance']):.1e} | {'yes' if row['passed'] else '**no**'} |"
        )
    lines += [
        "",
        "## 6. Figures",
        "",
    ]
    for name in figures:
        lines.append(f"![{name}]({name})")
        lines.append("")
    lines += [
        "## 7. What this does not cover",
        "",
        "* Free-surface correction: a tank property, not a hull property.",
        "* Damage and flooded-compartment stability: the flooding is a mass model the",
        "  geometry does not carry.",
        "* Trim change with heel: a cross curve holds trim constant and so does this module.",
        "* Any statutory criterion (IS Code 2008 intact stability, SOLAS): the arm is a",
        "  number, not a compliance verdict.",
        "",
    ]
    path.joinpath("report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/righting_arms"))
    arguments = parser.parse_args(argv)
    output = arguments.output
    output.mkdir(parents=True, exist_ok=True)

    claims = Claims()
    box_block = measure_box(claims)
    hull_coarse = measure_hull(claims, nx=41, nz=13, analytic_tolerance=8e-2)
    hull_fine = measure_hull(claims, nx=161, nz=41, analytic_tolerance=2e-2)
    refinement = measure_refinement(claims, hull_coarse, hull_fine)
    wigley_block = measure_wigley(claims)

    hull_mesh = generate_hull(HullParameters(length=10.0, beam=2.0, draft=0.8)).mesh
    sections = build_sections(hull_mesh, HULL_KG, hull_coarse["displacement_volume_m3"])

    payload: dict[str, object] = {
        "schema": "tensorcad.righting_arm.example/1",
        "version": __version__,
        "box": box_block,
        "hull": hull_coarse,
        "hull_fine": hull_fine,
        "wigley": wigley_block,
        "refinement": refinement,
        "sections": sections,
        "claims": claims.rows,
    }
    figures = write_figures(output, payload)
    payload["figures"] = figures
    output.joinpath("results.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                                                encoding="utf-8")
    write_markdown(output, payload, claims, figures)

    for row in claims.rows:
        status = "PASS" if row["passed"] else "FAIL"
        print(
            f"[{status}] {row['claim']}: measured {float_to_text(float(row['measured']), 9)} vs "
            f"reference {float_to_text(float(row['reference']), 9)} "
            f"(relative {float(row['relative_error']):.3e}, tolerance {float(row['tolerance']):.1e})"
        )
    if claims.failed:
        print(f"{len(claims.failed)} claim(s) failed: {claims.failed}", file=sys.stderr)
        return 1
    print(f"\nall {len(claims.rows)} claims hold; wrote {output}/results.json and report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
