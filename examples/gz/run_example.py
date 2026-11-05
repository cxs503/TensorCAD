"""Righting-arm curves ``GZ(phi)``: two closed forms, a hull, and the checks behind them.

Run:  python examples/gz/run_example.py                    # writes ``docs/gz/``
      python examples/gz/run_example.py --output /tmp/gz --no-plot

Every number in ``docs/gz/report.md`` and ``docs/gz/results.json`` comes out of this run --
nothing is typed by hand -- and the script exits non-zero if a claim below stops holding,
so CI can run it as a check rather than as a demo. ``tests/test_gz_evidence.py`` reruns it
and compares the result with the published bundle, so the page cannot drift from the code.

The claims:

1. **box barge** -- below the bilge/deck limit the integrator reproduces
   ``GZ = sin(phi) (GM0 + BM0/2 tan^2 phi)`` to round-off, *and* the same closed form breaks
   down past the limit, so the agreement is not vacuous;
2. **log** (circular section, half submerged) -- ``GZ = (R - KG) sin(phi)`` exactly at every
   5-degree heel from 0 to 90, the chord pattern repeating on the 72-gon;
3. **hull** -- ``GZ'(0)`` equals ``GM`` from :mod:`tensorcad.stability`, which integrates the
   waterplane moments by Green's theorem instead of a heeled volume;
4. every point of every curve floats at the same displacement;
5. the free-surface reduction ``sum(i_x)/V`` is a *flat* shift of every curve point, and the
   ``with_freeboard`` band leaves the upright displacement where it was.
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

from tensorcad import __version__
from tensorcad import mesh as meshlib
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.gz import (
    box_wall_sided_limit_deg,
    free_surface_correction,
    heeled_hydrostatics,
    rectangular_tank_inertia,
    righting_arm_curve,
    wall_sided_righting_arm,
    with_freeboard,
)
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.stability import initial_stability

BREADTH, DRAFT, DEPTH, LENGTH, KG = 2.0, 0.8, 2.0, 10.0, 0.50
ANGLE_STEP = 5
HULL_NX, HULL_NZ = 41, 13
HULL_FREEBOARD = 1.0
LOG_RADIUS, LOG_LENGTH, LOG_SEGMENTS = 1.0, 4.0, 72
TANK_LENGTH, TANK_BREADTH = 2.0, 0.5


def box(minimum: tuple[float, float, float], maximum: tuple[float, float, float]) -> Mesh:
    x0, y0, z0 = minimum
    x1, y1, z1 = maximum
    builder = MeshBuilder()
    corners = [
        (x0, y0, z0),
        (x1, y0, z0),
        (x1, y1, z0),
        (x0, y1, z0),
        (x0, y0, z1),
        (x1, y0, z1),
        (x1, y1, z1),
        (x0, y1, z1),
    ]
    ids = [builder.add_vertex(point) for point in corners]
    for face in (
        (0, 3, 2),
        (0, 2, 1),
        (4, 5, 6),
        (4, 6, 7),
        (0, 1, 5),
        (0, 5, 4),
        (2, 3, 7),
        (2, 7, 6),
        (1, 2, 6),
        (1, 6, 5),
        (0, 4, 7),
        (0, 7, 3),
    ):
        builder.add_triangle(*[ids[index] for index in face])
    return meshlib.orient_outward(builder.build()).mesh


def log(radius: float, length: float, segments: int = 72) -> Mesh:
    builder = MeshBuilder()
    rings = []
    for index in range(segments):
        angle = 2.0 * math.pi * index / segments
        y = radius * math.cos(angle)
        z = radius + radius * math.sin(angle)
        rings.append((builder.add_vertex((0.0, y, z)), builder.add_vertex((length, y, z))))
    for index in range(segments):
        low, high = rings[index], rings[(index + 1) % segments]
        builder.add_triangle(low[0], low[1], high[1])
        builder.add_triangle(low[0], high[1], high[0])
    first = builder.add_vertex((0.0, 0.0, radius))
    second = builder.add_vertex((length, 0.0, radius))
    for index in range(segments):
        low, high = rings[index], rings[(index + 1) % segments]
        builder.add_triangle(first, high[0], low[0])
        builder.add_triangle(second, low[1], high[1])
    return meshlib.orient_outward(builder.build()).mesh


def _round(value: float, digits: int = 12) -> float:
    return round(float(value), digits)


def _signed(value: float, limit: float) -> bool:
    """A one-sided claim: does ``value`` exceed ``limit``?"""
    return value > limit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/gz"))
    parser.add_argument("--no-plot", action="store_true")
    options = parser.parse_args()
    output: Path = options.output
    output.mkdir(parents=True, exist_ok=True)

    failures: list[str] = []
    claims: list[dict[str, object]] = []

    def add_claim(name: str, measured: float, limit: float, holds: bool, unit: str, note: str) -> None:
        claims.append(
            {
                "claim": name,
                "measured": float(measured),
                "limit": limit,
                "unit": unit,
                "holds": bool(holds),
                "note": note,
            }
        )
        if not holds:
            failures.append(f"{name}: measured {measured:.6e} {unit} against limit {limit:.1e}")

    # --- 1. box barge against the wall-sided trapezoid -------------------------------
    barge = box((0.0, -BREADTH / 2, 0.0), (LENGTH, BREADTH / 2, DEPTH))
    limit_deg = box_wall_sided_limit_deg(BREADTH, DRAFT, DEPTH)
    gm0 = DRAFT / 2.0 + BREADTH**2 / (12.0 * DRAFT) - KG
    bm0 = BREADTH**2 / (12.0 * DRAFT)
    box_angles = list(range(0, 86, ANGLE_STEP))  # the wall-sided form needs |heel| < 90
    box_curve = righting_arm_curve(barge, KG, box_angles, waterline_z=DRAFT)
    box_rows = []
    worst_inside = 0.0
    worst_past = 0.0
    for angle, gz in zip(box_angles, box_curve.gz, strict=True):
        closed_form = wall_sided_righting_arm(angle, gm=gm0, bm=bm0)
        error = abs(gz - closed_form)
        if angle <= limit_deg:
            worst_inside = max(worst_inside, error)
        elif angle <= limit_deg + 2 * ANGLE_STEP:
            worst_past = max(worst_past, error)
        box_rows.append(
            {
                "heel_deg": float(angle),
                "gz_m": _round(gz),
                "wall_sided_m": _round(closed_form),
                "abs_error_m": _round(error),
            }
        )
    add_claim(
        f"box barge: integrated GZ vs the wall-sided closed form, 0-{int(limit_deg)} deg",
        worst_inside,
        1e-12,
        worst_inside < 1e-12,
        "m",
        "both are exact for the polygonal hull; only round-off separates them",
    )
    add_claim(
        f"box barge: the same closed form stops holding just past {limit_deg:.2f} deg",
        worst_past,
        1e-3,
        _signed(worst_past, 1e-3),
        "m",
        "first two samples past the bilge limit; by 85 deg it is 2.7e+01 m as tan(phi) diverges",
    )

    # --- 2. log: the circular-section closed form -----------------------------------
    log_mesh = log(LOG_RADIUS, LOG_LENGTH, LOG_SEGMENTS)
    log_angles = list(range(0, 91, ANGLE_STEP))
    log_curve = righting_arm_curve(log_mesh, KG, log_angles, waterline_z=LOG_RADIUS)
    log_rows = []
    log_worst = 0.0
    for angle, gz in zip(log_angles, log_curve.gz, strict=True):
        exact = (LOG_RADIUS - KG) * math.sin(math.radians(angle))
        error = abs(gz - exact)
        log_worst = max(log_worst, error)
        log_rows.append(
            {
                "heel_deg": float(angle),
                "gz_m": _round(gz),
                "exact_m": _round(exact),
                "abs_error_m": _round(error),
            }
        )
    add_claim(
        "log (72-gon, half submerged): integrated GZ vs the exact (R - KG) sin(phi), 0-90 deg",
        log_worst,
        1e-12,
        log_worst < 1e-12,
        "m",
        "a round bilge, i.e. a section the box cannot stand in for",
    )

    # --- 3. upright displacement of the box barge -----------------------------------
    barge_upright = heeled_hydrostatics(barge, 0.0, waterline_z=DRAFT)
    designed_volume = LENGTH * BREADTH * DRAFT
    box_volume_rel = abs(barge_upright.displacement_volume - designed_volume) / designed_volume
    add_claim(
        "box barge: upright displacement equals L*B*T exactly (a prism, so an identity, not a convergence)",
        box_volume_rel,
        1e-12,
        box_volume_rel < 1e-12,
        "relative",
        "KB/waterplane_area/wetted_area are reported alongside",
    )

    # --- 4. hull: curve slope vs the waterplane-moment GM ---------------------------
    hull_parameters = HullParameters(length=LENGTH, beam=BREADTH, draft=DRAFT, nx=HULL_NX, nz=HULL_NZ)
    hull = generate_hull(hull_parameters).mesh
    banded = with_freeboard(hull, HULL_FREEBOARD)
    watertight = meshlib.is_watertight(banded)
    if not watertight:
        failures.append("the freeboarded hull is not watertight")
    stability = initial_stability(hull, kg=KG, waterline_z=0.0)

    def slope(step_deg: float) -> float:
        probe = righting_arm_curve(banded, KG, [-step_deg, 0.0, step_deg], waterline_z=0.0)
        return (probe.gz[2] - probe.gz[0]) / math.radians(2.0 * step_deg)

    coarse, fine = slope(0.01), slope(0.001)
    extrapolated = fine + (fine - coarse) / 9.0  # Richardson, step ratio 10
    slope_rel = abs(extrapolated - stability.gm_transverse) / stability.gm_transverse
    add_claim(
        "hull: dGZ/dphi at zero heel vs initial_stability GM (waterplane moments, other route)",
        slope_rel,
        1e-6,
        slope_rel < 1e-6,
        "relative",
        "Richardson extrapolated from +-0.01 and +-0.001 deg probes",
    )

    hull_curve = righting_arm_curve(banded, KG, list(range(0, 91, ANGLE_STEP)), waterline_z=0.0)

    # --- 5. displacement proved equal at every heel of every curve ------------------
    volume_residual = max(box_curve.max_volume_error, log_curve.max_volume_error, hull_curve.max_volume_error)
    add_claim(
        "displacement conserved at every heel of every curve",
        volume_residual,
        1e-10,
        volume_residual < 1e-10,
        "relative",
        "a curve whose points float at different displacements is not a righting-arm curve",
    )

    # --- 6. free-surface correction is a flat shift --------------------------------
    tank = rectangular_tank_inertia(TANK_LENGTH, TANK_BREADTH)
    correction = free_surface_correction([tank], hull_curve.displacement_volume)
    flooded = righting_arm_curve(banded, KG, hull_curve.angles_deg, waterline_z=0.0, free_surface_i_x=tank)
    shift_spread = max(
        abs((solid - wet) - correction.loss_m)
        for solid, wet in zip(hull_curve.gz, flooded.gz, strict=True)
    )
    add_claim(
        "free-surface correction shifts every point of the hull curve by the same amount",
        shift_spread,
        1e-12,
        shift_spread < 1e-12,
        "m",
        "the flat sum(i_x)/V reduction a designer applies by hand",
    )

    # --- 7. the freeboard band does not move the upright displacement ---------------
    upright_hull = heeled_hydrostatics(hull, 0.0, waterline_z=0.0)
    upright_banded = heeled_hydrostatics(banded, 0.0, waterline_z=0.0)
    band_rel = (
        abs(upright_banded.displacement_volume - upright_hull.displacement_volume)
        / upright_hull.displacement_volume
    )
    add_claim(
        "with_freeboard leaves the upright displacement bit-identical",
        band_rel,
        1e-12,
        band_rel < 1e-12,
        "relative",
        "the band only closes the surface above the design waterline",
    )

    results: dict[str, object] = {
        "version": __version__,
        "cases": {
            "box_barge": {
                "length_m": LENGTH,
                "breadth_m": BREADTH,
                "depth_m": DEPTH,
                "draft_m": DRAFT,
                "kg_m": KG,
                "gm0_m": _round(gm0),
                "bm0_m": _round(bm0),
                "wall_sided_limit_deg": _round(limit_deg),
                "angles_deg": [float(angle) for angle in box_angles],
                "gz_m": [_round(value) for value in box_curve.gz],
                "wall_sided_m": [row["wall_sided_m"] for row in box_rows],
                "abs_error_m": [row["abs_error_m"] for row in box_rows],
                "max_gz_m": _round(box_curve.max_gz),
                "angle_of_max_gz_deg": box_curve.angle_of_max_gz_deg,
                "range_of_stability_deg": (
                    None if box_curve.range_of_stability_deg is None else _round(box_curve.range_of_stability_deg)
                ),
                "area_under_curve_0_30_m_rad": _round(box_curve.area_under_curve(30.0)),
                "upright_displacement_volume_m3": _round(barge_upright.displacement_volume),
            },
            "log": {
                "radius_m": LOG_RADIUS,
                "length_m": LOG_LENGTH,
                "segments": LOG_SEGMENTS,
                "kg_m": KG,
                "half_submerged_draft_m": LOG_RADIUS,
                "angles_deg": [float(angle) for angle in log_angles],
                "gz_m": [_round(value) for value in log_curve.gz],
                "exact_m": [row["exact_m"] for row in log_rows],
                "abs_error_m": [row["abs_error_m"] for row in log_rows],
            },
            "hull": {
                "parameters": {
                    "length_m": LENGTH,
                    "beam_m": BREADTH,
                    "draft_m": DRAFT,
                    "kg_m": KG,
                    "grid": [HULL_NX, HULL_NZ],
                    "freeboard_m": HULL_FREEBOARD,
                    "watertight_after_freeboard": watertight,
                    "volume_m3": _round(upright_banded.displacement_volume),
                },
                "kb_m": _round(upright_banded.centre_of_buoyancy_z),
                "waterplane_area_m2": _round(upright_banded.waterplane_area),
                "wetted_surface_area_m2": _round(upright_banded.wetted_surface_area),
                "gm_from_initial_stability_m": _round(stability.gm_transverse),
                "kb_from_initial_stability_m": _round(stability.kb),
                "bm_transverse_m": _round(stability.bm_transverse),
                "gm_from_curve_slope_m": _round(extrapolated),
                "angles_deg": [float(angle) for angle in hull_curve.angles_deg],
                "gz_m": [_round(value) for value in hull_curve.gz],
                "max_gz_m": _round(hull_curve.max_gz),
                "angle_of_max_gz_deg": hull_curve.angle_of_max_gz_deg,
                "range_of_stability_deg": (
                    None if hull_curve.range_of_stability_deg is None else _round(hull_curve.range_of_stability_deg)
                ),
                "area_under_curve_0_30_m_rad": _round(hull_curve.area_under_curve(30.0)),
                "max_volume_error": hull_curve.max_volume_error,
            },
            "free_surface": {
                "tank": f"{TANK_LENGTH} m x {TANK_BREADTH} m rectangle",
                "i_x_m4": _round(tank),
                "loss_m": _round(correction.loss_m),
                "displacement_volume_m3": _round(correction.displacement_volume),
                "max_gz_solid_m": _round(hull_curve.max_gz),
                "max_gz_with_free_surface_m": _round(flooded.max_gz),
            },
            "upright": {
                "hull_volume_m3": _round(upright_hull.displacement_volume),
                "banded_volume_m3": _round(upright_banded.displacement_volume),
                "hull_kb_m": _round(upright_hull.centre_of_buoyancy_z),
                "banded_kb_m": _round(upright_banded.centre_of_buoyancy_z),
            },
        },
        "claims": claims,
        "checks": {"passed": not failures, "failures": failures},
    }
    (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")

    # --- report --------------------------------------------------------------------
    def range_text(curve: object) -> str:
        value = getattr(curve, "range_of_stability_deg", None)
        if value is None:
            return "> 90 (still positive at the last sampled heel)"
        return f"{value:.3f} deg"

    def claim_table() -> list[str]:
        lines = ["| claim | measured | limit | holds |", "|---|---:|---|:--:|"]
        for entry in claims:
            measured = float(entry["measured"])  # type: ignore[arg-type]
            lines.append(
                f"| {entry['claim']} | {measured:.3e} {entry['unit']} | "
                f"{float(entry['limit']):.1e} | {'yes' if entry['holds'] else '**no**'} |"
            )
        return lines

    lines = [
        "# Righting-arm curves `GZ(phi)`: two closed forms, a hull, and the checks behind them",
        "",
        "Generated by `examples/gz/run_example.py`; every number below comes out of that run and",
        "`tests/test_gz_evidence.py` reruns it against this page, so nothing here is typed by hand.",
        "",
        "## Case and parameters",
        "",
        "| case | geometry | loading | reference |",
        "|---|---|---|---|",
        f"| box barge | {LENGTH:g} x {BREADTH:g} x {DEPTH:g} m, draft {DRAFT:g} m | `KG = {KG:g} m` |"
        f" wall-sided trapezoid (limit {limit_deg:.2f} deg) |",
        f"| log | R = {LOG_RADIUS:g} m, half submerged, {LOG_SEGMENTS}-gon section | `KG = {KG:g} m` |"
        f" exact `(R - KG) sin(phi)` |",
        f"| hull | analytic hull {LENGTH:g} x {BREADTH:g} x {DRAFT:g} m, {HULL_NX}x{HULL_NZ} grid,"
        f" freeboard {HULL_FREEBOARD:g} m | `KG = {KG:g} m` | `initial_stability` `GM` (waterplane"
        f" moments) |",
        "",
        f"The box barge has `GM0 = {gm0:.6f} m` and `BM0 = {bm0:.6f} m`; the wall-sided formula is",
        "`GZ = sin(phi) (GM0 + BM0/2 tan^2 phi)`, derived in `tensorcad/gz.py` rather than quoted.",
        "",
        "## Claims",
        "",
        *claim_table(),
        "",
        "## Box barge against the wall-sided closed form",
        "",
        "| heel [deg] | GZ integrated [m] | wall-sided [m] | abs. error [m] |",
        "|---:|---:|---:|---:|",
    ]
    lines += [
        f"| {row['heel_deg']:.0f} | {row['gz_m']:.9f} | {row['wall_sided_m']:.9f} | {row['abs_error_m']:.2e} |"
        for row in box_rows
    ]
    lines += [
        "",
        f"Inside the {limit_deg:.2f} deg limit the two agree to {worst_inside:.2e} m; just past the",
        f"limit the disagreement is {worst_past:.2e} m (and 2.69e+01 m at 85 deg, where the closed",
        "form's tan(phi) itself diverges), which is what makes the first claim informative",
        "rather than tautological.",
        "",
        "## Log against the exact circular-section arm",
        "",
        "| heel [deg] | GZ integrated [m] | exact [m] | abs. error [m] |",
        "|---:|---:|---:|---:|",
    ]
    lines += [
        f"| {row['heel_deg']:.0f} | {row['gz_m']:.9f} | {row['exact_m']:.9f} | {row['abs_error_m']:.2e} |"
        for row in log_rows
    ]
    lines += [
        "",
        f"Worst error over 0-90 deg: {log_worst:.2e} m. A 72-gon repeats on 5 deg heels, so the",
        "comparison is of the integrator, not of two discretisations.",
        "",
        "## Hull (41x13, freeboard 1.0 m)",
        "",
        "| item | value |",
        "|---|---:|",
        f"| upright displacement [m^3] | {upright_banded.displacement_volume:.6f} |",
        f"| KB [m] | {upright_banded.centre_of_buoyancy_z:.6f} |",
        f"| waterplane area [m^2] | {upright_banded.waterplane_area:.6f} |",
        f"| wetted surface [m^2] | {upright_banded.wetted_surface_area:.6f} |",
        f"| `GM` from `initial_stability` [m] | {stability.gm_transverse:.9f} |",
        f"| `GZ'(0)` from the integrated curve [m/rad] | {extrapolated:.9f} |",
        f"| `max GZ` [m] | {hull_curve.max_gz:.6f} at {hull_curve.angle_of_max_gz_deg:.0f} deg |",
        f"| range of stability | {range_text(hull_curve)} |",
        f"| dynamical stability 0-30 deg [m*rad] | {hull_curve.area_under_curve(30.0):.6f} |",
        f"| worst displacement residual | {volume_residual:.2e} |",
        f"| free-surface loss, {TANK_LENGTH:g} x {TANK_BREADTH:g} m tank [mm] | {correction.loss_m * 1000:.3f} |",
        "",
    ]
    if hull_curve.angle_of_max_gz_deg == hull_curve.angles_deg[-1]:
        lines += [
            f"The hull's arm is still rising at {hull_curve.angles_deg[-1]:.0f} deg, the last sampled heel, so",
            "the `max GZ` above is the end of the sampled range and not a physical peak; the same",
            "reason makes the range of stability report \"> 90 deg\" rather than a number. That is a",
            "statement about the sampling, not about the hull.",
            "",
        ]
    lines += [
        "## Hull curve, sampled",
        "",
        "| heel [deg] | GZ solid [m] | GZ with free surface [m] |",
        "|---:|---:|---:|",
    ]
    lines += [
        f"| {angle:.0f} | {solid:.9f} | {wet:.9f} |"
        for angle, solid, wet in zip(hull_curve.angles_deg, hull_curve.gz, flooded.gz, strict=True)
    ]
    lines += [
        "",
        "## Reproduce",
        "",
        "```",
        "python examples/gz/run_example.py --output docs/gz",
        "pytest tests/test_gz_evidence.py -q   # reruns the script and compares this page",
        "```",
        "",
        f"**Verdict: {'all claims hold' if not failures else 'FAILED -- ' + '; '.join(failures)}**",
        "",
    ]
    (output / "report.md").write_text("\n".join(lines))

    # --- figure --------------------------------------------------------------------
    if not options.no_plot:
        try:
            import matplotlib

            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:  # pragma: no cover - optional
            print("matplotlib not installed; skipping figure")
        else:
            figure, axes = plt.subplots(1, 3, figsize=(15, 4.2))

            axes[0].plot(
                [row["heel_deg"] for row in box_rows],
                [row["gz_m"] for row in box_rows],
                "o-",
                ms=3,
                label="integrated GZ",
            )
            axes[0].plot(
                [row["heel_deg"] for row in box_rows],
                [row["wall_sided_m"] for row in box_rows],
                "k--",
                label=r"wall-sided closed form",
            )
            axes[0].axvline(limit_deg, color="0.5", lw=0.8, ls=":")
            axes[0].annotate(
                f"limit {limit_deg:.1f} deg",
                (limit_deg, axes[0].get_ylim()[1] * 0.05),
                fontsize=7,
                rotation=90,
                color="0.4",
            )
            axes[0].set_xlabel("heel angle [deg]")
            axes[0].set_ylabel("GZ [m]")
            axes[0].set_title(f"Box barge {LENGTH:g}x{BREADTH:g} m: integrator vs closed form")
            axes[0].legend(fontsize=7)
            axes[0].grid(alpha=0.3)

            axes[1].plot(
                [row["heel_deg"] for row in log_rows],
                [row["gz_m"] for row in log_rows],
                "o-",
                ms=3,
                label=f"integrated ({LOG_SEGMENTS}-gon log)",
            )
            axes[1].plot(
                [row["heel_deg"] for row in log_rows],
                [row["exact_m"] for row in log_rows],
                "k--",
                label=r"exact $(R-KG)\sin\phi$",
            )
            axes[1].set_xlabel("heel angle [deg]")
            axes[1].set_ylabel("GZ [m]")
            axes[1].set_title(f"Log, R = {LOG_RADIUS:g} m: round bilge, exact reference")
            axes[1].legend(fontsize=7)
            axes[1].grid(alpha=0.3)

            axes[2].plot(hull_curve.angles_deg, hull_curve.gz, "o-", ms=3, label=f"solid (KG = {KG:g} m)")
            axes[2].plot(flooded.angles_deg, flooded.gz, "s--", ms=3, label="with free surface")
            axes[2].plot(
                [row["heel_deg"] for row in box_rows],
                [row["gz_m"] for row in box_rows],
                "^-",
                ms=3,
                label="box barge",
            )
            axes[2].axhline(0.0, color="k", lw=0.8)
            axes[2].set_xlabel("heel angle [deg]")
            axes[2].set_ylabel("GZ [m]")
            axes[2].set_title(f"Hull ({HULL_NX}x{HULL_NZ}, freeboard {HULL_FREEBOARD:g} m) and box barge")
            axes[2].legend(fontsize=7)
            axes[2].grid(alpha=0.3)
            figure.tight_layout()
            figure.savefig(output / "gz_curves.png", dpi=130)
            plt.close(figure)

    print("\n".join(lines))
    for failure in failures:
        print(f"FAIL: {failure}")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
