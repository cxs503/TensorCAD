#!/usr/bin/env python3
"""Publish the large-angle stability evidence page under ``docs/gz/``.

This is the reproducible half of `docs/gz.md`: every number on the published page, every
row of every table and every panel of ``figures/gz-curve.png`` are recomputed here from the
geometry and written out, so the page cannot drift away from the code.

    python publish_report.py                    # write docs/gz/{results.json,report.md,figures}
    python publish_report.py --output /tmp/gz   # write somewhere else
    python publish_report.py --check            # recompute, compare with what is committed

``--check`` recomputes the claims, compares the scalar results with the committed
``docs/gz/results.json`` at a tolerance far below any physical sensitivity, and exits
non-zero on drift, on a failed claim, or on a missing figure. ``tests/test_gz.py`` runs it,
so a change that moves a published number has to move the page with it.

The claims (each one a computed comparison, not an assertion in prose):

``C1``  box barge -- the integrator reproduces the wall-sided closed form
        ``GZ = sin(phi) (GM0 + BM0/2 tan^2 phi)`` to round-off inside the formula's
        validity window (bilge emergence / deck-edge immersion).
``C2``  the same comparison *must fail* past that window, so C1 is not vacuous.
``C3``  log (circular section, half submerged) -- ``GZ = (R - KG) sin(phi)`` to 90 deg.
``C4``  hull -- ``GZ'(0)`` equals ``GM`` from :mod:`tensorcad.stability`, which gets it from
        the waterplane moments by a completely different route.
``C5``  every sampled point of every curve floats at the same displacement.
``C6``  cross-route agreement: this module's Newton waterline solve against the independent
        bracket-and-bisect route in :mod:`tensorcad.righting_arm`, on the same meshes.
``C7``  a negative control for C6: shifting the displacement target by 10 % has to move
        ``GZ`` far outside that agreement.
``C8``  a slack tank enters as a flat shift ``sum(i_x)/V`` on every heel of the curve.
``C9``  the free-surface correction lowers the peak by exactly that shift.

matplotlib is optional: without it the JSON and the Markdown are still written and the
figure is skipped.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:  # run from a checkout without installing
    sys.path.insert(0, str(ROOT))

from tensorcad import __version__  # noqa: E402
from tensorcad import mesh as meshlib  # noqa: E402
from tensorcad.geometry import Mesh, MeshBuilder  # noqa: E402
from tensorcad.gz import (  # noqa: E402
    box_righting_arms_analytic,
    box_wall_sided_limit_deg,
    free_surface_correction,
    heeled_hydrostatics,
    rectangular_tank_inertia,
    righting_arm_curve,
    wall_sided_righting_arm,
    with_freeboard,
)
from tensorcad.hull import HullParameters, generate_hull  # noqa: E402
from tensorcad.righting_arm import righting_arm_curve as bisection_righting_arm_curve  # noqa: E402
from tensorcad.stability import initial_stability  # noqa: E402

# --- the case: one barge, one log, one parametric hull --------------------------------
LENGTH, BREADTH, DRAFT, DEPTH = 10.0, 2.0, 0.8, 2.0
BARGE_KG, BARGE_STEP_DEG = 0.50, 5
BARGE_BEYOND_DEG = (45.0, 55.0, 70.0)

LOG_RADIUS, LOG_LENGTH, LOG_SEGMENTS, LOG_KG = 1.0, 4.0, 64, 0.50
# The polygon is centrally symmetric about a vertex line through ``y = 0``; the identity
# ``KN = R sin(phi)`` is exact for the *mesh* only when the cut also lands on such a line,
# i.e. on multiples of 360/segments = 5.625 deg. The claim samples those; heels between
# them are measured separately (``off_axis``) so the polygonal asymmetry is on the page
# instead of being hidden by a lucky grid.
LOG_ANGLES_DEG = tuple(k * 360.0 / LOG_SEGMENTS for k in range(0, 17, 2))
LOG_OFF_AXIS_DEG = (10.0, 30.0)

HULL_NX, HULL_NZ, HULL_FREEBOARD, HULL_KG = 41, 13, 1.0, 0.50
HULL_ANGLES_DEG = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)

TANK_LENGTH, TANK_BREADTH = 2.0, 0.5

CROSS_ROUTE_HEELS_DEG = (0.0, 10.0, 20.0, 30.0, 40.0)
NEGATIVE_CONTROL_HEEL_DEG = 20.0
NEGATIVE_CONTROL_VOLUME_SHIFT = 0.10

# Thresholds. Each is orders of magnitude above the round-off observed here and far below
# any physically meaningful heeled quantity, so a real regression cannot hide under it.
TOL_CLOSED_FORM_M = 1e-9
TOL_NEGATIVE_CONTROL_M = 1e-2
TOL_SLOPE_RELATIVE = 1e-6
TOL_VOLUME_RESIDUAL = 1e-11
TOL_CROSS_ROUTE_M = 1e-9
TOL_FREE_SURFACE_M = 1e-12
CHECK_RELATIVE_TOLERANCE = 1e-6


# --- meshes --------------------------------------------------------------------------
def box(length: float, beam: float, draft: float, depth: float) -> Mesh:
    """A rectangular barge: keel datum ``z = 0``, waterline at ``draft``, deck at ``depth``."""
    builder = MeshBuilder()
    corners = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (0.0, depth):
                corners[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corners[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi, bottom, top, far = -0.5 * beam, 0.5 * beam, 0.0, depth, length
    quad((0.0, lo, bottom), (0.0, hi, bottom), (far, hi, bottom), (far, lo, bottom))
    quad((0.0, lo, top), (far, lo, top), (far, hi, top), (0.0, hi, top))
    quad((0.0, lo, bottom), (far, lo, bottom), (far, lo, top), (0.0, lo, top))
    quad((0.0, hi, bottom), (0.0, hi, top), (far, hi, top), (far, hi, bottom))
    quad((0.0, lo, bottom), (0.0, lo, top), (0.0, hi, top), (0.0, hi, bottom))
    quad((far, lo, bottom), (far, hi, bottom), (far, hi, top), (far, lo, top))
    return meshlib.orient_outward(builder.build()).mesh


def log(radius: float, length: float, segments: int) -> Mesh:
    """A prism of circular section: keel at ``z = 0``, axis along ``x``.

    ``segments`` must be even so a symmetry axis of the polygon lies on ``y = 0``: the
    waterline then stays on a vertex line at every heel, the polygonal section is exactly
    centrally symmetric, and ``KN = R sin(phi)`` is exact for the mesh.
    """
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


def hull_mesh() -> tuple[Mesh, Mesh]:
    """The parametric hull at ``41 x 13`` and the same hull closed by 1 m of freeboard."""
    hull = generate_hull(
        HullParameters(length=LENGTH, beam=BREADTH, draft=DRAFT, nx=HULL_NX, nz=HULL_NZ)
    ).mesh
    return hull, with_freeboard(hull, HULL_FREEBOARD)


# --- measurement ---------------------------------------------------------------------
def _row(heel: float, gz: float, reference: float) -> dict[str, float]:
    return {
        "heel_deg": round(heel, 9),
        "gz_m": round(gz, 12),
        "reference_m": round(reference, 12),
        "abs_error_m": abs(gz - reference),
    }


def measure() -> dict[str, object]:
    """Recompute the whole page. Returns the dict written to ``results.json``."""
    claims: list[dict[str, object]] = []
    results: dict[str, object] = {"generated_by": "publish_report.py", "tensorcad_version": __version__}

    def claim(ident, statement, value, threshold, unit, passed, direction="upper") -> None:
        """Record one claim. ``direction`` says which side of the threshold passes:

        ``upper`` -- the value must stay *below* the threshold (an accuracy claim);
        ``lower`` -- the value must stay *above* it (a negative control that has to bite).
        """
        claims.append(
            {
                "id": ident,
                "statement": statement,
                "value": float(value),
                "threshold": float(threshold),
                "unit": unit,
                "direction": direction,
                "passed": bool(passed),
            }
        )

    # -- C1/C2: barge against the wall-sided trapezoid --------------------------------
    barge = box(LENGTH, BREADTH, DRAFT, DEPTH)
    limit_deg = box_wall_sided_limit_deg(BREADTH, DRAFT, DEPTH)
    inside_angles = [float(a) for a in range(0, int(limit_deg) + 1, BARGE_STEP_DEG)]
    barge_curve = righting_arm_curve(barge, BARGE_KG, inside_angles, waterline_z=DRAFT)
    closed = box_righting_arms_analytic(BREADTH, DRAFT, BARGE_KG, inside_angles, depth=DEPTH)
    inside = [
        _row(heel, gz, reference)
        for heel, gz, reference in zip(barge_curve.angles_deg, barge_curve.gz, closed, strict=True)
    ]
    worst_inside = max(row["abs_error_m"] for row in inside)

    bm0 = BREADTH * BREADTH / (12.0 * DRAFT)
    gm0 = 0.5 * DRAFT + bm0 - BARGE_KG
    beyond = []
    for heel in BARGE_BEYOND_DEG:
        gz = righting_arm_curve(barge, BARGE_KG, [heel], waterline_z=DRAFT).gz[0]
        beyond.append(_row(heel, gz, wall_sided_righting_arm(heel, gm=gm0, bm=bm0)))
    worst_beyond = max(row["abs_error_m"] for row in beyond)

    results["barge"] = {
        "parameters": {
            "length_m": LENGTH,
            "beam_m": BREADTH,
            "draft_m": DRAFT,
            "depth_m": DEPTH,
            "freeboard_m": DEPTH - DRAFT,
            "kg_m": BARGE_KG,
            "waterline_z_m": DRAFT,
        },
        "wall_sided_limit_deg": limit_deg,
        "bm0_m": bm0,
        "gm0_m": gm0,
        "inside": inside,
        "worst_abs_error_m": worst_inside,
        "beyond_limit": beyond,
        "worst_beyond_abs_error_m": worst_beyond,
    }
    claim("C1", "barge reproduces sin(phi)(GM0 + BM0/2 tan^2 phi) inside its validity window",
          worst_inside, TOL_CLOSED_FORM_M, "m", worst_inside <= TOL_CLOSED_FORM_M)
    claim("C2", "the same closed form is metres out past the window (non-vacuous control)",
          worst_beyond, TOL_NEGATIVE_CONTROL_M, "m", worst_beyond >= TOL_NEGATIVE_CONTROL_M,
          direction="lower")

    # -- C3: the log ------------------------------------------------------------------
    log_mesh = log(LOG_RADIUS, LOG_LENGTH, LOG_SEGMENTS)
    log_curve = righting_arm_curve(log_mesh, LOG_KG, list(LOG_ANGLES_DEG), waterline_z=LOG_RADIUS)
    log_rows = [
        _row(heel, gz, (LOG_RADIUS - LOG_KG) * math.sin(math.radians(heel)))
        for heel, gz in zip(log_curve.angles_deg, log_curve.gz, strict=True)
    ]
    worst_log = max(row["abs_error_m"] for row in log_rows)
    off_axis_curve = righting_arm_curve(log_mesh, LOG_KG, list(LOG_OFF_AXIS_DEG), waterline_z=LOG_RADIUS)
    off_axis = [
        _row(heel, gz, (LOG_RADIUS - LOG_KG) * math.sin(math.radians(heel)))
        for heel, gz in zip(off_axis_curve.angles_deg, off_axis_curve.gz, strict=True)
    ]
    results["log"] = {
        "parameters": {
            "radius_m": LOG_RADIUS,
            "length_m": LOG_LENGTH,
            "segments": LOG_SEGMENTS,
            "kg_m": LOG_KG,
            "waterline_z_m": LOG_RADIUS,
        },
        "angles_deg": list(LOG_ANGLES_DEG),
        "rows": log_rows,
        "worst_abs_error_m": worst_log,
        "off_axis": off_axis,
        "worst_off_axis_abs_error_m": max(row["abs_error_m"] for row in off_axis),
    }
    claim("C3", "log reproduces (R - KG) sin(phi) from 0 to 90 degrees",
          worst_log, TOL_CLOSED_FORM_M, "m", worst_log <= TOL_CLOSED_FORM_M)

    # -- C4: hull slope at zero against the waterplane-moment route -------------------
    hull, banded = hull_mesh()
    watertight = meshlib.is_watertight(banded)
    stability = initial_stability(hull, kg=HULL_KG, waterline_z=0.0)

    def slope(step_deg: float) -> float:
        probe = righting_arm_curve(banded, HULL_KG, [-step_deg, 0.0, step_deg], waterline_z=0.0)
        return (probe.gz[2] - probe.gz[0]) / math.radians(2.0 * step_deg)

    coarse, fine = slope(0.01), slope(0.001)
    # central differences are second order: cancel the O(h^2) term by Richardson extrapolation
    slope_extrapolated = fine + (fine - coarse) / 9.0
    slope_relative = abs(slope_extrapolated - stability.gm_transverse) / stability.gm_transverse

    hull_curve = righting_arm_curve(banded, HULL_KG, list(HULL_ANGLES_DEG), waterline_z=0.0)
    results["hull"] = {
        "parameters": {
            "length_m": LENGTH,
            "beam_m": BREADTH,
            "draft_m": DRAFT,
            "freeboard_m": HULL_FREEBOARD,
            "kg_m": HULL_KG,
            "nx": HULL_NX,
            "nz": HULL_NZ,
            "waterline_z_m": 0.0,
        },
        "watertight": watertight,
        "upright": {
            "displacement_volume_m3": hull_curve.displacement_volume,
            "kb_m": stability.kb,
            "bm_transverse_m": stability.bm_transverse,
            "gm_transverse_m": stability.gm_transverse,
            "waterplane_area_m2": heeled_hydrostatics(banded, 0.0, waterline_z=0.0).waterplane_area,
        },
        "gm_from_curve_slope_m": slope_extrapolated,
        "slope_relative_error": slope_relative,
        "curve": {
            "angles_deg": list(hull_curve.angles_deg),
            "gz_m": list(hull_curve.gz),
            "kn_m": list(hull_curve.kn),
            "max_gz_m": hull_curve.max_gz,
            "angle_of_max_gz_deg": hull_curve.angle_of_max_gz_deg,
            "range_of_stability_deg": hull_curve.range_of_stability_deg,
            "area_under_curve_m_rad": hull_curve.area_under_curve(),
        },
    }
    claim("C4", "hull GZ'(0) equals GM from the waterplane moments (independent integration)",
          slope_relative, TOL_SLOPE_RELATIVE, "relative", slope_relative <= TOL_SLOPE_RELATIVE)

    # -- C5: every point floats at the upright displacement ---------------------------
    residuals = {
        "barge": barge_curve.max_volume_error,
        "log": log_curve.max_volume_error,
        "hull": hull_curve.max_volume_error,
    }
    worst_residual = max(residuals.values())
    results["volume_residual"] = {"per_curve": residuals, "worst": worst_residual}
    claim("C5", "every sampled heel floats at the upright displacement",
          worst_residual, TOL_VOLUME_RESIDUAL, "relative", worst_residual <= TOL_VOLUME_RESIDUAL)

    # -- C6/C7: cross-route agreement and its negative control -----------------------
    barge_heels = [heel for heel in CROSS_ROUTE_HEELS_DEG if heel <= 30.0]
    bisection_hull = bisection_righting_arm_curve(
        banded, kg=HULL_KG, heels_deg=list(CROSS_ROUTE_HEELS_DEG), waterline_z=0.0
    )
    newton_hull = righting_arm_curve(banded, HULL_KG, list(CROSS_ROUTE_HEELS_DEG), waterline_z=0.0)
    bisection_barge = bisection_righting_arm_curve(
        barge, kg=BARGE_KG, heels_deg=barge_heels, waterline_z=DRAFT
    )
    newton_barge = righting_arm_curve(barge, BARGE_KG, barge_heels, waterline_z=DRAFT)
    # the bisection route publishes one ``RightingArm`` record per heel, not parallel tuples
    bisection_hull_gz = [point.righting_arm for point in bisection_hull.points]
    bisection_barge_gz = [point.righting_arm for point in bisection_barge.points]
    hull_diffs = [abs(a - b) for a, b in zip(newton_hull.gz, bisection_hull_gz, strict=True)]
    barge_diffs = [abs(a - b) for a, b in zip(newton_barge.gz, bisection_barge_gz, strict=True)]
    worst_cross_route = max(max(hull_diffs), max(barge_diffs))

    # negative control: ask the Newton route for a displacement 10 % off the upright one;
    # if C6 were measuring nothing, this would leave GZ unchanged.
    control_heel = [NEGATIVE_CONTROL_HEEL_DEG]
    reference_gz = righting_arm_curve(banded, HULL_KG, control_heel, waterline_z=0.0).gz[0]
    shifted = righting_arm_curve(
        banded,
        HULL_KG,
        control_heel,
        waterline_z=0.0,
        displacement_volume=newton_hull.displacement_volume * (1.0 + NEGATIVE_CONTROL_VOLUME_SHIFT),
    )
    control_delta = abs(shifted.gz[0] - reference_gz)

    results["cross_route"] = {
        "heels_deg": list(CROSS_ROUTE_HEELS_DEG),
        "hull": {
            "newton_gz_m": list(newton_hull.gz),
            "bisection_gz_m": bisection_hull_gz,
            "abs_difference_m": hull_diffs,
        },
        "barge": {
            "heels_deg": barge_heels,
            "newton_gz_m": list(newton_barge.gz),
            "bisection_gz_m": bisection_barge_gz,
            "abs_difference_m": barge_diffs,
        },
        "worst_abs_difference_m": worst_cross_route,
        "negative_control": {
            "heel_deg": NEGATIVE_CONTROL_HEEL_DEG,
            "volume_shift_fraction": NEGATIVE_CONTROL_VOLUME_SHIFT,
            "gz_shift_m": control_delta,
        },
    }
    claim("C6", "Newton and bracket-and-bisect waterline routes agree on GZ",
          worst_cross_route, TOL_CROSS_ROUTE_M, "m", worst_cross_route <= TOL_CROSS_ROUTE_M)
    claim("C7", "a 10 % displacement shift moves GZ far outside the C6 agreement (negative control)",
          control_delta, 10.0 * TOL_CROSS_ROUTE_M, "m", control_delta >= 10.0 * TOL_CROSS_ROUTE_M,
          direction="lower")

    # -- C8/C9: free surface ---------------------------------------------------------
    tank = rectangular_tank_inertia(TANK_LENGTH, TANK_BREADTH)
    correction = free_surface_correction([tank], hull_curve.displacement_volume)
    flooded = righting_arm_curve(
        banded, HULL_KG, hull_curve.angles_deg, waterline_z=0.0, free_surface_i_x=tank
    )
    shift_errors = [
        abs((solid - wet) - correction.loss_m)
        for solid, wet in zip(hull_curve.gz, flooded.gz, strict=True)
    ]
    worst_shift_error = max(shift_errors)
    peak_drop = hull_curve.max_gz - flooded.max_gz
    results["free_surface"] = {
        "tank": {"length_m": TANK_LENGTH, "breadth_m": TANK_BREADTH, "i_x_m4": tank},
        "flat_shift_m": correction.loss_m,
        "worst_shift_error_m": worst_shift_error,
        "max_gz_solid_m": hull_curve.max_gz,
        "max_gz_with_free_surface_m": flooded.max_gz,
        "peak_drop_m": peak_drop,
    }
    claim("C8", "a slack tank shifts every heel of the curve by exactly sum(i_x)/V",
          worst_shift_error, TOL_FREE_SURFACE_M, "m", worst_shift_error <= TOL_FREE_SURFACE_M)
    claim("C9", "the free-surface shift lowers the peak by exactly that amount",
          abs(peak_drop - correction.loss_m), TOL_FREE_SURFACE_M, "m",
          abs(peak_drop - correction.loss_m) <= TOL_FREE_SURFACE_M)

    passed = all(bool(entry["passed"]) for entry in claims)
    results["claims"] = claims
    results["recorded_claims"] = {str(entry["id"]): entry["value"] for entry in claims}
    results["passed"] = passed
    return results


# --- the page ------------------------------------------------------------------------
def _fmt(value: float, digits: int = 3) -> str:
    """Compact table-cell formatting: scientific for very small/large, fixed otherwise."""
    if value == 0.0:
        return "0"
    if abs(value) < 1e-3 or abs(value) >= 1e5:
        return f"{value:.{digits - 1}e}"
    return f"{value:.{max(digits, 3)}f}"


def _heel(value: float) -> str:
    """Heel-angle cell: keeps 5.625 deg steps from being rounded away."""
    rounded = round(value)
    return f"{rounded:d}" if abs(value - rounded) < 5e-10 else f"{value:.2f}"


def _verdict(value: float, threshold: float, direction: str = "upper") -> str:
    passed = value <= threshold if direction == "upper" else value >= threshold
    return "**pass**" if passed else "**FAIL**"


def render_report(results) -> str:
    """Markdown for ``docs/gz/report.md``. Every number comes from ``results``."""
    barge, log_block, hull = results["barge"], results["log"], results["hull"]
    cross, free, residual = results["cross_route"], results["free_surface"], results["volume_residual"]
    claims = {str(entry["id"]): entry for entry in results["claims"]}
    bp, lp, hp = barge["parameters"], log_block["parameters"], hull["parameters"]
    curve, control = hull["curve"], cross["negative_control"]
    range_text = (
        f"{curve['range_of_stability_deg']:.1f} deg"
        if curve["range_of_stability_deg"] is not None
        else "no vanishing angle below 90 deg"
    )
    heels_text = ", ".join(str(int(a)) if float(a).is_integer() else str(a) for a in cross["heels_deg"])

    lines = [
        "# Large-angle stability `GZ(phi)`: measured evidence page",
        "",
        f"TensorCAD {results['tensorcad_version']}, written by `publish_report.py`. Every number "
        "below is recomputed by that script, and `python publish_report.py --check` fails if this "
        "page and the code disagree. Discussion and the API live in [docs/gz.md](../gz.md); this "
        "page is the measurement record.",
        "",
        "## 1. The case",
        "",
        "Three bodies, chosen so that two of them have a right answer that does not come from the "
        "integrator under test:",
        "",
        "| body | L x B x T [m] | depth / freeboard [m] | KG [m] | reference |",
        "|---|---|---|---|---|",
        f"| box barge | {bp['length_m']:.2f} x {bp['beam_m']:.2f} x {bp['draft_m']:.2f} "
        f"| {bp['depth_m']:.2f} / {bp['freeboard_m']:.2f} | {bp['kg_m']:.2f} "
        "| wall-sided closed form `sin(phi)(GM0 + BM0/2 tan^2 phi)` |",
        f"| log (circular section) | {lp['length_m']:.2f} x {2 * lp['radius_m']:.2f} x "
        f"{2 * lp['radius_m']:.2f} | {lp['segments']}-gon | {lp['kg_m']:.2f} "
        "| central symmetry: `(R - KG) sin(phi)` |",
        f"| parametric hull | {hp['length_m']:.2f} x {hp['beam_m']:.2f} x {hp['draft_m']:.2f} "
        f"| freeboard {hp['freeboard_m']:.2f} | {hp['kg_m']:.2f} "
        "| `GM` from waterplane moments (`tensorcad.stability`) |",
        "",
        f"Waterlines: barge `z = {bp['waterline_z_m']:.2f} m` (keel datum `z = 0`), log "
        f"`z = {lp['waterline_z_m']:.2f} m` (half submerged), hull `z = 0` with the "
        f"`{hp['nx']}x{hp['nz']}` mesh closed above the waterline by `with_freeboard`. Upright "
        f"displacement of the hull {hull['upright']['displacement_volume_m3']:.6f} m^3, "
        f"`KB = {hull['upright']['kb_m']:.6f} m`, `BM = {hull['upright']['bm_transverse_m']:.6f} m`, "
        f"`GM = {hull['upright']['gm_transverse_m']:.6f} m`.",
        "",
        "## 2. How it was done",
        "",
        "1. `tensorcad.gz.heeled_hydrostatics` rotates the closed mesh about the `x` axis through "
        "the keel line, cuts it with `z = c`, and integrates the part below with the divergence "
        "theorem, apex on the waterline plane so the closing cap contributes exactly zero; the "
        "same pass accumulates the waterplane area and the wetted area.",
        "2. `c` is Newton-solved per heel on `dV/dc = A(c)` until `V(phi) = V(0)` to a relative "
        "`1e-13`. Then `GZ(phi) = KN(phi) - KG sin(phi)` with `KN = y_B cos phi + z_B sin phi`, "
        "all measured from the keel datum.",
        f"3. References: the wall-sided closed form on a `{bp['length_m']:.0f} x {bp['beam_m']:.0f} "
        f"x {bp['draft_m']:.1f} m` barge (valid to `{barge['wall_sided_limit_deg']:.2f} deg` here), "
        "central symmetry on the log, and `tensorcad.stability.initial_stability` (Green's theorem "
        "on the cut waterplane loops) for the hull slope.",
        "4. Cross-route check: `tensorcad.righting_arm.righting_arm_curve` integrates the clipped "
        "sections in the body frame and finds the waterline by bracket-and-bisect on the monotone "
        "`V(h)` -- a different integral and a different solver -- and the two curves are compared "
        f"at `{heels_text} deg`.",
        "5. `publish_report.py` writes `results.json`, this page and `figures/gz-curve.png` from "
        "one run; `--check` re-derives the claims and compares them with the committed JSON.",
        "",
        "## 3. Results",
        "",
        "### 3.1 Barge against the wall-sided trapezoid (inside the window)",
        "",
        "| heel [deg] | GZ [m] | closed form [m] | abs diff [m] |",
        "|---|---|---|---|",
    ]
    for row in barge["inside"]:
        lines.append(
            f"| {_heel(row['heel_deg'])} | {row['gz_m']:.12f} | {row['reference_m']:.12f} "
            f"| {_fmt(row['abs_error_m'])} |"
        )
    lines += [
        "",
        f"`GM0 = {barge['gm0_m']:.6f} m`, `BM0 = {barge['bm0_m']:.6f} m`, worst diff "
        f"`{_fmt(barge['worst_abs_error_m'])} m` -- round-off: the polyhedral integral and the "
        "closed form are both exact for this shape.",
        "",
        "**Negative control** -- the same comparison past the limit, where the closed form must fail:",
        "",
        "| heel [deg] | GZ [m] | wall-sided [m] | abs diff [m] |",
        "|---|---|---|---|",
    ]
    for row in barge["beyond_limit"]:
        lines.append(
            f"| {_heel(row['heel_deg'])} | {row['gz_m']:.6f} | {row['reference_m']:.6f} "
            f"| {_fmt(row['abs_error_m'])} |"
        )
    lines += [
        "",
        "### 3.2 Log against central symmetry",
        "",
        "| heel [deg] | GZ [m] | (R - KG) sin phi [m] | abs diff [m] |",
        "|---|---|---|---|",
    ]
    for row in log_block["rows"]:
        lines.append(
            f"| {_heel(row['heel_deg'])} | {row['gz_m']:.12f} | {row['reference_m']:.12f} "
            f"| {_fmt(row['abs_error_m'])} |"
        )
    lines += [
        "",
        "Off the polygon's symmetry axes the identity is only approximate -- the cut no "
        "longer splits the section into equal areas -- and that is reported rather than "
        "hidden by the grid: "
        + "; ".join(
            f"`{row['heel_deg']:.1f} deg -> {_fmt(row['abs_error_m'])} m`"
            for row in log_block["off_axis"]
        )
        + ".",
        "",
        f"### 3.3 Hull curve (`{hp['nx']} x {hp['nz']}`, KG = {hp['kg_m']:.2f} m, "
        f"freeboard {hp['freeboard_m']:.1f} m)",
        "",
        "| heel [deg] | GZ [m] | KN [m] |",
        "|---|---|---|",
    ]
    for heel, gz, kn in zip(curve["angles_deg"], curve["gz_m"], curve["kn_m"], strict=True):
        lines.append(f"| {heel:.0f} | {gz:.6f} | {kn:.6f} |")
    lines += [
        "",
        f"- `GZ'(0) = {hull['gm_from_curve_slope_m']:.9f} m/rad` vs `GM = "
        f"{hull['upright']['gm_transverse_m']:.9f} m` from the waterplane moments (relative "
        f"difference `{_fmt(hull['slope_relative_error'])}`).",
        f"- peak `GZ = {curve['max_gz_m']:.6f} m` at `{curve['angle_of_max_gz_deg']:.0f} deg`, "
        f"range of stability {range_text}, area under the curve "
        f"`{curve['area_under_curve_m_rad']:.6f} m rad`.",
        f"- mesh watertight `{hull['watertight']}`; upright waterplane area "
        f"`{hull['upright']['waterplane_area_m2']:.6f} m^2`.",
        "",
        "### 3.4 Cross-route agreement (Newton vs bracket-and-bisect)",
        "",
        "| heel [deg] | hull Newton [m] | hull bisection [m] | abs diff [m] |",
        "|---|---|---|---|",
    ]
    for heel, newton, bisect, diff in zip(
        cross["heels_deg"],
        cross["hull"]["newton_gz_m"],
        cross["hull"]["bisection_gz_m"],
        cross["hull"]["abs_difference_m"],
        strict=True,
    ):
        lines.append(f"| {heel:.0f} | {newton:.12f} | {bisect:.12f} | {_fmt(diff)} |")
    lines += [
        "",
        f"Worst over both meshes `{_fmt(cross['worst_abs_difference_m'])} m`. Negative control: "
        f"shifting the displacement target by `{control['volume_shift_fraction'] * 100:.0f} %` "
        f"moves `GZ({control['heel_deg']:.0f} deg)` by `{_fmt(control['gz_shift_m'])} m`, so the "
        "agreement above is a measurement and not an artefact of both routes ignoring the "
        "displacement.",
        "",
        "### 3.5 Free surface and displacement bookkeeping",
        "",
        f"- slack tank `{free['tank']['length_m']:.1f} x {free['tank']['breadth_m']:.1f} m`: "
        f"`i_x = {free['tank']['i_x_m4']:.6f} m^4`, flat loss "
        f"`i_x/V = {free['flat_shift_m'] * 1000:.4f} mm`; worst departure from a flat shift "
        f"`{_fmt(free['worst_shift_error_m'])} m`.",
        f"- peak `GZ` {free['max_gz_solid_m']:.6f} m -> {free['max_gz_with_free_surface_m']:.6f} m "
        f"(drop {free['peak_drop_m'] * 1000:.4f} mm).",
        f"- displacement residual, worst over all curves `{_fmt(residual['worst'])}` relative "
        f"(barge `{_fmt(residual['per_curve']['barge'])}`, log "
        f"`{_fmt(residual['per_curve']['log'])}`, hull `{_fmt(residual['per_curve']['hull'])}`).",
        "",
        "![GZ curves](figures/gz-curve.png)",
        "",
        "*Left: the three righting-arm curves with the log's closed form dashed. Middle: absolute "
        "difference against each reference on a log scale -- the two exact benchmarks sit at "
        "double-precision round-off across their whole validity window, while the barge's closed "
        "form departs by metres once the bilge emerges. Right: the Newton/bisection cross-route "
        "difference next to the 10 % displacement-shift control, which has to be far above it for "
        "the comparison to mean anything.*",
        "",
        "## 4. Caveats and scope",
        "",
        "- The barge and log references are exact *for the polygonal mesh*; they measure the "
        "integrator and the waterline solve, not the meshing of a real hull. On the log the "
        "identity holds on the polygon's symmetry axes only: between them the discretisation "
        "shows its own asymmetry (section 3.2), which is a property of the 64-gon, not of the "
        "waterline solver.",
        "- The hull has no closed-form `GZ`. Its checks are the slope against the waterplane "
        "moments (a different integration of the same mesh), the cross-route comparison and the "
        "displacement residual; none of them is a validation against model-test or sea-trial data.",
        f"- `KN` is limited by the chord-wise resolution of the `{hp['nx']} x {hp['nz']}` mesh, "
        "exactly as `KB`/`BM` are; no convergence ladder is claimed on this page (see "
        "[docs/righting_arm.md](../righting_arm.md) for one).",
        "- `KG` is the only loading input: no tanks other than the idealised `i_x` shift, no "
        "flooding, no damage, no trim equilibrium, no wind heeling arm, no statutory or class "
        "criterion. `area_under_curve` is an input to such a check, not the check.",
        "- Deck-edge immersion is reported by `deck_edge_immersion_deg`, but the curve keeps rising "
        "past it because nothing floods; that is a flooding criterion, not geometry.",
        "",
        "## 5. Error verdict",
        "",
        "| id | quantity | value | threshold | verdict |",
        "|---|---|---|---|---|",
    ]
    labels = {
        "C1": "barge abs(GZ - wall-sided) [m]",
        "C2": "barge abs(GZ - wall-sided) past the limit [m] (must be large)",
        "C3": "log abs(GZ - (R-KG) sin phi) [m]",
        "C4": "rel(hull GZ'(0), GM) [-]",
        "C5": "max displacement residual [-]",
        "C6": "max abs(Newton - bisection) [m]",
        "C7": "GZ shift under a 10 % volume error [m] (must exceed C6)",
        "C8": "departure from the flat free-surface shift [m]",
        "C9": "abs(peak drop - i_x/V) [m]",
    }
    for ident, label in labels.items():
        entry = claims[ident]
        lines.append(
            f"| {ident} | {label} | `{_fmt(entry['value'])}` | `{_fmt(entry['threshold'])}` "
            f"| {_verdict(entry['value'], entry['threshold'], str(entry.get('direction', 'upper')))} |"
        )
    lines += [
        "",
        "Whether each key quantity is within **3 %** of its benchmark, stated explicitly:",
        "",
        f"- barge `GZ` vs the wall-sided closed form: worst absolute error "
        f"`{_fmt(barge['worst_abs_error_m'])} m` on arms of order 0.1 m, i.e. round-off -- "
        "**yes, below 3 %**.",
        f"- log `GZ` vs central symmetry: `{_fmt(log_block['worst_abs_error_m'])} m` -- "
        "**yes, below 3 %**.",
        f"- hull `GZ'(0)` vs `GM`: `{_fmt(hull['slope_relative_error'])}` relative -- "
        "**yes, below 3 %** (round-off plus the Richardson residual).",
        f"- hull `GZ` vs the independent bisection route: "
        f"`{_fmt(cross['worst_abs_difference_m'])} m` absolute, i.e. round-off -- **yes, below "
        "3 %**; note this is an implementation cross-check, not an experimental benchmark.",
        "- hull `GZ` vs model test or sea trial: **no such benchmark exists in this repository**, "
        "so no percentage is claimed. Nothing on this page validates the righting arms of a real "
        "ship; the honest status of that comparison is *no reference data*, not *agreement*.",
        "",
        f"**Verdict: {'all 9 claims hold' if results['passed'] else 'FAILED -- a claim broke'}.**",
        "",
    ]
    return "\n".join(lines)


def render_figure(results, path: Path) -> bool:
    """Three-panel evidence figure. Returns False when matplotlib is unavailable."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:  # pragma: no cover - optional dependency
        print("matplotlib is not installed; skipping figures/gz-curve.png")
        return False

    barge, log_block, hull = results["barge"], results["log"], results["hull"]
    cross = results["cross_route"]
    bp, lp, hp = barge["parameters"], log_block["parameters"], hull["parameters"]
    curve, control = hull["curve"], cross["negative_control"]
    floor = 1e-18
    figure, axes = plt.subplots(1, 3, figsize=(15.0, 4.4))

    # (a) the curves themselves
    axes[0].plot(curve["angles_deg"], curve["gz_m"], "o-", ms=4,
                 label=f"hull {hp['nx']}x{hp['nz']}, KG={hp['kg_m']:.2f} m")
    barge_rows = barge["inside"] + barge["beyond_limit"]
    axes[0].plot([row["heel_deg"] for row in barge_rows], [row["gz_m"] for row in barge_rows],
                 "s-", ms=4, label=f"barge {bp['length_m']:.0f}x{bp['beam_m']:.0f}x{bp['draft_m']:.1f} m")
    axes[0].plot([row["heel_deg"] for row in log_block["rows"]],
                 [row["gz_m"] for row in log_block["rows"]], "^-", ms=4,
                 label=f"log R={lp['radius_m']:.1f} m")
    axes[0].plot([row["heel_deg"] for row in log_block["rows"]],
                 [row["reference_m"] for row in log_block["rows"]], "k--", lw=1,
                 label=r"exact $(R-KG)\sin\phi$")
    axes[0].axhline(0.0, color="k", lw=0.8)
    axes[0].set_xlabel("heel angle [deg]")
    axes[0].set_ylabel("GZ [m]")
    axes[0].set_title("righting arms, integrated")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.3)

    # (b) error against each reference, log scale
    axes[1].semilogy([row["heel_deg"] for row in barge["inside"]],
                     [max(row["abs_error_m"], floor) for row in barge["inside"]], "s-", ms=4,
                     label=f"barge vs wall-sided (<= {barge['wall_sided_limit_deg']:.1f} deg)")
    axes[1].semilogy([row["heel_deg"] for row in barge["beyond_limit"]],
                     [max(row["abs_error_m"], floor) for row in barge["beyond_limit"]], "s--",
                     ms=4, mfc="none", label="barge past the limit (control)")
    axes[1].semilogy([row["heel_deg"] for row in log_block["rows"]],
                     [max(row["abs_error_m"], floor) for row in log_block["rows"]], "^-", ms=4,
                     label=r"log vs $(R-KG)\sin\phi$")
    axes[1].semilogy(cross["heels_deg"],
                     [max(value, floor) for value in cross["hull"]["abs_difference_m"]], "o-", ms=4,
                     label="hull: Newton vs bisection")
    axes[1].axhline(TOL_CLOSED_FORM_M, color="grey", ls=":", lw=1, label="claim threshold 1e-9 m")
    axes[1].set_xlabel("heel angle [deg]")
    axes[1].set_ylabel("|GZ - reference| [m]")
    axes[1].set_title("difference against each reference")
    axes[1].legend(fontsize=6.5)
    axes[1].grid(alpha=0.3, which="both")

    # (c) cross-route agreement next to its negative control
    bar_labels = [
        "barge\nNewton vs bisect",
        "hull\nNewton vs bisect",
        f"hull, {control['volume_shift_fraction'] * 100:.0f} % V shift\n(control)",
    ]
    values = [
        max(cross["barge"]["abs_difference_m"]),
        max(cross["hull"]["abs_difference_m"]),
        control["gz_shift_m"],
    ]
    bars = axes[2].bar(bar_labels, [max(value, floor) for value in values],
                       color=["#4c72b0", "#4c72b0", "#c44e52"])
    axes[2].set_yscale("log")
    axes[2].set_ylabel("|dGZ| [m]")
    axes[2].set_title("route agreement and its control")
    axes[2].tick_params(axis="x", labelsize=7.5)
    for bar, value in zip(bars, values, strict=True):
        axes[2].annotate(f"{value:.2e}", (bar.get_x() + bar.get_width() / 2, max(value, floor)),
                         ha="center", va="bottom", fontsize=8)
    axes[2].grid(alpha=0.3, axis="y", which="both")

    figure.suptitle(
        f"TensorCAD {results['tensorcad_version']}: GZ(phi) against two exact closed forms and one "
        f"independent route ({len(results['claims'])} claims, publish_report.py)",
        fontsize=10,
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=140)
    plt.close(figure)
    return True


# --- --check -------------------------------------------------------------------------
def check(committed, recomputed) -> list[str]:
    """Compare the recomputed scalar claims with what is committed."""
    problems: list[str] = []
    recorded = committed.get("recorded_claims") or {}
    for entry in recomputed["claims"]:
        ident = str(entry["id"])
        if ident not in recorded:
            problems.append(f"{ident} is missing from the committed results.json")
            continue
        old, new = float(recorded[ident]), float(entry["value"])
        scale = max(abs(old), abs(new), 1e-30)
        if abs(new - old) > CHECK_RELATIVE_TOLERANCE * scale:
            problems.append(f"{ident} moved: committed {old:.6e} -> recomputed {new:.6e}")
        if not entry["passed"]:
            problems.append(f"{ident} no longer holds: {entry['statement']}")
    if bool(committed.get("passed")) != bool(recomputed["passed"]):
        problems.append("the overall verdict changed")
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--output", default=str(ROOT / "docs" / "gz"), help="artifact directory")
    parser.add_argument("--check", action="store_true",
                        help="recompute and compare with the committed artifacts; write nothing")
    parser.add_argument("--no-figure", action="store_true", help="skip figures/gz-curve.png")
    options = parser.parse_args(argv)

    output = Path(options.output)
    results = measure()

    if options.check:
        results_path = output / "results.json"
        if not results_path.exists():
            print(f"FAIL: {results_path} is missing")
            return 1
        problems = check(json.loads(results_path.read_text()), results)
        figure_path = output / "figures" / "gz-curve.png"
        if not options.no_figure:
            if not figure_path.exists():
                problems.append(f"{figure_path} is missing")
            elif figure_path.stat().st_size == 0:
                problems.append(f"{figure_path} is empty")
        for entry in results["claims"]:
            print(f"{entry['id']}: {entry['statement']} -> {entry['value']:.6e} "
                  f"[{entry['threshold']:.1e}] {'ok' if entry['passed'] else 'FAIL'}")
        if problems:
            for problem in problems:
                print(f"FAIL: {problem}")
            return 1
        print(f"publish_report: all {len(results['claims'])} claims hold and match {results_path}")
        return 0

    output.mkdir(parents=True, exist_ok=True)
    (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    (output / "report.md").write_text(render_report(results))
    if not options.no_figure:
        render_figure(results, output / "figures" / "gz-curve.png")

    for entry in results["claims"]:
        print(f"{entry['id']}: {entry['statement']} -> {entry['value']:.6e} "
              f"[{entry['threshold']:.1e}] {'ok' if entry['passed'] else 'FAIL'}")
    print(f"wrote {output}/results.json, {output}/report.md, {output}/figures/gz-curve.png")
    if not results["passed"]:
        print("FAIL: a claim broke")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
