"""Large-angle stability measured against a closed form and an independent section solve.

Nothing on this page is quoted from another run: the script writes ``results.json``,
``report.md`` and the figures, and **exits non-zero if one of its own claims stops
holding**, so ``docs/righting_arm.md`` cannot drift away from the code.

The verification ladder, in the order it is run:

1. **A barge**, where the righting arm has a closed form. For a rectangular body whose
   waterline still spans the full beam, ``GZ = sin(phi)(GM + BM tan(phi)**2 / 2)`` is
   exact, so the mesh integration has to reproduce it to round-off. That pins the
   geometry, the draft solve and the lever at once.
2. **The same barge past the closed form's limit**, where the waterline has left one of
   the sides. The reference there is an independent two-dimensional computation: the
   section of a prism is a rectangle, so the heeled problem reduces exactly to clipping a
   polygon and taking its shoelace centroid -- a different dimension and a different
   integral from the surface flux the module uses.
3. **The loading identities**, which hold for any body: raising ``KG`` by ``d`` lowers the
   arm by exactly ``d sin(phi)``, ``GZ(0) = 0``, and ``GZ`` is odd in the heel (the body is
   symmetric about the centreline).
4. **A hull with topsides**, where no closed form exists. Two things are checked instead:
   the slope of the arm at zero heel against ``KB + I_T / V - KG``, assembled from the
   *waterplane* integrator of ``tensorcad.stability`` (Green's theorem on a plane section,
   no shared code with the heeled-displacement integration), and the order of convergence
   of the arm under mesh refinement.
5. **The shipped hull family**, which is capped at its design waterline and therefore has
   no freeboard. The honest answer is an exception with the reason, and the script asserts
   that it is raised rather than printing a table of zeros.

Run it with ``python examples/righting_arm/run_example.py --output docs/righting_arm``.
matplotlib is optional: the JSON and the report are written either way, the PNGs are
skipped.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad import __version__
from tensorcad.geometry import Mesh, MeshBuilder
from tensorcad.hull import HullParameters, generate_hull
from tensorcad.righting_arm import (
    RightingArmError,
    box_reference_arm,
    immersed_volume,
    righting_arm,
    righting_arms,
    wall_sided_validity,
    wall_sided_vanishing_angle,
    zero_crossings,
)
from tensorcad.stability import waterplane_moments

#: The barge. ``depth = 3 * draft`` so the closed form's limit is bilge emergence at
#: ``atan(2 * draft / beam)`` = 45 deg, with the deck edge still 45 deg of heel away.
LENGTH, BEAM, DRAFT, DEPTH = 10.0, 4.0, 2.0, 6.0
BARGE_KG = 0.9
BARGE_UNSTABLE_KG = 1.9
BARGE_CURVE_STEP = 1.0

#: The hull with topsides: generated with a 1.6 m draft and floated at z = -0.9, so it
#: has 0.7 m immersed and 0.9 m of freeboard at the centreline.
HULL_LENGTH, HULL_BEAM, HULL_GENERATED_DRAFT, HULL_WATERLINE = 10.0, 2.0, 1.6, -0.9
HULL_KG = 0.30
HULL_KG_RAISED = 0.45
HULL_ANGLES = tuple(float(angle) for angle in range(0, 91, 5))
LADDER = ((41, 13), (81, 25), (161, 41))
LADDER_HEEL = 20.0

#: Selected heels for the tables: inside the closed form's regime, at its edge, and past it.
TABLE_ANGLES = (2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 44.0, 45.0, 46.0, 55.0, 60.0, 75.0, 90.0)


# ---------------------------------------------------------------------------
# reference geometries
# ---------------------------------------------------------------------------
def box(length: float, beam: float, depth: float, draft: float) -> Mesh:
    """Closed rectangular prism: keel on ``z = -draft``, top on ``z = depth - draft``."""
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * beam, 0.5 * beam):
            for z in (-draft, depth - draft):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0, p1, p2, p3) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi = -0.5 * beam, 0.5 * beam
    bottom, top = -draft, depth - draft
    quad((0.0, lo, bottom), (0.0, hi, bottom), (length, hi, bottom), (length, lo, bottom))
    quad((0.0, lo, top), (length, lo, top), (length, hi, top), (0.0, hi, top))
    quad((0.0, lo, bottom), (length, lo, bottom), (length, lo, top), (0.0, lo, top))
    quad((0.0, hi, bottom), (0.0, hi, top), (length, hi, top), (length, hi, bottom))
    quad((0.0, lo, bottom), (0.0, lo, top), (0.0, hi, top), (0.0, hi, bottom))
    quad((length, lo, bottom), (length, hi, bottom), (length, hi, top), (length, lo, top))
    return builder.build()


def flared_hull(nx: int, nz: int) -> Mesh:
    parameters = HullParameters(
        length=HULL_LENGTH, beam=HULL_BEAM, draft=HULL_GENERATED_DRAFT, nx=nx, nz=nz
    )
    return generate_hull(parameters).mesh


# ---------------------------------------------------------------------------
# independent two-dimensional reference for a prism
# ---------------------------------------------------------------------------
def _clip_below(polygon, waterline):
    clipped = []
    count = len(polygon)
    for index in range(count):
        p = polygon[index]
        q = polygon[(index + 1) % count]
        dp = p[1] - waterline
        dq = q[1] - waterline
        if dp <= 0.0:
            clipped.append(p)
        if (dp < 0.0 < dq) or (dq < 0.0 < dp):
            t = dp / (dp - dq)
            clipped.append((p[0] + t * (q[0] - p[0]), waterline))
    return clipped


def _polygon_area_and_centroid(polygon):
    area2 = cx = cy = 0.0
    count = len(polygon)
    for index in range(count):
        x1, y1 = polygon[index]
        x2, y2 = polygon[(index + 1) % count]
        cross = x1 * y2 - x2 * y1
        area2 += cross
        cx += (x1 + x2) * cross
        cy += (y1 + y2) * cross
    area = 0.5 * area2
    return abs(area), cx / (6.0 * area), cy / (6.0 * area)


def section_reference_arm(heel_deg, *, beam, draft, depth, kg):
    """``GZ`` of a rectangular prism by clipping its section in 2D (shoelace moments).

    Independent of the module under test in every step that matters: the region is a
    polygon in the plane of the section rather than a surface in space, the moments come
    from Green's theorem in 2D rather than from a divergence theorem in 3D, and the draft
    is solved here again rather than borrowed.
    """
    angle = math.radians(heel_deg)
    cos, sin = math.cos(angle), math.sin(angle)
    corners = [
        (-0.5 * beam, -draft),
        (0.5 * beam, -draft),
        (0.5 * beam, depth - draft),
        (-0.5 * beam, depth - draft),
    ]
    rotated = [(y * cos + z * sin, -y * sin + z * cos) for y, z in corners]
    target = beam * draft
    lower = min(point[1] for point in rotated)
    upper = max(point[1] for point in rotated)
    polygon = rotated
    for _ in range(200):
        waterline = 0.5 * (lower + upper)
        clipped = _clip_below(rotated, waterline)
        if len(clipped) < 3:
            break
        area, _, _ = _polygon_area_and_centroid(clipped)
        polygon = clipped
        if abs(area - target) <= 1e-15 * target:
            break
        if area > target:
            upper = waterline
        else:
            lower = waterline
    _, y_bar, _ = _polygon_area_and_centroid(polygon)
    return y_bar - (-draft + kg) * sin


# ---------------------------------------------------------------------------
# the two bodies
# ---------------------------------------------------------------------------
def measure_barge() -> dict[str, object]:
    mesh = box(LENGTH, BEAM, DEPTH, DRAFT)
    keel = min(vertex[2] for vertex in mesh.vertices)
    kb = 0.5 * DRAFT
    bm = (LENGTH * BEAM**3 / 12.0) / (LENGTH * BEAM * DRAFT)
    validity = wall_sided_validity(beam=BEAM, draft=DRAFT, depth=DEPTH)

    steps = tuple(BARGE_CURVE_STEP * index for index in range(int(180 / BARGE_CURVE_STEP) + 1))
    curve = righting_arms(mesh, kg=BARGE_KG, heel_angles=steps)

    table = []
    for heel in TABLE_ANGLES:
        arm = righting_arm(mesh, kg=BARGE_KG, heel_deg=heel)
        closed = box_reference_arm(heel, length=LENGTH, beam=BEAM, draft=DRAFT, kg=BARGE_KG)
        section = section_reference_arm(heel, beam=BEAM, draft=DRAFT, depth=DEPTH, kg=BARGE_KG)
        table.append(
            {
                "heel_deg": heel,
                "gz_m": arm.gz,
                "wall_sided_m": closed,
                "section_2d_m": section,
                "delta_wall_sided_m": arm.gz - closed,
                "delta_section_2d_m": arm.gz - section,
                "inside_wall_sided_regime": heel <= validity.binding_deg,
                "waterline_z_m": arm.waterline_z,
                "relative_volume_error": arm.relative_volume_error,
            }
        )

    # The two neutral loading identities, on a body where they can be checked exactly.
    mirrored = righting_arms(
        mesh, kg=BARGE_KG, heel_angles=(-60.0, -30.0, -10.0), displacement=curve.target_volume
    )
    symmetric = righting_arms(
        mesh, kg=BARGE_KG, heel_angles=(60.0, 30.0, 10.0), displacement=curve.target_volume
    )
    kg_identity = []
    identity_angles = (10.0, 30.0, 60.0)
    base = righting_arms(mesh, kg=BARGE_KG, heel_angles=identity_angles, displacement=curve.target_volume)
    raised = righting_arms(
        mesh, kg=BARGE_KG + 0.4, heel_angles=identity_angles, displacement=curve.target_volume
    )
    for base_arm, raised_arm in zip(base.arms, raised.arms):
        shift = 0.4 * math.sin(math.radians(base_arm.heel_deg))
        kg_identity.append(
            {
                "heel_deg": base_arm.heel_deg,
                "gz_base_m": base_arm.gz,
                "gz_raised_m": raised_arm.gz,
                "predicted_shift_m": -shift,
                "residual_m": raised_arm.gz - (base_arm.gz - shift),
            }
        )

    # An unstable body (GM < 0) comes to rest where the arm returns to zero.
    analytic_loll = wall_sided_vanishing_angle(kb=kb, bm=bm, kg=BARGE_UNSTABLE_KG)
    unstable = righting_arms(
        mesh, kg=BARGE_UNSTABLE_KG, heel_angles=tuple(0.5 * index for index in range(361))
    )
    crossings = zero_crossings(list(unstable.heel_deg), list(unstable.gz))
    measured_loll = next((angle for angle in crossings if angle > 0.0), None)

    # Dynamic stability: the same quadrature rule on the measured and the closed-form arm
    # separates the integration error from the geometry error.
    limit = validity.binding_deg
    fine = tuple(0.25 * index for index in range(int(limit / 0.25) + 1))
    area_curve = righting_arms(mesh, kg=BARGE_KG, heel_angles=fine)
    gm = kb + bm - BARGE_KG
    phi = math.radians(limit)
    area_closed = gm * (1.0 - math.cos(phi)) + 0.5 * bm * (
        math.cos(phi) + 1.0 / math.cos(phi) - 2.0
    )
    trapezoid_closed = 0.0
    for lower, upper in zip(fine, fine[1:]):
        trapezoid_closed += 0.25 * math.pi / 180.0 * 0.5 * (
            box_reference_arm(lower, length=LENGTH, beam=BEAM, draft=DRAFT, kg=BARGE_KG)
            + box_reference_arm(upper, length=LENGTH, beam=BEAM, draft=DRAFT, kg=BARGE_KG)
        )

    return {
        "geometry": {
            "length_m": LENGTH,
            "beam_m": BEAM,
            "draft_m": DRAFT,
            "depth_m": DEPTH,
            "freeboard_m": DEPTH - DRAFT,
            "keel_z_m": keel,
            "kg_m": BARGE_KG,
            "displaced_volume_m3": curve.target_volume,
            "kb_m": kb,
            "bm_m": bm,
            "gm_m": gm,
            "bilge_emergence_deg": validity.bilge_emergence_deg,
            "deck_immersion_deg": validity.deck_immersion_deg,
            "binding_deg": validity.binding_deg,
            "binding": validity.binding,
            "faces": len(mesh.faces),
        },
        "curve": curve.to_dict(),
        "table": table,
        "kg_identity": kg_identity,
        "mirror_residual_m": max(
            abs(mirror.gz + sym.gz) for mirror, sym in zip(mirrored.arms, symmetric.arms)
        ),
        "max_relative_volume_error": max(abs(arm.relative_volume_error) for arm in curve.arms),
        "max_delta_wall_sided_m": max(
            abs(row["delta_wall_sided_m"]) for row in table if row["inside_wall_sided_regime"]
        ),
        "max_delta_section_2d_m": max(abs(row["delta_section_2d_m"]) for row in table),
        "departure_at_edge_m": min(
            abs(row["delta_wall_sided_m"]) for row in table if row["heel_deg"] > validity.binding_deg
        ),
        "departure_past_regime_m": abs(
            righting_arm(mesh, kg=BARGE_KG, heel_deg=validity.binding_deg + 15.0).gz
            - box_reference_arm(
                validity.binding_deg + 15.0, length=LENGTH, beam=BEAM, draft=DRAFT, kg=BARGE_KG
            )
        ),
        "departure_past_regime_relative": abs(
            righting_arm(mesh, kg=BARGE_KG, heel_deg=validity.binding_deg + 15.0).gz
            - box_reference_arm(
                validity.binding_deg + 15.0, length=LENGTH, beam=BEAM, draft=DRAFT, kg=BARGE_KG
            )
        )
        / abs(box_reference_arm(validity.binding_deg + 15.0, length=LENGTH, beam=BEAM, draft=DRAFT, kg=BARGE_KG)),
        "unstable_kg_m": BARGE_UNSTABLE_KG,
        "unstable_gm_m": kb + bm - BARGE_UNSTABLE_KG,
        "analytic_loll_deg": analytic_loll,
        "measured_loll_deg": measured_loll,
        "loll_residual_deg": None
        if measured_loll is None or analytic_loll is None
        else measured_loll - analytic_loll,
        "area_to_limit_m_rad": area_curve.area(limit),
        "area_closed_form_m_rad": area_closed,
        "area_trapezoid_of_closed_form_m_rad": trapezoid_closed,
        "area_residual_m_rad": area_curve.area(limit) - area_closed,
        "area_quadrature_only_m_rad": trapezoid_closed - area_closed,
    }


def measure_hull() -> dict[str, object]:
    mesh = flared_hull(*LADDER[0])
    keel = min(vertex[2] for vertex in mesh.vertices)
    immersed = immersed_volume(mesh, HULL_WATERLINE)
    section = waterplane_moments(mesh, HULL_WATERLINE)
    kb = immersed.centre_of_buoyancy[2] - keel
    bm = section.transverse_second_moment / immersed.volume
    gm = kb + bm - HULL_KG

    curve = righting_arms(mesh, kg=HULL_KG, heel_angles=HULL_ANGLES, waterline_z=HULL_WATERLINE)
    raised = righting_arms(
        mesh, kg=HULL_KG_RAISED, heel_angles=HULL_ANGLES, waterline_z=HULL_WATERLINE
    )

    # The slope at zero heel, from the ratio GZ/sin(phi) at two small angles, with the
    # O(phi**2) term of the expansion removed by Richardson.
    slope = {}
    for heel in (2.0, 4.0):
        arm = righting_arm(mesh, kg=HULL_KG, heel_deg=heel, waterline_z=HULL_WATERLINE)
        slope[heel] = {
            "heel_deg": heel,
            "gz_m": arm.gz,
            "ratio_gz_over_sin_m": arm.gz / math.sin(math.radians(heel)),
        }
    ratio_two, ratio_four = slope[2.0]["ratio_gz_over_sin_m"], slope[4.0]["ratio_gz_over_sin_m"]
    richardson = (4.0**2 * ratio_two - 2.0**2 * ratio_four) / (4.0**2 - 2.0**2)

    ladder = []
    for nx, nz in LADDER:
        coarse = flared_hull(nx, nz) if (nx, nz) != LADDER[0] else mesh
        step_curve = righting_arms(
            coarse, kg=HULL_KG, heel_angles=(LADDER_HEEL,), waterline_z=HULL_WATERLINE
        )
        ladder.append(
            {
                "nx": nx,
                "nz": nz,
                "faces": len(coarse.faces),
                "gz_m": step_curve.arms[0].gz,
            }
        )
    coarse_arm, medium_arm, fine_arm = (row["gz_m"] for row in ladder)
    order = (
        math.log2(abs(coarse_arm - fine_arm) / abs(medium_arm - fine_arm))
        if medium_arm != fine_arm
        else float("inf")
    )

    kg_identity = [
        {
            "heel_deg": base.heel_deg,
            "residual_m": other.gz - (base.gz - (HULL_KG_RAISED - HULL_KG) * math.sin(base.heel_rad)),
        }
        for base, other in zip(curve.arms, raised.arms)
    ]

    # The shipped family: capped at its waterline, so no freeboard and no arm.
    shipped = generate_hull().mesh
    shipped_error = ""
    try:
        righting_arms(shipped, kg=0.3, heel_angles=(10.0,))
    except RightingArmError as error:
        shipped_error = str(error)

    return {
        "geometry": {
            "length_m": HULL_LENGTH,
            "beam_m": HULL_BEAM,
            "generated_draft_m": HULL_GENERATED_DRAFT,
            "design_waterline_z_m": HULL_WATERLINE,
            "immersed_draft_m": HULL_WATERLINE - keel,
            "freeboard_m": -HULL_WATERLINE,
            "kernel": "HullParameters(L, B, draft=2 t) floated at waterline_z = -t",
            "keel_z_m": keel,
            "kg_m": HULL_KG,
            "displaced_volume_m3": immersed.volume,
            "waterplane_area_m2": section.area,
            "transverse_second_moment_m4": section.transverse_second_moment,
            "kb_m": kb,
            "bm_m": bm,
            "gm_m": gm,
            "faces": len(mesh.faces),
        },
        "curve": curve.to_dict(),
        "raised_curve": raised.to_dict(),
        "slope": {
            "angles": list(slope.values()),
            "richardson_gm_m": richardson,
            "waterplane_gm_m": gm,
            "relative_residual": abs(richardson - gm) / gm,
        },
        "ladder": {
            "heel_deg": LADDER_HEEL,
            "rows": ladder,
            "observed_order": order,
            "coarse_relative_error": abs(coarse_arm - fine_arm) / abs(fine_arm),
        },
        "kg_identity_max_residual_m": max(abs(row["residual_m"]) for row in kg_identity),
        "max_relative_volume_error": max(abs(arm.relative_volume_error) for arm in curve.arms),
        "shipped_hull_error": shipped_error,
    }


# ---------------------------------------------------------------------------
# claims
# ---------------------------------------------------------------------------
def evaluate(results: dict[str, object]) -> list[dict[str, object]]:
    barge = results["barge"]
    hull = results["hull"]
    claims: list[dict[str, object]] = []

    def add(name: str, value: float, threshold: float, *, above: bool = False, detail: str = "") -> None:
        holds = value >= threshold if above else value <= threshold
        claims.append(
            {
                "claim": name,
                "value": value,
                "threshold": threshold,
                "holds": bool(holds),
                "detail": detail,
            }
        )

    add(
        "barge: the arm reproduces the wall-sided formula to round-off inside its regime",
        barge["max_delta_wall_sided_m"],
        1e-12,
    )
    add(
        "barge: the arm matches an independent 2D section computation at every heel "
        "(inside and outside the wall-sided regime)",
        barge["max_delta_section_2d_m"],
        1e-12,
    )
    add(
        "barge: the displacement is held at every heel angle",
        barge["max_relative_volume_error"],
        1e-12,
    )
    add("barge: GZ(0) = 0 and the curve is odd in the heel", barge["mirror_residual_m"], 1e-13)
    add(
        "barge: raising KG by d lowers the arm by exactly d sin(phi)",
        max(abs(row["residual_m"]) for row in barge["kg_identity"]),
        1e-13,
    )
    add(
        "barge with GM < 0: the resting heel matches tan(phi) = sqrt(-2 GM / BM)",
        abs(barge["loll_residual_deg"]),
        0.05,
        detail=f"analytic {barge['analytic_loll_deg']:.4f} deg",
    )
    add(
        "barge: the area under the arm matches the closed-form integral",
        abs(barge["area_residual_m_rad"]),
        5e-6,
        detail=f"quadrature alone contributes {barge['area_quadrature_only_m_rad']:.2e}",
    )
    add(
        "barge: 15 deg past its regime edge the closed form is off by more than 10 % of the arm "
        "(so the limit is not cosmetic)",
        barge["departure_past_regime_relative"],
        0.1,
        above=True,
        detail=f"absolute departure {barge['departure_past_regime_m']:.3e} m; "
        f"only {barge['departure_at_edge_m']:.1e} m at the first degree past the edge",
    )
    add(
        "hull: the slope of the arm at zero heel reproduces KB + I_T/V - KG from the "
        "waterplane integrator",
        hull["slope"]["relative_residual"],
        1e-4,
    )
    add(
        "hull: the arm is second order in the grid step (fitted order > 1.8)",
        hull["ladder"]["observed_order"],
        1.8,
        above=True,
        detail=f"coarse grid is {100 * hull['ladder']['coarse_relative_error']:.2f} % off the fine one",
    )
    add(
        "hull: the shipped default grid is within 5 % of the finest at 20 deg of heel",
        hull["ladder"]["coarse_relative_error"],
        0.05,
    )
    add(
        "hull: raising KG by d lowers the arm by exactly d sin(phi)",
        hull["kg_identity_max_residual_m"],
        1e-12,
    )
    add(
        "hull: the displacement is held at every heel angle",
        hull["max_relative_volume_error"],
        1e-12,
    )
    claims.append(
        {
            "claim": "shipped hull (capped at the waterline, no freeboard): the arm is refused, "
            "not returned as zeros",
            "value": 1.0 if "freeboard" in str(hull["shipped_hull_error"]) else 0.0,
            "threshold": 1.0,
            "holds": "freeboard" in str(hull["shipped_hull_error"]),
            "detail": str(hull["shipped_hull_error"])[:160],
        }
    )
    return claims


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def figure_curves(results: dict[str, object], output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    barge, hull = results["barge"], results["hull"]
    figure, (left, right) = plt.subplots(1, 2, figsize=(12.0, 4.6))

    table = barge["table"]
    heels = [row["heel_deg"] for row in table]
    left.plot(heels, [row["gz_m"] for row in table], "o-", label="mesh integration")
    left.plot(
        [row["heel_deg"] for row in table if row["inside_wall_sided_regime"]],
        [row["wall_sided_m"] for row in table if row["inside_wall_sided_regime"]],
        "--",
        label="wall-sided closed form (exact in its regime)",
    )
    left.plot(
        [row["heel_deg"] for row in table if not row["inside_wall_sided_regime"]],
        [row["wall_sided_m"] for row in table if not row["inside_wall_sided_regime"]],
        ":",
        color="tab:red",
        label="same formula past its regime (no claim)",
    )
    left.plot(
        heels, [row["section_2d_m"] for row in table], "x", color="black", label="independent 2D section"
    )
    left.axvline(
        barge["geometry"]["binding_deg"],
        color="grey",
        lw=0.8,
        label=f"regime edge ({barge['geometry']['binding_deg']:.1f} deg, {barge['geometry']['binding']})",
    )
    left.set_xlabel("heel [deg]")
    left.set_ylabel("GZ [m]")
    left.set_title(f"barge {LENGTH:g}x{BEAM:g}x{DEPTH:g} m, T={DRAFT:g} m, KG={BARGE_KG:g} m")
    left.grid(alpha=0.3)
    left.legend(fontsize=7, loc="upper left")

    arms = hull["curve"]["arms"]
    raised = hull["raised_curve"]["arms"]
    right.plot([a["heel_deg"] for a in arms], [a["gz_m"] for a in arms], "o-", label=f"KG = {HULL_KG:g} m")
    right.plot(
        [a["heel_deg"] for a in raised], [a["gz_m"] for a in raised], "s-", label=f"KG = {HULL_KG_RAISED:g} m"
    )
    gm = hull["geometry"]["gm_m"]
    tangent = [gm * math.sin(math.radians(a["heel_deg"])) for a in arms[:5]]
    right.plot([a["heel_deg"] for a in arms[:5]], tangent, "k:", label="GM sin(phi), GM from waterplane")
    right.set_xlabel("heel [deg]")
    right.set_ylabel("GZ [m]")
    right.set_title(
        f"hull with topsides (41x13, {hull['geometry']['immersed_draft_m']:.2f} m immersed, "
        f"{hull['geometry']['freeboard_m']:.2f} m freeboard)"
    )
    right.grid(alpha=0.3)
    right.legend(fontsize=8)

    figure.tight_layout()
    figure.savefig(output / "gz_curves.png", dpi=140)
    plt.close(figure)


def figure_deviation(results: dict[str, object], output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    barge = results["barge"]
    table = barge["table"]
    figure, (left, right) = plt.subplots(1, 2, figsize=(12.0, 4.4))

    left.semilogy(
        [row["heel_deg"] for row in table],
        [max(abs(row["delta_wall_sided_m"]), 1e-18) for row in table],
        "o-",
        label="|module - wall-sided|",
    )
    left.semilogy(
        [row["heel_deg"] for row in table],
        [max(abs(row["delta_section_2d_m"]), 1e-18) for row in table],
        "x-",
        label="|module - independent 2D section|",
    )
    left.axvline(barge["geometry"]["binding_deg"], color="grey", lw=0.8)
    left.set_xlabel("heel [deg]")
    left.set_ylabel("absolute difference [m]")
    left.set_title("a closed form has a regime; an integrator has an error")
    left.grid(alpha=0.3, which="both")
    left.legend(fontsize=8)

    ladder = results["hull"]["ladder"]
    steps = [1.0 / (row["nx"] - 1) for row in ladder["rows"]]
    errors = [abs(row["gz_m"] - ladder["rows"][-1]["gz_m"]) for row in ladder["rows"]]
    right.loglog(steps, [max(value, 1e-12) for value in errors], "o-", label="|GZ(coarse) - GZ(fine)|")
    reference = [errors[0] * (step / steps[0]) ** 2 for step in steps]
    right.loglog(steps, reference, "k--", label="second order reference")
    right.set_xlabel("grid step 1/(nx-1) [1/m]")
    right.set_ylabel("arm error at 20 deg [m]")
    right.set_title(f"hull arm convergence, fitted order {ladder['observed_order']:.2f}")
    right.grid(alpha=0.3, which="both")
    right.legend(fontsize=8)

    figure.tight_layout()
    figure.savefig(output / "deviation.png", dpi=140)
    plt.close(figure)


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def render_report(results: dict[str, object], claims: list[dict[str, object]]) -> str:
    barge, hull = results["barge"], results["hull"]
    geometry = barge["geometry"]
    lines = [
        "# Righting arms: measured tables",
        "",
        f"`tensorcad {results['version']}`, generated by "
        "`python examples/righting_arm/run_example.py`.",
        "",
        "## Claims",
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
        "## Barge: the arm against a closed form and an independent section solve",
        "",
        f"`L={geometry['length_m']:g}`, `B={geometry['beam_m']:g}`, `T={geometry['draft_m']:g}`, "
        f"`D={geometry['depth_m']:g}` m, `KG={geometry['kg_m']:g}` m, "
        f"`V={geometry['displaced_volume_m3']:.9f} m3`, `KB={geometry['kb_m']:g}`, "
        f"`BM={geometry['bm_m']:.9f}`, `GM={geometry['gm_m']:.9f}` m.",
        "",
        f"The closed form is exact only while the waterline spans the full beam; for this body "
        f"that ends at {geometry['binding_deg']:.4f} deg ({geometry['binding']}).",
        "",
        "| heel [deg] | GZ mesh [m] | wall-sided [m] | 2D section [m] | mesh - wall [m] | mesh - 2D [m] | in regime | waterline z [m] |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in barge["table"]:
        lines.append(
            f"| {row['heel_deg']:g} | {row['gz_m']:.12f} | {row['wall_sided_m']:.12f} | "
            f"{row['section_2d_m']:.12f} | {row['delta_wall_sided_m']:.2e} | "
            f"{row['delta_section_2d_m']:.2e} | "
            f"{'yes' if row['inside_wall_sided_regime'] else 'no'} | {row['waterline_z_m']:.2e} |"
        )

    lines += [
        "",
        "## Barge: loading identities",
        "",
        "| heel [deg] | GZ at KG | GZ at KG + 0.4 m | predicted shift [m] | residual [m] |",
        "|---|---|---|---|---|",
    ]
    for row in barge["kg_identity"]:
        lines.append(
            f"| {row['heel_deg']:g} | {row['gz_base_m']:.12f} | {row['gz_raised_m']:.12f} | "
            f"{row['predicted_shift_m']:.12f} | {row['residual_m']:.2e} |"
        )

    lines += [
        "",
        "## Barge: unstable body, dynamic stability",
        "",
        f"* `GM = {barge['unstable_gm_m']:.9f} m < 0` at `KG = {barge['unstable_kg_m']:g} m`: the arm "
        f"crosses zero at {barge['measured_loll_deg']:.6f} deg, the closed form predicts "
        f"{barge['analytic_loll_deg']:.6f} deg (residual {barge['loll_residual_deg']:+.2e} deg, "
        "grid-limited: the crossing is interpolated between 0.5 deg samples).",
        f"* Area under the arm to {geometry['binding_deg']:.2f} deg: "
        f"{barge['area_to_limit_m_rad']:.12f} m rad measured, "
        f"{barge['area_closed_form_m_rad']:.12f} m rad from the closed-form integral, "
        f"residual {barge['area_residual_m_rad']:.2e} (of which "
        f"{barge['area_quadrature_only_m_rad']:.2e} is the trapezoid on a 0.25 deg grid).",
        "",
        "## Hull with topsides",
        "",
        f"Generated with `draft={hull['geometry']['generated_draft_m']:g}` m and floated at "
        f"`z={hull['geometry']['design_waterline_z_m']:g}` m, so "
        f"{hull['geometry']['immersed_draft_m']:.3f} m immersed and "
        f"{hull['geometry']['freeboard_m']:.3f} m of freeboard at the centreline; "
        f"`V={hull['geometry']['displaced_volume_m3']:.9f} m3`, "
        f"`KB={hull['geometry']['kb_m']:.9f}`, `BM={hull['geometry']['bm_m']:.9f}`, "
        f"`GM={hull['geometry']['gm_m']:.9f}` m at `KG={hull['geometry']['kg_m']:g}` m.",
        "",
        "| heel [deg] | GZ [m] | displaced volume [m3] | relative volume error |",
        "|---|---|---|---|",
    ]
    for arm in hull["curve"]["arms"]:
        lines.append(
            f"| {arm['heel_deg']:g} | {arm['gz_m']:.12f} | {arm['displaced_volume_m3']:.9f} | "
            f"{arm['relative_volume_error']:.1e} |"
        )

    slope = hull["slope"]
    lines += [
        "",
        "### Slope at zero heel against the waterplane metacentric height",
        "",
        "| heel [deg] | GZ [m] | GZ / sin(phi) [m] |",
        "|---|---|---|",
    ]
    for row in slope["angles"]:
        lines.append(
            f"| {row['heel_deg']:g} | {row['gz_m']:.12f} | {row['ratio_gz_over_sin_m']:.12f} |"
        )
    lines += [
        "",
        f"Richardson limit of the ratio: `{slope['richardson_gm_m']:.12f}` m against "
        f"`KB + I_T/V - KG = {slope['waterplane_gm_m']:.12f}` m "
        f"(relative residual {slope['relative_residual']:.2e}).",
        "",
        "### Convergence of the arm at 20 deg of heel",
        "",
        "| nx | nz | faces | GZ [m] |",
        "|---|---|---|---|",
    ]
    for row in hull["ladder"]["rows"]:
        lines.append(f"| {row['nx']} | {row['nz']} | {row['faces']} | {row['gz_m']:.12f} |")
    lines += [
        "",
        f"Observed order {hull['ladder']['observed_order']:.3f}; the shipped default grid "
        f"(41x13) is {100 * hull['ladder']['coarse_relative_error']:.2f} % above the 161x41 value.",
        "",
        "## The shipped hull family has no large-angle curve",
        "",
        "`generate_hull` caps the surface at the design waterline, so the displacement is the "
        "whole body volume and no heeled waterline can hold it. The module raises:",
        "",
        "```",
        str(hull["shipped_hull_error"]),
        "```",
        "",
        "## Figures",
        "",
        "- `gz_curves.png` -- the barge's arm against the closed form and the 2D section solve, "
        "and the hull's arm at two loadings with the `GM sin(phi)` tangent.",
        "- `deviation.png` -- the two differences at 1e-18 resolution, and the arm's second-order "
        "convergence under refinement.",
        "",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("docs/righting_arm"))
    arguments = parser.parse_args(argv)
    output = arguments.output
    output.mkdir(parents=True, exist_ok=True)

    results: dict[str, object] = {
        "version": __version__,
        "barge": measure_barge(),
        "hull": measure_hull(),
    }
    claims = evaluate(results)
    results["claims"] = claims

    figures = []
    try:
        figure_curves(results, output)
        figure_deviation(results, output)
        figures = ["gz_curves.png", "deviation.png"]
    except ImportError:
        print("[righting_arm] matplotlib not available: figures skipped", file=sys.stderr)

    (output / "results.json").write_text(
        json.dumps(results, indent=2, sort_keys=True, default=str) + "\n"
    )
    (output / "report.md").write_text(render_report(results, claims))

    failed = [claim for claim in claims if not claim["holds"]]
    for claim in claims:
        status = "PASS" if claim["holds"] else "FAIL"
        print(
            f"[righting_arm] {status}  {claim['claim']}  "
            f"({claim['value']:.3e} vs {claim['threshold']:.0e})"
        )
    print(f"[righting_arm] wrote {output}/results.json, report.md, {', '.join(figures) or 'no figures'}")
    if failed:
        print(f"[righting_arm] {len(failed)} claim(s) do not hold", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())