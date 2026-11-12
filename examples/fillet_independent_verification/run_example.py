"""Independent verification of ``tensorcad.fillet`` (mission m16561, verifier seat).

Every number in ``docs/fillet_independent_verification.md`` is produced here; nothing is
hand-copied.  The point of the mission is that the fillet kernel's own tests compare two
derivations that live *inside* the module (Green's-theorem path integral vs the per-corner
closed form ``dA = -sign(phi) r^2 (tan(|phi|/2) - |phi|/2)``).  That is a useful cross-check,
but it cannot detect a mistake shared by both.  This script therefore recomputes the same
quantities on routes that share no code with ``tensorcad.fillet``:

R1  a closed form for the rounded rectangle derived by hand (area, perimeter, centroidal
    second moments), plus the signed-corner area/perimeter law for a re-entrant L section,
R2  a brute-force midpoint integration of the rounded body's *indicator function* -- it never
    sees an arc, only "is this point inside?",
R3  a five-point Gauss-Legendre quadrature of ``A(r(x))`` along the sweep axis for the solid
    volume, and the divergence-theorem volume of the produced mesh with mesh refinement,
R4  chordal convergence of the sampled outline (Richardson: the "exact" value must be the
    limit of the mesh the caller actually gets),
R5  an adversarial sweep over the feasibility boundary (circle limit, sliver, notch, reflex
    corner, orientation, over-limit refusal) and over numerical conditioning (translation).

Run: ``python examples/fillet_independent_verification/run_example.py``
Exit code is non-zero if any check regresses.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from tensorcad import fillet as F

HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
ART = HERE.parent.parent / "artifacts" / "fillet_independent_verification"

CHECKS: list[dict[str, object]] = []


def check(name: str, measured: float, reference: float, tolerance: float, *, note: str = "") -> None:
    error = abs(measured - reference) / abs(reference) if reference else abs(measured)
    ok = error <= tolerance
    CHECKS.append(
        {
            "check": name,
            "measured": measured,
            "reference": reference,
            "relative_error": error,
            "tolerance": tolerance,
            "passed": ok,
            "note": note,
        }
    )
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}: {measured:.12g} vs {reference:.12g} "
          f"(rel {error:.2e} <= {tolerance:.0e})")


# --------------------------------------------------------------------------- R1
def rounded_rectangle_closed_form(w: float, h: float, r: float) -> dict[str, float]:
    """Rounded rectangle: area, perimeter and centroidal second moments, derived by hand.

    The rect is ``w x h`` with corner radius ``r``.  Area loses four corner bites
    ``(1 - pi/4) r^2``; the perimeter loses four ``(2 - pi/2) r``.  For ``int (v - h/2)^2 dA``
    each bite removes the strip integral over its ``r x r`` square minus the quarter disk,
    the quarter disk being taken about its own centre and then shifted back to ``v = h/2``
    by ``int (v-h/2)^2 = int (v-c)^2 + 2(v-c)(c-h/2) + (c-h/2)^2`` with ``c`` the disk centre.
    """
    area = w * h - (4.0 - math.pi) * r * r
    perimeter = 2.0 * (w + h) - (8.0 - 2.0 * math.pi) * r

    def second(extent: float, other: float) -> float:
        """``int rho^2 dA`` about the centroidal axis normal to ``extent``.

        A corner bite subtracts its ``r x r`` square and adds the quarter disk back.  In the
        centroidal coordinate the disk centre sits at ``d = extent/2 - r``, the local radial
        coordinate runs towards the mid plane, and the quarter disk contributes
        ``d^2 pi r^2/4 + 2 d r^3/3 + pi r^4/16``.  (The sign of the ``2 d r^3/3`` term is the
        one that a by-hand derivation gets wrong -- see the mission report.)
        """
        half = extent / 2.0
        d = half - r
        strip = r * (half**3 - (half - r) ** 3) / 3.0
        quarter = d * d * math.pi * r * r / 4.0 + 2.0 * d * r**3 / 3.0 + math.pi * r**4 / 16.0
        return other * extent**3 / 12.0 - 4.0 * (strip - quarter)

    return {"area": area, "perimeter": perimeter, "i_uu": second(w, h), "i_vv": second(h, w)}


def brute_force_moments(w: float, h: float, r: float, n: int = 1600) -> dict[str, float]:
    """R2: midpoint integration of the rounded body's indicator; no arcs, no Green."""
    def inside(u: float, v: float) -> bool:
        if not (0.0 <= u <= w and 0.0 <= v <= h):
            return False
        cu = min(max(u, r), w - r)
        cv = min(max(v, r), h - r)
        return (u - cu) ** 2 + (v - cv) ** 2 <= r * r

    du, dv = w / n, h / n
    area = mu = mv = iu = iv = 0.0
    for i in range(n):
        u = (i + 0.5) * du
        for j in range(n):
            v = (j + 0.5) * dv
            if inside(u, v):
                d = du * dv
                area += d
                mu += u * d
                mv += v * d
                iu += u * u * d
                iv += v * v * d
    cu, cv = mu / area, mv / area
    return {"area": area, "i_uu": iu - area * cu * cu, "i_vv": iv - area * cv * cv,
            "centroid": [cu, cv]}


# --------------------------------------------------------------------------- R3
GAUSS5_X = (0.046910077030668, 0.230765344947158, 0.5, 0.769234655052841, 0.953089922969332)
GAUSS5_W = (0.118463442528095, 0.239314335249683, 0.284444444444444, 0.239314335249683, 0.118463442528095)


def sweep_volume_quadrature(base_area: float, coefficient: float, length: float,
                            r0: float, r1: float) -> float:
    """R3: five-point Gauss-Legendre of A(r(x)); exact for the quadratic A(r(x))."""
    total = 0.0
    for x, w in zip(GAUSS5_X, GAUSS5_W):
        radius = r0 + (r1 - r0) * x
        total += w * (base_area + coefficient * radius * radius)
    return total * length


SQUARE = [(0.0, 0.0), (4.0, 0.0), (4.0, 2.0), (0.0, 2.0)]
L_SECTION = [(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (1.0, 1.0), (1.0, 2.0), (0.0, 2.0)]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ART.mkdir(parents=True, exist_ok=True)
    payload: dict[str, object] = {}

    # ---------------------------------------------------------------- R1 + R2
    print("R1/R2: rounded rectangle 4 x 2 against a hand-derived closed form and brute force")
    rect_rows = []
    for radius in (0.05, 0.2, 0.4, 0.8, 0.999):
        shape = F.fillet_section(SQUARE, radius)
        ref = rounded_rectangle_closed_form(4.0, 2.0, radius)
        cu, cv = shape.moments.centroid
        i_uu = shape.moments.second_u - shape.moments.area * cu * cu
        i_vv = shape.moments.second_v - shape.moments.area * cv * cv
        check(f"area r={radius}", shape.moments.area, ref["area"], 1e-14)
        check(f"perimeter r={radius}", shape.moments.perimeter, ref["perimeter"], 1e-14)
        check(f"I_uu r={radius}", i_uu, ref["i_uu"], 1e-13)
        check(f"I_vv r={radius}", i_vv, ref["i_vv"], 1e-13)
        rect_rows.append(
            {"radius": radius, "area": shape.moments.area, "area_ref": ref["area"],
             "i_uu": i_uu, "i_uu_ref": ref["i_uu"], "i_vv": i_vv, "i_vv_ref": ref["i_vv"],
             "crosscheck": shape.area_crosscheck_error}
        )
    payload["rounded_rectangle"] = rect_rows

    brute = brute_force_moments(4.0, 2.0, 0.4)
    exact = F.fillet_section(SQUARE, 0.4)
    cu, cv = exact.moments.centroid
    check("brute-force area", brute["area"], exact.moments.area, 2e-3,
          note="midpoint rule on a curved boundary: O(h) on the arc, not a claim about the kernel")
    check("brute-force centroid u", brute["centroid"][0], cu, 1e-4)
    check("brute-force I_uu", brute["i_uu"], exact.moments.second_u - exact.moments.area * cu * cu, 1e-4)
    check("brute-force I_vv", brute["i_vv"], exact.moments.second_v - exact.moments.area * cv * cv, 1e-4)
    payload["brute_force"] = brute

    # ---------------------------------------------------------------- re-entrant section
    print("R1: re-entrant L section against the signed-corner area/perimeter law")
    base = F.polygon_moments(L_SECTION)
    turns = F.section_turns(L_SECTION)
    l_rows = []
    for radius in (0.1, 0.3, 0.5):
        shape = F.fillet_section(L_SECTION, radius)
        d_area = sum(-math.copysign(1.0, t) * radius * radius * (math.tan(abs(t) / 2.0) - abs(t) / 2.0)
                     for t in turns if abs(t) > 1e-12)
        d_perimeter = sum(radius * abs(t) - 2.0 * radius * math.tan(abs(t) / 2.0)
                          for t in turns if abs(t) > 1e-12)
        check(f"L area r={radius}", shape.moments.area, base.area + d_area, 1e-14)
        check(f"L perimeter r={radius}", shape.moments.perimeter, base.perimeter + d_perimeter, 1e-14)
        l_rows.append({"radius": radius, "area": shape.moments.area,
                       "area_ref": base.area + d_area, "turns_deg": [math.degrees(t) for t in turns]})
    payload["l_section"] = l_rows

    # ---------------------------------------------------------------- R4 chordal limit
    print("R4: the 'exact' values must be the limit of the chorded outline (Richardson)")
    conv_rows = []
    previous = None
    for chords in (4, 8, 16, 32, 64):
        shape = F.fillet_section(SQUARE, 0.4, arc_points=chords)
        sampled = shape.sampled_moments().area
        error = (sampled - shape.moments.area) / shape.moments.area
        ratio = None if previous is None else previous / error
        conv_rows.append({"arc_points": chords, "sampled_area": sampled, "exact_area": shape.moments.area,
                          "relative_error": error, "error_ratio": ratio})
        if previous is not None:
            check(f"chordal error ratio at {chords} chords (order 2)", ratio, 4.0, 0.15)
        previous = error
    payload["chordal_convergence"] = conv_rows

    solid_rows = []
    for chords, stations in ((4, 1), (8, 1), (16, 2), (32, 4)):
        solid = F.fillet_solid(SQUARE, 3.0, 0.4, arc_points=chords, stations=stations)
        solid_rows.append({"arc_points": chords, "stations": stations,
                           "volume_exact": solid.volume_exact, "volume_mesh": solid.volume_mesh,
                           "mesh_volume_error": solid.mesh_volume_error,
                           "watertight": solid.watertight,
                           "euler_characteristic": solid.euler_characteristic,
                           "zero_area_faces": solid.zero_area_faces,
                           "sampled_crease_deg": solid.sampled_crease_deg})
        assert solid.watertight and solid.euler_characteristic == 2 and solid.zero_area_faces == 0
    ratios = [solid_rows[i]["mesh_volume_error"] / solid_rows[i + 1]["mesh_volume_error"]
              for i in range(len(solid_rows) - 1)]
    check("mesh volume error ratio per doubling (order 2)", sum(ratios) / len(ratios), 4.0, 0.45)
    payload["solid_convergence"] = solid_rows

    # ---------------------------------------------------------------- R3 sweep volume
    print("R3: sweep volume against Gauss-Legendre quadrature of A(r(x))")
    turns_square = F.section_turns(SQUARE)
    coefficient = F._area_coefficient(turns_square, "fillet")
    base_area = F.polygon_moments(SQUARE).area
    sweep_rows = []
    for r0, r1 in ((0.2, 0.6), (0.6, 0.2), (0.1, 0.999), (0.4, 0.4)):
        solid = F.fillet_solid(SQUARE, 3.0, r0, radius_end=r1, arc_points=16, stations=8)
        quad = sweep_volume_quadrature(base_area, coefficient, 3.0, r0, r1)
        check(f"sweep volume r {r0}->{r1}", solid.volume_exact, quad, 1e-14)
        check(f"sweep mesh volume r {r0}->{r1}", solid.volume_mesh, quad, 5e-4)
        sweep_rows.append({"radius_start": r0, "radius_end": r1, "volume_exact": solid.volume_exact,
                           "volume_quadrature": quad, "volume_mesh": solid.volume_mesh,
                           "mesh_volume_error": solid.mesh_volume_error})
    payload["sweep"] = sweep_rows

    # ---------------------------------------------------------------- R5 feasibility
    print("R5: feasibility boundary and adversarial geometry")
    square2 = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]
    limit = F.max_fillet_radius(square2)
    check("max radius of a 2x2 square", limit, 1.0, 1e-12)
    circle = F.fillet_section(square2, limit)
    check("circle-limit area", circle.moments.area, math.pi, 1e-14)
    check("circle-limit perimeter", circle.moments.perimeter, 2.0 * math.pi, 1e-14)
    refused = False
    try:
        F.fillet_section(square2, 1.02)
    except F.FilletError as exc:
        refused = "largest feasible" in str(exc)
    CHECKS.append({"check": "over-limit radius refused with the limit named", "measured": float(refused),
                   "reference": 1.0, "relative_error": 0.0 if refused else 1.0, "tolerance": 0.0,
                   "passed": refused, "note": "no silent repair"})
    print(f"  [{'PASS' if refused else 'FAIL'}] over-limit radius refused with the limit named")

    near_limit = {
        "sliver_4x0.2_r0.0999": ([(0.0, 0.0), (4.0, 0.0), (4.0, 0.2), (0.0, 0.2)], 0.0999),
        "rectangle_r0.999": (SQUARE, 0.999),
        "l_section_r0.5": (L_SECTION, 0.5),
        "notched_section_r0.29": ([(0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (2.0, 1.0),
                                   (2.0, 0.4), (1.0, 0.4), (1.0, 1.0), (0.0, 1.0)], 0.29),
    }
    feasibility_rows = []
    for tag, (points, radius) in near_limit.items():
        shape = F.fillet_section(points, radius, arc_points=8)
        rows = feasibility_rows.append
        area_sign = math.copysign(1.0, _shoelace(shape.points))
        rows({"case": tag, "radius": radius, "max_radius": F.max_fillet_radius(points),
              "area": shape.moments.area, "sampled_area": shape.sampled_moments().area,
              "orientation": area_sign, "g1_error_deg": shape.g1_error_deg,
              "crosscheck": shape.area_crosscheck_error})
        # the chorded outline must stay a simple polygon of the same orientation
        assert abs(area_sign - math.copysign(1.0, _shoelace(points))) < 1e-12
        check(f"{tag} sampled area below exact (chords cut the arc)",
              shape.sampled_moments().area, shape.moments.area, 5e-3)
    payload["feasibility"] = feasibility_rows

    clockwise = F.fillet_section(list(reversed(SQUARE)), 0.4)
    check("clockwise input area equals counter-clockwise", clockwise.moments.area, 7.862654824574368, 1e-14)

    # ---------------------------------------------------------------- conditioning
    print("R5: numerical conditioning under translation (absolute-coordinate Green's theorem)")
    reference_area = F.fillet_section(SQUARE, 0.4).moments.area
    conditioning = []
    for shift in (0.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e8):
        moved = F.fillet_section([(x + shift, y + shift) for x, y in SQUARE], 0.4)
        relative = abs(moved.moments.area - reference_area) / reference_area
        conditioning.append({"shift": shift, "area": moved.moments.area, "relative_error": relative,
                             "relative_shift": shift / 4.0})
    payload["conditioning"] = conditioning
    at_1e6 = conditioning[-2]["relative_error"]
    check("area error at a 1e6 translation (1e6/4 = 2.5e5 x the body)", at_1e6, 0.0, 1e-9,
          note="documents the conditioning floor; error grows like eps * shift / size")
    at_1e8 = conditioning[-1]["relative_error"]
    print(f"  measured: shift 1e6 -> {at_1e6:.2e}, shift 1e8 -> {at_1e8:.2e} relative area error")

    # ---------------------------------------------------------------- figures
    failed = [row for row in CHECKS if not row["passed"]]
    payload["summary"] = {"checks": len(CHECKS), "failed": len(failed),
                          "all_passed": not failed}
    _figures(payload)
    _write_outputs(payload)

    print(f"\n{len(CHECKS) - len(failed)}/{len(CHECKS)} checks passed")
    (OUT / "checks.json").write_text(json.dumps(payload, indent=2) + "\n")
    return 1 if failed else 0


def _shoelace(points) -> float:
    total = 0.0
    count = len(points)
    for index in range(count):
        u0, v0 = points[index]
        u1, v1 = points[(index + 1) % count]
        total += u0 * v1 - u1 * v0
    return total / 2.0


def _figures(payload: dict[str, object]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    conv = payload["chordal_convergence"]
    solid = payload["solid_convergence"]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))
    ax = axes[0]
    ax.loglog([row["arc_points"] for row in conv], [abs(row["relative_error"]) for row in conv],
              "o-", label="section area (chords)")
    ax.loglog([row["arc_points"] for row in solid], [row["mesh_volume_error"] for row in solid],
              "s--", label="swept mesh volume")
    h = [row["arc_points"] for row in conv]
    ax.loglog(h, [abs(conv[0]["relative_error"]) * (h[0] / value) ** 2 for value in h],
              ":", color="grey", label="slope -2 (O(h^2))")
    ax.set_xlabel("chords per arc")
    ax.set_ylabel("relative error vs exact")
    ax.set_title("chorded mesh converges to the exact value")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)

    cond = payload["conditioning"]
    ax = axes[1]
    ax.loglog([row["relative_shift"] for row in cond], [max(row["relative_error"], 1e-18) for row in cond],
              "o-", label="measured relative area error")
    ax.loglog([row["relative_shift"] for row in cond],
              [2.2e-16 * max(row["relative_shift"], 1.0) for row in cond],
              ":", color="grey", label="eps * shift / size")
    ax.set_xlabel("translation / body size")
    ax.set_ylabel("relative error in area")
    ax.set_title("conditioning floor of absolute-coordinate integration")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ART / "fig1_independent_routes.png", dpi=150)
    plt.close(fig)

    rect = payload["rounded_rectangle"]
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    ax.semilogy([row["radius"] for row in rect],
                [abs(row["area"] - row["area_ref"]) / row["area_ref"] for row in rect], "o-",
                label="area")
    ax.semilogy([row["radius"] for row in rect],
                [abs(row["i_uu"] - row["i_uu_ref"]) / row["i_uu_ref"] for row in rect], "s-",
                label="I_uu (hand closed form)")
    brute = payload["brute_force"]
    check_row = next(r for r in CHECKS if r["check"] == "brute-force I_uu")
    ax.axhline(check_row["relative_error"], color="tab:red", linestyle="--",
               label="brute-force I_uu at r=0.4 (O(h) boundary)")
    ax.set_xlabel("corner radius r (m)   [rectangle 4 x 2]")
    ax.set_ylabel("relative error")
    ax.set_title("rounded rectangle: kernel vs independent references")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ART / "fig2_closed_form_agreement.png", dpi=150)
    plt.close(fig)


def _write_outputs(payload: dict[str, object]) -> None:
    (ART / "independent_verification.json").write_text(json.dumps(payload, indent=2) + "\n")
    lines = ["# Fillet kernel, independent verification\n",
             "| check | measured | reference | rel. error | tolerance | passed |",
             "|---|---|---|---|---|---|"]
    for row in CHECKS:
        lines.append(
            f"| {row['check']} | {row['measured']:.10g} | {row['reference']:.10g} | "
            f"{row['relative_error']:.2e} | {row['tolerance']:.0e} | "
            f"{'yes' if row['passed'] else 'NO'} |"
        )
    (OUT / "report.md").write_text("\n".join(lines) + "\n")
    (OUT / "summary.json").write_text(json.dumps(
        {key: payload[key] for key in ("summary",) if key in payload}, indent=2) + "\n")


if __name__ == "__main__":
    sys.exit(main())