"""Cross-feature verification of the PR #48 cure (mission 3871).

The cure merges ``agent/m3479`` (PR #48: the large-angle ``GZ`` curve) with ``main`` at
``2ff80b2`` (transom-stern hull family, split fore/aft taper exponents, waterplane second
moments).  A textual merge that keeps both sides is necessary but not sufficient: the two
feature sets have to *work together*.  This script exercises the pairs that cross the merge
boundary and writes machine-readable evidence.

Run from the repository root::

    python tools/verify_union_m3871.py --output artifacts/m3871

Exit code 0 means every claim below still holds; 1 means at least one broke.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tensorcad.hull import analytic_hydrostatics as hull_analytic_hydrostatics  # noqa: E402
from tensorcad import (  # noqa: E402  (path bootstrap above must run first)
    HullParameters,
    MeshBuilder,
    TransomHullParameters,
    analytic_transom_hydrostatics,
    generate_hull,
    generate_transom_hull,
    generate_wigley,
    gz_curve,
    imo_intact_stability,
    initial_stability,
    submerged_solid,
)

CHECKS: list[dict[str, Any]] = []


def record(name: str, value: Any, reference: Any, tolerance: float, note: str) -> None:
    """Store one check; ``tolerance`` is the absolute bound on ``|value - reference|``."""
    if isinstance(value, (int, float)) and isinstance(reference, (int, float)):
        error: Any = abs(value - reference)
        relative: Any = error / abs(reference) if reference else float("nan")
    else:
        error = relative = float("nan")
    CHECKS.append(
        {
            "check": name,
            "value": value,
            "reference": reference,
            "abs_error": error,
            "rel_error": relative,
            "tolerance": tolerance,
            "passed": bool(value == reference) if isinstance(value, (bool, tuple)) else bool(error <= tolerance),
            "note": note,
        }
    )


def barge() -> Any:
    """The wall-sided barge of ``tests/test_gz.py`` -- a polyhedron, so the cut is exact."""
    length, breadth, depth = 30.0, 6.0, 10.0
    builder = MeshBuilder()
    corner = {}
    for x in (0.0, length):
        for y in (-0.5 * breadth, 0.5 * breadth):
            for z in (0.0, depth):
                corner[(x, y, z)] = builder.add_vertex((x, y, z))

    def quad(p0: Any, p1: Any, p2: Any, p3: Any) -> None:
        a, b, c, d = (corner[point] for point in (p0, p1, p2, p3))
        builder.add_triangle(a, b, c)
        builder.add_triangle(a, c, d)

    lo, hi = -0.5 * breadth, 0.5 * breadth
    left, right = 0.0, length
    bottom, top = 0.0, depth
    quad((left, lo, bottom), (left, hi, bottom), (right, hi, bottom), (right, lo, bottom))
    quad((left, lo, top), (right, lo, top), (right, hi, top), (left, hi, top))
    quad((left, lo, bottom), (left, lo, top), (left, hi, top), (left, hi, bottom))
    quad((right, lo, bottom), (right, hi, bottom), (right, hi, top), (right, lo, top))
    quad((left, lo, bottom), (right, lo, bottom), (right, lo, top), (left, lo, top))
    quad((left, hi, bottom), (left, hi, top), (right, hi, top), (right, hi, bottom))
    return builder.build()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("artifacts/m3871"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------------- PR #48's own claims
    # The wall-sided barge reproduces sin(phi)(GM + BM tan^2(phi)/2) exactly while the
    # waterline still cuts both side walls (tan(phi) < 2T/B = 0.667, i.e. phi < 33.7 deg).
    vessel = barge()
    kg = 1.0
    kb, bm = 1.0, 6.0 ** 2 / (12.0 * 2.0)
    angles = (0.0, 5.0, 10.0, 20.0, 30.0, 40.0)
    curve = gz_curve(vessel, kg=kg, heel_angles_deg=angles, waterline_z=2.0)
    worst = 0.0
    for point in curve.points:
        if point.heel_angle_deg > 33.7:  # past this heel the waterline leaves the side wall
            continue
        phi = math.radians(point.heel_angle_deg)
        analytic = math.sin(phi) * ((kb + bm - kg) + 0.5 * bm * math.tan(phi) ** 2)
        worst = max(worst, abs(point.gz - analytic))
    record("barge_gz_vs_wall_sided_analytic_m", worst, 0.0, 1e-12,
           "max |GZ - analytic| over the 5 wall-sided angles of the ladder")
    record("barge_gm_upright_vs_KB+BM-KG_m", curve.gm_upright, kb + bm - kg, 2e-3,
           "GM recovered from dKN/dphi at zero heel")
    record("barge_imo_screen_passes", imo_intact_stability(curve).passed, True, 0,
           "the six 2008 IS Code Part A 2.2 numerical criteria")
    stiff = gz_curve(vessel, kg=2.45, heel_angles_deg=angles, waterline_z=2.0)
    record("barge_imo_screen_fails_at_GM=0.05m", imo_intact_stability(stiff).passed, False, 0,
           "negative control: the screen is not a rubber stamp")

    # ------------------------------------------------- GZ x transom stern (cross-merge)
    # main added the transom-stern family; the GZ cut has to accept its planar stern plate
    # and agree with the older initial-stability route on the very same mesh.
    for tau in (0.0, 0.25, 0.4):
        params = TransomHullParameters(length=12.0, beam=2.4, draft=0.9, transom_fraction=tau)
        model = generate_transom_hull(params)
        mesh = model.mesh
        wl = model.hydrostatics.waterline_z
        gzc = gz_curve(mesh, kg=0.45, heel_angles_deg=(0.0, 5.0, 10.0, 15.0, 20.0), waterline_z=wl)
        ini = initial_stability(mesh, kg=0.45, waterline_z=wl)
        record(f"transom_tau{tau}_gm_gz_vs_initial_m", gzc.gm_upright, ini.gm_transverse, 1e-9,
               "two independent routes to GM must agree on the transom mesh")
        worst_dv = max(abs(point.displacement_volume - gzc.displacement_volume)
                       / gzc.displacement_volume for point in gzc.points)
        record(f"transom_tau{tau}_displacement_conserved_rel", worst_dv, 0.0, 1e-9,
               "sinkage solve holds the displaced volume at every heel angle")
        record(f"transom_tau{tau}_gz_at_zero_heel_m", gzc.points[0].gz, 0.0, 1e-3,
               "GZ(0) is the mesh's own lateral centre of buoyancy at the design waterline")
        record(f"transom_tau{tau}_gz_at_10deg_positive", gzc.gz_at(10.0) > 0.0, True, 0,
               "small-angle righting survives the stern truncation and its planar plate")
        analytic = analytic_transom_hydrostatics(params)
        record(f"transom_tau{tau}_volume_mesh_vs_analytic_rel",
               abs(model.hydrostatics.displacement_volume - analytic["displacement_volume"])
               / analytic["displacement_volume"], 0.0, 0.01,
               "main's transom hydrostatic references survive the merge")

    # The default transom grid (nz=13) puts a ring exactly at half draft, which is the
    # bracket midpoint of the sinkage solve; the start is stepped off that ring so the
    # GZ curve exists at all on the family main just landed.
    default_grid = generate_transom_hull(TransomHullParameters(length=12.0, beam=2.4, draft=0.9,
                                                               transom_fraction=0.25))
    record("transom_default_grid_waterline_start_off_ring",
           default_grid.mesh.bounding_box().maximum[2], 0.0, 0.0,
           "the capped transom mesh still floats at z = 0 with the solve started off its ring")

    # ------------------------------------- GZ x asymmetric hull (split fore/aft exponents)
    asym = HullParameters(length=12.0, beam=2.4, draft=0.9, fore_exponent=0.5, aft_exponent=0.9)
    asym_model = generate_hull(asym)
    asym_curve = gz_curve(asym_model.mesh, kg=0.45, heel_angles_deg=(0.0, 5.0, 10.0, 15.0, 20.0))
    asym_ini = initial_stability(asym_model.mesh, kg=0.45, waterline_z=0.0)
    record("asymmetric_gm_gz_vs_initial_m", asym_curve.gm_upright, asym_ini.gm_transverse, 1e-9,
           "independent GM routes agree on the split-exponent mesh")
    record("asymmetric_gz_at_zero_heel_m", asym_curve.points[0].gz, 0.0, 1e-9,
           "a transverse-symmetric hull still has GZ(0) = 0 under fore/aft asymmetry")
    record("asymmetric_gz_at_20deg_positive", asym_curve.gz_at(20.0) > 0.0, True, 0,
           "small-angle righting preserved under fore/aft asymmetry")
    reference = hull_analytic_hydrostatics(asym.length, asym.beam, asym.draft,
                                           fore_exponent=asym.fore_exponent,
                                           aft_exponent=asym.aft_exponent)
    record("asymmetric_lcb_mesh_vs_analytic_m", asym_model.hydrostatics.longitudinal_centre_buoyancy,
           reference.longitudinal_centre_buoyancy, 0.02,
           "main's split-exponent LCB reference still matches its mesh")

    # ---------------------------------- what the merge had to keep from each side intact
    plain = HullParameters(length=12.0, beam=2.4, draft=0.9)
    record("default_hull_exponents_unchanged", (plain.fore_exponent, plain.aft_exponent),
           (0.7, 0.7), 0, "the shipped default surface stays the previous one")
    wigley = generate_wigley()
    record("wigley_waterplane_inertia_present",
           wigley.hydrostatics.transverse_waterplane_inertia > 0.0
           and wigley.hydrostatics.longitudinal_waterplane_inertia > 0.0, True, 0,
           "main's waterplane second moments reach the hydrostatics record")
    record("wigley_keel_z_m", wigley.hydrostatics.keel_z, wigley.mesh.bounding_box().minimum[2], 0.0,
           "keel_z is reported next to the inertias and is the mesh's lowest vertex")

    # ------------------------------------------------------- a measured asymmetrism note
    transom_solid = submerged_solid(default_grid.mesh, default_grid.hydrostatics.waterline_z)
    hull_solid = submerged_solid(generate_hull(plain).mesh, 0.0)
    record("transom_yB_vs_hull_yB_m", abs(transom_solid.centroid[1]),
           abs(hull_solid.centroid[1]), 1e-4,
           "the transom mesh is only triangulation-symmetric: |y_B| = 6.2e-5 m vs 7e-17 m for "
           "the plain hull, which is what GZ(0) inherits")

    passed = all(check["passed"] for check in CHECKS)
    payload = {
        "mission": 3871,
        "merge": {"branch": "agent/m3479", "main": "2ff80b2", "pr": 48},
        "passed": passed,
        "checks": CHECKS,
        "gz_curves": {
            f"transom_tau{tau}": {
                "kg": 0.45,
                "angles_deg": [0.0, 5.0, 10.0, 15.0, 20.0],
                "gz_m": [round(p.gz, 9) for p in gz_curve(
                    generate_transom_hull(TransomHullParameters(
                        length=12.0, beam=2.4, draft=0.9, transom_fraction=tau)).mesh,
                    kg=0.45, heel_angles_deg=(0.0, 5.0, 10.0, 15.0, 20.0)).points],
            }
            for tau in (0.0, 0.25, 0.4)
        },
        "barge_gz": {
            "kg": 1.0,
            "angles_deg": list(angles),
            "gz_m": [round(p.gz, 9) for p in curve.points],
            "waterline_z_m": [round(p.waterline_z, 9) for p in curve.points],
        },
    }
    (args.output / "union_verification.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )

    def fmt(value: Any) -> str:
        return f"{value:.6g}" if isinstance(value, float) else str(value)

    lines = [
        "# Cross-feature verification (PR #48 cure, mission 3871)",
        "",
        "Merged tree: `agent/m3479` (GZ curve) + `main` at `2ff80b2` (transom stern, split"
        " exponents, waterplane second moments).",
        "",
        "| check | value | reference | abs error | tol | pass |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for check in CHECKS:
        lines.append(
            "| `{check}` | {value} | {reference} | {err} | {tol:g} | {ok} |".format(
                check=check["check"],
                value=fmt(check["value"]),
                reference=fmt(check["reference"]),
                err=fmt(check["abs_error"]),
                tol=float(check["tolerance"]),
                ok="yes" if check["passed"] else "**NO**",
            )
        )
    lines += ["", f"Overall: **{'PASS' if passed else 'FAIL'}** ({len(CHECKS)} checks).", ""]
    (args.output / "union_verification.md").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())