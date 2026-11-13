#!/usr/bin/env python3
"""OpenCASCADE parity harness for ``tensorcad.nurbs_intersect`` -- gated, never faked.

The task anchor asks for the NURBS x NURBS intersection to be cross-checked against
OpenCASCADE / pythonocc on at least three cases. OCCT is a large optional binary wheel
(``cadquery-ocp`` for the OCP bindings used elsewhere in this package, or
``pythonocc-core`` from conda) and it is **not** installed in every environment, so two
things are separated here:

* this harness, which, when the wheels are present, builds the same analytic surfaces in
  OCCT, intersects them with ``GeomAPI_IntSS``, discretises the resulting curves and
  reports branch count, length ratio and the two-way deviation from our curves;
* the policy, which is to write ``occt_available: false`` and exit ``2`` when they are
  absent. **No reference numbers are invented in that case** -- the parity claim is then
  WITHHELD, not approximated by the in-repo cross-check (which is a different, weaker
  claim: it shares no code, but it is not an industry kernel).

Usage::

    python tools/xcheck_occt_nurbs_intersect.py [--output artifacts/occt_xcheck.json]

Exit codes: ``0`` parity within tolerance, ``1`` parity broken, ``2`` OCCT not available.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad.adapters.ocp import have_ocp
from tensorcad.nurbs import cylinder_surface, sphere_surface
from tensorcad.nurbs_intersect import nurbs_surface_intersection
from tensorcad.surface_intersection import PlanePatch, SpherePatch

#: Relative tolerance on the length and on the two-way deviation.
LENGTH_TOLERANCE = 1e-4
DEVIATION_TOLERANCE = 1e-4
SAMPLES_PER_CURVE = 4000


def occ_surfaces():
    """The OCCT surface objects, or ``None`` when the bindings are unavailable."""

    if not have_ocp():
        return None
    try:  # pragma: no cover - depends on the optional wheels
        from OCP.Geom import Geom_CylindricalSurface, Geom_Plane, Geom_SphericalSurface
        from OCP.gp import gp_Ax3, gp_Dir, gp_Pnt, gp_Pln
    except Exception:  # pragma: no cover
        return None

    def sphere(radius: float, centre=(0.0, 0.0, 0.0)):
        axis = gp_Ax3(gp_Pnt(*centre), gp_Dir(0.0, 0.0, 1.0))
        return Geom_SphericalSurface(axis, float(radius))

    def plane(point=(0.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)):
        return Geom_Plane(gp_Pln(gp_Pnt(*point), gp_Dir(*normal)))

    def cylinder(radius: float, centre=(0.0, 0.0, 0.0), direction=(0.0, 0.0, 1.0)):
        axis = gp_Ax3(gp_Pnt(*centre), gp_Dir(*direction))
        return Geom_CylindricalSurface(axis, float(radius))

    return {"sphere": sphere, "plane": plane, "cylinder": cylinder}


def occ_divergence(surface_one, surface_two, tolerance: float = 1e-10):
    """Intersect two OCCT surfaces: ``[(points, length), ...]`` per resulting curve."""

    from OCP.GCPnts import GCPnts_QuasiUniformAbscissa
    from OCP.GeomAdaptor import GeomAdaptor_Curve
    from OCP.GeomAPI import GeomAPI_IntSS

    solver = GeomAPI_IntSS(surface_one, surface_two, float(tolerance))
    if not solver.IsDone():  # pragma: no cover - only with a live OCCT
        raise RuntimeError("GeomAPI_IntSS failed to compute an intersection")
    curves = []
    for index in range(1, solver.NbLines() + 1):
        handle = solver.Line(index)
        adaptor = GeomAdaptor_Curve(handle)
        sampled = GCPnts_QuasiUniformAbscissa(adaptor, SAMPLES_PER_CURVE)
        if not sampled.IsDone():  # pragma: no cover
            continue
        points = []
        for k in range(1, sampled.NbPoints() + 1):
            point = adaptor.Value(sampled.Parameter(k))
            points.append((point.X(), point.Y(), point.Z()))
        length = sum(
            math.dist(a, b) for a, b in zip(points, points[1:])
        )
        curves.append({"points": points, "length": length})
    return curves


def polyline_deviation(points, reference) -> float:
    """Largest distance from ``points`` to the polyline ``reference``."""

    worst = 0.0
    segments = list(zip(reference, reference[1:]))
    for point in points:
        best = math.inf
        for a, b in segments:
            ab = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
            ap = (point[0] - a[0], point[1] - a[1], point[2] - a[2])
            denominator = ab[0] ** 2 + ab[1] ** 2 + ab[2] ** 2
            if denominator == 0.0:
                best = min(best, math.dist(point, a))
                continue
            t = (ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2]) / denominator
            t = max(0.0, min(1.0, t))
            foot = (a[0] + t * ab[0], a[1] + t * ab[1], a[2] + t * ab[2])
            best = min(best, math.dist(point, foot))
        worst = max(worst, best)
    return worst


def cases(occ) -> list[dict]:
    """Four pairs, all with a closed-form reference so a third opinion exists."""

    return [
        {
            "name": "sphere R=2 x plane z=0.5",
            "occ": (occ["sphere"](2.0), occ["plane"]((0.0, 0.0, 0.5), (0.0, 0.0, 1.0))),
            "ours": lambda: nurbs_surface_intersection(
                sphere_surface(2.0),
                PlanePatch(origin=(-3.0, -3.0, 0.5), du=(6.0, 0.0, 0.0), dv=(0.0, 6.0, 0.0)),
            ),
            "scale": 2.0,
        },
        {
            "name": "sphere R=2 x sphere R=1.2 at z=1",
            "occ": (occ["sphere"](2.0), occ["sphere"](1.2, (0.0, 0.0, 1.0))),
            "ours": lambda: nurbs_surface_intersection(
                sphere_surface(2.0), SpherePatch(radius=1.2, centre=(0.0, 0.0, 1.0))
            ),
            "scale": 2.0,
        },
        {
            "name": "cylinder R=1 x plane z=0.3",
            "occ": (occ["cylinder"](1.0), occ["plane"]((0.0, 0.0, 0.3), (0.0, 0.0, 1.0))),
            "ours": lambda: nurbs_surface_intersection(
                cylinder_surface(1.0, 4.0, base_z=-2.0),
                PlanePatch(origin=(-2.0, -2.0, 0.3), du=(4.0, 0.0, 0.0), dv=(0.0, 4.0, 0.0)),
            ),
            "scale": 1.0,
        },
        {
            "name": "cylinder R=1 x sphere R=1.5 at z=0",
            "occ": (occ["cylinder"](1.0), occ["sphere"](1.5)),
            "ours": lambda: nurbs_surface_intersection(
                cylinder_surface(1.0, 4.0, base_z=-2.0), SpherePatch(radius=1.5)
            ),
            "scale": 1.5,
        },
    ]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/occt_xcheck.json")
    arguments = parser.parse_args(argv)

    output = Path(arguments.output)
    if not output.is_absolute():
        output = REPO_ROOT / output
    output.parent.mkdir(parents=True, exist_ok=True)

    occ = occ_surfaces()
    if occ is None:
        payload = {
            "occt_available": False,
            "parity": "WITHHELD",
            "reason": "neither the OCP (cadquery-ocp) nor the pythonocc bindings import "
                      "in this environment; no OCCT numbers were generated and none are "
                      "reported. The in-repo cross-check against the implicit-field "
                      "marcher (see docs/nurbs_intersect.md section 3.4) is a different, "
                      "weaker claim and is not a substitute for it.",
            "cases": [],
        }
        output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print("OCCT bindings are not installed: parity WITHHELD, no numbers invented.")
        print(f"recorded: {output}")
        return 2

    rows = []
    broken = []
    for case in cases(occ):  # pragma: no cover - needs the optional wheels
        reference_curves = occ_divergence(*case["occ"])
        ours = case["ours"]()
        ours_points = [point for curve in ours.curves for point in curve.points]
        occ_points = [point for curve in reference_curves for point in curve["points"]]
        occ_length = sum(curve["length"] for curve in reference_curves)
        length_ratio = ours.total_length / occ_length if occ_length else math.inf
        deviation_ours = polyline_deviation(occ_points, ours_points) if ours_points else math.inf
        deviation_theirs = polyline_deviation(ours_points, occ_points) if occ_points else math.inf
        deviation = max(deviation_ours, deviation_theirs) / case["scale"]
        row = {
            "case": case["name"],
            "occ_curves": len(reference_curves),
            "our_curves": len(ours.curves),
            "occ_length": occ_length,
            "our_length": ours.total_length,
            "length_ratio": length_ratio,
            "two_way_deviation_over_scale": deviation,
            "agree": bool(
                len(reference_curves) == len(ours.curves)
                and abs(length_ratio - 1.0) <= LENGTH_TOLERANCE
                and deviation <= DEVIATION_TOLERANCE
            ),
        }
        rows.append(row)
        if not row["agree"]:
            broken.append(case["name"])
        print(f"{case['name']}: OCC {occ_length:.9f} vs ours {ours.total_length:.9f} "
              f"(ratio {length_ratio:.9f}), deviation {deviation:.2e} of scale")

    payload = {
        "occt_available": True,
        "parity": "BROKEN" if broken else "AGREED",
        "tolerance": {"length_relative": LENGTH_TOLERANCE,
                      "deviation_over_scale": DEVIATION_TOLERANCE},
        "cases": rows,
    }
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"recorded: {output}")
    return 1 if broken else 0


if __name__ == "__main__":
    raise SystemExit(main())