"""Adjudicate the two contradictory KVLCC2 displacement measurements.

Background
----------
``GAP-CAD-T1-IGES-PIPE`` needs KVLCC2 (the second publisher IGES) to pass the physical
3 % gate.  Two routes through the *same file* currently disagree:

* ``tensorcad.fleet_baseline.build_ship`` (dense in-trim sampling, watertight mesh)
  reported ``nabla = 312 568.0 m3`` vs the published 312 622 m3 (**-0.017 %**,
  ``docs/gap_cad_t1_iges_pipe_kvlcc2_density.md``);
* an exact divergence-theorem integral over the *unrepaired* patch set
  (``work/m22903 tools/baseline_fleet/kvlcc2_divergence_probe.py``) reported
  ``256 657 m3`` (**-17.90 %**) with a closure residual ``|∮n dA| = 11 272 m2``, and
  concluded KVLCC2 is ``WITHHELD`` because the patch set is not a closed,
  consistently oriented half-shell.

An open half-breadth hull *skin* cannot have ``∮n dA = 0``: the centreplane cut and the
design-waterline cap are simply not part of a hull-surface file (KCS behaves the same
way).  So a non-zero closure residual alone does **not** prove the geometry is short.
This script measures, on one import, the quantities that decide between the two
readings:

1. the *closure-free* displacement from :func:`hull_surface.hull_characteristics`
   (envelope-breadth integral; orientation is irrelevant to it), with a grid ladder;
2. the divergence-theorem integral over the mesh **as stored**, and again after
   :func:`mesh_repair.orient_consistently` plus a global outward flip;
3. the **cap correction**: for the half-solid, the centreplane cap at ``y = 0``
   contributes ``r·n = y = 0`` and therefore nothing, so the only missing cap is the
   design waterline, whose term is ``(1/3) * T * A_wl``.  If
   ``V_skin + (1/3) T A_wl`` lands on the published volume, then the -17.9 % is a
   missing-cap artefact of the measurement route, not a defect of the file.

Usage
-----
    python scripts/kvlcc2_route_adjudication.py --iges <kvlcc2.igs> \
        --out artifacts/gap_t1_kvlcc2_routes/evidence.json [--u-samples 96 --v-samples 48]

Everything is read-only with respect to the publisher bytes; the JSON records the file
name, its sha256 and the grid settings so the numbers can be re-derived by anyone who
holds the same file (the bytes themselves do not enter this repository).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

from tensorcad.fleet.frames import normalise
from tensorcad.geometry import Mesh
from tensorcad.geomquality.hull_surface import (
    PARTICULARS,
    _Slabs,
    classify_model,
    envelope_breadth,
    hull_characteristics,
    validate_against_particulars,
)
from tensorcad.geomquality.trim_mesh import import_cad
from tensorcad.mesh_repair import measure_mesh, orient_consistently


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(u, v):
    return (u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0])


def _dot(u, v):
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def divergence_integrals(vertices, faces) -> dict[str, float]:
    """``∮ n dA`` (vector) and ``(1/3) ∮ r·n dA`` (signed volume) for the faces as given."""
    nx = ny = nz = 0.0
    volume = 0.0
    for face in faces:
        a, b, c = (vertices[i] for i in face)
        cross = _cross(_sub(b, a), _sub(c, a))
        nx += cross[0] * 0.5
        ny += cross[1] * 0.5
        nz += cross[2] * 0.5
        # r·(n dA) with r = a: exact for a planar triangle
        volume += _dot(a, cross) / 6.0
    return {
        "closure_x": nx, "closure_y": ny, "closure_z": nz,
        "closure_norm": math.sqrt(nx * nx + ny * ny + nz * nz),
        "volume_signed": volume,
    }


def flip_if_inward(vertices, faces):
    """Global outward flip: the half-solid's skin must give a positive volume."""
    integrals = divergence_integrals(vertices, faces)
    if integrals["volume_signed"] >= 0.0:
        return faces, False
    return [tuple(reversed(face)) for face in faces], True


def waterplane_area(vertices, faces, waterline: float, stations: int = 201) -> float:
    """``∫ b(x, T) dx`` on the imported model (half-breadth: one side only)."""
    slabs = _Slabs.build(vertices, faces, 256)
    characteristics = classify_model(vertices)
    full = not characteristics["half_breadth"]
    xs = [min(p[0] for p in tri) for tri in slabs.triangles]
    x_lo = min(xs)
    x_hi = max(max(p[0] for p in tri) for tri in slabs.triangles)
    step = (x_hi - x_lo) / stations
    total = 0.0
    for index in range(stations):
        station = x_lo + (index + 0.5) * step
        total += envelope_breadth(slabs, station, waterline, full=full) * step
    return total


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iges", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--u-samples", type=int, default=96)
    parser.add_argument("--v-samples", type=int, default=48)
    parser.add_argument("--cache", type=Path, default=None,
                        help="optional JSON mesh cache so the heavy import runs once")
    args = parser.parse_args(argv)

    payload: dict[str, object] = {
        "file": args.iges.name,
        "sha256": sha256_of(args.iges),
        "samples": {"u": args.u_samples, "v": args.v_samples},
        "waterline": PARTICULARS["KVLCC2"]["draft"],
    }

    mesh_payload = None
    if args.cache and args.cache.exists():
        mesh_payload = json.loads(args.cache.read_text())
    if mesh_payload is None:
        started = time.time()
        imported = import_cad(args.iges, weld_tolerance=1e-3,
                              u_samples=args.u_samples, v_samples=args.v_samples)
        payload["import_seconds"] = round(time.time() - started, 1)
        payload["import_report"] = {key: value for key, value in imported.report.items()
                                    if key in ("route", "gates", "after_repair", "before_repair",
                                               "faces_in_transfer_band", "max_inversion_residual")}
        mesh_payload = {"vertices": [list(v) for v in imported.mesh.vertices],
                        "faces": [list(f) for f in imported.mesh.faces]}
        if args.cache:
            args.cache.parent.mkdir(parents=True, exist_ok=True)
            args.cache.write_text(json.dumps(mesh_payload))
    vertices = [tuple(v) for v in mesh_payload["vertices"]]
    faces = [tuple(f) for f in mesh_payload["faces"]]
    payload["triangles"] = len(faces)
    payload["vertices"] = len(vertices)

    # ---- frame + physical characteristics (route C: closure-free) -------------- #
    frame_points, frame_report = normalise(vertices)
    payload["frame"] = {"z_flip": frame_report.z_flip, "x_shift": frame_report.x_shift}
    payload["classify"] = {k: v for k, v in classify_model(frame_points).items()
                           if k in ("half_breadth", "y_min", "y_max", "symmetry_factor")}
    particulars = PARTICULARS["KVLCC2"]
    ladder = {}
    for stations, levels in ((41, 12), (61, 24), (61, 48), (101, 48)):
        started = time.time()
        got = hull_characteristics(frame_points, faces, waterline=particulars["draft"],
                                   stations=stations, levels=levels)
        report = validate_against_particulars(got, particulars, tolerance_percent=3.0)
        ladder[f"{stations}x{levels}"] = {
            "volume": got.volume,
            "cb": got.cb,
            "beam": got.beam,
            "s_wetted": got.s_wetted,
            "s_wetted_naive": got.s_wetted_naive,
            "s_wetted_discarded": got.s_wetted_discarded,
            "lwl": got.lwl,
            "errors_percent": {k: row["error_percent"] for k, row in report["quantities"].items()},
            "passed": report["passed"],
            "seconds": round(time.time() - started, 1),
        }
        payload["characteristics_ladder"] = ladder
        print(f"[{stations}x{levels}] volume={got.volume:.1f} m3 "
              f"({ladder[f'{stations}x{levels}']['errors_percent'].get('volume', float('nan')):+.3f} %) "
              f"cb={got.cb:.4f} S={got.s_wetted:.1f} m2", flush=True)

    # ---- divergence theorem on the mesh as stored, then repaired --------------- #
    as_stored = divergence_integrals(vertices, faces)
    payload["divergence_as_stored"] = as_stored
    repaired, step = orient_consistently(Mesh(vertices=tuple(vertices), faces=tuple(faces),
                                              units="m"))
    rep_vertices = list(repaired.vertices)
    rep_faces = list(repaired.faces)
    rep_flipped, global_flip = flip_if_inward(rep_vertices, rep_faces)
    if global_flip:
        rep_faces = rep_flipped
    payload["orientation_repair"] = {
        "step": getattr(step, "kind", str(step)),
        "faces_reversed": getattr(step, "faces_reversed", None),
        "global_flip_applied": global_flip,
        "divergence_after_repair": divergence_integrals(rep_vertices, rep_faces),
    }
    try:
        metrics = measure_mesh(Mesh(vertices=tuple(rep_vertices), faces=tuple(rep_faces), units="m"))
        payload["mesh_metrics"] = {name: getattr(metrics, name)
                                   for name in ("watertight", "manifold", "winding_conflicts",
                                                "components", "boundary_edges")
                                   if hasattr(metrics, name)}
    except Exception as exc:  # pragma: no cover - diagnostics only
        payload["mesh_metrics"] = {"error": str(exc)}

    # ---- cap correction: the only missing cap of a half-solid is the waterline -- #
    waterline = particulars["draft"]
    area_wl = waterplane_area(frame_points, faces, waterline)
    half_volume = as_stored["volume_signed"] if as_stored["volume_signed"] < 0.0 else -as_stored["volume_signed"]
    # Take the magnitude of the stored-orientation integral as the skin term; the sign
    # bookkeeping of an open shell is not meaningful, its magnitude is.
    skin = abs(as_stored["volume_signed"])
    sign = 1.0 if as_stored["volume_signed"] > 0.0 else -1.0
    cap_term = waterline * area_wl / 3.0
    payload["cap_correction"] = {
        "waterplane_area_half": area_wl,
        "cap_term_m3": cap_term,
        "skin_integral_m3": as_stored["volume_signed"],
        "skin_magnitude_m3": skin,
        "corrected_half_m3_same_sign": skin + cap_term * (1.0 if sign > 0 else -1.0),
        "corrected_full_m3": 2.0 * abs(skin + cap_term * (1.0 if sign > 0 else -1.0)),
        "published_full_m3": particulars["volume"],
        "skin_only_full_m3": 2.0 * skin,
    }
    payload["cap_correction"]["corrected_error_percent"] = (
        100.0 * (payload["cap_correction"]["corrected_full_m3"] - particulars["volume"])
        / particulars["volume"])
    payload["cap_correction"]["skin_only_error_percent"] = (
        100.0 * (payload["cap_correction"]["skin_only_full_m3"] - particulars["volume"])
        / particulars["volume"])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(json.dumps(payload["cap_correction"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())