#!/usr/bin/env python3
"""Scratch: which vertical slab of the KVLCC2 IGES is the immersed hull?

The published KVLCC2 IGES spans z from -27.07 m to 0 in a waterline-origin, z-down
frame.  The frame route assumes z = 0 *is* the design waterline and keeps the slab
0 ... draft below it, which reads 3.2x too little volume.  This probe measures the
volume of the bottom ``draft`` of the model instead (keel at the model's own extreme)
and scans the slab offset so the assumption is tested, not asserted.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from tensorcad.fleet.frames import normalise
from tensorcad.fleet.trimfield import sample_trimmed_patches, trim_polygons
from tensorcad.fleet_baseline import build_solid, field_from_points
from tensorcad.iges import read_iges

IGES = "/srv/tensor/agent/assets/geometry/hulls/KVLCC2/kvlcc2.igs"
DRAFT = 20.8
REF = 312622.0


def mesh_volume(mesh) -> float:
    volume = 0.0
    for a, b, c in mesh.faces:
        pa, pb, pc = mesh.vertices[a], mesh.vertices[b], mesh.vertices[c]
        ux, uy, uz = pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]
        vx, vy, vz = pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]
        nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        volume += (pa[0] * nx + pa[1] * ny + pa[2] * nz) / 6.0
    return abs(volume)


def main() -> int:
    iges = read_iges(IGES)
    records = trim_polygons(iges)
    points, report = sample_trimmed_patches(records, spacing=1.0 / 4.0)
    points, frame_report = normalise(points, {"z_axis": "auto", "x_origin": "auto"})
    zs = [p[2] for p in points]
    depth = max(zs) - min(zs)
    print(json.dumps({"raw_samples": len(points), "z_min": round(min(zs), 3),
                      "z_max": round(max(zs), 3), "depth": round(depth, 3),
                      "frame": frame_report.as_dict()}, ensure_ascii=False))

    for offset in (0.0, 3.0, 5.0, 6.27, 6.4, 7.0, depth - DRAFT):
        shifted = [(x, y, z - offset) for x, y, z in points]
        field = field_from_points(shifted, draft=DRAFT, dx=1.0, dz=0.1)
        mesh, info = build_solid(field, ring_points=240, draft=DRAFT)
        volume = mesh_volume(mesh)
        print(json.dumps({
            "slab_offset_from_model_top_m": round(offset, 2),
            "volume_m3": round(volume, 1),
            "relative_error": round(volume / REF - 1.0, 5),
            "cb": round(volume / (320.0 * 58.0 * DRAFT), 4),
            "slices": info.get("slices"),
            "max_slice_area_m2": round(info.get("max_slice_area_m2", 0.0), 1),
            "first_slice_area_m2": round(info.get("first_slice_area_m2", 0.0), 1),
            "last_slice_area_m2": round(info.get("last_slice_area_m2", 0.0), 1),
            "cover": field.__dict__.get("cover") if hasattr(field, "__dict__") else None,
        }, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())