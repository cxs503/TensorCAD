#!/usr/bin/env python3
"""Scratch: where does the half-breadth field lose the KVLCC2 hull?

Reports (a) the fraction of (x, z) field nodes with data, (b) the stations with no data at
all, (c) the reconstructed midship section area against Cm*B*T and (d) the waterplane area
against the published Cwp band, for the slab that reads the model's own bottom ``draft``.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from tensorcad.fleet.frames import normalise
from tensorcad.fleet.trimfield import sample_trimmed_patches, trim_polygons
from tensorcad.fleet_baseline import field_from_points
from tensorcad.iges import read_iges

IGES = "/root/tp-agent/assets/geometry/hulls/KVLCC2/kvlcc2.igs"
DRAFT, LPP, BEAM = 20.8, 320.0, 58.0


def main() -> int:
    iges = read_iges(IGES)
    points, report = sample_trimmed_patches(trim_polygons(iges), spacing=1.0 / 8.0)
    points, frame_report = normalise(points, {"z_axis": "auto", "x_origin": "auto"})
    zs = [p[2] for p in points]
    depth = max(zs) - min(zs)
    print(json.dumps({"samples": len(points), "depth": round(depth, 3),
                      "trim": report.get("mode")}), flush=True)

    slab_top = depth - DRAFT          # keel = model bottom; immersed slab is the bottom draft
    shifted = [(x, y, z - slab_top) for x, y, z in points]
    field = field_from_points(shifted, draft=DRAFT, dx=1.0, dz=0.5)
    nx, nz = field.nx, field.nz
    xs = [field.x(i) for i in range(nx)]
    # occupancy per station
    stations = []
    for i in range(nx):
        filled = sum(1 for k in range(nz) if field.b[i][k] > 0.0)
        stations.append(filled)
    empty = [round(xs[i], 1) for i in range(nx) if stations[i] == 0]
    print(json.dumps({
        "nx": nx, "nz": nz, "x_lo": round(xs[0], 2), "x_hi": round(xs[-1], 2),
        "cover": round(sum(stations) / (nx * nz), 4),
        "empty_stations": len(empty),
        "empty_station_x_first20": empty[:20],
        "median_filled_per_station": sorted(stations)[len(stations) // 2],
    }, ensure_ascii=False), flush=True)

    # midship and waterline slices
    midship = next(i for i in range(nx) if abs(xs[i] - 160.0) < 0.51)
    area_mid = 0.0
    for k in range(nz):
        b = field.b[midship][k]
        area_mid += 2.0 * b * (field.z(1) - field.z(0))
    waterline_area = 0.0
    for i in range(nx - 1):
        waterline_area += 2.0 * field.b[i][nz - 1] * (xs[i + 1] - xs[i])
    print(json.dumps({
        "midship_station_x": round(xs[midship], 1),
        "midship_area_m2": round(area_mid, 1),
        "Cm_BT_m2": round(0.998 * BEAM * DRAFT, 1),
        "waterplane_area_m2": round(waterline_area, 1),
        "published_waterplane_band_m2": [13000.0, 15500.0],
    }, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())