#!/usr/bin/env python3
"""Scratch: does the KVLCC2 iges_field reconstruction depend on patch sampling density?"""
import json
import sys
from dataclasses import replace
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from tensorcad.fleet.build import build_ship_geometry
from tensorcad.fleet.manifest import load_manifest

IGES = "/root/tp-agent/assets/geometry/hulls/KVLCC2/kvlcc2.igs"
FRAME = {"z_axis": "auto", "x_origin": "auto", "x_extent": "geometry"}
REF = 312622.0

manifest = load_manifest()
base = next(s for s in manifest.ships if s.name == "KVLCC2")
record = replace(base, geometry_route="iges_field", geometry_path=IGES, geometry_frame=dict(FRAME))

for spacing in (4.0, 8.0):
    built = build_ship_geometry(record, iges_path=IGES, dx=2.0, dz=0.1,
                                ring_points=240, samples_per_metre=spacing)
    vol = abs(built.details.get("signed_volume", 0.0)) or None
    details = built.details
    print(json.dumps({
        "spacing": spacing,
        "dx": 2.0,
        "route": built.route,
        "field": details.get("field"),
        "loft": details.get("loft"),
        "triangles": len(built.mesh.faces),
    }, ensure_ascii=False), flush=True)