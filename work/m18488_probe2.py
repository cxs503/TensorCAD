#!/usr/bin/env python3
"""Probe per-surface control-net bounding boxes to spot local frames / transforms."""
import sys, json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from tensorcad.iges import read_iges  # noqa: E402

FILES = {
    "JBC": "artifacts/baseline_fleet/raw/JBC/jbc.igs",
    "KVLCC2": "artifacts/baseline_fleet/raw/KVLCC2/kvlcc2.igs",
    "ONRT": "artifacts/baseline_fleet/raw/ONRT/onr5613-full_appendage_deckhous_mergedR.igs",
}

for name, rel in FILES.items():
    model = read_iges(REPO / rel)
    print("=" * 20, name, "surfaces", len(model.surfaces()), flush=True)
    rows = []
    for i, s in enumerate(model.surfaces()):
        pts = getattr(s, "control_points", None)
        if pts is None:
            rows.append((i, "no-control-points", {}))
            continue
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]; zs = [p[2] for p in pts]
        rows.append((i, len(pts), {
            "x": [round(min(xs), 2), round(max(xs), 2)],
            "y": [round(min(ys), 2), round(max(ys), 2)],
            "z": [round(min(zs), 2), round(max(zs), 2)],
            "layout": getattr(s, "layout", "?"),
            "amb": getattr(s, "layout_ambiguous", None),
        }))
    for r in rows:
        print(r, flush=True)
    attrs = [a for a in dir(model) if not a.startswith("_")]
    print("model attrs:", attrs, flush=True)