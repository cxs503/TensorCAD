#!/usr/bin/env python3
"""Diagnostic probe: per-hull IGES read / patch census / sample cost (m18488)."""
import sys, time, json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from tensorcad.iges import read_iges
from tensorcad.fleet_baseline import unit_scale, partition_surfaces

FILES = {
    "KCS": "artifacts/baseline_fleet/raw/KCS/KCS_hull_SVA.igs",
    "JBC": "artifacts/baseline_fleet/raw/JBC/jbc.igs",
    "KVLCC2": "artifacts/baseline_fleet/raw/KVLCC2/kvlcc2.igs",
    "ONRT": "artifacts/baseline_fleet/raw/ONRT/onr5613-full_appendage_deckhous_mergedR.igs",
}

out = {}
for name, rel in FILES.items():
    p = REPO / rel
    t0 = time.time()
    model = read_iges(p)
    t1 = time.time()
    surfaces = model.surfaces()
    t2 = time.time()
    usable, quarantined = partition_surfaces(surfaces)
    t3 = time.time()
    rec = {
        "bytes": p.stat().st_size,
        "read_s": round(t1 - t0, 2),
        "surfaces_s": round(t2 - t1, 2),
        "partition_s": round(t3 - t2, 2),
        "n_surfaces": len(surfaces),
        "n_usable": len(usable),
        "n_quarantined": len(quarantined),
        "unit_flag": model.unit_flag,
        "unit_name": model.unit_name,
        "unit_scale": unit_scale(model.unit_flag),
        "model_scale": model.model_scale,
        "units": model.units,
        "warnings": list(model.unit_warnings)[:5],
    }
    # bounding box of a coarse sample of the usable patches
    try:
        from tensorcad.fleet_baseline import sample_patches
        t4 = time.time()
        pts = list(sample_patches(usable, spacing=1.0))
        t5 = time.time()
        xs = [q[0] for q in pts]; ys = [q[1] for q in pts]; zs = [q[2] for q in pts]
        rec["sample_s"] = round(t5 - t4, 2)
        rec["n_points"] = len(pts)
        rec["bbox"] = [round(min(xs), 3), round(max(xs), 3),
                       round(min(ys), 3), round(max(ys), 3),
                       round(min(zs), 3), round(max(zs), 3)]
    except Exception as exc:
        rec["sample_error"] = f"{type(exc).__name__}: {exc}"
    out[name] = rec
    print(name, json.dumps(rec), flush=True)

(REPO / "work" / "m18488_probe.json").write_text(json.dumps(out, indent=1) + "\n")