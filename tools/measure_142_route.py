"""Before/after measurement of the type-142 route (field position vs token scan).

Usage (from the repository root)::

    python tools/measure_142_route.py before.json
    python tools/measure_142_route.py after.json

The script writes a JSON fingerprint of

* every 142 record in the archives that are present (fields actually used, loops per face),
* the tessellation outcome of the propeller fixture (triangles, refused faces),
* the KCS hull import (faces, gates, volume, free edges) when the asset is present.

It is deliberately read-only and does not touch the platform; the point is that the two
JSON files differ only by the reader under test, so a claim like "the phantom loop is
gone" is a diff, not an assertion.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from tensorcad.geomquality.iges_trim import read_iges_trimmed  # noqa: E402
from tensorcad.geomquality.trim_mesh import tessellate_trimmed_model  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ARCHIVES = {
    "KP505": ROOT / "benchmarks/anchors/sources/KP505_modelScale_250mm.igs",
    "KVLCC2": Path("/srv/tensor/agent/assets/geometry/hulls/KVLCC2/kvlcc2.igs"),
    "KCS": ROOT / "work/assets/kcs.igs",
}


def fingerprint_model(model, elapsed: float) -> dict:
    return {
        "elapsed_s": round(elapsed, 2),
        "faces": len(model.faces),
        "loops_per_face": [len(face.loops) for face in model.faces],
        "routes": [face.trim_routes() for face in model.faces],
        "cross_checks": [
            [None if value is None else round(value, 6) for value in face.cross_checks]
            for face in model.faces
        ],
        "refused_entries": len(model.refused),
    }


def main() -> int:
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "fingerprint.json")
    report: dict[str, object] = {}
    for name, path in ARCHIVES.items():
        if not path.exists():
            report[name] = {"skipped": "asset not present"}
            continue
        started = time.time()
        try:
            model = read_iges_trimmed(str(path))
        except Exception as exc:  # noqa: BLE001 - the failure itself is the measurement
            report[name] = {"error": f"{type(exc).__name__}: {exc}",
                            "elapsed_s": round(time.time() - started, 2)}
            continue
        entry = fingerprint_model(model, time.time() - started)
        if name == "KP505":
            tessellation = tessellate_trimmed_model(model)
            entry["tessellation"] = {
                "triangles": len(tessellation.faces),
                "refused_faces": sorted(tessellation.refused),
                "refusal_reasons": {str(k): str(v)[:200] for k, v in tessellation.refused.items()},
            }
        report[name] = entry
    out.write_text(json.dumps(report, indent=2, sort_keys=True))
    print(f"wrote {out}")
    for name, entry in report.items():
        if "error" in entry:
            print(f"  {name}: ERROR {entry['error']}")
        elif "skipped" in entry:
            print(f"  {name}: skipped")
        else:
            loops = entry["loops_per_face"]
            print(f"  {name}: {entry['faces']} faces, loops/face={loops}, "
                  f"{entry['elapsed_s']}s")
            if "tessellation" in entry:
                tess = entry["tessellation"]
                print(f"      triangles={tess['triangles']} refused={tess['refused_faces']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())