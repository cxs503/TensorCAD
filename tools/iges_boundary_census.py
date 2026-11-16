"""Where does the residual open boundary of a publisher hull IGES actually lie?

Given the KCS patch-set report, the surviving free edges either lie on the model's own
cut planes (deck edge, transom, centre plane) or on genuine defect gaps between patches.
This tool reports the census in model coordinates, so the documented verdict -- "the file
is an open surface, not a leaky solid" -- is a measurement rather than an opinion.

Usage::

    python tools/iges_boundary_census.py work/assets/kcs.igs --seam-tolerance 0.02
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from tensorcad.geomquality.stitch import free_edges  # noqa: E402
from tensorcad.geomquality.trim_mesh import import_iges  # noqa: E402


def census(vertices, faces, *, band: float) -> dict:
    edges = free_edges(faces)
    xs = [p[0] for p in vertices]
    ys = [p[1] for p in vertices]
    zs = [p[2] for p in vertices]
    lo = (min(xs), min(ys), min(zs))
    hi = (max(xs), max(ys), max(zs))
    counts = {"centre_plane_y": 0, "deck_z": 0, "keel_z": 0, "transom_x": 0, "bow_x": 0,
              "elsewhere": 0}
    lengths = {"centre_plane_y": 0.0, "deck_z": 0.0, "keel_z": 0.0, "transom_x": 0.0,
               "bow_x": 0.0, "elsewhere": 0.0}
    for a, b in edges:
        mid = tuple((vertices[a][i] + vertices[b][i]) / 2.0 for i in range(3))
        length = sum((vertices[a][i] - vertices[b][i]) ** 2 for i in range(3)) ** 0.5
        if abs(mid[1] - lo[1]) <= band:
            key = "centre_plane_y"
        elif abs(mid[2] - hi[2]) <= band:
            key = "deck_z"
        elif abs(mid[2] - lo[2]) <= band:
            key = "keel_z"
        elif abs(mid[0] - lo[0]) <= band:
            key = "transom_x"
        elif abs(mid[0] - hi[0]) <= band:
            key = "bow_x"
        else:
            key = "elsewhere"
        counts[key] += 1
        lengths[key] += length
    total = len(edges)
    return {
        "free_edges": total,
        "bounding_box": {"lo": lo, "hi": hi},
        "band": band,
        "counts": counts,
        "fraction": {key: (value / total if total else 0.0) for key, value in counts.items()},
        "length": {key: round(value, 3) for key, value in lengths.items()},
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("path")
    parser.add_argument("--seam-tolerance", type=float, default=None)
    parser.add_argument("--band", type=float, default=None,
                        help="bounding-box band in model units (default 0.5%% of extent)")
    parser.add_argument("--json", type=Path, default=None)
    args = parser.parse_args(argv)

    result = import_iges(args.path, seam_tolerance=args.seam_tolerance)
    report = result.report
    box = report["after_repair"]["bounding_box"]
    extent = max(box[1][i] - box[0][i] for i in range(3))
    band = args.band if args.band is not None else 0.005 * extent
    payload = {
        "source": str(args.path),
        "seam_tolerance": args.seam_tolerance,
        "extent": extent,
        "census": census(result.mesh.vertices, result.mesh.faces, band=band),
        "gates": {name: bool(gate.get("passed"))
                  for name, gate in report.get("gates", {}).items()},
    }
    print(json.dumps(payload, indent=1))
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())