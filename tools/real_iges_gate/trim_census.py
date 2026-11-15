"""Per-face census of the trimmed reader: does each trim loop actually lie on its surface?

The import pipeline refuses a face whose trim loop cannot be inverted onto the face's
surface within tolerance (``tensorcad.geomquality.trim_mesh.tessellate_trimmed_model``),
and that refusal is what decides whether a publisher IGES is meshed as a *trimmed B-rep*
or whether the pipeline falls back to an untrimmed patch set and silently loses the trim.
When a whole file falls back, the interesting question is which faces were refused and by
how much, so this tool reports the per-face inversion residual and the wall time, without
pretending to judge the geometry.

    python tools/real_iges_gate/trim_census.py \
        --files /root/tp-agent/assets/geometry/hulls/JBC/jbc.igs \
        --out artifacts/real_iges_gate/trim_census.json

The tolerance is the pipeline's own default (``1e-6 * model extent``, file units) unless
``--residual-tolerance`` is given, so the census speaks about the same numbers the gate
does.  Nothing is hard-coded per file: the reader decides the faces.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.geomquality.iges_trim import read_iges_trimmed  # noqa: E402
from tensorcad.geomquality.trim_mesh import (  # noqa: E402
    _loop_to_uv,
    _model_extent_from_faces,
)


def census_one(path: Path, residual_tolerance: float | None) -> dict[str, Any]:
    started = time.time()
    model = read_iges_trimmed(path)
    read_seconds = time.time() - started
    extent = _model_extent_from_faces(model.faces)
    tolerance = (residual_tolerance if residual_tolerance is not None
                 else 1e-6 * max(extent, 1e-9))
    faces: list[dict[str, Any]] = []
    # The pipeline's wall clock is dominated by loop inversion, so it is measured per face.
    for index, face in enumerate(model.faces):
        face_started = time.time()
        worst = 0.0
        lengths: list[int] = []
        for loop in face.loops:
            parameters, residual = _loop_to_uv(face.surface, loop, tolerance)
            worst = max(worst, residual)
            lengths.append(len(parameters))
        faces.append({
            "face": index,
            "surface": type(face.surface).__name__,
            "domain": [list(face.surface.domain[0]), list(face.surface.domain[1])],
            "loops": len(face.loops),
            "loop_points": lengths,
            "max_residual": worst,
            "residual_ratio": worst / tolerance if tolerance else None,
            "accepted": worst <= tolerance,
            "seconds": round(time.time() - face_started, 3),
        })
    refused = [row["face"] for row in faces if not row["accepted"]]
    return {
        "path": str(path),
        "bytes": path.stat().st_size,
        "faces": len(model.faces),
        "extent": extent,
        "residual_tolerance": tolerance,
        "read_seconds": round(read_seconds, 2),
        "loop_seconds": round(sum(row["seconds"] for row in faces), 2),
        "refused_faces": refused,
        "refused_count": len(refused),
        "worst_residual": max((row["max_residual"] for row in faces), default=0.0),
        "faces_detail": faces,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--files", nargs="+", required=True)
    parser.add_argument("--out", default="artifacts/real_iges_gate/trim_census.json")
    parser.add_argument("--residual-tolerance", type=float, default=None)
    args = parser.parse_args(argv)

    records = []
    for name in args.files:
        path = Path(name)
        if not path.exists():
            records.append({"path": str(path), "status": "absent"})
            continue
        record = census_one(path, args.residual_tolerance)
        record["status"] = "measured"
        records.append(record)
        print(f"{path.name}: faces {record['faces']} refused {record['refused_count']} "
              f"worst {record['worst_residual']:.3e} "
              f"(tol {record['residual_tolerance']:.3e}) "
              f"loop {record['loop_seconds']}s")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"records": records}, indent=1))
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())