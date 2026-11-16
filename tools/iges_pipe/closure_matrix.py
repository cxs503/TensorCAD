"""Import the three real publisher IGES hulls and emit the opening/gate matrix.

Run this on a compute node -- the KVLCC2 and ONRT files take minutes, not seconds::

    python tools/iges_pipe/closure_matrix.py --out artifacts/iges_closure_matrix

Every number in ``matrix.json`` comes from the pipeline's own report; the published
principal dimensions are read from the same ``fleet_specs.json`` the platform's
baseline-fleet gate uses, so nothing here is typed in by hand.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

ASSET_ROOT = Path("/root/tp-agent/assets/geometry/hulls")

FILES = {
    "KCS": ASSET_ROOT / "KCS" / "KCS_hull_SVA.igs",
    "KVLCC2": ASSET_ROOT / "KVLCC2" / "kvlcc2.igs",
    "ONRT": ASSET_ROOT / "ONRT" / "onrt_wo_gap" / "onr5613-full_appendage_deckhous_mergedR.igs",
}


def _published_specs() -> dict:
    path = ROOT / "benchmarks" / "verified" / "baseline_fleet" / "fleet_specs.json"
    if not path.exists():
        return {}
    payload = json.loads(path.read_text())
    ships = payload.get("ships", payload)
    return ships if isinstance(ships, dict) else {}


def _dimensions(vertices) -> dict:
    if not vertices:
        return {}
    lo = [min(point[axis] for point in vertices) for axis in range(3)]
    hi = [max(point[axis] for point in vertices) for axis in range(3)]
    return {
        "bbox_min": lo,
        "bbox_max": hi,
        "length_bbox": hi[0] - lo[0],
        "beam_bbox": hi[1] - lo[1],
        "half_beam_bbox": max(abs(lo[1]), abs(hi[1])),
        "height_bbox": hi[2] - lo[2],
    }


def one(name: str, path: Path, published: dict) -> dict:
    from tensorcad.geomquality.trim_mesh import import_cad

    started = time.time()
    result = import_cad(path, seal_openings=True)
    report = result.report
    openings = report.get("openings", {})
    analysis = openings.get("analysis_before", {})
    after_seal = openings.get("analysis_after", {})
    dims = _dimensions([list(vertex) for vertex in result.mesh.vertices])
    record = {
        "file": str(path),
        "bytes": path.stat().st_size,
        "seconds": round(time.time() - started, 1),
        "triangles_out": len(result.mesh.faces),
        "gates": {key: verdict["passed"] for key, verdict in report["gates"].items()},
        "gate_values": {key: verdict["value"] for key, verdict in report["gates"].items()},
        "signed_volume": report["after_repair"]["signed_volume"],
        "surface_area": report["after_repair"]["area"],
        "dimensions": dims,
        "openings": {
            "verdict": analysis.get("verdict"),
            "free_edges": analysis.get("free_edges"),
            "planar_edges": analysis.get("planar_edges"),
            "defect_edges": analysis.get("defect_edges"),
            "loops": analysis.get("loops", []),
            "loops_after_seal": after_seal.get("loops", []),
            "free_edges_after_seal": after_seal.get("free_edges"),
            "triangles_added": openings.get("seal", {}).get("triangles_added"),
            "capped_area": openings.get("seal", {}).get("capped_area"),
        },
    }
    spec = published.get(name, {}) if isinstance(published, dict) else {}
    if spec.get("beam"):
        record["beam_error_vs_published"] = (abs(dims["beam_bbox"] - spec["beam"]) / spec["beam"])
    if spec.get("lpp"):
        record["length_bbox_over_lpp_error"] = ((dims["length_bbox"] - spec["lpp"]) / spec["lpp"])
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="artifacts/iges_closure_matrix")
    parser.add_argument("--only", nargs="*", default=None)
    args = parser.parse_args(argv)

    published = _published_specs()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix: dict[str, object] = {"published_reference": published, "ships": {}}
    for name, path in FILES.items():
        if args.only and name not in args.only:
            continue
        if not path.exists():
            matrix["ships"][name] = {"file": str(path), "status": "asset_missing"}
            continue
        try:
            matrix["ships"][name] = one(name, path, published)
        except Exception as exc:  # a refusal is a result, not a crash
            matrix["ships"][name] = {"file": str(path),
                                     "status": "refused",
                                     "error": f"{type(exc).__name__}: {exc}"}
        print(name, matrix["ships"][name].get("openings", {}).get("verdict"),
              matrix["ships"][name].get("gates"), flush=True)
    (out_dir / "matrix.json").write_text(json.dumps(matrix, indent=1))
    print(f"wrote {out_dir / 'matrix.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())