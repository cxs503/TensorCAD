"""Run the mesh-derived voxelisation ladders and commit the artifact.

Two sweeps are merged into one evidence record so both gates can be fitted from
a single artifact:

A  triangle axis at one fixed grid (reaches the 10^7 band by subdivision)
B  cell axis at one fixed mesh

The complete artifact JSON is printed to stdout as well as written to disk, so
the numbers survive even if only the job log is retrievable.
"""
from __future__ import annotations

import json
from pathlib import Path

from tensorcad.perf_voxel_mesh import SCHEMA, run_voxelize_mesh_benchmark, scaling_gates

OUT = Path("benchmarks/perf_scale/artifacts/ladder_voxelize_mesh.json")


def main() -> int:
    a = run_voxelize_mesh_benchmark(
        bases=((41, 13),),
        levels_ladder=(0, 3, 4, 5, 6),
        grids=((16, 8, 4),),
        band_factor=1.0,
    )
    b = run_voxelize_mesh_benchmark(
        bases=((41, 13),),
        levels_ladder=(4,),
        grids=((8, 4, 2), (24, 12, 6), (32, 16, 8)),
        band_factor=1.0,
    )
    rows = a["results"] + b["results"]
    gates = scaling_gates(rows)
    warnings = sorted(set(a["cluster_warnings"]) | set(b["cluster_warnings"]))
    artifact = {
        "schema": SCHEMA,
        "passed": all(g["passed"] for g in gates if g["status"] == "measured"),
        "verdict": "pass" if all(g["passed"] for g in gates if g["status"] == "measured") else "fail",
        "gates_checked": sum(1 for g in gates if g["status"] == "measured"),
        "gates_not_run": sum(1 for g in gates if g["status"] == "not_run"),
        "seconds": round(a["seconds"] + b["seconds"], 3),
        "band_factor": 1.0,
        "solver": a["solver"],
        "environment": a["environment"],
        "ladder_a": {"bases": a["bases"], "levels": a["levels"], "grids": a["grids"]},
        "ladder_b": {"bases": b["bases"], "levels": b["levels"], "grids": b["grids"]},
        "cluster_warnings": warnings,
        "gates": gates,
        "results": rows,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n")
    print("=== ARTIFACT BEGIN ===")
    print(json.dumps(artifact, indent=2, ensure_ascii=False))
    print("=== ARTIFACT END ===")
    print(f"wrote {OUT} rows={len(rows)} verdict={artifact['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())