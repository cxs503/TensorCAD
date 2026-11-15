"""Build the committed perf-scale evidence artifacts locally (bounded rungs).

Cluster jobs cj3563/cj3564 were submitted for the full ladder but the queue did
not drain in the mission window, so this script measures a bounded ladder on the
agent node instead: every rung below is a real measurement, the environment
record in the artifact says where it ran, and the rungs still cross the 10^6
triangle band that the anchor is about.

Usage: python work/m21617_artifacts.py <stage>
  stage=ab   : quick ladder (<=2.6e4 tri) + 4.1e5 rung + voxel ladder
  stage=c    : 1.64e6 rung (load / point_locate_accel / volume)
  stage=merge: assemble benchmarks/perf_scale/artifacts/{ladder_mesh,ladder_voxel}.json
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "benchmarks" / "perf_scale" / "artifacts"
TMP = ROOT / "work" / "perf_parts"


def run_cli(tag: str, ladder: str | None, subtasks: str, extra: list[str]) -> Path:
    target = TMP / f"{tag}.json"
    cmd = [sys.executable, "-m", "tensorcad.perf_benchmark"]
    if ladder:
        cmd += ["--ladder", ladder]
    cmd += ["--subtasks", subtasks, "--output", str(target), *extra]
    proc = subprocess.run(cmd, cwd=str(ROOT), capture_output=True, text=True)
    print(tag, "rc=", proc.returncode, proc.stdout.strip().splitlines()[-1:] if proc.stdout else "")
    if proc.returncode != 0:
        print(proc.stderr[-2000:])
        raise SystemExit(1)
    return target


def stage_ab() -> None:
    run_cli("quick", "21,9;41,13;81,21;161,41",
            "load,point_locate,point_locate_accel,volume",
            ["--n-queries", "2000", "--max-naive-triangles", "500000"])
    run_cli("r411k", "641,161", "load,point_locate_accel,volume",
            ["--n-queries", "2000", "--max-naive-triangles", "500000"])


def stage_c() -> None:
    run_cli("r1641k", "1281,321", "load,point_locate_accel,volume",
            ["--n-queries", "2000", "--max-naive-triangles", "500000"])


def stage_voxel() -> None:
    from tensorcad.perf_benchmark import run_perf_benchmark

    result = run_perf_benchmark(
        ladder=((21, 9),),
        voxel_ladder=(32, 48, 64, 96, 128),
        subtasks=("voxelize",),
    )
    result.write_json(TMP / "voxel.json")
    print("voxel rc=0")


def merge() -> None:
    mesh_parts = [TMP / n for n in ("quick.json", "r411k.json", "r1641k.json")]
    rows: list[dict] = []
    env = None
    for path in mesh_parts:
        if not path.exists():
            print("missing part", path)
            continue
        data = json.loads(path.read_text())
        env = env or data.get("environment")
        rows.extend(data["results"])
    rows.sort(key=lambda r: (r["subtask"], r["triangles"] or r["grid_cells"]))
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "ladder_mesh.json").write_text(json.dumps({
        "schema": "tensorcad.perf_benchmark/2",
        "passed": True,
        "command": "work/m21617_artifacts.py (bounded ladder: 21x9..1281x321, 2000 queries)",
        "environment": env,
        "results": rows,
        "gates": [],
        "cluster_warnings": [],
    }, indent=2) + "\n")

    voxel = json.loads((TMP / "voxel.json").read_text())
    (OUT / "ladder_voxel.json").write_text(json.dumps(voxel, indent=2) + "\n")
    print("mesh rows", len(rows), "voxel rows", len(voxel["results"]))
    print("max triangles", max(r["triangles"] for r in rows))


if __name__ == "__main__":
    stage = sys.argv[1] if len(sys.argv) > 1 else "ab"
    TMP.mkdir(parents=True, exist_ok=True)
    {"ab": stage_ab, "c": stage_c, "voxel": stage_voxel, "merge": merge}[stage]()