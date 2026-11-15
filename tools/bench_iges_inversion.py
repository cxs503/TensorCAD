"""Regenerate the GAP-CAD-T1-IGES-PIPE inversion-speedup measurements.

Runs the same workload with the fix stashed (baseline) and applied (candidate) so the
report table is generated, never hand-copied.  Writes artifacts/gap_t1_iges_inversion/
perf.json.

Usage:  python tools/bench_iges_inversion.py [--quick]
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "work" / "assets"
OUT = ROOT / "artifacts" / "gap_t1_iges_inversion"


def micro(kvlcc2: Path, faces: int = 2) -> dict:
    sys.path.insert(0, str(ROOT))
    from tensorcad.geomquality.iges_trim import read_iges_trimmed
    from tensorcad.geomquality.trim_mesh import _loop_to_uv, _model_extent_from_faces

    model = read_iges_trimmed(str(kvlcc2))
    extent = _model_extent_from_faces(model.faces)
    tolerance = 1e-6 * max(extent, 1e-9)
    loops = [(f.surface, l) for f in model.faces[:faces] for l in f.loops]
    points = 0
    worst = 0.0
    t0 = time.perf_counter()
    for surface, loop in loops:
        params, residual = _loop_to_uv(surface, loop, tolerance)
        points += len(params)
        worst = max(worst, residual)
    seconds = time.perf_counter() - t0
    return {
        "faces": faces,
        "loops": len(loops),
        "loop_points": points,
        "inversion_s": round(seconds, 3),
        "us_per_point": round(1e6 * seconds / max(points, 1)),
        "worst_residual": float(f"{worst:.6e}"),
        "tolerance": float(f"{tolerance:.6e}"),
    }


def kcs_import(kcs: Path) -> dict:
    sys.path.insert(0, str(ROOT))
    from tensorcad.geomquality.trim_mesh import import_cad

    t0 = time.perf_counter()
    result = import_cad(str(kcs), weld_tolerance=1e-3, seam_tolerance=1e-3)
    seconds = time.perf_counter() - t0
    gates = {name: bool(v["passed"]) for name, v in result.report["gates"].items()}
    return {
        "seconds": round(seconds, 2),
        "triangles": len(result.mesh.faces),
        "vertices": len(result.mesh.vertices),
        "refused": len(result.report["import"].get("refused", {})),
        "gates": gates,
    }


def run_both(label: str) -> dict:
    micro_result = micro(ASSETS / "kvlcc2.igs")
    kcs_result = kcs_import(ASSETS / "kcs.igs")
    return {"label": label, "inversion_micro": micro_result, "kcs_import": kcs_result}


def main() -> int:
    if not (ASSETS / "kvlcc2.igs").exists():
        print(f"assets missing under {ASSETS}", file=sys.stderr)
        return 2
    candidate = run_both("with_derivative_net_cache")
    patched = {}
    for rel in ("tensorcad/nurbs.py", "tensorcad/geomquality/trim_mesh.py"):
        patched[rel] = (ROOT / rel).read_bytes()
    try:
        subprocess.run(["git", "checkout", "--", *patched], check=True, cwd=ROOT)
        baseline = run_both("baseline_from_head")
    finally:
        for rel, blob in patched.items():
            (ROOT / rel).write_bytes(blob)
    payload = {
        "baseline": baseline,
        "candidate": candidate,
        "speedup_us_per_point": round(
            baseline["inversion_micro"]["us_per_point"]
            / max(candidate["inversion_micro"]["us_per_point"], 1), 2),
        "speedup_kcs_import": round(
            baseline["kcs_import"]["seconds"] / max(candidate["kcs_import"]["seconds"], 1e-9), 2),
        "identical_kcs_mesh": (
            baseline["kcs_import"]["triangles"] == candidate["kcs_import"]["triangles"]
            and baseline["kcs_import"]["vertices"] == candidate["kcs_import"]["vertices"]
            and baseline["kcs_import"]["gates"] == candidate["kcs_import"]["gates"]
        ),
        "identical_worst_residual": (
            baseline["inversion_micro"]["worst_residual"]
            == candidate["inversion_micro"]["worst_residual"]
        ),
    }
    if payload["speedup_us_per_point"] < 1.2:
        (OUT / "perf.json").unlink(missing_ok=True)
        print("BASELINE NOT REVERTED (speedup %.2f): artifact withheld" % payload["speedup_us_per_point"])
        return 3
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "perf.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())