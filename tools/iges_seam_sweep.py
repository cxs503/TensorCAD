"""Seam-tolerance sweep for the IGES/STEP → watertight mesh import pipeline.

The seam stage is the only repair that can close the T-junction gaps of a publisher
patch set, but its budget trades closure against surface movement: a tolerance wide
enough to close the seam moves the hull by that much.  This tool measures that trade
instead of asserting it: it imports the same file at a ladder of seam tolerances and
records, per rung, the gate verdicts, the surviving open boundary and the gap
classification of :mod:`tools.iges_gap_diag`.

Numbers in ``docs/iges_pipe_seam_study.md`` come from this tool's JSON, not from a
report transcribed by hand.

Usage::

    python tools/iges_seam_sweep.py work/assets/kcs.igs \
        --tolerances none 0.001 0.02 0.2 --out artifacts/iges_seam/kcs.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tools"))

from iges_gap_diag import classify  # noqa: E402
from tensorcad.geomquality.trim_mesh import import_iges  # noqa: E402


def _rung(path: Path, tolerance: float | None, *, seal: bool, eps: float,
          symmetry: bool = False) -> dict:
    started = time.time()
    try:
        result = import_iges(path, seam_tolerance=tolerance, seal_openings=seal,
                             symmetry_closure=symmetry)
    except Exception as exc:  # a rung that cannot run is a result, recorded as one
        return {"seam_tolerance": tolerance, "seal_openings": seal,
                "error": f"{type(exc).__name__}: {exc}", "elapsed_seconds": time.time() - started}
    report = result.report
    gates = {name: {"passed": bool(gate.get("passed")), "value": gate.get("value"),
                    "limit": gate.get("limit")}
             for name, gate in report.get("gates", {}).items()}
    seam = report.get("seam", {})
    stitch = seam.get("stitch", {})
    return {
        "seam_tolerance": tolerance,
        "seal_openings": seal,
        "symmetry_closure": symmetry,
        "elapsed_seconds": round(time.time() - started, 2),
        "gates": gates,
        "gates_passed": sum(1 for gate in gates.values() if gate["passed"]),
        "seam": {
            "applied": seam.get("applied"),
            "splits_total": stitch.get("splits_total"),
            "free_edges_before": stitch.get("free_edges_before"),
            "free_edges_after": stitch.get("free_edges_after"),
            "manifold_degraded": seam.get("manifold_degraded"),
            "aborted": stitch.get("aborted"),
        },
        "seal": report.get("openings", {}).get("seal"),
        "symmetry": {key: value for key, value in report.get("symmetry", {}).items()
                     if key != "detected"} | {"detected_axis":
                     (report.get("symmetry", {}).get("detected") or {}).get("axis_name")},
        "after_repair": report.get("after_repair"),
        "classification": classify(result.mesh.vertices, result.mesh.faces,
                                   eps=eps, seam=tolerance or 0.0),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("path")
    parser.add_argument("--tolerances", nargs="+", default=["none", "0.001", "0.02", "0.2"])
    parser.add_argument("--seal-openings", action="store_true")
    parser.add_argument("--symmetry-closure", action="store_true")
    parser.add_argument("--eps", type=float, default=1e-6)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv)

    rungs = []
    for token in args.tolerances:
        tolerance = None if token.lower() in ("none", "off", "measure") else float(token)
        rung = _rung(Path(args.path), tolerance, seal=args.seal_openings, eps=args.eps,
                     symmetry=args.symmetry_closure)
        rung["free_edges"] = (rung.get("after_repair") or {}).get("free_edges")
        rung["components"] = (rung.get("after_repair") or {}).get("components")
        print(json.dumps({key: rung.get(key) for key in
                          ("seam_tolerance", "gates_passed", "free_edges", "components",
                           "elapsed_seconds", "error")}), flush=True)
        rungs.append(rung)

    payload = {"source": str(args.path), "rungs": rungs}
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=1) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())