#!/usr/bin/env python3
"""Board audit harness for m17476 (GAP-CAD-B3-SOLVER-INTERFACE) -- machine re-checkable.

Audited subject: TensorCAD worktree /root/tp-agent/work/m17476 @ agent/m17476
(base 16eefd3, deliverables uncommitted at audit time).

What this script does (all read-only with respect to the subject worktree):

  default      static evidence re-check: file/line pointers, placeholder scan, frame-contract
               refusal behaviour, artifact-vs-reference arithmetic recomputed from
               artifacts/m17476/solver_interface/results.json, caller grep.
  --pytest     additionally runs the subject's test file and records the outcome.
  --repro      additionally re-runs the shipped repro command into a scratch directory and
               compares every comparison against the shipped artifact (takes ~250 s; run it
               in the background, the interactive shell window is 120 s).

Emits one JSON document on stdout so the verdict can be re-derived from raw evidence.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import subprocess
import sys
from pathlib import Path

SUBJECT = Path("/root/tp-agent/work/m17476")
ARTIFACTS = SUBJECT / "artifacts/m17476/solver_interface"
IMPL = SUBJECT / "tensorcad/solver_interface.py"
TESTS = SUBJECT / "tests/test_solver_interface.py"
DOCS = SUBJECT / "docs/solver_interface.md"
BENCH = SUBJECT / "tensorcad/benchmark.py"
INIT = SUBJECT / "tensorcad/__init__.py"

PLACEHOLDER = re.compile(r"NotImplementedError|raise NotImplemented|TODO|FIXME|placeholder|stub", re.I)


def line_of(path: Path | str, needle: str) -> int | None:
    text = path if isinstance(path, str) else path.read_text(encoding="utf-8")
    try:
        for number, line in enumerate(text.splitlines(), start=1):
            if needle in line:
                return number
    except (OSError, TypeError):
        return None
    return None


def check_impl() -> dict:
    text = IMPL.read_text(encoding="utf-8")
    hits = [(n, line.strip()[:110]) for n, line in enumerate(text.splitlines(), 1) if PLACEHOLDER.search(line)]
    return {
        "path": str(IMPL.relative_to(SUBJECT)),
        "lines": text.count("\n") + 1,
        "bytes": IMPL.stat().st_size,
        "placeholder_hits": hits,
        "entry_main_line": line_of(IMPL, "def main(argv"),
        "run_functions": {
            name: line_of(IMPL, f"def run_{name}(")
            for name in ("bem_added_mass", "fem_plate", "lbm_channel")
        },
        "contract_line": line_of(IMPL, "FRAME_CONTRACT_ID ="),
        "solver_api_calls": {
            "tensorbem.compute_added_mass_matrix_3d": line_of(IMPL, "compute_added_mass_matrix_3d"),
            "tensorbem.TriangleMesh": line_of(IMPL, "bem.TriangleMesh("),
            "tensorfem.plate.q4_mindlin_stiffness": line_of(IMPL, "plate_mod.q4_mindlin_stiffness"),
            "tensorfem.plate.kirchhoff_sine_center_deflection": line_of(
                IMPL, "kirchhoff_sine_center_deflection"
            ),
            "tensorlbm.collide_bgk": line_of(IMPL, "lbm.collide_bgk"),
            "tensorlbm.macroscopic": line_of(IMPL, "lbm.macroscopic"),
        },
    }


def check_entry() -> dict:
    text = INIT.read_text(encoding="utf-8")
    return {
        "module_entry": f"tensorcad/solver_interface.py:{line_of(IMPL, 'def main(argv')} "
        "-> PYTHONPATH=. python -m tensorcad.solver_interface --out <dir>",
        "package_export_line": line_of(INIT, "from .solver_interface import"),
        "export_names": sorted(
            name.strip()
            for name in re.findall(
                r"from \.solver_interface import \(([^)]*)\)", text, flags=re.S
            )[0].replace("\n", " ").split(",")
            if name.strip()
        ),
    }


def check_contract() -> dict:
    """The z-up contract must be enforced, not merely declared."""
    sys.path.insert(0, str(SUBJECT))
    import importlib  # noqa: PLC0415

    mod = importlib.import_module("tensorcad.solver_interface")
    from tensorcad.mesh import Mesh  # noqa: PLC0415

    good = Mesh(
        vertices=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
        faces=[(0, 2, 1), (0, 1, 3), (1, 2, 3), (0, 3, 2)],
        coordinate_frame=mod.CONTRACT_FRAME,
    )
    out = {"contract_id": mod.FRAME_CONTRACT_ID, "contract_frame": mod.CONTRACT_FRAME}
    try:
        mod.frame_of(
            Mesh(vertices=good.vertices, faces=good.faces, coordinate_frame="y-up")
        )
        out["refuses_wrong_frame"] = False
    except mod.FrameContractError:
        out["refuses_wrong_frame"] = True
    try:
        mod.frame_of(
            Mesh(vertices=good.vertices, faces=good.faces,
                 coordinate_frame=mod.CONTRACT_FRAME, units="mm")
        )
        out["refuses_wrong_units"] = False
    except mod.FrameContractError:
        out["refuses_wrong_units"] = True
    try:
        out["accepts_contract_frame"] = mod.frame_of(good) == mod.CONTRACT_FRAME
    except Exception as exc:  # pragma: no cover - reported, not hidden
        out["accepts_contract_frame"] = f"{type(exc).__name__}: {exc}"
    return out


def _rel(a: float, b: float) -> float:
    return abs(a - b) / abs(b) if b else abs(a - b)


def check_artifacts() -> dict:
    payload = json.loads((ARTIFACTS / "results.json").read_text(encoding="utf-8"))
    rows, gate_fail, recomputed_bad = [], [], []
    for run in payload["runs"]:
        for cmp in run.get("comparisons", []):
            measured, reference = cmp.get("measured"), cmp.get("reference")
            row = {
                "solver": run["solver"],
                "name": cmp["name"],
                "measured": measured,
                "reference": reference,
                "rel_error": cmp.get("rel_error"),
                "tolerance": cmp.get("tolerance"),
                "passed": cmp.get("passed"),
            }
            rows.append(row)
            if cmp.get("passed") is False:
                gate_fail.append(row)
            # independent arithmetic re-check of the claimed relative error (Unity: |meas/ref-1|)
            if measured is not None and reference not in (None, 0):
                claimed = cmp.get("rel_error")
                actual = measured / reference - 1.0
                if claimed is not None and _rel(claimed, actual) > 1e-9:
                    recomputed_bad.append(
                        {"name": cmp["name"], "claimed": claimed, "recomputed": actual}
                    )
    ran = [run["solver"] for run in payload["runs"] if run.get("status") == "ran"]
    files = sorted(str(p.relative_to(ARTIFACTS)) for p in ARTIFACTS.rglob("*") if p.is_file())
    return {
        "statuses": {run["solver"]: run.get("status") for run in payload["runs"]},
        "ran": ran,
        "gates_failed": payload.get("gates_failed"),
        "comparisons": len(rows),
        "gated_rows": sum(1 for r in rows if r["tolerance"] is not None),
        "failed_gated_rows": gate_fail,
        "rel_error_recomputed_mismatches": recomputed_bad,
        "artifact_files": files,
        "solver_files_recorded": {
            run["solver"]: run.get("quantities", {}).get("solver_file") for run in payload["runs"]
        },
        "solver_versions_recorded": {
            run["solver"]: run.get("quantities", {}).get("solver_version") for run in payload["runs"]
        },
        "environment": payload.get("environment", {}),
    }


def check_callers() -> dict:
    bench = BENCH.read_text(encoding="utf-8")
    init = INIT.read_text(encoding="utf-8")
    platform = subprocess.run(
        ["grep", "-rl", "solver_interface", "/root/tensor-ecosystem/TensorPlatform/backend"],
        capture_output=True, text=True, check=False,
    ).stdout.split()
    return {
        "in_repo_benchmark_case": line_of(BENCH, "def solver_handoff_case"),
        "benchmark_calls": {
            name: line_of(BENCH, f"solver_interface.run_{name}") for name in
            ("bem_added_mass", "fem_plate", "lbm_channel")
        },
        "package_export": line_of(INIT, "from .solver_interface import"),
        "platform_side_callers": platform,
    }


def run_pytest() -> dict:
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "tests/test_solver_interface.py", "-q"],
        cwd=SUBJECT, capture_output=True, text=True, check=False,
        env={"PYTHONPATH": ".", "PATH": "/usr/bin:/bin"},
    )
    tail = [line for line in proc.stdout.strip().splitlines() if line.strip()][-2:]
    return {"exit": proc.returncode, "tail": tail}


def run_repro(scratch: Path) -> dict:
    proc = subprocess.run(
        [sys.executable, "-m", "tensorcad.solver_interface", "--out", str(scratch)],
        cwd=SUBJECT, capture_output=True, text=True, check=False,
        env={"PYTHONPATH": ".", "PATH": "/usr/bin:/bin"},
    )
    fresh_path = scratch / "results.json"
    if not fresh_path.exists():
        return {"exit": proc.returncode, "error": "no results.json written", "stderr": proc.stderr[-800:]}
    shipped = json.loads((ARTIFACTS / "results.json").read_text(encoding="utf-8"))
    fresh = json.loads(fresh_path.read_text(encoding="utf-8"))

    def table(payload):
        out = {}
        for run in payload["runs"]:
            for cmp in run.get("comparisons", []):
                out[(run["solver"], cmp["name"])] = cmp
        return out

    a, b = table(shipped), table(fresh)
    worst, flips = 0.0, []
    for key in set(a) & set(b):
        x, y = a[key].get("measured"), b[key].get("measured")
        if x is not None and y is not None:
            worst = max(worst, _rel(y, x))
        if a[key].get("passed") != b[key].get("passed"):
            flips.append(key)
    return {
        "exit": proc.returncode,
        "wall_seconds": None,
        "comparisons_shipped": len(a),
        "comparisons_fresh": len(b),
        "keys_only_in_shipped": sorted(":".join(k) for k in set(a) - set(b)),
        "keys_only_in_fresh": sorted(":".join(k) for k in set(b) - set(a)),
        "max_relative_difference": worst,
        "gate_flips": flips,
        "fresh_gates_failed": fresh.get("gates_failed"),
        "summary_tail": [line for line in proc.stdout.strip().splitlines() if line.strip()][-1:],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pytest", action="store_true")
    parser.add_argument("--repro", action="store_true")
    parser.add_argument("--scratch", default="/tmp/m17576_audit_repro")
    args = parser.parse_args()
    if not SUBJECT.exists():
        print(json.dumps({"error": f"subject worktree missing: {SUBJECT}"}))
        return 2
    out = {
        "anchor": "GAP-CAD-B3-SOLVER-INTERFACE",
        "audited_mission": 17476,
        "audit_mission": 17576,
        "subject": str(SUBJECT),
        "impl": check_impl(),
        "entry": check_entry(),
        "contract": check_contract(),
        "artifacts": check_artifacts(),
        "callers": check_callers(),
        "docs_present": DOCS.exists(),
    }
    if args.pytest:
        out["pytest"] = run_pytest()
    if args.repro:
        out["repro"] = run_repro(Path(args.scratch))
    print(json.dumps(out, indent=2, sort_keys=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())