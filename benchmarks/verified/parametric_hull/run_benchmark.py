"""Benchmark: parametric hull geometry and hydrostatics against closed-form references.

Runs the verification campaign (``examples/parametric_hull_dag/run_verification.py``)
and, on top of it, the end-to-end command-line handoff
(``python -m tensorcad.parametric_hull_cli wigley --gmsh ...``), then writes
``result.json``: one row per benchmarked quantity with the *reference* value, the
*measured* value, the error and the acceptance threshold. That table is what a reader
(or a downstream gate) can compare without re-running anything.

    python benchmarks/verified/parametric_hull/run_benchmark.py            # full run
    python benchmarks/verified/parametric_hull/run_benchmark.py --reuse     # reuse record

References, and why they are references and not self-consistency:

* Wigley III / generalised Wigley: the closed-form hydrostatics
  ``V = p q L B T / ((p+1)(q+1))``, ``Aw = p L B / (p+1)``, ``I_T``, ``I_L``, ``VCB``.
  These are the textbook values (``V = 4/9 L B T``, ``Aw = 2/3 L B``, ``VCB = -0.375 T``
  for ``p = q = 2``); ``tests/test_parametric_hull_dag.py`` pins the closed form to
  those numbers independently of the mesh code.
* DARPA SUBOFF: the four-segment radius table of Groves, Huang & Chang (1989),
  transcribed in :mod:`tensorcad.suboff` and pinned by ``tests/test_suboff.py``. The
  hull is measured against that table, not against a model of it.
* Mesh quality: the scaled-Jacobian / aspect-ratio definitions in
  :func:`tensorcad.parametric_hull.mesh_quality_report`, thresholds from the node's
  acceptance criteria (min scaled Jacobian > 0.3, max aspect ratio < 5).

Scope limit, stated up front: no experimental hull-form data is used, because none of
the three quantities benchmarked here (geometry, displacement, waterplane) is an
experimental observable for these bodies. The external anchor is the published SUBOFF
offset table; the rest is closed-form mathematics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:  # the checkout, not an installed wheel, is what runs
    sys.path.insert(0, str(ROOT))
HERE = pathlib.Path(__file__).resolve().parent
RECORD = HERE / "verification.json"
RESULT = HERE / "result.json"
DRIVER = ROOT / "examples" / "parametric_hull_dag" / "run_verification.py"

REPRO_CMD = "python benchmarks/verified/parametric_hull/run_benchmark.py"


def _sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def _git_commit() -> str:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True
        )
        return out.stdout.strip()
    except Exception:
        return "unknown"


def run_driver(*, figures: bool) -> dict:
    """Run the verification campaign and return its record."""
    cmd = [sys.executable, str(DRIVER), "--out", str(HERE)]
    if not figures:
        cmd.append("--no-figures")
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(
            f"verification driver failed (rc={proc.returncode})\n"
            f"{proc.stdout[-4000:]}\n{proc.stderr[-4000:]}"
        )
    return json.loads(RECORD.read_text())


def run_cli_handoff() -> dict:
    """End-to-end: CLI -> gmsh 2.2 file, measured back from the exported file."""
    msh = HERE / "wigley_L10.msh"
    module = "tensorcad.parametric_hull_cli"
    proc = subprocess.run(
        [sys.executable, "-m", module, "wigley", "--length", "10", "--beam", "1",
         "--draft", "0.625", "--stations", "11", "--points", "12",
         "--mesh-n-u", "97", "--mesh-n-v", "16", "--gmsh", str(msh),
         "--json", str(HERE / "wigley_L10.json")],
        cwd=ROOT, capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise SystemExit(
            f"CLI handoff failed (rc={proc.returncode})\n"
            f"{proc.stdout[-3000:]}\n{proc.stderr[-3000:]}"
        )
    record = json.loads(proc.stdout)
    # The exported file is the evidence: re-read it instead of trusting the record.
    from tensorcad.adapters.fem import read_gmsh

    parsed = read_gmsh(msh)
    return {
        "repro": "python -m tensorcad.parametric_hull_cli wigley ... --gmsh wigley_L10.msh",
        "mesh_file": str(msh.relative_to(ROOT)),
        "mesh_sha256": _sha256(msh),
        "gmsh_lines": sum(1 for _ in msh.open()),
        "nodes_in_file": int(parsed["nodes"]),
        "triangles_in_file": int(parsed["triangles"]),
        "watertight": record["gmsh"]["watertight"],
        "edge_manifold": record["gmsh"]["edge_manifold"],
        "euler_characteristic": record["gmsh"]["euler_characteristic"],
        "volume_from_file_m3": record["gmsh"]["volume_m3"],
        "measured_volume_m3": record["hydrostatics"]["displaced_volume_m3"],
        "analytic_volume_m3": record["analytic_reference"]["displaced_volume"],
        "gate_failures": record["gate_failures"],
    }


def build_cases(record: dict, cli: dict) -> list[dict]:
    """One row per benchmarked quantity: reference, measured, error, threshold."""
    a, b, c, d, e, f = (
        record["criterion_a"], record["criterion_b"], record["criterion_c"],
        record["criterion_d"], record["criterion_e"], record["criterion_f"],
    )
    cases = [
        {
            "id": "a_wigley_ordinates",
            "quantity": (
                f"max relative error of hull ordinates y(x,z) at {a['n_stations']} stations "
                f"x {a['n_waterlines']} waterlines"
            ),
            "reference": "y = (B/2)(1-(2x-L)^2/L^2)(1-(z/T)^2), Wigley III closed form",
            "measured": a["value"],
            "error": a["value"],          # the metric *is* the relative error
            "threshold": a["threshold"],
            "passed": bool(record["verdicts"]["a"]),
        },
        {
            "id": "b_wigley_displacement",
            "quantity": "displaced volume [m^3]",
            "reference": f"V = p q L B T/((p+1)(q+1)) = {b['analytic_volume']:.6f}",
            "reference_value": b["analytic_volume"],
            "measured": b["measured_volume"],
            "error": (b["measured_volume"] - b["analytic_volume"]) / b["analytic_volume"],
            "error_measure": "relative",
            "threshold": b["threshold"],
            "passed": bool(record["verdicts"]["b"]),
        },
        {
            "id": "c_wigley_waterplane",
            "quantity": "waterplane area [m^2]",
            "reference": f"Aw = p L B/(p+1) = {c['analytic_waterplane_area']:.6f}",
            "reference_value": c["analytic_waterplane_area"],
            "measured": c["measured_waterplane_area"],
            "error": (c["measured_waterplane_area"] - c["analytic_waterplane_area"])
                     / c["analytic_waterplane_area"],
            "error_measure": "relative",
            "threshold": c["threshold"],
            "passed": bool(record["verdicts"]["c"]),
        },
        {
            "id": "d_suboff_offsets",
            "quantity": f"max relative radius error over {d['n_stations']} SUBOFF stations",
            "reference": "DARPA SUBOFF four-segment radius table (Groves, Huang & Chang 1989)",
            "measured": d["value"],
            "error": d["value"],
            "threshold": d["threshold"],
            "passed": bool(record["verdicts"]["d"]),
        },
        {
            "id": "e_min_scaled_jacobian",
            "quantity": "minimum scaled Jacobian over the hull triangles",
            "reference": "criterion floor 0.30 (equilateral = 0.866)",
            "measured": e["value"]["min_scaled_jacobian"],
            "error": e["value"]["min_scaled_jacobian"] - e["threshold"]["min_scaled_jacobian"],
            "error_measure": "absolute margin above the floor",
            "threshold": e["threshold"]["min_scaled_jacobian"],
            "passed": bool(e["passed"]),
        },
        {
            "id": "e_max_aspect_ratio",
            "quantity": "maximum triangle aspect ratio",
            "reference": "criterion ceiling 5.0",
            "measured": e["value"]["max_aspect_ratio"],
            "error": e["threshold"]["max_aspect_ratio"] - e["value"]["max_aspect_ratio"],
            "error_measure": "absolute margin below the ceiling",
            "threshold": e["threshold"]["max_aspect_ratio"],
            "passed": bool(e["passed"]),
        },
        {
            "id": "f_parametric_cb_sweep",
            "quantity": "max relative volume error over the Cb sweep (Cb 0.35 -> 0.55)",
            "reference": "closed-form volume at each Cb",
            "measured": f["value"],
            "error": f["value"],
            "threshold": f["threshold"],
            "passed": bool(record["verdicts"]["f"]),
        },
        {
            "id": "cli_gmsh_volume",
            "quantity": "displaced volume of the exported gmsh file [m^3]",
            "reference": (
                f"V = 4/9 L B T = {cli['analytic_volume_m3']:.6f} for L=10, B=1, T=0.625"
            ),
            "reference_value": cli["analytic_volume_m3"],
            "measured": cli["volume_from_file_m3"],
            "error": (cli["volume_from_file_m3"] - cli["analytic_volume_m3"])
                     / cli["analytic_volume_m3"],
            "error_measure": "relative",
            "threshold": 5e-3,
            "passed": bool(cli["watertight"] and cli["edge_manifold"]),
        },
    ]
    return cases


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="parametric hull benchmark (ref vs got)")
    parser.add_argument("--reuse", action="store_true", help="reuse an existing verification.json")
    parser.add_argument("--no-figures", action="store_true", help="skip figure regeneration")
    args = parser.parse_args(argv)

    if args.reuse and RECORD.is_file():
        record = json.loads(RECORD.read_text())
    else:
        record = run_driver(figures=not args.no_figures)
    cli = run_cli_handoff()
    cases = build_cases(record, cli)
    payload = {
        "benchmark": "parametric_hull_geometry_and_hydrostatics",
        "node": record["node"],
        "repository": record["repository"],
        "repro_cmd": REPRO_CMD,
        "commit": _git_commit(),
        "reference_hull": record["reference_hull"],
        "mesh": record["mesh"],
        "cases": cases,
        "all_passed": all(case["passed"] for case in cases),
        "cli_handoff": cli,
        "convergence": {
            "observed_order": record["convergence"]["observed_order"],
            "richardson_extrapolated_rel_error": record["convergence"][
                "richardson_extrapolated_rel_error"],
        },
        "artifacts": {
            "verification_record": str(RECORD.relative_to(ROOT)),
            "verification_sha256": _sha256(RECORD),
            "figures": [
                {"name": figure.get("name"), "path": figure.get("path")}
                for figure in record.get("figures", [])
            ],
        },
    }
    RESULT.write_text(json.dumps(payload, indent=2) + "\n")

    width = max(len(case["id"]) for case in cases)
    print(f"{'case'.ljust(width)}  {'measured':>12}  {'error':>12}  {'threshold':>10}  verdict")
    for case in cases:
        print(
            f"{case['id'].ljust(width)}  {case['measured']:12.6g}  {case['error']:12.4e}  "
            f"{case['threshold']:10.4g}  {'PASS' if case['passed'] else 'FAIL'}"
        )
    print(f"\nall_passed = {payload['all_passed']}   record -> {RESULT.relative_to(ROOT)}")
    return 0 if payload["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
