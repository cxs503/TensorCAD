#!/usr/bin/env python3
"""Resolution-controlled, machine-gated hull export, verified end to end.

The example is a *measurement*, not a demo: it generates the parametric hull on a
ladder of grids, compares every grid against the analytic displacement volume
(an integral the mesh never sees), shows the machine gate rejecting the coarsest
one, refines a coarse request until the analytic error target is met, writes a
binary STL, reads the STL back, and re-runs the gate on the recovered triangle
soup. It exits non-zero if any of those claims does not hold.

Usage::

    python examples/hull_stl_export/run_example.py [--output results/hull-export]

Outputs (in ``--output``):
    resolution_study.json   per-grid volume, analytic error and facet count
    coarse_gate.json        the gate report for an under-resolved export
    resolved/               hull.obj, hull.stl, offsets.csv, metadata.json
    summary.json            the resolution plan and the STL read-back evidence
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:  # allow running from a checkout without installing
    sys.path.insert(0, str(REPO_ROOT))

from tensorcad import (
    HullParameters,
    analytic_volume,
    export_hull,
    gate_mesh,
    generate_hull,
    read_stl,
    resolve_resolution,
)
from tensorcad import mesh as meshlib

LADDER: tuple[tuple[int, int], ...] = ((11, 5), (21, 9), (41, 13), (81, 21), (161, 41))
TARGET_VOLUME_ERROR = 0.002


def resolution_study(length: float, beam: float, draft: float) -> list[dict[str, object]]:
    exact = analytic_volume(length, beam, draft)
    rows: list[dict[str, object]] = []
    for nx, nz in LADDER:
        model = generate_hull(HullParameters(length, beam, draft, nx, nz))
        volume = model.hydrostatics.displacement_volume if model.hydrostatics else float("nan")
        rows.append(
            {
                "nx": nx,
                "nz": nz,
                "triangles": len(model.mesh.faces),
                "vertices": len(model.mesh.vertices),
                "volume_m3": volume,
                "analytic_volume_m3": exact,
                "relative_error": volume / exact - 1.0,
            }
        )
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="results/hull-export")
    args = parser.parse_args(argv)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    length, beam, draft = 10.0, 2.0, 0.8
    exact = analytic_volume(length, beam, draft)

    print("## Resolution ladder against the analytic displacement volume")
    rows = resolution_study(length, beam, draft)
    (output / "resolution_study.json").write_text(json.dumps(rows, indent=2) + "\n")
    print("| nx | nz | triangles | volume [m^3] | error vs analytic |")
    print("|---:|---:|---:|---:|---:|")
    for row in rows:
        print(
            f"| {row['nx']} | {row['nz']} | {row['triangles']} | "
            f"{float(row['volume_m3']):.6f} | {100.0 * float(row['relative_error']):+.3f} % |"
        )

    # 1. The gate must reject the coarsest grid, even though it is closed and manifold.
    coarse = generate_hull(HullParameters(length, beam, draft, 11, 5))
    coarse_report = gate_mesh(coarse.mesh, reference_volume=exact, volume_tolerance=0.02)
    (output / "coarse_gate.json").write_text(json.dumps(coarse_report.to_dict(), indent=2) + "\n")
    assert not coarse_report.passed, "the gate accepted an under-resolved mesh"
    assert "reference_volume" in {check.name for check in coarse_report.failures}
    print(
        "\ncoarse 11x5 gate: FAILED as required -> "
        f"{[check.name for check in coarse_report.failures]} "
        f"(volume error {100.0 * coarse_report.measurements['relative_volume_error']:+.3f} %)"
    )

    # 2. Refine until the analytic target is met, then export and verify.
    plan = resolve_resolution(HullParameters(length, beam, draft, 11, 5),
                              target_volume_error=TARGET_VOLUME_ERROR)
    assert plan.accepted, plan.to_dict()
    artifacts = export_hull(
        coarse, output / "resolved", stl=True, target_volume_error=TARGET_VOLUME_ERROR
    )
    metadata = json.loads(Path(artifacts["metadata"]).read_text())
    resolved = metadata["resolution"]
    stl_record = metadata["stl"]
    recovered = read_stl(artifacts["stl"])
    recovered_volume = meshlib.total_signed_volume(recovered)

    summary = {
        "analytic_volume_m3": exact,
        "target_volume_error": TARGET_VOLUME_ERROR,
        "plan": plan.to_dict(),
        "metadata_resolution": resolved,
        "resolved_mesh_relative_error": metadata["gate"]["measurements"]["relative_volume_error"],
        "stl": {
            "path": artifacts["stl"],
            "triangles": len(recovered.faces),
            "vertices": len(recovered.vertices),
            "watertight": meshlib.is_watertight(recovered),
            "edge_manifold": meshlib.is_edge_manifold(recovered),
            "signed_volume_m3": recovered_volume,
            "relative_volume_error_vs_analytic": recovered_volume / exact - 1.0,
            "read_back_gate_passed": stl_record["gate"]["passed"],
        },
        "artifacts": artifacts,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(
        f"\nresolved grid {resolved['parameters']['nx']}x{resolved['parameters']['nz']} "
        f"({resolved['triangles']} triangles): "
        f"{100.0 * float(resolved['relative_error']):+.4f} % vs analytic "
        f"(target {100.0 * TARGET_VOLUME_ERROR:g} %)"
    )
    print(
        f"STL read back: {summary['stl']['triangles']} triangles, "
        f"watertight={summary['stl']['watertight']}, "
        f"euler={len(recovered.vertices) - len(meshlib.edge_face_counts(recovered)) + len(recovered.faces)}, "
        f"volume error vs analytic {100.0 * float(summary['stl']['relative_volume_error_vs_analytic']):+.5f} %"
    )

    checks = [
        ("|plan error| <= target", abs(plan.relative_error) <= TARGET_VOLUME_ERROR),
        ("coarse gate failed on the analytic volume",
         "reference_volume" in {check.name for check in coarse_report.failures}),
        ("resolved gate passed", metadata["gate"]["passed"]),
        ("STL read back watertight", bool(summary["stl"]["watertight"])),
        ("STL read back edge-manifold", bool(summary["stl"]["edge_manifold"])),
        ("STL volume within target of analytic",
         abs(float(summary["stl"]["relative_volume_error_vs_analytic"])) <= TARGET_VOLUME_ERROR),
    ]
    failed = [name for name, ok in checks if not ok]
    for name, ok in checks:
        print(f"- {'PASS' if ok else 'FAIL'} {name}")
    print(f"\nwrote {output}/resolution_study.json, summary.json and resolved/ artifacts")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
