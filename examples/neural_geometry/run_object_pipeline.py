#!/usr/bin/env python3
"""Run the neural geometry object pipeline and *verify* it numerically.

This is the end-to-end benchmark for the 0.3 neural geometry layer: it fits an
object, edits it, derives a mesh, re-opens it, and deliberately attempts an edit
that must be rejected. Every claim it prints is checked, and the script exits
non-zero if any of them stops holding.

What it proves, in order:

1. a field can be persisted, reloaded and queried (round trip through the bundle);
2. a local edit is committed only when it meets measured gates, and the result is
   a new revision whose parent hash matches what was on disk;
3. the derived mesh is watertight, edge-manifold and genus 0 — the property a
   solver actually needs;
4. a failing edit leaves the active revision byte-identical and appends the
   failure (with finite-safe metrics) to the history.

Usage::

    python examples/neural_geometry/run_object_pipeline.py --output results/object
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:  # allow running from a checkout without installing
    sys.path.insert(0, str(ROOT))

from tensorcad.neural.objects import (
    Tolerances,
    apply_dent_edit,
    extract_object,
    fit_sphere,
)
from tensorcad.neural_geometry import ObjectStore, read_manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", default="results/object", help="directory for the object and reports")
    parser.add_argument("--radius", type=float, default=0.5)
    parser.add_argument("--fit-steps", type=int, default=1200)
    parser.add_argument("--edit-steps", type=int, default=1000)
    parser.add_argument("--resolution", type=int, default=96)
    parser.add_argument("--seed", type=int, default=503)
    parser.add_argument("--reject-depth", type=float, default=0.25,
                        help="dent depth that must fail validation")
    args = parser.parse_args(argv)

    failures: list[str] = []

    def expect(label: str, condition: bool, detail: str = "") -> None:
        print(f"  [{'PASS' if condition else 'FAIL'}] {label}{': ' + detail if detail else ''}")
        if not condition:
            failures.append(label)

    output = Path(args.output)
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)
    target = output / "sphere.ngo"

    print("1. fit revision 1")
    started = time.perf_counter()
    fit = fit_sphere(
        target,
        radius=args.radius,
        steps=args.fit_steps,
        seed=args.seed,
        tolerance=Tolerances(),
        command=["run_object_pipeline", "--fit-steps", str(args.fit_steps)],
    )
    fit_seconds = time.perf_counter() - started
    store = ObjectStore(target)
    manifest = store.manifest()
    print(
        f"     field MAE {fit.metrics['field_mae_m']:.6f} m, "
        f"surface error mean {fit.metrics['boundary_error_mean_m']:.6f} m, "
        f"coverage {fit.metrics['ray_coverage']:.3f}"
    )
    expect("revision 1 published", manifest.revision == 1)
    expect("parent is absent on a root revision", manifest.lineage.parent_revision is None)
    expect("gates passed", fit.passed)
    expect("weights hash recorded", manifest.weights.sha256 == fit.weights_sha256)
    expect("history has one committed attempt", [e["status"] for e in store.history()] == ["committed"])
    parent_hash = store.active_weights_sha256()

    print("2. local dent edit (transaction)")
    edit = apply_dent_edit(
        target,
        depth=0.06,
        sigma=0.22,
        steps=args.edit_steps,
        seed=args.seed,
        extract_resolution=args.resolution,
    )
    metrics = edit.metrics
    print(
        f"     dent peak {metrics['dent_peak_measured_m']:.5f} m (target "
        f"{metrics['dent_peak_target_m']:.5f}), region error mean "
        f"{metrics['dent_region_error_mean_m']:.5f} m"
    )
    print(
        f"     protected mean {metrics['protected_mean_abs_displacement_m']:.5f} m, "
        f"P95 {metrics['protected_p95_abs_displacement_m']:.5f} m, "
        f"max {metrics['protected_max_abs_displacement_m']:.5f} m"
    )
    expect("edit committed", edit.committed)
    expect("new revision is 2", edit.active_revision == 2)
    current = store.manifest()
    expect("lineage points at revision 1", current.lineage.parent_revision == 1 and current.lineage.root_revision == 1)
    expect(
        "training record names the parent weights",
        current.raw["training"]["parent_weights_sha256"] == parent_hash,
    )
    expect("dent reached at least 80% of its target depth",
           metrics["dent_peak_measured_m"] >= 0.8 * 0.06)
    expect("protected region stayed inside its gates", not edit.failed_checks)

    print("3. derived mesh")
    _mesh, report, records = extract_object(target, resolution=args.resolution, formats=("obj", "stl"))
    print(
        f"     {report.vertices} vertices, {report.triangles} triangles, "
        f"volume {report.signed_volume:.6f} m^3, area {report.surface_area:.6f} m^2"
    )
    expect("derived mesh is watertight", report.watertight)
    expect("derived mesh is edge-manifold", report.edge_manifold)
    expect("derived mesh is genus zero", report.euler_characteristic == 2)
    expect("derived mesh is positively oriented", report.signed_volume > 0.0)
    expect("every derived file exists", all((target / record["path"]).is_file() for record in records))
    expect(
        "committed manifest carries the derived artifacts",
        {artifact["format"] for artifact in store.manifest().derived} == {"obj", "stl"},
    )
    analytic_sphere = 4.0 / 3.0 * 3.141592653589793 * args.radius**3
    expect(
        "dented volume is below the sphere it started from",
        report.signed_volume < analytic_sphere,
        f"{report.signed_volume:.6f} < {analytic_sphere:.6f}",
    )

    print("4. failure injection: an edit that must be rejected")
    before_weights = store.active_weights_sha256()
    before_manifest = (store.root / "manifest.json").read_bytes()
    rejected = apply_dent_edit(
        target, depth=args.reject_depth, sigma=0.15, steps=args.edit_steps, seed=args.seed
    )
    print(f"     status {rejected.status}; failed gates: {', '.join(rejected.failed_checks)}")
    expect("rejected edit did not commit", not rejected.committed)
    expect("rejected edit names the gates it broke", bool(rejected.failed_checks))
    expect("active revision unchanged", store.active_revision() == 2)
    expect("weights bytes unchanged", store.active_weights_sha256() == before_weights)
    expect("manifest bytes unchanged", (store.root / "manifest.json").read_bytes() == before_manifest)
    history = store.history()
    expect("failure is recorded in history", history[-1]["status"] == "rolled_back")
    expect(
        "history keeps the metrics that justified the rejection",
        isinstance(history[-1].get("metrics"), dict) and bool(history[-1].get("failed_checks")),
    )
    expect(
        "history is valid standard JSON",
        _is_standard_json((store.root / "history.json").read_text(encoding="utf-8")),
    )
    expect(
        "re-opened object still verifies",
        read_manifest(store.root).revision == 2,
    )

    summary = {
        "object": str(target),
        "revision": store.active_revision(),
        "fit_seconds": fit_seconds,
        "edit_seconds": edit.elapsed_seconds,
        "metrics": {
            "fit": fit.metrics,
            "edit": metrics,
            "mesh": report.to_dict(),
        },
        "history": [
            {"revision": entry.get("revision"), "status": entry.get("status"),
             "operation": entry.get("operation"), "failed_checks": entry.get("failed_checks", [])}
            for entry in history
        ],
        "checks_failed": failures,
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print()
    if failures:
        print(f"FAILED {len(failures)} check(s): {', '.join(failures)}")
        return 1
    print("all checks passed")
    return 0


def _is_standard_json(text: str) -> bool:
    try:
        json.loads(text, parse_constant=lambda token: (_ for _ in ()).throw(ValueError(token)))
    except ValueError:
        return False
    return True


if __name__ == "__main__":
    raise SystemExit(main())
