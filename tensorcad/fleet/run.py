"""Run the baseline fleet through the intake gate and write reproducible artefacts.

``run_fleet`` is the single call behind both the library entry point and the CLI. It
does four things and records all of them:

1. build every registered hull through its declared geometry route;
2. run the full gate set on each mesh, keeping the number behind every verdict;
3. cross-check the two independent geometry routes (form-family loft vs lofted offset
   table) on one ship, so a gate that only ever sees one route cannot hide a bug;
4. run the negative controls, which must be *rejected* -- a gate that passes everything
   is indistinguishable from no gate at all.

Artefacts land in ``<out_dir>``: ``evidence.json`` (everything), ``ref_vs_got.csv``
(the benchmark table), ``manifests/<ship>.json`` (per-ship provenance manifest) and
``settings.json`` (so the run can be repeated bit-for-bit).
"""

from __future__ import annotations

import csv
import json
import math
import platform
import sys
from pathlib import Path
from typing import Any

from tensorcad import mesh as meshlib
from tensorcad.fleet.build import (
    RouteRefusal,
    build_ship_geometry,
    loft_station_table,
    sample_offset_table,
)
from tensorcad.fleet.controls import (
    clipped_cb_record,
    duplicate_manifest,
    duplicate_stl_pair,
    placeholder_box,
    wrong_volume_record,
)
from tensorcad.fleet.gates import (
    VOLUME_TOLERANCE,
    evaluate_gates,
    gate_particulars_consistent,
    measure,
    placeholder_reasons,
    refused,
)
from tensorcad.fleet.manifest import FleetManifest, fleet_sha_index, load_manifest, sha256_file

SCHEMA = "tensorcad.fleet.run.v1"


def _hydrostatics(mesh) -> dict[str, Any]:
    """Independent volume/inertia pass through the hydrostatics solver.

    This is the platform-side caller: the geometry the gate admitted is handed to
    :func:`tensorcad.hydrostatics.compute_hydrostatics` rather than to the same code
    path that built it, so a meshing error that fooled the builder is caught here.
    """
    from tensorcad.hydrostatics import compute_hydrostatics

    result = compute_hydrostatics(mesh, waterline_z=0.0, water_density=1025.0)
    return {
        "displacement_volume": result.displacement_volume,
        "displacement_mass": result.displacement_mass,
        "keel_z": result.keel_z,
        "lcb": result.longitudinal_centre_buoyancy,
        "vcb": result.vertical_centre_buoyancy,
        "block_coefficient": result.block_coefficient,
        "waterplane_area": result.waterplane_area,
        "wetted_surface_area": result.wetted_surface_area,
        "watertight": meshlib.is_watertight(mesh),
    }


def _convergence(record, *, nx: int, nz: int) -> dict[str, Any] | None:
    """Grid convergence of the enclosed volume on three successively halved grids.

    The geometry routes interpolate the section shape linearly between waterlines, so the
    volume has a first-order deficit that must be shown to vanish. Three grids give the
    observed order ``p`` and a Richardson estimate of the exact volume; the gate then
    reports both the fine-grid volume (what the artefact actually is) and the extrapolated
    one (what the law meant), so a discretisation deficit can never be tuned away by
    choosing a friendly tolerance.
    """
    grids = [(max(31, nx // 4), max(9, nz // 4)), (max(41, nx // 2), max(13, nz // 2)), (nx, nz)]
    volumes: list[float] = []
    for gx, gz in grids:
        try:
            built = build_ship_geometry(record, nx=gx, nz=gz)
        except RouteRefusal as exc:
            return {"refused": str(exc), "grids": [list(g) for g in grids], "volumes": volumes}
        volumes.append(meshlib.total_signed_volume(built.mesh))
    reference = record.reference_volume()
    errors = [100.0 * (volume - reference) / reference for volume in volumes]
    if len(volumes) == 3:
        delta_fine = volumes[2] - volumes[1]
        delta_coarse = volumes[1] - volumes[0]
        ratio = abs(delta_coarse) / abs(delta_fine) if delta_fine else float("nan")
        order = math.log2(ratio) if ratio and ratio > 0 else float("nan")
        extrapolated = volumes[2] + delta_fine / (2.0 ** order - 1.0) if order == order and order > 0.1 else float("nan")
        gci = 1.25 * abs(delta_fine) / abs(volumes[2]) / (2.0 ** order - 1.0) if order == order and order > 0.1 else float("nan")
    else:
        ratio = order = extrapolated = gci = float("nan")
    return {
        "grids": [list(g) for g in grids],
        "volumes": volumes,
        "relative_errors_percent": errors,
        "observed_order": order,
        "error_ratio_coarse_over_fine": ratio,
        "richardson_volume": extrapolated,
        "richardson_relative_error_percent": 100.0 * (extrapolated - reference) / reference
        if extrapolated == extrapolated else None,
        "gci_fine": gci if gci == gci else None,
        "monotone": volumes[0] <= volumes[1] <= volumes[2],
    }


def _ship_block(record, manifest: FleetManifest, *, nx: int, nz: int,
                volume_tolerance: float, convergence: bool = True) -> dict[str, Any]:
    try:
        built = build_ship_geometry(record, nx=nx, nz=nz)
    except RouteRefusal as exc:
        block = refused(record, str(exc))
        block["record"] = record.as_dict()
        return block
    verdict = evaluate_gates(record, built.mesh, manifest, volume_tolerance=volume_tolerance)
    verdict["record"] = record.as_dict()
    verdict["build"] = built.as_dict()
    verdict["hydrostatics"] = _hydrostatics(built.mesh)
    # the hydrostatics pass must agree with the topology measurement, else one of the
    # two volume integrations is wrong
    mesh_volume = verdict["topology"]["signed_volume"]
    hydro_volume = verdict["hydrostatics"]["displacement_volume"]
    scale = max(abs(mesh_volume), abs(hydro_volume), 1e-12)
    verdict["volume_integration_agreement"] = abs(mesh_volume - hydro_volume) / scale
    verdict["volume_integration_ok"] = verdict["volume_integration_agreement"] < 1e-6
    if convergence:
        verdict["convergence"] = _convergence(record, nx=nx, nz=nz)
    return verdict


def _negative_controls(manifest: FleetManifest, out_dir: Path) -> list[dict[str, Any]]:
    """Every control here must be rejected; the run records the measured reason."""
    controls: list[dict[str, Any]] = []

    # 1. the 304-triangle placeholder signature
    placeholder = placeholder_box(304)
    facts = measure(placeholder)
    reasons = placeholder_reasons(facts, mesh=placeholder)
    controls.append(
        {
            "control": "placeholder_304_triangles",
            "expected": "rejected",
            "rejected": bool(reasons),
            "reasons": reasons,
            "measurements": {"triangle_count": facts.triangle_count, "vertices": facts.vertex_count,
                             "signed_volume": facts.signed_volume, "boundary_edges": facts.boundary_edge_count},
        }
    )

    # 2. two ships sharing one byte-identical geometry file
    control_dir = out_dir / "controls"
    first, second = duplicate_stl_pair(control_dir)
    sha_first, sha_second = sha256_file(first), sha256_file(second)
    duplicate = duplicate_manifest(manifest, first, second)
    collisions = duplicate.duplicate_geometry()
    verdicts = [evaluate_gates(ship, placeholder, duplicate) for ship in duplicate.ships]
    controls.append(
        {
            "control": "byte_identical_pair",
            "expected": "rejected",
            "rejected": bool(collisions) and all(not v["passed"] for v in verdicts),
            "sha256_first": sha_first,
            "sha256_second": sha_second,
            "files_identical": sha_first == sha_second,
            "collisions": {sha[:16]: list(names) for sha, names in collisions.items()},
            "failed_gates": [v["failed_gates"] for v in verdicts],
        }
    )

    # 3. Cb clipped into the form family's window: builds, but the volume gate must fail
    base = next((ship for ship in manifest.ships if ship.name == "KCS"), manifest.ships[0])
    clipped = clipped_cb_record(base)
    try:
        built = build_ship_geometry(clipped, nx=41, nz=13)
        verdict = evaluate_gates(clipped, built.mesh, manifest)
        failed = verdict["failed_gates"]
        controls.append(
            {
                "control": "clipped_cb",
                "expected": "rejected",
                "rejected": bool(failed),
                "failed_gates": failed,
                "volume_relative_error": verdict["volume_relative_error"],
                "detail": next(g["detail"] for g in verdict["gates"] if g["name"] == "G7_volume_vs_reference"),
            }
        )
    except RouteRefusal as exc:  # pragma: no cover - only if the control itself is invalid
        controls.append({"control": "clipped_cb", "expected": "rejected", "rejected": True,
                         "detail": f"route refused the control outright: {exc}", "failed_gates": ["G_build"]})

    # 4. a record whose declared displacement contradicts its own particulars
    wrong = wrong_volume_record(base)
    gate = gate_particulars_consistent(wrong)
    controls.append(
        {
            "control": "incoherent_particulars",
            "expected": "rejected",
            "rejected": not gate.passed,
            "failed_gates": [] if gate.passed else [gate.name],
            "detail": gate.detail,
        }
    )
    return controls


def _cross_route_check(manifest: FleetManifest, *, nx: int, nz: int) -> dict[str, Any]:
    """Loft one ship by the *other* route and compare the two volumes."""
    base = next((ship for ship in manifest.ships if ship.name == "KCS"), manifest.ships[0])
    primary = build_ship_geometry(base, nx=nx, nz=nz)
    import dataclasses as _dc
    other = _dc.replace(base, geometry_route="parametric_sections")
    lofted = build_ship_geometry(other, nx=nx, nz=nz).mesh
    v_primary = meshlib.total_signed_volume(primary.mesh)
    v_loft = meshlib.total_signed_volume(lofted)
    scale = max(abs(v_primary), abs(v_loft), 1e-12)
    return {
        "ship": base.name,
        "primary_route": primary.route,
        "secondary_route": "parametric_sections (independent Beta section-area law)",
        "volume_primary": v_primary,
        "volume_lofted": v_loft,
        "relative_difference": abs(v_primary - v_loft) / scale,
        "agree_within_1pct": abs(v_primary - v_loft) / scale < 0.01,
        "lofted_watertight": meshlib.is_watertight(lofted),
        "lofted_boundary_edges": len(meshlib.boundary_edges(lofted)),
        "lofted_euler": len(lofted.vertices) - len(meshlib.edge_face_counts(lofted)) + len(lofted.faces),
    }


def run_fleet(
    out_dir: str | Path,
    *,
    manifest_path: str | Path | None = None,
    nx: int = 481,
    nz: int = 97,
    volume_tolerance: float = VOLUME_TOLERANCE,
    write: bool = True,
) -> dict[str, Any]:
    """Build, gate and record the whole fleet. Returns the evidence dictionary."""
    manifest = load_manifest(manifest_path)
    root = Path(out_dir)
    ships = [
        _ship_block(record, manifest, nx=nx, nz=nz, volume_tolerance=volume_tolerance)
        for record in manifest.ships
    ]
    controls = _negative_controls(manifest, root)
    cross_route = _cross_route_check(manifest, nx=max(41, nx // 3), nz=nz)

    passed = [block["ship"] for block in ships if block["passed"]]
    failed = [block["ship"] for block in ships if not block["passed"]]
    measured = [block for block in ships if block["topology"] is not None]
    worst_error = max((block["volume_relative_error"] for block in measured), default=None)
    worst_ship = max(measured, key=lambda b: b["volume_relative_error"])["ship"] if measured else None
    external = [block for block in measured if block["reference_is_external"]]
    extrapolated = [
        block["convergence"]["richardson_relative_error_percent"] / 100.0
        for block in measured
        if block.get("convergence") and block["convergence"].get("richardson_relative_error_percent") is not None
    ]
    orders = [
        block["convergence"]["observed_order"]
        for block in measured
        if block.get("convergence") and block["convergence"].get("observed_order") is not None
    ]

    evidence: dict[str, Any] = {
        "schema": SCHEMA,
        "anchor": "GAP-CAD-T1-BASELINE-FLEET",
        "mission": "m17335",
        "settings": {
            "nx": nx,
            "nz": nz,
            "volume_tolerance": volume_tolerance,
            "manifest_path": None if manifest.path is None else str(manifest.path),
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        },
        "provenance": fleet_sha_index(manifest),
        "ships": ships,
        "cross_route_check": cross_route,
        "negative_controls": controls,
        "summary": {
            "ships_registered": len(manifest.ships),
            "ships_admitted": len(passed),
            "ships_refused_or_failed": len(failed),
            "admitted": passed,
            "failed": failed,
            "worst_volume_relative_error": worst_error,
            "worst_volume_ship": worst_ship,
            "ships_with_external_reference": [block["ship"] for block in external],
            "controls_all_rejected": all(control["rejected"] for control in controls),
            "negative_control_count": len(controls),
            "cross_route_agrees": cross_route["agree_within_1pct"],
            "worst_richardson_relative_error": max(extrapolated) if extrapolated else None,
            "observed_convergence_orders": orders,
            "volume_integration_all_agree": all(
                block.get("volume_integration_ok", False) for block in measured
            ),
            "hydrostatics_solver_calls": len(measured),
        },
    }

    if write:
        root.mkdir(parents=True, exist_ok=True)
        (root / "manifests").mkdir(exist_ok=True)
        (root / "evidence.json").write_text(json.dumps(evidence, indent=2, default=_json_default) + "\n")
        (root / "settings.json").write_text(json.dumps(evidence["settings"], indent=2) + "\n")
        _write_csv(root / "ref_vs_got.csv", ships)
        for block in ships:
            (root / "manifests" / f"{block['ship']}.json").write_text(
                json.dumps(
                    {
                        "ship": block["ship"],
                        "record": block["record"],
                        "verdict": {
                            "passed": block["passed"],
                            "failed_gates": block["failed_gates"],
                            "reference_basis": block["reference_basis"],
                            "volume_relative_error": block["volume_relative_error"],
                        },
                        "topology": block["topology"],
                        "gates": block["gates"],
                    },
                    indent=2,
                    default=_json_default,
                )
                + "\n"
            )
    return evidence


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    raise TypeError(f"not JSON serialisable: {type(value).__name__}")


def _write_csv(path: Path, ships: list[dict[str, Any]]) -> None:
    """The benchmark table: reference vs gotten, with the basis of the reference stated."""
    columns = [
        "ship",
        "reference_basis",
        "reference_is_external",
        "reference_volume_m3",
        "mesh_volume_m3",
        "volume_rel_error",
        "tolerance",
        "within_tolerance",
        "triangles",
        "boundary_edges",
        "euler_characteristic",
        "passed",
        "failed_gates",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for block in ships:
            record = block["record"]
            topology = block["topology"] or {}
            volume_gate = next(
                (gate for gate in block["gates"] if gate["name"] == "G7_volume_vs_reference"), None
            )
            writer.writerow(
                [
                    block["ship"],
                    block["reference_basis"],
                    block["reference_is_external"],
                    f"{record['reference_volume']:.6f}",
                    "" if not topology else f"{topology['signed_volume']:.6f}",
                    "" if block["volume_relative_error"] is None else f"{block['volume_relative_error']:.8f}",
                    0.01,
                    "" if volume_gate is None else volume_gate["passed"],
                    topology.get("triangle_count", ""),
                    topology.get("boundary_edge_count", ""),
                    topology.get("euler_characteristic", ""),
                    block["passed"],
                    ";".join(block["failed_gates"]),
                ]
            )


def main(argv: list[str] | None = None) -> int:  # pragma: no cover - thin wrapper
    from tensorcad.fleet.cli import main as cli_main

    return cli_main(argv)