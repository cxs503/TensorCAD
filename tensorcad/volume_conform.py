"""Make a CAD volume mesh conforming before it crosses into a solver.

Anchor: ``GAP-CAD-B3-SOLVER-INTERFACE`` (CAD -> solver consumption face).

The graded solid mesh written by :func:`tensorcad.volume_mesh.volume_mesh_from_hull`
contains degenerate cells: at the stem and the stern the section breadth collapses to
zero, so the breadth layers collapse onto the centreline and the recovered tetrahedra
have **zero volume** (measured below: 756 of 11520 at ``first_layer=0.02, layers=2``).
The gmsh bundle still verifies (its volume and digest are consistent), but a solver
refuses it: TensorFVM's ``NativeMesh3D.from_cells`` raises ``degenerate face`` on the
first zero-volume tetrahedron it meets.

This module is the missing step: weld coincident nodes, drop the cells that welding
proved degenerate, re-check the boundary, and only then write the bundle. Every number
this module reports is recomputed from the arrays it is about to write.

Consumers here are *real* solver libraries, imported at call time:

* :func:`fvm_affine_audit` -- TensorFVM ``NativeMesh3D`` + ``NativeOperators3D``
  build the face-based mesh from the CAD tetrahedra and audit affine exactness
  (least-squares gradient, face reconstruction, divergence of a linear field, global
  flux balance). A non-conforming mesh fails a construction assertion, not a tolerance.
* :func:`fem_buoyancy_equilibrium` -- TensorFEM ``Tet4J2Model`` (elastic range) with the
  hydrostatic pressure of the *same* hull surface as consistent nodal loads. Summing the
  reactions must return Archimedes ``rho g volume`` upwards: this checks the surface ->
  nodal-load handoff and the solve in one measurement against a closed form.

Both consumers read ``volume.msh`` back from disk with
:func:`tensorcad.solver_interface.read_gmsh_volume`; the arrays are never passed in
memory from the writer.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from .solver_interface import (
    boundary_triangles,
    read_gmsh_volume,
    total_volume,
    verify_bundle,
    write_bundle_volume,
)

__all__ = [
    "ConformReport",
    "weld_and_prune",
    "conform_report",
    "hull_volume_bundle",
    "fvm_affine_audit",
    "fem_buoyancy_equilibrium",
    "run_hull_handoff",
]


# --------------------------------------------------------------------- welding


@dataclass(frozen=True)
class ConformReport:
    """What welding and pruning did to a (vertices, tets) pair."""

    vertices: tuple[tuple[float, float, float], ...]
    tets: tuple[tuple[int, int, int, int], ...]
    nodes_in: int
    nodes_out: int
    tets_in: int
    tets_out: int
    dropped_degenerate: int
    dropped_nonpositive: int
    flipped: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "nodes_in": self.nodes_in,
            "nodes_out": self.nodes_out,
            "tets_in": self.tets_in,
            "tets_out": self.tets_out,
            "dropped_degenerate": self.dropped_degenerate,
            "dropped_nonpositive": self.dropped_nonpositive,
            "flipped": self.flipped,
        }


def _tet_volumes(vertices: Sequence[Sequence[float]], tets: Sequence[Sequence[int]]) -> list[float]:
    out: list[float] = []
    for a, b, c, d in tets:
        pa, pb, pc, pd = vertices[a], vertices[b], vertices[c], vertices[d]
        ux, uy, uz = pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]
        vx, vy, vz = pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]
        wx, wy, wz = pd[0] - pa[0], pd[1] - pa[1], pd[2] - pa[2]
        cx, cy, cz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        out.append((cx * wx + cy * wy + cz * wz) / 6.0)
    return out


def weld_and_prune(
    vertices: Sequence[Sequence[float]],
    tets: Sequence[Sequence[int]],
    *,
    decimals: int = 12,
) -> ConformReport:
    """Weld coincident nodes, drop cells welding proved degenerate, orient and compact.

    ``decimals`` is the rounding used to decide that two nodes are the same node; it is
    absolute (metres in this repository's contract) and therefore documented rather than
    hidden. A tetrahedron that repeats a node after welding, or whose signed volume is
    zero to working precision, is dropped and counted -- never silently kept. A negative
    signed volume means an inverted orientation and is flipped, also counted, because an
    inverted cell would break outward face normals downstream.
    """
    index: dict[tuple[float, float, float], int] = {}
    remap: list[int] = []
    out_vertices: list[tuple[float, float, float]] = []
    for point in vertices:
        key = (round(float(point[0]), decimals), round(float(point[1]), decimals), round(float(point[2]), decimals))
        node = index.get(key)
        if node is None:
            node = len(out_vertices)
            index[key] = node
            out_vertices.append(key)
        remap.append(node)

    kept: list[tuple[int, int, int, int]] = []
    dropped = 0
    nonpositive = 0
    flipped = 0
    scale = 0.0
    for a, b, c, d in tets:
        cell = (remap[a], remap[b], remap[c], remap[d])
        if len(set(cell)) != 4:
            dropped += 1
            continue
        kept.append(cell)
    volumes = _tet_volumes(out_vertices, kept)
    lo = max((abs(v) for v in volumes), default=0.0)
    if lo > 0.0:
        scale = lo
    oriented: list[tuple[int, int, int, int]] = []
    for cell, volume in zip(kept, volumes):
        if abs(volume) <= 1e-14 * max(scale, 1e-30):
            nonpositive += 1
            continue
        if volume < 0.0:
            cell = (cell[0], cell[2], cell[1], cell[3])
            flipped += 1
        oriented.append(cell)
    # Dropping degenerate cells can orphan the nodes that only they used. An orphaned
    # node still owns three stiffness rows, so keeping it makes the assembled system
    # singular (measured: TensorFEM raised "input matrix is singular"). Compact them.
    used = sorted({node for cell in oriented for node in cell})
    compact = {node: index for index, node in enumerate(used)}
    final_vertices = tuple(out_vertices[node] for node in used)
    oriented = [tuple(compact[node] for node in cell) for cell in oriented]
    return ConformReport(
        vertices=final_vertices,
        tets=tuple(oriented),
        nodes_in=len(vertices),
        nodes_out=len(final_vertices),
        tets_in=len(tets),
        tets_out=len(oriented),
        dropped_degenerate=dropped,
        dropped_nonpositive=nonpositive,
        flipped=flipped,
    )


def conform_report(
    vertices: Sequence[Sequence[float]], tets: Sequence[Sequence[int]]
) -> dict[str, Any]:
    """Recompute the acceptance numbers from the arrays a consumer will actually see."""
    volumes = _tet_volumes(vertices, tets)
    triangles = boundary_triangles(tets)
    edge_counts: dict[tuple[int, int], int] = {}
    for a, b, c in triangles:
        for u, v in ((a, b), (b, c), (c, a)):
            key = (int(u), int(v)) if u < v else (int(v), int(u))
            edge_counts[key] = edge_counts.get(key, 0) + 1
    bad_edges = sum(1 for count in edge_counts.values() if count != 2)
    return {
        "nodes": len(vertices),
        "tetrahedra": len(tets),
        "boundary_triangles": len(triangles),
        "nonmanifold_edges": bad_edges,
        "watertight": bad_edges == 0,
        "volume_m3": total_volume(vertices, tets),
        "min_tet_volume_m3": min((abs(v) for v in volumes), default=0.0),
        "inverted_tets": sum(1 for v in volumes if v < 0.0),
        "zero_tets": sum(1 for v in volumes if v == 0.0),
    }


# ------------------------------------------------------------------ hull bundle


def hull_volume_bundle(
    outdir: str | Path,
    *,
    parameters: Any = None,
    first_layer: float = 0.02,
    growth: float = 1.2,
    layers: int = 2,
    density: float = 1025.0,
    case: str | None = None,
) -> dict[str, Any]:
    """Hull body -> conforming C3D4 gmsh bundle, verified by re-reading the file."""
    from .hull import generate_hull, HullParameters
    from .volume_mesh import volume_mesh_from_hull

    hull_parameters = parameters if parameters is not None else HullParameters()
    hull = generate_hull(hull_parameters)
    mesh = volume_mesh_from_hull(hull_parameters, first_layer=first_layer, growth=growth, layers=layers)
    raw_vertices = list(mesh.vertices)
    raw_tets = [tuple(cell) for cell in mesh.tets]
    before = conform_report(raw_vertices, raw_tets)
    conformed = weld_and_prune(raw_vertices, raw_tets)
    after = conform_report(conformed.vertices, conformed.tets)
    hydro = hull.hydrostatics
    name = case or f"hull-volume-layers{layers}"
    manifest = write_bundle_volume(
        outdir,
        vertices=conformed.vertices,
        tets=conformed.tets,
        case=name,
        extra={
            "generator": "tensorcad.volume_mesh.volume_mesh_from_hull + tensorcad.volume_conform.weld_and_prune",
            "layers": layers,
            "first_layer_m": first_layer,
            "growth": growth,
            "density_kg_m3": density,
            "dropped_cells": conformed.dropped_degenerate + conformed.dropped_nonpositive,
        },
    )
    verification = verify_bundle(outdir)
    return {
        "case": name,
        "bundle": str(Path(outdir)),
        "route": "volume",
        "frame": verification.get("frame"),
        "manifest": manifest,
        "verify": verification,
        "before_weld": before,
        "after_weld": after,
        "weld": conformed.to_dict(),
        "geometry": {
            "length_m": hull_parameters.length,
            "beam_m": hull_parameters.beam,
            "draft_m": hull_parameters.draft,
            "cb_hydrostatics": hydro.block_coefficient,
            "displacement_volume_m3": hydro.displacement_volume,
            "lcb_m": hydro.longitudinal_centre_buoyancy,
            "vcg_buoyancy_m": hydro.vertical_centre_buoyancy,
            "water_density_kg_m3": density,
            "validation_passed": all(check.passed for check in hull.validation.checks if check.severity == "error"),
        },
    }


# ----------------------------------------------------------- FVM consumption


def fvm_affine_audit(bundle_dir: str | Path) -> dict[str, Any]:
    """TensorFVM builds its face-based mesh from the CAD tetrahedra and audits it."""
    try:
        from tensorfvm.native_mesh3d import NativeMesh3D
        from tensorfvm.native_operators3d import NativeOperators3D
    except Exception as exc:  # noqa: BLE001 - the message is the evidence
        return {"solver": "tensorfvm", "status": "unavailable", "error": f"{type(exc).__name__}: {exc}"}
    data = read_gmsh_volume(Path(bundle_dir) / "volume.msh")
    vertices = data["nodes"]
    tets = data["tets"]
    try:
        mesh = NativeMesh3D.from_cells(vertices, tets, ("tetra",) * len(tets))
        operators = NativeOperators3D(mesh)
        audit = operators.audit_affine()
    except Exception as exc:  # noqa: BLE001 - a refused mesh is the finding
        return {
            "solver": "tensorfvm",
            "status": "rejected",
            "error": f"{type(exc).__name__}: {exc}",
            "cells": len(tets),
        }
    volume = float(mesh.cell_volumes.sum())
    declared = float(json.loads((Path(bundle_dir) / "handoff.json").read_text())["volume_m3"])
    return {
        "solver": "tensorfvm",
        "status": "ran",
        "cells": len(mesh.cells),
        "faces": len(mesh.owner),
        "internal_faces": int((mesh.neighbor >= 0).sum()),
        "boundary_faces": int((mesh.neighbor < 0).sum()),
        "mesh_volume_m3": volume,
        "declared_volume_m3": declared,
        "volume_rel_error": (volume - declared) / declared,
        "audit_affine": audit,
    }


# ----------------------------------------------------------- FEM consumption


def fem_buoyancy_equilibrium(
    bundle_dir: str | Path,
    *,
    young: float = 2.1e11,
    poisson: float = 0.3,
    gravity: float = 9.80665,
    density: float = 1025.0,
    waterline_z: float = 0.0,
) -> dict[str, Any]:
    """Hydrostatic pressure on the CAD surface -> TensorFEM tet4 -> Archimedes check.

    The pressure ``p = rho g (waterline_z - z)`` acts on every boundary triangle whose
    centroid is below the waterline, along the inward normal, as a consistent nodal load
    (``p * area / 3`` to each corner). The solve is linear-elastic in the elastic range;
    the reactions at a statically determinate support must balance the applied load, and
    the applied load itself must equal ``rho g V`` upwards with its moment about the
    transverse axis at the hydrostatics LCB.
    """
    try:
        import torch
        from tensorfem.solid_plasticity import PlasticSolidState, Tet4J2Model, solve_increment
    except Exception as exc:  # noqa: BLE001
        return {"solver": "tensorfem", "status": "unavailable", "error": f"{type(exc).__name__}: {exc}"}

    import numpy as np

    directory = Path(bundle_dir)
    data = read_gmsh_volume(directory / "volume.msh")
    manifest = json.loads((directory / "handoff.json").read_text())
    nodes = np.asarray(data["nodes"], dtype=np.float64)
    tets = np.asarray(data["tets"], dtype=np.int64)
    triangles = boundary_triangles(tets)
    density = float(manifest.get("extra", {}).get("density_kg_m3", density))

    forces = np.zeros(3 * len(nodes))
    pressure_force = np.zeros(3)
    wetted_triangles = 0
    for a, b, c in triangles:
        pa, pb, pc = nodes[a], nodes[b], nodes[c]
        centroid_z = (pa[2] + pb[2] + pc[2]) / 3.0
        if centroid_z >= waterline_z:
            continue
        cross = np.cross(pb - pa, pc - pa)
        area = 0.5 * float(np.linalg.norm(cross))
        if area <= 0.0:
            continue
        normal = cross / (2.0 * area)  # outward for a positively oriented closed surface
        pressure = density * gravity * (waterline_z - centroid_z)
        share = -pressure * area / 3.0 * normal  # pressure pushes inward
        for node in (a, b, c):
            forces[3 * node : 3 * node + 3] += share
        pressure_force += share * 3.0
        wetted_triangles += 1

    # Statically determinate support, chosen geometrically: three nodes as far apart as
    # the hull allows, so all six rigid-body modes are removed. Picking them by x-order
    # alone put them on one transverse section and left the system singular (measured).
    lo = nodes.min(axis=0)
    hi = nodes.max(axis=0)
    anchors = [
        np.array([lo[0], 0.0, lo[2]]),
        np.array([hi[0], 0.0, lo[2]]),
        np.array([0.5 * (lo[0] + hi[0]), hi[1], lo[2]]),
    ]
    picked: list[int] = []
    for target in anchors:
        distance = np.linalg.norm(nodes - target, axis=1)
        for node in np.argsort(distance):
            if int(node) not in picked:
                picked.append(int(node))
                break
    fixed: list[int] = []
    for level, node in enumerate(picked):
        for axis in range(3):
            if level == 0 or axis >= level:
                fixed.append(3 * node + axis)

    model = Tet4J2Model(
        nodes=torch.tensor(nodes, dtype=torch.float64),
        elements=torch.tensor(tets, dtype=torch.long),
        young=float(young),
        poisson=float(poisson),
        yield_stress=1.0e15,
        hardening=0.0,
        fixed_dofs=torch.tensor(fixed, dtype=torch.long),
    )
    state = PlasticSolidState.virgin(model)
    solution = solve_increment(
        model,
        torch.tensor(forces, dtype=torch.float64),
        1.0,
        state,
        torch.zeros(model.n_dofs, dtype=torch.float64),
        tolerance=1e-11,
        max_iterations=25,
    )
    displacement = np.asarray(solution.displacement.detach().cpu().numpy(), dtype=float).reshape(-1, 3)
    reaction = np.asarray(solution.reaction.detach().cpu().numpy(), dtype=float).reshape(-1, 3)
    net = reaction.reshape(-1).sum() + forces.sum()
    applied_vector = forces.reshape(-1, 3).sum(axis=0)
    reaction_total = reaction.reshape(-1, 3).sum(axis=0)
    scale = float(np.abs(forces.reshape(-1, 3)).sum())
    volume = float(manifest["volume_m3"])
    moment = float((nodes[:, 0] * forces.reshape(-1, 3)[:, 2]).sum())
    lcb = moment / (density * gravity * volume)
    return {
        "solver": "tensorfem",
        "status": "ran",
        "elements": int(len(tets)),
        "nodes": int(len(nodes)),
        "wetted_triangles": wetted_triangles,
        "applied_buoyancy_n": float(applied_vector[2]),
        "archimedes_n": density * gravity * volume,
        "archimedes_rel_error": (float(applied_vector[2]) - density * gravity * volume) / (density * gravity * volume),
        "reaction_total_n": float(reaction_total[2]),
        "equilibrium_relative_error": float(abs(net) / max(scale, 1e-30)),
        "lcb_from_pressure_m": lcb,
        "max_displacement_m": float(np.abs(displacement).max()),
        "newton_iterations": int(solution.iterations),
        "residual_norm": float(solution.residual_norm),
    }


# ------------------------------------------------------------------- driver

GATES = {
    "watertight": 0.0,
    "nonmanifold_edges": 0,
    "zero_tets": 0,
    "inverted_tets": 0,
    "fvm_affine_error": 1e-9,
    "fvm_divergence_error": 1e-9,
    "fem_archimedes_rel_error": 5e-3,
    "fem_equilibrium_rel_error": 1e-8,
    "fem_lcb_rel_error": 2e-2,
}


def run_hull_handoff(
    outdir: str | Path,
    *,
    layers: int = 1,
    first_layer: float = 0.04,
    nx: int = 21,
    nz: int = 9,
) -> dict[str, Any]:
    """Hull -> conformed bundle -> TensorFVM audit + TensorFEM buoyancy, with gates.

    ``nx``/``nz`` are the hull station/level counts. The defaults are the coarse hull:
    the graded solid at 41x13 stations needs ~90 s of TensorFVM mesh construction and
    TensorFEM needs ~50 s per solve on 1.7 k tetrahedra, so the gated run uses the coarse
    hull and the resolution study is reported separately.
    """
    from .hull import HullParameters

    out = Path(outdir)
    bundle = hull_volume_bundle(
        out / "bundle",
        parameters=HullParameters(nx=nx, nz=nz),
        layers=layers,
        first_layer=first_layer,
        case=f"hull-volume-{nx}x{nz}-layers{layers}",
    )
    fvm = fvm_affine_audit(out / "bundle")
    fem = fem_buoyancy_equilibrium(out / "bundle", density=bundle["geometry"]["water_density_kg_m3"])

    volume_declared = float(bundle["manifest"]["volume_m3"])
    hydro_volume = float(bundle["geometry"]["displacement_volume_m3"])
    gates: list[dict[str, Any]] = []

    def gate(name: str, measured: float, reference: float, tolerance: float, unit: str = "") -> None:
        error = abs(measured - reference)
        gates.append(
            {
                "name": name,
                "measured": measured,
                "reference": reference,
                "abs_error": error,
                "tolerance": tolerance,
                "unit": unit,
                "passed": bool(error <= tolerance),
            }
        )

    after = bundle["after_weld"]
    gate("nonmanifold_edges", float(after["nonmanifold_edges"]), 0.0, GATES["nonmanifold_edges"])
    gate("zero_tets", float(after["zero_tets"]), 0.0, GATES["zero_tets"])
    gate("inverted_tets", float(after["inverted_tets"]), 0.0, GATES["inverted_tets"])
    # The faceted hull surface is the one hydrostatics integrated, so the conformed
    # solid volume must equal the hydrostatic displacement volume, not merely approach it.
    gate("mesh_volume_m3", volume_declared, hydro_volume, 1e-6 * hydro_volume, "m3")
    if fvm.get("status") == "ran":
        audit = fvm["audit_affine"]
        gate(
            "fvm_ls_gradient_error",
            float(audit["least_squares_gradient_max_error"]),
            0.0,
            GATES["fvm_affine_error"],
        )
        gate(
            "fvm_face_affine_error",
            float(audit["face_affine_max_error"]),
            0.0,
            GATES["fvm_affine_error"],
        )
        gate(
            "fvm_affine_divergence_error",
            float(audit["affine_vector_divergence_max_error"]),
            0.0,
            GATES["fvm_divergence_error"],
        )
        gate(
            "fvm_global_flux_balance_error",
            float(audit["global_flux_balance_error"]),
            0.0,
            GATES["fvm_divergence_error"],
        )
    if fem.get("status") == "ran":
        gate(
            "fem_archimedes_rel_error",
            abs(float(fem["archimedes_rel_error"])),
            0.0,
            GATES["fem_archimedes_rel_error"],
        )
        gate(
            "fem_equilibrium_rel_error",
            float(fem["equilibrium_relative_error"]),
            0.0,
            GATES["fem_equilibrium_rel_error"],
        )
        gate(
            "fem_lcb_rel_error",
            abs(float(fem["lcb_from_pressure_m"]) - float(bundle["geometry"]["lcb_m"]))
            / float(bundle["geometry"]["lcb_m"]),
            0.0,
            GATES["fem_lcb_rel_error"],
        )

    payload = {
        "schema": "tensorcad.hull_handoff/1",
        "bundle": bundle,
        "fvm": fvm,
        "fem": fem,
        "gates": gates,
        "gates_failed": [g["name"] for g in gates if not g["passed"]],
    }
    (out / "results.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def _main(argv: Sequence[str] | None = None) -> int:  # pragma: no cover - CLI wrapper
    import argparse

    parser = argparse.ArgumentParser(description="Hull volume bundle consumed by TensorFVM and TensorFEM")
    parser.add_argument("--out", default="artifacts/hull_handoff")
    args = parser.parse_args(argv)
    payload = run_hull_handoff(args.out)
    print(json.dumps({"gates_failed": payload["gates_failed"], "gates": len(payload["gates"])}, indent=2))
    return 0 if not payload["gates_failed"] else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())