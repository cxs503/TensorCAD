"""Consume the watertight-solid handoff the way a finite-volume solver does.

Why this module exists
----------------------
:data:`tensorcad.solver_interface.COORDINATE_CONTRACT` describes four routes, and
three of them already had a consumer (:mod:`tensorcad.solver_consumers`: a surface
table for BEM, a 2-D grid for a FEM plate, an occupancy lattice for LBM).  The
fourth -- **watertight surface -> tetrahedral volume mesh -> FVM** -- was written
and re-read by the contract verifier but never solved on: the FVM leg of the
earlier round was a 2-D planar channel, so "a watertight CAD solid reaches a
volume solver" was an assertion with no number behind it.

Here the chain is closed and each link is machine-checked:

1. A **CAD solid** comes out of this repository's own volume mesher
   (:func:`tensorcad.volume_mesh.build_volume_mesh`, a hull section grid) and is
   written by :func:`tensorcad.solver_interface.write_bundle_volume` as gmsh 2.2
   ASCII: one ``C3D4`` tetrahedron block plus one physically tagged ``wall``
   triangle block, in the contract frame (z-up, metres, waterline ``z = 0``,
   keel ``z = -T``).
2. :func:`tensorcad.solver_interface.verify_bundle` re-reads that file with an
   independent parser and refuses a moved node, a wrong digest or a non-``z-up``
   frame.  It runs first and is a hard gate.
3. The consumer parses the file **from disk** and hands the vertices and
   tetrahedra to TensorFVM's own mesh object,
   ``tensorfvm.native_mesh3d.NativeMesh3D.from_cells`` -- the same constructor
   ``from_gmsh`` uses once the gmsh SDK has read the file.  TensorFVM then owns
   the topology: it rebuilds faces, orients them, rejects non-manifold or
   non-convex cells, and audits closure.
4. The solve is steady scalar conduction ``-div(kappa grad phi) = s`` with
   Dirichlet data on the tagged wall faces, by
   ``tensorfvm.native_hybrid_scalar3d.NativeHybridScalar3D`` -- conservative
   cell-centred finite volume.

The references are closed-form and are computed here, never taken from the
solver:

``patch``
    A linear field ``phi = a + b.x + c.y + d.z`` is harmonic, so conduction with
    Dirichlet data taken from it must be reproduced **exactly** on any conforming
    mesh.  This is a real patch test, in contrast with the earlier round's
    ``patch_test_relative_linf`` headline, which an independent audit showed was
    measuring a boundary-layer incompatibility (a Poisson's-ratio control
    experiment moved it by 985x) rather than discretisation error.
``harmonic``
    ``phi = sin(pi x / L) sin(pi y / B) sinh(gamma (z + T))`` with
    ``gamma = pi sqrt(1/L^2 + 1/B^2)`` is an exact harmonic function on the
    barge box ``[0, L] x [0, B] x [-T, 0]``, so the same Dirichlet problem has a
    non-trivial exact solution to compare against, and the error can be followed
    under refinement.
``frame``
    The same patch field applied to the *same* table after a y-up/z-up mix-up,
    which is the mistake the 5415 rounds paid for.  The measured error says how
    large that mistake is on this geometry -- a number, not a warning.

Findings recorded here rather than hidden
-----------------------------------------
* A *collapsed* column in a hull section grid (breadth -> 0 at the keel line or
  at a pointed stem) is accepted by :mod:`tensorcad.volume_mesh` but its tet split
  then contains topologically degenerate faces, which TensorFVM's auditor
  refuses with ``degenerate face``.  The consequence for the handoff is that a
  CAD solid destined for a volume solver must not be emitted with collapsed
  cells; the hull used here therefore has a flat keel width and a transom (both
  ordinary features of a real hull), and the point is stated in
  ``docs/fvm_handoff.md`` instead of being papered over by dropping cells (which
  would leave the surface open).
* ``NativeMesh3D.from_gmsh`` needs the gmsh SDK, which has no loongarch64 wheel
  in this environment, so the file is parsed by this repository's reader and the
  arrays enter TensorFVM through ``from_cells``.  What is *not* claimed: that
  ''the solver opened the gmsh file itself''.  What is claimed: the numbers the
  solver used are the ones on disk (digest-verified) and the mesh topology is the
  solver's own.
"""
from __future__ import annotations

import importlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

from .solver_consumers import (
    CaseResult,
    SolverUnavailableError,
    gate,
    solver_fingerprint,
    watertight_report,
)
from .solver_interface import (
    COORDINATE_CONTRACT,
    HandoffError,
    read_gmsh_volume,
    verify_bundle,
    write_bundle_volume,
)
from .volume_mesh import SectionGrid, VolumeMesh, build_volume_mesh, box_section_grid

__all__ = [
    "FLAT_KEEL_FRACTION",
    "TRANSOM_FRACTION",
    "barge_bundle",
    "exact_harmonic",
    "exact_linear",
    "fvm_harmonic_from_bundle",
    "fvm_patch_from_bundle",
    "frame_counterexample",
    "hull_section_grid",
    "hull_volume_bundle",
    "run_all_cases",
    "summary_markdown",
]

#: Fraction of the maximum half breadth still present at the keel line.  A real
#: hull has a flat of the bottom; a value of ``0`` collapses the column and its
#: tet split then contains degenerate faces (see the module docstring).
FLAT_KEEL_FRACTION = 0.10
#: Fraction of the maximum half breadth still present at the end stations, i.e.
#: a transom rather than a knife-edge stem.
TRANSOM_FRACTION = 0.15


# --------------------------------------------------------------------------- fields


def exact_linear(
    points: Sequence[Sequence[float]],
    coefficients: Sequence[float] = (1.0, 0.3, -0.7, 0.2),
) -> Any:
    """``a + b.x + c.y + d.z``, harmonic by construction."""
    import numpy as np

    a, b, c, d = (float(v) for v in coefficients)
    p = np.atleast_2d(np.asarray(points, dtype=float))
    return a + b * p[:, 0] + c * p[:, 1] + d * p[:, 2]


def exact_harmonic(points: Sequence[Sequence[float]], *, length: float, beam: float, draft: float) -> Any:
    """Exact harmonic field on the box ``[0,L] x [0,B] x [-T,0]`` (separation of variables).

    ``gamma`` follows from ``sin(pi x/L) sin(pi y/B)`` satisfying the 2-D
    Laplace equation, so the ``sinh`` factor makes the 3-D field harmonic.
    """
    import numpy as np

    p = np.atleast_2d(np.asarray(points, dtype=float))
    gamma = math.pi * math.sqrt(1.0 / length**2 + 1.0 / beam**2)
    return (
        np.sin(math.pi * p[:, 0] / length)
        * np.sin(math.pi * p[:, 1] / beam)
        * np.sinh(gamma * (p[:, 2] + draft))
    )


# --------------------------------------------------------------------------- geometry


def hull_section_grid(
    *,
    length: float = 10.0,
    beam: float = 1.0,
    draft: float = 0.6,
    nx: int = 13,
    nz: int = 7,
) -> SectionGrid:
    """A Wigley-derived hull with a flat keel of the bottom and a transom.

    The half breadth is the Wigley paraboloid times two shrinkage factors that
    stop it reaching zero at the keel line and at the end stations, which keeps
    every mesh column non-degenerate (see the module docstring).  The body is
    still a hull: waterline at ``z = 0``, keel at ``z = -draft``, single-valued
    half breadth in ``y``.
    """
    xs = tuple(length * i / (nx - 1) for i in range(nx))
    zs = tuple(-draft + draft * j / (nz - 1) for j in range(nz))

    def breadth(i: int, j: int) -> float:
        x, z = xs[i], zs[j]
        end = TRANSOM_FRACTION + (1.0 - TRANSOM_FRACTION) * (1.0 - (2.0 * x / length - 1.0) ** 2)
        keel = FLAT_KEEL_FRACTION + (1.0 - FLAT_KEEL_FRACTION) * (1.0 - (z / draft) ** 2)
        return 0.5 * beam * end * keel

    return SectionGrid(xs=xs, zs=zs, breadth=breadth)


def _write_bundle(outdir: Path, mesh: VolumeMesh, *, case: str, extra: dict[str, Any]) -> dict[str, Any]:
    return write_bundle_volume(
        outdir,
        vertices=mesh.vertices,
        tets=mesh.tets,
        case=case,
        extra={"mesh_route": "tensorcad.volume_mesh.build_volume_mesh(split_tets=True)", **extra},
    )


def barge_bundle(outdir: str | Path, *, length: float, beam: float, draft: float, nx: int, nz: int, layers: int) -> dict[str, Any]:
    """Write a watertight box-section hull solid (the analytic-control geometry)."""
    grid = box_section_grid(length=length, beam=beam, draft=draft, nx=nx, nz=nz)
    mesh = build_volume_mesh(grid, first_layer=0.006, growth=1.2, layers=layers, split_tets=True)
    return _write_bundle(
        Path(outdir),
        mesh,
        case=f"barge_{length:g}x{beam:g}x{draft:g}_{nx}x{nz}x{layers}",
        extra={"analytic": "linear patch test and separation-of-variables harmonic field"},
    )


def hull_volume_bundle(outdir: str | Path, *, length: float, beam: float, draft: float, nx: int, nz: int, layers: int) -> dict[str, Any]:
    """Write a watertight *curved* hull solid -- the ship-shaped FVM input."""
    grid = hull_section_grid(length=length, beam=beam, draft=draft, nx=nx, nz=nz)
    mesh = build_volume_mesh(grid, first_layer=0.006, growth=1.2, layers=layers, split_tets=True)
    return _write_bundle(
        Path(outdir),
        mesh,
        case=f"hull_{length:g}x{beam:g}x{draft:g}_{nx}x{nz}x{layers}",
        extra={"analytic": "linear patch test only: a curved solid has no closed-form harmonic field here"},
    )


# --------------------------------------------------------------------------- consumption


def _pin_solver_path() -> str | None:
    """Make the TensorFVM checkout explicit, not whatever a ``.pth`` file injected.

    This host has several agent worktrees of TensorFVM on ``sys.path`` through
    editable installs, so an unprefixed ``import tensorfvm`` can resolve to a
    colleague's half-finished copy.  ``TENSORFVM_SRC`` wins, then the ecosystem
    checkout, and the resolved path and git sha are recorded in every result by
    :func:`tensorcad.solver_consumers.solver_fingerprint`.
    """
    import os
    import sys

    for candidate in (os.environ.get("TENSORFVM_SRC"), "/srv/tensor/workspace/TensorFVM/src"):
        if candidate and (Path(candidate) / "tensorfvm" / "__init__.py").is_file():
            if sys.path[:1] != [candidate]:
                sys.path.insert(0, candidate)
            return candidate
    return None


def _load(module_name: str) -> Any:
    try:
        return importlib.import_module(module_name)
    except Exception as exc:  # noqa: BLE001 - the text is the evidence
        raise SolverUnavailableError(f"{module_name}: {type(exc).__name__}: {exc}") from exc


def limit_gate(name: str, measured: float, limit: float, unit: str = "") -> dict[str, Any]:
    """An upper-bound check: ``measured <= limit`` passes.

    :func:`tensorcad.solver_consumers.gate` is a ref-vs-got ratio gate, which
    needs a non-zero physical reference; an error bound has none, so the bound is
    stated as such instead of being dressed up as a reference value.
    """
    value = float(measured)
    bound = float(limit)
    return {
        "name": name,
        "measured": value,
        "reference": bound,
        "rel_error": value,
        "tolerance": bound,
        "unit": unit,
        "kind": "upper_bound",
        "passed": bool(value <= bound),
    }


@dataclass(frozen=True)
class _OpenBundle:
    """What was actually read off disk, and the solver-side mesh built from it."""

    manifest: dict[str, Any]
    verdict: dict[str, Any]
    nodes: Any
    tets: Any
    mesh: Any
    module: Any


def _open_bundle(bundle_dir: str | Path) -> _OpenBundle:
    """Hard gate: verify the bundle, parse the table, let the solver build the mesh."""
    bundle = Path(bundle_dir)
    _pin_solver_path()
    module = _load("tensorfvm.native_mesh3d")
    _load("tensorfvm.native_hybrid_scalar3d")

    verdict = verify_bundle(bundle)
    if not verdict.get("ok", False):
        raise HandoffError(f"bundle at {bundle} did not survive verification: {verdict}")
    if verdict.get("route") != "volume":
        raise HandoffError(f"bundle at {bundle} re-read as route {verdict.get('route')!r}, not 'volume'")

    table = read_gmsh_volume(bundle / "volume.msh")  # second reader, not the writer's state
    import numpy as np

    nodes = np.asarray(table["nodes"], dtype=float)
    tets = np.asarray(table["tets"], dtype=int)
    if tets.size == 0:
        raise HandoffError(f"bundle at {bundle} carries no tetrahedra")
    cells = [tuple(int(v) for v in tet) for tet in tets]
    mesh = module.NativeMesh3D.from_cells(nodes, cells, ["tetra"] * len(cells))
    return _OpenBundle(
        manifest=json.loads((bundle / "handoff.json").read_text()),
        verdict=verdict,
        nodes=nodes,
        tets=tets,
        mesh=mesh,
        module=module,
    )


def _solve(
    opened: _OpenBundle,
    field: Callable[[Any], Any],
    *,
    conductivity: float = 1.0,
) -> dict[str, Any]:
    """Dirichlet conduction with wall data from ``field``; return errors and ledgers."""
    import numpy as np

    _pin_solver_path()
    scalar = _load("tensorfvm.native_hybrid_scalar3d")
    solver = scalar.NativeHybridScalar3D(opened.mesh, conductivity)
    boundary = opened.mesh.neighbor < 0
    values = np.zeros(len(opened.mesh.owner), dtype=float)
    values[boundary] = field(opened.mesh.face_centers[boundary])
    result = solver.solve(boundary_values=values)

    cells = np.asarray(result["cell_values"], dtype=float)
    centres = opened.mesh.cell_centers
    reference = field(centres)
    volume = opened.mesh.cell_volumes
    error = np.abs(cells - reference)
    scale = float(np.max(np.abs(reference))) or 1.0
    l2 = math.sqrt(float((volume * error**2).sum()) / max(float((volume * reference**2).sum()), 1e-300))
    return {
        "cell_values": cells,
        "cell_centers": centres,
        "reference_values": reference,
        "relative_linf": float(error.max() / scale),
        "relative_l2": l2,
        "linear_relative_residual": float(result["linear_relative_residual"]),
        "global_relative_balance": float(result["global_relative_balance"]),
        "cell_balance_max": float(result["cell_balance_max"]),
    }


def _mesh_facts(opened: _OpenBundle, report: dict[str, Any]) -> dict[str, Any]:
    audit = opened.mesh.audit()
    return {
        "cells": int(len(opened.mesh.cells)),
        "faces": int(len(opened.mesh.owner)),
        "boundary_faces": int(sum(1 for n in opened.mesh.neighbor if n < 0)),
        "mesh_volume_m3": float(opened.mesh.cell_volumes.sum()),
        "declared_volume_m3": float(opened.manifest["volume_m3"]),
        "volume_rel_error_vs_declared": (
            float(opened.mesh.cell_volumes.sum()) - float(opened.manifest["volume_m3"])
        )
        / float(opened.manifest["volume_m3"]),
        "maximum_relative_closure": float(audit["maximum_relative_closure"]),
        "minimum_face_projection": float(audit["minimum_face_projection"]),
        "maximum_nonorthogonality_degrees": float(audit["maximum_nonorthogonality_degrees"]),
        "watertight_bad_edges": int(report["bad_edges"]),
        "watertight_euler_characteristic": int(report["euler_characteristic"]),
        "frame": str(opened.manifest["frame"].get("name")),
        "frame_up_axis": str(opened.manifest["frame"].get("up_axis")),
    }


def _watertight(bundle_dir: str | Path) -> dict[str, Any]:
    """Edge census of the wall triangles, plus the Euler characteristic.

    A closed genus-0 surface has ``V - E + F = 2``; a hole or a duplicated use
    moves it, and so does a non-manifold wall.  This is measured on the file,
    not on the mesher's intent.
    """
    table = read_gmsh_volume(Path(bundle_dir) / "volume.msh")
    report = watertight_report(table["nodes"], table["triangles"])
    used = {v for face in table["triangles"] for v in face}
    report["euler_characteristic"] = len(used) - report["edges"] + len(table["triangles"])
    report["boundary_triangles"] = len(table["triangles"])
    return report


def fvm_patch_from_bundle(
    bundle_dir: str | Path,
    *,
    coefficients: Sequence[float] = (1.0, 0.3, -0.7, 0.2),
    tolerance: float = 1e-9,
) -> CaseResult:
    """Linear-field patch test: a consistent FVM operator must reproduce it exactly."""
    try:
        opened = _open_bundle(bundle_dir)
    except SolverUnavailableError as exc:
        return CaseResult(case=str(bundle_dir), route="volume_msh", solver="tensorfvm", status="unavailable", metrics={}, error=str(exc))
    field = lambda p: exact_linear(p, coefficients)  # noqa: E731
    solved = _solve(opened, field)
    metrics = {
        **_mesh_facts(opened, _watertight(bundle_dir)),
        "field": "linear (harmonic)",
        "coefficients": list(coefficients),
        "relative_linf": solved["relative_linf"],
        "relative_l2": solved["relative_l2"],
        "linear_relative_residual": solved["linear_relative_residual"],
        "global_relative_balance": solved["global_relative_balance"],
        "cell_balance_max": solved["cell_balance_max"],
    }
    gates = (
        limit_gate("patch_relative_linf", solved["relative_linf"], tolerance),
        gate("mesh_volume_vs_declared", metrics["mesh_volume_m3"], metrics["declared_volume_m3"], 1e-10, "m3"),
    )
    return CaseResult(
        case=str(opened.manifest.get("case", bundle_dir)),
        route="volume_msh",
        solver="tensorfvm.native_hybrid_scalar3d",
        status="ran",
        metrics=metrics,
        gates=gates,
        solver_info=solver_fingerprint(opened.module),
        note="Dirichlet wall data taken from an exact linear field; the solve must return the field itself",
    )


def fvm_harmonic_from_bundle(
    bundle_dir: str | Path,
    *,
    length: float,
    beam: float,
    draft: float,
    tolerance: float = 0.15,
) -> CaseResult:
    """Steady conduction against the separation-of-variables harmonic field."""
    try:
        opened = _open_bundle(bundle_dir)
    except SolverUnavailableError as exc:
        return CaseResult(case=str(bundle_dir), route="volume_msh", solver="tensorfvm", status="unavailable", metrics={}, error=str(exc))
    field = lambda p: exact_harmonic(p, length=length, beam=beam, draft=draft)  # noqa: E731
    solved = _solve(opened, field)
    metrics = {
        **_mesh_facts(opened, _watertight(bundle_dir)),
        "field": "sin(pi x/L) sin(pi y/B) sinh(gamma (z+T))",
        "gamma_1_per_m": math.pi * math.sqrt(1.0 / length**2 + 1.0 / beam**2),
        "relative_linf": solved["relative_linf"],
        "relative_l2": solved["relative_l2"],
        "linear_relative_residual": solved["linear_relative_residual"],
        "global_relative_balance": solved["global_relative_balance"],
        "cell_balance_max": solved["cell_balance_max"],
    }
    return CaseResult(
        case=str(opened.manifest.get("case", bundle_dir)),
        route="volume_msh",
        solver="tensorfvm.native_hybrid_scalar3d",
        status="ran",
        metrics=metrics,
        gates=(limit_gate("harmonic_relative_l2", solved["relative_l2"], tolerance),),
        solver_info=solver_fingerprint(opened.module),
        note="closed-form harmonic field on the box; error is discretisation, not handoff",
    )


def frame_counterexample(
    bundle_dir: str | Path,
    *,
    coefficients: Sequence[float] = (1.0, 0.3, -0.7, 0.2),
) -> CaseResult:
    """Price the z-up/y-up mistake on this geometry: same table, wrong frame.

    The field is defined in the contract frame.  Reading the *same* vertices as
    if ``y`` were up is exactly the mix-up the 5415 rounds paid for; applying
    that field to the swapped coordinates measures what the mistake costs.
    """
    try:
        opened = _open_bundle(bundle_dir)
    except SolverUnavailableError as exc:
        return CaseResult(case=str(bundle_dir), route="volume_msh", solver="tensorfvm", status="unavailable", metrics={}, error=str(exc))
    import numpy as np

    right = _solve(opened, lambda p: exact_linear(p, coefficients))

    _pin_solver_path()
    scalar = _load("tensorfvm.native_hybrid_scalar3d")
    solver = scalar.NativeHybridScalar3D(opened.mesh, 1.0)
    boundary = opened.mesh.neighbor < 0
    values = np.zeros(len(opened.mesh.owner), dtype=float)
    values[boundary] = exact_linear(opened.mesh.face_centers[boundary], coefficients)

    swapped = opened.mesh.cell_centers[:, [0, 2, 1]]  # the y-up reading of the same z-up table
    reference = exact_linear(swapped, coefficients)
    cells = np.asarray(right["cell_values"], dtype=float)
    error = np.abs(cells - reference)
    wrong = float(error.max() / max(float(np.abs(reference).max()), 1e-30))
    return CaseResult(
        case=str(opened.manifest.get("case", bundle_dir)),
        route="volume_msh",
        solver="tensorfvm.native_hybrid_scalar3d",
        status="ran",
        metrics={
            "zup_relative_linf": right["relative_linf"],
            "yup_relative_linf": wrong,
            "frame_price_ratio": (wrong / right["relative_linf"]) if right["relative_linf"] else float("inf"),
            "frame": str(opened.manifest["frame"].get("name")),
        },
        gates=(limit_gate("zup_relative_linf", right["relative_linf"], 1e-9),),
        solver_info=solver_fingerprint(opened.module),
        note="the same solved field compared with a z-up and a y-up rendering of the same coordinates",
    )


# --------------------------------------------------------------------------- driver


def run_all_cases(
    outdir: str | Path,
    *,
    barge: tuple[float, float, float] = (2.0, 0.4, 0.3),
    barge_resolutions: Sequence[tuple[int, int, int]] = ((9, 5, 2), (17, 9, 3)),
    hull: tuple[float, float, float] = (10.0, 1.0, 0.6),
    hull_resolutions: Sequence[tuple[int, int, int]] = ((13, 7, 2),),
) -> dict[str, Any]:
    """Write the CAD solids, consume each one through FVM, gate every number."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    cases: list[CaseResult] = []
    length, beam, draft = barge

    for nx, nz, layers in barge_resolutions:
        sub = out / f"barge_{nx}x{nz}x{layers}"
        barge_bundle(sub, length=length, beam=beam, draft=draft, nx=nx, nz=nz, layers=layers)
        cases.append(fvm_patch_from_bundle(sub))
        cases.append(fvm_harmonic_from_bundle(sub, length=length, beam=beam, draft=draft))
        if (nx, nz, layers) == tuple(barge_resolutions)[0]:
            cases.append(frame_counterexample(sub))

    hull_length, hull_beam, hull_draft = hull
    for nx, nz, layers in hull_resolutions:
        sub = out / f"hull_{nx}x{nz}x{layers}"
        hull_volume_bundle(sub, length=hull_length, beam=hull_beam, draft=hull_draft, nx=nx, nz=nz, layers=layers)
        result = fvm_patch_from_bundle(sub)
        cases.append(result)

    payload = {
        "frame": dict(COORDINATE_CONTRACT),
        "schema": "tensorcad.fvm_handoff/1",
        "cases": [c.to_dict() for c in cases],
        "all_gates_passed": all(not c.gates_failed for c in cases if c.status == "ran"),
        "statuses": sorted({c.status for c in cases}),
    }
    (out / "results.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (out / "summary.md").write_text(summary_markdown(payload), encoding="utf-8")
    return payload


def summary_markdown(payload: dict[str, Any]) -> str:
    """Render the result payload as the markdown tables the docs quote verbatim."""
    lines = ["| case | cells | check | value | bound / reference | verdict |", "|---|---|---|---|---|---|"]
    for case in payload["cases"]:
        if case["status"] != "ran":
            lines.append(f"| {case['case']} | - | solver | - | - | {case['status']} |")
            continue
        m = case["metrics"]
        for g in case["gates"]:
            if g.get("kind") == "upper_bound":
                bound = f"<= {g['tolerance']:.3g}"
            else:
                bound = f"ref {g['reference']:.6g} (rel {g['rel_error']:+.2e})"
            lines.append(
                f"| {case['case']} | {m.get('cells', '-')} | {g['name']} | {g['measured']:.6g} | {bound} | "
                f"{'PASS' if g['passed'] else 'FAIL'} |"
            )
    lines.append("")
    lines.append("| case | cells | quantity | value |")
    lines.append("|---|---|---|---|")
    for case in payload["cases"]:
        if case["status"] != "ran":
            continue
        m = case["metrics"]
        for key in (
            "relative_linf",
            "relative_l2",
            "linear_relative_residual",
            "global_relative_balance",
            "maximum_nonorthogonality_degrees",
            "watertight_bad_edges",
            "watertight_euler_characteristic",
            "volume_rel_error_vs_declared",
        ):
            if key in m:
                lines.append(f"| {case['case']} | {m.get('cells', '-')} | {key} | {m[key]:.6g} |")
        for key in ("zup_relative_linf", "yup_relative_linf", "frame_price_ratio"):
            if key in m:
                lines.append(f"| {case['case']} | - | {key} | {m[key]:.6g} |")
    lines.append("")
    lines.append(f"all_gates_passed = {payload['all_gates_passed']}; statuses = {payload['statuses']}")
    return "\n".join(lines) + "\n"


def _main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="artifacts/fvm_handoff")
    parser.add_argument("--quick", action="store_true", help="one barge and one hull resolution")
    args = parser.parse_args(argv)
    if args.quick:
        payload = run_all_cases(args.out, barge_resolutions=((9, 5, 2),), hull_resolutions=((13, 7, 2),))
    else:
        payload = run_all_cases(args.out)
    print(summary_markdown(payload))
    for case in payload["cases"]:
        if case["status"] != "ran":
            print(f"{case['case']}: {case['status']} {case['error']}")
    return 0 if payload["all_gates_passed"] else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())