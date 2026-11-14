"""Consume a CAD handoff the way a solver does: read the files, run the solver, gate the number.

Why this module exists
----------------------
:mod:`tensorcad.solver_interface` defines the handoff *contract*: a z-up frame,
outward normals, solver-facing tables on disk, and :func:`~tensorcad.solver_interface.verify_bundle`
which re-reads those tables with an independent parser and refuses a bundle whose
declared numbers do not survive the re-read.  What that module did not do is
**consume** the bundle: no code read ``surface.json`` back off disk and handed it
to a solver, so "exported" was still not "accepted", and a bundle could be
shipped with an in-memory object silently standing in for the file.

Here the files *are* the interface, and the claim is machine-checked:

1. :func:`tensorcad.solver_interface.verify_bundle` runs first and is a hard gate.
   A tampered table, a wrong digest, or a manifest whose frame is not ``z-up``
   raises instead of producing a number.
2. The consumer then parses the table itself (its own reader, not the writer's
   objects), rebuilds arrays, and only then calls the solver library.
3. A solver library that cannot be imported (or does not expose the expected API)
   is reported as ``status="unavailable"`` with the exception text.  Missing
   evidence is a status, never a zero and never a pass.
4. Each measured quantity is compared with a **closed-form** reference computed
   here, not supplied by the solver, and the relative error is written out.

Routes and their consumers (one end-to-end case each):

=================  ==========================  ======================================
bundle route       solver library              case / closed-form reference
=================  ==========================  ======================================
``surface_json``   TensorBEM                   submerged sphere, added mass
                                               ``2/3 pi rho a^3`` (Lamb)
``grid_json``      TensorFEM                   simply-supported plate under a sine
                                               load, ``q0 a^4/(4 pi^4 D)`` (Navier)
``lattice_json``   TensorLBM                   body-force channel, ``F H^2/(8 nu)``
=================  ==========================  ======================================

Cross-repository finding recorded here rather than hidden: the TensorLBM release
measured in this environment streams populations with the sign of its own lattice
direction table reversed (``stream`` rolls by ``-c`` instead of ``+c`` for a
``(9, ny, nx)`` table).  :func:`stream_probe` measures that and the result is
written into ``results.json``; the channel case therefore uses TensorLBM's
equilibrium, BGK collision and macroscopic moments with a local half-way
bounce-back streaming step -- a hybrid, stated as such, not a solver that "just
worked".
"""
from __future__ import annotations

import importlib
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from .solver_interface import (
    COORDINATE_CONTRACT,
    HandoffError,
    channel_grid_2d,
    ellipsoid_mesh,
    read_grid_table,
    read_lattice_table,
    signed_volume,
    verify_bundle,
    write_bundle_grid,
    write_bundle_lattice,
    write_bundle_surface,
)

__all__ = [
    "RESULT_SCHEMA",
    "SOLVER_MODULES",
    "CaseResult",
    "SolverUnavailableError",
    "bem_added_mass_from_bundle",
    "fem_plate_deflection_from_bundle",
    "gate",
    "lbm_channel_from_bundle",
    "run_all_cases",
    "solver_fingerprint",
    "stream_probe",
    "watertight_report",
]

#: The solver families the handoff contract is tested against.
SOLVER_MODULES = ("tensorbem", "tensorfem", "tensorlbm")

#: Bumped on any incompatible change to ``results.json``.
RESULT_SCHEMA = "tensorcad.solver_consumers/1"


class SolverUnavailableError(RuntimeError):
    """A solver library is absent or does not expose the expected API.

    This is a *status*, not a failure of the handoff: the bundle may be perfect
    and the host may simply have no TensorBEM installed.
    """


@dataclass(frozen=True)
class CaseResult:
    """One end-to-end consumption: what was read, what was measured, what passed."""

    case: str
    route: str
    solver: str
    status: str
    metrics: dict[str, Any]
    gates: tuple[dict[str, Any], ...] = ()
    solver_info: dict[str, Any] | None = None
    note: str = ""
    error: str = ""

    @property
    def gates_failed(self) -> tuple[str, ...]:
        return tuple(g["name"] for g in self.gates if not g["passed"])

    def to_dict(self) -> dict[str, Any]:
        return {
            "case": self.case,
            "route": self.route,
            "solver": self.solver,
            "status": self.status,
            "metrics": self.metrics,
            "gates": [dict(g) for g in self.gates],
            "gates_failed": list(self.gates_failed),
            "solver_info": dict(self.solver_info or {}),
            "note": self.note,
            "error": self.error,
        }


def _load(module_name: str) -> Any:
    try:
        return importlib.import_module(module_name)
    except Exception as exc:  # noqa: BLE001 - the message is the evidence
        raise SolverUnavailableError(f"{module_name}: {type(exc).__name__}: {exc}") from exc


def solver_fingerprint(module: Any) -> dict[str, Any]:
    """Record *which* copy of a solver library produced a number.

    Repository checkouts of the same library differ across agent worktrees, so a
    result naming only ``module.__version__`` (often ``None``) is not
    reproducible.  The git commit of the checkout the import resolved to is part
    of the evidence.
    """
    path = str(getattr(module, "__file__", "") or "")
    root = None
    for parent in Path(path).resolve().parents:
        if (parent / ".git").exists():
            root = parent
            break
    sha = None
    if root is not None:
        try:
            sha = (
                subprocess.run(
                    ["git", "-C", str(root), "rev-parse", "HEAD"],
                    capture_output=True,
                    text=True,
                    timeout=10,
                    check=False,
                ).stdout.strip()
                or None
            )
        except (OSError, subprocess.SubprocessError):  # pragma: no cover - host dependent
            sha = None
    return {
        "module": getattr(module, "__name__", "?"),
        "path": path,
        "version": getattr(module, "__version__", None),
        "repo_root": str(root) if root is not None else None,
        "git_sha": sha,
    }


def gate(name: str, measured: float, reference: float, tolerance: float, unit: str = "") -> dict[str, Any]:
    """A ref-vs-got comparison whose error is recomputed, not asserted."""
    measured_f = float(measured)
    reference_f = float(reference)
    rel = (measured_f - reference_f) / reference_f if reference_f else float("nan")
    return {
        "name": name,
        "measured": measured_f,
        "reference": reference_f,
        "rel_error": rel,
        "tolerance": tolerance,
        "unit": unit,
        "passed": bool(abs(rel) <= tolerance),
    }


def _numpy() -> Any:
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover - exercised only without NumPy
        raise SolverUnavailableError("the consumers need NumPy for their tables") from exc
    return np


# --------------------------------------------------------------- mesh integrity


def watertight_report(
    vertices: Sequence[Sequence[float]], faces: Sequence[Sequence[int]]
) -> dict[str, Any]:
    """Edge census of a triangle surface: closed, consistently oriented, outward.

    A BEM panel code silently produces nonsense on a surface with a hole or a
    flipped panel, so the consumer measures all three instead of trusting the
    exporter.  Every undirected edge must appear exactly twice, and every
    *directed* edge exactly once -- that is what makes the two uses opposite.
    """
    directed: dict[tuple[int, int], int] = {}
    for a, b, c in faces:
        for u, v in ((a, b), (b, c), (c, a)):
            directed[(u, v)] = directed.get((u, v), 0) + 1
    undirected: dict[tuple[int, int], int] = {}
    for (u, v), count in directed.items():
        key = (u, v) if u < v else (v, u)
        undirected[key] = undirected.get(key, 0) + count
    bad_edges = [key for key, count in undirected.items() if count != 2]
    duplicated = [key for key, count in directed.items() if count > 1]
    volume = signed_volume(vertices, faces)
    return {
        "edges": len(undirected),
        "bad_edges": len(bad_edges),
        "duplicate_directed_edges": len(duplicated),
        "closed": not bad_edges,
        "orientation_consistent": not duplicated,
        "signed_volume_m3": volume,
        "outward": volume > 0.0,
    }


def _read_surface_table(path: Path) -> tuple[list[tuple[float, float, float]], list[tuple[int, int, int]]]:
    table = json.loads(path.read_text(encoding="utf-8"))
    vertices = [tuple(float(c) for c in v) for v in table["vertices"]]
    faces = [tuple(int(i) for i in f) for f in table["faces"]]
    if any(len(face) != 3 for face in faces):
        raise HandoffError("a BEM surface table must hold triangles")
    if any(not (0 <= i < len(vertices)) for face in faces for i in face):
        raise HandoffError("a surface triangle references a node outside the table")
    return vertices, faces


def _verified_manifest(bundle_dir: str | Path, *, route: str, table: str) -> dict[str, Any]:
    """Run the contract verifier, then insist on the route the caller needs."""
    manifest = verify_bundle(bundle_dir)
    if manifest.get("route") != route:
        raise HandoffError(
            f"bundle at {bundle_dir} re-read as route {manifest.get('route')!r}, not {route!r}"
        )
    if manifest.get("frame") != COORDINATE_CONTRACT["name"]:
        raise HandoffError(f"bundle frame {manifest.get('frame')!r} is not the z-up contract")
    declared = json.loads((Path(bundle_dir) / "handoff.json").read_text(encoding="utf-8"))["tables"]
    if table not in declared:
        raise HandoffError(f"bundle does not name {table!r} among its tables")
    return manifest


# ------------------------------------------------------------------- consumers


def bem_added_mass_from_bundle(
    bundle_dir: str | Path,
    *,
    density: float = 1000.0,
    center: Sequence[float] = (0.0, 0.0, 0.0),
    quadrature_order: int = 6,
) -> CaseResult:
    """Read ``surface.json`` off disk and hand the mesh to :mod:`tensorbem`.

    The bundle must re-read as the ``surface`` route.  The triangle table is
    parsed here, checked for closure and outward orientation, and only then given
    to ``TriangleMesh`` -- a hole or a flipped panel is refused before the solver
    sees it.
    """
    bundle = Path(bundle_dir)
    manifest = _verified_manifest(bundle, route="surface", table="surface_json")
    vertices, faces = _read_surface_table(bundle / "surface.json")
    integrity = watertight_report(vertices, faces)
    if not integrity["closed"]:
        raise HandoffError(f"surface has {integrity['bad_edges']} non-manifold edges")
    if not integrity["orientation_consistent"]:
        raise HandoffError("surface has duplicated directed edges: panels are inconsistent")
    if not integrity["outward"]:
        raise HandoffError("surface normals point inward (negative signed volume)")

    np = _numpy()
    module = _load("tensorbem")
    info = solver_fingerprint(module)
    mesh = module.TriangleMesh(np.asarray(vertices, dtype=float), np.asarray(faces, dtype=int))
    result = module.compute_added_mass_matrix_3d(
        mesh,
        density=density,
        center=tuple(float(c) for c in center),
        quadrature_order=quadrature_order,
    )
    raw = getattr(result, "added_mass_matrix", None)
    if raw is None:
        raw = getattr(result, "added_mass", None)
    if raw is None:
        raw = getattr(result, "matrix", None)
    if raw is None:  # pragma: no cover - API guard
        raise SolverUnavailableError(
            "tensorbem result exposes no added-mass matrix (looked for "
            ".added_mass_matrix, .added_mass, .matrix)"
        )
    matrix = np.asarray(raw, dtype=float)
    metrics: dict[str, Any] = {
        "nodes": len(vertices),
        "panels": len(faces),
        "signed_volume_m3": integrity["signed_volume_m3"],
        "edges": integrity["edges"],
        "bundle_digest": manifest.get("digest_sha256"),
        "added_mass_matrix": [[float(v) for v in row] for row in matrix],
    }
    for label, index in (("surge_kg", 0), ("sway_kg", 1), ("heave_kg", 2)):
        metrics[label] = float(matrix[index][index])
    for label in ("max_relative_residual_l2", "symmetry_error", "double_layer_row_sum_error", "condition_number"):
        value = getattr(result, label, None)
        if value is not None:
            metrics[label] = float(value)
    return CaseResult(
        case=str(manifest.get("case", "")),
        route="surface_json",
        solver="tensorbem",
        status="ran",
        metrics=metrics,
        solver_info=info,
        note="panels parsed from surface.json after verify_bundle(); no in-memory object is used",
    )


def _plate_mesh_from_table(bundle_dir: Path) -> tuple[Any, Any, dict[str, Any]]:
    """``grid.json`` -> torch nodes/connectivity, winding and manifoldness checked."""
    grid = read_grid_table(bundle_dir / "grid.json")
    vertices = grid["vertices"]
    quads = grid["quads"]
    if any(abs(v[2]) > 1e-12 for v in vertices):
        raise HandoffError("a plate mid-surface must lie on z = 0 in the contract frame")
    directed: dict[tuple[int, int], int] = {}
    for quad in quads:
        for u, v in zip(quad, quad[1:] + quad[:1]):
            directed[(u, v)] = directed.get((u, v), 0) + 1
    if any(count != 1 for count in directed.values()):
        raise HandoffError("grid quads share a directed edge: connectivity is not a manifold mesh")
    area = 0.0
    for quad in quads:
        a, b, c, d = (vertices[i] for i in quad)
        area += 0.5 * abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))
        area += 0.5 * abs((c[0] - a[0]) * (d[1] - a[1]) - (d[0] - a[0]) * (c[1] - a[1]))
    import torch

    nodes = torch.tensor([[v[0], v[1]] for v in vertices], dtype=torch.float64)
    elems = torch.tensor([list(q) for q in quads], dtype=torch.long)
    return nodes, elems, {"nodes": len(vertices), "quads": len(quads), "area_m2": area}


def fem_plate_deflection_from_bundle(
    bundle_dir: str | Path,
    *,
    side: float,
    young: float = 7.0e10,
    poisson: float = 0.3,
    thickness: float = 0.05,
    q0: float = 1.0e4,
) -> CaseResult:
    """Read ``grid.json`` off disk and assemble a Mindlin Q4 plate in :mod:`tensorfem`.

    The element stiffness, the shear treatment and the simply-supported
    convention come from ``tensorfem.plate``; the *mesh* comes from the file.
    That is the point of the route: the CAD side decides the discretisation and
    the solver does not re-generate it from its own parameters.
    """
    bundle = Path(bundle_dir)
    manifest = _verified_manifest(bundle, route="grid", table="grid_json")
    nodes, elems, grid_info = _plate_mesh_from_table(bundle)

    import torch

    plate = _load("tensorfem.plate")
    info = solver_fingerprint(plate)
    ndof = 3 * nodes.shape[0]
    stiffness = torch.zeros((ndof, ndof), dtype=torch.float64)
    load = torch.zeros(ndof, dtype=torch.float64)
    gp = 1.0 / 3.0**0.5
    for conn in elems:
        xy = nodes[conn]
        ke = plate.q4_mindlin_stiffness(xy, young, poisson, thickness)
        fe = torch.zeros(12, dtype=torch.float64)
        for xif in (-gp, gp):
            for etaf in (-gp, gp):
                xi = torch.tensor(xif, dtype=torch.float64)
                eta = torch.tensor(etaf, dtype=torch.float64)
                shape, dshape = plate._shape(xi, eta)
                jac = dshape.T @ xy
                point = shape @ xy
                q = q0 * torch.sin(torch.pi * point[0] / side) * torch.sin(torch.pi * point[1] / side)
                fe[0::3] += shape * q * torch.linalg.det(jac)
        dofs = torch.stack((3 * conn, 3 * conn + 1, 3 * conn + 2), dim=1).reshape(-1)
        stiffness[dofs[:, None], dofs] += ke
        load[dofs] += fe

    tol = side * 1e-12
    xedge = (nodes[:, 0].abs() < tol) | ((nodes[:, 0] - side).abs() < tol)
    yedge = (nodes[:, 1].abs() < tol) | ((nodes[:, 1] - side).abs() < tol)
    edge = xedge | yedge
    fixed = torch.cat(
        (
            3 * torch.nonzero(edge, as_tuple=False).flatten(),
            3 * torch.nonzero(yedge, as_tuple=False).flatten() + 1,
            3 * torch.nonzero(xedge, as_tuple=False).flatten() + 2,
        )
    ).unique()
    all_dofs = torch.arange(ndof)
    free = all_dofs[~torch.isin(all_dofs, fixed)]
    solution = torch.zeros(ndof, dtype=torch.float64)
    solution[free] = torch.linalg.solve(stiffness[free[:, None], free], load[free])
    residual = stiffness @ solution - load
    w = solution[0::3]
    centre_index = int(torch.argmin((nodes[:, 0] - side / 2) ** 2 + (nodes[:, 1] - side / 2) ** 2))
    metrics: dict[str, Any] = {
        **grid_info,
        "centre_deflection_m": float(w[centre_index]),
        "max_deflection_m": float(w.abs().max()),
        "free_dof_residual_norm": float(torch.linalg.norm(residual[free])),
        "load_norm": float(torch.linalg.norm(load)),
        "bundle_digest": manifest.get("digest_sha256"),
    }
    return CaseResult(
        case=str(manifest.get("case", "")),
        route="grid_json",
        solver="tensorfem.plate",
        status="ran",
        metrics=metrics,
        solver_info=info,
        note="mesh from grid.json; element stiffness and SS convention from tensorfem.plate",
    )


def lbm_channel_from_bundle(
    bundle_dir: str | Path,
    *,
    tau: float = 0.8,
    body_force: float = 2.0e-5,
    max_steps: int = 20000,
    tol: float = 1.0e-12,
) -> CaseResult:
    """Read ``lattice.json`` off disk and run a body-force channel with :mod:`tensorlbm`.

    The occupancy table decides which cells are fluid; the solver never sees the
    CAD representation.  TensorLBM supplies the equilibrium, the BGK collision and
    the macroscopic moments; streaming is local, with half-way bounce-back on the
    ``k = 0`` and ``k = nk - 1`` cell faces, because the measured release of
    ``stream`` transports populations against its own lattice directions (see
    :func:`stream_probe`).
    """
    bundle = Path(bundle_dir)
    manifest = _verified_manifest(bundle, route="lattice", table="lattice_json")
    data = read_lattice_table(bundle / "lattice.json")
    shape = data["shape"]
    spacing = data["spacing"]
    cells = data["fluid_cells"]
    if shape[1] != 1:
        raise HandoffError("the channel case is the x-z plane (one cell through y)")
    expected = shape[0] * shape[2]
    if len(cells) != expected:
        raise HandoffError(f"the channel lattice expects every cell to be fluid, got {len(cells)} of {expected}")
    walls = data.get("walls", {})
    if walls.get("normal_axis") != "z":
        raise HandoffError("the channel lattice must declare walls normal to z")

    import torch

    lbm = _load("tensorlbm")
    info = solver_fingerprint(lbm)
    nz, nx = shape[2], shape[0]
    c = torch.as_tensor(lbm.d2q9.C, dtype=torch.long)
    weights = torch.as_tensor(lbm.d2q9.W, dtype=torch.float64)
    opposite = {i: int(torch.nonzero((c == -c[i]).all(dim=1)).flatten()[0]) for i in range(9)}

    rho0 = torch.ones((nz, nx), dtype=torch.float64)
    zeros = torch.zeros((nz, nx), dtype=torch.float64)
    f = lbm.equilibrium(rho0, zeros, zeros)

    force = torch.zeros(2, dtype=torch.float64)
    force[0] = body_force

    def forcing(u_x: Any, u_y: Any) -> Any:
        """Guo source for a constant body force: S_i = w_i [3 c.F + 9 (c.u)(c.F) - 3 u.F].

        The plain source (no ``1 - 1/(2 tau)`` split) is deliberate: with the split
        factor and the half-force velocity the measured steady amplitude came out at
        ``(1 - 1/(2 tau)) F`` instead of ``F`` -- a 62.5 % shortfall at ``tau = 0.8``
        that the closed-form gate catches.  The residual of this form is the known
        constant ``O(F)`` offset, reported as ``forcing_offset_lattice`` rather than
        hidden.
        """
        cx = c[:, 0].to(torch.float64).view(9, 1, 1)
        cz = c[:, 1].to(torch.float64).view(9, 1, 1)
        cu = cx * u_x + cz * u_y
        cf = cx * force[0] + cz * force[1]
        uf = u_x * force[0] + u_y * force[1]
        w = weights.view(9, 1, 1)
        return w * (3.0 * cf + 9.0 * cu * cf - 3.0 * uf)

    def stream(f_post: Any) -> Any:
        out = torch.empty_like(f_post)
        for i in range(9):
            cx, cz = int(c[i, 0]), int(c[i, 1])
            out[i] = torch.roll(f_post[i], shifts=(cz, cx), dims=(0, 1))
        for i in range(9):
            cz = int(c[i, 1])
            j = opposite[i]
            if cz > 0:
                out[i][0, :] = f_post[j][0, :]
            elif cz < 0:
                out[i][nz - 1, :] = f_post[j][nz - 1, :]
        return out

    steps = 0
    change = float("inf")
    while steps < max_steps:
        rho, ux, uy = lbm.macroscopic(f)
        f_post = lbm.collide_bgk(f, tau) + forcing(ux, uy)
        f_next = stream(f_post)
        change = float((f_next - f).abs().max())
        f = f_next
        steps += 1
        if change < tol:
            break

    rho, ux, uy = lbm.macroscopic(f)
    column = ux[:, nx // 2]
    centre_index = int(torch.argmax(column).item())
    metrics: dict[str, Any] = {
        "shape": list(shape),
        "dx_m": spacing[0],
        "dz_m": spacing[2],
        "height_m": spacing[2] * nz,
        "fluid_cells": len(cells),
        "steps": steps,
        "steady_change": change,
        "u_max_lattice": float(column.max()),
        "u_max_z_m": (centre_index + 0.5) * spacing[2],
        "centreline_profile_lattice": [float(v) for v in column.tolist()],
        "mass_drift": float(abs(rho.sum() - rho.numel())),
        "curvature_lattice": float(
            max(
                abs(column[k + 1] - 2.0 * column[k] + column[k - 1])
                for k in range(1, nz - 1)
            )
        ),
        "bundle_digest": manifest.get("digest_sha256"),
    }
    return CaseResult(
        case=str(manifest.get("case", "")),
        route="lattice_json",
        solver="tensorlbm",
        status="ran",
        metrics=metrics,
        solver_info=info,
        note="occupancy from lattice.json; TensorLBM equilibrium/BGK/moments + local half-way bounce-back",
    )


def stream_probe() -> dict[str, Any]:
    """Measure which way ``tensorlbm.stream`` actually transports a population.

    A unit population on the ``+x`` direction is placed at one node; after one
    streaming step it must appear one node further in ``+x`` for a
    ``(9, ny, nx)`` table.  The reported ``dx_lattice`` is measured, not claimed.
    """
    try:
        import torch

        lbm = _load("tensorlbm")
        nx = ny = 6
        f = torch.zeros((9, ny, nx), dtype=torch.float64)
        f[1, 2, 2] = 1.0  # population 1 -> +x
        moved = lbm.stream(f)
        index = torch.nonzero(moved[1] > 0.5)
        if index.numel() == 0:
            return {"status": "indeterminate", "dx_lattice": None}
        dx = int(index[0][1].item()) - 2
        return {
            "status": "measured",
            "dx_lattice": dx,
            "expected_dx_lattice": 1,
            "agrees_with_sign_convention": dx == 1,
            "note": "stream(f) must move a +x population to x+1 for a (9, ny, nx) table",
        }
    except Exception as exc:  # noqa: BLE001 - a probe must not fail the run
        return {"status": "unavailable", "error": f"{type(exc).__name__}: {exc}"}


# ------------------------------------------------------------------ case driver


def _ellipsoid_bundle(
    outdir: Path, *, semi_axes: Sequence[float], n_theta: int, n_phi: int, case: str
) -> dict[str, Any]:
    mesh = ellipsoid_mesh(semi_axes, n_theta=n_theta, n_phi=n_phi)
    return write_bundle_surface(
        outdir,
        mesh=mesh,
        case=case,
        extra={
            "generator": "tensorcad.solver_interface.ellipsoid_mesh",
            "n_theta": n_theta,
            "n_phi": n_phi,
            "aging": "surface handed over as a written triangle table",
        },
    )


def _plate_bundle(outdir: Path, *, side: float, n: int, case: str) -> dict[str, Any]:
    grid = channel_grid_2d(length=side, height=side, nx=n, ny=n)
    xs = [point[0] for point in grid[0]]
    ys = [row[0][1] for row in grid]
    return write_bundle_grid(
        outdir,
        xs=xs,
        ys=ys,
        case=case,
        extra={"generator": "tensorcad.solver_interface.channel_grid_2d"},
    )


def _channel_bundle(
    outdir: Path, *, length: float, height: float, nx: int, nz: int, case: str
) -> dict[str, Any]:
    dx, dz = length / nx, height / nz
    cells = [(i, 0, k) for i in range(nx) for k in range(nz)]
    return write_bundle_lattice(
        outdir,
        shape=(nx, 1, nz),
        origin=(0.0, 0.0, 0.0),
        spacing=(dx, length, dz),
        fluid_cells=cells,
        walls={
            "normal_axis": "z",
            "lower_k": 0,
            "upper_k": nz - 1,
            "kind": "no-slip half-way bounce-back",
        },
        case=case,
        extra={
            "generator": "tensorcad.solver_interface.channel_grid_2d",
            "note": "x-z plane, one cell through y",
        },
    )


def _voxelisation_gate() -> dict[str, Any]:
    """CAD -> voxel occupancy: a curved solid must converge as O(dx^2)."""
    from .geometry_field import EllipsoidField
    from .voxel import to_voxel_grid

    semi = (0.5, 0.35, 0.25)
    analytic = 4.0 / 3.0 * math.pi * semi[0] * semi[1] * semi[2]
    field = EllipsoidField(center=(0.0, 0.0, 0.0), semi_axes=semi)
    out: dict[str, Any] = {"analytic_volume_m3": analytic, "semi_axes_m": list(semi)}
    for resolution in (16, 32):
        grid = to_voxel_grid(field, resolution=resolution)
        volume = float(grid.occupied_volume())
        out[f"voxel_volume_{resolution}"] = volume
        out[f"voxel_rel_error_{resolution}"] = (volume - analytic) / analytic
        out[f"voxel_count_{resolution}"] = int(grid.inside_count())
    return out


def _attach(
    result: CaseResult,
    gates: Sequence[dict[str, Any]],
    *,
    extra: dict[str, Any],
) -> CaseResult:
    """Attach gates and driver-side reference bookkeeping to a measurement."""
    metrics = dict(result.metrics)
    metrics.update(extra)
    return CaseResult(
        case=result.case,
        route=result.route,
        solver=result.solver,
        status=result.status,
        metrics=metrics,
        gates=tuple(gates),
        solver_info=result.solver_info,
        note=result.note,
    )


def _unavailable(case: str, route: str, solver: str, exc: Exception) -> CaseResult:
    return CaseResult(case=case, route=route, solver=solver, status="unavailable", metrics={}, error=str(exc))


def _rejected(case: str, route: str, solver: str, exc: Exception) -> CaseResult:
    return CaseResult(case=case, route=route, solver=solver, status="rejected", metrics={}, error=str(exc))


def run_all_cases(
    outdir: str | Path,
    *,
    bem_resolutions: Sequence[tuple[int, int]] = ((16, 32), (24, 48)),
    fem_resolutions: Sequence[int] = (8, 16),
    lbm_resolutions: Sequence[int] = (8, 16),
    sphere_radius: float = 0.5,
    density: float = 1000.0,
    plate_side: float = 2.0,
    plate_thickness: float = 0.05,
    young: float = 7.0e10,
    poisson: float = 0.3,
    q0: float = 1.0e4,
    tau: float = 0.8,
    body_force: float = 2.0e-5,
    max_steps: int = 20000,
) -> dict[str, Any]:
    """Write bundles from the geometry library, consume each one, gate the numbers."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    cases: list[CaseResult] = []

    # --- BEM: submerged sphere, added mass reference 2/3 pi rho a^3 (Lamb)
    reference_added_mass = 2.0 / 3.0 * math.pi * density * sphere_radius**3
    analytic_sphere_volume = 4.0 / 3.0 * math.pi * sphere_radius**3
    for n_theta, n_phi in bem_resolutions:
        sub = out / f"bem_{n_theta}x{n_phi}"
        try:
            manifest = _ellipsoid_bundle(
                sub,
                semi_axes=(sphere_radius, sphere_radius, sphere_radius),
                n_theta=n_theta,
                n_phi=n_phi,
                case=f"sphere_added_mass_{n_theta}x{n_phi}",
            )
            result = bem_added_mass_from_bundle(
                sub, density=density, center=(0.0, 0.0, -0.75), quadrature_order=6
            )
            gates = (
                gate("enclosed_volume", result.metrics["signed_volume_m3"], analytic_sphere_volume, 0.02, "m3"),
                gate("added_mass_heave", result.metrics["heave_kg"], reference_added_mass, 0.05, "kg"),
                gate("added_mass_surge", result.metrics["surge_kg"], reference_added_mass, 0.05, "kg"),
            )
            cases.append(_attach(result, gates, extra={"declared_volume_m3": manifest["signed_volume_m3"]}))
        except SolverUnavailableError as exc:
            cases.append(_unavailable(f"sphere_added_mass_{n_theta}x{n_phi}", "surface_json", "tensorbem", exc))
        except HandoffError as exc:
            cases.append(_rejected(f"sphere_added_mass_{n_theta}x{n_phi}", "surface_json", "tensorbem", exc))

    # --- FEM: simply supported plate, sine load, reference q0 a^4 / (4 pi^4 D)
    rigidity = young * plate_thickness**3 / (12.0 * (1.0 - poisson**2))
    reference_deflection = q0 * plate_side**4 / (4.0 * math.pi**4 * rigidity)
    for n in fem_resolutions:
        sub = out / f"fem_{n}x{n}"
        try:
            _plate_bundle(sub, side=plate_side, n=n, case=f"plate_sine_{n}x{n}")
            result = fem_plate_deflection_from_bundle(
                sub,
                side=plate_side,
                young=young,
                poisson=poisson,
                thickness=plate_thickness,
                q0=q0,
            )
            gates = (
                gate("centre_deflection", result.metrics["centre_deflection_m"], reference_deflection, 0.05, "m"),
                gate("plate_area", result.metrics["area_m2"], plate_side**2, 1e-12, "m2"),
            )
            cases.append(_attach(result, gates, extra={"reference_deflection_m": reference_deflection, "rigidity_Nm": rigidity}))
        except SolverUnavailableError as exc:
            cases.append(_unavailable(f"plate_sine_{n}x{n}", "grid_json", "tensorfem", exc))
        except HandoffError as exc:
            cases.append(_rejected(f"plate_sine_{n}x{n}", "grid_json", "tensorfem", exc))

    # --- LBM: body-force channel, reference F H^2 / (8 nu)
    viscosity = (tau - 0.5) / 3.0
    for nz in lbm_resolutions:
        sub = out / f"lbm_{nz}"
        height = 0.4
        try:
            _channel_bundle(sub, length=0.8, height=height, nx=2 * nz, nz=nz, case=f"channel_poiseuille_{nz}")
            result = lbm_channel_from_bundle(sub, tau=tau, body_force=body_force, max_steps=max_steps)
            reference_umax = body_force * nz**2 / (8.0 * viscosity)
            # The plain Guo source leaves a constant O(F) offset, so the *shape* is
            # asserted through the discrete curvature (a parabola has constant second
            # difference -F/nu) instead of through an absolute profile comparison.
            curvature_ref = body_force / viscosity
            shape_error = abs(result.metrics["curvature_lattice"] - curvature_ref) / curvature_ref
            gates = (
                gate("peak_velocity_lattice", result.metrics["u_max_lattice"], reference_umax, 0.05, "lattice"),
                gate("centreline_curvature", result.metrics["curvature_lattice"], curvature_ref, 0.01, "lattice"),
                {
                    "name": "mass_conservation",
                    "measured": result.metrics["mass_drift"],
                    "reference": 0.0,
                    "rel_error": 0.0,
                    "tolerance": 1e-9,
                    "unit": "lattice",
                    "passed": bool(result.metrics["mass_drift"] < 1e-9),
                },
            )
            cases.append(
                _attach(
                    result,
                    gates,
                    extra={
                        "reference_u_max_lattice": reference_umax,
                        "height_m": height,
                        "analytic_umax_SI": body_force * height**2 / (8.0 * viscosity),
                        "forcing_offset_lattice": result.metrics["u_max_lattice"] - reference_umax,
                        "shape_error_rel": shape_error,
                    },
                )
            )
        except SolverUnavailableError as exc:
            cases.append(_unavailable(f"channel_poiseuille_{nz}", "lattice_json", "tensorlbm", exc))
        except HandoffError as exc:
            cases.append(_rejected(f"channel_poiseuille_{nz}", "lattice_json", "tensorlbm", exc))

    results = {
        "schema": RESULT_SCHEMA,
        "frame": dict(COORDINATE_CONTRACT),
        "cases": [case.to_dict() for case in cases],
        "voxelisation": _voxelisation_gate(),
        "stream_probe": stream_probe(),
        "gates_failed": [f"{case.case}:{name}" for case in cases for name in case.gates_failed],
    }
    results["statuses"] = {case.case: case.status for case in cases}
    (out / "results.json").write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return results


def _main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Consume CAD bundles with the three solver families")
    parser.add_argument("--out", default="artifacts/solver_consumers", help="output directory")
    parser.add_argument("--quick", action="store_true", help="coarse resolutions only")
    args = parser.parse_args(argv)
    if args.quick:
        results = run_all_cases(
            args.out,
            bem_resolutions=((12, 24),),
            fem_resolutions=(8,),
            lbm_resolutions=(8,),
            max_steps=4000,
        )
    else:
        results = run_all_cases(args.out)
    print(json.dumps({"statuses": results["statuses"], "gates_failed": results["gates_failed"]}, indent=2))
    return 0 if not results["gates_failed"] else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
