"""Three solvers, one geometry source, one frame contract -- end to end.

Each stage is the same: ``tensorcad`` writes a bundle (tables + declared
numbers), :func:`tensorcad.solver_interface.verify_bundle` re-reads it with an
independent parser, and only then does a solver package consume the files.

1. ``fem``  -- box solid -> ``C3D4`` gmsh -> TensorFEM ``Tet4J2Model`` in the
   elastic range -> uniaxial tension patch test against the exact linear field.
2. ``bem``  -- CAD ellipsoid surface -> triangle table -> TensorBEM
   ``compute_added_mass_matrix_3d`` -> analytic sphere added mass ``0.5 rho V``.
3. ``fvm``  -- CAD planar channel grid -> TensorFVM ``Mesh2D`` backend ->
   collocated SIMPLE -> plane-Poiseuille profile ``6 ubar y/h (1 - y/h)``.

Run with::

    PYTHONPATH=.:/root/tensor-ecosystem/TensorFVM/src python \
        examples/solver_interface/run_three_solvers.py [--only fem,bem,fvm]

Numbers land in ``artifacts/solver_interfaces/summary.json``; every entry
records the bundle it came from so a reader can re-derive it.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tensorcad import solver_interface as si  # noqa: E402
from tensorcad.volume_mesh import box_section_grid, build_volume_mesh  # noqa: E402

ART = ROOT / "artifacts" / "solver_interfaces"


def _write_result(case: str, payload: dict) -> Path:
    out = ART / case
    out.mkdir(parents=True, exist_ok=True)
    path = out / "results.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return path


# --------------------------------------------------------------------------- FEM


def run_fem_case(length=2.0, beam=0.4, draft=0.3, *, nx=9, nz=5, layers=2,
                 young=2.1e11, poisson=0.3, stress=1.0e6) -> dict:
    import torch
    from tensorfem.solid_plasticity import Tet4J2Model, PlasticSolidState, solve_increment

    started = time.time()
    grid = box_section_grid(length=length, beam=beam, draft=draft, nx=nx, nz=nz)
    volume = build_volume_mesh(grid, first_layer=0.02 * draft, growth=1.2, layers=layers,
                               split_tets=True)
    case_dir = ART / "fem_tension"
    manifest = si.write_bundle_volume(
        case_dir,
        vertices=volume.vertices,
        tets=volume.tets,
        wall_name="wall",
        case=f"box-uniaxial-tension nx={nx} nz={nz} layers={layers}",
        extra={
            "solver": "TensorFEM Tet4J2Model (elastic range, yield_stress=1e9*stress)",
            "geometry": {"length_m": length, "beam_m": beam, "draft_m": draft},
            "mesh_route": "tensorcad.volume_mesh.build_volume_mesh(split_tets=True)",
        },
    )
    check = si.verify_bundle(case_dir)

    data = si.read_gmsh_volume(case_dir / "volume.msh")
    nodes = np.asarray(data["nodes"], dtype=np.float64)
    tets = np.asarray(data["tets"], dtype=np.int64)
    nnodes = len(nodes)

    area = beam * draft  # x = length face
    x_min, x_max = nodes[:, 0].min(), nodes[:, 0].max()
    face = np.isclose(nodes[:, 0], x_max, atol=1e-12)
    external = np.zeros(3 * nnodes)
    # Consistent nodal loads for a uniform traction: each loaded boundary triangle
    # carries stress*area and splits it equally over its three corners. Equal force
    # per node would be inconsistent on a graded face and would pollute the patch test.
    loaded = 0
    for a, b, c in si.boundary_triangles(tets):
        if face[a] and face[b] and face[c]:
            pa, pb, pc = nodes[a], nodes[b], nodes[c]
            tri_area = 0.5 * float(np.linalg.norm(np.cross(pb - pa, pc - pa)))
            share = stress * tri_area / 3.0
            for node in (a, b, c):
                external[3 * node] += share
            loaded += 1

    fixed = []
    base = np.isclose(nodes[:, 0], x_min, atol=1e-12)
    for node in np.nonzero(base)[0]:
        fixed.append(3 * int(node) + 0)  # u_x = 0 on the clamped end
    pin = int(np.nonzero(base)[0][0])
    fixed += [3 * pin + 1, 3 * pin + 2]  # remove rigid translation in y, z

    model = Tet4J2Model(
        nodes=torch.tensor(nodes, dtype=torch.float64),
        elements=torch.tensor(tets, dtype=torch.long),
        young=float(young),
        poisson=float(poisson),
        yield_stress=1.0e9 * stress,  # elastic: no element can yield
        hardening=0.0,
        fixed_dofs=torch.tensor(sorted(set(fixed)), dtype=torch.long),
    )
    state = PlasticSolidState.virgin(model)
    solution = solve_increment(
        model,
        torch.tensor(external, dtype=torch.float64),
        1.0,
        state,
        torch.zeros(model.n_dofs, dtype=torch.float64),
        tolerance=1e-10,
        max_iterations=25,
    )
    u = solution.displacement.numpy().reshape(-1, 3)

    exact = np.stack(
        (
            stress * nodes[:, 0] / young,
            -poisson * stress * nodes[:, 1] / young,
            -poisson * stress * nodes[:, 2] / young,
        ),
        axis=1,
    )
    offset = u[pin] - exact[pin]  # the pinned node absorbs the free translation
    residual = (u - exact) - offset
    scale = np.abs(exact - exact[pin]).max()
    linf = float(np.abs(residual).max() / scale)
    internal = np.asarray(solution.reaction.detach().cpu().numpy(), dtype=float).reshape(-1, 3)
    reaction = float(np.abs(internal[base, 0]).sum())
    applied = float(external.reshape(-1, 3)[:, 0].sum())

    payload = {
        "case": "fem",
        "solver": "TensorFEM",
        "bundle": str((case_dir / "handoff.json").relative_to(ROOT)),
        "verify": check,
        "parameters": {
            "young_pa": young,
            "poisson": poisson,
            "applied_stress_pa": stress,
            "length_m": length,
            "beam_m": beam,
            "draft_m": draft,
            "nx": nx,
            "nz": nz,
            "layers": layers,
        },
        "declared_volume_m3": manifest["volume_m3"],
        "analytic_volume_m3": length * beam * draft,
        "loaded_triangles": loaded,
        "external_tension_n": float(external.reshape(-1, 3)[:, 0].sum()),
        "reaction_at_clamped_end_n": reaction,
        "equilibrium_relative_error": abs(reaction - applied) / applied,
        "patch_test_relative_linf": linf,
        "newton_iterations": int(solution.iterations),
        "newton_residual": float(solution.residual_norm),
        "wall_time_s": time.time() - started,
    }
    _write_result("fem_tension", payload)
    return payload


# --------------------------------------------------------------------------- BEM


def run_bem_case(*, radius=1.0, density=1000.0, levels=((10, 16), (16, 24), (22, 32))) -> dict:
    from tensorbem.geometry3d import TriangleMesh
    from tensorbem.solver3d import compute_added_mass_matrix_3d

    started = time.time()
    rows = []
    for n_theta, n_phi in levels:
        mesh = si.ellipsoid_mesh((radius, radius, radius), n_theta=n_theta, n_phi=n_phi)
        case_dir = ART / "bem_sphere"
        manifest = si.write_bundle_surface(
            case_dir, mesh,
            case=f"sphere added mass n_theta={n_theta} n_phi={n_phi}",
            extra={"solver": "TensorBEM compute_added_mass_matrix_3d",
                   "semi_axes_m": [radius, radius, radius]},
        )
        check = si.verify_bundle(case_dir)
        table = json.loads((case_dir / "surface.json").read_text())
        panel_mesh = TriangleMesh(
            np.asarray(table["vertices"], dtype=float), np.asarray(table["faces"], dtype=np.int64)
        )
        result = compute_added_mass_matrix_3d(panel_mesh, density=density, quadrature_order=3)
        matrix = np.asarray(result.added_mass_matrix, dtype=float)
        volume = 4.0 / 3.0 * math.pi * radius ** 3
        analytic = 0.5 * density * volume  # sphere: A_ii = 0.5 rho V, i = 1,2,3
        diagonals = [float(matrix[i, i]) for i in range(3)]
        rows.append({
            "n_theta": n_theta,
            "n_phi": n_phi,
            "panels": int(panel_mesh.panel_count),
            "digest_verified": check["digest_matches"],
            "discretised_volume_m3": check["signed_volume_m3"],
            "analytic_volume_m3": volume,
            "added_mass_surge_kg": diagonals[0],
            "added_mass_sway_kg": diagonals[1],
            "added_mass_heave_kg": diagonals[2],
            "analytic_added_mass_kg": analytic,
            "relative_error_surge": abs(diagonals[0] - analytic) / analytic,
            "relative_error_heave": abs(diagonals[2] - analytic) / analytic,
            "condition_number": float(result.condition_number),
            "max_relative_residual_l2": float(result.max_relative_residual_l2),
            "double_layer_row_sum_error": float(result.double_layer_row_sum_error),
        })
    payload = {
        "case": "bem",
        "solver": "TensorBEM",
        "bundle": str((ART / "bem_sphere" / "handoff.json").relative_to(ROOT)),
        "density_kg_m3": density,
        "radius_m": radius,
        "levels": rows,
        "wall_time_s": time.time() - started,
    }
    _write_result("bem_sphere", payload)
    return payload


# --------------------------------------------------------------------------- FVM


def _channel_mesh(config):
    """Build a TensorFVM ``Mesh2D`` from the CAD channel grid file."""
    import torch

    bundle = json.loads(Path(config.mesh_file).read_text())
    nodes = bundle["vertices"]  # [(x, y), ...] row-major by y
    nx, ny = bundle["nx"], bundle["ny"]
    opts = dict(dtype=torch.float64, device=config.device)
    x = torch.tensor([point[0] for point in nodes[0]], **opts)
    y = torch.tensor([row[0][1] for row in nodes], **opts)
    grid_y, grid_x = torch.meshgrid(y, x, indexing="ij")
    vertices = torch.stack((grid_x, grid_y), -1)
    polygon = torch.stack((vertices[:-1, :-1], vertices[:-1, 1:],
                           vertices[1:, 1:], vertices[1:, :-1]), -2)
    nxt = polygon.roll(-1, -2)
    cross = polygon[..., 0] * nxt[..., 1] - nxt[..., 0] * polygon[..., 1]
    volumes = cross.sum(-1) / 2
    centers = ((polygon + nxt) * cross[..., None]).sum(-2) / (6 * volumes[..., None])
    if not torch.all(volumes > 0):
        raise ValueError("CAD channel grid produced non-positive cell volumes")
    owners, neighbors, labels, ends = [], [], [], []
    known = {}
    for j in range(ny):
        for i in range(nx):
            ids = [(j, i), (j, i + 1), (j + 1, i + 1), (j + 1, i)]
            for a, b in zip(ids, ids[1:] + ids[:1]):
                key = tuple(sorted((a, b)))
                if key in known:
                    neighbors[known[key]] = j * nx + i
                    continue
                known[key] = len(owners)
                owners.append(j * nx + i)
                neighbors.append(-1)
                ends.append(torch.stack((vertices[a], vertices[b])))
                if a[0] == b[0] and a[0] in (0, ny):
                    label = "wall"
                elif a[1] == b[1] == 0:
                    label = "inlet"
                elif a[1] == b[1] == nx:
                    label = "outlet"
                else:
                    label = "interior"
                labels.append(label)
    owner = torch.tensor(owners, dtype=torch.long, device=config.device)
    neighbor = torch.tensor(neighbors, dtype=torch.long, device=config.device)
    face_vertices = torch.stack(ends)
    face_centers = face_vertices.mean(1)
    edge = face_vertices[:, 1] - face_vertices[:, 0]
    face_area_vectors = torch.stack((edge[:, 1], -edge[:, 0]), -1)
    face_lengths = torch.linalg.vector_norm(edge, dim=-1)
    interior = neighbor >= 0
    return type("CadChannelMesh", (), {
        "vertices": vertices,
        "centers": centers,
        "volumes": volumes,
        "field_shape": (ny, nx),
        "owner": owner,
        "neighbor": neighbor,
        "face_vertices": face_vertices,
        "face_centers": face_centers,
        "face_area_vectors": face_area_vectors,
        "face_lengths": face_lengths,
        "face_normals": face_area_vectors / face_lengths[:, None],
        "interior": interior,
        "boundary": ~interior,
        "boundary_labels": tuple(labels),
        "masks": {name: torch.tensor([s == name for s in labels], dtype=torch.bool,
                                     device=config.device)
                  for name in ("inlet", "outlet", "wall", "far-field", "cylinder", "airfoil")},
    })()


def run_fvm_case(*, length=20.0, height=1.0, nx=40, ny=12, ubar=1.0, re_h=5.0,
                 max_iterations=300, tolerance=1e-6) -> dict:
    from tensorfvm.backend_registry import SolverBackend, register_backend
    from tensorfvm.solver import SolverConfig
    from tensorfvm.body_fitted import BodyFittedSolver

    started = time.time()
    case_dir = ART / "fvm_channel"
    case_dir.mkdir(parents=True, exist_ok=True)
    grid = si.channel_grid_2d(length=length, height=height, nx=nx, ny=ny)
    bundle = {
        "schema": si.BUNDLE_SCHEMA,
        "frame": si.COORDINATE_CONTRACT,
        "solver_family": "planar Mesh2D (FVM)",
        "nx": nx,
        "ny": ny,
        "vertices": [[[float(cx), float(cy)] for cx, cy in row] for row in grid],
        "analytic": {"u_max_over_bar": 1.5, "dpdx": -12.0 * (ubar * height / re_h) * ubar
                     / height ** 2},
    }
    mesh_file = case_dir / "channel.json"
    mesh_file.write_text(json.dumps(bundle) + "\n")
    frames = {
        "schema": si.BUNDLE_SCHEMA,
        "case": f"plane channel nx={nx} ny={ny}",
        "frame": dict(si.COORDINATE_CONTRACT),
        "solver_family": "planar Mesh2D (FVM)",
        "tables": {"channel_json": "channel.json"},
        "counts": {"cells": nx * ny},
        "extra": {"length_m": length, "height_m": height,
                  "solver": "TensorFVM BodyFittedSolver (collocated SIMPLE)"},
    }
    (case_dir / "handoff.json").write_text(json.dumps(frames, indent=2, sort_keys=True) + "\n")

    register_backend(SolverBackend(
        name="tensorcad-channel",
        solver_factory=lambda c: BodyFittedSolver(c),
        mesh_factory=_channel_mesh,
        structured=True,
        requires_mesh_file=True,
        required_mesh_masks=("inlet", "outlet", "wall"),
        reference_length="height",
    ))
    config = SolverConfig(
        nx=nx, ny=ny, length=length, height=height, mesh_type="tensorcad-channel",
        mesh_file=str(mesh_file), inlet_velocity=ubar, density=1.0, reynolds=re_h,
        max_iterations=max_iterations, tolerance=tolerance, outer_boundary="channel",
    )
    solver = BodyFittedSolver(config)
    result = solver.solve()
    mesh = solver.mesh
    u_field, v_field = result.cell_center_velocity()
    u = np.asarray(u_field.detach().cpu().numpy(), dtype=float).ravel()
    centers = mesh.centers.reshape(-1, 2).detach().cpu().numpy()
    xc, yc = centers[:, 0], centers[:, 1]
    band = (xc > 0.80 * length) & (xc < 0.90 * length)
    exact = 6.0 * ubar * yc / height * (1.0 - yc / height)
    l2 = float(np.linalg.norm(u[band] - exact[band]) / np.linalg.norm(exact[band]))
    linf = float(np.abs(u[band] - exact[band]).max() / (1.5 * ubar))
    payload = {
        "case": "fvm",
        "solver": "TensorFVM",
        "bundle": str((case_dir / "handoff.json").relative_to(ROOT)),
        "parameters": {"length_m": length, "height_m": height, "nx": nx, "ny": ny,
                       "re_h": re_h, "bar_velocity_m_s": ubar,
                       "nu_m2_s": ubar * height / re_h, "max_iterations": max_iterations},
        "converged": bool(result.converged),
        "iterations": len(result.history) if getattr(result, "history", None) is not None else None,
        "profile_station": [0.80 * length, 0.90 * length],
        "profile_relative_l2": l2,
        "profile_relative_linf_over_umax": linf,
        "mean_velocity_in_band_m_s": float(u[band].mean()),
        "analytic_u_max_m_s": 1.5 * ubar,
        "wall_time_s": time.time() - started,
        "profile": {"y_over_h": (yc[band] / height).tolist(), "u_over_ubar": (u[band] / ubar).tolist()},
    }
    _write_result("fvm_channel", payload)
    return payload


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", default="fem,bem,fvm")
    args = parser.parse_args(argv)
    wanted = {part.strip() for part in args.only.split(",") if part.strip()}
    summary: dict[str, object] = {"schema": si.BUNDLE_SCHEMA, "frame": si.COORDINATE_CONTRACT}
    if "fem" in wanted:
        summary["fem_coarse"] = run_fem_case(nx=7, nz=5, layers=2)
        summary["fem_fine"] = run_fem_case(nx=11, nz=7, layers=2)
    if "bem" in wanted:
        summary["bem"] = run_bem_case()
    if "fvm" in wanted:
        summary["fvm"] = run_fvm_case()
    ART.mkdir(parents=True, exist_ok=True)
    (ART / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: (v if not isinstance(v, dict) else
                          {kk: vv for kk, vv in v.items() if kk in
                           ("patch_test_relative_linf", "equilibrium_relative_error",
                            "levels", "profile_relative_l2", "profile_relative_linf_over_umax",
                            "converged", "wall_time_s")})
                      for k, v in summary.items() if k != "frame"}, indent=2)[:4000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())