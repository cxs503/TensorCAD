"""CAD -> occupancy lattice -> 3-D LBM: the lattice route *consumed*, not just exported.

Why this module exists
----------------------
:mod:`tensorcad.voxel_solid` produces a lattice from a closed CAD surface and
verifies it against closed forms -- but it says of itself, correctly, *"this is the
geometry face of the LBM route, not an LBM result"*.  The only lattice any solver
had ever consumed was the box-aligned 2-D channel whose fluid cells the caller typed
in (``solver_consumers.lbm_channel_from_bundle``).  So the lattice route repeated,
one face further out, the defect this repository keeps finding: **exported is not
accepted**.  Nothing had ever run a flow on a lattice produced from CAD geometry.

Here a closed CAD solid is voxelised, written to ``lattice.json``, re-read from disk
by an independent reader, and then handed to a 3-D LBM solve whose result is gated
against a closed form the geometry never saw.

The case: Hagen--Poiseuille flow in a duct of the CAD solid's own interior
-----------------------------------------------------------------------
A closed cylinder (``nurbs.cylinder_solid_mesh``, watertight by construction) is
voxelised.  Its *interior* cells are the fluid; every other cell is wall.  The
cylinder's axis is the periodic direction, so the lattice is one periodic unit of an
infinitely long duct driven by a constant body force ``g`` -- the standard
fully-developed setup, with no inlet/outlet to model and no entry-length error.  The
closed form is exact for that configuration:

.. math::  u(r) = \\frac{g (R^2 - r^2)}{4 \\nu}, \\qquad
           Q = \\frac{\\pi g R^4}{8 \\nu}, \\qquad
           u_{max} = \\frac{g R^2}{4 \\nu}

The point of the case is not that LBM can do Poiseuille flow; it is that **the
staircase boundary a CAD lattice imposes has a measurable, first-order cost**, and
this module measures it rather than asserting it.  Two references are therefore
reported for every run:

* :math:`R` — the true CAD radius.  The gap contains both the staircase *area*
  deficit and the staircase *shape* error.
* :math:`R_{eff} = \\sqrt{A_{vox}/\\pi}` — the area-equivalent radius of the lattice
  actually written.  The gap to this reference isolates the shape error, i.e. the
  part that a finer lattice must reduce and an area correction cannot.

Contract
--------
Everything is the existing handoff contract, unchanged: ``z-up`` frame, metres, one
cell ordering ``occupancy[i + nx*(j + ny*k)]`` with ``i`` along ``+x`` fastest, and
``verify_bundle`` as a hard gate before the solver reads anything.  The duct axis is
``+z`` (``z-up``), so the flow is vertical -- this is a verification geometry, not a
ship case, and no drag or resistance claim is made.

Cross-repository finding (measured, not inherited by assumption)
----------------------------------------------------------------
:func:`tensorcad.solver_consumers.stream_probe` measured that the 2-D
``tensorlbm.stream`` transports populations against its own ``d2q9.C``.  :func:`stream3d_probe`
here measures the 3-D ``stream3d`` and finds the **same reversal** on D3Q19
(``shifts = (-cz, -cy, -cx)`` where a pull formulation needs ``(+cz, +cy, +cx)``), so
the defect is in both code paths, not in the 2-D one.  TensorLBM's equilibrium,
BGK collision and moments are used unchanged; streaming is local here.  That makes
this case a hybrid and it is labelled as one.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any, Sequence

from .nurbs import cylinder_solid_mesh
from .solver_consumers import gate, solver_fingerprint
from .solver_interface import (
    COORDINATE_CONTRACT,
    HandoffError,
    read_lattice_table,
    verify_bundle,
    write_bundle_lattice,
)
from .voxel_solid import solid_lattice

__all__ = [
    "FLOW_SCHEMA",
    "duct_bundle",
    "duct_lattice",
    "gate_order",
    "lattice_flow_from_bundle",
    "poiseuille_reference",
    "run_lattice_flow_case",
    "stream3d_probe",
]

#: Bumped on any incompatible change to ``results.json`` / the metrics record.
FLOW_SCHEMA = "tensorcad.lattice_flow/1"


# ------------------------------------------------------------------- references
def poiseuille_reference(body_force: float, radius: float, nu: float) -> dict[str, float]:
    """Closed form for a body-force-driven duct, in whatever units the caller uses.

    Lattice units are the convention here (``rho = 1`` so ``mu = nu``); pass
    ``radius`` in cells to get the reference in lattice units.
    """
    if radius <= 0.0:
        raise ValueError(f"radius must be positive, got {radius!r}")
    if nu <= 0.0:
        raise ValueError(f"kinematic viscosity must be positive, got {nu!r}")
    return {
        "radius": radius,
        "u_max": body_force * radius**2 / (4.0 * nu),
        "u_mean": body_force * radius**2 / (8.0 * nu),
        "flow_rate": math.pi * body_force * radius**4 / (8.0 * nu),
    }


def gate_order(errors: Sequence[float], spacings: Sequence[float]) -> dict[str, Any]:
    """Observed order of convergence of ``|error|`` against the cell size.

    Reported as a per-pair order, because a single number over three levels hides a
    level that did not participate.
    """
    pairs = list(zip(errors, errors[1:], spacings, spacings[1:]))
    orders = []
    for e_coarse, e_fine, h_coarse, h_fine in pairs:
        if e_coarse <= 0.0 or e_fine <= 0.0 or h_fine <= 0.0 or h_coarse <= h_fine:
            orders.append(None)
            continue
        orders.append(math.log(e_coarse / e_fine) / math.log(h_coarse / h_fine))
    return {"orders": orders, "scheme": "p = ln(e_coarse/e_fine) / ln(h_coarse/h_fine)"}


# ------------------------------------------------------------------- geometry
def duct_lattice(
    *,
    radius: float,
    length: float,
    radial_cells: int,
    pad_cells: int = 2,
    theta_samples: int = 192,
    v_samples: int = 4,
) -> dict[str, Any]:
    """Voxelise a closed CAD cylinder and return its *interior* as the fluid region.

    ``radial_cells`` is the radius expressed in cells, so the cell size is
    ``radius / radial_cells`` and the lattice is the same geometry at every
    resolution.  ``pad_cells`` cells of empty box are added on the ``x``/``y`` sides
    so the fluid can never touch a lattice face that is not periodic; that makes the
    "outside the lattice bounds" case unreachable for the flow, rather than relying
    on a bounce-back rule at a face that does not exist.
    """
    if radial_cells < 2:
        raise ValueError(f"radial_cells must be at least 2, got {radial_cells!r}")
    h = radius / radial_cells
    pad = pad_cells * h
    lo = (-radius - pad, -radius - pad, 0.0)
    hi = (radius + pad, radius + pad, length)
    mesh = cylinder_solid_mesh(
        radius, length, base_z=0.0, u_samples=theta_samples, v_samples=v_samples
    )
    lattice = solid_lattice(mesh, spacing=h, bounds=[lo, hi], ray_axis="x")
    shape = tuple(int(v) for v in lattice.shape)
    nx, ny, nz = shape
    spacing = tuple(float(v) for v in lattice.spacing)
    solid = lattice.occupancy

    fluid: list[tuple[int, int, int]] = []
    per_plane = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                if solid[i + nx * (j + ny * k)]:
                    fluid.append((i, j, k))
                    per_plane += 1
    if not fluid:
        raise HandoffError("the CAD solid voxelised to no interior cells")
    area = len(fluid) * spacing[0] * spacing[1] / nz
    return {
        "mesh": mesh,
        "shape": shape,
        "origin": lo,
        "spacing": spacing,
        "cell_size_m": h,
        "fluid_cells": fluid,
        "cell_area_m2": spacing[0] * spacing[1],
        "cross_section_area_m2": area,
        "radius_equivalent_m": math.sqrt(area / math.pi),
        "solid_cells": len(fluid),
        "build_seconds": lattice.seconds,
    }


def duct_bundle(
    outdir: str | Path,
    *,
    radius: float = 0.05,
    radial_cells: int = 8,
    axial_cells: int = 4,
    pad_cells: int = 2,
    theta_samples: int = 192,
) -> dict[str, Any]:
    """Write the duct lattice bundle (``lattice.json`` + ``handoff.json``).

    The table is the interface: the solver below never sees ``duct_lattice``'s
    in-memory arrays, it re-reads this file.
    """
    if axial_cells < 2:
        raise ValueError(f"axial_cells must be at least 2, got {axial_cells!r}")
    h = radius / radial_cells
    length = axial_cells * h
    built = duct_lattice(
        radius=radius,
        length=length,
        radial_cells=radial_cells,
        pad_cells=pad_cells,
        theta_samples=theta_samples,
    )
    case = f"duct_poiseuille_r{radial_cells}"
    manifest = write_bundle_lattice(
        outdir,
        shape=built["shape"],
        origin=built["origin"],
        spacing=built["spacing"],
        fluid_cells=built["fluid_cells"],
        walls={
            "kind": "staircase_cad_surface",
            "axial_axis": "z",
            "axial_boundary": "periodic",
            "transverse_boundary": "solid (outside the CAD interior)",
        },
        case=case,
        extra={
            "geometry": {
                "source": "tensorcad.nurbs.cylinder_solid_mesh",
                "radius_m": radius,
                "length_m": length,
                "radial_cells": radial_cells,
                "axial_cells": axial_cells,
                "pad_cells": pad_cells,
                "triangles": len(built["mesh"].faces),
                "vertices": len(built["mesh"].vertices),
            },
            "cross_section": {
                "analytic_area_m2": math.pi * radius**2,
                "voxel_area_m2": built["cross_section_area_m2"],
                "radius_equivalent_m": built["radius_equivalent_m"],
            },
            "frame": COORDINATE_CONTRACT["name"],
        },
    )
    return {"manifest": manifest, "lattice": built, "case": case, "outdir": str(outdir)}


def _absolute_gate(name: str, measured: float, tolerance: float, unit: str) -> dict[str, Any]:
    """A gate whose reference is exact zero, so the error is absolute, not relative."""
    measured_f = float(measured)
    return {
        "name": name,
        "measured": measured_f,
        "reference": 0.0,
        "rel_error": None,
        "abs_error": measured_f,
        "tolerance": float(tolerance),
        "unit": unit,
        "passed": bool(abs(measured_f) <= float(tolerance)),
    }


# ------------------------------------------------------------------ 3-D solver
def _fluid_mask(bundle_dir: str | Path) -> dict[str, Any]:
    """Run the contract verifier, then rebuild the occupancy mask from the file."""
    manifest = verify_bundle(bundle_dir)
    if manifest.get("route") != "lattice":
        raise HandoffError(
            f"the duct flow needs a lattice bundle, got route {manifest.get('route')!r}"
        )
    if manifest.get("frame") != COORDINATE_CONTRACT["name"]:
        raise HandoffError(f"bundle frame {manifest.get('frame')!r} is not the z-up contract")
    data = read_lattice_table(Path(bundle_dir) / "lattice.json")
    nx, ny, nz = data["shape"]
    fluid = [[[False] * nx for _ in range(ny)] for _ in range(nz)]
    for i, j, k in data["fluid_cells"]:
        fluid[k][j][i] = True
    return {"manifest": manifest, "data": data, "fluid": fluid}


def _shifted_wall_mask(mask: Any, c: tuple[int, int, int]) -> Any:
    """Wall mask at ``x + c`` (numpy bool), with periodic ``z`` and wall ``y``/``x``."""
    import numpy as np

    nz, ny, nx = mask.shape
    cz, cy, cx = c
    padded = np.pad(mask, ((0, 0), (1, 1), (1, 1)), constant_values=True)
    z_index = (np.arange(nz) + cz) % nz
    y_index = np.arange(ny) + 1 + cy
    x_index = np.arange(nx) + 1 + cx
    return np.asarray(padded[np.ix_(z_index, y_index, x_index)], dtype=bool)


def lattice_flow_from_bundle(
    bundle_dir: str | Path,
    *,
    tau: float = 0.8,
    body_force: float = 2.0e-5,
    max_steps: int = 40000,
    tol: float = 1.0e-13,
) -> dict[str, Any]:
    """Read a CAD lattice off disk and solve the Poisson--Poiseuille duct flow in 3-D.

    TensorLBM supplies ``equilibrium3d``, ``collide_bgk``, ``macroscopic3d`` and the
    D3Q19 lattice constants.  Streaming is local (pull form ``+c``, half-way
    bounce-back on every non-fluid face, periodic along ``+z``) because the release's
    ``stream3d`` transports populations against its own direction table -- see
    :func:`stream3d_probe`.

    The body force enters as the *plain* Guo source
    ``S_i = w_i [3 c.F + 9 (c.u)(c.F) - 3 u.F]``: it sums to zero (mass conserving)
    and leaves the known constant ``O(F)`` momentum offset that the split form would
    remove -- the same choice, and the same stated offset, as the 2-D channel case.
    """
    import numpy as np
    import torch
    from .solver_consumers import _load  # local: keeps the core import stdlib-only

    # one thread on purpose: the arrays are small (10^3--10^4 cells) and BLAS/LAPACK
    # thread pools on them cost more than they give -- measured 59 CPU-minutes for a
    # 2-minute serial solve before this line existed.
    torch.set_num_threads(1)

    reading = _fluid_mask(bundle_dir)
    data = reading["data"]
    nx, ny, nz = data["shape"]
    fluid = torch.as_tensor(np.asarray(reading["fluid"], dtype=bool))
    wall = ~fluid
    fluid_count = int(fluid.sum().item())
    if nz < 2:
        raise HandoffError("the duct needs at least two cells along the periodic axis")

    lbm = _load("tensorlbm")
    info = solver_fingerprint(lbm)
    c_table = torch.as_tensor(lbm.d3q19.C, dtype=torch.long)
    weights = torch.as_tensor(lbm.d3q19.W, dtype=torch.float64)
    opposite = [int(v) for v in lbm.d3q19.OPPOSITE]
    nu = (tau - 0.5) / 3.0
    if nu <= 0.0:
        raise ValueError(f"tau must exceed 0.5, got {tau!r}")

    rolls = [(int(c_table[i, 2]), int(c_table[i, 1]), int(c_table[i, 0])) for i in range(19)]
    wall_at = [torch.as_tensor(_shifted_wall_mask(wall.numpy(), rolls[i])) for i in range(19)]

    rho0 = torch.ones((nz, ny, nx), dtype=torch.float64)
    zero = torch.zeros((nz, ny, nx), dtype=torch.float64)
    at_rest = lbm.equilibrium3d(rho0, zero, zero, zero)
    fluid_b = fluid.unsqueeze(0)
    f = at_rest.clone()

    force = torch.zeros(3, dtype=torch.float64)
    force[2] = float(body_force)  # along +z, the duct axis
    fw = weights.view(19, 1, 1, 1)

    def forcing(ux: Any, uy: Any, uz: Any) -> Any:
        cx = c_table[:, 0].to(torch.float64).view(19, 1, 1, 1)
        cy = c_table[:, 1].to(torch.float64).view(19, 1, 1, 1)
        cz = c_table[:, 2].to(torch.float64).view(19, 1, 1, 1)
        cf = cx * force[0] + cy * force[1] + cz * force[2]
        cu = cx * ux + cy * uy + cz * uz
        uf = ux * force[0] + uy * force[1] + uz * force[2]
        return fw * (3.0 * cf + 9.0 * cu * cf - 3.0 * uf)

    def collide(f_now: Any) -> Any:
        """Single-relaxation-time collision, written here for D3Q19.

        ``tensorlbm.collide_bgk`` and ``tensorlbm.collide_mrt`` both resolve to the
        D2Q9 implementation in this release (``collide_mrt`` literally calls
        ``collide_bgk``, which calls ``d2q9.macroscopic``), so there is no usable 3-D
        collision to call.  Equilibrium and moments are still TensorLBM's; the
        relaxation ``f - (f - feq)/tau`` is one line and is stated as local.
        """
        rho_l, ux_l, uy_l, uz_l = lbm.macroscopic3d(f_now)
        return f_now - (f_now - lbm.equilibrium3d(rho_l, ux_l, uy_l, uz_l)) / tau

    def stream(f_post: Any) -> Any:
        out = torch.empty_like(f_post)
        for i in range(19):
            cz, cy, cx = rolls[i]
            out[i] = torch.roll(f_post[i], shifts=(cz, cy, cx), dims=(0, 1, 2))
        for i in range(1, 19):
            j = opposite[i]
            if j == i:
                continue
            out[j] = torch.where(wall_at[i], f_post[i], out[j])
        return out

    started = time.perf_counter()
    steps = 0
    change = float("inf")
    while steps < max_steps:
        _, ux, uy, uz = lbm.macroscopic3d(f)
        f_post = collide(f) + forcing(ux, uy, uz)
        f_next = torch.where(fluid_b, stream(f_post), at_rest)
        change = float(torch.where(fluid, f_next - f, torch.zeros_like(f)).abs().max())
        f = f_next
        steps += 1
        if change < tol:
            break
    solve_seconds = time.perf_counter() - started

    rho, ux, uy, uz = lbm.macroscopic3d(f)
    axial = torch.where(fluid, uz, torch.zeros_like(uz))
    transverse = torch.where(fluid, torch.sqrt(ux**2 + uy**2), torch.zeros_like(uz))
    fluid_area_cells = fluid_count / nz
    r_eff = math.sqrt(fluid_area_cells / math.pi)
    reference = poiseuille_reference(body_force, r_eff, nu)
    # the flow rate is a flux *through one cross-section*: summing every plane of a
    # developed duct over-counts it by the number of periodic planes.
    measured_rate = float(axial.sum().item()) / nz
    # radial profile for the shape evidence: bin the fluid cells by distance from the axis
    idx = torch.nonzero(fluid)  # (n, 3) as (k, j, i)
    mid = torch.tensor([(nz - 1) / 2.0, (ny - 1) / 2.0, (nx - 1) / 2.0], dtype=torch.float64)
    radial = torch.sqrt(((idx[:, 1:].to(torch.float64) - mid[1:]) ** 2).sum(dim=1))
    values = axial[fluid]
    bins = torch.unique(torch.round(radial * 2.0) / 2.0)
    profile = []
    for radius_bin in bins.tolist():
        selected = torch.abs(radial - radius_bin) <= 0.25
        if bool(selected.any()):
            profile.append([radius_bin, float(values[selected].mean().item()), int(selected.sum().item())])
    return {
        "case": str(reading["manifest"].get("case", "")),
        "route": "lattice_json",
        "solver": "tensorlbm",
        "schema": FLOW_SCHEMA,
        "status": "ran",
        "solver_info": info,
        "note": (
            "CAD interior lattice from lattice.json; TensorLBM equilibrium3d/BGK/macroscopic3d "
            "+ local D3Q19 pull streaming with half-way bounce-back (stream3d sign reversed, measured)"
        ),
        "metrics": {
            "shape": [nx, ny, nz],
            "flow_axis": "+z",
            "tau": tau,
            "nu_lattice": nu,
            "body_force": body_force,
            "fluid_cells": fluid_count,
            "fluid_area_cells": fluid_area_cells,
            "radius_equivalent_cells": r_eff,
            "steps": steps,
            "converged": change < tol,
            "steady_change": change,
            "solve_seconds": solve_seconds,
            "u_max_lattice": float(axial.max().item()),
            "u_mean_lattice": measured_rate / fluid_area_cells,
            "flow_rate_lattice": measured_rate,
            "reference_u_max_lattice": reference["u_max"],
            "reference_u_mean_lattice": reference["u_mean"],
            "reference_flow_rate_lattice": reference["flow_rate"],
            "mass_drift": float(abs(rho[fluid].sum().item() - fluid_count)),
            "transverse_speed_max": float(transverse.max().item()),
            "axial_profile": profile,
            "bundle_digest": reading["manifest"].get("digest_sha256"),
        },
    }


def stream3d_probe() -> dict[str, Any]:
    """Measure which way ``tensorlbm.stream3d`` actually transports a population.

    One population on the ``+x`` direction at one node; after a single streaming step
    it must appear one node further in ``+x`` for a ``(19, nz, ny, nx)`` table.  The
    reported ``dx_lattice`` is measured.
    """
    try:
        import torch
        from .solver_consumers import _load

        lbm = _load("tensorlbm")
        shape = (19, 5, 5, 5)
        f = torch.zeros(shape, dtype=torch.float64)
        f[1, 2, 2, 2] = 1.0
        moved = lbm.stream3d(f)
        index = torch.nonzero(moved[1] > 0.5)
        if index.numel() == 0:
            return {"status": "indeterminate", "dx_lattice": None}
        dx = int(index[0][2].item()) - 2
        return {
            "status": "measured",
            "dx_lattice": dx,
            "expected_dx_lattice": 1,
            "agrees_with_pull_convention": dx == 1,
            "note": "stream3d(f) must move a +x population to x+1 for a (19, nz, ny, nx) table",
        }
    except Exception as exc:  # noqa: BLE001 - a probe must not fail the run
        return {"status": "unavailable", "error": f"{type(exc).__name__}: {exc}"}


# --------------------------------------------------------------------- driver
def run_lattice_flow_case(
    outdir: str | Path,
    *,
    resolutions: Sequence[int] = (6, 10, 14),
    radius: float = 0.05,
    axial_cells: int = 4,
    tau: float = 0.8,
    body_force: float = 2.0e-5,
    max_steps: int = 40000,
    tol: float = 1.0e-13,
    theta_samples: int = 128,
) -> dict[str, Any]:
    """Voxelise the CAD duct at several resolutions, solve each, and gate the numbers."""
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    nu = (tau - 0.5) / 3.0
    cases: list[dict[str, Any]] = []
    for radial_cells in resolutions:
        bundle_dir = out / f"duct_r{radial_cells}"
        written = duct_bundle(
            bundle_dir,
            radius=radius,
            radial_cells=radial_cells,
            axial_cells=axial_cells,
            theta_samples=theta_samples,
        )
        result = lattice_flow_from_bundle(
            bundle_dir, tau=tau, body_force=body_force, max_steps=max_steps, tol=tol
        )
        metrics = result["metrics"]
        geometry = written["lattice"]
        r_equiv_m = geometry["radius_equivalent_m"]
        r_equiv_cells = metrics["radius_equivalent_cells"]
        true_reference = poiseuille_reference(body_force, float(radial_cells), nu)
        area_error = geometry["cross_section_area_m2"] / (math.pi * radius**2) - 1.0
        gates = [
            gate("voxel_cross_section_area_m2", geometry["cross_section_area_m2"],
                 math.pi * radius**2, 0.10, "m2"),
            gate("flow_rate_lattice", metrics["flow_rate_lattice"],
                 true_reference["flow_rate"], 0.25, "cell^2/step"),
            gate("flow_rate_vs_area_equivalent", metrics["flow_rate_lattice"],
                 metrics["reference_flow_rate_lattice"], 0.10, "cell^2/step"),
            gate("u_max_lattice", metrics["u_max_lattice"], true_reference["u_max"], 0.25, "cell/step"),
            gate("u_max_vs_area_equivalent", metrics["u_max_lattice"],
                 metrics["reference_u_max_lattice"], 0.10, "cell/step"),
            _absolute_gate(
                "mass_drift", metrics["mass_drift"], 1e-8 * max(1.0, metrics["fluid_cells"]), "mass"
            ),
            _absolute_gate(
                "transverse_speed_max",
                metrics["transverse_speed_max"],
                5.0 * body_force * float(radial_cells) ** 2 / (4.0 * nu),
                "cell/step",
            ),
        ]
        cases.append(
            {
                "radial_cells": radial_cells,
                "bundle": str(bundle_dir),
                "frame": written["manifest"].get("frame"),
                "triangles": len(geometry["mesh"].faces),
                "cell_size_m": geometry["cell_size_m"],
                "cross_section_area_m2": geometry["cross_section_area_m2"],
                "radius_equivalent_m": r_equiv_m,
                "radius_equivalent_cells": r_equiv_cells,
                "area_rel_error": area_error,
                "result": result,
                "reference": true_reference,
                "gates": gates,
                "gates_failed": [g["name"] for g in gates if not g["passed"]],
            }
        )

    errors = [abs(c["result"]["metrics"]["flow_rate_lattice"] / c["reference"]["flow_rate"] - 1.0) for c in cases]
    equiv_errors = [
        abs(c["result"]["metrics"]["flow_rate_lattice"] / c["result"]["metrics"]["reference_flow_rate_lattice"] - 1.0)
        for c in cases
    ]
    spacings = [c["cell_size_m"] / radius for c in cases]
    cone = gate_order(errors, spacings)
    cone_equiv = gate_order(equiv_errors, spacings)
    payload = {
        "schema": FLOW_SCHEMA,
        "case_family": "duct_poiseuille_3d",
        "frame": COORDINATE_CONTRACT["name"],
        "radius_m": radius,
        "axial_cells": axial_cells,
        "tau": tau,
        "nu_lattice": nu,
        "body_force": body_force,
        "stream3d_probe": stream3d_probe(),
        "cases": cases,
        "convergence": {
            "resolutions": list(resolutions),
            "cell_size_over_radius": spacings,
            "flow_rate_rel_error_vs_cad_radius": errors,
            "flow_rate_rel_error_vs_area_equivalent": equiv_errors,
            "order_vs_cad_radius": cone,
            "order_vs_area_equivalent": cone_equiv,
        },
        "gates_failed": [
            f"r{c['radial_cells']}:{name}" for c in cases for name in c["gates_failed"]
        ],
    }
    (out / "results.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


def _main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default="artifacts/lattice_flow")
    parser.add_argument("--resolutions", default="6,10,14")
    parser.add_argument("--tau", type=float, default=0.8)
    parser.add_argument("--body-force", type=float, default=2.0e-5)
    parser.add_argument("--max-steps", type=int, default=40000)
    args = parser.parse_args(argv)
    resolutions = tuple(int(v) for v in args.resolutions.split(","))
    payload = run_lattice_flow_case(
        args.out,
        resolutions=resolutions,
        tau=args.tau,
        body_force=args.body_force,
        max_steps=args.max_steps,
    )
    for case in payload["cases"]:
        metrics = case["result"]["metrics"]
        print(
            f"r={case['radial_cells']:>2}  h={case['cell_size_m']*1e3:7.3f} mm  "
            f"Q={metrics['flow_rate_lattice']:.6e}  ref={case['reference']['flow_rate']:.6e}  "
            f"area={case['area_rel_error']*100:+.3f}%  steps={metrics['steps']}"
        )
    print("gates_failed:", payload["gates_failed"] or "none")
    print("stream3d_probe:", payload["stream3d_probe"])
    return 0 if not payload["gates_failed"] else 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())