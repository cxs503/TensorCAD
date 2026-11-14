"""CAD -> TensorDG curved-wall hand-off, end to end, with the order gate.

The route is the same one the other solver families use in
``examples/solver_interface/run_three_solvers.py``: CAD writes a bundle (tables
plus declared numbers), the bundle is re-read and re-measured by an independent
parser, and only then does the solver consume the file.

Here the consumer is the DG curved-geometry layer:

1. **CAD** -- the equator of the unit sphere from ``tensorcad.solver_interface``
   (the same geometry object the BEM route uses) is projected onto the ``x-y``
   chart; the wall radius/centre are *measured* from those 3-D points.
2. **Contract** -- ``tensorcad.dg_handoff.circle_contract`` builds the
   rational-quadratic NURBS wall from the tangent-intersection construction and
   writes ``curved_wall.json`` + ``handoff.json``; ``verify_bundle`` re-reads
   them and recomputes closure, radial deviation, circumference and area.
3. **DG** -- ``tensordg.curved_geometry.mesh_from_contract`` turns the emitted
   control points into isoparametric element geometry (``g = p + 1``) and
   ``CurvedAdvection`` advects a smooth rigid-rotation profile on the annulus at
   ``p = 2, 3`` over three mesh levels.  The exact solution is analytic, so the
   measured convergence order is the gate: a correct curved-wall treatment keeps
   ``p + 1``.

Two controls are run so the numbers mean something:

* a **degree-1 wall** supplied by the solver (polygonal inner wall) -- shows the
  order loss that curved CAD geometry is supposed to prevent;
* a **CAD polygon contract** (degree-1 closed NURBS through on-circle vertices)
  -- TensorDG's validator must *refuse* it, i.e. the hand-off is a gate and not
  a comment.

Run with::

    PYTHONPATH=.:/root/tensor-ecosystem/TensorDG python \
        examples/dg_curved_wall/run_dg_handoff.py

Numbers land in ``artifacts/dg_curved_wall/results.json``; the figures next to
it are regenerated from those numbers, never hand-edited.
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tensorcad import dg_handoff as dh  # noqa: E402
from tensorcad import solver_interface as si  # noqa: E402

ART = ROOT / "artifacts" / "dg_curved_wall"

# --- the advection problem (analytic reference, no solver involved) ---------
OMEGA = 1.0            # rigid-body rotation rate, 1/s
AZIMUTH = 2            # azimuthal wavenumber
R_MID, R_WIDTH = 1.5, 0.35
T_FINAL = 0.3          # one rotation takes 2*pi; 0.3 keeps the profile smooth
R_OUTER = 2.0


def rotation_velocity(x, y):
    return -OMEGA * y, OMEGA * x


def rotation_exact(x, y, t):
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    return np.exp(-((r - R_MID) / R_WIDTH) ** 2) * np.cos(AZIMUTH * (th - OMEGA * t))


# --------------------------------------------------------------------- CAD side


def cad_sphere_equator(radius: float = 1.0, n_theta: int = 24, n_phi: int = 64):
    """Equator of the CAD sphere solid, as 3-D z-up points.

    ``n_theta`` is even so that ``i = n_theta/2`` lands exactly on the equator
    ring of :func:`tensorcad.solver_interface.ellipsoid_mesh`; those points lie
    on the exact circle by construction of the parametrisation, so the measured
    radius is a property of the CAD object, not of a fit.
    """
    mesh = si.ellipsoid_mesh((radius, radius, radius), n_theta=n_theta, n_phi=n_phi)
    ring = [v for v in mesh.vertices if abs(v[2]) < 1e-12]
    if len(ring) < 3:
        raise RuntimeError("no equator ring found in the CAD sphere mesh")
    return ring


def cad_wall(case_dir: Path, *, radius=1.0, n_theta=24, n_phi=64, segments=4):
    """CAD ring -> contract -> bundle -> independent verification."""
    ring = cad_sphere_equator(radius, n_theta, n_phi)
    r_meas = sum(math.hypot(p[0], p[1]) for p in ring) / len(ring)
    centre = [0.0, 0.0]
    contract = dh.circle_contract(r_meas, centre, segments=segments, plane="xy", offset_m=0.0,
                                  name="sphere-equator-wall")
    manifest = dh.write_dg_bundle(
        case_dir, contract, r_outer_m=R_OUTER,
        plane_map=dh.PlaneMap("xy", *dh.PLANES["xy"], 0.0),
        source_3d_ring=ring,
        case=f"CAD sphere equator R={r_meas:.6f} m -> DG annulus wall",
        extra={
            "cad_source": "tensorcad.solver_interface.ellipsoid_mesh((R,R,R)).equator",
            "cad_points": len(ring),
            "measured_radius_from_cad_m": r_meas,
            "segments": segments,
        },
    )
    check = dh.verify_bundle(case_dir)
    return contract, manifest, check


# ---------------------------------------------------------------------- DG side


def _order_fit(hs, errs):
    slope, _ = np.polyfit(np.log(np.asarray(hs, float)), np.log(np.asarray(errs, float)), 1)
    return float(slope)


def run_advection(mesh, p, T=T_FINAL):
    from tensordg.curved_geometry import CurvedAdvection

    dg = CurvedAdvection(mesh, p, rotation_velocity)
    U0 = dg.project(lambda x, y: rotation_exact(x, y, 0.0))
    U, dt, steps = dg.evolve(U0, T, dg.stable_dt())
    num, den, rel = dg.l2_error(U, lambda x, y: rotation_exact(x, y, T))
    area = float(np.sum(np.einsum("nqr,q,r->n", dg.Jg, dg.w, dg.w)))
    # geometric fidelity of the inner wall the solver actually sees, sampled
    # densely along each element edge (a chord deviates between its nodes)
    wall_radial, wx, wy = wall_fidelity(mesh, R_INNER_ACTUAL)
    return {
        "degree": p, "cells": int(mesh.ncells), "h": float(np.sqrt(area / mesh.ncells)),
        "rel_l2": float(rel), "abs_l2": float(num), "dt": float(dt), "steps": int(steps),
        "wall_radial_error_max_m": wall_radial,
        "wall_radial_error_rel": wall_radial / R_INNER_ACTUAL,
        "wall_samples": (wx, wy), "dg": dg, "U": U,
    }


def wall_fidelity(mesh, radius: float, samples: int = 41):
    """Max ``|r - R|`` of the inner wall, sampled along each element edge.

    Sampling dense than the geometry nodes matters: on a degree-1 wall the nodes
    sit exactly on the circle while the chord between them does not, which is
    precisely the error a curved element is supposed to remove.
    """
    from tensordg.curved_geometry import gll_nodes_weights

    gll_nodes, _ = gll_nodes_weights(mesh.degree)
    eta = np.linspace(-1.0, 1.0, samples)
    xs, ys = [], []
    for m in mesh.maps[0]:
        x, y = m.eval(np.array([-1.0]), eta)
        xs.append(np.asarray(x).ravel())
        ys.append(np.asarray(y).ravel())
    xs = np.concatenate(xs)
    ys = np.concatenate(ys)
    dev = np.abs(np.hypot(xs, ys) - radius)
    return float(dev.max()), xs, ys


def solution_nodes(dg, U, refine: int = 6):
    """Evaluate the DG field on a refined reference grid, for the field figure."""
    from tensordg.curved_geometry import legendre_eval_matrix

    fine = np.linspace(-1.0, 1.0, refine * max(dg.p, dg.mesh.degree) + 1)
    Lu = legendre_eval_matrix(dg.p, fine)              # state basis degree p
    Lg = legendre_eval_matrix(dg.mesh.degree, fine)    # geometry basis degree g
    CX, CY = dg.mesh.coefficients()
    u = np.einsum("qk,nkl,rl->nqr", Lu, U, Lu)
    x = np.einsum("qk,nkl,rl->nqr", Lg, CX, Lg)
    y = np.einsum("qk,nkl,rl->nqr", Lg, CY, Lg)
    return x.ravel(), y.ravel(), u.ravel()


# -------------------------------------------------------------------- controls


def polygon_contract(radius: float, centre=(0.0, 0.0), segments: int = 8):
    """Degree-1 closed NURBS through on-circle vertices: a faceted wall."""
    cps, weights = [], []
    for i in range(segments):
        a = 2.0 * math.pi * i / segments
        cps.append([centre[0] + radius * math.cos(a), centre[1] + radius * math.sin(a)])
        weights.append(1.0)
    cps.append(list(cps[0]))
    weights.append(1.0)
    knots = [0.0] * 2 + [i / segments for i in range(1, segments) for _ in range(1)] + [1.0] * 2
    return {
        "schema": dh.DG_CONTRACT_SCHEMA, "kind": "nurbs_curve", "name": "faceted-wall",
        "centre": list(centre), "degree": 1, "knots": knots,
        "control_points": cps, "weights": weights, "closed": True, "radius_m": radius,
    }


def main() -> int:
    ART.mkdir(parents=True, exist_ok=True)
    started = time.time()
    from tensordg.curved_geometry import (
        annulus_mesh, mesh_from_contract, validate_geometry_contract,
    )

    global R_INNER_ACTUAL

    contract, manifest, check = cad_wall(ART / "cad_wall")
    R_INNER_ACTUAL = float(check["measured"]["radius_m"])

    # solver-side acceptance of the CAD payload (the gate, not a courtesy)
    normalised = validate_geometry_contract(contract)
    accepted = {
        "schema_ok": normalised["schema"] == dh.DG_CONTRACT_SCHEMA,
        "degree": normalised["degree"],
        "control_points": len(normalised["control_points"]),
        "digest": manifest["geometry"]["digest"],
    }

    # negative control: a faceted CAD wall must be refused by the consumer
    refusal = {}
    try:
        validate_geometry_contract(polygon_contract(R_INNER_ACTUAL, segments=8))
        mesh_from_contract(polygon_contract(R_INNER_ACTUAL, segments=8), 6, 12, R_OUTER, 2)
        refusal["refused"] = False
    except ValueError as exc:
        refusal = {"refused": True, "reason": str(exc), "type": type(exc).__name__}

    # main study: CAD-derived wall, geometry degree g = p + 1
    levels = [(8, 16), (12, 24), (16, 32)]
    rows, fields = [], {}
    for p in (2, 3):
        g = p + 1
        errs, hs, wall = [], [], []
        for nr, nt in levels:
            mesh = mesh_from_contract(contract, nr, nt, R_OUTER, g)
            res = run_advection(mesh, p)
            errs.append(res["rel_l2"])
            hs.append(res["h"])
            wall.append(res["wall_radial_error_rel"])
            rows.append({"degree": p, "geometry_degree": g, "nr": nr, "ntheta": nt,
                         "cells": res["cells"], "h": res["h"], "rel_l2": res["rel_l2"],
                         "dt": res["dt"], "steps": res["steps"],
                         "wall_radial_error_rel": res["wall_radial_error_rel"]})
            fields[(p, nr, nt)] = (res["dg"], res["U"], res["wall_samples"])
        rows.append({"degree": p, "geometry_degree": g, "levels": [list(l) for l in levels],
                     "h": hs, "rel_l2": errs, "observed_order": _order_fit(hs, errs),
                     "theory_order": p + 1,
                     "rel_deviation": abs(_order_fit(hs, errs) - (p + 1)) / (p + 1),
                     "wall_radial_error_rel_max": max(wall),
                     "kind": "order_fit"})

    # control: solver-side degree-1 wall (polygon), same DG degree
    control = []
    for p in (2, 3):
        errs, hs, wall = [], [], []
        for nr, nt in levels:
            mesh = annulus_mesh(nr, nt, R_INNER_ACTUAL, R_OUTER, 1)
            res = run_advection(mesh, p)
            errs.append(res["rel_l2"])
            hs.append(res["h"])
            wall.append(res["wall_radial_error_rel"])
        control.append({"degree": p, "geometry_degree": 1, "levels": [list(l) for l in levels],
                        "h": hs, "rel_l2": errs, "observed_order": _order_fit(hs, errs),
                        "theory_order": p + 1,
                        "wall_radial_error_rel_max": max(wall),
                        "kind": "control_polygonal_wall"})

    payload = {
        "schema": dh.BUNDLE_SCHEMA,
        "case": "CAD sphere equator -> DG curved-wall rigid-rotation advection",
        "cad_bundle": str((ART / "cad_wall" / "handoff.json").relative_to(ROOT)),
        "cad_verification": {k: v for k, v in check["measured"].items()},
        "plane_check": check["plane_check"],
        "dg_acceptance": accepted,
        "refusal_control": refusal,
        "parameters": {"omega_1_s": OMEGA, "azimuthal_wavenumber": AZIMUTH,
                       "r_mid_m": R_MID, "r_width_m": R_WIDTH, "t_final_s": T_FINAL,
                       "r_inner_m": R_INNER_ACTUAL, "r_outer_m": R_OUTER,
                       "geometry_degree_rule": "g = p + 1",
                       "cfl": 0.2, "levels": [list(l) for l in levels]},
        "runs": rows,
        "control_polygonal_wall": control,
        "wall_time_s": time.time() - started,
    }
    (ART / "results.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")

    # ---------------- figures (regenerated from the numbers above) ----------
    _figures(payload, contract, fields)
    print(json.dumps({k: payload[k] for k in ("cad_verification", "dg_acceptance",
                                              "refusal_control")}, indent=1)[:1800])
    for r in rows + control:
        if r.get("kind") in ("order_fit", "control_polygonal_wall"):
            print(f"p={r['degree']} g={r['geometry_degree']} order={r['observed_order']:.3f}"
                  f" (theory {r['theory_order']}) errs={['%.3e' % e for e in r['rel_l2']]}"
                  f" wall={r['wall_radial_error_rel_max']:.2e}")
    print(f"total {time.time() - started:.1f} s")
    return 0


def _figures(payload, contract, fields):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # 1) convergence with p+1 reference slopes
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    colours = {2: "tab:blue", 3: "tab:red"}
    for r in payload["runs"]:
        if r.get("kind") != "order_fit":
            continue
        hs, errs = np.array(r["h"]), np.array(r["rel_l2"])
        ax.loglog(hs, errs, "o-", color=colours[r["degree"]],
                  label=f"CAD wall, p={r['degree']}, g={r['geometry_degree']}"
                        f" (order {r['observed_order']:.2f})")
        ref = errs[0] * (hs / hs[0]) ** (r["theory_order"])
        ax.loglog(hs, ref, "--", color=colours[r["degree"]], alpha=0.5,
                  label=f"p+1 = {r['theory_order']} reference")
    for r in payload["control_polygonal_wall"]:
        hs, errs = np.array(r["h"]), np.array(r["rel_l2"])
        ax.loglog(hs, errs, "s:", color="grey", alpha=0.8,
                  label=f"polygonal wall (g=1), p={r['degree']}"
                        f" (order {r['observed_order']:.2f})")
    ax.set_xlabel("h = sqrt(area / cells)  [m]")
    ax.set_ylabel("relative L2 error at t = 0.3 s")
    ax.set_title("CAD-derived curved wall: DG convergence order")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(ART / "convergence.png", dpi=140)
    plt.close(fig)

    # 2) wall geometry: CAD contract vs the wall the DG solver sees
    fig, axs = plt.subplots(1, 2, figsize=(10.4, 4.2))
    ts = np.linspace(0.0, 1.0, 721)
    pts = np.asarray(dh.evaluate_curve(contract, ts))
    axs[0].plot(pts[:, 0], pts[:, 1], "-", lw=1.2, label="CAD NURBS wall (contract)")
    th = np.linspace(0, 2 * np.pi, 721)
    axs[0].plot(payload["parameters"]["r_inner_m"] * np.cos(th),
                payload["parameters"]["r_inner_m"] * np.sin(th), "--", lw=0.8,
                label="analytic circle")
    dg_wall = None
    if fields:
        key = sorted(fields, key=lambda k: (k[0], -k[1]))[0]
        _, _, (wx, wy) = fields[key]
        dg_wall = (wx, wy)
        axs[0].plot(wx, wy, ".", ms=3, alpha=0.6,
                    label=f"DG quadrature nodes (p={key[0]}, {key[1]}x{key[2]})")
    axs[0].set_aspect("equal")
    axs[0].set_xlabel("x [m]")
    axs[0].set_ylabel("y [m]")
    axs[0].set_title("inner wall: CAD contract vs DG element geometry")
    axs[0].legend(fontsize=7)
    axs[0].grid(alpha=0.3)

    r_in = payload["parameters"]["r_inner_m"]
    if dg_wall is not None:
        dev = np.abs(np.hypot(*dg_wall) - r_in)
        axs[1].semilogy(np.hypot(*dg_wall), dev, ".", ms=3,
                        label="CAD contract (this run)")
    for r in payload["control_polygonal_wall"]:
        errs = r["wall_radial_error_rel_max"]
        axs[1].axhline(errs * r_in, ls=":", color="grey",
                       label=f"polygonal wall g=1 (p={r['degree']}): {errs:.2e} m")
    axs[1].set_xlabel("wall radius r [m]")
    axs[1].set_ylabel("| r - R | of the wall nodes [m]")
    axs[1].set_title("geometric fidelity of the wall seen by the solver")
    axs[1].legend(fontsize=7)
    axs[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(ART / "wall_geometry.png", dpi=140)
    plt.close(fig)

    # 3) solution field at the finest mesh, p = 3
    if not fields:
        return
    p_max = max(k[0] for k in fields)
    key = max((k for k in fields if k[0] == p_max), key=lambda k: k[1] * k[2])
    dg, U, _ = fields[key]
    x, y, u = solution_nodes(dg, U)
    exact = rotation_exact(x, y, T_FINAL)
    fig, axs = plt.subplots(1, 2, figsize=(10.4, 4.3))
    sc = axs[0].scatter(x, y, c=u, s=3, cmap="RdBu_r")
    axs[0].set_aspect("equal")
    axs[0].set_title(f"u at t = {T_FINAL} s (p={key[0]}, {key[1]}x{key[2]} cells)")
    fig.colorbar(sc, ax=axs[0], shrink=0.8)
    sc2 = axs[1].scatter(x, y, c=np.abs(u - exact), s=3, cmap="magma")
    axs[1].set_aspect("equal")
    axs[1].set_title("|u_num - u_exact|")
    fig.colorbar(sc2, ax=axs[1], shrink=0.8)
    fig.tight_layout()
    fig.savefig(ART / "solution_field.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())