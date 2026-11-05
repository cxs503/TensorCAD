"""Cure-PR-#39 verification: the union of main and agent/m2925, measured.

Three questions, each answered with a number rather than an assertion of intent:

1. m2925's route is alive on the merged tree -- ``compute_hydrostatics`` reports the
   waterplane second moments and ``keel_z``, and ``analytic_hydrostatics`` reproduces
   them from closed forms (rectangle barge exact, hull mesh converging).
2. main's route is alive on the merged tree -- the asymmetric fore/aft exponents
   (LCB/LCF away from midship, Richardson-free quadrature) and the transom and Wigley
   families still build and still agree with their own references.
3. The two routes compose honestly -- the inertias are withheld (nan) exactly where
   their closed forms do not hold, i.e. on a fore/aft-asymmetric waterplane, and the
   ``basis`` string says so.

Writes ``artifacts/m3768/merge_verification.json`` and the report figure.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from tensorcad import transom, wigley
from tensorcad.geometry import MeshBuilder
from tensorcad.hull import HullParameters, analytic_hydrostatics, generate_hull
from tensorcad.hydrostatics import compute_hydrostatics

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "m3768"

L, B, T = 10.0, 2.0, 0.8


def barge(length: float, beam: float, depth: float):
    """Axis-aligned box: ``I_T = L B^3/12``, ``I_L = B L^3/12``, ``keel_z = -T``.

    Built face by face so the reference is independent of the hull loft: only the
    waterplane rectangle (at ``z = 0``) and the keel plane (``z = -depth``) matter.
    """
    builder = MeshBuilder()
    hx, hy = length / 2.0, beam / 2.0
    x0, x1, y0, y1, z0, z1 = -hx, hx, -hy, hy, -depth, 0.0
    v = {
        (ix, iy, iz): builder.add_vertex((x, y, z))
        for ix, x in ((0, x0), (1, x1))
        for iy, y in ((0, y0), (1, y1))
        for iz, z in ((0, z0), (1, z1))
    }
    quads = (
        ((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0)),  # keel
        ((0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)),  # waterplane
        ((0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)),
        ((0, 1, 0), (0, 1, 1), (1, 1, 1), (1, 1, 0)),
        ((0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0)),
        ((1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)),
    )
    for quad in quads:
        builder.add_quad(*(v[key] for key in quad))
    return builder.build()


def main() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    data: dict[str, object] = {"geometry": {"L": L, "B": B, "T": T}}

    # (1) m2925 route: barge against the textbook second moments.
    box = compute_hydrostatics(barge(L, B, T))
    ref_it, ref_il = L * B**3 / 12.0, B * L**3 / 12.0
    data["barge"] = {
        "I_T": box.transverse_waterplane_inertia,
        "I_T_ref": ref_it,
        "I_T_rel_err": (box.transverse_waterplane_inertia - ref_it) / ref_it,
        "I_L": box.longitudinal_waterplane_inertia,
        "I_L_ref": ref_il,
        "I_L_rel_err": (box.longitudinal_waterplane_inertia - ref_il) / ref_il,
        "keel_z": box.keel_z,
        "keel_z_ref": -T,
    }

    # (2) mesh vs closed forms on the parametric family, several grids.
    grids = [(21, 9), (41, 13), (81, 25), (161, 41)]
    rows = []
    for nx, nz in grids:
        model = generate_hull(HullParameters(L, B, T, nx, nz))
        mesh = compute_hydrostatics(model.mesh)
        ref = analytic_hydrostatics(L, B, T)
        rows.append(
            {
                "nx": nx,
                "nz": nz,
                "I_T_mesh": mesh.transverse_waterplane_inertia,
                "I_T_ref": ref.transverse_waterplane_inertia,
                "I_T_rel_err": (mesh.transverse_waterplane_inertia - ref.transverse_waterplane_inertia)
                / ref.transverse_waterplane_inertia,
                "I_L_mesh": mesh.longitudinal_waterplane_inertia,
                "I_L_ref": ref.longitudinal_waterplane_inertia,
                "I_L_rel_err": (mesh.longitudinal_waterplane_inertia - ref.longitudinal_waterplane_inertia)
                / ref.longitudinal_waterplane_inertia,
                "keel_z_mesh": mesh.keel_z,
                "basis_analytic": ref.basis,
            }
        )
    data["hull_vs_closed_form"] = rows

    # (3) main's route: asymmetric exponents + the two newer hull families.
    asym = analytic_hydrostatics(L, B, T, fore_exponent=0.5, aft_exponent=0.9)
    sym = analytic_hydrostatics(L, B, T)
    data["asymmetric"] = {
        "LCB": asym.longitudinal_centre_buoyancy,
        "LCF": asym.longitudinal_centre_flotation,
        "LCB_symmetric": sym.longitudinal_centre_buoyancy,
        "mid": 0.5 * L,
        "I_T_is_nan": math.isnan(asym.transverse_waterplane_inertia),
        "I_L_is_nan": math.isnan(asym.longitudinal_waterplane_inertia),
        "keel_z": asym.keel_z,
        "basis": asym.basis,
        "Cb": asym.block_coefficient,
    }

    wparams = wigley.WigleyParameters()
    wmodel = wigley.generate_wigley(wparams)
    wmesh = wigley.compute_hydrostatics(wmodel.mesh)
    wref = wigley.wigley_analytic_hydrostatics(wparams)
    data["wigley"] = {
        "Cb_mesh": wmesh.block_coefficient,
        "Cb_ref": wref.block_coefficient,
        "Cb_rel_err": (wmesh.block_coefficient - wref.block_coefficient) / wref.block_coefficient,
        "volume_mesh": wmesh.displacement_volume,
        "volume_ref": wref.displacement_volume,
        "volume_rel_err": (wmesh.displacement_volume - wref.displacement_volume)
        / wref.displacement_volume,
        "analytic_kind": type(wref).__name__,
        "record_has_inertia": hasattr(wref, "transverse_waterplane_inertia"),
    }

    tparams = transom.TransomHullParameters()
    tmodel = transom.generate_transom_hull(tparams)
    tmesh = transom.compute_hydrostatics(tmodel.mesh)
    tref = transom.analytic_transom_hydrostatics(tparams)
    data["transom"] = {
        "tau": getattr(tparams, "transom_fraction", None),
        "Cb_mesh": tmesh.block_coefficient,
        "Cb_ref": tref["block_coefficient"],
        "Cb_rel_err": (tmesh.block_coefficient - tref["block_coefficient"])
        / tref["block_coefficient"],
        "record_kind": "dict (transom reference ladder)",
        "volume_rel_err": (tmesh.displacement_volume - tref["displacement_volume"])
        / tref["displacement_volume"],
    }

    # Figure: signed mesh error of the two waterplane second moments vs grid, plus
    # main's asymmetric LCB/LCF route as an inset bar.
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10.0, 3.8))
    h = [r["nx"] * r["nz"] for r in rows]
    ax.plot(h, [100 * r["I_T_rel_err"] for r in rows], "o-", label="I_T (∫y²dA)")
    ax.plot(h, [100 * r["I_L_rel_err"] for r in rows], "s-", label="I_L (∫(x−LCF)²dA)")
    ax.axhline(0.0, color="k", lw=0.8, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("grid resolution  nx·nz  (stations × depth rings)")
    ax.set_ylabel("mesh vs closed form  [%]")
    ax.set_title("m2925 route after the merge: waterplane second moments")
    ax.legend()
    ax.grid(alpha=0.3)

    labels = ["LCB", "LCF"]
    ax2.bar(labels, [asym.longitudinal_centre_buoyancy, asym.longitudinal_centre_flotation],
            color="#4c72b0", width=0.5)
    ax2.axhline(0.5 * L, color="r", ls="--", lw=1.0, label="midship L/2")
    ax2.set_ylabel("x from aft end  [m]")
    ax2.set_title("main route: asymmetric fore=0.5 / aft=0.9")
    ax2.legend()
    ax2.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(OUT / "merge_verification.png", dpi=100)
    plt.close(fig)

    (OUT / "merge_verification.json").write_text(json.dumps(data, indent=2, default=str))
    print(json.dumps({k: v for k, v in data.items() if k != "hull_vs_closed_form"}, indent=2, default=str))
    for r in rows:
        print(
            f"nx={r['nx']:4d} nz={r['nz']:3d}  I_T {100*r['I_T_rel_err']:+8.4f}%  "
            f"I_L {100*r['I_L_rel_err']:+8.4f}%  keel_z={r['keel_z_mesh']:.6f}"
        )
    return data


if __name__ == "__main__":
    main()
