#!/usr/bin/env python3
"""Evidence figure for the PR #37 (agent/m2873) semantic-union cure.

Four panels, every number read from a CSV written by the branch's own example
(`examples/transom_stern/run_example.py`) or from the package itself -- nothing
is typed in by hand:

  A  geometry: waterplane plan + body plan of the transom hull (tau = 0.20) and
     the two analytic reference sections it is truncated at (tau = 0)
  B  discretisation: signed volume error vs the analytic reference, log-log,
     for the four truncations of the example's convergence table
  C  reference ladder: the sinc^p power series against two independent
     evaluations (Gamma closed form at tau=0, tanh-sinh quadrature at tau>0)
  D  design space: Cb and LCB/LOA against the immersed transom-area ratio

Run: PYTHONPATH=. python3 docs/m3736/make_figure.py
"""
from __future__ import annotations

import csv
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import tensorcad as tc  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[1]
DATA = REPO / "docs/m2873/outputs/transom_stern"


def rows(name: str) -> list[dict[str, str]]:
    with open(DATA / name, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def panel_geometry(ax):
    tau = 0.20
    mesh, outline = tc.build_transom_mesh(tc.TransomHullParameters(transom_fraction=tau))
    v = np.asarray(mesh.vertices, dtype=float)
    # waterline ring: vertices at z = 0 (the mesh is capped at the design waterline)
    wl = v[np.isclose(v[:, 2], 0.0, atol=1e-9)]
    half = wl[np.argsort(wl[:, 0])]
    starboard = half[half[:, 1] >= 0.0]
    ax.plot(starboard[:, 0], starboard[:, 1], color="#1f4e79", lw=1.6,
            label=f"waterline, tau={tau:.2f}")
    ax.plot(starboard[:, 0], -starboard[:, 1], color="#1f4e79", lw=1.6)
    # transom plate: the planar section at x = 0
    tz = [p[0] for p in outline]
    ty = [p[1] for p in outline]
    ax.plot([0.0] * len(ty), ty, color="#c00000", lw=2.0, label="transom plate (x=0)")
    ax.plot([0.0] * len(ty), [-y for y in ty], color="#c00000", lw=2.0)
    # the pointed hull the family degenerates to at tau = 0, for scale
    p0 = tc.TransomHullParameters(transom_fraction=0.0)
    m0, _ = tc.build_transom_mesh(p0)
    v0 = np.asarray(m0.vertices, dtype=float)
    wl0 = v0[np.isclose(v0[:, 2], 0.0, atol=1e-9)]
    s0 = wl0[(wl0[:, 1] >= 0.0)]
    s0 = s0[np.argsort(s0[:, 0])]
    ax.plot(s0[:, 0], s0[:, 1], color="#7f7f7f", lw=1.0, ls="--",
            label="tau=0 (pointed) for reference")
    ax.axhline(0.0, color="k", lw=0.4)
    ax.set_xlabel("x / m  (aft end → forward)")
    ax.set_ylabel("y / m  (half-breadth)")
    ax.set_title("A  waterplane plan + immersed transom plate")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc="lower right")
    return {"tau": tau, "faces": len(mesh.faces), "outline_points": len(outline),
            "v_mesh": float(np.asarray(tc.analytic_transom_hydrostatics(
                tc.TransomHullParameters(transom_fraction=tau))["displacement_volume"]))}


def panel_convergence(ax):
    data = rows("convergence.csv")
    out = {}
    for tau in sorted({float(r["tau"]) for r in data}):
        rs = [r for r in data if float(r["tau"]) == tau]
        n = np.array([float(r["nx"]) * float(r["nz"]) for r in rs])
        err = np.array([abs(float(r["volume_error"])) for r in rs])
        order = np.polyfit(np.log(n), np.log(err), 1)[0]
        ax.loglog(n, err, "o-", lw=1.2, ms=4,
                  label=f"tau={tau:.2f}  (slope {order:.2f} vs nx·nz)")
        out[f"{tau:.2f}"] = {"order": float(order),
                             "err_finest": float(err[-1])}
    ax.set_xlabel("nx · nz  (mesh size)")
    ax.set_ylabel("|V_mesh − V_analytic| / V_analytic")
    ax.set_title("B  volume error vs analytic reference")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=7)
    return out


def panel_ladder(ax):
    data = rows("reference_ladder.csv")
    ref = []
    for r in data:
        if r["series_vs_closed_form"]:
            ref.append((0.0, float(r["series_at_tau0"]),
                        float(r["gamma_closed_form"])))
        if r["series_vs_quadrature"]:
            ref.append((float(r["tau"]), float(r["series_incomplete"]),
                        float(r["tanh_sinh_incomplete"])))
    ref = np.array(ref)
    rel = np.abs(ref[:, 1] - ref[:, 2]) / np.abs(ref[:, 2])
    tau0 = ref[:, 0] == 0.0
    ax.semilogy(np.where(tau0, 1, 2) + np.random.default_rng(0).uniform(-0.05, 0.05, len(rel)),
                np.maximum(rel, 1e-17), "o", ms=5,
                color="#2e7d32")
    ax.axhline(np.finfo(float).eps, color="k", ls="--", lw=0.8, label="float64 eps")
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"tau=0\nvs Γ closed form\n({int(tau0.sum())} pts)",
                        f"tau>0\nvs tanh-sinh\n({int((~tau0).sum())} pts)"])
    ax.set_ylabel("relative disagreement")
    ax.set_title("C  truncated taper integral: series vs independent")
    ax.grid(alpha=0.3, axis="y")
    ax.legend(fontsize=7)
    return {"points": int(len(rel)), "worst": float(rel.max()),
            "worst_over_eps": float(rel.max() / np.finfo(float).eps)}


def panel_design(ax):
    data = rows("design_table.csv")
    out = {}
    for p in sorted({float(r["exponent"]) for r in data}):
        rs = [r for r in data if float(r["exponent"]) == p]
        ar = [float(r["transom_area_ratio"]) for r in rs]
        cb = [float(r["block_coefficient"]) for r in rs]
        lcb = [float(r["lcb_over_length"]) for r in rs]
        ax.plot(ar, cb, "o-", ms=4, lw=1.2, label=f"Cb, p={p}")
        ax.plot(ar, lcb, "s--", ms=3.5, lw=1.0, alpha=0.75, label=f"LCB/LOA, p={p}")
        out[f"{p}"] = {"cb_range": [min(cb), max(cb)],
                       "lcb_range": [min(lcb), max(lcb)]}
    ax.set_xlabel("immersed transom area ratio  A_t/A_m")
    ax.set_ylabel("Cb  /  LCB/LOA")
    ax.set_title("D  design space dialled by transom area")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, ncol=2)
    return out


def main() -> int:
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.0))
    summary = {
        "geometry": panel_geometry(axes[0][0]),
        "convergence": panel_convergence(axes[0][1]),
        "ladder": panel_ladder(axes[1][0]),
        "design": panel_design(axes[1][1]),
    }
    fig.suptitle("TensorCAD PR #37 cure — transom-stern family after the semantic union "
                 f"(v{tc.__version__})", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = HERE / "merge_evidence.png"
    fig.savefig(out, dpi=150)
    (HERE / "figure_data.json").write_text(json.dumps(summary, indent=2) + "\n",
                                          encoding="utf-8")
    print(json.dumps(summary, indent=2)[:2000])
    print("wrote", out, out.stat().st_size, "bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
