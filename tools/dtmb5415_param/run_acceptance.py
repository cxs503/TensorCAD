#!/usr/bin/env python3
"""Acceptance run for GAP-CAD-C1-5415-PARAM-REBUILD (config 0 vs exact hull, 20-config gate,
gradient chain status).

Predeclared before looking at any result (docs/gap_cad_c1_5415_param_rebuild.md):
  * mesh resolution := median triangle edge length of the config-0 mesh;
  * acceptance 1   := >95 % of the DTMB5415_Hull_exact vertices lie within that resolution of
                      the config-0 mesh (nearest point-to-triangle distance);
  * acceptance 2   := 20 sampled configurations (seeded) pass the geometry gate below;
  * acceptance 3   := analytic dV/dv agrees with a central difference of the same functional
                      to < 1 % relative for all 27 variables (reported honestly otherwise).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from tensorcad.dtmb5415 import DATASET, REFERENCE, WATERLINE_Z, DTMB5415Space

OUT = Path("artifacts/m17472_5415_param")
def _out_dir() -> Path:
    import sys
    if "--out" in sys.argv:
        return Path(sys.argv[sys.argv.index("--out") + 1])
    return OUT
DATASET_MAT = "/root/ship-data-official/dtmb5415/DTMB5415-database.mat"
HOLDOUT_MODULUS = 5          # same split as the fit
N_SAMPLED = 20
SEED = 20261113


def read_stl(path: str) -> np.ndarray:
    raw = open(path, "rb").read()
    n = int(np.frombuffer(raw[80:84], dtype="<u4")[0])
    rec = np.frombuffer(raw[84:84 + 50 * n],
                        dtype=np.dtype([("n", "<3f4"), ("v", "<9f4"), ("a", "<u2")]))
    return rec["v"].reshape(n, 3, 3).astype(np.float64)


def point_tri_dist(p, a, b, c):
    ab, ac, ap = b - a, c - a, p - a
    d1 = (ab * ap).sum(-1); d2 = (ac * ap).sum(-1)
    bp = p - b; d3 = (ab * bp).sum(-1); d4 = (ac * bp).sum(-1)
    cp = p - c; d5 = (ab * cp).sum(-1); d6 = (ac * cp).sum(-1)
    vc = d1 * d4 - d3 * d2; vb = d5 * d2 - d1 * d6; va = d3 * d6 - d5 * d4
    den = va + vb + vc
    with np.errstate(divide="ignore", invalid="ignore"):
        v = vb / den; w = vc / den
    proj = a + v[..., None] * ab + w[..., None] * ac
    e = np.linalg.norm(ap, axis=-1)
    e = np.minimum(e, np.where((d3 >= 0) & (d4 <= d3), np.linalg.norm(bp, axis=-1), np.inf))
    e = np.minimum(e, np.where((d6 >= 0) & (d5 <= d6), np.linalg.norm(cp, axis=-1), np.inf))
    inside = (va >= 0) & (vb >= 0) & (vc >= 0)
    e = np.minimum(e, np.where(inside, np.linalg.norm(p - proj, axis=-1), np.inf))

    def seg(p, a, b):
        ab = b - a
        t = np.clip(((p - a) * ab).sum(-1) / (ab * ab).sum(-1), 0, 1)
        return np.linalg.norm(p - (a + t[..., None] * ab), axis=-1)
    e = np.minimum(e, seg(p, a, b)); e = np.minimum(e, seg(p, a, c)); e = np.minimum(e, seg(p, b, c))
    return e


def dist_to_mesh(verts, tris, points, k=12):
    from scipy.spatial import cKDTree
    tri = verts[tris]
    tree = cKDTree(tri.mean(axis=1))
    _, idx = tree.query(points, k=min(k, len(tri)))
    d = np.full(len(points), np.inf)
    for kk in range(idx.shape[1]):
        f = tri[idx[:, kk]]
        d = np.minimum(d, point_tri_dist(points, f[:, 0], f[:, 1], f[:, 2]))
    return d


def gate(q: dict, base: dict) -> list[dict]:
    checks = [
        ("finite_vertices", q["non_finite_vertices"] == 0, q["non_finite_vertices"], 0),
        ("no_degenerate_faces", q["degenerate_faces"] == 0, q["degenerate_faces"], 0),
        ("no_non_manifold_edges", q["non_manifold_edges"] == 0, q["non_manifold_edges"], 0),
        ("boundary_edges_stable", abs(q["boundary_edges"] - base["boundary_edges"])
         <= 0.25 * base["boundary_edges"], q["boundary_edges"], base["boundary_edges"]),
        ("min_angle_deg", q["min_angle_deg"] >= 20.0, round(q["min_angle_deg"], 2), 20.0),
        ("frac_angle_below_5deg", q["frac_angle_below_5deg"] <= 0.001,
         round(q["frac_angle_below_5deg"], 6), 0.001),
        ("aspect_p99", q["aspect_p99"] <= 10.0, round(q["aspect_p99"], 2), 10.0),
        ("stations_monotone", bool(q["stations_monotone"]), q["stations_monotone"], True),
        ("half_breadth_non_negative", q["half_breadth_min_m"] >= -1e-6,
         round(q["half_breadth_min_m"], 9), -1e-6),
        ("edge_max_m", q["edge_max_m"] <= 8.0, round(q["edge_max_m"], 3), 8.0),
    ]
    return [{"check": c, "value": v, "threshold": t, "ok": bool(ok)} for c, ok, v, t in checks]


def main() -> None:
    global OUT
    OUT = _out_dir()
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    space = DTMB5415Space()
    v0 = space.zero()
    q0 = space.quality(v0)
    verts0 = space.vertices(v0)
    tris = space.faces()
    resolution = float(q0["edge_median_m"])          # predeclared mesh resolution
    hydro0 = space.hydrostatics(v0)
    vol0 = space.volume(v0)

    # ---- acceptance 1: config 0 vs DTMB5415_Hull_exact --------------------------
    exact = read_stl(REFERENCE["path"])
    exact_pts = np.unique(exact.reshape(-1, 3), axis=0)
    d0 = dist_to_mesh(verts0, tris, exact_pts)
    acc1 = {
        "reference": REFERENCE["path"],
        "reference_sha256": REFERENCE["sha256"],
        "n_reference_points": int(len(exact_pts)),
        "mesh_resolution_m": resolution,
        "median_m": float(np.median(d0)), "p95_m": float(np.percentile(d0, 95)),
        "max_m": float(d0.max()),
        "frac_within_resolution": float((d0 < resolution).mean()),
        "frac_within_half_resolution": float((d0 < 0.5 * resolution).mean()),
        "threshold": 0.95,
        "pass": bool((d0 < resolution).mean() > 0.95),
    }

    # ---- acceptance 2: sampled configurations through the gate ------------------
    import scipy.io as sio
    d = sio.loadmat(DATASET_MAT, variable_names=["variablesMatrix", "xMatrix", "yMatrix",
                                                 "zMatrix"])
    V = np.asarray(d["variablesMatrix"]).T
    clouds = np.stack([np.asarray(d["xMatrix"]).T, np.asarray(d["yMatrix"]).T,
                       np.asarray(d["zMatrix"]).T], axis=2)
    hold = np.flatnonzero(np.arange(V.shape[0]) % HOLDOUT_MODULUS == 0)
    rng = np.random.default_rng(SEED)
    pick = np.sort(rng.choice(hold, size=N_SAMPLED, replace=False))
    rows, gate_rows = [], []
    for k in pick:
        v = V[k]
        pred = space.cloud(v)
        true = clouds[k]
        err = np.linalg.norm(pred - true, axis=1)
        q = space.quality(v)
        h = space.hydrostatics(v)
        g = gate(q, q0)
        ok = all(c["ok"] for c in g)
        rows.append({
            "config_index": int(k),
            "rms_point_m": float(np.sqrt((err ** 2).mean())),
            "max_point_m": float(err.max()),
            "rms_pct_beam": float(100 * np.sqrt((err ** 2).mean()) / space.beam_full_m),
            "rms_pct_half_beam": float(100 * np.sqrt((err ** 2).mean()) / space.beam_m),
            "half_breadth_raw_min_m": q["half_breadth_raw_min_m"],
            "half_breadth_min_m": q["half_breadth_min_m"],
            "half_breadth_clamped_vertices": q["half_breadth_clamped_vertices"],
            "edge_median_m": q["edge_median_m"], "min_angle_deg": q["min_angle_deg"],
            "aspect_p99": q["aspect_p99"], "boundary_edges": q["boundary_edges"],
            "degenerate_faces": q["degenerate_faces"],
            "non_manifold_edges": q["non_manifold_edges"],
            "stations_monotone": q["stations_monotone"],
            "volume_m3": space.volume(v),
            "displacement_m3": h["displacement_volume_m3"], "cb": h["cb"],
            "beam_wl_m": h["beam_waterline_m"], "gate_pass": ok,
            "checks_failed": [c["check"] for c in g if not c["ok"]],
        })
        gate_rows.append({"config_index": int(k), "checks": g, "pass": ok})
    # a second, harder sample: corners of the box [-1, 1]^27 sampled uniformly
    uni_rows = []
    for n in range(N_SAMPLED):
        v = rng.uniform(-1.0, 1.0, size=27)
        q = space.quality(v)
        h = space.hydrostatics(v)
        g = gate(q, q0)
        uni_rows.append({"sample": n, "gate_pass": all(c["ok"] for c in g),
                         "volume_m3": space.volume(v), "cb": h["cb"],
                         "displacement_m3": h["displacement_volume_m3"],
                         "min_angle_deg": q["min_angle_deg"],
                         "checks_failed": [c["check"] for c in g if not c["ok"]]})

    # ---- acceptance 3: gradient chain ------------------------------------------
    hstep = 1e-5
    fd_vol = np.zeros(27)
    for k in range(27):
        e = np.zeros(27); e[k] = hstep
        fd_vol[k] = (space.volume(e) - space.volume(-e)) / (2 * hstep)
    ana_vol = space.volume_gradient(v0)
    rel = np.abs(ana_vol - fd_vol) / np.maximum(np.abs(fd_vol), 1e-12)
    jac = space.vertex_jacobian(v0)
    jac_err = []
    for k in (0, 5, 13, 26):
        e = np.zeros(27); e[k] = hstep
        fd_v = (space.vertices(e) - space.vertices(-e)) / (2 * hstep)
        jac_err.append(float(np.abs(fd_v - jac[:, :, k]).max()))
    grad = {
        "vertex_jacobian_shape": list(jac.shape),
        "vertex_jacobian_fd_max_abs_err_by_mode": jac_err,
        "vertex_jacobian_exact": bool(max(jac_err) < 1e-6),
        "volume_gradient_max_abs_diff": float(np.abs(ana_vol - fd_vol).max()),
        "volume_gradient_max_rel_err": float(rel.max()),
        "volume_gradient_pass": bool(rel.max() < 0.01),
        "volume_gradient_analytic": [float(x) for x in ana_vol],
        "volume_gradient_fd": [float(x) for x in fd_vol],
        "waterline_clip_differentiated": False,
    }

    result = {
        "anchor": "GAP-CAD-C1-5415-PARAM-REBUILD",
        "dataset": {"source_sha256": DATASET["source_sha256"],
                    "n_configurations": DATASET["n_configurations"],
                    "n_design_variables": DATASET["n_design_variables"]},
        "config0": {"volume_shell_m3": vol0, "hydrostatics": hydro0, "quality": q0},
        "acceptance": {"config0_vs_exact": acc1,
                       "sampled_configs_gate": {
                           "n": len(rows), "seed": SEED, "holdout_modulus": HOLDOUT_MODULUS,
                           "pass": bool(all(r["gate_pass"] for r in rows)),
                           "n_pass": int(sum(r["gate_pass"] for r in rows)), "rows": rows},
                       "uniform_samples_gate": {
                           "n": len(uni_rows), "seed": SEED,
                           "n_pass": int(sum(r["gate_pass"] for r in uni_rows)),
                           "rows": uni_rows},
                       "gradient_chain": grad},
        "criteria": {
            "acceptance_1_config0_vs_exact": {
                "criterion": ">95% of DTMB5415_Hull_exact vertices within the config-0 mesh "
                             "resolution",
                "value": acc1["frac_within_resolution"], "threshold": 0.95,
                "pass": acc1["pass"]},
            "acceptance_2_sampled_gate": {
                "criterion": "20 sampled dataset configurations pass the geometry quality gate",
                "value": f'{sum(r["gate_pass"] for r in rows)}/{len(rows)}',
                "threshold": "20/20",
                "pass": bool(all(r["gate_pass"] for r in rows))},
            "acceptance_3_gradient_chain": {
                "criterion": "analytic dV/dv within 1% of central difference, all 27 variables",
                "value": grad["volume_gradient_max_rel_err"], "threshold": 0.01,
                "pass": grad["volume_gradient_pass"],
                "vertex_jacobian_exact": grad["vertex_jacobian_exact"]},
        },
        "model": {
            "shipped_holdout_rms_m": space.metrics.get(
                "M3_pca_quadratic_shipped", {}).get("holdout", {}).get("rms_point_m"),
            "beam_full_m": space.beam_full_m, "beam_half_m": space.beam_m,
            "note": "geometry quality gate uses the full-beam percentage; the round-0 report "
                    "divided by the half breadth and double counted the linear modal block",
        },
        "wall_clock_s": None,
    }
    (OUT / "acceptance.json").write_text(json.dumps(result, indent=2))
    with open(OUT / "sampled_configs.csv", "w") as fh:
        keys = list(rows[0])
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(f"{r[k]:.6g}" if isinstance(r[k], float)
                              else "|".join(r[k]) if isinstance(r[k], list) else str(r[k])
                              for k in keys) + "\n")
    with open(OUT / "gate_samples.csv", "w") as fh:
        fh.write("sample,gate_pass,volume_m3,cb,displacement_m3,min_angle_deg,failed\n")
        for r in uni_rows:
            fh.write(f"{r['sample']},{r['gate_pass']},{r['volume_m3']:.6g},{r['cb']:.6g},"
                     f"{r['displacement_m3']:.6g},{r['min_angle_deg']:.6g},"
                     f"{'|'.join(r['checks_failed'])}\n")

    # ---- figures ---------------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].hist(d0, bins=80, color="#2b6cb0")
    ax[0].axvline(resolution, color="crimson", label=f"mesh resolution {resolution:.3f} m")
    ax[0].set_yscale("log"); ax[0].set_xlabel("distance to config-0 mesh (m)")
    ax[0].set_ylabel("count (log)"); ax[0].legend()
    ax[0].set_title("config 0 vs DTMB5415_Hull_exact")
    ax[1].plot(np.sort(d0), np.linspace(0, 1, d0.size))
    ax[1].axvline(resolution, color="crimson"); ax[1].axhline(0.95, ls=":", color="k")
    ax[1].set_xlabel("distance (m)"); ax[1].set_ylabel("cumulative fraction")
    ax[1].set_title(f"{(d0 < resolution).mean() * 100:.2f}% within resolution")
    fig.tight_layout(); fig.savefig(OUT / "fig_config0_distance.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    idxs = [0, 40, 80, 120]
    g0 = space.station_grid(v0)
    for i in idxs:
        ax[0].plot(g0[0][:, i, 1], g0[0][:, i, 2], "-", label=f"fore station {i}")
    ax[0].set_xlabel("y (m)"); ax[0].set_ylabel("z (m)"); ax[0].legend(fontsize=7)
    ax[0].set_title("parametric sections at v = 0")
    vols = [r["volume_m3"] for r in rows]
    ax[1].plot([r["config_index"] for r in rows], vols, "o-")
    ax[1].axhline(vol0, color="crimson", ls="--", label="config 0")
    ax[1].set_xlabel("dataset configuration index"); ax[1].set_ylabel("shell volume (m$^3$)")
    ax[1].legend(); ax[1].set_title("20 sampled configurations")
    fig.tight_layout(); fig.savefig(OUT / "fig_sections_and_volume.png", dpi=140)
    plt.close(fig)

    # symmetry-plane restoration: the raw surrogate profile vs the feasible geometry
    worst = min(rows, key=lambda r: r["half_breadth_raw_min_m"])
    vw = V[worst["config_index"]]
    raw_w = space.raw_station_grid(vw)
    fea_w = space.station_grid(vw)
    patch, _, col = np.unravel_index(np.argmin(raw_w[..., 1]), raw_w[..., 1].shape)
    prof_raw = raw_w[patch, :, col, 1]
    prof_fea = fea_w[patch, :, col, 1]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(prof_raw, "o-", ms=3, label=f"raw surrogate (min {prof_raw.min():.3f} m)")
    ax[0].plot(prof_fea, "s-", ms=3, label="after y >= 0 projection")
    ax[0].axhline(0.0, color="k", ls=":", lw=1)
    ax[0].set_xlabel(f"girth index (patch {patch}, station {col})")
    ax[0].set_ylabel("half breadth y (m)")
    ax[0].set_title(f"config {worst['config_index']}: symmetry-plane restoration")
    ax[0].legend(fontsize=8)
    mins = [r["half_breadth_raw_min_m"] for r in rows]
    ax[1].bar(range(len(rows)), mins, color="#c05621")
    ax[1].axhline(0.0, color="k", ls=":")
    ax[1].axhline(-1e-6, color="crimson", ls="--", label="gate tolerance -1e-6 m")
    ax[1].set_xlabel("sampled configuration (20)"); ax[1].set_ylabel("raw min half breadth (m)")
    ax[1].set_title("surrogate overshoot, repaired to 0 by construction")
    ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "fig_symmetry_plane_restoration.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.6, 5.0))
    lim = max(np.abs(ana_vol).max(), np.abs(fd_vol).max()) * 1.1
    ax.plot([-lim, lim], [-lim, lim], "k:", lw=1)
    ax.plot(fd_vol, ana_vol, "o", color="#c53030")
    ax.set_xlabel("central difference dV/dv (m$^3$ per unit v)")
    ax.set_ylabel("analytic dV/dv")
    ax.set_title(f"volume gradient: max rel err {rel.max():.2f} (NOT exact -> reported)")
    fig.tight_layout(); fig.savefig(OUT / "fig_gradient_parity.png", dpi=140)
    plt.close(fig)

    result["wall_clock_s"] = round(time.time() - t0, 1)
    (OUT / "acceptance.json").write_text(json.dumps(result, indent=2))
    print(json.dumps({"config0": acc1["frac_within_resolution"],
                      "median_dist": acc1["median_m"], "resolution": resolution,
                      "gate_pass": result["acceptance"]["sampled_configs_gate"]["n_pass"],
                      "gate_n": len(rows),
                      "uniform_pass": result["acceptance"]["uniform_samples_gate"]["n_pass"],
                      "grad_rel": grad["volume_gradient_max_rel_err"],
                      "vertex_jac_exact": grad["vertex_jacobian_exact"],
                      "hydro0": hydro0, "wall_s": result["wall_clock_s"]}, indent=2))


if __name__ == "__main__":
    main()
