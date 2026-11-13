#!/usr/bin/env python3
"""Fit the 27-DV design space of the Zenodo DTMB 5415 database (9001 x 27).

The Zenodo dataset ("Design-Space Dimensionality Reduction Benchmark Dataset - DTMB 5415
Hull", Serani et al.) ships 9001 hull configurations, each as 2250 surface points
(6750 geometric DOF) plus the 27 normalised design variables that generated it.  What it
does *not* ship is the parametric model itself, so the map ``v -> geometry`` has to be
recovered.  This script recovers it from the data and stores a compact, differentiable
surrogate that the TensorCAD module :mod:`tensorcad.dtmb5415` consumes.

Model (predeclared, see docs/gap_cad_c1_5415_param_rebuild.md):

    d(v) = Phi . ( C1 v + C2 (v (x) v) )        # displacement of the 2250 base points

* ``d`` is the lagrangian displacement of the baseline point cloud (config 0), i.e.
  ``d_k = P_k - P_0`` for configuration ``k``; ``v = 0`` maps to the exact baseline hull.
* ``Phi`` (6750 x K) are the leading displacement PCA modes of the *fit* split.
* ``C1`` (K x 27) / ``C2`` (K x 378) are least-squares coefficients on the training split.
* Rationale for the quadratic lift: the dataset README states the geometry comes from a
  *nonlinear* parametric model, so a purely linear map is expected to leave a residual;
  the quadratic term measures how much of that residual is captured.

Split (predeclared): configuration index ``k % 5 == 0`` -> holdout, the rest -> fit.
Configurations are independent samples (Monte Carlo design of experiments), so this is a
plain holdout, not a time series split.

Outputs (into --out):
  design_space.npz        Phi, C1, C2, baseline cloud, feature/DOF layout, split spec
  design_space_metrics.json   training / holdout residuals for M0/M1/M2, PCA spectrum
  design_space_holdout.csv    per-configuration residuals for the holdout split
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np
import scipy.io as sio

#: Dataset matrices we need (the .mat is ~1 GiB, so never load the rest).
MAT_KEYS = ["variablesMatrix", "xMatrix", "yMatrix", "zMatrix"]

#: Predeclared PCA rank of the shipped displacement basis.
DEFAULT_RANK = 40

#: Predeclared holdout rule: every 5th configuration (index % 5 == 0).
HOLDOUT_MODULUS = 5

#: Reference model-scale table quoted by the dataset README (Zenodo record), used for the
#: geometry sanity check (Lpp, beam, draft, displacement at model scale).
README_MODEL_SCALE = {
    "Lpp_m": 5.720,
    "beam_m": 0.760,
    "draft_m": 0.248,
    "displacement_t": 0.549,
    "rho_kg_m3": 998.5,
}
README_FULL_SCALE_LPP_M = 142.0


def sha256(path: str | os.PathLike, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk)
            if not block:
                break
            h.update(block)
    return h.hexdigest()


def load_clouds(mat_path: str):
    """Return (variables (9001,27), clouds (9001,2250,3))."""
    d = sio.loadmat(mat_path, variable_names=MAT_KEYS)
    v = np.asarray(d["variablesMatrix"].T, dtype=np.float64)          # (9001, 27)
    x = np.asarray(d["xMatrix"], dtype=np.float64)                    # (2250, 9001)
    y = np.asarray(d["yMatrix"], dtype=np.float64)
    z = np.asarray(d["zMatrix"], dtype=np.float64)
    clouds = np.stack([x.T, y.T, z.T], axis=2)                        # (9001, 2250, 3)
    return v, clouds


def quadratic_features(v: np.ndarray) -> np.ndarray:
    """[v, vech(v v^T)] with the diagonal kept once (27 + 378 = 405 columns)."""
    n, m = v.shape
    iu, ju = np.triu_indices(m)
    prod = (v[:, :, None] * v[:, None, :])[:, iu, ju]
    return np.concatenate([v, prod], axis=1)


def _fit(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Least squares without intercept (v = 0 must map to the baseline exactly)."""
    XtX = X.T @ X
    XtY = X.T @ Y
    return np.linalg.solve(XtX + 1e-12 * np.eye(XtX.shape[0]), XtY)


def _residual_stats(pred: np.ndarray, true: np.ndarray, beam: float) -> dict:
    err = np.linalg.norm(pred - true, axis=2)             # (n_cfg, 2250)
    per = err.max(axis=1)
    stats = {
        "n_configs": int(err.shape[0]),
        "rms_point_m": float(np.sqrt((err ** 2).mean())),
        "median_config_max_m": float(np.median(per)),
        "p95_config_max_m": float(np.percentile(per, 95)),
        "worst_config_max_m": float(per.max()),
        "rms_point_pct_beam": float(100.0 * np.sqrt((err ** 2).mean()) / beam),
        "p95_config_max_pct_beam": float(100.0 * np.percentile(per, 95) / beam),
    }
    return stats


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mat", default="/root/ship-data-official/dtmb5415/DTMB5415-database.mat")
    ap.add_argument("--out", default="runs/m17352_5415_param")
    ap.add_argument("--rank", type=int, default=DEFAULT_RANK)
    args = ap.parse_args()

    t0 = time.time()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    v, clouds = load_clouds(args.mat)
    n_cfg, n_pts, _ = clouds.shape
    base = clouds[0]
    disp = (clouds - base).reshape(n_cfg, n_pts * 3)      # (9001, 6750)

    idx = np.arange(n_cfg)
    hold = idx % HOLDOUT_MODULUS == 0
    fit = ~hold
    X1 = v
    X2 = quadratic_features(v)

    print(f"[fit] configs={n_cfg} fit={fit.sum()} holdout={hold.sum()} rank={args.rank}",
          flush=True)

    # ---- baseline identity: what does v = 0 give? --------------------------------
    base_identity = float(np.abs(v[0]).max())

    # ---- full-space regression (upper bound of the family) -----------------------
    A1 = _fit(X1[fit], disp[fit])
    A2 = _fit(X2[fit], disp[fit])
    pred1_hold = (X1[hold] @ A1).reshape(-1, n_pts, 3) + base
    pred2_hold = (X2[hold] @ A2).reshape(-1, n_pts, 3) + base
    pred1_fit = (X1[fit] @ A1).reshape(-1, n_pts, 3) + base
    pred2_fit = (X2[fit] @ A2).reshape(-1, n_pts, 3) + base
    beam = float(base[:, 1].max() - base[:, 1].min())

    # ---- PCA of the displacement field (fit split only) --------------------------
    Df = disp[fit]
    rng = np.random.default_rng(20261113)
    Omega = rng.standard_normal((Df.shape[1], args.rank + 20))
    Yr = Df @ Omega
    Q, _ = np.linalg.qr(Yr)                              # sample space (7200, k)
    for _ in range(2):                                   # power iterations
        Q, _ = np.linalg.qr(Df.T @ Q)                    # -> DOF space (6750, k)
        Q, _ = np.linalg.qr(Df @ Q)                      # -> sample space (7200, k)
    Q, _ = np.linalg.qr(Df.T @ Q)                        # finish in DOF space
    B = Q.T @ Df.T                                        # (k, n_samples)
    Ub, s, _ = np.linalg.svd(B, full_matrices=False)
    Phi = Q @ Ub
    if Phi.shape[1] > args.rank:
        Phi = Phi[:, :args.rank]
        s = s[:args.rank]
    var = float((Df ** 2).sum())
    var_captured = float((s ** 2).sum() / var)
    # deterministic sign convention: largest-|value| entry of each mode positive
    for j in range(Phi.shape[1]):
        if Phi[np.argmax(np.abs(Phi[:, j])), j] < 0:
            Phi[:, j] = -Phi[:, j]

    # ---- reduced model: project, regress in modal coordinates -------------------
    a_fit = Df @ Phi
    C1 = np.linalg.lstsq(X1[fit], a_fit, rcond=None)[0]
    C2 = np.linalg.lstsq(X2[fit], a_fit, rcond=None)[0]
    pred_red_hold = ((X2[hold] @ C2) @ Phi.T).reshape(-1, n_pts, 3) + base
    pred_red_fit = ((X2[fit] @ C2) @ Phi.T).reshape(-1, n_pts, 3) + base

    stats = {
        "baseline_v_max_abs": base_identity,
        "pca_rank": int(Phi.shape[1]),
        "pca_variance_captured": var_captured,
        "pca_singular_values": [float(x) for x in s],
        "M0_baseline_only": {
            "holdout": _residual_stats(np.broadcast_to(base, pred1_hold.shape), clouds[hold], beam),
        },
        "M1_linear_full": {
            "fit": _residual_stats(pred1_fit, clouds[fit], beam),
            "holdout": _residual_stats(pred1_hold, clouds[hold], beam),
        },
        "M2_quadratic_full": {
            "fit": _residual_stats(pred2_fit, clouds[fit], beam),
            "holdout": _residual_stats(pred2_hold, clouds[hold], beam),
        },
        "M3_pca_quadratic_shipped": {
            "fit": _residual_stats(pred_red_fit, clouds[fit], beam),
            "holdout": _residual_stats(pred_red_hold, clouds[hold], beam),
        },
        "beam_m": beam,
        "split": {"holdout_modulus": HOLDOUT_MODULUS,
                  "n_fit": int(fit.sum()), "n_holdout": int(hold.sum())},
        "elapsed_s": None,
    }

    # ---- per-config holdout residuals ------------------------------------------
    err = np.linalg.norm(pred_red_hold - clouds[hold], axis=2)
    with open(out / "design_space_holdout.csv", "w") as fh:
        fh.write("config_index,rms_point_m,max_point_m\n")
        for k, e in zip(idx[hold], err):
            fh.write(f"{int(k)},{float(np.sqrt((e ** 2).mean())):.6f},{float(e.max()):.6f}\n")

    np.savez_compressed(
        out / "design_space.npz",
        phi=Phi.astype(np.float32),
        c1=C1.astype(np.float64),
        c2=C2.astype(np.float64),
        base_cloud=base.astype(np.float64),
        beam_m=np.float64(beam),
        n_points=np.int64(n_pts),
    )
    (out / "design_space_metrics.json").write_text(json.dumps(stats, indent=2))

    stats["elapsed_s"] = round(time.time() - t0, 1)
    stats["source_mat"] = {"path": args.mat,
                           "size_bytes": os.path.getsize(args.mat),
                           "sha256": sha256(args.mat)}
    (out / "design_space_metrics.json").write_text(json.dumps(stats, indent=2))
    print(json.dumps({k: stats[k] for k in
                      ("pca_rank", "pca_variance_captured", "M1_linear_full",
                       "M2_quadratic_full", "M3_pca_quadratic_shipped")}, indent=2))


if __name__ == "__main__":
    main()