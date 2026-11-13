#!/usr/bin/env python3
"""Cross-check the *shipped* reduced model (``docs/data/dtmb5415/design_space.npz``)
against the design database itself, and re-derive the holdout accuracy independently.

Why this exists: ``fit_design_space.py`` reports residuals for the model *it just fitted*
in memory; ``run_acceptance.py`` measures ``DTMB5415Space.cloud`` on the *serialised*
artifact.  Round 0 of this anchor reported two different numbers for the same quantity
(0.252 m holdout RMS from the fit, 0.786 m per-configuration RMS from the shipped model),
so the artifact -> module round trip is audited here on its own evidence, not on trust.

Outputs (``--out``):
  * ``shipped_model_check.json``  — round-trip identity, holdout residuals, worst configs
  * ``shipped_model_check.csv``   — per-configuration residuals of the shipped model

Usage::

    PYTHONPATH=. python3 tools/dtmb5415_param/verify_shipped_model.py \
        --out artifacts/m17472_5415_param --n-configs 400
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import scipy.io as sio

from tensorcad.dtmb5415 import (DATASET, GRID, MESH_RESOLUTION, _DATA_DIR,
                                DTMB5415Space)

MAT_KEYS = ["variablesMatrix", "xMatrix", "yMatrix", "zMatrix"]
HOLDOUT_MODULUS = 5          # predeclared split, identical to fit_design_space.py


def load_clouds(mat_path: str):
    d = sio.loadmat(mat_path, variable_names=MAT_KEYS)
    v = np.asarray(d["variablesMatrix"].T, dtype=np.float64)
    clouds = np.stack([np.asarray(d["xMatrix"]).T, np.asarray(d["yMatrix"]).T,
                       np.asarray(d["zMatrix"]).T], axis=2)
    return v, clouds


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mat", default=DATASET["source_file"])
    ap.add_argument("--out", default="artifacts/m17472_5415_param")
    ap.add_argument("--n-configs", type=int, default=400,
                    help="holdout configurations re-measured (0 = all 1801)")
    ap.add_argument("--seed", type=int, default=20261113)
    ap.add_argument("--fixture", default=None,
                    help="write a small (design vector, true cloud) fixture for the unit tests "
                         "so CI never needs the 1.0 GB design database")
    ap.add_argument("--rank", type=int, default=40,
                    help="PCA rank used for the independent refit (must match the fit)")
    args = ap.parse_args()

    t0 = time.time()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    space = DTMB5415Space()
    v, clouds = load_clouds(args.mat)
    n_cfg, n_pts, _ = clouds.shape
    base = clouds[0]
    disp = (clouds - base).reshape(n_cfg, n_pts * 3)
    idx = np.arange(n_cfg)
    hold = idx % HOLDOUT_MODULUS == 0
    fit = ~hold

    # ---- 1. round trip: shipped arrays vs the module's own prediction -------------
    # The module is a pure function of the npz; check that the npz really encodes the
    # quoted (rank-40, quadratic) model rather than something else.
    npz = np.load(Path(_DATA_DIR) / "design_space.npz")
    round_trip = {
        "npz_keys": sorted(npz.files),
        "phi_shape": list(npz["phi"].shape),
        "phi_dtype": str(npz["phi"].dtype),
        "c1_shape": list(npz["c1"].shape),
        "c2_shape": list(npz["c2"].shape),
        "base_cloud_shape": list(npz["base_cloud"].shape),
        "rank_matches_metrics": int(npz["phi"].shape[1]) == int(space.metrics.get("pca_rank", -1)),
        "n_points_matches_grid": int(npz["base_cloud"].shape[0]) ==
        GRID["n_girth"] * GRID["n_station"] * GRID["n_patches"],
        "beam_m": float(npz["beam_m"]),
        "base_cloud_equals_dataset_column_0": float(np.abs(npz["base_cloud"] - base).max()),
        "phi_float32_rounding_m": None,      # filled below
    }
    # how much error the float32 storage of phi alone can inject, measured not assumed
    phi64 = np.asarray(npz["phi"], dtype=np.float64)
    phi32 = np.asarray(npz["phi"], dtype=np.float64).astype(np.float32).astype(np.float64)
    amp_scale = max(float(np.abs(space._amplitudes(v[k])).max()) for k in (1, 2, 3))
    round_trip["phi_float32_rounding_m"] = float(np.abs(phi64 - phi32).max() * amp_scale)

    # ---- 2. holdout accuracy of the *shipped* model ------------------------------
    rng = np.random.default_rng(args.seed)
    hold_idx = idx[hold]
    pick = (np.sort(rng.choice(hold_idx, size=min(args.n_configs, hold_idx.size),
                               replace=False))
            if args.n_configs else hold_idx)
    pred = np.stack([space.cloud(v[k]) for k in pick])
    err = np.linalg.norm(pred - clouds[pick], axis=2)
    per_cfg_rms = np.sqrt((err ** 2).mean(axis=1))
    per_cfg_max = err.max(axis=1)
    beam = float(base[:, 1].max() - base[:, 1].min())

    # ---- 3. independent refit of the same model family (predeclared) -------------
    # Identical estimator to fit_design_space.py, recomputed here so the ship-side number
    # can never be an accident of a stale array.
    def quad_feat(x):
        iu, ju = np.triu_indices(x.shape[1])
        prod = (x[:, :, None] * x[:, None, :])[:, iu, ju]
        return np.concatenate([x, prod], axis=1)

    X2 = quad_feat(v)
    Df = disp[fit]
    Om = rng.standard_normal((Df.shape[1], args.rank + 20))
    Q, _ = np.linalg.qr(Df @ Om)
    for _ in range(2):
        Q, _ = np.linalg.qr(Df.T @ Q)
        Q, _ = np.linalg.qr(Df @ Q)
    Q, _ = np.linalg.qr(Df.T @ Q)
    Ub, s, _ = np.linalg.svd(Q.T @ Df.T, full_matrices=False)
    Phi = Q @ Ub
    Phi = Phi[:, :args.rank]
    C2 = np.linalg.lstsq(X2[fit], Df @ Phi, rcond=None)[0]
    refit_pred = np.stack([(X2[k] @ C2) @ Phi.T for k in pick]).reshape(-1, n_pts, 3) + base
    refit_err = np.linalg.norm(refit_pred - clouds[pick], axis=2)

    # ---- 4. round-0 defect baseline: amplitudes with the linear block counted twice --
    # The shipped module used to return ``c1 v + c2 [v, vech(vv^T)]``; ``c2`` already
    # contains the linear columns, so this doubles the linear displacement.  Measured here
    # so the before/after numbers in the report are reproducible from the audit alone.
    double = np.stack([space.base_cloud + (space.phi @ (
        space.c1 @ v[k] + space.c2 @ np.concatenate(
            [v[k], space._quadratic(v[k])]))).reshape(-1, 3) for k in pick])
    dbl_err = np.linalg.norm(double - clouds[pick], axis=2)
    dbl_rms = float(np.sqrt((dbl_err ** 2).mean()))

    worst = np.argsort(-per_cfg_rms)[:10]
    result = {
        "anchor": "GAP-CAD-C1-5415-PARAM-REBUILD",
        "purpose": "audit of the shipped reduced model against the design database",
        "dataset": {"source_sha256": DATASET["source_sha256"],
                    "source_file": DATASET["source_file"],
                    "n_configurations": n_cfg, "n_points_per_configuration": n_pts},
        "shipped_artifact": {
            "npz": "docs/data/dtmb5415/design_space.npz",
            "metrics": space.metrics.get("split", {}),
            "beam_half_breadth_m": float(space.beam_m),
            "beam_full_m": beam,
            "note": "DTMB5415Space.beam_m is the HALF breadth of the base cloud; the full "
                    "beam is 2x that (see 'beam_full_m')",
            "mesh_resolution": MESH_RESOLUTION,
        },
        "round_trip": round_trip,
        "holdout_shipped": {
            "n_configs": int(len(pick)), "seed": args.seed,
            "holdout_modulus": HOLDOUT_MODULUS,
            "rms_point_m": float(np.sqrt((err ** 2).mean())),
            "rms_point_pct_full_beam": float(100 * np.sqrt((err ** 2).mean()) / beam),
            "median_config_rms_m": float(np.median(per_cfg_rms)),
            "p95_config_rms_m": float(np.percentile(per_cfg_rms, 95)),
            "worst_config_rms_m": float(per_cfg_rms.max()),
            "median_config_max_m": float(np.median(per_cfg_max)),
            "p95_config_max_m": float(np.percentile(per_cfg_max, 95)),
            "worst_config_max_m": float(per_cfg_max.max()),
        },
        "holdout_independent_refit": {
            "rank": args.rank,
            "rms_point_m": float(np.sqrt((refit_err ** 2).mean())),
            "rms_point_pct_full_beam": float(100 * np.sqrt((refit_err ** 2).mean()) / beam),
            "agrees_with_shipped_pct": float(
                100 * abs(np.sqrt((refit_err ** 2).mean()) - np.sqrt((err ** 2).mean()))
                / max(np.sqrt((err ** 2).mean()), 1e-12)),
        },
        "round0_defect_baseline": {
            "description": "amplitudes with the linear block counted twice "
                           "(c1 v + c2 [v, vech]) as shipped in round 0",
            "rms_point_m": dbl_rms,
            "rms_point_pct_full_beam": float(100 * dbl_rms / beam),
            "ratio_vs_fixed": float(dbl_rms / max(np.sqrt((err ** 2).mean()), 1e-12)),
        },
        "worst_configs": [
            {"config_index": int(pick[j]), "rms_m": float(per_cfg_rms[j]),
             "max_m": float(per_cfg_max[j]),
             "rms_pct_full_beam": float(100 * per_cfg_rms[j] / beam)} for j in worst],
        "wall_clock_s": None,
    }

    with open(out / "shipped_model_check.csv", "w") as fh:
        fh.write("config_index,rms_point_m,max_point_m,rms_pct_full_beam\n")
        for k, r, mx in zip(pick, per_cfg_rms, per_cfg_max):
            fh.write(f"{int(k)},{r:.6f},{mx:.6f},{100 * r / beam:.4f}\n")
    if args.fixture:
        fixed = np.argsort(per_cfg_rms)          # sorted so the fixture is fingerprintable
        fx = Path(args.fixture)
        fx.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            fx,
            config_index=pick[fixed].astype(np.int64),
            design_vector=v[pick[fixed]].astype(np.float32),
            true_cloud=clouds[pick[fixed]].astype(np.float32),
            per_config_rms_m=per_cfg_rms[fixed].astype(np.float64),
            source_sha256=np.array(DATASET["source_sha256"]),
        )
        result["fixture"] = {"path": str(fx), "n_configs": int(len(pick)),
                             "bytes": fx.stat().st_size,
                             "note": "true clouds + design vectors resampled from the design "
                                     "database, float32; lets the unit tests check the shipped "
                                     "model against the database without the 1.0 GB file"}

    result["wall_clock_s"] = round(time.time() - t0, 1)
    (out / "shipped_model_check.json").write_text(json.dumps(result, indent=2))
    print(json.dumps({k: result[k] for k in ("round_trip", "holdout_shipped",
                                             "holdout_independent_refit",
                                             "round0_defect_baseline", "wall_clock_s")},
                     indent=2)[:4000])


if __name__ == "__main__":
    main()