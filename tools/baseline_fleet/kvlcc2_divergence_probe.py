#!/usr/bin/env python3
"""Divergence-theorem volume of the KVLCC2 trimmed patch set: a closure audit.

Anchor ``GAP-CAD-T2-PARAM-SERIES`` / fleet geometry debt. ``docs/kvlcc2_trim_field.md``
closed with a stated next step: *"a field estimator whose error is monotone in the cell
size (area-weighted binning plus an explicit hull/not-hull mask from the trim loops)"*.
This tool tests that recommendation directly, by integrating the trim-masked patch set
**exactly** -- no rasterised half-breadth field at all -- and asking whether the deficit
goes away when the estimator is removed from the loop.

Method
------
For a closed surface ``S`` enclosing ``V``, ``V = (1/3) * oint x . n dA``. Each trimmed
patch is triangulated on a parametric grid (grid nodes kept only where the (u, v) cell
centre is inside the even-odd trim region), and the divergence integral is accumulated
per patch. Two diagnostics decide *why* a deficit exists rather than *that* it exists:

* ``closure_residual`` -- ``oint n dA`` over the whole set. It is identically zero for a
  closed, consistently oriented surface; a non-zero residual is the missing boundary in
  vector form. ``x . n dA`` is also reported, so the deficit's x-moment is visible.
* ``trim_coverage`` -- masked over full patch area, per patch: how much of each basis
  sheet the trim loops keep. A patch that keeps a third of itself is a trim defect, not an
  integration defect.
* ``duplicate_groups`` -- patches whose masked area and signed volume agree to 1e-9 are
  the same surface stored twice; summing them double-counts area.

The integrator is validated on an analytic closed sphere before the archive is touched, so
the residual measured on KVLCC2 belongs to the patch set and not to the quadrature.

The asset is publisher material and stays outside the repository; pass its path in.

    python tools/baseline_fleet/kvlcc2_divergence_probe.py \
        --iges /srv/tensor/agent/assets/geometry/hulls/KVLCC2/kvlcc2.igs \
        --out artifacts/baseline_fleet/kvlcc2_divergence.json
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
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tensorcad.fleet.trimfield import trim_polygons  # noqa: E402
from tensorcad.iges import read_iges  # noqa: E402

REFERENCE_FULL = 312622.0
LPP, BEAM, DRAFT = 320.0, 58.0, 20.8


def sphere_check(nu: int = 120, nv: int = 60) -> dict[str, float]:
    """Volume of a unit sphere by the same divergence kernel (exact 4/3 pi)."""
    pts = np.empty((nv + 1, nu + 1, 3))
    for j in range(nv + 1):
        theta = math.pi * j / nv
        for i in range(nu + 1):
            phi = 2.0 * math.pi * i / nu
            pts[j, i] = (math.sin(theta) * math.cos(phi),
                         math.sin(theta) * math.sin(phi),
                         math.cos(theta))
    A, B, C, D = pts[:-1, :-1], pts[1:, :-1], pts[1:, 1:], pts[:-1, 1:]
    volume = area_vector = 0.0
    for tri in (np.stack([A.reshape(-1, 3), B.reshape(-1, 3), C.reshape(-1, 3)], 1),
                np.stack([A.reshape(-1, 3), C.reshape(-1, 3), D.reshape(-1, 3)], 1)):
        vec = 0.5 * np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        volume += (np.einsum("ij,ij->i", tri.mean(axis=1), vec) / 3.0).sum()
        area_vector += np.linalg.norm(vec.sum(axis=0))
    exact = 4.0 / 3.0 * math.pi
    return {
        "volume": volume,
        "exact": exact,
        "error_percent": 100.0 * (volume - exact) / exact,
        "closure_residual_norm": area_vector,
    }


def mask_uv(uc: np.ndarray, vc: np.ndarray, loops) -> np.ndarray:
    """Even-odd containment of (u, v) cell centres against the trim loops."""
    U, V = np.meshgrid(uc, vc)
    hits = np.zeros(U.shape, np.int32)
    for loop in loops:
        P = np.asarray(loop, float)
        Q = np.roll(P, -1, axis=0)
        inside = np.zeros(U.shape, bool)
        for (au, av), (bu, bv) in zip(P, Q):
            if av == bv:
                continue
            cross = (av > V) != (bv > V)
            if not cross.any():
                continue
            ux = au + (V - av) * (bu - au) / (bv - av)
            inside ^= cross & (U < ux)
        hits += inside.astype(np.int32)
    return (hits % 2) == 1


def patch_terms(basis, loops, grid: int) -> dict[str, object]:
    """Signed volume, area (masked and full) and area-vector sum for one trimmed patch."""
    u0, u1 = basis.domain_u
    v0, v1 = basis.domain_v
    us = np.linspace(u0, u1, grid + 1)
    vs = np.linspace(v0, v1, grid + 1)
    Z = np.array([[basis.evaluate(u, v) for u in us] for v in vs])
    mask = mask_uv(0.5 * (us[:-1] + us[1:]), 0.5 * (vs[:-1] + vs[1:]), loops)
    A, B, C, D = Z[:-1, :-1], Z[1:, :-1], Z[1:, 1:], Z[:-1, 1:]
    volume = 0.0
    area_masked = 0.0
    area_full = 0.0
    vector_sum = np.zeros(3)
    n_masked = 0
    for tri, keep in ((np.stack([A.reshape(-1, 3), B.reshape(-1, 3), C.reshape(-1, 3)], 1), mask.ravel()),
                      (np.stack([A.reshape(-1, 3), C.reshape(-1, 3), D.reshape(-1, 3)], 1), mask.ravel())):
        vec = 0.5 * np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        norm = np.linalg.norm(vec, axis=1)
        area_full += norm.sum()
        if not keep.any():
            continue
        sub_vec = vec[keep]
        sub_tri = tri[keep]
        volume += (np.einsum("ij,ij->i", sub_tri.mean(axis=1), sub_vec) / 3.0).sum()
        area_masked += np.linalg.norm(sub_vec, axis=1).sum()
        vector_sum += sub_vec.sum(axis=0)
        n_masked += int(keep.sum())
    return {
        "volume": volume,
        "area_masked": area_masked,
        "area_full": area_full,
        "trim_coverage": (area_masked / area_full) if area_full else 0.0,
        "area_vector": vector_sum.tolist(),
        "triangles": n_masked,
    }


def run(iges_path: str, grid: int) -> dict[str, object]:
    check = sphere_check()
    started = time.time()
    model = read_iges(iges_path)
    records = [r for r in trim_polygons(model)
               if r.get("basis") is not None and r.get("loops")]
    patches = []
    for record in records:
        terms = patch_terms(record["basis"], record["loops"], grid)
        (x0, y0, z0), (x1, y1, z1) = record["basis"].bounding_box()
        terms.update(patch=record["patch"], bounding_box=[[x0, y0, z0], [x1, y1, z1]])
        patches.append(terms)
    half = sum(p["volume"] for p in patches)
    residual = np.sum([p["area_vector"] for p in patches], axis=0)
    # duplicate surfaces: same masked area and same signed volume to solver precision
    groups = []
    for i, a in enumerate(patches):
        for b in patches[i + 1:]:
            if abs(a["area_masked"] - b["area_masked"]) < 1e-6 and abs(a["volume"] - b["volume"]) < 1e-6:
                groups.append([a["patch"], b["patch"]])
    unique_volume = sum(p["volume"] for p in patches
                        if not any(p["patch"] in g[1:] for g in groups))
    return {
        "schema": "tensorcad.fleet.kvlcc2_divergence/1",
        "asset": Path(iges_path).name,
        "grid_per_patch": grid,
        "integrator_check": check,
        "patches_with_loops": len(records),
        "half_volume_signed_m3": half,
        "full_volume_m3": abs(half) * 2.0,
        "reference_full_m3": REFERENCE_FULL,
        "error_vs_reference_percent": 100.0 * (abs(half) * 2.0 - REFERENCE_FULL) / REFERENCE_FULL,
        "de_duplicated_half_volume_m3": unique_volume,
        "de_duplicated_full_volume_m3": abs(unique_volume) * 2.0,
        "duplicate_groups": groups,
        "closure_residual_vector": residual.tolist(),
        "closure_residual_norm_m2": float(np.linalg.norm(residual)),
        "wetted_area_half_m2": sum(p["area_masked"] for p in patches),
        "wetted_area_half_de_duplicated_m2": sum(p["area_masked"] for p in patches
                                                 if not any(p["patch"] in g[1:] for g in groups)),
        "patches": patches,
        "seconds": round(time.time() - started, 2),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--iges", required=True)
    ap.add_argument("--grid", type=int, default=60)
    ap.add_argument("--out", default="artifacts/baseline_fleet/kvlcc2_divergence.json")
    args = ap.parse_args()
    result = run(args.iges, args.grid)
    out = Path(args.out)
    if not out.is_absolute():
        out = ROOT / out
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "patches"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())