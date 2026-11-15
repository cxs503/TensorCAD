"""Measure the cost of the type-128/126 span-index defect, legacy vs fixed.

Reproduce:
    python tools/baseline_fleet/iges_span_defect.py --out runs/m21680_span_defect/evidence.json

The legacy arm re-installs the pre-fix :func:`tensorcad.iges._find_span` (the one that
compared against ``knots[n]`` instead of ``knots[n+1]`` and returned ``n-1`` on the
early-out), so both arms run on the *same bytes* in the *same process* and the only
difference is the defect.

Every number in ``evidence.json`` is computed here; nothing is transcribed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path

import tensorcad.iges as I
from tensorcad.iges import read_iges

#: bytes of the official archives, cached outside the repository (publishers grant no
#: redistribution, so the files are not committed).  sha256 is recorded in the evidence
#: so a reviewer can tell whether the same bytes were measured.
CANDIDATES = {
    "KCS": [
        "/root/sim-assets/m21253-379b12515423/artifacts/baseline_fleet/raw/closure/kcs/KCS_hull_SVA.igs",
        "/root/sim-assets/m20045-285e6d47d6f1/.cache/baseline_fleet/iges/KCS_hull_SVA.igs",
    ],
    "KVLCC2": [
        "/root/sim-assets/m20163-fedab06121c6/artifacts/baseline_fleet/raw/KVLCC2/kvlcc2.igs",
        "/root/sim-assets/m20045-285e6d47d6f1/.cache/baseline_fleet/iges/kvlcc2.igs",
    ],
    "JBC": [
        "/root/sim-assets/m21253-379b12515423/artifacts/baseline_fleet/raw/closure/jbc/jbc.igs",
        "/root/sim-assets/m20045-285e6d47d6f1/.cache/baseline_fleet/iges/jbc.igs",
    ],
    "ONRT": [
        "/root/sim-assets/m21253-379b12515423/artifacts/baseline_fleet/raw/closure/onrt/onr5613-full_appendage_deckhous_mergedR.igs",
        "/root/sim-assets/m20045-285e6d47d6f1/.cache/baseline_fleet/iges/onrt_wo_gap/onr5613-full_appendage_deckhous_mergedR.igs",
    ],
}

def _fixed_find_span(n, degree, knots, value):  # noqa: ANN001 - post-fix code, restated here
    lo, hi = degree, n + 1
    if value >= knots[n + 1]:
        return n
    if value <= knots[degree]:
        return degree
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if value < knots[mid]:
            hi = mid
        else:
            lo = mid
    return lo


def _legacy_find_span(n, degree, knots, value):  # noqa: ANN001 - byte-for-byte pre-fix code
    lo, hi = degree, n
    if value >= knots[hi]:
        return hi - 1
    if value <= knots[lo]:
        return lo
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if value < knots[mid]:
            hi = mid
        else:
            lo = mid
    return lo


def install(legacy: bool) -> None:
    I._find_span = _legacy_find_span if legacy else _fixed_find_span
    assert I._find_span is _legacy_find_span if legacy else I._find_span is _fixed_find_span


def surface_health(surfaces, samples: int = 13) -> dict:
    """Per-surface evaluation health: degenerate (all-zero) and out-of-hull samples."""
    zero = 0
    outside = 0
    worst = 0.0
    span = 0.0
    for s in surfaces:
        (x0, y0, z0), (x1, y1, z1) = s.bounding_box()
        u0, u1 = s.domain_u
        v0, v1 = s.domain_v
        patch_zero = 0
        patch_out = 0
        for a in range(samples):
            for b in range(samples):
                p = s.evaluate(u0 + (u1 - u0) * a / (samples - 1), v0 + (v1 - v0) * b / (samples - 1))
                if abs(p[0]) < 1e-12 and abs(p[1]) < 1e-12 and abs(p[2]) < 1e-12:
                    patch_zero += 1
                ok = (x0 - 1e-6 <= p[0] <= x1 + 1e-6 and y0 - 1e-6 <= p[1] <= y1 + 1e-6
                      and z0 - 1e-6 <= p[2] <= z1 + 1e-6)
                if not ok:
                    patch_out += 1
        zero += patch_zero
        outside += patch_out
        if patch_zero or patch_out:
            span += 1
        worst = max(worst, patch_zero + patch_out)
    return {
        "surfaces": len(surfaces),
        "samples_per_surface": samples * samples,
        "degenerate_samples": zero,
        "outside_control_net_samples": outside,
        "surfaces_with_defect": int(span),
        "worst_surface_defect_samples": int(worst),
    }


# --------------------------------------------------------------------------- #
# displacement volume of a half hull by a height field (no topology needed)
# --------------------------------------------------------------------------- #
def tris_from_patches(patches, grid: int) -> list:
    tris = []
    for basis, loops, point_in_polygon in patches:
        u0, u1 = basis.domain_u
        v0, v1 = basis.domain_v
        pts = [[basis.evaluate(u0 + (u1 - u0) * i / grid, v0 + (v1 - v0) * j / grid)
                for i in range(grid + 1)] for j in range(grid + 1)]
        for j in range(grid):
            for i in range(grid):
                if loops is not None:
                    uc = u0 + (u1 - u0) * (i + 0.5) / grid
                    vc = v0 + (v1 - v0) * (j + 0.5) / grid
                    if sum(1 for lp in loops if point_in_polygon(uc, vc, lp)) % 2 == 0:
                        continue
                a, b, c, d = pts[j][i], pts[j][i + 1], pts[j + 1][i + 1], pts[j + 1][i]
                tris.append((a, b, c))
                tris.append((a, c, d))
    return tris


def height_field(tris, cell: float):
    if not tris:
        return {}, cell, cell
    xs = [p[0] for t in tris for p in t]
    ys = [p[1] for t in tris for p in t]
    xlo, xhi = min(xs) - cell, max(xs) + cell
    ylo, yhi = min(ys) - cell, max(ys) + cell
    nx = max(1, int(math.ceil((xhi - xlo) / cell)))
    ny = max(1, int(math.ceil((yhi - ylo) / cell)))
    dx = (xhi - xlo) / nx
    dy = (yhi - ylo) / ny
    field: dict[tuple[int, int], float] = {}
    for tri in tris:
        (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) = tri
        det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
        if abs(det) < 1e-14:
            continue
        i0 = max(0, int((min(x1, x2, x3) - xlo) / dx))
        i1 = min(nx - 1, int((max(x1, x2, x3) - xlo) / dx))
        j0 = max(0, int((min(y1, y2, y3) - ylo) / dy))
        j1 = min(ny - 1, int((max(y1, y2, y3) - ylo) / dy))
        for i in range(i0, i1 + 1):
            cx = xlo + (i + 0.5) * dx
            for j in range(j0, j1 + 1):
                cy = ylo + (j + 0.5) * dy
                l1 = ((y2 - y3) * (cx - x3) + (x3 - x2) * (cy - y3)) / det
                l2 = ((y3 - y1) * (cx - x3) + (x1 - x3) * (cy - y3)) / det
                l3 = 1.0 - l1 - l2
                if l1 < -1e-9 or l2 < -1e-9 or l3 < -1e-9:
                    continue
                z = l1 * z1 + l2 * z2 + l3 * z3
                key = (i, j)
                if key not in field or z < field[key]:
                    field[key] = z
    return field, dx, dy


def volume_below(field, dx, dy, waterline: float, half: bool) -> tuple[float, int]:
    depth = [max(waterline - z, 0.0) for z in field.values()]
    wet = sum(1 for d in depth if d > 0.0)
    total = sum(depth) * dx * dy
    return (2.0 * total if half else total), wet


def kcs_volume(path: str, cell: float, grid: int, waterline: float = 10.8) -> dict:
    f = read_iges(path)
    surfaces = f.surfaces()
    patches = [(s, None, None) for s in surfaces]
    tris = tris_from_patches(patches, grid)
    field, dx, dy = height_field(tris, cell)
    vol, wet = volume_below(field, dx, dy, waterline, half=True)
    return {
        "surfaces": len(surfaces),
        "triangles": len(tris),
        "grid": grid,
        "cell_m": cell,
        "waterline_m": waterline,
        "wet_columns": wet,
        "displacement_volume_m3": round(vol, 1),
        "public_volume_m3": 52030.0,
        "error_pct": round(100.0 * (vol - 52030.0) / 52030.0, 3),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="runs/m21680_span_defect/evidence.json")
    args = ap.parse_args()

    evidence: dict = {
        "schema": "tensorcad.iges.span_defect.v1",
        "anchor": "GAP-CAD-T1-BASELINE-FLEET",
        "defect": {
            "site": "tensorcad/iges.py::_find_span",
            "legacy": "lo, hi = degree, n; if value >= knots[hi]: return hi - 1",
            "fixed": "lo, hi = degree, n + 1; if value >= knots[n + 1]: return n",
            "consequence": (
                "the early-out fires from the second-to-last support knot onwards, so the "
                "returned span is n-1 instead of n; the basis-function index "
                "span-degree+i then starts at -1, Python wraps it to the END of the knot "
                "vector and every basis value comes back 0.0 -- evaluate() divides 0/0 and "
                "returns (0,0,0) silently. A surface whose u or v direction has a single "
                "span (a Bezier patch, count == degree+1) hits this over its whole domain."
            ),
        },
        "ships": {},
        "kcs_volume": {},
    }

    for ship, paths in CANDIDATES.items():
        path = next((p for p in paths if os.path.exists(p)), None)
        if path is None:
            evidence["ships"][ship] = {"status": "bytes_not_cached_locally"}
            continue
        raw = Path(path).read_bytes()
        entry: dict = {
            "path": path,
            "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest(),
        }
        for arm in ("legacy", "fixed"):
            install(legacy=(arm == "legacy"))
            f = read_iges(path)
            health = surface_health(f.surfaces())
            entry[arm] = health
            if ship == "KVLCC2":
                from tensorcad.fleet.trimfield import point_in_polygon, trim_polygons

                records = trim_polygons(f)
                entry[arm]["trim_patches_decided"] = sum(1 for r in records if r["verdict"] == "decided")
                entry[arm]["trim_patches_total"] = len(records)
                patches = [
                    (r["basis"], r["loops"], point_in_polygon)
                    for r in records
                    if r.get("basis") is not None and r.get("loops")
                ]
            else:
                patches = [(s, None, None) for s in f.surfaces()]
            tris = tris_from_patches(patches, 40)
            field, dx, dy = height_field(tris, 2.0)
            zs = [p[2] for t in tris for p in t] or [0.0]
            entry[arm]["tessellated_triangles"] = len(tris)
            entry[arm]["height_field_columns"] = len(field)
            entry[arm]["z_range_m"] = [round(min(zs), 3), round(max(zs), 3)]
        install(legacy=False)
        if ship == "KCS":
            evidence["kcs_volume"]["fixed_cell_0.5_grid_24"] = kcs_volume(path, 0.5, 24)
            evidence["kcs_volume"]["fixed_cell_0.25_grid_32"] = kcs_volume(path, 0.25, 32)
        evidence["ships"][ship] = entry

    install(legacy=False)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(evidence, indent=1))
    print(json.dumps({k: v for k, v in evidence.items() if k != "ships"}, indent=1))
    for ship, entry in evidence["ships"].items():
        if "legacy" not in entry:
            print(ship, entry)
            continue
        print(f"{ship}: legacy {entry['legacy']}\n      fixed  {entry['fixed']}")
    print("written", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())