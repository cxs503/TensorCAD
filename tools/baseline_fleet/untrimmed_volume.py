#!/usr/bin/env python3
"""Divergence-theorem volume of the type-128 basis surfaces, *without* their trims.

This is the measurement behind "the patch set does not close the hull": it answers
what a consumer that ignores 141/142/143/144 integrates.  It is deliberately not a
verdict -- an untrimmed patch set is not a hull.

    python tools/baseline_fleet/untrimmed_volume.py FILE.igs --ref 312622 --grid 40
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from tensorcad.iges import read_iges, unit_scale  # noqa: E402


def integrate(path: str, grid: int = 40, scale: float | None = None) -> dict:
    model = read_iges(path)
    factor = unit_scale(model.unit_flag) * model.model_scale if scale is None else scale
    total = 0.0
    per_patch = []
    for index, surface in enumerate(model.surfaces()):
        u0, u1 = surface.domain_u
        v0, v1 = surface.domain_v
        volume = 0.0
        for i in range(grid):
            for j in range(grid):
                a = surface.evaluate(u0 + (u1 - u0) * i / grid, v0 + (v1 - v0) * j / grid)
                b = surface.evaluate(u0 + (u1 - u0) * (i + 1) / grid, v0 + (v1 - v0) * j / grid)
                c = surface.evaluate(u0 + (u1 - u0) * (i + 1) / grid, v0 + (v1 - v0) * (j + 1) / grid)
                d = surface.evaluate(u0 + (u1 - u0) * i / grid, v0 + (v1 - v0) * (j + 1) / grid)
                for p, q, r in ((a, b, c), (a, c, d)):
                    volume += (
                        p[0] * (q[1] * r[2] - q[2] * r[1])
                        - p[1] * (q[0] * r[2] - q[2] * r[0])
                        + p[2] * (q[0] * r[1] - q[1] * r[0])
                    ) / 6.0
        volume *= factor ** 3
        total += volume
        per_patch.append(round(volume, 2))
    return {
        "file": Path(path).name,
        "grid": grid,
        "unit_factor_to_metres": factor,
        "surfaces": len(per_patch),
        "signed_volume_m3": round(total, 2),
        "magnitude_m3": round(abs(total), 2),
        "per_patch_m3": per_patch,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file")
    parser.add_argument("--grid", type=int, default=40)
    parser.add_argument("--ref", type=float, default=None, help="published displacement (m3)")
    parser.add_argument("--scale", type=float, default=None, help="override unit factor")
    parser.add_argument("--out", default=None)
    args = parser.parse_args(argv)
    row = integrate(args.file, args.grid, args.scale)
    if args.ref:
        row["reference_m3"] = args.ref
        row["rel_error_percent"] = round((row["magnitude_m3"] - args.ref) / args.ref * 100.0, 3)
    print(json.dumps({k: v for k, v in row.items() if k != "per_patch_m3"}, indent=2))
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(row, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
