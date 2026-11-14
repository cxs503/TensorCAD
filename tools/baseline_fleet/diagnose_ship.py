"""Per-ship diagnostic for the baseline fleet chain: where does a hull stop?

The fleet producer (:mod:`tensorcad.fleet_baseline`) can refuse a hull at several
distinct stages -- units, scale, frame, sampling, field coverage, slice count.
The benchmark only records the *terminal* message, so a refusal hides the stage
that actually broke.  This script walks the same stages as
:func:`tensorcad.fleet_baseline.build_ship` and prints each intermediate fact,
so a WITHHELD verdict can be attributed to one cause instead of guessed.

    python tools/baseline_fleet/diagnose_ship.py KVLCC2
    python tools/baseline_fleet/diagnose_ship.py --all --specs <fleet_specs.json>

Read-only: nothing is written, nothing is uploaded.  Standard library only plus
the package itself.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

from tensorcad.fleet_baseline import (
    DEFAULT_DX,
    DEFAULT_DZ,
    DEFAULT_SAMPLES_PER_METRE,
    field_from_points,
    partition_surfaces,
    sample_patches,
    specs_from_json,
    slice_profile,
)
from tensorcad.iges import read_iges, unit_scale


def _span(values):
    return (min(values), max(values)) if values else (float("nan"), float("nan"))


def diagnose(spec, samples_per_metre: float) -> dict:
    out: dict = {"ship": spec.name, "iges_path": spec.iges_path}
    model = read_iges(spec.iges_path)
    surfaces = model.surfaces()
    factor = unit_scale(model.unit_flag) * model.model_scale
    out["unit_flag"] = model.unit_flag
    out["unit_name"] = model.unit_name
    out["model_scale"] = model.model_scale
    out["factor_to_metres"] = factor
    out["patches_total"] = len(surfaces)
    usable, quarantined = partition_surfaces(surfaces)
    out["patches_usable"] = len(usable)
    out["patches_quarantined"] = len(quarantined)
    if not usable:
        out["stage_failed"] = "partition_surfaces"
        return out

    raw = sample_patches(usable, spacing=1.0 / max(samples_per_metre, 1e-6))
    out["samples"] = len(raw)
    xs = [p[0] * factor for p in raw]
    zs = [p[2] * factor for p in raw]
    ys = [p[1] * factor for p in raw]
    x0, x1 = _span(xs)
    z0, z1 = _span(zs)
    y0, y1 = _span(ys)
    out["bbox_m"] = {"x": [x0, x1], "y": [y0, y1], "z": [z0, z1]}
    out["extent_x_m"] = x1 - x0
    out["extent_x_over_lpp"] = (x1 - x0) / spec.length_pp if spec.length_pp else float("nan")

    # same frame guard as build_ship
    shift = 0.0
    if z0 < -0.5 * spec.draft and -0.05 * spec.draft <= z1 <= 0.05 * spec.draft:
        shift = -z1
    out["frame_shift_z"] = shift
    points = [(x, y, z + shift) for x, y, z in zip(xs, ys, zs)]
    z0 += shift
    z1 += shift
    out["z_after_shift"] = [z0, z1]
    out["draft"] = spec.draft
    below = sum(1 for _, _, z in points if 0.0 <= z <= spec.draft)
    out["samples_below_draft"] = below
    out["samples_below_draft_fraction"] = below / len(points) if points else 0.0

    if not below:
        out["stage_failed"] = "no samples in [0, draft] after frame guard"
        return out

    field = field_from_points(points, spec.draft, dx=DEFAULT_DX, dz=DEFAULT_DZ)
    levels = [
        k for k in range(field.nz) if any(field.b[i][k] > 0.0 for i in range(field.nx))
    ]
    out["field"] = field.to_dict()
    out["field_usable_levels"] = len(levels)
    out["field_usable_level_x_range"] = [
        [field.x(0), field.x(field.nx - 1)] if levels else None
    ][0]
    # slice_profile on each level with data: a slice is "usable" only if it has
    # at least two distinct x positions.
    usable_slices = 0
    slice_report = []
    for k in levels:
        profile = slice_profile(field, k)
        if profile is None:
            slice_report.append([round(field.z(k), 3), "no profile"])
            continue
        xs_k = profile[0]
        if len(xs_k) >= 2 and abs(xs_k[-1] - xs_k[0]) > 1e-9:
            usable_slices += 1
            slice_report.append([round(field.z(k), 3), round(xs_k[-1] - xs_k[0], 3)])
        else:
            slice_report.append([round(field.z(k), 3), "degenerate x-range"])
    out["usable_slices"] = usable_slices
    out["slice_report"] = slice_report[:8] + ([["...", len(slice_report)]] if len(slice_report) > 8 else [])
    out["stage_failed"] = None if usable_slices >= 2 else "build_solid (<2 usable slices)"
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ships", nargs="*", help="ship names (default: all in --specs)")
    parser.add_argument("--specs", default="benchmarks/verified/baseline_fleet/fleet_specs.json")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--samples-per-metre", type=float, default=DEFAULT_SAMPLES_PER_METRE)
    parser.add_argument("--out", default=None, help="write the JSON report here")
    args = parser.parse_args(argv)

    specs = specs_from_json(args.specs)
    wanted = None if (args.all or not args.ships) else set(args.ships)
    reports = []
    for spec in specs:
        if wanted is not None and spec.name not in wanted:
            continue
        try:
            reports.append(diagnose(spec, args.samples_per_metre))
        except Exception as exc:  # a refusal is evidence, not a crash
            reports.append({"ship": spec.name, "stage_failed": f"{type(exc).__name__}: {exc}"})
    payload = {"schema": "tensorcad.fleet.diagnose/1.0", "ships": reports}
    text = json.dumps(payload, indent=1)
    print(text)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(text)
    return 0 if all(r.get("stage_failed") is None for r in reports) else 1


if __name__ == "__main__":
    sys.exit(main())