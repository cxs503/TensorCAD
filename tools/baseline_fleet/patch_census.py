#!/usr/bin/env python3
"""Census of how each fleet hull's IGES file binds its trimmed patches.

    python tools/baseline_fleet/patch_census.py --out docs/fleet_patch_census.md \
        --json benchmarks/verified/fleet_refs/patch_census.json

For every hull in the fleet manifest this reports, per file:

* the entity histogram (which types the hull arrives in);
* how many type-144/143 wrappers resolve their **outer trim boundary** and how many
  of those curves resolve their basis surface -- the two bindings that a reader has
  to get right before trimming means anything (see ``tensorcad.iges``);
* the length of the sampled (u, v) trim polygon per patch;
* whether the patch's basis surface evaluates to a *degenerate net* (every sample at
  the model origin), which is a decode failure, not a flat piece of hull.

The point of the census is attribution.  A hull whose benchmark verdict reads
"field produced fewer than two usable waterline slices" has been recorded WITHHELD,
and that message names a downstream symptom.  The census names the upstream cause
per patch, so the next round of work starts from a count instead of a guess.

Everything here is read-only over the asset store; publisher bytes are not copied
into the repository.  Numbers reach the markdown only through this program.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from tensorcad.iges import (  # noqa: E402
    IgesFile,
    read_iges,
)

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "benchmarks" / "verified" / "baseline_fleet" / "fleet_specs.json"


def _inside(polygon: list[tuple[float, float]], u: float, v: float) -> bool:
    """Even-odd point-in-polygon test on the trim loop."""
    inside = False
    count = len(polygon)
    for index in range(count):
        u0, v0 = polygon[index]
        u1, v1 = polygon[(index + 1) % count]
        if (v0 > v) != (v1 > v):
            crossing = u0 + (v - v0) * (u1 - u0) / (v1 - v0)
            if u < crossing:
                inside = not inside
    return inside


def patch_census(path: Path, grid: int = 12) -> dict[str, Any]:
    """Per-patch binding + degeneracy census for one IGES file."""
    file: IgesFile = read_iges(str(path))
    histogram = file.stats()
    rows: list[dict[str, Any]] = []
    for wrapper, basis, curves in file.trimmed_surfaces():
        row: dict[str, Any] = {
            "wrapper_de": wrapper.entity_type,
            "basis_surface_bound": basis is not None,
            "boundary_count": wrapper.boundary_count,
            "boundary_bound": wrapper.boundary is not None,
            "curve_on_surface": len(curves),
            "uv_loop_points": 0,
            "sample_points": 0,
            "degenerate_samples": 0,
            "bbox": None,
        }
        if curves and curves[0].surface:
            row["surface_pointer_resolved"] = True
        else:
            row["surface_pointer_resolved"] = False
        if curves and curves[0].curve_uv and basis is not None:
            try:
                polygon = file.trim_polygon_uv(curves[0].curve_uv, 200)
            except Exception as exc:  # noqa: BLE001 - the census reports, it does not raise
                row["trim_error"] = f"{type(exc).__name__}: {exc}"
                polygon = []
            row["uv_loop_points"] = len(polygon)
            if polygon and basis is not None:
                u0, u1 = basis.domain_u
                v0, v1 = basis.domain_v
                xs, ys, zs = [], [], []
                for i in range(grid):
                    u = u0 + (u1 - u0) * (i + 0.5) / grid
                    for j in range(grid):
                        v = v0 + (v1 - v0) * (j + 0.5) / grid
                        if not _inside(polygon, u, v):
                            continue
                        point = basis.evaluate(u, v)
                        row["sample_points"] += 1
                        if all(abs(q) < 1e-9 for q in point):
                            row["degenerate_samples"] += 1
                            continue
                        xs.append(point[0])
                        ys.append(point[1])
                        zs.append(point[2])
                if xs:
                    row["bbox"] = [
                        [round(min(xs), 2), round(min(ys), 2), round(min(zs), 2)],
                        [round(max(xs), 2), round(max(ys), 2), round(max(zs), 2)],
                    ]
        rows.append(row)
    usable = [r for r in rows if r["sample_points"] and r["degenerate_samples"] < r["sample_points"]]
    return {
        "file": str(path.name),
        "histogram": histogram,
        "patches": len(rows),
        "trim_bound": sum(1 for r in rows if r["boundary_bound"]),
        "surface_pointer_resolved": sum(1 for r in rows if r["surface_pointer_resolved"]),
        "uv_loops": sum(1 for r in rows if r["uv_loop_points"]),
        "degenerate_patches": sum(
            1 for r in rows if r["sample_points"] and r["degenerate_samples"] == r["sample_points"]
        ),
        "partly_degenerate_patches": sum(
            1 for r in rows if 0 < r["degenerate_samples"] < r["sample_points"]
        ),
        "usable_patches": len(usable),
        "per_patch": rows,
    }


def _iges_strings(node: Any) -> list[str]:
    """Every string in a manifest record that names an ``.igs`` file, at any depth."""
    found: list[str] = []
    if isinstance(node, str):
        if node.lower().endswith(".igs"):
            found.append(node)
    elif isinstance(node, dict):
        for value in node.values():
            found.extend(_iges_strings(value))
    elif isinstance(node, list):
        for value in node:
            found.extend(_iges_strings(value))
    return found


def _resolve(path: str) -> Path | None:
    """Resolve a manifest geometry path (``asset://...``) against the asset store."""
    if not str(path).strip():
        return None
    if path.startswith("asset://"):
        candidate = Path("/root/tp-agent/assets") / path[len("asset://"):]
        return candidate if candidate.exists() else None
    candidate = Path(path)
    return candidate if candidate.exists() else None


def census_manifest(manifest: Path = MANIFEST) -> dict[str, Any]:
    data = json.loads(manifest.read_text())
    hulls: dict[str, Any] = {}
    for record in data["ships"]:
        name = record["name"]
        candidates = _iges_strings(record)
        raw = candidates[0] if candidates else ""
        path = None
        for candidate in candidates or [raw]:
            path = _resolve(candidate)
            if path is not None:
                raw = candidate
                break
        if path is None:
            hulls[name] = {"file": raw, "status": "asset-absent"}
            continue
        entry = patch_census(path)
        entry["status"] = "read"
        entry["reference_volume_m3"] = record.get("published_volume_m3")
        hulls[name] = entry
    return {"schema": "tensorcad.fleet.patch_census.v1", "hulls": hulls}


def render(payload: dict[str, Any]) -> str:
    lines: list[str] = []
    add = lines.append
    add("# Fleet patch census — what each IGES file binds, per trimmed patch\n")
    add(
        "Generated by `tools/baseline_fleet/patch_census.py`; every number below is injected "
        "from the JSON it writes (`benchmarks/verified/fleet_refs/patch_census.json`), none is typed by hand.\n"
    )
    add("## Summary\n")
    add("| hull | file | #144 | trim bound | uv loops | surface ptr | usable patches | degenerate |")
    add("|---|---|---|---|---|---|---|---|")
    for name, entry in sorted(payload["hulls"].items()):
        if entry.get("status") != "read":
            add(f"| {name} | `{entry.get('file')}` | — | — | — | — | — | asset absent |")
            continue
        deg = f"{entry['degenerate_patches']} full / {entry['partly_degenerate_patches']} partial"
        add(
            f"| {name} | `{entry['file']}` | {entry['patches']} | {entry['trim_bound']} | "
            f"{entry['uv_loops']} | {entry['surface_pointer_resolved']} | {entry['usable_patches']} | {deg} |"
        )
    add("")
    add("## Entity histograms\n")
    for name, entry in sorted(payload["hulls"].items()):
        if entry.get("status") != "read":
            continue
        hist = ", ".join(f"type {key}: {value}" for key, value in entry["histogram"].items())
        add(f"* **{name}** — {hist}")
    add("")
    add("## Per-patch detail\n")
    for name, entry in sorted(payload["hulls"].items()):
        if entry.get("status") != "read":
            continue
        add(f"### {name} (`{entry['file']}`)\n")
        add("| # | basis bound | outer boundary | uv loop pts | samples | degenerate | bbox x | bbox y | bbox z |")
        add("|---|---|---|---|---|---|---|---|---|")
        for index, row in enumerate(entry["per_patch"]):
            bbox = row["bbox"]
            bx = f"{bbox[0][0]} … {bbox[1][0]}" if bbox else "—"
            by = f"{bbox[0][1]} … {bbox[1][1]}" if bbox else "—"
            bz = f"{bbox[0][2]} … {bbox[1][2]}" if bbox else "—"
            add(
                f"| {index} | {'yes' if row['basis_surface_bound'] else 'NO'} | "
                f"{'yes' if row['boundary_bound'] else 'NO'} | {row['uv_loop_points']} | "
                f"{row['sample_points']} | {row['degenerate_samples']} | {bx} | {by} | {bz} |"
            )
        add("")
    add("## Reading\n")
    add(
        "* **trim bound / #144** is the binding `tensorcad.iges` has to resolve before a trimmed patch "
        "means anything. Before the fix behind this census the count was `0` on every hull that ships "
        "type-144 wrappers, because the outer boundary was read out of parameter 4 and the "
        "curve-on-surface basis surface out of parameter 2 — both off by one field. A mis-bound patch "
        "is silently treated as its whole basis surface, which inflates `max|y|` in the half-breadth "
        "field and — on a file whose patches are mostly untrimmed overhang — moves the volume."
    )
    add(
        "* **degenerate** counts patches whose basis surface evaluates to the model origin at every "
        "sample inside its own trim loop. Those patches contribute nothing to a reconstruction, so a "
        "volume built from this file is a *partial* hull, and any verdict on it must say so."
    )
    add(
        "* The census is deliberately blind to the reference volume: it reports reader coverage, not "
        "accuracy. A hull can have every patch bound and still be `WITHHELD` for a frame reason."
    )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", default=str(MANIFEST))
    parser.add_argument("--json", default=str(ROOT / "benchmarks/verified/fleet_refs/patch_census.json"))
    parser.add_argument("--out", default=str(ROOT / "docs/fleet_patch_census.md"))
    args = parser.parse_args(argv)
    payload = census_manifest(Path(args.manifest))
    Path(args.json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.json).write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n")
    Path(args.out).write_text(render(payload))
    print(f"wrote {args.json} and {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())