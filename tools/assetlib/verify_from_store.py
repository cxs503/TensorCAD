#!/usr/bin/env python3
"""Prove the durable store's bytes are *usable*, not merely present.

    python tools/assetlib/verify_from_store.py --out artifacts/assetlib/physics_from_store.json
    python tools/assetlib/verify_from_store.py --ships KCS --dx 0.5 0.25 0.125

For every hull in ``fleet_specs.json`` this script

1. resolves the spec pointer through :mod:`tensorcad.assetlib.consumer` (``asset://``
   URI -> store path) and records the SHA-256 it actually read, so the run is tied to
   bytes rather than to a path;
2. rebuilds the hull with :func:`tensorcad.fleet_baseline.build_ship` at each requested
   grid resolution (the resolution study doubles as the convergence check);
3. compares displacement volume, ``Cb`` and ``LCB`` against the published particulars
   the spec carries, with the relative error and the acceptance verdict.

The point of step 1 is that a registry with no consumer is a file nobody reads: a
number produced here can only exist if the store resolved and hashed.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from tensorcad.assetlib.consumer import resolve_pointer  # noqa: E402
from tensorcad.assetlib.registry import load_registry  # noqa: E402
from tensorcad.assetlib.store import sha256_file  # noqa: E402
from tensorcad.fleet_baseline import build_ship, specs_from_json  # noqa: E402

DEFAULT_SPECS = REPO / "benchmarks" / "verified" / "baseline_fleet" / "fleet_specs.json"
DEFAULT_REGISTRY = REPO / "assets" / "geometry_registry.json"
#: acceptance tolerance on the displacement volume, from the verified benchmark case
TOLERANCE_PERCENT = 1.0


def rel_error(ref: float | None, got: float | None) -> float | None:
    if ref in (None, 0) or got is None:
        return None
    return (float(got) - float(ref)) / float(ref) * 100.0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--specs", default=str(DEFAULT_SPECS))
    parser.add_argument("--registry", default=str(DEFAULT_REGISTRY))
    parser.add_argument("--ships", nargs="*", default=None)
    parser.add_argument("--dx", nargs="+", type=float, default=[0.5, 0.25])
    parser.add_argument("--out", default=str(REPO / "artifacts" / "assetlib" / "physics_from_store.json"))
    args = parser.parse_args(argv)

    registry = load_registry(args.registry)
    specs = specs_from_json(args.specs, asset_registry=registry)
    ships = [s for s in specs if args.ships is None or s.name in set(args.ships)]
    report: dict[str, object] = {
        "schema": "tensorcad.assetlib.physics_from_store.v1",
        "anchor": "GAP-GEOM-ASSET-LIB",
        "registry": str(args.registry),
        "specs": str(args.specs),
        "tolerance_percent": TOLERANCE_PERCENT,
        "resolutions_m": args.dx,
        "hulls": {},
    }
    failures = 0
    for spec in ships:
        entry: dict[str, object] = {
            "spec_pointer": spec.iges_path,
            "iges_sha256_expected": spec.iges_sha256,
            "published_volume_m3": spec.published_volume_m3,
            "published_cb": spec.published_cb,
            "published_lcb_percent": spec.published_lcb_percent,
            "resolutions": {},
        }
        binding = registry.bindings.get(spec.name)
        entry["asset_id"] = binding.asset_id if binding else ""
        try:
            resolved = resolve_pointer(spec.iges_path, spec.iges_sha256, registry=registry)
            entry["resolved_path"] = resolved
            entry["iges_sha256_read"] = sha256_file(resolved)
            entry["sha_match"] = entry["iges_sha256_read"] == spec.iges_sha256
        except Exception as exc:  # noqa: BLE001 - recorded, not swallowed
            entry["resolve_error"] = f"{type(exc).__name__}: {exc}"
            report["hulls"][spec.name] = entry
            failures += 1
            continue
        for dx in args.dx:
            started = time.time()
            try:
                _mesh, info = build_ship(spec, dx=dx)
                derived = info.get("derived", {}) or {}
                mesh_info = info.get("mesh", {}) or {}
                checks = info.get("checks", {}) or {}
                row = {
                    "verdict": info.get("verdict"),
                    "volume_m3": derived.get("volume_m3"),
                    "field_volume_m3": derived.get("field_volume_m3"),
                    "mesh_vs_field_volume_percent": derived.get("mesh_vs_field_volume_percent"),
                    "cb": derived.get("cb"),
                    "lcb_percent": derived.get("lcb_percent_of_lpp"),
                    "length_waterline_m": derived.get("length_waterline_m"),
                    "watertight": mesh_info.get("is_watertight"),
                    "edge_manifold": mesh_info.get("is_edge_manifold"),
                    "triangles": mesh_info.get("triangle_count"),
                    "checks": checks,
                    "reasons": info.get("reasons", []),
                    "seconds": round(time.time() - started, 2),
                }
                row["volume_rel_error_percent"] = rel_error(spec.published_volume_m3, row["volume_m3"])
                row["cb_rel_error_percent"] = rel_error(spec.published_cb, row["cb"])
                row["lcb_abs_error_percent_lpp"] = (
                    None
                    if spec.published_lcb_percent is None or row["lcb_percent"] is None
                    else row["lcb_percent"] - spec.published_lcb_percent
                )
                row["within_tolerance"] = (
                    None
                    if row["volume_rel_error_percent"] is None
                    else abs(row["volume_rel_error_percent"]) <= TOLERANCE_PERCENT
                )
                entry["resolutions"][f"dx={dx}"] = row
            except Exception as exc:  # noqa: BLE001
                entry["resolutions"][f"dx={dx}"] = {"error": f"{type(exc).__name__}: {exc}"}
        report["hulls"][spec.name] = entry
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps(report, indent=1))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())