#!/usr/bin/env python3
"""Growth metrics for the geometry asset library, plus the summary figure.

    python tools/assetlib/growth_report.py --out-dir artifacts/assetlib

Two families of metric, both machine-computed from artifacts already in git:

*demand side* (can downstream still work?)
    ``resolvable_geometry_pointer_fraction`` — dead geometry pointers in the fleet spec
    before/after repair; ``spec_ship_coverage`` — share of spec ships bound to a
    registered asset; ``store_managed_bytes`` — bytes reachable through the store.

*supply side* (is the library actually growing?)
    ``assets`` / ``bytes`` / ``by_hull`` from the registry, differentiated against the
    last line of ``assets/geometry_registry_history.jsonl`` when one exists.

The baseline is the state of the *committed* spec file before the repair
(``artifacts/assetlib/pointer_scan_before.json``, produced from ``git show
HEAD:...fleet_specs.json``), not a remembered number.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from tensorcad.assetlib.registry import load_registry  # noqa: E402


def _read(path: Path) -> dict:
    return json.loads(path.read_text()) if path.is_file() else {}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", default=str(REPO / "artifacts" / "assetlib"))
    parser.add_argument("--registry", default=str(REPO / "assets" / "geometry_registry.json"))
    parser.add_argument("--history", default=str(REPO / "assets" / "geometry_registry_history.jsonl"))
    parser.add_argument("--figure", default=None)
    args = parser.parse_args(argv)
    out = Path(args.out_dir)

    registry = load_registry(args.registry)
    current = registry.snapshot_dict()
    history = Path(args.history)
    snapshots = [json.loads(line) for line in history.read_text().splitlines() if line.strip()] if history.is_file() else []
    previous = snapshots[-2] if len(snapshots) > 1 else None
    growth = registry.growth(previous) if previous else None

    before = _read(out / "pointer_scan_before.json")
    after = _read(out / "pointer_scan_after.json")
    physics = _read(out / "physics_from_store.json")

    def geom_fraction(scan: dict) -> float:
        detail = scan.get("detail", [])
        geometry = [p for p in detail if p.get("kind") == "geometry"]
        if not geometry:
            return 0.0
        return round(sum(1 for p in geometry if p["status"] == "resolved") / len(geometry), 6)

    row_before = {
        "label": "HEAD spec file (dead-pointer era)",
        "store_managed_bytes": 0,
        "store_assets": 0,
        "resolvable_geometry_pointer_fraction": geom_fraction(before),
        "dead_pointers": before.get("dead", 0),
        "withheld_pointers": before.get("withheld", 0),
        "spec_ship_coverage": 0.0,
    }
    row_after = {
        "label": "after asset-library repair",
        "store_managed_bytes": current["bytes"],
        "store_assets": current["assets"],
        "resolvable_geometry_pointer_fraction": geom_fraction(after),
        "dead_pointers": after.get("dead", 0),
        "withheld_pointers": after.get("withheld", 0),
        "spec_ship_coverage": current["coverage"],
    }
    report = {
        "schema": "tensorcad.assetlib.growth_report.v1",
        "anchor": "GAP-GEOM-ASSET-LIB",
        "baseline": row_before,
        "current": row_after,
        "delta": {
            key: (round(row_after[key] - row_before[key], 6)
                  if isinstance(row_before[key], (int, float)) and not isinstance(row_before[key], bool)
                  else None)
            for key in row_after
            if key != "label"
        },
        "registry_snapshot": {k: current[k] for k in ("assets", "bytes", "verified", "verified_fraction",
                                                      "spec_bindings", "spec_ships", "coverage", "by_hull")},
        "history_lines": len(snapshots),
        "snapshot_growth": None if growth is None else growth["delta"],
        "physics": {
            ship: {
                "sha_match": entry.get("sha_match"),
                "asset_id": entry.get("asset_id"),
                "resolutions": {
                    key: {
                        "volume_rel_error_percent": row.get("volume_rel_error_percent"),
                        "verdict": row.get("verdict"),
                    }
                    for key, row in entry.get("resolutions", {}).items()
                },
            }
            for ship, entry in physics.get("hulls", {}).items()
        },
    }
    out.mkdir(parents=True, exist_ok=True)
    target = out / "growth_report.json"
    target.write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps(report, indent=1))

    figure = Path(args.figure) if args.figure else out / "assetlib_growth.png"
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - figure is evidence, not the product
        print(f"figure skipped: {exc}")
        return 0
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    axes[0].bar(["HEAD spec\n(dead pointers)", "after repair"],
                [row_before["resolvable_geometry_pointer_fraction"],
                 row_after["resolvable_geometry_pointer_fraction"]],
                color=["#b03a2e", "#1e8449"])
    axes[0].set_ylim(0, 1.05)
    axes[0].set_ylabel("resolvable geometry pointers / all geometry pointers")
    axes[0].set_title("(a) dead-pointer removal in fleet_specs.json")
    for index, value in enumerate([row_before["resolvable_geometry_pointer_fraction"],
                                   row_after["resolvable_geometry_pointer_fraction"]]):
        axes[0].text(index, value + 0.03, f"{value:.2f}", ha="center")

    hulls = sorted(current["by_hull"])
    axes[1].bar(hulls, [current["by_hull"][h]["bytes"] / 1e6 for h in hulls], color="#2874a6")
    axes[1].set_ylabel("store-managed bytes [MB]")
    axes[1].set_title(f"(b) library: {current['assets']} assets, {current['bytes']/1e6:.1f} MB, "
                      f"coverage {current['coverage']:.0%}")

    kcs = physics.get("hulls", {}).get("KCS", {}).get("resolutions", {})
    dxs = sorted((float(k.split("=")[1]), v) for k, v in kcs.items() if v.get("volume_rel_error_percent") is not None)
    if dxs:
        axes[2].plot([d for d, _ in dxs], [v["volume_rel_error_percent"] for _, v in dxs], "o-", label="got vs published")
        axes[2].axhline(0.0, color="#1e8449", label="published 52030 m3 (ref)")
        axes[2].axhline(1.0, color="#b03a2e", ls="--", label="+1% acceptance")
        axes[2].set_xlabel("grid spacing dx [m] (finer to the left)")
        axes[2].set_ylabel("displacement volume error [%]")
        axes[2].set_title("(c) KCS rebuilt from store bytes")
        axes[2].invert_xaxis()
        axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(figure, dpi=130)
    print(f"figure: {figure}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())