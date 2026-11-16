#!/usr/bin/env python3
"""Reachability report for the geometry asset library (round-3 hull-primary work).

Answers two questions with data instead of prose:

1. *Which hulls can a consumer actually ask for by name?* — before this round only the
   four ships a fleet spec binds could be resolved; the store held six hulls, so two were
   unreachable while ``coverage`` still read 1.0.
2. *Is the byte the resolver returns the byte that was measured?* — the DTC hull's
   SHA-256 is compared with the fingerprint recorded by
   ``benchmarks/verified/fleet_refs/dtc_hull_measured.json``, and that record's error is
   plotted against the platform's 3 % line.

    python tools/assetlib/resolution_report.py --out artifacts/assetlib_hull_primary.png

Everything printed is computed from ``assets/geometry_registry.json``, the durable store
and the measurement record; nothing is typed in by hand.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from tensorcad.assetlib.registry import RegistryError, load_registry  # noqa: E402
from tensorcad.assetlib.store import AssetStore  # noqa: E402

REGISTRY = REPO / "assets" / "geometry_registry.json"
DTC_RECORD = REPO / "benchmarks" / "verified" / "fleet_refs" / "dtc_hull_measured.json"


def before_after(registry, store) -> list[dict]:
    """Per hull: reachable by the old rule (spec binding only) and by the new one."""
    rows = []
    for hull in sorted(registry.hulls):
        bound = registry.bindings.get(hull)
        old = bool(bound and bound.asset_id)
        try:
            entry, rule = registry.hull_primary(hull)
            new, asset_id = True, entry.asset_id
        except RegistryError as exc:
            new, asset_id, rule = False, "", str(exc)[:60]
        path = ""
        if new:
            path = store.verify_asset(asset_id, registry.entry(asset_id).sha256)
        rows.append(
            {
                "hull": hull,
                "old_reachable": old,
                "new_reachable": new,
                "rule": rule,
                "asset_id": asset_id,
                "path": path,
                "manifest_status": str(registry.hulls.get(hull, {}).get("status", "")),
            }
        )
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", default=str(REGISTRY))
    parser.add_argument("--out", default=str(REPO / "artifacts" / "assetlib_hull_primary.png"))
    parser.add_argument("--json", default="")
    args = parser.parse_args(argv)

    registry = load_registry(args.registry)
    store = AssetStore(registry.store_root or None)
    rows = before_after(registry, store)
    snapshot = registry.snapshot_dict()

    dtc = json.loads(DTC_RECORD.read_text())
    dtc_sha = dtc["source"]["hull_sha256"]
    dtc_row = next((row for row in rows if row["hull"] == "DTC"), None)
    dtc_sha_match = bool(dtc_row and dtc_row["asset_id"] == "geometry/hulls/DTC/DTC_hull.igs")
    if dtc_sha_match:
        dtc_sha_match = registry.entry(dtc_row["asset_id"]).sha256 == dtc_sha

    summary = {
        "hulls": len(rows),
        "old_reachable": sum(row["old_reachable"] for row in rows),
        "new_reachable": sum(row["new_reachable"] for row in rows),
        "coverage_spec_ships": snapshot["coverage"],
        # the new denominator applied to the *old* rules: what hull_coverage would have
        # reported before this round, so the delta is not a change of definition only
        "hull_coverage_before": round(sum(row["old_reachable"] for row in rows) / len(rows), 6),
        "hull_coverage": snapshot["hull_coverage"],
        "rows": rows,
        "dtc": {
            "measured_sha256": dtc_sha,
            "resolved_sha256": registry.entry(dtc_row["asset_id"]).sha256 if dtc_row else "",
            "sha_match": dtc_sha_match,
            "checks": dtc["measured"]["checks"],
        },
    }
    print(json.dumps(summary, indent=1))
    if args.json:
        Path(args.json).write_text(json.dumps(summary, indent=1) + "\n")

    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:  # pragma: no cover - optional dependency
        print(f"matplotlib unavailable ({exc}); JSON only", file=sys.stderr)
        return 0

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(12.4, 4.6))

    hulls = [row["hull"] for row in rows]
    old = [1 if row["old_reachable"] else 0 for row in rows]
    new = [1 if row["new_reachable"] else 0 for row in rows]
    ypos = range(len(hulls))
    ax0.barh([y + 0.2 for y in ypos], old, height=0.38, color="#9aa5b1", label="round 1/2: spec binding only")
    ax0.barh([y - 0.2 for y in ypos], new, height=0.38, color="#2f7ed8", label="round 3: + hull-primary rules")
    for y, row in zip(ypos, rows):
        if not row["old_reachable"] and row["new_reachable"]:
            ax0.text(0.02, y - 0.2, f"  now reachable via {row['rule']}", va="center", fontsize=8, color="#12467b")
    ax0.set_yticks(list(ypos))
    ax0.set_yticklabels(hulls)
    ax0.set_xlim(0, 1.6)
    ax0.set_xlabel("reachable by name (1 = yes)")
    ax0.set_title(f"hull reachability: {sum(old)}/{len(rows)} -> {sum(new)}/{len(rows)}")
    ax0.legend(loc="lower right", fontsize=8)
    ax0.grid(axis="x", alpha=0.3)

    checks = dtc["measured"]["checks"]
    labels = {"displacement_volume": "displacement nabla", "block_coefficient": "Cb", "beam": "beam B"}
    names, errors = [], []
    for key, label in labels.items():
        names.append(label)
        errors.append(checks[key]["rel_error_percent"])
    colors = ["#2f7ed8" if abs(err) < 1.0 else "#d9534f" for err in errors]
    ax1.bar(names, errors, color=colors, width=0.5)
    ax1.axhline(3.0, color="#d9534f", ls="--", lw=1.2, label="platform line 3 %")
    ax1.axhline(-3.0, color="#d9534f", ls="--", lw=1.2)
    ax1.axhline(1.0, color="#888", ls=":", lw=1.0, label="acceptance 1 %")
    ax1.axhline(-1.0, color="#888", ls=":", lw=1.0)
    for i, err in enumerate(errors):
        ax1.annotate(f"{err:+.4f} %", (i, err), textcoords="offset points", xytext=(0, 6 if err >= 0 else -14), ha="center", fontsize=8)
    ax1.set_ylim(-2.0, 2.0)
    ax1.set_ylabel("rel. error vs published [%]")
    ax1.set_title("DTC: the byte now resolvable is the measured byte\n(registry sha == measurement record sha)")
    ax1.legend(fontsize=8, loc="upper right")
    ax1.grid(axis="y", alpha=0.3)

    fig.suptitle(
        f"geometry asset library: hull coverage {sum(old) / len(rows):.2f} -> {snapshot['hull_coverage']:.2f} "
        f"({snapshot['hulls_resolvable']}/{snapshot['hulls']} hulls), spec coverage {snapshot['coverage']:.2f}"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"figure -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())