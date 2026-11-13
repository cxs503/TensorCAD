"""Run the geometry contract validator on all 8 STL assets and produce the compliance table.

This script knows the *declared* units and coordinate frame for each asset
(because STL files do not carry this metadata) and passes them to the validator.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

# Add repo root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tensorcad.geom_contract import validate_stl_file, ContractReport

# Declared metadata for each asset (STL does not carry this)
ASSET_METADATA: dict[str, dict[str, str]] = {
    "wigley": {
        "units": "m",
        "coordinate_frame": "x=longitudinal, from the aft end (x=0) toward the forward end (x=length); y=transverse, positive to starboard; z=vertical, positive up; z=0 is the design waterline, the keel lies at z=-draft",
    },
    "series60_cb60": {
        "units": "m",
        "coordinate_frame": "x aft->fore, y port->starboard, z up; right-handed",
    },
    "series60_cb75": {
        "units": "m",
        "coordinate_frame": "x aft->fore, y port->starboard, z up; right-handed",
    },
    "box_barge": {
        "units": "m",
        "coordinate_frame": "x aft->fore, y port->starboard, z up; right-handed",
    },
    "sphere_r1": {
        "units": "m",
        "coordinate_frame": "x fore, y starboard, z up; right-handed",
    },
    # Red-team assets
    "inverted_wigley": {
        "units": "m",
        "coordinate_frame": "x=longitudinal, from the aft end (x=0) toward the forward end (x=length); y=transverse, positive to starboard; z=vertical, positive up; z=0 is the design waterline, the keel lies at z=-draft",
    },
    "wrong_units": {
        "units": "mm",  # deliberately wrong: geometry is in metres but declared as mm
        "coordinate_frame": "x=longitudinal, from the aft end (x=0) toward the forward end (x=length); y=transverse, positive to starboard; z=vertical, positive up; z=0 is the design waterline, the keel lies at z=-draft",
    },
    "left_handed": {
        "units": "m",
        "coordinate_frame": "x fore, y starboard, z down; left-handed",  # deliberately wrong
    },
}


def run_compliance(assets_dir: Path) -> list[ContractReport]:
    """Validate all 8 assets and return reports."""
    reports: list[ContractReport] = []
    for name in sorted(ASSET_METADATA.keys()):
        stl_path = assets_dir / f"{name}.stl"
        if not stl_path.exists():
            print(f"  SKIP: {stl_path} not found")
            continue
        meta = ASSET_METADATA[name]
        report = validate_stl_file(
            stl_path,
            units=meta["units"],
            coordinate_frame=meta["coordinate_frame"],
        )
        reports.append(report)
    return reports


def print_table(reports: list[ContractReport]) -> None:
    """Print a human-readable compliance table."""
    print("=" * 90)
    print(f"{'Asset':<20} {'Status':<8} {'Units':<6} {'Diagonal (m)':<14} {'Volume':<14} {'Watertight':<10} {'Violations'}")
    print("=" * 90)
    for r in reports:
        status = "PASS" if r.passed else "FAIL"
        info = r.mesh_info
        diag = info.get("diagonal_m", "N/A")
        vol = info.get("signed_volume", "N/A")
        wt = "yes" if info.get("watertight") else "no"
        units = info.get("units", "?")
        viol = "; ".join(f"[{v.convention}] {v.message}" for v in r.violations) if r.violations else "—"
        print(f"{r.source.split('/')[-1]:<20} {status:<8} {units:<6} {str(diag):<14} {str(vol):<14} {wt:<10} {viol}")
    print("=" * 90)
    passed = sum(1 for r in reports if r.passed)
    print(f"\n{passed}/{len(reports)} assets passed the geometry contract.")


def main() -> None:
    assets_dir = Path(__file__).resolve().parent.parent / "artifacts" / "geom_contract_assets"
    reports = run_compliance(assets_dir)

    # Print human-readable table
    print_table(reports)

    # Save JSON report
    report_path = assets_dir / "compliance_report.json"
    payload = [r.to_dict() for r in reports]
    report_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"\nJSON report saved to {report_path}")


if __name__ == "__main__":
    main()
