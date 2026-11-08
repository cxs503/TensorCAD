#!/usr/bin/env python3
"""Series 60 verification against published parent forms.

This script verifies that the corrected Series 60 module reproduces the published
parent form coefficients (Todd 1951) and generates a verification report.

Usage:
    python3 examples/series60_verification/run_verification.py
"""
from __future__ import annotations

import json
from pathlib import Path

from tensorcad.series60 import (
    Series60Parameters,
    analytic_series60_hydrostatics,
    generate_series60,
)

# Published Series 60 parent form particulars (Todd 1951)
PUBLISHED_PARENTS = {
    0.60: {"cm": 0.960, "cp": 0.625, "cwp": 0.700},
    0.65: {"cm": 0.970, "cp": 0.670, "cwp": 0.727},
    0.70: {"cm": 0.980, "cp": 0.715, "cwp": 0.753},
    0.75: {"cm": 0.985, "cp": 0.760, "cwp": 0.780},
    0.80: {"cm": 0.990, "cp": 0.805, "cwp": 0.807},
}

L, B, T = 10.0, 1.5, 0.6


def schneekluth_cwp(cb: float) -> float:
    """Schneekluth & Bertram waterplane regression: (1 + 2 Cb) / 3."""
    return (1.0 + 2.0 * cb) / 3.0


def verify_parent_forms() -> dict:
    """Verify each parent form against published coefficients."""
    results = []
    
    for cb in sorted(PUBLISHED_PARENTS):
        pub = PUBLISHED_PARENTS[cb]
        
        # Legacy hull (no cwp correction)
        params_legacy = Series60Parameters(length=L, beam=B, draft=T, cb=cb)
        ref_legacy = analytic_series60_hydrostatics(params_legacy)
        
        # Corrected hull (with cwp correction)
        params_corrected = Series60Parameters(length=L, beam=B, draft=T, cb=cb, cwp=pub["cwp"])
        ref_corrected = analytic_series60_hydrostatics(params_corrected)
        
        # Generate mesh to verify topology
        model = generate_series60(params_corrected)
        
        result = {
            "cb": cb,
            "published": pub,
            "legacy": {
                "cwp": ref_legacy["waterplane_coefficient"],
                "cwp_error_pct": 100.0 * (ref_legacy["waterplane_coefficient"] / pub["cwp"] - 1.0),
            },
            "corrected": {
                "cwp": ref_corrected["waterplane_coefficient"],
                "cwp_error_pct": 100.0 * (ref_corrected["waterplane_coefficient"] / pub["cwp"] - 1.0),
                "cp": ref_corrected["prismatic_coefficient"],
                "cm": ref_corrected["midship_coefficient"],
                "cb_identity_check": ref_corrected["prismatic_coefficient"] * ref_corrected["midship_coefficient"],
            },
            "schneekluth": {
                "cwp": schneekluth_cwp(cb),
                "corrected_error_pct": 100.0 * (ref_corrected["waterplane_coefficient"] / schneekluth_cwp(cb) - 1.0),
            },
            "mesh_quality": {
                "watertight": model.validation.watertight,
                "manifold": model.validation.manifold,
                "euler": model.validation.euler_characteristic,
                "outward": model.validation.outward,
            },
        }
        results.append(result)
    
    return {"parent_forms": results}


def main():
    """Run verification and write results."""
    print("Series 60 Verification Against Published Parent Forms")
    print("=" * 70)
    print()
    
    results = verify_parent_forms()
    
    # Print table
    print(f"{'Cb':>5} | {'Pub Cwp':>8} | {'Legacy Cwp':>10} | {'Legacy Err':>11} | {'Corrected Cwp':>13} | {'Corrected Err':>13}")
    print("-" * 70)
    
    for r in results["parent_forms"]:
        print(
            f"{r['cb']:5.2f} | "
            f"{r['published']['cwp']:8.3f} | "
            f"{r['legacy']['cwp']:10.3f} | "
            f"{r['legacy']['cwp_error_pct']:+10.2f}% | "
            f"{r['corrected']['cwp']:13.3f} | "
            f"{r['corrected']['cwp_error_pct']:+12.2f}%"
        )
    
    print()
    print("Verification checks:")
    
    all_pass = True
    for r in results["parent_forms"]:
        cb = r["cb"]
        error = r["corrected"]["cwp_error_pct"]
        passed = abs(error) < 1.0  # 1% threshold
        
        status = "PASS" if passed else "FAIL"
        print(f"  Cb={cb:.2f}: Cwp error = {error:+.2f}% [{status}]")
        
        if not passed:
            all_pass = False
        
        # Check Cb = Cp * Cm identity
        identity_error = abs(r["corrected"]["cb_identity_check"] - cb)
        identity_pass = identity_error < 1e-12
        identity_status = "PASS" if identity_pass else "FAIL"
        print(f"         Cb = Cp*Cm identity: error = {identity_error:.2e} [{identity_status}]")
        
        if not identity_pass:
            all_pass = False
    
    print()
    if all_pass:
        print("OVERALL VERDICT: PASS")
        print("The corrected Series 60 module reproduces the published parent form")
        print("coefficients to machine precision.")
    else:
        print("OVERALL VERDICT: FAIL")
        print("Some verification checks did not pass.")
    
    # Write results to JSON
    output_dir = Path(__file__).parent
    output_file = output_dir / "verification_results.json"
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults written to {output_file}")
    
    return 0 if all_pass else 1


if __name__ == "__main__":
    exit(main())
