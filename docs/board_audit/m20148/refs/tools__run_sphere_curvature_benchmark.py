#!/usr/bin/env python3
"""Integration caller: platform-side entry for sphere curvature benchmark.

This script is called by the TensorPlatform to run the CASE-CAD-CURVATURE-SPHERE
benchmark and return machine-readable results.

Usage:
    python tools/run_sphere_curvature_benchmark.py [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from benchmarks.verified.sphere_mean_curvature.run_benchmark import run as run_benchmark


def main() -> None:
    parser = argparse.ArgumentParser(description="Run sphere curvature benchmark")
    parser.add_argument("--out", type=str, default=None, help="Output directory")
    parser.add_argument("--json", action="store_true", help="Output result as JSON to stdout")
    args = parser.parse_args()

    out_dir = Path(args.out) if args.out else None
    result = run_benchmark(out_dir)

    if args.json:
        print(json.dumps(result, indent=2))
    else:
        s = result["summary"]
        print(f"CASE: {result['case']}")
        print(f"Gates: {s['passed']}/{s['total_gates']} passed")
        print(f"All hold: {result['all_gates_hold']}")
        print(f"Worst error: {s['worst_error_pct']:.4f}%")

    sys.exit(0 if result["all_gates_hold"] else 1)


if __name__ == "__main__":
    main()
