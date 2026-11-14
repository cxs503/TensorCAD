#!/usr/bin/env python3
"""Entry point for CASE-CAD-CURVATURE-SPHERE benchmark.

Usage:
    python -m tensorcad.benchmarks.sphere_curvature [--out DIR]

This is the canonical entry point for the sphere mean curvature verification case.
It runs all gates, produces result.json and gates_verified.json, and exits
non-zero if any gate fails.
"""
from __future__ import annotations

import sys
from pathlib import Path

# Ensure repo root is on path
_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from benchmarks.verified.sphere_mean_curvature.run_benchmark import main as benchmark_main

if __name__ == "__main__":
    benchmark_main()
