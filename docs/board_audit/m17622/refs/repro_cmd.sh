#!/usr/bin/env bash
# One-command reproduction of GAP-CAD-A3-CURVATURE-FAIRNESS (m17465 rework round).
#
#   bash runs/m17465/repro_cmd.sh
#
# Writes, under runs/m17465/:
#   benchmark/result.json + gates_verified.json + figures/   (reference vs measured)
#   cli_{sphere,wigley,hull}.json                            (product CLI records)
#   cli_via_package_entry.json                               (tensorcad-hull fairness ...)
#   crosscheck_three_cases.json                              (3-case third-party cross-check)
#   pytest.log, manifest.json                                (tests + sha256 of every artifact)
# Exit code is non-zero if any step fails; the benchmark's own exit code already
# encodes "all gates pass", so a broken claim cannot be reported as a green run.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="$ROOT/runs/m17465"
cd "$ROOT"
mkdir -p "$OUT"

echo "== benchmark (gates G1-G6 + diagnostic G7)"
python benchmarks/verified/curvature_fairness/run_benchmark.py --out "$OUT/benchmark" \
    | tee "$OUT/benchmark.log"

echo "== product CLI: python -m tensorcad.fairness_cli"
python -m tensorcad.fairness_cli sphere --radius 1 --level 24 --out "$OUT/cli_sphere.json" \
    > /dev/null
python -m tensorcad.fairness_cli wigley --nx 121 --nz 37 --out "$OUT/cli_wigley.json" \
    > /dev/null
python -m tensorcad.fairness_cli hull --out "$OUT/cli_hull.json" > /dev/null

echo "== package console entry: tensorcad-hull fairness (tensorcad/cli.py dispatch)"
python -m tensorcad.cli fairness wigley --nx 81 --nz 25 --out "$OUT/cli_via_package_entry.json" \
    > /dev/null

echo "== three-case cross-check (meshmagick WITHHELD -> trimesh + closed forms)"
python tools/crosscheck_curvature_thirdparty.py --out "$OUT/crosscheck_three_cases.json" \
    | tee "$OUT/crosscheck.log"

echo "== tests"
python -m pytest tests/test_fairness_energy.py tests/test_fairness_entry.py -q \
    | tee "$OUT/pytest.log"

echo "== manifest"
python "$OUT/make_manifest.py"

echo "done: $(ls "$OUT")"
