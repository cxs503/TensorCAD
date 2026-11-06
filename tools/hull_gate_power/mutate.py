#!/usr/bin/env python3
"""Gate-power audit of the parametric hull generator's closed forms.

Continues the 2026-11-04 claim-gate-power audit
(``docs/analysis/2026-11-04-claim-gate-power.md``), which explicitly left
``hull.py`` uncovered: "the generator's own closed forms are shared constants --
a mutant of the exponent would have to be paired with the analytic side to be
meaningful". This harness performs that pairing: one single-defect mutant at a
time, injected by exact-string replacement into ``tensorcad/hull.py`` (the
generator) or ``tensorcad/hydrostatics_analytic.py`` (the analytic reference),
with the test files that guard the documented claim run against each mutant.

A mutant the suite does not notice is a documented number with no teeth.

Method (same discipline as the 2026-11-04 audit):

* one defect at a time, patched by exact-string replacement (a pattern that is
  not unique is refused, not patched);
* a clean-tree guard aborts if ``git diff -- tensorcad`` is non-empty, because
  an interrupted run would leave one mutant in place;
* the file is restored from its original bytes with a SHA-256 check before the
  next mutant starts;
* each mutant runs only the test files that guard its claim, so a kill is
  attributable; baselines (unmutated) are recorded for the same file sets.

Usage::

    python3 tools/hull_gate_power/mutate.py --out results.json [--only HUL-1,ANA-2]

The run is resumable: mutants already present in ``--out`` are skipped.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HULL = ROOT / "tensorcad" / "hull.py"
ANA = ROOT / "tensorcad" / "hydrostatics_analytic.py"

GEN_TESTS = ["tests/test_hull.py", "tests/test_analytic_hydrostatics.py"]
ANA_TESTS = ["tests/test_hydrostatics_analytic.py", "tests/test_analytic_hydrostatics.py", "tests/test_hull.py"]

MUTANTS: list[dict] = [
    {
        "id": "HUL-1-taper-exponent",
        "kind": "model",
        "file": "hull",
        "old": "        _breadth_exponent(parameters, index),\n    )",
        "new": "        0.8,  # mutant HUL-1\n    )",
        "claim": "mesh volume 7.446554 m^3 at 41x13 vs analytic 7.512698 (-0.88 %): the sampled surface must be the surface the reference integrates",
        "tests": GEN_TESTS,
    },
    {
        "id": "HUL-2-vertical-taper",
        "kind": "model",
        "file": "hull",
        "old": "    return math.sqrt(level / (parameters.nz - 1))",
        "new": "    return level / (parameters.nz - 1)  # mutant HUL-2",
        "claim": "the sqrt vertical profile (the 2/3 volume factor of B*L*T*I*2/3)",
        "tests": GEN_TESTS,
    },
    {
        "id": "HUL-3-end-zeroing",
        "kind": "gate",
        "file": "hull",
        "old": "    if index == 0 or index == parameters.nx - 1:\n        return 0.0",
        "new": "    if False:  # mutant HUL-3\n        return 0.0",
        "claim": "the end stations close on the centreline (math.sin(pi) is not exactly zero in floating point)",
        "tests": GEN_TESTS,
    },
    {
        "id": "HUL-4-midship-split",
        "kind": "gate",
        "file": "hull",
        "old": "    if 2 * index >= parameters.nx - 1:",
        "new": "    if 2 * index > parameters.nx - 1:  # mutant HUL-4",
        "claim": "midship belongs to both halves; the tie must not bend the surface",
        "tests": GEN_TESTS,
    },
    {
        "id": "ANA-1-volume-factor",
        "kind": "value",
        "file": "ana",
        "old": "    volume = beam * length * draft * taper * (2.0 / 3.0)",
        "new": "    volume = beam * length * draft * taper * (1.0 / 2.0)  # mutant ANA-1",
        "claim": "displaced volume 7.512698 m^3 = B*L*T*I*2/3",
        "tests": ANA_TESTS,
    },
    {
        "id": "ANA-2-vcb",
        "kind": "value",
        "file": "ana",
        "old": "        vertical_centre_buoyancy=-0.4 * draft,",
        "new": "        vertical_centre_buoyancy=-0.5 * draft,  # mutant ANA-2",
        "claim": "VCB = -0.4*T (0.6*T above the keel)",
        "claim_note": "the z-weighted mean of zeta**(1/2) on [-T, 0]",
        "tests": ANA_TESTS,
    },
    {
        "id": "ANA-3-waterplane-half-percent",
        "kind": "tolerance",
        "file": "ana",
        "old": "    waterplane = beam * length * taper",
        "new": "    waterplane = beam * length * taper * 0.995  # mutant ANA-3",
        "claim": "waterplane area 14.061669 m^2 at 41x13; the fine grid is within half a percent on every quantity",
        "tests": ANA_TESTS,
    },
    {
        "id": "ANA-4-taper-arg",
        "kind": "value",
        "file": "ana",
        "old": "math.gamma(0.5 * (exponent + 1.0))",
        "new": "math.gamma(0.5 * (exponent + 1.5))  # mutant ANA-4",
        "claim": "taper_integral is a Beta integral: e=0 -> 1, e=1 -> 2/pi, e=2 -> 1/2",
        "tests": ANA_TESTS,
    },
    {
        "id": "ANA-5-volume-half-percent",
        "kind": "tolerance",
        "file": "ana",
        "old": "    volume = beam * length * draft * taper * (2.0 / 3.0)",
        "new": "    volume = beam * length * draft * taper * (2.0 / 3.0) * 1.005  # mutant ANA-5",
        "claim": "the documented -0.880 % volume grid error, pinned to 5e-5",
        "tests": ANA_TESTS,
    },
    {
        "id": "HUL-5-taper-integral-const",
        "kind": "value",
        "file": "hull",
        "old": "TAPER_INTEGRAL = (2.0 / math.pi) * math.sqrt(math.pi) * math.gamma(0.85) / (2.0 * math.gamma(1.35))",
        "new": "TAPER_INTEGRAL = (1.0 / math.pi) * math.sqrt(math.pi) * math.gamma(0.85) / (2.0 * math.gamma(1.35))  # mutant HUL-5",
        "claim": "TAPER_INTEGRAL = int_0^1 sin(pi*x)^0.7 dx, cross-derived from the volume identity",
        "tests": ANA_TESTS,
    },
]


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def clean_tree() -> bool:
    out = subprocess.run(
        ["git", "diff", "--quiet", "--", "tensorcad"], cwd=ROOT, capture_output=True
    )
    return out.returncode == 0


def run_pytest(tests: list[str]) -> dict:
    started = time.time()
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", *tests, "-q", "--no-header", "-p", "no:cacheprovider"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    tail = (proc.stdout + proc.stderr).strip().splitlines()
    summary = next((line for line in reversed(tail) if "passed" in line or "error" in line), "")
    failed = [
        line.strip()
        for line in (proc.stdout + proc.stderr).splitlines()
        if line.startswith("FAILED") or line.startswith("ERROR")
    ]
    return {
        "returncode": proc.returncode,
        "killed": proc.returncode != 0,
        "summary": summary,
        "failed": failed,
        "wall_seconds": round(time.time() - started, 2),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    parser.add_argument("--only", default="")
    args = parser.parse_args()

    out_path = ROOT / args.out
    results: dict = {}
    if out_path.exists():
        results = json.loads(out_path.read_text())

    originals = {
        "hull": HULL.read_bytes(),
        "ana": ANA.read_bytes(),
    }
    digests = {key: sha256(data) for key, data in originals.items()}

    selected = [m for m in MUTANTS if not args.only or m["id"] in args.only.split(",")]
    for mutant in selected:
        mid = mutant["id"]
        if mid in results:
            print(f"[skip] {mid} already recorded", flush=True)
            continue
        path = HULL if mutant["file"] == "hull" else ANA
        original = originals[mutant["file"]]
        if not clean_tree():
            print(f"[abort] dirty tree before {mid}; restoring and stopping", flush=True)
            HULL.write_bytes(originals["hull"])
            ANA.write_bytes(originals["ana"])
            return 1
        text = original.decode()
        if text.count(mutant["old"]) != 1:
            results[mid] = {
                "status": "pattern-not-unique",
                "occurrences": text.count(mutant["old"]),
                "mutant": mutant,
            }
            out_path.write_text(json.dumps(results, indent=2))
            print(f"[skip] {mid}: pattern occurs {text.count(mutant['old'])} times", flush=True)
            continue
        path.write_text(text.replace(mutant["old"], mutant["new"]))
        if sha256(path.read_bytes()) == digests[mutant["file"]]:
            print(f"[abort] {mid}: patch had no effect", flush=True)
            path.write_bytes(original)
            return 1
        run = run_pytest(mutant["tests"])
        path.write_bytes(original)
        if sha256(path.read_bytes()) != digests[mutant["file"]]:
            print(f"[abort] {mid}: restore failed", flush=True)
            return 1
        results[mid] = {"status": "ok", "mutant": mutant, "run": run}
        out_path.write_text(json.dumps(results, indent=2))
        verdict = "KILLED" if run["killed"] else "SURVIVED"
        print(f"[done] {mid}: {verdict} ({run['wall_seconds']}s) {run['summary']}", flush=True)

    # Baselines for the two test sets (recorded once).
    for name, tests in (("gen", GEN_TESTS), ("ana", ANA_TESTS)):
        if name not in results:
            results[f"baseline-{name}"] = {"status": "baseline", "run": run_pytest(tests)}
            out_path.write_text(json.dumps(results, indent=2))
            print(f"[done] baseline-{name}: {results[f'baseline-{name}']['run']['summary']}", flush=True)

    killed = sum(1 for r in results.values() if r.get("status") == "ok" and r["run"]["killed"])
    total = sum(1 for r in results.values() if r.get("status") == "ok")
    print(f"gate power: {killed}/{total} killed", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
