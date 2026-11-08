#!/usr/bin/env python3
"""Gate-power audit of TensorCAD's *regulatory* numbers.

Continues the gate-power line of `docs/analysis/2026-11-04-claim-gate-power.md` and
`docs/analysis/2026-11-06-hull-gate-power.md`. Those two audits covered the hydrostatic
closed forms and the parametric hull generator. Neither touched the three modules whose
numbers are *verdicts*: `tensorcad.statutory` (IS Code 2008 part A s. 2.2), the weather
criterion in `tensorcad.intact_stability` (s. 2.3) and `tensorcad.damage` (lost buoyancy).

That is the highest-stakes verification debt in the repository, because the failure mode is
not a number being slightly off -- it is a ship being signed off. A test that re-asserts
whatever the implementation does cannot detect a threshold typo, and a doc that quotes
`0.055 m rad` is not evidence that the code compares against `0.055`.

So the claims are attacked instead of read: one single-defect mutant at a time, injected by
exact-string replacement (a pattern that occurs zero or several times is refused, not
patched), with the test files that are supposed to protect the claim run against it. A
mutant the suite does not notice is a rule with no teeth.

Method (same discipline as the two earlier audits):

* one defect at a time, exact-string replacement, pattern must be unique;
* a clean-tree guard aborts if `git diff -- tensorcad` is non-empty, because an interrupted
  run would leave one mutant in place and silently contaminate the next;
* the file is restored from its original bytes with a SHA-256 check before the next mutant;
* only the test files that guard the claim are run, so a kill is attributable;
* baselines (unmutated) are recorded for the same file sets, so "killed" always means
  "red *because of* the mutation";
* a `--phase` label keeps the before/after of a remediation in one JSON instead of
  overwriting the finding it fixed.

Usage::

    python3 tools/regulatory_gate_power/mutate.py --out runs/regulatory_gate_power.json \
        --phase before [--only STA-1,IS-2]

The run is resumable: mutants already present in `--out` under the same phase are skipped.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

STAT = ROOT / "tensorcad" / "statutory.py"
IS = ROOT / "tensorcad" / "intact_stability.py"
DMG = ROOT / "tensorcad" / "damage.py"

FILES = {"stat": STAT, "is": IS, "dmg": DMG}

STAT_TESTS = ["tests/test_statutory.py"]
IS_TESTS = ["tests/test_intact_stability.py", "tests/test_intact_stability_validation.py"]
DMG_TESTS = ["tests/test_damage.py"]

#: Every mutant is a *plausible* engineering error, not random noise: a drifted threshold, a
#: dropped clamp, a dropped balance equation, a digit transposition in a quoted constant, a
#: tolerance that quietly turns a converged number into an approximate one.
MUTANTS = [
    # --- IS Code 2008 part A, section 2.2 (tensorcad.statutory) -------------------------
    {
        "id": "STA-1-threshold-2.2.1a",
        "kind": "value",
        "file": "stat",
        "old": "AREA_TO_30_MIN_M_RAD = 0.055",
        "new": "AREA_TO_30_MIN_M_RAD = 0.035  # mutant STA-1",
        "claim": "2.2.1a requires 0.055 m rad under GZ to 30 deg",
        "tests": STAT_TESTS,
    },
    {
        "id": "STA-2-flooding-angle-clamp",
        "kind": "gate",
        "file": "stat",
        "old": (
            "    if flooding_angle_deg is not None:\n"
            "        ceiling = min(ceiling, flooding_angle_deg)"
        ),
        "new": (
            "    if False:  # mutant STA-2\n"
            "        ceiling = min(ceiling, flooding_angle_deg)"
        ),
        "claim": "2.2.1b/c integrate to theta_f when it is below 40 deg, not always to 40",
        "tests": STAT_TESTS,
    },
    {
        "id": "STA-3-indeterminate-band",
        "kind": "gate",
        "file": "stat",
        "old": "    if uncertainty > floor and abs(margin) <= uncertainty:",
        "new": "    if False and uncertainty > floor and abs(margin) <= uncertainty:  # mutant STA-3",
        "claim": "a margin inside the stated numerical uncertainty is indeterminate, not pass",
        "tests": STAT_TESTS,
    },
    {
        "id": "STA-4-peak-angle-25",
        "kind": "value",
        "file": "stat",
        "old": "ANGLE_OF_MAX_GZ_MIN_DEG = 25.0",
        "new": "ANGLE_OF_MAX_GZ_MIN_DEG = 20.0  # mutant STA-4",
        "claim": "2.2.3 requires the heel of maximum GZ to be at least 25 deg",
        "tests": STAT_TESTS,
    },
    {
        "id": "STA-5-area-quadrature",
        "kind": "tolerance",
        "file": "stat",
        "old": "DEFAULT_AREA_TOLERANCE = 1e-11",
        "new": "DEFAULT_AREA_TOLERANCE = 1e-3  # mutant STA-5",
        "claim": "the areas are adaptive-Simpson converged to 1e-11 m rad, so they do not decide a verdict",
        "tests": STAT_TESTS,
    },
    # --- IS Code 2008 part A, section 2.3 weather criterion (tensorcad.intact_stability) --
    {
        "id": "IS-1-gust-factor",
        "kind": "value",
        "file": "is",
        "old": "WIND_GUST_FACTOR = 1.5",
        "new": "WIND_GUST_FACTOR = 1.0  # mutant IS-1",
        "claim": "2.3.2 applies a gust factor of 1.5 to the steady wind heeling lever",
        "tests": IS_TESTS,
    },
    {
        "id": "IS-2-wind-pressure-504",
        "kind": "value",
        "file": "is",
        "old": "WIND_PRESSURE_PA = 504.0",
        "new": "WIND_PRESSURE_PA = 540.0  # mutant IS-2 (digit transposition)",
        "claim": "the steady wind heeling lever uses 504 Pa (lw1 = P A Z / (1000 g Delta))",
        "tests": IS_TESTS,
    },
    {
        "id": "IS-3-deck-edge-reduction",
        "kind": "gate",
        "file": "is",
        "old": "        limit = min(limit, 0.8 * deck_edge)",
        "new": "        limit = limit  # mutant IS-3",
        "claim": "the 2.3.1 heel limit is min(16 deg, 0.8 x deck-edge immersion)",
        "tests": IS_TESTS,
    },
    {
        "id": "IS-4-steady-heel-16",
        "kind": "value",
        "file": "is",
        "old": "WEATHER_MAX_STEADY_HEEL_DEG = 16.0",
        "new": "WEATHER_MAX_STEADY_HEEL_DEG = 20.0  # mutant IS-4",
        "claim": "the steady wind heel shall not exceed 16 deg",
        "tests": IS_TESTS,
    },
    {
        "id": "IS-5-theta2-cap-50",
        "kind": "value",
        "file": "is",
        "old": "WEATHER_THETA2_CAP_DEG = 50.0",
        "new": "WEATHER_THETA2_CAP_DEG = 40.0  # mutant IS-5",
        "claim": "2.3.2 caps theta2 at 50 deg in addition to down-flooding / second intercept",
        "tests": IS_TESTS,
    },
    # --- lost-buoyancy damage stability (tensorcad.damage) -------------------------------
    {
        "id": "DAM-1-heel-balance",
        "kind": "gate",
        "file": "dmg",
        "old": (
            "    if tcb is not None:\n"
            "        residuals.append((integrals.tcb - tcb) / sections.beam)"
        ),
        "new": (
            "    if False:  # mutant DAM-1\n"
            "        residuals.append((integrals.tcb - tcb) / sections.beam)"
        ),
        "claim": "heel is an output: the net volume's centroid y must equal the ship's y_G",
        "tests": DMG_TESTS,
    },
    {
        "id": "DAM-2-permeability-waterplane",
        "kind": "value",
        "file": "dmg",
        "old": "                            cmp_area += mu * weight * part",
        "new": "                            cmp_area += weight * part  # mutant DAM-2",
        "claim": "the net waterplane subtracts mu x flooded area (and no separate free-surface "
        "correction is added, which would double count)",
        "tests": DMG_TESTS,
    },
    {
        "id": "DAM-3-moment-tolerance",
        "kind": "tolerance",
        "file": "dmg",
        "old": "MOMENT_TOLERANCE = 1e-10",
        "new": "MOMENT_TOLERANCE = 1e-4  # mutant DAM-3",
        "claim": "the equilibrium solve is converged to 1e-10 relative on every balance",
        "tests": DMG_TESTS,
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
    combined = proc.stdout + proc.stderr
    tail = combined.strip().splitlines()
    summary = next((line for line in reversed(tail) if "passed" in line or "error" in line), "")
    failed = [
        line.strip()
        for line in combined.splitlines()
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
    parser.add_argument("--phase", default="before")
    parser.add_argument("--only", default="")
    parser.add_argument("--skip-baseline", action="store_true")
    args = parser.parse_args()

    out_path = ROOT / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    document: dict = {}
    if out_path.exists():
        document = json.loads(out_path.read_text())
    phase = document.setdefault(args.phase, {})

    originals = {key: path.read_bytes() for key, path in FILES.items()}
    digests = {key: sha256(data) for key, data in originals.items()}

    def restore_all() -> None:
        for key, data in originals.items():
            FILES[key].write_bytes(data)
            assert sha256(FILES[key].read_bytes()) == digests[key], f"restore failed for {key}"

    if not args.skip_baseline:
        for label, tests in (
            ("stat", STAT_TESTS),
            ("is", IS_TESTS),
            ("dmg", DMG_TESTS),
        ):
            if label in phase.get("_baseline", {}):
                continue
            result = run_pytest(tests)
            phase.setdefault("_baseline", {})[label] = {**result, "tests": tests}
            print(f"[baseline] {label}: {result['summary']} ({result['wall_seconds']}s)", flush=True)
            out_path.write_text(json.dumps(document, indent=2))
            if result["killed"]:
                print(f"[abort] baseline for {label} is red; a kill would not be attributable")
                return 1

    selected = [m for m in MUTANTS if not args.only or m["id"] in args.only.split(",")]
    for mutant in selected:
        mid = mutant["id"]
        if mid in phase:
            print(f"[skip] {mid} already recorded in phase {args.phase}", flush=True)
            continue
        if not clean_tree():
            restore_all()
            print(f"[abort] tree dirty before {mid}; restored and stopping", flush=True)
            return 1
        path = FILES[mutant["file"]]
        text = originals[mutant["file"]].decode()
        occurrences = text.count(mutant["old"])
        if occurrences != 1:
            restore_all()
            phase[mid] = {
                "status": "pattern-not-unique",
                "occurrences": occurrences,
                "mutant": {k: v for k, v in mutant.items() if k != "tests"},
            }
            out_path.write_text(json.dumps(document, indent=2))
            print(f"[skip] {mid}: pattern occurs {occurrences} times", flush=True)
            continue
        path.write_text(text.replace(mutant["old"], mutant["new"]))
        try:
            result = run_pytest(mutant["tests"])
        finally:
            restore_all()
        phase[mid] = {
            "status": "killed" if result["killed"] else "survived",
            "kind": mutant["kind"],
            "file": mutant["file"],
            "claim": mutant["claim"],
            "tests": mutant["tests"],
            "old": mutant["old"],
            "new": mutant["new"],
            "summary": result["summary"],
            "failed": result["failed"],
            "wall_seconds": result["wall_seconds"],
        }
        out_path.write_text(json.dumps(document, indent=2))
        print(f"[{phase[mid]['status'].upper():8s}] {mid}: {result['summary']}", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())