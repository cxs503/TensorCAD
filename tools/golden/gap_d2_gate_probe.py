#!/usr/bin/env python3
"""Exercise the pre-commit gate in its four decision modes and record the evidence.

The hook's contract has four outcomes, and three of them are about *not* paying for a
35-second suite: a documentation-only commit must pass through untouched, a geometry
edit must run the dual-derived suite, a named-but-missing suite must be refused, and a
checkout with no TensorPapers beside it must warn rather than block (CI is the arm that
fails closed).  This script reproduces each one with ``GOLDEN_STAGED_FILE`` -- the same
test seam the unit tests use -- and writes ``artifacts/gap-cad-d2-golden/hook_probe.json``
so the gate document quotes a run instead of a memory.

    python3 tools/golden/gap_d2_gate_probe.py --golden-dir /path/to/papers/gap-cad-d2-golden-regression
    python3 tools/golden/gap_d2_gate_probe.py            # geometry arm skipped if not found
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HOOK = os.path.join(REPO_ROOT, ".githooks", "pre-commit-geometry-golden.sh")
CHAIN = os.path.join(REPO_ROOT, ".githooks", "pre-commit")
OUT = os.path.join(REPO_ROOT, "artifacts", "gap-cad-d2-golden", "hook_probe.json")

DOC_ONLY = ["docs/baseline_fleet.md", "README.md"]
GEOMETRY = ["tensorcad/hull.py"]

DEFAULT_CANDIDATES = (
    os.path.join(REPO_ROOT, "..", "TensorPapers", "papers", "gap-cad-d2-golden-regression"),
    os.path.join(os.path.expanduser("~"), "TensorPapers", "papers", "gap-cad-d2-golden-regression"),
)


def _run(staged: list[str], golden_dir: str | None, entry: str = HOOK,
         tensorcad_path: str | None = None) -> dict:
    env = dict(os.environ, PYTHONPATH=tensorcad_path or REPO_ROOT)
    if tensorcad_path is not None:
        env["TENSORCAD_PATH"] = tensorcad_path
    env.pop("GOLDEN_DIR", None)
    if golden_dir is not None:
        env["GOLDEN_DIR"] = golden_dir
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as handle:
        handle.write("\n".join(staged) + "\n")
        staged_file = handle.name
    env["GOLDEN_STAGED_FILE"] = staged_file
    if entry == CHAIN:
        # the chain's first arm (tensorcad.golden) reads its own staging seam
        env["TENSORCAD_GOLDEN_STAGED_FILE"] = staged_file
    started = time.time()
    proc = subprocess.run([entry], cwd=REPO_ROOT, env=env, capture_output=True, text=True,
                          timeout=1800)
    elapsed = time.time() - started
    os.unlink(staged_file)
    tail = (proc.stdout + proc.stderr).strip().splitlines()
    return {
        "entry": os.path.relpath(entry, REPO_ROOT),
        "tensorcad_path": tensorcad_path or REPO_ROOT,
        "staged": staged,
        "golden_dir": golden_dir,
        "exit_code": proc.returncode,
        "wall_s": round(elapsed, 2),
        "last_line": tail[-1] if tail else "",
        # "[golden] ok" is printed only after both suite.py and mutants.py returned 0.
        "ran_suite": "[golden] ok" in proc.stdout,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--golden-dir", default=os.environ.get("GOLDEN_DIR", ""))
    parser.add_argument("--skip-suite", action="store_true",
                        help="skip the geometry arm (the ~35 s one)")
    parser.add_argument("--regression-path", default="",
                        help="a perturbed TensorCAD checkout; the hook must block on it")
    args = parser.parse_args(argv)

    suite_dir = args.golden_dir
    if not suite_dir:
        for candidate in DEFAULT_CANDIDATES:
            if os.path.isfile(os.path.join(candidate, "suite.py")):
                suite_dir = candidate
                break

    probe: dict[str, object] = {
        "hook": os.path.relpath(HOOK, REPO_ROOT),
        "chain": os.path.relpath(CHAIN, REPO_ROOT),
        "modes": [],
    }
    # (mode, staged paths, GOLDEN_DIR, entry point)
    modes = [
        ("chain_docs_only", DOC_ONLY, "/nonexistent-suite", CHAIN),
        ("hook_docs_only", DOC_ONLY, "/nonexistent-suite", HOOK),
    ]
    if suite_dir and not args.skip_suite:
        modes.append(("hook_geometry_edit", GEOMETRY, suite_dir, HOOK))
        modes.append(("chain_geometry_edit", GEOMETRY, suite_dir, CHAIN))
    if args.regression_path and suite_dir:
        modes.append(("hook_regression_blocked", GEOMETRY, suite_dir, HOOK))
    modes.append(("hook_named_suite_missing", GEOMETRY, "/nonexistent-suite", HOOK))
    if not any(os.path.isfile(os.path.join(c, "suite.py")) for c in DEFAULT_CANDIDATES):
        modes.append(("hook_no_suite_beside_checkout", GEOMETRY, None, HOOK))

    for name, staged, golden_dir, entry in modes:
        tensorcad_path = args.regression_path if name == "hook_regression_blocked" else None
        record = _run(staged, golden_dir, entry, tensorcad_path)
        record["mode"] = name
        probe["modes"].append(record)
        print(f"{name}: exit={record['exit_code']} wall={record['wall_s']}s :: {record['last_line']}")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as handle:
        json.dump(probe, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())