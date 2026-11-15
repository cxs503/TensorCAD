#!/usr/bin/env python3
"""Before/after evidence for the recorded-criterion check (GAP-CAD-D2-GOLDEN-REGRESSION).

The gate stores, for every asserted metric, the tolerance it asserts against. In
``--quick`` mode -- the mode the git hook runs -- there is no coarse level, so that stored
band *is* the criterion. Nothing used to re-derive it, so a committed JSON edit could set
the criterion while every measured number, reference and fingerprint stayed honest.

This script reproduces that in two trees and runs the real gate on each:

``before``  ``tensorcad/golden.py`` as of ``--before-rev`` (the revision the fix is diffed
            against, so the claim "this used to pass" is bound to a commit, not a memory);
``after``   the worktree's ``tensorcad/golden.py``.

Each tree gets the same hand edit -- the Wigley block-coefficient band x100 -- and the gate
is run in quick mode (the hook's mode) and, for ``after``, in full mode (CI's mode). A
control run with the honest baseline is included per tree so a harness that simply fails
cannot be mistaken for a detection.

Usage::

    python tools/golden/criterion_provenance.py --before-rev HEAD~1 \
        --json docs/golden/criterion_provenance.json --png docs/golden/criterion_provenance.png
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASELINE_REL = Path("docs/golden/tensorcad_geometry_golden.json")
DATA_REL = Path("docs/data/dtmb5415")
CASE, METRIC, FACTOR = "wigley", "block_coefficient", 100.0
SCHEMA = "tensorcad.golden.criterion_provenance/1"


def _git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout


def _source(revision: str | None) -> tuple[bytes, dict]:
    if revision is None:
        text = (ROOT / "tensorcad" / "golden.py").read_bytes()
        ident = {"kind": "worktree", "dirty": bool(_git("status", "--porcelain").strip())}
    else:
        text = subprocess.run(
            ["git", "show", f"{revision}:tensorcad/golden.py"],
            cwd=ROOT, capture_output=True, check=True,
        ).stdout
        ident = {"kind": "git", "rev": revision, "commit": _git("rev-parse", revision).strip()}
    ident["sha256"] = hashlib.sha256(text).hexdigest()
    return text, ident


def _stage(destination: Path, source: bytes) -> None:
    shutil.copytree(
        ROOT / "tensorcad", destination / "tensorcad",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    shutil.copytree(ROOT / DATA_REL, destination / DATA_REL)
    target = destination / BASELINE_REL
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / BASELINE_REL, target)
    (destination / "tensorcad" / "golden.py").write_bytes(source)


def _inflate(root: Path) -> dict:
    path = root / BASELINE_REL
    data = json.loads(path.read_text(encoding="utf-8"))
    entry = data["cases"][CASE]["tolerances"][METRIC]
    honest = float(entry["value"])
    entry["value"] = honest * FACTOR
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"metric": f"{CASE}/{METRIC}", "derived": honest, "recorded": entry["value"]}


def _run_gate(root: Path, *, quick: bool) -> dict:
    report = root / "gate.json"
    command = [sys.executable, "-m", "tensorcad.golden", "--check", "--json", str(report)]
    if quick:
        command.append("--quick")
    completed = subprocess.run(command, cwd=root, capture_output=True, text=True, check=False)
    payload = json.loads(report.read_text(encoding="utf-8")) if report.is_file() else {}
    findings = list(payload.get("findings") or [])
    cases = payload.get("cases") or {}
    for case in (cases.values() if isinstance(cases, dict) else cases):
        findings.extend((case or {}).get("findings") or [])
    codes = sorted({str(f.get("code")) for f in findings})
    return {
        "quick": quick,
        "returncode": completed.returncode,
        "passed": payload.get("passed") is True and completed.returncode == 0,
        "codes": codes,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before-rev", required=True, help="revision to run as 'before'")
    parser.add_argument("--json", help="write the evidence here")
    parser.add_argument("--png", help="write a figure here")
    args = parser.parse_args(argv)

    measurements: dict[str, dict] = {}
    workspace = Path(tempfile.mkdtemp(prefix="tensorcad-criterion-"))
    try:
        for label, revision in (("before", args.before_rev), ("after", None)):
            source, ident = _source(revision)
            root = workspace / label
            root.mkdir()
            _stage(root, source)
            control = _run_gate(root, quick=True)
            edit = _inflate(root)
            runs = {"control_quick": control, "inflated_quick": _run_gate(root, quick=True)}
            if label == "after":
                runs["inflated_full"] = _run_gate(root, quick=False)
            measurements[label] = {"source": ident, "baseline_edit": edit, "runs": runs}
            print(
                f"{label:<7}{ident['sha256'][:12]}  control_quick="
                f"{'PASS' if control['passed'] else 'FAIL'}  "
                f"inflated_quick={'PASS' if runs['inflated_quick']['passed'] else 'FAIL'} "
                f"{runs['inflated_quick']['codes']}  "
                + (
                    f"inflated_full={'PASS' if runs['inflated_full']['passed'] else 'FAIL'} "
                    f"{runs['inflated_full']['codes']}"
                    if "inflated_full" in runs else ""
                ),
                flush=True,
            )
    finally:
        shutil.rmtree(workspace, ignore_errors=True)

    before = measurements["before"]["runs"]["inflated_quick"]
    after = measurements["after"]["runs"]["inflated_quick"]
    payload = {
        "schema": SCHEMA,
        "question": (
            "can a hand edit of the recorded acceptance band decide the golden gate while "
            "every measurement, reference and geometry fingerprint stays honest?"
        ),
        "baseline_edit": measurements["after"]["baseline_edit"],
        "measurements": measurements,
        "verdict": {
            "before_admits_the_edit": bool(before["passed"]),
            "after_refuses_the_edit": not bool(after["passed"]),
            "refusal_code_present": "tolerance_not_derived" in after["codes"],
        },
    }
    payload["passed"] = all(payload["verdict"].values()) and (
        measurements["before"]["runs"]["control_quick"]["passed"]
        and measurements["after"]["runs"]["control_quick"]["passed"]
        and not measurements["after"]["runs"]["inflated_full"]["passed"]
    )
    if args.json:
        Path(args.json).write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    if args.png:
        _figure(payload, Path(args.png))
    print("verdict:", json.dumps(payload["verdict"], sort_keys=True), "passed:", payload["passed"])
    return 0 if payload["passed"] else 1


def _figure(payload: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    edit = payload["baseline_edit"]
    quick = {label: payload["measurements"][label]["runs"]["inflated_quick"]["passed"]
             for label in ("before", "after")}
    full = payload["measurements"]["after"]["runs"]["inflated_full"]["passed"]
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.0, 4.2))
    left.bar(
        ["derived\n(2-grid)", f"recorded\n(hand edit x{FACTOR:g})"],
        [edit["derived"], edit["recorded"]], color=["#1f77b4", "#d62728"], width=0.55,
    )
    left.set_yscale("log")
    left.set_ylabel(f"acceptance band, {edit['metric']} (relative)")
    left.set_title("the criterion the gate asserts against")
    for index, value in enumerate((edit["derived"], edit["recorded"])):
        left.annotate(f"{value:.2e}", (index, value), textcoords="offset points",
                      xytext=(0, 6), ha="center", fontsize=9)
    left.grid(alpha=0.3, axis="y")

    labels = ["quick\ngit hook", "quick\ngit hook", "full\nCI"]
    verdicts = [quick["before"], quick["after"], full]
    names = ["before fix", "after fix", "after fix"]
    right.bar(labels, [1, 1, 1], color=["#d62728" if v else "#2ca02c" for v in verdicts], width=0.55)
    for index, (passed, name) in enumerate(zip(verdicts, names)):
        right.annotate(f"{'PASS -> hole' if passed else 'FAIL -> refused'}\n({name})",
                       (index, 1), textcoords="offset points", xytext=(0, -22),
                       ha="center", fontsize=9, color="white")
    right.set_yticks([])
    right.set_title("gate verdict on the inflated-band tree")
    right.grid(alpha=0.3, axis="y")
    figure.suptitle("golden gate: a recorded band must be its own derivation", fontsize=11)
    figure.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=150)
    plt.close(figure)


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
