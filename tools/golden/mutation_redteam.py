#!/usr/bin/env python3
"""Mutation red-team of the golden geometry gate (GAP-CAD-D2-GOLDEN-REGRESSION).

A regression gate is only worth its runtime if it can fail. This harness takes the
committed tree, breaks one thing at a time in a scratch copy, and runs the *real* gate
(``python -m tensorcad.golden --check``) against it. Every mutant must be killed, i.e.
the gate must exit non-zero **and** report the finding code the mutation is supposed to
trip. A null control (the unmodified copy) must pass, otherwise the harness itself is
what fails.

The mutants are chosen to cover the three ways a golden suite can be vacuous:

``reference``
    the mutation changes a *reference source* (the closed form), so the stored
    baseline reference no longer equals what its declared source computes. If the gate
    only compared numbers it would never notice a stale/hand-copied criterion.
``numeric``
    the mutation changes the geometry that the integrated quantities are taken on
    (grid, case definition, mirror axis); the drift lock on the measurements must fire.
``identity``
    the mutation leaves every measured quantity bit-identical and only changes the
    mesh itself (inverted winding: the surface area takes the norm and the
    sphere/ellipsoid volume the absolute value). Only the geometry fingerprint can
    catch it -- this is the mutant that justifies the fingerprint assertion.

Usage::

    python tools/golden/mutation_redteam.py               # quick gate (default, CI)
    python tools/golden/mutation_redteam.py --full        # full two-grid gate
    python tools/golden/mutation_redteam.py --only orientation-winding --json out.json
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASELINE_REL = Path("docs/golden/tensorcad_geometry_golden.json")
DATA_REL = Path("docs/data/dtmb5415")
SCHEMA = "tensorcad.golden.redteam/1"


class HarnessError(RuntimeError):
    """The harness could not apply a mutation: not a verdict about the gate."""


@dataclass(frozen=True)
class Mutant:
    """One deliberate defect: a unique textual patch plus the codes it must trip."""

    id: str
    kind: str  # "reference" | "numeric" | "identity"
    target: str
    find: str
    replace: str
    expected: tuple[str, ...]
    rationale: str
    expected_full: tuple[str, ...] | None = None


MUTANTS: tuple[Mutant, ...] = (
    Mutant(
        id="reference-wigley-waterplane",
        kind="reference",
        target="tensorcad/wigley.py",
        find="    return length * beam * n / (n + 1.0)",
        replace="    return length * beam * n / (n + 1.001)",
        expected=("reference_not_reproducible",),
        rationale=(
            "the closed form of the Wigley waterplane area is perturbed by 1e-3; the "
            "stored reference no longer equals what its declared source computes, so a "
            "hand-copied criterion cannot survive a source change"
        ),
    ),
    Mutant(
        id="reference-series60-waterplane",
        kind="reference",
        target="tensorcad/series60.py",
        find="    return L * B * n / (n + 1.0)",
        replace="    return L * B * n / (n + 1.002)",
        expected=("reference_not_reproducible",),
        rationale="same proof on the second analytic family (Series 60 exponents)",
    ),
    Mutant(
        id="numeric-wigley-grid-silently-coarsened",
        kind="numeric",
        target="tensorcad/golden.py",
        find='{"coarse": (41, 13), "fine": (161, 41)}',
        replace='{"coarse": (41, 13), "fine": (121, 31)}',
        expected=("geometry_changed",),
        expected_full=("drift", "geometry_changed"),
        rationale=(
            "the fine grid is silently coarsened without re-blessing: in full mode the "
            "fine measurements move (drift) and in quick mode, where only the coarse "
            "level is re-measured, the fine mesh fingerprint is what fails"
        ),
    ),
    Mutant(
        id="identity-orientation-winding",
        kind="identity",
        target="tensorcad/neural/extract.py",
        find="    faces[flip] = faces[flip][:, [0, 2, 1]]",
        replace="    faces[flip] = faces[flip]  # redteam: orientation repair disabled",
        expected=("geometry_changed",),
        rationale=(
            "marching tetrahedra stops repairing inverted winding: the surface area "
            "takes the norm, the sphere/ellipsoid volume the absolute value and the "
            "extents are unchanged, so every measured number is bit-identical and only "
            "the geometry fingerprint can fail"
        ),
    ),
    Mutant(
        id="numeric-sphere-radius",
        kind="numeric",
        target="tensorcad/golden.py",
        find="    RADIUS = 1.0",
        replace="    RADIUS = 1.0005",
        expected=("drift", "geometry_changed"),
        rationale=(
            "the case definition itself changes (sphere radius 0.05 % larger): the "
            "reference moves with it, so only the drift lock and the fingerprint fire"
        ),
    ),
    Mutant(
        id="numeric-5415-mirror-axis",
        kind="numeric",
        target="tensorcad/golden.py",
        find="    MIRROR_AXIS = 1",
        replace="    MIRROR_AXIS = 0",
        expected=("invariant_violated",),
        rationale=(
            "the case declares the mirror axis as x although the shipped shell mirrors "
            "about y: the measured extents do not change at all, so only the invariant "
            "identity can fail -- a check that is neither drift nor fingerprint"
        ),
    ),
    Mutant(
        id="baseline-reference-hand-edited",
        kind="reference",
        target=str(BASELINE_REL),
        find='"displacement_volume_m3": 2.7777777777777777',
        replace='"displacement_volume_m3": 2.79',
        expected=("reference_not_reproducible",),
        rationale=(
            "someone edits the blessed Wigley displacement by +0.4 % by hand; the gate "
            "recomputes the reference from the source and refuses the pass"
        ),
    ),
)


def stage(destination: Path) -> None:
    """Copy the parts of the tree the gate reads into ``destination``."""
    shutil.copytree(
        ROOT / "tensorcad",
        destination / "tensorcad",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    shutil.copytree(ROOT / DATA_REL, destination / DATA_REL)
    target = destination / BASELINE_REL
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ROOT / BASELINE_REL, target)


def apply_mutant(mutant: Mutant, root: Path) -> None:
    path = root / mutant.target
    if not path.is_file():
        raise HarnessError(f"{mutant.id}: {mutant.target} not present in the staged tree")
    text = path.read_text(encoding="utf-8")
    hits = text.count(mutant.find)
    if hits != 1:
        raise HarnessError(
            f"{mutant.id}: the patch anchor occurs {hits} times in {mutant.target}; "
            "a mutation must be unambiguous, so the harness refuses to guess"
        )
    path.write_text(text.replace(mutant.find, mutant.replace), encoding="utf-8")


def _finding_codes(report: dict) -> list[str]:
    codes = [str(finding.get("code")) for finding in report.get("findings") or []]
    for case in report.get("cases") or []:
        codes.extend(str(finding.get("code")) for finding in case.get("findings") or [])
    return sorted(set(codes))


def run_gate(root: Path, *, quick: bool) -> dict:
    """Run the committed gate inside ``root`` and return its verdict."""
    report_path = root / "gate.json"
    if report_path.exists():
        report_path.unlink()
    command = [sys.executable, "-m", "tensorcad.golden", "--check", "--json", str(report_path)]
    if quick:
        command.append("--quick")
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(root) + (
        os.pathsep + environment["PYTHONPATH"] if environment.get("PYTHONPATH") else ""
    )
    environment.pop("TENSORCAD_GOLDEN_SKIP", None)
    started = time.perf_counter()
    completed = subprocess.run(
        command, cwd=root, env=environment, capture_output=True, text=True, check=False
    )
    seconds = time.perf_counter() - started
    report: dict = {}
    if report_path.is_file():
        try:
            report = json.loads(report_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as failure:  # pragma: no cover - defensive
            report = {"unreadable": str(failure)}
    return {
        "returncode": completed.returncode,
        "seconds": round(seconds, 2),
        "codes": _finding_codes(report),
        "passed": bool(report.get("passed")),
        "output": (completed.stdout + completed.stderr)[-4000:],
    }


def expectations(mutant: Mutant, *, quick: bool) -> tuple[str, ...]:
    if not quick and mutant.expected_full is not None:
        return mutant.expected_full
    return mutant.expected


def run_redteam(
    *,
    quick: bool = True,
    only: tuple[str, ...] | None = None,
    keep: Path | None = None,
) -> dict:
    mutants = [mutant for mutant in MUTANTS if not only or mutant.id in only]
    started = time.perf_counter()
    workspace = Path(tempfile.mkdtemp(prefix="tensorcad-golden-redteam-"))
    results: list[dict] = []
    try:
        control_root = workspace / "control"
        control_root.mkdir()
        stage(control_root)
        control = run_gate(control_root, quick=quick)
        control["passed"] = control["returncode"] == 0 and control["passed"] is not False
        for mutant in mutants:
            case_root = workspace / mutant.id
            case_root.mkdir()
            stage(case_root)
            apply_mutant(mutant, case_root)
            outcome = run_gate(case_root, quick=quick)
            observed = outcome["codes"]
            expected = expectations(mutant, quick=quick)
            killed = outcome["returncode"] != 0 and set(expected) <= set(observed)
            results.append(
                {
                    "id": mutant.id,
                    "kind": mutant.kind,
                    "target": mutant.target,
                    "expected": list(expected),
                    "observed": observed,
                    "killed": killed,
                    "returncode": outcome["returncode"],
                    "seconds": outcome["seconds"],
                    "rationale": mutant.rationale,
                    "output_tail": "" if killed else outcome["output"],
                }
            )
            print(
                f"{'KILLED ' if killed else 'SURVIVED'} {mutant.id:<38} "
                f"expected={','.join(expected):<44} observed={','.join(observed) or '-'}",
                flush=True,
            )
        if keep is not None and keep.is_dir():
            shutil.copytree(workspace, keep, dirs_exist_ok=True)
    finally:
        if keep is None:
            shutil.rmtree(workspace, ignore_errors=True)
    killed = sum(1 for result in results if result["killed"])
    return {
        "schema": SCHEMA,
        "quick": quick,
        "seconds": round(time.perf_counter() - started, 2),
        "control": {
            "passed": bool(control.get("passed")),
            "returncode": control.get("returncode"),
            "seconds": control.get("seconds"),
            "codes": control.get("codes"),
        },
        "mutants": results,
        "killed": killed,
        "total": len(results),
        "kill_rate": (killed / len(results)) if results else 0.0,
        "passed": bool(control.get("passed")) and killed == len(results) and bool(results),
    }


def format_markdown(report: dict) -> str:
    lines = [
        "# Golden geometry gate: mutation red-team",
        "",
        f"- schema: `{report['schema']}`",
        f"- gate mode: {'quick (coarse level, as the git hook runs it)' if report['quick'] else 'full two-grid'}",
        f"- null control (unmodified copy): {'PASS' if report['control']['passed'] else 'FAIL'}",
        f"- mutants killed: **{report['killed']}/{report['total']}** "
        f"(kill rate {report['kill_rate']:.0%}), {report['seconds']} s wall clock",
        "",
        "| mutant | kind | patch site | expected finding | observed findings | verdict |",
        "|---|---|---|---|---|---|",
    ]
    for mutant in report["mutants"]:
        lines.append(
            f"| `{mutant['id']}` | {mutant['kind']} | `{mutant['target']}` | "
            f"{', '.join(f'`{code}`' for code in mutant['expected']) or '-'} | "
            f"{', '.join(f'`{code}`' for code in mutant['observed']) or '-'} | "
            f"{'KILLED' if mutant['killed'] else '**SURVIVED**'} |"
        )
    lines += ["", "## What each mutant proves", ""]
    for mutant in report["mutants"]:
        lines.append(f"- **{mutant['id']}** ({mutant['kind']}): {mutant['rationale']}")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--full", action="store_true", help="run the full two-grid gate")
    parser.add_argument("--only", default="", help="comma separated mutant ids")
    parser.add_argument("--json", metavar="PATH", help="write the machine-readable report here")
    parser.add_argument("--markdown", metavar="PATH", help="write the markdown report here")
    parser.add_argument(
        "--keep", metavar="DIR", help="keep the staged mutant trees in DIR (debugging)"
    )
    arguments = parser.parse_args(argv)
    only = tuple(item.strip() for item in arguments.only.split(",") if item.strip())
    keep = Path(arguments.keep) if arguments.keep else None
    if keep is not None:
        keep.mkdir(parents=True, exist_ok=True)
    report = run_redteam(quick=not arguments.full, only=only or None, keep=keep)
    if arguments.json:
        Path(arguments.json).write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    if arguments.markdown:
        Path(arguments.markdown).write_text(format_markdown(report), encoding="utf-8")
    print(
        f"red-team {'PASS' if report['passed'] else 'FAIL'}: "
        f"{report['killed']}/{report['total']} mutants killed, "
        f"control {'passed' if report['control']['passed'] else 'FAILED'}"
    )
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())