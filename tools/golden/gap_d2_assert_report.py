#!/usr/bin/env python3
"""Assert the acceptance criteria of the GAP-CAD-D2 golden run (GAP-CAD-D2-GOLDEN-REGRESSION).

Both carriers of the gate -- the CI workflow and the pre-commit hook -- have to decide
"is this green?" from the same reading of the run, so the reading lives in one file
instead of in two shell heredocs that can drift.

    python3 tools/golden/gap_d2_assert_report.py <suite-dir>

``<suite-dir>`` is the TensorPapers ``papers/gap-cad-d2-golden-regression`` checkout; the
script reads the ``outputs/golden_suite.json`` and ``outputs/mutation_report.json`` that
the run just wrote.  Exit 0 = accepted, exit 1 = rejected with a one-line reason.

Why not just trust the exit codes of ``suite.py`` / ``mutants.py``:

* ``suite.py`` exits non-zero for a FAIL verdict only.  A case whose geometric gate fails
  (the generated hull no longer sits on its own analytic surface) is reported as
  SKIP/UNAVAILABLE and the script still exits 0 -- measured: +3 % half-breadth in
  ``tensorcad/wigley.py`` gives ``suite.py`` exit 0 with two cases SKIP
  (``artifacts/gap-cad-d2-golden/regression_control.json``).
* ``mutants.py`` exits non-zero, but the *contract* is about rates (control arm 8/8,
  detection power 1.0, WITHHELD still refused), and a rate can fall while the exit code
  stays 0 if the arm's own bookkeeping changes.

So the gate reads the report, not just the exit status.  That is the whole point of the
acceptance criteria being numbers.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

EXPECTED_CASES = 8


class Rejected(Exception):
    """A run that must not be allowed through."""


def _load(path: str) -> dict:
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def check_suite(report: dict) -> str:
    cases = report.get("results", [])
    if len(cases) != EXPECTED_CASES:
        raise Rejected(
            f"expected {EXPECTED_CASES} golden cases, found {len(cases)}: "
            f"{[c.get('case') for c in cases]}"
        )
    unmeasured = [
        (case["case"], case["verdict"])
        for case in cases
        if not str(case["verdict"]).startswith("PASS")
    ]
    if unmeasured:
        raise Rejected(f"golden case(s) not measured or failed: {unmeasured}")
    verdicts = [f["verdict"] for case in cases for f in case["findings"]]
    bad = sorted({v for v in verdicts if v not in ("PASS", "WITHHELD")})
    if bad:
        raise Rejected(f"unexpected finding verdicts: {bad}")
    return (
        f"{len(cases)} cases PASS; findings "
        f"{verdicts.count('PASS')} PASS / {verdicts.count('WITHHELD')} WITHHELD"
    )


def check_mutants(report: dict) -> str:
    summary = report.get("summary")
    if summary is None:
        raise Rejected("mutation report has no summary section")

    control = summary["control_arm"]
    if len(control) != EXPECTED_CASES or not all(
        str(v).startswith("PASS") for v in control.values()
    ):
        raise Rejected(f"control arm is not {EXPECTED_CASES}/{EXPECTED_CASES} PASS: {control}")
    if summary["control_arm_all_pass"] is not True:
        raise Rejected("control_arm_all_pass is not true")

    for key in ("strong_detection_rate", "gate_detection_rate"):
        if summary[key] != 1.0:
            raise Rejected(f"{key} = {summary[key]} (must be 1.0)")
    if summary["unexpected_misses"] != 0:
        raise Rejected(f"unexpected misses: {summary['misses']}")
    if summary["false_alarms"] != 0:
        raise Rejected(f"false alarms: {summary['false_alarm_examples']}")
    if summary["withheld_refused"] != summary["withheld_probes"]:
        raise Rejected(
            f"WITHHELD probes not refused: {summary['withheld_refused']}"
            f"/{summary['withheld_probes']}"
        )
    return (
        f"control arm {len(control)}/{len(control)} PASS, "
        f"detection power {summary['strong_detection_rate']:.2f} "
        f"({summary['strong_caught']}/{summary['must_catch_mutants']}), "
        f"WITHHELD refused {summary['withheld_refused']}/{summary['withheld_probes']}, "
        f"false alarms {summary['false_alarms']}"
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite_dir", help="papers/gap-cad-d2-golden-regression checkout")
    parser.add_argument("--suite-only", action="store_true",
                        help="judge only outputs/golden_suite.json (the sweep is next)")
    parser.add_argument("--mutants-only", action="store_true",
                        help="judge only outputs/mutation_report.json")
    args = parser.parse_args(argv)

    outputs = os.path.join(args.suite_dir, "outputs")
    if args.suite_only and args.mutants_only:
        parser.error("--suite-only and --mutants-only are mutually exclusive")
    try:
        suite_line = None if args.mutants_only else check_suite(
            _load(os.path.join(outputs, "golden_suite.json")))
        mutant_line = None if args.suite_only else check_mutants(
            _load(os.path.join(outputs, "mutation_report.json")))
    except Rejected as rejection:
        print(f"REJECTED: {rejection}", file=sys.stderr)
        return 1
    except (OSError, KeyError, TypeError) as error:
        print(f"REJECTED: cannot judge the run ({error})", file=sys.stderr)
        return 1

    if suite_line is not None:
        print(f"accepted: {suite_line}")
    if mutant_line is not None:
        print(f"accepted: {mutant_line}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())