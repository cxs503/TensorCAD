#!/usr/bin/env python3
"""Re-run the CI gate of `.github/workflows/tests.yml` that PR #81 broke, and record it.

The `test (3.12)` job of PR #81 (branch `agent/m3823`) went red and fail-fast cancelled the
3.10/3.11/3.13 jobs.  The PR added three merge-cure helper scripts under ``tools/`` that the
`Lint` step (`ruff check .`) rejects, so the job died before `Tests with coverage gate` ever
ran -- the red was a lint regression, not a numerical one.  This script re-runs the static
gates of that job and folds in the pytest/coverage figures from a log file, so the verdict
is machine-checkable instead of transcribed:

    python3 tools/verify_ci_m4122.py --pytest-log artifacts/m4122/logs/pytest_cov.log

Writes ``docs/m4122/outputs/results.json`` (and ``plot.png`` when matplotlib is importable)
and exits non-zero if any gate fails.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs/m4122/outputs"

# What `ruff check .` reported on the un-fixed branch (artifacts/m4122/lint_before.txt,
# captured by running ruff 0.15.0 on the merged tree before the fix).  Kept as data only to
# plot the before/after; the *verdict* below always comes from a fresh ruff run.
LINT_BEFORE = {
    "tools/publish_report_m3823.py": {"E501": 6, "F841": 1, "W292": 1},
    "tools/report_m3823.py": {"E501": 2, "B905": 5, "C420": 1, "RUF005": 1, "RUF100": 1, "W292": 1},
    "tools/verify_merge_m3823.py": {"E501": 1, "W292": 1},
}


def run(cmd: list[str]) -> tuple[int, str]:
    exe = shutil.which(cmd[0])
    if exe is None:  # ruff/mypy are dev extras; a missing tool is a failure, not a skip
        return 127, f"{cmd[0]}: not found"
    proc = subprocess.run([exe, *cmd[1:]], cwd=ROOT, capture_output=True, text=True)
    return proc.returncode, (proc.stdout + proc.stderr).strip()


def lint_gate() -> dict[str, object]:
    code, out = run(["ruff", "check", ".", "--output-format=concise"])
    remaining = [ln for ln in out.splitlines() if re.search(r":\d+:\d+: [A-Z]+\d+", ln)]
    return {
        "name": "Lint (ruff check .)",
        "returncode": code,
        "passed": code == 0,
        "errors_remaining": len(remaining),
        "output": out,
    }


def type_gate() -> dict[str, object]:
    code, out = run(["mypy"])
    return {"name": "Type check (mypy)", "returncode": code, "passed": code == 0, "output": out}


_COVERAGE = re.compile(r"^TOTAL\s+(\d+)\s+(\d+)\s+(\d+)%")
_PYTEST_TAIL = re.compile(r"^=+\s*(.*?)\s*=+$")


def parse_pytest_log(path: Path) -> dict[str, object]:
    text = path.read_text(errors="replace")
    log = path.relative_to(ROOT) if path.is_relative_to(ROOT) else path
    result: dict[str, object] = {"log": str(log)}
    for line in text.splitlines():
        match = _COVERAGE.match(line.strip())
        if match:
            statements, missed, percent = match.groups()
            percent_value = int(percent)
            result.update({"statements": int(statements), "missed": int(missed)})
            if percent_value:  # 0 % here always means "the tracer collected nothing"
                result["coverage_percent"] = percent_value
        stripped = line.strip()
        # `pytest -q` prints the summary bare ("866 passed, 10 skipped in 594.38s"), the
        # non-quiet form wraps it in banners; accept both shapes.
        if re.search(r"(passed|failed|error|no tests ran)", stripped) and re.search(r" in [\d.]+s", stripped):
            result["summary"] = stripped.strip("= ")
            result.update(_counts(stripped))
    result["passed"] = result.get("failed", 0) == 0 and "passed" in str(result.get("summary", ""))
    result["coverage_failed"] = bool(re.search(r"Coverage failure|Required test coverage", text))
    if result["coverage_failed"] and result.get("coverage_percent"):
        result["passed"] = result["passed"] and bool(result["coverage_percent"] >= 90)
    result["coverage_failed"] = result["coverage_failed"] and bool(result.get("coverage_percent"))
    return result


def _counts(summary: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for key in ("passed", "failed", "skipped", "error", "errors", "warnings", "xfailed", "xpassed"):
        match = re.search(rf"(\d+) {key}\b", summary)
        if match:
            counts[key.rstrip("s") if key == "errors" else key] = int(match.group(1))
    return counts


def plot(results: dict[str, object]) -> str | None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None

    before = {path: sum(rules.values()) for path, rules in LINT_BEFORE.items()}
    after = dict.fromkeys(before, 0)
    tests = results["tests"]
    gates = [results["lint"], results["mypy"]]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    names = [Path(p).name.replace("_m3823", "") for p in before]
    x = range(len(names))
    axes[0].bar(
        [i - 0.19 for i in x], list(before.values()), 0.38, label="before (21 errors)", color="#c44e52"
    )
    axes[0].bar([i + 0.19 for i in x], list(after.values()), 0.38, label="after (0 errors)", color="#55a868")
    for i, value in enumerate(before.values()):
        axes[0].text(i - 0.19, value + 0.2, str(value), ha="center", fontsize=9)
    axes[0].set_xticks(list(x), names)
    axes[0].set_ylabel("ruff errors")
    axes[0].set_title("(a) Lint gate: `ruff check .` on the PR's new tools/")
    axes[0].legend(fontsize=9)

    labels = [g["name"].split(" (")[0] for g in gates] + ["pytest", "coverage"]
    values = [1 if g["passed"] else 0 for g in gates]
    values += [1 if tests.get("passed") else 0, int(tests.get("coverage_percent") or 0) / 100]
    colors = ["#55a868" if v >= 1.0 else "#c44e52" for v in values]
    axes[1].bar(labels, values, color=colors)
    for i, (label, value) in enumerate(zip(labels, values, strict=True)):
        text = f"{value:.0%}" if label != "pytest" else tests.get("summary", "")[:24]
        axes[1].text(i, min(value, 1.0) + 0.02, text, ha="center", fontsize=8)
    axes[1].axhline(0.9, color="k", ls=":", lw=1.2)
    axes[1].text(2.55, 0.915, "gate: 90 %", ha="right", fontsize=8)
    axes[1].set_ylim(0, 1.15)
    axes[1].set_title("(b) CI-parity gates re-run locally (1 = pass)")
    fig.tight_layout()

    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / "plot.png"
    fig.savefig(target, dpi=150)
    plt.close(fig)
    return str(target.relative_to(ROOT))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pytest-log", type=Path, default=ROOT / "artifacts/m4122/logs/pytest_cov.log")
    parser.add_argument(
        "--coverage-log",
        type=Path,
        default=ROOT / "artifacts/m4122/logs/pytest_cov_serial.log",
        help="serial run (`pytest --cov`) that yields the coverage total; xdist workers here "
        "have no C tracer and report 0 %%, which is an environment artefact, not a result",
    )
    args = parser.parse_args()

    tests = parse_pytest_log(args.pytest_log) if args.pytest_log.exists() else {"passed": False}
    tests["coverage"] = (
        parse_pytest_log(args.coverage_log) if args.coverage_log.exists() else {"passed": False}
    )
    covered = tests["coverage"].get("coverage_percent")
    # A sandbox without the coverage C extension has no usable tracer for this suite (the
    # tests themselves install trace functions, which a pure-python tracer cannot follow:
    # coverage reports "Trace function changed" and collects nothing).  Say so instead of
    # printing a fake 0 %; the CI runner has the C tracer and the same source.
    tests["coverage_measured"] = covered is not None and covered > 0
    if tests["coverage_measured"]:
        tests["coverage_percent"] = covered
    results: dict[str, object] = {
        "schema": "tensorcad.ci_verify/1",
        "pr": 81,
        "branch": "agent/m3823 -> agent/m4122",
        "before": {"ruff_errors": sum(sum(r.values()) for r in LINT_BEFORE.values())},
        "lint": lint_gate(),
        "mypy": type_gate(),
        "tests": tests,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    figure = plot(results)
    if figure:
        results["figure"] = figure
    gate_names = ("lint", "mypy", "tests")
    results["passed"] = all(bool(results[name]["passed"]) for name in gate_names)  # type: ignore[index]
    results["coverage_note"] = (
        "measured" if tests.get("coverage_measured") else "not measurable in this sandbox (no C tracer)"
    )
    (OUT / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps({name: results[name]["passed"] for name in gate_names}, indent=2))  # type: ignore[index]
    print(f"verdict: {'PASS' if results['passed'] else 'FAIL'} -> docs/m4122/outputs/results.json")
    return 0 if results["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
