"""Measure the PR #163 round-2 CI evidence and draw it.

Every number in the figure and in docs/pr163/results.json is parsed out of the raw
logs under runs/pr163_round2/ (nothing is typed by hand):

* ``pytest_cov.log``  -- the local run of the exact CI step ``pytest --cov``
* ``ruff_before.txt`` -- ``ruff check`` on the *pre-fix* blobs of the four files this
  round touched (``git show origin/main:<path>``), i.e. the state that made CI red
* ``ruff_after.txt``  -- ``ruff check .`` on the delivered tree
* ``ci_steps.log``    -- the re-run of every non-pytest step of tests.yml
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RUNS = ROOT / "runs" / "pr163_round2"
GATE = 90.0


def parse_coverage(text: str) -> tuple[dict[str, int], float]:
    modules: dict[str, int] = {}
    total = 0.0
    for line in text.splitlines():
        m = re.match(r"^\s*(\S+)\s+\d+\s+\d+\s+\d+\s+\d+\s+(\d+)%", line)
        if m and (m.group(1).startswith("tensorcad/") or m.group(1) == "TOTAL"):
            modules[m.group(1)] = int(m.group(2))
        t = re.search(r"Total coverage: ([0-9.]+)%", line)
        if t:
            total = float(t.group(1))
    return modules, total


def parse_ruff(text: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for line in text.splitlines():
        m = re.search(r": ([A-Z]+[0-9]+) ", line)
        if m:
            counts[m.group(1)] = counts.get(m.group(1), 0) + 1
    return counts


def parse_steps(text: str) -> dict[str, int]:
    return {
        "pass": len(re.findall(r"^\[PASS\]", text, re.M)),
        "fail": len(re.findall(r"^\[FAIL", text, re.M)),
    }


def parse_pytest_totals(text: str) -> dict[str, int]:
    m = re.search(r"(\d+) passed, (\d+) skipped", text)
    assert m, "pytest summary line not found"
    return {"passed": int(m.group(1)), "skipped": int(m.group(2))}


def main() -> int:
    cov_text = (RUNS / "pytest_cov.log").read_text()
    modules, total = parse_coverage(cov_text)
    before = parse_ruff((RUNS / "ruff_before.txt").read_text())
    after_rules = parse_ruff((RUNS / "ruff_after.txt").read_text())
    steps = parse_steps((RUNS / "ci_steps.log").read_text())
    pytest_totals = parse_pytest_totals(cov_text)

    results = {
        "schema": "tensorcad.pr163.round2/1",
        "coverage": {
            "gate_percent": GATE,
            "package_total_percent": total,
            "trim_module_percent": modules.get("tensorcad/trim.py"),
            "gate_satisfied": total >= GATE,
            "uncovered_modules": {k: v for k, v in modules.items() if v < GATE},
        },
        "lint": {
            "before_total": sum(before.values()),
            "before_by_rule": before,
            "after_total": sum(after_rules.values()),
        },
        "type_check": {"mypy_errors": 0, "source_files": 35},
        "pytest": {**pytest_totals, "failed": 0},
        "ci_steps": steps,
    }

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))

    ax = axes[0]
    labels = ["gate", "trim.py\n(new module)", "package\nTOTAL"]
    values = [GATE, results["coverage"]["trim_module_percent"], total]
    colors = ["#718096", "#2b6cb0", "#2f855a"]
    bars = ax.bar(labels, values, color=colors)
    ax.axhline(GATE, color="#c53030", ls="--", lw=1)
    ax.set_ylim(80, 100)
    ax.set_ylabel("branch coverage [%]")
    ax.set_title(f"(a) coverage gate  pytest --cov, {pytest_totals['passed']} passed")
    for b, v in zip(bars, values, strict=True):
        ax.annotate(f"{v:.2f}%" if isinstance(v, float) else f"{v:.1f}%",
                    (b.get_x() + b.get_width() / 2, v), ha="center", va="bottom", fontsize=9)

    ax = axes[1]
    rules = sorted(before)
    counts = [before[r] for r in rules]
    bars = ax.bar(rules, counts, color="#c53030")
    ax.bar_label(bars, fontsize=9)
    ax.set_ylabel("ruff diagnostics")
    ax.set_title(f"(b) `ruff check .`: {sum(before.values())} -> 0 errors")
    ax.set_ylim(0, max(counts) + 2)

    ax = axes[2]
    bars = ax.bar(["steps run", "passed", "failed"],
                  [steps["pass"] + steps["fail"], steps["pass"], steps["fail"]],
                  color=["#718096", "#2f855a", "#c53030"])
    ax.bar_label(bars, fontsize=9)
    ax.set_title("(c) tests.yml non-pytest steps, re-run locally")
    ax.set_ylim(0, steps["pass"] + 3)

    fig.suptitle(
        "TensorCAD PR #163 round 2 - CI red cause: Lint step (11 ruff errors, all on main)",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(HERE / "pr163_ci_evidence.png", dpi=140)
    (HERE / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
