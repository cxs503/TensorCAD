#!/usr/bin/env python3
"""Publish the Series-60 gate-power record as documentation, from the record itself.

The anchor ``GAP-CAD-T2-PARAM-SERIES`` carries a clause "gate-power 突变测试沿用
>=80 %".  The numbers behind that clause used to live in prose: a hand-written
table in ``docs/series60_param_series_completion.md`` claimed 31 mutants, five
executed and an 80.0 % rate, while a second page claimed 14/14 -- and the seven
entry/family-comparison mutants (``S60-15`` .. ``S60-21``) had no recorded run at
all.  Nothing here is typed by hand: every count, rate, failing-test name and
wall time is read out of

    artifacts/series60_gate_power/gate_power_series60.json

which is produced only by ``tools/hull_gate_power/mutate.py`` (it patches one
defect into the sources, runs the guarding tests, restores the original bytes and
records the measured exit status).

Outputs::

    docs/series60_gate_power_record.md          full page, tables + verdict
    artifacts/series60_gate_power/gate_power_summary.json   machine-readable
    docs/series60_param_series_completion.md    the "Gate-Power 突变测试" section

A survivor is printed as a survivor: the page states the rate and marks the entry
``SURVIVED``.  An incomplete record (fewer mutants than the harness defines) is
reported as incomplete, never as a rate over whatever happens to be present.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RECORD = ROOT / "artifacts" / "series60_gate_power" / "gate_power_series60.json"
SUMMARY = ROOT / "artifacts" / "series60_gate_power" / "gate_power_summary.json"
PAGE = ROOT / "docs" / "series60_gate_power_record.md"
COMPLETION = ROOT / "docs" / "series60_param_series_completion.md"
HARNESS = ROOT / "tools" / "hull_gate_power" / "mutate.py"
SECTION_HEAD = "### Gate-Power 突变测试"
SECTION_END = "## 已知局限"
#: The anchor's bar.  Not moved after the fact.
BAR = 0.80


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def harness_ids() -> list[str]:
    spec = importlib.util.spec_from_file_location("hull_gate_power_mutate", HARNESS)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return [m["id"] for m in module.MUTANTS if m["id"].startswith("S60-")]


def summarise(record: dict, defined: list[str]) -> dict:
    rows, survivors = [], []
    for mutant_id in sorted(record):
        if mutant_id.startswith("baseline"):
            continue  # a baseline is the unmutated control, not a mutant
        entry = record[mutant_id]
        run = entry.get("run", {})
        mutant = entry.get("mutant", {})
        killed = bool(run.get("killed"))
        rows.append({
            "id": mutant_id,
            "claim": mutant.get("claim", ""),
            "file": mutant.get("file", ""),
            "tests": mutant.get("tests", []),
            "killed": killed,
            "returncode": run.get("returncode"),
            "summary": run.get("summary", ""),
            "failed": run.get("failed", []),
            "wall_seconds": run.get("wall_seconds"),
            "status": entry.get("status"),
        })
        if not killed:
            survivors.append(mutant_id)
    executed = [row for row in rows if row["status"] == "ok"]
    killed = sum(1 for row in executed if row["killed"])
    total = len(executed)
    mutated = {key for key in record if not key.startswith("baseline")}
    missing = sorted(set(defined) - mutated)
    rate = killed / total if total else None
    return {
        "record": str(RECORD.relative_to(ROOT)),
        "record_sha256": sha256_of(RECORD) if RECORD.exists() else None,
        "harness_sha256": sha256_of(HARNESS),
        "chain": "series 60 (S60-*)",
        "defined_by_harness": len(defined),
        "recorded": len(mutated),
        "executed": total,
        "killed": killed,
        "survived": len(survivors),
        "kill_rate": rate,
        "bar": BAR,
        "meets_bar": bool(rate is not None and rate >= BAR),
        "complete": not missing,
        "missing": missing,
        "survivors": survivors,
        "rows": rows,
        "wall_seconds_total": round(sum(row["wall_seconds"] or 0.0 for row in executed), 1),
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "environment": {"python": platform.python_version(), "platform": platform.platform()},
    }


def page(summary: dict) -> str:
    rate = "n/a" if summary["kill_rate"] is None else f"{100 * summary['kill_rate']:.1f} %"
    verdict = ("PASS" if summary["meets_bar"] and summary["complete"]
               else "FAIL" if not summary["meets_bar"] else "INCOMPLETE")
    lines = [
        "# Series 60 gate power — the executed record",
        "",
        f"Generated from `{summary['record']}` by "
        "`tools/series60_family/gate_power_record_page.py`; no number below is hand-typed.",
        "",
        f"* verdict: **{verdict}** (bar {100 * BAR:.0f} %, measured rate {rate})",
        f"* mutants defined by the harness: {summary['defined_by_harness']}",
        f"* mutants recorded: {summary['recorded']} (missing: "
        f"{summary['missing'] if summary['missing'] else 'none'})",
        f"* killed / survived: {summary['killed']} / {summary['survived']}",
        f"* wall time in the mutated runs: {summary['wall_seconds_total']} s",
        f"* record sha256: `{summary['record_sha256']}`",
        f"* harness sha256: `{summary['harness_sha256']}`",
        f"* generated {summary['generated_utc']} on python {summary['environment']['python']}",
        "",
        "## What a kill means here",
        "",
        "One single-defect mutant at a time is injected into `tensorcad/series60.py` or "
        "`tensorcad/series60_cli.py` by exact-string replacement, only the test files that "
        "guard the mutated claim are run, and the file is restored from its original bytes "
        "with a SHA-256 check before the next mutant starts. A mutant survives when the "
        "guarding tests still pass, which is a gap in the tests, not in the hull.",
        "",
        "| mutant | claim it breaks | guarding tests | verdict | failing tests | wall [s] |",
        "|---|---|---|---|---|---|",
    ]
    for row in summary["rows"]:
        tests = ", ".join(Path(t).name for t in row["tests"])
        failed = ", ".join(f.split("::")[-1] for f in row["failed"]) or "—"
        mark = "KILLED" if row["killed"] else "**SURVIVED**"
        lines.append(
            f"| `{row['id']}` | {row['claim']} | {tests} | {mark} | {failed} "
            f"| {row['wall_seconds']:.1f} |"
        )
    lines += ["", "## Survivors", ""]
    if summary["survivors"]:
        lines.append("The following mutants were **not** killed by the current tests. They "
                     "name a test gap and are listed rather than smoothed over:")
        lines.append("")
        for mutant_id in summary["survivors"]:
            row = next(r for r in summary["rows"] if r["id"] == mutant_id)
            lines.append(f"* `{mutant_id}` — {row['claim']} (summary: {row['summary']})")
    else:
        lines.append("None: every mutant in the record was killed by at least one named test.")
    lines += [
        "",
        "## Reproduce",
        "",
        "```",
        "python3 tools/hull_gate_power/mutate.py \\",
        "  --out artifacts/series60_gate_power/gate_power_series60.json --only "
        + ",".join(sorted(r["id"] for r in summary["rows"])),
        "python3 tools/series60_family/gate_power_record_page.py",
        "python3 tools/series60_family/anchor_acceptance.py",
        "python3 -m pytest tests/test_series60_gate_power_record.py -q",
        "```",
        "",
        "## Scope and limits",
        "",
        f"* This record covers the Series-60 chain only (`S60-*`, {summary['executed']} "
        "mutants). The legacy hull/analytic mutants (`HUL-*`, `ANA-*`) are audited "
        "separately in `docs/hull_gate_power.md`; mixing the two chains in one rate is "
        "what made the earlier 80.0 % figure unreadable.",
        "* A kill proves the tests notice the defect, not that the closed form is right: "
        "the external shape anchor (Todd's tabulated offsets) is still WITHHELD — see "
        "`docs/series60_anchor_acceptance.md`.",
        "* Mutation runs edit the working tree. They are only safe with a clean "
        "`git diff -- tensorcad`; the harness refuses to start otherwise and aborts on a "
        "leftover injection marker.",
        "",
    ]
    return "\n".join(lines)


def section(summary: dict) -> str:
    rate = "n/a" if summary["kill_rate"] is None else f"{100 * summary['kill_rate']:.1f}%"
    verdict = ("pass" if summary["meets_bar"] and summary["complete"]
               else "fail" if not summary["meets_bar"] else "incomplete")
    out = [
        SECTION_HEAD,
        "",
        "本节的数字由 `tools/series60_family/gate_power_record_page.py` 从 "
        f"`{summary['record']}` 程序化注入，不手抄。",
        "",
        f"- **家族突变体总数**（harness 定义）: {summary['defined_by_harness']}（`S60-*`）",
        f"- **记录中的突变体**: {summary['recorded']}（缺: "
        f"{summary['missing'] if summary['missing'] else '无'}）",
        f"- **已执行/杀死/幸存**: {summary['executed']} / {summary['killed']} / "
        f"{summary['survived']}",
        f"- **杀死率**: {rate}（验收线 {100 * BAR:.0f}%）→ **{verdict}**",
        f"- **记录 sha256**: `{summary['record_sha256']}`",
        f"- **完整表格**: `docs/series60_gate_power_record.md`",
        "",
    ]
    if summary["survivors"]:
        out.append("> 幸存突变体（测试缺口，如实列出）: "
                   + ", ".join(f"`{m}`" for m in summary["survivors"]))
        out.append("")
    else:
        out.append("> 无幸存突变体：每个突变都被具名测试杀死。")
        out.append("")
    return "\n".join(out)


def patch_completion(summary: dict) -> bool:
    text = COMPLETION.read_text()
    start = text.find(SECTION_HEAD)
    end = text.find(SECTION_END)
    if start < 0 or end < 0 or end < start:
        return False
    COMPLETION.write_text(text[:start] + section(summary) + "\n" + text[end:])
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--record", default=str(RECORD))
    args = parser.parse_args()
    record_path = Path(args.record)
    if not record_path.is_absolute():
        record_path = ROOT / record_path
    if not record_path.exists():
        print(f"[abort] no mutation record at {record_path}; run the harness first "
              "(tools/hull_gate_power/mutate.py). No page is written from an absent run.")
        return 2
    record = json.loads(record_path.read_text())
    summary = summarise(record, harness_ids())
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n")
    PAGE.write_text(page(summary))
    patched = patch_completion(summary)
    print(f"[written] {PAGE.relative_to(ROOT)}")
    print(f"[written] {SUMMARY.relative_to(ROOT)}")
    print(f"[patched] {COMPLETION.relative_to(ROOT)}: {patched}")
    print(f"[killed]  {summary['killed']}/{summary['executed']} "
          f"({100 * (summary['kill_rate'] or 0):.1f}%) complete={summary['complete']} "
          f"meets_bar={summary['meets_bar']}")
    return 0 if summary["meets_bar"] and summary["complete"] else 1


if __name__ == "__main__":
    raise SystemExit(main())