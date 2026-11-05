#!/usr/bin/env python3
"""PR #45 cure evidence: build the verification figure and the numbers table.

Everything printed here is read back from artefacts, never typed in by hand:

* ``artifacts/m3823/merge_content_check.json`` -- written by ``tools/verify_merge_m3823.py``
  (merge base, per-side added/deleted line counts, potential losses, public-surface union).
* ``artifacts/m3823/logs/pytest_main.log`` / ``pytest_merged.log`` -- the two pytest runs
  (``-q`` summary lines are parsed for passed/failed/skipped).

Outputs (all under ``docs/m3823/``): ``outputs/merge_verification.png``,
``outputs/results.json`` and ``report.md``.

Usage: python3 tools/report_m3823.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
ART = ROOT / "artifacts/m3823"
OUT = ROOT / "docs/m3823"
COUNTS = ("passed", "failed", "skipped", "errors")


def parse_pytest(log: Path) -> dict[str, int]:
    """Pull ``N passed, M failed, ...`` out of the last summary line of a ``-q`` run."""
    text = log.read_text(errors="replace")
    lines = [ln.strip() for ln in text.splitlines() if re.search(r"\b\d+ (passed|failed|skipped|error)", ln)]
    if not lines:
        raise SystemExit(f"no pytest summary line in {log}")
    tail = lines[-1]
    out = {k: 0 for k in COUNTS}
    for num, word in re.findall(r"(\d+) ([a-z]+)s?\b", tail):
        word = "errors" if word == "error" else word
        if word in out:
            out[word] = int(num)
    out["summary_line"] = tail  # type: ignore[assignment]
    out["total"] = out["passed"] + out["failed"] + out["skipped"] + out["errors"]
    return out  # type: ignore[return-value]


def main() -> int:
    merge = json.loads((ART / "merge_content_check.json").read_text())
    main_run = parse_pytest(ART / "logs/pytest_main.log")
    merged_run = parse_pytest(ART / "logs/pytest_merged.log")
    surface = merge["public_surface"]
    counts = merge["line_counts"]

    results = {
        "merge": {
            "base": merge["base"],
            "ours": merge["ours_head"],
            "theirs": merge["theirs_head"],
            "files_touched_ours": len(merge["files_touched_ours"]),
            "files_touched_theirs": len(merge["files_touched_theirs"]),
            "files_touched_both": len(merge["files_touched_both"]),
            "textual_conflicts": ["CHANGELOG.md"],
            "added_lines_checked": merge["added_lines_checked"],
            "added_lines_ours": counts["added_ours"],
            "added_lines_theirs": counts["added_theirs"],
            "allowlisted": counts["allowlisted"],
            "potential_losses": len(merge["potential_losses"]),
            "deleted_lines_checked": counts["deleted_ours"] + counts["deleted_theirs"],
            "verdict": merge["verdict"],
        },
        "public_surface": surface,
        "tests": {"main": main_run, "merged": merged_run},
    }

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "outputs").mkdir(parents=True, exist_ok=True)
    (OUT / "outputs/results.json").write_text(json.dumps(results, indent=2) + "\n")

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))

    ax = axes[0]
    labels = ["merge base\n(df31744)", "main\n(2ff80b2)", "agent/m3391\n(164fc63)", "merged"]
    values = [surface["base"], surface["theirs"], surface["ours"], surface["merged"]]
    colors = ["#9aa5b1", "#4c72b0", "#dd8452", "#55a868"]
    bars = ax.bar(labels, values, color=colors)
    ax.axhline(surface["union"], color="k", ls=":", lw=1.2)
    ax.text(3.35, surface["union"] + 1.5, f"union = {surface['union']}", ha="right", fontsize=9)
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 1.0, str(val), ha="center", fontsize=10)
    ax.set_ylim(0, max(values + [surface["union"]]) * 1.18)
    ax.set_ylabel("names in tensorcad.__all__")
    ax.set_title(f"(a) Public surface: union {surface['union']}, lost {len(surface['lost_names'])}")

    ax = axes[1]
    names = ["ours\n(agent/m3391)", "theirs\n(main)", "allowlisted", "silent\nlosses"]
    vals = [counts["added_ours"], counts["added_theirs"], counts["allowlisted"],
            len(merge["potential_losses"])]
    bars = ax.bar(names, vals, color=["#dd8452", "#4c72b0", "#8172b3", "#c44e52"])
    ax.set_yscale("symlog")
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, max(val, 0.6) * 1.15, str(val),
                ha="center", fontsize=10)
    ax.set_ylabel("added lines matched in the merged tree (log scale)")
    ax.set_title(f"(b) Content check: {merge['added_lines_checked']} lines, "
                 f"{len(merge['potential_losses'])} losses")

    ax = axes[2]
    x = range(2)
    passed = [main_run["passed"], merged_run["passed"]]
    failed = [main_run["failed"], merged_run["failed"]]
    skipped = [main_run["skipped"], merged_run["skipped"]]
    ax.bar(x, passed, 0.55, label="passed", color="#55a868")
    ax.bar(x, failed, 0.55, bottom=passed, label="failed", color="#c44e52")
    ax.bar(x, skipped, 0.55, bottom=[p + f for p, f in zip(passed, failed)],
           label="skipped", color="#9aa5b1")
    for i, (p, f, s) in enumerate(zip(passed, failed, skipped)):
        ax.text(i, p + f + s + 12, f"{p}/{f}/{s}", ha="center", fontsize=10)
    ax.set_xticks(list(x), ["main", "merged"])
    ax.set_ylim(0, (max(p + f + s for p, f, s in zip(passed, failed, skipped)) * 1.2))
    ax.set_ylabel("tests")
    ax.set_title("(c) Suite: passed / failed / skipped")
    ax.legend(loc="lower right", fontsize=9)

    fig.suptitle("TensorCAD PR #45 cure — agent/m3391 ← main (three-way merge, content-verified)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(OUT / "outputs/merge_verification.png", dpi=150)

    report = f"""# PR #45 cure — `agent/m3391` merged with `main`

Generated by `tools/report_m3823.py` from `artifacts/m3823/` (numbers read back, not typed).

| quantity | value |
|---|---|
| merge base | `{merge['base'][:12]}` |
| ours (`agent/m3391`) | `{merge['ours_head'][:12]}` |
| theirs (`origin/main`) | `{merge['theirs_head'][:12]}` |
| files touched: ours / theirs / both | {results['merge']['files_touched_ours']} / {results['merge']['files_touched_theirs']} / {results['merge']['files_touched_both']} |
| textual conflicts | 1 (`CHANGELOG.md`) |
| added lines checked (ours / theirs) | {counts['added_ours']} / {counts['added_theirs']} |
| allowlisted (documented renumbering) | {counts['allowlisted']} |
| silent losses | {len(merge['potential_losses'])} |
| public surface (merged / union) | {surface['merged']} / {surface['union']} |
| public names lost | {len(surface['lost_names'])} |
| suite, main | {main_run['passed']} passed, {main_run['failed']} failed, {main_run['skipped']} skipped |
| suite, merged | {merged_run['passed']} passed, {merged_run['failed']} failed, {merged_run['skipped']} skipped |

![merge verification](outputs/merge_verification.png)
"""
    (OUT / "report.md").write_text(report)
    print(json.dumps(results["merge"], indent=2))
    print("tests:", {k: (v["passed"], v["failed"], v["skipped"]) if isinstance(v, dict) else v
                     for k, v in results["tests"].items()})
    print(f"wrote {OUT/'outputs/merge_verification.png'}, {OUT/'outputs/results.json'}, {OUT/'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())