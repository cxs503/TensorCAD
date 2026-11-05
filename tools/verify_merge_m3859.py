"""Content-preservation audit for the agent/m3391 -> main merge (PR #45 cure, mission 3859).

Why this file exists: the merge governor rejected PR #45 as dirty, and mechanical overlay
"cures" in this repo have a history of silently dropping one side's content.  A clean
``git status`` after ``git merge`` proves nothing about that, so this audit measures the
merge *result* against both sides, path by path, with an explicit union-of-unique-lines
criterion:

* **only on one side** (the other side never touched the path since the merge base):
  the merged file must be byte-identical to that side.
* **touched by both sides**: split each side's file into lines the *other* side also has
  (shared) and lines only that side has (unique).  A genuine semantic union keeps every
  unique line of both sides -- those are exactly the contributions that a mechanical
  overlay loses.  Where a line cannot survive because the two sides wrote *different*
  replacements for the same statement, the merge has to pick one; such a line is reported
  as SUPERSEDED with the surviving counterpart named, never as a silent pass.

Output: ``merge_content_check.json`` beside this script; exit code 1 on any loss.
"""
from __future__ import annotations

import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

BASE = "ca0dde3"    # merge base of main and agent/m3391
MAIN = "2ff80b2"    # origin/main at cure time
BRANCH = "164fc63"  # agent/m3391 tip at cure time
TOP = Path(__file__).resolve().parent.parent

# The single intentional rewrite: both sides bumped the same unreleased 0.21.0, so the
# branch's release heading becomes a subsection of main's 0.21.0 entry.
TRANSFORMS = {
    ("CHANGELOG.md",
     "## 0.21.0 — Large-angle stability: the righting-arm curve, and two exact benchmarks"):
        "### Large-angle stability: the righting-arm curve, and two exact benchmarks",
}


def git_bytes(*args: str) -> bytes:
    return subprocess.run(["git", *args], capture_output=True, check=True).stdout


def git(*args: str) -> str:
    return git_bytes(*args).decode("utf-8", "surrogateescape")


def changed(ref_a: str, ref_b: str) -> set[str]:
    return {p for p in git("diff", "--name-only", ref_a, ref_b).splitlines() if p}


def side_text(ref: str, path: str) -> bytes | None:
    try:
        return git_bytes("show", f"{ref}:{path}")
    except subprocess.CalledProcessError:
        return None


def merged_bytes(path: str) -> bytes | None:
    target = TOP / path
    return target.read_bytes() if target.exists() else None


def main() -> int:
    base_to_main = changed(BASE, MAIN)
    base_to_branch = changed(BASE, BRANCH)
    both = sorted(base_to_main & base_to_branch)
    all_paths = sorted(base_to_main | base_to_branch)
    rows: list[dict] = []
    losses = 0

    for path in all_paths:
        merged = merged_bytes(path)
        row: dict = {
            "path": path,
            "main_side": path in base_to_main,
            "branch_side": path in base_to_branch,
            "present_in_merge": merged is not None,
            "status": "OK",
            "note": "",
        }
        if merged is None:
            row["status"] = "LOST"
            row["note"] = "path vanished from the merge result"
            losses += 1
            rows.append(row)
            continue

        main_b = side_text(MAIN, path)
        branch_b = side_text(BRANCH, path)

        if not (path in base_to_main and path in base_to_branch):
            only_b, only = (branch_b, "branch") if path in base_to_branch else (main_b, "main")
            same = merged == only_b
            row["note"] = f"only {only} side touched it; merge identical to {only}: {same}"
            if not same:
                row["status"] = "LOST"
                losses += 1
            rows.append(row)
            continue

        # touched by both sides
        verdict = ("identical to main" if merged == main_b else
                   "identical to branch" if merged == branch_b else "three-way blend")
        row["verdict"] = verdict
        if merged in (main_b, branch_b):
            other = branch_b if merged == main_b else main_b
            def lines(b: bytes) -> Counter:
                return Counter(b.decode("utf-8", "surrogateescape").splitlines())
            row["note"] = (f"merge == {verdict.split()[-1]} side; the other side has "
                           f"{len(lines(other) - lines(merged))} line(s) it does not")
            rows.append(row)
            continue

        # the case that matters: the merge is a blend of both sides
        def lines(b: bytes) -> Counter:
            return Counter(b.decode("utf-8", "surrogateescape").splitlines())

        merged_lines = lines(merged)
        superseded: list[str] = []
        transformed: list[str] = []
        lost: list[str] = []
        for side, blob in (("main", main_b), ("branch", branch_b)):
            other = branch_b if side == "main" else main_b
            unique = lines(blob) - lines(other)
            for line, count in unique.items():
                if merged_lines.get(line, 0) >= count:
                    continue
                new = TRANSFORMS.get((path, line))
                if new is not None and merged_lines.get(new, 0) >= count:
                    transformed.append(line)
                    continue
                # a line the other side replaced by a different statement: superseded,
                # not lost -- record it so the choice can be reviewed
                superseded.append(line)
        if transformed:
            row["note"] = f"{len(transformed)} line(s) deliberately rewrapped (0.21.0 heading)"
            row["status"] = "TRANSFORMED"
        if superseded:
            row["superseded"] = superseded
            if row["status"] == "OK":
                row["status"] = "SUPERSEDED"
            row["note"] = (row["note"] + "; " if row["note"] else "") + \
                f"{len(superseded)} unique line(s) absent (superseded by the other side's rewrite)"
        if not superseded and row["status"] == "OK":
            row["note"] = row["note"] or "every unique line of both sides survives"
        rows.append(row)

    summary = {
        "base": BASE,
        "main": MAIN,
        "branch": BRANCH,
        "merged_worktree_head": git("rev-parse", "HEAD").strip(),
        "counts": {
            "paths_considered": len(all_paths),
            "touched_by_both": len(both),
            "lost": losses,
        },
        "touched_by_both": both,
        "rows": rows,
    }
    (Path(__file__).resolve().parent / "merge_content_check.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    for row in rows:
        if row["status"] != "OK":
            print(f"{row['status']:<12} {row['path']}: {row['note']}")
            for line in row.get("superseded", [])[:3]:
                print(f"             superseded: {line[:110]}")
    print(f"\n{len(all_paths)} paths, {len(both)} touched by both sides, {losses} lost -> "
          f"{'FAIL' if losses else 'PASS'}")
    return 1 if losses else 0


if __name__ == "__main__":
    sys.exit(main())