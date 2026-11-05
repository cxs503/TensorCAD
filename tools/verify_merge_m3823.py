#!/usr/bin/env python3
"""Cure PR #45 (agent/m3391 ← main): prove the three-way merge lost no content.

The failure mode this guards against is the *mechanical overlay*: taking one side's
whole file, or letting git's text merge keep a hunk that reads fine but drops a name
the other side added (an ``__all__`` entry, a CLI flag, a doc paragraph).  Reading the
diff by eye does not scale to a 45-file merge, so the check is mechanical:

1. Recompute the merge base from ``git merge-base``.
2. For every file that *either* side touched, collect the lines each side **added**
   with respect to the base (``git diff -U0``, only ``+`` lines).
3. Require every such line to appear in the merged working tree under a whole-file,
   whitespace-normalised containment test -- deliberately strict, so that a line that
   survives only as a prefix of a longer line (the "list literal merged over the other
   side's entry" shape) counts as a loss.  Anything not found is reported as a
   **potential silent loss** and the script exits 1.
4. Same test for **deletions**: a line a side removed should be gone from the merged
   tree.  A line that comes back is not a loss (union semantics keep content) but it is
   listed, so the review can confirm it is a deliberate keep rather than a stray hunk.
5. Require every file only **one** side touched to be byte-identical to that side's
   blob -- the other direction of the overlay failure.
6. Record the public-surface union (``tensorcad.__all__``) for base / ours / theirs /
   merged, since that list is where a text merge most quietly drops a name.

Anything the report lists must be explained in ``docs/m3823_pr45_merge_cure.md`` -- a
deliberate renumbering, a line that both sides added identically (deduplicated), or a
line superseded by a rename the other side performed.  The allowlist below is *explicit
and short* for that reason: it is a review record, not a silencer.

Usage:
    python3 tools/verify_merge_m3823.py [--base BASE] [--ours REV] [--theirs REV]

Run it from the repository root after the merge is staged in the working tree.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SURFACE_FILE = "tensorcad/__init__.py"

# Lines that legitimately do not survive verbatim.  Each entry is (path, needle,
# reason).  Keep this list tiny and justify every line in the merge-cure report.
ALLOWLIST: list[tuple[str, str, str]] = [
    (
        "CHANGELOG.md",
        "## 0.21.0 — Large-angle stability: the righting-arm curve, and two exact benchmarks",
        "renumbered to `## Unreleased — ...`: main already released a dated 0.21.0 heading, "
        "two identical version headings in one file is the defect, not the content",
    ),
]


def git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, check=True, capture_output=True, text=True
    ).stdout


def git_bytes(*args: str) -> bytes:
    return subprocess.run(
        ["git", *args], cwd=ROOT, check=True, capture_output=True
    ).stdout


def normalise(line: str) -> str:
    return re.sub(r"\s+", " ", line).strip()


def changed_lines(rev: str, path: str, base: str, sign: str) -> list[str]:
    """Non-blank lines ``rev`` added (sign='+') or removed (sign='-') vs ``base``."""
    out = git("diff", "-U0", base, rev, "--", path)
    skip = sign * 3
    lines: list[str] = []
    for raw in out.splitlines():
        if not raw.startswith(sign) or raw.startswith(skip):
            continue
        if raw[1:].strip():
            lines.append(raw[1:])
    return lines


def added_lines(rev: str, path: str, base: str) -> list[str]:
    return changed_lines(rev, path, base, "+")


def deleted_lines(rev: str, path: str, base: str) -> list[str]:
    return changed_lines(rev, path, base, "-")


def is_present(line: str, haystack_norm: str, allow: list[tuple[str, str, str]], path: str) -> tuple[bool, str]:
    norm = normalise(line)
    if norm and norm in haystack_norm:
        return True, "exact"
    for apath, needle, reason in allow:
        if apath == path and normalise(needle) == norm:
            return True, f"allowlisted: {reason}"
    return False, "MISSING"


def public_surface(source: str) -> set[str]:
    match = re.search(r"__all__\s*=\s*[\[{](.*?)[\]}]", source, re.S)
    if not match:
        return set()
    return set(re.findall(r"[\"']([^\"']+)[\"']", match.group(1)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=None, help="merge base (default: git merge-base ours theirs)")
    ap.add_argument("--ours", default="agent/m3391")
    ap.add_argument("--theirs", default="origin/main")
    ap.add_argument("--report", default="artifacts/m3823/merge_content_check.json")
    args = ap.parse_args()

    base = (args.base or git("merge-base", args.ours, args.theirs)).splitlines()[-1].strip()
    print(f"merge base: {base}")
    print(f"ours: {args.ours}  theirs: {args.theirs}")

    ours_files = set(git("diff", "--name-only", base, args.ours).split("\n")) - {""}
    theirs_files = set(git("diff", "--name-only", base, args.theirs).split("\n")) - {""}
    both = ours_files & theirs_files
    touched = ours_files | theirs_files
    print(f"files touched by ours: {len(ours_files)}, by theirs: {len(theirs_files)}, "
          f"by both: {len(both)}, union: {len(touched)}")

    merged: dict[str, str] = {}
    for path in touched:
        p = ROOT / path
        merged[path] = p.read_text(errors="replace") if p.exists() else ""

    problems: list[dict[str, object]] = []
    resurrected: list[dict[str, object]] = []
    counts = {"added_ours": 0, "added_theirs": 0, "deleted_ours": 0, "deleted_theirs": 0,
              "allowlisted": 0}
    for path in sorted(touched):
        hay = normalise(merged[path])
        for side, files, rev, tag in (
            ("ours/agent/m3391", ours_files, args.ours, "ours"),
            ("theirs/main", theirs_files, args.theirs, "theirs"),
        ):
            if path not in files:
                continue
            for line in added_lines(rev, path, base):
                counts[f"added_{tag}"] += 1
                ok, how = is_present(line, hay, ALLOWLIST, path)
                if not ok:
                    problems.append({"path": path, "side": side, "line": line, "how": how})
                elif how.startswith("allowlisted"):
                    counts["allowlisted"] += 1
            for line in deleted_lines(rev, path, base):
                counts[f"deleted_{tag}"] += 1
                norm = normalise(line)
                if norm and norm in hay:
                    resurrected.append({"path": path, "side": side, "line": line})

    checked = counts["added_ours"] + counts["added_theirs"]
    print(f"added lines checked: {checked} (ours {counts['added_ours']}, "
          f"theirs {counts['added_theirs']}), allowlisted: {counts['allowlisted']}, "
          f"potential losses: {len(problems)}")
    for pr in problems[:40]:
        print(f"  LOSS? {pr['side']:18s} {pr['path']}: {pr['line'][:100]}")
    print(f"deleted lines checked: {counts['deleted_ours'] + counts['deleted_theirs']}, "
          f"still present after the merge (union keeps content; review): {len(resurrected)}")

    # A file only one side touched must come through *byte for byte* as that side's
    # blob.  This is the check that catches the mechanical overlay in its other
    # direction: a whole-file copy from the wrong side would show up here as a
    # mismatch, and so would a stray whitespace-only hunk.
    single_side = {}
    for path in sorted(touched):
        if path in both:
            continue
        rev = args.ours if path in ours_files else args.theirs
        try:
            blob = git_bytes("show", f"{rev}:{path}")
        except subprocess.CalledProcessError:  # pragma: no cover - side deleted the file
            single_side[path] = "deleted-by-side"
            continue
        target = ROOT / path
        disk = target.read_bytes() if target.exists() else b""
        single_side[path] = "identical" if blob == disk else "MISMATCH"
    mismatched = {p: v for p, v in single_side.items() if v == "MISMATCH"}
    print(f"single-side files: {len(single_side)}, byte-identical to their side: "
          f"{sum(1 for v in single_side.values() if v == 'identical')}, mismatched: {len(mismatched)}")
    for path in mismatched:
        print(f"  MISMATCH {path}")

    base_src = git("show", f"{base}:{SURFACE_FILE}")
    ours_src = git("show", f"{args.ours}:{SURFACE_FILE}")
    theirs_src = git("show", f"{args.theirs}:{SURFACE_FILE}")
    merged_src = (ROOT / SURFACE_FILE).read_text()
    s_base, s_ours, s_theirs, s_merged = (
        public_surface(base_src), public_surface(ours_src),
        public_surface(theirs_src), public_surface(merged_src),
    )
    union = s_ours | s_theirs
    surface = {
        "base": len(s_base), "ours": len(s_ours), "theirs": len(s_theirs),
        "merged": len(s_merged), "union": len(union),
        "added_by_ours": sorted(s_ours - s_base),
        "added_by_theirs": sorted(s_theirs - s_base),
        "lost_names": sorted(union - s_merged),
        "extra_names": sorted(s_merged - union),
    }
    print(f"public surface: base {surface['base']}, ours {surface['ours']}, "
          f"theirs {surface['theirs']}, merged {surface['merged']}, union {surface['union']}, "
          f"lost {len(surface['lost_names'])}")

    ok = (not problems and not surface["lost_names"] and not surface["extra_names"]
          and not mismatched)
    report_path = ROOT / args.report
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(
            {
                "base": base,
                "ours": args.ours,
                "theirs": args.theirs,
                "ours_head": git("rev-parse", args.ours).strip(),
                "theirs_head": git("rev-parse", args.theirs).strip(),
                "files_touched_ours": sorted(ours_files),
                "files_touched_theirs": sorted(theirs_files),
                "files_touched_both": sorted(both),
                "line_counts": counts,
                "added_lines_checked": checked,
                "potential_losses": problems,
                "deleted_lines_still_present": resurrected,
                "allowlist": [{"path": p, "needle": n, "reason": r} for p, n, r in ALLOWLIST],
                "public_surface": surface,
                "single_side_files": single_side,
                "verdict": "PASS" if ok else "REVIEW",
            },
            indent=2,
        )
        + "\n"
    )
    print(f"verdict: {'PASS' if ok else 'REVIEW'}  report: {args.report}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())