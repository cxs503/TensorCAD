#!/usr/bin/env python3
"""Anti-silent-loss audit for the cure of PR #37 (agent/m2873 x main).

The merge that cures a dirty PR must be a *semantic union*: nothing either side
added may disappear, and nothing either side deleted for a reason may come back
by accident.  Mechanical overlays have silently dropped content in this repo
before, so the property is checked instead of asserted.

Two independent checks:

``files``
    every path added or modified by either side (relative to the merge base)
    exists in the merged tree, unless a side deleted it on purpose;
``lines``
    every non-trivial line either side added (relative to the merge base) is
    present in the merged file.  Lines are normalised (stripped, blank and
    conflict markers/version digits dropped) so pure renumbering -- e.g. the
    changelog header moving from ``0.14.0`` to the merged ``0.21.0`` -- does not
    raise a false positive; every such exemption is listed explicitly in
    ``EXPECTED_RELABELS`` and is fatal if it stops matching.

Usage (inside the worktree, after ``git add`` of the resolution)::

    python3 docs/merge_pr37/verify_union.py --base <sha> --ours <sha> --theirs <sha> --out report.json
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

# Lines that legitimately differ between an input side and the merged tree:
# (path, side, exact line in that side) -> the line that replaces it in the merge.
EXPECTED_RELABELS = {
    ("CHANGELOG.md", "ours", "## 0.14.0 \u2014 2026-11-03"): "## 0.21.0 \u2014 2026-11-04",
    ("pyproject.toml", "ours", 'version = "0.14.0"'): 'version = "0.21.0"',
    ("pyproject.toml", "theirs", 'version = "0.20.0"'): 'version = "0.21.0"',
    ("tensorcad/__init__.py", "ours", '__version__ = "0.14.0"'): '__version__ = "0.21.0"',
    ("tensorcad/__init__.py", "theirs", '__version__ = "0.20.0"'): '__version__ = "0.21.0"',
    (
        "tests/test_public_surface.py",
        "ours",
        'assert tensorcad.__version__ == "0.14.0"',
    ): 'assert tensorcad.__version__ == "0.21.0"',
    (
        "tests/test_public_surface.py",
        "theirs",
        'assert tensorcad.__version__ == "0.20.0"',
    ): 'assert tensorcad.__version__ == "0.21.0"',
    ("tests/test_public_surface.py", "ours", "# transom-stern family (0.14.0)"): None,
    # the module tuple was folded onto one line so both new modules stay listed
    (
        "tests/test_public_surface.py",
        "ours",
        '"transom",',
    ): '"wigley", "transom", "hull", "design", "suboff", "suboff_appendages", "cli",',
    (
        "tests/test_public_surface.py",
        "theirs",
        '"wigley", "hull", "design", "suboff", "suboff_appendages", "cli",',
    ): '"wigley", "transom", "hull", "design", "suboff", "suboff_appendages", "cli",',
}

TRIVIAL = re.compile(r"^[\s{}\[\]()]*$")


def git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], check=True, capture_output=True, text=True
    ).stdout


def names(ref: str, base: str) -> dict[str, str]:
    """path -> status letter, for paths touched between base and ref."""
    out = git("diff", "--name-status", "-z", f"{base}", f"{ref}", "--")
    fields = [f for f in out.split("\0") if f]
    result: dict[str, str] = {}
    it = iter(fields)
    for status in it:
        path = next(it)
        result[path] = status[0]
    return result


def added_lines(base: str, ref: str, path: str) -> list[str]:
    diff = git("diff", "-U0", f"{base}", f"{ref}", "--", path)
    lines = []
    for line in diff.splitlines():
        if line.startswith("+") and not line.startswith("+++"):
            text = line[1:].strip()
            if text and not TRIVIAL.match(text) and not text.startswith(("<<<<<<<", "=======", ">>>>>>>")):
                lines.append(text)
    return lines


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--ours", required=True, help="mission branch tip (agent/m2873)")
    ap.add_argument("--theirs", required=True, help="main tip merged in")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    sides = {"ours": args.ours, "theirs": args.theirs}
    touched = {side: names(ref, args.base) for side, ref in sides.items()}
    merged = set(Path(".").rglob("*"))
    report: dict = {"base": args.base, "ours": args.ours, "theirs": args.theirs}

    # -- file set ---------------------------------------------------------
    missing: list[str] = []
    intended_deletions: dict[str, list[str]] = {"ours": [], "theirs": []}
    for side, table in touched.items():
        for path, status in table.items():
            exists = Path(path).exists()
            if status == "D":
                if exists:
                    intended_deletions[side].append(path)  # informational only
                continue
            if not exists:
                # dropped unless the other side deleted it on purpose
                other = "theirs" if side == "ours" else "ours"
                if touched[other].get(path) == "D":
                    intended_deletions[other].append(path)
                else:
                    missing.append(f"{side}::{path}")

    # -- line set ---------------------------------------------------------
    merged_cache: dict[str, set[str]] = {}
    dropped: list[dict[str, str]] = []
    relabelled: list[dict[str, str]] = []
    files = sorted(set().union(*[set(t) for t in touched.values()]))
    for path in files:
        if not Path(path).exists() or Path(path).is_dir():
            continue
        if path not in merged_cache:
            merged_cache[path] = {
                ln.strip() for ln in Path(path).read_text(errors="replace").splitlines()
            }
        present = merged_cache[path]
        for side, ref in sides.items():
            if touched[side].get(path) in (None, "D"):
                continue
            for line in added_lines(args.base, ref, path):
                if line in present:
                    continue
                replacement = EXPECTED_RELABELS.get((path, side, line), "sentinel")
                if replacement == "sentinel":
                    dropped.append({"side": side, "path": path, "line": line})
                elif replacement is not None and replacement not in present:
                    dropped.append(
                        {"side": side, "path": path, "line": line,
                         "note": f"expected relabel -> {replacement!r} also absent"}
                    )
                else:
                    relabelled.append({"side": side, "path": path, "line": line,
                                       "merged_as": replacement})

    # -- registry / tuple sanity -----------------------------------------
    import tensorcad  # noqa: PLC0415  (import after the text audit, on purpose)

    dup_all = [n for n, c in Counter(tensorcad.__all__).items() if c > 1]

    report.update(
        files_checked=len(files),
        files_missing=missing,
        intended_deletions=intended_deletions,
        lines_dropped=dropped,
        lines_relabelled=relabelled,
        public_surface={
            "version": tensorcad.__version__,
            "len___all__": len(tensorcad.__all__),
            "duplicates": dup_all,
        },
        ok=not missing and not dropped and not dup_all,
    )
    text = json.dumps(report, indent=2, ensure_ascii=False)
    if args.out:
        Path(args.out).write_text(text + "\n")
    print(text)
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())