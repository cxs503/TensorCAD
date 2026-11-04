#!/usr/bin/env python3
"""Audit the semantic union of PR #37 (agent/m2873) into main.

For every path either side touched since the merge base, compare the merged
working tree against that side's blob.  Nothing introduced by one side may be
absent from the result unless the OTHER side also touched it (a genuine
three-way conflict, resolved by hand -> listed explicitly).

Usage: python3 docs/m3736/verify_union.py [--json out.json]
Exit code 1 on any silent loss.
"""
from __future__ import annotations

import json
import subprocess
import sys

BASE = "6c568c3"            # merge base of origin/main and origin/agent/m2873
OURS = "origin/main"
THEIRS = "origin/agent/m2873"


def git(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True, check=True).stdout


def tree(rev: str) -> dict[str, str]:
    """path -> blob sha for a revision."""
    out = {}
    for line in git("ls-tree", "-r", rev).splitlines():
        meta, path = line.split("\t", 1)
        _mode, _type, sha = meta.split()
        out[path] = sha
    return out


def changed(rev: str) -> dict[str, str]:
    """path -> status for base..rev."""
    res = {}
    for line in git("diff", "--name-status", BASE, rev).splitlines():
        parts = line.split("\t")
        res[parts[-1]] = parts[0][0]
    return res


def worktree_blob(path: str) -> str | None:
    try:
        return git("hash-object", path).strip()
    except subprocess.CalledProcessError:
        return None


def main() -> int:
    ours, theirs = changed(OURS), changed(THEIRS)
    both = sorted(set(ours) & set(theirs))
    only_ours = sorted(set(ours) - set(theirs))
    only_theirs = sorted(set(theirs) - set(ours))

    report = {"base": BASE, "ours": OURS, "theirs": THEIRS,
              "counts": {"ours_only": len(only_ours), "theirs_only": len(only_theirs),
                         "both": len(both)},
              "conflict_files": [], "ours_only": {}, "theirs_only": {},
              "silent_loss": [], "ok": True}

    their_tree = tree(THEIRS)
    # their-side paths must be present byte-identically, unless we also touched it
    for path in only_theirs + both:
        want = their_tree.get(path)
        if theirs[path] == "D":
            got = worktree_blob(path)
            if got is not None:
                report["silent_loss"].append({"path": path, "side": "theirs",
                                              "why": "deletion not applied"})
            continue
        got = worktree_blob(path)
        same = got == want
        if path in both:
            report["conflict_files"].append({"path": path, "theirs_kept_verbatim": same,
                                             "ours_status": ours[path],
                                             "theirs_status": theirs[path]})
        elif not same:
            report["silent_loss"].append({"path": path, "side": "theirs",
                                          "why": "content differs from branch blob"})
        else:
            report["theirs_only"][path] = "kept verbatim"

    our_tree = tree(OURS)
    for path in only_ours:
        if ours[path] == "D":
            if worktree_blob(path) is not None:
                report["silent_loss"].append({"path": path, "side": "ours",
                                              "why": "deletion not applied"})
            continue
        got = worktree_blob(path)
        if got != our_tree.get(path):
            report["silent_loss"].append({"path": path, "side": "ours",
                                          "why": "content differs from main blob"})
        else:
            report["ours_only"][path] = "kept verbatim"

    report["ok"] = not report["silent_loss"]
    text = json.dumps(report, indent=2, ensure_ascii=False)
    if "--json" in sys.argv:
        open(sys.argv[sys.argv.index("--json") + 1], "w", encoding="utf-8").write(text + "\n")
    print(text)
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
