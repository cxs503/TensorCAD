"""Verify that a merge kept *both* sides' content (guard against silent overlay loss).

Reads the three commits (merge base, ours, theirs, merge result) straight out of git and
compares file-by-file.  Three checks, each of which a mechanical `git checkout --theirs`
style overlay would fail:

1. **Additions survive.** Every path added on either side relative to the merge base
   exists in the merge result, with byte-identical content -- unless the path is one of
   the explicitly declared *conflict* paths, which is checked separately.
2. **One-sided edits survive.** A path modified on exactly one side must be byte-identical
   to that side in the result.
3. **Conflict resolution is a union, not a choice.** For each declared conflict path the
   result must contain (as a multiset of lines) every non-trivial line that either side
   contributed relative to the base, and must carry no unresolved conflict markers.

Also asserts that ``tensorcad.__all__`` in the result is a superset of both sides' and that
no file from either side disappeared.
"""

from __future__ import annotations

import json
import subprocess
import sys
from collections import Counter

REPO = "/root/tp-agent/work/m3832"

# Declared conflict paths and the justification recorded for each one.
CONFLICTS = {
    "CHANGELOG.md": "union: righting-arm bullets (ours) + '### Added' block (theirs) under 0.21.0",
    "tests/test_public_surface.py": "union: module list = {stability, righting_arm, bonjean} + {transom}",
}


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", REPO, *args], capture_output=True, text=True, check=True).stdout


def ls_tree(rev: str) -> dict[str, str]:
    """path -> blob sha for every file in ``rev``."""
    out = {}
    for line in git("ls-tree", "-r", rev).splitlines():
        meta, path = line.split("\t", 1)
        mode, kind, sha = meta.split()
        if kind == "blob":
            out[path] = sha
    return out


def blob(sha: str) -> str:
    return subprocess.run(["git", "-C", REPO, "cat-file", "-p", sha],
                          capture_output=True, text=True, check=True).stdout


def nontrivial_lines(text: str) -> Counter:
    return Counter(ln.strip() for ln in text.splitlines() if ln.strip() and ln.strip() not in ("{", "}", "(", ")"))


def main() -> int:
    base, ours, theirs, merged = sys.argv[1:5]
    trees = {name: ls_tree(rev) for name, rev in
             (("base", base), ("ours", ours), ("theirs", theirs), ("merged", merged))}
    failures: list[str] = []
    report = {"commits": {"base": base, "ours": ours, "theirs": theirs, "merged": merged},
              "conflicts": dict(CONFLICTS), "checks": {}, "failures": failures}

    # ---- 0. nothing vanished ----
    for side in ("ours", "theirs"):
        for path in trees[side]:
            if path not in trees["merged"]:
                failures.append(f"{side}: file dropped by merge: {path}")

    # ---- 1 & 2. additions / one-sided edits are byte-identical ----
    added, onesided, both = [], [], []
    for path in sorted(set(trees["ours"]) | set(trees["theirs"])):
        b = trees["base"].get(path)
        o = trees["ours"].get(path)
        t = trees["theirs"].get(path)
        m = trees["merged"].get(path)
        if path in CONFLICTS:
            both.append(path)
            continue
        if o != b and t == b:      # ours only
            onesided.append(path)
            if m != o:
                failures.append(f"ours-only edit lost: {path}")
        elif t != b and o == b:    # theirs only
            onesided.append(path)
            if m != t:
                failures.append(f"theirs-only edit lost: {path}")
        elif o == t:               # both sides made the same edit
            added.append(path)
            if m != o:
                failures.append(f"identical-both-sides edit changed: {path}")
        elif o != b and t != b:    # both sides edited differently -> git either merged or conflicted
            both.append(path)
            if m is None:
                failures.append(f"both-sides edit vanished: {path}")
    # files added on either side (base absent) must appear verbatim
    for path in trees["ours"]:
        if path not in trees["base"] and path not in CONFLICTS and trees["merged"].get(path) != trees["ours"][path]:
            failures.append(f"ours-added file altered: {path}")
    for path in trees["theirs"]:
        if path not in trees["base"] and path not in CONFLICTS and trees["merged"].get(path) != trees["theirs"][path]:
            failures.append(f"theirs-added file altered: {path}")

    # ---- 3. conflict paths: union of both sides' lines, no markers ----
    conflict_report = {}
    for path, why in CONFLICTS.items():
        mtext = blob(trees["merged"][path])
        for marker in ("<<<<<<<", "=======", ">>>>>>>"):
            if any(ln.startswith(marker) for ln in mtext.splitlines()):
                failures.append(f"{path}: unresolved conflict marker {marker}")
        mlines = nontrivial_lines(mtext)
        per_side = {}
        for side in ("ours", "theirs"):
            side_lines = nontrivial_lines(blob(trees[side][path]))
            base_lines = nontrivial_lines(blob(trees["base"][path]))
            contributed = side_lines - base_lines
            missing = contributed - mlines
            per_side[side] = {"contributed": sum(contributed.values()),
                              "missing": sorted(missing.elements())}
            if missing:
                failures.append(f"{path}: {side} lines missing from resolution: {sorted(missing)}")
        conflict_report[path] = {"why": why, **per_side,
                                 "merged_lines": len(mlines), "identical_to_ours": mtext in (
                                     blob(trees["ours"][path]),),
                                 "identical_to_theirs": mtext in (blob(trees["theirs"][path]),)}
    report["checks"]["conflict_union"] = conflict_report
    report["checks"]["one_sided_edits_verified"] = len(onesided)
    report["checks"]["identical_both_sides"] = len(added)
    report["checks"]["both_sides_edited"] = len(both)

    # ---- 4. public surface is a superset of both sides ----
    surfaces = {}
    for side in ("ours", "theirs", "merged"):
        if "tensorcad/__init__.py" not in trees[side]:
            continue
        text = blob(trees[side]["tensorcad/__init__.py"])
        start = text.index("__all__ = [")
        end = text.index("]", start)
        names = {s.strip().strip('",') for s in text[start:end].replace("__all__ = [", "").split(",") if s.strip()}
        surfaces[side] = names
    for side in ("ours", "theirs"):
        missing = surfaces[side] - surfaces["merged"]
        if missing:
            failures.append(f"__all__ lost names present on {side}: {sorted(missing)}")
    report["checks"]["public_surface"] = {k: len(v) for k, v in surfaces.items()}
    report["checks"]["public_surface_union_ok"] = not any(
        f"__all__ lost" in f for f in failures)

    print(json.dumps(report, indent=2)[:6000])
    with open(REPO + "/artifacts/m3832/merge_verification.json", "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())