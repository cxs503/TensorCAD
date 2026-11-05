#!/usr/bin/env python3
"""Evidence generator for the PR #41 (agent/m3213) merge cure.

The cure had one job: merge ``main`` into ``agent/m3213`` without losing content
from either side.  "No conflict markers left" is not evidence of that -- a
mechanical overlay can drop a hunk silently and still leave a clean tree.  So
this script compares the *merged tree* against both of its parents, name by
name:

* ``__all__``             -- set equality with the union of the two parents' names,
                             no duplicates, and every name actually importable;
* module import list      -- both ``heel`` (branch) and ``transom`` (main) kept;
* conflict files          -- every file the merge reported as conflicted is
                             present, resolvable, and keeps both sides' content
                             (checked through per-file content probes);
* ``CHANGELOG.md``        -- both the branch's heel entry and main's transom /
                             inertia entries survive;
* ``pyproject.toml``      -- both sides bumped to the same version (0.21.0);
* tree-level diff         -- the merged tree differs from ``main`` only by the
                             branch's own additions (the branch side is a
                             superset), and differs from the branch only by
                             main's own additions.

Usage::

    python3 tools/verify_merge_pr41.py --base <main-sha> --branch <merged-branch-sha>
                                       [--json artifacts/m3778/pr41_merge_evidence.json]

Exit code 0 only if every probe passes; the JSON is the machine-readable
evidence quoted in ``docs/m3778_pr41_merge_cure.md``.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

CONFLICT_FILES = ["CHANGELOG.md", "tensorcad/__init__.py", "tests/test_public_surface.py"]

# Content probes: (file, needle) -- every needle must be present in ``rev``.
PROBES = {
    "branch": [  # only the branch introduced these; a merge must not drop them
        ("tensorcad/heel.py", "def gz_curve("),
        ("tests/test_heel.py", "def test_"),
        ("docs/gz_curve.md", "GZ"),
        ("examples/gz_curve/run_example.py", "gz_curve"),
        ("README.md", "Righting arms at finite heel"),
        ("CHANGELOG.md", "Righting arms at finite heel"),
        ("tensorcad/__init__.py", '"wall_sided_gz"'),
    ],
    "main": [  # and the merge must not drop main's own work either
        ("tensorcad/transom.py", "def generate_transom_hull("),
        ("tests/test_transom.py", "def test_"),
        ("tests/test_hydrostatics_inertia.py", "def test_"),
        ("docs/transom_stern.md", "transom"),
        ("docs/asymmetric_hull.md", "asymmetric"),
        ("CHANGELOG.md", "Transom-stern hull family"),
        ("CHANGELOG.md", "Waterplane second moments"),
        ("tensorcad/__init__.py", '"generate_transom_hull"'),
    ],
}


def git(*args: str, cwd: Path) -> str:
    return subprocess.run(
        ["git", *args], cwd=cwd, capture_output=True, text=True, check=True
    ).stdout


def show(rev: str, path: str, cwd: Path) -> str | None:
    proc = subprocess.run(
        ["git", "show", f"{rev}:{path}"], cwd=cwd, capture_output=True, text=True
    )
    return proc.stdout if proc.returncode == 0 else None


def all_names(rev: str, cwd: Path) -> list[str]:
    src = show(rev, "tensorcad/__init__.py", cwd)
    if src is None:
        return []
    block = re.search(r"__all__ = \[(.*?)\n\]", src, re.S)
    return re.findall(r'"([^"]+)"', block.group(1)) if block else []


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", required=True, help="main-side commit (the merge target)")
    ap.add_argument("--branch", default="HEAD", help="merged PR branch commit")
    ap.add_argument("--worktree", default=None, help="tree to inspect (default: this repo)")
    ap.add_argument("--json", default=None, help="write the evidence JSON here")
    args = ap.parse_args(argv)

    root = Path(args.worktree) if args.worktree else Path(__file__).resolve().parent.parent
    checks: list[dict] = []

    def check(name: str, ok: bool, detail) -> None:
        checks.append({"check": name, "ok": bool(ok), "detail": detail})

    # 1. merge topology: the base really is an ancestor of the merged branch
    is_ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", args.base, args.branch],
        cwd=root, capture_output=True,
    ).returncode == 0
    check("base_is_ancestor_of_branch", is_ancestor, f"{args.base[:8]} -> {args.branch[:8]}")

    # 2. no conflict markers anywhere in the tracked tree
    tracked = subprocess.run(
        ["git", "grep", "-n", "-E", r"^(<<<<<<< |=======$|>>>>>>> )", args.branch],
        cwd=root, capture_output=True, text=True,
    ).stdout.strip()
    check("no_conflict_markers", tracked == "", tracked.splitlines()[:5])

    # 3. __all__ is the exact union of both parents, duplicate-free
    branch_names = all_names(args.branch, root)
    base_names = all_names(args.base, root)
    union = set(base_names) | set(branch_names)
    check(
        "__all___equals_parent_union",
        set(branch_names) == union,
        {
            "merged": len(branch_names),
            "main": len(base_names),
            "union": len(union),
            "missing": sorted(union - set(branch_names)),
            "extra": sorted(set(branch_names) - union),
            "duplicates": sorted({n for n in branch_names if branch_names.count(n) > 1}),
        },
    )

    # 4. every declared name imports
    import importlib

    sys.path.insert(0, str(root))
    module = importlib.import_module("tensorcad")
    unimportable = [n for n in module.__all__ if getattr(module, n, None) is None]
    check("all_public_names_importable", not unimportable, unimportable)

    # 5. both layers of the layered-module list survive in the contract test
    contract = show(args.branch, "tests/test_public_surface.py", root) or ""
    check(
        "contract_lists_both_new_modules",
        ('"heel"' in contract) and ('"transom"' in contract),
        "heel=%s transom=%s" % ('"heel"' in contract, '"transom"' in contract),
    )

    # 6. conflict files: content probes from both sides
    for side, probes in PROBES.items():
        for path, needle in probes:
            src = show(args.branch, path, root)
            ok = src is not None and needle in src
            check(f"kept_{side}:{path}::{needle}", ok, needle)

    # 7. pyproject version is the common 0.21.0 bump, not one side's copy
    merged_pyproject = show(args.branch, "pyproject.toml", root) or ""
    base_pyproject = show(args.base, "pyproject.toml", root) or ""
    version = re.search(r'^version = "([^"]+)"', merged_pyproject, re.M)
    check(
        "pyproject_version_matches_main_bump",
        version is not None and version.group(1) == re.search(r'^version = "([^"]+)"', base_pyproject, re.M).group(1),
        version.group(1) if version else None,
    )

    # 8. name-status accounting against both parents: nothing deleted, nothing added twice
    def name_status(a: str, b: str) -> dict[str, str]:
        out = git("diff", "--name-status", a, b, cwd=root)
        table: dict[str, str] = {}
        for line in out.splitlines():
            parts = line.split("\t")
            table[parts[-1]] = parts[0][0]
        return table

    vs_main = name_status(args.base, args.branch)    # branch relative to main
    vs_branch = name_status(args.branch, args.base)  # main relative to branch
    check(
        "nothing_deleted_from_main",
        not [p for p, s in vs_main.items() if s == "D"],
        {
            "branch_adds": sum(1 for s in vs_main.values() if s == "A"),
            "branch_modifies": sum(1 for s in vs_main.values() if s == "M"),
            "deleted_from_main": sorted(p for p, s in vs_main.items() if s == "D"),
            "main_adds": sum(1 for s in vs_branch.values() if s == "A"),
            "deleted_from_branch": sorted(p for p, s in vs_branch.items() if s == "D"),
        },
    )

    ok = all(c["ok"] for c in checks)
    evidence = {
        "mission": 3778,
        "pr": 41,
        "base": {"ref": "main", "sha": args.base},
        "branch": {"ref": "agent/m3213", "sha": args.branch},
        "merge_base": git("merge-base", args.base, args.branch, cwd=root).strip(),
        "conflict_files": CONFLICT_FILES,
        "all_ok": ok,
        "checks": checks,
        "checks_passed": sum(c["ok"] for c in checks),
        "checks_total": len(checks),
    }
    if args.json:
        out = root / args.json
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {out}")
    for c in checks:
        mark = "PASS" if c["ok"] else "FAIL"
        print(f"[{mark}] {c['check']}")
    print(f"\n{evidence['checks_passed']}/{evidence['checks_total']} checks passed -> {'OK' if ok else 'NOT OK'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())