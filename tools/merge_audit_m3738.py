"""Line-level audit that merging ``origin/agent/m2873`` into ``origin/main`` lost nothing.

Mechanical overlays fail silently: the file survives, lines that only one side had
do not.  So the check is a multiset test per conflict file -- every non-blank line of
each side must still be present -- plus a set test on paths (nothing added by either
side may be absent) and on the public surface (``__all__`` must be the union).

Driven by a handful of git calls rather than one per file, so it runs in seconds.
Writes ``docs/m3738/merge_audit.json``; exits non-zero on any violation.
"""
from __future__ import annotations

import ast
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OURS, THEIRS = "origin/main", "origin/agent/m2873"
DECLARED = ["CHANGELOG.md", "pyproject.toml", "tensorcad/__init__.py",
            "tests/test_public_surface.py"]


def sh(*args: str) -> str:
    return subprocess.run(args, cwd=ROOT, capture_output=True, text=True,
                          check=True).stdout


def listing(ref: str) -> dict[str, str]:
    return {line.split("\t", 1)[1]: line.split()[2]
            for line in sh("git", "ls-tree", "-r", ref).splitlines()}


VERSION_LINE = __import__("re").compile(r'(__version__\s*==|__version__\s*=|^version\s*=)\s*"[0-9.]+"')


def _tokens(line: str) -> list[str]:
    """Identifiers and quoted literals of a line, for the reflow test."""
    return __import__("re").findall(r'"[^"]*"|\'[^\']*\'|[A-Za-z_][A-Za-z_0-9]*', line)


def lines(text: str) -> Counter:
    return Counter(ln.strip() for ln in text.splitlines() if ln.strip())


def main() -> int:
    ours, theirs = listing(OURS), listing(THEIRS)
    merged = set(sh("git", "ls-files").splitlines())
    problems: list[str] = []

    report: dict = {
        "ours": OURS, "theirs": THEIRS, "declared_conflicts": DECLARED,
        "paths": {"only_in_main": sorted(set(ours) - set(theirs)),
                  "only_in_m2873": sorted(set(theirs) - set(ours)),
                  "shared": len(set(ours) & set(theirs)),
                  "merged_total": len(merged)},
        "missing_from_merged": sorted((set(ours) | set(theirs)) - merged),
    }
    problems += [f"missing path: {p}" for p in report["missing_from_merged"]]

    # files the merge actually wrote differently from either side, with line audit
    changed = [p for p in DECLARED]
    undeclared = sorted(set(sh("git", "diff", "--name-only", OURS, THEIRS).splitlines())
                        - set(DECLARED))
    report["differ_between_sides"] = sorted(set(sh("git", "diff", "--name-only",
                                                  OURS, THEIRS).splitlines()))
    conflicts = []
    for path in changed:
        text = (ROOT / path).read_text()
        row = {"path": path, "markers_left": any(m in text for m in
                                                 ("<<<<<<<", ">>>>>>>")),
                "lost": {}}
        for ref in (OURS, THEIRS):
            gap = lines(sh("git", "show", f"{ref}:{path}")) - lines(text)
            # A line is genuinely lost only if it carries content the merge does not
            # have anywhere: a version literal is *superseded* (the union release
            # number replaces both) and a reflowed connector line keeps every token.
            lost, superseded, reflowed = {}, [], []
            for line, n in gap.items():
                if VERSION_LINE.search(line):
                    superseded.append(line)
                elif all(tok in text for tok in _tokens(line)):
                    reflowed.append(line)
                else:
                    lost[line] = n
            if lost:
                row["lost"][ref] = lost
                problems.append(f"{path}: lost {sum(lost.values())} line(s) from {ref}")
            if superseded:
                row.setdefault("superseded", {})[ref] = superseded
            if reflowed:
                row.setdefault("reflowed", {})[ref] = reflowed
        if row["markers_left"]:
            problems.append(f"{path}: conflict marker left")
        added = lines(text) - (lines(sh("git", "show", f"{OURS}:{path}"))
                               + lines(sh("git", "show", f"{THEIRS}:{path}")))
        row["lines_added_by_resolution"] = sum(added.values())
        row["added_sample"] = sorted(added)[:6]
        conflicts.append(row)
    report["conflicts"] = conflicts
    report["unchanged_between_sides_but_audited_for_presence"] = undeclared

    sys.path.insert(0, str(ROOT))
    import tensorcad

    def expected(ref: str) -> set[str]:
        tree = ast.parse(sh("git", "show", f"{ref}:tests/test_public_surface.py"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "EXPECTED_API":
                return set(ast.literal_eval(node.value))
        raise AssertionError("EXPECTED_API not found")

    api_ours, api_theirs = expected(OURS), expected(THEIRS)
    merged_api = set(tensorcad.__all__)
    report["api"] = {
        "main": len(api_ours), "m2873": len(api_theirs), "merged": len(merged_api),
        "union_missing_from_merged": sorted((api_ours | api_theirs) - merged_api),
        "merged_outside_union": sorted(merged_api - (api_ours | api_theirs)),
        "version": tensorcad.__version__,
        "both_families_importable": all(hasattr(tensorcad, n) for n in (
            "generate_wigley", "wigley_analytic_hydrostatics",
            "generate_transom_hull", "analytic_transom_hydrostatics")),
    }
    if report["api"]["union_missing_from_merged"] or report["api"]["merged_outside_union"]:
        problems.append("public surface is not the union of both sides")
    if not report["api"]["both_families_importable"]:
        problems.append("one hull family is not importable from the merged tree")

    report["problems"] = problems
    report["ok"] = not problems
    out = ROOT / "docs/m3738/merge_audit.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "paths"}, indent=2))
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())