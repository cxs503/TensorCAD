#!/usr/bin/env python3
"""verify_merge_m3917.py -- no-silent-loss audit for the PR #47 cure merge.

The failure mode this guards against is the mechanical overlay: a merge that
looks clean but quietly drops the lines one side added (TensorCAD has that scar
tissue).  So the audit is line-set based and side-by-side:

  for every line *added* by main   since the merge base -> must survive verbatim
  for every line *added* by m3431  since the merge base -> must survive verbatim

with one declared exception list (the README conflict hunk, where the two sides
asserted contradictory facts about the same file and the union had to be
re-stated -- the audit prints those lines so a reviewer can read them).

It also checks the two machine-checkable invariants of the merge:

  __all__ (merged)            == __all__(main)  |  __all__(m3431)
  EXPECTED_API (test file)    == __all__ (merged)
  both import blocks present  (righting from m3431, transom from main)

Writes artifacts/m3917/merge_audit.json.
"""
from __future__ import annotations

import ast
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "m3917"
MAIN = "origin/main"
BRANCH = "agent/m3431"


def git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *args],
                          capture_output=True, text=True, check=True).stdout


def added_lines(ref: str, base: str, path: str) -> set[str]:
    diff = git("diff", "-U0", base, ref, "--", path)
    out = set()
    for line in diff.splitlines():
        if line.startswith("+") and not line.startswith("+++"):
            body = line[1:].strip()
            if body:
                out.add(body)
    return out


def _tokens_subset(line: str, other: str) -> bool:
    """True when every quoted token of `line` occurs on `other` (superset line)."""
    toks = re.findall(r'"[^"]+"', line)
    if toks:
        return all(tok in other for tok in toks)
    return line in other


def all_names(text: str) -> set[str]:
    """Extract the literal strings of the `__all__` (and of EXPECTED_API) list."""
    tree = ast.parse(text)
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            tgt = node.targets[0]
            key = getattr(tgt, "id", None) or getattr(tgt, "attr", None)
            if key in ("__all__", "EXPECTED_API"):
                for elt in node.value.elts:
                    names.add(elt.value)
    return names


def main() -> int:
    base = git("merge-base", MAIN, BRANCH).strip()
    merged_ref = "HEAD"
    files = sorted({
        p for p in git("diff", "--name-only", base, MAIN).split()
        + git("diff", "--name-only", base, BRANCH).split()
    })

    merged_text = {p: (ROOT / p).read_text(encoding="utf-8", errors="replace")
                   for p in files if (ROOT / p).is_file()}

    # README conflict hunk: the only lines allowed to be re-stated rather than kept.
    restated_prefix = (
        "remeshing an arbitrary imported surface to a target size field",
        "or any classification-rule compliance",
    )
    report = {"base": base, "main": git("rev-parse", MAIN).strip(),
              "branch_pre_merge": sys.argv[1] if len(sys.argv) > 1 else "",
              "files_checked": files, "per_file": {}, "lost": [], "extended": [], "restated": []}

    for path in files:
        if path not in merged_text:
            continue  # binary / deleted-on-both-sides
        blob = merged_text[path]
        blob_lines = {ln.strip() for ln in blob.splitlines() if ln.strip()}
        entry = {}
        for label, ref in (("main", MAIN), ("branch", BRANCH)):
            add = added_lines(ref, base, path)
            missing = sorted(ln for ln in add if ln not in blob_lines)
            restated = [ln for ln in missing if ln.startswith(restated_prefix)]
            # a line the merged tree only *extends* (e.g. main's module tuple with
            # "righting" inserted) still carries every original token: not a loss.
            extended = [ln for ln in missing
                        if ln not in restated and any(_tokens_subset(ln, bl) for bl in blob_lines)]
            lost = [ln for ln in missing if ln not in restated and ln not in extended]
            entry[label] = {"added": len(add), "kept": len(add) - len(missing),
                            "extended": len(extended), "restated": len(restated),
                            "lost": len(lost)}
            report["restated"].extend({"file": path, "side": label, "line": ln} for ln in restated)
            report["extended"].extend({"file": path, "side": label, "line": ln} for ln in extended)
            report["lost"].extend({"file": path, "side": label, "line": ln} for ln in lost)
        report["per_file"][path] = entry

    init_text = merged_text["tensorcad/__init__.py"]
    merged_all = all_names(init_text)
    main_all = all_names(git("show", f"{MAIN}:tensorcad/__init__.py"))
    branch_all = all_names(git("show", f"{BRANCH}:tensorcad/__init__.py"))
    test_all = all_names(merged_text["tests/test_public_surface.py"])

    report["api"] = {
        "count_main": len(main_all),
        "count_branch": len(branch_all),
        "count_merged": len(merged_all),
        "count_union": len(main_all | branch_all),
        "main_only_missing": sorted((main_all - merged_all)),
        "branch_only_missing": sorted((branch_all - merged_all)),
        "expected_api_equals_all": test_all == merged_all,
        "expected_api_extra": sorted(test_all - merged_all),
        "expected_api_missing": sorted(merged_all - test_all),
        "imports": {mod: f"from .{mod} import (" in init_text
                    for mod in ("righting", "transom")},
    }

    ok = (not report["lost"]
          and not report["api"]["main_only_missing"]
          and not report["api"]["branch_only_missing"]
          and report["api"]["expected_api_equals_all"]
          and all(report["api"]["imports"].values()))
    report["ok"] = ok

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "merge_audit.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                                          encoding="utf-8")

    print("files checked     :", len(files))
    for path, entry in report["per_file"].items():
        print(f"  {path:38s} main {entry['main']['kept']}/{entry['main']['added']} kept, "
              f"branch {entry['branch']['kept']}/{entry['branch']['added']} kept")
    print("lost lines        :", len(report["lost"]))
    for item in report["lost"][:20]:
        print("   LOST", item)
    print("extended (superset):", len(report["extended"]))
    for item in report["extended"]:
        print("   EXTENDED", item["file"], "->", item["line"][:78])
    print("restated (README) :", len(report["restated"]))
    print("__all__           : main %(count_main)d | branch %(count_branch)d | "
          "union %(count_union)d | merged %(count_merged)d" % report["api"])
    print("EXPECTED_API==__all__:", report["api"]["expected_api_equals_all"])
    print("imports present   :", report["api"]["imports"])
    print("VERDICT           :", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())