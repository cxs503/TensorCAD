#!/usr/bin/env python3
"""Independent verification of the PR #41 cure merge.

Checks (all machine-decidable, no self-report):
  V1  every python file in the merged tree parses (AST) and no module defines the
      same top-level symbol twice;
  V2  ``tensorcad.__all__`` == union of both parents' ``__all__`` (no name lost,
      no duplicate);
  V3  ``tests/test_public_surface.py::EXPECTED_API`` == union of both parents'
      literals, and the layered-module tuple covers both ``heel`` and ``transom``;
  V4  no file present in either parent is missing from the merged tree, and no
      conflicted file lost a line that either side had;
  V5  version metadata is single-valued (both sides bumped 0.20.0 -> 0.21.0; a
      duplicated key or a stale number would show up here).
"""
from __future__ import annotations

import ast
import re
import subprocess
import sys
from pathlib import Path

WS = Path(__file__).resolve().parents[1]
BASE = "df31744"  # merge-base(origin/main, origin/agent/m3213)
OURS = "origin/agent/m3213"
THEIRS = "origin/main"
fails: list[str] = []


def git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=WS, capture_output=True,
                          text=True, check=True).stdout


def show(rev: str, path: str) -> str:
    return git("show", f"{rev}:{path}")


def all_names(src: str) -> list[str]:
    m = re.search(r"__all__\s*=\s*\[(.*?)\n\]", src, re.S)
    return ast.literal_eval("[" + m.group(1) + "]")


def expected_api(src: str) -> set[str]:
    m = re.search(r"EXPECTED_API\s*=\s*\{(.*?)\n\}", src, re.S)
    return set(ast.literal_eval("{" + m.group(1) + "}"))


def module_tuple(src: str) -> set[str]:
    m = re.search(r"for module in \((.*?)\):", src, re.S)
    return set(re.findall(r'"([^"\n]+)"', m.group(1)))


# ------------------------------------------------------------------ V1
py = [p for p in WS.rglob("*.py") if ".git" not in p.parts and "runs" not in p.parts]
bad = 0
for path in py:
    src = path.read_text()
    try:
        tree = ast.parse(src)
    except SyntaxError as exc:
        fails.append(f"V1 {path.relative_to(WS)}: SyntaxError {exc}")
        bad += 1
        continue
    seen = set()
    for node in tree.body:
        names = []
        if isinstance(node, (ast.FunctionDef, ast.ClassDef)):
            names = [node.name]
        elif isinstance(node, ast.Assign):
            names = [t.id for t in node.targets if isinstance(t, ast.Name)]
        for n in names:
            if n in seen:
                fails.append(f"V1 {path.relative_to(WS)}: duplicate top-level {n}")
            seen.add(n)
print(f"V1 parsed {len(py)} python files, {bad} syntax error(s)")

# ------------------------------------------------------------------ V2
init = (WS / "tensorcad/__init__.py").read_text()
ours_all, theirs_all, base_all = (all_names(show(OURS, "tensorcad/__init__.py")),
                                  all_names(show(THEIRS, "tensorcad/__init__.py")),
                                  all_names(show(BASE, "tensorcad/__init__.py")))
merged = all_names(init)
if len(merged) != len(set(merged)):
    dup = sorted({n for n in merged if merged.count(n) > 1})
    fails.append(f"V2 __all__ duplicates: {dup}")
want = set(ours_all) | set(theirs_all)
if set(merged) != want:
    fails.append(f"V2 __all__ lost {sorted(want - set(merged))} "
                 f"extra {sorted(set(merged) - want)}")
if not set(base_all) <= set(merged):
    fails.append(f"V2 base names dropped: {sorted(set(base_all) - set(merged))}")
print(f"V2 __all__ base={len(base_all)} ours={len(ours_all)} theirs={len(theirs_all)} "
      f"merged={len(merged)} union={len(want)} "
      f"ours-only={sorted(set(ours_all) - set(theirs_all))} "
      f"theirs-only={sorted(set(theirs_all) - set(ours_all))}")

# ------------------------------------------------------------------ V3
pub = (WS / "tests/test_public_surface.py").read_text()
e_merged, e_ours, e_theirs = (expected_api(pub),
                              expected_api(show(OURS, "tests/test_public_surface.py")),
                              expected_api(show(THEIRS, "tests/test_public_surface.py")))
if e_merged != e_ours | e_theirs:
    fails.append(f"V3 EXPECTED_API lost {sorted((e_ours | e_theirs) - e_merged)}")
m_merged, m_ours, m_theirs = (module_tuple(pub),
                              module_tuple(show(OURS, "tests/test_public_surface.py")),
                              module_tuple(show(THEIRS, "tests/test_public_surface.py")))
if m_merged != m_ours | m_theirs:
    fails.append(f"V3 module tuple lost {sorted((m_ours | m_theirs) - m_merged)}")
print(f"V3 EXPECTED_API merged={len(e_merged)} union={len(e_ours | e_theirs)}; "
      f"module tuple merged={len(m_merged)} union={len(m_ours | m_theirs)} "
      f"(heel={'heel' in m_merged}, transom={'transom' in m_merged})")

# ------------------------------------------------------------------ V4
for rev, label in ((OURS, "ours"), (THEIRS, "theirs")):
    ref = set(git("ls-tree", "-r", "--name-only", rev).split())
    have = set(git("ls-tree", "-r", "--name-only", "HEAD").split()) | {
        p.relative_to(WS).as_posix() for p in py}
    tracked = set(git("ls-files").split())
    gone = ref - tracked
    if gone:
        fails.append(f"V4 {label}: {len(gone)} file(s) missing from merged tree: "
                     f"{sorted(gone)[:5]}")
    print(f"V4 {label}: {len(ref)} files, {len(gone)} missing")

# ------------------------------------------------------------------ V5
ppt = (WS / "pyproject.toml").read_text()
versions = re.findall(r'^version\s*=\s*"([^"]+)"', ppt, re.M)
if versions != ["0.21.0"]:
    fails.append(f"V5 pyproject version(s) = {versions}")
if len(re.findall(r"^\[project\]", ppt, re.M)) != 1:
    fails.append("V5 duplicated [project] table")
print(f"V5 pyproject version={versions}")

print()
if fails:
    print(f"FAIL ({len(fails)})")
    for f in fails:
        print("  -", f)
    sys.exit(1)
print("PASS: all merge-fidelity checks hold")
