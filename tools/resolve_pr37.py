"""Resolve the four PR #37 conflict files as a semantic union (no side silently dropped).

Run from the repository root while ``git merge origin/agent/m2873`` is in progress::

    python3 tools/resolve_pr37.py

Each conflict is handled by an explicit rule, not by "take ours/theirs":

* ``CHANGELOG.md``      -- keep *both* release sections verbatim, add a 0.21.0 header
                           that explains the union, newest version first.
* ``pyproject.toml``    -- the two sides disagree on a scalar; bump to 0.21.0.
* ``tensorcad/__init__.py`` -- ``__all__`` is the sorted union, ``__version__`` 0.21.0.
* ``tests/test_public_surface.py`` -- EXPECTED_API is the union of both sides, the
                           module list is the union, the version assertion follows the
                           bump.  The union is then checked against the real
                           ``tensorcad.__all__`` so a missed name cannot pass silently.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION = "0.21.0"

MERGE_NOTE = """\
## 0.21.0 — 2026-11-03

The transom-stern family (`tensorcad/transom.py`, authored on `agent/m2873` where it was
versioned 0.19.0 against a 0.16.0 `main`) merged onto the 0.20.0 Wigley line. Both hull
families are present and neither release section below is edited: 0.19.0 was never
released on its own, so nothing is renumbered and no entry is dropped.
"""

BLOCK = re.compile(
    r"<<<<<<< [^\n]*\n(?P<ours>.*?)\n=======\n(?P<theirs>.*?)\n>>>>>>> [^\n]*\n",
    re.S,
)
NAMES = re.compile(r'"([A-Za-z_][A-Za-z0-9_.]*)"')


def wrap(names: list[str], indent: str, width: int) -> str:
    lines, line = [], indent
    for name in names:
        piece = f'"{name}", '
        if len(line) + len(piece) - 1 > width and line.strip():
            lines.append(line.rstrip())
            line = indent
        line += piece
    lines.append(line.rstrip())
    return "\n".join(lines)


def changelog(ours: str, theirs: str) -> str:
    return MERGE_NOTE + "\n" + ours + "\n" + theirs


def pyproject(ours: str, theirs: str) -> str:
    return f'version = "{VERSION}"'


def init_py(ours: str, theirs: str) -> str:
    if "generate_" in ours or "generate_" in theirs:
        names = sorted(set(NAMES.findall(ours)) | set(NAMES.findall(theirs)))
        return wrap(names, "    ", 99)
    return f'__version__ = "{VERSION}"'


def public_surface(ours: str, theirs: str) -> str:
    if "assert" in ours:
        return f'    assert tensorcad.__version__ == "{VERSION}"'
    found = set(NAMES.findall(ours)) | set(NAMES.findall(theirs))
    if "adapters." in ours or "adapters." in theirs:
        ordered = sorted(found, key=lambda s: s.replace("adapters.", "~"))
        return wrap(ordered, "                  ", 99) + "):"
    return wrap(sorted(found), "    ", 95)


HANDLERS = {
    "CHANGELOG.md": changelog,
    "pyproject.toml": pyproject,
    "__init__.py": init_py,
    "test_public_surface.py": public_surface,
}

for rel in ("CHANGELOG.md", "pyproject.toml", "tensorcad/__init__.py",
            "tests/test_public_surface.py"):
    path = ROOT / rel
    text = path.read_text(encoding="utf-8")
    handler = HANDLERS[path.name]
    count = len(BLOCK.findall(text))
    if not count:
        raise SystemExit(f"no conflict markers in {rel}")
    new = BLOCK.sub(
        lambda m, resolve=handler: resolve(m.group("ours"), m.group("theirs")) + "\n",
        text,
    )
    if "<<<<<<<" in new or ">>>>>>>" in new:
        raise SystemExit(f"markers survived in {rel}")
    path.write_text(new, encoding="utf-8")
    print(f"resolved {rel}: {count} conflict block(s)")
