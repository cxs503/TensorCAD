#!/usr/bin/env python3
"""Resolve the three conflicted files of PR #41 (agent/m3213 <- main).

Policy: **semantic union**.  Every line/token added by either side must survive,
nothing may be silently dropped (a mechanical ``--ours``/``--theirs``/overlay
would drop one side), and nothing may be duplicated.

Per-file strategy
-----------------
``tensorcad/__init__.py``   both sides appended names to the ``__all__`` literal;
                           keep both insertions (union), ordered within the
                           conflict hunk so the literal stays readable.
``CHANGELOG.md``            both sides opened the same ``## 0.21.0`` section;
                           one section, theirs' ``### Added`` heading, both
                           bullet lists underneath.
``tests/test_public_surface.py``
                           the module tuple gained ``heel`` (ours) and
                           ``transom`` (theirs); token-level union so both
                           modules stay import-covered.

Every resolution is verified: for each conflict hunk the set of non-blank lines
(and, for the token case, quoted names) of *both* sides must be contained in the
resolved file.  A lossy resolution fails the script.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

WS = Path(__file__).resolve().parents[1]
OURS_MARK, SEP, THEIRS_MARK = "<<<<<<<", "=======", ">>>>>>>"


def conflict_blocks(text: str):
    """Return [(start, end, ours_lines, theirs_lines)] for each conflict."""
    lines = text.splitlines(keepends=True)
    out, i = [], 0
    while i < len(lines):
        if lines[i].startswith(OURS_MARK):
            start, ours, theirs, mode = i, [], [], "ours"
            i += 1
            while i < len(lines) and not lines[i].startswith(THEIRS_MARK):
                if lines[i].startswith(SEP):
                    mode = "theirs"
                elif mode == "ours":
                    ours.append(lines[i])
                else:
                    theirs.append(lines[i])
                i += 1
            out.append((start, i, ours, theirs))
        i += 1
    return out


def strip_markers(text: str, resolve) -> str:
    lines = text.splitlines(keepends=True)
    blocks = {b[0]: b for b in conflict_blocks(text)}
    out, i = [], 0
    while i < len(lines):
        if i in blocks:
            _, end, ours, theirs = blocks[i]
            out.extend(resolve(ours, theirs))
            i = end + 1
        else:
            out.append(lines[i])
            i += 1
    return "".join(out)


# ------------------------------------------------------------ __init__.py
def resolve_init(ours, theirs):
    """Union the two insertions, ordered by the name inside the hunk."""
    return sorted(ours + theirs, key=lambda line: line.strip().strip('",'))


# ------------------------------------------------------------- CHANGELOG.md
def resolve_changelog(ours, theirs):
    """Theirs carries the ``### Added`` heading; put both bullet lists under it."""
    return theirs + ["\n"] + ours


# --------------------------------------------------- tests/test_public_surface
def resolve_public_surface(ours, theirs):
    """Token-level union of the module tuple: keep theirs' wrapping, insert the
    names only ours declares right after the token that precedes them in ours."""
    ours_txt, theirs_txt = "".join(ours), "".join(theirs)
    o_tokens = re.findall(r'"([^"\n]+)"', ours_txt)
    t_tokens = re.findall(r'"([^"\n]+)"', theirs_txt)
    out = theirs_txt
    for idx, name in enumerate(o_tokens):
        if name in t_tokens:
            continue
        anchor = next((t for t in reversed(o_tokens[:idx]) if t in t_tokens), None)
        if anchor is None:
            raise SystemExit(f"no anchor for {name!r}")
        m = re.search(r'"%s",\s*' % re.escape(anchor), out)
        if not m:
            raise SystemExit(f"anchor {anchor!r} not found in theirs")
        out = out[: m.end()] + f'"{name}", ' + out[m.end():]
    return [out]


RESOLUTIONS = {
    "tensorcad/__init__.py": resolve_init,
    "CHANGELOG.md": resolve_changelog,
    "tests/test_public_surface.py": resolve_public_surface,
}


def verify(rel: str, resolved: str, blocks) -> list[str]:
    """Both sides' content must be present; no markers left."""
    lost = []
    for marker in (OURS_MARK + " ", ">>>>>>> "):
        if marker in resolved:
            lost.append(f"marker left: {marker!r}")
    res_names = set(re.findall(r'"([^"\n]+)"', resolved))
    res_lines = {line.rstrip("\n").strip() for line in resolved.splitlines()}
    for _, _, ours, theirs in blocks:
        for side in (ours, theirs):
            for line in side:
                s = line.rstrip("\n").strip()
                if not s:
                    continue
                if s not in res_lines:
                    names = re.findall(r'"([^"\n]+)"', line)
                    if names and all(n in res_names for n in names):
                        continue  # token survived, line was re-wrapped
                    lost.append(s[:90])
    return lost


def main() -> int:
    rc = 0
    for rel, fn in RESOLUTIONS.items():
        path = WS / rel
        raw = path.read_text()
        if OURS_MARK not in raw:
            print(f"SKIP {rel} (no conflict)")
            continue
        blocks = conflict_blocks(raw)
        resolved = strip_markers(raw, fn)
        lost = verify(rel, resolved, blocks)
        if lost:
            print(f"FAIL {rel}: {len(lost)} item(s) lost -> {lost[:4]}")
            rc = 1
            continue
        path.write_text(resolved)
        print(f"OK   {rel}: {len(blocks)} hunk(s) unioned, nothing lost")
    return rc


if __name__ == "__main__":
    sys.exit(main())