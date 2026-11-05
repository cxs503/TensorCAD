#!/usr/bin/env python3
"""resolve_merge_m3917.py -- semantic-union resolution of the PR #47 merge.

m3917 cure of TensorCAD PR #47 (branch agent/m3431 onto main).

The merge base is e604f02; both sides independently published 0.21.0:

  main  (since base) : 0.21.0 = waterplane second moments + transom-stern family
  m3431 (one commit) : 0.21.0 = large-angle righting arms (GZ curves)

Three files collide.  Every resolution below is a *union*, i.e. no name, bullet
or sentence added on either side is dropped; where the two sides make
contradictory factual claims about the same release, the claim is re-stated so
that it is true of the merged tree (README "what it does not provide" list).

Run from the worktree root after `git merge origin/main` left conflicts.
"""
from __future__ import annotations

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]

CONFLICT_RE = re.compile(
    r"<<<<<<< [^\n]*\n(.*?)=======\n(.*?)>>>>>>> [^\n]*\n", re.DOTALL
)


def blocks(text: str):
    return CONFLICT_RE.findall(text)


def resolve_init_union() -> str:
    """tensorcad/__init__.py: both sides appended names to ``__all__``.

    Union, keeping the list's rough alphabetical order (main's names first where
    they sort earlier).  Both import blocks were already auto-merged.
    """
    path = ROOT / "tensorcad/__init__.py"
    text = path.read_text(encoding="utf-8")

    # conflict 1: ours = righting helpers, theirs = transom sinc-power helper
    pair1 = (
        '    "solve_plane_offset",\n    "submerged_properties",\n',
        '    "sinc_power_coefficients",\n',
    )
    union1 = (
        '    "sinc_power_coefficients",\n'
        '    "solve_plane_offset",\n'
        '    "submerged_properties",\n'
    )
    # conflict 2: ours = righting helper, theirs = transom/taper helpers
    pair2 = (
        '    "upright_displacement",\n',
        '    "transom_area_ratio",\n    "transom_fraction_for_area_ratio",\n'
        '    "truncated_taper_integral",\n    "truncated_taper_moment",\n',
    )
    union2 = (
        '    "transom_area_ratio",\n'
        '    "transom_fraction_for_area_ratio",\n'
        '    "truncated_taper_integral",\n'
        '    "truncated_taper_moment",\n'
        '    "upright_displacement",\n'
    )

    found = blocks(text)
    if len(found) != 2:
        raise SystemExit("__init__.py: expected 2 conflicts, got %d" % len(found))
    if [f for f in found] != [pair1, pair2]:
        raise SystemExit("__init__.py: conflicts are not the expected pairs:\n%r" % (found,))

    text = CONFLICT_RE.sub(lambda m: union1 if m.group(1) == pair1[0] else union2, text)
    assert "<<<<<<<" not in text and ">>>>>>>" not in text

    # both import blocks must survive (union of *imports*, not only of __all__)
    for block in ("from .righting import (", "from .transom import ("):
        assert block in text, block
    path.write_text(text, encoding="utf-8")
    return "tensorcad/__init__.py"


def resolve_changelog_union() -> str:
    """CHANGELOG.md: one 0.21.0 section per side -> one section, both bullet sets.

    Both sections carry the same heading (``## 0.21.0 — 2026-11-04``), so the
    union is a single release note: main's ``### Added`` block first, then
    m3431's righting bullets appended to the same list.
    """
    path = ROOT / "CHANGELOG.md"
    text = path.read_text(encoding="utf-8")
    found = blocks(text)
    if len(found) != 1:
        raise SystemExit("CHANGELOG.md: expected 1 conflict, got %d" % len(found))
    ours, theirs = found[0]

    merged = theirs.rstrip("\n") + "\n" + ours.rstrip("\n") + "\n"
    text = CONFLICT_RE.sub(lambda m: merged, text, count=1)
    assert "<<<<<<<" not in text and ">>>>>>>" not in text
    path.write_text(text, encoding="utf-8")
    return "CHANGELOG.md"


README_MERGED = (
    "remeshing an arbitrary imported surface to a target size field (mesh conditioning only coarsens; `resolve_resolution` instead regenerates the *parametric* hull at a grid resolved to a stated analytic error), a hull form that is a ship rather than a benchmark (three families exist -- `sin(πx/L)^0.7 · √(z/draft)`, the parabolic Wigley, and the truncation of that hull at a station, closed by a transom plate, which is how real displacement hulls end their aft body -- all parametric, the first two fore/aft-symmetric by default; still no series-60 station data, parallel middle body or appendage *layout*), curvature-driven remeshing or an automatic fairing loop (the fairness audit reports where the surface is unfair, it does not repair it), lightship weight, free-surface or damage stability, stability *criteria* (`GZ` curves exist since 0.21; the rule limits are the caller's), resistance, or any classification-rule compliance (initial GM and the righting-arm curve are reported, the rules that judge them are not). "
    "`metadata.json` carries that as `engineering_status`, and no number in this repository is reported without a check that produced it."
)


def resolve_readme_union() -> str:
    """README.md: both sides rewrote the same 'what it does not provide' list.

    main says the righting-arm curve / GZ does not exist (true of main); m3431
    says it does since 0.21 (true of m3431) and still lists no transom (true of
    m3431).  Both statements are false for the merged tree, so the union keeps
    every *item* from both lists and re-states the two corrected facts: three
    hull families (one with a transom) exist, and GZ exists while the stability
    *criteria* do not.
    """
    path = ROOT / "README.md"
    text = path.read_text(encoding="utf-8")
    found = blocks(text)
    if len(found) != 1:
        raise SystemExit("README.md: expected 1 conflict, got %d" % len(found))
    text = CONFLICT_RE.sub(lambda m: README_MERGED + "\n", text, count=1)
    assert "<<<<<<<" not in text and ">>>>>>>" not in text

    # the merged tree also has a `righting` module to import; main added
    # `transom` to this tuple, m3431 forgot `righting`.
    path.write_text(text, encoding="utf-8")
    return "README.md"


def resolve_public_surface_union() -> str:
    """tests/test_public_surface.py: git auto-merged the name lists; make sure the
    layered-module tuple carries *both* new modules (main: transom, m3431: righting).
    """
    path = ROOT / "tests/test_public_surface.py"
    text = path.read_text(encoding="utf-8")
    old = '"wigley", "transom", "hull",'
    new = '"wigley", "transom", "righting", "hull",'
    assert old in text, "transom entry not found in layered-module tuple"
    text = text.replace(old, new, 1)
    path.write_text(text, encoding="utf-8")
    return "tests/test_public_surface.py"


def main() -> int:
    touched = [
        resolve_init_union(),
        resolve_changelog_union(),
        resolve_readme_union(),
        resolve_public_surface_union(),
    ]
    print("resolved:", ", ".join(touched))
    return 0


if __name__ == "__main__":
    sys.exit(main())