#!/usr/bin/env python3
"""Resolve the CHANGELOG.md conflict left by the m2925 cure merge (verbatim union).

Both sides shipped a 0.14.0 section for the same feature developed independently
(main: plane-section route; m2925: flux-identity route + analytic reference). The
resolution keeps *both* bullet lists verbatim, merged into a single 0.14.0 section,
with only the three qualifiers that the merge itself made necessary:

  * each side's test bullet is labelled with its route, so the two counts (19 / 18)
    are not read as contradictory -- the merged suite is 31 functions / 37 checks;
  * m2925's example bullet names the artifacts it really writes after the merge
    (`results_hydrostatics.json`, `report_hydrostatics.md`, `plot.png`), because the
    merge gave the flux route its own files instead of overwriting main's.
"""
import pathlib
import subprocess

def blob(rev, path="CHANGELOG.md"):
    return subprocess.run(["git", "show", f"{rev}:{path}"],
                          capture_output=True, text=True, check=True).stdout

def section(text, start="## 0.14.0 — 2026-11-03", stop="## 0.13.0"):
    return text.split(start, 1)[1].split(stop, 1)[0]

def bullets(block):
    """Split a changelog body into top-level '- ' bullets (continuations kept)."""
    out, cur = [], None
    for line in block.splitlines():
        if line.startswith("- "):
            if cur is not None:
                out.append(cur)
            cur = [line]
        elif cur is not None:
            if line.strip():
                cur.append(line)
            elif cur:
                out.append(cur)
                cur = None
    if cur:
        out.append(cur)
    return ["\n".join(b).rstrip() for b in out]

MAIN, BRANCH = "311f45a", "ad86f4f"
main_sec = section(blob(MAIN))
branch_sec = section(blob(BRANCH))

# sub-sections per side
def sub(block, name):
    head = f"### {name}"
    part = block.split(head, 1)[1]
    for nxt in ("### Added", "### Measured", "### Changed"):
        if nxt in part:
            part = part.split(nxt, 1)[0]
    return part

main_added = bullets(sub(main_sec, "Added"))
main_meas = bullets(sub(main_sec, "Measured"))
br_added = bullets(sub(branch_sec, "Added"))
br_meas = bullets(sub(branch_sec, "Measured"))

# --- the only edits the merge forces (see module docstring) -------------------
repl = {
    "- `tests/test_stability.py` — 18 tests.":
        "- `tests/test_stability.py` — 18 tests for this route "
        "(the merged file carries both suites, 31 tests).",
    "- `tests/test_stability.py` adds 19 tests; `docs/stability.md` derives the formulas, reports\n"
    "  the ladder and names the limits.":
        "- `tests/test_stability.py` adds 19 tests for the plane-section route "
        "(31 with the flux-identity suite); `docs/stability.md` derives the formulas, reports\n"
        "  the ladder and names the limits.",
    "- `examples/stability/run_example.py` — ladder 21×9 … 321×81, writes\n"
    "  `docs/stability/{results.json,report.md,plot.png}`, exits 1 if a claim fails.":
        "- `examples/stability/run_example_hydrostatics.py` — ladder 21×9 … 321×81 through the\n"
        "  flux-identity route, writes `docs/stability/{results_hydrostatics.json,\n"
        "  report_hydrostatics.md,plot.png}`, exits 1 if a claim fails.",
}
def apply_edits(items):
    return [repl.get(b, b) for b in items]

br_added = apply_edits(br_added)
main_added = apply_edits(main_added)

def dedup(*groups):
    seen, out = set(), []
    for g in groups:
        for b in g:
            key = " ".join(b.split())
            if key in seen:
                continue
            seen.add(key)
            out.append(b)
    return out

merged = ["## 0.14.0 — 2026-11-03",
          "",
          "Two independently developed routes to the same metacentric heights landed together;",
          "the merge keeps both, because the second route is the cheapest cross-check the data",
          "structure has (see `docs/stability.md`, \"Second route\").",
          "",
          "### Added"]
merged += dedup(main_added, br_added)
merged += ["", "### Measured"]
merged += dedup(main_meas, br_meas)

# --- splice into the conflicted file -----------------------------------------
path = pathlib.Path("CHANGELOG.md")
lines = path.read_text().splitlines()
i_head = lines.index("<<<<<<< HEAD")
i_eq = lines.index("=======")
i_tail = lines.index(">>>>>>> origin/main")
header = lines[:i_head]
main_block = lines[i_eq + 1:i_tail]
footer = lines[i_tail + 1:]
i_014 = main_block.index("## 0.14.0 — 2026-11-03")

text = "\n".join(header + main_block[:i_014] + merged + footer).rstrip("\n") + "\n"
assert not any(m in text for m in ("<<<<<<<", "=======", ">>>>>>>"))
assert text.count("## 0.14.0 — 2026-11-03") == 1
path.write_text(text)
print("wrote CHANGELOG.md:", len(text.splitlines()), "lines,",
      len(dedup(main_added, br_added)), "Added bullets,", len(dedup(main_meas, br_meas)),
      "Measured bullets")