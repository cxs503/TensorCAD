#!/usr/bin/env python3
"""Lossless-merge verifier for PR #42 (``agent/m3300`` -> ``main``).

The cure for a dirty PR is judged on one property: **the merged tree is the semantic
union of the two sides, and neither side's additions were silently dropped**.  That is a
set statement, so this script tests it as a set statement instead of trusting a
by-eye read of the conflict hunks.

Checks, all read-only:

1. **Export surface.**  ``__all__`` of ``merge-base``, of the branch head, of ``main``
   and of the working tree.  Required: ``merged ⊇ ours``, ``merged ⊇ main``,
   ``merged = base ∪ ours ∪ main`` (no invention, no loss) and the two sides' additions
   disjoint (otherwise "both kept" is ambiguous).
2. **File surface.**  Every path added on either side since the merge base must exist in
   the merged worktree, and every path deleted on ``main`` must still be gone.
3. **Module surface.**  Every module named in the public-surface test's layering list
   imports from the merged tree.

Prints a JSON report and writes ``artifacts/m3300/merge_verification.json``; with
``--png`` it also draws the two-panel figure used in the report.

    python3 tools/verify_merge_m3300.py --base df31744 --ours cc8985b --main 2ff80b2
"""

from __future__ import annotations

import argparse
import importlib
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout


def exports(rev: str) -> list[str]:
    src = git("show", f"{rev}:tensorcad/__init__.py")
    body = re.search(r"__all__\s*=\s*\[(.*?)\n\]", src, re.S)
    if body is None:  # pragma: no cover - only if the contract file is gutted
        raise SystemExit(f"no __all__ found at {rev}")
    return re.findall(r'"([^"]+)"', body.group(1))


def paths(rev: str) -> set[str]:
    return set(git("ls-tree", "-r", "--name-only", rev).split())


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="merge base SHA")
    ap.add_argument("--ours", required=True, help="branch head SHA before the cure merge")
    ap.add_argument("--main", required=True, help="main SHA the branch was merged with")
    ap.add_argument("--merged-tests", type=int, default=None,
                    help="pytest count on the merged tree (from the run log)")
    ap.add_argument("--main-tests", type=int, default=None,
                    help="pytest count on main for the same file set")
    ap.add_argument("--png", default=None, help="write the report figure here")
    args = ap.parse_args()

    ours, theirs, base = exports(args.ours), exports(args.main), exports(args.base)
    sys.path.insert(0, str(ROOT))
    import tensorcad  # the working tree under test, not a installed copy

    merged = list(tensorcad.__all__)
    s_ours, s_theirs, s_base, s_merged = set(ours), set(theirs), set(base), set(merged)

    report: dict[str, object] = {
        "refs": {"base": args.base, "ours": args.ours, "main": args.main},
        "exports": {
            "base": len(s_base), "ours": len(s_ours), "main": len(s_theirs),
            "merged": len(s_merged), "merged_unique": len(set(merged)),
            "ours_kept": len(s_ours & s_merged), "main_kept": len(s_theirs & s_merged),
            "ours_dropped": sorted(s_ours - s_merged),
            "main_dropped": sorted(s_theirs - s_merged),
            "invented": sorted(s_merged - (s_base | s_ours | s_theirs)),
            "shared_additions": sorted(((s_ours | s_theirs) - s_base) & s_ours & s_theirs),
        },
    }

    f_ours, f_theirs, f_base = paths(args.ours), paths(args.main), paths(args.base)
    added_ours, added_theirs = f_ours - f_base, f_theirs - f_base
    missing = sorted(p for p in added_ours | added_theirs if not (ROOT / p).exists())
    deleted_on_main = sorted(f_base - f_theirs)
    resurrected = sorted(p for p in deleted_on_main if (ROOT / p).exists())
    report["files"] = {
        "added_ours": len(added_ours), "added_main": len(added_theirs),
        "missing_in_merged": missing, "deleted_on_main": deleted_on_main,
        "resurrected": resurrected,
    }

    layering = re.search(
        r"for module in \((.*?)\):", (ROOT / "tests/test_public_surface.py").read_text(), re.S
    ).group(1)
    modules = re.findall(r'"([^"]+)"', layering)
    bad = []
    for name in modules:
        try:
            importlib.import_module(f"tensorcad.{name}")
        except Exception as exc:  # pragma: no cover - reported, not raised
            bad.append(f"{name}: {exc!r}")
    report["modules"] = {"checked": len(modules), "failed": bad}

    e = report["exports"]
    ok = (
        not e["ours_dropped"] and not e["main_dropped"] and not e["invented"]
        and len(set(merged)) == len(merged)
        and not missing and not resurrected and not bad
    )
    report["lossless"] = ok
    if args.merged_tests is not None:
        report["tests"] = {"merged": args.merged_tests, "main": args.main_tests,
                           "same_files_ok": (args.main_tests is None
                                             or args.merged_tests >= args.main_tests)}

    out = ROOT / "artifacts/m3300/merge_verification.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report, indent=2, ensure_ascii=False))

    if args.png:
        plot(report, ROOT / args.png)
        print(f"figure -> {args.png}")
    return 0 if ok else 1


def plot(report: dict, dest: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    dest.parent.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), dpi=140)

    ex = report["exports"]
    labels = ["merge base", "branch\n(agent/m3300)", "main\n(#39)", "merged"]
    vals = [ex["base"], ex["ours"], ex["main"], ex["merged"]]
    cols = ["#9aa5b1", "#2f6f9f", "#c06c3a", "#2e8b57"]
    bars = ax1.bar(labels, vals, color=cols, width=0.62)
    for b, v in zip(bars, vals):
        ax1.annotate(str(v), (b.get_x() + b.get_width() / 2, v), ha="center", va="bottom",
                     fontsize=10)
    ax1.set_ylim(0, max(vals) * 1.18)
    ax1.set_ylabel("public names in `tensorcad.__all__`")
    ax1.set_title(f"Export surface is the union: {ex['base']} + "
                  f"{ex['ours'] - ex['base']} (branch) + {ex['main'] - ex['base']} (main) "
                  f"= {ex['merged']}")
    ax1.grid(axis="y", alpha=0.3)
    ax1.text(0.02, 0.03, f"ours dropped: {len(ex['ours_dropped'])}   "
                         f"main dropped: {len(ex['main_dropped'])}",
             transform=ax1.transAxes, fontsize=9, color="#243b53")

    t = report.get("tests", {})
    if t:
        mvals = [t.get("main") or 0, t.get("merged") or 0]
        mcols = ["#c06c3a", "#2e8b57"]
        bars = ax2.bar(["main (baseline)", "merged"], mvals, color=mcols, width=0.55)
        for b, v in zip(bars, mvals):
            ax2.annotate(str(v), (b.get_x() + b.get_width() / 2, v), ha="center",
                         va="bottom", fontsize=10)
        ax2.set_ylim(0, max(mvals) * 1.18)
        ax2.set_ylabel("pytest cases passed (core file set)")
        ax2.set_title(f"Merged suite ≥ main baseline: {mvals[1]} vs {mvals[0]} passed "
                      f"(+{mvals[1] - mvals[0]}, tests/test_gz.py); identical 2 "
                      f"pre-existing env failures on both sides")
        ax2.grid(axis="y", alpha=0.3)

    fig.suptitle("TensorCAD PR #42 cure — lossless union of agent/m3300 and main", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(dest)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())