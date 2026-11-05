#!/usr/bin/env python3
"""PR #41 cure gate: the merge of `main` into `agent/m3213` keeps **both** sides, and the
branch stays under the merge-governor's 200 KB artifact cap.

Why the check is written this way
---------------------------------
The local clone is shallow (`git rev-parse --is-shallow-repository` -> true, history
grafted at the boundary), so `git merge-base 292d794 origin/main` does not resolve and a
textbook three-way comparison is impossible here.  That property -- nothing silently
dropped -- is still decidable from whole trees, so it is decided that way:

  A1/A2  tree level: no path of `main` disappears from the merge, and every path the
         branch owns (present in the branch tip, absent from main) survives byte-exact.
  B      main side: for every path whose merged blob differs from main's, main's lines
         are a sub-multiset of the merged lines, up to the documented allowances below.
         (Multiset, not set: a line cannot be "kept" by one surviving copy when the
         other side also needed it.)
  C      branch side: an explicit inventory of what this branch contributes -- module
         blob, tests, example, docs, CHANGELOG bullet, exported names, module tuple.
         Line-level containment is *not* used for files main has since rewritten (a
         superseded line is supposed to vanish there and would only produce noise).
  D      artifact cap: sum of artifact-class blob sizes in the PR delta against main,
         using the same classifier as `merge_governor.pr_artifact_bytes`.

Deliberate, documented reductions are listed in DELIBERATE; anything undocumented fails
the gate.  Exit status 0 == PR #41 is safe to push.
"""
from __future__ import annotations

import io
import json
import math
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
MERGE = "HEAD"
MAIN = "origin/main"
BRANCH_TIP = "292d794"            # agent/m3213 as GitHub currently sees it (pre-cure)
CAP_BYTES = 200 * 1024
ART_RE = re.compile(r"(outputs|figures|artifacts|__pycache__)/|\.png$|\.npz$|\.pt$|\.vti$|\.pdf$|\.pyc$")

#: path -> (lines allowed to disappear from main's side, reason)
DELIBERATE = {
    "tests/test_public_surface.py": (
        1,
        "the layered-module tuple was re-wrapped onto fewer lines when `heel` was added; "
        "the member list is unchanged and asserted member-by-member in check C",
    ),
    "tensorcad/__init__.py": (
        0,
        "`righting_arm` moved to its alphabetical slot; no name removed",
    ),
    "CHANGELOG.md": (
        0,
        "one heel bullet block kept (both parents carried a copy); text asserted in check C",
    ),
}

#: re-encoded figures: same pixels to within the stated PSNR, re-written to fit the cap
REENCODED = {
    "docs/gz_curve/gz_barge.png": 40.0,
    "docs/gz_curve/gz_sphere.png": 40.0,
    "docs/gz_curve/sphere_ladder.png": 40.0,
}


def git(*args: str, binary: bool = False):
    out = subprocess.run(("git",) + args, cwd=REPO, check=True, capture_output=True)
    return out.stdout if binary else out.stdout.decode()


def ls_tree(rev: str) -> dict[str, str]:
    files = {}
    for line in git("ls-tree", "-r", rev).splitlines():
        meta, path = line.split("\t", 1)
        _mode, kind, sha = meta.split()
        if kind == "blob":
            files[path] = sha
    return files


def lines_of(sha: str) -> Counter:
    return Counter(git("cat-file", "blob", sha, binary=True).decode("utf-8", "replace").splitlines())


def psnr(path: str, sha_ref: str) -> float:
    import numpy as np
    from PIL import Image

    a = np.asarray(Image.open(io.BytesIO(git("cat-file", "blob", sha_ref, binary=True))).convert("RGB"), float)
    b = np.asarray(Image.open(REPO / path).convert("RGB"), float)
    if a.shape != b.shape:
        return -1.0
    mse = ((a - b) ** 2).mean()
    return 99.0 if mse == 0 else 10 * math.log10(255.0 ** 2 / mse)


def main() -> int:
    merge, main_tree, tip = ls_tree(MERGE), ls_tree(MAIN), ls_tree(BRANCH_TIP)
    failures: list[str] = []
    rows: list[tuple] = []

    def record(area: str, detail: str, ok: bool, note: str = "") -> None:
        rows.append((area, detail, "ok" if ok else "FAIL", note))
        if not ok:
            failures.append(f"[{area}] {detail} {note}".rstrip())

    # ---- A. tree level -----------------------------------------------------------
    for path in main_tree:
        if path not in merge:
            record("A1", path, False, "present on main, missing from the merge")
    record("A1", f"{len(main_tree)} main paths present in merge",
           all(p in merge for p in main_tree))
    ours_only = sorted(set(tip) - set(main_tree))
    for path in ours_only:
        if path in REENCODED:      # deliberately re-encoded; judged by PSNR in A2b/B
            continue
        record("A2", path, merge.get(path) == tip[path],
               "branch-owned path changed or lost by the merge")
    for path in sorted(set(ours_only) & set(REENCODED)):
        value = psnr(path, tip[path])
        record("A2b", path, value >= REENCODED[path],
               f"re-encoded vs branch original, PSNR {value:.1f} dB (>= {REENCODED[path]} dB)")
    record("A2", f"{len(ours_only)} branch-owned paths intact", True,
           ", ".join(Path(p).name for p in ours_only[:6]) + " ...")

    # ---- B. main side lines ------------------------------------------------------
    checked = 0
    for path, sha in main_tree.items():
        msha = merge.get(path)
        if msha is None or msha == sha:
            continue
        if path in REENCODED:
            checked += 1
            value = psnr(path, sha)
            record("B", path, value >= REENCODED[path],
                   f"re-encoded figure, PSNR {value:.1f} dB vs main's blob (>= {REENCODED[path]} dB)")
            continue
        want, have = lines_of(sha), lines_of(msha)
        lost = sum(max(0, n - have[ln]) for ln, n in want.items())
        allowed, why = DELIBERATE.get(path, (0, ""))
        checked += 1
        record("B", path, lost <= allowed,
               f"{lost} main line(s) missing, allowance {allowed}"
               + (f" ({why[:70]})" if lost else ""))

    # ---- C. branch inventory -----------------------------------------------------
    def blob_equal(path: str) -> bool:
        return merge.get(path) == tip.get(path)

    for path in ("tensorcad/heel.py", "tests/test_heel.py", "examples/gz_curve/run_example.py",
                 "tensorcad/heel.py", "docs/gz_curve.md", "docs/gz_curve/report.md",
                 "docs/gz_curve/results.json"):
        record("C", f"{path} byte-identical to branch tip", blob_equal(path))

    for path, floor in REENCODED.items():
        record("C", f"{path} re-encoded", merge.get(path) is not None,
               f"PSNR floor {floor} dB asserted in check B")

    def exports(rev: str) -> set[str]:
        text = git("show", f"{rev}:tensorcad/__init__.py")
        block = text.split("__all__ = [", 1)[1].split("]", 1)[0]
        return set(re.findall(r'"([^"]+)"', block))

    e_tip, e_merge = exports(BRANCH_TIP), exports(MERGE)
    missing = sorted(e_tip - e_merge)
    record("C", f"{len(e_tip)} branch __all__ names survive", not missing, f"missing {missing}")

    def modules(rev: str) -> set[str]:
        text = git("show", f"{rev}:tests/test_public_surface.py")
        block = text.split("for module in (", 1)[1].split("):", 1)[0]
        return set(re.findall(r'"([^"]+)"', block))

    m_tip, m_merge = modules(BRANCH_TIP), modules(MERGE)
    record("C", f"{len(m_tip)} branch module-tuple members survive",
           not (m_tip - m_merge), f"missing {sorted(m_tip - m_merge)}")

    changelog = git("show", f"{MERGE}:CHANGELOG.md")
    record("C", "CHANGELOG carries the heel feature bullet",
           changelog.count("Righting arms at finite heel") == 1,
           f"count={changelog.count('Righting arms at finite heel')}")

    # this branch never touched the CLI: assert that main's cli.py is carried over
    # verbatim rather than inventing a subcommand that the branch does not ship.
    record("C", "cli.py is main's blob, untouched by this branch",
           merge.get("tensorcad/cli.py") == main_tree.get("tensorcad/cli.py"))

    # ---- D. artifact cap ---------------------------------------------------------
    delta = [p for p in git("diff", "--name-only", MAIN, MERGE).split() if p]
    artifacts = {p: int(git("cat-file", "-s", f"{MERGE}:{p}")) for p in delta
                 if ART_RE.search(p) and p in merge}
    total = sum(artifacts.values())
    record("D", f"artifact-class bytes in PR delta {total/1024:.1f} KB < {CAP_BYTES/1024:.0f} KB",
           total < CAP_BYTES, ", ".join(f"{Path(p).name}={s/1024:.0f}KB"
                                        for p, s in sorted(artifacts.items())))

    # ---- write evidence ----------------------------------------------------------
    out = REPO / "runs/pr41"
    out.mkdir(parents=True, exist_ok=True)
    with (out / "union_check.tsv").open("w") as fh:
        fh.write("area\tcheck\tstatus\tnote\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    summary = {
        "merge_commit": git("rev-parse", MERGE).strip(),
        "merge_parents": git("rev-list", "--parents", "-n1", MERGE).split()[1:],
        "main": git("rev-parse", MAIN).strip(),
        "branch_tip_github_sees": git("rev-parse", BRANCH_TIP).strip(),
        "shallow_clone": git("rev-parse", "--is-shallow-repository").strip() == "true",
        "method": "base-free whole-tree containment (merge-base unavailable: grafted history)",
        "main_paths": len(main_tree),
        "branch_owned_paths": len(ours_only),
        "main_side_files_checked": checked,
        "artifact_bytes": total,
        "artifact_cap_bytes": CAP_BYTES,
        "artifact_files": {p: s for p, s in sorted(artifacts.items())},
        "checks": len(rows),
        "failures": failures,
        "verdict": "union holds on both sides; PR ready to push" if not failures
                   else "LOSS/DEFECT DETECTED",
    }
    (out / "union_check.json").write_text(json.dumps(summary, indent=2) + "\n")

    for r in rows:
        print(f"  {r[2]:4s} {r[0]}  {r[1]}" + (f"   [{r[3]}]" if r[2] != "ok" and r[3] else ""))
    print(f"\n{len(rows)} checks, {len(failures)} failures -- {summary['verdict']}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())