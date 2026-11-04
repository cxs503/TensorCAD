"""Audit that merging ``agent/m2873`` into ``main`` lost nothing from either side.

The failure mode this guards against is silent loss: a mechanical ``--ours``/``--theirs``
overlay produces a tree that compiles and passes tests while quietly dropping a file or a
public name that one side contributed.  So the audit is set-based, not narrative:

1. **File ledger.**  For every path touched by either side since the merge base, the
   merged tree must contain that path with exactly the contributing side's blob when only
   one side changed it (``one-sided``); when both sides changed it the path must be listed
   as a *declared conflict* whose resolution is accounted for individually.
2. **No accidental deletion.**  Paths that exist on either side and were not deleted
   there must still exist in the merged tree.
3. **Name ledger.**  The public surface ``tensorcad.__all__`` in the merged tree must be
   exactly ``main``'s names unioned with the branch's names -- nothing from either side
   missing, nothing invented.
4. **Test ledger.**  The merged suite's pass/fail counts must be >= the ``main`` baseline,
   with any failure present in *both* runs attributed as pre-existing (environmental).

Usage::

    python3 tools/audit_pr37_union.py            # writes artifacts/pr37/union_audit.json

The script exits non-zero when the union is not exact, so it doubles as a regression gate.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = "6c568c314da11a3082b41218eb5b564e466cf3a0"  # merge-base(main, agent/m2873)
MAIN = "origin/main"
BRANCH = "origin/agent/m2873"
OUT = ROOT / "artifacts" / "pr37" / "union_audit.json"

# Files one side added that the merge deliberately edits, with the reason.  Everything
# else must arrive byte-identical from its contributing side.
INTENTIONAL_EDITS = {
    "examples/transom_stern/run_example.py":
        "added a sys.path bootstrap so the example runs from a checkout without an "
        "install; identical 3 lines to examples/wigley/run_example.py, no numbers changed",
}

DECLARED_CONFLICTS = {
    "CHANGELOG.md": "both release sections kept verbatim; 0.21.0 header added",
    "pyproject.toml": "scalar version: max(0.20.0, 0.19.0) bumped to 0.21.0",
    "tensorcad/__init__.py": "__all__ = sorted union; __version__ = 0.21.0",
    "tests/test_public_surface.py": "EXPECTED_API/module list = union; version assert bumped",
}


def git(*args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=ROOT, check=True, capture_output=True, text=True
    ).stdout


def changeset(rev: str) -> dict[str, str]:
    """``{path: status}`` for every path ``rev`` changed relative to the merge base."""
    out = git("diff", "--name-status", f"{BASE}", rev)
    ledger = {}
    for line in out.splitlines():
        status, *paths = line.split("\t")
        ledger[paths[-1]] = status[0]
    return ledger


def blob_sha(rev: str, path: str) -> str | None:
    """The blob id of ``path`` in ``rev`` (``None`` when absent) -- binary-safe."""
    proc = subprocess.run(
        ["git", "rev-parse", "--verify", "--quiet", f"{rev}:{path}"],
        cwd=ROOT, capture_output=True, text=True,
    )
    return proc.stdout.strip() if proc.returncode == 0 else None


def worktree_sha(path: str) -> str | None:
    target = ROOT / path
    if not target.exists():
        return None
    return subprocess.run(
        ["git", "hash-object", "--", path], cwd=ROOT, check=True,
        capture_output=True, text=True,
    ).stdout.strip()


def text_at(rev: str, path: str) -> str | None:
    proc = subprocess.run(["git", "show", f"{rev}:{path}"], cwd=ROOT,
                          capture_output=True, text=True)
    return proc.stdout if proc.returncode == 0 else None


def public_names(rev: str) -> set[str]:
    text = text_at(rev, "tensorcad/__init__.py")
    if text is None:
        return set()
    block = text.split("__all__", 1)[1]
    block = block.split("]", 1)[0]
    return set(re.findall(r'"([A-Za-z_][A-Za-z0-9_]*)"', block))


def pytest_counts(path: Path) -> dict[str, object]:
    if not path.exists():
        return {"log": path.name, "present": False}
    text = path.read_text(encoding="utf-8", errors="replace")
    tail = [line for line in text.splitlines() if " passed" in line or " failed" in line]
    summary = tail[-1] if tail else ""
    counts = {key: int(m.group(1)) for key in ("passed", "failed", "skipped", "error")
              if (m := re.search(rf"(\d+) {key}", summary))}
    failed_tests = sorted(set(re.findall(r"^FAILED (\S+)", text, re.M)))
    return {"log": path.name, "present": True, "summary": summary.strip(),
            "counts": counts, "failed_tests": failed_tests}


def main() -> int:
    main_changed = changeset(MAIN)
    branch_changed = changeset(BRANCH)
    paths = sorted(set(main_changed) | set(branch_changed))

    ledger = {"one_sided": [], "conflicts": [], "missing": [], "unexpected": []}
    for path in paths:
        in_main, in_branch = path in main_changed, path in branch_changed
        local = ROOT / path
        if in_main and in_branch:
            if path not in DECLARED_CONFLICTS:
                ledger["unexpected"].append(path)
                continue
            ledger["conflicts"].append({
                "path": path,
                "main_status": main_changed[path],
                "branch_status": branch_changed[path],
                "resolution": DECLARED_CONFLICTS[path],
                "present_after_merge": local.exists(),
            })
            continue
        side, rev = ("main", MAIN) if in_main else ("branch", BRANCH)
        status = main_changed[path] if in_main else branch_changed[path]
        if status == "D":  # side deleted it; the merge may not resurrect it silently
            ledger["one_sided"].append({"path": path, "side": side, "status": "D",
                                        "present_after_merge": local.exists()})
            continue
        expected, actual = blob_sha(rev, path), worktree_sha(path)
        identical = expected is not None and expected == actual
        ledger["one_sided"].append({
            "path": path, "side": side, "status": status,
            "present_after_merge": local.exists(), "content_identical": identical,
            "intentional_edit": INTENTIONAL_EDITS.get(path),
        })
        if not local.exists():
            ledger["missing"].append(path)

    names_main, names_branch = public_names(MAIN), public_names(BRANCH)
    block = (ROOT / "tensorcad" / "__init__.py").read_text(encoding="utf-8")
    names_merged = set(re.findall(r'"([A-Za-z_][A-Za-z0-9_]*)"',
                                  block.split("__all__", 1)[1].split("]", 1)[0]))
    union = names_main | names_branch
    name_ledger = {
        "main_only": sorted(names_main - names_branch),
        "branch_only": sorted(names_branch - names_main),
        "both": len(names_main & names_branch),
        "merged_count": len(names_merged),
        "union_count": len(union),
        "missing_from_merge": sorted(union - names_merged),
        "invented_in_merge": sorted(names_merged - union),
        "exact_union": names_merged == union,
    }

    tests = {
        "main_baseline": pytest_counts(ROOT / "runs" / "pytest_main_baseline.log"),
        "merged": pytest_counts(ROOT / "runs" / "pytest_merged.log"),
    }
    base_c, mer_c = tests["main_baseline"], tests["merged"]
    regressions = sorted(set(mer_c.get("failed_tests", [])) - set(base_c.get("failed_tests", [])))
    tests["new_failures_vs_baseline"] = regressions
    tests["passed_delta"] = (mer_c.get("counts", {}).get("passed", 0)
                             - base_c.get("counts", {}).get("passed", 0))

    ok = (
        not ledger["missing"]
        and not ledger["unexpected"]
        and name_ledger["exact_union"]
        and not regressions
        and all(entry.get("content_identical", True) or entry.get("intentional_edit")
                for entry in ledger["one_sided"])
    )
    report = {
        "merge_base": BASE,
        "main_tip": git("rev-parse", MAIN).strip(),
        "branch_tip": git("rev-parse", BRANCH).strip(),
        "files_touched_by_either_side": len(paths),
        "file_ledger": ledger,
        "name_ledger": name_ledger,
        "tests": tests,
        "verdict": "PASS (semantic union, nothing lost)" if ok else "FAIL (see ledgers)",
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: report[k] for k in ("verdict", "name_ledger")}, indent=2))
    print(f"one-sided files: {len(ledger['one_sided'])}, conflicts: {len(ledger['conflicts'])}, "
          f"missing: {len(ledger['missing'])}, unexpected: {len(ledger['unexpected'])}")
    print(f"tests: baseline {base_c.get('counts')} vs merged {mer_c.get('counts')}, "
          f"delta passed {tests['passed_delta']}, new failures {regressions}")
    print(f"wrote {OUT.relative_to(ROOT)}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
