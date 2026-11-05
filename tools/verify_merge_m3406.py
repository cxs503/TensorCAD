"""Verify that merging ``main`` into ``agent/m3406`` (PR #46) lost nothing.

The cure for PR #46 is a *semantic union*: ``main`` brought the transom-stern and
asymmetric-hull families (plus the waterplane-inertia fields on ``Hydrostatics``) and
the branch brought the righting-arm / GZ route.  Two files were edited on both sides
(``tensorcad/__init__.py`` and ``tests/test_public_surface.py``); everything else was
touched by a single side and must therefore be byte-identical to that side after the
merge.

This script turns the merge claim into numbers, so it can be re-checked on any later
tree instead of being taken on trust:

1. **No loss on single-side files** -- for every path changed by exactly one side,
   the working tree must equal that side's blob byte for byte.
2. **Public surface is the exact union** -- ``tensorcad.__all__`` must equal
   ``base | main | branch`` with no duplicates, and ``tests/test_public_surface.py``
   must assert exactly that set.
3. **Branch numbers reproduce** -- ``examples/stability/run_gz.py`` is re-run and its
   JSON compared, leaf by leaf, against the branch's committed
   ``docs/stability/gz_results.json`` (the 0.20.0 baseline).
4. **Main numbers reproduce** -- ``examples/transom_stern/run_example.py`` is re-run
   and its three CSVs compared against ``docs/m2873/outputs/transom_stern``.

Run with ``--full`` to re-run both examples (about 3 minutes); without it only the
static checks (1, 2) are done.  Exit code is non-zero if any check fails.

    python tools/verify_merge_m3406.py --full --json artifacts/m3924/merge_verification.json
"""

from __future__ import annotations

import argparse
import ast
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
#: refs the PR #46 cure merged: 0900b21 = old agent/m3406 head, 2ff80b2 = main head
DEFAULT_BRANCH = "0900b21"
DEFAULT_MAIN = "2ff80b2"
NUMERIC_REL_TOL = 1e-12


def _git(*args: str) -> bytes:
    return subprocess.run(("git",) + args, cwd=REPO, capture_output=True, check=True).stdout


def _git_text(*args: str) -> str:
    return _git(*args).decode()


def _show(ref: str, path: str) -> bytes:
    return _git("show", f"{ref}:{path}")


def _api(ref: str) -> list[str]:
    src = _show(ref, "tensorcad/__init__.py").decode()
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", None) == "__all__":
            return [ast.literal_eval(elt) for elt in node.value.elts]
    raise SystemExit(f"no __all__ in {ref}:tensorcad/__init__.py")


def _expected_api(ref: str) -> set[str]:
    src = _show(ref, "tests/test_public_surface.py").decode()
    for node in ast.parse(src).body:
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", None) == "EXPECTED_API":
            return set(ast.literal_eval(node.value))
    raise SystemExit(f"no EXPECTED_API in {ref}:tests/test_public_surface.py")


def check_single_side_files(base: str, main: str, branch: str) -> dict[str, Any]:
    main_files = set(_git_text("diff", "--name-only", base, main).split())
    branch_files = set(_git_text("diff", "--name-only", base, branch).split())
    report: dict[str, Any] = {"overlap": sorted(main_files & branch_files), "sides": {}}
    for tag, files, ref in (
        ("main_only", main_files - branch_files, main),
        ("branch_only", branch_files - main_files, branch),
    ):
        mismatched = []
        for path in sorted(files):
            on_disk = REPO / path
            if not on_disk.exists():
                mismatched.append([path, "missing from worktree"])
            elif on_disk.read_bytes() != _show(ref, path):
                mismatched.append([path, "differs from " + ref])
        report["sides"][tag] = {"files": len(files), "mismatched": mismatched}
    report["ok"] = not any(side["mismatched"] for side in report["sides"].values())
    return report


def check_public_surface(base: str, main: str, branch: str) -> dict[str, Any]:
    union = set(_api(base)) | set(_api(main)) | set(_api(branch))
    merged = _api("HEAD")
    expected_now = _expected_api("HEAD")
    expected_union = _expected_api(main) | _expected_api(branch)
    duplicates = sorted({name for name in merged if merged.count(name) > 1})
    info = {
        "counts": {
            "base": len(set(_api(base))),
            "main": len(set(_api(main))),
            "branch": len(set(_api(branch))),
            "union": len(union),
            "merged_all": len(merged),
            "expected_api_merged": len(expected_now),
        },
        "missing_from_all": sorted(union - set(merged)),
        "extra_in_all": sorted(set(merged) - union),
        "duplicates": duplicates,
        "expected_api_missing": sorted(expected_union - expected_now),
        "expected_api_extra": sorted(expected_now - expected_union),
        "both_sides_symbols_present": sorted(
            (set(_api(main)) - set(_api(base))) & (set(_api(branch)) - set(_api(base)))
        )
        or "none (the two sides added disjoint symbols)",
    }
    info["ok"] = (
        union == set(merged)
        and not duplicates
        and expected_union == expected_now
        and expected_now == set(merged)
    )
    return info


def _compare_json(branch_doc: Path, regenerated: Path) -> dict[str, Any]:
    want = json.loads(branch_doc.read_text())
    got = json.loads(regenerated.read_text())
    worst = {"rel": 0.0, "path": "", "branch": None, "merged": None}
    count = 0

    versions: dict[str, Any] = {}

    def walk(a: Any, b: Any, path: str) -> None:
        nonlocal count, worst
        if path.endswith("/version") or path == "/version":
            # The branch's committed artifact carries the version stamp of the tree it
            # was generated on (0.20.0); the merged tree reports 0.21.0 from main. The
            # stamp is provenance, not a measured quantity, so record it instead of
            # failing on it -- every number below is still compared exactly.
            versions["baseline"] = a
            versions["regenerated"] = b
            return
        if isinstance(a, dict):
            if set(a) != set(b):
                raise AssertionError(f"key set differs at {path}: {set(a) ^ set(b)}")
            for key in a:
                walk(a[key], b[key], f"{path}/{key}")
        elif isinstance(a, list):
            if len(a) != len(b):
                raise AssertionError(f"length differs at {path}: {len(a)} vs {len(b)}")
            for index, (left, right) in enumerate(zip(a, b)):
                walk(left, right, f"{path}[{index}]")
        elif isinstance(a, bool) or isinstance(a, str):
            if a != b:
                raise AssertionError(f"value differs at {path}: {a!r} vs {b!r}")
        elif isinstance(a, (int, float)):
            count += 1
            scale = max(abs(a), 1e-30)
            rel = abs(a - b) / scale if a != 0 else abs(b)
            if not math.isfinite(rel):
                raise AssertionError(f"non-finite difference at {path}: {a!r} vs {b!r}")
            if rel > worst["rel"]:
                worst = {"rel": rel, "path": path, "branch": a, "merged": b}

    walk(want, got, "")
    return {"numeric_leaves": count, "max_rel_diff": worst["rel"], "worst_path": worst["path"],
            "version_stamp": versions, "ok": worst["rel"] <= NUMERIC_REL_TOL}


def _run(script: Path, output: Path) -> None:
    subprocess.run(
        [sys.executable, str(script), "--output", str(output)],
        cwd=REPO,
        check=True,
        capture_output=True,
        env={**__import__("os").environ, "PYTHONPATH": str(REPO)},
    )


def check_gz_reproduction(full: bool) -> dict[str, Any]:
    if not full:
        return {"skipped": "pass --full"}
    with tempfile.TemporaryDirectory() as tmp:
        _run(REPO / "examples/stability/run_gz.py", Path(tmp))
        result = _compare_json(REPO / "docs/stability/gz_results.json", Path(tmp) / "gz_results.json")
    result["baseline"] = "docs/stability/gz_results.json"
    return result


def check_transom_reproduction(full: bool) -> dict[str, Any]:
    if not full:
        return {"skipped": "pass --full"}
    baseline = REPO / "docs/m2873/outputs/transom_stern"
    with tempfile.TemporaryDirectory() as tmp:
        _run(REPO / "examples/transom_stern/run_example.py", Path(tmp))
        files = {}
        for name in ("design_table.csv", "convergence.csv", "reference_ladder.csv"):
            files[name] = (Path(tmp) / name).read_bytes() == (baseline / name).read_bytes()
    return {"baseline": "docs/m2873/outputs/transom_stern", "byte_identical": files,
            "ok": all(files.values())}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--base", default=None, help="merge base (default: git merge-base main branch)")
    parser.add_argument("--main", default=DEFAULT_MAIN)
    parser.add_argument("--branch", default=DEFAULT_BRANCH)
    parser.add_argument("--full", action="store_true", help="also re-run the two examples")
    parser.add_argument("--json", type=Path, default=None, help="write the evidence here")
    args = parser.parse_args(argv)

    base = args.base or _git_text("merge-base", args.main, args.branch).strip()
    report: dict[str, Any] = {
        "commits": {"base": base, "main": args.main, "branch": args.branch,
                    "merged_head": _git_text("rev-parse", "HEAD").strip()},
        "single_side_files": check_single_side_files(base, args.main, args.branch),
        "public_surface": check_public_surface(base, args.main, args.branch),
        "gz_reproduction": check_gz_reproduction(args.full),
        "transom_reproduction": check_transom_reproduction(args.full),
    }
    keys = ("single_side_files", "public_surface", "gz_reproduction", "transom_reproduction")
    report["all_ok"] = all(report[key].get("ok", True) for key in keys)

    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text)
    print(text)
    return 0 if report["all_ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())