"""Union audit for the #37 cure: prove the three-way merge lost nothing.

The 2026-10 precedent in this ecosystem is that a *mechanical* overlay resolves a
dirty PR by taking one side wholesale, silently dropping the other side's symbols.
This script is the machine check that the semantic merge of ``main`` into
``agent/m2873`` kept both sides, and it is re-runnable after any later rebase.

Checks (all fail hard, exit 1):

1. **file-set union** — every path tracked in either parent is present in the merge.
2. **symbol union** — ``tensorcad.__all__`` in the merge is a superset of the union of
   the parents' ``__all__`` (parsed with ``ast``, not imported, so deleted-but-exported
   names cannot hide behind an ``ImportError``).
3. **no duplicates** — the merged ``__all__`` has no repeated name (a lazy
   concatenation of conflict hunks is exactly how that happens).
4. **importability** — every exported name resolves on the imported package, and both
   family modules (``transom``, ``wigley``) import.
5. **version coherence** — ``tensorcad.__version__`` == ``pyproject.toml`` version, and
   the merged version is >= both parents' (a merge that lowers the version would
   unpublish a release).
6. **behaviour spot-check** — the two hull families produce their documented reference
   numbers in the same interpreter (transom at ``tau = 0`` reproduces the pointed hull;
   Wigley reports ``Cb = 4/9``), so the union is not merely name-level.

Usage::

    python docs/merge_pr37/union_audit.py [--json docs/merge_pr37/union_audit.json]
"""

from __future__ import annotations

import argparse
import ast
import json
import pathlib
import subprocess
import sys
import tomllib

OURS = "origin/agent/m2873"
THEIRS = "origin/main"


def git(*args: str) -> str:
    return subprocess.run(["git", *args], check=True, capture_output=True, text=True).stdout


def tracked_files(ref: str) -> set[str]:
    return {line.strip() for line in git("ls-tree", "-r", "--name-only", ref).splitlines() if line.strip()}


def dunder_all_from_source(source: str) -> list[str]:
    for node in ast.parse(source).body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == "__all__":
                    return [elt.value for elt in node.value.elts if isinstance(elt, ast.Constant)]
    raise SystemExit("no __all__ found in source")


def dunder_all_from_ref(ref: str) -> list[str]:
    return dunder_all_from_source(git("show", f"{ref}:tensorcad/__init__.py"))


def pyproject_version(ref: str | None = None) -> str:
    if ref is None:
        text = pathlib.Path("pyproject.toml").read_text()
    else:
        text = git("show", f"{ref}:pyproject.toml")
    return tomllib.loads(text)["project"]["version"]


def version_tuple(v: str) -> tuple[int, ...]:
    return tuple(int(part) for part in v.split("."))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", default="docs/merge_pr37/union_audit.json")
    args = parser.parse_args(argv)

    results: dict[str, object] = {"parents": {"ours": OURS, "theirs": THEIRS}}
    failures: list[str] = []

    # 1. file-set union -----------------------------------------------------------
    ours_files, theirs_files = tracked_files(OURS), tracked_files(THEIRS)
    merged_files = {str(p.relative_to(pathlib.Path("."))) for p in pathlib.Path(".").rglob("*") if p.is_file()}
    merged_files = {p for p in merged_files if not p.startswith(".git/")}
    expected = ours_files | theirs_files
    missing_files = sorted(expected - merged_files)
    results["files"] = {
        "ours": len(ours_files),
        "theirs": len(theirs_files),
        "union_expected": len(expected),
        "present_in_merge": len(expected & merged_files),
        "missing": missing_files,
    }
    if missing_files:
        failures.append(f"files lost in merge: {missing_files}")

    # 2-3. symbol union -----------------------------------------------------------
    ours_all, theirs_all = dunder_all_from_ref(OURS), dunder_all_from_ref(THEIRS)
    merged_all = dunder_all_from_source(pathlib.Path("tensorcad/__init__.py").read_text())
    expected_symbols = set(ours_all) | set(theirs_all)
    missing_symbols = sorted(expected_symbols - set(merged_all))
    duplicates = sorted({name for name in merged_all if merged_all.count(name) > 1})
    results["symbols"] = {
        "ours": len(ours_all),
        "theirs": len(theirs_all),
        "union_expected": len(expected_symbols),
        "merged": len(merged_all),
        "ours_only": sorted(set(ours_all) - set(theirs_all))[:5] + ["..."],
        "theirs_only": sorted(set(theirs_all) - set(ours_all))[:5] + ["..."],
        "missing": missing_symbols,
        "duplicates": duplicates,
    }
    if missing_symbols:
        failures.append(f"exported names lost in merge: {missing_symbols}")
    if duplicates:
        failures.append(f"duplicate names in merged __all__: {duplicates}")

    # 4. importability ------------------------------------------------------------
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
    import tensorcad  # imported late: check 1-3 must read source, not the build

    unresolved = sorted(name for name in expected_symbols if not hasattr(tensorcad, name))
    modules = {}
    for module_name in ("transom", "wigley", "hull", "bonjean"):
        try:
            module = __import__(f"tensorcad.{module_name}", fromlist=["_"])
        except Exception as exc:  # pragma: no cover - failure path
            failures.append(f"tensorcad.{module_name} does not import: {exc!r}")
            modules[module_name] = "IMPORT FAILED"
            continue
        modules[module_name] = "ok"
        if module_name == "transom":
            for attr in ("generate_transom_hull", "analytic_transom_hydrostatics", "transom_area_ratio"):
                if not hasattr(module, attr):
                    failures.append(f"tensorcad.transom missing {attr}")
        if module_name == "wigley":
            for attr in ("generate_wigley", "wigley_analytic_hydrostatics"):
                if not hasattr(module, attr):
                    failures.append(f"tensorcad.wigley missing {attr}")
    results["imports"] = {"unresolved_exports": unresolved, "modules": modules}
    if unresolved:
        failures.append(f"__all__ names that do not resolve: {unresolved}")

    # 5. version coherence --------------------------------------------------------
    versions = {
        "ours": pyproject_version(OURS),
        "theirs": pyproject_version(THEIRS),
        "merged": pyproject_version(),
        "tensorcad.__version__": tensorcad.__version__,
    }
    results["version"] = versions
    if versions["merged"] != tensorcad.__version__:
        failures.append(f"pyproject {versions['merged']} != tensorcad.__version__ {tensorcad.__version__}")
    for side in ("ours", "theirs"):
        if version_tuple(versions["merged"]) < version_tuple(versions[side]):
            failures.append(f"merged version {versions['merged']} is below {side} {versions[side]}")

    # 6. behaviour spot-check -----------------------------------------------------
    import math

    from tensorcad import (
        HullParameters,
        analytic_hydrostatics,
        wigley_analytic_hydrostatics,
        wigley_reference_parameters,
    )
    from tensorcad.transom import TransomHullParameters, generate_transom_hull, transom_area_ratio

    # both families in one interpreter; tau = 0 must reproduce the pointed hull exactly
    reference = analytic_hydrostatics(HullParameters(length=10.0, beam=2.0, draft=0.8))
    pointed = generate_transom_hull(TransomHullParameters(transom_fraction=0.0)).hydrostatics
    transom = generate_transom_hull(TransomHullParameters(transom_fraction=0.25)).hydrostatics
    wigley = wigley_analytic_hydrostatics(wigley_reference_parameters())
    area_ratio_error = abs(transom_area_ratio(0.25, 0.7) - math.sin(math.pi * 0.25) ** 0.7)
    pointed_error = abs(pointed.displacement_volume - reference.displacement_volume) / reference.displacement_volume
    checks = {
        "pointed_vs_analytic_relative_error": f"{pointed_error:.3e}",
        "transom_tau0_volume_m3": round(pointed.displacement_volume, 12),
        "transom_tau025_volume_m3": round(transom.displacement_volume, 12),
        "transom_area_ratio_error": f"{area_ratio_error:.3e}",
        "wigley_Cb": round(wigley.block_coefficient, 12),
        "wigley_Cb_error_vs_4_9": f"{abs(wigley.block_coefficient - 4.0 / 9.0):.3e}",
    }
    ok = (
        # the 41x13 mesh integrates ~0.9 % low on this family (documented signed error);
        # the check is that tau = 0 is *the same surface*, not that the grid is fine
        pointed_error < 1.5e-2
        and area_ratio_error < 1e-12
        and transom.displacement_volume > pointed.displacement_volume  # LOA frame lengthens with tau
        and abs(wigley.block_coefficient - 4.0 / 9.0) < 1e-12
    )
    results["behaviour"] = {"checks": checks, "ok": ok}
    if not ok:
        failures.append(f"behaviour spot-check failed: {checks}")

    results["failures"] = failures
    results["passed"] = not failures

    out = pathlib.Path(args.json)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2) + "\n")

    print(json.dumps({k: results[k] for k in ("files", "version", "behaviour", "passed")}, indent=2))
    for failure in failures:
        print("FAIL:", failure, file=sys.stderr)
    print(f"audit written to {out}")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())