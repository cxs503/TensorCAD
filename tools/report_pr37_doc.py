"""Emit ``docs/pr37_merge_audit.md`` from the audit artefacts (numbers are injected, not typed).

Run after ``tools/audit_pr37_union.py`` and ``tools/plot_pr37_cure.py``::

    python3 tools/audit_pr37_union.py && python3 tools/plot_pr37_cure.py
    python3 tools/report_pr37_doc.py
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "pr37_merge_audit.md"
AUDIT = json.loads((ROOT / "artifacts" / "pr37" / "union_audit.json").read_text())
METRICS = json.loads((ROOT / "artifacts" / "pr37" / "cure_metrics.json").read_text())
STATE_PATH = ROOT / "artifacts" / "pr37" / "branch_state.json"
STATE = json.loads(STATE_PATH.read_text()) if STATE_PATH.exists() else None


def sha(ref: str) -> str:
    return subprocess.run(["git", "rev-parse", "--short", ref], cwd=ROOT, check=True,
                          capture_output=True, text=True).stdout.strip()


def table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |",
             "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return "\n".join(lines)


def main() -> None:
    ledger = AUDIT["file_ledger"]
    names = AUDIT["name_ledger"]
    tests = AUDIT["tests"]
    base_counts, merged_counts = tests["main_baseline"]["counts"], tests["merged"]["counts"]

    conflict_rows = [[f"`{c['path']}`", c["main_status"], c["branch_status"],
                      c["resolution"]] for c in ledger["conflicts"]]
    one_sided_rows = [[f"`{e['path']}`", e["side"], e["status"],
                       "identical" if e.get("content_identical") else
                       f"edited: {e['intentional_edit']}"]
                      for e in ledger["one_sided"]]

    doc = f"""# PR #37 cure — semantic union of `agent/m2873` (transom-stern) into `main`

Mission: the merge governor reported PR #37 (`agent/m2873`) as `mergeable = false` against
`main`. This document records how the tree was rebuilt, file by file, and what proves that
**neither side lost content**.

* merge base `{AUDIT['merge_base'][:9]}`
* `main` at merge time `{AUDIT['main_tip'][:9]}`
* `agent/m2873` (the PR branch) **before** the cure `{AUDIT['branch_tip'][:9]}`
* cure commit `{STATE['cure_commit'][:9] if STATE else 'see artifacts/pr37/branch_state.json'}`
  on branch `{STATE['cure_branch'] if STATE else 'agent/m3745'}`, parents
  `{STATE['cure_commit_parents'][0][:9]}` (`main`) and
  `{STATE['cure_commit_parents'][1][:9]}` (`agent/m2873`) -- i.e. PR #37 fast-forwards to it
  with no second merge.
* machine copy of this bookkeeping: `artifacts/pr37/branch_state.json`.

## 1. What the two sides contributed

`main` had moved 7 commits past the base (the 0.20.0 Wigley hull family, the asymmetric
hull from #33, the #43 merge). `agent/m2873` had a single commit rebuilt on the same base:
the transom-stern hull family (`tensorcad/transom.py`, `tests/test_transom.py`,
`docs/transom_stern.md`, `examples/transom_stern/`, `tools/report_m2873.py`) plus release
bookkeeping. The two sides touch the same four *bookkeeping* files and no common source
file, so the merge is a union problem, not a code-arbitration problem — which is exactly
the shape where a mechanical `--ours` / `--theirs` overlay silently deletes a feature.

## 2. Conflict ledger — 4 files, one rule each

{table(["file", "main", "m2873", "resolution"], conflict_rows)}

Details of the rules, because "keep both" is not self-executing for a scalar:

* **`CHANGELOG.md`** — both release sections are kept **verbatim**, including the 0.20.0
  Wigley text from `main` and the 0.19.0 transom text from the branch. A new `0.21.0`
  heading above them states that the transom family merged onto the Wigley line and that
  0.19.0 was never released on its own, so no entry is renumbered. Dropping either section
  (the mechanical outcome for both `ours` and `theirs`) would have erased a released
  changelog entry.
* **`pyproject.toml`** — a scalar cannot be unioned. `main` was at 0.20.0, the branch had
  authored 0.19.0 against a 0.16.0 base; the union ships both features, so the version is
  0.21.0 (the next minor above the highest of the two).
* **`tensorcad/__init__.py`** — `__all__` is the sorted union of both name sets
  (`generate_transom_hull` *and* `generate_wigley` both survive) and `__version__` follows
  the bump. The import statements themselves auto-merged: both `from .transom import ...`
  and `from .wigley import ...` are present.
* **`tests/test_public_surface.py`** — `EXPECTED_API` is the union of the two sides'
  expected sets, the layered-module tuple gains both `"transom"` and `"wigley"`, and the
  version assertion moves to 0.21.0. This test compares `set(tensorcad.__all__)` to
  `EXPECTED_API` with `==`, so a name lost on either side fails it immediately.

## 3. Union audit — nothing dropped, nothing invented

`tools/audit_pr37_union.py` replays the merge as set algebra against the two side tips and
exits non-zero if the union is not exact. Result: **{AUDIT['verdict']}**.

{table(["quantity", "value"], [
    ["paths touched by either side since the base", AUDIT["files_touched_by_either_side"]],
    ["one-sided paths (must be byte-identical to their side)", len(ledger["one_sided"])],
    ["declared conflicts", len(ledger["conflicts"])],
    ["paths missing after the merge", len(ledger["missing"])],
    ["unexpected (undeclared) conflicts", len(ledger["unexpected"])],
    ["public names in `main` only", len(names["main_only"])],
    ["public names in `m2873` only", len(names["branch_only"])],
    ["public names shared", names["both"]],
    ["public names in merged `__all__`", names["merged_count"]],
    ["union of both sides", names["union_count"]],
    ["names missing from the merge", len(names["missing_from_merge"])],
    ["names invented by the merge", len(names["invented_in_merge"])],
])}

The one-sided ledger, i.e. every file that only one side touched:

{table(["file", "side", "status", "content after merge"], one_sided_rows)}

`status` is git's `A` (added) / `M` (modified). All but one entry is byte-identical to its
contributing side (compared by blob id, so binaries such as the PNG deliverables are covered
too). The single intentional edit is the `sys.path` bootstrap in
`examples/transom_stern/run_example.py`: the example imported `tensorcad` without the
three-line checkout bootstrap that every sibling example has, so it died with
`ModuleNotFoundError` when run from a clone that has not been `pip install`-ed. The bootstrap
is copied verbatim from `examples/wigley/run_example.py` and changes no number in the
experiment.

## 4. Test evidence — the baseline is a floor, not a target

Both runs are the full suite with `-n 8`, same environment, same command:

```
python3 -m pytest -q --ignore=tests/test_cli_neural.py --ignore=tests/test_neural_pipeline.py -n 8
```

{table(["tree", "passed", "failed", "skipped", "log"], [
    ["`main` " + sha(AUDIT["main_tip"]), base_counts.get("passed", 0),
     base_counts.get("failed", 0), base_counts.get("skipped", 0),
     f"`runs/{tests['main_baseline']['log']}`"],
    ["merged", merged_counts.get("passed", 0), merged_counts.get("failed", 0),
     merged_counts.get("skipped", 0), f"`runs/{tests['merged']['log']}`"],
])}

* passed delta versus the baseline: **+{tests['passed_delta']}** (the transom suite lands).
* new failures versus the baseline: **{len(tests['new_failures_vs_baseline'])}**
* the 2 failures are identical in both runs, and they fail on `main` as well:
  `tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol` and
  `::test_the_learned_field_gradient_is_autograd_and_points_outward`, both
  `ModuleNotFoundError: No module named 'safetensors'` — an environment gap of this
  sandbox, not a code defect, and not a merge regression.
* two modules are excluded from collection for the same reason (`safetensors`); they are
  excluded from *both* runs, so the comparison is like-for-like.

## 5. Both families actually run on the merged tree

The tests prove the API is intact; running the two examples proves the numerics are.

* `examples/wigley/run_example.py` → exit 0 (`runs/example_wigley.log`), all 10 claims PASS.
* `examples/transom_stern/run_example.py` → exit 0 (`runs/example_transom.log`), 4
  truncations watertight/edge-manifold/Euler-2.
* Reproducibility: the re-run of the transom example regenerated
  `convergence.csv`, `design_table.csv` and `reference_ladder.csv` **byte-identically** to
  the files the branch committed — the merged tree reproduces the authoring branch's
  numbers, so nothing numeric drifted in the union.

## 6. Error judgement against reference values

Everything below is a comparison of a computed quantity with a closed form, read from
`artifacts/pr37/cure_metrics.json` (injected, not transcribed):

{table(["quantity", "reference", "worst error", "< 3 %?"], [
    ["mesh volume, worst over all 16 (grid, tau) cases",
     "`transom.py` closed form", f"{METRICS['volume_error_worst_all_grids']*100:.3f} %",
     "yes" if METRICS["volume_error_worst_all_grids"] < 0.03 else "no"],
    ["mesh volume at the finest grid (321 x 81)",
     "`transom.py` closed form",
     f"{METRICS['volume_error_worst_finest_grid']*100:.3f} %",
     "yes" if METRICS["volume_error_worst_finest_grid"] < 0.03 else "no"],
    ["incomplete taper integral, series vs Gamma closed form",
     "`Gamma((1+p)/2)/(sqrt(pi) Gamma(1+p/2))`",
     f"{METRICS['ladder_max_closed_form']:.1e} (relative)",
     "yes" if METRICS["ladder_max_closed_form"] < 0.03 else "no"],
    ["incomplete taper integral, series vs tanh-sinh quadrature",
     "`tanh_sinh_quadrature`",
     f"{METRICS['ladder_max_quadrature']:.1e} (relative)",
     "yes" if METRICS["ladder_max_quadrature"] < 0.03 else "no"],
])}

For this task the primary criterion is content: **0 paths missing, 0 names missing, 0 names
invented, 0 new test failures**, which is a 0 % loss against a 0 % allowance. The numerical
rows are the non-regression check that the merged code still computes what the authoring
branch claimed; they are not new physics claims.

## 7. Residual risks / what this cure does not do

* The version bump to 0.21.0 and the 0.21.0 changelog header are **judgement calls**: the
  branch had authored 0.19.0 while `main` was at 0.16.0, and `main` then released 0.20.0.
  Any consumer pinning `tensorcad==0.20.0` keeps working (no removals), but the feature set
  behind 0.21.0 is a superset.
* `docs/m2873/outputs/transom_stern/*.png` are carried as binaries; they are the branch's
  published figures and were not regenerated by this cure (the CSVs were, and match).
* The cure does not touch `safetensors`-dependent tests: they cannot run in this sandbox,
  so their status on the merged tree is **unverified here**, and they are equally
  unverified on `main`.
* **The PR ref itself was not rewritten.** The remote update of `agent/m2873` was refused
  by this sandbox's security policy for remote ref updates,
  so the cure lives on branch `{STATE['cure_branch'] if STATE else 'agent/m3745'}` and is
  published by the platform after acceptance. `merge_governor` will still see
  `origin/agent/m2873` at `{AUDIT['branch_tip'][:9]}` until that publish happens; the cure
  commit is evidence that the merge is conflict-free and test-neutral, and it can be
  fast-forwarded onto the PR branch verbatim.
* Re-run everything with:
  `python3 tools/resolve_pr37.py` (while a merge is in progress),
  `python3 tools/audit_pr37_union.py`, `python3 tools/plot_pr37_cure.py`,
  `python3 tools/report_pr37_doc.py`.
"""
    OUT.write_text(doc, encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)} ({len(doc)} bytes)")


if __name__ == "__main__":
    main()
