"""Board audit m20064 -- CASE-CAD-CURVATURE-SPHERE (mission 20002), closure verification.

Nothing here trusts the audited mission's own report: every claim in
``evidence.json`` is produced by re-reading the audited work tree, re-running the
committed runner in a clean process, diffing the fresh numbers against the
committed ``result.json``, querying the benchmark ledgers, and re-running the
pinning test module.

    python docs/board_audit/m20064/reproduce_and_evidence.py

Outputs (all inside this worktree):

* ``runs/m20064_sphere_curvature/{result.json,gates_verified.json,figures/,repro.log}``
  -- the fresh, independent measurement;
* ``docs/board_audit/m20064/evidence.json`` -- the machine-checked evidence bundle;
* ``docs/board_audit/m20064/gates_table.md`` -- reference vs measured per gate.

The audited delivery lives on ``agent/m20002`` (commit ``4558ef6``); this branch
does not carry it, so the runner is executed *in the audited tree* and bound to
the audited commit.  Its identity is pinned by sha256 of every input file below.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

AUDITED_TREE = Path("/root/tp-agent/work/m20002")
AUDITED_COMMIT = "4558ef6d3d18e91090d777ab05987528a85fe72e"
CASE_REL = Path("benchmarks/verified/sphere_mean_curvature")
CASE = AUDITED_TREE / CASE_REL
ANCHOR = "CASE-CAD-CURVATURE-SPHERE"

HERE = Path(__file__).resolve().parent
WORKTREE = HERE.parents[2]
RUNS = WORKTREE / "runs" / "m20064_sphere_curvature"
DB = os.environ.get("TP_PG_DB", "tensorplatform")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def run(argv: list[str], *, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run(argv, cwd=cwd, capture_output=True, text=True)


def git(*args: str, cwd: Path = AUDITED_TREE) -> str:
    return run(["git", *args], cwd=cwd).stdout.strip()


def psql(query: str) -> list[list[str]]:
    """Read-only query against the platform ledger (psql -At, pipe separated)."""
    proc = run(["psql", "-U", os.environ.get("PGUSER", "postgres"), "-d", DB, "-Atc", query])
    if proc.returncode != 0:
        raise RuntimeError(f"psql failed: {proc.stderr.strip()}")
    return [line.split("|") for line in proc.stdout.splitlines() if line.strip()]


def numeric_diff(committed: object, fresh: object, *, skip: set[str], path: str = "") -> list[dict]:
    """Every leaf that differs between the committed and the fresh result."""
    out: list[dict] = []
    if isinstance(committed, dict) and isinstance(fresh, dict):
        for key in committed:
            if key in skip:
                continue
            if key not in fresh:
                out.append({"path": f"{path}/{key}", "committed": "present", "fresh": "absent"})
                continue
            out += numeric_diff(committed[key], fresh[key], skip=skip, path=f"{path}/{key}")
    elif isinstance(committed, list) and isinstance(fresh, list):
        if len(committed) != len(fresh):
            out.append({"path": path, "committed": f"len {len(committed)}", "fresh": f"len {len(fresh)}"})
        else:
            for index, (a, b) in enumerate(zip(committed, fresh)):
                out += numeric_diff(a, b, skip=skip, path=f"{path}[{index}]")
    elif isinstance(committed, float) and isinstance(fresh, float):
        if abs(committed - fresh) > 1e-12 * max(1.0, abs(committed)):
            out.append({"path": path, "committed": committed, "fresh": fresh})
    elif committed != fresh:
        out.append({"path": path, "committed": committed, "fresh": fresh})
    return out


def gate_rows(gates: dict) -> list[dict]:
    rows = []
    for name, gate in gates["gates"].items():
        rows.append(
            {
                "gate": name.split(":")[0],
                "criterion": name,
                "reference": gate["reference"],
                "measured": gate["measured"],
                "error_percent": gate["error_percent"],
                "tolerance": gate["tolerance"],
                "verdict": gate["verdict"],
            }
        )
    for name, gate in (gates.get("diagnostics") or {}).items():
        rows.append(
            {
                "gate": name.split(":")[0],
                "criterion": name,
                "reference": gate.get("reference"),
                "measured": gate.get("measured"),
                "error_percent": gate.get("error_percent"),
                "tolerance": gate.get("tolerance"),
                "verdict": f"diagnostic-{gate['verdict']}",
            }
        )
    return rows


def tolerance_check(row: dict) -> dict:
    """Is the gate's own got-vs-ref inside its own tolerance?  Two conventions live
    in this file: relative-error gates (tolerance is a fraction) and absolute gates
    (G3, order vs 2 +- 0.5).  Both are re-derived here instead of trusted."""
    if row["measured"] is None or row["reference"] is None:
        return {"pass": None, "relative_error": None, "note": "not a numeric gate (diagnostic without a value)"}
    if row["gate"] == "G3":
        delta = abs(row["measured"] - row["reference"])
        return {"pass": delta <= row["tolerance"], "absolute_delta": delta, "tolerance": row["tolerance"]}
    if row["reference"] == 0.0:
        # G2 convention: the reference is "zero error" and the measured value *is* the
        # error, so the comparison is absolute, not relative.
        return {"pass": abs(row["measured"]) <= row["tolerance"], "absolute": abs(row["measured"]),
                "tolerance": row["tolerance"], "convention": "absolute (reference = 0 error)"}
    relative = abs(row["measured"] - row["reference"]) / abs(row["reference"])
    return {"pass": relative <= row["tolerance"], "relative_error": relative, "tolerance": row["tolerance"],
            "convention": "relative"}


def main() -> int:
    RUNS.mkdir(parents=True, exist_ok=True)
    facts: dict[str, object] = {}
    facts["audited_mission"] = 20002
    facts["audited_tree"] = str(AUDITED_TREE)
    facts["generated_at"] = datetime.now(timezone.utc).isoformat()

    # -- 1. identity of the audited code --------------------------------------
    head = git("rev-parse", "HEAD")
    facts["audited_head"] = head
    facts["head_matches_expected"] = head == AUDITED_COMMIT
    facts["dirty"] = git("status", "--porcelain") != ""
    facts["on_main"] = run(["git", "merge-base", "--is-ancestor", head, "origin/main"], cwd=AUDITED_TREE).returncode == 0
    facts["branches_containing_head"] = git("branch", "--contains", head).replace("*", "").split()
    facts["commit_subject"] = git("log", "-1", "--format=%s", head)
    facts["commit_files"] = git("show", "--stat=1000", "--format=", "--name-only", head).split()

    # -- 2. do the three deliverables exist? ----------------------------------
    artifacts = {}
    for rel in [
        CASE_REL / "run_benchmark.py",
        CASE_REL / "result.json",
        CASE_REL / "gates_verified.json",
        CASE_REL / "README.md",
        Path("docs/sphere_mean_curvature.md"),
        Path("tests/test_case_cad_curvature_sphere.py"),
    ] + sorted(str(p.relative_to(AUDITED_TREE)) for p in (CASE / "figures").glob("*.png")):
        path = AUDITED_TREE / rel
        artifacts[str(rel)] = {"exists": path.is_file(), "bytes": path.stat().st_size if path.is_file() else 0,
                               "sha256": sha256(path) if path.is_file() else None}
    facts["artifacts"] = artifacts
    facts["figures_present"] = sorted(p.name for p in (CASE / "figures").glob("*.png"))

    # -- 3. the benchmark ledger is the sole judge: what does it say? ---------
    def count(query: str) -> str:
        return psql(query)[0][0]

    facts["ledger"] = {
        "tensor_benchmark_cases_rows_for_anchor": count(
            f"select count(*) from tensor_benchmark_cases where anchor_slug = '{ANCHOR}'"),
        "tensor_benchmark_cases_rows_cad_any": count(
            "select count(*) from tensor_benchmark_cases where anchor_slug ilike '%CAD%'"),
        "tensor_benchmark_cases_total": count("select count(*) from tensor_benchmark_cases"),
        "tensorbench_results_rows_metric_curvature": count(
            "select count(*) from tensorbench_results where metric ilike '%curvature%'"),
        "tensorbench_results_rows_total": count("select count(*) from tensorbench_results"),
        "tensorbench_evaluations_rows_total": count("select count(*) from tensorbench_evaluations"),
        "cad_anchor_rows": psql(
            "select anchor_slug, case_name, status, coalesce(passed_at::text,'-') "
            "from tensor_benchmark_cases where repo = 'TensorCAD' order by id"),
    }

    # -- 4. independent re-run -------------------------------------------------
    committed = json.loads((CASE / "result.json").read_text())
    out_dir = RUNS
    repro = run([sys.executable, str(CASE_REL / "run_benchmark.py"), "--out", str(out_dir)], cwd=AUDITED_TREE)
    (out_dir / "repro.log").write_text(
        f"$ cd {AUDITED_TREE} && {sys.executable} {CASE_REL}/run_benchmark.py --out {out_dir}\n"
        f"commit {head}\nreturncode {repro.returncode}\n--- stdout ---\n{repro.stdout}\n--- stderr ---\n{repro.stderr}\n"
    )
    facts["repro_returncode"] = repro.returncode
    facts["repro_stdout"] = repro.stdout.strip().splitlines()

    fresh = json.loads((out_dir / "result.json").read_text())
    facts["fresh_result_sha256"] = sha256(out_dir / "result.json")
    facts["committed_result_sha256"] = sha256(CASE / "result.json")
    mismatches = numeric_diff(committed, fresh, skip={"seconds", "environment", "figures", "generated_at"})
    facts["numeric_mismatches_vs_committed"] = mismatches
    facts["numeric_leaves_compared"] = None  # filled below
    facts["headline"] = fresh["headline"]
    facts["fresh_seconds"] = fresh["seconds"]
    facts["environment"] = fresh["environment"]

    # count leaves for the "N leaves identical" statement
    def count_leaves(node: object, skip: set[str]) -> int:
        if isinstance(node, dict):
            return sum(count_leaves(v, skip) for k, v in node.items() if k not in skip)
        if isinstance(node, list):
            return sum(count_leaves(v, skip) for v in node)
        return 1

    facts["numeric_leaves_compared"] = count_leaves(committed, {"seconds", "environment", "figures", "generated_at"})
    facts["identical_leaves"] = facts["numeric_leaves_compared"] - sum(1 for _ in mismatches)

    # -- 5. gate table, re-derived --------------------------------------------
    gates = json.loads((out_dir / "gates_verified.json").read_text())
    rows = gate_rows(gates)
    for row in rows:
        row["independent_check"] = tolerance_check(row)
    facts["gates"] = rows
    facts["gates_all_pass_claim"] = gates["all_gates_pass"]
    facts["gates_failing_claim"] = gates["failing"]
    facts["gates_numeric_fail_after_recheck"] = [
        row["gate"] for row in rows if row["independent_check"]["pass"] is False
    ]
    facts["diagnostics"] = gates.get("diagnostics")

    lines = ["| gate | criterion | reference | measured | error | tolerance | mission verdict | audit re-check |",
             "|---|---|---|---|---|---|---|---|"]
    for row in rows:
        err = "—" if row["error_percent"] is None else f"{row['error_percent']:+.4g} %"
        ref = "—" if row["reference"] is None else f"{row['reference']:.10g}"
        got = "—" if row["measured"] is None else f"{row['measured']:.10g}"
        check = row["independent_check"]
        verdict = "pass" if check["pass"] else ("fail" if check["pass"] is False else "n/a")
        lines.append(
            f"| {row['gate']} | {row['criterion'].split(': ', 1)[-1]} | {ref} | {got} | {err} | "
            f"{row['tolerance']} | {row['verdict']} | {verdict} |")
    (HERE / "gates_table.md").write_text("\n".join(lines) + "\n")

    # -- 6. the pinned test module, re-run ------------------------------------
    test = run([sys.executable, "-m", "pytest", "tests/test_case_cad_curvature_sphere.py", "-q"], cwd=AUDITED_TREE)
    tail = [line for line in test.stdout.splitlines() if line.strip()][-3:]
    facts["tests"] = {"returncode": test.returncode, "tail": tail,
                      "summary": next((line for line in reversed(test.stdout.splitlines()) if "passed" in line or "failed" in line), "")}
    facts["test_names"] = re.findall(r"^def (test_\w+)", (AUDITED_TREE / "tests/test_case_cad_curvature_sphere.py").read_text(), re.M)

    # -- 7. integration callers (L3 evidence), grepped not asserted -----------
    callers = []
    for pattern in ["vertex_mean_curvature", "measure("]:
        proc = run(["grep", "-rn", "--include=*.py", pattern,
                    "tensorcad", "tools", "examples", "scripts"], cwd=AUDITED_TREE)
        for line in proc.stdout.splitlines():
            if line.startswith("tensorcad/curvature.py:"):
                continue  # definition site
            callers.append(line)
    facts["integration_callers_repo"] = callers
    facts["integration_callers_platform"] = run(
        ["bash", "-lc",
         "grep -rn 'vertex_mean_curvature\\|tensorcad.curvature' --include=*.py /root/tensor-ecosystem/TensorPlatform/backend "
         "/root/tensor-ecosystem/TensorPlatform/api 2>/dev/null | head -5"],
    ).stdout.strip().splitlines()
    facts["console_scripts"] = [
        line for line in (AUDITED_TREE / "pyproject.toml").read_text().splitlines()
        if re.match(r"^tensorcad[-a-z0-9]* = ", line)
    ]

    # -- 8. report shape ------------------------------------------------------
    doc = (AUDITED_TREE / "docs/sphere_mean_curvature.md").read_text()
    facts["doc_sections"] = re.findall(r"^## (.+)$", doc, re.M)
    facts["doc_bytes"] = len(doc.encode())
    facts["result_schema"] = committed.get("schema")
    facts["result_sections"] = sorted(committed.keys())

    (HERE / "evidence.json").write_text(json.dumps(facts, indent=2, ensure_ascii=False, default=str) + "\n")

    print(json.dumps({k: facts[k] for k in (
        "audited_head", "head_matches_expected", "on_main", "branches_containing_head",
        "repro_returncode", "numeric_leaves_compared", "identical_leaves",
        "gates_numeric_fail_after_recheck", "tests", "ledger")}, indent=2, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())