#!/usr/bin/env python3
"""Board-audit m20544 of mission m20512 (anchor XVAL-GEOM-VOLCROSS, repo TensorCAD).

Independence contract
---------------------
The audited worktree /srv/tensor/agent/work/m20512 is treated read-only: everything
below either (a) reads its committed artifacts, (b) recomputes the published
numbers with arithmetic written here (no import of ``tensorcad``), or (c) re-runs
the subject's own declared ``repro_cmd`` in a disposable copy under /tmp and
records the exit code and stdout.  Nothing in the subject tree is written.

What this script answers
------------------------
1. do the claimed run artifacts exist, and what are their sha256?
2. does the platform referee table ``tensor_benchmark_cases`` carry a case for
   this anchor (the audit mandate says that table is the sole referee)?
3. do the artifact numbers reproduce under *independent* arithmetic
   (analytic Wigley volume, per-rung deviations, measured deficit order,
   Richardson limit, KCS block coefficient / calibration-order gap)?
4. do the in-repo benchmark gates equal the gates recomputed from ``result.json``?
5. does the posted mission report (``tensor_reports``) match the artifacts, or
   does it carry hand-copied / truncated numbers?
6. is ``repro_cmd`` executable and does it return the same gate values?

Run: python3 audit/m20544_audit.py            (writes audit/m20544_evidence.json)
"""
from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import psycopg

SUBJ = Path("/srv/tensor/agent/work/m20512")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DSN = "host=/var/run/postgresql port=5432 dbname=tensorplatform user=postgres"
ANCHOR = "XVAL-GEOM-VOLCROSS"
MISSION = 20512
AUDIT_MISSION = 20544

ART = SUBJ / "artifacts/xval_geom_volcross"
BENCH = SUBJ / "benchmarks/verified/volume_cross_engine"
KCS = {"length": 230.0, "beam": 32.2, "draft": 10.8, "cb_published": 0.6505,
       "volume_published": 52030.0}
WIGLEY = {"length": 10.0, "beam": 1.0, "draft": 0.625}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def q(sql: str, params: tuple = ()) -> list[tuple]:
    with psycopg.connect(DSN) as cx, cx.cursor() as cur:
        cur.execute(sql, params)
        return cur.fetchall()


# ---------------------------------------------------------------------------
# 1. artifact inventory
# ---------------------------------------------------------------------------
def inventory() -> dict:
    files = [
        "tensorcad/volume_xval.py", "tensorcad/volume_xval_cli.py",
        "tensorcad/fleet_gate.py", "tensorcad/cli.py", "tensorcad/__init__.py",
        "tests/test_volume_xval.py", "tests/test_volume_xval_entry.py",
        "artifacts/xval_geom_volcross/results.json",
        "artifacts/xval_geom_volcross/results.md",
        "artifacts/xval_geom_volcross/figures/convergence.png",
        "artifacts/xval_geom_volcross/figures/deviation_vs_reference.png",
        "artifacts/xval_geom_volcross/figures/hull_geometry.png",
        "benchmarks/verified/volume_cross_engine/run_benchmark.py",
        "benchmarks/verified/volume_cross_engine/result.json",
        "benchmarks/verified/volume_cross_engine/gates_verified.json",
        "benchmarks/verified/volume_cross_engine/README.md",
        "benchmarks/verified/volume_cross_engine/figures/volume_cross_engine.png",
        "tools/volume_xval_report.py", "docs/volume_xval.md",
    ]
    out = {}
    for rel in files:
        p = SUBJ / rel
        out[rel] = {"exists": p.is_file(), "bytes": p.stat().st_size if p.is_file() else 0,
                    "sha256": sha(p) if p.is_file() else None}
    head = subprocess.run(["git", "rev-parse", "-s", "HEAD"], cwd=SUBJ,
                          capture_output=True, text=True).stdout.strip()
    status = subprocess.run(["git", "status", "--porcelain"], cwd=SUBJ,
                            capture_output=True, text=True).stdout.strip()
    out["_worktree"] = {"head_commit": head, "dirty_entries": len(status.splitlines()),
                        "branch": subprocess.run(
                            ["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=SUBJ,
                            capture_output=True, text=True).stdout.strip()}
    return out


# ---------------------------------------------------------------------------
# 2. referee table
# ---------------------------------------------------------------------------
def referee() -> dict:
    total, anchors = q("SELECT count(*), count(DISTINCT anchor_slug) FROM tensor_benchmark_cases")[0]
    rows = q("SELECT case_name, metric, tolerance, repo, status FROM tensor_benchmark_cases "
             "WHERE anchor_slug=%s ORDER BY case_name", (ANCHOR,))
    related = q("SELECT anchor_slug, case_name, status FROM tensor_benchmark_cases "
                "WHERE repo='TensorCAD' ORDER BY anchor_slug, case_name")
    prev = q("SELECT scope, level, gate_state, gate_level, verdict_by, created_at "
             "FROM tensor_maturity_verdicts WHERE scope LIKE %s ORDER BY created_at", ("%VOLCROSS%",))
    return {
        "anchor": ANCHOR,
        "table_cases": total,
        "table_distinct_anchors": anchors,
        "rows_for_anchor": len(rows),
        "rows": rows,
        "tensorcad_anchor_cases": related,
        "previous_verdicts_same_anchor": [
            {"scope": r[0], "level": r[1], "gate_state": r[2], "gate_level": r[3],
             "verdict_by": r[4], "created_at": str(r[5])} for r in prev],
    }


# ---------------------------------------------------------------------------
# 3. independent recomputation from the committed artifacts
# ---------------------------------------------------------------------------
def fit_order_own(samples: list[tuple[float, float]]) -> float | None:
    """Own log-log least squares; no tensorcad import."""
    pts = [(h, d) for h, d in samples if h > 0 and d > 0]
    if len(pts) < 2:
        return None
    xs = [math.log(h) for h, _ in pts]
    ys = [math.log(d) for _, d in pts]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    den = sum((x - mx) ** 2 for x in xs)
    if den == 0:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / den


def richardson_own(fine: float, coarse: float, p: float) -> float:
    return fine + (fine - coarse) / (2 ** p - 1)


def recompute() -> dict:
    payload = json.loads((ART / "results.json").read_text())
    cases = {c["key"]: c for c in payload["cases"]}
    out: dict = {"cases": {}}

    # --- Wigley: analytic reference, own arithmetic -------------------------
    wig = cases["wigley"]
    exact = 4.0 * WIGLEY["length"] * WIGLEY["beam"] * WIGLEY["draft"] / 9.0
    rows = wig["rows"]
    dev = [{"stations": r["stations"],
            "engine_a": r["engine_a_trimesh_m3"],
            "engine_b": r["engine_b_numpy_stl_m3"],
            "my_vs_ref": (r["engine_a_trimesh_m3"] - exact) / exact,
            "artifact_vs_ref": r["engine_a_vs_reference"],
            "my_engine_spread": abs(r["engine_a_trimesh_m3"] - r["engine_b_numpy_stl_m3"])
            / abs(r["engine_a_trimesh_m3"]),
            "artifact_engine_spread": r.get("engine_agreement")} for r in rows]
    for d in dev:
        d["spread_delta_abs"] = abs(d["my_engine_spread"] - (d["artifact_engine_spread"] or 0.0))
    ladder = [(r.get("h_m", r.get("h")), exact - r["engine_a_trimesh_m3"]) for r in rows]
    p = fit_order_own(ladder)
    fine, coarse = rows[-1]["engine_a_trimesh_m3"], rows[-2]["engine_a_trimesh_m3"]
    out["cases"]["wigley"] = {
        "analytic_volume_own_4LBT_over_9": exact,
        "artifact_reference_volume": wig["reference_volume_m3"],
        "reference_matches_own": abs(wig["reference_volume_m3"] - exact) <= 1e-12 * exact,
        "rows": dev,
        "max_abs_delta_between_my_deviation_and_artifact": max(
            abs(d["my_vs_ref"] - d["artifact_vs_ref"]) for d in dev),
        "my_measured_order": p,
        "artifact_measured_order": wig.get("attribution", {}).get("observed_order"),
        "my_richardson_limit": richardson_own(fine, coarse, p),
        "my_richardson_vs_reference": (richardson_own(fine, coarse, p) - exact) / exact,
        "finest_engine_a_vs_reference_own": (fine - exact) / exact,
    }

    # --- KCS full scale: geometry must be checked against Cb / L / B / T ----
    kcs = cases["kcs-fullscale"]
    lbt = KCS["length"] * KCS["beam"] * KCS["draft"]
    v_pipeline = kcs["pipeline_volume_m3"]
    out["cases"]["kcs-fullscale"] = {
        "LBT_m3_own": lbt,
        "artifact_LBT_m3": kcs["principal_particulars"]["LBT_m3"],
        "pipeline_volume_m3": v_pipeline,
        "cb_from_pipeline_volume_own": v_pipeline / lbt,
        "cb_published": KCS["cb_published"],
        "cb_excess_relative": v_pipeline / lbt / KCS["cb_published"] - 1.0,
        "vs_published_volume_own": (v_pipeline - KCS["volume_published"]) / KCS["volume_published"],
        "artifact_pipeline_vs_reference": kcs["rows"][-1]["pipeline_vs_reference"],
        "finest_engine_a_vs_published_own": (
            kcs["rows"][-1]["engine_a_trimesh_m3"] - KCS["volume_published"]) / KCS["volume_published"],
        "finest_engine_a_vs_published_artifact": kcs["rows"][-1]["engine_a_vs_reference"],
        "topology_finest": kcs["rows"][-1]["topology"],
    }

    # --- reference-free rows must say so -----------------------------------
    out["cases"]["kcs-modelscale"] = {
        "reference_kind": cases["kcs-modelscale"]["reference_kind"],
        "reference_volume": cases["kcs-modelscale"].get("reference_volume_m3"),
        "cb_from_pipeline_own": cases["kcs-modelscale"]["pipeline_volume_m3"] / (3.0 * 0.4 * 0.15),
    }
    battery = payload.get("integrator_battery") or payload.get("battery")
    out["integrator_battery_present"] = battery is not None
    if isinstance(battery, dict):
        analytic = battery["reference_volume_m3"]
        out["integrator_battery"] = {
            "analytic": analytic,
            "engine_a_vs_analytic_own": (battery["engine_a_trimesh_m3"] - analytic) / analytic,
            "self_vs_engine_a_artifact": battery.get("self_vs_engine_a"),
            "watertight": battery["topology"]["watertight"],
            "genus": battery["topology"]["closed_genus"],
        }
    return out


# ---------------------------------------------------------------------------
# 4. gates: committed vs recomputed from the live result.json
# ---------------------------------------------------------------------------
def gates() -> dict:
    res = json.loads((BENCH / "result.json").read_text())
    committed = json.loads((BENCH / "gates_verified.json").read_text())
    live = {g["name"]: g for g in res["gates"]}
    rows = res["live"]["rows"]
    ctx = res["live"]["context"]
    mine = {
        "G1_engine_agreement": max(r["engine_spread"] for r in rows),
        "G2_reference_at_finest_wigley": abs(rows[-1]["engine_a_vs_reference"]),
        "G3_quadrature_is_exact_for_the_reference_hull": abs(ctx["sectional_vs_reference"]),
        "G4_measured_order": ctx["attribution"]["observed_order"],
    }
    table = {}
    for g in committed["gates"]:
        name = g["name"]
        table[name] = {"committed_value": g["value"], "limit": g["limit"],
                       "passed": g["passed"],
                       "my_recomputed": mine.get(name),
                       "matches": (mine.get(name) is None or
                                   (isinstance(g["value"], (int, float)) and mine.get(name) is not None
                                    and abs(mine[name] - g["value"]) <= 1e-12 * max(1.0, abs(g["value"]))))}
    return {"gates": table, "committed_passed": committed["passed"],
            "repro_cmd_declared": res.get("repro_cmd"),
            "environment": res.get("environment"),
            "numpy_stl_reported": res.get("environment", {}).get("numpy_stl")}


# ---------------------------------------------------------------------------
# 5. posted report vs artifacts
# ---------------------------------------------------------------------------
def report_vs_artifacts() -> dict:
    rows = q("SELECT id, title, markdown, jsonb_array_length(figures) FROM tensor_reports "
             "WHERE mission_id=%s ORDER BY id DESC", (MISSION,))
    if not rows:
        return {"report_rows": 0}
    rid, title, md, nfig = rows[0]
    payload = json.loads((ART / "results.json").read_text())
    cases = {c["key"]: c for c in payload["cases"]}
    checks = []

    # (a) wigley ladder table in report §3 must carry the artifact numbers
    for r in cases["wigley"]["rows"]:
        for token, expect in ((f"{r['engine_a_trimesh_m3']:.6f}", r["engine_a_vs_reference"]),
                              (f"{100 * r['engine_a_vs_reference']:+.4f}", None)):
            checks.append({"token": token, "in_report": token in md, "expects": expect})

    # (b) committed multi-hull table: are artifact volumes printed at full precision?
    truncations = []
    for r in cases["wigley"]["rows"]:
        v = r["engine_a_trimesh_m3"]
        one_dp = f"| {v:.1f} |"
        if one_dp in md:
            truncations.append({"artifact_value": v, "rendered_in_report": f"{v:.1f}",
                                "relative_error_of_rendering": abs(v - float(f"{v:.1f}")) / v})
    for key in ("kcs-fullscale",):
        r = cases[key]["rows"][-1]
        token = f"{r['engine_a_trimesh_m3']:.1f}"
        if token in md:
            truncations.append({"artifact_value": r["engine_a_trimesh_m3"],
                                "rendered_in_report": token,
                                "relative_error_of_rendering":
                                    abs(r["engine_a_trimesh_m3"] - float(token)) / r["engine_a_trimesh_m3"]})

    # (c) the report's honest-FAIL sentence must match an actual test
    honest_fail = "61×13" in md or "61x13" in md
    return {
        "report_rows": len(rows), "report_id": rid, "title": title,
        "figures": nfig, "markdown_chars": len(md),
        "token_checks": checks,
        "tokens_present": sum(1 for c in checks if c["in_report"]), "tokens_checked": len(checks),
        "precision_truncations": truncations,
        "honest_fail_sentence_present": honest_fail,
        "sections_present": {k: (k in md) for k in ("①", "②", "③", "④", "⑤", "⑥")},
    }


# ---------------------------------------------------------------------------
# 6. repro: run the subject's declared command in a disposable copy
# ---------------------------------------------------------------------------
def repro_run() -> dict:
    tmp = Path("/tmp/m20544_repro")
    if tmp.exists():
        shutil.rmtree(tmp)
    shutil.copytree(SUBJ, tmp, ignore=shutil.ignore_patterns(".git", ".pytest_cache", "__pycache__"))
    log: dict = {"copy": str(tmp)}
    bench = subprocess.run([sys.executable, "benchmarks/verified/volume_cross_engine/run_benchmark.py"],
                           cwd=tmp, capture_output=True, text=True, timeout=3600)
    log["run_benchmark"] = {"returncode": bench.returncode,
                            "stdout_tail": bench.stdout.strip().splitlines()[-3:],
                            "stderr_tail": bench.stderr.strip().splitlines()[-2:]}
    shipped = json.loads((BENCH / "result.json").read_text())
    rerun = json.loads((tmp / "benchmarks/verified/volume_cross_engine/result.json").read_text())
    log["gates_identical"] = (shipped["gates"] == rerun["gates"])
    log["live_rows_identical"] = (shipped["live"]["rows"] == rerun["live"]["rows"])
    log["committed_problems"] = rerun["committed"]["problems"]
    tests = subprocess.run([sys.executable, "-m", "pytest",
                            "tests/test_volume_xval.py", "tests/test_volume_xval_entry.py",
                            "-q"], cwd=tmp, capture_output=True, text=True, timeout=3600)
    tail = tests.stdout.strip().splitlines()[-2:]
    log["pytest"] = {"returncode": tests.returncode, "tail": tail}
    # the report tool regenerates the shipped artifact (results.md unchanged except stamp)
    import os
    env = dict(os.environ, PYTHONPATH=str(tmp))
    gen = subprocess.run([sys.executable, "tools/volume_xval_report.py",
                          "--out", str(tmp / "artifacts/xval_geom_volcross")],
                         cwd=tmp, capture_output=True, text=True, timeout=3600, env=env)
    log["report_tool"] = {"returncode": gen.returncode,
                          "stdout_tail": gen.stdout.strip().splitlines()[-2:]}
    md_shipped = (ART / "results.md").read_text().splitlines()
    md_rerun = (tmp / "artifacts/xval_geom_volcross/results.md").read_text().splitlines()
    diffs = [i for i, (a, b) in enumerate(zip(md_shipped, md_rerun)) if a != b]
    log["results_md_diff_lines"] = diffs
    log["results_md_identical_except_stamp"] = (
        len(diffs) == 1 and md_shipped[diffs[0]].startswith("Generated:"))
    j_shipped = json.loads((ART / "results.json").read_text())
    j_rerun = json.loads((tmp / "artifacts/xval_geom_volcross/results.json").read_text())
    log["results_json_identical"] = j_shipped == j_rerun
    if not log["results_json_identical"]:
        log["results_json_keys_differing"] = [
            k for k in set(j_shipped) | set(j_rerun) if j_shipped.get(k) != j_rerun.get(k)]

    # CLI entry: exit codes 0 / 1 / 2 are part of the contract
    def cli(*args):
        r = subprocess.run([sys.executable, "-m", "tensorcad.cli", "volume-xval", *args],
                           cwd=tmp, capture_output=True, text=True, timeout=1800, env=env)
        return {"args": list(args), "returncode": r.returncode,
                "last_line": (r.stdout.strip().splitlines() or [""])[-1]}

    log["cli"] = [cli("wigley"),
                  cli("wigley", "--stations", "61", "--levels", "13"),
                  cli("mesh", "/tmp/does_not_exist.stl")]
    return log


def main() -> int:
    ev = {"audit_mission": AUDIT_MISSION, "subject_mission": MISSION, "anchor": ANCHOR,
          "repo": "TensorCAD", "subject_worktree": str(SUBJ),
          "inventory": inventory(), "referee": referee(), "recompute": recompute(),
          "gates": gates(), "report": report_vs_artifacts(), "repro": repro_run()}
    (HERE / "m20544_evidence_raw.json").write_text(
        json.dumps(ev, indent=1, ensure_ascii=False, default=str))
    print(json.dumps({"referee_rows": ev["referee"]["rows_for_anchor"],
                      "wigley_order_own": ev["recompute"]["cases"]["wigley"]["my_measured_order"],
                      "wigley_order_artifact":
                          ev["recompute"]["cases"]["wigley"]["artifact_measured_order"],
                      "wigley_finest_vs_ref_own":
                          ev["recompute"]["cases"]["wigley"]["finest_engine_a_vs_reference_own"],
                      "kcs_cb_actual": ev["recompute"]["cases"]["kcs-fullscale"]["cb_from_pipeline_volume_own"],
                      "gates_identical": ev["repro"]["gates_identical"],
                      "pytest": ev["repro"]["pytest"],
                      "cli": ev["repro"]["cli"],
                      "report_truncations": ev["report"]["precision_truncations"]},
                     indent=1, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())