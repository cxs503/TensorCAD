#!/usr/bin/env python3
"""Cross-check audit of mission m17333 (anchor GAP-CAD-T1-IGES-PIPE).

Read-only: runs the subject mission's own reproduce command, re-measures the
artifacts with independent tooling (trimesh), re-reads the only third-party
trimmed IGES (KP505, outside the subject tree) and records git/DB facts.
Every number in docs/audit_m17343_m17333_iges_pipe.md comes from this JSON.

Usage:
    python tools/audit_m17343_crosscheck.py \
        --subject /root/tp-agent/work/m17333 \
        --kp505 /root/tp-agent/work/m15295/benchmarks/anchors/sources/KP505_modelScale_250mm.igs \
        --out artifacts/audit_m17343

Nothing here writes to the subject tree except its own `reproduce.py` output
(runs/gap_t1_iges, already the committed evidence path; bytes are checked
against the committed sha256).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import time
from collections import Counter
from pathlib import Path

ROOTS = Path(__file__).resolve().parents[1]


def run(cmd, cwd=None, timeout=300):
    t = time.time()
    p = subprocess.run(cmd, shell=isinstance(cmd, str), cwd=cwd, capture_output=True,
                       text=True, timeout=timeout)
    return {"cmd": cmd if isinstance(cmd, str) else " ".join(cmd), "rc": p.returncode,
            "secs": round(time.time() - t, 2),
            "stdout": p.stdout[-4000:], "stderr": p.stderr[-2000:]}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def census(path: Path) -> dict:
    types = Counter()
    for line in path.read_text(errors="replace").replace("\r", "").splitlines():
        if len(line) == 80 and line[72] == "D":
            try:
                types[int(line[0:8])] += 1
            except ValueError:
                pass
    return dict(sorted(types.items()))


def kp505_loop_residuals(kp505: Path, budget_s: float = 25.0) -> dict:
    """Invert the 3-D trim loops of the third-party trimmed file; bounded in time."""
    from tensorcad.geomquality.iges_trim import read_iges_trimmed
    from tensorcad.geomquality.trim_mesh import _loop_to_uv
    m = read_iges_trimmed(kp505)
    extent = max(abs(c) for f in m.faces for loop in f.loops for p in loop for c in p)
    tol = 1e-6 * extent  # same rule as tessellate_trimmed_model
    start = time.time()
    loops, worst = [], 0.0
    for f in m.faces:
        for k, loop in enumerate(f.loops):
            if time.time() - start > budget_s:
                loops.append({"face": f.index, "loop": k, "points": len(loop),
                              "status": "not_run_time_budget"})
                continue
            uv, res = _loop_to_uv(f.surface, loop, tol)
            worst = max(worst, res)
            loops.append({"face": f.index, "loop": k, "points": len(loop),
                          "residual": res, "ratio_to_tolerance": res / tol})
    return {"faces_read": len(m.faces), "metres_per_unit": m.metres_per_unit,
            "model_extent_file_units": extent, "residual_tolerance": tol,
            "worst_residual": worst, "loops": loops,
            "loops_over_tolerance": sum(1 for l in loops
                                        if l.get("residual", 0.0) > tol),
            "loops_measured": sum(1 for l in loops if "residual" in l)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--subject", default="/root/tp-agent/work/m17333")
    ap.add_argument("--kp505", default="/root/tp-agent/work/m15295/benchmarks/anchors/"
                                      "sources/KP505_modelScale_250mm.igs")
    ap.add_argument("--kp-report", default="/root/tp-agent/work/m17250/runs/gap_t1_iges/"
                                           "KP505_modelScale_250mm.report.json")
    ap.add_argument("--out", default="artifacts/audit_m17343")
    args = ap.parse_args()
    subj, out = Path(args.subject), (ROOTS / args.out)
    out.mkdir(parents=True, exist_ok=True)
    checks: dict = {"audit_mission": 17343, "subject_mission": 17333,
                    "anchor": "GAP-CAD-T1-IGES-PIPE", "subject_tree": str(subj),
                    "ts": time.strftime("%FT%TZ")}

    # 1. git facts: is the delivery in a commit, in main, did it touch code?
    checks["git"] = {
        "head": run("git rev-parse HEAD", cwd=subj)["stdout"].strip(),
        "head_subject": run("git log --oneline -1", cwd=subj)["stdout"].strip(),
        "branch": run("git branch --show-current", cwd=subj)["stdout"].strip(),
        "files_in_head_commit": run("git show --stat --name-only --format= HEAD", cwd=subj)["stdout"].split(),
        "status_porcelain": run("git status --short", cwd=subj)["stdout"].splitlines(),
        "code_files_touched_by_delivery": [],
        "main_head_of_audit_tree": run("git rev-parse HEAD", cwd=ROOTS)["stdout"].strip(),
        "delivery_in_audit_main": run(
            f"git merge-base --is-ancestor {run('git rev-parse HEAD', cwd=subj)['stdout'].strip()} HEAD "
            f"&& echo YES || echo NO", cwd=ROOTS)["stdout"].strip(),
    }
    touched = [f for f in checks["git"]["files_in_head_commit"]
               if f.startswith("tensorcad/") or f.startswith("tests/")]
    checks["git"]["code_files_touched_by_delivery"] = touched
    checks["git"]["touched_tensorcad_package"] = [f for f in touched if f.startswith("tensorcad/")]

    # 2. inputs: provenance sha256 + entity census (is it trimmed or a patch set?)
    prov = json.loads((subj / "artifacts/gap_t1_iges/provenance.json").read_text())
    prov_by_name = {f["path"].split("/")[-1]: f for f in prov["files"]}
    inputs = {}
    for p in sorted((subj / "artifacts/gap_t1_iges/inputs").glob("*.iges")):
        inputs[p.name] = {"sha256": sha256(p), "bytes": p.stat().st_size,
                          "entity_census": census(p),
                          "in_provenance": p.name in prov_by_name,
                          "sha_matches_provenance": prov_by_name.get(p.name, {}).get("sha256") == sha256(p)}
    checks["inputs"] = inputs

    # 3. reproduce the subject run, compare with the committed artifacts
    rec = {p.name: json.loads(p.read_text())
           for p in (subj / "runs/gap_t1_iges").glob("*.report.json")}
    before = {k: {"volume": v["after_repair"]["signed_volume"],
                  "triangles": v["after_repair"]["triangles"],
                  "sha256": sha256(subj / "runs/gap_t1_iges" / k.replace(".report.json", ".stl"))}
              for k, v in rec.items()}
    checks["reproduce"] = run("python artifacts/gap_t1_iges/reproduce.py", cwd=subj)
    after = {}
    for k in before:
        r = json.loads((subj / "runs/gap_t1_iges" / k).read_text())
        after[k] = {"volume": r["after_repair"]["signed_volume"],
                    "triangles": r["after_repair"]["triangles"],
                    "sha256": sha256(subj / "runs/gap_t1_iges" / k.replace(".report.json", ".stl"))}
    checks["reproduce"]["bit_identical"] = {
        k: {"volume_equal": before[k]["volume"] == after[k]["volume"],
            "triangles_equal": before[k]["triangles"] == after[k]["triangles"],
            "stl_sha_equal": before[k]["sha256"] == after[k]["sha256"]} for k in before}

    ev = json.loads((subj / "runs/gap_t1_iges/evidence.json").read_text())
    checks["benchmark_self_consistency"] = [
        {"file": f["path"].split("/")[-1], "route": f["route"],
         "iges_volume": f["signed_volume"],
         "analytic_volume": f["benchmark"]["analytic_volume"],
         "relative_error_pct": f["benchmark"]["relative_error_pct"]}
        for f in ev["files"]]
    checks["evidence_gold_standard_label"] = ev.get("gold_standard")

    # 4. independent mesh audit (trimesh, not the subject's own audit code)
    import trimesh
    ind = {}
    for stl in sorted((subj / "runs/gap_t1_iges").glob("*.stl")):
        m = trimesh.load(str(stl), process=True)
        m.merge_vertices()
        ind[stl.name] = {"faces": len(m.faces), "vertices": len(m.vertices),
                         "watertight": bool(m.is_watertight),
                         "winding_consistent": bool(m.is_winding_consistent),
                         "volume": float(m.volume), "euler": int(m.euler_number)}
    checks["independent_trimesh"] = ind

    # 5. tests (subject's committed version; and the audit tree's version)
    checks["tests_subject"] = run(["python", "-m", "pytest", "tests/test_iges_trim_import.py", "-q"], cwd=subj)
    checks["tests_subject"]["summary"] = checks["tests_subject"]["stdout"].strip().splitlines()[-1]
    checks["tests_audit_tree"] = run(["python", "-m", "pytest", "tests/test_iges_trim_import.py", "-q"], cwd=ROOTS)
    checks["tests_audit_tree"]["summary"] = checks["tests_audit_tree"]["stdout"].strip().splitlines()[-1]
    checks["tests_exercising_trimmed_end_to_end"] = run(
        "grep -rn 'type.*144\\|144/143\\|read_iges_trimmed' tests/ | wc -l", cwd=subj)["stdout"].strip()

    # 6. integration callers inside the package
    checks["integration_callers"] = {
        "grep_import_cad_or_trim_mesh_outside_geomquality": run(
            "grep -rn 'import_cad\\|trim_mesh' --include=*.py tensorcad/ | "
            "grep -v 'geomquality/trim_mesh.py' | grep -v 'geomquality/import_cli.py' | wc -l", cwd=subj)["stdout"].strip(),
        "exported_in_geomquality_init": run(
            "grep -c 'trim_mesh\\|import_cli\\|import_cad' tensorcad/geomquality/__init__.py", cwd=subj)["stdout"].strip(),
        "main_cli_reference": run("grep -c 'iges\\|geomquality' tensorcad/cli.py", cwd=subj)["stdout"].strip(),
    }

    # 7. benchmark referee table
    q = ("SELECT count(*) FROM tensor_benchmark_cases WHERE anchor_slug="
         "'GAP-CAD-T1-IGES-PIPE'")
    qt = "SELECT count(*) FROM tensor_benchmark_cases"
    checks["benchmark_referee"] = {
        "rows_for_anchor": run(f"psql -U postgres -d tensorplatform -tAc \"{q}\"")["stdout"].strip(),
        "rows_total": run(f"psql -U postgres -d tensorplatform -tAc \"{qt}\"")["stdout"].strip(),
    }

    # 8. the only third-party trimmed IGES: does the trimmed route work?
    kp = Path(args.kp505)
    checks["kp505"] = {"path": str(kp), "exists": kp.exists()}
    if kp.exists():
        checks["kp505"]["sha256"] = sha256(kp)
        checks["kp505"]["entity_census"] = census(kp)
        t0 = time.time()
        checks["kp505"]["loop_inversion"] = kp505_loop_residuals(kp)
        checks["kp505"]["loop_inversion"]["wall_seconds"] = round(time.time() - t0, 1)
    kpr = Path(args.kp_report)
    if kpr.exists():
        d = json.loads(kpr.read_text())
        checks["kp505"]["subject_prior_run"] = {
            "report": str(kpr), "route": d.get("route"),
            "faces_in_file": d["import"]["faces_in_file"],
            "faces_tessellated": d["import"]["faces_tessellated"],
            "faces_refused": d["import"]["faces_refused"],
            "refused": d["import"]["refused"],
            "triangles": d["after_repair"]["triangles"],
            "watertight": d["after_repair"]["watertight"]}

    # 9. the committed orphan input (not in provenance / manifest)
    wig = subj / "artifacts/gap_t1_iges/inputs/wigley_iii.iges"
    if wig.exists():
        from tensorcad.geomquality.trim_mesh import import_cad
        r = import_cad(wig)
        checks["wigley_orphan"] = {
            "in_provenance": "wigley_iii.iges" in prov_by_name,
            "entity_census": census(wig),
            "route": r.report["route"],
            "gates": {k: v["passed"] for k, v in r.report["gates"].items()},
            "triangles": r.report["after_repair"]["triangles"],
            "watertight": r.report["after_repair"]["watertight"]}

    # 10. subject self-claim (mission result JSON)
    got = run(["psql", "-U", "postgres", "-d", "tensorplatform", "-tAc",
               "SELECT json_build_object('status',status,'head_commit',"
               "result->>'head_commit','checks',result->>'checks','git_status',"
               "result->>'git_status') FROM tensor_dev_missions WHERE id=17333"])
    checks["mission_17333_self_claim"] = got["stdout"].strip()

    (out / "checks.json").write_text(json.dumps(checks, indent=2, default=str))
    print(json.dumps({k: checks[k] for k in (
        "benchmark_referee", "tests_subject", "tests_audit_tree",
        "integration_callers")}, indent=1, default=str)[:2000])
    print("wrote", out / "checks.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())