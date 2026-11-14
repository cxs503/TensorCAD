"""Land the board-audit verdict for m20416 (XVAL-GEOM-VOLCROSS) in TensorCAD.

Verdict = L1 (reopened), not L2: the *algorithm* is L2-grade and proved by a
disposable repair probe, but the delivered artefact contains no runnable,
repository-resident implementation, so the capability has not landed.

Every pointer below is an absolute path that exists today, so the maturity gate's
own re-verification (files_exist / impl_real / callers_verified) reads the same
facts the text claims.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

WORK = Path("/root/tp-agent/work/m20416")
HERE = Path(__file__).resolve().parent
EVID = HERE / "evidence_m20416.json"
PROBE = Path("/tmp/m20464_probe")
MY_MISSION = 20464
ANCHOR = "XVAL-GEOM-VOLCROSS"
LEVEL = 1

ev = json.loads(EVID.read_text())
impl = WORK / "tensorcad/volume_xval.py"
tests = WORK / "tests/test_volume_xval.py"
tool = WORK / "tools/volume_xval_study.py"
art = WORK / "artifacts/volume_xval_study.json"

for p in (impl, tests, tool, art):
    assert p.exists(), f"pointer check failed: {p}"
assert not ev["compile"]["impl_compiles"], "as-delivered tree compiles now; re-audit"
assert ev["referee"]["rows_for_anchor"] == 0, "referee rows appeared; re-audit"
assert ev["probe"]["pytest_rc"] == 0 and "13 passed" in " ".join(ev["probe"]["pytest_after_repair"])
assert (PROBE / "tensorcad/volume_xval.py").exists()

S = ev["numbers"]["ships"]


def dev(ship: str, prefix: str):
    for r in S[ship]["routes"]:
        if r["engine"].startswith(prefix):
            return r.get("abs_err_pct")
    return None


loft = ev["numbers"]["loft_convergence"]["KCS"]
rich = loft["richardson_volume_m3"]
ref_kcs = S["KCS"]["reference"]["volume_m3"]

evidence = {
    "impl_files": str(impl),
    "entry": str(tool),
    "run_artifacts": str(art),
    "tests": [
        "test_ray_engine_reproduces_exact_box",
        "test_ray_engine_reproduces_inscribed_sphere_volume",
        "test_ray_uncertainty_is_a_measured_bar",
        "test_trimesh_route_agrees_with_the_repository_summation",
        "test_brep_engine_status_is_never_a_silent_pass",
        "test_nurbs_quadrature_route_converges_in_order",
        "test_wigley_row_is_inside_the_band_and_its_routes_are_exact",
        "test_fleet_rows_agree_between_engines_and_disagree_with_reference",
        "test_closure_reference_is_labelled_withheld",
        "test_loft_convergence_is_second_order_ish_and_extrapolates_inside_the_band",
        "test_direction_sweep_separates_stations_from_levels",
        "test_locate_discrepancy_names_the_geometry_when_engines_agree",
        "test_report_table_is_markdown_and_carries_evidence_per_row",
    ],
    "integration_callers": "NONE: grep -rn 'volume_xval|cross_validate_(hull|standard_hulls)' over the "
                           "whole TensorCAD tree excluding the three deliverables returns no caller; "
                           "tensorcad/api/ has no reference. Not wired into solver or REST.",
    "benchmark": (
        "tensor_benchmark_cases WHERE anchor_slug='XVAL-GEOM-VOLCROSS' -> 0 rows (table has 131 "
        "rows): NO machine referee. Judged on the module's own declared 0.3 % band + in-artefact "
        f"references: Wigley analytic ref {S['Wigley']['reference']['volume_m3']} -> "
        f"{dev('Wigley','self-'):.4f}% (self divergence) / {dev('Wigley','ray'):.4f}% (ray); "
        f"KCS ref {S['KCS']['reference']['volume_m3']} -> {dev('KCS','self-'):.4f}%; "
        f"JBC ref {S['JBC']['reference']['volume_m3']} -> {dev('JBC','self-'):.4f}%; "
        f"Series60 ref {S['Series60']['reference']['volume_m3']} -> "
        f"{dev('Series60','self-'):.4f}%; engine spread self-vs-trimesh 1.0e-8; KCS Richardson "
        f"limit {rich:.2f} m3 = {(rich-ref_kcs)/ref_kcs*100:+.4f}% vs ref. All < 3 % (owner "
        "standard), but the 0.3 % band is exceeded on 3 of 4 hulls and those 3 references are "
        "self-declared non-external (recalled literature verified=false / closure)."
    ),
    "repro_cmd": (
        "PYTHONPATH=. python -m pytest tests/test_volume_xval.py -q  # AS DELIVERED: 1 collection "
        "error, IndentationError volume_xval.py:1053; PYTHONPATH=. python tools/volume_xval_study.py"
        " # AS DELIVERED: ImportError. In the disposable copy /tmp/m20464_probe after a 6-line "
        "mechanical repair: 13 passed in 22.18 s and the study regenerates the artefact."
    ),
    "external_anchor": (
        "Wigley closed form L*B*T*n*m/((n+1)(m+1)) evaluated in code (exact, REAL) + third-party "
        "trimesh 5.0.0 through a binary-STL round trip (independent code base and parser; 1.0e-8 "
        "agreement with the in-memory divergence sum = float32 STL precision). B-rep solid route "
        "(gmsh/OCP) honestly WITHHELD: no loongarch64 wheel on the configured index."
    ),
    # ---- below: the audit record the gate does not score, kept for the human reviewer ----
    "verdict_level": {
        "level": LEVEL,
        "why_not_L2": (
            "L2 = 真代码 输入→计算→输出. The code exists (1072 lines) and the numbers reproduce, "
            "but as delivered the module is not importable (IndentationError volume_xval.py:1053), "
            "all six deliverables are untracked (git status ??, head_commit=null, branch=null, "
            "HEAD:tensorcad/volume_xval.py absent), so the TensorCAD repository today contains no "
            "runnable volume cross-validation. A capability that cannot be imported has not "
            "landed; L2 on potential would write off a broken delivery."
        ),
        "substance_measured": (
            "A 6-line mechanical repair (remove 1 stray token, 1 duplicated guard, dedent 20 "
            "lines; no logic added or removed) in /tmp/m20464_probe yields 13/13 tests green in "
            "22.18 s and regenerates the committed artefact with 179/179 numeric leaves identical "
            "and 3 explanatory strings changed. So the repair cost is trivial and the numbers are "
            "real; promotion to L2 needs commit + fix, and to L4 needs a referee row."
        ),
        "L3": "FAIL: 0 callers in the whole tree, none in tensorcad/api/.",
        "L4": "FAIL: no referee row (0/131 anchors), 3 of 4 references non-external, and not "
              "reproducible as shipped.",
        "L5": "FAIL: nothing beyond the above.",
    },
    "defect": {
        "where": "tensorcad/volume_xval.py:1043-1073 (_append_sweeps)",
        "what": "botched manual patch: stray bare token `station` (l.1049), duplicated "
                "`if not sweeps: return` guard (l.1050-1051), body over-indented by 4 "
                "(l.1053-1073), file ends without trailing newline",
        "class": "packaging, not algorithmic",
    },
    "overstated_claims": [
        "tensor_dev_missions.evidence.finished_status='succeeded' and accepted=true while the same "
        "result text records the module as broken ('文件当前 IndentationError ... 如实记录为 FAIL', "
        "fail_stage=budget_exhausted+forced_finish, server_quality_audit.gate_valid=false).",
        "The agent's self-reported '13 passed / 19.78 s' refers to a code state that no longer "
        "exists in the tree; as shipped, pytest cannot even collect.",
        "All six deliverables are untracked; nothing of m20416 reached the TensorCAD repository.",
        "tensor_reports WHERE mission_id=20416 -> 0 rows: the mission report was never posted.",
    ],
    "audit_scope_limits": (
        "Cross-check only: the audited worktree was never modified; the repair probe ran on a "
        "disposable copy under /tmp/m20464_probe. No judgement on trimesh's own correctness "
        "beyond the 1.0e-8 STL round trip, nor on the recalled literature values behind KCS/JBC."
    ),
    "referee_query": "SELECT case_name,metric,tolerance FROM tensor_benchmark_cases "
                      "WHERE anchor_slug='XVAL-GEOM-VOLCROSS' -> 0 rows",
    "audit_script": str(HERE / "audit_m20416.py"),
    "evidence_bundle": str(EVID),
}

row = {
    "repo": "TensorCAD",
    "scope": f"m20416 {ANCHOR}",
    "mission_id": MY_MISSION,
    "level": LEVEL,
    "verdict_by": "board-audit",
}
ev_json = json.dumps(evidence, ensure_ascii=False)


def psql(sql: str) -> str:
    return subprocess.run(
        ["psql", "-U", "postgres", "-d", "tensorplatform", "-t", "-A", "-c", sql],
        capture_output=True, text=True, check=False,
    ).stdout.strip()


SQL = (
    "UPDATE tensor_maturity_verdicts SET level={lv}, evidence='{ev}', verdict_by='board-audit' "
    "WHERE repo='TensorCAD' AND scope='m20416 {anchor}' AND mission_id={mid} "
    "AND coalesce(gate_state,'pending')='pending' RETURNING id, level, gate_state"
).format(lv=row["level"], ev=ev_json.replace("'", "''"), anchor=ANCHOR, mid=MY_MISSION)

if __name__ == "__main__":
    (HERE / "verdict_row_m20416.json").write_text(json.dumps({**row, "evidence": evidence},
                                                             indent=1, ensure_ascii=False))
    print(psql(SQL) or "UPDATE affected no pending row")
    print(psql("SELECT id, repo, scope, level, verdict_by, gate_state FROM tensor_maturity_verdicts "
               "ORDER BY id DESC LIMIT 2"))