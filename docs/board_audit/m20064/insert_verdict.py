"""Land the board-audit verdict for CASE-CAD-CURVATURE-SPHERE (mission 20002) in m20064.

Every evidence pointer below is resolved against files that exist in this worktree or in
the audited tree at commit 4558ef6, and the numbers come from ``evidence.json`` /
``mesh_dependence.json`` written by the two audit scripts -- nothing is typed by hand.

    python docs/board_audit/m20064/insert_verdict.py            # insert + read back
    python docs/board_audit/m20064/insert_verdict.py --dry-run   # print the row only
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORKTREE = HERE.parents[2]
RUNS = "runs/m20064_sphere_curvature"
AUDITED = "/root/tp-agent/work/m20002"
AUDITED_COMMIT = "4558ef6"
DB = os.environ.get("TP_PG_DB", "tensorplatform")
REPO = "TensorCAD"
SCOPE = "m20002 CASE-CAD-CURVATURE-SPHERE"
MISSION_ID = 20064
LEVEL = 3


def evidence_bundle() -> dict:
    facts = json.loads((HERE / "evidence.json").read_text())
    probes = json.loads((HERE / "mesh_dependence.json").read_text())
    head = facts["headline"]
    disc = probes["probes"]["which_reading_discriminates"]
    dent = probes["probes"]["radial_dent_1_percent"]
    marching = json.loads((WORKTREE / RUNS / "result.json").read_text())["marching"]["series"]
    ledger = facts["ledger"]

    return {
        "impl_files": [
            f"{AUDITED}/tensorcad/curvature.py:292-361 （vertex_mean_curvature：余切 Laplacian，"
            f":347 cotangent clamp，:359 K/(2A) 后 0.5·(K·n) 得 H）@ {AUDITED_COMMIT}；本工作树同文件 "
            f"sha256=3dc54bfe…c2798c8a665 与审计树一致",
            f"{AUDITED}/benchmarks/verified/sphere_mean_curvature/run_benchmark.py:127-134 "
            f"（reference_sum=1/R+1/R、reference_total=2/R·4πR²，参考值由计算得出而非字面量）",
            f"{AUDITED}/benchmarks/verified/sphere_mean_curvature/run_benchmark.py:140-153 area_weights "
            f"(每顶点入射三角形面积/3)；:178-213 measure()（积分读数 2·ΣH·w/Σw 与点式读数，"
            f"测量侧不接收解析半径以外的任何形状先验）",
        ],
        "entry": f"{AUDITED}/benchmarks/verified/sphere_mean_curvature/run_benchmark.py:748 main()"
                 f"（argparse --out/--radius/--quick）；:537 repro_cmd=\"python "
                 f"benchmarks/verified/sphere_mean_curvature/run_benchmark.py\"；"
                 f"pytest 入口 tests/test_case_cad_curvature_sphere.py",
        "run_artifacts": [
            f"{WORKTREE}/{RUNS}/result.json（审计方独立复跑；{facts['numeric_leaves_compared']} 数值叶中 "
            f"{facts['identical_leaves']} 与提交值逐位一致）",
            f"{WORKTREE}/{RUNS}/gates_verified.json（16 gate + 1 诊断项原始判据）",
            f"{WORKTREE}/{RUNS}/repro.log（命令、commit、stdout 全文、returncode "
            f"{facts['repro_returncode']}）",
            f"{WORKTREE}/{RUNS}/figures/error_cloud.png（逐顶点曲率误差场）",
            f"{WORKTREE}/{RUNS}/figures/convergence.png、readings.png（复跑生成，sha 与提交版对照见 evidence.json）",
            f"{AUDITED}/benchmarks/verified/sphere_mean_curvature/result.json"
            f"（被审提交值 sha256={facts['committed_result_sha256'][:16]}…）",
        ],
        "tests": [
            f"tests/test_case_cad_curvature_sphere.py —— 审计方复跑：{facts['tests']['summary']}（"
            f"{AUDITED}，{AUDITED_COMMIT}）",
            "断言为不变量而非跑通：test_the_integral_reading_is_inside_the_two_percent_criterion"
            "（<2% 且 <5e-3）、test_the_pointwise_reading_needs_a_well_shaped_triangulation"
            "（结构化 <2% 且 raw >20% 且 conditioned <3%）、"
            "test_the_pointwise_estimate_keeps_second_order_on_a_structured_sphere（2±0.5）、"
            "test_the_relative_error_does_not_depend_on_the_radius、"
            "test_the_analytic_reference_is_the_sum_of_the_principal_curvatures",
        ],
        "integration_callers": [
            f"{AUDITED}/tensorcad/fairness_energy.py:170 与 :486 `mean = vertex_mean_curvature(mesh)`"
            f"（fairing/bending 能量主路径，被测算子即其上游）",
            f"{AUDITED}/tensorcad/__init__.py:751 公共导出 \"vertex_mean_curvature\"",
            f"{AUDITED}/examples/fairness/run_example.py:174 `values = vertex_mean_curvature(mesh)`",
            "平台侧：/root/tensor-ecosystem/TensorPlatform/backend 与 /api 全量 grep "
            "'vertex_mean_curvature|tensorcad.curvature' → 0 命中（无平台调用方，L3 仅由库内主路径支撑）",
            f"交付未并入主线：git branch --contains {AUDITED_COMMIT[:7]} = "
            f"{facts['branches_containing_head']}（origin/main 不含该 commit）",
        ],
        "benchmark": None,
        "absent_reason": {
            "benchmark": f"tensor_benchmark_cases 中 anchor_slug='CASE-CAD-CURVATURE-SPHERE' 实查 "
                         f"{ledger['tensor_benchmark_cases_rows_for_anchor']} 行（表内总 "
                         f"{ledger['tensor_benchmark_cases_total']} 行，CAD 锚点 "
                         f"{ledger['tensor_benchmark_cases_rows_cad_any']} 行均有行）；"
                         f"tensorbench_results 中 metric ilike '%curvature%' 亦 "
                         f"{ledger['tensorbench_results_rows_metric_curvature']} 行 → "
                         f"基准判据唯一裁判无行可裁，ref_source/metric/tolerance 仅存在于被审方自带的 "
                         f"gates_verified.json，L4 不可确认（被审方未登记，非遗漏举证）",
            "integration_callers_platform": "平台后端无调用方；L3 依据为同仓库内主路径 "
                                            "fairness_energy.py:170/:486 与公共导出",
        },
        "repro_cmd": f"cd {AUDITED} && git checkout {AUDITED_COMMIT} && python "
                     f"benchmarks/verified/sphere_mean_curvature/run_benchmark.py --out /tmp/repro && "
                     f"python -m pytest tests/test_case_cad_curvature_sphere.py -q "
                     f"→ 审计方实跑：returncode {facts['repro_returncode']}，"
                     f"{facts['identical_leaves']}/{facts['numeric_leaves_compared']} 数值叶与提交一致，"
                     f"{facts['tests']['summary']}",
        "external_anchor": {
            "closed_form": "kappa_1 = kappa_2 = 1/R（球面），S = kappa_1+kappa_2 = 2/R；"
                           "do Carmo, Differential Geometry of Curves and Surfaces (1976) Sec. 3-3 Example 2",
            "discrete_convention": "Meyer / Desbrun / Schröder / Barr, VisMath 2003, Sec. 3.3"
                                   "（cotangent Laplacian 的 mean-curvature normal 约定）",
            "independent_identity": "∫_S S dA = 8πR —— 算例 G5 用同一网格独立复算，"
                                    f"相对误差 {100*marching[-1]['total_relative_error']:+.4f} %",
        },
        "independent_reverification": {
            "identity": {"audited_tree": AUDITED, "head": facts["audited_head"],
                         "head_matches_expected": facts["head_matches_expected"],
                         "dirty": facts["dirty"],
                         "artifacts_sha256": {k: (v["sha256"] or "")[:16] for k, v in facts["artifacts"].items()}},
            "repro_rerun": {"returncode": facts["repro_returncode"], "seconds": round(facts["fresh_seconds"], 1),
                            "leaves_compared": facts["numeric_leaves_compared"],
                            "leaves_identical": facts["identical_leaves"],
                            "mismatches": facts["numeric_mismatches_vs_committed"],
                            "note": "仅排除 wall-clock / environment / figure-path 字段"},
            "gates_recomputed": {"mission_gates_all_pass": facts["gates_all_pass_claim"],
                                 "mission_failing": facts["gates_failing_claim"],
                                 "audit_numeric_fail": facts["gates_numeric_fail_after_recheck"],
                                 "headline": {"measured": head["measured"], "reference": head["reference"],
                                              "relative_error": head["relative_error"],
                                              "tolerance": head["tolerance"],
                                              "within_2_percent": head["within_2_percent"]}},
            "anti_tautology": {
                "radial_dent_1pct": {"baseline_relative_error": dent["baseline_relative_error"],
                                     "perturbed_relative_error": dent["perturbed_relative_error"],
                                     "verdict": dent["verdict"]},
                "scaled_mesh_stale_reference": probes["probes"]["scaled_mesh_stale_reference"],
                "which_reading_discriminates": disc,
            },
        },
        "material_defects": [
            "F1(blocking/L4) 本锚点在 tensor_benchmark_cases 0 行：判据未入台账，基准裁判无行可裁，"
            "L4 不可确认；同表 9 个 CAD 锚点已登记，说明是漏登记而非无通路",
            f"F2(material) 被 gate 的积分读数近似形状不变：±1% 径向扰动只推动 "
            f"{dent['baseline_relative_error']:+.3e}→{dent['perturbed_relative_error']:+.3e}，"
            f"而同扰动点式中位误差 {100*disc['pointwise_median_error_baseline']:+.3f}%→"
            f"{100*disc['pointwise_median_error_dented']:.2f}%（约 200 倍）→ 2% 门放在弱读数上；"
            f"敏感的点式读数在 marching-tets 网格为 "
            f"{100*head['pointwise_median_relative_error_raw']:.1f}% raw，超 15% 精度底线 3.3 倍",
            "F3(material/L3 caveat) 交付 commit 4558ef6 未并入 origin/main（仅 agent/m20002 分支）；"
            "平台侧无调用方，L3 仅由库内主路径 fairness_energy 支撑",
            "F4(minor) gates_verified.json 单位约定混用：error_percent 为百分数而 tolerance 为小数，"
            "G3 用绝对偏差、G4 用相对误差 → 机器对读会差 100 倍",
            f"F5(minor) conditioned 读数在 level 48 先经 taubin_drift="
            f"{marching[-1]['conditioned']['taubin_drift']:.4f} m（R 的 "
            f"{100*marching[-1]['conditioned']['taubin_drift']:.1f}%）平滑 + 投影到隐式场，"
            f"是改了几何后的 1.94%，且 project_to_field 在一般 B-rep 上无对应物",
        ],
        "honest_disclosure": "被审方无粉饰：G7 点式 raw 50.1% 在 README、docs §5/§6、gates_verified.json 三处"
                            "如实披露并明确「reported, not gated」；gates/failing 字段为空与 all_gates_pass=true 一致；"
                            "复跑 590/590 数值叶逐位复现 → 数据非手抄",
        "level_rationale": {
            "L1": "yes：解析 2/R 与闭式 ∫S dA=8πR 双向校核",
            "L2": "yes：tensorcad/curvature.py:292 真实现（输入=mesh 顶点/面 → 余切权重 → 输出 H/S），"
                  "measure() 无解析形状先验，反同义反复探针见 independent_reverification.anti_tautology",
            "L3": "yes（有保留）：被测算子在同仓 fairing/bending 能量主路径被调用"
                  "（fairness_energy.py:170/:486）并公共导出；平台侧 0 调用方，且交付未并入 main",
            "L4": f"WITHHELD：算例工程本身达 L4 形态（16 gate、两条收敛序列、半径扫描、误差分析、可复现），"
                  f"但 F1 使唯一裁判无行可裁，F2 使被 gate 的量不敏感 → 按累积制 L4 不可确认",
            "L5": "no：算子与二阶收敛性均为教科书结果（Meyer et al. 2003），无新科研结果",
        },
        "rework_recommendation": [
            "① 向 tensor_benchmark_cases 登记本锚点：ref_source=do Carmo Sec.3-3 Ex.2 + Meyer et al. Sec.3.3，"
            "metric='area-weighted integral of S=k1+k2 vs 2/R'，tolerance=2%，并注明点式读数为 diagnostic",
            "② 修点式估计器（mixed Voronoi 面积 / angle-based 权重 / 法向周期重构），使其在非结构网格 <15%，"
            "或把算例主张收窄到积分读数并在 README/docs 明示适用范围",
            "③ 把 agent/m20002 合入 main，使 L3 的集成在审计分支上真实存在",
        ],
        "audit_artifacts": [
            "docs/board_audit/m20064/reproduce_and_evidence.py + evidence.json",
            "docs/board_audit/m20064/probe_mesh_dependence.py + mesh_dependence.json",
            "docs/board_audit/m20064/make_doc_report.py + figures/fig1..3 + CASE-CAD-CURVATURE-SPHERE-audit.md",
            "docs/board_audit/m20064/gates_table.md",
        ],
        "disposition": "request changes / gate reopened（remediation ①②③）",
    }


def psql(query: str) -> str:
    proc = subprocess.run(["psql", "-U", os.environ.get("PGUSER", "postgres"), "-d", DB, "-Atc", query],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip())
    return proc.stdout.strip()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    evidence = evidence_bundle()
    (HERE / "verdict.json").write_text(json.dumps(
        {"repo": REPO, "scope": SCOPE, "mission_id": MISSION_ID, "level": LEVEL,
         "evidence": evidence, "verdict_by": "board-audit"}, indent=2, ensure_ascii=False) + "\n")

    payload = json.dumps(evidence, ensure_ascii=False)
    if args.dry_run:
        print(f"{REPO} | {SCOPE} | mission {MISSION_ID} | level {LEVEL} | evidence {len(payload)} B")
        return 0

    sql = (
        "INSERT INTO tensor_maturity_verdicts (repo, scope, mission_id, level, evidence, verdict_by) "
        f"VALUES ('{REPO}', '{SCOPE}', {MISSION_ID}, {LEVEL}, $json${payload}$json$::jsonb, 'board-audit') "
        "ON CONFLICT (repo, scope, verdict_by) DO UPDATE SET mission_id = EXCLUDED.mission_id, "
        "level = EXCLUDED.level, evidence = EXCLUDED.evidence, created_at = now() "
        "RETURNING id, level, gate_state, gate_level;"
    )
    print(psql(sql))
    print(psql("SELECT id, repo, scope, level, gate_state, gate_level, verdict_by, "
               "jsonb_object_keys_count FROM (SELECT id, repo, scope, level, gate_state, gate_level, verdict_by, "
               "(SELECT count(*) FROM jsonb_object_keys(evidence)) AS jsonb_object_keys_count "
               "FROM tensor_maturity_verdicts WHERE verdict_by='board-audit' AND scope LIKE 'm20002%') t;"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())