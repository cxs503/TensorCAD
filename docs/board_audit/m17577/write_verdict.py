#!/usr/bin/env python3
"""m17577（board-audit）落表脚本：把 m17472 的裁定写入 tensor_maturity_verdicts。

设计要点（沿用 m15207 审计纪律）
  1. 证据指针落在本审计工作树 docs/board_audit/m17577/refs/ 内（被审工作树
     /root/tp-agent/work/m17472 随时可能被回收）；gate 的 files_exist 要求纯路径。
  2. 缺失的证据键显式置 null（integration_callers / benchmark），不用文字凑数：
     L3/L4 缺就是缺，让机器门把级别压回 L2。
  3. INSERT … ON CONFLICT (repo, scope, verdict_by) DO UPDATE，重复跑幂等。

裁定：level=2（L2 独立功能实现，独立复跑验证）
  - L2 满足：impl/entry/run_artifacts/tests 全部真实存在且复跑数值逐位一致
  - L3 不满足：无主求解器/平台侧调用方（仅同仓 CLI + 验收工具）
  - L4 不可确认：tensor_benchmark_cases 无本锚点行（基准判据唯一裁判无行可裁）；
    且验收判据 2 = 19/20（配置 4895 half_breadth_non_negative，被审方已如实披露）

用法：python3 docs/board_audit/m17577/write_verdict.py [--dry-run]
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

DB = "tensorplatform"
HERE = Path(__file__).resolve().parent
WT = HERE.parent.parent.parent  # /root/tp-agent/work/m17577
AUDIT_MID = 17577
TARGET_MID = 17472
REPO = "TensorCAD"
SCOPE = "m17472 GAP-CAD-C1-5415-PARAM-REBUILD"

REFS = "docs/board_audit/m17577/refs"
ART = f"{REFS}/artifacts/m17472_5415_param"
IMPL = f"{REFS}/impl"


def build_evidence() -> dict:
    return {
        # —— gate files_exist 要求纯路径 ——
        "impl_files": [
            f"{IMPL}/dtmb5415.py",
            f"{IMPL}/dtmb5415_cli.py",
            f"{IMPL}/tools/run_acceptance.py",
            f"{IMPL}/tools/verify_geometry.py",
            f"{IMPL}/tools/verify_shipped_model.py",
            f"{IMPL}/tools/fit_design_space.py",
        ],
        "impl_files_detail": {
            f"{IMPL}/dtmb5415.py": "746 行真实现：27 设计变量→位移场 d(v)=Φ·(C2·[v,vech(v⊗v)])→(girth,station) 样条细化→水密 mesh；解析 dV/dv；DATASET sha256 钉桩",
            f"{IMPL}/dtmb5415_cli.py": "162 行 CLI（--config-index/--design-vector/--stl/--gradient/--gate）",
            f"{IMPL}/tools/run_acceptance.py": "311 行验收入口（config0 逐点对照+20 配置几何门+梯度链）",
            f"{IMPL}/tools/verify_geometry.py": "195 行 Cb/L/B/T/∇ 对标（3% 门预置）",
            f"{IMPL}/tools/verify_shipped_model.py": "222 行 交付物 vs 独立重拟合一致性",
            f"{IMPL}/tools/fit_design_space.py": "238 行 PCA-40+二次回归拟合器（87.7 s，npz 已入库）",
        },
        "entry": f"{IMPL}/dtmb5415_cli.py（console script tensorcad-dtmb5415，pyproject.toml:43）",
        "run_artifacts": [
            f"{ART}/acceptance.json",
            f"{ART}/geometry_validation.json",
            f"{ART}/shipped_model_check.json",
            f"{ART}/shipped_model_check.csv",
            f"{ART}/sampled_configs.csv",
            f"{ART}/gate_samples.csv",
            f"{ART}/fig_config0_distance.png",
            f"{ART}/fig_geometry_validation.png",
            f"{ART}/fig_gradient_parity.png",
            f"{ART}/fig_sections_and_volume.png",
            f"{ART}/fit/design_space.npz",
            f"{ART}/fit/design_space_metrics.json",
        ],
        "tests": [
            "tests/test_dtmb5415_param.py（m17577 独立复跑：17 passed, 6.41s）",
            "test_point_grid_layout / test_cloud_map_matches_shipped_coefficients / test_shipped_model_holdout_accuracy / test_no_double_counted_linear_block / test_volume_gradient_matches_central_difference / test_principal_particulars_against_readme / test_cli_entry_point 等 17 项",
        ],
        # —— 以下两项缺失是实情，显式 null，不用文字凑数 ——
        "integration_callers": None,
        "benchmark": None,
        "repro_cmd": ("PYTHONPATH=. python3 tools/dtmb5415_param/run_acceptance.py && "
                      "PYTHONPATH=. python3 tools/dtmb5415_param/verify_geometry.py --out artifacts/m17472_5415_param && "
                      "PYTHONPATH=. python3 tools/dtmb5415_param/verify_shipped_model.py --n-configs 400 --out artifacts/m17472_5415_param && "
                      "PYTHONPATH=. python3 -m pytest tests/test_dtmb5415_param.py -q"
                      " → m17577 独立复跑：3 个 JSON 数值逐位一致（仅 wall-clock 字段不同），2 个 CSV 逐字节一致（refs/*.sha + refs/repro_*.log）"),
        "external_anchor": {
            "source": "Zenodo DTMB 5415 design-space dataset（Serani et al. 2025, Engineering with Computers）：9001 配置×27 变量，DTMB5415-database.mat 1086721690 B sha256=16b14ba5…364bee2（DATASET.source_sha256 钉桩）",
            "reference_mesh": "DTMB5415_Hull_exact.stl sha256=d40b17b1…3977fd（82794 唯一顶点）",
        },
        "absent_reason": {
            "integration_callers": ("TensorPlatform backend/api 全量 grep：DTMB5415 仅出现在 benchmark_reference*.json 数据文件，"
                                    "无主求解器/平台侧代码调用方；同仓调用方仅 CLI console script 与验收工具 → L3 不满足"),
            "benchmark": ("tensor_benchmark_cases 无 anchor_slug='GAP-CAD-C1-5415-PARAM-REBUILD' 行（m17577 实查 0 行；"
                          "全表 55 行中 CAD 锚点仅 GAP-CAD-B2-MESH-FORMAT-MATRIX / GAP-CAD-T1-BASELINE-FLEET）→ "
                          "基准判据唯一裁判无行可裁，L4 不可确认；被审方 evidence.json 的 benchmark 节为自声明，未入台账"),
        },
        "independent_reverification": {
            "repro_rerun": "m17577 于被审工作树 /root/tp-agent/work/m17472 复跑完整 repro_cmd：acceptance/geometry/shipped 三 JSON 除 wall-clock 外数值逐位一致；sampled_configs.csv 与 gate_samples.csv 逐字节一致",
            "tests_rerun": "17 passed, 6.41s（本机）",
            "report": "tensor_reports id=16450，六段齐全，19/20 几何门在标题/§3.3/§6 三处如实披露",
            "acceptance": {
                "criterion_1_config0_vs_exact": {"got": 0.9990216682368288, "threshold": 0.95, "pass": True},
                "criterion_2_sampled_gate": {"got": "19/20", "threshold": "20/20", "pass": False,
                                             "failing": "config 4895 half_breadth_non_negative（艉部中线 ~1e-6 m 越界，PCA-40 模型族限制）"},
                "criterion_3_gradient_chain": {"got": 5.2795e-08, "threshold": 0.01, "pass": True},
            },
            "key_numbers_vs_ref": {
                "Lpp_m": {"ref": 142.0, "got": 143.264947, "diff_pct": 0.891},
                "beam_m": {"ref": 18.867133, "got": 19.090662, "diff_pct": 1.185},
                "draft_m": {"ref": 6.156643, "got": 6.156643, "diff_pct": 0.0},
                "displacement_volume_m3": {"ref": 8412.037925, "got": 8452.75781, "diff_pct": 0.484},
                "cb": {"ref": 0.509992, "got": 0.501988, "diff_pct": -1.569},
                "holdout_rms_m": {"ref_independent_refit": 0.25188, "got": 0.25190, "diff_pct": 0.0099},
            },
        },
        "material_defects": [
            "tensor_benchmark_cases 无本锚点行：ref/tolerance 判据未入台账，L4 无法经基准台账确认（被审方未登记）",
            "验收判据 2「抽样 20 配置过几何质量门」= 19/20 未达标（config 4895，被审方已如实披露，非隐瞒）",
            "无主求解器/平台侧调用方，L3 不满足（仅同仓 CLI+验收工具）",
            "被审分支 agent/m17472（commit 979fb83）尚未并入 origin/main（mission 状态 succeeded，合并待门禁）",
        ],
        "honest_disclosure": "被审方报告对 19/20 与模型有效域（数据集流形邻域）均如实披露，无粉饰；缺陷在台账登记与判据 2 未闭环",
        "rework_recommendation": ("① 向 tensor_benchmark_cases 登记本锚点行（ref=Zenodo README×λ 主尺度参数 + "
                                  "DTMB5415_Hull_exact.stl sha256，tolerance 按报告 §6 口径 3%）；"
                                  "② 修复或正式收窄判据 2（提高 PCA rank / 流形投影约束，或把判据限定在数据集流形邻域）；"
                                  "③ 完成后重审 L4"),
    }


def insert(evidence: dict, dry_run: bool = False) -> str:
    sql = f"""
INSERT INTO tensor_maturity_verdicts (repo, scope, mission_id, level, evidence, verdict_by)
VALUES ('{REPO}', '{SCOPE}', {AUDIT_MID}, 2, $$ {json.dumps(evidence, ensure_ascii=False)} $$::jsonb, 'board-audit')
ON CONFLICT (repo, scope, verdict_by) DO UPDATE
  SET level = EXCLUDED.level, evidence = EXCLUDED.evidence,
      mission_id = EXCLUDED.mission_id
RETURNING id, coalesce(gate_state,''), coalesce(gate_level,-1);
"""
    sql_path = "/tmp/m17577_verdict.sql"
    Path(sql_path).write_text(sql, encoding="utf-8")
    if dry_run:
        print(sql[:400], "...")
        return "dry-run"
    r = subprocess.run(["su", "postgres", "-c",
                        f"psql -d {DB} -t -A -F '|' -f {sql_path}"],
                       capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        raise SystemExit(f"insert failed: {r.stderr[:600]}")
    line = [ln for ln in r.stdout.strip().splitlines() if ln.count("|") >= 2][-1]
    vid, state, glevel = line.split("|")[:3]
    return f"verdict_id={vid} gate_state={state or 'pending'} gate_level={glevel}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    ev = build_evidence()
    (HERE / "verdict_evidence.json").write_text(
        json.dumps(ev, ensure_ascii=False, indent=2), encoding="utf-8")
    print(insert(ev, dry_run=a.dry_run))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
