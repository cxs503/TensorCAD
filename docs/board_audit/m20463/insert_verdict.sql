-- Board audit verdict for m20389 REPRO-GEOM-FLEET-REFS
-- Audited by: board-audit (m20463)
-- Date: 2026-11-14

INSERT INTO tensor_maturity_verdicts (repo, scope, mission_id, level, evidence, verdict_by) VALUES (
  'TensorCAD',
  'm20389 REPRO-GEOM-FLEET-REFS',
  20463,
  2,
  $json${
  "anchor": "REPRO-GEOM-FLEET-REFS",
  "audited_mission": 20389,
  "audited_worktree": "/root/tp-agent/work/m20389 branch=agent/m20389 commit=1f5763d",
  "impl_files": [
    "tensorcad/fleet/crossref.py:1-454 (Cross-source reference library: parse_source/cross_reference/closure_check/coverage_summary; 真代码：HTML表格解析→家族谱系分类→交叉验证判定→闭锁自检)",
    "tools/baseline_fleet/crossref_refs.py:1-370 (fetch子命令：从出版商页面抓取HTML并SHA256钉扎)",
    "tools/baseline_fleet/crossref_report.py:1-318 (报告渲染：从crossref.json生成markdown+matplotlib图)",
    "benchmarks/verified/fleet_refs/crossref.json:1-1588 (完整交叉参考库：6船×7量=42格，含sources/report/coverage/comparisons/closure/probes)"
  ],
  "entry": "python tools/baseline_fleet/crossref_refs.py build --raw artifacts/fleet_refs/raw && python tools/baseline_fleet/crossref_report.py --library benchmarks/verified/fleet_refs/crossref.json --doc docs/fleet_refs_crossref.md",
  "run_artifacts": [
    "benchmarks/verified/fleet_refs/crossref.json (1588行，schema=tensorcad.fleet_refs.crossref/1.0，anchor=REPRO-GEOM-FLEET-REFS)",
    "docs/fleet_refs_crossref.md (115行，六段报告：介绍/方法/参数表/交叉验证/灵敏度/失败模式)",
    "docs/figs/fleet_refs/closure_check.png (Cb闭锁自检图)",
    "docs/figs/fleet_refs/coverage_matrix.png (42格覆盖矩阵)",
    "docs/figs/fleet_refs/displacement_ref_vs_got.png (排水量ref vs got对比)"
  ],
  "tests": {
    "files": ["tests/test_fleet_refs_crossref.py:1-201 (14 tests)"],
    "audit_run": "python -m pytest tests/test_fleet_refs_crossref.py -q → 12 passed, 2 failed, 0.46s",
    "failures": [
      "test_fallback_column_is_flagged_and_does_not_vote: 期望fallback列标记后verdict=MISSING，实际得SINGLE_SOURCE——fallback列的cell虽被标记unit_basis='fallback_column:...'但cross_reference()未排除其投票",
      "test_conflict_is_reported_not_averaged: 期望family_values键为family='F'，实际得source key='a'——cross_reference()用source key而非family作为family_values的键"
    ],
    "passing": [
      "test_label_comes_from_the_first_cell_only",
      "test_mass_reference_is_converted_with_the_stated_density",
      "test_two_documents_from_one_lineage_are_not_two_sources",
      "test_every_source_is_digest_pinned_the_same_way_after_a_fetch",
      "test_cross_verified_cells_have_two_lineages_and_agree",
      "test_kcs_is_fully_cross_verified",
      "test_kvlcc2_collapses_to_one_lineage_and_exposes_the_bwl_typo",
      "test_honest_gaps_are_recorded",
      "test_missing_reconstructions_are_not_quoted_as_errors",
      "test_closure_uses_each_sources_own_digits",
      "test_report_is_rendered_from_the_library",
      "test_registry_covers_all_six_hulls"
    ]
  },
  "integration_callers": "none beyond tests: grep -rn 'fleet.crossref\\|fleet_refs' --include='*.py' 仅见 tests/test_fleet_refs_crossref.py 和 tools/baseline_fleet/crossref_report.py 调用 → L3 不达",
  "benchmark": {
    "anchor": "REPRO-GEOM-FLEET-REFS",
    "registered": "SELF_DECLARED — anchor在crossref.json:3和crossref.py:3自声明",
    "coverage_total_cells": 42,
    "coverage_cross_verified": 11,
    "coverage_single_source": 13,
    "coverage_conflict": 2,
    "coverage_missing": 16,
    "ships_fully_cross_verified": ["KCS"],
    "ships_failing_min_two_sources": ["DTC", "JBC", "KVLCC2", "Series60"],
    "key_metrics": {
      "KCS_volume_error": "+0.30% (VERIFIED, within 1% tolerance)",
      "KVLCC2_volume_error": "WITHHELD (geometry reconstruction failed: fewer than two usable waterline slices)",
      "JBC_volume_error": "-18.28% (WITHHELD, outside 1% tolerance)",
      "ONRT_volume_error": "WITHHELD (model-scale bytes, not comparable to full-scale displacement)",
      "DTC": "MISSING (no reachable numeric source)",
      "Series60": "MISSING (no reachable numeric source)"
    },
    "honest_gaps_documented": true,
    "kvlcc2_bwl_typo_detected": true,
    "onrt_mass_vs_volume_conflict_detected": true
  },
  "repro_cmd": {
    "documented_command": "python tools/baseline_fleet/crossref_refs.py fetch --out artifacts/fleet_refs && python tools/baseline_fleet/crossref_refs.py build --raw artifacts/fleet_refs/raw && python tools/baseline_fleet/crossref_report.py --library benchmarks/verified/fleet_refs/crossref.json --doc docs/fleet_refs_crossref.md",
    "audit_verification": "crossref.json已存在且结构完整(42格/6船/7量)；fetch需网络未重跑；build+report从已有JSON可复现文档",
    "tests_reproducible": "12/14 pass; 2 fail是实现bug非环境问题"
  },
  "external_anchor": "出版商原始网页(NMRI T2015/SIMMAN 2014/2020/ITTC)的SHA256钉扎字节；Cb·Lpp·Bwl·T闭锁自检为独立数学验证",
  "report_six_sections": {
    "section_1_intro": "✅ 物理背景+六标模+参数表",
    "section_2_method": "✅ 方法/流程/工具链",
    "section_3_results": "✅ 数字表(42格verdict+4船ref vs got+闭锁自检)",
    "section_4_figures": "✅ 3张matplotlib图(closure_check/coverage_matrix/displacement_ref_vs_got)",
    "section_5_caveats": "✅ 灵敏度分析+失败模式+诚实差距",
    "section_6_error_judgment": "✅ KCS +0.30%<3% ✅; KVLCC2/JBC/ONRT/DTC/Series60 明确声明WITHHELD/MISSING原因"
  },
  "level_rationale": "L2：独立功能实现(真代码：HTML解析→谱系分类→交叉验证→闭锁自检→报告渲染，1142行核心+工具代码)；42格参考库结构完整；报告六段齐全；3图程序化生成。不达L3：无求解器/平台集成调用者(仅tests和tools调用)。不达L4：2/14测试失败(实现bug：fallback列投票未排除、family_values键用source key而非family)；仅1/6船(KCS)满足≥2独立来源要求；3/6船benchmark WITHHELD。",
  "verdict_summary": "三件套核验：run_artifacts ✅(crossref.json 1588行+report.md 115行+3图全存在)；ref vs got 部分✅(KCS +0.30%在容差内，但KVLCC2/JBC/ONRT WITHHELD，DTC/Series60 MISSING)；repro_cmd 部分✅(文档命令完整，但2/14测试失败暴露实现bug)。报告六段✅齐全。诚实度高：主动暴露KVLCC2 Bwl=580排版错误、ONRT质量/体积混淆、DTC/Series60无源。判级L2：功能实质为独立数据管线，测试失败和覆盖率不足阻止L4。",
  "gate_notes": "mission声称succeeded但实际仅KCS满足≥2独立来源要求(1/6)；2测试失败未修复即交差；maturity_gate应reopened或confirmed-with-notes"
}$json$::jsonb,
  'board-audit'
) ON CONFLICT (repo, scope, verdict_by) DO UPDATE
  SET mission_id = EXCLUDED.mission_id,
      level = EXCLUDED.level,
      evidence = EXCLUDED.evidence,
      gate_state = 'pending',
      gate_level = NULL,
      gate_note = NULL
RETURNING id, repo, scope, level, gate_state;
