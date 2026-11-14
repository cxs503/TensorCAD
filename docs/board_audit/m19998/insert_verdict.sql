INSERT INTO tensor_maturity_verdicts (repo, scope, mission_id, level, evidence, verdict_by) VALUES ('TensorCAD', 'm19974 CASE-CAD-VOLUME-CLOSED', 19998, 4, $json${
  "anchor": "CASE-CAD-VOLUME-CLOSED",
  "audited_mission": 19974,
  "audited_worktree": "/root/tp-agent/work/m19974 branch=agent/m19974",
  "impl_files": [
    "tensorcad/mesh_volume.py:1-450 (analytic_sphere_volume/analytic_box_volume/check_closed_body/convergence_study/observed_order/defect_detection; 真代码：解析闭式→闭体检查→收敛阶梯→负控，无硬编码结果)",
    "tensorcad/geometry_audit.py:346-500 (audit_watertightness: 自由边/非流形边/重复面/绕向冲突/有符号体积/缝距)",
    "benchmarks/verified/mesh_volume_watertight/run_benchmark.py:1-450 (判据组装+工件输出+环境记录)"
  ],
  "entry": "python benchmarks/verified/mesh_volume_watertight/run_benchmark.py --out runs/mesh_volume_watertight",
  "run_artifacts": [
    "benchmarks/verified/mesh_volume_watertight/result.json (25KB, 32 gates, anchor=CASE-CAD-VOLUME-CLOSED, passed=True)",
    "benchmarks/verified/mesh_volume_watertight/convergence.csv (7 levels, h/V/error/open_edges/euler)",
    "benchmarks/verified/mesh_volume_watertight/convergence.png (log-log error vs h, slope-2 reference, 0.1% line)",
    "benchmarks/verified/mesh_volume_watertight/gates.png (32 gates, green/red pass/fail)",
    "benchmarks/verified/mesh_volume_watertight/box_stl.stl (float32 round-trip sample)",
    "runs/mesh_volume_watertight/repro.log (full run output, 32/32 ✓, 16.074s)"
  ],
  "tests": "tests/test_mesh_volume.py: 12/12 passed (24.78s): 断言均为不变量（rel_tol/exact_tol/route_tol），非仅执行",
  "integration_callers": "none：grep -rn 'mesh_volume' --include='*.py'（排除模块自身与其测试）仅见 benchmarks/verified/mesh_volume_watertight/run_benchmark.py 与 tests/test_mesh_volume.py 调用 → L3 不达",
  "benchmark": {
    "anchor": "CASE-CAD-VOLUME-CLOSED",
    "registered": "SELF_DECLARED — 案号在 run_benchmark.py:419 与 result.json:3 自声明；tolerances 预声明于 tensorcad/mesh_volume.py 模块常量",
    "gates_total": 32,
    "gates_passed": 32,
    "gates_failed": 0,
    "key_metrics": {
      "box_volume_error": "0.000e+00 (exact)",
      "sphere_sub6_volume_error": "-1.352e-04 (-0.0135%, within 0.1%)",
      "observed_order_volume": 1.9888,
      "observed_order_area": 2.0051,
      "negative_controls_detected": 4
    }
  },
  "repro_cmd": "cd /root/tp-agent/work/m19974 && python benchmarks/verified/mesh_volume_watertight/run_benchmark.py (本单复算: 32/32 ✓, 16.074s); python -m pytest tests/test_mesh_volume.py -v (12/12 passed)",
  "external_anchor": "解析闭式（Euclid《几何原本》XII.18；立体几何 V=4/3πr³, A=4πr², V=abc）；trimesh 5.0.0 第三方交叉",
  "report": "docs/mesh_volume_watertight_bench.md 六段齐全；所有数字由 publish_report.py 从 result.json 程序化注入，零手抄",
  "level_rationale": "L4：完整基准算例（32 条判据全过）+ 误差分析（收敛阶梯 p=1.989/2.005、尺度不变性、刚体运动不变性）+ 可复现（repro_cmd 验证通过）。不达 L3：无求解器/平台集成调用者。不达 L5：验证基准，非科研工程成果。",
  "verdict_summary": "三件套：run_artifacts ✅（result.json/convergence.csv/convergence.png/gates.png/box_stl.stl/repro.log 全存在）；ref vs got ✅（32/32 判据全过，最差 -0.0135% < 0.1%）；repro_cmd ✅（本单实跑 32/32 ✓, 12/12 tests ✓）。报告六段 ✅齐全。判级 L4。"
}$json$::jsonb, 'board-audit') ON CONFLICT (repo, scope, verdict_by) DO UPDATE SET level=EXCLUDED.level, evidence=EXCLUDED.evidence, mission_id=EXCLUDED.mission_id, gate_state='pending', gate_level=NULL, gate_note=NULL RETURNING id, level, gate_state;
