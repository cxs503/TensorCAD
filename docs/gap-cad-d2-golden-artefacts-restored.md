# GAP-CAD-D2 金标回归门禁：committed 运行工件恢复（m22631）

## 问题

锚点 `GAP-CAD-D2-GOLDEN-REGRESSION`（几何金标回归套件+CI）的**门禁自测套件**
`tests/test_geometry_golden_gate.py` 在本次 HEAD（`8023726`）上是红的：

```
FAILED test_assert_report_accepts_the_committed_run
FAILED test_assert_report_rejects_a_case_that_stopped_being_measured[...]  (×4)
FAILED test_assert_report_rejects_lost_detection_power[...]              (×5)
FAILED test_assert_report_runs_as_the_gate_invokes_it
FAILED test_gate_document_is_generated_from_the_committed_artefacts
11 failed, 11 passed
```

根因：`tools/golden/gap_d2_assert_report.py`（CI 与 pre-commit 钩子共用的**唯一判据**）与
`tools/golden/gap_d2_gate_docs.py --check`（门禁文档生成器）都从
`artifacts/gap-cad-d2-golden/*.json` 读取**本次运行的 committed 工件**（`golden_suite.json`、
`mutation_report.json`、`run_meta.json`、`series60_convergence.json`、`hook_probe.json`、
`regression_control.json`）。该目录在安装门禁的提交 `3a79b8a` 里存在，但**不在当前 HEAD 的
树里**（`artifacts/` 只剩 `m3768`），于是：

- 判据脚本读不到工件 → `FileNotFoundError` → 门禁自测全红；
- 文档生成器 `--check` 因工件缺失而判陈旧 → 文档与工件脱钩。

这不是"门禁没装"，而是**门禁的自我验证欠账**：门禁本身（CI 工作流 + 钩子 + golden.py）都在，
但"门禁判据所依赖的 committed 运行工件"没随当前分支落地，导致门禁自己的测试套件无法证明
自己是绿的。

## 修复

从安装门禁的提交 `3a79b8a`（`geometry-golden: install the GAP-CAD-D2 dual-derived golden gate`）
把 `artifacts/gap-cad-d2-golden/` 的 12 个 committed 运行工件原样恢复到当前工作树：

```
artifacts/gap-cad-d2-golden/
  golden_suite.json                          # 8 金标 case 的 dual-derived 判定（EXACT 1e-9 / QUAD 1e-6）
  mutation_report.json                       # 突变审计：控制臂 8/8 PASS、检测力 1.0、WITHHELD 仍拒判
  run_meta.json                              # 运行元数据（退出码/墙钟/快照）
  series60_convergence.json                  # Series60 默认网格 Cb 欠填收敛记录
  hook_probe.json                            # 钩子七种决策模式实测
  regression_control.json                    # +3% Wigley 半宽 / +2% 水线面回归控制
  *.log                                      # suite.py / mutants.py 运行日志
```

这些工件是**第三方判据**（TensorPapers `papers/gap-cad-d2-golden-regression` 的 dual-derived
suite）的 committed 快照——判据代码与金标数值都不住在本仓，本仓只暴露 `TENSORCAD_PATH`，
因此恢复它们不会把门禁"bless 成与自己一致"。

## 验证

- `python -m pytest tests/test_geometry_golden_gate.py -q` → **22 passed**（此前 11 failed）；
- `python -m tensorcad.golden --check --json /tmp/golden_full.json` → `passed: True`，
  4 个金标 case（wigley/series60/sphere_ellipsoid/dtmb5415）全部在 committed 基线内；
- `python -m pytest tests/test_wigley.py tests/test_series60.py -q` → **50 passed**；
- 钩子行为实测：
  - 仅暂存 `docs/*.md` → `no staged geometry change, check skipped`（exit 0，不启动套件）；
  - 暂存 `tensorcad/hull.py` + `GOLDEN_DIR=/nonexistent-suite` → `refusing to guess`（拒绝，不静默降级）。

## 误差判定（对拍判据零手抄）

本单是**门禁自测套件修复**，不是新增几何量。门禁判据本身由 `gap_d2_assert_report.py` 从
committed 工件读取（零手抄），恢复后自测全绿。金标 case 的实测误差（如 Wigley Cb −2.36e-03、
Series60 Cb −1.77e-02、球/椭球体积 −3.61e-03）均落在各自**由两档网格变化推导**的容差带内
（`<3%` 判据由容差推导机制保证，非手敲），详见 `docs/golden/tensorcad_geometry_golden.json`。
