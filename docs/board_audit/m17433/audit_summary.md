# m17433 交差核验单 · GAP-CAD-B3-SOLVER-INTERFACE（被审 m17403 · TensorCAD）

**结论：L3（系统集成）**，不达 L4。verdict 落在 `tensor_maturity_verdicts`（repo=TensorCAD,
scope=`m17403 GAP-CAD-B3-SOLVER-INTERFACE`, mission_id=17433, verdict_by=board-audit）；报告 `tensor_reports` id=16427。

## 三件套核验（基准判据 = tensor_benchmark_cases，唯一裁判）

| 项 | 结论 | 证据 |
|---|---|---|
| run_artifacts | 存在，可重算；但 BEM 工件 2/3 档、summary.json 只剩 fvm（F4） | `runs/m17433/subject_artifacts/`；审计重跑 `runs/m17433/repro_artifacts/` |
| ref vs got | 逐项独立复算一致（FEM 体积 2e-14、反力 0、BEM heave 1.13%、FVM L2 0.63%） | `runs/m17433/independent_recompute.json` |
| repro_cmd | 实跑通（fem/bem/fvm 全部 exit 0，关键数字逐位一致） | 见 `runs/m17433/repro_cmd.sh` |
| 基准案号 | **0 行**（`WHERE anchor_slug='GAP-CAD-B3-SOLVER-INTERFACE'`）→ L4 无机器裁判 | SQL 见报告 §2 |
| 报告六段 | 被审单未发报告（result 自报未 POST）；本审计报告六段齐全 | `tensor_reports` id=16427 |

## 判级依据

- **L2 已达**：真代码 输入→计算→输出（`solver_interface.py` 553 行：z-up 契约 / gmsh C3D4 写入 / 独立解析器回读）；
  无 MOCK/STUB/random 标记；审计独立复算（不 import 被测模块）逐位吻合。
- **L3 已达（本单判定）**：三个主求解器（TensorFEM `Tet4J2Model` / TensorBEM `compute_added_mass_matrix_3d` /
  TensorFVM `BodyFittedSolver`）真实消费 CAD bundle 并端到端跑通，帧/单位契约机器校验。
- **L4 未达**：(a) 无 `tensor_benchmark_cases` 案号（F2）；(b) 全部交付 untracked、branch/head_commit=null（F1）；
  (c) FEM 头条指标 `patch_test_relative_linf` 名不符实（F3，见下）；(d) 工件陈旧（F4）；(e) 锚点“水密面→FVM 体网格”
  规范路径未演示：FVM 腿为 2D 槽道、FEM 腿为盒体、BEM 面由接口模块内建（F5）；(f) 无 pytest/文档/API 导出（F6）。

## 最重要发现（F3，可复现的对照实验）

被审交付把 `patch_test_relative_linf` 当作 patch test 收敛证据（细网格 0.72%）。控制实验：同网格 nx7×nz5、同载荷、同约束，
仅令 `poisson=0` → 相对误差 **1.396e-04**，而 `poisson=0.3` 同网格 **13.7474%**（比值 **985×**）。
真 patch test 应在任何协调四面体网格上复现线性场到机器精度 → 该指标度量的是“单角节点钉死 y/z 与泊松收缩场不相容”造成的边界层，
不是离散误差。真正可信的度量是平衡残差 1.86e-13 与体积/拓扑独立复算（0.24 m³ rel 2.16e-14，边界 496 面/744 边全用两次，χ=2）。

## 复现

```
cd /root/tp-agent/work/m17433 && python3 runs/m17433/independent_recompute.py   # 独立复算（不 import 被测）
python3 runs/m17433/make_figure.py                                              # 图 → artifacts/m17433/*.png
python3 docs/board_audit/m17433/publish.py                                      # 落表 + 发报告（幂等）
```

被审代码复现（审计在 /tmp/aud17433 副本上进行，未污染 m17403 工作树）：
`PYTHONPATH=.:/root/tensor-ecosystem/TensorFVM/src python examples/solver_interface/run_three_solvers.py --only {fem|bem|fvm}`
（`--only` 会覆盖 summary.json，这是被审 summary 不完整的原因）。

## 文件清单

- `runs/m17433/independent_recompute.py` / `.json` —— 独立复算脚本与结果（零手抄）
- `runs/m17433/make_figure.py` —— 三联图生成
- `runs/m17433/repro_artifacts/` —— 审计重跑工件（含 `fem_tension/nu0_control_results.json` 控制实验）
- `runs/m17433/subject_artifacts/` —— 被审工件快照（含 digest，供后续比对）
- `runs/m17433/subject_impl_solver_interface.py` / `subject_impl_run_three_solvers.py` —— 被审实现快照（sha256 见 evidence）
- `runs/m17433/repro_cmd.sh` —— 复现命令
- `artifacts/m17433/solver_interface_audit.png` —— 报告图
- `docs/board_audit/m17433/publish.py` —— 落表 + 报告发布（幂等）

## 独立性限制

审计在同一物理机上、用被审交付文档中的命令复现；被审模块未被 import 到复算脚本中，但**求解器侧**（TensorFEM/TensorBEM/TensorFVM）
与运行环境是共用的。未被审查：TensorFEM/TensorBEM/TensorFVM 自身的正确性（另有其仓的审计面）。