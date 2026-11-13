# Board audit — m17476 / GAP-CAD-B3-SOLVER-INTERFACE（交差核验单，非重做）

- **被审单**：m17476（TensorCAD，返修轮1，status=succeeded，agent branch `agent/m17476`，base `16eefd3`）
- **审计单**：m17576（repo=TensorCAD，board-audit）
- **审计对象**：worktree `/root/tp-agent/work/m17476`（交付物在审计时为**未提交**工作树：`?? tensorcad/solver_interface.py`, `?? tests/test_solver_interface.py`, `?? docs/solver_interface.md`, `?? artifacts/m17476/`, `M pyproject.toml tensorcad/__init__.py tensorcad/benchmark.py`）
- **锚点验收判据**（m17476 goal 原文）：三求解器各一个端到端算例真跑通（坐标系契约统一 z-up）；输入=统一几何库输出
- **审计结论**：**L4（基准验证）**，附 6 条 finding（无阻断级缺陷，2 条需返修叙事/工程收口）

## 1. 取证与复核方法（全部可机器复跑）

审计器：`docs/board_audit/m17576/verify_m17576.py`（静态证据重算 + 契约行为探针 + 工件算术复算；`--pytest` 复跑测试，`--repro` 复跑 repro_cmd 并与出厂工件逐项比对）。
本次原始输出：`docs/board_audit/m17576/audit_m17576.json`；复跑日志：`docs/board_audit/m17576/repro_auditor_run.log`；被审单 result/report 快照：`mission_17476_result.json` / `report_16459.md`。

| 核验项 | 判定 | 真实指针 |
|---|---|---|
| impl_files（非占位） | ✅ | `tensorcad/solver_interface.py` 1854 行 / 75759 B；占位扫描（TODO/FIXME/NotImplemented/stub）**0 命中**；`run_bem_added_mass:929` / `run_fem_plate:1160` / `run_lbm_channel:1393`；契约 `FRAME_CONTRACT_ID:76` |
| entry（入口可执行） | ✅ | `tensorcad/solver_interface.py:1824 main()` → `PYTHONPATH=. python -m tensorcad.solver_interface --out <dir>`；包导出 `tensorcad/__init__.py:74`（11 个符号） |
| run_artifacts（真实工件） | ✅ | `artifacts/m17476/solver_interface/` 16 个文件：`results.json`（29 comparisons / statuses: 3×`ran` / `gates_failed=[]`）、`bem/{panels,manifest}_{16,24}.json`、`fem/plate_{8x8,16x16}.json`、`lbm/{lattice,voxels}_{8,16}.json`、`convergence.png`、`lbm_profile.png` |
| tests | ✅ | `tests/test_solver_interface.py`，审计员复跑 **11 passed / 21.77 s**；断言的是不变量（错坐标系拒收、错单位拒收、闭合外向球、加密单调、digest 稳定性、体素盒精确、面板闭合残差<1e-12），不是「跑过即过」 |
| integration_callers（真调用） | ✅（仓内）/ ⚠️（平台侧 0） | `tensorcad/benchmark.py:2260 solver_handoff_case`（case C）在 `:2283-2285` 调用三个 run_*；真求解器 API 调用：`tensorbem.TriangleMesh:986` + `compute_added_mass_matrix_3d:989`、`tensorfem.plate.q4_mindlin_stiffness:1234` + `kirchhoff_sine_center_deflection:1185`、`tensorlbm.equilibrium/collide_bgk/macroscopic:1470-1491`。`grep -r solver_interface TensorPlatform/backend` 仅命中报告文本载荷，**无平台代码调用方** |
| benchmark（ref vs got） | ✅ 数值在容差内 / ⚠️ 无台账登记 | 24 条带容差门禁 **全部 passed**（`gates_failed=[]`）；artifacts 里 rel_error 经审计器逐条用 `measured/reference-1` 独立复算，**0 条不符**；但 `SELECT ... FROM tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-B3-SOLVER-INTERFACE'` = **0 行**（全表 54 行，TensorCAD 仅 2 行别的锚点）→ 见 F2 |
| repro_cmd（可执行） | ✅ 位对齐复现 | 审计员 2026-11-13 07:01–07:03 真跑 `PYTHONPATH=. python -m tensorcad.solver_interface --out /tmp/repro_m17576/full`：**WALL 245 s**，`ran=['tensorbem','tensorfem','tensorlbm']`、`unavailable=[]`、`gates_failed=[]`；与出厂 `results.json` 逐项比对 **29/29 相同、最大相对差 0.0、门禁翻转 0** |
| external_anchor | ✅ | 闭式解：球附加质量 `2/3 π ρ a³`（Lamb）、`V=4/3πa³`、`A=4πa²`；Kirchhoff 板中心挠度 `q0 a⁴/(4π⁴D)`；体力 Poiseuille `u_max = F H²/(8ν)` |

## 2. 基准裁判逐项（ref vs got，<3% 判定）

| 量（细网格） | ref | got | 误差 | 容差 | 判 |
|---|---|---|---|---|---|
| BEM 24×48 体积 | 0.523599 m³ | 0.519871 | −0.712 % | 1 % | ✅ |
| BEM 24×48 湿面积 | 3.141593 m² | 3.130391 | −0.357 % | 1 % | ✅ |
| BEM 24×48 附加质量(垂荡) | 261.799 kg | 260.296 | −0.574 % | 5 % | ✅ |
| BEM 24×48 附加质量(纵荡) | 261.799 kg | 260.933 | −0.331 % | 5 % | ✅ |
| FEM 16×16 中心挠度 | 5.12478e-4 m | 5.13181e-4 | +0.137 % | 5 % | ✅ |
| FEM 8×8 中心挠度 | 5.12478e-4 m | 5.09827e-4 | −0.517 % | 5 % | ✅ |
| LBM 16 峰值速度 | 6.375e-3 | 6.352e-3 | −0.361 % | 5 % | ✅ |
| LBM 8 峰值速度 | 1.575e-3 | 1.552e-3 | −1.460 % | 5 % | ✅ |
| LBM 剖面形状（细） | 抛物线 | 抛物线 | 3.76e-9（绝对） | 1e-6 | ✅ |
| CAD→体素单元数 16 | 4096 | 4096 | 0 | 0 | ✅ |

**是否 <3%：是**。所有被判量最大误差 1.460 %（出现在最粗的 LBM 8 单元，属收敛阶梯上级）；细网格最大 0.712 %。
辅助不变量：BEM 闭合残差 3.9e-17、双层面行和误差 8.4e-15；LBM 质量守恒 3.7e-13、对称性 1.5e-14。

## 3. 判级：L4

| 级 | 判据 | 结论 |
|---|---|---|
| L0 文档 | — | 超额 |
| L1 数学/单元验证 | 闭合/拓扑/收敛不变量 | ✅ 11 测试断言不变量 |
| L2 独立功能实现 | 真代码 输入→计算→输出 | ✅ 1854 行非占位 + CLI 入口 + 工件 |
| L3 系统集成 | 接入主求解器/平台 | ✅ 三个真求解器库 API 被真调用（三个不同仓），仓内 benchmark case C 为主调用方；⚠️ 无平台/船型上层调用方 |
| **L4 基准验证** | 完整算例+误差分析+可复现 | ✅ 三个端到端算例 × 闭式解参考 × 两分辨率收敛 × 24 门禁全绿 × 审计员位对齐复现（245 s） |
| L5 科研工程成果 | 新方法/新物理结论 | ✗ 不主张：算例为标准闭式解，属接口工程而非科研结论 |

不支持 L5 的理由：三个算例（球、方板、直通道）是标准验证题，不是船舶阻力/耐波新结论；LBM 为 2D x-z 平面，报告自己声明不宣称阻力预报。

## 4. Findings（按严重度）

- **F1 叙事夸大（medium，要求文档更正）**：`docs/solver_interface.md` 写「The files are the interface: the runs read what a downstream consumer would read」，`run.json` note 写「the geometry reaches the panel code as a written JSON table, not as an in-memory object」；但代码实际把**内存**数组交给求解器：`solver_interface.py:986` `bem.TriangleMesh(np.asarray(bundle.mesh.vertices), ...)`。JSON 是**留痕**不是**接口**（除 `results.json:1773` 增量合并外，全文无回读 `panels_*.json` / `plate_*.json` 的代码）。功能不受损（求解器真跑），但契约叙事强于实现；且 m17403 版本里真正的回读校验 `verify_bundle()`（main `tensorcad/solver_interface.py:495`）在本次重写中**丢失**。
- **F2 基准无台账登记（medium）**：`tensor_benchmark_cases` 对本锚点 0 行 → 机器裁判只能退回仓内 `tensorcad/benchmark.py:case_c_solver_interface` 的门禁与闭式解参考。与 m17414（TensorBEM L4）同款 finding 同款处置：登记 case_c 三行（case_name/metric/tolerance）可闭合。
- **F3 门禁分辨率事后选取（low-medium）**：16×32 体积误差 −1.596 % 超 1 % 容差，最终把 gated 行定为**最细分辨率**后通过（m17476 result checks 自认「已按 24×48 实测 −0.712 % 调至通过」）。缓解证据：粗网格失败数与 O(h²) 阶梯**照实公布**，且另有一条尺度无关的 `refinement_reduces_the_error` 门禁（不看绝对值）独立看收敛方向 → 不构成造假，但属「看结果调门」痕迹，后续若改容差须留 before/after。
- **F4 与 main 的同路径冲突（medium，工程收口）**：m17476 base=`16eefd3`，**早于** `b33a3a2`（把 m17403 的另一个 `tensorcad/solver_interface.py` 并进 main，553 行，API=`write_bundle_volume/verify_bundle/to_z_up`）。两份同名模块 API 不相交、无合并决议；直接合入会冲突并静默丢掉 m17403 的写盘+回读校验能力。且交付物**未提交**（工作树 `??`），自动推分支将推不出任何东西 → 需 commit 并 rebase 到含 `b33a3a2` 的 main 后重解冲突。
- **F5 环境未钉住（medium，可复现性边界）**：三个求解器由三个**别的 agent worktree** 的 editable 安装提供（`solver_file` 记录 m17241 / m17361 / m17457），results.json 只记 `__file__`+version（TensorFEM/TensorLBM 的 version 为 null）；审计时同一主机默认导入已漂到 `/root/tp-agent/work/m17591/TensorLBM`（本模块自身的 `import tensorcad` 也漂到 m17333，故 repro_cmd 必须带 `PYTHONPATH=.`）。**位对齐复现成立，但属主机绑定**：干净机器上会 `status=unavailable`、退出码 2（诚实状态，非假通过）。建议记录三个求解器仓的 commit SHA。
- **F6 跨仓发现已复核（informational）**：交付物称 `tensorlbm.stream` 方向与自身 `d2q9.C` 反号（+x 布居被搬向 x−1），因此自写半程反弹 streaming、只借 TensorLBM 的平衡/碰撞/宏观量。审计员独立探针（另一 checkout `m17591`）：`f[1,4,4]=1 → stream → 位置 (4,3)`，**x 递减，复现其说法**。后果：LBM 链路是「TensorLBM 内核 + 本地 streaming」混合，「三求解器真跑通」在这条链上是**部分真跑通**，报告已在 ⑥ 段如实披露，未粉饰。

## 5. 未跑项（如实记 not_run，不记为通过）

1. `case_c_solver_interface` 的逐门禁明细：case C 曾完整跑通（`repro_auditor_run.log`：3/3 solvers ran，约 148 s），但该单的 result 里 `gates_failed=[]` 来自 CLI 主路径；审计员捕获 `CaseResult` 无 `passed` 属性（类定义 `benchmark.py:232` 只有 `gates/failed_gates`），逐门禁明细**未捕获**。
2. 仓库全量 pytest（只复跑 `tests/test_solver_interface.py`；被审单自报全量未跑）。
3. README/CHANGELOG 的 Unreleased 条目（被审单自报 patch 处 AssertionError，未写入）——审计确认 `docs/solver_interface.md` 存在但 README 未提及该能力。
4. 跨架构（x86/GitHub CI）复现：本审计仅在 loongarch64 主机同环境进行。
5. 船舶域算例（船体→求解器）：本锚点交付的是接口层，三个算例为球/板/通道，**未**用 KCS/KVLCC2 等船体几何；这属于 GAP 的下一层而非本单验收线。

## 6. 处置（disposition）

**Approve with notes（L4 成立，不返修功能；三项收口建议）**：①改正 F1 的「文件即接口」叙事（或补真正的回读路径，把 m17403 的 `verify_bundle` 合并回来）；②F2 把 case_c 三行登记进 `tensor_benchmark_cases`；③F4 提交交付物并 rebase 解决与 `b33a3a2` 的同路径冲突（显式保留/合并两版能力），F5 补记三仓 commit SHA。

## 7. 复跑方式

```bash
python3 docs/board_audit/m17576/verify_m17576.py                # 静态重算，秒级
python3 docs/board_audit/m17576/verify_m17576.py --pytest       # + 测试复跑
python3 docs/board_audit/m17576/verify_m17576.py --repro        # + repro_cmd 位对齐复现（~250 s，后台跑）
```