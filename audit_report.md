# AUDIT-ASSET-12REPO · TensorCAD 质量审计报告

- 审计任务：m14491（锚点 AUDIT-ASSET-12REPO，A+B 双面：资产盘点 + 功能真实性 + 成熟度判级）
- 审计对象：`cxs503/TensorCAD` main @ `8c5131b`（含近期已并 agent PR：#295–#311 等）
- 审计者：audit2-ab（verdict_by），工作树 `/root/tp-agent/work/m14491`，分支 `agent/m14491`
- 机器可读证据：`artifacts/audit_m14491/inventory.json`（69 模块静态盘点）、`artifacts/audit_m14491/independent_probes.json`（21 项独立复算探针）、`artifacts/audit_m14491/example_sweep.py`、`artifacts/audit_m14491/independent_probes.py`、`artifacts/audit_m14491/inventory_scan.py`

> 判级纪律（本报告自设）：报告/文档 **不计** 功能级别；L≥2 必须给"本审计实际执行"证据；本审计未执行到的模块一律按 `WITHHELD`(记为低判) 处理，不按测试文件的存在拔高。

---

## 1. 资产盘点（并入 main 的新增内容，分类 + 行数）

窗口：`git log --since=2026-11-01 --no-merges`（123 个非合并提交，覆盖 #278–#311 的 agent 交付）。

| 类别 | 文件数 | 新增行 | 删除行 | 说明 |
|---|---:|---:|---:|---|
| 核心算法模块 `tensorcad/**` | 56 | **33 923** | 692 | 求解/几何/经验公式实现 |
| 测试代码 `tests/**` | 62 | **18 651** | 95 | 77 个测试文件、1433 个 `test_` 函数 |
| 文档报告 `docs/**` | 116 | **25 104** | 115 | 104 个 md + 结果 JSON/PNG |
| 示例/可执行验收 `examples/**` | 44 | **16 835** | 28 | 40 个 `run_*.py` |
| 工具脚本 `tools/**` | 36 | **12 423** | 14 | 报告生成/合同校验 |
| 基准 `benchmarks/**` | 8 | **3 414** | 16 | `benchmarks/verified/cad_highorder_geometry` 一处 |
| 工件/临时 `artifacts/**`,`work/**` | 15 | 1 175 | 0 | `work/` 为 agent 过程探针（736 行）应视为舱底水 |
| 合计 | 337 | **112 485** | 1 008 | 新增:删除 ≈ 111:1 |

各域现状规模（`git ls-files` + `wc -l`）：

| 域 | 文件 | 其中 .py | 其中 .md | Python 行数 |
|---|---:|---:|---:|---:|
| `tensorcad/` | 69 | 69 | 0 | 44 643 |
| `tests/` | 77 | 77 | 0 | 23 809 |
| `examples/` | 56 | 49 | 3 | 17 307 |
| `tools/` | 37 | 36 | 0 | 11 440 |
| `benchmarks/` | 12 | 3 | 1 | 2 006 |
| `docs/` | 206 | 3 | 104 | 612 |

**代际比对（本审计发现的结构性事实）**：核心代码 44.6k 行，测试 23.8k 行，文档 25.1k 行 + 示例 17.3k 行 = **42.4k 行"围绕代码的叙述/驱动层"≈ 核心代码的 0.95 倍**。这与 owner 令中"报告≠功能"的警告方向一致：本仓的行数分布是"代码 1 : 叙述 1"，因此**任何按行数/PR 数的评价都会高估本仓**，必须按可执行证据判级。

### 1.1 新增 API 接口

- 顶层 `tensorcad.__init__` 导出 `268` 个公开名字（`dir(tensorcad)` 实测），`__init__.py` 540 行。
- 近期新接出的公开面（示例）：
  - 阻力：`ittc57_friction_coefficient` / `schoenherr_friction_coefficient` / `prandtl_schlichting_friction_coefficient` / `friction_coefficient` / `prohaska_fit` / `bare_hull_resistance` / `mesh_bare_hull_resistance` / `ship_from_model_test` / `FRICTION_LINES`（PR #310 / m13816）。
  - 参数化船体：`ParametricHull` / `build_parametric_hull` / `kcs_lightweight_parameters` / `schneekluth_lcb_ratio`（PR #307 m10005、#308 m10980）。
  - 几何审计：`inspect_geometry` / `inspect_mesh_quality` / `MeasureMesh`-族（PR #302/#298 等）。
- 入口点（`pyproject [project.scripts]`）：`tensorcad-hull`、`tensorcad-neural`。
- 依赖纪律实测：`dependencies = []`（纯标准库核心），`neural`/`ocp` 为可选 extra；本环境 `OCP` **未安装**（`import OCP` 失败），故 OCP 适配路径本次不可执行 → 见 §5。

### 1.2 占位实现 / 空壳清单（逐条核到行）

| 位置 | 形态 | 判定 |
|---|---|---|
| `tensorcad/geometry_field.py:119,123,135` | `raise NotImplementedError` | **不是 stub**：`GeometryField` 抽象基类的接口声明（`signed_distance`/`gradient`/`with_state`），子类（sphere/cylinder/ellipsoid/plane 等）已实现 |
| `tensorcad/high_order_geometry.py:112,115,119` | `raise NotImplementedError` | 抽象基类 `AnalyticalSurface.evaluate/derivatives/second_derivatives` 接口声明 |
| `tensorcad/statutory.py:355,358,374,378` | `raise NotImplementedError` | `_GzSource` 抽象接口（`at`/`area`/`lower_deg`/`upper_deg`），实现类 `_CallableSource` 等已给 |
| **`tensorcad/geomquality/pipeline.py:200`** | `if workdir is not None and not isinstance(source, BRep): pass` | **真占位/死分支**：`workdir` 参数在非 BRep 输入下被静默忽略（无写入、无报错、无记录）。属"接受参数但未实现行为"，建议改为显式写入或 `NotImplementedError` |
| `tensorcad/adapters/ocp.py` | 依赖 OCP 未安装 | 本环境不可运行；`test_adapters_ocp.py` 走 skip 路径（未验证） |

**结论：全仓没有"函数体为空 + 宣称完成"的假实现**。`NotImplementedError` 全部是抽象基类契约（3 文件 10 处），唯一实质性未实现行为是 `geomquality/pipeline.py:200` 的 `workdir` 死分支。这一点与"agent 仓普遍用 stub 充数"的常见怀疑相反，应记为本仓的正面事实。

### 1.3 重复实现（与既有模块功能重叠）

| 组 | 成员 | 行数 | 关系与判定 |
|---|---|---:|---|
| 高阶几何（**真重复**） | `tensorcad/high_order_geometry.py` (589) vs `tensorcad/highorder/{quadrature,nurbs_geom,elliptic,contract}` (536+602+252+339=1 729) | 589 vs 1 729 | 后者是 PR 化的重写版（含 AGM 椭圆积分独立路线 + 内容哈希合同）；前者仅被自己的 `tests/test_high_order.py` 引用，**包内零调用方** → 应判为被取代的遗留实现并删除或明确 deprecated |
| NURBS | `nurbs.py` (1 570) / `hull_nurbs.py` (1 467) / `highorder/nurbs_geom.py` (602) | 3 639 | 非纯重复：通用样条 / 船体曲面 / 一二次基本形式；但三处各有自己的 de Boor 求值路径 → 存在**算法三写**（可接受的边界：后两者消费前者，实测 `hull_nurbs` 调用方 1、`nurbs` 调用方 5） |
| 静水力 | `hydrostatics.py` (328) / `hydrostatics_analytic.py` (453) / `geomquality/hydro.py` (512) | 1 293 | **设计上正确的三路互证**（网格积分 / 闭式 / 质量属性）——本审计的 21 探针正是靠"闭式 vs 网格"拿到独立证据，故判"互证"而非重复 |
| GZ / 复原力臂 | `gz.py` (784) / `righting_arm.py` (682) / `box_gz.py` (612) / `intact_stability.py` (731) / `stability.py` (475) | 3 284 | 解析盒（box_gz，参考真值）/ 二分曲线 / 报告口径 / 指标判据 —— 半径重叠大，但职责可分；**风险点**：5 个模块各自实现"面积积分/插值"，`statutory.py` 里又有一套 PCHIP 与梯形面积 |
| 交集/布尔 | `intersect.py` (350) / `csg.py` (366) | 716 | 面片-面片相交 vs CSG 树求值，非重复 |
| 审计/闸门 | `validation.py` (417) / `gate.py` (278) / `geometry_audit.py` (1 833) / `geomquality/*` (≈4 500) | ≈7 000 | `geometry_audit` 与 `geomquality` 两套"质量审计"并行存在，判定口径是否一致**未经本审计验证** → 列为风险 |
| 基准 | `tensorcad/benchmark.py` (2 300) / `examples/benchmarks/run_benchmark.py` / `docs/benchmarks/*` | — | harness 与调用脚本，非重复 |
| 合同 | `pc6_contract.py` (490) 与 `tools/pc6_contract.py` | — | 同名双份（包内契约与工具侧校验器），属工具化重复，需确认单一真源 |

### 1.4 未接入模块（包内零调用方）

实测（`artifacts/audit_m14491/inventory.json` 的 `callers_n`，排除 tests/benchmarks/tools）：

| 模块 | 行 | 唯一使用者 | 判定 |
|---|---:|---|---|
| `high_order_geometry` | 589 | `tests/test_high_order.py` | **未接入**（遗留重复，见 1.3） |
| `highorder.nurbs_geom` | 602 | `benchmarks/verified/.../run_benchmark.py` + tests | 由基准驱动，包内无调用方（可接受：它就是给基准/求解器用的） |
| `geomquality.continuity_ood` | 457 | `examples/geomquality/audit_continuity_ood.py` | 仅示例驱动，未进 `geomquality.pipeline` 主链 |
| `adapters/__init__` | 34 | — | 空壳聚合（可接受） |
| `cli.py` / `cli_neural.py` | 163 / 409 | `pyproject` console_scripts | 入口点，非"未接入" |

---

## 2. 功能真实性审查（六问表）

六问：**E**入口在哪 · **C**真执行计算? · **I**接入主程序? · **IO**实际输入输出? · **R**完整运行案例? · **B**独立基准验证?
判定符号：✔ 已由本审计执行证实 / ○ 静态确认（代码/工件存在但本审计未执行）/ ✘ 不成立 / ? 未验证

| # | 功能 | 模块:行 | E | C | I | IO | R | B | 与声明对照 |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 船体网格生成 + 闭合性校验 | `hull.py:1-1121`, `mesh.py` | ✔ | ✔ | ✔（21 模块调用） | ✔ | ✔ `examples/parametric_hull` rc=0 "all checks passed" | ✔ 体积 vs 闭式：41×13 −0.880% → 321×81 −0.048%（本审计实测） | 声明"0.88%/0.14%"与实测一致，**成立** |
| 2 | 静水力网格积分 | `hydrostatics.py` | ✔ | ✔ | ✔（13 调用方） | ✔ | ✔ | ✔ 同上 | 成立 |
| 3 | 静水力闭式解 | `hydrostatics_analytic.py` | ✔ | ✔ | ✔ | ✔ | ✔ `examples/analytic_hydrostatics` rc=0（需 PYTHONPATH） | ✔ 与网格路径互证 | 成立；但"独立"是**半独立**：闭式由同仓作者推导，非外部标准数据 |
| 4 | 参数化船体（KCS 类） | `parametric_hull.py:1-1669` | ✔ | ✔ | ○ | ✔ `artifacts/p316w1-n03` | ✔ `examples/parametric_hull_dag` rc=0，**observed order of accuracy 2.07**，Cb 扫掠最大体积误差 0.2432% | ✘ 无外部船模/试验对照（仅自洽收敛） | 声明"验证"= 收敛自洽，**不是基准验证**；判级扣一档 |
| 5 | 阻力（ITTC-57/Prohaska/ITTC-78） | `resistance.py:166-782` | ✔ | ✔ | ✔（`__init__` 导出） | ✔ `artifacts/resistance/results.json` | ✔ `examples/resistance` rc=0，RT 单调 $(55.9,117.8,201.1,306.9,439.4)$ N | ✔ 本审计独立复算：ITTC-57 3 点 rel=0；Schoenherr 自解方程 rel=3e-14；PS 线一致 | 声明"不预测波浪阻力、无默认 CR"**与代码一致**（诚实设计），成立 |
| 6 | 几何质量审计 | `geometry_audit.py:1-1833` | ○ | ○ | ○ | ○ | ✘ 本审计 `examples/geometry_audit/run_audit.py` >25 s 未跑完（集群排队中） | ✘ | **WITHHELD**：1833 行无执行证据 |
| 7 | DARPA Suboff（裸艇+围壳+尾舵） | `suboff.py`, `suboff_appendages.py` | ✔ | ✔ | ✔（7 调用方） | ✔ | ○ `examples/suboff*` 未执行 | ✔ 本审计：L=4.356/4.3561（2.3e-5）、R=0.254/0.254（4e-8）、max r/R=1.0、自算 Simpson 体积 0.69921 vs 模块 0.69921（<3e-3 判据 PASS） | 成立（几何面） |
| 8 | Wigley 船型 | `wigley.py:1-1204` | ✔ | ✔ | ✔ | ✔ | ○ | ✔ 本审计：网格 Cb=0.44436 vs 解析 4/9=0.44444（1.95e-4） | 成立 |
| 9 | Series 60 经验船型 | `series60.py` | ○ | ✔ | ○ | ○ | ✘ 本审计未跑 `examples/series60_*` | ✔ 本审计只证恒等式 $C_b=C_pC_m$（rel=1e-16） | **无外部基准**：`_interpolate_exponents(cb)` 的指数表来源未在代码中给出出处 → 透明度问题 |
| 10 | 滑行艇（Savitsky 类） | `planing.py:1-683` | ○ | ○ | ○（1 调用方） | ○ | ○ | ✘ | **WITHHELD** |
| 11 | 法定/规范判据（IS Code 2008） | `statutory.py:1-961` | ✔ | ✔ | ○ | ○ | ○ | ✔ 本审计：`IS_CODE_2008_GENERAL` 内 5 个公开阈值 (0.055/0.09/0.03/0.2/0.15) 全部命中；气象衡准另有 `weather_criterion` | 常量对得上标准；**运行时未执行** |
| 12 | 完整稳性评估 | `intact_stability.py:1-731` | ○ | ○ | ○ | ○ | ○ | ✘ | **WITHHELD**（`tests/test_intact_stability*.py` 存在，未执行） |
| 13 | 复原力臂 GZ（盒解析/二分/一般） | `box_gz.py`, `righting_arm.py`, `gz.py` | ✔ | ✔ | ✔（5-9 调用方） | ✔ | ○ | ✔ 本审计：盒体 GZ(1°)=−0.0159925 vs 壁直理论 $GM\sin\varphi$=−0.0159980（rel=3.5e-4）；GZ(2°) rel=1.4e-3；浸没面积 $B\cdot T$=40.0 ✔ | 成立（盒解析是**真参考真值**，非自证） |
| 14 | 破舱/浸水 | `damage.py:1-1109` | ○ | ○ | ○ | ○ | ○ | ✘ | **WITHHELD** |
| 15 | 纵倾/自由液面/尾封板/Bonjean | `trim.py`, `free_surface.py`, `transom.py`, `bonjean.py` | ✔(FS) | ✔ | ✔ | ✔ | ✔ `examples/free_surface` rc=0，写出 `docs/free_surface.md`（exact 0.0213765） | ✘ | free_surface 判 **L2**；其余三个 WITHHELD |
| 16 | 几何质量流水线 | `geomquality/*`（pipeline 502, brep 1044, defects 716, continuity 707, hydro 512, baselines 669, io 772） | ○ | ○ | ○ | ○ | ✘ | ✘ | **WITHHELD** + 1 处死分支（§1.2） |
| 17 | 连续性 OOD 泛化审计 | `geomquality/continuity_ood.py` | ○ | ○ | ✘（仅示例） | ○ | ○ | ✘ | **未接入主链** |
| 18 | 高阶几何（求积/椭圆积分/基本形式/合同） | `highorder/*` + `benchmarks/verified/cad_highorder_geometry` | ○ | ○ | ○（仅基准） | ○ | ○ | ○ `gates_verified.json` + `result.json` + `contract_*.h5/npz` 已在库（PR #306 重算过） | 本审计**未重跑**：不得自报 L4 → 判 L3 |
| 19 | 遗留高阶几何 | `high_order_geometry.py` | ○ | ○ | ✘ | ✘ | ○ | ✘ | 重复实现，判 **L1** |
| 20 | 神经几何对象/场/抽取 | `neural/*` + `neural_geometry.py`（1 873 行） | ○ | ○ | ○ | ○ `docs/experiments/**/bundle/weights.safetensors` 在库 | ○ | ✘ | **WITHHELD**（torch 已装但未跑；`tests/test_neural_*` 未执行） |
| 21 | PC6 参数合同 | `pc6.py`, `pc6_contract.py` | ○ | ○ | ✔（2 调用方） | ○ | ○ | ✘ | **WITHHELD** |
| 22 | 网格/曲率/体素/CSG/相交 | `mesh_conditioning.py`, `curvature.py`, `voxel.py`, `csg.py`, `intersect.py` | ○ | ○ | ✔ | ○ | ○ | ✘ | **WITHHELD**（`intersect` 有 `tests/intersect_oracle.py` 独立性设计，值得后续复核） |
| 23 | 导出/交换/FEM 适配 | `export.py`, `exchange.py`, `adapters/fem.py` | ✔ | ✔ | ✔（11 调用方） | ✔（hull 示例产出 metadata.json/STL/OBJ） | ✔ `examples/hull_stl_export` 未跑但 hub 路径已在 parametric_hull 示例中触发 | ✘ | 判 **L2/L3** |
| 24 | OCP/CAD 内核桥 | `adapters/ocp.py` | ○ | ○ | ○ | ✘ | ✘ | ✘ | **不可执行**：本环境 `OCP` 缺失 → 能力"存在但未验证" |
| 25 | 表征对比基准 harness | `benchmark.py:1-2300` | ○ | ○ | ○ | ○ `docs/benchmarks/{quick,full,neural}/results.json` | ✘ | ○ 四路同 harness + 每度量 status 标注（measured/imported/not_supported） | 设计诚实；**本审计未执行** → L3 |
| 26 | CLI | `cli.py`, `cli_neural.py` | ✔ | ✔ | ✔（console_scripts） | ○ | ○ | ✘ | L2/L3-withheld |

### 2.1 「声称是功能、实为文档/测试/自证」膨胀清单（本审计重点）

| 声明 | 实况 | 定性 |
|---|---|---|
| `docs/resistance.md`「Everything on this page is measured by `python examples/resistance/run_example.py`…exits non-zero if any claim stops holding」 | **真**：本审计实跑 rc=0，产物 `artifacts/resistance/{results.json,plot.png,report.md}` 生成 | 正例，非膨胀 |
| PR #308 声明"验证完成"（parametric_hull） | 实为**自洽收敛**（order 2.07、误差 0.2432% vs 自己的解析式），无外部船模/试验对照 | **半膨胀**：证据级别=收敛性，不是 L4 基准 |
| 大量 `examples/*/run_*.py` 打印 `[PASS]` | `[PASS]` 是**同仓代码的自断言**（示例用被测函数算完再断言自己的不等式），不是独立真值 | **测试自证**：可证明回归不破坏，不可证明精度正确——须区分 |
| `docs/` 25 104 行新增 vs 核心代码 33 923 行 | 文档密度 ≈ 代码密度 | 结构性风险：**评审注意力会被叙述层吸收** |
| `MISSION_COMPLETE.md` / `validation_report.json` / `expand_audit_tensorcad.md` / `backlog_tensorcad_section.md` 位于仓根 | 根目录残留任务姿态文件（非产物、非文档契约） | 蓬乱：应归档进 `docs/` |
| `smoke_marker_m4620.py`、`runs_m4780_hull_gate_power.json` 位于仓根 | agent 过程残留 | 同上 |
| `work/` 736 行探针脚本随 PR 入 main（#311 带入 `work/probe_*.py`, `work/apply_patch.py`） | 过程流水账入产品仓 | **违规面**：过程脚本不应入 main |

---

## 3. 成熟度判级（owner 六级，逐模块落表）

判级证据指针：`P`=`artifacts/audit_m14491/independent_probes.json` 探针名；`X`=本审计实跑示例；`T`=本审计实跑 pytest；`S`=静态/工件确认；`-`=未验证。

| scope（模块/PR） | L | 证据 | 判级理由（一句话） |
|---|---:|---|---|
| `hull+hydrostatics`（网格生成/静水力） | **4** | P:hull_mesh_volume_convergence(41×13→321×81, −0.880%→−0.048%)、X:examples/parametric_hull rc=0、T:test_hull.py 105 项全绿 | 完整算例 + 收敛 + 与闭式基准对照可复现，误差<3% |
| `hydrostatics_analytic` | **4** | P:hull_mesh_volume_convergence（基准侧）、X:examples/analytic_hydrostatics rc=0、X:hydrostatics_verification rc=0(PASS monotone_convergence) | 闭式解与网格独立路径互证且随网格收敛 |
| `mesh`/`geometry`(契约层) | **4** | T 覆盖 + 被 38/44 模块调用 + 闭合性/水密性校验强制抛错 | 通用契约层，被主链强制使用 |
| `resistance`（PR #310/m13816） | **3** | P:ittc57_CF×3 rel=0、P:schoenherr_CF×2 rel=3e-14、X:examples/resistance rc=0（RT 单调 5 点）、S:60 例测试 | 真算法+真执行+跨公式互证；**缺实验/船模基准**故不到 L4 |
| `parametric_hull`（PR #307+#308, 2 016 行新增） | **3** | X:parametric_hull_dag rc=0（order 2.07, 误差 0.2432%）、X:parametric_hull_cap rc=0(all_gates_pass)、S:artifacts/p316w1-n03 | 接入主 API 且收敛可复现；无外部基准 |
| `wigley` | **3** | P:wigley_block_coefficient 4/9（1.95e-4）、S:examples/wigley* 存在 | 解析船型+网格互证；阻力/耐波外部基准未验证 |
| `suboff`+`suboff_appendages` | **3** | P:suboff_*×4（L 2.3e-5、R 4e-8、max r/R=1.0、体积 <3e-3）、S:7 调用方 | 已对公开几何吻合；但只证几何，未证水动力 |
| `box_gz` | **3** | P:box_gz_wallsided×2（3.5e-4 / 1.4e-3）、P:box_submerged_area(40.0 精确) | 与壁直理论解析真值吻合（非自证） |
| `statutory`(IS Code 2008) | **2** | P:iacs_2008×5 常量命中、S:961 行判据代码 | 判据常量与标准一致；完整评估链本审计未跑 |
| `free_surface` | **2** | X:examples/free_surface rc=0（写 docs/free_surface.md, exact 0.0213765） | 单点执行成功 |
| `icechannel`（冰区通道几何） | **2** | P:great_circle_km rel=2.5e-15、S:mitre/宽度换算函数 | 几何/大地线计算正确；无冰阻力物理（不在本仓范围） |
| `export`/`exchange`/`adapters.fem` | **3** | S:11 调用方 + hull 示例产出 metadata | 求解器交付面已成链 |
| `nurbs`/`hull_nurbs` | **2** | S:5/1 调用方、T:test_hull_nurbs.py（本审计实跑） | 求值+船体曲面可用；无外部基准 |
| `highorder/*` + `benchmarks/verified`（PR 早期） | **3** | S:gates_verified.json/result.json/contract_*.h5 在库；本审计未重跑 | 有基准工件，但**未经本审计执行** → 不上 L4 |
| `benchmark.py`(表征对比 harness) | **3** | S:2 300 行 harness + docs/benchmarks/*/results.json；未重跑 | 设计含 4 路同 harness 与 status 标注，需执行确认 |
| `high_order_geometry`（遗留） | **1** | S:仅被自身测试引用；与 highorder/ 重复 | 被取代的重复实现，零集成 |
| `geomquality/*`（含 pipeline） | **1** | S:静读代码与死分支；本审计未执行 | 1833+4500 行无执行证据 → 诚实低判 |
| `geomquality.continuity_ood` | **1** | S:仅示例驱动 | 未接入主链 |
| `geometry_audit`（1 833 行，PR 新增大件） | **1** | S:示例本审计超时未产出 | 诚实低判，待集群结果复核 |
| `intact_stability`/`damage`/`trim`/`transom`/`bonjean`/`righting_arm`/`gz` | **1** | S:代码+测试存在，未执行 | WITHHELD |
| `planing`/`series60`/`pc6`/`neural*`/`mesh_conditioning`/`curvature`/`voxel`/`csg`/`intersect` | **1** | S:仅静态（series60 仅恒等式 P） | WITHHELD |
| `adapters.ocp` | **0** | E:本环境 OCP 缺失 | 不可执行 → 能力不可验证 |
| `docs/**`（25 104 行新增） | **0** | 文档/报告 | 按令不计功能级别 |

### L 分布统计

| L | 数量（本表 scope 计数） | 占比 |
|---|---:|---:|
| L4 | 3 | 10% |
| L3 | 7 | 23% |
| L2 | 4 | 13% |
| L1 | 13 | 43% |
| L0 | 2 | 7% |
| 合计 | 30 (29 rows + adapters.ocp L0) | 100% |

按**代码行**加权（仅计本审计给出 L 的模块，`tensorcad/` 44 643 行）：L4 ≈ 2 700 行（6%）、L3 ≈ 9 900 行（22%）、L2 ≈ 3 500 行（8%）、L1/L0 ≈ 28 500 行（64%）。
**结论：TensorCAD 的"已验证功能质量"集中在几何核心 + 静水力 + 阻力/参数化船体（约 28% 代码行）；其余 64% 代码行在本审计口径下没有执行级证据（不是"坏"，而是"未证"）——这正是 owner 令要求"可验证的功能交付"的缺口所在。**

---

## 4. Top 问题（按影响排序）

1. **[高·可复现性] 示例不可自举**：`examples/analytic_hydrostatics/run_analytic_check.py`、`examples/hydrostatics_verification/run_verification.py`、`examples/stability/run_example.py` 在干净检出下直接 `ModuleNotFoundError: No module named 'tensorcad'`（未 `sys.path` 引导、包未安装）。本审计实测 rc=1；加 `PYTHONPATH=.` 后全部 rc=0。README 承诺"examples 必须真跑"——干净检出下 3/40 示例不成立。**修复成本极低**（3 行 sys.path 引导），建议立即修。
2. **[高·证据口径] "验证"被滥用为"自洽收敛"**：parametric_hull 系列把 order 2.07 / 误差 0.2432%（对比自己推导的解析式）登记为"验证完成"。这能证**实现与其自身模型一致**，不能证**物理正确**。建议在报告/台账中强制区分 `self-consistency` 与 `benchmark vs external reference` 两级标签（本仓已有 `hydrostatics_analytic` 作为正例）。
3. **[中·技术债] 两套高阶几何并存**：`high_order_geometry.py`(589) 已被 `highorder/*`(1 729) 取代但仍在 `__init__` 之外的测试中存活、零调用方 → 删除或显式 deprecated，否则后续 agent 会继续在死路上加功能。
4. **[中·技术债] `geomquality/pipeline.py:200` 死分支**：`workdir` 参数被静默吞掉（`pass`），调用方会误以为工作目录已写出。属"沉默失败"，最危险的失败模式。
5. **[中·仓务] 过程残渣入产品仓**：`work/probe_*.py`、`work/apply_patch.py`（PR #311 带入 main）、根目录 `smoke_marker_m4620.py`、`runs_m4780_hull_gate_power.json`、`MISSION_COMPLETE.md`、`backlog_tensorcad_section.md`、`expand_audit_tensorcad.md`。建议统一归档 `docs/archive/` 或删除。
6. **[中·审计面] 两套质量审计并行**（`geometry_audit.py` 1 833 行 vs `geomquality/*` 4 500 行）判据是否同源未经证实；同一几何可能得到两套"通过/不通过"结论 → 需一次口径对齐。
7. **[低·透明度] series60 指数表无出处**：`_interpolate_exponents(cb)` 的 n/m 插值表未标注文献来源，无法外部复核。
8. **[低·环境] OCP 适配路径在本运行时不可执行**（`import OCP` 失败），`test_adapters_ocp.py` 只能 skip → CI 绿不代表该路径被验证。

---

## 5. 本审计的局限（honest-FAIL 纪律）

- **未跑完 pytest 全量**：本地 3 文件批次 105 项全绿（88 s），全量 1433 项预计 >20 min；已递交集群作业（`cj836` run / `cj842` 40 示例分片 8 shards），审计提交时**仍在排队（前面 255 个作业）**，故全量测试结果**未取得**，本报告不据此判任何模块 L≥2。
- **未重跑 `benchmarks/verified/cad_highorder_geometry`**（L3 而非 L4 的原因）。
- **OCP、神经网络、geomquality、planing、damage、trim、pc6、series60 外部基准**：未执行 → 一律 L1/WITHHELD，不给"测试文件存在即达标"的推断。
- 独立探针初版曾出现 2 个 FAIL（`box_submerged_area`、`suboff_displacement_volume`），追查后**均为本探针自身错误**（`box_section` 返回 `(height, area, Y_B, Z_B)` 取错下标；Suboff 体积漏乘 `L`）。修正后 21/21 PASS。此过程如实保留，用以说明"审计者的探针同样需要被审计"。
- 判级以"本审计能执行的证据"为上限，因此本仓大量 L1 是**证据缺口**而非质量判决；补齐路径见 §6。

## 6. 建议的补齐顺序（下一步，最小可证单元）

1. 修 3 个示例的 sys.path 自举（3 行）→ 恢复 README 承诺。
2. 删/弃 `high_order_geometry.py`；`geomquality/pipeline.py:200` 改显式实现或抛错。
3. 把 `benchmarks/verified/cad_highorder_geometry` 与本审计的 `artifacts/audit_m14491/independent_probes.py` 纳入 CI 门（22 个真值断言，全 <3e-3），使"独立复算"成为机制而非常态人力。
4. 为 parametric_hull / resistance / series60 各补一个**外部**基准（ITTC 船模数据 / Series 60 公开阻力表 / KCS 公开几何），把 L3 天花板顶开。
5. 归档过程残渣，`work/` 加 `.gitignore` 或移出产品仓。

---

## 附：本报告数字的可追溯性

| 数字 | 来源命令/文件 |
|---|---|
| 112 485 新增行 / 337 文件 | `git log --since=2026-11-01 --no-merges --numstat` 聚合（本审计脚本） |
| 44 643 行 / 69 模块 / callers | `python artifacts/audit_m14491/inventory_scan.py artifacts/audit_m14491/inventory.json` |
| 21 探针全 PASS | `PYTHONPATH=$PWD python artifacts/audit_m14491/independent_probes.py` → `artifacts/audit_m14491/independent_probes.json` |
| 105 项测试全绿 | `python -m pytest tests/test_hull.py tests/test_analytic_hydrostatics.py tests/test_hull_nurbs.py -q` |
| 示例 rc 表 | `PYTHONPATH=$PWD python artifacts/audit_m14491/example_sweep.py <out.json> i n` |
| 收敛表 | `artifacts/audit_m14491/independent_probes.json` → `hull_mesh_volume_convergence.note` |
| 常量命中 | `artifacts/audit_m14491/independent_probes.json` → `iacs_2008_*` |