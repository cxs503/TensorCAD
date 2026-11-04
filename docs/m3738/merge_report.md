# PR #37 治愈：`agent/m2873` 与 `main` 的三方合并（语义并集，零内容丢失）

- 仓库 / 分支：`TensorCAD`，治理分支 `agent/m2873`；本单执行分支 `agent/m3738`
- 旧 SHA（合并前 `agent/m2873`）：`206698a4f723fc70706744cbd2cedd11f46b0d39`
- `main` 基线：`8389c7e883554a391931fd1419c7441f33f97236`（"Merge #33 (m2575) via merge-governor"）
- 合并共同祖先（merge base）：`6c568c314da11a3082b41218eb5b564e466cf3a0`
- 新 SHA（合并提交）：`a0a7a39b1c0c88b0e63536769a6a1f053cf3f9d8`
- 合并提交父节点：`8389c7e`（main） + `206698a`（agent/m2873）

## 1. 为什么变脏

两条线在同一个基线上各自长出了一个船型族，并且**各自都动了版本号**，于是三方合并后
`git` 只能把同一个 hunk 判成冲突：

| 侧 | 新增能力 | 版本 |
|---|---|---|
| `main` | `tensorcad/wigley.py`（抛物 Wigley 族，精确解析静水力）+ `hull.py` 非对称族大幅演进 | 0.20.0 |
| `agent/m2873` | `tensorcad/transom.py`（艉板截断族，不完整 taper 闭式解） | 0.19.0 |

冲突不在算法，而在**两套公开面/目录/版本的登记口**。这正是机械 overlay 最容易静默丢内容的位置。

## 2. 冲突文件与逐文件解法

`git merge --no-commit --no-ff origin/agent/m2873` 报 4 个 `UU` 冲突（其余 27 个文件自动合并成功）：

| # | 文件 | 冲突性质 | 解法（语义并集） |
|---|---|---|---|
| 1 | `tensorcad/__init__.py` | `__all__` 中同一插槽 `generate_wigley` vs `generate_transom_hull`；`__version__` 两侧各写一次 | `__all__` 按字母序**两者都保留**（`generate_transom_hull` 在前、`generate_wigley` 在后）；两侧的 `from .transom import (...)` 与 5 段 `from .wigley import (...)` 全部并存；`__version__ = "0.21.0"`（承载两族的发布号） |
| 2 | `tests/test_public_surface.py` | `EXPECTED_API` 两侧各自追加自己的名字块；版本断言；分层模块元组 | 名字块**并集**（107 ∪ 110 = 120 个名字，`set(__all__) == EXPECTED_API` 仍精确成立）；分层元组含 `transom` 与 `wigley`；版本断言 `0.21.0` |
| 3 | `pyproject.toml` | `version` 一行 | `version = "0.21.0"` |
| 4 | `CHANGELOG.md` | 两侧在同一位置插入各自 release 小节 | **两节全留**，按版本降序排：新增 `0.21.0`（并集说明）→ `0.20.0`（Wigley）→ `0.19.0`（艉板）→ `0.17.0`（非对称，main 原有）。无小节被改写或删除 |

版本号取 **0.21.0** 而非任一侧的旧号，理由：0.20.0 与 0.19.0 都被对方分支"占位"，
并集发布必须是新号；否则等于声称艉板族是 0.20.0 的一部分（错误）或把 Wigley 退回 0.19.0（倒退）。
这是本次**唯一的"加一行/替代一行"**，其余全部是纯并集。

## 3. 无内容丢失的机器证明

`tools/merge_audit_m3738.py`（`python3 tools/merge_audit_m3738.py`，退出码 0）做三层审计，
结果落在 `docs/m3738/merge_audit.json`：

1. **路径集合层**：两侧文件并集中**没有任何路径**在合并树里缺失（`missing_from_merged: []`）；
   main 独有 27 个路径、m2873 独有 17 个路径，全部在树内。
2. **行多重集层**：对 4 个冲突文件，分别取 `origin/main` 与 `origin/agent/m2873` 的
   **非空行多重集**减去合并结果的多重集，逐行分类为
   - 真实丢失（判 FAIL）——**0 行**；
   - 版本字面量被并集发布号取代（`0.20.0`/`0.19.0` → `0.21.0`，已声明）—— 4 行（每侧 2 行）；
   - 重排的连接行（其所有 token 都仍在文件内，例如 `test_public_surface.py` 的分层模块元组）—— 3 行。
   逐文件：`CHANGELOG.md` 加 14 行、`pyproject.toml` 加 1 行、`__init__.py` 加 1 行、
   `test_public_surface.py` 加 3 行，且均无残留冲突标记。
3. **公开面集合层**：`merged = 120` 名字，`(main ∪ m2873) − merged = ∅`、`merged − (main ∪ m2873) = ∅`，
   两族符号 `generate_wigley / wigley_analytic_hydrostatics / generate_transom_hull /
   analytic_transom_hydrostatics` 均可从合并树 `import`。

## 4. 测试证据（本地）

环境：CPython 3.11.6、pytest 9.2.0、numpy 可用。`test_cli_neural.py` 与
`test_neural_pipeline.py` 在**两侧都**因缺少可选依赖 `safetensors` 而收集失败（`neural` extra 未装），
与本次合并无关，两侧同样排除。逐块结果（`docs/m3738/test_results.json`，由运行摘要程序化写入）：

| 测试块 | merged 通过 | main 通过 | skip | fail |
|---|---|---|---|---|
| 基础层（adapters/analytic/export/extract/gate/hydrostatics）10 文件 | 96 | 96 | 10 / 10 | 0 |
| 非对称族 + bonjean/cli/curvature/design/parameters/resolution | 159 | 159 | 0 | 0 |
| `test_hull.py` | 44 | 44 | 0 | 0 |
| mesh / conditioning / mirror / validation / voxel | 97 | 97 | 0 | 0 |
| `test_hydrostatics_analytic.py` | 47 | 47 | 0 | 0 |
| `test_benchmark.py` | 32 | 32 | 1 / 1 | 0 |
| stability / suboff / appendages / public_surface | 78 | 78 | 0 | 0 |
| `test_geometry_field.py` | 47 | 47 | 0 | **2 / 2** |
| `test_transom.py`（合并新增） | **96** | — | 0 | 0 |
| `test_wigley.py`（合并新增） | **30** | — | 0 | 0 |
| **合计** | **726** | **600** | 11 / 11 | 2 / 2 |

**结论：合并树 726 通过 ≥ main 基线 600 通过（+126 = 96 艉板 + 30 Wigley），
失败集合与基线逐条相同、未新增任何失败。**

两处 2 个失败是 `test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol`
与 `::test_the_learned_field_gradient_is_autograd_and_points_outward`，在纯 `origin/main`
工作树（`git worktree add /tmp/baseline_m3738 origin/main`）上**同样 2 failed**，
属 torch 学习场环境问题，非合并引入。

> 计算资源说明：本机负载约 220（其他智能体并发），单次全量套件超过工具 120 s 上限，
> 故按块顺序执行并逐块记录摘要；集群节点缺 numpy，`--only-binary` 无 loong64 轮子，
> 无法作为等价环境，未采用（在此如实声明）。

## 5. 注意点 / 局限

- 版本号 0.21.0 是**策略性决定**，两个上游分支的 CHANGELOG 行文未被改写；若治理侧要求保持
  0.20.0，只需改 `pyproject.toml`、`__init__.py`、`tests/test_public_surface.py` 三处字面量
  （审计脚本仍会并行全绿，因为版本行被显式归入"superseded"类别）。
- `CHANGELOG.md` 中 `0.21.0` 小节的定量陈述（Wigley 体积 `2.777777777777778 m³` =
  `4/9 × 10 × 1 × 0.625`、`tau = 0` 逐位复原尖艉）取自合并树实测
  （`wigley_analytic_hydrostatics` 与 `test_transom.py`），未手抄上游文档。
- 未修改任何一侧的算法文件（`wigley.py` / `transom.py` / `hull.py` 均字节不变），
  未改 `.git` 配置，未自行 push；`agent/m2873` 的更新交由平台推送 +
  `merge_governor` 下一周期复验。
- 局限：全量套件在本机负载下无法一次性跑完，采用分块等价覆盖（同一 pytest 进程参数、
  同一环境），分块之间无共享状态依赖的已知风险，但严格意义上不是单次单进程执行。

## 6. 误差判定

本单是合并治愈，无物理量计算，判定量是**离散的集合/计数**而非连续数值：

| 指标 | 基准 | 实测 | 误差 | <3%? |
|---|---|---|---|---|
| 合并后丢失内容行（两测） | 0 | 0 | 0 | **是（0/0，精确为 0）** |
| 公开面缺失名字数 | 0（并集应全部可导入） | 0 | 0 | **是** |
| 通过测试数 | main 600 | merged 726 | +126（+21.0%），无一项下降 | **是（不低于基线）** |
| 新增失败数 | 0 | 0 | 0 | **是** |

## 7. 产物

- `tools/merge_audit_m3738.py` —— 三层（路径/行/公开面）丢失审计，退出码可作 CI 门
- `tools/figures_m3738.py` —— 从 JSON 读数据绘图（零手抄）
- `docs/m3738/merge_audit.json`、`docs/m3738/test_results.json`
- `docs/m3738/test_comparison.png`、`docs/m3738/line_preservation.png`
- 本文件 `docs/m3738/merge_report.md`