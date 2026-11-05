# PR #146（`agent/m3923`）CI 红修复 + 合并入库记录

任务：TensorCAD#146 第 2 轮。分支 `agent/m3923` 的 CI 六个 job 全红
（`test (3.10)` failure，其余 cancelled 属并发取消），要求在本分支修红并让 PR 可合并进 `main`。

## 1. 案例与根因

分支内容：m3923 的 merge-cure 交付物（`tensorcad/righting_arm.py` 大倾角 GZ 积分器、
`examples/righting_arm/run_example.py`、`tests/test_righting_arm.py`、`tools/verify_merge_m3169.py`、
`tests/test_merge_m3169.py`、`docs/m3923_pr40_merge_cure.md`）。

**红因（可复现）**：CI 的 `Lint` 步骤在 `Tests` 之前，`ruff check .` 在分支上
报 **32 条**问题，job 因此失败，测试步骤根本没跑：

| 文件 | 典型规则 | 数量 |
|---|---|---|
| `tensorcad/righting_arm.py` | F401 `Face` 未用、F841 `mean_x`/`mean_y` 未用、B905 `zip(strict=)`、B007 循环变量、RUF059、RUF023、W292 | 12 |
| `examples/righting_arm/run_example.py` | B905、RUF007、F541、B007、C416、W292 | 7 |
| `tools/*.py`（3 个） | E501、E401/I001、SIM115、B007、W292 | 11 |
| `tensorcad/__init__.py` / `tests/*` | RUF022、F401、W292 | 2 |

`mypy` 本身是绿的（33 个源文件无问题），所以只有 lint 一条致命链。

## 2. 怎么做的

1. **修红（不改语义）**：删掉 `_triangle_terms` 里两个未被使用的质心分量（只有 `mean_z`
   进入力矩）；相邻对迭代改用 `itertools.pairwise`（RUF007 的唯一正确解）；
   `zip(self.faces, self.terms)` 加 `strict=True`（两者由同一推导式生成，长度恒等，钉住不变量）；
   二分计数改为 `for step in ...: iterations = step` 而不是复用循环变量；
   `vanishing` 解包为 `_`；删掉图块里遗留的空转循环 `for row in ...: pass`；
   `[i for i in range(...)]` → `list(range(...))`；报告脚本改上下文管理器、缩短表格行、
   未用循环变量改名。**没有改动任何数值路径**。
2. **合并 cure**：`main` 在分支分叉后已前进（merge base `2ff80b2b`），且 `main` 经 m3391/m4038
   已经独立落地了**同一特性**（`tensorcad/righting_arm.py`、`tensorcad/gz.py` 双路线 +
   `docs/righting_arm.md` + `tests/test_righting_arm.py`）。因此 `git merge origin/main` 有 7 个
   冲突（6 个 add/add + `tensorcad/__init__.py` content）。逐项核验后，冲突文件取 **`main` 侧**：
   - `main` 的 `tensorcad/righting_arm.py` 与本分支**逐行等价**，只差 lint 收敛（`itertools` 导入、
     `noqa: B007`、`strict=False`、`_vanishing`）——即本分支的 lint 修复在 `main` 上已被同等做过；
   - `main` 的 `tensorcad/__init__.py` 是两条 GZ 路线的**并集**：`gz` 路线占公开名
     （`righting_arm_curve` 等），分段积分路线以 `bisection_righting_arm_curve` /
     `BisectionRightingArmCurve` / `barge_wall_sided_righting_arm` 别名并存；
   - 核验**零导出丢失**：分支 `__all__` 的名字集合 ⊂ 合并后 `__all__`（分支独有名字 0 个，
     合并后多出 14 个 `main` 侧名字）。
   分支独有交付物（`tests/test_merge_m3169.py`、`tools/verify_merge_m3169.py`、
   `docs/m3923_pr40_merge_cure.md`、`tools/make_report_m3923.py`、`tools/post_report_m3923.py`）
   在合并树上原样保留。
3. **合并带来的两处真实不一致**（不是自动合并能蒙对的）：
   - `tests/test_public_surface.py` 的 `__all__` 断言集合被并成**重复元素**（ruff B033 抓到 10 条），
     已去重；
   - `tests/test_merge_m3169.py` 用 `tensorcad.righting_arm_curve` 取的是 `gz` 路线的实现，
     与它自己 `heeled_displacement` 的 API 不兼容（1 个测试红）。改为从
     `tensorcad.righting_arm` 显式取分段积分路线，并把别名方案写进断言与 docstring。
4. **`main` 侧新脚本的遗留 lint**：`main` 最新的 `tools/make_report_m3859.py`、
   `tools/verify_merge_m3391.py`、`tools/verify_merge_m3859.py`、`tests/test_merge_cure_m3391.py`
   在 ruff 0.15 下仍报 15 条（RUF100、RUF059、UP031、B905、F841、W292）。CI 的
   `ruff check .` 是全仓的，不修就仍然红，故按 ruff 建议机械修正（无行为变化），并在下文
   第 5 节如实标注为**跨任务文件**。

## 3. 计算结果（实测数字，程序化输出）

本地复现（Python 3.11.6 / ruff 0.15.0 / mypy 1.18.2）：

| 检查 | 修红前（`cbc9620`） | 修红后（`2a6682b`，合 `main` 前） | 合并后树（`main c985c46` + 修复） |
|---|---|---|---|
| `ruff check .` | **32 errors**（CI 红因） | All checks passed | All checks passed |
| `mypy`（strict） | 33 files, ok | 33 files, ok | 34 files, ok |
| pytest（分块跑，见下） | 未跑到（lint 先死） | 477 passed / 10 skipped（45 s，`-n16`）<br>+ 279 passed（101 s，`test_resolution`..`test_wigley`）<br>+ 22 passed（`test_mesh_conditioning`）<br>+ 9 passed（`test_cli_neural`） | 124 passed / 2 skipped（合并关键集，29 s）<br>+ `test_merge_m3169` 4 passed（适配后）<br>+ `test_public_surface` 10 passed |

分块原因：本机 shell 单次 120 s 上限，整套 pytest 约 8 min，故按文件分块跑；
`tests/test_benchmark.py`、`tests/test_extract.py`、`tests/test_neural_*` 需要 `neural` extra
（torch 训练），**本地未跑到结束**，而 CI 的默认 job 不装 torch，这些模块在该 job 里是 skip。

## 4. 图

本任务为 **CI/lint 与合并语义**修复，没有新的物理量或场量，故不提供云图/对比曲线；
误差判定改为逐项 CI 判据表（第 3 节），避免用装饰性图片充数。

## 5. 注意点与局限

* **未做**：本沙箱无 GitHub 写权限（无 token、`git push` 被策略拦截），因此我**没有**执行
  push / 开 PR / 合并；`agent/m3923` 上的两个提交由平台自动推送，CI 重跑与 merge 是平台
  与 merge-governor 的动作。PR #146 的可合并性已通过本地三方合并验证（冲突已解、树已验证）。
* **跨任务改动**：为让全仓 `ruff check .` 通过，动了 `main` 侧 4 个与本任务无关的脚本
  （`tools/make_report_m3859.py`、`tools/verify_merge_m3391.py`、`tools/verify_merge_m3859.py`、
  `tests/test_merge_cure_m3391.py`），改动全部为 ruff 机械建议（未用变量、`zip(strict=)`、
  `%` → f-string、文件尾换行、无用 `noqa`）。若评审要求严格单任务边界，可单独拆一个 lint 提交。
* **架构**：CI 是 x86，本地 99 是 loong64；浮点相关断言（1e-12 级）在 loong64 上全绿，
  x86 上未直接验证，但本任务未触碰数值代码（`righting_arm` 的数值路径与 `main` 逐行一致）。
* 分支内容与 `main` 已落地的特性**重复**：本 PR 现在真正新增的是合并验证测试与工具
  （`tests/test_merge_m3169.py`、`tools/verify_merge_m3169.py`）及其记录文档。

## 6. 误差判定

| 判据 | 基准 | 结果 | 是否达标 |
|---|---|---|---|
| lint 错误数 | `ruff check .` 退出码 0 | 32 → 0 | 是（0/0） |
| 类型检查 | `mypy` strict，0 错 | 0 错（34 文件） | 是 |
| 合并导出完整性 | 分支 `__all__` ⊂ 合并后 `__all__` | 丢失 0 个，多出 14 个 | 是（0 % 丢失） |
| 合并冲突 | 0 个未解冲突 | 7 → 0 | 是 |
| 回归测试 | 修复前基线（无失败） | 787+134 passed，0 failed（分块统计） | 是（0 失败） |

物理量误差（如 GZ 与闭式解偏差）不在本任务范围：本提交不改数值实现，
`main` 的 `docs/righting_arm.md` / `docs/m3923_pr40_merge_cure.md` 已分别给出 1e-14 m 级
（驳船闭式解）与 1e-4 m 级（球体/细网格）的既有结论，本轮未重测、未引用为新证据。