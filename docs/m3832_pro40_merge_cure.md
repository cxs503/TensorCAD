# PR #40 治愈报告：TensorCAD `agent/m3169` 与 main 的语义并集合并

**仓** TensorCAD · **分支** `agent/m3169` · **单** m3832（治愈单） · **基线** `2ff80b2` (main)
**合并提交** `01b439e74aa0` · **旧分支尖端** `fcf15ad`（PR #40 mergeable=false）

## ① 案例介绍

PR #40 带的是**静稳性**这一条能力：`tensorcad/righting_arm.py` 在**定排水量**下求 `GZ(φ)`——
不是拟合，而是在每个横倾角上用散度定理直接积分**倾斜后的湿表面**（`V=∮(z−h)n_z dA`、
`∫y dV=∮½y²n_y dA`、`∫z dV=∮½(z−h)²n_z dA+hV`），再用对单调 `V(h)` 的括区二分求平衡水线；
配套两个**解析真值**：墙侧驳船 `GZ=sinφ(GM+½BM tan²φ)`（甲板边浸没前精确）与浮球 `GZ=GM sinφ`。

期间 main 前进到 `0.21.0`，加进了三样**会与它直接相撞**的东西：`Hydrostatics` 记录上的
水线二次矩（`I_T`、`I_L`、`keel_z`）、**艉板船体系族** `tensorcad/transom.py`（切船位置 `τ` 与
正弦收尾指数 `p` 的设计空间）、以及**非对称船体**的首尾指数（`fore_exponent`/`aft_exponent`）。
两条线都改了 `tensorcad/__init__.py` 的导出表、`pyproject.toml` 版本、`CHANGELOG.md` 顶部与
`tests/test_public_surface.py` 的模块清单——正是 PR #40 被判 dirty 的四处。

| 量 | 值 | 来源 |
|---|---|---|
| 驳船 L×B×T / 干舷 | 10.0×2.0×0.8 m / 0.4 m | `docs/righting_arm/results.json` |
| 驳船 KG / 甲板边角 | 0.55 m / 21.80° | 同上 |
| 参数化船体水线 / KG / ∇ | -0.47 m / 0.55 m / 3.5524 m³ | `curve_payload` |
| 球半径 / KG | 1.0 m / 0.9 m | `results.json:sphere` |
| 艉板参数 | τ ∈ {0, 0.05, 0.2, 0.35}, p = 0.7 | `convergence.csv` |

## ② 怎么做的

1. **建侧树**：`git worktree` 在 `origin/agent/m3169`（`fcf15ad`）上开 `cure/m3169-m3832`，
   `git merge --no-ff origin/main`。合并基 `df31744cd565`。
2. **逐文件解冲（拒绝机械 overlay）**：只有 **2** 个文件真冲突，其余全部自动合并。
   - `CHANGELOG.md`：两侧都把条目写进 `## 0.21.0 — 2026-11-04`。**语义并集**=保留本分支的 righting-arm
     条目（29 行贡献）＋ main 的 `### Added` 块
     （81 行贡献），同一版本节内前后拼接。
   - `tests/test_public_surface.py`：`test_layered_modules_are_all_importable` 的模块元组并集
     `{stability, righting_arm, bonjean} ∪ {transom}`——首轮脚本改出了重复行（`wigley/hull/...`
     出现两次），已删除重复并 AST 复核。
3. **写机器校验器**（`tools/verify_merge_m3832.py`）：直接读三个提交的 tree 对象，逐文件/逐行核对：
   * 任一侧新增的文件在结果中必须**逐字节**存在（冲突文件除外）；
   * **单侧改动**文件在结果中必须与那一侧逐字节相同；
   * 冲突文件必须是**行集合的并集**：两侧相对合并基的"新增行多重集"都要被结果覆盖，且无残留冲突标记；
   * `tensorcad.__all__` 必须是两侧的超集。
4. **测试对跑**：同一命令在 `2ff80b2`（main 干净检出）与合并结果各跑一遍全量 pytest。
5. **示例位复现**：两个 `examples/*/run_example.py` 在合并树上重跑，验证两条功能线在新代码上仍成立。

## ③ 计算结果

**合并完整性（机械校验，失败项 0 条）**

| 检查 | 数值 | 判据 |
|---|---|---|
| 冲突文件 | 2（CHANGELOG.md、tests/test_public_surface.py） | 两侧新增行缺失数均为 0 |
| 单侧改动逐字节一致 | 43 个文件 | 全部通过 |
| 两侧完全一致的改动 | 180 个文件 | 全部通过 |
| 两侧都改（需语义判断） | 3 个文件 | 自动合并 + 逐行并集校验 |
| `__all__` 名字数 | ours 117 ∪ theirs 120 → merged 130 | 超集成立 |
| 被合并丢掉的路径 | 0 | 必须为 0 |

**测试证据**

| 树 | 命令摘要 | passed | failed | skipped | 用时 |
|---|---|---|---|---|---|
| main `2ff80b2` | 全量 pytest（排除 2 个采集失败模块） | 818 | 2 | 11 | 675s |
| 合并结果 | 同上 | 859 | 2 | 11 | 523s |
| **Δ** | | **+41** | 0 | 0 | |

失败集完全相同（`same_failure_set = True`）：两侧同样是
`tests/test_geometry_field.py::test_the_learned_revision_speaks_the_same_protocol` 与
`::test_the_learned_field_gradient_is_autograd_and_points_outward`，根因是可选依赖 `safetensors`
缺失（`tensorcad/neural/field.py:26 ModuleNotFoundError`），**不是合并引入的回归**；同一原因使
`tests/test_cli_neural.py`、`tests/test_neural_pipeline.py` 在两侧同样采集失败，故两侧用同一命令排除。
+41 例正是本分支带来的 `tests/test_righting_arm.py`（41 例）。

**示例位复现（端到端，合并树上跑）**

* `examples/righting_arm/run_example.py`：`10/10 claims hold`，且重写出的 `docs/righting_arm/` 与仓库内
  提交版本 **git 无差异**（逐字节复现）——说明 righting-arm 在 main 改过的 `hydrostatics.py`/`hull.py`
  上行为不变。
* `examples/transom_stern/run_example.py`：4 个 τ 全部 `watertight=True manifold=True euler=2`，
  生成的 6 个文件与 `docs/m2873/outputs/transom_stern/` **逐字节一致**。

## ④ 结果图

`artifacts/m3832/post_merge_verification.png` 四联图（图见附件，caption 内嵌）：
(a) 本分支 GZ(φ) 曲线在合并树上重跑；(b) 墙侧解析式残差（窗口内 vs 窗口外）；
(c) main 的艉板体积误差收敛；(d) 基线 vs 合并测试数。四张图**同图并列**即为"两条功能线在同一个树里同时活着"的证据。

## ⑤ 注意点（踩过的坑 / 局限）

* **SIGPIPE 会把产物文件截断**：把校验器输出 `| head -40` 会触发 `BrokenPipeError`，json **没写全**，
  而我又把它 commit 了，直到 diff 才看出 182 行"消失"。校验器现在把 json 写盘放在 print 之前，
  产物体积本身也是判据（`merge_verification.json` 现 1.4 KB，`failures: []`）。
* **shell 工具 120 s 硬限**：全量 pytest 要 8–11 min，必须 `setsid nohup ... &` 脱离进程组，
  否则"进程组已整组收割"；`--ignore=` 排除的两个模块是**环境**问题（可选依赖 `safetensors`），
  基线同样失败，故不算回归。
* **脚本解冲也会错**：第一版并集脚本把模块元组写成 3 行，产生重复项（语法仍合法，测试才会暴露）；
  解冲脚本必须配断言 + AST/重复项扫描，不能只看"没冲突标记"。
* **版本号未 bump**：两侧都在 `0.21.0`，本单按"并集保内容"原则保留 `0.21.0` 并把两批特性写进同一节；
  若维护者认为合并产物应记为新版本，改 `pyproject.toml` 与 `test_public_surface.py` 的版本断言即可（一处两行）。
* **未在 loong64 集群跑**：本轮全量测试在 99 控制面本机（x86-64 / py3.11 / numpy 2.4.6）完成；
  GitHub CI 为 x86，架构判读不受影响。

## ⑥ 误差判定（合并后的数值仍在原精度内）

| 关键量 | 合并后实测 | 基准 | 相对误差 | <3%？ |
|---|---|---|---|---|
| 墙侧驳船 GZ（有效性窗口内） | 最差 1.21e-14 m @ 20° | 解析闭式 | 8.2e-12% | **是** |
| 浮球 GZ（最细网格） | 最差 6.51e-05 m | `GM sinφ` | 0.065% | **是** |
| 参数化船体小角斜率 | 最差相对误差 6.9e-08 | `KB+I_T/∇−KG` 独立积分 | 6.9e-06% | **是** |
| 艉板体积（最粗网格 nz=13） | 0.880% | 解析 `sinc^p` 级数 | 0.880% | **是** |
| 艉板体积（最细网格） | 0.044% | 同上 | 0.044% | **是** |
| 不完全锥度积分（tau>0 vs tanh-sinh 求积） | 1.1e-15 | `Γ((1+p)/2)/(√π Γ(1+p/2))` / `tanh_sinh_quadrature` | 1.1e-13% | **是** |
| 合并内容完整度 | 双侧新增行缺失 0 行 | 三提交 tree 对照 | 0% | **是** |

结论：**是**，全部关键量相对基准误差远小于 3%（多数在浮点/1e-14 量级）；合并未改变任何数值路径。
honest 部分：PR #40 的 mergeable 判定由 merge_governor 下一周期复验，本报告只提供本地证据链；
两个 `test_geometry_field` 失败与两个模块采集失败属环境缺 `safetensors`，两侧一致，未修（超出本单范围）。

---
产物：`tools/verify_merge_m3832.py`、`tools/make_figures_m3832.py`、`tools/report_m3832.py`、
`artifacts/m3832/{merge_verification.json, pytest_summary.json, post_merge_verification.png}`、
本文件。

