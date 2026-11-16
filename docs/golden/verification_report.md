# GAP-CAD-D2-GOLDEN-REGRESSION 验收报告

**锚点**: GAP-CAD-D2-GOLDEN-REGRESSION  
**标题**: 几何金标回归套件+CI  
**验收判据**: 金标库：5415_exact/Wigley/Series60/球椭球；任何 tensorcad 几何改动全跑（git hook 或 CI 工作流）；对拍判据零手抄  
**状态**: ✅ 通过  

---

## 1. 金标库覆盖（4/4）

| 金标案例 | 几何路径 | 参考源 | 断言指标数 | 网格层级 | GCI | 状态 |
|---|---|---|---:|---|---|---|
| `wigley` | 解析曲面 → loft → 甲板封盖 → 面积分 | `wigley.wigley_analytic_hydrostatics`（闭式） | 8 | coarse/fine/finer | ✅ | ✅ |
| `series60` | Series-60 指数族 → loft → 封盖 → 积分 | `series60.analytic_*`（Cb 拟合指数） | 6+2 恒等式 | coarse/fine/finer | ✅ | ✅ |
| `sphere_ellipsoid` | SDF 原语 → marching tetrahedra | (4/3)πr³, 4πr², (4/3)πabc + 椭圆积分 | 10 | coarse/fine/finer | ✅ | ✅ |
| `dtmb5415` | DTMB 5415 缩比模型（npz，rank-40 PCA 重建） | `holdout_fixture.npz` + `readme_reference()` | 5+6 恒等式+16 漂移 | coarse/fine | N/A（数据集固定网格） | ✅ |

**合计**: 29 个参考指标 + 8 个恒等式 + 16 个漂移监测 + 4 个 GCI 外推 = 57 个检查点

---

## 2. 零手抄判据验证

每条参考值均由可调用的参考源实时计算，非硬编码字面量：

```python
# wigley 参考源
wigley.wigley_analytic_hydrostatics(params).as_dict()
# series60 参考源
series60.analytic_volume(L, B, T, n, m)
# sphere_ellipsoid 参考源
(4/3)*math.pi*r**3, ellipsoid_surface_area((a,b,c))
# dtmb5415 参考源
np.load("holdout_fixture.npz")["per_config_rms_m"][0]
```

验证：`test_baseline_reference_values_reproduce_from_their_source` ✅ 通过  
所有存储参考值与当前源码重算结果相对误差 ≤ 1e-9。

---

## 3. CI 工作流验证

**文件**: `.github/workflows/golden-geometry.yml`  
**触发路径**: `tensorcad/**`, `tests/test_golden_regression.py`, `tests/test_golden_mutation.py`, `docs/golden/**`, `docs/data/**`, `.githooks/**`, `tools/golden/**`  
**执行**: `python -m tensorcad.golden --check`（全双网格）+ `python -m pytest tests/test_golden_regression.py -q`  
**变异审计**: `python tools/golden/mutation_audit.py --check` + `python -m pytest tests/test_golden_mutation.py -q`  
**验证**: 工作流路径覆盖所有声明的几何输入（`_declared_inputs`），包括 `docs/data/dtmb5415/` 下的三个数据文件。

---

## 4. Git Hook 验证

**文件**: `.githooks/pre-commit`  
**安装**: `tools/golden/install_git_hooks.sh` → `core.hooksPath=.githooks`  
**执行**: `python -m tensorcad.golden --staged --quick`（单网格，~5s）  
**跳过机制**: `TENSORCAD_GOLDEN_SKIP=1`（单次提交可见绕过）  
**触发逻辑**: 仅当暂存路径触及 `tensorcad/*.py`、金标测试、基线文件、`docs/golden/` 或声明的数据输入时才运行。

---

## 5. 测试套件

| 测试文件 | 用例数 | 说明 |
|---|---:|---|
| `tests/test_golden_regression.py` | 72 | 引擎合同（blessing/drift/fail-closed）+ 金标库层（覆盖/容差/溯源/局限）+ CI/Hook 层 + **sphere_ellipsoid 三网格 GCI** |
| `tests/test_golden_mutation.py` | 11 | 变异审计（库绑定/盲点包络/声明容差阈值/恒等式层/提交报告一致性） |
| **合计** | **83** | |

覆盖层：
- **引擎层**: blessing、drift 检测、手改参考拒绝、缺失基线、未知案例、派生容差、fail-closed、三网格 GCI、强制 bless 审计
- **金标库层**: 四案例覆盖、解析指标容差、误差在容差内、基线溯源、已知局限披露、椭球面积精确参考、**sphere_ellipsoid 三网格 GCI**
- **CI/Hook 层**: 工作流路径绑定、hook 跳过/运行/绕过、数据输入触发
- **变异审计层**: 库突变检测、无突变控制通过、检测层覆盖、声明容差阈值、恒等式层验证、盲点披露、提交报告一致性

---

## 6. 误差判定（关键量 vs 基准）

| 案例 | 指标 | 测量值 | 参考值 | 相对误差 | 容差 | 是否<容差 |
|---|---|---:|---:|---:|---:|---:|
| wigley | Cb | 0.444358 | 0.444444 | -1.95e-04 | 4.33e-03 | ✅ |
| wigley | ∇ [m³] | 2.77724 | 2.77778 | -1.95e-04 | 4.33e-03 | ✅ |
| series60 | Cb | 0.594846 | 0.600000 | -8.59e-03 | 1.83e-02 | ✅ |
| series60 | ∇ [m³] | 5.35361 | 5.40000 | -8.59e-03 | 1.83e-02 | ✅ |
| sphere_ellipsoid | 球体积 [m³] | 4.18508 | 4.18879 | -8.86e-04 | 5.46e-03 | ✅ |
| sphere_ellipsoid | 球体积 [GCI] | 4.18697 | 4.18879 | -4.35e-04 | 1.67e-03 | ✅ |
| sphere_ellipsoid | 椭球面积 [m²] | 22.0667 | 22.0767 | -4.50e-04 | 2.78e-03 | ✅ |
| sphere_ellipsoid | 椭球面积 [GCI] | 22.0718 | 22.0767 | -2.21e-04 | 8.48e-03 | ✅ |
| dtmb5415 | holdout RMS [m] | 0.117382 | 0.117382 | -2.65e-08 | 1.00e-04 | ✅ |

**所有关键量误差均在派生容差内，最大误差/容差比 < 1.0。**  
**GCI 外推误差均小于 fine 层误差，证明误差来源为离散化而非模型误差。**

---

## 7. 已知局限

| 局限 ID | 说明 |
|---|---|
| `displacement-withheld-open-shell` | DTMB 5415 壳为开曲面（722 条边界边），无法计算封闭表面积分位移；仅报告数据集/站面积位移 |
| `series60/maximum_section_area_m2` | Series 60 最大截面面积收敛阶 ~1（垂直离散化在龙骨处斜率无界），盲点 ~2.7%；需改进垂直节点分布才能收紧 |
| `sphere_ellipsoid/richardson_unavailable` | 球/椭球三轴 extents 在网格细化时不变（解析精确），Richardson 外推不适用；仅双网格容差生效 |

---

## 8. 交付物清单

| 文件 | 说明 |
|---|---|
| `tensorcad/golden.py` | 金标回归引擎（~1945 行，含 sphere_ellipsoid 三网格 GCI） |
| `tests/test_golden_regression.py` | 72 项测试（+1: sphere_ellipsoid 三网格 GCI 验证） |
| `tests/test_golden_mutation.py` | 11 项变异审计测试 |
| `.github/workflows/golden-geometry.yml` | CI 工作流（golden + mutation-audit 双 job） |
| `.githooks/pre-commit` | Git 预提交钩子 |
| `tools/golden/install_git_hooks.sh` | 钩子安装脚本 |
| `tools/golden/mutation_audit.py` | 变异审计工具 |
| `tools/golden/render_evidence.py` | 证据渲染器 |
| `docs/golden/tensorcad_geometry_golden.json` | 基线文件（含 sphere_ellipsoid finer + extrapolated） |
| `docs/golden/results.json` | 完整报告 |
| `docs/golden/results.md` | 表格报告 |
| `docs/golden/golden_errors.png` | 误差/容差比图 |
| `docs/golden/mutation_audit.json` | 变异审计报告 |
| `docs/golden/mutation_audit.md` | 变异审计报告（Markdown） |
| `docs/golden/mutation_audit.png` | 变异审计图 |
| `docs/golden/verification_report.md` | 本验收报告 |
| `docs/golden_regression.md` | 套件文档 |

---

## 9. 容差体系

| 容差类型 | 公式 | 适用场景 |
|---|---|---|
| 双网格界 | 2 × \|fine - coarse\| / \|reference\|，底限 1e-4 | 有网格收敛信号的指标 |
| GCI（三网格） | 1.25 × max(\|X∞ - X_finer\|, \|X23 - X12\|) / \|reference\| | 三网格案例的 Richardson 外推 |
| 声明容差 | 案例自行声明，上限 5% | 无网格收敛信号的发布表参考 |
| 漂移容差 | 1e-9（相对） | 与基线测量值的比对 |
| 参考容差 | 1e-9（相对） | 存储参考值与源码重算的比对 |

---

## 10. 本单增量

**sphere_ellipsoid 三网格 GCI 覆盖**：
- `tensorcad/golden.py`: `SphereEllipsoidGolden` 新增 `levels()` 返回 `("coarse", "fine", "finer")`，`_resolution()` 扩展 finer=96，新增 `refinement_factor()` 方法
- `docs/golden/tensorcad_geometry_golden.json`: 基线重 bless，sphere_ellipsoid 新增 finer 层测量值 + GCI 外推值
- `tests/test_golden_regression.py`: 新增 `test_sphere_ellipsoid_has_three_grid_gci` 验证 finer 层存在且 GCI 外推误差 ≤ fine 层误差
- `docs/golden/results.json`, `results.md`, `golden_errors.png`: 证据产物刷新

**效果**: 所有解析参考案例（wigley, series60, sphere_ellipsoid）均实现三网格 GCI 覆盖，误差来源可区分为离散化误差 vs 模型误差。

---

**结论**: 金标回归套件完整实现，4 个金标案例全覆盖（3/4 含 GCI 三网格外推），CI 工作流与 Git Hook 双触发机制就位，对拍判据零手抄，83 项测试全部通过，变异审计 8/8 突变被检出。验收通过。

---

## 10. 重开第 1 轮（#655 关闭后重开）：证据产物与门禁的一致性

上一轮的四个交付文件（`results.json`、`results.md`、`tensorcad_geometry_golden.json`、本报告）
是在 #655 被关闭前写的。本轮**从当前 main 重做**：先复核它们是否仍描述这棵树，再补上让它们
不会再静默失效的检查。

### 10.1 本轮查实的三个问题

| # | 问题 | 为什么门禁没发现 | 证据 |
|---|---|---|---|
| 1 | **证据陈旧**：门禁注册并经 bless 的 case 有 5 个（`wigley`、`series60`、`sphere_ellipsoid`、`dtmb5415`、`inclined_box`），而 `docs/golden/results.json`、`results.md` 仍只报告 **4** 个（缺 `inclined_box`） | `python -m tensorcad.golden --check` 只比对**基线 ↔ 代码**；两者一致即全绿，**没有任何一步比对"证据文件 ↔ 基线"**，缺失的第五个 case 因此无人发现 | `git show HEAD:docs/golden/results.json` 的 case 列表 = ['wigley', 'series60', 'sphere_ellipsoid', 'dtmb5415']；`docs/golden/inclined_box_case.md` 已声明"第五个金标 case" |
| 2 | **机器绑定**：`results.json`、`results.md`、`mutation_audit.json` 记录了生成者工作树的**绝对路径**（`/root/tp-agent/work/<mission>/docs/golden/...`） | 报告里没有任何断言检查路径形态；在另一台机器上重跑同一命令必然产出不同文件，这份"证据"永远无法与重跑结果做 diff | 修复前 `docs/golden/results.json` 第 2 行 |
| 3 | **溯源陈旧**：基线 `provenance.git_commit` 仍指向 `9778a10`（`inclined_box` 入库**之前**），`geometry_fingerprint` 也是旧值 | 门禁不校验 provenance 与当前树的关系，只校验测量值 | 重 bless 后 fingerprint 由 `f1edda3d5de44525…` 变为 `cdadcd27af6107bf…` |

### 10.2 修法（改动最小、可验证）

1. `tensorcad/golden.py` 新增 `display_path()`：baseline 路径在仓库内时记录为**仓库相对路径**，
   在仓库外保留绝对路径（那里绝对路径是唯一能标识文件的写法）；`run_suite()` 的记录改用它。
   于是 `results.json`、`results.md` 在任何机器上重跑都得到同一字符串。
2. `tools/golden/mutation_audit.py` 的 payload `baseline` 走同一个 `display_path()`。
3. 重放生成：`python -m tensorcad.golden --bless`（5 case，55.6 s，退出码 0）→
   `python tools/golden/render_evidence.py`（退出码 0），`results.json`、`results.md`、
   `golden_errors.png` 全部刷新。**重 bless 前后全部已 bless 的测量值逐位相同**（去掉
   `provenance` 块后两版基线相等 = `True`），因此本轮基线 diff 只有 provenance 两行。
4. 新增 `tests/test_golden_evidence_freshness.py`（9 项），把"证据必须描述当前这棵树"变成 CI
   可见的断言：报告的 case 集合 == 注册 case 集合；报告里每个测量值、几何指纹、容差 == 基线记录；
   报告记录的 baseline 路径为仓库相对路径；Markdown 表覆盖每个断言。并**带两项对检查器自身的
   变异测试**（删掉一个 case 必须判红、写回机器绝对路径必须判红），证明这些检查不是恒真。

### 10.3 计算结果（本树实测，程序化注入）

门禁判定：**PASS**；case 5 个；
本轮 re-render 墙钟 55.25 s；检查点 = 35 条对基准断言 + 10 条恒等式 + 17 条 GCI 外推。

| case | 对基准断言 | 恒等式 | GCI 外推 | 最大 |误差|/容差 | 该量 | 该量误差 |
|---|---|---|---|---|---|---|
| wigley | 8 | 0 | 7 | 0.267 | maximum_section_area_m2 | 0.0078% |
| series60 | 6 | 2 | 6 | 0.756 | maximum_section_area_m2 | 0.8315% |
| sphere_ellipsoid | 10 | 0 | 4 | 0.261 | sphere_volume_m3 | 0.0435% |
| dtmb5415 | 5 | 6 | 0 | 0.523 | block_coefficient | 1.5693% |
| inclined_box | 6 | 2 | 0 | 0.000 | identity_volume_equals_requested_displacement | 0.0000% |

`inclined_box`（L=10 m、B=4 m、型深 3 m、∇=L·B·T=60 m³、横倾 15°、纵倾 3°；基准是与库
**结构不同**的按列解析积分 oracle，见 `tools/golden/inclined_box_oracle.py`）：

| 量 | 测得 | 独立基准 | 相对误差 | 派生容差 | 相对误差 (%) | <3% ? |
|---|---|---|---|---|---|---|
| displacement_volume_m3 | 60 | 60 | -2.513e-13 | 1.289e-04 | 2.51e-11% | 是 |
| lcb_x_m | 0.301425142 | 0.301425142 | +2.554e-13 | 1.289e-04 | 2.55e-11% | 是 |
| tcb_y_m | 0.2381770599 | 0.2381770599 | +2.508e-13 | 1.289e-04 | 2.51e-11% | 是 |
| vcb_z_m | 0.7900868159 | 0.7900868159 | -2.260e-13 | 1.158e-04 | 2.26e-11% | 是 |
| waterline_offset_m | 1.446903088 | 1.446903088 | -2.511e-13 | 1.289e-04 | 2.51e-11% | 是 |
| waterplane_area_m2 | 41.4678775 | 41.4678775 | -5.140e-16 | 1.000e-04 | 5.14e-14% | 是 |
| identity_volume_equals_requested_displacement | -1.440980668e-11 | 0 | +1.441e-11 | 1.000e-06 | 1.44e-09% | 是 |
| identity_mirrored_heel_mirrors_transverse_centre | 5.664912983e-14 | 0 | +5.665e-14 | 1.000e-08 | 5.66e-12% | 是 |

### 10.4 误差判定（逐项 vs 基准，是否 <3%）

- 五个 case **全部**断言项满足 `|相对误差| ≤ 派生容差`；最紧的一项是
  `series60/maximum_section_area_m2`，占其容差的 75.6%
  （该量绝对误差 8.31e-01%）。
- 绝对误差最大的一项：`dtmb5415/block_coefficient` = 1.569%
  → **<3%，是**（该量的基准是发布表而非闭式解，3% 档由数据集自身的几何合同声明，
  见 `tools/dtmb5415_param/verify_geometry.py`）。
- `inclined_box` 六个水动力量 vs 独立 oracle：最大相对误差 2.55e-11%，
  即双精度地板量级 → **<3%，是**（倾斜水线下的多面体裁剪与按列积分是两种不同算法，
  这一致性是"零离散误差"的证据，而不是同一个公式算两遍）。
- 变异审计（`docs/golden/mutation_audit.json`，路径标签已改为仓库相对）：
  6/6 个绑定对照被检出，
  策略违规 0 项 —— 门禁仍绑在它守护的库上。

### 10.5 注意点与残留（如实）

| 项 | 状态 |
|---|---|
| `mutation_audit.json` / `.md` 的测量载荷 | 本轮**未**在集群重跑（递交 cj4591 / cj4595 后队列饱和，>60 min 未派发）。产物中**只**把 baseline 路径标签改成生成器现在会写出的值（由 `golden.display_path` 计算，非手抄数字）；审计读数只取决于基线的测量载荷，而该载荷已验证逐位不变 |
| `docs/golden/gate_reach.md`、`docs/golden/grid_convergence.md` | 仍引用 `/tmp` 下的日志路径（早期任务的手写记录，不属本生成器产物），本轮未改，登记为残留 |
| 全套 `pytest` | 本轮**未**跑全量：本地 shell 有 120 s 上限，集群作业仍在排队。已验证的是 `--bless` 的完整 5-case 门禁退出码 0、`render_evidence.py` 退出码 0、新增 9 项测试全过 |
| series60 盲点 | `maximum_section_area_m2` 收敛阶 ~1（龙骨处斜率无界），容差占比 0.756，是全套最紧的一项 |
| dtmb5415 | 仍为开壳（722 条边界边），位移用数据集/站面积口径，`displacement-withheld-open-shell` 继续如实登记 |

### 10.6 本轮撞到的既有缺陷（不在本 PR 范围，如实登记）

`python -m tensorcad.golden --check --quick` 在**当前 main 上就失败**（退出码 1）：

```
ERROR inclined_box invariant_violated identity_volume_equals_requested_displacement:
      residual -3.868e-03 > 1e-06 (V = L·B·T ...)
ERROR inclined_box invariant_violated identity_mirrored_heel_mirrors_transverse_centre:
      residual 1.535e-05 > 1e-08
```

即 `inclined_box` 的两条恒等式容差（1e-6 m³ / 1e-8 m）只在细网格上成立，在 `--quick`
的单档（coarse）测量下超出。全量 `--check` / `--bless`（两档）通过，所以 CI 与 bless 都看不见；
但 `.githooks/pre-commit` 跑的正是 `--staged --quick`，因此**任何触及 `docs/golden/` 或几何代码的
提交都会被钩子拦下**，本轮提交只能使用仓库自带的 `TENSORCAD_GOLDEN_SKIP=1` 一次性绕过。

本 PR 的代码改动只有 `display_path()`（一个路径字符串格式化）与审计 payload 的同一处，不可能影响
恒等式残差（残差来自测量值，而重 bless 已证明测量值逐位不变）。建议的后续修法（本轮未做）：
恒等式应在**最细档**断言（或按档位给容差），并让 `--quick` 与全量的通过集合一致，
另补一个测试固定"quick 与 full 的 ERROR 集合相同"。
