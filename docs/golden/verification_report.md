# GAP-CAD-D2-GOLDEN-REGRESSION 验收报告

**锚点**: GAP-CAD-D2-GOLDEN-REGRESSION  
**标题**: 几何金标回归套件+CI  
**验收判据**: 金标库：5415_exact/Wigley/Series60/球椭球；任何 tensorcad 几何改动全跑（git hook 或 CI 工作流）；对拍判据零手抄  
**状态**: ✅ 通过  

---

## 1. 金标库覆盖（4/4）

| 金标案例 | 几何路径 | 参考源 | 断言指标数 | 网格层级 | 状态 |
|---|---|---|---:|---|---|
| `wigley` | 解析曲面 → loft → 甲板封盖 → 面积分 | `wigley.wigley_analytic_hydrostatics`（闭式） | 8 | coarse/fine/finer | ✅ |
| `series60` | Series-60 指数族 → loft → 封盖 → 积分 | `series60.analytic_*`（Cb 拟合指数） | 6+2 恒等式 | coarse/fine/finer | ✅ |
| `sphere_ellipsoid` | SDF 原语 →  marching tetrahedra | (4/3)πr³, 4πr², (4/3)πabc + 椭圆积分 | 10 | coarse/fine | ✅ |
| `dtmb5415` | DTMB 5415 缩比模型（npz，rank-40 PCA 重建） | `holdout_fixture.npz` + `readme_reference()` | 5+6 恒等式+16 漂移 | coarse/fine | ✅ |

**合计**: 29 个参考指标 + 8 个恒等式 + 16 个漂移监测 = 53 个检查点

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
| `tests/test_golden_regression.py` | 71 | 引擎合同（blessing/drift/fail-closed）+ 金标库层（覆盖/容差/溯源/局限）+ CI/Hook 层 |
| `tests/test_golden_mutation.py` | 11 | 变异审计（库绑定/盲点包络/声明容差阈值/恒等式层/提交报告一致性） |
| **合计** | **82** | |

覆盖层：
- **引擎层**: blessing、drift 检测、手改参考拒绝、缺失基线、未知案例、派生容差、fail-closed、三网格 GCI、强制 bless 审计
- **金标库层**: 四案例覆盖、解析指标容差、误差在容差内、基线溯源、已知局限披露、椭球面积精确参考
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
| sphere_ellipsoid | 椭球面积 [m²] | 22.0667 | 22.0767 | -4.50e-04 | 2.78e-03 | ✅ |
| dtmb5415 | holdout RMS [m] | 0.117382 | 0.117382 | -2.65e-08 | 1.00e-04 | ✅ |

**所有关键量误差均在派生容差内，最大误差/容差比 < 1.0。**

---

## 7. 已知局限

| 局限 ID | 说明 |
|---|---|
| `displacement-withheld-open-shell` | DTMB 5415 壳为开曲面（722 条边界边），无法计算封闭表面积分位移；仅报告数据集/站面积位移 |
| `series60/maximum_section_area_m2` | Series 60 最大截面面积收敛阶 ~1（垂直离散化在龙骨处斜率无界），盲点 ~2.7%；需改进垂直节点分布才能收紧 |

---

## 8. 交付物清单

| 文件 | 说明 |
|---|---|
| `tensorcad/golden.py` | 金标回归引擎（1931 行） |
| `tests/test_golden_regression.py` | 71 项测试 |
| `tests/test_golden_mutation.py` | 11 项变异审计测试 |
| `.github/workflows/golden-geometry.yml` | CI 工作流（golden + mutation-audit 双 job） |
| `.githooks/pre-commit` | Git 预提交钩子 |
| `tools/golden/install_git_hooks.sh` | 钩子安装脚本 |
| `tools/golden/mutation_audit.py` | 变异审计工具 |
| `tools/golden/mutation_redteam.py` | 红队变异工具 |
| `tools/golden/forced_bless_audit.py` | 强制 bless 审计工具 |
| `tools/golden/convergence_study.py` | 网格收敛研究工具 |
| `tools/golden/ellipsoid_area_evidence.py` | 椭球面积证据工具 |
| `tools/golden/criterion_provenance.py` | 判据溯源工具 |
| `tools/golden/render_evidence.py` | 证据渲染器 |
| `docs/golden/tensorcad_geometry_golden.json` | 基线文件 |
| `docs/golden/results.json` | 完整报告 |
| `docs/golden/results.md` | 表格报告 |
| `docs/golden/mutation_audit.json` | 变异审计报告 |
| `docs/golden/mutation_audit.md` | 变异审计报告（Markdown） |
| `docs/golden/mutation_audit.png` | 变异审计图 |
| `docs/golden/grid_convergence.md` | 网格收敛分析 |
| `docs/golden/grid_convergence.png` | 网格收敛图 |
| `docs/golden/criterion_provenance.md` | 判据溯源文档 |
| `docs/golden/criterion_provenance.json` | 判据溯源数据 |
| `docs/golden/criterion_provenance.png` | 判据溯源图 |
| `docs/golden/forced_bless_audit.md` | 强制 bless 审计 |
| `docs/golden/forced_bless_audit.json` | 强制 bless 审计数据 |
| `docs/golden/forced_bless_audit.png` | 强制 bless 审计图 |
| `docs/golden/redteam.md` | 红队变异报告 |
| `docs/golden/redteam_quick.json` | 红队变异数据 |
| `docs/golden/analytic_reference.md` | 椭球面积解析参考 |
| `docs/golden/gate_reach.md` | 门禁覆盖范围文档 |
| `docs/golden/golden_errors.png` | 误差/容差比图 |
| `docs/golden/dtmb5415_particulars.png` | 5415 主尺度图 |
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

**结论**: 金标回归套件完整实现，4 个金标案例全覆盖，CI 工作流与 Git Hook 双触发机制就位，对拍判据零手抄，82 项测试全部通过，变异审计 8/8 突变被检出。验收通过。
