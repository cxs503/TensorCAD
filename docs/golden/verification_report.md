# GAP-CAD-D2-GOLDEN-REGRESSION 验收报告

**锚点**: GAP-CAD-D2-GOLDEN-REGRESSION  
**标题**: 几何金标回归套件+CI  
**验收判据**: 金标库：5415_exact/Wigley/Series60/球椭球；任何 tensorcad 几何改动全跑（git hook 或 CI 工作流）；对拍判据零手抄  
**状态**: ✅ 通过  

---

## 1. 金标库覆盖（4/4）

| 金标案例 | 几何路径 | 参考源 | 断言指标数 | 状态 |
|---|---|---|---:|---|
| `wigley` | 解析曲面 → loft → 甲板封盖 → 面积分 | `wigley_analytic_hydrostatics`（闭式） | 8 | ✅ |
| `series60` | Series-60 指数族 → loft → 封盖 → 积分 | `series60.analytic_*`（Cb 拟合指数） | 8 | ✅ |
| `sphere_ellipsoid` | SDF 原语 →  marching tetrahedra | (4/3)πr³, 4πr², (4/3)πabc + 椭圆积分 | 10 | ✅ |
| `dtmb5415` | DTMB 5415 缩比模型（npz） | `holdout_fixture.npz:per_config_rms_m[0]` | 7 | ✅ |

**合计**: 33 个断言指标 + 8 个恒等式 + 20 个漂移监测 = 61 个检查点

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
**触发路径**: `tensorcad/**`, `tests/test_golden_regression.py`, `docs/golden/**`, `.githooks/**`  
**执行**: `python -m tensorcad.golden --check`（全双网格）  
**验证**: `test_ci_workflow_is_wired_to_geometry_paths` ✅ 通过

---

## 4. Git Hook 验证

**文件**: `.githooks/pre-commit`  
**安装**: `tools/golden/install_git_hooks.sh` → `core.hooksPath=.githooks`  
**执行**: `python -m tensorcad.golden --staged --quick`（单网格，~5s）  
**跳过机制**: `TENSORCAD_GOLDEN_SKIP=1`（单次提交可见绕过）  
**验证**:
- `test_git_hook_skips_when_nothing_geometric_is_staged` ✅
- `test_git_hook_runs_the_gate_on_a_geometry_change` ✅
- `test_git_hook_honours_the_skip_switch` ✅

---

## 5. 测试套件

| 测试文件 | 用例数 | 通过 | 失败 |
|---|---:|---:|---:|
| `tests/test_golden_regression.py` | 40 | 40 | 0 |

覆盖层：
- **引擎层**（22 项）: blessing、drift 检测、手改参考拒绝、缺失基线、未知案例、派生容差、fail-closed
- **金标库层**（10 项）: 四案例覆盖、解析指标容差、误差在容差内、基线溯源、已知局限披露
- **CI/Hook 层**（8 项）: 工作流路径绑定、hook 跳过/运行/绕过、CLI 报告

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

---

## 8. 交付物清单

| 文件 | 说明 |
|---|---|
| `tensorcad/golden.py` | 金标回归引擎（1142 行） |
| `tests/test_golden_regression.py` | 40 项测试 |
| `.github/workflows/golden-geometry.yml` | CI 工作流 |
| `.githooks/pre-commit` | Git 预提交钩子 |
| `tools/golden/install_git_hooks.sh` | 钩子安装脚本 |
| `tools/golden/render_evidence.py` | 证据渲染器 |
| `docs/golden/tensorcad_geometry_golden.json` | 基线文件 |
| `docs/golden/results.json` | 完整报告 |
| `docs/golden/results.md` | 表格报告 |
| `docs/golden/golden_errors.png` | 误差/容差比图 |
| `docs/golden_regression.md` | 套件文档 |
| `docs/golden/verification_report.md` | 本验收报告 |

---

**结论**: 金标回归套件完整实现，4 个金标案例全覆盖，CI 工作流与 Git Hook 双触发机制就位，对拍判据零手抄，40 项测试全部通过。验收通过。
