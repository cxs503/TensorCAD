# 几何金标回归套件（Golden Geometry Regression Suite）

**锚点**: GAP-CAD-D2-GOLDEN-REGRESSION  
**版本**: tensorcad.golden/1  

---

## 概述

几何金标回归套件是 TensorCAD 的几何管线质量门禁。任何对 `tensorcad/` 下几何代码的改动，
都会自动触发对四个金标船型的重新测量，并与提交的基线对拍。

### 三条铁律

1. **零手抄判据**：参考值从不硬编码，每次运行都由声明的参考源实时重算。编辑基线文件中的数字无法通过检查。
2. **容差派生，非人为选定**：接受容差 = 2 × |fine - coarse| / |reference|（双网格界），对不动的指标设 1e-4 相对底限。
3. **Fail Closed**：缺失基线、未知案例、缺失数据文件、无 geometry_sha256 的基线条目——全部报错，永不"跳过"。

---

## 金标库（4 案例）

| 案例 | 几何路径 | 参考源 | 网格层级 |
|---|---|---|---|
| `wigley` | 解析曲面 → loft → 甲板封盖 → 面积分 | `wigley.wigley_analytic_hydrostatics`（闭式） | coarse/fine/finer |
| `series60` | Series-60 指数族 → loft → 封盖 → 积分 | `series60.analytic_*`（Cb 拟合指数） | coarse/fine/finer |
| `sphere_ellipsoid` | SDF 原语 → marching tetrahedra | (4/3)πr³, 4πr², (4/3)πabc + 椭圆积分 | coarse/fine |
| `dtmb5415` | DTMB 5415 缩比模型（npz，rank-40 PCA 重建） | `holdout_fixture.npz` + `readme_reference()` | coarse/fine |

### 断言指标汇总

- **Wigley**: 8 个静水力指标（Cb, ∇, Awp, Cp, Cm, VCB, LCB, Sw）
- **Series 60**: 6 个静水力指标 + 2 个恒等式（Cb=Cp·Cm, V=Cb·L·B·T）
- **Sphere/Ellipsoid**: 10 个指标（体积、表面积、三轴 extents × 2）
- **DTMB 5415**: 5 个参考指标 + 6 个恒等式（Cb=∇/LBT, 镜像对称, 开曲面拓扑）

**合计**: 33 个断言指标 + 8 个恒等式 + 20 个漂移监测 = 61 个检查点

---

## 使用方式

```bash
# CI / Git Hook 门禁（双网格全跑）
python -m tensorcad.golden --check

# 快速检查（仅 coarse 网格，~5s，适合 pre-commit hook）
python -m tensorcad.golden --check --quick

# 记录新基线（拒绝在有 error 时执行）
python -m tensorcad.golden --bless

# 列出金标案例
python -m tensorcad.golden --list

# 输出 JSON 报告
python -m tensorcad.golden --check --json /tmp/golden_report.json
```

### Git Hook 安装

```bash
tools/golden/install_git_hooks.sh
# 设置 core.hooksPath = .githooks
# 单次绕过: TENSORCAD_GOLDEN_SKIP=1 git commit ...
```

---

## 容差体系

| 容差类型 | 公式 | 适用场景 |
|---|---|---|
| 双网格界 | 2 × |fine - coarse| / |reference|，底限 1e-4 | 有网格收敛信号的指标 |
| GCI（三网格） | 1.25 × max(|X∞ - X_finer|, |X23 - X12|) / |reference| | 三网格案例的 Richardson 外推 |
| 声明容差 | 案例自行声明，上限 5% | 无网格收敛信号的发布表参考 |
| 漂移容差 | 1e-9（相对） | 与基线测量值的比对 |
| 参考容差 | 1e-9（相对） | 存储参考值与源码重算的比对 |

---

## CI 工作流

**文件**: `.github/workflows/golden-geometry.yml`

**触发路径**:
- `tensorcad/**`
- `tests/test_golden_regression.py`
- `tests/test_golden_mutation.py`
- `tools/golden/**`
- `docs/golden/**`
- `docs/data/**`
- `.githooks/**`
- `.github/workflows/golden-geometry.yml`

**Jobs**:
1. `golden`: 运行金标检查 + 验证报告 schema + 运行测试套件
2. `mutation-audit`: 变异审计（证明门禁有牙齿）

---

## 已知局限

| 局限 ID | 说明 |
|---|---|
| `displacement-withheld-open-shell` | DTMB 5415 壳为开曲面（722 条边界边），无法计算封闭表面积分位移；仅报告数据集/站面积位移 |

---

## 相关文件

| 文件 | 说明 |
|---|---|
| `tensorcad/golden.py` | 金标回归引擎 |
| `tests/test_golden_regression.py` | 引擎合同 + 金标库测试 |
| `tests/test_golden_mutation.py` | 变异审计测试 |
| `.github/workflows/golden-geometry.yml` | CI 工作流 |
| `.githooks/pre-commit` | Git 预提交钩子 |
| `tools/golden/install_git_hooks.sh` | 钩子安装脚本 |
| `tools/golden/mutation_audit.py` | 变异审计工具 |
| `docs/golden/tensorcad_geometry_golden.json` | 基线文件 |
| `docs/golden/mutation_audit.md` | 变异审计报告 |
| `docs/golden/grid_convergence.md` | 网格收敛分析 |
| `docs/golden/criterion_provenance.md` | 判据溯源文档 |
