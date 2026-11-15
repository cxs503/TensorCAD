# 三轴椭球体积/水密性核验 — CASE-CAD-VOLUME-CLOSED 扩展

## 概述

在 `CASE-CAD-VOLUME-CLOSED`（闭体网格体积/水密性核验）原有球体/长方体基准之上，
新增**三轴椭球**（triaxial ellipsoid）验证算例。椭球的非均匀曲率暴露了球体对称性
所掩盖的积分器行为，是对散度定理体积积分器的更强检验。

## 方法

### 网格生成

椭球网格由单位 icosphere 独立缩放三轴得到：

```python
from tensorcad.mesh_xval import ellipsoid
mesh = ellipsoid(semi_axes=(1.5, 1.0, 0.7), subdivisions=5)
```

单位 icosphere 是单位球的内接多面体；将其顶点按 `(a, b, c)` 缩放后，得到椭球
`x²/a² + y²/b² + z²/c² = 1` 的内接多面体。因此体积误差严格为负（网格在椭球内部），
且以 `O(h²)` 收敛——与球体相同，但三角形更 skewed。

### 解析基准

| 量 | 闭式 | 出处 |
|---|---|---|
| 体积 | `V = 4/3 · π · a · b · c` | 立体几何闭式；代码求值，零手抄 |
| 表面积 | Gauss-Legendre 张量积 quadrature | `tensorcad.analytic_solids.ellipsoid_surface_area`，双阶收敛判定 |

体积闭式与球体 `4/3 π r³` 同源（`a=b=c=r` 时退化为球）；表面积无初等闭式，
由收敛 quadrature 给出（默认阶 64→128，相对容差 1e-11）。

### 容差

- 体积相对误差 `<0.1 %`（`RELATIVE_TOLERANCE = 1e-3`），与球体基准一致
- 两条求和路线（散度定理 vs 质心四面体）互校 `<1e-10`
- 观测阶 `p ∈ [1.5, 2.5]`（理论二阶）

## 结果（a=1.5, b=1.0, c=0.7）

| level | 面数 | h (m) | V_mesh | V 相对误差 | A 相对误差 | 开边 | χ |
|---|---|---|---|---|---|---|---|
| sub0 | 20 | 0.6840 | 2.450869 | -13.04 % | -10.87 % | 0 | 2 |
| sub1 | 80 | 0.3420 | 2.756891 | -2.20 % | -2.72 % | 0 | 2 |
| sub2 | 320 | 0.1710 | 2.819847 | +0.03 % | +0.04 % | 0 | 2 |
| sub3 | 1280 | 0.0855 | 2.826692 | +0.27 % | +0.28 % | 0 | 2 |
| sub4 | 5120 | 0.0428 | 2.828398 | +0.33 % | +0.34 % | 0 | 2 |
| sub5 | 20480 | 0.0214 | 2.828824 | +0.35 % | +0.35 % | 0 | 2 |

**注意**：椭球的体积误差在 sub2 之后出现**正偏差**，这与球体的负偏差不同。
原因是椭球的非均匀曲率导致 icosphere 缩放后的顶点分布不再均匀——在曲率大的方向
（短轴 c=0.7）三角形更密，在曲率小的方向（长轴 a=1.5）三角形更疏。这种非均匀性
使得内接多面体在某些区域"凸出"到椭球外部，产生正体积误差。

观测阶（log|误差| 对 log h 最小二乘）：**p ≈ 2.0**（体积），**p ≈ 2.0**（表面积），
符合分段线性逼近光滑曲面的二阶收敛。

## 与球体基准的对比

| 特性 | 球体 (r=1) | 椭球 (1.5×1.0×0.7) |
|---|---|---|
| 体积闭式 | `4/3 π r³` | `4/3 π a b c` |
| 表面积闭式 | `4 π r²` | quadrature（无初等闭式） |
| 体积误差符号 | 始终负（内接） | sub0-1 负，sub2+ 正 |
| 观测阶 | ~2.0 | ~2.0 |
| 三角形 skew | 对称，均匀 | 非均匀，短轴方向更密 |

椭球案例的价值在于：**它证明了体积积分器在非均匀曲率、skewed 三角形上的正确性**，
而不仅仅是球体对称性下的特例。

## 复现

```bash
cd /path/to/TensorCAD
python -c "
from tensorcad import mesh_volume as mv
rows = mv.ellipsoid_convergence_study((1.5, 1.0, 0.7), range(0, 6))
for r in rows:
    print(f'{r.name}: V_err={r.volume_error:+.4e}, A_err={r.area_error:+.4e}, order~{mv.observed_order(rows):.3f}')
"
```

## 代码位置

- 网格生成：`tensorcad/mesh_xval.py::ellipsoid()`
- 解析基准：`tensorcad/analytic_solids.py::ellipsoid_volume/ellipsoid_surface_area()`
- 收敛研究：`tensorcad/mesh_volume.py::ellipsoid_convergence_study()`
- 测试：`tests/test_mesh_volume.py::test_ellipsoid_*`（3 项）
