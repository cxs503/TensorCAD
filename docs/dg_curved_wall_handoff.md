# CAD → DG curved-wall hand-off (`tensordg.curved_geometry.v1`)

锚点：`GAP-CAD-B3-SOLVER-INTERFACE`（CAD→求解器消费面，其中「DG 高阶曲边」一条此前无证据）。

## 1. 问题设定

TensorDG 的 `tensordg/curved_geometry.py` 消费一份 `tensordg.curved_geometry.v1` 合同
（NURBS 控制点 + 节点向量 + 权重），用 `mesh_from_contract` 生成等参 DG 单元几何。
其 docstring 自称这是「TensorCAD(N03) → TensorDG 的交接」，但**合同在本仓没有任何生产者**：
唯一能产出该 payload 的是 TensorDG 自己的 `circle_contract`。即求解器把几何发给自己，
CAD 误差无法进入链路，几何错误也无法被判定为几何缺陷而非求解器缺陷。

本轮补的就是这一段：**CAD 侧成为几何的权威来源，并且可机器复核**。

## 2. 方法（关键代码指针）

| 环节 | 实现 | 说明 |
|---|---|---|
| 构造 | `tensorcad/dg_handoff.py::circle_contract` | 用 CAD 的**切点相交**构造有理二次圆弧：段中控制点 = 两端切线的交点（距心 `R/cos(Δθ/2)`），权重 `cos(Δθ/2)`。由 `radius/centre/segments/plane` 参数化生成，不抄任何常数表 |
| 独立求值 | `evaluate_curve` / `evaluate_curve_derivative` | 自写有理 de Boor（Piegl & Tiller A3.2/A4.2，齐次坐标 + 商法则），不与消费端共用代码 |
| 度量 | `curve_geometry_report` | 闭合间隙、最大径向偏差、周长（按节点区间分片 20 点 Gauss 求积 `∫|C'(t)|dt`）、面积（Green 定理）；基准 `2πR`、`πR²` 为解析值 |
| 坐标系 | `PlaneMap` | 2D 图表放进 z-up 三维契约：记录 `u/v/normal` 轴、平面偏移、**图表手性**（`x-y` 为 +1，`x-z` 为 −1）|
| 交接 | `write_dg_bundle` → `curved_wall.json` + `handoff.json` | 合同是**文件**不是私有 import；三维源环随包落盘 |
| 复核 | `verify_bundle` | 独立解析器回读：重算全部申报数字、digest、平面投影；不一致即 `HandoffError` 拒收 |

端到端算例：`examples/dg_curved_wall/run_dg_handoff.py`（CAD 球赤道 → 合同 → TensorDG 环形
区域等参几何 → 刚体旋转对流 DG 求解 → 收敛阶门禁）。

## 3. 算例与参数

- CAD 源：`solver_interface.ellipsoid_mesh((1,1,1), 24, 64)` 的赤道环（与 BEM 路线同一个几何对象），
  三维 z-up，投影到 `x-y` 水位面图（`plane=xy`, `offset_m=0`）。
- 域：环形 `r_inner = R`（CAD 内壁）、`r_outer = 2.0 m`；几何阶 `g = p + 1`（超参数化）。
- 物理：`u_t + ∇·(A u) = 0`，刚体旋转 `A = (-Ωy, Ωx)`，`Ω = 1 s⁻¹`；
  解析解 `exp(-((r-1.5)/0.35)²)·cos(2(θ-Ωt))`，在环上光滑、壁面通量为零。
- 终时 `t = 0.3 s`；CFL 0.2；网格档 (8×16), (12×24), (16×32)；p = 2, 3。

## 4. 结果

### 4.1 CAD 几何合同复核（`verify_bundle`）

| 量 | 值 | 解析基准 | 相对误差 |
|---|---:|---:|---:|
| 最大径向偏差 | 2.22e-16 m | 0（正圆） | 2.22e-16 |
| 闭合间隙 | 0.0 m | 0 | — |
| 周长 | 6.283185307179591 m | 2πR = 6.283185307179586 | 7.07e-16 |
| 包围面积 | 3.1415926535897953 m² | πR² = 3.141592653589793 | 7.07e-16 |
| 平面检查 | 面外 0 m，壁面径向误差 0 | 赤道环在 z=0 | 0 |

合同被 TensorDG 自己的 `validate_geometry_contract` 接受（digest `e232f15f…f36f`），
且本仓求值与 TensorDG `nurbs_curve_points` 在 401 个参数点上的最大差 **3.33e-16 m**（跨实现互校）。

### 4.2 DG 收敛阶（相对 L2 误差 vs 解析解）

| p | g | 8×16 | 12×24 | 16×32 | 观测阶 | 理论 p+1 | 偏差 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | 3 | 2.441e-03 | 7.130e-04 | 3.048e-04 | **3.004** | 3 | 0.13 % |
| 3 | 4 | 1.242e-04 | 2.251e-05 | 8.050e-06 | **3.965** | 4 | 0.88 % |

壁面几何保真度（内壁沿单元边密采样的 `max|r−R|/R`）：p=2 → **1.24e-05**，p=3 → **2.33e-08**。

### 4.3 对照（曲边几何到底值不值）

| 情形 | p | 观测阶 | 壁面几何误差 (rel) |
|---|---:|---:|---:|
| CAD 合同壁（g = p+1） | 2 | 3.004 | 1.24e-05 |
| CAD 合同壁（g = p+1） | 3 | 3.965 | 2.33e-08 |
| 求解器侧折线壁（g = 1） | 2 | **1.251** | 1.92e-02 |
| 求解器侧折线壁（g = 1） | 3 | **1.171** | 1.92e-02 |

折线壁的 1.92e-02 恰等于 `1 − cos(π/16) = 0.0192`，即弦与圆的矢高；阶由此塌到 1 阶——
高阶 DG 在低阶几何上不涨阶，这正是 CAD 必须把曲边交给求解器的量化理由。

**负向控制（门是真的）**：把内壁做成过圆顶点的 1 次闭合 NURBS（折线壁）交合同时，
TensorDG 拒绝并给出原因
`inner wall is not circular; the annulus builder needs a circular wall (got radius spread 5.027e-02)`。
交接不是注释，是闸门。

图：`artifacts/dg_curved_wall/convergence.png`（误差-步长与 p+1 参考斜率 + 折线壁对照）、
`wall_geometry.png`（CAD 合同曲线 / DG 壁面节点 / 径向偏差）、`solution_field.png`（t=0.3 s 数值场与逐点误差）。
数值：`artifacts/dg_curved_wall/results.json`、`artifacts/dg_curved_wall/cad_wall/handoff.json`。

## 5. 注意点与局限（如实）

- **只支持圆形内壁**：TensorDG 的 annulus builder 要求 `max|r−R| ≤ 1e-8 R`，通用船体剖面/肋位曲线
  暂时无法通过该入口；本模块只提供圆（`circle_contract`），不做测量几何的 NURBS 拟合。
- **合同是曲线（壁面），不是面/体**：走通的是 DG 的曲壁路线，不是曲边体网格。
- 平面图 `x-z`（中纵剖面）手性为 −1：2D 逆时针积分对应 +y 反向的物理法向。本算例用 `x-y`（+1），
  该字段已在 `handoff.json` 声明，消费端不得默认。
- 2D 平面内的「旋转 Ω」是教学/验证用的解析流动，**不宣称**任何阻力、推进或耐波预报能力。
- 折线壁对照用的是求解器内部 `annulus_mesh(gdeg=1)`，不是 CAD 产物（CAD 折线壁会被门拒绝，
  见上），两者用途不同，不能混为一谈。

## 6. 误差判定（是否 <3%）

| 关键量 | 基准 | 误差 | <3%？ |
|---|---|---:|---|
| 合同几何：最大径向偏差 | 正圆（解析） | 2.22e-16 rel | ✅ |
| 合同几何：周长 | 2πR | 7.07e-16 rel | ✅ |
| 合同几何：面积 | πR² | 7.07e-16 rel | ✅ |
| 跨实现互校（本仓 vs TensorDG 求值） | 401 点最大差 | 3.33e-16 m | ✅ |
| p=2 收敛阶 | p+1 = 3 | 3.004（+0.13 %） | ✅ |
| p=3 收敛阶 | p+1 = 4 | 3.965（−0.88 %） | ✅ |
| 折线壁对照 | p+1 = 3 / 4 | 1.251 / 1.171（几何限阶，**如实记录为失败对照**） | ❌（预期失败，用于证明曲边必要） |

## 7. 复现

```bash
PYTHONPATH=.:/root/tensor-ecosystem/TensorDG python examples/dg_curved_wall/run_dg_handoff.py
PYTHONPATH=.:/root/tensor-ecosystem/TensorDG python -m pytest tests/test_dg_handoff.py -q   # 15 passed
```
平台报告：id 18722（scenario_key `cad_dg_curved_wall_handoff`）。
