# 网格质量指标库 · 求解器容器入口（OpenFOAM polyMesh）

**锚点**：`GAP-CAD-D1-MESH-QUALITY-LIB` —— 网格质量指标库（aspect/skew/nonorth/jacobian/水密/流形/取向）
**代码**：`tensorcad/meshquality/polymesh.py`（读取 + 度量 + 门），`tests/test_meshquality_polymesh.py`（14 项）
**基准与数字**：`examples/mesh_quality/run_polymesh_benchmark.py` → `artifacts/mesh_quality_polymesh/results.json` + `plot.png`

---

## 1. 为什么补这一块

库在 `meshquality.metrics` / `.topology` / `.report` 里已经能算四面体/六面体/棱柱/棱锥的
Verdict/Gmsh 指标，也能算 **OpenFOAM 定义**的内部面非正交角——但那个入口吃的是本仓库自己的
`(points, hexes)` 连通表。求解器交接过来的是 `constant/polyMesh`（`points`/`faces`/`owner`/
`neighbour`/`boundary`）。此前仓库**没有任何代码能读它**，于是：

* 文档里唯一对标 `checkMesh` 的量（non-orthogonality）无法在真实算例上测量；
* 质量门只能判本仓库自己生成的网格，不能判求解器手里的网格 —— “可作 CI 门用”只对本仓几何成立。

本模块把这一半补齐：读 polyMesh → 用 **OpenFOAM 自己的定义**算 `checkMesh` 打印的三个数
→ 复用同一个 `GateReport` 出门。

## 2. 定义（这是实现之间真正会分歧的地方，故逐条写明）

| 量 | 定义 | 理想值 |
|---|---|---|
| 面心 / 面面积向量 | 绕首点的三角扇形分解（与 `Foam::face::centre`、`Foam::face::area` 同法）；四边形面心得其面积质心，**不是**四个顶点的平均 | — |
| 单元体心 / 体积 | `cEst` = 面心按面积加权平均；再自 `cEst` 向各面切金字塔，累计有符号体积与金字塔质心，体心取体积加权平均（`primitiveMesh::calcCellCentresAndVols`）。平面面封闭多面体下**体积精确、体心为真实体积质心** | — |
| `max_non_orthogonality` | `acos( Sf·d / (|Sf||d|) )`，度；`Sf` 为 owner 外法向面积向量，`d` = 邻元体心 − owner 体心 | 0°（笛卡尔） |
| `max_skewness` | `t = ((Cf−C)·Sf)/(d·Sf)`；过体心线 `C + t d` 与面平面交于 `Pi`；报告 `|Pi−Cf|/|d|` | 0 |
| `cell_volume` | 同上金字塔分解，取绝对值；`<= 0` 计为非正单元 | > 0 |

`checkMesh` 的告警带被取作默认阈值（**政策，不是测量**）：non-orth ≤ 70°、skewness ≤ 4。
两者都可用 `limits=` 覆盖，基准脚本会打印测量值与阈值的余量。

## 3. 解析基准（差异 vs 闭式解，逐项）

| 算例 | 量 | 实现值 | 解析值 | 差异 | 是否 < 1e-12 |
|---|---|---:|---:|---:|---|
| 均匀剪切 0.5（x += 0.5 y，3×3×3 单元） | non-orth | 26.565051177078022° | `atan(0.5)` = 26.56505117707799° | 3.20e-14 | 是 |
| 同上 | skewness | 5.62e-16 | 0 | 5.62e-16 | 是 |
| 上下两单元，上层顶面在 x 向平移 2 | non-orth | 45.000000000000010° | 45° | 1.4e-14 | 是 |
| 同上 | skewness | 0.35355339059327373 | `0.5/√2` = 0.35355339059327373 | 0 | 是 |
| 笛卡尔 4×3×2 | non-orth / skewness | 0 / 0 | 0 / 0 | 0 | 是 |

上面“均匀剪切”那一行是有意义的：剪切网格的**面法向仍沿坐标轴**，而体心连线倾斜
`atan(shear)`，所以非正交角必须非零、skewness 必须为零（平行四边形面心就落在体心连线上）。
一个把 skewness 写成“面心偏移”的实现会在这行失败；一个把非正交角写成体心连线与面法向
夹角的实现会在这行通过、在下面“平移 2”那行失败（那里两者分别是 45° 与 0.354）。

## 4. 与仓库内独立路径对拍

`topology.non_orthogonality(points, hexes, "hex")` 是同一物理量的**另一条代码路径**，但它用
单元体心的**顶点平均**近似。二者只在仿射网格（平行六面体，顶点平均=质心）上应当逐位一致：

| 网格 | polyMesh 路径 max non-orth | hex 路径 max non-orth | 差 | 说明 |
|---|---:|---:|---:|---|
| cartesian_box_4x3x2 | 0.0000° | 0.0000° | 0 | 仿射，精确一致 |
| sheared_box | 26.5651° | 26.5651° | 2.84e-14 | 仿射，浮点噪声 |
| hull_nx17_nz9_L1 | 48.3781° | 51.1508° | −2.77° | 非仿射：体心定义不同 |
| hull_nx33_nz17_L3 | 94.4652° | 73.2772° | +21.19° | 同上；且已越过 70°/90° |

船体网格上的两个差值**不是 bug**，而是两个定义的差：求解器（`checkMesh`）用金字塔分解体心，
`topology` 用顶点平均。谁对？在对 `checkMesh` 说话时是前者。这也说明**顶点平均体心在细船体
网格上会把 94° 的面报成 73°**，是原来那条路径的一个已量化的近似局限。

## 5. 真实船体几何（本仓库自己的产物过求解器容器）

`volume_mesh_from_hull(HullParameters(nx, nz), first_layer=0.05, layers=L)` → `polymesh_from_hexes`
→ `write_polymesh` → `read_polymesh` → `check_polymesh`（往返后测量逐位相同）：

| 网格 | 单元 | 内部面 | max non-orth | mean non-orth | max skew | 最小单元体积 | 总体积 | 非正单元 | 退化面 | 门 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| cartesian_box_4x3x2 | 24 | 46 | 0.0000° | 0.0000° | 0.0000 | 1.000e+00 | 24.000000 | 0 | 0 | PASS |
| sheared_box | 27 | 54 | 26.5651° | 17.7100° | 0.0000 | 1.000e+00 | 27.000000 | 0 | 0 | PASS |
| stacked_shear | 2 | 1 | 45.0000° | 45.0000° | 0.3536 | 1.000e+00 | 2.000000 | 0 | 0 | PASS |
| hull_nx17_nz9_L1 | 256 | 592 | 48.3781° | 9.0686° | 1.9258 | 1.564e-03 | 7.354892 | 0 | 0 | PASS |
| hull_nx33_nz17_L3 | 3072 | 8416 | 94.4652° | 20.4325° | 11.0631 | 9.647e-06 | 7.458713 | 0 | 0 | **FAIL**（non-orth、skew） |

读法（这是本单最有用的一条结论）：

1. 船体结构化体网格在**粗档（nx17/nz9，1 层）**上交给 OpenFOAM 是**可用**的：max non-orth
   48.4°、max skewness 1.93，都在 `checkMesh` 告警带内。
2. **细档（nx33/nz17，3 层、first_layer 0.05、growth 1.2）**上 max non-orth 94.5°、max skew
   11.06 —— 已经越过 `checkMesh` 的 70°/4 告警带，门判 FAIL。来源是边界层内外尺度比过大
   （最内层 0.05 m 对最外层近 1 m，同时站距/层距在各向异性方向拉长），不是拓扑错误：无零/负
   体积单元、无退化面，总体积 7.458713 m³ 与粗档 7.354892 m³ 的一致性由加密解释（+1.4%）。
3. 因此这条链给下游（TensorFEM/TensorBEM/OpenFOAM 消费方）一个**可执行的判据**：给定层数/首层
   厚度，先过这道门再交接。这也是“可作 CI 门用”在真实船型上的第一次落地。

图 `artifacts/mesh_quality_polymesh/plot.png`：左=船体 polyMesh 非正交角分布（红线 70°）；
中=skewness 分布（红线 4）；右=两条路径的 max non-orth 对照 + `atan(0.5)` 解析线。

## 6. 注意点 / 局限（不得含糊）

* **没有与运行中的 `checkMesh` 对拍**：本机 loongarch64，OpenFOAM 与 VTK 均无可用二进制
  （与 `docs/mesh_quality.md` 第 5 节同一限制），也不能联网取源码。因此本模块的对标是
  **闭式解（1e-12）+ 仓库内独立路径（仿射网格 3e-14）**，而非外部程序输出。命令已备：
  在有 OpenFOAM 的机器上 `checkMesh -case <case>` 的 `Max non-orthogonality` /
  `Max skewness` 应与本模块的 `max_non_orthogonality` / `max_skewness` 对到 1e-6 内。
* **只读 ASCII polyMesh**：`format binary;` 直接报错并给 `foamFormatConvert -ascii` 提示，
  不猜、不静默。
* **只处理 hexahedral 连通表 → polyMesh**（`polymesh_from_hexes`）；任意多面体单元可读入并度量
  （读路径不限制面数），但本仓库没有多面体生成器，故未验证多面体网格。
* **patch 只记录、不分类**：`empty`/`cyclic`/`symmetry` 等类型保留在 `Patch.kind` 里但不改变
  度量；`empty` 边界的面参与 boundary 计数。
* **无内部面的网格判 FAIL**（`n_internal_faces == 0`）：与 `topology.gate_non_orthogonality`
  的既有约定一致——不用“空通过”掩盖一个没有内部面的网格。
* **零面积面不拒收**：塌陷列（龙骨/艏柱）会留下零面积面，此时面法向不可定向、该面非正交角记
  NaN（不进统计），单元体积仍计入并接受 `cell_volume` 检查。拒绝读入会让“测量缺陷”变成不可能。

## 7. 复现

```
python examples/mesh_quality/run_polymesh_benchmark.py            # 数字 + 图 + 解析判据（不达标 exit 1）
python -m pytest tests/test_meshquality_polymesh.py -q            # 14 项
python -m tensorcad.meshquality.cli --polymesh case/constant/polyMesh --json q.json   # CI 门
```

CLI 退出码：`0` 全过、`1` 有检查失败、`2` 读不了输入（工具坏了 vs 网格坏了可区分）。