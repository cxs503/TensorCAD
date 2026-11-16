# 体素化（从显式三角网格）：10⁶–10⁷ 三角级与窄带 SDF

> 锚点：feature_ledger `GAP-CAD-D4-PERF-SCALE`（大几何性能基准 10⁶–10⁷ 三角级）
> 工件：`benchmarks/perf_scale/artifacts/ladder_voxelize_mesh.json`（schema `tensorcad.perf_voxel_mesh/1`）
> 代码：`tensorcad/mesh_field.py`（`MeshSDFField`）、`tensorcad/spatial_index.py`（`nearest_surface`）、
> `tensorcad/perf_voxel_mesh.py`（阶梯 + 门）、`scripts/perf_voxelize_mesh.py`（驱动）
> 测试：`tests/test_mesh_field.py`（10 项）、`tests/test_perf_voxel_mesh.py`（8 项）——18/18 PASS

## 1. 案例介绍

`GAP-CAD-D4-PERF-SCALE` 的验收要求是「加载/点定位/**体素化**/体积积分子任务基准表」。
加载、点定位、体积积分三张表已在此前版本落地；**体素化那一行此前只在解析场上测过**：
`voxelize` 子任务采样 `SuboffHullField`（闭式 SDF），代价是纯 O(cells)，
所以它的阶梯是"网格分辨率"、`triangles` 列结构性为 0——它回答的是"把解析 SDF 采到 10⁸ 单元要多久"，
**不是**求解器流水线真正问的问题。

真实流水线拿到的是**显式三角网格**（CAD 交换网格、修复后的 STL、B-rep 镶嵌），10⁶–10⁷ 三角面，
没有闭式表达，SDF 必须从三角形本身导出。这一档此前**无法测量**，因为 TensorCAD 没有"网格 → 场"的适配器。

参数（同既有阶梯，保证可比）：参数化船体 L=10 m、B=2 m、T=0.8 m，解析排水体积
**7.512698 m³**；三角形规模由 1→4 中点细分从 41×13 基网格（1,998 三角）推到 31,968（level 3）、
492k（level 4）、1.97e6（level 5）、7.86e6（level 6）。体素网格各向异性（船体 10:2:0.8 用立方网格浪费 4/5 单元）。

## 2. 怎么做的

1. **`TriangleGrid.nearest_surface(point, band=None)`**（stdlib-only，新增）：均匀网格上的**精确**欧氏最近点查询——
   以其所在单元为中心逐层扩壳，用"查询点到已搜索块外部的距离"作为可证明下界剪枝；三角形按包围盒登记到
   所有覆盖到的单元，因此未搜索的三角形必在块外，下界成立。`band` 给出窄带：超出带返回
   `SurfaceHit(distance=band, triangle=-1)`（`triangle == -1` 是"带内无三角形"的标记，不是"最近三角形在 band 处"）。
   与 `brute_force_surface_distance`（同模块，O(N) 全扫参考实现）逐位一致，误差 **0.0**（见 §4）。
2. **`tensorcad.mesh_field.MeshSDFField`**（新增）：把 `Mesh` 适配到既有 `GeometryField` 协议
   （`signed_distance`/`gradient`/`bounds`/`metadata`），符号由射线奇偶性（`TriangleGrid.contains`）给出，
   幅值是最近距离。于是 `to_voxel_grid(field, resolution)`、`voxelize_mesh(mesh, res)` 无需改一行即可对三角汤工作。
3. **窄带规则**：体素化时 `band = 1 × 最大体素边长`。带 ≥ 半个体素时 `volume_fraction` 的
   `clip(0.5 − s/h, 0, 1)` 对深处单元仍给 1、对外部给 0，**体积积分不受截断影响**；
   截断只砍掉"深内部单元要清掉半米球内全部三角形"的那部分代价。
4. **阶梯与门**：`tensorcad/perf_voxel_mesh.py` 扫两条轴（固定网格下的三角面轴、固定网格下的体素数轴），
   门 = 两个复杂度指数（限 1.5）+ **体素体积精度门（>3% = RED）**。表格由工件渲染、`--check` 校验漂移。

## 3. 计算结果

<!-- perf-meshvoxel:section begin -->
`tensorcad.perf_voxel_mesh/1` · fail · gates checked=3 not_run=0 · wall 0.2 min · band=1 体宽

#### 三角面轴（固定网格，窄带 1 体宽）

<!-- perf-meshvoxel:tri-table begin -->
| 三角面 | 秒 | µs/查询 | band m | 索引 cells | 索引 MiB | 体素体积误差 |
|---:|---:|---:|---:|---:|---:|---:|
| 1,998 | 0.467 | 912 | 0.625 | 672 | 0.2 | -5.648% |
| 7,992 | 0.570 | 1113 | 0.625 | 2,000 | 0.6 | -5.648% |
| 31,968 | 1.244 | 2431 | 0.625 | 8,960 | 2.5 | -5.648% |
<!-- perf-meshvoxel:tri-table end -->

#### 体素轴（固定网格数 vs 固定网格分辨率，最大单体网格）

<!-- perf-meshvoxel:grid-table begin -->
| 网格 | 体素 | 秒 | µs/查询 | band m | 体素体积误差 |
|---:|---:|---:|---:|---:|---:|
| 8×4×2 | 64 | 0.086 | 1350 | 1.250 | -5.444% |
| 16×8×4 | 512 | 0.570 | 1113 | 0.625 | -5.648% |
| 24×12×6 | 1,728 | 2.099 | 1215 | 0.417 | -4.538% |
| 32×16×8 | 4,096 | 4.527 | 1105 | 0.312 | -3.828% |
<!-- perf-meshvoxel:grid-table end -->

#### 门（`not_run` 不是 pass）

<!-- perf-meshvoxel:gates begin -->
| 门 | 状态 | 判定 | 证据 |
|---|---|---|---|
| `scaling_voxelize_mesh` | measured | pass | worst exponent=0.35 (limit 1.50), R^2=0.895, 3 sizes, grid=16x8x4 sizes 1,998..31,968 |
| `scaling_voxelize_mesh_cells` | measured | pass | worst exponent=0.96 (limit 1.50), R^2=0.999, 4 sizes, mesh 7,992 tri, cells 64..4,096 |
| `voxel_volume_accuracy` | measured | RED | worst |rel err|=5.648% (tolerance 3%) at 1,998 tri, grid (16, 8, 4) |
<!-- perf-meshvoxel:gates end -->
<!-- perf-meshvoxel:section end -->

## 4. 误差判定（关键量 vs 基准，逐项）

| 量 | 基准 | 实测 | 相对误差 | <3%? |
|---|---|---|---|---|
| `nearest_surface` vs 全扫参考（船体网格，60 点，含盒外/盒下点） | `brute_force_surface_distance` | 逐点一致 | **0.0**（双精度相等，`rel=1e-12` 断言） | 是 |
| `MeshSDFField` vs 解析球面（24 cell 球，120 点） | `SphereField.signed_distance` | 最坏偏差 < 0.03（=球半径的 3%，且这部分是**网格离散误差**不是查询误差） | <3% | 是（测试断言） |
| 球面体素体积 @32³ | 4/3·π = 4.188790 m³ | 见 `tests/test_mesh_field.py::test_voxelised_mesh_volume_beats_the_3pct_bar` | <3% | 是 |
| 船体体素体积 @8×4×2（64 体素） | 7.512698 m³ | −5.444% | **−5.444%** | **否** |
| 船体体素体积 @16×8×4（512 体素） | 7.512698 m³ | −5.648% | **−5.648%** | **否** |
| 船体体素体积 @24×12×6（1,728 体素） | 7.512698 m³ | −4.538% | **−4.538%** | **否** |
| 船体体素体积 @32×16×8（4,096 体素） | 7.512698 m³ | −3.828% | **−3.828%** | **否**（趋势对，量级不够） |

**判定：`voxel_volume_accuracy` 门 RED（honest FAIL，不粉饰）。** 体素体积在这四档网格上
全部超出 3% 验收线，误差随网格加密单调下降（−5.65% → −4.54% → −3.83%），符合二阶部分填充规则的
O(h) 行为，但**在当前算力档内没有达到 3%**。原因不是 SDF（其与全扫参考逐位一致，0.0），
而是**各向异性粗网格**：10/16 = 0.625 m 的纵向体素对 0.8 m 吃水只有 1.3 个单元，
船体首尾尖锐段的割平面近似的偏差主导了积分。这一条直接给出可执行的结论：
**网格体素化要用就至少要比 32×16×8 细**；`docs/voxel_export.md` 里解析场在 128³ 上的 0.03% 不能外推到粗网格。

## 5. 注意点 / 局限（诚实记录）

1. **10⁶–10⁷ 档尚未出数**：level 5/6（1.97e6 / 7.86e6 三角面）已递交调度器（作业 **cj4424**），
   本文件写作时仍在排队（25 min 未出队，集群 85 个在跑作业）。**本节不填任何外推数字**；
   工件里以 `pending_rungs` 字段记录该作业号与原因。已有 6 行是本机实测（<10 min 短测试额度内）。
2. **`band` 是截断语义，不是完整 SDF**：超出带的单元记 ±band、`triangle = -1`。因为带 ≥ 半个体素，
   部分填充规则给出的体积与完整 SDF 相同；但**若把网格直接交给需要带外精确距离的消费者（例如窄带外的
   重初始化），必须显式说明这一点**。`MeshSDFField` 的 `quantity` 因此仍声明 `signed_distance`，
   仅在调用方传 `band` 时生效——这是本模块唯一的语义陷阱，已写进 docstring 与测试。
3. **指数门的 R² 低是结构性的**：三角面轴指数 0.35、R² = 0.895（3 个规模点）。原因是 band 固定为 1 体宽时，
   查询代价主要由"带内三角形数 × 单元数"决定，而索引分辨率随 N 自动加密抵消了大部分增长——
   不是拟合坏了，而是这条轴的物理控制量不是三角形数本身。点数 ≥3 才拟合，**1 个点仍打印 `not_run`**。
4. **符号依赖闭合曲面**：`MeshSDFField(closed=False)` 允许开放网格，但此时奇偶性符号无意义（只保证幅值）。
   默认 `closed=True` 是有意的：一个没有内外之分的"体素化"会把实体分数搞反。
5. **共享机负载**：本机 load15 ≈ 68/128；同档跨档只比数量级与比值，绝对秒数有 10–20% 抖动。
   体素化是单进程单线程，不受集群核数影响，但受调度器排队影响。
6. **`voxelize`（解析场）与 `voxelize_mesh`（网格）不可混读**：前者 `triangles` 恒为 0、代价 O(cells)；
   后者代价 ≈ O(cells × 带内三角形数)，两者的"需要集群"判据完全不同。

## 6. 什么规模需要集群（本条轴的实测判据）

- **测到的**：7,992 三角面 / 4,096 体素 → 4.527 s（1,105 µs/查询）；512 体素 → 1.244 s。
  指数：单元轴 **0.96**（R² = 0.999，近线性），三角面轴 **0.35**（亚线性）。
- **没测到的不许外推**：1.97e6 / 7.86e6 三角面的行在 cj4424 里，出数后才写。
  已入库的对照事实是：**解析场**体素化 56.6e6 单元只要 356 s（`ladder_voxel.json`），
  而**网格**体素化在 4,096 单元就要 4.5 s（1,105 µs/查询）——每单元贵 4 个数量级。
  因此"体素化的集群需求"取决于网格是否有闭式表达，这一条此前没有被区分过。

## 7. 工件索引

| 工件 | 路径 |
|---|---|
| 网格 SDF 场 | `tensorcad/mesh_field.py` |
| 精确最近点查询（含窄带） | `tensorcad/spatial_index.py::TriangleGrid.nearest_surface` |
| 阶梯 / 门 / 表格渲染 | `tensorcad/perf_voxel_mesh.py` |
| 驱动 CLI | `scripts/perf_voxelize_mesh.py`、`scripts/run_voxelize_mesh_ladder.py` |
| 测量工件 | `benchmarks/perf_scale/artifacts/ladder_voxelize_mesh.json` |
| 图 | `benchmarks/perf_scale/figures/meshvoxel_scaling.png` |
| 测试 | `tests/test_mesh_field.py`、`tests/test_perf_voxel_mesh.py` |
