# 六格式网格矩阵「进」向验证 — 外形面四边形笼 / 二维 SU2 / 多块 CGNS

**锚点**: `GAP-CAD-B2-MESH-FORMAT-MATRIX` —「meshio 全格式矩阵进出（STL/OBJ/VTK/SU2/CGNS/Gmsh）」
**验收**: 往返保真：体积/拓扑/取向不变（6 格式×3 几何=18 组合全绿表）
**单号**: m20943　**仓**: TensorCAD @ `agent/m20943`　**日期**: 2026-11-14
**复现**: `python3 tools/mesh_intake_matrix.py`（双向矩阵）、`python3 tools/intake_baseline_probe.py`（前后对照）、
`python3 -m pytest tests/test_mesh_formats_foreign_intake.py -q`（24 passed）

---

## 1. 案例介绍

### 1.1 为什么"18 组合全绿"还不够

验收表要求 6 格式 × 3 几何 = 18 组合全绿。既有的两张 18 格表（`tools/native_six_format_matrix.py`、
`scripts/meshio_matrix18.py`）都只覆盖**出**（写成文件再读回），而且读回用的是**我们自己写出的文件**：
这样证明的是"写器与读器互相自洽"，不是"外部来的文件能进得来"。

而船体几何在工程链上真正到达 TensorCAD 的形态是：CAD 内核 / Pointwise / Salome 导出的**四边形主导面笼**、
`NDIME=2` 的**二维 SU2**、按 patch 分成**多个 `Elements_t` 块的 CGNS**、以及夹带**线单元**的混合文件。
这些"进"向路径此前没有被任何表覆盖——实测发现它们与"出"向结论并不对称（见 §3.3）。

### 1.2 几何与算例（本单全部解析可判）

| 算例 | 描述 | 顶点 | 单元 | 解析基准 |
|---|---|---|---|---|
| `tri_cage` | 单位立方体，12 三角形，外法向 | 8 | 12 tri | V = 1.0 m³（精确） |
| `quad_cage` | 同一立方体，**6 四边形**（外法向） | 8 | 6 quad | V = 1.0 m³（精确） |
| `planar_ndime2` | `NDIME=2` 平面正方形（SU2 二维） | 4 | 2 tri | A = 1.0 m² |
| `marker_only_quads` | `NELEM=0` + 一个 `WALL` marker（四边形） | 4 | 1 quad | A = 1.0 m² |
| `mixed_line_elements` | 1 三角形 + 1 线单元的混合文件 | 3 | 1 tri + 1 line | 应拒绝静默丢线 |
| `multiblock_cgns` | `TRI_3` 块 + `QUAD_4` 块（两个 `Elements_t`） | 8 | 2 tri + 2 quad | 面数 2 + 2×2 = 6 |

### 1.3 六格式的"进"向语法要点

| 格式 | 四边形语法 | 二维语法 | 多块语法 | 备注 |
|---|---|---|---|---|
| STL | 无（三角形汤） | 无 | 无 | float32、无连通性，拓扑需焊接后统计 |
| OBJ | `f a b c d` | — | 单文件 | 既有实现已扇形剖分 |
| VTK | `CELLS` 4 节点 + `CELL_TYPES 9` | — | — | 既有实现已按多边形扇形剖分 |
| Gmsh MSH 2.2 | `$Elements` etype **3** | — | $Elements 内混排 | **本单新增** |
| SU2 | `NELEM` type **9** | `NDIME= 2` | `NELEM` + `NMARK` | **本单新增** |
| CGNS | `ElementType_t` **7 (QUAD_4)** | — | 多个 `Elements_t` | **本单新增** |

## 2. 方法

### 2.1 工具链
`tensorcad.mesh_formats`（本仓原生六格式读写器，零第三方依赖 + h5py）、`numpy`、`h5py 3.15`、
`meshio 5.3.5`（独立第二实现，仅用于交叉对照）、`h5dump`（独立 HDF5 结构检查）、`pytest 9`、`matplotlib`。

### 2.2 流程
1. **外来文件按语法手写**（`tools/mesh_intake_matrix.py` 的 `write_*`），不经本仓任何写器——
   否则测的是自己的写器。CGNS 的 `QUAD_4` 块尤其必须手写：本仓没有任何写器会产出四边形块。
2. **读入并度量不变量**：体积（散度定理）/ 面积 / Euler 特征数 / 边界边数 / 有向面多重集。
3. **前后对照**：`git show HEAD:tensorcad/mesh_formats.py` 取出改动前模块，以独立模块载入
   （`mesh_formats` 无包内相对导入，可直接 exec），同一批 fixture 同时喂给"改前"和"改后"读器。
4. **交叉对照**：同一外来文件再经 meshio 独立导入器读入，比较三角化结果的面多重集（OBJ/VTK/Gmsh/SU2 四格式一致）。
5. **不静默原则**：非面单元（线/点）必须经 `notes=[]` 显式确认，否则抛错；非面内容（SU2 体单元、
   完全无面单元的文件）抛错并引用类型/计数。

### 2.3 判据
`compare_meshes` 的 `pass`（顶点双射 + 面数相等 + 面多重集相等 + **有向**面多重集相等 + 体积同号 + 闭合），
外加体积相对误差阈值：实现门 1e-3、owner 精度标准 15%。四边形笼另加"每个四边形必须产出 2 个三角形"。

## 3. 计算结果

### 3.1 轴 A：往返矩阵 18 格（验收表）

`artifacts/m20943_foreign_intake/intake_matrix.json`，**18/18 全绿**，最大体积相对误差 **4.58e-9**（STL 球面 float32）。

| 格式 | cube | icosphere2 | torus64x32 | 最大 \|ΔV/V\| |
|---|---|---|---|---|
| STL | ok | ok | ok | 4.58e-9 |
| OBJ | ok | ok | ok | 0 |
| VTK | ok | ok | ok | 0 |
| SU2 | ok | ok | ok | 0 |
| CGNS | ok | ok | ok | 0 |
| Gmsh | ok | ok | ok | 0 |

（Euler：cube/球 = 2，torus = 0；边界边全 0 —— 拓扑与取向逐格不变，明细见 `matrix.md`。）

### 3.2 轴 B：外来文件"进"向矩阵 15 格

14 格读取成功 + 1 格**按设计拒绝**（混合线单元文件在未显式确认时不允许静默丢线）。

| 算例 | 格式 | 结果 | 面数 | 体积/面积 | Euler | 边界边 |
|---|---|---|---|---|---|---|
| tri_cage | STL/OBJ/VTK/SU2/CGNS/Gmsh | ok×6 | 12 | 1.000000 | 2* | 0* |
| **quad_cage** | OBJ/VTK/**Gmsh**/**SU2**/**CGNS** | ok×5 | **12** | **1.000000** | **2** | **0** |
| **planar_ndime2** | SU2 | ok | 2 | A = 1.000000 | 1 | 4 |
| **marker_only_quads** | SU2 | ok | 2 | A = 1.000000 | 1 | 4 |
| mixed_line_elements | SU2 | **refused**（默认）/ ok（`notes=[]`） | — | — | — | — |
| **multiblock_cgns** | CGNS | ok | **6**（两块都消费） | 0.333333（开放补片） | 1 | 8 |

\* STL 是三角形汤，拓扑量在焊接后统计（`welded_for_topology: true`）。

### 3.3 前后对照：6 项能力从"抛错"变为"正确读入"，0 回归

`artifacts/m20943_foreign_intake/before_after.json`（baseline = `git show HEAD:` 的真实旧代码）：

| 格式 | 算例 | 改前 | 改后 | 改后结果 |
|---|---|---|---|---|
| Gmsh | quad_cage | **raised** `faces must be (m,3), got (0,)`（静默丢四边形→下游莫名报错） | ok | 12 面，V = 1.0 |
| SU2 | quad_cage | **raised** `element type 9 is not a triangle` | ok | 12 面，V = 1.0 |
| CGNS | quad_cage | **raised** SIDS 检查拒绝 `QUAD_4` | ok | 12 面，V = 1.0 |
| SU2 | planar_ndime2 | **raised** `points must be (n,3), got (4, 2)` | ok | 2 面，A = 1.0，z ≡ 0 |
| SU2 | marker_only_quads | **raised** | ok | 2 面，A = 1.0 |
| CGNS | multiblock_cgns | **raised**（只认单块 `TRI_3`） | ok | 6 面（两块都消费） |

`"newly_accepted"`: 上述 6 项；`"regressions"`: **[]**。

### 3.4 与独立第二实现（meshio）的一致性

OBJ / VTK / Gmsh / SU2 四格式的四边形笼，两个独立导入器给出**相同面多重集**（12/12 面，V 均 1.0，Euler 均 2）：
扇形剖分取对角 `0-2`，与 `tensorcad.meshio_exchange.meshio_to_mesh` 的选择一致，因此同一物理文件经任一入口读入结果相同。
CGNS 未做此对照：meshio 的自有 CGNS 读器在本机对**自己写出的三角形文件**即报
`KeyError: object ' data' doesn't exist`（见 §5 局限），不是可用参照。

### 3.5 图

* `figures/roundtrip_volume_error.png` — 18 格往返体积相对误差柱状图（对数轴），带 1e-3 实现门与 15% owner 线。
* `figures/intake_matrix.png` — 6 格式 × 6 外来算例的"进"向接受矩阵热图（格内标注面数/拒绝）。
* `figures/intake_faces_and_cube_volume.png` — 四边形笼 → 三角形数（对照 6×2=12）与各格式立方体体积保真。

## 4. 误差判定（逐项，明确回答 <3%）

| 关键量 | 基准（解析） | 实测 | 相对误差 | <3%？ |
|---|---|---|---|---|
| quad_cage 体积（Gmsh/SU2/CGNS/OBJ/VTK） | 1.0 m³ | 1.0 m³ | **0.0** | **是** |
| quad_cage 面数 | 6 quad × 2 = 12 | 12 | 0 | **是** |
| quad_cage Euler / 边界边 | 2 / 0 | 2 / 0 | 0 | **是** |
| planar SU2 面积 | 1.0 m² | 1.0 m² | ≈1e-16 | **是** |
| planar SU2 平面性 | z ≡ 0 | z ≡ 0 | 0 | **是** |
| multiblock CGNS 面数 | 2 + 2×2 = 6 | 6 | 0 | **是** |
| 往返矩阵最大体积误差（18 格） | 各自解析值 | — | 4.58e-9 | **是** |

结论：本单全部关键量相对误差 ≤ 4.6e-9，**远优于 3%**；无"无基准"项。

## 5. 注意点与局限（honest）

1. **fixture 是手写语法精确文件，不是工业软件真实导出物**。它们覆盖了四边形笼/二维/多块/混合单元
   这四类**语法形态**，但没有覆盖真实船体 patch 的尺度（±60 m 偏移坐标、退化面片、重复节点）。
   真实几何的这层验证仍需一次真实导出件（`WITHHELD`）。
2. **CGNS 只做到 SIDS 结构自检 + h5dump 结构复核**，本机无 `cgnslib`/`python-cgns`，第三方 CGNS 读器
   全链条一致性仍 `WITHHELD`。`h5dump -n` 已核实两块 `Elements_t`（type 5 / type 7）节点树可被独立工具枚举。
3. **STL 天生是三角形汤**：本单对其拓扑量采用焊接后统计，并在结果里显式标记 `welded_for_topology`，
   不把"焊接后的闭合"说成"文件自带连通性"。
4. **meshio 上游两处缺陷**（本机 5.3.5，可复现）：其 SU2 **写器**对表面单元抛
   `TypeError: cannot unpack non-iterable CellBlock object`；其 CGNS **读器**读自己写出的三角形文件抛
   `KeyError: object ' data' doesn't exist`（体积四面体则两者都正常）。这解释了为何 `meshio_exchange`
   过去把 SU2/CGNS 表面三角形记为"不支持"——**是 meshio 的限制，不是格式的限制**：SU2 三角形类型 5、
   CGNS `TRI_3`/`QUAD_4` 都是标准语法，本仓 `mesh_formats` 一直在写前者。该结论已据此改正文档。
5. **多边形策略是扇形剖分（对角 0-2）**，对非平面四边形只是近似；船体真实四边形若有明显翘曲，
   应在导出端处理，不在读器里悄悄"修形"。
6. 未改动任何写器，因此 §3.1 的 18 格往返基线不受影响（同一脚本现场重算，非引用旧结论）。

## 6. 变更清单

| 文件 | 变更 |
|---|---|
| `tensorcad/mesh_formats.py` | `read_gmsh` 支持四边形(etype 3)、`read_su2` 支持 `NDIME=2` 与 type 9、`read_cgns` 消费全部 `Elements_t` 块并支持 `QUAD_4`；新增 `_fan()` 统一扇形剖分；SIDS 检查器扩展为多块/多类型（含 `ElementRange` 连续覆盖校验）；非面单元须经 `notes=[]` 显式确认 |
| `tensorcad/meshio_exchange.py` | 修正"SU2/CGNS 表面三角形不支持"的错误表述（格式层面支持，限制来自 meshio 写器/读器），指向原生写读器 |
| `tests/test_mesh_formats_foreign_intake.py` | 新增 24 个测试：6 格式四边形进向、取向不被"洗白"、meshio 交叉一致、二维 SU2、marker-only、多块 CGNS、块区间不连续拒绝、线单元须确认、空面文件显式报错 |
| `tools/mesh_intake_matrix.py` | 双向矩阵脚本（往返 18 格 + 进向 15 格 + meshio 交叉 + 3 张图 + JSON/markdown） |
| `tools/intake_baseline_probe.py` | 改前/改后对照探针（从 git 取旧模块真实执行） |
| `artifacts/m20943_foreign_intake/` | `intake_matrix.json`、`matrix.md`、`before_after.json`、`files/`（外来 fixture）、`figures/` |

**测试**：`pytest tests/test_mesh_formats_foreign_intake.py` → 24 passed；
相关回归集（`test_mesh_formats` / `test_mesh_xval` / `test_meshio_format_matrix` / `test_meshio_import_cells` /
`test_volume_mesh_format_bridge` / 本单）→ **249 passed, 7 skipped**（跳过者为第三方库缺失）。