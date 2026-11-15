# TensorCAD 六格式网格往返保真矩阵 — 综合验证报告

**报告编号**: m22384  
**日期**: 2026-11-15  
**锚点**: `GAP-CAD-B2-MESH-FORMAT-MATRIX`  
**验收标准**: 往返保真：体积/拓扑/取向不变（6 格式×3 几何=18 组合全绿表）  
**仓库**: TensorCAD @ `agent/m22384`  
**状态**: ✅ 验证通过，全量测试绿（249 passed, 0 failed）

---

## 1. 案例介绍

### 1.1 物理背景

在船舶总体性能计算管线中，船体几何模型需要在不同系统之间可靠交换：
- **CAD 系统** → 网格生成器 → **CFD 求解器**（SU2/OpenFOAM）→ **结构分析工具**
- 每个环节支持的网格格式不同，格式转换必须保持几何保真性

本工作验证 TensorCAD 原生六格式读写器的往返保真性：同一几何体写入某格式后读回，
**体积、拓扑、取向**必须保持不变。这是船舶设计管线中几何交换的基础能力。

### 1.2 测试几何体

| 几何体 | 类型 | 拓扑 | 解析体积 | 顶点数 | 面数 |
|--------|------|------|----------|--------|------|
| cube | 凸多面体 | genus 0, Euler=2 | 8.0 m³（精确） | 8 | 12 |
| icosphere(1.0, sub=2) | 球面近似 | genus 0, Euler=2 | 4.047005 m³ | 162 | 320 |
| torus(1.0, 0.3, 64×32) | 环面 | genus 1, Euler=0 | 1.762302 m³ | 2048 | 4096 |

**选择理由**：
- **cube**：最简单的封闭多面体，解析体积精确，是格式保真的基准
- **icosphere**：曲面离散化代表，顶点分布均匀，检验浮点精度
- **torus**：非平凡拓扑（genus 1），检验格式对孔洞几何的处理

### 1.3 六格式特性

| 格式 | 扩展名 | 精度 | 拓扑存储 | 特殊约束 |
|------|--------|------|----------|----------|
| **STL** | .stl | float32（格式定义） | facet soup，无拓扑 | 法向量需与顶点绕向一致 |
| **OBJ** | .obj | float64（17位有效数字） | 显式索引 | 1-based 索引 |
| **VTK** | .vtk | float64 | UNSTRUCTURED_GRID | meshio 只读此 dialect |
| **SU2** | .su2 | float64 | NELEM 块 | 三角形类型 ID=5 |
| **CGNS** | .cgns | float64 | SIDS/HDF5 | 需 CGNSLibraryVersion 节点 |
| **Gmsh** | .msh | float64 | MSH 2.2 ASCII | 三角形类型 ID=2 |

---

## 2. 方法

### 2.1 工具链

| 组件 | 模块 | 说明 |
|------|------|------|
| 写入器 | `tensorcad.mesh_formats.write()` | 原生六格式写入器，不依赖 meshio |
| 读取器 | `tensorcad.mesh_formats.read()` | 原生六格式读取器，从格式规范编写 |
| 比较器 | `tensorcad.mesh_formats.compare_meshes()` | 几何+拓扑+取向比较 |
| 验证器 | `tensorcad.mesh_formats.sids_violations()` | CGNS SIDS 合规性检查 |

### 2.2 往返流程

```
原始 TriMesh → write(fmt, path) → 磁盘文件 → read(path, fmt) → 读回 TriMesh → compare_meshes()
```

### 2.3 比较不变量

`compare_meshes()` 验证以下不变量：
1. **顶点双射**：读回顶点集与原始顶点集存在一一映射（float64 格式要求精确匹配）
2. **面数相等**：输入/输出面数相同
3. **面多重集相等**：面的索引集合相同
4. **有向面多重集相等**：面的取向（顶点顺序）相同
5. **体积符号相等**：体积正负号不变（外向法向保持）
6. **封闭性**：读回网格无边界边

### 2.4 容差设置

- **STL**：`tol=1e-5`（float32 精度限制，格式本身定义为 float32）
- **其他格式**：`tol=0.0`（要求精确双射，float64 往返无损）

### 2.5 三重验证管线

本锚点的验证覆盖三条独立路径：

| 路径 | 测试文件 | 组合数 | 说明 |
|------|----------|--------|------|
| **原生路径** | `test_mesh_formats.py` | 18 | 自有读写器，不经第三方库 |
| **meshio 路径** | `test_meshio_format_matrix.py` | 24 | 经 meshio 库的往返 |
| **规格级审计** | `test_format_conformance.py` | 36 | 规格级解析器，独立验证 |
| **外来摄入** | `test_mesh_formats_foreign_intake.py` | 24 | 四边形/多边形网格摄入 |

---

## 3. 计算结果

### 3.1 全矩阵结果（18 组合）

| 几何体 | 格式 | 通过 | 体积(参考) | 体积(读回) | 相对误差 | 顶点(入→出) | 面(入→出) | Euler | 文件大小 |
|--------|------|------|------------|------------|----------|-------------|-----------|-------|----------|
| cube | stl | ✅ | 8.000000 | 8.000000 | 0.00e+00 | 8→8 | 12→12 | 2→2 | 684 B |
| cube | obj | ✅ | 8.000000 | 8.000000 | 0.00e+00 | 8→8 | 12→12 | 2→2 | 200 B |
| cube | vtk | ✅ | 8.000000 | 8.000000 | 0.00e+00 | 8→8 | 12→12 | 2→2 | 321 B |
| cube | su2 | ✅ | 8.000000 | 8.000000 | 0.00e+00 | 8→8 | 12→12 | 2→2 | 245 B |
| cube | cgns | ✅ | 8.000000 | 8.000000 | 0.00e+00 | 8→8 | 12→12 | 2→2 | 20,040 B |
| cube | gmsh | ✅ | 8.000000 | 8.000000 | 0.00e+00 | 8→8 | 12→12 | 2→2 | 303 B |
| icosphere2 | stl | ✅ | 4.047005 | 4.047005 | 4.58e-09 | 162→162 | 320→320 | 2→2 | 16,084 B |
| icosphere2 | obj | ✅ | 4.047005 | 4.047005 | 0.00e+00 | 162→162 | 320→320 | 2→2 | 12,722 B |
| icosphere2 | vtk | ✅ | 4.047005 | 4.047005 | 0.00e+00 | 162→162 | 320→320 | 2→2 | 13,146 B |
| icosphere2 | su2 | ✅ | 4.047005 | 4.047005 | 0.00e+00 | 162→162 | 320→320 | 2→2 | 13,067 B |
| icosphere2 | cgns | ✅ | 4.047005 | 4.047005 | 0.00e+00 | 162→162 | 320→320 | 2→2 | 27,992 B |
| icosphere2 | gmsh | ✅ | 4.047005 | 4.047005 | 0.00e+00 | 162→162 | 320→320 | 2→2 | 14,799 B |
| torus64x32 | stl | ✅ | 1.762302 | 1.762302 | 6.91e-10 | 2048→2048 | 4096→4096 | 0→0 | 204,884 B |
| torus64x32 | obj | ✅ | 1.762302 | 1.762302 | 0.00e+00 | 2048→2048 | 4096→4096 | 0→0 | 190,608 B |
| torus64x32 | vtk | ✅ | 1.762302 | 1.762302 | 0.00e+00 | 2048→2048 | 4096→4096 | 0→0 | 194,809 B |
| torus64x32 | su2 | ✅ | 1.762302 | 1.762302 | 0.00e+00 | 2048→2048 | 4096→4096 | 0→0 | 194,728 B |
| torus64x32 | cgns | ✅ | 1.762302 | 1.762302 | 0.00e+00 | 2048→2048 | 4096→4096 | 0→0 | 118,376 B |
| torus64x32 | gmsh | ✅ | 1.762302 | 1.762302 | 0.00e+00 | 2048→2048 | 4096→4096 | 0→0 | 223,261 B |

### 3.2 汇总统计

| 指标 | 值 |
|------|-----|
| 总组合数 | 18 |
| 通过数 | **18** |
| 失败数 | **0** |
| 最坏相对体积误差 | **4.58×10⁻⁹**（icosphere2/STL） |
| 拓扑保持（全部 18） | ✅ |
| 全部在验收阈值内（<1e-5） | ✅ |

### 3.3 测试套件汇总

| 测试文件 | 测试数 | 通过 | 失败 | 跳过 |
|----------|--------|------|------|------|
| `test_mesh_formats.py` | 38 | 36 | 0 | 2 |
| `test_meshio_format_matrix.py` | 145 | 141 | 0 | 4 |
| `test_format_conformance.py` | 48 | 48 | 0 | 0 |
| `test_mesh_formats_foreign_intake.py` | 24 | 24 | 0 | 0 |
| **合计** | **255** | **249** | **0** | **6** |

---

## 4. 结果图

### 4.1 体积误差热图

![Volume Relative Error Heatmap](../artifacts/m22384/figures/matrix_heatmap.png)

*图 1：六格式×三几何体往返体积相对误差热图（symlog 刻度）。STL 因 float32 精度限制显示 ~10⁻⁹ 量级误差，其余格式为零误差。*

### 4.2 拓扑不变量验证

![Topology Verification](../artifacts/m22384/figures/topology_verification.png)

*图 2：四项拓扑不变量（顶点数、面数、Euler 示性数、封闭性）在 18 组合中全部保持。*

### 4.3 文件大小比较

![File Size Comparison](../artifacts/m22384/figures/file_size_comparison.png)

*图 3：各格式输出文件大小（对数刻度）。CGNS 因 HDF5 元数据开销最大，OBJ 最紧凑。*

### 4.4 体积对比

![Volume Comparison](../artifacts/m22384/figures/volume_comparison.png)

*图 4：各格式往返体积与参考值对比（红色虚线为解析参考值）。所有格式均在参考值附近，STL 有微小 float32 偏差。*

---

## 5. 注意点

### 5.1 STL 精度限制

STL 格式本身定义为 float32，因此往返后顶点坐标有 ~7 位有效数字的精度损失。
对于 icosphere2（半径 1.0），最坏相对误差为 4.58×10⁻⁹，远低于验收阈值 1e-5。
对于大尺度模型（如 100 m 船体），float32 的绝对误差约为 10⁻⁵ m，在工程可接受范围内。

### 5.2 STL 是 facet soup

STL 不存储拓扑信息，每个面独立存储三个顶点坐标。读回后需要焊接（weld）才能恢复拓扑。
本模块的 `weld()` 函数通过舍入不变量焊接，decimals=12 时能正确恢复原始拓扑。

### 5.3 CGNS SIDS 合规性

本模块写入的 CGNS 文件包含完整的 SIDS 元数据：
- `/CGNSLibraryVersion` 节点（R4 类型，值 4.3）
- `/Base` 节点（CGNSBase_t 标签）
- `/Base/Zone1/GridCoordinates` 节点（CoordinateX/Y/Z）
- `/Base/Zone1/GridElements` 节点（ElementType=5 for TRI_3）

`sids_violations()` 检查确认零违规。但本环境无 cgnslib/ParaView，第三方 CGNS 库的兼容性
按 **WITHHELD** 记。

### 5.4 SU2 格式选择

SU2 写入器默认使用 `nelem_mode="volume"`，将三角形放在 NELEM 块中。这是 meshio 可读的 dialect。
`nelem_mode="boundary"` 也是合法 SU2，但 meshio 的 SU2 读取器有缺陷（`StopIteration`），
因此不推荐。

### 5.5 VTK dialect 选择

VTK 写入器默认使用 `UNSTRUCTURED_GRID` 而非 `POLYDATA`。虽然两者都是合法 VTK，
但 meshio 的 legacy VTK 读取器只接受 `UNSTRUCTURED_GRID`/`STRUCTURED_*`。
`POLYDATA` 可通过 `dataset="POLYDATA"` 参数使用。

### 5.6 四边形/多边形摄入

外来 CAD 文件通常是四边形主导的网格。本模块的读取器支持：
- 四边形扇形三角化（0-2 对角线，与 meshio 路径一致）
- 多边形扇形三角化（从第一个顶点出发）
- 平面 SU2（NDIME=2）自动补零 Z 坐标

### 5.7 适用范围

- 本验证针对**表面三角形网格**（TriMesh）
- 体积网格（四面体）的 SU2/CGNS 往返由 `meshio_exchange` 模块单独验证
- 四边形/多边形网格的导入由 `format_conformance` 模块处理（扇形三角化）

---

## 6. 误差判定

### 6.1 关键量 vs 基准

| 关键量 | 基准值 | 实测最坏值 | 相对误差 | 是否 <3% |
|--------|--------|------------|----------|----------|
| cube 体积 | 8.0 m³ | 8.0 m³ | 0.00e+00 | ✅ 是 |
| icosphere2 体积 | 4.047005 m³ | 4.047005 m³ | 4.58e-09 | ✅ 是 |
| torus64x32 体积 | 1.762302 m³ | 1.762302 m³ | 6.91e-10 | ✅ 是 |
| 顶点数保持 | 精确 | 精确 | 0 | ✅ 是 |
| 面数保持 | 精确 | 精确 | 0 | ✅ 是 |
| Euler 示性数 | 精确 | 精确 | 0 | ✅ 是 |
| 封闭性 | 全封闭 | 全封闭 | — | ✅ 是 |

### 6.2 结论

**所有 18 组合的体积相对误差均远小于 3%（最坏 4.58×10⁻⁹ ≈ 4.58×10⁻⁷%），
拓扑和取向完全保持。验收标准达成。**

---

## 附录：工件索引

| 工件 | 路径 |
|------|------|
| 综合验证报告 | `docs/mesh_format_matrix_comprehensive_report.md` |
| 结果 JSON | `artifacts/m22384/results.json` |
| 误差热图 | `artifacts/m22384/figures/matrix_heatmap.png` |
| 拓扑验证图 | `artifacts/m22384/figures/topology_verification.png` |
| 文件大小图 | `artifacts/m22384/figures/file_size_comparison.png` |
| 体积对比图 | `artifacts/m22384/figures/volume_comparison.png` |
| 生成脚本 | `scripts/generate_matrix_report_figures.py` |
| 测试用例 | `tests/test_mesh_formats.py` |

## 复现

```bash
# 全部相关测试（249 passed）
python -m pytest tests/test_mesh_formats.py tests/test_meshio_format_matrix.py \
  tests/test_format_conformance.py tests/test_mesh_formats_foreign_intake.py -q

# 生成报告工件和图表
python scripts/generate_matrix_report_figures.py
```
