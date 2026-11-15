# GAP-CAD-B2-MESH-FORMAT-MATRIX — 扩展验证报告（m21535）

**锚点**: `GAP-CAD-B2-MESH-FORMAT-MATRIX`
**验收**: 往返保真：体积/拓扑/取向不变（6 格式×3 几何=18 组合全绿表）
**Mission**: m21535
**日期**: 2026-11-15
**仓库**: TensorCAD @ agent/m21535
**状态**: ✅ 验证通过，22/22 适用组合全绿 + 2 N/A（正确拒绝）

---

## ① 案例介绍

### 物理背景
船舶总体性能设计中，船体几何需要在 CAD 系统、CFD 求解器（SU2/OpenFOAM）、结构分析工具
之间可靠交换。meshio 库提供统一的 Python 接口读写多种网格格式。本单验证 TensorCAD 表面网格
经 meshio 往返转换后不失真，并扩展验证到第 4 种几何——环面（torus，genus 1，非星形体）。

### 几何基准（4 几何）

| 几何 | 类型 | 顶点/单元 | 解析体积 [m³] | 来源 |
|---|---|---|---|---|
| cube | 立方体 (edge=2.0) | 8 / 12 | 8.0 | 解析精确 |
| sphere | 球面离散 (icosphere, r=1.0, subdiv=2) | 162 / 320 | 4π/3 ≈ 4.18879 | 解析精确 |
| torus | 环面 (genus 1, R=1.0, r=0.3) | 2048 / 4096 | 2π²Rr² ≈ 1.77653 | 解析精确 |
| wigley_hull | Wigley 抛物船型 (L=10, B=1, T=0.625) | 1001 / 1998 | 4·L·B·T/9 ≈ 2.77778 | `tensorcad.wigley` 解析闭式 |

### 格式矩阵（6 格式）

| 格式 | 扩展名 | 表面三角 | 四面体 | 精度 |
|---|---|---|---|---|
| STL | .stl | ✅ REAL | — | float32（格式定义）/ ASCII（meshio 实测） |
| OBJ | .obj | ✅ REAL | — | float64 |
| VTK | .vtk | ✅ REAL | — | float64 |
| Gmsh | .msh | ✅ REAL | — | float64 |
| SU2 | .su2 | NOT_SUPPORTED | ✅ REAL | float64 |
| CGNS | .cgns | NOT_SUPPORTED | ✅ REAL | float64 |

体积格式（SU2/CGNS）先将几何四面体化（star tessellation）再导出。
**torus 是非星形体**，star tessellation 正确拒绝它——这是预期行为，不是缺陷。

---

## ② 方法

### 验证管线

**pytest 单元测试**（`tests/test_meshio_format_matrix.py`）
- 6 格式 × 4 几何 = 24 组合（其中 2 个 N/A），验证体积/拓扑/取向
- 141 项测试，全部通过

**规格级独立审计**（`tests/test_format_conformance.py`）
- 用 `tensorcad.format_conformance` 的规格级解析器（不经 meshio）
- 48 项测试，全部通过

**扩展审计工具**（`tools/meshio_matrix_expanded_audit.py`）
- 新增 torus 几何，覆盖 genus-1 拓扑
- 生成 matplotlib 可视化证据

### 判据
- 体积相对误差 < 1e-12（STL 因 float32 存储放宽到 1e-3）
- 顶点数/单元数一致
- 有向体积为正向（outward）
- 表面闭合（边价数全为 2）

---

## ③ 计算结果

### 24 组合审计表

| Format | Geometry | Topology | V_orig [m³] | V_rt [m³] | rel. err | Topology | Orientation | Status |
|--------|----------|----------|-------------|-----------|----------|----------|-------------|--------|
| stl | cube | triangle (8v/12c) | 8.000000000 | 8.000000000 | 0.00e+00 | ok | ok | **PASS** |
| stl | sphere | triangle (162v/320c) | 4.047005367 | 4.047005367 | 0.00e+00 | ok | ok | **PASS** |
| stl | torus | triangle (2048v/4096c) | 1.762301535 | 1.762301535 | 0.00e+00 | ok | ok | **PASS** |
| stl | wigley_hull | triangle (1001v/1998c) | 2.771222150 | 2.771222150 | 0.00e+00 | ok | ok | **PASS** |
| obj | cube | triangle (8v/12c) | 8.000000000 | 8.000000000 | 0.00e+00 | ok | ok | **PASS** |
| obj | sphere | triangle (162v/320c) | 4.047005367 | 4.047005367 | 0.00e+00 | ok | ok | **PASS** |
| obj | torus | triangle (2048v/4096c) | 1.762301535 | 1.762301535 | 0.00e+00 | ok | ok | **PASS** |
| obj | wigley_hull | triangle (1001v/1998c) | 2.771222150 | 2.771222150 | 0.00e+00 | ok | ok | **PASS** |
| vtk | cube | triangle (8v/12c) | 8.000000000 | 8.000000000 | 0.00e+00 | ok | ok | **PASS** |
| vtk | sphere | triangle (162v/320c) | 4.047005367 | 4.047005367 | 0.00e+00 | ok | ok | **PASS** |
| vtk | torus | triangle (2048v/4096c) | 1.762301535 | 1.762301535 | 0.00e+00 | ok | ok | **PASS** |
| vtk | wigley_hull | triangle (1001v/1998c) | 2.771222150 | 2.771222150 | 0.00e+00 | ok | ok | **PASS** |
| gmsh | cube | triangle (8v/12c) | 8.000000000 | 8.000000000 | 0.00e+00 | ok | ok | **PASS** |
| gmsh | sphere | triangle (162v/320c) | 4.047005367 | 4.047005367 | 0.00e+00 | ok | ok | **PASS** |
| gmsh | torus | triangle (2048v/4096c) | 1.762301535 | 1.762301535 | 0.00e+00 | ok | ok | **PASS** |
| gmsh | wigley_hull | triangle (1001v/1998c) | 2.771222150 | 2.771222150 | 0.00e+00 | ok | ok | **PASS** |
| su2 | cube | tetra (9v/12c) | 8.000000000 | 8.000000000 | 0.00e+00 | ok | ok | **PASS** |
| su2 | sphere | tetra (163v/320c) | 4.047005367 | 4.047005367 | 0.00e+00 | ok | ok | **PASS** |
| su2 | torus | N/A | 1.762301535 | — | — | N/A | N/A | **N/A** |
| su2 | wigley_hull | tetra (1002v/1996c) | 2.771222150 | 2.771222150 | 0.00e+00 | ok | ok | **PASS** |
| cgns | cube | tetra (9v/12c) | 8.000000000 | 8.000000000 | 0.00e+00 | ok | ok | **PASS** |
| cgns | sphere | tetra (163v/320c) | 4.047005367 | 4.047005367 | 0.00e+00 | ok | ok | **PASS** |
| cgns | torus | N/A | 1.762301535 | — | — | N/A | N/A | **N/A** |
| cgns | wigley_hull | tetra (1002v/1996c) | 2.771222150 | 2.771222150 | 0.00e+00 | ok | ok | **PASS** |

**22/22 适用组合全绿，2 N/A（torus + 体积格式，正确拒绝）。**

### 解析体积对比

| 几何 | 解析 V [m³] | 网格 V [m³] | 离散误差 % |
|---|---|---|---|
| cube | 8.000000 | 8.000000 | 0.0000 % |
| sphere | 4.188790 | 4.047005 | 3.3849 % |
| torus | 1.776529 | 1.762302 | 0.8008 % |
| wigley_hull | 2.777778 | 2.771222 | 0.2360 % |

sphere 的 3.38% 离散误差来自 icosphere subdiv=2 的粗糙离散化（320 面片逼近球面），
不是格式损失——提高 subdiv 即可收敛。

---

## ④ 结果图

### 体积误差热力图
![volume_error_heatmap](../runs/m21535_matrix24/figures/volume_error_heatmap.png)

### 状态网格
![status_grid](../runs/m21535_matrix24/figures/status_grid.png)

### 解析 vs 离散体积
![analytic_vs_mesh_volume](../runs/m21535_matrix24/figures/analytic_vs_mesh_volume.png)

### 环面往返可视化（genus 1）
![torus_roundtrip](../runs/m21535_matrix24/figures/torus_roundtrip.png)

### Wigley 船体往返可视化
![wigley_hull_roundtrip](../runs/m21535_matrix24/figures/wigley_hull_roundtrip.png)

---

## ⑤ 注意点与局限

1. **torus 是非星形体**：从任何内部点都有部分面片被"遮挡"（孔洞在 apex 后面），
   star tessellation 正确拒绝它。这是预期行为，不是缺陷。
   表面格式（STL/OBJ/VTK/Gmsh）直接处理三角面，不受此限制。

2. **STL 二进制按规范是 float32**（体积容差取格式下限 1e-3）；
   meshio 默认写 ASCII STL（文本保留 double），故实测体积误差为 0。

3. **sphere 离散误差 3.38%** 来自 icosphere subdiv=2 的粗糙离散化，
   不是格式损失。提高 subdiv 可收敛到解析值。

4. **CGNS SIDS 合规**：`tensorcad.mesh_cgns` 写完整 SIDS 元数据（0 违规）；
   meshio 自带 CGNS writer 缺失 14 个必需节点/属性。

5. **SU2/CGNS 表面三角 NOT_SUPPORTED** 是 meshio 上游限制（meshio 5.3.5），
   不是格式本身限制。

---

## ⑥ 误差判定（是否 <3%）

| 关键量 | 基准 | 实测 | 误差 | 是否 <3% |
|---|---|---|---|---|
| 格式保真体积误差（out） | 源网格体积 | 最大 0.00e+00 | 0 % | ✅ 是 |
| 格式保真体积误差（in） | 源网格体积 | 0.00e+00 | 0 % | ✅ 是 |
| 文件体积 vs 解析体积（sphere） | 4π/3 ≈ 4.18879 | 4.04701 | 3.38 % | ⚠️ 离散误差（非格式损失） |
| 文件体积 vs 解析体积（torus） | 2π²Rr² ≈ 1.77653 | 1.76230 | 0.80 % | ✅ 是 |
| 文件体积 vs 解析体积（wigley） | 4LBT/9 ≈ 2.77778 | 2.77122 | 0.24 % | ✅ 是 |
| 拓扑（顶点/单元数） | 源网格 | 22 组合精确一致 | 0 | ✅ 是 |
| 取向（有向体积符号） | outward | 22 组合全部 outward | — | ✅ 是 |

**结论：22/22 适用组合全绿，验收通过。**

---

## 测试套件汇总

| 测试文件 | 测试数 | 通过 | 失败 | 跳过 |
|---|---|---|---|---|
| `test_meshio_format_matrix.py` | 141 | 141 | 0 | 4 |
| `test_format_conformance.py` | 48 | 48 | 0 | 0 |
| `test_mesh_formats_foreign_intake.py` | 24 | 24 | 0 | 0 |
| **合计** | **213** | **213** | **0** | **4** |

---

## 复现

```bash
# 扩展审计（24 格）
python tools/meshio_matrix_expanded_audit.py

# 全部相关测试
python -m pytest tests/test_meshio_format_matrix.py tests/test_format_conformance.py \
  tests/test_mesh_formats_foreign_intake.py -q
```

## 工件索引

| 工件 | 路径 |
|---|---|
| 扩展验证报告 | `docs/meshio_format_matrix_expanded_m21535.md` |
| 审计 JSON | `runs/m21535_matrix24/matrix24.json` |
| 审计 Markdown | `runs/m21535_matrix24/matrix24.md` |
| 体积误差热力图 | `runs/m21535_matrix24/figures/volume_error_heatmap.png` |
| 状态网格图 | `runs/m21535_matrix24/figures/status_grid.png` |
| 解析 vs 离散体积图 | `runs/m21535_matrix24/figures/analytic_vs_mesh_volume.png` |
| 环面往返图 | `runs/m21535_matrix24/figures/torus_roundtrip.png` |
| Wigley 船体往返图 | `runs/m21535_matrix24/figures/wigley_hull_roundtrip.png` |
| 审计工具 | `tools/meshio_matrix_expanded_audit.py` |
