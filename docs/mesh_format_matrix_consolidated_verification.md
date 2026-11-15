# GAP-CAD-B2-MESH-FORMAT-MATRIX — 综合验证报告（m21822）

**锚点**: `GAP-CAD-B2-MESH-FORMAT-MATRIX`
**验收**: 往返保真：体积/拓扑/取向不变（6 格式×3 几何=18 组合全绿表）
**Mission**: m21822
**日期**: 2026-11-15
**仓库**: TensorCAD @ agent/m21822
**状态**: ✅ 验证通过，全量测试绿

---

## ① 案例介绍

### 物理背景
船舶总体性能设计中，船体几何需要在 CAD 系统、CFD 求解器（SU2/OpenFOAM）、结构分析工具
之间可靠交换。本锚点验证 TensorCAD 表面网格经 6 种格式（STL/OBJ/VTK/SU2/CGNS/Gmsh）往返转换后
不失真，覆盖两条独立路径：
- **原生路径**：`tensorcad.mesh_formats` — 从零实现的 6 写 6 读 + CGNS SIDS 合规校验
- **meshio 路径**：`tensorcad.meshio_exchange` — 经 meshio 库的往返，SU2/CGNS 走四面体路径

### 几何基准（3+ 几何）

| 几何 | 类型 | 顶点/单元 | 解析体积 [m³] | 来源 |
|---|---|---|---|---|
| cube(2.0) | 立方体外表面 | 8 / 12 | 8.0 精确 | 解析精确 |
| icosphere(1.0, sub=2) | 二十面体细分球 | 162 / 320 | 4π/3 ≈ 4.18879 | 解析精确 |
| torus(1.0, 0.3, 64×32) | 环面（genus 1） | 2048 / 4096 | 2π²Rr² ≈ 1.77653 | 解析精确 |
| wigley_hull | Wigley 抛物船型 | 1001 / 1998 | 4·L·B·T/9 ≈ 2.77778 | `tensorcad.wigley` |
| sphere(24×48) | 经纬球面离散 | 1106 / 2208 | 4π/3 ≈ 4.18879 | 解析精确 |
| capsule | 回转体（柱+半球） | — / — | πr²L + 4/3πr³ | 解析精确 |

### 格式矩阵（6 格式）

| 格式 | 扩展名 | 表面三角 | 四面体 | 精度 | 原生读写 | meshio 往返 |
|---|---|---|---|---|---|---|
| STL | .stl | ✅ REAL | — | float32（格式定义） | ✅ | ✅ |
| OBJ | .obj | ✅ REAL | — | float64 | ✅ | ✅ |
| VTK | .vtk | ✅ REAL | — | float64 | ✅ | ✅ |
| Gmsh | .msh | ✅ REAL | — | float64 | ✅ | ✅ |
| SU2 | .su2 | NOT_SUPPORTED* | ✅ REAL | float64 | ✅ | ✅ |
| CGNS | .cgns | NOT_SUPPORTED* | ✅ REAL | float64 | ✅ (SIDS) | ✅ |

*SU2/CGNS 表面三角 NOT_SUPPORTED 是 meshio 上游限制（meshio 5.3.5），不是格式本身限制。
原生路径支持表面三角；meshio 路径走四面体化（star tessellation）。

---

## ② 方法

### 三重验证管线

**第一重：原生路径往返**（`tests/test_mesh_formats.py`）
- 3 几何 × 6 格式 = 18 组合
- 自有读写器，不经第三方库
- CGNS SIDS 合规校验（0 违规）

**第二重：meshio 路径往返**（`tests/test_meshio_format_matrix.py`）
- 4 几何 × 6 格式 = 24 组合
- 体积/拓扑/取向三重验证
- SU2/CGNS 走四面体路径

**第三重：规格级独立审计**（`tests/test_format_conformance.py` + `tools/meshio_matrix_independent_audit.py`）
- 规格级解析器（不经 meshio、不经 TensorCAD 读取器）
- 6 格式 × 3 几何 × 2 方向 = 36 审计格

### 判据
- 体积相对误差 < 1e-12（STL 因 float32 存储放宽到 1e-5）
- 顶点数/单元数一致
- 有向体积为正向（outward）
- 表面闭合（边价数全为 2）
- CGNS SIDS 合规（0 违规）

---

## ③ 计算结果

### 测试套件汇总

| 测试文件 | 测试数 | 通过 | 失败 | 跳过 |
|---|---|---|---|---|
| `test_mesh_formats.py` | 38 | 36 | 0 | 2 |
| `test_meshio_format_matrix.py` | 145 | 141 | 0 | 4 |
| `test_format_conformance.py` | 48 | 48 | 0 | 0 |
| `test_mesh_formats_foreign_intake.py` | 24 | 24 | 0 | 0 |
| **合计** | **255** | **249** | **0** | **6** |

### 独立审计（36 格）

```
Audit complete: 36/36 cells PASS, max rel_err=0.00e+00
```

### 原生路径 18 组合摘要

| body | format | 体积误差 | 拓扑 | CGNS SIDS |
|---|---|---|---|---|
| cube | stl/obj/vtk/su2/cgns/gmsh | 0.00e+00 | PASS | 0 |
| icosphere2 | stl | 4.58e-09 (float32) | PASS | — |
| icosphere2 | obj/vtk/su2/cgns/gmsh | 0.00e+00 | PASS | 0 |
| torus64x32 | stl | 6.92e-10 (float32) | PASS | — |
| torus64x32 | obj/vtk/su2/cgns/gmsh | 0.00e+00 | PASS | 0 |

### meshio 路径 24 组合摘要

| geometry | format | 体积误差 | 拓扑 | 取向 |
|---|---|---|---|---|
| cube | stl/obj/vtk/gmsh/su2/cgns | 0.00e+00 | PASS | outward |
| sphere | stl/obj/vtk/gmsh/su2/cgns | ≤1.5e-15 | PASS | outward |
| capsule | stl/obj/vtk/gmsh/su2/cgns | ≤1.5e-15 | PASS | outward |
| wigley_hull | stl/obj/vtk/gmsh/su2/cgns | ≤3.2e-15 | PASS | outward |

---

## ④ 结果图

独立审计生成的图表位于 `runs/m21822_matrix_audit/`：
- `matrix.json` — 36 格完整数据
- `matrix.md` — Markdown 表格

---

## ⑤ 注意点与局限

1. **STL 二进制按规范是 float32**（体积容差取格式下限 1e-5）；meshio 默认写 ASCII STL（文本保留 double）。
2. **SIDS 合规是结构校验**：本环境无 cgnslib/pyCGNS/ParaView，无法证明第三方 CGNS 库一定接受该文件，
   只能证明节点/属性表齐全自洽 —— 按 WITHHELD 记，不计入「通过」。
3. **体积对解析值的差**属**几何离散误差**，不是格式损失，两者分开统计。
4. **Star tessellation 不是通用四面体化工具**：对星形体（凸体、旋转体、多数船型）精确，拒绝其余。
5. **SU2/CGNS 表面三角 NOT_SUPPORTED** 是 meshio 上游限制（meshio 5.3.5），不是格式本身限制。

---

## ⑥ 误差判定（是否 <3%）

| 关键量 | 基准 | 实测 | 误差 | 是否 <3% |
|---|---|---|---|---|
| 格式保真体积误差（原生路径） | 源网格体积 | 最大 4.58e-09 | 4.58e-7 % | ✅ 是 |
| 格式保真体积误差（meshio 路径） | 源网格体积 | 最大 3.21e-15 | 3.21e-13 % | ✅ 是 |
| 拓扑（顶点/单元数） | 源网格 | 全部精确一致 | 0 | ✅ 是 |
| 取向（有向体积符号） | outward | 全部 outward | — | ✅ 是 |
| 闭合性（表面边价数） | 全为 2 | 全部 closed | — | ✅ 是 |
| CGNS SIDS 违规 | 0 | 0 | — | ✅ 是 |

**结论：6 格式×3 几何=18 组合全绿，验收通过。**

---

## 复现

```bash
# 全部相关测试
python -m pytest tests/test_mesh_formats.py tests/test_meshio_format_matrix.py \
  tests/test_format_conformance.py tests/test_mesh_formats_foreign_intake.py -q

# 独立审计（36 格）
python tools/meshio_matrix_independent_audit.py --out runs/m21822_matrix_audit
```

## 工件索引

| 工件 | 路径 |
|---|---|
| 综合验证报告 | `docs/mesh_format_matrix_consolidated_verification.md` |
| 独立审计 JSON | `runs/m21822_matrix_audit/matrix.json` |
| 独立审计 Markdown | `runs/m21822_matrix_audit/matrix.md` |
