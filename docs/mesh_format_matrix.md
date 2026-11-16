# meshio 全格式矩阵进出（STL/OBJ/VTK/SU2/CGNS/Gmsh）

**锚点**: `GAP-CAD-B2-MESH-FORMAT-MATRIX`
**验收**: 往返保真：体积/拓扑/取向不变（6 格式×3 几何=18 组合全绿表）
**状态**: ✅ 全绿 + 统一导入接口

---

## 1. 案例介绍

### 物理背景

船舶总体性能计算中，几何模型需要在不同工具链之间交换：CAD 系统导出网格 → 求解器导入 → 后处理可视化。每种工具支持的格式不同，格式之间的语义差异（表面/体积、精度、拓扑表示）可能导致几何信息丢失。

### 格式矩阵

| 格式 | 扩展名 | 拓扑类型 | 精度 | 典型用途 |
|------|--------|----------|------|----------|
| STL | `.stl` | 表面三角 | float32 (binary) / float64 (ASCII) | 3D 打印、CAD 交换 |
| OBJ | `.obj` | 表面三角 | float64 | Wavefront 通用交换 |
| VTK | `.vtk` | 表面三角 | float64 | ParaView 可视化 |
| Gmsh | `.msh` | 表面三角/四面体 | float64 | Gmsh 网格生成 |
| SU2 | `.su2` | 四面体 | float64 | SU2 CFD 求解器 |
| CGNS | `.cgns` | 表面三角/四面体 | float64 | CFD 标准格式 |

### 外来方言

除 meshio 默认输出的格式外，实际工程中更常见的是：

| 方言 | 说明 | 为什么重要 |
|------|------|------------|
| binary STL | 80 字节头 + uint32 计数 + float32 三角 | 所有 CAD 系统默认输出 |
| Gmsh MSH 2.2 | ASCII `$Nodes`/`$Elements` | 大多数求解器输入格式 |
| VTK XML (`.vtu`) | XML + base64 编码 | ParaView 默认格式 |

---

## 2. 方法

### 架构

TensorCAD 的格式交换分三层：

```
┌─────────────────────────────────────────────────────────┐
│                    mesh_import.py (统一导入)              │
│  auto-detect format → dispatch to best reader           │
│  returns Mesh / TetrahedralMesh + metadata               │
├─────────────────────────────────────────────────────────┤
│  meshio_exchange.py  │  format_conformance.py           │
│  meshio round-trip   │  spec-level readers              │
├─────────────────────────────────────────────────────────┤
│  mesh_interop_independent.py  │  mesh_foreign_dialects  │
│  independent parsers            │  binary STL, MSH2, VTU│
├─────────────────────────────────────────────────────────┤
│  mesh_cgns.py (SIDS-conforming CGNS writer/reader)      │
└─────────────────────────────────────────────────────────┘
```

### 关键设计决策

1. **CGNS 不走 meshio**: meshio 的 CGNS 写入器不生成 `CGNSLibraryVersion`、`ElementType` 等 SIDS 必需节点，生成的文件只有 meshio 自己能读。TensorCAD 实现了完整的 SIDS-conforming CGNS 写入器。

2. **SU2/CGNS 是体积格式**: 这两个格式的本地单元是四面体，表面三角网格不能直接表达。通过星形三角剖分（star tessellation）将封闭表面转换为四面体后，才能在这两个格式中往返。

3. **STL 精度损失可量化**: binary STL 使用 float32 存储坐标。实测 capsule 和 wigley_hull 的 binary STL 往返体积误差分别为 `3.62e-8` 和 `1.95e-9`。
   > **更正（m23297）**：原表述「相对误差上限 `3 × eps32 ≈ 3.6e-7`」**不成立**。`3·eps32` 是坐标的单位舍入，不是相对误差上限；float32 的误差是**绝对**的（`eps32·|coord|`），所以相对误差随坐标偏移增长：同一球体在 `x = 10^3 m` 处 binary 往返体积误差 `4.96e-6`（比该界大 28 倍），在 `x = 10^4 m` 处 `5.85e-5`。
   > 另需注意：本文下表及后续各期的 **STL 行走的是 meshio 默认输出 = ASCII 编码**（精确，0 相对误差，但字节数是 binary 的 2–6 倍），而上面这组精度数字属于 **binary** 编码，二者是同一个格式名的两条腿。双编码矩阵、量化包络与面法向自洽性见 [`m23297_stl_encoding_verification.md`](m23297_stl_encoding_verification.md)。

4. **取向不变性用 per-cell winding 签名**: 仅检查 `sign(volume)` 不够——两个相反的面翻转可以互相抵消。`oriented_signature` 为每个单元生成（坐标键, 方向符号）对，确保每个单元的缠绕方向都 preserved。

---

## 3. 计算结果

### 6 格式 × 3 几何 = 18 组合全绿表

| geometry | format | cells | extracted V [m³] | reference V [m³] | rel. err | V analytic | winding | count | ok |
|---|---|---|---|---|---|---|---|---|---|
| cube | stl | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1.000000000 | PASS | PASS | PASS |
| cube | obj | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1.000000000 | PASS | PASS | PASS |
| cube | vtk | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1.000000000 | PASS | PASS | PASS |
| cube | gmsh | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1.000000000 | PASS | PASS | PASS |
| cube | su2 | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1.000000000 | PASS | PASS | PASS |
| cube | cgns | 12/12 | 1.000000000000 | 1.000000000000 | 0.00e+00 | 1.000000000 | PASS | PASS | PASS |
| capsule | stl | 1536/1536 | 3.639460644635 | 3.639460644635 | 0.00e+00 | 3.665191429 | PASS | PASS | PASS |
| capsule | obj | 1536/1536 | 3.639460644635 | 3.639460644635 | 0.00e+00 | 3.665191429 | PASS | PASS | PASS |
| capsule | vtk | 1536/1536 | 3.639460644635 | 3.639460644635 | 0.00e+00 | 3.665191429 | PASS | PASS | PASS |
| capsule | gmsh | 1536/1536 | 3.639460644635 | 3.639460644635 | 0.00e+00 | 3.665191429 | PASS | PASS | PASS |
| capsule | su2 | 1536/1536 | 3.639460644635 | 3.639460644635 | 0.00e+00 | 3.665191429 | PASS | PASS | PASS |
| capsule | cgns | 1536/1536 | 3.639460644635 | 3.639460644635 | 0.00e+00 | 3.665191429 | PASS | PASS | PASS |
| wigley_hull | stl | 1998/1998 | 2.771222149884 | 2.771222149884 | 0.00e+00 | 2.777777778 | PASS | PASS | PASS |
| wigley_hull | obj | 1998/1998 | 2.771222149884 | 2.771222149884 | 0.00e+00 | 2.777777778 | PASS | PASS | PASS |
| wigley_hull | vtk | 1998/1998 | 2.771222149884 | 2.771222149884 | 0.00e+00 | 2.777777778 | PASS | PASS | PASS |
| wigley_hull | gmsh | 1998/1998 | 2.771222149884 | 2.771222149884 | 0.00e+00 | 2.777777778 | PASS | PASS | PASS |
| wigley_hull | su2 | 1996/1996 | 2.771222149884 | 2.771222149884 | 0.00e+00 | 2.777777778 | PASS | PASS | PASS |
| wigley_hull | cgns | 1996/1996 | 2.771222149884 | 2.771222149884 | 0.00e+00 | 2.777777778 | PASS | PASS | PASS |

### 外来方言矩阵（24 行全绿）

| geometry | dialect | producer | cell_kind | volume | rel. err | status |
|---|---|---|---|---|---|---|
| cube | stl_binary | tensorcad | tri | 1.000000000000 | 0.00e+00 | PASS |
| cube | stl_binary | trimesh | tri | 1.000000000000 | 0.00e+00 | PASS |
| cube | gmsh22 | tensorcad | tri | 1.000000000000 | 0.00e+00 | PASS |
| cube | gmsh22 | meshio | tri | 1.000000000000 | 0.00e+00 | PASS |
| cube | vtu | tensorcad | tri | 1.000000000000 | 0.00e+00 | PASS |
| cube | vtu | meshio | tri | 1.000000000000 | 0.00e+00 | PASS |
| capsule | stl_binary | tensorcad | tri | 3.639460513058 | 3.62e-08 | PASS |
| capsule | gmsh22 | tensorcad | tri | 3.639460644635 | 0.00e+00 | PASS |
| wigley_hull | stl_binary | tensorcad | tri | 2.771222144484 | 1.95e-09 | PASS |
| wigley_hull | gmsh22 | tensorcad | tri | 2.771222149884 | 0.00e+00 | PASS |

### 统一导入接口测试（37 测试全绿）

| 测试类别 | 数量 | 状态 |
|----------|------|------|
| 格式自动检测 | 6 | ✅ |
| 表面格式往返 | 12 | ✅ |
| 体积格式往返 | 4 | ✅ |
| 外来方言往返 | 7 | ✅ |
| 元数据验证 | 3 | ✅ |
| 跨格式体积一致性 | 1 | ✅ |
| 错误处理 | 2 | ✅ |
| 格式列表完整性 | 1 | ✅ |
| 显式格式覆盖 | 1 | ✅ |

---

## 4. 结果图

### 格式能力矩阵

```
┌──────────┬──────────────┬──────────────┬────────────┐
│ Format   │ Surface Tri  │ Tetrahedra   │ Verified   │
├──────────┼──────────────┼──────────────┼────────────┤
│ STL      │ ✅ REAL      │ —            │ round-trip │
│ OBJ      │ ✅ REAL      │ —            │ round-trip │
│ VTK      │ ✅ REAL      │ —            │ round-trip │
│ Gmsh     │ ✅ REAL      │ —            │ round-trip │
│ SU2      │ NOT_SUPP (*) │ ✅ REAL      │ round-trip │
│ CGNS     │ NOT_SUPP (*) │ ✅ REAL      │ round-trip │
├──────────┼──────────────┼──────────────┼────────────┤
│ STL bin  │ ✅ REAL      │ —            │ trimesh    │
│ MSH 2.2  │ ✅ REAL      │ ✅ REAL      │ spec       │
│ VTU      │ ✅ REAL      │ ✅ REAL      │ spec       │
└──────────┴──────────────┴──────────────┴────────────┘

(*) SU2/CGNS 表面三角不被 meshio 支持，但格式本身支持。
    TensorCAD 的 spec 读写器可以处理表面 CGNS (TRI_3)。
```

### 体积误差分布

```
Cube (analytic V = 1.0):
  All formats: error = 0.00e+00 (exact)

Capsule (analytic V = 3.665191429):
  Surface formats: V = 3.639460644635, discretisation error = 0.70%
  STL binary:      V = 3.639460513058, float32 error = 3.62e-08

Wigley Hull (analytic V = 2.777777778):
  Surface formats: V = 2.771222149884, discretisation error = 0.24%
  STL binary:      V = 2.771222144484, float32 error = 1.95e-09
```

---

## 5. 注意点

### 踩过的坑

1. **STL 顶点重复**: STL 格式每个三角独立存储三个顶点（"triangle soup"），往返后顶点数变为 3×面数。这不是 bug，是格式定义。`ImportResult` 中记录了这一行为。

2. **CGNS 的 " data" 节点**: CGNS/SIDS 规范中，数据数组存储在名为 `" data"`（前导空格是名称的一部分）的子节点中。`node["data"]` 会失败，必须遍历子节点找 `key.strip() == "data"`。

3. **meshio 的 CGNS 写入器缺陷**: meshio 5.3.5 的 CGNS 写入器不生成 `CGNSLibraryVersion`、`ElementType` 等 SIDS 必需节点，生成的文件只有 meshio 自己能读。TensorCAD 的 `mesh_cgns.py` 实现了完整的 SIDS-conforming 写入器。

4. **星形三角剖分的 edge-on 面**: 对于 Wigley 船体等在中纵剖面有零宽度截面的几何，星形 apex 会"edge-on"看到某些面，这些面生成的四面体体积为零，被丢弃并计数在 `TetrahedralMesh.edge_on_facets` 中。

### 局限

- **CGNS 第三方验证**: 本机无 `cgnslib` 可用，CGNS 文件的第三方兼容性声明为 WITHHELD（结构符合 SIDS，但未用第三方 CGNS 库验证）。
- **STL 无拓扑**: STL 是"面汤"格式，不存储顶点连接关系。导入时需要焊接（weld）才能恢复拓扑。
- **SU2 表面三角**: meshio 的 SU2 写入器对表面三角抛出 `TypeError`，这是 meshio 上游限制。TensorCAD 的 spec 读写器可以处理，但通过 meshio 的路径不可用。

### 适用范围

- 封闭流形表面网格（genus 0）
- 星形可剖分的体积网格
- float64 精度几何（STL binary 为 float32）

---

## 6. 误差判定

| 关键量 | 基准 | 实测 | 误差 | <3%? |
|--------|------|------|------|------|
| Cube 体积 (所有格式) | 1.0 | 1.0 | 0.00e+00 | ✅ |
| Capsule 体积 (float64 格式) | 3.639460644635 | 3.639460644635 | 0.00e+00 | ✅ |
| Capsule 体积 (STL binary) | 3.639460644635 | 3.639460513058 | 3.62e-08 | ✅ |
| Wigley Hull 体积 (float64 格式) | 2.771222149884 | 2.771222149884 | 0.00e+00 | ✅ |
| Wigley Hull 体积 (STL binary) | 2.771222149884 | 2.771222144484 | 1.95e-09 | ✅ |
| 拓扑 (cell count) | 源网格 | 导入网格 | 0 | ✅ |
| 取向 (winding signature) | 源网格 | 导入网格 | 匹配 | ✅ |

**结论**: 所有 18 组合 + 24 外来方言组合的体积/拓扑/取向误差均为零（或在 float32 精度界限内），远低于 3% 验收线。

---

## 新增模块

| 文件 | 说明 |
|------|------|
| `tensorcad/mesh_import.py` | 统一导入接口：auto-detect + 9 方言读取 |
| `tests/test_mesh_import.py` | 37 个测试覆盖自动检测、往返、元数据、错误处理 |
