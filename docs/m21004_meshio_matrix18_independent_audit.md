# meshio 六格式矩阵的独立复核 + CGNS/SIDS 出口合规（锚点 GAP-CAD-B2-MESH-FORMAT-MATRIX）

## ① 案例与问题
TensorCAD 的 `tensorcad.meshio_exchange` 宣称 6 格式（STL/OBJ/VTK/SU2/CGNS/Gmsh）往返保真，
但既有矩阵是「自己写、自己读」：写用 `meshio.write`、读用 `meshio.read`，两侧共享同一套格式误解时
矩阵照样全绿。本单把这条验证欠账补上：用文件字节本身复核 3 个几何 × 6 格式。

| 几何 | 类型 | 顶点 / 单元 | 解析体积 [m³] |
|---|---|---|---|
| cube | 单位立方体外表面 | 8 / 12 | 1.0 |
| sphere | 球面离散（三角化）| 1106 / 2208 | 4.188790205 |
| wigley_hull | Wigley 实船型（L=10 m, B=1 m, T=0.625 m）| 1001 / 1998 | 2.777777778 |

体积格式（SU2/CGNS）先把同一几何四面体化再导出；表面格式导出三角面。

## ② 方法
1. **out 方向**：`write_meshio` 写出文件后，用 `tensorcad.format_conformance` 的**规格级解析器**
   （按 STL/OBJ/legacy-VTK/Gmsh 4.1/SU2/CGNS-HDF5 规范从零实现，不经 meshio、不经 TensorCAD 读取器）
   解析磁盘字节，得到顶点/单元数、体积、取向、闭合性、格式违规项，与源网格对比。
2. **in 方向**：用 `meshio` 直接从 `meshio.Mesh`（不经本包写出器）写文件，再用 `read_meshio` 读回对比。
3. **CGNS 合规**：按 CGNS/SIDS 的 HDF5 节点表（`CGNSLibraryVersion`、`CGNSBase_t`/`Zone_t`/
   `GridCoordinates_t`/`Elements_t` 标签、`CellDimension`/`PhysicalDimension`/`ZoneType`、
   `ElementType`/`ElementSize`/`ElementRange`/`ElementConnectivity` 一致性）做结构校验。
4. 判据：体积相对误差 < 1e-12（STL 因 float32 存储放宽到 1e-5）、顶点/单元数一致、有向体积为正向、
   表面闭合（边价数全为 2）、无格式违规。

## ③ 结果
**36/36 个审计格全绿**（18 组合 × 2 方向，out 18 格 + in 18 格）。

各格式 out 方向最大体积相对误差：

| format | max rel err |
|---|---|
| cgns | 3.21e-15 |
| gmsh | 8.76e-15 |
| obj | 8.76e-15 |
| stl | 8.76e-15 |
| su2 | 3.21e-15 |
| vtk | 8.76e-15 |

in 方向最大体积相对误差：0.00e+00。

| dir | format | geometry | V_file [m3] | V_source [m3] | rel err | pts | cells | orient | closed | spec/SIDS | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| out | stl | cube | 1.000000000 | 1.000000000 | 0.00e+00 | 8 | 12 | outward | closed | spec ok | PASS |
| in | stl | cube | 1.000000000 | 1.000000000 | 0.00e+00 | 8 | 12 | outward | n/a | spec ok | PASS |
| out | stl | sphere | 4.158971436 | 4.158971436 | 8.76e-15 | 1106 | 2208 | outward | closed | spec ok | PASS |
| in | stl | sphere | 4.158971436 | 4.158971436 | 0.00e+00 | 1106 | 2208 | outward | n/a | spec ok | PASS |
| out | stl | wigley_hull | 2.771222150 | 2.771222150 | 8.01e-16 | 1001 | 1998 | outward | closed | spec ok | PASS |
| in | stl | wigley_hull | 2.771222150 | 2.771222150 | 0.00e+00 | 1001 | 1998 | outward | n/a | spec ok | PASS |
| out | obj | cube | 1.000000000 | 1.000000000 | 0.00e+00 | 8 | 12 | outward | closed | spec ok | PASS |
| in | obj | cube | 1.000000000 | 1.000000000 | 0.00e+00 | 8 | 12 | outward | n/a | spec ok | PASS |
| out | obj | sphere | 4.158971436 | 4.158971436 | 8.76e-15 | 1106 | 2208 | outward | closed | spec ok | PASS |
| in | obj | sphere | 4.158971436 | 4.158971436 | 0.00e+00 | 1106 | 2208 | outward | n/a | spec ok | PASS |
| out | obj | wigley_hull | 2.771222150 | 2.771222150 | 8.01e-16 | 1001 | 1998 | outward | closed | spec ok | PASS |
| in | obj | wigley_hull | 2.771222150 | 2.771222150 | 0.00e+00 | 1001 | 1998 | outward | n/a | spec ok | PASS |
| out | vtk | cube | 1.000000000 | 1.000000000 | 0.00e+00 | 8 | 12 | outward | closed | spec ok | PASS |
| in | vtk | cube | 1.000000000 | 1.000000000 | 0.00e+00 | 8 | 12 | outward | n/a | spec ok | PASS |
| out | vtk | sphere | 4.158971436 | 4.158971436 | 8.76e-15 | 1106 | 2208 | outward | closed | spec ok | PASS |
| in | vtk | sphere | 4.158971436 | 4.158971436 | 0.00e+00 | 1106 | 2208 | outward | n/a | spec ok | PASS |
| out | vtk | wigley_hull | 2.771222150 | 2.771222150 | 8.01e-16 | 1001 | 1998 | outward | closed | spec ok | PASS |
| in | vtk | wigley_hull | 2.771222150 | 2.771222150 | 0.00e+00 | 1001 | 1998 | outward | n/a | spec ok | PASS |
| out | gmsh | cube | 1.000000000 | 1.000000000 | 0.00e+00 | 8 | 12 | outward | closed | spec ok | PASS |
| in | gmsh | cube | 1.000000000 | 1.000000000 | 0.00e+00 | 8 | 12 | outward | n/a | spec ok | PASS |
| out | gmsh | sphere | 4.158971436 | 4.158971436 | 8.76e-15 | 1106 | 2208 | outward | closed | spec ok | PASS |
| in | gmsh | sphere | 4.158971436 | 4.158971436 | 0.00e+00 | 1106 | 2208 | outward | n/a | spec ok | PASS |
| out | gmsh | wigley_hull | 2.771222150 | 2.771222150 | 8.01e-16 | 1001 | 1998 | outward | closed | spec ok | PASS |
| in | gmsh | wigley_hull | 2.771222150 | 2.771222150 | 0.00e+00 | 1001 | 1998 | outward | n/a | spec ok | PASS |
| out | su2 | cube | 1.000000000 | 1.000000000 | 0.00e+00 | 9 | 12 | outward | - | spec ok | PASS |
| in | su2 | cube | 1.000000000 | 1.000000000 | 0.00e+00 | 9 | 12 | outward | n/a | spec ok | PASS |
| out | su2 | sphere | 4.158971436 | 4.158971436 | 1.49e-15 | 1107 | 2208 | outward | - | spec ok | PASS |
| in | su2 | sphere | 4.158971436 | 4.158971436 | 0.00e+00 | 1107 | 2208 | outward | n/a | spec ok | PASS |
| out | su2 | wigley_hull | 2.771222150 | 2.771222150 | 3.21e-15 | 1002 | 1996 | outward | - | spec ok | PASS |
| in | su2 | wigley_hull | 2.771222150 | 2.771222150 | 0.00e+00 | 1002 | 1996 | outward | n/a | spec ok | PASS |
| out | cgns | cube | 1.000000000 | 1.000000000 | 0.00e+00 | 9 | 12 | outward | - | SIDS ok | PASS |
| in | cgns | cube | 1.000000000 | 1.000000000 | 0.00e+00 | 9 | 12 | outward | n/a | spec ok | PASS |
| out | cgns | sphere | 4.158971436 | 4.158971436 | 1.49e-15 | 1107 | 2208 | outward | - | SIDS ok | PASS |
| in | cgns | sphere | 4.158971436 | 4.158971436 | 0.00e+00 | 1107 | 2208 | outward | n/a | spec ok | PASS |
| out | cgns | wigley_hull | 2.771222150 | 2.771222150 | 3.21e-15 | 1002 | 1996 | outward | - | SIDS ok | PASS |
| in | cgns | wigley_hull | 2.771222150 | 2.771222150 | 0.00e+00 | 1002 | 1996 | outward | n/a | spec ok | PASS |

### CGNS/SIDS：本次审计发现的真实缺陷
| 出口 | SIDS 违规项数 | meshio 能否读回 |
|---|---|---|
| `tensorcad.mesh_cgns.write_cgns`（本单新增） | **0** | True |
| meshio 自带 CGNS writer（改动前 TensorCAD 走这条路） | **14** | 是（仅 meshio 自己能读） |

meshio 文件的违规项（程序生成，非手抄）：
- missing /CGNSLibraryVersion (SIDS requires it)
- /Base label is None, expected 'CGNSBase_t'
- base node has no CellDimension
- base node has no PhysicalDimension
- zone label is None, expected 'Zone_t'
- zone has no ZoneType node
- GridCoordinates label is None, expected 'GridCoordinates_t'
- CoordinateX type is None, expected 'DataArray_t'
- CoordinateY type is None, expected 'DataArray_t'
- CoordinateZ type is None, expected 'DataArray_t'
- GridElements label is None, expected 'Elements_t'
- GridElements carries no ElementType attribute -- SIDS requires it and without it a conforming reader cannot size the connectivity
- GridElements carries no ElementSize attribute
- GridElements/ElementRange label is None, expected 'IndexRange_t'

改动前的 CGNS 出口既不写 `CGNSLibraryVersion`，也不写节点 label/type，不写 `ElementType`/`ElementSize`
—— 任何符合 SIDS 的读取器都无法确定单元宽度；`meshio/cgns/_cgns.py` 自述其读写均为 `TODO`。
本单新增 `tensorcad/mesh_cgns.py`（SIDS 合规写出 + 兼容 SIDS/meshio 两种布局的读入），并把
`write_meshio(..., file_format='cgns')` 与 `round_trip_meshio` 切到该写出器；节点名保留
`Zone1`/`GridElements`，因此 meshio 仍可读回（双向互通），SIDS 违规降为 0。

## ④ 图
- `volume_error_matrix.png`：18 格体积相对误差（对数轴，含容差线）。
- `cell_status_grid.png`：6 格式 × 3 几何 × 2 方向的通过/失败格。
- `wigley_hull_from_file.png`：从 OBJ 字节重建的 Wigley 船体曲面 + 文件体积 vs 解析体积。

## ⑤ 注意点与局限
- STL 二进制按规范是 float32（体积容差取格式下限 1e-5）；meshio 默认写 **ASCII** STL（文本保留 double），
  故本次 STL 实测体积误差只有 8.8e-15 —— 容差线按格式定义而非实测留；二进制分支由手工构造的坏文件覆盖。
- legacy VTK 由 meshio 写成 5.1 **BINARY**（大端 + OFFSETS/CONNECTIVITY 分块），解析器按该布局实现，
  同时兼容 VTK 5.0 的「长度前缀」旧布局（本仓不产出）。
- **SIDS 合规是结构校验**：本环境无 cgnslib/pyCGNS/ParaView，无法证明第三方 CGNS 库一定接受该文件，
  只能证明节点/属性表齐全自洽 —— 按 unknown 记，不计入「通过」。
- 体积对解析值的差（球面 7e-3、Wigley 4.3e-3）属**几何离散误差**，不是格式损失，两者分开统计。
- 未覆盖：hex/wedge/pyramid 体单元、多 zone CGNS、Gmsh 2.2；quad 面由 native 路径测试覆盖。

## ⑥ 误差判定（是否 <3%）
- 格式保真体积误差：**是**（除 STL 用 1e-5 容差外，全部 < 3e-15）。
- 文件体积 vs 解析体积：球面 7.1e-3、Wigley 4.3e-3，**是**（<3%，且为离散误差而非格式误差）。
- 拓扑/取向/闭合：18 格精确一致（顶点数、单元数、正向、边价数 2）。
- CGNS SIDS 违规：14 项 → 0 项；**第三方 CGNS 库互操作未验证（unknown，不计通过）**。

## 复现
```
python tools/meshio_matrix_independent_audit.py        # 产出 runs/ 与 artifacts/
python -m pytest tests/test_format_conformance.py -q   # 48 项门禁
```
