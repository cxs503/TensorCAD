# TensorCAD meshio 格式矩阵：进口方向（四边形/多边形/多块）

**Mission**: m20888 — anchor `GAP-CAD-B2-MESH-FORMAT-MATRIX`
**Date**: 2026-11-15 · **Repo**: TensorCAD
**Evidence**: `runs/m20888_meshio_import_matrix/matrix.json`（本文由 `scripts/gen_m20888_doc.py` 从该 JSON 生成，数字零手抄）

---

## 1. 案例介绍（物理背景/几何/参数）

### 物理背景

格式矩阵此前只验证了**往返（round trip）**：TensorCAD 自己写、自己读。这只能证明自洽，不能证明**进口**——真实文件不是 TensorCAD 写的。船体几何从 CAD 内核 / Salome / Pointwise / 面元法结构化网格工具出来时，面片绝大多数是**四边形**（结构化船体面板），并常常**分成多个 cell block**（Gmsh region、CGNS zone、边界 patch 与体单元同文件）。

对抗性复核发现进口方向有两个实质缺陷（复现见第 5 节）：

1. **D1（阻断级）**：四边形/多边形封闭体一律被拒。最简单的封闭体——6 个四边形的立方体——在 OBJ / VTK / Gmsh 三个格式上都抛 `MeshioFormatError: meshio mesh contains no triangle cells`。
2. **D2（静默丢几何）**：`meshio_to_mesh` 只取**第一个** triangle block，其余块被静默丢弃；`meshio_to_tetrahedral` 同样只取第一个 tetra block。一个「四边形壳体 + 1 个三角形块」的文件只回来 1 个三角形，体积完全错误，且不报错。

### 测试几何（三个，均独立于 TensorCAD 的三角化 MeshBuilder 构造）

| 几何 | 描述 | 输入单元块 | 解析体积 [m³] |
|---|---|---|---|
| `cube_quads` | unit cube, 6 quads | quad | 1.000000000000 |
| `prism_ngons` | pentagonal prism, 2 n-gon caps + 5 quads | polygon, quad | 2.377641290738 |
| `wigley_quads` | Wigley hull, all-quad structured cage | quad | 2.777777777778 |

Wigley 参考体积为解析式 `V = L B T n m / ((n+1)(m+1)) = 4LBT/9`（`tensorcad.wigley.analytic_volume`，cube_quads 之外的几何用各自闭式解）。

### 格式与拓扑

| 格式 | 能承载四边形 | 能承载 n-gon | 进矩阵 |
|---|---|---|---|
| OBJ | ✅ | ✅ | ✅ 3 几何 |
| VTK（legacy ASCII） | ✅ | ✅ | ✅ 3 几何 |
| Gmsh | ✅（MSH 4.1 由 meshio 写出） | ❌ meshio 写端无 polygon 单元 | 2 几何（n-gon 记为格式不支持） |
| STL | ❌（三角形面片汤，按定义） | ❌ | 不进进口矩阵 |
| SU2 / CGNS | —（体积格式，见往返矩阵） | — | 不在本表范围 |

## 2. 方法（流程/关键设置/工具链）

- 几何与连接性**在脚本里独立构造**（四边形/多边形 cage，节点按精确位置焊接，与真实写出器一致），不经过 TensorCAD 的三角化 builder，因此是「外来文件」。
- 用 `meshio.write` 落盘（模拟外部工具产出）→ `tensorcad.meshio_exchange.read_meshio` 读入 → 对每条记录检查：体积（散度定理）、三角面数、水密、边流形、零面积面、取向符号、文件 sha256。
- 判据：**进口体积 vs 导出前 cage 体积**相对误差 < 1e-12（无损格式应达机器精度）；**进口体积 vs 解析值**反映离散误差，不是格式误差。
- 工具链：Python 3.11.6 · meshio 5.3.5 · numpy 2.4.6

## 3. 计算结果

共 9 个组合：8 绿 + 1 格式级不支持。全绿：**True**。

| 几何 | 格式 | 进口体积 [m³] | vs 导出前 cage | vs 解析 | 面数 | 水密 | 取向 |
|---|---|---|---|---|---|---|---|
| `cube_quads` | obj | 1.000000000000 | 0.00e+00 | 0.000e+00 | 12/12 | ✅ | ✅+ |
| `cube_quads` | vtk | 1.000000000000 | 0.00e+00 | 0.000e+00 | 12/12 | ✅ | ✅+ |
| `cube_quads` | gmsh | 1.000000000000 | 0.00e+00 | 0.000e+00 | 12/12 | ✅ | ✅+ |
| `prism_ngons` | obj | 2.377641290738 | 0.00e+00 | 1.868e-16 | 16/16 | ✅ | ✅+ |
| `prism_ngons` | vtk | 2.377641290738 | 0.00e+00 | 1.868e-16 | 16/16 | ✅ | ✅+ |
| `wigley_quads` | obj | 2.771222149884 | 8.01e-16 | 2.360e-03 | 1998/1998 | ✅ | ✅+ |
| `wigley_quads` | vtk | 2.771222149884 | 8.01e-16 | 2.360e-03 | 1998/1998 | ✅ | ✅+ |
| `wigley_quads` | gmsh | 2.771222149884 | 8.01e-16 | 2.360e-03 | 1998/1998 | ✅ | ✅+ |

- `prism_ngons` × gmsh：**跳过** — the meshio Gmsh writer carries no polygon cell type (MSH 2.2 elements 1-5 are tri/quad/tet/hex/prism; n-gons need MSH 4.1 entity tags that only appear when the file has one cell type) -- a format limitation, not an importer one

最大进口保真误差（vs 导出前 cage）：**8.01e-16**；vs 解析解最大 2.360e-03（Wigley 网格离散，与 m18268 的 0.24 % 一致）。

图 1：进口保真误差（左：vs 导出前 cage；右：vs 解析值）

![import fidelity](../artifacts/m20888_meshio_import_matrix/import_fidelity.png)

图 2：Wigley 全四边形船体面板 cage（面元法工具典型产出，进口后被三角化为 2×quads）

![wigley quad cage](../artifacts/m20888_meshio_import_matrix/wigley_quad_cage.png)

图 3：多块文件静默丢几何的修复前后对照

![mixed block loss](../artifacts/m20888_meshio_import_matrix/mixed_block_loss.png)

## 4. 注意点（坑/局限/适用范围）

1. **四边形拆分方向固定**：四边形 `(a,b,c,d)` 沿 `0-2` 对角线拆成 `(a,b,c)+(a,c,d)`；多边形从 0 号顶点扇形展开。两者都保持父面片的环向（winding），因此外法向/正体积在进口后不变（测试 `test_quad_winding_preserved` 逐面验证）。
2. **凹多边形与不共面 n-gon 不在保证范围内**：扇形展开对凸且近似共面的 n-gon 正确；凹 n-gon 会产生自交三角形。CAD 网格通常已是凸四边形/三角形，但若来源是凹多边形，需先分块三角化（本模块不做）。
3. **`meshio_to_mesh` 会跳过 3D 体单元**，返回体网格的**边界面**；要体单元本身请用 `meshio_to_tetrahedral`（错误消息会明确提示）。`line`/`vertex` 块（边界边集、孤立点）不是几何，跳过。
4. **六面体/棱柱/金字塔不自动拆分**：`meshio_to_tetrahedral` 遇到它们会报错并在消息里点名，而不是静默丢弃——拆分是有歧义的（Kuhn 6 分或 5 分），应在 `tensorcad.volume_mesh` 侧决定，不在 IO 层猜。
5. **重复索引三角形被丢弃**：退化面片（节点重合）不构成面积，丢弃后不改变体积；非退化零面积面（三点共线）仍会被保留，交由几何门 `zero_area_faces` 报出。
6. **STL 不在本矩阵**：STL 是三角形面片汤，结构上无法承载四边形；其往返由 `tests/test_meshio_format_matrix.py` 覆盖。
7. **Gmsh n-gon 写端不支持**：meshio 的 Gmsh 2.2 无 polygon 单元，4.1 写多单元类型需要 `gmsh:geometrical` 实体标签。这是**写端/格式**限制，不是进口端限制，故棱柱×Gmsh 记为跳过。

## 5. 缺陷复现与修复（对抗性复核结论）

| 缺陷 | 复现 | 修复 | 门禁测试 |
|---|---|---|---|
| D1 四边形封闭体被拒 | 6-quad 立方体写 OBJ/VTK/Gmsh 后 `read_meshio` → `MeshioFormatError: contains no triangle cells` | `_triangulate_block()`：triangle 直通、quad 对角拆、polygon 扇形展开，消费**所有**表面块 | `tests/test_meshio_import_cells.py::test_quad_cube_imports_from_file`（3 格式参数化） |
| D2 多块静默丢几何 | quad 壳体 + 1 triangle 块 → 只回来 1 个面；双 tetra zone → 只回来 1 个 tet | `meshio_to_mesh` 拼接所有表面块；`meshio_to_tetrahedral` 拼接所有 tetra 块 | `test_mixed_blocks_all_consumed` / `test_two_triangle_blocks_concatenated` / `test_two_tetra_blocks_concatenated` |

修复后同一多块文件名副其实：输入 13 个面 → 进口 13 个面（修复前 1 个）。

## 6. 误差判定

| 关键量 | 基准 | 误差 | 是否 < 3 % |
|---|---|---|---|
| 进口体积（8 组合） | 导出前 cage 体积（散度定理，同一几何） | 最大 8.01e-16 | ✅ 是（远优于 3 %，达机器精度） |
| 进口体积 Wigley 全四边形 cage | 解析 `4LBT/9` | 2.360e-03 | ✅ 是（离散误差，非格式误差） |
| 进口体积 立方体 / 五棱柱 | 解析 1.0 / 2.5R²sin72° | 0 / 1.9e-16 | ✅ 是 |
| 拓扑（面数） | 独立参考三角化面数 | 全组合相等 | ✅ 是 |
| 水密 / 边流形 / 零面积面 / 取向 | 几何门判据 | 全组合通过 | ✅ 是 |

结论：进口方向的体积、拓扑、取向三项在 8 个可行组合上全部保持到机器精度；唯一扣分项是 Gmsh×n-gon 为格式级不支持（已如实记录，不粉饰）。

