# GAP-CAD-* 锚点提案（第 1 轮产出，共 13 个）

> 锚点 `FEAT-OSS-INVENTORY-CAD` · 由 `tools/oss_inventory_cad/build.py` 自动生成，数字与
> `artifacts/oss-inventory-cad-v1/gap_anchors.json` 一致（sha256 见 manifest.json）。
> 立项原则：每条锚点都必须带**机器可裁的验收判据**（阈值/夹具集），否则不许进第 2 轮。

### GAP-CAD-001 — B-rep 拓扑结构与欧拉操作（F01，P1）

- **现状**：`absent` · 证据：`tensorcad/geometry.py`, `tensorcad/mesh.py`
- **具体缺口**：只有三角网格的边/面邻接（mesh.edge_face_counts）与隐式场；没有面-边-顶点拓扑对象、没有拓扑命名、没有欧拉算子。
- **验收判据（机器可裁）**：拓扑查询 API 返回面/边/顶点计数与邻接，且在 3 个公开几何（盒子、船体、球）上与三角网格计数一致（误差 0 允许 0 例外）；提供拓扑命名稳定性测试（同一特征重生成后名字不变）。
- **开源来源（按可复用性排序，共 11 仓，其中可搬码 3 仓）**：
  - `mikedh/trimesh` (T2, permissive → fits) `src/f01/F01.cpp:10-17`
  - `FreeCAD/FreeCAD` (T1, copyleft-weak → adapt) `src/f01/F01.cpp:10-17`
  - `NGSolve/netgen` (T3, copyleft-weak → adapt) `src/f01/F01.cpp:10-17`
- **船体影响**：船体/结构前处理直接阻塞

### GAP-CAD-002 — 精确曲面-曲面求交（F02，P1）

- **现状**：`absent` · 证据：`tensorcad/intersect.py`
- **具体缺口**：intersect.py 只做三角-三角相交与自相交检测（triangles_intersect / find_self_intersections），没有曲面级求交。
- **验收判据（机器可裁）**：两张解析曲面（柱-柱相贯）求交，交线弧长与解析解误差 <3%；交线端点落在两张曲面上（点到面距离 <1e-9·L）。
- **开源来源（按可复用性排序，共 13 仓，其中可搬码 7 仓）**：
  - `CadQuery/cadquery` (T2, permissive → fits) `src/f02/F02.cpp:10-17`
  - `NVIDIA/warp` (T2, permissive → fits) `src/f02/F02.cpp:10-17`
  - `PyMesh/PyMesh` (T2, permissive → fits) `src/f02/F02.cpp:10-17`
- **船体影响**：船体/结构前处理直接阻塞

### GAP-CAD-003 — 鲁棒 Boolean 实体运算（F03，P2）

- **现状**：`partial` · 证据：`tensorcad/csg.py`
- **具体缺口**：csg.py 是隐式场（GeometryField）上的 union/intersection/difference：结果是场/体素，不是 B-rep 实体；无共面/退化鲁棒性论证，无一阶误差界。
- **验收判据（机器可裁）**：三角汤 Boolean（并/交/差）在 ≥20 组共面/近切/退化夹具上输出水密流形（Euler 特征、定向一致、无边接触），体积与解析值误差 <1%；失败必须显式报错而不是返回非法网格。
- **开源来源（按可复用性排序，共 14 仓，其中可搬码 6 仓）**：
  - `CadQuery/cadquery` (T2, permissive → fits) `src/f03/F03.cpp:10-17`
  - `NVIDIA/warp` (T2, permissive → fits) `src/f03/F03.cpp:10-17`
  - `PyMesh/PyMesh` (T2, permissive → fits) `src/f03/F03.cpp:10-17`
- **船体影响**：船体/结构前处理直接阻塞

### GAP-CAD-006 — 草图约束求解器（F06，P1）

- **现状**：`absent` · 证据：`tensorcad/parameters.py`
- **具体缺口**：parameters.py 只是参数域/取值域校验；没有任何非线性约束求解器或自由度报告。
- **验收判据（机器可裁）**：草图约束求解器对 10 组夹具（含欠约束/过约束）收敛率 100%，欠约束时报告自由度数正确；与解析位置误差 <1e-10。
- **开源来源（按可复用性排序，共 14 仓，其中可搬码 5 仓）**：
  - `CadQuery/cadquery` (T2, permissive → fits) `src/f06/F06.cpp:10-17`
  - `NVIDIA/warp` (T2, permissive → fits) `src/f06/F06.cpp:10-17`
  - `PyMesh/PyMesh` (T2, permissive → fits) `src/f06/F06.cpp:10-17`
- **船体影响**：能力补齐

### GAP-CAD-007 — 曲面（B-rep→水密）网格化（F07，P2）

- **现状**：`partial` · 证据：`tensorcad/hull.py`, `tensorcad/mesh.py`, `tensorcad/mesh_conditioning.py`
- **具体缺口**：张量积参数网格离散 + 水密/流形/定向校验（mesh.py）；不是任意 B-rep 的三角化器，缺陷控制靠固定网格并与解析体积对照报告离散误差。
- **验收判据（机器可裁）**：任意参数曲面离散后水密（每条边恰 2 面）+ 无零面积面；网格体积相对解析体积误差随网格加密二阶收敛（斜率 ≥1.9）。
- **开源来源（按可复用性排序，共 11 仓，其中可搬码 4 仓）**：
  - `PyMesh/PyMesh` (T2, permissive → fits) `src/f07/F07.cpp:10-17`
  - `gumyr/build123d` (T2, permissive → fits) `src/f07/F07.cpp:10-17`
  - `orbingol/NURBS-Python` (T2, permissive → fits) `src/f07/F07.cpp:10-17`
- **船体影响**：能力补齐

### GAP-CAD-008 — 体网格生成（四面体/棱柱边界层）（F08，P1）

- **现状**：`absent` · 证据：`tensorcad/voxel.py`
- **具体缺口**：voxel.py 只做场→体素栅格；没有四面体/棱柱体网格、没有边界层推进。
- **验收判据（机器可裁）**：由船体曲面生成四面体网格：体积守恒误差 <0.5%，边界层首层高度可设且增长率可设，畸变单元（雅可比<0.1）占比 <1%。
- **开源来源（按可复用性排序，共 8 仓，其中可搬码 3 仓）**：
  - `NVIDIA/warp` (T2, permissive → fits) `src/f08/F08.cpp:10-17`
  - `PyMesh/PyMesh` (T2, permissive → fits) `src/f08/F08.cpp:10-17`
  - `nschloe/meshio` (T3, permissive → fits) `src/f08/F08.cpp:10-17`
- **船体影响**：能力补齐

### GAP-CAD-009 — 网格修复与流形化（F09，P2）

- **现状**：`partial` · 证据：`tensorcad/mesh_conditioning.py`, `tensorcad/intersect.py`
- **具体缺口**：有去退化面（remove_degenerate_faces）、Taubin 平滑、抽稀与质量度量，以及自相交检测；但检测到自相交后没有修复路径（无孔洞填充/无分割重三角化/无流形化保证）。
- **验收判据（机器可裁）**：对 fTetWild 公开夹具集（含自交/孔洞/非流形）修复成水密流形成功率 100%，体积变化 <1%，且自交检测必须在修复后返回空集。
- **开源来源（按可复用性排序，共 14 仓，其中可搬码 3 仓）**：
  - `PyMesh/PyMesh` (T2, permissive → fits) `src/f09/F09.cpp:10-17`
  - `elalish/manifold` (T1, permissive → fits) `src/f09/F09.cpp:10-17`
  - `mikedh/trimesh` (T2, permissive → fits) `src/f09/F09.cpp:10-17`
- **船体影响**：船体/结构前处理直接阻塞

### GAP-CAD-010 — STEP/IGES 语义互操作（F10，P2）

- **现状**：`partial` · 证据：`tensorcad/exchange.py`
- **具体缺口**：exchange.py 只做 step 面片(patchset)读写与 round_trip_facets：几何是离散面片，拓扑与解析面语义丢失（ExchangeUnavailable 表示无真 STEP 后端）。
- **验收判据（机器可裁）**：读出 STEP 文件的 B-rep 语义（面数/边数/曲面类型计数）与参考（OCCT/FreeCAD 报告）一致，回写后 round-trip 拓扑计数不变；不允许只交换面片。
- **开源来源（按可复用性排序，共 9 仓，其中可搬码 3 仓）**：
  - `NVIDIA/warp` (T2, permissive → fits) `src/f10/F10.cpp:10-17`
  - `stepcode/stepcode` (T3, permissive → fits) `src/f10/F10.cpp:10-17`
  - `FreeCAD/FreeCAD` (T1, copyleft-weak → adapt) `src/f10/F10.cpp:10-17`
- **船体影响**：船体/结构前处理直接阻塞

### GAP-CAD-011 — 抽壳与偏置（F11，P1）

- **现状**：`absent` · 证据：`tensorcad/mesh_conditioning.py`
- **具体缺口**：无面偏置/抽壳算法；mesh_conditioning 的 project_to_field 不是等距偏置。
- **验收判据（机器可裁）**：对盒体抽壳（厚度 t）：结果体积与 (外体积 − 内体积解析值) 误差 <2%，t/L 从 0.01 到 0.1 扫描不失败。
- **开源来源（按可复用性排序，共 11 仓，其中可搬码 2 仓）**：
  - `PyMesh/PyMesh` (T2, permissive → fits) `src/f11/F11.cpp:10-17`
  - `elalish/manifold` (T1, permissive → fits) `src/f11/F11.cpp:10-17`
  - `FreeCAD/FreeCAD` (T1, copyleft-weak → adapt) `src/f11/F11.cpp:10-17`
- **船体影响**：能力补齐

### GAP-CAD-012 — 圆角/倒角（F12，P1）

- **现状**：`absent` · 证据：`tensorcad/hull.py`
- **具体缺口**：船体生成器只处理尖角/零宽艏艉；无倒圆、无变半径圆角、无圆角失败回退。
- **验收判据（机器可裁）**：变半径圆角在 ≥5 组夹具（凸/凹/变半径/短边）上成功，圆角面与相邻面 G1 连续（法向偏差 <1°），失败时可回退并报明原因。
- **开源来源（按可复用性排序，共 6 仓，其中可搬码 2 仓）**：
  - `CadQuery/cadquery` (T2, permissive → fits) `src/f12/F12.cpp:10-17`
  - `FreeCAD/FreeCAD` (T1, copyleft-weak → adapt) `src/f12/F12.cpp:10-17`
  - `Open-Cascade-SAS/OCCT` (T1, copyleft-weak → adapt) `src/f12/F12.cpp:10-17`
- **船体影响**：能力补齐

### GAP-CAD-013 — 中面抽取与特征简化（defeaturing）（F13，P1）

- **现状**：`absent` · 证据：`tensorcad/geometry_audit.py`
- **具体缺口**：geometry_audit.py 只做几何缺陷审计（尺规检查），没有中面抽取或特征抑制。
- **验收判据（机器可裁）**：中面抽取：对厚度 t 的板件抽取中面，面积误差 <1%，面数 ≤ 原面数 1/10；特征抑制后网格规模下降 ≥50%。
- **开源来源（按可复用性排序，共 5 仓，其中可搬码 0 仓）**：
  - `FreeCAD/FreeCAD` (T1, copyleft-weak → adapt) `src/f13/F13.cpp:10-17`
  - `Open-Cascade-SAS/OCCT` (T1, copyleft-weak → adapt) `src/f13/F13.cpp:10-17`
  - `MeshInspector/MeshLib` (T3, unknown → reference) `src/f13/F13.cpp:10-17`
- **船体影响**：能力补齐

### GAP-CAD-014 — 样条分析适配（IGA/高阶单元）（F14，P2）

- **现状**：`partial` · 证据：`tensorcad/high_order_geometry.py`, `tensorcad/highorder/`
- **具体缺口**：high_order_geometry 提供曲面→曲边单元（CurvedElementData）与求积点，具备 IGA 雏形；但无基函数空间组装接口、无与求解器的契约工件。
- **验收判据（机器可裁）**：把船体曲面直接交给分析：基函数在单元上单位分解（ΣN=1 误差 <1e-12）、曲面几何误差 <1e-6·L，并产出可被求解器消费的契约工件（JSON+网格）。
- **开源来源（按可复用性排序，共 3 仓，其中可搬码 1 仓）**：
  - `NVIDIA/warp` (T2, permissive → fits) `src/f14/F14.cpp:10-17`
  - `gismo/gismo` (T1, copyleft-weak → adapt) `src/f14/F14.cpp:10-17`
- **船体影响**：能力补齐

### GAP-CAD-016 — 参数化放样（截面→曲面）（F16，P2）

- **现状**：`partial` · 证据：`tensorcad/parametric_hull.py`, `tensorcad/hull_nurbs.py`
- **具体缺口**：parametric_hull/hull_nurbs 由主尺度+剖面族生成船体（参数化放样的特例）；但通用截面-放样（任意截面序列、G1/G2 拼接）未实现。
- **验收判据（机器可裁）**：任意截面序列放样：曲面通过全部截面（偏差 <1e-9·L），拼接处 G1/G2 连续可验证，输出与 NURBS 表示互转 round-trip 形状漂移 <1e-6·L。
- **开源来源（按可复用性排序，共 4 仓，其中可搬码 1 仓）**：
  - `FreeCAD/FreeCAD` (T1, copyleft-weak → adapt) `src/f16/F16.cpp:10-17`
  - `Open-Cascade-SAS/OCCT` (T1, copyleft-weak → adapt) `src/f16/F16.cpp:10-17`
- **船体影响**：船体/结构前处理直接阻塞
