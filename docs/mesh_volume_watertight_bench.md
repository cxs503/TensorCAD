# 闭体网格体积/水密性核验（TensorCAD）— CASE-CAD-VOLUME-CLOSED

任务：算例池第 1 轮一例。判据 `<0.1 %`，ref = 解析几何（球/长方体）。本文所有数字由
`benchmarks/verified/mesh_volume_watertight/publish_report.py` 从同目录 `result.json` 程序化注入，
无手抄；一键复现命令见末节。

## ① 问题与 ref 出处

**要判定的问题**：TensorCAD 的闭体三角网格（a）体积是否等于该体的解析体积，（b）是否真的水密
（闭合、流形、无重复面、绕向一致且外向、亏格 0）。两者是**两个独立结论**，不能互相替代——
下文 `unwelded_soup` 控制组即为此设：它的体积误差严格为 0，拓扑却是敞开的。

**ref = 解析闭式，全部在代码里求值（零手抄）**，每一条都有明确出处：

| 量 | 闭式 | 出处 |
|---|---|---|
| 球体积 | `V = 4/3·π·r³` | 立体几何闭式；本报告不写十进制常数，由 `math.pi` 求值 |
| 球体积（第二推导） | `V = (2/3)·(π·r²·2r)` | Euclid《几何原本》XII.18（球 = 外切圆柱的 2/3，穷竭法），代码内独立求值，r=1 处两式相对差 0.0e+00（逐位相同） |
| 球表面积 | `A = 4·π·r²` | 同上闭式 |
| 长方体体积/表面积 | `V = a·b·c`，`A = 2(ab+bc+ca)` | 同上闭式 |
| 正二十面体体积（R 外接球） | `V = (5/12)(3+√5)·a³`，`a = 4R/√(10+2√5)` | 正多面体闭式（Euclid XIII）；这正是 `sub0` icosphere 的自身闭式 |

**为什么 ref 不是"手册上某个数"**：球/长方体的解析体积是无量纲闭式而非表格数值，任何查表都
只能查到 10 位有效数字并引入手抄风险；本任务因此把 ref 定义为**代码求值的闭式**，并用两条
独立推导（`4/3πr³` 与 Euclid 外切圆柱式）互相校验出处，差值 0.0e+00（r=1 处逐位相同；r=1.234 时为 2.3e-16，即只差舍入）。
正二十面体一行是第三重校验：它是 20 个**非轴对齐**三角形的闭式，能把体积积分器从"只会算六个平面"里拉出来。

**参考对照**：仓库内既有 icosphere 基准值 `4.047044679978`（`docs/board_audit_m17634/`，sub2, r=1，
即本报告阶梯的 sub2 行 4.047044679979）——两者一致到 1e-12，说明本次阶梯与
既有证据链同源。

## ② 方法与实现（复用现有能力，不另起炉灶）

| 环节 | 复用的 TensorCAD 能力 |
|---|---|
| 生成闭体 | `tensorcad.mesh_xval.cube`（长方体）、`icosphere`（球，inscribed 细分） |
| 体积（主路线） | `tensorcad.mesh_xval.TriMesh.volume()`：散度定理 `V = 1/6 Σ a·(b×c)` |
| 体积（第二路线） | 新增 `tensorcad.mesh_volume.volume_by_tet_centroid`：以网格质心为顶点的四面体有符号求和（不同求和次序/不同量级，用于捕捉符号与索引错误） |
| 水密性判定 | `tensorcad.geometry_audit.audit_watertightness`（自由边、非流形边、重复面、绕向冲突、有符号体积、缝距） |
| 缺陷注入 | `geometry_audit.inject_hole / inject_flip_faces / inject_duplicate_faces`（负控） |
| 第三方交叉 | `trimesh`（OSS）对同一批三角形独立算体积 |
| 文件级核验 | `mesh_xval.write_stl` 写二进制 STL → `trimesh` 读回 → `mesh_xval.weld` 焊接后复测 |
| 新增薄封装 | `tensorcad/mesh_volume.py`：解析闭式 + 闭体检查记录 + 收敛阶梯 + 负控；`benchmarks/verified/mesh_volume_watertight/run_benchmark.py`：判据与工件 |

**设计要点**：把误差拆成**积分器误差**与**离散化误差**。长方体是三角网格能精确表示的体，容差因此定在
1e-12（只许浮点舍入）；球的三角化是**内接**多面体，体积必然偏小，误差随 h² 下降，所以只对阶梯
最细一档（sub6）判 `<0.1 %`，并把整条阶梯（含未达标档）如实列出。

**预声明阈值（跑之前固定，写在模块常量里，不看结果再调）**：

- `EXACT_TOLERANCE = 1e-12`：长方体族（精确可表示，只许舍入）；
- `RELATIVE_TOLERANCE = 0.001`（= 0.1 %）：曲面体最细档与缩放档；
- `ROUTE_TOLERANCE = 1e-10`：两条求和路线互校；
- `FLOAT32_TOLERANCE = 1e-06`：二进制 STL（float32，相对精度 ~1.2e-7）回读，误差源是文件格式而非积分器；
- 观测阶带 `[1.5, 2.5]`：分段线性逼近光滑曲面的二阶；
- **反证判据**：sub4 一档**必须**超出 0.1 %（若粗糙内接多面体"给出"球体积，则测量与几何脱钩）。

## ③ 参数表（无量纲/量纲完整）

| 项 | 取值 |
|---|---|
| 长方体算例 | 2×2×2 m³；1.5×0.75×3 m³；2×2×2 m³ 绕 Z-Y-X 旋转 (31°, −17°, 73°) 再平移 (13.7, −42.1, 7.3) m |
| 球半径 | r = 1.0 m（阶梯）；r = 2.5 m（无量纲/缩放校验） |
| 球离散 | icosphere 细分 level = 0…6，面数 20·4^level（20 … 81920） |
| 判据档位 | volume/area @ sub6；反证档 sub4；第二路线一致性全档 |
| STL 回读算例 | 1.7×1.1×2.9 m³（半边长 0.85/0.55/1.45 在 float32 下不可精确表示，故格式损失可见） |
| 审计参考长度 | 球的直径 2r；长方体取 max(a,b,c) |
| 容差 | 见 ② 预声明；无量纲量：体积/表面积相对误差、观测阶 p |

## ④ 对照表（ref vs got + 误差逐项）

| 算例 | 判据量 | ref（解析闭式，代码求值） | got（TensorCAD 测得） | 相对误差 | 容差 | ✓/✗ |
|---|---|---|---|---|---|---|
| `box_2x2x2` | volume | 8 | 8 | +0.0000 % | 1e-12 | ✓ |
| `box_2x2x2` | area | 24 | 24 | +0.0000 % | 1e-12 | ✓ |
| `box_2x2x2` | route_agreement | 8 | 8 | +0.0000 % | 1e-10 | ✓ |
| `box_2x2x2` | topology_closed | 1 | 1 | — | 精确判定 | ✓ |
| `box_2x2x2` | trimesh_volume | 8 | 8 | +0.0000 % | 1e-09 | ✓ |
| `box_1.5x0.75x3` | volume | 3.375 | 3.375 | +0.0000 % | 1e-12 | ✓ |
| `box_1.5x0.75x3` | area | 15.75 | 15.75 | +0.0000 % | 1e-12 | ✓ |
| `box_1.5x0.75x3` | route_agreement | 3.375 | 3.375 | +0.0000 % | 1e-10 | ✓ |
| `box_1.5x0.75x3` | topology_closed | 1 | 1 | — | 精确判定 | ✓ |
| `box_1.5x0.75x3` | trimesh_volume | 3.375 | 3.375 | +0.0000 % | 1e-09 | ✓ |
| `box_2x2x2_rot_trans` | volume | 8 | 8 | -0.0000 % | 1e-12 | ✓ |
| `box_2x2x2_rot_trans` | area | 24 | 24 | +0.0000 % | 1e-12 | ✓ |
| `box_2x2x2_rot_trans` | route_agreement | 8 | 8 | +0.0000 % | 1e-10 | ✓ |
| `box_2x2x2_rot_trans` | topology_closed | 1 | 1 | — | 精确判定 | ✓ |
| `box_2x2x2_rot_trans` | trimesh_volume | 8 | 8 | +0.0000 % | 1e-09 | ✓ |
| `sphere_r1_sub0` | icosahedron_closed_form | 2.53615071012 | 2.53615071012 | +0.0000 % | 1e-12 | ✓ |
| `sphere_r1_sub6` | volume | 4.18879020479 | 4.18822373818 | -0.0135 % | +0.1000 % | ✓ |
| `sphere_r1_sub6` | area | 12.5663706144 | 12.5654311425 | -0.0075 % | +0.1000 % | ✓ |
| `sphere_r1_sub6` | route_agreement | 4.18822373818 | 4.18822373818 | +0.0000 % | 1e-10 | ✓ |
| `sphere_r1_sub6` | topology_closed | 1 | 1 | — | 精确判定 | ✓ |
| `reference_provenance` | sphere_volume_two_derivations | 4.18879020479 | 4.18879020479 | +0.0000 % | 1e-15 | ✓ |
| `sphere_r1_sub4` | coarse_must_exceed_tolerance | 0.001 | 0.00216083 | — | 精确判定 | ✓ |
| `sphere_ladder` | observed_order_volume | 2 | 1.98876 | — | 精确判定 | ✓ |
| `sphere_r2.5_sub6` | volume | 65.4498469498 | 65.4409959091 | -0.0135 % | +0.1000 % | ✓ |
| `scaling_invariance` | volume_error_difference | -0.000135233940842 | -0.000135233940842 | -2.019e-16 | +1.000e-12 | ✓ |
| `box_stl_roundtrip` | volume | 5.423 | 5.42300044799 | +0.0000 % | 1e-06 | ✓ |
| `control:hole` | defect_detected | 1 | 1 | — | 精确判定 | ✓ |
| `control:flipped_facet` | defect_detected | 1 | 1 | — | 精确判定 | ✓ |
| `control:duplicate_facet` | defect_detected | 1 | 1 | — | 精确判定 | ✓ |
| `control:unwelded_soup` | defect_detected | 1 | 1 | — | 精确判定 | ✓ |
| `blind_spot:duplicate_single` | defect_detected_by_audit | 1 | 1 | — | 精确判定 | ✓ |
| `blind_spot:flipped_single` | defect_detected_by_audit | 1 | 1 | — | 精确判定 | ✓ |

**结论**：全部 32 条判据 通过；体积/表面积关键量相对误差
最差为曲面体最细档 -0.0135 %（|·| = 1.35e-04），
**是否 <0.1 %：是**；**是否 <3 %：是**。
长方体族体积误差 +0.000e+00（= 0，精确可表示）；旋转平移后 -8.327e-15（1e-14 量级，舍入底噪）。
正二十面体闭式（2.536150710120）与散度定理实测（2.536150710120）一致到 1.8e-16（相对），
即在 20 个非轴对齐三角形上体积积分器已被独立验证。

## ⑤ 收敛性与参数敏感性（≥2 点）

**（1）离散化阶梯**（主敏感轴，h = 平均边长）：

| level | 面数 | h（平均边长, m） | V_mesh | V 相对误差 | A 相对误差 | 开边 | χ |
|---|---|---|---|---|---|---|---|
| sphere_r1_sub0 | 20 | 1.0515e+00 | 2.536150710120 | -39.4539 % | -23.8082 % | 0 | 2 |
| sphere_r1_sub1 | 80 | 5.8228e-01 | 3.658712208512 | -12.6547 % | -7.1655 % | 0 | 2 |
| sphere_r1_sub2 | 320 | 2.9933e-01 | 4.047044679979 | -3.3839 % | -1.8822 % | 0 | 2 |
| sphere_r1_sub3 | 1280 | 1.5073e-01 | 4.152740817093 | -0.8606 % | -0.4765 % | 0 | 2 |
| sphere_r1_sub4 | 5120 | 7.5499e-02 | 4.179738947995 | -0.2161 % | -0.1195 % | 0 | 2 |
| sphere_r1_sub5 | 20480 | 3.7766e-02 | 4.186524949279 | -0.0541 % | -0.0299 % | 0 | 2 |
| sphere_r1_sub6 | 81920 | 1.8885e-02 | 4.188223738180 | -0.0135 % | -0.0075 % | 0 | 2 |

观测阶（log|误差| 对 log h 最小二乘）：体积 **p = 1.989**，表面积 **p = 2.005**
（理论二阶，分段线性内接）。每细分一档误差降 ~4 倍，符合二阶；首个进入 0.1 % 的档位 = **`sphere_r1_sub5`**。

**（2）尺度不变性**（误差应是无量纲量）：r = 1.0 与 r = 2.5 在 sub6 的体积相对误差之差
-2.019e-16（容差 1e-12）——r 改变 2.5 倍，相对误差不动，说明误差来自拓扑/角度而非量级。

**（3）刚体运动不变性**：长方体绕三轴旋转（31°,−17°,73°）并平移到远离原点处，体积误差
-8.327e-15、表面积 +5.921e-16（相对 28.2 的量级均为 1e-15 级），说明体积路线对坐标系无依赖。

**（4）容差档位敏感性（重要）**：若把 `<0.1 %` 判据施加到**所有**档位，sub0…sub4 会全部判负
（sub0 -39.4539 %、sub4 -0.2161 %）。这不是缺陷而是内接多面体的
离散化性质；因此判据档位与反证档同时预声明，避免"看到结果再选分辨率"。

**（5）分辨率与"看不见的缺陷"**：在最细档 sub6（4.188223738 m³）上只注入**一个**重复面，
体积效应 +1.153e-05（远小于 0.1 % 判据，`invisible_to_volume_gate = True`），审计却报
`nonmanifold_edge, duplicate_facets`（非流形边 3 条）；只翻**一个**面，体积效应 -2.307e-05，审计报
`normal_inconsistency`（绕向冲突 3 条）。也就是说，判据越"体量宽松"，水密性审计越不可省。

**（6）"精确 vs 文件格式"**：同一长方体在内存中误差为 0，写成二进制 STL（float32）回读后误差
+8.261e-08（0.7 × float32 eps）。0.1 % 判据仍以两个数量级余量通过，但
说明"文件经手"会引入与积分器无关的误差项，故单列 `FLOAT32_TOLERANCE`。

## ⑥ 失败模式与工件

**负控（检出能力证明）**——干净球（sub2, 320 面）注入四类缺陷：

| 注入缺陷 | 做法 | 审计判定 | 命名 defect | 体积相对误差 | 检出 |
|---|---|---|---|---|---|
| `hole` | facets removed inside an x/z window | open=16, nonmanifold=0, dup=0, winding=0 | seam_gap, seam_gap | -2.3666 % | ✓ |
| `flipped_facet` | facet winding reversed on 4 facets | open=0, nonmanifold=0, dup=0, winding=12 | normal_inconsistency | -2.3766 % | ✓ |
| `duplicate_facet` | 3 facets repeated | open=0, nonmanifold=9, dup=3, winding=0 | nonmanifold_edge, duplicate_facets | +0.8958 % | ✓ |
| `unwelded_soup` | connectivity destroyed (STL before welding) | open=960, nonmanifold=0, dup=0, winding=0 | component_count | +0.0000 % | ✓ |

四类全部检出且缺陷名词与 `geometry_audit` 的既有词表一致（`seam_gap` / `normal_inconsistency` /
`nonmanifold_edge` + `duplicate_facets` / `component_count`）。关键发现：`unwelded_soup` 的体积误差严格为 0
（+0.0000 %）却敞开 960 条边——**体积核验与水密性核验不可互替**；
反过来 `flipped_facet` 只翻 4/320 个面就让体积偏 -2.3766 %，说明体积数值对绕向敏感，可作水密性的旁证但不能替代审计。

**失败模式（如实记录）**：

1. **粗分辨率内接球不满足 0.1 %**：sub≤4 一律超出（sub4 -0.2161 %）。任何把三角化球的体积当"解析球体积"用的下游流程，在该分辨率下会带 0.2 % 系统偏低。
2. **STL 是三角形汤**：未焊接的 STL 在体积上完全正确、在拓扑上完全敞开；`mesh_xval.weld` 是必要前置步骤，本算例把它写进流程。
3. **低分辨率下体积判据有盲区**：同一类缺陷（1 个重复面 / 1 个翻面）在粗网格 sub2 上把体积推动
   +0.8958 % / -2.3766 %（超判据），在 sub6 上只有
   +1.153e-05 / -2.307e-05（**在**判据内）。因此"体积误差 <0.1 %"**不能**当作水密性证明：
   网格越细，拓扑缺陷在体积上越隐形，而审计的检出不受分辨率影响（两档均检出）。
4. **生态摩擦（如实记）**：`/tensor/knowledge/search` 后端 engine 离线（`127.0.0.1:8021 connection refused`），`/tensor/knowledge/assets` 报 `asset_dir_missing`（`/root/tensor-ecosystem/TensorKnowledge/artifacts/ke_export`），故开工三查的知识库/OSS 检索两路**不可用**，只能退化为：(a) 仓库内既有基准与文档检索；(b) negative_knowledge_ledger 检索。台账检索到的同域负知识均为"result 未提供 checks.reason"的流程性失败（m19360/m19040/m18843/m18671/m18606/m17355），本单据此在 finish 的 checks 中显式给出判据与证据位置。
5. **既有知识台账一处预先存在的缺口（非本单引入）**：`python -m pytest tests/test_knowledge_governance.py` 的
   `test_every_entry_reproduces_today` 在 `nurbs_intersect.sphere_plane_circle` 一条上失败，原因是它 pin 的工件
   `artifacts/m17346/nurbs_intersect/results.json` 不在当前检出中（MISSING）。本单新增的 3 条 `mesh_volume.*` 台账
   条目逐条 `verify_entry` 均为 PASS（box 精确/二十面体闭式/二阶观测阶），未触碰该条。如实记录，未代为修复。
6. **`trimesh` 体积与 TensorCAD 同源不同实现**：长方体族一致到 ≤1.3e-15（相对），可作第三方交叉；但它对同一批三角形复算，不能替代解析 ref。

**工件（均为真实文件，`git` 内可查）**：

| 工件 | 内容 |
|---|---|
| `benchmarks/verified/mesh_volume_watertight/result.json` | 全部判据、体记录、阶梯、负控、环境版本 |
| `benchmarks/verified/mesh_volume_watertight/convergence.csv` | 阶梯逐档 h / V / 误差 / 开边 / χ |
| `benchmarks/verified/mesh_volume_watertight/convergence.png` | log-log 误差-步长 + 斜率 2 参考线 + 0.1 % 线 |
| `benchmarks/verified/mesh_volume_watertight/gates.png` | 每条判据的 |相对误差| 与 0.1 % 线（绿/红） |
| `benchmarks/verified/mesh_volume_watertight/box_stl.stl` | 经 writer 写出的二进制 STL（float32）回读样本 |
| `tensorcad/mesh_volume.py` | 新增能力：解析闭式、闭体检查、收敛/观测阶、负控 |
| `tests/test_mesh_volume.py` | 12 项测试（含整算例跑一遍并断言全部判据） |

## 一键复现

```bash
python benchmarks/verified/mesh_volume_watertight/run_benchmark.py --out benchmarks/verified/mesh_volume_watertight   # 约 18 s，判据不过则退出码 1
python -m pytest tests/test_mesh_volume.py -q             # 12 passed
python benchmarks/verified/mesh_volume_watertight/publish_report.py   # 重新生成本文（数字仍来自 result.json）
```

运行环境：Python 3.11.6，numpy 2.4.6，
trimesh 5.0.0，matplotlib 3.10.8；
平台 Linux-6.6.52-1.7.lns23.loongarch64-loongarch64-with-glibc2.38；算例耗时 17.88 s（单机、单线程 numpy，<8 核，无需集群）。
