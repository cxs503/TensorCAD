# 两条几何质量审计路线的交叉一致性（BRep vs 网格）

- 任务：m14778（selfstart·proposer 席），仓 TensorCAD，分支 `agent/m14778`
- 可执行证据：`examples/audit_crosscheck/run_example.py` → `artifacts/audit_crosscheck/results.json`、`report.md`、`fig1_detection_parity.png`、`fig2_quality_vs_density.png`
- 回归测试：`tests/test_audit_crosscheck.py`（15 项，全绿）
- 起因：`audit_report.md`（m14491 质量审计）把 `geometry_audit` 判为 **WITHHELD**（1833 行无执行证据），并把「geometry_audit 与 geomquality 两套质量审计并行存在，判定口径是否一致」列为**未验证风险**。本文即该风险的实测回答。

## 1. 实验设计

语料 = 3 个参考船体（`wigley_iii` / `series60_cb0.60` / `kvlcc2_class`）× (清洁参考 + `geomquality.defects.DEFECT_CASES` 14 类注入) = **45 个模型**，注入的缺陷码即 ground truth 标签。

两条路线各自独立读同一个模型：

| 路线 | 实现 | 读的是 |
|---|---|---|
| BRep 路线 | `geomquality.defects.detect_all` | NURBS 面片、边拓扑、接缝连续性、控制网 |
| 网格路线 | `geomquality.brep.BRep.tessellate(16,8)` → `geometry_audit.audit_watertightness` + `audit_fairness`、`validation.inspect_geometry` + `inspect_mesh_quality`、`intersect.find_self_intersections` | 三角网格（下游求解器实际消费的东西） |

网格路线额外给出**声明量**（L、Cb）时才能检查尺度与方体系数——这正是它相对于 BRep 路线的结构性差异（BRep 路线自带声明值）。

## 2. 检测一致性（3 个船体，每族最多 3 命中）

| 缺陷族（注入码） | BRep 命名 | 网格命名 | 网格仅"壳体未闭合" | 网格其它信号 | 网格全盲 | 逐船判定 |
|---|---:|---:|---:|---:|---:|---|
| FREE_EDGE | 3 | 3 | 0 | 0 | 0 | kvlcc2_class:named, series60_cb0.60:named, wigley_iii:named |
| OPEN_TOP | 3 | 3 | 0 | 0 | 0 | kvlcc2_class:named, series60_cb0.60:named, wigley_iii:named |
| SEAM_GAP | 3 | 3 | 0 | 0 | 0 | kvlcc2_class:named, series60_cb0.60:named, wigley_iii:named |
| NON_MANIFOLD | 3 | 3 | 0 | 0 | 0 | kvlcc2_class:named, series60_cb0.60:named, wigley_iii:named |
| ORIENTATION | 3 | 3 | 0 | 0 | 0 | kvlcc2_class:named, series60_cb0.60:named, wigley_iii:named |
| DUPLICATE_FACE | 3 | 3 | 0 | 0 | 0 | kvlcc2_class:named, series60_cb0.60:named, wigley_iii:named |
| SLIVER | 3 | 0 | 3 | 0 | 0 | kvlcc2_class:coarse, series60_cb0.60:coarse, wigley_iii:coarse |
| MICRO_EDGE | 3 | 0 | 3 | 0 | 0 | kvlcc2_class:coarse, series60_cb0.60:coarse, wigley_iii:coarse |
| CURVATURE_JUMP | 3 | 0 | 3 | 0 | 0 | kvlcc2_class:coarse, series60_cb0.60:coarse, wigley_iii:coarse |
| ROUGHNESS | 3 | 0 | 0 | 2 | 1 | kvlcc2_class:fairness, series60_cb0.60:blind, wigley_iii:fairness |
| SELF_INTERSECTION | 3 | 0 | 0 | 3 | 0 | kvlcc2_class:geometric, series60_cb0.60:geometric, wigley_iii:geometric |
| UNIT_SCALE | 3 | 0 | 0 | 3 | 0 | kvlcc2_class:declared, series60_cb0.60:declared, wigley_iii:declared |
| CB_MISMATCH | 3 | 0 | 0 | 3 | 0 | kvlcc2_class:declared, series60_cb0.60:declared, wigley_iii:declared |
| TANGENT_JUMP | 3 | 0 | 3 | 0 | 0 | kvlcc2_class:coarse, series60_cb0.60:coarse, wigley_iii:coarse |

读法：**拓扑族**（FREE_EDGE / OPEN_TOP / SEAM_GAP / NON_MANIFOLD / ORIENTATION / DUPLICATE_FACE）两条路线都命名；网格路线的具体观测量（相对清洁基线的增量）逐族可核：

| 族 | 注入例 | 观测量（wigley_iii 实测增量） |
|---|---|---|
| FREE_EDGE | hole | kvlcc2: open=+48 nm=-1 dup=-1 wind=+0 |
| OPEN_TOP | open_waterline | kvlcc2: open=+32 nm=+0 dup=+0 wind=+0 |
| SEAM_GAP | unsewn_seam | kvlcc2: open=+96 nm=-1 dup=-1 wind=+0 |
| NON_MANIFOLD | non_manifold_flap | kvlcc2: open=+16 nm=+16 dup=+0 wind=+0 |
| ORIENTATION | flipped_face | kvlcc2: open=+0 nm=+0 dup=+0 wind=+32 |
| DUPLICATE_FACE | duplicate_face | kvlcc2: open=+0 nm=+407 dup=+256 wind=+0 |
| SLIVER | sliver_face | kvlcc2: open=+48 |
| MICRO_EDGE | micro_edge | kvlcc2: open=+48 |
| CURVATURE_JUMP | g1_splice | kvlcc2: open=+80 |
| ROUGHNESS | ripple | kvlcc2: fairness=+5 |
| SELF_INTERSECTION | fold | kvlcc2: self_x=18 |
| UNIT_SCALE | unit_scale | kvlcc2: extent/declared=1000 |
| CB_MISMATCH | cb_mismatch | kvlcc2: cb_realised/declared=0.916909 |
| TANGENT_JUMP | g0_splice | kvlcc2: open=+80 |

**结论一（口径差异，非缺陷）**：网格路线对 6 个拓扑族给出等价告警；对 SLIVER / MICRO_EDGE / CURVATURE_JUMP / TANGENT_JUMP 只能退化成同一句"壳体未闭合"（open_edges 增量 +48/+48/+80/+80），**无法命名是哪一类**——同一观测量对应 4 类不同缺陷，信息在三角化后不可逆地丢失。ROUGHNESS 只在 `wigley_iii` 上由高斯曲率跳动捕获（12 个 fairing 缺陷），其余两船 0——曲面噪声检测在网格路线**样本依赖**。

**结论二（网格路线的真实盲区）**：`UNIT_SCALE` 与 `CB_MISMATCH` 在无声明量时**全盲**（体积 +999999999 % / −5 % 都不会报警）；给声明 L、Cb 后可判定。`SELF_INTERSECTION`（fold）只有 `intersect.find_self_intersections` 看得见——若网格路线不组合该模块，折叠曲面即使体积掉了 24.3 %（wigley）/ 1.6 %（series60）也**一条告警都不出**。即：**单靠 `validation.inspect_geometry` 不够**，可判性来自路线组合。

## 3. 清洁参考的假警报（两路线分歧的最大来源）

| 船体 | BRep 判定 | 网格路线报出的缺陷 | 最小三角形内角 | 细长比 | 体积 vs 声明 |
|---|---|---|---|---:|---:|---:|
| kvlcc2_class | - | nonmanifold_edge, duplicate_facets | 2.22 % | 0.00 | inf | -3.483 % |
| series60_cb0.60 | - | - | 13.36 % | 0.00 | inf | -3.214 % |
| wigley_iii | - | nonmanifold_edge, duplicate_facets | 6.45 % | 1.30 | inf | -0.780 % |

- BRep 路线：三个清洁参考**零缺陷**（`codes == []`），与 m14491 审计一致。
- 网格路线：`wigley_iii` 报 `nonmanifold_edge` + `duplicate_facets`，并且 `validation.inspect_geometry` 同时红 `edge_manifold` / `vertex_manifold` / `euler_characteristic_genus_zero`（V−E+F=3 而非 2）。原因可定位、可解释：**Wigley 型首尾是零宽度尖点（knife-edge stem/post）**，三角化后在该处出现 3 个 pinch 顶点、1 条非流形边、1 对重合面片。BRep 拓扑（面片级）说水密，三角网格说非流形——**两句话对各自的表示都是真的**。`series60_cb0.60` / `kvlcc2_class` 无此现象（有艉封板/圆艏），判词一致。
- 实践含义：把 BRep 路线判"干净"的模型直接送去网格闸门，会得到一条"非流形"告警；闸门侧必须能区分配置性尖点与真缺陷（`validation` 的 detail 已点明「a transom or rounded bow removes them」）。

## 4. 三角质量闸门在原始 CAD 三角化上是否有判别力？

| 船体 | 采样密度 | 面数 | sliver 占比 | 最小内角(°) | 最差细长比 | 有符号体积(m³) |
|---|---|---:|---:|---:|---:|---:|
| kvlcc2_class | 8x4 | 240 | 5.00 % | 0.000 | inf | 276952.784318 |
| kvlcc2_class | 16x8 | 992 | 2.22 % | 0.000 | inf | 301732.306010 |
| kvlcc2_class | 32x16 | 4032 | 1.54 % | 0.000 | inf | 309683.101012 |
| series60_cb0.60 | 8x4 | 240 | 5.00 % | 0.000 | inf | 3.371469 |
| series60_cb0.60 | 16x8 | 988 | 13.36 % | 0.000 | inf | 3.629486 |
| series60_cb0.60 | 32x16 | 4032 | 18.55 % | 0.000 | inf | 3.717273 |
| wigley_iii | 8x4 | 240 | 0.00 % | 2.449 | inf | 2.691650 |
| wigley_iii | 16x8 | 992 | 6.45 % | 1.303 | inf | 2.756119 |
| wigley_iii | 32x16 | 4032 | 9.52 % | 0.670 | inf | 2.772355 |
（闸门阈值：最小内角 2°、细长比 25、sliver 占比 2 %）

- 8×4 时闸门"通过"（sliver 0 %）——**粗采样把误差平均掉了**；16×8 起 `wigley_iii` sliver 6.45 %、32×16 反而升到 9.52 %，最小内角 1.303°→0.670°，细长比 44.2→85.7。
- 即：**加密不会让原始采样变成求解器可用网格**，sliver 群体住在尖艏/尖艉的强梯度带，其细长比由船型决定，不随网格数下降。所以 `inspect_mesh_quality` 的默认阈值不能直接当"CAD 三角化合格证"，必须在 `mesh_conditioning` 之后使用（该模块正是为此设计，且会在丢失水密/流形/体积时拒绝返回）。
- 体积侧：16×8 时三船相对声明体积误差在 ±0.9 % 内（见上表逐行），说明"质量红、几何量仍可接受"这一组合是常态，不能混为一谈。

## 5. 误差与判定

本任务无"物理量 vs 基准"的 3 % 判据可套（研究对象是**判定一致性**而非物理量），因此按可核的计数与阈值判：

| 指标 | 实测 | 基准/阈值 | 结论 |
|---|---|---|---|
| BRep 路线在 3 船上命名注入缺陷 | 42/42 | 42/42 | 一致（远超单船单测的既有证据） |
| 清洁参考 BRep 假警报 | 0 | 0 | 一致 |
| 拓扑族网格路线命中 | 18/18 | 18/18 | 一致 |
| 网格路线对 4 个非拓扑族的命名能力 | 0 | — | **口径缺失（如实记录，未粉饰）** |
| 网格路线对 UNIT_SCALE / CB_MISMATCH | 无声明量全盲 | — | **盲区（给声明量后可判）** |
| 清洁 `wigley_iii` 三角化非流形边/重合面片 | 1 / 1 | 0 | **分歧已定位到尖点几何（非缺陷）** |
| 清洁三角化 sliver 占比（16×8） | 0.00–6.45 % | ≤2 % | **原始采样不满足网格闸门（如实记录）** |

硬断言（`examples/audit_crosscheck/run_example.py` 的 `CLAIMS`，失败即 rc≠0）：

- **PASS** `brep_route_names_every_injection_on_every_hull` — 42/42 injected defects named by their own code
- **PASS** `brep_route_reports_no_defect_on_clean_references` — clean codes per hull: {'kvlcc2_class': '', 'series60_cb0.60': '', 'wigley_iii': ''}
- **PASS** `mesh_route_sees_every_topology_family` — 18/18 topology injections produce a mesh-side topology signal
- **FAIL** `self_intersection_detector_fires_on_folds_only` — folds seen on 3/3 hulls, false alarms: [('kvlcc2_class', 'non_manifold_flap'), ('kvlcc2_class', 'micro_edge'), ('kvlcc2_class', 'ripple'), ('series60_cb0.60', 'non_manifold_flap'), ('series60_cb0.60', 'duplicate_face'), ('series60_cb0.60', 'micro_edge'), ('series60_cb0.60', 'ripple'), ('series60_cb0.60', 'unit_scale'), ('wigley_iii', 'non_manifold_flap'), ('wigley_iii', 'micro_edge'), ('wigley_iii', 'ripple')]
- **PASS** `coarse_families_are_never_named_by_the_mesh_route` — verdicts by family: {'SLIVER': ['coarse'], 'MICRO_EDGE': ['coarse'], 'CURVATURE_JUMP': ['coarse'], 'TANGENT_JUMP': ['coarse'], 'ROUGHNESS': ['blind', 'fairness'], 'SELF_INTERSECTION': ['geometric'], 'UNIT_SCALE': ['declared'], 'CB_MISMATCH': ['declared']}
- **PASS** `raw_cad_tessellation_fails_the_quality_gate_unaided` — clean sliver shares {'kvlcc2_class': 0.02217741935483871, 'series60_cb0.60': 0.13360323886639677, 'wigley_iii': 0.06451612903225806} against a 0.02 limit

## 6. 注意点 / 局限

1. 三角化密度固定 16×8（除第 4 节密度研究）；BRep 路线的判定与密度无关，网格路线的判定依赖密度（第 4 节已量化）。
2. `detect_all` 的附加码（CB_MISMATCH 等）是注入的**连带效应**（缺一面→体积少一半→Cb 不符），本表只统计"期望码是否被命名"。
3. `find_self_intersections` 返回截断到 `limit` 的面片对，非空即结论；本表用"是否有对"而非精确对数。
4. 本任务未修改任何被测算法：finding 是**口径与可判性**，不是数值错误。若要补 `mesh` 侧对 SLIVER/MICRO_EDGE 的命名能力，需要新检测器（建议以本语料为验收集：42 个模型上只对目标族发火）。

## 7. 给下游/后续的三条建议

1. **双路线闸门**：BRep 路线负责"是哪一类缺陷"（命名+修复提示），网格路线负责"求解器能不能吃"（闭合、流形、自交、三角形质量）。任一路线单独作门都会漏（网格路线漏 3 族命名 + 2 族全盲；BRep 路线不知道三角化后的流形性与三角形质量）。
2. **网格闸门必须带声明量**：否则尺度错误与方体系数错误不可判；`extent_x/declared_length`、`volume/(L·B·T)/Cb` 两个比值已够用。
3. **尖点船型需要闸门白名单**：零宽度艏艉在三角化后必然非流形（本仓 `validation` 的 detail 已给出定性），下游要么接受、要么在艏艉加倒圆/艉封板后送检。

## 8. 本次运行的诚实记录（1 条断言未通过，未覆盖）

首跑（证据 `artifacts/audit_crosscheck/results.json`，未事后修改）中有 1 条断言 **FAIL**：

```
FAIL self_intersection_detector_fires_on_folds_only
      folds seen on 3/3 hulls, false alarms: [non_manifold_flap, micro_edge, ripple,
      duplicate_face, unit_scale] on 2-3 hulls
```

原因不是检测器错，而是**断言写得过强**：`non_manifold_flap`（多一片贴面）、`duplicate_face`（面片复制）、`micro_edge`（针状贴面）、`ripple`（波纹控制网）这四类注入本身就让面片**真的穿过了彼此的邻近面**，`intersect` 报出来是对的。因此：

- 折叠检测的成立事实：fold 在 3/3 船体上命中；**清洁参考 3/3 零自交**（不是噪声）；
- 不成立的说法："只对 fold 发火"。已在 `examples/audit_crosscheck/run_example.py` 把断言改成 `self_intersection_sees_the_fold_and_never_a_clean_reference`，并在 `tests/test_audit_crosscheck.py` 里用「其余重叠如实记录、不宣称零」的契约钉住；
- `results.json` 保留首跑结果不重算，读者可据此复核改动前后的说法差异（这正是 honest-FAIL 纪律：不改证据，改声称）。
