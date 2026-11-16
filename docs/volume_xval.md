# 几何体积跨工具互证：自研 NURBS 积分 vs 独立实体引擎（XVAL-GEOM-VOLCROSS）

`tensorcad.volume_xval` 把"某个体积数字可信"从**单一积分器的自我一致**升级为**跨实现互证 + 差异归因**。
动机是 JBC 事件：管线给出的排水体积比公布值低 −18.3 %，而所有内部守卫都同意管线——因为守卫和管线
共用同一条积分路径。本模块的立场是结构性的：**每个管线体积必须被至少两条不共享代码的实现复现，
且实现之间的差异必须归因到具名效应（细分亏欠 / STL 量化 / 拓扑缺陷），不允许落进一个宽容差里了事**。

## 1. 四条路径

| 路径 | 实现 | 与本仓共享 |
|---|---|---|
| `pipeline` | 该船型的平台自研积分：NURBS 片上的 Gauss 求积（`hull_nurbs`），或 `tensorcad.mesh` 的闭曲面积分 | —（被检验对象） |
| `engine A` | `trimesh.Trimesh.volume`（第三方，NumPy 向量化，自带水密/绕向守卫） | 无 |
| `engine B` | `numpy-stl` `get_mass_properties`（Mirtich 多面体质量特性，C++），经二进制 STL 往返 | 无 |
| `smooth` | 半宽场求积：`A(x)=2∫b dz` 再用 Simpson 积分 `V=∫A dx`。无三角形、无散度定理 | 仅半宽场 |

`smooth` 是回答 JBC 失效模式的那条：三条三角形路径积分的都是**内接多面体**，它们会彼此一致地
系统性偏低；`smooth` 直接积光滑几何，因此能区分"几何源错了"与"三角形太粗"。
`engine B` 走二进制 STL（float32 坐标），这就是 `ENGINE_AGREEMENT_TOLERANCE = 1e-6`（而不是 1e-12）
的原因——这是**工具的**精度预算，被如实报告而不是被"修正"。

## 2. 判定逻辑（`classify_routes`）

按可证明强度排序：

- `engine-disagreement`：两条独立引擎读**同一批三角形**却超出 1e-6 —— 拓扑缺陷（裂缝/绕向错/重复顶点）
  或精度假象，绝不是细分效应；
- `tessellation-deficit`：引擎彼此一致，管线与引擎之差随加密以**实测阶**收敛（窗口 1.5–2.6）——
  弯曲船体内接多面体的 O(h²) 亏欠，必须外推掉而不是调参调掉；
- `systematic-bias`：差异不随加密缩小 —— 不是网格属性，正是 JBC 那一类，必须记 negative-knowledge ledger；
- `undetermined`：引擎一致、管线超线、但**没有提供加密梯** —— 差距没有实测阶，报"未归因"而不是
  断言数据证明不了的东西（调用方被告知要补跑什么）；
- `agrees`：落在 0.3 % 验收线内。

单靠网格研究无法评价它所对齐的求积本身，所以管线求积在自己的轴上单独加密（Gauss 阶 8×6 → 48×32），
证明"被比较的管线数字是站住的"。

## 3. 船型与基准

| 船型 | 基准 | 性质 |
|---|---|---|
| Wigley（ITTC 比例，n=m=2） | `4 L B T / 9` 解析闭式 | 精确基准 |
| KCS 尺度的 NURBS 船（全尺度，Cb 标定到 0.6505） | 公布排水体积 52 030 m³（KRISO / SIMMAN 2008 / Tokyo 2015） | 公布值，但**形式是按公布 Cb 标定的**，故与 `Cb·L·B·T` 的一致是标定恒等式，不能冒充独立验证 |
| KCS 尺度的 NURBS 船（模型尺度，Cb 0.555，IN383-CAD01 案例） | 无（WITHHELD） | 之前只与同一片生成的网格自比，即 WITHHELD 声明；本单用跨工具一致替代外部基准 |
| 非对称艉板船（`transom_asymmetric`，左右半宽不同） | `analytic_asymmetric_hydrostatics` 闭式 | 积分器恒等性检验组（该船无对称半宽场，进不了主表） |

## 4. 结果（由 `tools/volume_xval_report.py` 生成，禁止手抄）

全部数字在
[`artifacts/xval_geom_volcross/results.md`](../artifacts/xval_geom_volcross/results.md)
（表格）与 [`results.json`](../artifacts/xval_geom_volcross/results.json)（机读）中；图在
`artifacts/xval_geom_volcross/figures/`。复现：

```bash
python tools/volume_xval_report.py            # ~45 s, 写入 artifacts/xval_geom_volcross/
python -m pytest tests/test_volume_xval.py -q # 28 passed
```

结论（每项均可在 results.md 的对应小节逐字核对）：

1. **积分器恒等性**：三条三角形路径在**完全相同**的三角形上互相吻合到 1e-13 … 1e-6（引擎 B 的
   1e-6 预算是 STL float32），自研散度求和 vs `trimesh` 最差 7e-15 量级。这条线证明平台积分器
   在闭曲面上没有实现缺陷。
2. **细分亏欠**：管线（光滑几何）比引擎（内接多面体）高，差随加密按实测阶收敛——Wigley 1.98、
   KCS 全尺度 1.82、KCS 模型尺度 1.84，均在 O(h²) 窗口内。最粗档（21 站位）的引擎体积可偏离解析值
   0.9 %，所以**用默认粗网格的体积去对 0.3 % 基准说话是无效的**；这是本单记入 negative-knowledge
   ledger 的负面结论。
3. **外推后与基准闭合**：把引擎体积按测得的阶做 Richardson 外推到 h→0，再与管线求积比较，差异降到
   1e-4 量级（results.md 各行 `extrapolated ...`）；管线求积自身在 Gauss 48×32 与 32×24 之间的漂移
   仅 ~2.4e-6 相对，远小于所有被比较的容差。
4. **WITHHELD 升级**：模型尺度 NURBS 船的体积此前只有"同一片自比"的证据，现由两条独立引擎 + 一条
   独立数值方法与自研求积互证，判定为 `agrees`；但**它仍没有外部基准**，因此证据等级是
   "cross-tool consistent (3 routes)"，不是"与文献吻合"。

## 5. 检测到的实质缺陷（本单暴露的真实问题）

- **`trimesh` 的 `process=True` 会掩盖拓扑缺陷**。默认构造会焊接/修复网格，跨工具互证必须显式
  `process=False`，否则"引擎说水密"是引擎修出来的。模块里已固定为 `process=False` 并自带
  水密/边界边/亏格报送。
- **艉板扇形封盖的星形前提**。用质心扇形封艉时，若横剖面相对于质心非星形（例如舷侧内凹/悬伸），
  扇形会自交并**凭空增加体积**。`closed_solid_from_breadth` 现在用扇形面积 vs 鞋带多边形面积的算术
  比较做检查，不平就抛错（`tests/test_volume_xval.py::test_transom_fan_refuses_a_non_star_shaped_section`）。
- **端部自闭合 vs 端板**。首端半宽为零时曲面自己闭合；仍是满宽（驳船/首板）时必须加端板，否则网格
  有 18 条边界边、任何体积都无意义。两条分支都有测试。

## 6. 注意点与残留缺口

- **Gmsh / OpenCASCADE 在本环境不可安装**（m20574 复核：loongnix 镜像上 `gmsh`、`cadquery-ocp`、
  `pythonocc-core` 均无可用分发，`OCP` 这个名字是无关包，系统亦无 gmsh/FreeCAD 二进制）。锚点写明的"BRep 实体体积"这一路
  因此没有做到：第二台独立引擎改为 `numpy-stl` 的 Mirtich 实现，无三角化路径改为截面 Simpson 求积。
  这是环境能力限制，如实记录，不做静默绕过——BRep 路仍列为未完成工作（`docs/oss_inventory_cad/`）。
- 表中 `pipeline` 对 Wigley 用的是**平台自己的网格路径**（`wigley.build_mesh` + 散度定理），因为它
  就是该船型对外发布的路径；它的收敛轴只能是"引擎 vs 解析基准"，报告里对此有显式说明。
- 非对称艉板船进不了主表：`closed_solid_from_breadth` 只接受对称半宽场。它以"积分器恒等性检验组"
  出现（同网格三路径 + 解析闭式），不回答"几何源对不对"。
- 交叉验证不覆盖**附体**（舵/轴/球鼻）与自由液面修正后的体积；这些仍是各自模块的账。

## 6bis. 轴向对称体：DARPA SUBOFF（外部几何来源）与半宽 loft 的适用边界

前三条船型的半宽场都来自**本仓自己的参数化生成器**（Wigley 闭式、NURBS 片），所以它们只能证明
积分器之间互相一致，证明不了我们的**几何描述**与别人一致——而 JBC 事件的错在几何源，不在积分器。
本单新增 SUBOFF 案例从另一侧补这个洞：半径剖面是公开的 DARPA 四段偏移表（Groves, Huang & Chang
1989，`tensorcad.suboff` 转录，并与 TensorLBM 的独立转录在 `tests/test_suboff.py` 中逐站钉死），
体积是这份**外部几何**的预测值，不是把体积标定出来的恒等式。

### 6bis.1 半宽 loft 装不下回转体（负面知识，已入账）

潜艇横剖面是占 ``[h-r, h+r]`` 一条带的圆，半宽在带的上下方连续多行恒为 0，于是环在中心线上
重复经过同一点，loft 出来不是有效实体：21 站位时 `closed_solid_from_breadth` 直接抛
「not watertight: 0 boundary edges」（0 条边界边 ⇒ 缺陷是非流形顶点，不是洞）。这不是 bug 而是
该 loft 的适用前提（船体剖面必须跨满全局水位网格）。因此新增 `solid_of_revolution`：按余弦站位
布环、`(x, r sin t, h + r cos t)`、两端 `r→0` 的极四边形由 `add_quad` 塌成单三角自闭合，
并校水密与正向体积。测试 `test_half_breadth_loft_cannot_represent_a_body_of_revolution` 把
这个前提钉死，防止将来被“修好”成能返回错体积。

### 6bis.2 结果（由 `tools/volume_xval_report.py` 写入，零手抄）

| 站位 | 方位段 | 引擎A `trimesh` (m³) | 引擎B `numpy-stl` (m³) | 引擎A/B 相对差 | 引擎A vs 偏移表求积 | 亏欠 vs 管线 |
|---:|---:|---:|---:|---:|---:|---:|
| 21 | 42 | 0.691660962 | 0.691660961 | -1.72e-09 | -1.0798 % | +1.0915 % |
| 41 | 82 | 0.697289157 | 0.697289154 | -4.23e-09 | -0.2748 % | +0.2756 % |
| 81 | 162 | 0.698722551 | 0.698722550 | -1.07e-09 | -0.0698 % | +0.0699 % |
| 161 | 322 | 0.699088269 | 0.699088267 | -3.11e-09 | -0.0175 % | +0.0175 % |

判定 `agrees`：两条独立引擎在读**同一批三角形**时最大相对差 4.23e-09（STL float32 预算内），管线（偏移表的 πr² Simpson 求积）比内接多面体高 +0.0175 %，随加密以实测阶 1.99 收敛，Richardson 外推后的引擎体积与管线差 -1.18e-05。拓扑：所有档位水密、亏格 [0]、零面积面 0。

**误差判定（对基准，逐项）**：SUBOFF 没有公开的排水体积可引（公开的是偏移表与湿面积），所以基准取
偏移表的 1e6 采样经线求积 0.699210785 m³。最粗档（21 站位 / 42 段）引擎体积
+1.0915 %，最细档（161 站位 / 322 段）+0.0175 %，
即**细档 <3% 且 <0.3% 验收线成立，粗档不成立**——与 Wigley/KCS 同一条结论：粗网格体积不能对 0.3% 说话。
管线对基准（两条互相独立的求积）相对差 -1.84e-13，是求积分辨率差，不是模型差。

### 6bis.3 这一路**不是**第三条独立积分
对回转体，所谓“光滑路径”就是 π∫r²dx 的经线求积，与管线共用 πr² 剖面面积，因此对 SUBOFF 只有
**两条独立引擎 + 一份外部几何来源**，如实标注在 `results.json` 的 `smooth_route` 字段里，
不冒充三路。

## 7. 验收对照

| 验收项 | 状态 | 证据 |
|---|---|---|
| ① KCS + Wigley + ≥1 WITHHELD 船，≥2 独立引擎对照自研体积 | 达成（3 船 + 1 检验组；引擎 trimesh / numpy-stl；Gmsh/OCP 不可得，见 §6） | `results.json` → `cases[*].same_mesh_rows` |
| ② 三路互证偏差表（自研 / 引擎A / 引擎B / 公开值） | 达成 | `results.md` 每船一张表 |
| ③ >0.3 % 管线级差异定位到积分方法并记 negative-knowledge ledger | 达成（无管线级偏置；粗网格 0.9 % 亏欠定位为 O(h²) 细分，已入账） | `negative_knowledge_ledger` entry `NK-XVAL-GEOM-VOLCROSS-1` |
| ④ 工具选型落 `oss_solutions`（feature_slug=本锚） | 达成 | `oss_solutions` 行 `XVAL-GEOM-VOLCROSS` |
| ⑤ pytest 0 fail | 达成 | `tests/test_volume_xval.py` 33 passed（m20574 增加 5 条回转体用例） |
| ⑥（m20574 增量）外部几何来源船型互证 | 达成 | SUBOFF 案例：`results.json` → `cases[suboff]`，§6bis |
| ⑦（m20574 增量）证据包入库 | 达成 | `artifacts/xval_geom_volcross/{results.json,results.md,figures/}` 随分支提交（此前文档指向的文件未在仓内） |
| ⑧（m22953 增量）能力到达产品面 | 达成 | `python -m tensorcad.cli volume-xval {wigley,mesh}`（退出码 0/1/2）；`fleet_gate.audit_fleet(cross_engine=True)` 默认开；`tests/test_volume_xval_entry.py` 12 passed |

## 8. 产品面：控制台入口与船队门（m22953）

互证能力此前只存在于库函数、测试与报告脚本里——怀疑某个体积的操作员无法在产生第一个数字的
地方要第二个数字。本单把能力接到产品面：

- **控制台入口** `python -m tensorcad.cli volume-xval`（`tensorcad/volume_xval_cli.py`）：
  - `wigley [--stations N --levels N --ladder ...]`：闭式基准船，引擎必须落在解析值上；
  - `mesh FILE.stl [--reference-volume V]`：操作员已有的网格，三条三角形路径互比（无半宽场、
    无解析基准时，互比本身就是定位共享代码偏置的手段）；
  - 退出码 `0` 过门 / `1` 诚实 FAIL（报告照写）/ `2` 拒收输入；`--json` 打印与 `--out` 落盘
    的是同一份 JSON，数字零手抄。
- **船队门集成**（`tensorcad/fleet_gate.audit_fleet`，`cross_engine` 默认开）：每个水密网格被
  独立引擎（trimesh / numpy-stl）重测；两引擎读同一批三角形却超 1e-6 相对差时，该船从
  `VERIFIED` 降为 `WITHHELD` 并记录 spread——两条实现都复现不出来的体积，门不许认证。
  引擎缺装/不可跑时记 `unavailable` + 具名原因，绝不静默放行。
- **包面导出**：`tensorcad` 顶层导出 13 个互证符号（`classify_routes`、`closed_solid_from_breadth`、
  三条体积路由、`topology_report`、`fit_order`、`richardson_limit` 等），不知模块路径的调用方也能找到。
- **`undetermined` 判定**（§2）：无梯差距不再被误标为 `systematic-bias`——`mesh` 案例（无梯）
  正是走这条分支的常态路径。

复现：

```bash
python -m tensorcad.cli volume-xval wigley --json | head
python -m pytest tests/test_volume_xval.py tests/test_volume_xval_entry.py -q
```