# GAP-CAD-T1-IGES-PIPE：进口管道的门控记账与反演精度（m21127）

**仓** TensorCAD · **分支** agent/m21127 · **锚点** GAP-CAD-T1-IGES-PIPE
**验收原文**：管道 CLI+库函数；3 个真实 IGES 全过门；金标回归测试绿；依赖锁定

本单不新增 CAD 读入能力，而是**把这条锚点的三处记账错误改对**——门控误判（把合格样本判成不合格）、
依赖锁测试误判（把内部相对导入当成第三方包）、以及真实文件中一处**被误拒的面**（反演精度不足导致
整条 trim 链退化为无 trim 的 patch-set 路线）。三处都有修前/修后的实测数字。

## 1. 案例

输入是出版方原样交付的 IGES（资产库 `asset://`，非本仓自往返写出）：

| 样本 | 字节 | D 段实体 | 出版值（`fleet_specs.json`） |
|---|---|---|---|
| KCS | 246 164 | `{128: 55}` | Lpp 230.0 m、B 32.2 m、∇ 52 030 m³、Cb 0.6505 |
| JBC（`jbc.igs`） | 1 713 554 | `{126:172, 128:45, 141:19, 143:19, 314:1, 406:1}` | Lpp 280.0 m、B 45.0 m、∇ 178 369.9 m³ |
| KVLCC2 | 6 270 129 | `{102:38, 126:304, 128:19, 142:19, 144:19}` | Lpp 320.0 m、B 58.0 m、∇ 312 622 m³ |

三条样本**都没有 type-124 变换实体**，故"忽略变换矩阵"不是本轮任何失效的原因。

## 2. 怎么做的

工具链不变（`import_cad` → 三角面 → `trimesh` 独立量测 → 与出版值比），本单改动集中三处：

1. `tools/real_iges_gate/run_gate.py`：把评分从 `run_one` 里抽成纯函数
   `score_against_published(measured, spec, tolerance)`，并区分**可比的精度量**与**单侧包含量**。
2. `tests/test_iges_pipeline_gate.py`：依赖锁测试的 AST 判定修正相对导入。
3. `tensorcad/geomquality/iges_trim.py` + `trim_mesh.py`：`BoundedParametricSurface.partials()`
   暴露样条**精确一阶偏导**；`_gauss_newton` 增加 `exact` 开关，仅在**回退栅格**上使用精确切矢，
   种子走法（多数点）仍用便宜的差分。

新增诊断工具 `tools/real_iges_gate/trim_census.py`：逐面报告反演残差与耗时，回答"哪个面被拒、
差多少"。

## 3. 计算结果

### 3.1 门控误判（修前/修后，同一 KCS 输入）

| 项 | 修前 | 修后 | 判据 |
|---|---|---|---|
| KCS `all_within_tolerance` | **false** | **true** | 门控退出码 1 → 0 |
| 记为精度误差的量 | `length_bbox 6.019 %` + `beam_bbox 0.072 %` | 仅 `beam_bbox 0.072 %` | ≤3 % |
| 长度判定 | 双侧相对误差（把包围盒当 Lpp） | 单侧包含 `bbox_x ≥ Lpp` → true | 尺度/单位错才失败 |

`length_bbox` 与 Lpp 差 6.019 % **不是误差**：KCS 的船壳 x 跨度本身就是 −6.000…237.844 m
（艏球与艉外伸都在两柱之外），Lpp 是设计水线上的水静力基准，不是包围盒长度。旧判据把一个
尺度完全正确的样本判成不合格；新判据只用这个数字抓它真正抓得住的失效（mm/m 混用，如 ONRT）。

### 3.2 被拒面的反演精度（JBC boss，逐面实测）

| 面 | 修前残差（差分切矢） | 修后残差（精确切矢） | 容差 | 修前 | 修后 |
|---:|---|---|---|---|---|
| 0 | 5.95e-13 | 9.53e-13 | 3.856e-03 | 接受 | 接受 |
| 1 | 9.94e-13 | 9.11e-13 | 3.856e-03 | 接受 | 接受 |
| **2** | **3.50e-02** | **9.67e-13** | 3.856e-03 | **拒绝** | **接受** |
| 3 | 9.54e-13 | 9.68e-13 | 3.856e-03 | 接受 | 接受 |
| 4 | 3.02e-03 | 3.02e-03 | 3.856e-03 | 接受（临界） | 接受（临界） |
| 5 | 9.82e-13 | 9.57e-13 | 3.856e-03 | 接受 | 接受 |

面 2 修前超容差 9.1×，会被 `tessellate_trimmed_model` 拒绝，整条 trim 路线随即退化成无 trim 的
patch-set（正是 JBC 体积 −18.28 % 那条记录的来路）。精确切矢把它降到 1e-12 量级，面不再被拒。
**代价**：面 0–5 的反演总耗时 19.6 s → 27.1 s（+38 %），因为样条导数网求值比三次普通求值更贵。
故精确切矢只放在回退栅格上，不放在种子走法上。

### 3.3 依赖锁（验收第 4 项）

| 项 | 修前 | 修后 |
|---|---|---|
| `TestDependencyLock::test_iges_core_reader_is_stdlib_only` | **FAIL**：`iges.py imports third-party: ['geomquality.iges_trim']` | **PASS** |
| `TestDependencyLock::test_core_has_no_dependencies` | PASS | PASS |

根因：`iges.py` 里是 `from .geomquality.iges_trim import _tokenize_parameters`（**相对**导入），
而 `ast` 把前导点存在 `node.level` 里、`node.module` 只留 `geomquality.iges_trim`；测试只查
`module.startswith(".")`，于是把本包的一个子模块当成外部第三方包。核心包 `dependencies = []` 的
声明本身从未被违反。

### 3.4 金标回归（验收第 3 项）

`tests/test_iges_gold_regression.py`：**9 passed**（DTMB5415_Hull_exact，26.6 s），改动前后一致。

## 4. 结果图

`artifacts/real_iges_gate/inversion_before_after.png`：左=逐面反演残差（对数轴，含容差线，标出面 2
从 3.50e-02 降到 9.67e-13 的跨容差跳变）；右=逐面耗时前后对比（+38 %）。数字由
`artifacts/real_iges_gate/inversion_before_after.json` 程序化注入，无手抄。

## 5. 注意点 / 局限

1. **一次被否定的优化**（负面结果，已入账）：先按 1 次求值给 49 个栅格种子排序、只从最近 6 个起解，
   期望 5× 加速。实测面 27 从 9.9e-13 崩到 1.4e+01——最近节点与控制点**不在同一吸引域**，
   GN 是局部精化，采样距离不是收敛性的代理。该写法已回退，不做"性能换正确性"。
2. **JBC/KVLCC2 的整文件门控本单仍未取得**：JBC 45 面全量反演在本机 120 s 硬顶下跑不完，
   已递交集群（cj3215，4 核，本单收尾时仍 `queued`），与上一轮 cj2960 同样卡在排队上。
   故"3 个真实 IGES 全过门"这一条**本单只兑现 KCS 一项**，JBC/KVLCC2 记 WITHHELD，不粉饰。
3. 面 4 的残差 3.02e-03 距容差 3.856e-03 只剩 22 % 余量，属临界接受，换子午精度它会翻面。
4. KCS 样本是**水线下开壳**，`displacement` 仍然 WITHHELD（水线/艉切口非平面可封闭），
   门控通过靠的是 beam（0.072 %）+ 单侧包含 + 0 拒绝面，不含任何体积结论。

## 6. 误差判定

| 关键量 | 计算值 | 基准 | 误差 | <3 % ? |
|---|---|---|---|---|
| KCS 型宽 B | 32.2234 m | 32.2 m（NMRI T2015） | **0.072 %** | **是** |
| KCS 排水量 / Cb | WITHHELD | 52 030 m³ / 0.6505 | 不估 | 无基准值可判（开壳） |
| KCS bbox 长 | 243.844 m | Lpp 230.0 m | +6.019 %（**非误差**，单侧包含通过） | 不适用 |
| JBC boss 面 2 反演残差 | 9.67e-13 | 容差 3.856e-03 | 超容差倍数 1（修前 9.1） | **是**（残差口径） |
| JBC 排水量 / KVLCC2 门控 | WITHHELD | — | 未取得 | 见 §5.2 |

## 7. 工件

- `tools/real_iges_gate/run_gate.py`（评分抽出为 `score_against_published`，长度改为单侧包含）
- `tools/real_iges_gate/trim_census.py`（新增：逐面反演残差/耗时普查）
- `tensorcad/geomquality/iges_trim.py`（`BoundedParametricSurface.partials`）
- `tensorcad/geomquality/trim_mesh.py`（`_gauss_newton(..., exact=)`，回退栅格用精确切矢）
- `tests/test_iges_pipeline_gate.py`（相对导入判定修正）
- `tests/test_real_iges_gate_scoring.py`（新增 5 条：包围盒长不记精度误差、短于 Lpp 必须失败、
  beam 仍记误差、体积/Cb 有值即记、WITHHELD 不得记 0 %）
- `artifacts/real_iges_gate/`（summary.json、inversion_before_after.{json,png}）