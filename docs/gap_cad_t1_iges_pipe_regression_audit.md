# GAP-CAD-T1-IGES-PIPE 回归审计：ShipSpec 签名 + KVLCC2 单位门/密度修复丢失（m23349）

## 1. 案例与触发
锚点 `GAP-CAD-T1-IGES-PIPE`（IGES/STEP→水密网格进口管道，金标回归=DTMB5415_Hull_exact）。
本单核对 feature_ledger 时该锚为 `open`，故直接认领做实质增量。开工先跑锚点自带的验收测试
`tests/test_iges_pipeline_gate.py`（"3 个真实 IGES 全过门"），发现**收集即失败**：

```
ERROR collecting tests/test_iges_pipeline_gate.py
TypeError: ShipSpec.__init__() got an unexpected keyword argument 'geometry_scale'
```

即锚点自己的验收测试连收集都过不去——这是"验证欠账"类（已合并功能对验收基准复核）的典型。

## 2. 缺陷 A（已修复）：`ShipSpec` 缺 `geometry_scale`/`scale_basis` 字段
`tensorcad/fleet_baseline.py::ShipSpec`（dataclass）只有 `name/draft/length_pp/beam/…/notes`，
**没有** `geometry_scale`/`scale_basis`。但两处调用方都传这两个关键字：
- 锚点验收测试 `tests/test_iges_pipeline_gate.py`（KCS/KVLCC2/JBC 三个 `ShipSpec(...)` 均传 `geometry_scale=…`）；
- `tensorcad/fleet/official_geometry.py::OfficialShip.spec()`（第 302 行 `spec = ship.spec()` 构造 `ShipSpec` 时传 `geometry_scale=self.geometry_scale, scale_basis=self.scale_basis`）。

后果：锚点验收测试**无法收集**（TypeError），`OfficialShip.spec()` 是**潜在运行时崩溃**（任何走
official 几何路线的船都会炸）。`scale_calibration` 模块文档明确 `geometry_scale` 是"模型→全尺度因子"
（KCS/KVLCC2 全尺度=1.0，JBC/ONRT 是模型），语义已定，只是 dataclass 没跟上。

**修复**（最小）：在 `ShipSpec` 加两个字段（放在 `notes` 前，与 `official_geometry.py` 字段序一致）：
```python
geometry_scale: float = 1.0
scale_basis: str = "geometry file is full scale (publisher archive)"
```
`to_dict()` 自动带上（`dict(self.__dict__)`），报告 `spec` 节随之携带该 provenance。

**验证**：`test_iges_pipeline_gate.py` 从"收集失败"变为"8 个测试收集成功"；
`OfficialShip.spec()` 不再崩溃；KCS 路线 VERIFIED、TrimMesh 路线、依赖锁定 5 项全绿。

## 3. 缺陷 B（未修复，诚实记录）：KVLCC2 单位门/密度修复在合并中丢失
`test_kvlcc2_fleet_baseline_status` 断言**修后**行为（`unit_gate.verdict == "declared_contradicted_by_geometry"`、
`density_resampled_for_adopted_factor`、∇ 312 568 m³ = −0.017 %、VERIFIED），但当前 `build_ship` **没有**这些：
```
tensorcad/fleet_baseline.py:489  factor = unit_scale(model.unit_flag)   # KVLCC2 → 0.001（毫米）
```
KVLCC2 全局段声明毫米（flag 2）但坐标是米制，`unit_scale(2)=0.001` 把 320 m 船缩成 0.32 m，
半宽场 `filled_nodes=1/209`（空场），`build_solid` 抛 `ValueError: field produced fewer than two usable waterline slices`。

**根因**：`git log -S "measured_unit_factor"` 显示该逻辑在 `b352a8c2`(m21995)/`c8b07393`(m21706) 实现，
但 `d0da2b85`(Merge #698, m23054) 的合并结果已**不含** `measured_unit_factor`——合并冲突把整个
"trimmed-patch-loops 路线 + 单位门 + 密度跟随"架构**丢弃**，当前 `build_ship` 回退到简单的
"untrimmed 采样 + 声明单位"版本。`docs/gap_cad_t1_iges_pipe_kvlcc2_density.md` 与测试都描述的是
被丢掉的架构。

**为何不在此单修复**：KVLCC2 的修复不是加一个因子——它需要 `sample_trimmed_patches`（19/19 全裁剪面，
untrimmed 域是 trim 域的 6× 且低于龙骨 7 m）、`solid_from_points`、`partition_surfaces`、`surface_census`
等一整条被合并丢弃的路线。这是**结构性回归**，超出"最小修复"范围，且当前代码已朝不同方向演化，
强行回填会破坏 KCS 等已绿路线。**honest FAIL 不覆盖**：本单如实记录该回归，不把 KVLCC2 测试改成
"假装通过"。

**JBC 同理**：`test_jbc_fleet_baseline_status` 断言 `ValueError("fewer than two usable waterline slices")`
（修后拒绝语义），当前 `build_ship` 因空场在 `build_solid` 抛同名错误——但那是**空场**（采样没进船体），
不是"帧/单位"拒绝，语义不同。该测试同样描述被丢掉的架构。

## 4. 误差判定（本单可验证量）
| 关键量 | 基准 | 实测 | 误差 | < 3 % ? |
|---|---|---|---|---|
| 锚点验收测试收集 | 应收集 | 8 个收集成功（修前 TypeError） | — | 是（修复） |
| `OfficialShip.spec()` | 应可构造 | 构造成功（修前崩溃） | — | 是（修复） |
| KCS 路线 | VERIFIED | VERIFIED（∇ +0.30 %） | — | 是 |
| DTMB5415 金标回归 | 水密+单壳 | 9 项全绿（test_iges_gold_regression） | — | 是 |
| KVLCC2 单位门/密度 | 修后 −0.017 % | **未实现**（空场） | — | **否（回归，未修）** |

## 5. 注意点与局限
- 本单只做"最小、可验证"的 `ShipSpec` 缺陷修复 + 诚实回归审计；KVLCC2/JBC 结构性回归**未修**，
  已如实记录，不粉饰。
- `test_gap_cad_t1_iges_pipe.py` 在本机 120 s 墙钟内跑不完（含 DTMB5415 loft + 真实 IGES 重采样），
  记 not_run；金标回归 `test_iges_gold_regression.py` 已单独全绿。
- 建议后续：派单恢复 `c8b07393` 的 trimmed-patch-loops 路线（KVLCC2 单位门+密度），
  或把 `test_kvlcc2_fleet_baseline_status`/`test_jbc_fleet_baseline_status` 改为对当前简单路线的
  诚实断言（WITHHELD + 空场原因），二选一，不得两者都留。

## 6. 文件
`tensorcad/fleet_baseline.py`（ShipSpec 加字段）· `docs/gap_cad_t1_iges_pipe_regression_audit.md`（本审计）。
