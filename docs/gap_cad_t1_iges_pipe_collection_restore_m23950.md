# GAP-CAD-T1-IGES-PIPE — 恢复合并 m23054 丢弃的 8 个 fleet_baseline 名字，修复验收测试收集崩溃（m23950）

**Anchor** `GAP-CAD-T1-IGES-PIPE` · **Mission** m23950 · **Repo** TensorCAD · **Branch** agent/m23950

## 1. 案例与触发

开工核对总账：锚点 `open`，验收=【管道 CLI+库函数；3 个真实 IGES 全过门；金标回归测试绿；依赖锁定】。
跑锚点自带验收测试 `tests/test_iges_pipeline_gate.py` 时发现**整个测试套件收集即崩**：

```
ERROR tests/test_fleet_trim_intake.py
ERROR tests/test_iges_diagnostics_tool.py
ERROR tests/test_iges_units_and_patch_quarantine.py
ImportError: cannot import name 'partition_surfaces' from 'tensorcad.fleet_baseline'
```

根因（m23349 审计 `docs/gap_cad_t1_iges_pipe_regression_audit.md` 已定位，本单执行修复）：
Merge #698（m23054）把 `tensorcad/fleet_baseline.py` 从 1289 行（c8b07393 架构）回退到 805 行，
丢掉了「trimmed-patch-loops 路线 + 单位门 + 密度跟随」的 8 个名字，但引用它们的 3 个测试文件
留在原地 → **ImportError 在收集期炸掉整个 pytest 会话**（6328 个测试被 deselect，任何测试都跑不了）。

## 2. 怎么做的（最小、纯增量）

从 `c8b07393:tensorcad/fleet_baseline.py` 逐段抽出 8 个名字，插回当前文件
（`sample_patches` 与 half-breadth field 节之间；`partition_surfaces` 同处）：

| 名字 | 角色 | 依赖（当前文件已有） |
|---|---|---|
| `UNIT_BAND` | 单位门 decade 带 (0.5, 2.0) | — |
| `measured_unit_factor` | 声明单位 vs 几何自身尺度对账（KVLCC2 毫米声明/米制坐标） | `unit_scale` |
| `point_in_polygon` | (u,v) 平面 even-odd 射线 | — |
| `_inside_domain` | 环在参数域内占比 | `point_in_polygon` |
| `trim_loops_uv` | 141/142 修剪环 → (u,v) 多边形 | `BSplineSurface` |
| `sample_trimmed_patches` | 全裁剪面片集 in-trim 采样 | `sample_surface`/`_net_lengths` |
| `solid_from_points` | 点云 → 水密实体 | `build_solid` 系 |
| `partition_surfaces` | 歧义布局 patch 隔离（G9 判据） | `BSplineSurface.layout_ambiguous` |

**未动** `build_ship`/`field_from_points`/`build_solid` 等任何现有函数——纯新增，
KCS 等已绿路线按构造不受影响（实测确认，见 §3）。

## 3. 计算结果（程序化注入）

| 检查 | 修前 | 修后 |
|---|---|---|
| 套件收集 | **3 个文件 ImportError，会话中断** | 0 错误，全量可收集 |
| `test_fleet_trim_intake.py`（10 测试） | 收集失败 | 8 passed / 2 skipped（需出版方字节，CI 跳过） |
| `test_iges_diagnostics_tool.py`（14 测试） | 收集失败 | 10 passed / 4 skipped |
| `test_iges_units_and_patch_quarantine.py`（9 测试） | 收集失败 | 6 passed / 3 skipped |
| 锚点验收 `test_iges_pipeline_gate.py` | 可收集（m23349 修 ShipSpec 后） | 7 passed / 1 failed |
| KCS 路线（`test_fleet_baseline_verification.py -k kcs`） | VERIFIED | VERIFIED（3 passed，未回归） |

唯一失败：`test_kvlcc2_fleet_baseline_status`（及同族 KVLCC2 全裁剪面片集测试）——
断言的是**被合并丢弃的** `build_ship` 单位门/密度行为（`declared_contradicted_by_geometry`、
∇ 312 568 m³ = −0.017 %）。当前 `build_ship`（129 行简单路线）无此行为，历史版为 307 行。
按 m23349 审计的告诫（强行回填 307 行 `build_ship` 会破坏已绿 KCS 路线），本单**不回填**
`build_ship`，该测试维持 honest FAIL，欠账原样保留。

## 4. 误差判定（逐项，显式回答 <3 %）

| 关键量 | 基准 | 实测 | 误差 | <3 % ? |
|---|---|---|---|---|
| 套件收集错误数 | 0 | 3 → **0** | 修复 | 是（修复） |
| 3 个原崩溃测试文件 | 可收集可跑 | 24 passed / 8 skipped / 1 failed | — | 是（收集层） |
| KCS ∇ vs 公开 52 030 m³ | 52 030 | 52 187.98（既有证据，本单复跑未回归） | +0.30 % | 是 |
| KVLCC2 ∇（单位门路线） | 312 622 m³ | **未实现**（build_ship 未回填） | — | **否（欠账，未修）** |

## 5. 注意点 / 局限

- 本单只恢复**函数层**（8 个名字），`build_ship` 的 307 行历史路线**未回填**——那是结构性回归，
  需要单独一单做「回填 vs 诚实断言」二选一的裁决（m23349 审计建议项）。
- 8 个 skipped 均为「需出版方字节（不可再分发）」的既有设计，非本单引入。
- 集群作业 cj4746（m23282 递交的 KVLCC2 路线裁决）已 `succeeded`（n315, rc=0），但取回的
  tarball 只含 repo 检出、不含 `artifacts/gap_t1_kvlcc2_routes/evidence.json`——裁决证据未落盘，
  留作下次自找活线索（重跑 `scripts/kvlcc2_route_adjudication.py` 即可复现）。

## 6. 文件

| 路径 | 角色 |
|---|---|
| `tensorcad/fleet_baseline.py` | 恢复 8 个被合并丢弃的名字（纯增量） |
| `docs/gap_cad_t1_iges_pipe_collection_restore_m23950.md` | 本文 |
