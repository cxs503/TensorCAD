# GAP-CAD-T1-IGES-PIPE — m24043 路线感知 build_ship 回填（criterion 2 欠账关闭）

**Anchor** `GAP-CAD-T1-IGES-PIPE` · **Mission** m24043 · **Repo** TensorCAD · **Branch** agent/m24043
**日期** 2026-11-16 · **上游基线** c70c6401（origin/main）

## 1. 本单干了什么

按 m23992 状态核对（`docs/gap_cad_t1_iges_pipe_m23992_status.md` §2.1）定位的可执行修复，
把 m21995（commit `625bf3f7`）`fleet_baseline.build_ship` 的「单位门 + 密度跟随 +
fully_trimmed→trimmed_patch_loops + dz 一次细化」逻辑**路线感知**地回填进
Merge #698（m23054）回退后的 `tensorcad/fleet_baseline.py`：

- `fully_trimmed`（19/19 裁剪，KVLCC2）→ 裁剪环路线（unit gate 采用因子重采样 +
  `sample_trimmed_patches` + `solid_from_points` + dz 0.5× 一次细化并报告）；
- 非 fully_trimmed（KCS 0/55 裁剪、JBC 等）→ 保持现路线不动，未触碰已绿行为；
- 同步恢复该路线依赖的报告字段与规格字段（均为 m23054 回退时丢失、m23950 只恢复了
  辅助函数未恢复本体）：
  - `ShipSpec`: `geometry_scale: float | None`（None=按 L_PP 推断并降级为 Cb 门）、
    `waterline_z_file` / `waterline_basis`（声明水线只被佐证不被信任）；
  - `solid_from_points`: `waterline_z_file` 帧移 + `frame_source` 记录；
  - `measure_mesh`: `extent_m`（尺度门用）；
  - `specs_from_json`: geometry_scale/scale_basis/waterline_z_file/waterline_basis 透传
    （此前 JSON 驱动的 fleet 运行把每个档案静默当全尺度）；
  - 报告新增 `source.unit_gate` / `source.route` / `source.trim_sampling` /
    `source.surface_census` / `frame_check` / `volume_gate_basis` /
    `checks.scale_consistent_with_published_lpp` 等；
  - 常量 `FRAME_CHECK_TOL_PERCENT = 2.0`。

改动文件：`tensorcad/fleet_baseline.py`（单文件，+330/−120 行量级，全部为上述恢复项）。

## 2. 验证（本机，资产权限齐全）

| 测试 | 回填前 | 回填后 |
|---|---|---|
| `tests/test_iges_pipeline_gate.py::TestFleetBaselineGate::test_kvlcc2_fleet_baseline_status` | **FAIL**（honest，m23950 起） | **PASS** |
| `tests/test_iges_pipeline_gate.py::TestFleetBaselineGate::test_kcs_fleet_baseline_passes` | PASS | PASS（不回归） |
| `tests/test_iges_pipeline_gate.py::TestFleetBaselineGate::test_jbc_fleet_baseline_status` | PASS | PASS（拒绝语义不变） |
| `tests/test_iges_surface_census.py::test_refusal_names_the_curve_entities_without_blaming_the_hull` | **FAIL**（同根因：回退丢了 126 曲线点名拒绝语） | **PASS** |
| `tests/test_iges_pipeline_gate.py` 全文件 | 1 failed / 1 passed | **2 passed** |
| `tests/test_iges_surface_census.py` + `tests/test_iges_baseline.py` | 1 failed / 其余过 | **全过** |
| `tests/test_fleet_frame_scale_closure.py` | 5 skipped（本地无 workshop 档案） | 5 skipped（同因，非本单范围） |
| `tests/test_fleet_iges_field_integration.py` `test_fleet_trim_intake.py` `test_assetlib.py` `test_fleet_frame_autodetect.py` `test_iges_diagnostics_tool.py` | — | **55 passed, 5 skipped** |
| `tests/test_fleet_baseline.py` | — | **29 passed** |
| `tests/test_fleet_baseline_verification.py` | — | **28 passed** |

KVLCC2 过门数字（测试断言内）：单位门 `declared_contradicted_by_geometry`（声明 mm、
坐标实为 m）；密度重采样 declared 125.0 → adopted 0.125 file units；
体积 312 568 m³ vs 公开 312 622 m³ = **−0.017 %**（<1 % 门）；watertight + edge_manifold；
verdict **VERIFIED**。

金标回归（criterion 3）：`tests/test_iges_gold_regression.py` 走 `dtmb5415`/`gold_loft`
路径，不经过 `fleet_baseline`，本单改动不触及；其绿状态维持 m23992 核对结论
（金标 STL 在控制面 `/root/ship-performance-platform/backend/data/geometry/ships/DTMB5415_Hull_exact.stl`）。

## 3. 未做 / 残留

- **KVLCC2/ONRT trim_mesh 路线墙钟**（m23992 §2.2）：142 字段位解析修复与 NURBS 导数
  memoize 已合并，但 120 s 门内的墙钟实测仍未完成（需在有资产权限的节点跑
  `test_real_iges_gate.py` 计时）。criterion 2 的「3 个真实 IGES 全过门」中 KVLCC2 的
  fleet_baseline 路线已绿；trim_mesh 路线的墙钟判定留作下次自找活。
- `test_fleet_frame_scale_closure.py` 5 项 skip：依赖 `artifacts/baseline_fleet/iges/`
  下的 JBC/ONRT workshop 档案（不随仓分发），非代码缺陷。

## 4. 给下次自找活的建议线索

1. **首选**：在有资产权限的节点实测 KVLCC2/ONRT trim_mesh 墙钟（`test_real_iges_gate.py`
   + `test_iges_pipeline_gate.py::TestTrimMeshGate` 计时），判定 criterion 2 是否可整体关闭。
2. 建议新总账条目（若 trim_mesh 墙钟也过）：`GAP-CAD-T1-IGES-PIPE-TRIM-CLOCK`，
   验收 = KVLCC2/ONRT trim_mesh 路线 <120 s 且 G1/G2/G5 过门，所属域 CAD 进口。
