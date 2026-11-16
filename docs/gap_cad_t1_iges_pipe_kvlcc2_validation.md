# GAP-CAD-T1-IGES-PIPE — KVLCC2 物理验证收口（验证欠账）

**Anchor**: GAP-CAD-T1-IGES-PIPE · **Mission**: m22886 · **Repo**: TensorCAD · **Branch**: agent/m22886

## 1. 为什么这一轮存在（验证欠账）

锚点验收条款「3 个真实 IGES 全过门」中，KCS 已物理验证绿（`test_kcs_workshop_geometry_matches_published_particulars`：
体积 52 030 m³ / 湿表面积 9 424 m² / Cb 0.651 全部 <3%）。但 **KVLCC2**（第二个真实 IGES）从未被物理验证：

- round6（m22618）记录 KVLCC2 为 `NOT_RUN`（系统 load > 100，导入超时）；
- KVLCC2 资产 MANIFEST 自述「ref 312622m3 未接 runner；waterline-slice bug open」。

根因不是形状缺陷而是**坐标系缺陷**：KVLCC2 导出把 `z=0` 放在设计水线、`z` 向下增长
（实测 `0 ... -28.067 m`），`x=0` 在船首（船尾 `x ~= -328 m`）。喂给「保留 `0<=z<=draft`」的
路由会**什么都不留**——水线切片消失。这是 `tensorcad.fleet.frames` 已测量并修复的帧缺陷
（`normalise` 纯函数，`test_fleet_frames.py` 已绿）。

## 2. 怎么做的

对 KVLCC2 真实 IGES 文件（`kvlcc2.igs`，6.3 MB，38 面片）：

1. `import_cad(path, weld_tolerance=1e-3)` 导入（G1/G2 全读 + 修剪在面上）；
2. `normalise(vertices)` 应用声明式帧变换（`z_axis=auto` 翻 z、`x_origin=auto` 平移 x），
   记录 `z_flip` / `x_shift` 证据；
3. `hull_characteristics(pts, faces, waterline=20.8, stations=61, levels=48)` 测体积/Cb/S；
4. `validate_against_particulars(got, PARTICULARS["KVLCC2"], tolerance_percent=3.0)` 判 3% 门。

## 3. 验收对照

| 验收项 | 要求 | 状态 |
|--------|------|------|
| 管道 CLI+库函数 | 可用 | ✅（既有） |
| 3 个真实 IGES 全过门 | KCS ✅ / **KVLCC2 本轮** / ONRT（scale 待定） | KCS 已绿；KVLCC2 本轮补 |
| 金标回归测试绿 | DTMB5415 全 7 门 | ✅（既有，4 测试绿） |
| 依赖锁定 | `requirements-iges-pipe.lock` | ✅（既有） |

## 4. 新增测试

`tests/test_hull_surface_validation.py::test_kvlcc2_workshop_geometry_matches_published_particulars`
——断言 KVLCC2 帧归一化后体积/Cb/S 对 312 622 m³ / 27 194 m² / Cb 0.8098 全部 <3%，
且 `z_flip=True`、`x_shift<0`（帧证据）、`half_breadth=True`、`watertight=False`（诚实拓扑）。

## 5. 注意点 / 局限

- KVLCC2 导入重（6.3 MB / 38 面片，>120 s），本机工具超时，已递交集群（cj4465/cj4466）。
- ONRT 为 model-scale 坐标（MANIFEST：「scale factor hunt open」），需先定标尺才能物理验证，
  本单未做（留作下次自找活线索）。
- 帧归一化是纯函数（点云），不依赖 6 MB 文件即可测——`test_fleet_frames.py` 已覆盖。
