# GAP-CAD-T1-BASELINE-FLEET 几何缺口分析

**锚点**: `GAP-CAD-T1-BASELINE-FLEET`  
**日期**: 2026-11-15  
**仓库**: TensorCAD

---

## 现状总结

| 验收标准 | 状态 | 说明 |
|----------|------|------|
| 6 船全部水密+流形 | ✅ 通过 | 全部 watertight=True, edge_manifold=True, Euler=2 |
| 体积 vs 公开排水量误差 <1% | ✅ 通过 | 最差 KVLCC2 @ 0.54% |
| 每船 manifest（源URL/license/sha/验证结论） | ✅ 通过 | FleetManifest + SourceRegistry 完整 |
| 占位 STL 替换方案 | ✅ 通过 | 304 三角形占位符被检测并拒绝 |

## 核心缺口：型族代理 ≠ 真实 IGES 几何

当前 6 船均使用**参数化型族代理**（form-family surrogate），而非官方 workshop 发布的 IGES 几何：

| 船名 | 当前路线 | 真实几何来源 | IGES 文件状态 |
|------|----------|-------------|--------------|
| KCS | series60_form | SIMMAN 2014 | ✅ 已获取 (SHA: b50adace...) |
| KVLCC2 | parametric_sections | SIMMAN 2014 | ✅ 已获取 (SHA: dc032647...) |
| ONRT | parametric_sections | SIMMAN 2020 / Gothenburg 2010 | ✅ 已获取 (SHA: bd747eff...) |
| JBC | parametric_sections | Tokyo 2015 | ✅ 已获取 (SHA: e19800ee...) |
| DTC | series60_form | El Moctar et al. (2012) | ❌ URL 404，需新链接 |
| Series60 | series60_form | DTMB Report 1710 | N/A（使用型值表） |

### 型族代理的局限

1. **无球鼻艏**（bulbous bow）：KCS、KVLCC2、JBC 均有显著球鼻艏，型族无法表示
2. **无艉封板细节**（transom stern）：DTC、KCS 有艉封板，型族使用对称艉部
3. **无局部光顺**：真实船型的局部曲率变化无法用全局截面律捕捉
4. **体积系统偏差**：型族体积基于 Cb×L×B×T，与真实 IGES 几何的体积可能有 2-5% 偏差

### IGES 路线就绪状态

`iges_field` 路线（IGES → 半宽场 → 水密网格）已实现并集成：

| 组件 | 状态 | 文件 |
|------|------|------|
| IGES BSplineSurface 解析 | ✅ | `tensorcad/iges.py` |
| 表面采样 + 半宽场构建 | ✅ | `tensorcad/fleet_baseline.py::field_from_points` |
| 水线切片 → 水密网格 | ✅ | `tensorcad/fleet_baseline.py::build_solid` |
| 路线调度 | ✅ | `tensorcad/fleet/build.py::build_ship_geometry` |
| 门控集成 | ✅ | `tensorcad/fleet/gates.py` |
| 测试覆盖 | ✅ | `tests/test_fleet_iges_field_integration.py` (12 tests) |

## 下一步行动

### 立即可做（无需网络）
1. ✅ 完成（本单）：验证报告 + IGES 路线集成测试 + 缺口分析文档
2. 更新 acquisition_plan.json 中 DTC 的 URL（需 custodian 提供）

### 需要网络/文件访问
1. 将 KCS/KVLCC2/ONRT/JBC 的 IGES 文件复制到本环境
2. 更新 fleet_manifest.json 中对应船的 `geometry_route = "iges_field"` 和 `geometry_path`
3. 运行 `run_fleet` 验证 IGES 几何的体积误差
4. 更新 source_registry.json 中对应船的 `verified = true`

### 建议新条目（供下次自找活）
- **slug**: `GAP-CAD-T2-IGES-FLEET-INTEGRATION`
- **标题**: 将 4 船（KCS/KVLCC2/ONRT/JBC）的 IGES 几何集成到基准舰队
- **量化验收**: 4 船 `geometry_route = "iges_field"`，体积 vs 官方排水量误差 <3%，manifest 中 `verified = true`
- **所属域**: TensorCAD / fleet
