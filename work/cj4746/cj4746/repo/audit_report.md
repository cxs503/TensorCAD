# 交差审计报告：m17249 GAP-CAD-T1-BASELINE-FLEET

**审计任务**：m17269
**审计对象**：m17249（GAP-CAD-T1-BASELINE-FLEET·标准船型权威几何集）
**审计日期**：2026-11-13
**审计人**：board-audit

---

## 执行摘要

m17249 交付了 TensorCAD 的 fleet_gate 工具（检测占位网格缺陷），7/7 单元测试通过，4/6 权威 IGES 源已获取并记录出处。但：
- **基准用例未注册**：`tensor_benchmark_cases` 中无 `GAP-CAD-T1-BASELINE-FLEET` 条目
- **6 船几何未转换**：0/6 网格获得 VERIFIED（TensorCAD 无 CAD 内核）
- **主目标 honest-FAIL**：体积误差 <1% 的验收线未达成

**判级：L2**（独立功能实现）

---

## ① 取证核验

### 三件套

| 项目 | 状态 | 证据 |
|---|---|---|
| run_artifacts | ✅ | `artifacts/baseline_fleet/{sources.json, published.json, placeholder_audit.json, iges_inventory.json, figs/}` |
| benchmark ref vs got | ❌ | `SELECT ... FROM tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-T1-BASELINE-FLEET'` → 0 行 |
| repro_cmd | ✅ | `pytest tests/test_fleet_gate.py -v` → 7/7 pass；`audit.py audit` → 正确识别 6 PLACEHOLDER |

### 报告六段

| 段落 | 状态 |
|---|---|
| ① 案例介绍 | ✅ |
| ② 方法 | ✅ |
| ③ 结果 | ✅ |
| ④ 结果图 | ✅ 2 张 |
| ⑤ 注意点 | ✅ |
| ⑥ 误差判定 | ✅ honest-FAIL |

---

## ② 代码审查

### 实现文件
- `tensorcad/fleet_gate.py`：MeshAudit, judge_mesh, audit_fleet, ShipVerdict（~200 行）
- `tools/baseline_fleet/audit.py`：CLI 包装
- `tools/baseline_fleet/fetch_sources.py`：HTTPS 出处获取
- `tools/baseline_fleet/make_figures.py`：图生成

### 测试
- `tests/test_fleet_gate.py`：7 测试全部通过
  - test_mesh_helpers_agree_on_synthetic_cube
  - test_binary_and_ascii_stl_round_trip
  - test_placeholder_box_is_rejected_by_facet_count
  - test_duplicate_bytes_across_two_ships_is_placeholder
  - test_watertight_form_passes_volume_gate_and_open_form_withholds
  - test_provenance_validator_requires_official_host_licence_and_sha
  - test_audit_fleet_writes_json_report

### 复现验证
```
$ pytest tests/test_fleet_gate.py -v
7 passed, 1 warning in 1.45s

$ audit.py audit /root/ship-m2056/backend/data/geometry/ships --published published.json
DTC        PLACEHOLDER  tris=304
DTMB5415   WITHHELD     tris=127624
JBC        PLACEHOLDER  tris=304
KCS        PLACEHOLDER  tris=304
KVLCC1     PLACEHOLDER  tris=304 (byte-identical to KVLCC2)
KVLCC2     PLACEHOLDER  tris=304 (byte-identical to KVLCC1)
ONRT       PLACEHOLDER  tris=304
pass: False
```

---

## ③ 判级

| 级别 | 标准 | 判定 | 理由 |
|---|---|---|---|
| L0 | 文档成果 | ✅ | 有报告、README |
| L1 | 数学/单元验证 | ✅ | 7 个单元测试 |
| L2 | 独立功能实现 | ✅ | fleet_gate.py：输入→计算→输出 |
| L3 | 系统集成 | ❌ | 未接入主求解器/平台 API |
| L4 | 基准验证 | ❌ | 无基准用例注册、无完整算例 |
| L5 | 科研工程成果 | ❌ | |

**最终判级：L2**

---

## ④ verdict 落表

```sql
INSERT INTO tensor_maturity_verdicts (repo, scope, mission_id, level, evidence, verdict_by)
VALUES ('TensorCAD', 'm17249 GAP-CAD-T1-BASELINE-FLEET', 17269, 2, '<JSON>', 'board-audit');
```

已执行，verdict id=504。

---

## ⑤ 关键发现

1. **基准用例缺失**：验收标准未注册到 `tensor_benchmark_cases`，机器门无法自动裁定
2. **IGES 不支持**：TensorCAD 声明 IGES 为 NOT_SUPPORTED，无法执行几何转换
3. **占位网格检测正确**：6 个 304 三角面占位网格被正确识别
4. **KVLCC1/KVLCC2 字节级重复**：同一文件两个名字，被正确检测
5. **4/6 源获取**：DTC 和 Series60 未找到（404/410 Gone）

---

## ⑥ 建议

1. 注册 `GAP-CAD-T1-BASELINE-FLEET` 基准用例到 `tensor_benchmark_cases`
2. 添加 CAD 内核（OpenCASCADE）以支持 IGES→mesh 转换
3. 将 fleet_gate 集成到 ship-platform CI
