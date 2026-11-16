# DTC LCB/Lwl cross-reference gap: defect record and closure plan

Anchor `REPRO-GEOM-FLEET-REFS` — 六标模文献交叉参考值库：每船 ≥2 独立来源。

## 1. 现状（从 `benchmarks/verified/fleet_refs/crossref.json` 读，非手抄）

DTC 在四个命名量（∇ / Cb / LCB / Lwl）上的交叉参考状态：

| 量 | 独立 lineage 数 | 判定 | 文档 |
|---|---|---|---|
| ∇ (m³) | 2 (JSR-2017, DTC-PRIMARY) | CROSS_VERIFIED | dtc_ugent, dtc_strath, dtc_ije |
| Cb | 2 (JSR-2017, DTC-PRIMARY) | CROSS_VERIFIED | dtc_ugent, dtc_strath, dtc_ije |
| LCB (% Lpp) | 0 | **MISSING**（但见 §2 缺陷） | — |
| Lwl (m) | 0 | **MISSING** | — |

`acceptance.json` 判 DTC 在 `lcb_percent` 与 `lwl` 上不满足「≥2 独立来源」。

## 2. 已定位的缺陷：DTC `lcb_percent` 交叉污染

`crossref.json` 的 `report.DTC.quantities.lcb_percent` 含一条**派生** cell：

```
source = "kvlcc2_ugent"   (ship=KVLCC2, JSR-2017)
raw    = "LCG (m) 171.1 169.0 174.059 170.460  (KVLCC2 | heavy ballast | design | light ballast)"
value  = -0.969296 %Lpp
derivation = "LCG 174.059 m, ap datum ... (174.0590 - 177.500) / 355.000 = -0.9693 %Lpp"
```

**缺陷**：`kvlcc2_ugent` 是 JSR-2017 论文 Table 1 的**多船表**（KVLCC2 | DTC heavy | DTC design | DTC light 四列）。其 `design` 列 LCG=174.059 m 是 **DTC 设计吃水** 的纵向浮心，但该 source 被标记为 `ship="KVLCC2"`，于是这条 DTC 值被**错误归因**到 KVLCC2 源键下、并作为 DTC 的 `lcb_percent` 派生 cell 出现。

后果：
- DTC `lcb_percent` 被 `coverage_summary` 计为「有 1 文档」，但该文档是 KVLCC2 源——**不是** DTC 的独立来源；
- 若按此计数，DTC `lcb_percent` 会被误判为「SINGLE_DOCUMENT」而非「MISSING」，掩盖真实缺口；
- 该派生值本身（-0.9693 %Lpp）是**正确**的 DTC 设计吃水 LCB（JSR-2017 表），但归因错误。

**修复方向**：`kvlcc2_ugent` 应识别为多船源（KVLCC2|DTC），其 DTC 列派生 cell 应归到 DTC 源键（如 `dtc_ugent` 同源 JSR-2017），或显式标记为「DTC 派生、非独立」。

## 3. DTC Lwl 缺口

DTC `lwl` 在库中 **0 文档**（`report.DTC.quantities.lwl.entries` 为空）。DTC 主参考（Moctar/Shigunov/Zorn 2012, *Ship Technology Research* 59(3)）给出 **Lwl = 360.0 m**（Lpp=355.5 m）。现有三个 DTC 源（dtc_ugent/dtc_strath/dtc_ije）的 quote 均未覆盖 `lwl` 行。

## 4. 关闭计划（下一轮）

1. **修 §2 归因**：把 `kvlcc2_ugent` 的 DTC 设计吃水 LCG 派生 cell 从 DTC `lcb_percent` 中剔除（或改归 `dtc_ugent` 同源），使 DTC `lcb_percent` 如实回到 `MISSING`；
2. **补 DTC `lwl`**：新增一个含 DTC `Lwl=360.0` 行的独立文档（DTC 主论文或复现其全表的开放获取论文），pin SHA-256 后进 `SOURCES`；
3. **补 DTC `lcb_percent` 第二 lineage**：DTC 主论文 LCB≈1.5%Lpp fwd 与 JSR-2017 派生值（-0.969%Lpp）构成两独立 lineage，使 DTC `lcb_percent` 达 CROSS_VERIFIED。

## 5. 验证

本文件为缺陷记录与计划；`tests/test_fleet_refs_crossref.py` 已绿（13 passed）。下一轮按 §4 执行后重跑 `tools/baseline_fleet/refs_acceptance.py check` 复核 DTC 四量。
