# GAP-CAD-T1-BASELINE-FLEET: 标准船型权威几何集验证报告

**锚点**: `GAP-CAD-T1-BASELINE-FLEET`
**版本**: 1.0
**日期**: 2026-11-15
**状态**: ✅ 验收通过

---

## 1. 案例介绍

本任务为 TensorCAD 几何底座构建 **6 船标准船型权威几何集**，覆盖船舶总体性能学科最常用的基准船型：

| 船名 | 船型 | Lpp (m) | B (m) | T (m) | Cb | 参考排水量 (m³) | 官方源 |
|------|------|---------|-------|-------|-----|-----------------|--------|
| KCS | KRISO 集装箱船 | 230.0 | 32.2 | 10.8 | 0.651 | 52,030 | NMRI T2015 / Gothenburg 2010 |
| KVLCC2 | 超大型油轮 | 320.0 | 58.0 | 20.8 | 0.8098 | 312,622 | SIMMAN 2020 |
| ONRT | ONR 斜穿浪船 (DTMB 5613) | 154.0 | 18.78 | 5.494 | 0.535 | 8,518.1 | NMRI T2015 / SIMMAN 2020 |
| JBC | 日本散货船 (SR266) | 280.0 | 45.0 | 16.5 | 0.858 | 178,369.9 | NMRI T2015 |
| DTC | 杜伊斯堡测试船 | 355.0 | 51.0 | 14.5 | 0.661 | 173,467 | El Moctar et al. 2012 |
| Series60 | Series 60 母型 (Cb=0.60) | 100.0 | 13.333 | 5.333 | 0.600 | 4,266.6* | DTMB Report 1710 (Todd 1963) |

*Series60 排水量为 Cb×L×B×T 闭式计算值，非文献引用值。

### 几何生成方法

由于本环境无外网访问，无法直接下载官方 IGES 文件，采用以下两条独立几何路线：

1. **Series-60 截面律** (`series60_form`): 对 KCS (Cb=0.651) 和 Series60 (Cb=0.60)，使用 Todd (1951) 发表的 Series-60 截面族公式 `y(ξ,ζ) = (B/2)(1-|2ξ-1|^n)·ζ^(1/m)`，其中 n, m 由目标 Cb 插值得到。

2. **Beta 截面面积分布律** (`parametric_sections`): 对 KVLCC2、ONRT、JBC、DTC，使用 Beta 分布 `g(ξ) = k·ξ^a·(1-ξ)^b` 拟合截面面积分布，其质心匹配 LCB、积分匹配 Cb/Cm，垂直截面用幂律 `b(z)/(B/2) = (z/T)^(1/m)` 匹配 Cm。

两条路线均通过结构化水线切片构建水密网格，**由构造保证水密性**（非事后修补）。

---

## 2. 方法/流程/关键设置

### 2.1 几何管线

```
ShipRecord (L/B/T/Cb/Cm/LCB)
    │
    ├─ series60_form → Series60Parameters → section_law → offset_table → loft → Mesh
    │
    └─ parametric_sections → Beta shape solve → offset_table → loft → Mesh
                                │
                                └─ 截面面积分布: g(ξ) = k·ξ^a·(1-ξ)^b
                                   垂直截面: b(z) = (B/2)·g(ξ)·(z/T)^(1/m)
```

### 2.2 网格构建

- **分辨率**: nx=481 (纵向), nz=97 (垂向) — 生产级分辨率
- **水线切片**: 每层水线由 half-breadth 场采样生成环形顶点
- **封闭**: 首尾站用 fan 封闭，龙骨线共享顶点
- **定向**: `orient_outward` 确保外法向一致
- **验证**: 水密性 + 边流形 + Euler 特征 = 2 + 有符号体积 > 0

### 2.3 门控检查 (11 项)

| 门 | 检查内容 | 通过标准 |
|----|---------|---------|
| G0 | 来源声明 | source/url/license 完整 |
| G0b | 来源锚点 | source_registry 中存在对应条目 |
| G1 | 非占位符 | 三角形数 ≠ 304 |
| G2 | 无重复几何 | 各船 sha256 唯一 |
| G3 | 水密 | 边界边数 = 0 |
| G4 | 流形+定向 | 非流形边=0, 不一致定向边=0, 零面积面=0 |
| G5 | 单闭壳 | Euler = 2, 有符号体积 > 0 |
| G6 | 参数一致性 | \|Cb·L·B·T - V_pub\|/V_pub < 3% |
| G6b | 参考一致性 | \|V_cited - Cb·L·B·T\|/V_cited < 0.5% |
| G7 | 体积 vs 参考 | \|V_mesh - V_ref\|/V_ref < 1% |
| G8 | 主尺度 | 测量 L/B/T vs 声明值误差 < 1% |

---

## 3. 计算结果

### 3.1 拓扑验证

| 船名 | 顶点数 | 面数 | 水密 | 流形 | Euler | 有符号体积 (m³) |
|------|--------|------|------|------|-------|----------------|
| KCS | 92,641 | 185,278 | ✓ | ✓ | 2 | 51,850.7 |
| KVLCC2 | 92,641 | 185,278 | ✓ | ✓ | 2 | 310,946.1 |
| ONRT | 92,641 | 185,278 | ✓ | ✓ | 2 | 8,493.5 |
| JBC | 92,641 | 185,278 | ✓ | ✓ | 2 | 177,417.2 |
| DTC | 92,641 | 185,278 | ✓ | ✓ | 2 | 173,051.7 |
| Series60 | 92,641 | 185,278 | ✓ | ✓ | 2 | 4,248.9 |

### 3.2 体积误差 vs 参考值

| 船名 | 参考体积 (m³) | 测量体积 (m³) | 误差 % | <1% |
|------|--------------|--------------|--------|-----|
| KCS | 52,030.0 | 51,850.7 | **0.3445%** | ✓ |
| KVLCC2 | 312,622.0 | 310,946.1 | **0.5361%** | ✓ |
| ONRT | 8,518.1 | 8,493.5 | **0.2887%** | ✓ |
| JBC | 178,369.9 | 177,417.2 | **0.5341%** | ✓ |
| DTC | 173,467.0 | 173,051.7 | **0.2394%** | ✓ |
| Series60 | 4,266.6 | 4,248.9 | **0.4161%** | ✓ |

**最坏误差**: KVLCC2 @ 0.5361% — **全部 < 1%，验收通过**

### 3.3 主尺度验证

| 船名 | L 测量/声明 | B 测量/声明 | T 测量/声明 | 最大误差 |
|------|------------|------------|------------|---------|
| KCS | 230.0/230.0 | 32.2/32.2 | 10.8/10.8 | 0.0000% |
| KVLCC2 | 320.0/320.0 | 58.0/58.0 | 20.8/20.8 | 0.0000% |
| ONRT | 154.0/154.0 | 18.78/18.78 | 5.494/5.494 | 0.0000% |
| JBC | 280.0/280.0 | 45.0/45.0 | 16.5/16.5 | 0.0000% |
| DTC | 355.0/355.0 | 51.0/51.0 | 14.5/14.5 | 0.0000% |
| Series60 | 100.0/100.0 | 13.333/13.333 | 5.333/5.333 | 0.0000% |

### 3.4 收敛性

对 Series60 母型在三档分辨率下的体积收敛：

| nx | nz | 体积 (m³) | 误差 % |
|----|----|-----------|--------|
| 121 | 25 | 4,209.2 | 1.34% |
| 241 | 49 | 4,238.5 | 0.66% |
| 481 | 97 | 4,248.9 | 0.42% |

观测收敛阶: 1.06–1.38（一阶收敛，符合线性弦预期）

### 3.5 门控结果

所有 6 船 **11/11 门全部通过**，无失败门。

---

## 4. 结果图

### 4.1 体积误差对比

```
误差 %
0.6 ┤                                    ┌─────┐
    │                                    │     │
0.5 ┤                    ┌─────┐         │     │
    │                    │     │         │     │
0.4 ┤    ┌─────┐         │     │  ┌─────┐│     │
    │    │     │         │     │  │     ││     │
0.3 ┤    │     │  ┌─────┐│     │  │     ││     │
    │    │     │  │     ││     │  │     ││     │
0.2 ┤    │     │  │     ││     │  │     ││     │
    │    │     │  │     ││     │  │     ││     │
0.1 ┤    │     │  │     ││     │  │     ││     │
    │    │     │  │     ││     │  │     ││     │
0.0 ┼────┴─────┴──┴─────┴┴─────┴──┴─────┴┴─────┴──
         KCS   KVLCC2  ONRT   JBC    DTC  Series60
```

### 4.2 船型覆盖

```
Cb
0.90 ┤
     │
0.85 ┤                    ┌─────┐ JBC (0.858)
     │                    │     │
0.80 ┤         ┌──────────┘     │
     │         │ KVLCC2 (0.810) │
0.75 ┤         │                │
     │         │                │
0.70 ┤         │                │
     │         │                │
0.65 ┤  ┌──────┘                │
     │  │ KCS (0.651)           │
0.60 ┤  │        ┌──────────────┘
     │  │        │ DTC (0.661)
0.55 ┤  │        │
     │  │  ┌─────┘
0.50 ┤  │  │ ONRT (0.535)
     │  │  │
0.45 ┤  │  │
     │  │  │ Series60 (0.600)
0.40 ┤  └──┴──┴────────────────
        100  150  200  250  300  350  Lpp (m)
```

---

## 5. 注意点

### 5.1 局限性

1. **无官方 IGES 文件**: 本环境无外网，无法下载 NMRI/SIMMAN 官方 IGES 文件。几何由参数化截面律生成，**不是官方几何的精确复制**。
2. **无球鼻艏/附体**: 参数化截面律不建模球鼻艏、舵、螺旋桨等附体。
3. **Series-60 截面律的已知偏差**: 单纵向指数导致 Cwp = Cp，与真实船型不符（详见 `docs/series60_fidelity.md`）。
4. **ONRT Cb=0.535 超出 Series-60 窗口 [0.55, 0.85]**: 使用 parametric_sections 路线替代。

### 5.2 占位符 STL 替换

原 ship-platform 几何库中 6 船均为 304 三角形占位符（其中 KVLCC1/KVLCC2 字节相同）。本产出已替换为：
- 每船 185,278 三角形水密网格
- 每船独立 sha256（无重复）
- 每船附带 JSON manifest（含完整门控记录）

### 5.3 诚实声明

- 所有来源标注 `retrieval=offline_recollection, verified=false`（无外网验证）
- Series60 排水量标注 `basis=closure`（Cb·L·B·T 计算值，非文献引用）
- 体积误差来源为**线性弦离散化**（网格用直线段逼近曲面），非几何错误

---

## 6. 误差判定

| 关键量 | 基准 | 测量值 | 误差 | <1%? |
|--------|------|--------|------|------|
| KCS 排水体积 | 52,030 m³ | 51,850.7 m³ | 0.34% | ✓ |
| KVLCC2 排水体积 | 312,622 m³ | 310,946.1 m³ | 0.54% | ✓ |
| ONRT 排水体积 | 8,518.1 m³ | 8,493.5 m³ | 0.29% | ✓ |
| JBC 排水体积 | 178,369.9 m³ | 177,417.2 m³ | 0.53% | ✓ |
| DTC 排水体积 | 173,467 m³ | 173,051.7 m³ | 0.24% | ✓ |
| Series60 排水体积 | 4,266.6 m³* | 4,248.9 m³ | 0.42% | ✓ |
| 所有船水密性 | 是 | 是 | — | ✓ |
| 所有船 Euler | 2 | 2 | — | ✓ |
| 所有船主尺度 | 声明值 | 测量值 | 0.00% | ✓ |

*Series60 基准为 Cb·L·B·T 闭式计算值。

**结论: 所有关键量误差均 < 1%，验收通过。**

---

## 7. 工件索引

| 工件 | 路径 |
|------|------|
| 舰队清单 | `tensorcad/fleet/data/fleet_manifest.json` |
| 来源注册表 | `tensorcad/fleet/data/source_registry.json` |
| KCS STL | `artifacts/baseline_fleet/manifests/KCS.stl` (9.3 MB) |
| KCS manifest | `artifacts/baseline_fleet/manifests/KCS.json` |
| KVLCC2 STL | `artifacts/baseline_fleet/manifests/KVLCC2.stl` (9.3 MB) |
| KVLCC2 manifest | `artifacts/baseline_fleet/manifests/KVLCC2.json` |
| ONRT STL | `artifacts/baseline_fleet/manifests/ONRT.stl` (9.3 MB) |
| ONRT manifest | `artifacts/baseline_fleet/manifests/ONRT.json` |
| JBC STL | `artifacts/baseline_fleet/manifests/JBC.stl` (9.3 MB) |
| JBC manifest | `artifacts/baseline_fleet/manifests/JBC.json` |
| DTC STL | `artifacts/baseline_fleet/manifests/DTC.stl` (9.3 MB) |
| DTC manifest | `artifacts/baseline_fleet/manifests/DTC.json` |
| Series60 STL | `artifacts/baseline_fleet/manifests/Series60.stl` (9.3 MB) |
| Series60 manifest | `artifacts/baseline_fleet/manifests/Series60.json` |
| 完整证据包 | `runs/gap_cad_t1_baseline_fleet/evidence.json` |
| 测试套件 | `tests/test_fleet_baseline.py` (29 tests, all pass) |

---

## 8. 测试验证

```
$ python -m pytest tests/test_fleet_baseline.py -v
======================== 29 passed in 75.57s =========================
```

关键测试项：
- `test_manifest_registers_the_six_standard_hulls_with_provenance` — 6 船注册 + 来源完整
- `test_every_hull_is_watertight_manifold_and_spans_its_dimensions` — 6 船水密+流形+主尺度
- `test_production_resolution_meets_one_percent_for_all_ships` — 6 船体积误差 < 1%
- `test_production_resolution_convergence_order` — 收敛阶验证
- `test_two_independent_geometry_laws_agree_on_the_same_hull` — 双路线交叉验证
- `test_placeholder_304_triangles_is_rejected_by_name` — 占位符拒绝
- `test_byte_identical_hull_pair_is_rejected` — 重复几何拒绝
- `test_run_fleet_writes_artefacts` — 完整管线 + 工件写入
