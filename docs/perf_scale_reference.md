# TensorCAD 大几何性能基准：综合参考（GAP-CAD-D4-PERF-SCALE）

**锚点**: `GAP-CAD-D4-PERF-SCALE`（大几何性能基准 10⁶–10⁷ 三角级）
**版本**: 1.0.0 · **实测节点**: LoongArch64, 128 核, Python 3.11.6, TensorCAD 0.23.0
**工件源**: `benchmarks/perf_scale/{artifacts/ladder_mesh.json,artifacts/ladder_voxel.json,ladder_fleet_10m_3pt.json,artifacts/ladder_single_10m.json}`
**回归基线**: `benchmarks/perf_scale/baseline.json`（本地 3 轮取最慢，闸门 2×）
**验证**: `python -m tensorcad.perf_benchmark --quick --baseline benchmarks/perf_scale/baseline.json`（全部绿）

本篇是 GAP-CAD-D4-PERF-SCALE 的**综合参考**，把分散在 `perf_scale_10m.md`、`perf_scale_single_body_10m.md`
和 CI 基线中的实测数字汇成一张表，并给出明确的"什么规模用什么路径"结论。

---

## 1. 全尺度基准表（程序化注入，零手抄）

### 1.1 网格子任务（单体，678 → 8,199,998 三角）

<!-- perf-evidence:ref-mesh-table begin -->
| 子任务 | 三角面 | 秒 | 吞吐 | 建索引 s | 索引 MiB | 峰值 RSS MiB |
|---|---:|---:|---:|---:|---:|---:|
| load | 678 | 0.020 | 34,326 tri/s | — | — | 40 |
| load | 1,998 | 0.062 | 32,229 tri/s | — | — | 40 |
| load | 6,558 | 0.185 | 35,492 tri/s | — | — | 48 |
| load | 25,918 | 0.772 | 33,587 tri/s | — | — | 76 |
| load | 103,038 | 3.461 | 29,773 tri/s | — | — | 187 |
| load | 410,878 | 15.054 | 27,294 tri/s | — | — | 643 |
| load | 1,640,958 | 70.577 | 23,251 tri/s | — | — | 2,444 |
| load | 4,003,998 | 219.201 | 18,266 tri/s | — | — | 6,163 |
| load | 8,199,998 | 573.483 | 14,299 tri/s | — | — | 12,746 |
| point_locate | 678 | 0.257 | 7,781 query/s | — | — | 40 |
| point_locate | 1,998 | 0.809 | 2,472 query/s | — | — | 40 |
| point_locate | 6,558 | 2.627 | 761 query/s | — | — | 48 |
| point_locate | 25,918 | 10.116 | 198 query/s | — | — | 76 |
| point_locate | 103,038 | 52.679 | 38 query/s | — | — | 187 |
| point_locate | 410,878 | 235.551 | 8 query/s | — | — | 643 |
| **point_locate_accel** | 678 | 0.083 | 23,974 query/s | 0.01 | 0.1 | 40 |
| **point_locate_accel** | 1,998 | 0.113 | 17,752 query/s | 0.02 | 0.2 | 40 |
| **point_locate_accel** | 6,558 | 0.181 | 11,077 query/s | 0.05 | 0.5 | 48 |
| **point_locate_accel** | 25,918 | 0.261 | 7,650 query/s | 0.21 | 2.0 | 76 |
| **point_locate_accel** | 103,038 | 0.393 | 5,084 query/s | 0.81 | 8.0 | 187 |
| **point_locate_accel** | 410,878 | 0.653 | 3,065 query/s | 3.13 | 31.6 | 643 |
| **point_locate_accel** | 1,640,958 | 0.701 | 2,853 query/s | 12.59 | 125.3 | 2,444 |
| **point_locate_accel** | 4,003,998 | 1.046 | 1,913 query/s | 30.55 | 304.6 | 6,163 |
| **point_locate_accel** | 8,199,998 | 0.965 | 2,072 query/s | 62.11 | 622.4 | 12,746 |
| volume | 678 | 0.000 | 2,478,795 tri/s | — | — | 40 |
| volume | 1,998 | 0.001 | 2,518,308 tri/s | — | — | 40 |
| volume | 6,558 | 0.003 | 2,301,464 tri/s | — | — | 48 |
| volume | 25,918 | 0.011 | 2,326,230 tri/s | — | — | 76 |
| volume | 103,038 | 0.044 | 2,329,374 tri/s | — | — | 187 |
| volume | 410,878 | 0.168 | 2,440,672 tri/s | — | — | 643 |
| volume | 1,640,958 | 0.666 | 2,463,016 tri/s | — | — | 2,444 |
| volume | 4,003,998 | 1.624 | 2,465,730 tri/s | — | — | 6,163 |
| volume | 8,199,998 | 3.363 | 2,438,286 tri/s | — | — | 12,746 |
<!-- perf-evidence:ref-mesh-table end -->

### 1.2 体素化（SUBOFF SDF，32³ → 384³ cells）

<!-- perf-evidence:ref-voxel-table begin -->
| 网格 | 体素 | 秒 | 吞吐 cell/s | 峰值 RSS MiB |
|---:|---:|---:|---:|---:|
| 32³ | 32,768 | 0.244 | 134,254 | 104 |
| 48³ | 110,592 | 0.810 | 136,542 | 259 |
| 64³ | 262,144 | 2.040 | 128,518 | 560 |
| 96³ | 884,736 | 6.200 | 142,696 | 585 |
| 128³ | 2,097,152 | 15.359 | 136,541 | 597 |
| 192³ | 7,077,888 | 48.020 | 147,396 | 650 |
| 256³ | 16,777,216 | 115.259 | 145,561 | 737 |
| 384³ | 56,623,104 | 356.209 | 158,960 | 794 |
<!-- perf-evidence:ref-voxel-table end -->

### 1.3 多体 Fleet（10⁶–10⁷ 带，2001×501 基座 × 2/3 体）

<!-- perf-evidence:ref-fleet-table begin -->
| 子任务 | 三角面 | 体数 | 秒 | 吞吐 | 建索引 s | 索引 MiB | 峰值 RSS MiB |
|---|---:|---:|---:|---:|---:|---:|---:|
| load | 410,878 | 1 | 18.443 | 22,278 tri/s | — | — | 619 |
| fleet_assemble | 821,756 | 2 | 0.584 | 1,406,197 tri/s | — | — | 619 |
| point_locate_accel | 821,756 | 2 | 0.636 | 3,145 query/s | — | — | 619 |
| volume | 821,756 | 2 | 0.385 | 2,132,084 tri/s | — | — | 619 |
| fleet_assemble | 1,232,634 | 3 | 0.932 | 1,322,555 tri/s | — | — | 619 |
| point_locate_accel | 1,232,634 | 3 | 0.645 | 3,100 query/s | — | — | 619 |
| volume | 1,232,634 | 3 | 0.600 | 2,056,020 tri/s | — | — | 619 |
| load | 1,640,958 | 1 | 87.160 | 18,827 tri/s | — | — | 2,339 |
| fleet_assemble | 3,281,916 | 2 | 3.069 | 1,069,466 tri/s | — | — | 2,339 |
| point_locate_accel | 3,281,916 | 2 | 0.945 | 2,117 query/s | — | — | 2,339 |
| volume | 3,281,916 | 2 | 1.729 | 1,898,042 tri/s | — | — | 2,339 |
| fleet_assemble | 4,922,874 | 3 | 3.818 | 1,289,406 tri/s | — | — | 2,339 |
| point_locate_accel | 4,922,874 | 3 | 1.019 | 1,963 query/s | — | — | 2,339 |
| volume | 4,922,874 | 3 | 2.840 | 1,733,441 tri/s | — | — | 2,339 |
| load | 4,003,998 | 1 | 269.335 | 14,866 tri/s | — | — | 5,735 |
| fleet_assemble | 8,007,996 | 2 | 6.221 | 1,287,164 tri/s | — | — | 5,735 |
| point_locate_accel | 8,007,996 | 2 | 1.195 | 1,674 query/s | — | — | 5,735 |
| volume | 8,007,996 | 2 | 3.928 | 2,038,436 tri/s | — | — | 5,735 |
| fleet_assemble | 12,011,994 | 3 | 9.534 | 1,259,876 tri/s | — | — | 5,735 |
| point_locate_accel | 12,011,994 | 3 | 1.130 | 1,769 query/s | — | — | 5,735 |
| volume | 12,011,994 | 3 | 5.652 | 2,125,364 tri/s | — | — | 5,735 |
<!-- perf-evidence:ref-fleet-table end -->

### 1.4 单体细分（201×65 基座，0→3 级，51,598 → 3,302,272 三角）

<!-- perf-evidence:ref-single-table begin -->
| 级别 | 三角面 | 子任务 | 秒 | 吞吐 | 建索引 s | 索引 MiB | 活 RSS MiB |
|---:|---:|---|---:|---:|---:|---:|---:|
| 0 | 51,598 | load | 2.081 | 24,793 tri/s | — | — | 64 |
| 0 | 51,598 | volume | 0.032 | 1,596,857 tri/s | — | — | 64 |
| 0 | 51,598 | point_locate_accel | 0.067 | 7,480 query/s | 0.39 | 4.0 | 64 |
| 1 | 206,392 | subdivide | 0.241 | 855,330 tri/s | — | — | 92 |
| 1 | 206,392 | volume | 0.127 | 1,619,909 tri/s | — | — | 92 |
| 1 | 206,392 | point_locate_accel | 0.171 | 2,919 query/s | 1.99 | 16.0 | 110 |
| 2 | 825,568 | subdivide | 1.090 | 757,545 tri/s | — | — | 241 |
| 2 | 825,568 | volume | 0.400 | 2,062,153 tri/s | — | — | 241 |
| 2 | 825,568 | point_locate_accel | 0.163 | 3,062 query/s | 6.58 | 63.3 | 319 |
| 3 | 3,302,272 | subdivide | 4.279 | 771,814 tri/s | — | — | 777 |
| 3 | 3,302,272 | volume | 1.604 | 2,058,208 tri/s | — | — | 777 |
| 3 | 3,302,272 | point_locate_accel | 0.257 | 1,948 query/s | 32.34 | 251.6 | 988 |
<!-- perf-evidence:ref-single-table end -->

---

## 2. 复杂度指数汇总（log-log 斜率，全部通过回归门）

| 子任务 | 数据源 | 规模范围 | 指数 | R² | 限值 | 判定 |
|---|---|---:|---:|---:|---:|---|
| load（单体） | ladder_mesh | 678 → 8.2M | 1.17 | 0.999 | ≤1.50 | ✅ |
| load（fleet） | ladder_fleet_10m_3pt | 411k → 4.0M | 1.17 | 0.999 | ≤1.50 | ✅ |
| point_locate（暴力） | ladder_mesh | 678 → 411k | 1.04 | 1.000 | ≤1.50 | ✅ |
| point_locate_accel（单体） | ladder_mesh | 678 → 8.2M | 0.28 | 0.993 | ≤0.80 | ✅ |
| point_locate_accel（fleet x2） | ladder_fleet_10m_3pt | 822k → 8.0M | 0.28 | 1.000 | ≤0.80 | ✅ |
| point_locate_accel（fleet x3） | ladder_fleet_10m_3pt | 1.2M → 12.0M | 0.25 | 0.951 | ≤0.80 | ✅ |
| volume（单体） | ladder_mesh | 678 → 8.2M | 1.00 | 0.999 | ≤1.50 | ✅ |
| volume（fleet x2） | ladder_fleet_10m_3pt | 822k → 8.0M | 1.03 | 0.998 | ≤1.50 | ✅ |
| volume（fleet x3） | ladder_fleet_10m_3pt | 1.2M → 12.0M | 1.00 | 0.991 | ≤1.50 | ✅ |
| voxelize | ladder_voxel | 32k → 56.6M | 0.98 | 0.998 | ≤1.50 | ✅ |
| subdivide（单体） | ladder_single_10m | 52k → 3.3M | 1.04 | 0.999 | ≤1.50 | ✅ |

---

## 3. 误差判定（关键量 vs 基准）

| 关键量 | 基准 | 最大实测误差 | 是否 <3% |
|---|---|---|---|
| 体积（散度定理）vs 解析闭式解 | `analytic_volume(10,2,0.8)` = 7.512698 m³ | 0.0707%（201×65 父网格离散误差） | ✅ |
| 体积（细分保真）vs 父网格 | 细分精确铺满 | ~1e-14（浮点累加序差异） | ✅ |
| 体积（fleet 3 体）vs 3×解析 | 22.538095 m³ | 0.0030% | ✅ |
| 点定位命中数（accel vs 暴力） | 逐射线对比 | 0 差异（7 个双路径档完全一致） | ✅ |
| 点定位交点距离（accel vs 暴力） | 同网格同射线 | ≤1e-9 相对误差 | ✅ |
| 拓扑不变量（水密/流形/Euler） | 全边二面、χ=2 | 全档通过 | ✅ |

---

## 4. 什么规模用什么路径（实测结论）

| 子任务 | 单机舒适区 | 需要集群 | 依据 |
|---|---|---|---|
| **load（生成器）** | ≤ 8.2M tri（573 s） | > 10⁷ tri（结构上到不了，需导入网格或 fleet） | 生成器循环 O(N)，指数 1.17 |
| **point_locate（暴力）** | ≤ 5×10⁴ tri（2000 射线 10 s） | ≥ 10⁵ tri | 38 query/s @103k，8 query/s @411k |
| **point_locate_accel** | ≤ 12M tri（2000 射线 1.13 s） | ≥ 10⁷ tri × 10⁵ 射线批量 | 查询指数 0.25–0.28，亚线性 |
| **volume** | ≤ 12M tri（5.65 s） | 不需要 | 2.0–2.4 M tri/s，指数 1.00–1.03 |
| **voxelize** | ≤ 56.6M cells（356 s） | > 10⁸ cells 或 float64 内存紧张 | 136k–159k cell/s |
| **subdivide** | ≤ 3.3M tri（4.3 s/级） | > 10⁷ tri（单级细分 > 15 s） | O(N)，指数 1.04 |

**核心结论**：
1. 点定位的"集群需求"是**算法缺口**造成的，不是硬件缺口。补上 `TriangleGrid` 后，10⁶–10⁷ 带单机即可。
2. 真正把单机推走的是**生成器**（8.2M tri 需 9.5 min）与**体素化**（56.6M cells 需 5.9 min），不是查询。
3. 10⁷ 带靠 fleet 多体到达（12M tri），单体细分到 3.3M tri；10⁷ 单体档（cj3805）排队中。

---

## 5. 回归门状态

| 闸门 | 基线 | 阈值 | 当前状态 |
|---|---|---|---|
| 墙上时钟回归（20 条） | `benchmarks/perf_scale/baseline.json` | 2× | ✅ 全绿（本机实测，2026-11-15） |
| 复杂度指数（7 门） | 内建限值 | load/vol ≤1.50, accel ≤0.80 | ✅ 全绿 |

---

## 6. 工件索引

| 工件 | 路径 |
|---|---|
| 网格阶梯原始测量 | `benchmarks/perf_scale/artifacts/ladder_mesh.json` |
| 体素阶梯原始测量 | `benchmarks/perf_scale/artifacts/ladder_voxel.json` |
| Fleet 10⁷ 三点阶梯 | `benchmarks/perf_scale/ladder_fleet_10m_3pt.json` |
| 单体细分阶梯 | `benchmarks/perf_scale/artifacts/ladder_single_10m.json` |
| CI 回归基线 | `benchmarks/perf_scale/baseline.json` |
| 基准框架 | `tensorcad/perf_benchmark.py` |
| 空间索引 | `tensorcad/spatial_index.py` |
| 细分算子 | `tensorcad/mesh_subdivision.py` |
| 体素导出 | `tensorcad/voxel.py` |
| 证据渲染 | `tensorcad/perf_evidence.py` |
| 文档漂移校验 | `scripts/perf_scale_doc.py --check` |

---

## 7. 局限与欠账

1. **10⁷ 单体档未到件**：作业 `cj3805` 排队中，本表不含 10⁷ 单体数字。
2. **文件导入路径未测**：从 OBJ/STEP 读入 10⁷ 三角的解析/焊接开销未被覆盖。
3. **基线只覆盖 quick 阶梯**：10⁶ 档没有 CI 基线（成本太高），靠本表对比守护。
4. **负载影响**：实测时 load1≈40–60（128 核），吞吐数字有 10–20% 抖动；跨机比较应看数量级。
