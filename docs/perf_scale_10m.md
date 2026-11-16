# TensorCAD 大几何性能基准：10⁶–10⁷ 三角级实测（GAP-CAD-D4-PERF-SCALE）

**版本**: 2.0.0（证据包受版本控制，表格程序化注入）· **锚点**: `GAP-CAD-D4-PERF-SCALE`
**实测节点**: LoongArch64, 128 核, Python 3.11.6, TensorCAD 0.23.0
**作业**: `cj2652`（网格阶梯）、`cj2653`（体素阶梯）、`cj2669`（超限档探针，干净失败）、`cj3081`（10⁷ 多体带）、
本版重算 `cj3563`（网格阶梯）/ `cj3564`（体素阶梯）
**工件（受版本控制、可直接复核）**: `benchmarks/perf_scale/artifacts/{ladder_mesh.json,ladder_voxel.json,ladder_import.json,ladder_single_10m.json}`、
`benchmarks/perf_scale/{baseline.json,ladder_fleet_10m_3pt.json,ladder_3pt_mechanism_smoke.json}`；
图 `benchmarks/perf_scale/figures/*.png`（12 张）；渲染/校验 `tensorcad/perf_evidence.py` + `scripts/perf_scale_doc.py`

本篇取代 `docs/perf_benchmark.md` 的 0.2.0 数据（那版最高只测到 103k 三角）——**本版第一次真正测到 10⁶ 带**。

---

## 1. 案例介绍

参数化船体（L=10 m, B=2 m, T=0.8 m）网格阶梯把三角面数从 6.4×10² 推到 4.0×10⁶（4 个数量级），四个子任务在同一阶梯上测：

| 子任务 | 实现路径 | 复杂度 | 备注 |
|---|---|---|---|
| `load` | `hull.generate_hull` | O(N) Python 循环 | 每次真生成（不复用缓存） |
| `point_locate` | `benchmark._RayCaster` 逐三角遍历 | O(N·Q) | 无加速结构，参考实现 |
| `point_locate_accel` | **新增** `spatial_index.TriangleGrid`（均匀网格 + DDA） | O(cells·Q) | 本单新增路径 |
| `volume` | `mesh.total_signed_volume` 散度定理 | O(N) | 解析值 7.512698 m³ 为基准 |
| `voxelize` | `voxel.to_voxel_grid`（SUBOFF SDF, NumPy 向量化） | O(G³) | 32³–384³ |

**为什么先补加速结构**：旧文档记录点定位在 10⁵ 三角面只有 13 条射线/秒，结论是"超过 5×10⁴ 三角面建议用集群"。10⁶ 带上按此推算单条射线就要 ~1 s，2000 条射线 ≈ 半小时——这说明瓶颈是**算法**而不是机器。因此本单先实现空间索引，再在同一阶梯上对比两条路径，用实测判定"什么规模真的需要集群"。

### 验收标准（feature_ledger）
1. 加载/点定位/体素化/体积积分子任务基准表 → 44 条测量记录（0.678k → 8,199,998 三角，体素 3.3e4 → 5.66e7 cells），另有 10⁶–10⁷ 带多体三点阶梯 21 条记录（→ 12,011,994 三角）
2. 性能回归门（超基线 2×=RED）→ 已落地 CI 基线 `benchmarks/perf_scale/baseline.json`（本地 3 轮取最慢，闸门 20 条全绿）
3. 记录什么规模需要集群 → 见 §5

---

## 2. 怎么做的

**新增 `tensorcad/spatial_index.py`（stdlib-only）**
- **CSR 桶布局**：三角面引用存在扁平 `array('i')`（32 位）里，按 cell 前缀和寻址（两趟计数排序）。10⁷ 三角面量级下是几百 MB，而不是 `dict[cell]→list` 的多 GB。
- **近立方单元**：各轴 cell 数与包围盒尺度成正比，细长船体盒不会退化成薄片。
- **超大三角面回退表**：包围盒覆盖 cell 数超阈值的三角面不进网格（否则一个大面能把内存打爆），而是每条射线都测——答案仍然精确。
- **双精度 Möller–Trumbore**，方向不必归一化；`cast_ray` 取最近正交点，`cast_all` 取全部交点（共享边重复命中按容差去重），`contains` 用奇偶性判定内外。

**基准框架改造（`tensorcad/perf_benchmark.py`）**
- 新子任务 `point_locate_accel`：分开报告 `build_seconds`（一次性建索引）与 `seconds`（2000 条射线查询）。
- **网格缓存（上限 2 档）**：旧实现每个子任务都重新生成网格，10⁶ 档会把全部墙钟花在生成器上；现在每档只生成一次。
- `max_naive_triangles`（默认 3×10⁵）：超过就跳过 O(N·Q) 参考实现并标 `not_run`（注：本次作业用 5×10⁵ 阈值，故 410,878 档仍有暴力参考数据）。
- 每条记录新增 `peak_rss_mb` / `index_cells` / `index_entries` / `index_bytes`；schema 升到 `tensorcad.perf_benchmark/2`。
- **回归基线迁移**：旧基线写在 `runs/perf_baseline.json`，而合并门的路径过滤会把 `runs/` 从分支里剔掉——也就是说回归闸门在 main 上一直是**空转**（拿不到基线），`tests/test_perf_benchmark.py::test_baseline_file_exists_and_valid` 在干净检出上必然失败。现迁移到受版本控制的 `benchmarks/perf_scale/baseline.json`，由 `scripts/build_perf_baseline.py`（3 轮取每键最慢，抑制共享机负载噪声）生成。

**执行方式**：全部重计算经 `run_compute` 递交（本机只跑 ≤2 min 的测试），节点由调度器指派。

---

## 3. 计算结果（全部来自受版本控制的 `benchmarks/perf_scale/artifacts/*.json`，程序化渲染）

### 3.1 网格子任务

下表由 `benchmarks/perf_scale/artifacts/ladder_mesh.json` **渲染**（`python scripts/perf_scale_doc.py --write`），
校验 `python scripts/perf_scale_doc.py --check` 在表与工件不一致时失败——不再手抄。

<!-- perf-evidence:mesh-table begin -->
| 子任务 | 三角面 | 秒 | 吞吐 | 建索引 s | 索引 MB | 峰值 RSS MiB |
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
| **point_locate** | 1,638,400 | not_run | brute-force reference skipped above 500,000 triangles (O(N) per ray); use point_locate_accel | — | — | — |
| **point_locate** | 4,000,000 | not_run | brute-force reference skipped above 500,000 triangles (O(N) per ray); use point_locate_accel | — | — | — |
| **point_locate** | 8,192,000 | not_run | brute-force reference skipped above 500,000 triangles (O(N) per ray); use point_locate_accel | — | — | — |
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
<!-- perf-evidence:mesh-table end -->

**同档加速比**（2000 条射线，同一起点/方向集）：25.9k 档 38.7×；103k 档 133.9×；411k 档 361.0×。10⁶ 档暴力路径已不可测（411k 档实测 235.6 s，吞吐 8 query/s；按 O(N) 线性外推 2000 条射线：1.64M 档 ≈ 16 min，4.0M 档 ≈ 38 min，8.2M 档 ≈ 78 min——这也解释了 `max_naive_triangles` 为何必须跳过这些档）。

### 3.2 体素化

同样由 `benchmarks/perf_scale/artifacts/ladder_voxel.json` 渲染。

<!-- perf-evidence:voxel-table begin -->
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
<!-- perf-evidence:voxel-table end -->

### 3.3 图

**基础阶梯图**（`ladder_mesh.json` + `ladder_voxel.json`）：
- `scaling.png`：墙钟 vs 规模 log-log（4 子任务，678 → 8.2M 三角）
- `point_locate.png`：暴力 vs 索引点定位吞吐对比
- `index_cost.png`：建索引时间/索引体积 vs 三角面
- `voxelize.png`：体素化缩放（32³ → 384³）

**10⁶–10⁷ 三点阶梯图**（`ladder_fleet_10m_3pt.json`）：
- `fig3pt_scaling.png`：每个子任务的拟合指数与门限（§8.2）
- `fig3pt_memory.png`：网格占用/索引/活内存 vs 规模

**OBJ 导入/导出阶梯图**（`ladder_import.json`）：
- `import_scaling.png`：load/obj_write/obj_read/volume 墙钟 vs 规模
- `import_throughput.png`：各子任务吞吐 vs 规模
- `import_volume_accuracy.png`：OBJ 往返体积精度（as-imported vs 解析）

**单体 10⁷ 细分阶梯图**（`ladder_single_10m.json`）：
- `single10m_scaling.png`：load/subdivide/volume/point_locate_accel 墙钟 vs 规模
- `single10m_memory.png`：网格/索引/RSS 内存 vs 规模
- `single10m_volume_accuracy.png`：细分体积精度（机器精度级守恒）

### 3.4 OBJ 导入/导出阶梯

由 `benchmarks/perf_scale/artifacts/ladder_import.json` 渲染。4 档（25,918 → 4,003,998 三角），
4 个子任务：load → obj_write → obj_read → volume 完整往返。

| 子任务 | 三角面 | 秒 | 吞吐 | 文件大小 | 文件吞吐 | 体积误差 |
|---|---:|---:|---:|---:|---:|---:|
| load | 25,918 | 0.751 | 34,492 tri/s | — | — | — |
| obj_write | 25,918 | 0.057 | 454,683 tri/s | 1.0 MiB | 18.0 MiB/s | — |
| obj_read | 25,918 | 0.093 | 278,628 tri/s | 1.0 MiB | 11.0 MiB/s | −0.136% |
| volume | 25,918 | 0.010 | 2,520,007 tri/s | — | — | −0.136% |
| load | 410,878 | 15.406 | 26,671 tri/s | — | — | — |
| obj_write | 410,878 | 1.187 | 346,078 tri/s | 18.2 MiB | 15.3 MiB/s | — |
| obj_read | 410,878 | 1.676 | 245,172 tri/s | 18.2 MiB | 10.9 MiB/s | −0.017% |
| volume | 410,878 | 0.176 | 2,340,416 tri/s | — | — | −0.017% |
| load | 1,640,958 | 110.528 | 14,847 tri/s | — | — | — |
| obj_write | 1,640,958 | 5.372 | 305,477 tri/s | 75.4 MiB | 14.0 MiB/s | — |
| obj_read | 1,640,958 | 7.586 | 216,328 tri/s | 75.4 MiB | 9.9 MiB/s | −0.006% |
| volume | 1,640,958 | 0.751 | 2,185,953 tri/s | — | — | −0.006% |
| load | 4,003,998 | 283.716 | 14,113 tri/s | — | — | — |
| obj_write | 4,003,998 | 11.013 | 363,585 tri/s | 206.3 MiB | 18.7 MiB/s | — |
| obj_read | 4,003,998 | 17.095 | 234,217 tri/s | 206.3 MiB | 12.1 MiB/s | −0.003% |
| volume | 4,003,998 | 1.877 | 2,133,446 tri/s | — | — | −0.003% |

**复杂度指数**（4 点拟合）：load 1.19、obj_write 1.06、obj_read 1.04、volume 1.03 —— 全部 ≤1.50 限值，无退化。

**OBJ 往返精度**：as-imported 体积 vs 解析基准，最大误差 −0.136%（25,918 档），4.0M 档仅 −0.003%。
OBJ 文本格式在 4.0M 三角面时文件 206 MiB，写入 18.7 MiB/s、读取 12.1 MiB/s。

### 3.5 单体 10⁷ 细分阶梯

由 `benchmarks/perf_scale/artifacts/ladder_single_10m.json` 渲染。
以 1→4 中点细分为几何源，单个**连通**闭合体直接做到 13.1M 三角面（非 fleet 多体装配）。

| 子任务 | 三角面 | 秒 | 吞吐 | 建索引 s | 索引 MiB | 活 RSS MiB |
|---|---:|---:|---:|---:|---:|---:|
| load (401×129) | 205,598 | 7.268 | 28,287 tri/s | — | — | 106 |
| subdivide → 822k | 822,392 | 0.874 | 941,178 tri/s | — | — | 242 |
| subdivide → 3.3M | 3,289,568 | 3.852 | 854,025 tri/s | — | — | 754 |
| subdivide → 13.2M | 13,158,272 | 17.531 | 750,587 tri/s | — | — | 2,530 |
| load (1601×513) | 3,279,998 | 186.747 | 17,564 tri/s | — | — | 2,953 |
| subdivide → 13.1M | 13,119,992 | 16.890 | 776,773 tri/s | — | — | 2,530 |
| point_locate_accel @13.2M | 13,158,272 | 1.549 | 1,291 query/s | 98.0 | 998 | 3,568 |
| volume @13.2M | 13,158,272 | 5.861 | 2,245,011 tri/s | — | — | 2,530 |

**体积精度**：细分不改变几何，体积误差继承自父网格（401×129 档 −0.0246%，1601×513 档 −0.0030%），
细分后各档体积漂移 ≤4.6×10⁻¹³（机器精度级），水密/流形/Euler 特征全部保持。

**内存**：13.2M 三角面网格对象 2.86 GB，空间索引 998 MiB，活 RSS 3.57 GB。

---

## 4. 误差判定（关键量 vs 基准）（关键量 vs 基准）

| 量 | 基准 | 实测 | 相对误差 | <3%? |
|---|---|---|---|---|
| 体积 @103,038 tri | 7.512698 m³（解析） | 7.509115 | −0.048% | 是 |
| 体积 @1,640,958 tri | 7.512698 m³ | 7.512257 | −0.0059% | 是 |
| 体积 @4,003,998 tri | 7.512698 m³ | 7.512474 | −0.0030% | 是 |
| 体积 @8,199,998 tri（单体生成上限） | 7.512698 m³ | 7.512493 | −0.0027% | 是 |
| 体积 @12,011,994 tri（3 体 fleet） | 22.538095 m³ | 22.537421 | −0.0030% | 是 |
| OBJ as-imported @25,918 tri | 7.512698 m³ | 7.502488 | −0.136% | 是 |
| OBJ as-imported @4,003,998 tri | 7.512698 m³ | 7.512474 | −0.003% | 是 |
| 细分体积 @13,158,272 tri（401×129 父） | 7.512698 m³ | 7.510847 | −0.0246% | 是 |
| 细分体积 @13,119,992 tri（1601×513 父） | 7.512698 m³ | 7.512473 | −0.0030% | 是 |
| 细分体积漂移（机器精度） | 0 | ≤4.6×10⁻¹³ | — | 是 |
| `cast_ray` 交点距离 | 暴力 Möller–Trumbore（同网格） | 一致 | ≤1e-9 相对 | 是 |
| `contains` 内外判定 | 闭网格几何（Cb/体积符号） | 全对 | 0 错判 | 是 |
| 命中数（accel vs 暴力，2000 射线） | 154/185/218/253/286/317（7 个双路径档） | 逐档完全一致 | 0 | 是 |
| 复杂度指数（10⁶–10⁷ 带，7 门） | 限值 load/assemble/volume ≤1.50、accel ≤0.80 | 1.17 / 1.05 / 1.02 / 1.03 / 1.00 / 0.28 / 0.25 | 全部低于限值 | 是（见 §9.4） |

**结论：关键几何量全部 <3%，点定位与参考实现逐位一致；加速结构不改变答案。**

---

## 5. 什么规模需要集群（实测判定，取代旧推测）

| 子任务 | 单机舒适区（实测） | 需要在集群上做 | 依据（全部为 `cj3563`/`cj3564` 实测） |
|---|---|---|---|
| `load` | 8,199,998 tri（573.5 s）——生成器上限 | > 10⁷ tri 单网格结构上做不到（§6 第 1 条），必须改用导入网格或 fleet 多体 | 22.3k→14.3k tri/s（生成器循环，指数 1.17） |
| `point_locate`（暴力） | ≤ 5×10⁴ tri（2000 射线 10 s） | ≥ 10⁵ tri | 200 query/s @25.9k，38 query/s @103k，1.2 query/s @411k |
| `point_locate_accel` | **8,199,998 tri（2000 射线 0.965 s；建索引 62.1 s / 622 MiB）** | 单机上 ≥10⁷ tri × 10⁵ 射线级批量查询 | 查询指数 0.28（x2）/0.25（x3），亚线性；索引 954 MiB @12.0e6 tri |
| `volume` | 12,011,994 tri（5.65 s） | 不需要 | 2.0–2.4 M tri/s，指数 1.03/1.00 |
| `voxelize` | 56,623,104 cells（356.2 s, 794 MiB 峰值） | >10⁸ cells 或 dtype=float64 且内存紧张时 | 136k–159k cell/s |

**核心结论**：点定位的"集群需求"是**算法缺口**造成的，不是硬件缺口。同一个 410,878 三角面的网格，暴力路径 235.6 s，加索引后 0.653 s（361×）——补上加速结构后，"10⁵ 三角面就该上集群"的判断被推翻，10⁶–10⁷ 带单机即可（12.0e6 tri 的 2000 条射线查询 1.13 s，索引 910 MiB）。真正把单机推走的是**生成器**（8.2e6 tri 需 9.5 min、网格对象 1.74 GB）与**体素化**（56.6e6 cells 需 5.9 min），不是查询。

---

## 6. 注意点 / 局限（诚实记录）

1. **10⁷ 带靠 fleet 多体到达，单网格到不了**：`HullParameters` 上限（`nz ≤ 513 levels`、`nx ≤ 4001 stations`）使**单体**网格最大 **8,199,998 三角面**（4001×513，`cj3563` 已实测，load 573.5 s）；10⁷ 带由 `build_fleet_mesh` 的 3 个互不相交副本达到 **12,011,994 三角面**（`cj3514` 实测）。因此"10⁷ 三角面"是**多体装配**规模，不是"一个 10⁷ 三角面网格"——导入式大网格仍未被实测，这条限制是结构性的，不是本单没跑。**部分已补（见 `docs/perf_scale_single_body_10m.md`）**：以 1→4 中点细分为几何源，单个**连通**闭合体直接做到 10⁷ 三角面（体积/面积按构造守恒，作为校验而非估计）；"从文件导入 10⁷ 三角"这段欠账由 **§11** 交付（`obj_write`/`obj_read` 一等子任务 + 指数门 + 2× 回归基线键）。
2. **越限档是干净失败**：`cj2669` 的两档返回 `status=failed` 并附生成器报错原话，没有静默跳过、也没有伪造数据；框架把异常记进结果表而非吞掉。
3. **暴力参考在 4.1×10⁵ 以上被主动跳过**（`max_naive_triangles`），标 `not_run` 而不是给一个不存在的秒数。
4. **回归基线的口径**：`benchmarks/perf_scale/baseline.json` 是共享机 3 轮取最慢，闸门阈值 2×；它只覆盖 quick 阶梯（≤2.6×10⁴ tri），10⁶ 档没有 CI 基线（成本太高），10⁶ 档的守护靠本单表格对比。
5. **负载影响**：n315 实测时 load1≈40（128 核），吞吐数字有 10–20% 抖动；同档之间比较应看数量级。
6. **`_RayCaster` 的角度前缀过滤**只对星形面（船体）成立，通用网格请用 `TriangleGrid` 或 `brute_force_nearest`。
7. **内存（实测口径已分开）**：4.0×10⁶ 档网格对象 871 MiB，索引 305 MiB，活 RSS 1.39 GB；8.2×10⁶ 档网格 1.74 GB、索引 622 MiB、活 RSS 2.95 GB。`peak_rss_mb` 是**进程高水位**（`ru_maxrss`，前面档位撑起来的峰值会留在后面），8.2e6 档显示 12.7 GB 是因为同一进程里先缓存过 4.0e6 网格；"这一档要多少 RAM"应读 `mesh_mb`+`index_bytes`+`rss_mb`。12.0e6 tri 的 3 体 fleet：网格 2.61 GB、索引 910 MiB、活 RSS 3.50 GB（`cj3514`）。

---

## 7. 工件索引

| 工件 | 路径 |
|---|---|
| 空间索引 | `tensorcad/spatial_index.py` |
| 索引测试（17 项，含与暴力实现逐位对比） | `tests/test_spatial_index.py` |
| 基准框架（含 `point_locate_accel`、网格缓存、峰值 RSS） | `tensorcad/perf_benchmark.py` |
| CI 回归基线（受版本控制） | `benchmarks/perf_scale/baseline.json` |
| 基线生成器 | `scripts/build_perf_baseline.py` |
| 原始测量（受版本控制） | `benchmarks/perf_scale/artifacts/{ladder_mesh,ladder_voxel}.json`、`benchmarks/perf_scale/ladder_fleet_10m_3pt.json` |
| 表格渲染/漂移校验 | `tensorcad/perf_evidence.py`、`scripts/perf_scale_doc.py` |
| 证据包测试 | `tests/test_perf_scale_evidence.py` |
| 10⁶–10⁷ 三点阶梯驱动 / 表格 / 作图 | `scripts/perf_scale_10m.py`、`scripts/perf_3pt_report.py`、`scripts/gen_perf_3pt_figures.py` |
| 图（基础阶梯） | `benchmarks/perf_scale/figures/{scaling,point_locate,index_cost,voxelize}.png` |
| 图（10⁶–10⁷ 三点） | `benchmarks/perf_scale/figures/{fig3pt_scaling,fig3pt_memory}.png` |
| 图（OBJ 导入/导出） | `benchmarks/perf_scale/figures/{import_scaling,import_throughput,import_volume_accuracy}.png` |
| 图（单体 10⁷ 细分） | `benchmarks/perf_scale/figures/{single10m_scaling,single10m_memory,single10m_volume_accuracy}.png` |
| 作图脚本 | `scripts/gen_perf_scale_figures.py`、`scripts/gen_perf_3pt_figures.py`、`scripts/gen_perf_import_figures.py`、`scripts/gen_perf_single10m_figures.py` |

---

## 8. 10⁷ 带：多体 (fleet) 几何源 + 复杂度回归门（新增，`tensorcad.perf_benchmark`）

### 8.1 为什么需要 fleet 源

单船参数化船体受 `HullParameters` 上限约束（`nx ≤ 4001` stations、`nz ≤ 513` levels），
最大只到 **8,192,000 三角面**，10⁷ 带用 `generate_hull` 结构上到不了（§6 第 1 条的欠账）。
新增 `build_fleet_mesh(model, copies)`：把 k 个互不相交的平移副本（默认间距 2 个型宽，留 1 倍型宽水隙）
组装成一个网格，**k 倍三角面以 O(N) 拷贝得到，不起生成器循环**；因为刚体彼此不接触，
体积积分保留闭式解 —— 总排水 = `copies × hull.analytic_volume(10,2,0.8)`，每个 rung 的
`note` 里直接写出 `analytic=` 与 `rel_err=`，快而不准的 rung 无法蒙混。

诚实边界（写进代码 docstring）：fleet 不是"单个大网格导入"的替代品；装配开销单列为
`fleet_assemble` 子任务与 `load`（生成）分开报，二者的工程问题不同（生成器规模 vs 纯拷贝规模）。
`_RayCaster` 是星形单体参考实现，`copies > 1` 时直接 `ValueError`，不静默给错答案。

CLI：`python -m tensorcad.perf_benchmark --fleet-tier --ladder 2001,501 --fleet-copies 2,3 ...`
（2001×501 单体 = 4.0e6 → 2 体 8.0e6 → 3 体 1.2e7，一次生成两档装配）。

### 8.2 复杂度（指数）回归门 —— 补上 10⁶ 带没有 CI 闸门的洞

§6 第 4 条承认：墙上时钟基线只覆盖 quick 阶梯（≤2.6×10⁴ tri），10⁶ 带没有 CI 基线。
新增 `scaling_gates()`：对每个子任务拟合 `log(秒) vs log(规模)` 的最小二乘指数，超限即 RED。

- 这个门是**尺度不变量**：把某个 O(N) 通道改成 O(N²)，即便是小阶梯也会给出指数 ~2 → RED；
  而 2× 墙上时钟门在小网格上完全看不见。
- 数据不足（不同规模 < `MIN_SCALING_POINTS=3`）时报 `status=not_run`，**不是**静默通过
  （`GateVerdict.status` 新字段，已进 JSON）。
- 上限：`load/fleet_assemble/point_locate/volume/voxelize ≤ 1.5`、`point_locate_accel ≤ 0.8`。

实测（本机 quick 阶梯，4 档，`python -m tensorcad.perf_benchmark --quick --n-queries 200 --baseline benchmarks/perf_scale/baseline.json`，26 个闸门全绿）：

| 门 | 拟合指数 | 上限 | 判定 |
|---|---:|---:|---|
| `load` | 1.08 | 1.50 | ✅ pass |
| `point_locate` | 1.08 | 1.50 | ✅ pass |
| `point_locate_accel` | 0.39 | 0.80 | ✅ pass |
| `volume` | 1.11 | 1.50 | ✅ pass |
| `voxelize` | 1.07 | 1.50 | ✅ pass |

结论：加速点定位路径指数 0.39（远低于 O(N)），体积/生成/暴力路径 ~1.08–1.11（线性，留 0.4 余量）。
`tests/test_perf_benchmark.py::test_scaling_gate_flags_quadratic_blowup` 是**反证测试**：
注入指数 ~2 的合成数据必须 RED；`test_scaling_gate_is_not_run_below_min_points` 守住 `not_run` 语义。

### 8.3 10⁷ 带实测：已拿到（`cj3081`，n330, rc=0）

上一版本的这一节是一份 honest FAIL，并把根因写成「取件口径：作业写在 worktree 里的新文件不会被 fetch 打包」。
**这个诊断是错的**，本单已把它推翻（见 §8.4 第 1 条）：两次失败作业的日志原文是
`bash: line 1: cd: /root/tp-agent/work/m20813: Permission denied`（rc=1）——作业根本没启动，
被测工作树属于另一张单（m20813）且节点侧无读权限。取件侧一切正常：本轮 `cj3081` 的
`artifacts/perf_scale/ladder_fleet_10m.json` 原样出现在 fetch 包里。

命令（经 `run_compute` 递交，cores=4，作业墙钟 523 s）：

```
python scripts/perf_scale_10m.py --base 2001:501 --copies 2,3 --n-queries 2000 \
    --output artifacts/perf_scale/ladder_fleet_10m.json
```

节点 n330（LoongArch64, 128 核），实测时 load1≈49/128。

| 子任务 | 体数 | 三角面 | 秒 | 吞吐 | 建索引 s | 索引 MB | 峰值 RSS MB |
|---|---:|---:|---:|---:|---:|---:|---:|
| `load`（生成单体） | 1 | 4,003,998 | 321.14 | 12,468 tri/s | – | – | 5,770 |
| `fleet_assemble` | 2 | 8,007,996 | 6.52 | 1.23 M tri/s | – | – | 5,770 |
| `point_locate_accel`（2000 射线） | 2 | 8,007,996 | 1.285 | 1,557 query/s | 68.31 | 637 | 5,770 |
| `volume` | 2 | 8,007,996 | 4.01 | 2.00 M tri/s | – | – | 5,770 |
| `fleet_assemble` | 3 | 12,011,994 | 9.50 | 1.26 M tri/s | – | – | 5,770 |
| `point_locate_accel`（2000 射线） | 3 | **12,011,994** | **1.247** | **1,604 query/s** | 100.26 | 954 | 5,770 |
| `volume` | 3 | **12,011,994** | **6.05** | 1.99 M tri/s | – | – | 5,770 |

索引规模：x2 档 2,030,644 cells、entries/tri=1.37；x3 档 3,030,075 cells、entries/tri=1.36
（近立方单元选址在细长多体包围盒上未爆桶）。原始数据：`benchmarks/perf_scale/ladder_fleet_10m.json`
（同一份也留在 `artifacts/perf_scale/`），作业 stdout 另有 gzip+base64 payload（§8.5）。

**精度（关键量 vs 解析基准，逐项）**

| 量 | 基准（`copies × analytic_volume(10,2,0.8)`） | 实测 | 相对误差 | <3%? |
|---|---|---:|---:|---|
| 体积 @8,007,996 tri（2 体） | 15.025396 m³ | 15.024947 m³ | −0.0030% | 是 |
| 体积 @12,011,994 tri（3 体） | 22.538095 m³ | 22.537421 m³ | −0.0030% | 是 |

**10⁷ 带的判定（取代 §5 的外推）**：12.0e6 三角面上，2000 条射线经网格加速只需 **1.25 s**。
按 §3.1 实测的暴力吞吐（7.9 query/s @411k，近似 ∝1/N）外推，同一查询集需
`2000 / (7.9 × 411k/12.0M) ≈ 7.4e3 s ≈ 2.1 h`，即加速结构在 10⁷ 带把点定位成本压低 **≈5.9e3 倍**。
结论更强：**10⁷ 三角面的点定位不需要集群**，单节点 1.25 s。真正的单机瓶颈是
(a) 参数化生成器（4.0e6 用 321 s，10⁷ 约 800 s）与 (b) 一次性建索引（12.0e6 用 100 s / 954 MB），
而不是查询本身。

**回归门（诚实口径）**：fleet 阶梯只有一个基础档 × 两个体数，每个 `scaling_*` 门都只有 1 个规模点，
因此全部 `not_run`，本次运行 `verdict=not_run`——**不是 pass**。§8.2 的 26 门全绿来自 quick 阶梯
（4 档、指数可拟合），两套数据不可混用。

### 8.4 本单修正的两个结论

1. **「10⁷ 取不回来」的根因在提交侧，不在取件侧**：`cj2967`/`cj2976` 的作业日志均为
   `cd: /root/tp-agent/work/m20813: Permission denied`、rc=1，作业一行都没执行。教训：递交前核对
   `stage_dir` 是否本工单工作树；失败要读 `jobs/cjN.log` 的 `rc` 与 `output_tail`，
   不能用未经证据支持的归因（原「取件口径」说法已被本轮 fetch 包直接证伪）。
2. **`not_run` 门曾被计入「全绿」**：`GateVerdict(passed=True, status=not_run)` 让
   `PerfBenchmarkResult.passed` 为 True、CLI 打印 "All gates passed."，而实际一个门都没算。
   现新增 `gates_checked` / `gates_not_run` / `verdict`（`pass`/`fail`/`not_run` 三态）并进 JSON；
   CLI 在无门可算时打印 "This is NOT a pass."。反证测试：
   `tests/test_perf_benchmark.py::test_verdict_distinguishes_not_checked_from_passed`、
   `::test_cli_does_not_claim_pass_when_no_gate_was_checked`（后者跑单档阶梯，断言输出含 "NOT a pass"）。

### 8.5 结果传输：三条独立通道（`tensorcad/perf_scale_artifact.py`）

一次 5–10 分钟的重计算不该死于打包约定，结果走三条互不依赖的通道：

1. `--output` JSON 落在**受版本控制**的 `benchmarks/perf_scale/ladder_fleet_10m.json`
   （`artifacts/` 部分子路径会被发布门禁滤掉，故基准产物不放那里）；
2. `PERF10M_JSON_BEGIN/END` 之间的 gzip+base64 payload（按 900 字符分块），可用
   `python -m tensorcad.perf_scale_artifact <job.log> -o recovered.json` 从日志还原；
3. 人读表格（子任务/体数/三角面/秒/建索引/索引 MB/RSS）。

`parse_emitted` 对**截断**报错（缺块=错误，绝不返回部分结果）；`tests/test_perf_scale_artifact.py`
（8 项）覆盖往返无损、分块、丢块、缺标记、坏 base64、日志恢复与失败退出码。

**局限（诚实记录）**
- `peak_rss_mb` 在 x2/x3 档与 `load` 档完全相同（5,770 MB），而 12.0e6 三角面的 Python 元组
  按体积推断不该停在这个量级；该字段取进程 `ru_maxrss` 高水位，在 loongarch64/n330 上疑似
  报告粒度或释放时机问题。**本单不基于 RSS 给出 10⁷ 内存结论**，需单独标定。
- 体数只有 2 和 3（无第三个规模点），所以复杂度门 `not_run`。要拿到 10⁷ 带的指数门，
  需再补一档，例如 `--base 1281:321 --copies 2,4,8`（1.6e6/3.3e6/6.6e6/1.3e7）。

---

## 9. 三点阶梯：把 10^6–10^7 带的指数门从 not_run 变成可判定（机制完成，10^7 重算排队中）

### 9.1 为什么 §8.3 的指数门必然全部 not_run

`scaling_gates()` 在**一个（子任务, 体数）分组内**拟合 log(秒) vs log(规模)。§8.3 的扶手只有
**一个基础档**（2001×501）× 两个体数（2、3），于是每个分组只有 **1 个规模点**：四个子任务的
not_run 不是数据坏了，而是阶梯形状本身不可能拟合指数。§8.5 末尾已把这条写成欠账。

本单落地的修法两件：

1. `scripts/perf_scale_10m.py --bases "641:161,1281:321,2001:501"`：**多基础档**，同一体数下
   得到 3 个规模点 → 分组内可拟合指数（`--base` 单档用法保持不变）。
2. `tensorcad.perf_benchmark.FLEET_SUBTASKS`：多体扶手只声明**适用**子任务
   （`point_locate_accel`、`volume`）。`point_locate`（暴力参考假定星形单体）与 `voxelize`
   （按网格分辨率而非网格规模成阶梯）改为每档一条 `status=not_run` + 原因，不再跑出 `failed`
   行——`failed` 读起来像测量事故，实际是"结构上不适用"。

### 9.2 机制验证（本机小阶梯：门真的会变 measured）

命令：`python -m tensorcad.perf_benchmark --fleet-tier --ladder "21,9;41,13;61,19"
--fleet-copies 2 --n-queries 200 --output benchmarks/perf_scale/ladder_3pt_mechanism_smoke.json`

| 门 | 判定 | 细节 |
|---|---|---|
| `scaling_fleet_assemble[x2]` | pass | exponent=0.94 (limit 1.50), R^2=1.000, 3 sizes 1,356..8,876 |
| `scaling_load` | pass | exponent=0.88 (limit 1.50), R^2=0.977, 3 sizes 678..4,438 |
| `scaling_point_locate_accel[x2]` | pass | exponent=0.22 (limit 0.80), R^2=0.512, 3 sizes 1,356..8,876 |
| `scaling_volume[x2]` | pass | exponent=0.91 (limit 1.50), R^2=0.954, 3 sizes 1,356..8,876 |

gates_checked=4、gates_not_run=0、verdict=**pass**。
这些指数来自 10^3–10^4 三角面的合成小阶梯，**只证明"三点阶梯 → 门可判定"这条机制，
不构成任何性能结论**（小网格计时被解释器开销主导）。

### 9.3 内存量纲：三个字段各自是什么（修正 §8.5 第一条局限）

§8.5 记录"x2/x3 档 `peak_rss_mb` 与 load 档完全相同（5,770 MB），疑似报告粒度问题"。本单把语义
讲清，并补两个新字段：

| 字段 | 来源 | 语义 |
|---|---|---|
| `peak_rss_mb` | `resource.getrusage(RUSAGE_SELF).ru_maxrss` | **进程级高水位，不能逐档复位**；前面档位达到过的峰值会一直显示在后面档上。本机实测该字段工作正常（分配 300 MB 后 20.6 → 315.4 MB，与 `VmHWM` 一致），**不是** loongarch 粒度缺陷。 |
| `rss_mb`（新增） | `/proc/self/status` `VmRSS` | 该档结束时的**活内存**，能区分档位。 |
| `mesh_mb`（新增） | 采样 `sys.getsizeof` 外推（`mesh_footprint_mb`） | 网格容器自身占多少内存，**与机器负载无关**；这才是"10^7 三角面要多少 RAM"的可复现答案。 |

于是 §8.3 表里 x2/x3 峰值相同，最简解释是"后续档位没有超过 load 档（4.0e6 生成）已达的进程高水位"；
给出逐档真实内存须用 `rss_mb`/`mesh_mb` 重测。**本单不下 10^7 内存结论**：`ladder_fleet_10m.json`
是旧版 JSON，没有这两个字段。

### 9.4 10^6–10^7 带实测：三点阶梯已到件，指数门欠账清零（本版）

§9.4 上一版是欠账（`cj3514` 递交 28 min 仍 `queued`，`10^7` 指数门全 `not_run`）。本版把三个作业
全部取回到件，**七个复杂度门从 `not_run` 变成 `measured`，verdict=pass**：

| 作业 | 内容 | 节点 / 返回码 | 环境（load1） | 落盘工件 |
|---|---|---|---|---|
| `cj3514` | 三点 fleet 阶梯 `--bases 641:161,1281:321,2001:501 --copies 2,3`（0.82e6→12.01e6 tri） | n330 / 0 | loongarch64, 128 核, load1 90.76 | `benchmarks/perf_scale/ladder_fleet_10m_3pt.json` |
| `cj3563` | 单体全阶梯到生成器上限 4001×513（8,199,998 tri） | n330 / 0 | load1 46.76 | `benchmarks/perf_scale/artifacts/ladder_mesh.json` |
| `cj3564` | 体素阶梯到 384³（56,623,104 cells） | n330 / 0 | load1 51.11 | `benchmarks/perf_scale/artifacts/ladder_voxel.json` |

两张表都由 `tensorcad.perf_evidence` 从提交的 JSON **渲染**（`python scripts/perf_scale_doc.py --write`），
`--check` 与 `tests/test_perf_scale_evidence.py` 在表与工件不一致时失败——与 §3 同一机制，没有手抄。

<!-- perf-evidence:fleet3pt-gates begin -->
| 门 | 判定 | 拟合指数 | 限值 | R² | 规模点数 | 细节 |
|---|---|---:|---:|---:|---:|---|
| `scaling_fleet_assemble[x2]` | pass | 1.05 | 1.50 | 0.989 | 3 | exponent=1.05 (limit 1.50), R^2=0.989, 3 sizes 821,756..8,007,996 |
| `scaling_fleet_assemble[x3]` | pass | 1.02 | 1.50 | 1.000 | 3 | exponent=1.02 (limit 1.50), R^2=1.000, 3 sizes 1,232,634..12,011,994 |
| `scaling_load` | pass | 1.17 | 1.50 | 0.999 | 3 | exponent=1.17 (limit 1.50), R^2=0.999, 3 sizes 410,878..4,003,998 |
| `scaling_point_locate_accel[x2]` | pass | 0.28 | 0.80 | 1.000 | 3 | exponent=0.28 (limit 0.80), R^2=1.000, 3 sizes 821,756..8,007,996 |
| `scaling_point_locate_accel[x3]` | pass | 0.25 | 0.80 | 0.951 | 3 | exponent=0.25 (limit 0.80), R^2=0.951, 3 sizes 1,232,634..12,011,994 |
| `scaling_volume[x2]` | pass | 1.03 | 1.50 | 0.998 | 3 | exponent=1.03 (limit 1.50), R^2=0.998, 3 sizes 821,756..8,007,996 |
| `scaling_volume[x3]` | pass | 1.00 | 1.50 | 0.991 | 3 | exponent=1.00 (limit 1.50), R^2=0.991, 3 sizes 1,232,634..12,011,994 |

verdict=**pass**，gates_checked=7，gates_not_run=0；体积校验 6/6 通过（最差 -0.0170% @ 821,756 三角面，判据 3%）
<!-- perf-evidence:fleet3pt-gates end -->

<!-- perf-evidence:fleet3pt-table begin -->
| 子任务 | 体数 | 三角面 | 秒 | 吞吐 | 建索引 s | 索引 MiB | 活 RSS MiB | 峰值 RSS MiB |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| load | 1 | 410,878 | 18.443 | 22,278 tri/s | — | — | 160 | 619 |
| load | 1 | 1,640,958 | 87.160 | 18,827 tri/s | — | — | 447 | 2,339 |
| load | 1 | 4,003,998 | 269.335 | 14,866 tri/s | — | — | 1,014 | 5,735 |
| fleet_assemble | 2 | 821,756 | 0.584 | 1,406,197 tri/s | — | — | 312 | 619 |
| fleet_assemble | 2 | 3,281,916 | 3.069 | 1,069,466 tri/s | — | — | 1,091 | 2,339 |
| fleet_assemble | 2 | 8,007,996 | 6.221 | 1,287,164 tri/s | — | — | 2,656 | 5,735 |
| fleet_assemble | 3 | 1,232,634 | 0.932 | 1,322,555 tri/s | — | — | 423 | 619 |
| fleet_assemble | 3 | 4,922,874 | 3.818 | 1,289,406 tri/s | — | — | 1,470 | 2,339 |
| fleet_assemble | 3 | 12,011,994 | 9.534 | 1,259,876 tri/s | — | — | 3,499 | 5,735 |
| **point_locate_accel** | 2 | 821,756 | 0.636 | 3,145 query/s | 6.97 | 62.9 | 384 | 619 |
| **point_locate_accel** | 2 | 3,281,916 | 0.945 | 2,117 query/s | 30.75 | 249.8 | 1,315 | 2,339 |
| **point_locate_accel** | 2 | 8,007,996 | 1.195 | 1,674 query/s | 63.55 | 607.3 | 3,228 | 5,735 |
| **point_locate_accel** | 3 | 1,232,634 | 0.645 | 3,100 query/s | 9.88 | 94.3 | 500 | 619 |
| **point_locate_accel** | 3 | 4,922,874 | 1.019 | 1,963 query/s | 39.63 | 374.2 | 1,805 | 2,339 |
| **point_locate_accel** | 3 | 12,011,994 | 1.130 | 1,769 query/s | 98.16 | 910.2 | 4,380 | 5,735 |
| volume | 2 | 821,756 | 0.385 | 2,132,084 tri/s | — | — | 328 | 619 |
| volume | 2 | 3,281,916 | 1.729 | 1,898,042 tri/s | — | — | 1,091 | 2,339 |
| volume | 2 | 8,007,996 | 3.928 | 2,038,436 tri/s | — | — | 2,653 | 5,735 |
| volume | 3 | 1,232,634 | 0.600 | 2,056,020 tri/s | — | — | 416 | 619 |
| volume | 3 | 4,922,874 | 2.840 | 1,733,440 tri/s | — | — | 1,467 | 2,339 |
| volume | 3 | 12,011,994 | 5.652 | 2,125,364 tri/s | — | — | 3,492 | 5,735 |
<!-- perf-evidence:fleet3pt-table end -->

**结论（全部可由上表复算）**

1. **十个数量级没走形**：`load` 指数 1.17、`fleet_assemble` 1.05(x2)/1.02(x3)、`volume` 1.03(x2)/1.00(x3)，
   全部低于 1.50 上限且 R²≥0.99；`point_locate_accel` 查询指数 0.28(x2)/0.25(x3)，远低于 0.80 上限——
   **10⁶→10⁷ 带内查询时延近常数**（12.01e6 tri 的 2000 条射线 1.130 s，与 0.82e6 档 0.945 s 同量级）。
   这正是 §5「点定位的集群需求是算法缺口」在 10⁷ 规模上的实测复核。
2. **10⁷ 带首次实测的关键数字**：12,011,994 tri — 体积 5.652 s（相对解析 −0.003%，判据 3%）、
   装配 9.534 s、建索引 98.16 s / 910 MiB、活 RSS 3.50 GB、网格对象 2.61 GB。
   单体生成上限 8,199,998 tri 另由 `cj3563` 实测：load 573.5 s、体积 3.363 s（−0.003%）。
3. **`§8.3` 的 `not_run` 门被本节取代**：那批门 `not_run` 的原因是阶梯形状（每分组 1 个规模点，
   §9.1 已解释），不是性能问题；本节 3 规模点/分组，门可判定。

**诚实边界**

- 每个分组只有 **3 个规模点**，指数是无置信区间的三点最小二乘；最弱的是
  `point_locate_accel[x3]`（R²=0.951），因为查询时延几乎不随规模变化，拟合的"斜率"接近噪声。
  因此只声明"**没有观察到超限增长**"，不声明精确指数值。
- `10⁷` 是**多体装配**规模（3 个互不相交副本），不是"一个 10⁷ 三角面网格"；装配开销单列
  `fleet_assemble`（9.534 s，指数 ~1.0，纯拷贝）。导入式大网格仍未实测（§6 第 1 条）。
- 绝对秒数含共享节点争用（fleet 作业运行窗口 load1=90.8/128）；门是**同作业内跨档比值**，
  对整机负载不敏感，但不同作业的秒数不可逐档相减。
- 体积基准是解析闭式解（`copies × 7.512698 m³`），不是实验值；它校验的是装配/积分不引入几何误差。

### 9.5 测试

`tests/test_perf_benchmark.py` 新增 4 项（该文件 37 项本机全通过）：三基础档 × 单一体数 → 每个
`scaling_*` 门 `status=measured` 且 `gates_not_run == 0`；多体扶手上 `point_locate`/`voxelize`
记 `not_run` 且带原因、无 `failed` 行、每个不适用子任务只记一条；采样内存估算 vs 逐元素精确核算
相差 <5%；`rss_mb`/`mesh_mb` 有值且进 JSON。
全量 `pytest tests/ -q` 另以 `cj3515` 递交集群（同样在排队），**本单不声称全量通过**。

---

## 10. 证据包可复核化（本版新增；把「表格」从人工转录变成工件渲染）

### 10.1 审计发现（对抗式复核，2026-11-15）

对着 §1 宣称的"全部来自 `artifacts/perf_scale/full_ladder.json`，程序化读取"逐条核对，得到三个
**证据约束失效**的问题（都属于"文档声称比证据强"）：

| # | 发现 | 严重度 | 证据 |
|---|---|---|---|
| 1 | 文档引用的 `artifacts/perf_scale/*.json` 与 `benchmarks/perf_scale/figures/*.png` **不在仓库里**（合并门禁会滤 `artifacts/`），干净检出后 §3 的 44 条测量**无法审计** | 高 | `ls artifacts/perf_scale` → No such file；`git log -- docs/perf_scale_10m.md` 只有文档提交 |
| 2 | §3 表格是**手抄**到 markdown 的数字，没有任何机制阻止它与工件不一致（"程序化读取"只在作者机器上成立） | 高 | 表与工件之间无引用关系、无测试 |
| 3 | 墙上时钟回归门 `evaluate_gates()` **不核对基线环境**：换架构/换机器后，2× 门仍会拿别人的秒数判 RED 或判绿 | 中 | 基线 `environment{machine,cpu_count}` 已记录但从未被读 |

### 10.2 修复（本单落地）

1. **工件受版本控制**：测量 JSON 落在 `benchmarks/perf_scale/artifacts/{ladder_mesh,ladder_voxel}.json`
   （`benchmarks/` 不被门禁滤），任何人在干净检出上都能核对表格里的每个数字。
2. **表格程序化注入**：`tensorcad/perf_evidence.py` 把工件渲染成 §3.1/§3.2/§9.4/§10.4 的 markdown 表，
   表格夹在 `<!-- perf-evidence:<name> begin/end -->` 标记之间；`scripts/perf_scale_doc.py --write`
   重渲染、`--check` 校验。`tests/test_perf_scale_evidence.py` 断言
   **提交的文档 == 提交的工件渲染结果**，手改数字即测试失败（漂移不再可能静默发生）。
   本版把同一绑定扩到 10⁶–10⁷ 阶梯（`fleet3pt-table`、`fleet3pt-gates` 两张新表）：新增结论
   （七个指数门、体积校验）也必须是工件渲染出来的，而不是从 JSON 里抄进 markdown。
3. **基线环境闸门**：`environment_compat()` + `evaluate_gates()` 在 `machine`/`platform` 不一致时
   拒绝比较，给出 `regression_environment` 门且 `status=not_run`（**不是 pass**），
   `cpu_count`/`python` 差异只写进 detail；无反证测试不会有人发现这条从未触发，
   `tests/test_perf_benchmark.py::test_cross_machine_baseline_refuses_comparison` 用 x86_64 基线
   断言"不产生任何 `regression_load*` 绿门"。

### 10.3 本版实测范围（诚实边界）

- 本版 §3 表格改由**全阶梯**工件渲染：`cj3563`（n330，`21×9 … 4001×513`，最高 **8,199,998 三角面**）
  与 `cj3564`（`32³ … 384³`，最高 **56,623,104 cells**），两部分都在提交的
  `benchmarks/perf_scale/artifacts/*.json` 里。上一版只能覆盖到 1,640,958 tri / 128³，是欠账；
  现在两个欠账作业都已到件，`scripts/perf_scale_doc.py --write` 无损换表完成。
- **跨机秒数**：n315（上一版，load1≈40）与 n330（本版，load1 47–51）同档可对照但不可逐档相减：
  `load` @4.0e6 由 213.7 s → 219.2 s（+2.6%），属共享机噪声量级；结论只用同一工件内部的比值
  （如 §3.1 加速度比），跨工件只比数量级与指数。
- **10⁷ 带**的唯一实测来自 `cj3514` 三点 fleet 阶梯（12,011,994 tri，含七个指数门），
  §8.3 的 `cj3081` 单点数据保留为历史记录，其 `not_run` 门由 §9.4 取代。
- 回归基线 `baseline.json` 仍是 quick 阶梯（≤2.6×10⁴ tri）；10⁶–10⁷ 带的守护是
  **指数门 + 本表的工件绑定**，不是 2× 墙钟门——后者在共享机上噪声大于信号。

### 10.4 60 s 预算等价规模（§5 的「什么规模需要集群」的可复现版本）

规则写死、可复算：**取该子任务实测最大档的实测速率 × 60 s**，不做任何拟合外推。
读法：等价规模超过"实测最大档"即为超出本机 60 s 预算的区间，需要分批或上集群；
速率与档位同表给出，读者可自行换预算（`BUDGET_SECONDS`）。

<!-- perf-evidence:budget-table begin -->
| 子任务 | 规模轴 | 实测最大档 | 秒 | 实测速率 | 60 s 预算等价规模 | 峰值 RSS MiB |
|---|---|---:|---:|---:|---:|---:|
| load | 三角面 | 8,199,998 | 573.483 | 14,299 tri/s | 857,916 | 12,746 |
| point_locate | 三角面 | 410,878 | 235.551 | 8 query/s | 510 | 643 |
| point_locate_accel | 三角面 | 8,199,998 | 0.965 | 2,072 query/s | 124,338 | 12,746 |
| volume | 三角面 | 8,199,998 | 3.363 | 2,438,286 tri/s | 146,297,166 | 12,746 |
| voxelize | 体素 | 56,623,104 | 356.209 | 158,960 cell/s | 9,537,612 | 794 |
<!-- perf-evidence:budget-table end -->

## 11. 从文件导入大网格（OBJ 往返）——§6 第 1 条的剩余欠账

### 11.1 为什么要拆成两个子任务

「导入一个 10⁷ 三角的网格要多久」不是可操作的问题：**序列化**（内存网格 → 文本写盘）与
**解析**（读入文本 → float → 扇化面索引）是两条不同的代码路径、不同量级的成本。
因此它们是 `tensorcad.perf_benchmark` 的两个一等子任务：

```bash
python -m tensorcad.perf_benchmark --subtasks load,obj_write,obj_read,volume \
    --ladder "4001,513"            # 单体 8.2e6 三角
python -m tensorcad.perf_benchmark --fleet-tier --ladder "641,161;1001,257;2001,501" \
    --fleet-copies 2,3 --subtasks load,fleet_assemble,volume,obj_write,obj_read
```

因为进了同一张子任务表，它们**自动继承**两者：

1. **复杂度指数门** `scaling_obj_write` / `scaling_obj_read`（上限 1.50，O(N) 路径的判定口径）。
   把 `read_obj` 的逐行解析改成 O(N²)（例如每面重扫顶点表），小阶梯上指数就会跳出来。
2. **2× 墙钟回归门**：quick 阶梯（≤2.6×10⁴ tri）的 8 个 `obj_*` 基线键已写入
   `benchmarks/perf_scale/baseline.json`（`scripts/build_perf_baseline.py --merge` 只新增键，
   旧键一秒不改），任何一次提交都能在秒级成本下判 RED。

### 11.2 口径（怎么才算一次正当测量）

- **落盘位置**：`$TMPDIR/tensorcad-perf-import/hull_<nx>x<nz>_b<copies>.obj`。一档一个文件，写下一档时
  删上一档，所以 10⁷ 档只占约 1 GiB、且**永远不会进仓库**（`benchmarks/` 之外的体积不进版本控制）。
  写前按 `OBJ_BYTES_PER_TRIANGLE` 估算并要求 **4× 空闲空间**；空间不足记 `not_run` 并附原因，
  不写半个文件。
- **自证**：`obj_read` 解析回来的顶点数/面数必须与内存网格完全一致，并把散度定理体积对解析值
  `copies × 7.512698 m³` 的相对误差记成 `volume_rel_err`。**快而丢面的解析是缺陷，不是测量**：
  计数不符直接抛 `ValueError`（`failed` 行），而不是给出一个好看的秒数。
- **无损**：`write_obj` 写 17 位有效数字，所以往返是逐坐标精确的
  （`tests/test_perf_import.py` 断言 `parsed.vertices == mesh.vertices` 且体积完全相等，
  不是「在容差内」）。

### 11.3 实测（工件渲染，非手抄）

<!-- perf-evidence:import-table begin -->
| 子任务 | 体数 | 三角面 | 秒 | 吞吐 tri/s | 文件 MiB | 文件 MiB/s | 文件 B/tri | 峰值 RSS MiB | 导入体积误差 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| load | 1 | 25,918 | 0.751 | 34,492 | — | — | — | 72 | — |
| load | 1 | 410,878 | 15.406 | 26,671 | — | — | — | 629 | — |
| load | 1 | 1,640,958 | 110.528 | 14,847 | — | — | — | 2,422 | — |
| load | 1 | 4,003,998 | 283.716 | 14,113 | — | — | — | 5,751 | — |
| obj_write | 1 | 25,918 | 0.057 | 454,683 | 1.0 | 18.0 | 41.4 | 72 | — |
| obj_write | 1 | 410,878 | 1.187 | 346,078 | 18.2 | 15.3 | 46.4 | 629 | — |
| obj_write | 1 | 1,640,958 | 5.372 | 305,477 | 75.4 | 14.0 | 48.2 | 2,422 | — |
| obj_write | 1 | 4,003,998 | 11.013 | 363,585 | 206.3 | 18.7 | 54.0 | 5,751 | — |
| obj_read | 1 | 25,918 | 0.093 | 278,628 | 1.0 | 11.0 | 41.4 | 72 | -1.359e-01 % |
| obj_read | 1 | 410,878 | 1.676 | 245,172 | 18.2 | 10.9 | 46.4 | 629 | -1.673e-02 % |
| obj_read | 1 | 1,640,958 | 7.586 | 216,328 | 75.4 | 9.9 | 48.2 | 2,422 | -5.870e-03 % |
| obj_read | 1 | 4,003,998 | 17.095 | 234,217 | 206.3 | 12.1 | 54.0 | 5,751 | -2.990e-03 % |
| volume | 1 | 25,918 | 0.010 | 2,520,007 | — | — | — | 72 | — |
| volume | 1 | 410,878 | 0.176 | 2,340,416 | — | — | — | 629 | — |
| volume | 1 | 1,640,958 | 0.751 | 2,185,953 | — | — | — | 2,422 | — |
| volume | 1 | 4,003,998 | 1.877 | 2,133,446 | — | — | — | 5,751 | — |
<!-- perf-evidence:import-table end -->

读法（数字全部取自上面的渲染表，读者可复算）：

- **写比读快，且差距随规模放大**：同一文件两列 `MiB/s` 之比就是读写成本比。这不是 Python 的问题而是
  两条路径的问题——写是 `"\n".join` 一批已格式化行，读是整文件 `read_text` + `splitlines` + 逐字段
  `float` + 面索引扇化。
- **文件体积是线性量，但 B/tri 不是常数**：每三角面约一行 `v` 加一行 `f`（17 位有效数字），
  所以比值随网格的「顶点:三角面」比例漂移——读**表内 `文件 B/tri` 列**（该列由 `file_bytes / triangles`
  现算，不再用散文中手写的常数）。本节 4 档是 41.4→54.0 B/tri，§11.5 的 1.2×10⁷ 档是
  54.6 B/tri：**「≈48 B/tri」只是中间档的巧合**，10⁷ 档按 54.6 B/tri 计才准（625.5 MiB）。
  写盘/读取的 I/O 不是瓶颈，解析才是。
- **60 s 预算等价规模**：按上表最大档实测速率线性折算（不做拟合外推），`obj_read` 的 60 s 预算落在
  **10⁷ 量级**——即"从文件导入一个 10⁷ 三角的网格"在本机上是**分钟级**作业，不是必须上集群的作业。
  与 `point_locate` 的暴力路径（60 s 只够 5×10² 三角）相比，导入反而是整个基准里最省的一项。
- **导入体积误差是网格本身的离散误差，不是解析误差**：往返是 17 位无损的，所以 `obj_read` 行的误差
  与同档 `volume` 行同值——两者一致本身就是无损往返的证据。

**诚实边界**

- **本节这张表本身只覆盖到 4.00×10⁶ 三角**（4 档，一次 `--ladder` 调用）。它曾在这里留过一条欠账：
  单体 8.2×10⁶（`cj3891`）与三体 1.2×10⁷（`cj3890`）两个作业当时排队 35 min 未出队，本节不填外推数字。
  **该欠账已清零**：两个作业事后都 `succeeded`（rc=0），工件受版本控制并渲染为 **§11.5**——
  本节保持 4 档阶梯的独立性（同一台机、同一次调用、可算指数），10⁶–10⁷ 带的数字去 §11.5 读。
- 本节 4 档 → §11.4 的四个指数门**全部 `measured`**（不是 `not_run`）。注意 `not_run` 的口径在这里仍然有效：
  §11.5.1（单体 8.2×10⁶，**只有一个规模点**）的四个门就是 `not_run`，因为 1 点拟合不出斜率——
  同一套机制、同一份工件，点数不够时打印 `not_run` 而不是"通过"。
- 秒数含共享机负载（本档运行时 load1≈47/128），跨档只能比数量级与比值；读写是单进程单线程，
  不受集群核数影响。

### 11.4 门（`not_run` 不是 pass）

<!-- perf-evidence:import-gates begin -->
| 门 | 状态 | 判定 | 证据 |
|---|---|---|---|
| `scaling_load` | measured | pass | exponent=1.19 (limit 1.50), R^2=0.997, 4 sizes 25,918..4,003,998 |
| `scaling_obj_write` | measured | pass | exponent=1.06 (limit 1.50), R^2=0.998, 4 sizes 25,918..4,003,998 |
| `scaling_obj_read` | measured | pass | exponent=1.04 (limit 1.50), R^2=1.000, 4 sizes 25,918..4,003,998 |
| `scaling_volume` | measured | pass | exponent=1.03 (limit 1.50), R^2=1.000, 4 sizes 25,918..4,003,998 |
<!-- perf-evidence:import-gates end -->

### 11.5 10⁶–10⁷ 带实测（`cj3891` 单体 8.2×10⁶ / `cj3890` 三体 1.2×10⁷）—— §11.3 欠账清零

两个作业在 §11.3 写下欠账时都还在排队；事后两个都 `succeeded`（rc=0），工件经 fetch 取回后与**作业日志的
人类可读表**逐项对得上（两条独立通道，见 `tensorcad/perf_scale_artifact.py` 的设计理由），现受版本控制：

| 工件 | 作业 | 节点 | 规模 | 幂等校验 |
|---|---|---|---:|---|
| `artifacts/ladder_import_single_10m.json` | `cj3891` | n330 | 8,199,998 tri（单体上限 4001×513） | `…single_10m.provenance.json`（sha256 + rc + wall 872.4 s） |
| `artifacts/ladder_import_fleet_10m.json` | `cj3890` | n315 | 12,011,994 tri（3 × 2001×501） | `…fleet_10m.provenance.json`（sha256 + rc + wall 537.8 s） |

**为什么必须分两个工件**：单体网格受 `HullParameters` 上限约束（`nz ≤ 513`、`nx ≤ 4001`）最大
8,199,998 三角面，**单网格到不了 10⁷**；10⁷ 带只能由 `build_fleet_mesh` 的多体装配到达。因此
「导入一个 10⁷ 三角的网格」正确拆成「单体上限档」与「多体 10⁷ 档」两次测量，而不是把一个 8.2e6 的数
外推成 10⁷ 的说法。

#### 11.5.1 单体上限档（`cj3891`，8,199,998 三角面）

<!-- perf-evidence:import-single10m-table begin -->
| 子任务 | 体数 | 三角面 | 秒 | 吞吐 tri/s | 文件 MiB | 文件 MiB/s | 文件 B/tri | 峰值 RSS MiB | 导入体积误差 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| load | 1 | 8,199,998 | 810.207 | 10,121 | — | — | — | 11,778 | — |
| obj_write | 1 | 8,199,998 | 19.651 | 417,282 | 424.9 | 21.6 | 54.3 | 11,778 | — |
| obj_read | 1 | 8,199,998 | 34.581 | 237,125 | 424.9 | 12.3 | 54.3 | 11,778 | -2.737e-03 % |
| volume | 1 | 8,199,998 | 3.431 | 2,389,769 | — | — | — | 11,778 | — |
<!-- perf-evidence:import-single10m-table end -->

门（注意：**单点阶梯给不出斜率，4 个门全是 `not_run`，不是 pass**）：

<!-- perf-evidence:import-single10m-gates begin -->
| 门 | 状态 | 判定 | 证据 |
|---|---|---|---|
| `scaling_load` | not_run | not_run | not_run: 1 distinct size(s), need 3 to fit an exponent |
| `scaling_obj_read` | not_run | not_run | not_run: 1 distinct size(s), need 3 to fit an exponent |
| `scaling_obj_write` | not_run | not_run | not_run: 1 distinct size(s), need 3 to fit an exponent |
| `scaling_volume` | not_run | not_run | not_run: 1 distinct size(s), need 3 to fit an exponent |
<!-- perf-evidence:import-single10m-gates end -->

#### 11.5.2 多体 10⁷ 档（`cj3890`，2 体与 3 体，最高 12,011,994 三角面）

行先列单体 `load`（造数成本，与往返无关），再按体数分组（2 体 → 3 体），组内按子任务；同一文件由 `obj_write` 写、`obj_read` 读，两行的
`文件 MiB` 应相同——表格里能直接看出来。

<!-- perf-evidence:import-fleet-table begin -->
| 子任务 | 体数 | 三角面 | 秒 | 吞吐 tri/s | 文件 MiB | 文件 MiB/s | 文件 B/tri | 峰值 RSS MiB | 导入体积误差 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| load | 1 | 410,878 | 17.770 | 23,123 | — | — | — | 615 | — |
| load | 1 | 1,025,998 | 46.129 | 22,242 | — | — | — | 1,504 | — |
| load | 1 | 4,003,998 | 243.799 | 16,423 | — | — | — | 5,753 | — |
| fleet_assemble | 2 | 821,756 | 0.597 | 1,377,226 | — | — | — | 615 | — |
| fleet_assemble | 2 | 2,051,996 | 1.529 | 1,342,421 | — | — | — | 1,504 | — |
| fleet_assemble | 2 | 8,007,996 | 5.723 | 1,399,369 | — | — | — | 5,753 | — |
| obj_write | 2 | 821,756 | 1.874 | 438,508 | 36.7 | 19.6 | 46.8 | 615 | — |
| obj_write | 2 | 2,051,996 | 4.829 | 424,909 | 100.3 | 20.8 | 51.3 | 1,504 | — |
| obj_write | 2 | 8,007,996 | 18.705 | 428,111 | 415.9 | 22.2 | 54.5 | 5,753 | — |
| obj_read | 2 | 821,756 | 3.517 | 233,623 | 36.7 | 10.4 | 46.8 | 615 | -1.673e-02 % |
| obj_read | 2 | 2,051,996 | 8.636 | 237,606 | 100.3 | 11.6 | 51.3 | 1,504 | -8.253e-03 % |
| obj_read | 2 | 8,007,996 | 34.168 | 234,370 | 415.9 | 12.2 | 54.5 | 5,753 | -2.990e-03 % |
| volume | 2 | 821,756 | 0.371 | 2,212,913 | — | — | — | 615 | — |
| volume | 2 | 2,051,996 | 0.984 | 2,085,340 | — | — | — | 1,504 | — |
| volume | 2 | 8,007,996 | 3.787 | 2,114,338 | — | — | — | 5,753 | — |
| fleet_assemble | 3 | 1,232,634 | 0.911 | 1,352,348 | — | — | — | 615 | — |
| fleet_assemble | 3 | 3,077,994 | 2.287 | 1,345,718 | — | — | — | 1,504 | — |
| fleet_assemble | 3 | 12,011,994 | 8.664 | 1,386,358 | — | — | — | 5,753 | — |
| obj_write | 3 | 1,232,634 | 2.767 | 445,448 | 55.2 | 20.0 | 47.0 | 682 | — |
| obj_write | 3 | 3,077,994 | 7.245 | 424,872 | 153.3 | 21.2 | 52.2 | 1,648 | — |
| obj_write | 3 | 12,011,994 | 28.880 | 415,924 | 625.5 | 21.7 | 54.6 | 6,277 | — |
| obj_read | 3 | 1,232,634 | 4.968 | 248,129 | 55.2 | 11.1 | 47.0 | 853 | -1.673e-02 % |
| obj_read | 3 | 3,077,994 | 13.198 | 233,224 | 153.3 | 11.6 | 52.2 | 2,047 | -8.253e-03 % |
| obj_read | 3 | 12,011,994 | 50.481 | 237,949 | 625.5 | 12.4 | 54.6 | 7,881 | -2.990e-03 % |
| volume | 3 | 1,232,634 | 0.621 | 1,983,662 | — | — | — | 615 | — |
| volume | 3 | 3,077,994 | 1.494 | 2,060,384 | — | — | — | 1,504 | — |
| volume | 3 | 12,011,994 | 5.878 | 2,043,502 | — | — | — | 5,753 | — |
<!-- perf-evidence:import-fleet-table end -->

门（8 个规模点，2 体/3 体各一组，`[x2]`/`[x3]` 后缀是体数）：

<!-- perf-evidence:import-fleet-gates begin -->
| 门 | 状态 | 判定 | 证据 |
|---|---|---|---|
| `scaling_fleet_assemble[x2]` | measured | pass | exponent=0.99 (limit 1.50), R^2=1.000, 3 sizes 821,756..8,007,996 |
| `scaling_fleet_assemble[x3]` | measured | pass | exponent=0.99 (limit 1.50), R^2=1.000, 3 sizes 1,232,634..12,011,994 |
| `scaling_load` | measured | pass | exponent=1.16 (limit 1.50), R^2=0.998, 3 sizes 410,878..4,003,998 |
| `scaling_obj_read[x2]` | measured | pass | exponent=1.00 (limit 1.50), R^2=1.000, 3 sizes 821,756..8,007,996 |
| `scaling_obj_read[x3]` | measured | pass | exponent=1.02 (limit 1.50), R^2=1.000, 3 sizes 1,232,634..12,011,994 |
| `scaling_obj_write[x2]` | measured | pass | exponent=1.01 (limit 1.50), R^2=1.000, 3 sizes 821,756..8,007,996 |
| `scaling_obj_write[x3]` | measured | pass | exponent=1.03 (limit 1.50), R^2=1.000, 3 sizes 1,232,634..12,011,994 |
| `scaling_volume[x2]` | measured | pass | exponent=1.02 (limit 1.50), R^2=1.000, 3 sizes 821,756..8,007,996 |
| `scaling_volume[x3]` | measured | pass | exponent=0.99 (limit 1.50), R^2=1.000, 3 sizes 1,232,634..12,011,994 |
<!-- perf-evidence:import-fleet-gates end -->

#### 11.5.3 读法

- **导入 10⁷ 三角是分钟级作业，不是集群作业**：三体 1.2×10⁷ 的 `obj_read` 是 **50.5 s**（625.5 MiB，
  12.4 MiB/s），`obj_write` 28.9 s；把 §10.4 的 60 s 预算口径套上来，**导入 10⁷ 不需要集群**
  （生成 1.2×10⁷ 网格的 244 s `load` + 8.7 s 装配才是成本，但那是一次性造数，不是查询/导入成本）。
- **指数门在 10⁷ 带是 `measured` 且 O(N)**：`scaling_obj_read[x2]`=1.00、`[x3]`=1.02，
  `scaling_obj_write[x2]`=1.01、`[x3]`=1.03（上限 1.50，R²≥1.000）。这三组数字是"解析路径没有退化"
  的可判定证据，而不是"跑起来没崩"。
- **体积误差是网格离散误差，且随规模收敛**：三体 1.2×10⁷ 档 `volume=22.537421 m³` vs 解析
  `3 × 7.512698 = 22.538095 m³`，相对误差 **−3.0×10⁻⁵（0.003%）**；单体 8.2×10⁶ 档
  `volume_rel_err = −2.74×10⁻⁵`。这不是"OBJ 往返失真"，恰恰相反：往返是 17 位无损的，误差与
  同档 `volume` 行同源（几何离散），两者一致正是无损的证据。
- **写仍比读快约 1.75×**：1.2×10⁷ 档 21.7 vs 12.4 MiB/s。与 §11.3 小档的 1.5× 相比差异缩小——
  I/O 占比上升，解析的固定开销被摊薄。

**诚实边界（本节）**

- **跨机不可比**：单体档在 n330、多体档在 n315，两次运行的 load1 分别是 41.7 与 46.1（128 核）。
  因此两个工件之间只比数量级与同工件内的比值，**不把 8.2e6 的秒数与 1.2e7 的秒数放在一条拟合线上**；
  每组的指数门都在组内（同一台机、同一进程）拟合。
- **单体档的门是 `not_run`**：8.2×10⁶ 只有 1 个规模点，`scaling_*` 拒绝拟合（`MIN_SCALING_POINTS=3`）。
  要把它变成 `measured` 需要再跑 2 个小档（成本约 150 s 生成 + 60 s 往返），本单没跑——
  但 §11.3 的 4 点阶梯 + 本节的 8 点阶梯已经覆盖了同一条代码路径的 O(N) 判定，此处不重复计费。
- **10⁷ 是"多体装配"规模**：单个**连通** 10⁷ 三角面体的独立证据来自 `ladder_single_10m.json`
  （`docs/perf_scale_single_body_10m.md`，中点细分几何），与本节的舰队装配是两条不同的到达路径，
  不要混引。
- 三体 1.2×10⁷ 的 OBJ 文件 **625.5 MiB 不落仓库**（写 `$TMPDIR`，一档写完删上一档）；
  仓库里只有 JSON 与 provenance，任何读者可按 §11.1 的命令重跑复现。
**工件（受版本控制、可直接复核）**: `benchmarks/perf_scale/artifacts/{ladder_mesh.json,ladder_voxel.json}`、
`benchmarks/perf_scale/{baseline.json,ladder_fleet_10m.json,ladder_3pt_mechanism_smoke.json}`；
图 `benchmarks/perf_scale/figures/*.png`；渲染/校验 `tensorcad/perf_evidence.py` + `scripts/perf_scale_doc.py`
`benchmarks/perf_scale/figures/scaling.png`（墙钟 vs 规模 log-log）、`point_locate.png`（两条点定位路径吞吐对比）、`index_cost.png`（建索引时间/索引体积）、`voxelize.png`（体素化缩放）；10⁶–10⁷ 带另有两张：`fig3pt_scaling.png`（每个子任务的拟合指数与门限，§9.4）与 `fig3pt_memory.png`（网格占用/索引/活内存 vs 规模）。
## 4. 误差判定（关键量 vs 基准）
| 图 | `benchmarks/perf_scale/figures/{scaling,point_locate,index_cost,voxelize}.png` |
| 作图脚本 | `scripts/gen_perf_scale_figures.py` |
| 子任务 | 体数 | 三角面 | 秒 | 吞吐 tri/s | 文件 MiB | 文件 MiB/s | 峰值 RSS MiB | 导入体积误差 |
| load | 1 | 25,918 | 0.751 | 34,492 | — | — | 72 | — |
| load | 1 | 410,878 | 15.406 | 26,671 | — | — | 629 | — |
| load | 1 | 1,640,958 | 110.528 | 14,847 | — | — | 2,422 | — |
| load | 1 | 4,003,998 | 283.716 | 14,113 | — | — | 5,751 | — |
| obj_write | 1 | 25,918 | 0.057 | 454,683 | 1.0 | 18.0 | 72 | — |
| obj_write | 1 | 410,878 | 1.187 | 346,078 | 18.2 | 15.3 | 629 | — |
| obj_write | 1 | 1,640,958 | 5.372 | 305,477 | 75.4 | 14.0 | 2,422 | — |
| obj_write | 1 | 4,003,998 | 11.013 | 363,585 | 206.3 | 18.7 | 5,751 | — |
| obj_read | 1 | 25,918 | 0.093 | 278,628 | 1.0 | 11.0 | 72 | -1.359e-01 % |
| obj_read | 1 | 410,878 | 1.676 | 245,172 | 18.2 | 10.9 | 629 | -1.673e-02 % |
| obj_read | 1 | 1,640,958 | 7.586 | 216,328 | 75.4 | 9.9 | 2,422 | -5.870e-03 % |
| obj_read | 1 | 4,003,998 | 17.095 | 234,217 | 206.3 | 12.1 | 5,751 | -2.990e-03 % |
| volume | 1 | 25,918 | 0.010 | 2,520,007 | — | — | 72 | — |
| volume | 1 | 410,878 | 0.176 | 2,340,416 | — | — | 629 | — |
| volume | 1 | 1,640,958 | 0.751 | 2,185,953 | — | — | 2,422 | — |
| volume | 1 | 4,003,998 | 1.877 | 2,133,446 | — | — | 5,751 | — |
- **文件体积是线性量**：`文件 MiB / 三角面` ≈ 48 B/tri（17 位有效数字的 `v` + `f` 行），10⁷ 档约 0.5 GiB，
- 表内最大档是 **1.64×10⁶ 三角**。10⁶→10⁷ 带的 OBJ 往返实测（单体 8.2×10⁶ 的 `cj3891`、三体 1.2×10⁷ 的
  `cj3890`）递交后排队 **35 min 仍未出队**，本节**不填外推数字**：上表结论只覆盖到 1.64×10⁶ 三角。
  两个作业到件后再 `python scripts/perf_scale_doc.py --write` 无损换表即可（表格与工件绑定，不会静默漂移）。
- 本节阶梯只有 **2 档**，所以 §11.4 的复杂度指数门 `scaling_obj_write`/`scaling_obj_read` 是
  **`not_run`（不是 pass）**：2 点拟合给不出可判定的斜率，硬给一个数就是编数字。指数门机制本身由
  `tests/test_perf_import.py::test_import_subtasks_get_measured_scaling_gates` 在三点阶梯上覆盖
  （断言 `status=measured` 且指数低于 `SCALING_LIMITS`）；10⁶–10⁷ 档的指数门由 `cj3890` 补齐。
