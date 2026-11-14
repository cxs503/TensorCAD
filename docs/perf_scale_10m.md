# TensorCAD 大几何性能基准：10⁶–10⁷ 三角级实测（GAP-CAD-D4-PERF-SCALE）

**版本**: 1.0.0 · **锚点**: `GAP-CAD-D4-PERF-SCALE` · **任务**: m20350
**实测节点**: n315（LoongArch64, 128 核, Python 3.11.6, TensorCAD 0.23.0）
**作业**: `cj2652`（网格阶梯 load/point_locate/point_locate_accel/volume）、`cj2653`（体素阶梯）、`cj2669`（超限档探针，干净失败）
**工件**: `artifacts/perf_scale/{full_ladder.json,ladder_mesh.json,ladder_voxel.json,ladder_10m.json}`；图 `benchmarks/perf_scale/figures/*.png`；CI 基线 `benchmarks/perf_scale/baseline.json`

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
1. 加载/点定位/体素化/体积积分子任务基准表 → 44 条测量记录（0.64k → 4.0e6 三角，体素 3.3e4 → 5.7e7 cells）
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

## 3. 计算结果（全部来自 `artifacts/perf_scale/full_ladder.json`，程序化读取）

### 3.1 网格子任务

| 子任务 | 三角面 | 秒 | 吞吐 | 建索引 s | 索引 MB | 峰值 RSS MB |
|---|---:|---:|---:|---:|---:|---:|
| load | 678 | 0.019 | 34,962 tri/s | – | – | 36 |
| load | 25,918 | 0.813 | 31,887 tri/s | – | – | 76 |
| load | 410,878 | 15.58 | 26,370 tri/s | – | – | 640 |
| load | 1,640,958 | 71.28 | 23,020 tri/s | – | – | 2,438 |
| load | 4,003,998 | 213.70 | 18,737 tri/s | – | – | 6,151 |
| point_locate（暴力） | 678 | 0.257 | 7,771 query/s | – | – | 40 |
| point_locate（暴力） | 25,918 | 10.81 | 185 query/s | – | – | 76 |
| point_locate（暴力） | 103,038 | 60.39 | 33.1 query/s | – | – | 191 |
| point_locate（暴力） | 410,878 | 253.04 | 7.9 query/s | – | – | 640 |
| point_locate（暴力） | 1,640,958 | not_run | 超 5×10⁵ 阈值跳过 | – | – | – |
| **point_locate_accel** | 678 | 0.080 | 25,077 query/s | 0.006 | 0.06 | 40 |
| **point_locate_accel** | 25,918 | 0.295 | 6,784 query/s | 0.24 | 2.1 | 76 |
| **point_locate_accel** | 103,038 | 0.385 | 5,191 query/s | 0.81 | 8.4 | 191 |
| **point_locate_accel** | 410,878 | 0.783 | 2,553 query/s | 3.32 | 33.2 | 640 |
| **point_locate_accel** | 1,640,958 | 0.701 | 2,853 query/s | 12.59 | 131.4 | 2,438 |
| **point_locate_accel** | 4,003,998 | 1.090 | 1,836 query/s | 31.11 | 319.4 | 6,151 |
| volume | 103,038 | 0.045 | 2.30 M tri/s | – | – | 191 |
| volume | 1,640,958 | 0.670 | 2.45 M tri/s | – | – | 2,438 |
| volume | 4,003,998 | 1.644 | 2.44 M tri/s | – | – | 6,151 |

**同档加速比**（2000 条射线，同一起点/方向集）：25.9k 档 36.7×；103k 档 157×；411k 档 323×。10⁶ 档暴力路径已不可测（按 7.9 query/s 外推 2000 条射线 ≈ 4.2 min，1.64M 档 ≈ 17 min，4.0M 档 ≈ 42 min）。

### 3.2 体素化

| 网格 | 体素 | 秒 | 吞吐 cell/s | 峰值 RSS MB |
|---:|---:|---:|---:|---:|
| 32³ | 32,768 | 0.287 | 114,029 | 103 |
| 64³ | 262,144 | 2.444 | 107,273 | 560 |
| 128³ | 2,097,152 | 15.15 | 138,433 | 601 |
| 192³ | 7,077,888 | 48.02 | 147,382 | 640 |
| 256³ | 16,777,216 | 115.62 | 145,105 | 735 |
| 384³ | 56,623,104 | 366.12 | 154,658 | 795 |

### 3.3 图

`benchmarks/perf_scale/figures/scaling.png`（墙钟 vs 规模 log-log）、`point_locate.png`（两条点定位路径吞吐对比）、`index_cost.png`（建索引时间/索引体积）、`voxelize.png`（体素化缩放）。

---

## 4. 误差判定（关键量 vs 基准）

| 量 | 基准 | 实测 | 相对误差 | <3%? |
|---|---|---|---|---|
| 体积 @103,038 tri | 7.512698 m³（解析） | 7.509115 | −0.048% | 是 |
| 体积 @1,640,958 tri | 7.512698 m³ | 7.512257 | −0.0059% | 是 |
| 体积 @4,003,998 tri | 7.512698 m³ | 7.512474 | −0.0030% | 是 |
| `cast_ray` 交点距离 | 暴力 Möller–Trumbore（同网格） | 一致 | ≤1e-9 相对 | 是 |
| `contains` 内外判定 | 闭网格几何（Cb/体积符号） | 全对 | 0 错判 | 是 |
| 命中数（accel vs 暴力，2000 射线） | 154/185/218/253/286/317/346/363 | 逐档完全一致 | 0 | 是 |

**结论：关键几何量全部 <3%，点定位与参考实现逐位一致；加速结构不改变答案。**

---

## 5. 什么规模需要集群（实测判定，取代旧推测）

| 子任务 | 单机舒适区（实测） | 需要在集群上做 | 依据 |
|---|---|---|---|
| `load` | ≤ 4.0×10⁶ tri（214 s, 6.2 GB） | > 10⁷ tri 需先解决生成器限制（见 §6），或改用导入网格 | 吞吐 18.7k–35k tri/s 线性 |
| `point_locate`（暴力） | ≤ 5×10⁴ tri（2000 射线 ≤ 20 s） | ≥ 10⁵ tri | 33 query/s @103k, 7.9 query/s @411k |
| `point_locate_accel` | **≤ 4×10⁶ tri（2000 射线 1.1 s）** | 单机上 ≥10⁷ × 10⁵ 射线级别的批量查询 | 建索引 31 s / 319 MB @4.0e6，查询线性 |
| `volume` | ≤ 10⁷ tri（~4 s） | 不需要 | 2.4 M tri/s |
| `voxelize` | ≤ 5.7×10⁷ cells（366 s, 795 MB） | >10⁸ cells 或 dtype=float64 且内存紧张时 | 145k–155k cell/s |

**核心结论**：点定位的"集群需求"是**算法缺口**造成的，不是硬件缺口。同一个 411k 三角面的网格，暴力路径需要 253 s，加索引后 0.78 s（323×）——补上加速结构后，"10⁵ 三角面就该上集群"的判断被推翻，10⁶ 带单机即可。

---

## 6. 注意点 / 局限（诚实记录）

1. **10⁷ 档未能实测（honest FAIL）**：原计划用 `nz=801` 与 `nx=4884` 两档冲 10⁷，都被 `HullParameters` 的网格上限拒绝（`nz ≤ 513 levels`、`nx ≤ 4001 stations`）。参数化船体能构造的最大档是 **8,192,000 三角面**（4001×513）。因此本单实测上限是 4.0×10⁶，10⁷ 带只有外推。阶梯末档已改为 (4001, 513)（可构造），但**本轮未测**（约需 10 min 生成 + 1 min 建索引），是明确的待办，不是已完成项。
2. **越限档是干净失败**：`cj2669` 的两档返回 `status=failed` 并附生成器报错原话，没有静默跳过、也没有伪造数据；框架把异常记进结果表而非吞掉。
3. **暴力参考在 4.1×10⁵ 以上被主动跳过**（`max_naive_triangles`），标 `not_run` 而不是给一个不存在的秒数。
4. **回归基线的口径**：`benchmarks/perf_scale/baseline.json` 是共享机 3 轮取最慢，闸门阈值 2×；它只覆盖 quick 阶梯（≤2.6×10⁴ tri），10⁶ 档没有 CI 基线（成本太高），10⁶ 档的守护靠本单表格对比。
5. **负载影响**：n315 实测时 load1≈40（128 核），吞吐数字有 10–20% 抖动；同档之间比较应看数量级。
6. **`_RayCaster` 的角度前缀过滤**只对星形面（船体）成立，通用网格请用 `TriangleGrid` 或 `brute_force_nearest`。
7. **内存**：4.0×10⁶ 档峰值 6.2 GB（Python 对象为主）；索引额外 319 MB。10⁷ 档按线性外推需 ~15 GB，单节点内存不是瓶颈。

---

## 7. 工件索引

| 工件 | 路径 |
|---|---|
| 空间索引 | `tensorcad/spatial_index.py` |
| 索引测试（17 项，含与暴力实现逐位对比） | `tests/test_spatial_index.py` |
| 基准框架（含 `point_locate_accel`、网格缓存、峰值 RSS） | `tensorcad/perf_benchmark.py` |
| CI 回归基线（受版本控制） | `benchmarks/perf_scale/baseline.json` |
| 基线生成器 | `scripts/build_perf_baseline.py` |
| 原始测量 | `artifacts/perf_scale/{full_ladder,ladder_mesh,ladder_voxel,ladder_10m}.json` |
| 图 | `benchmarks/perf_scale/figures/{scaling,point_locate,index_cost,voxelize}.png` |
| 作图脚本 | `scripts/gen_perf_scale_figures.py` |
