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

### 9.4 10^6–10^7 重算：已递交，收尾时仍在排队（诚实记录）

```
python scripts/perf_scale_10m.py --bases "641:161,1281:321,2001:501" --copies 2,3 \
    --n-queries 2000 --output benchmarks/perf_scale/ladder_fleet_10m_3pt.json
```

| 体数 | 三角面（3 点） |
|---:|---|
| 1（load） | 0.41e6 → 1.64e6 → 4.00e6 |
| 2 | 0.82e6 → 3.28e6 → 8.01e6 |
| 3 | 1.23e6 → 4.92e6 → 12.01e6 |

作业 `cj3514` 递交后 28 min 仍为 `queued`（当时集群 50 个作业在飞），**本单未取得 10^7 新数字**。因此：

- §8.3 的 10^7 表格**仍是唯一实测**，本单不覆盖、也不外推替代；
- 10^7 带的指数门**仍是欠账**（not_run），只因缺一次重算，不因缺机制；
- 重算到件后：`python scripts/perf_3pt_report.py benchmarks/perf_scale/ladder_fleet_10m_3pt.json`
  取表（零手抄）、`python scripts/gen_perf_3pt_figures.py` 出图，再把 9.1–9.3 的结论换成实测指数。

### 9.5 测试

`tests/test_perf_benchmark.py` 新增 4 项（该文件 37 项本机全通过）：三基础档 × 单一体数 → 每个
`scaling_*` 门 `status=measured` 且 `gates_not_run == 0`；多体扶手上 `point_locate`/`voxelize`
记 `not_run` 且带原因、无 `failed` 行、每个不适用子任务只记一条；采样内存估算 vs 逐元素精确核算
相差 <5%；`rss_mb`/`mesh_mb` 有值且进 JSON。
全量 `pytest tests/ -q` 另以 `cj3515` 递交集群（同样在排队），**本单不声称全量通过**。
