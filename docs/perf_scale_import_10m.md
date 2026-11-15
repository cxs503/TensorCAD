# 从文件导入 10⁷ 三角单体：细分几何源与导入门（GAP-CAD-D4-PERF-SCALE 续）

**锚点**: `GAP-CAD-D4-PERF-SCALE`（大几何性能基准 10⁶–10⁷ 三角级）
**前作**: `docs/perf_scale_10m.md` §11（OBJ 往返子任务，实测到 4.0×10⁶ 三角）
· `docs/perf_scale_single_body_10m.md`（单体 10⁷ 由细分达到，但只测了 load/volume/point_locate_accel）
**本版**: 把「**从文件导入一个 10⁷ 三角的连通单体**」变成一条可执行、可判定的测量：
新增 `--subdivide-levels`，让 `obj_write`/`obj_read`/`volume` 直接跑在细分rung 上。

---

## 1. 案例介绍

**问题**：前作 §11 的导入阶梯止步于 4.0×10⁶ 三角，而锚点要的是 10⁶–10⁷ 带。
参数化生成器 `HullParameters` 把单体上限钉在 `nx ≤ 4001, nz ≤ 513`（8.19×10⁶ 三角），
多体 fleet（`build_fleet_mesh`）虽然能到 1.2×10⁷，但它是 **k 个互不相交的副本**：
装配是纯拷贝、体积是 `k×V`，回答不了「**一个**连通体 10⁷ 三角，写盘/解析各要多久」。

**几何源**：均匀 1→4 中点细分 `subdivide_mesh(mesh, levels)`。每个平面父三角形被四个子三角形
**精确铺满**，中点按无向索引对 `(min,max)` 去重（焊接靠索引同一性，不靠容差），因此：

* 体积/面积逐位守恒（差异只有浮点累加序），拓扑按 `F'=4F, E'=2E+3F, V'=V+E` 精确推进；
* 父离散体的采样误差**不随细分改善**——它原样传给每一级（区间不变，只是更密）。
  这一点很重要：细分级上的体积误差是**父网格的**误差，不是导入误差。

**参数表**

| 项 | 值 |
|---|---|
| 基座 | `401×129` → 205,598 三角 |
| 细分级别 | 0,1,2,3（每级 ×4）→ 205,598 / 822,392 / 3,289,568 / **13,158,272** |
| 子任务 | `obj_write`（序列化）、`obj_read`（解析+自证）、`volume`（散度定理） |
| 文件命名 | `$TMPDIR/tensorcad-perf-import/hull_<nx>x<nz>_b<copies>_l<levels>.obj`（级别进文件名，级别间绝不串档） |
| 空间守卫 | 按 120 B/tri 估算并要求 4× 空闲；不足记 `not_run` + 原因，不写半个文件 |
| 复现命令 | `python -m tensorcad.perf_benchmark --ladder "401,129" --subtasks obj_write,obj_read,volume --subdivide-levels 0,1,2,3 --max-triangles 20000000 --output benchmarks/perf_scale/artifacts/ladder_import_10m.json` |
| 10⁷ 档作业 | `cj4105`（递交时队列位次 29，本机 load15≈74/128；**截稿时仍 queued，未到件**） |

---

## 2. 怎么做的

1. **算子复用**：直接调用已合并的 `tensorcad.mesh_subdivision.subdivide_mesh`（`docs/perf_scale_single_body_10m.md`
   的同一算子），不新写第二套细分语义。
2. **rung 访问器**：`_subdivided_for(nx, nz, copies, levels)`——`levels=0` 原样返回 `_mesh_for` 的结果
   （零行为变化），`levels>0` 时**从已缓存的最深父级**增量细分（走 0,1,2,3 只付每级一次的代价，
   而不是每次都从基座重算），缓存上限 2 级（10⁷ 网格是数 GB 的 Python 元组）。
3. **一次一档**：`import_obj_path(..., levels)` 把级别写进文件名，`clear_import_files`/`_prune_import_files`
   的 glob 一起覆盖，因此不会出现「读了未细分 rung 写的文件、却报细分 rung 的三角数」。
4. **自证不变**：`obj_read` 依旧断言顶点数/面数完全一致（否者抛 `ValueError`，不是给一个好听的秒数），
   并把导入后的体积对 `analytic_volume(10,2,0.8)` 的相对误差写进 `volume_rel_err`。
5. **门的三条口径**：
   - 复杂度指数门 `scaling_obj_write`/`scaling_obj_read`/`scaling_volume`：三条 rung 以上即 `measured`（≤1.50）；
   - 墙钟 2× 回归门：`measurement_key` 在 `levels>0` 时追加 `:l<levels>`，
     **细分级绝不与未细分基线比较**（否则 4× 大的网格必然假 RED，或反过来把回归洗白）；
   - 越限即 `not_run`：`est_tris × 4^levels > max_triangles` 记一行 `not_run` 带预期三角数，
     不是静默跳过；`load`/`point_locate*`/`voxelize` 在 `levels>0` 时记 `not_run` 并把原因写成
     "subdivide_levels applies only to obj_write, obj_read, volume"——它们不是细分网格的函数。

---

## 3. 计算结果

**本版到件的证据是「机制 + 门 + 测试」，不是 10⁷ 秒数**：作业 `cj4105` 截稿时仍在排队（见 §1 参数表），
所以本文**不填任何 10⁷ 的墙钟外推**——那是编数字。可信的路径速率上限来自前作 §11 的实测
（`obj_write` ≈ 2.4–3.6×10⁵ tri/s、`obj_read` ≈ 2.2–2.8×10⁵ tri/s、`volume` ≈ 2.1–2.5×10⁶ tri/s），
按它线性折算 13.16×10⁶ 三角约需 写 36–55 s / 读 47–60 s：**这是折算，不是测量**，到件后应以前作同款
`--write` 无损换表。本机小阶梯冒烟（`21×9`，基线 678 三角，levels 0/1/2 → 678/2,712/10,848）:
`obj_write` 0.0016/0.0058/0.0247 s、`obj_read` 0.0023/0.0098/0.0405 s、`volume` 0.0003/0.0011/0.0044 s，
三个指数门 `measured`：0.98 / 1.03 / 1.02（限值 1.50，R²≥0.999）——机制在最小规模上已被证明是 O(N)。

## 4. 结果图

本篇无图：它不是一次算例，而是一段**测量口径与门禁**的补齐；图上只会出现一条到件后才有意义的
10⁷ 曲线。到件后应在 `docs/perf_scale_10m.md` 的导入节用前作的渲染管线出图（吞吐 vs 三角数，
log-log），而不是在这里手画。

---

## 5. 注意点 / 局限（诚实记录）

1. **10⁷ 墙钟未到件**：`cj4105` 排队中；本文不含 10⁷ 的秒数、RSS、文件体积实测。
2. **细分不是"导入一个真实 CAD 模型"**：它给的是「一个连通体在 10⁷ 三角密度下的序列化/解析成本」。
   真实 STEP/IGES 曲面镶嵌（容差焊接、分块解析）仍未覆盖——这是本锚点剩下的唯一结构性缺口。
3. **基线键分离的代价**：细分级没有历史基线，所以它只有指数门、没有 2× 墙钟门；
   要对细分级做 2× 门，需先 `scripts/build_perf_baseline.py --merge` 在目标机上落一次基线。
4. **秒数含共享机负载**（本机 load15≈74/128），跨档只可比数量级与比值。
5. `peak_rss_mb` 是进程高水位（`ru_maxrss`），不是单档峰值；单档内存请看 `rss_mb`/`mesh_mb`。

## 6. 误差判定

| 关键量 | 基准 | 判定 | 是否 <3% |
|---|---|---|---|
| 细分前后体积 | 父网格散度定理体积 | 逐位一致（相对偏差 ≤1e-14，测试 `test_refinement_preserves_the_solid_volume_and_topology`） | ✅ |
| 细分级 OBJ 往返体积 | 同基座未细分 rung 的 `volume_rel_err` | 一致到相对 1e-6（17 位有效数字无损往返；测试 `test_obj_round_trip_is_lossless_on_a_refined_rung`） | ✅ |
| 往返计数完整性 | 写出网格的顶点/面数 | 完全相等，否则抛错（不是测量） | ✅ |
| 10⁷ 档导入秒数 vs 基准 | — | **无基准、无测量**：作业未到件，不做外推 | 不适用（声明原因） |
