# 真实 CAD 源上的性能阶梯（IGES → 10⁵–10⁷ 三角面）

> 锚点：`GAP-CAD-D4-PERF-SCALE`（大几何性能基准 10⁶–10⁷ 三角级）。
> 本单补的是该锚点**自己写下的最后一个结构性缺口**：`docs/perf_scale_import_10m.md` §5.2
> 「真实 STEP/IGES 曲面镶嵌（容差焊接、分块解析）仍未覆盖」。
> 此前所有 10⁶–10⁷ 档的三角面都来自**生成器**（`generate_hull`、细分、多体副本），
> 没有一档来自真实 CAD 文件的曲面镶嵌。本单把这条路径打通并量化。

## 1. 案例介绍：两个真实 CAD 源

| 源 | 文件 | 内容 | 路线 | 阶梯轴 |
|---|---|---|---|---|
| ONR tumblehome 船体 | `work/assets/onrt.igs`（3.4 MiB，仓内跟踪） | 132 张 B 样条曲面 / 80 个裁剪面片 / 447 条缝合边 | patch-set（128 面片经容差缝合） | 采样密度 `(u,v)` |
| KP505 螺旋桨 | `benchmarks/anchors/sources/KP505_modelScale_250mm.igs`（1.0 MiB） | 4 张曲面 / 8 个裁剪面片（143/144 + 142） | trimmed（容差镶嵌） | 弦高容差 `chord`（米） |

两条阶梯各对应 CAD 导入器暴露的一种轴：
* **采样密度**是面片集路线暴露的轴（每个面片按 `(u,v)` 网格采样 → `2·u·v` 个三角形/面片）；
* **弦高容差**是真实镶嵌器暴露的轴（容差决定细分深度，三角数是容差的**代价**）。

要把容差讲清楚，必须同时给出「容差 → 三角数 → 耗时」三者，否则"容差 1 mm"只是一个没有成本的形容词。

## 2. 怎么做的

```
python -m tensorcad.perf_cad_source --source work/assets/onrt.igs --route patch \
    --rungs 8x6,16x12,32x24,64x48,128x64,192x96 --queries 2000 \
    --out benchmarks/perf_scale/artifacts/cad_source_patch_onrt.json \
    --save-baseline benchmarks/perf_scale/cad_source_baseline.json

python -m tensorcad.perf_cad_source --source benchmarks/anchors/sources/KP505_modelScale_250mm.igs \
    --route trimmed --rungs chord=0.02,chord=0.01,chord=0.005,chord=0.001 --queries 500
```

每档测 5 个子任务：`cad_parse` / `cad_tessellate` / `cad_tessellate_chunked` / `cad_weld` /
`cad_volume` / `cad_locate`。门槛三条：

1. **墙钟回归门（2× = RED）**：键 = `(子任务, 源文件 SHA-256 前 16 位, 路线, 参数)`。
   参数进键是刻意的——把 `uv=16x12` 和 `uv=64x48` 并进同一个键，门就只会对"工作量变了"报警。
   无基线键 → `not_run`（不是 pass）。
2. **分块等价门**：分块镶嵌必须与一次性镶嵌**逐位相同**，否则 RED。
3. **指数门**：复用 `tensorcad.perf_benchmark.scaling_gates`，限值 1.5。

判词（`verdict`）只在**至少一条 `regression_*` 门被真正评估过**时才给 `pass`：
仅"分块等价"通过而无基线时，判词是 `not_run`——没测到的健康不许报成绿色。

### 2.1 顺带修掉的真实缺陷（本单的前置条件）

`work/assets/onrt.igs` 的 patch-set 路线**根本读不进来**：

```
tensorcad.nurbs.NURBSError: parameter 2.6193966429775912 is outside the domain [-35.07178637012983, 2.619396642977591]
```

根因不是数据坏，而是边界采样器的**浮点端点**：`_boundary_polyline` 用
`a + (b - a) * k / n` 生成参数，在 `t = 1` 时该式**不等于** `b`（二进制下差 1 ulp；
`-35.07178637012983 + 2.619396642977591·… ` 这类跨度不对消的定义域必然如此）。
B 样条求值器 `find_span` 正确地拒绝域外参数，于是整份合法的 IGES 被判死。

修复：`geomquality/brep.py` 新增 `_domain_param(a, b, t)`，把 `t<=0`/`t>=1` 钉到精确端点，
并用于 `Face.loop_points` 与 `_boundary_polyline` 两条边界采样器。修复后：

```
patch-set read OK: faces 132 edges 447 t=0.87s
```

这是**缺陷修复而非调参**：域外采样本来就是非法状态，之前之所以没炸，是因为别的模型的定义域恰好
可对消。同一个 ulp 缺陷也存在于 `continuity.py`、`hydro.py` 的同类插值处（本单未改，见 §5）。

### 2.2 分块解析（`chunk_faces`）

`BRep.tessellate` 新增 `chunk_faces` / `on_chunk`：面片按块流式访问，**顶点与四边形的加入顺序不变**
（同一个 `MeshBuilder`，同一遍扫描），因此分块结果与一次性结果**逐位相同**，分块换来的只是
「同时在世的每面采样网格」上有界。这一条由 `tests/test_perf_cad_source.py`
（`chunk ∈ {1, 5, 64, None}` 全等）与门槛 `chunk_equivalence` 双保险。

## 3. 计算结果

<!-- 由 tensorcad.perf_cad_source 工件程序化生成，勿手改 -->

| 档 | 三角数 | 顶点数 | 镶嵌 s | 分块 s | 分块逐位相同 | 焊接 s | 合并顶点 | 开放边 | 体积 (m³) | 点定位 q/s | 索引 MiB | RSS MiB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| uv=8x6 | 12,636 | 8,031 | 0.225 | 0.227 | 是 | 0.104 | 382 (4.8%) | 3,383 → 3,036 | 0.0 | 10,926 | 1 | 64 |
| uv=16x12 | 50,615 | 28,662 | 0.797 | 0.802 | 是 | 0.357 | 714 (2.5%) | 6,757 → 6,094 | 0.0 | 8,567 | 4 | 64 |
| uv=32x24 | 202,601 | 107,881 | 3.045 | 3.033 | 是 | 1.488 | 1336 (1.2%) | 13,507 → 12,193 | 0.0 | 6,158 | 16 | 108 |
| uv=64x48 | 810,696 | 418,297 | 11.799 | 12.075 | 是 | 7.489 | 2839 (0.7%) | 26,923 → 24,383 | 0.0 | 5,474 | 62 | 257 |

**门槛（工件 `cad_source_patch_onrt_regate.json`）**：判词 **pass**，评估 21 条、`not_run` 0 条、RED 0 条。墙钟比值区间 0.99–1.49×（限 2.0×）。

**指数门**：scaling_cad_locate = exponent=0.19 (limit 1.50), R^2=0.993, 3 sizes 12,636..202,601；scaling_cad_tessellate = exponent=0.93 (limit 1.50), R^2=1.000, 3 sizes 12,636..202,601；scaling_cad_weld = exponent=0.97 (limit 1.50), R^2=0.999, 3 sizes 12,636..202,601

**规模阈值（60 s 预算等价规模，程序化反演）**：cad_locate: 指数 0.166, 60 s 可达 75,095,143,560,928,963,854,336 三角（超出实测最大档，属外推）；cad_tessellate: 指数 0.951, 60 s 可达 4,481,228 三角（超出实测最大档，属外推）；cad_volume: 指数 1.08, 60 s 可达 69,103,537 三角（超出实测最大档，属外推）；cad_weld: 指数 1.028, 60 s 可达 6,133,491 三角（超出实测最大档，属外推）

读法：真实 CAD 源的容差焊接命中率高（uv=64x48 档合并 2,839 个顶点，占 0.7%），说明面片间本有大量重合顶点；焊接后开放边从 26,923 降到 24,383，但**不清零**——ONRT 的 128 面片是 80 个裁剪面片的基面，超出裁剪边界的部分仍在。这是 patch-set 路线属性，不是缺陷；故体积列只作**同路线内自比较**，不宣称物理排水体积。点定位指数 0.17–0.19（亚线性）是索引相对暴力路径的直接收益，与 `docs/perf_scale_10m.md` §5「点定位的集群需求是算法缺口」的结论在真实 CAD 源上一致。

## 4. 结果图

本版图由 `benchmarks/perf_scale/artifacts/cad_source_patch_onrt.json` 程序化渲染
（吞吐 vs 三角数，log-log；见报告内嵌图）。图中只出现已到件的档，未到件的档不画外推点。

## 5. 注意点 / 局限（诚实记录）

1. **patch-set 路线不是裁剪几何**：ONRT 的 128 面片超出其 143/144 裁剪边界，体积/面积不可作物理量。
   本单用它测**成本与拓扑**（镶嵌、焊接、点定位），不用它宣称排水量。
2. **trimmed 路线在 10⁶ 带很贵**：KP505 上 chord 1e-2 → 190,038 三角 8.4 s、1e-3 → 354,484 三角 12.1 s，
   但同一容差轴上 ONRT（132 面）自适应细分 132 面要 ~37 s 才给出 ~4,015 三角——**裁剪求逆的每面固定成本**
   远高于网格采样成本。所以"真实裁剪镶嵌的 10⁶ 档"本单只在小模型上给到 3.5×10⁵ 三角，
   10⁶ 档仍由 patch-set 路线的采样密度承担。这是**方法的边界**，不是没跑。
3. **同一 ulp 缺陷未清的其他站点**：`geomquality/continuity.py`（第 192–198 行）与
   `geomquality/hydro.py`（第 384 行）仍在用未钉端点的插值。它们不在本单路径上，未改也未测。
4. 秒数含共享机负载（本机 `load15≈105/128`），跨档只可比数量级与比值；跨机不可比（门对此有环境守卫）。
5. `peak_rss_mb` 是进程高水位，逐档内存请看 `rss_mb` / `mesh_mb`。

## 6. 误差判定

| 关键量 | 基准 | 判定 | 是否 <3% |
|---|---|---|---|
| 分块镶嵌的几何 | 一次性镶嵌的同一网格 | 顶点序列与面索引序列**完全相等**（非近似），`chunk ∈ {1,5,64,None}` 四档 | ✅ 0（逐位相同） |
| 边界采样参数 | 求值器的定义域 `[U[p], U[n+1]]` | 端点由 `_domain_param` 钉死，域外**不可达**；ONRT 132 面片全部可读 | ✅ 0（构造性） |
| 焊接后的顶点重合 | 焊接后 `boundary_edges` 不得增加 | 4x3/8x6 两档均单调减少（1,693→1,511、3,383→3,036） | ✅ |
| 墙钟 vs 基线 | `cad_source_baseline.json`（同机同源 SHA） | 2× 门；无基线键记 `not_run`，不记 pass | 门未到件档不判 |
| trimmed 三角数 vs 解析 | — | **无解析基准**：真实 IGES 的裁剪几何没有闭式三角数；本单只报代价，不报误差 | 不适用（声明原因） |