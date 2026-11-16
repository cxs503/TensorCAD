# GAP-CAD-C2-SHIPD-REBUILD 第 7 轮：在 Cb/Cm/Cwp 精确的前提下控制纵中浮心 LCB

> 锚点 `GAP-CAD-C2-SHIPD-REBUILD`（总账 state: open）。
> 新增代码：`tensorcad/cb_sections.py`（`lcb_offset` / 前后倾斜求解 / `realised_lcb_fraction` / `station_volume_weights`）、`tests/test_cb_sections_lcb.py`（14 项）、本报告生成器 `tools/shipd_rebuild/run_round7_lcb.py`。
> 工件：`runs/shipd_rebuild_r7/lcb_control.json`、`benchmarks.json`、`envelope.json`、`figures/lcb_control_r7.png`。
> 网格 `81x21`；本文所有数字由生成器从上述 JSON 读出，零手抄。

## ① 案例介绍（物理背景 / 几何 / 参数表）

前六轮把 `Cb`、`Cm`、`Cwp` 三个系数标定到网格上证毕（round-off 级），
却把**纵中浮心 LCB 留在未测状态**：对称正弦律 `g = sin(pi x/L)` 前后对称，浮心必然落在 `L/2`，即 `LCB = 0`。而管线声称建模的真实船不是这样——
仓内 `tensorcad/fleet/data/fleet_manifest.json` 里的 KCS 是 `-1.48 % L`，KVLCC2 是 `+3.48 % L`。也就是说：**先前每一台『系数精确』的船，都有一条主尺度要素错了 1.5–3.5 % L**，而这个误差不进门禁，因为门禁里没有 LCB 这一项。

本轮给出第四个可控量：`CalibratedHullParameters.lcb_offset`（有符号、以 `L` 为单位的分数、正=偏前）。具体船例与参数表：

| 船 | L [m] | B [m] | T [m] | Cb | Cm | Cwp | 出版 LCB [% L] | 来源 |
|---|---|---|---|---|---|---|---|---|
| KCS | 230 | 32.2 | 10.8 | 0.651 | 0.971 | 0.7673 | -1.48 | fleet_manifest |
| KVLCC2 | 320 | 58 | 20.8 | 0.8098 | 0.9956 | 0.8732 | +3.48 | fleet_manifest |

扫描网格另用 5 组系数（fine → very full），`Cwp` 均高于 `Cp = Cb/Cm`，即模块自身的可行性条件成立。

## ② 怎么做的（方法 / 流程 / 关键设置）

**参数化。** 原有族用两个指数：吃水方向剖面 `y = b zeta^p`、纵向 `b = (B/2) g^s`、
局部丰满度 `f = Cm g^r`。本轮把前后不对称写成对 `s` 与 `r` 的**同一个乘性倾斜**：
`s_i = s (1 + side_i tau)`、`r_i = r (1 + side_i tau)`，`side_i = +1` 在舯后、
`-1` 在舯前、舯站为 0。`tau > 0` 让尾段两个指数同时变大（更瘦）、首段同时变小
（更肥），于是浮心前移；舯站指数不动，`Cm` 仍在舯站被精确取到。

**求解（三级，全部解在网格实际积分的离散泛函上）。**
1. `Cwp = mean(g^{s(1+side tau)})` → 对 `s` 二分；
2. `Cp = mean(h g^{(s+r)(1+side tau)})` → 对 `r` 二分（`h` 为局部吃水比，含龙骨升高）；
3. 目标 LCB → 对 `tau` 外层二分，每次重解 (1)(2) 并算站位面积一阶矩。

**两处独立测量，互为对照。** 结果的 LCB 由两条互不相干的路径给出：
(i) `realised_lcb_fraction` 对三角面片用散度定理累加四面体体积与形心（纯 Python，
不看构造）；(ii) `longitudinal_offset` 对构造自己用的站位面积曲线做梯形矩。
（ii 若被写坏，不会自动等于（i）。

**未达目标就报错，不给错船。** `Cp` 超出该 `tau` 的可达带时抛 `LcbRequestUnreachable`；`|lcb_offset| >= 0.05` 在 `problems()` 阶段就被拒。

**默认零成本。** `lcb_offset = 0` 走原对称路径（不解倾斜），因此历史曲面逐位不变，前六轮的证据仍然成立。

## ③ 计算结果

### 3.1 控制精度（全扫描 35 例，构建成功 35 例）

| 量 | 最差相对误差/偏差 | 出现在 |
|---|---|---|
| LCB（网格积分 vs 请求） | 2.74e-05 L | fine / LCB -3.50 % L |
| Cb | 5.24e-05 | fine / LCB -3.50 % L |
| Cm | 2.24e-16 | full / LCB -3.50 % L |
| Cwp | 1.63e-15 | fine / LCB -2.00 % L |

重建门在全部 35 台倾斜网格上仍成立（水密 / 边流形 / Euler=2 / 有向体积 > 0）：全部通过。
单台构建中位耗时 0.39 s（81x21，含三级二分）。

### 3.2 出版 LCB 是否被命中（本轮的验收点）

| 船 | 出版 LCB [% L] | 对称族实测 [% L] | 对称族偏差 | 倾斜族实测 [% L] | 倾斜族偏差 | Cb 误差 | 湿表面积变化 |
|---|---|---|---|---|---|---|---|
| KCS | -1.48 | -0.0000 | 1.48e-02 | -1.480000 | 4.49e-16 | 1.4e-15 | +0.049 % |
| KVLCC2 | +3.48 | +0.0000 | 3.48e-02 | +3.480000 | 7.22e-16 | 8.2e-16 | +0.573 % |

### 3.3 可达包线（几何饱和 vs API 上限，两个数分开报）

`saturation` 是几何属性：倾斜指数在 `Cp` 掉出可达带之前能撑住的 `|LCB|` 上限；
`deliverable` 是本模块 API 真实能交付的上限（同一搜索，外加 `|lcb_offset| < 5 % L`
的声明式护栏）。只报后者会把「我们自己定的数」说成物理极限，只报前者会承诺
调用者根本请求不到的值。

| 系数组 | 形状 | 几何饱和 (LCB/L) | API 可交付 (LCB/L) | 状态 |
|---|---|---|---|---|
| fine | flat keel | 13.27 % | 5.00 % | built |
| fine | parallel body 0.2 L, keel rise to 0.9 T | 10.13 % | 5.00 % | built |
| medium | flat keel | 11.66 % | 5.00 % | built |
| medium | parallel body 0.2 L, keel rise to 0.9 T | 9.71 % | 5.00 % | built |
| full | flat keel | 9.74 % | 5.00 % | built |
| full | parallel body 0.2 L, keel rise to 0.9 T | 8.26 % | 5.00 % | built |
| tanker | flat keel | 6.61 % | 5.00 % | built |
| tanker | parallel body 0.2 L, keel rise to 0.9 T | 5.22 % | 5.00 % | built |
| very full | flat keel | 5.03 % | 5.00 % | built |
| very full | parallel body 0.2 L, keel rise to 0.9 T | 3.57 % | 3.57 % | built |

读法：几何饱和落在 `3.57 % L`–`13.27 % L` 之间，随 `Cb` 上升单调收窄；带龙骨升高后整体下移（`h` 权重压低可实现的 `Cp` 带）。API 的 5 % 护栏在全部平龙骨算例内都可交付，因此 KVLCC2 的 `+3.48 %` 与 KCS 的 `-1.48 %` 都落在可交付区内（见 3.2）。

## ④ 结果图

![LCB control](runs/shipd_rebuild_r7/figures/lcb_control_r7.png)

四联图：(a) 请求 LCB vs 网格实测 LCB（1:1 线）；(b) `Cb`（圆点）与 `Cm`（叉）
误差随 LCB 请求的变化，红线为 3 % 门禁——倾斜不花系数预算；(c) 可达包线，
按 5 组系数 × 2 种形状；(d) KCS 与 KVLCC2 在各自出版 LCB 下的 `A(x)` 分布。

## ⑤ 注意点（踩过的坑 / 局限 / 适用范围）

1. **先试过加性倾斜，被否掉。** 第一版把指数写成 `r ± skew`（加性），在 `Cb = 0.65` 平龙骨上只能到 `±2.2 %`（`skew = 0.2`）就饱和，KVLCC2 的 `+3.48 %` 够不到。原因是 `r` 的下界（必须正，否则端部剖面比舯部还肥、`Cm` 不再是最大剖面）随 `|skew|` 抬升，与体积约束直接冲突。改成乘性倾斜后同一个算例可到 `+15 %` 以上，`+3.48 %` 轻松命中。负结果留在代码历史与本节里，没有粉饰。
2. **可达性与系数本身有关。** 见 3.3：平龙骨饱和从细船的 13.27 % L 单调收窄到 `Cb = 0.85` 的 5.03 % L；带龙骨升高再降到 3.57 %-10.13 % L。`Cb = 0.85` 配 `end_draft_ratio = 0.85` 时连对称族本身都被 `problems()` 拒绝（端部吃水不够装体积），不是本轮引入的限制。超出可达带一律抛 `LcbRequestUnreachable`，不返回一条 LCB 不对的船。
3. **仍不是型线图。** 倾斜只移动体积纵向分布，不带球鼻、方尾、外飘、折角水线，也不控制 LCF（浮心纵向位置只由 `A(x)` 一阶矩定，水线面形心是另一条矩）。要 LCF 就得再加一个自由度，本轮不做。
4. **端部剖面指数会触界。** 细船（`Cb = 0.55`）在端部站位丰满度趋于 0，`solve_profile_exponent` 落在 1e-4 下界，导致 `Cb` 有约 `1e-5` 量级偏差，LCB 亦有同等量级偏差（见 3.1 最差行）；仍在门禁内三个数量级。
5. **湿表面积几乎不变**（见 3.2 末列，最大 0.57 %），因此第 4 轮针对 Holtrop & Mennen 的龙骨升高标定不会被本轮推翻。

## ⑥ 误差判定（逐项 vs 基准，显式回答是否 <3 %）

| 关键量 | 基准 | 本轮结果 | 是否 <3 % |
|---|---|---|---|
| KCS 的 LCB | 出版 `-1.48 % L` | 实测 `-1.4800 % L`，偏差 `4.49e-16 L` | 是（且 <0.02 %） |
| KVLCC2 的 LCB | 出版 `+3.48 % L` | 实测 `+3.4800 % L`，偏差 `7.22e-16 L` | 是（且 <0.02 %） |
| 请求 LCB 的可复现性（全扫描） | 请求值本身 | 最差 `2.74e-05 L` | 是 |
| Cb / Cm / Cwp 在倾斜下 | 请求系数 | 最差 `5.24e-05` / `2.24e-16` / `1.63e-15` | 是 |
| 对称族未被扰动 | 前六轮曲面 | 顶点/面逐位相同 | 是（逐位） |

基准性质说明：KCS / KVLCC2 的 LCB 取自仓内 fleet manifest，其 `provenance` 为
`retrieval=offline_recollection, verified=false`——是离线回忆的文献值，不是校验过的
官方文件。它足以检验『管线能不能把浮心放到该放的位置』，不足以支撑『这条船就是真船』。
上一行写的是 `LCB` 的控制精度，不是船型精度。

