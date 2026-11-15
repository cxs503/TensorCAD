# GAP-CAD-T2-PARAM-SERIES 验收裁定：可执行、带证据分级

**任务** m22471（selfstart, deepseek-v4-flash 席）· 生成于 2026-11-15T22:00:17Z
**锚点** `GAP-CAD-T2-PARAM-SERIES` · **裁定** `incomplete`

## ① 案例介绍

Series 60（Todd 1951）是单桨商船的标准系统船模系列，母型 Cb 阶梯 0.60/0.65/0.70/0.75/0.80。
TensorCAD 的 `tensorcad.series60` 以解析站型律
$y(\xi,\zeta)=(B/2)\,f_{\rm deck}(\xi)\,\zeta^{1/m}$ 生成闭式船体，再由 161×41 网格积分静水力。
本锚点的台账验收共三句：**五档 Cb 生成 + 静水力对照公开值；gate-power 突变测试 ≥80%；文档**。

主尺度：L = 10.0 m，B = 1.5 m，T = 0.6 m，
网格 161×41，水密度 1025 kg/m³。

**为什么本单不是"重做"**：逐句核对后，三句验收在前序任务里事实上已落地（五档生成 + 系数对照、
14/14 突变、五份文档齐备），台账状态 `open` 是滞后。所以本单做的是**验证欠账**：把这三句从
"只能靠读散文"变成**机器可裁定**，并修掉审计（m17308 F1/F2/F3/F5）指出的判级缺陷。

## ② 怎么做的

新增 `tools/series60_family/anchor_acceptance.py`：把锚点验收拆成 5 条子句，每条给出
`pass / fail / not_run / withheld` 与**证据等级**（`external` / `closure` / `withheld`，封闭词表）：

1. `geometry_precondition_cb_l_b` — 先验 Cb/L/B（平台纪律：几何未验不得引用下游系数）。
   L/B/T 按构造精确，门 0.1%；**Cb 由网格体积积分而来，带网格首阶亏损，门取 2%**
   （与模块自陈的离散化门一致——本单实测发现：把 Cb 也压到 0.1% 会让每一条收敛的船体都"失败"）。
2. `five_cb_generation_and_published_hydrostatics` — 五档生成 + 闭固性 + 五系数对照公开值，门 3%。
3. `gate_power_mutation_kill_rate` — 读**机器可读**突变记录；记录缺失即 `not_run`（=未知，绝不是"通过"）。
4. `documentation` — 五份参考页存在且非空。
5. `external_shape_anchor_todd_offsets` — 外部型值锚，**`withheld`**（见 ⑤）。

全部数字程序化注入，公开表 `docs/data/series60_published_parents.json`（sha256 `16ac564989a87d2b…`）
只读不抄写。核心自洽检查：`closure_where_external_should_be = []`，
即**没有任何闭合对照被报成外部验证**——这正是审计 F1 的根因，现在由测试而非注释守住。

## ③ 计算结果

### 3.1 静水力 vs 公开母型表（门 3%）

| Cb | 量 | 模型 | 公开值 | 相对误差 | 等级 | ≤3% |
|---|---|---|---|---|---|---|
| 0.60 | block_coefficient | 0.595651 | 0.600000 | -0.725% | closure | ✅ |
| 0.60 | lcb_over_length | 0.475084 | 0.475000 | +0.018% | closure | ✅ |
| 0.60 | midship_coefficient | 0.950068 | 0.960000 | -1.035% | closure | ✅ |
| 0.60 | prismatic_coefficient | 0.626957 | 0.625000 | +0.313% | closure | ✅ |
| 0.60 | waterplane_coefficient | 0.699968 | 0.700000 | -0.005% | closure | ✅ |
| 0.65 | block_coefficient | 0.644631 | 0.650000 | -0.826% | closure | ✅ |
| 0.65 | lcb_over_length | 0.483067 | 0.483000 | +0.014% | closure | ✅ |
| 0.65 | midship_coefficient | 0.959603 | 0.970000 | -1.072% | closure | ✅ |
| 0.65 | prismatic_coefficient | 0.671768 | 0.670000 | +0.264% | closure | ✅ |
| 0.65 | waterplane_coefficient | 0.726965 | 0.727000 | -0.005% | closure | ✅ |
| 0.70 | block_coefficient | 0.693424 | 0.700000 | -0.939% | closure | ✅ |
| 0.70 | lcb_over_length | 0.490041 | 0.490000 | +0.008% | closure | ✅ |
| 0.70 | midship_coefficient | 0.967917 | 0.980000 | -1.233% | closure | ✅ |
| 0.70 | prismatic_coefficient | 0.716409 | 0.715000 | +0.197% | closure | ✅ |
| 0.70 | waterplane_coefficient | 0.752960 | 0.753000 | -0.005% | closure | ✅ |
| 0.75 | block_coefficient | 0.742023 | 0.750000 | -1.064% | closure | ✅ |
| 0.75 | lcb_over_length | 0.496014 | 0.496000 | +0.003% | closure | ✅ |
| 0.75 | midship_coefficient | 0.975230 | 0.985000 | -0.992% | closure | ✅ |
| 0.75 | prismatic_coefficient | 0.760869 | 0.760000 | +0.114% | closure | ✅ |
| 0.75 | waterplane_coefficient | 0.779954 | 0.780000 | -0.006% | closure | ✅ |
| 0.80 | block_coefficient | 0.790335 | 0.800000 | -1.208% | closure | ✅ |
| 0.80 | lcb_over_length | 0.501998 | 0.502000 | -0.000% | closure | ✅ |
| 0.80 | midship_coefficient | 0.981713 | 0.990000 | -0.837% | closure | ✅ |
| 0.80 | prismatic_coefficient | 0.805057 | 0.805000 | +0.007% | closure | ✅ |
| 0.80 | waterplane_coefficient | 0.806946 | 0.807000 | -0.007% | closure | ✅ |

最坏一行：**+1.233%**（midship_coefficient @ Cb=0.70），
在 3% 门内，且非零——门是真在咬合的。

> **等级为 `closure`**：Cb/Cwp/LCB 是形状求解的**输入**，Cp/Cm 由公开 Cp 阶梯派生。
> 这组一致性度量的是代数，不是与真实船型的一致。注意 Cm 支路——Cb=Cp·Cm 在公开表里本身
> 不自洽（Cb=0.80 处 0.8/0.805=0.9938 vs 表列 0.990），所以"+0.38%"类偏差量的是公开表的
> 内部不自洽，不是实现精度。

### 3.2 几何前置门（Cb/L/B）

| Cb | 量 | 实测 | 请求 | 相对 | 门 | 结论 |
|---|---|---|---|---|---|---|
| 0.60 | length | 10.000000 | 10.000000 | +0.0000% | 0.0010 | ✅ |
| 0.60 | beam | 1.500000 | 1.500000 | +0.0000% | 0.0010 | ✅ |
| 0.60 | draft | 0.600000 | 0.600000 | +0.0000% | 0.0010 | ✅ |
| 0.60 | block_coefficient | 0.595651 | 0.600000 | -0.7248% | 0.0200 | ✅ |
| 0.65 | length | 10.000000 | 10.000000 | +0.0000% | 0.0010 | ✅ |
| 0.65 | beam | 1.500000 | 1.500000 | +0.0000% | 0.0010 | ✅ |
| 0.65 | draft | 0.600000 | 0.600000 | +0.0000% | 0.0010 | ✅ |
| 0.65 | block_coefficient | 0.644631 | 0.650000 | -0.8261% | 0.0200 | ✅ |
| 0.70 | length | 10.000000 | 10.000000 | +0.0000% | 0.0010 | ✅ |
| 0.70 | beam | 1.500000 | 1.500000 | +0.0000% | 0.0010 | ✅ |
| 0.70 | draft | 0.600000 | 0.600000 | +0.0000% | 0.0010 | ✅ |
| 0.70 | block_coefficient | 0.693424 | 0.700000 | -0.9394% | 0.0200 | ✅ |
| 0.75 | length | 10.000000 | 10.000000 | +0.0000% | 0.0010 | ✅ |
| 0.75 | beam | 1.500000 | 1.500000 | +0.0000% | 0.0010 | ✅ |
| 0.75 | draft | 0.600000 | 0.600000 | +0.0000% | 0.0010 | ✅ |
| 0.75 | block_coefficient | 0.742023 | 0.750000 | -1.0636% | 0.0200 | ✅ |
| 0.80 | length | 10.000000 | 10.000000 | +0.0000% | 0.0010 | ✅ |
| 0.80 | beam | 1.500000 | 1.500000 | +0.0000% | 0.0010 | ✅ |
| 0.80 | draft | 0.600000 | 0.600000 | +0.0000% | 0.0010 | ✅ |
| 0.80 | block_coefficient | 0.790335 | 0.800000 | -1.2081% | 0.0200 | ✅ |

### 3.3 突变门（gate-power）

- 声明判据：kill rate ≥ 80%
- 记录：`None` → 状态 **not_run**，kill_rate = None
- 说明：no mutation record at /srv/tensor/agent/work/m22471/artifacts/series60_gate_power/gate_power_series60.json; the kill rate is UNKNOWN, which is not the same as passed -- regenerate it with tools/hull_gate_power/mutate.py --only <S60-*>

## ④ 结果图

![anchor acceptance](artifacts/series60_anchor_acceptance/figures/anchor_acceptance.png)

左：五档 Cb 下五个系数相对公开值的偏差（红虚线为 ±3% 门；图例标注每条的证据等级）。
右：五条子句的状态裁定（绿 pass / 橙 not_run / 灰 withheld）——**橙色是未知，不是通过**。

## ⑤ 注意点（踩过的坑 / 局限 / 适用范围）

1. **外部型值锚 WITHHELD，这是诚实的边界，不是造假。** 沙箱内没有 Todd 原表，仓内从未落地过
   型值表：`tests/test_series60_stations_published_comparison.py` 比的是网格 vs 自身解析律，
   属于离散化自检。本族是"积分特性正确的解析替身"，**不重现 Todd 型值**，任何阻力结论不得
   建立在它上面。
2. **"输入目标"不可能被证伪。** 本单先写了一个"把公开 LCB 挪 0.02 应当失败"的负对照，
   结果它**通过**了——因为 LCB 是求解目标，表一动模型就跟着动。这恰好把 `closure` 这个等级
   钉死了，于是该"失败"的测试被改写成一条显式的**闭包性测试**（若哪天它开始失败，说明这些量
   真的变成了预测，等级必须升为 `external`）。
3. **能证伪的负对照选 Cm**：Cm 是公开表里有、而形状求解不消费的一列，篡改 10% 必然触发 3% 门。
4. **Cwp 不可达即拒绝**：把公开 Cwp 压 10% 会被 `Cb/(Cp-Cb)` 的可行域拒绝并抛 `ValueError`，
   不会静默夹取——这条也写成了测试。
5. **修掉审计 F5**：`tools/series60_family/audit.py --out /abs/path` 原来在算完全部数字后崩于
   `relative_to(ROOT)`；现改为 `display_path()`，仓内给相对路径、仓外给绝对路径，并有回归测试。
   审计 F1（role 误标）在当前代码中已修复，本单加测试防复发；F2（gate_power.log 与 json 相反）
   已由文档侧 14/14 + 合并 23/24 (96%) 收口。

## ⑥ 误差判定：是否 < 3%？

| 关键量 | 相对公开值最坏误差 | 门 | 是否 < 3% |
|---|---|---|---|
| Cb（五档） | +1.208% | 3% | ✅ |
| Cp（五档） | +0.313% | 3% | ✅ |
| Cm（五档） | +1.233% | 3% | ✅ |
| Cwp（五档） | +0.007% | 3% | ✅ |
| LCB/L（五档） | +0.018% | 3% | ✅ |
| 网格体积 vs 闭式 | +1.208% | 2% | ✅ |
| **站型型值 vs Todd 原表** | **无基准（WITHHELD）** | — | **无法回答** |

**回答**：五档系数与公开母型表的一致**全部 < 3%**（最坏 +1.233%），
但这是**闭合一致性**，不是外部基准验证。**形状对 Todd 型值这一关键量没有基准，因此不能声称
"< 3%"，等级 WITHHELD。** 整体裁定 `incomplete`：在突变记录到位前，gate-power 子句为
`not_run`（未知），故不签 pass。

## 复现

```bash
python3 tools/series60_family/anchor_acceptance.py
python3 tools/series60_family/anchor_acceptance.py --mutations <gate_power_series60.json>
python3 -m pytest tests/test_series60_anchor_acceptance.py -q   # 15 passed
```
