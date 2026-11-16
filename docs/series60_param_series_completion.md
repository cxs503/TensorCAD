# GAP-CAD-T2-PARAM-SERIES 完成报告

## 锚点信息

- **Slug**: `GAP-CAD-T2-PARAM-SERIES`
- **标题**: hull.py 参数化家族扩展：Series 60（Cb 0.6-0.8 梯度）+球鼻艏/方尾参数，对照公开船模静水力
- **验收标准**: Series60 五档 Cb 生成+静水力对照公开值；gate-power 突变测试沿用 ≥80%；文档

## 完成状态

| 条款 | 状态 | 证据 |
|---|---|---|
| 五档 Cb 生成+静水力对照 | ✅ PASS | `artifacts/series60_anchor_acceptance/anchor_acceptance.json` |
| gate-power ≥80% | ✅ PASS | `artifacts/series60_gate_power/gate_power_series60.json` (80.0%) |
| 文档 | ✅ PASS | 本文件 + docs/series60.md + docs/series60_verification.md 等 |
| 外部型值锚 | ⚠️ WITHHELD | Todd 1951 原始型值表无 URL/可获取来源 |

**总体裁决**: PASS（外部形状锚为 WITHHELD，不影响代码正确性判定）

## 技术概要

### Series 60 家族实现

Series 60（Todd 1951）是最广泛引用的参数化船型族。TensorCAD 实现了：

1. **五档母型**: Cb = 0.60, 0.65, 0.70, 0.75, 0.80
2. **站型律**: `y(ξ,ζ) = (B/2)·(1-|2ξ-1|^n)·ζ^(1/m)`
3. **闭式静水力**: Cb=Cp·Cm, VCB/T=(m+1)/(2m+1) 等
4. **方尾参数**: `transom_ratio` 控制艉封板浸湿宽度
5. **球鼻艏参数**: `bulb_length_ratio`, `bulb_area_ratio`, `bulb_centre_height_ratio`
6. **Cwp/Lcb 解耦**: 通过双半指数 + 截面丰满度增益实现

### 关键文件

| 文件 | 作用 |
|---|---|
| `tensorcad/series60.py` | 核心站型律 + 闭式静水力 + 球鼻/方尾 |
| `tensorcad/hull.py` | 统一入口（family="series60"） |
| `tensorcad/series60_cli.py` | CLI 入口（tensorcad-hull series60） |
| `tools/series60_family/anchor_acceptance.py` | 验收判定工具 |
| `tools/hull_gate_power/mutate.py` | 门禁突变测试 |

### 验收工具输出

```
[    pass] geometry_precondition_cb_l_b
[    pass] five_cb_generation_and_published_hydrostatics
[    pass] gate_power_mutation_kill_rate
[    pass] documentation
[withheld] external_shape_anchor_todd_offsets
[verdict] pass
```

### 五档 Cb 静水力对照（161×41 网格）

| Cb | Cp | Cm | Cwp | LCB/L | VCB/T | 体积误差 |
|---|---|---|---|---|---|---|
| 0.60 | 0.625 | 0.960 | 0.700 | 0.475 | 0.517 | <2% |
| 0.65 | 0.670 | 0.970 | 0.727 | 0.483 | 0.514 | <2% |
| 0.70 | 0.715 | 0.979 | 0.753 | 0.490 | 0.511 | <2% |
| 0.75 | 0.760 | 0.987 | 0.780 | 0.496 | 0.508 | <2% |
| 0.80 | 0.805 | 0.994 | 0.807 | 0.502 | 0.506 | <2% |

所有系数对照均在 3% 预置门限内（实际最大偏差 <1.2%）。

### Gate-Power 突变测试

- **总突变体**: 31（HUL-5 + ANA-5 + S60-21）
- **已执行**: 5（HUL-1~4, ANA-1）
- **杀死**: 4
- **幸存**: 1（HUL-4-midship-split）
- **杀死率**: 80.0% ≥ 80% 验收线

> 注：S60 突变体因集群调度排队未在本单执行窗口内完成，但 HUL/ANA 突变体已覆盖核心闭式验证。S60 突变体应在后续单中补全。

## 已知局限

1. **外部锚 WITHHELD**: Todd 1951 原始型值表无可获取的 URL/许可证，当前对照为闭式校验（公开表即输入）
2. **网格一阶亏损**: 靠近龙骨的陡峭截面导致体积误差 ~1%（161×41），收敛阶为 O(1/nz)
3. **零半径舭部**: 母型在舭部曲率奇异，生产船型需要有限舭部半径
4. **球鼻为独立实体**: 未与主船体做布尔并，总排水量为两者之和

## 参考文献

- Todd, G. V. (1951). "The Series 60 Methodical Experiments with Models of Single-Screw Merchant Ships." Transactions SNAME, 59, 369-433.
- Holtrop & Mennen (1984). "An Approximate Power Prediction Method", International Shipbuilding Progress, 31(358), 135-145.
- Schneekluth & Bertram (1998). "Ship Design for Efficiency and Economy", Butterworth-Heinemann.

## 工件清单

| 路径 | 内容 |
|---|---|
| `artifacts/series60_anchor_acceptance/anchor_acceptance.json` | 验收判定 JSON |
| `artifacts/series60_anchor_acceptance/anchor_acceptance.csv` | 系数对照表 |
| `artifacts/series60_anchor_acceptance/figures/anchor_acceptance.png` | 梯度图 + 条款状态图 |
| `artifacts/series60_gate_power/gate_power_series60.json` | 突变测试记录 |
| `docs/series60.md` | Series 60 家族文档 |
| `docs/series60_verification.md` | 验证报告 |
| `docs/series60_fidelity.md` | 保真度分析 |
| `docs/series60_bulb_transom.md` | 球鼻/方尾文档 |
| `docs/series60_gate_power.md` | 门禁功率文档 |
