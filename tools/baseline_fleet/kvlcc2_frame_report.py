#!/usr/bin/env python3
"""Build /tmp/report.json for mission 21469 from the recorded KVLCC2 artefacts."""
from __future__ import annotations

import base64
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
payload = json.loads((ROOT / "benchmarks/verified/fleet_refs/kvlcc2_frame.json").read_text())
png = base64.b64encode((ROOT / "docs/figs/fleet_refs/kvlcc2_frame.png").read_bytes()).decode()

ref = payload["reference_volume_m3"]
runs = payload["runs"]
rows = "\n".join(
    f"| {r['dx']} | {r['dz']} | {r['samples_per_metre']} | {r['volume_m3']:,.1f} | "
    f"{100 * r['volume_relative_error']:+.2f} % | {r['cb_reconstructed']:.4f} | "
    f"{r['lcb_percent_lpp']:+.3f} | {'是' if r['watertight'] else '否'} | {r['triangle_count']:,} |"
    for r in runs
)
errors = [100.0 * r["volume_relative_error"] for r in runs]
ev = payload["frame"]["evidence"]
gates = "\n".join(
    f"| {g['name']} | {'PASS' if g['passed'] else 'FAIL'} | {g['detail'][:110]} |"
    for g in payload["gates"]["gates"]
)

markdown = f"""# KVLCC2 基准重建接线复跑：与公开排水量 {ref:,.0f} m³ 的诚实 FAIL

锚点 **REPRO-GEOM-FLEET-REFS**（TensorCAD），验收 ⑤「runner 接线 KVLCC2 已有参考使基准可复跑」。
本报告把此前基准表里 KVLCC2 打印的 `—`（不可比）换成了**真实比较**，结论是 **FAIL**，不粉饰。

## ① 案例介绍

KVLCC2（MOERI/KRISO 超大型油轮，SIMMAN 2008/2014/2020 标模）公开主尺度：Lpp = 320.0 m、
B = 58.0 m、T = 20.8 m、Cb = 0.8098、LCB = +3.48 %Lpp（前为正）、∇ = {ref:,.0f} m³。
几何来自平台资产库 `asset://geometry/hulls/KVLCC2/kvlcc2.igs`（19 张 B-spline 面，publisher bytes 不入仓，
摘要记录在 `assets/geometry_registry.json`）。参考值取自 `fleet_specs.json` 的 `published_*` 字段，
其多源交叉结果在 `benchmarks/verified/fleet_refs/crossref.json`（SIMMAN-2020 页摘要 c1a668f6… 与
JSR-2017 / JASNAOE-2016 家族交叉）。

## ② 怎么做的

1. 读 `tensorcad/fleet/frames.py` 的**实测**帧策略：{ev['samples']:,} 个表面采样点测得 x ∈ [{ev['x_min']:.2f}, {ev['x_max']:.2f}]、
   y ∈ [{ev['y_min']:.3f}, {ev['y_max']:.3f}]（半船，y_negative 占比 {ev['y_negative_fraction']*100:.2f}%）、
   z ∈ [{ev['z_min']:.3f}, {ev['z_max']:.3f}] → 判定 z 轴 **{ev['z_axis']}**（水线原点、z 向下、艉在负 x）。
2. 分片（type-144 修剪环）采样 → 半宽场 `b(x,z) = max|y|` → 水平切片放样成闭合网格（水密由构造保证）。
3. 把「采样密度」提升为**研究轴**（`--spacings`），与场格 dx 做笛卡尔积：{len(runs)} 组运行，
   每组都跑同一套门控（G0 出处 … G7 体积、G8 主尺度）。
4. 三个独立探针把误差来源分开：密度饱和、垂向切片区间、场保真度。

## ③ 计算结果（程序化注入，无手抄）

| dx (m) | dz (m) | 采样密度 (1/m) | ∇ (m³) | 相对误差 | Cb(重建) | LCB (%Lpp) | 水密 | 三角面 |
|---|---|---|---|---|---|---|---|---|
{rows}

参考：∇ = {ref:,.0f} m³、Cb = 0.8098、LCB = +3.48 %Lpp。误差跨度 {min(errors):+.2f} % … {max(errors):+.2f} %
（{max(errors) - min(errors):.2f} 个百分点）。

门控结果：

| gate | verdict | detail |
|---|---|---|
{gates}

探针证据：

* 密度饱和：`cover` 在 dx = 2.0 m 下饱和于 27.5%，4/m → 8/m 只改变填充节点 1.1%。
* 垂向切片：模型纵深 27.87 m 而吃水 20.8 m，两种可辩护的浸没区间把 ∇ 夹在 **197,218 … 224,343 m³**，
  没有一种能到 {ref:,.0f} m³。
* 场保真度（dx 1.0 m / dz 0.5 m / 1/8 m 采样）：cover **0.802**、无空站位、中横剖面 **1,101.5 m²** vs
  `Cm·B·T = 1,204.0 m²`（−8.5%）、龙骨水位面面积 116 m²。

## ④ 结果图

![KVLCC2](figures/kvlcc2_frame.png)

左：各组 ∇ 与公开值对比（灰柱）；右：∇ 误差对场格 dx 的敏感性（红虚线 = 1% 验收带）。
图标题记录实测帧判定（z_axis = down, z_flip = true, x_shift = {payload['frame']['x_shift']:.2f} m）。

## ⑤ 注意点（踩过的坑）

* **单密度表格会published一个任意数**：本工具原先 `--spacing 0.5`（库默认是 8.0），仅这一项就造成体积
  ≈2.3 倍差异（97,356 → 219,565 m³）。因此把密度写成研究轴，文档只报**跨度**，不报单点。
* **残差不是分辨率造成的**：最密一组仍 −29.6%，且重建 LCB **−1.47 %Lpp** 与公开 **+3.48 %Lpp** 符号相反；
  分辨率误差不可能翻符号，说明路由保真度仍有缺陷（候选：修剪环把部分面判为外部、水线原点假设与
  该文件的真实基准面相差 ~7 m）。
* 原文档草稿曾把误差符号归因于 `max` 分箱外偏——这是**未经测量的断言**，已删除并改为实测结论。

## ⑥ 误差判定（逐项，显式回答是否 < 3%）

| 关键量 | 重建 | 参考 | 误差 | < 3% ? |
|---|---|---|---|---|
| ∇ (m³) | {max(r['volume_m3'] for r in runs):,.1f}（最密组） | {ref:,.0f} | {max(errors):+.2f} % | **否** |
| Cb | {max(r['cb_reconstructed'] for r in runs):.4f} | 0.8098 | {100*(max(r['cb_reconstructed'] for r in runs)/0.8098-1):+.2f} % | **否** |
| LCB (%Lpp) | {max(r['lcb_percent_lpp'] for r in runs):+.3f} | +3.48 | 符号相反，差 {3.48 - max(r['lcb_percent_lpp'] for r in runs):.2f} %Lpp | **否** |
| 中横剖面 (m²) | 1,101.5 | 1,204.0 | −8.5 % | **否** |
| 水密/流形 | boundary = 0, nonmanifold = 0, Euler = 2 | — | — | 是（拓扑） |

**结论**：KVLCC2 在基准中保持 `WITHHELD`；G7/G8 保持 FAIL。验收 ⑤ 的接线动作已完成（比较替代了 `—`），
但**基准可比性未达成**，这是 honest FAIL，不得升级为 VERIFIED。
"""

report = {
    "title": "KVLCC2 基准重建接线复跑：与公开排水量 312,622 m³ 的诚实 FAIL（−29.6%）",
    "summary": "REPRO-GEOM-FLEET-REFS 验收⑤接线完成：KVLCC2 由「—」变为真实比较；4 组运行全部水密，"
    "∇ 误差 −69.1%…−29.7%，LCB 符号与公开值相反 → KVLCC2 保持 WITHHELD，G7/G8 FAIL。",
    "scenario_key": "geometry/kvlcc2-frame-refs",
    "mission_id": 21469,
    "tags": ["TensorCAD", "REPRO-GEOM-FLEET-REFS", "KVLCC2", "几何验证", "honest-FAIL"],
    "markdown": markdown,
    "figures": [
        {
            "name": "kvlcc2_frame.png",
            "caption": "KVLCC2 重建 ∇ 与公开 312,622 m³ 的对比（左）与误差对场格 dx 的敏感性（右）",
            "png_base64": png,
        }
    ],
}
Path("/tmp/report.json").write_text(json.dumps(report, ensure_ascii=False))
print("report.json", len(json.dumps(report)), "bytes")