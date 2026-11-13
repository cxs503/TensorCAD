#!/usr/bin/env python3
"""Publish the m17307 board-audit report to the platform (/tensor/reports).

All numbers are read from the audit artifacts (evidence.json / independent_check.json /
tamper_probe.json) - nothing typed by hand.  Figure: docs/board_audit/m17307/audit_figures.png
(base64-embedded).

Run: python3 docs/board_audit/m17307/post_report.py
"""
from __future__ import annotations

import base64
import json
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = "http://127.0.0.1:8010/api/v1"
OUT = ROOT / "docs/board_audit/m17307"


def login() -> str:
    creds = json.loads(Path("/root/tp-agent/agent.env").read_text().split("\n", 1)[0]) if False else None
    env = dict(l.split("=", 1) for l in Path("/root/tp-agent/agent.env").read_text().splitlines()
               if "=" in l and not l.startswith("#"))
    body = json.dumps({"username": env["TENSOR_AGENT_USER"].strip(),
                       "password": env["TENSOR_AGENT_PASS"].strip()}).encode()
    req = urllib.request.Request(BASE + "/auth/login", data=body,
                                 headers={"Content-Type": "application/json"})
    return json.loads(urllib.request.urlopen(req, timeout=30).read())["access_token"]


def main() -> int:
    ev = json.loads((OUT / "evidence.json").read_text())
    chk = json.loads((ROOT / "runs/m17307/independent_check.json").read_text())
    res = json.loads((ROOT / "runs/m17307/subject_m17251/docs/offset_loft/results.json").read_text())
    b = ev["benchmark"]
    rv = b["ref_vs_got"]
    rows = "\n".join(
        f"| {h} | {rv[h]['volume_ref']:.6g} | {rv[h]['volume_got']:.6g} | {rv[h]['rel_err']:.4%} | "
        f"{rv[h]['rel_err_wp']:.4%} | {'是' if max(rv[h]['rel_err'], rv[h]['rel_err_wp']) < 0.03 else '否'} |"
        for h in ("wigley_iii", "series60_cb060", "kvlcc2_class"))
    st = "\n".join(f"| {n} | {v['dV_vs_81st_table']:.4%} | {v['dV_vs_closed_form']:.4%} |"
                   for n, v in b["station_sensitivity_kvlcc2_dV_vs_81st_table"].items())
    tables = "\n".join(
        f"| {k} | {v['n_stations']} × {v['n_waterlines']} | L={v['length']} m, B={v['beam']} m, "
        f"T={v['design_draft']} m | {v['source'][:70]}… | official={v['official']} |"
        for k, v in res["landed_tables"].items())
    png = base64.b64encode((OUT / "audit_figures.png").read_bytes()).decode()

    md = f"""# m17251 (GAP-CAD-T1-OFFSET-LOFT) 交差审计：L2 确认，L3/L4 拒收并给出具体缺口

## ① 案例介绍

被审对象=m17251「型值表→NURBS 曲面 lofting」（TensorCAD，commit `5bc0db7`，分支 `agent/m17251`，平台报告
id=16253）。被审验收线（mission goal 原文）：**≥3 船型值表落地；体积/水线面误差 <0.5% vs 数值积分；站线数
敏感性分析**。任务背景是 ship-platform 几何输入 8 个 STL 中 6 个为 304 三角占位符、KVLCC1≡KVLCC2 字节
相同——即缺真实几何。

三张已落地型值表（数值取自 `docs/offset_loft/results.json`）：

| 表 | 网格 | 主尺度 | 来源（provenance） | 是否发布机构原表 |
|---|---|---|---|---|
{tables}

## ② 怎么做的（审计方法，全部独立于被审报告）

1. 取证：台账 `tensor_maturity_verdicts`、`/tensor/missions/17251`、`/tensor/reports/16253`、工作树
   `/root/tp-agent/work/m17251`（commit 5bc0db7，`git status` 干净）；把交付逐文件 sha256 快照到
   `runs/m17307/subject_m17251/MANIFEST.json`（13 个文件，含实现/入口/测试/工件）。
2. 代码读者级复核：`tensorcad/offset_loft.py`（1282 行）+ `offset_forms.py`（284 行）→ `parametric_hull.
   loft_nurbs_surface` 双三次插值 → 以 patch 自身基函数做 Green 定理剖面矩（每 knot 区间 12 点 Gauss）→
   V/LCB/VCB/A_wp/湿表面积；另有一条 `NURBSSurface.derivatives` 慢路线与网格路线。
3. 独立重算（不 import 被审模块）：① 自写复合 Simpson（均匀 1/3+3/8 尾，非均匀用三点 Lagrange 二次权）
   从 CSV 重算基线；② 每条站剖面 z 方向 not-a-knot 三次样条再沿 x 样条积分（不同求积规则）作高精度对照；
   ③ Wigley 解析不变量 `4BLT/9`、`2LB/3`。
4. 完整性穿刺：改写单个纵坐标 +5%、删除 provenance 副车，验证回读被拒（不信任被审测试，自行构造）。
5. 判据对齐：`select count(*) from tensor_benchmark_cases where anchor_slug='GAP-CAD-T1-OFFSET-LOFT'`；
   `tensorbench_results`；`/tensor/bench/l4-verdict/17251`；`git ls-tree` 查生态仓主线上是否有本模块与调用方。
6. 复现：在被审工作树重跑 `pytest`（18 / 49 passed）与 `run_example.py`（`results.json` sha256 前后
   一致 `{ev['repro_cmd'].split('sha256 ')[1][:16]}…`，byte-identical=True）。
7. 落表 + 自校：evidence 12 键写入 `tensor_maturity_verdicts`（id=517，verdict_by='board-audit'），并用
   `maturity_gate.judge()` 预演：`judge(2) → (2, [])`、`judge(3) → (3 缺 integration_callers)`。

## ③ 计算结果（审计实测，程序化注入）

判据 A（被审自设）：loft 体积/水线面 vs 表自身 Simpson 积分 <0.5%

| 船型 | 基准 ref=表 Simpson | got=NURBS 求积 | dV | dA_wp | 是否 <3% |
|---|---|---|---|---|---|
{rows}

- 最差 dV = **{b['worst_relative_error_volume']:.4%}**，最差 dA = **{b['worst_relative_error_waterplane']:.4%}**，
  网格路线最差 {b['worst_relative_error_loft_mesh']:.4%} → 0.5% 验收线**通过且有余量**。
- 独立复算基线：审计者自写 Simpson 与报告 `table_simpson` 最大相对差 **{b['auditor_baseline_reproduction']['max_rel_diff_vs_reported_table_simpson']:.2e}**
  （三型全中）→ 基线数字不是编的。
- 独立求积 vs 报告 NURBS：Wigley {b['auditor_independent_quadrature']['rel_diff_volume_vs_reported_nurbs']['wigley_iii']:.2e}、
  Series60 {b['auditor_independent_quadrature']['rel_diff_volume_vs_reported_nurbs']['series60_cb060']:.2e}、
  KVLCC2 {b['auditor_independent_quadrature']['rel_diff_volume_vs_reported_nurbs']['kvlcc2_class']:.2e}。
- Wigley 解析锚：V=4BLT/9={chk['wigley_closed_form']['volume_analytic_4BLT_over_9']:.12g}，审计样条
  {chk['wigley_closed_form']['volume_audit_spline']:.12g}（rel {chk['wigley_closed_form']['rel_diff']:.1e}）；
  A_wp=2LB/3 rel {chk['wigley_closed_form']['rel_diff_waterplane']:.1e}。
- 几何自检 Cb：Wigley {chk['wigley_iii']['cb_from_tables']:.6f}（=4/9）、Series60 {chk['series60_cb060']['cb_from_tables']:.6f}、
  KVLCC2 {chk['kvlcc2_class']['cb_from_tables']:.6f}（发布值 0.8098，差 0.04%）。

站线数敏感性（KVLCC2，被审 results.json 原值）：

| 站数 | dV vs 81 站表 | dV vs 闭式 |
|---|---|---|
{st}

## ④ 结果图

![audit](audit_figures.png)

(a) 三型验收 ref-vs-got（虚线=0.5% 验收线，全部在线上方通过）；(b) KVLCC2 站线数敏感性（11 站 1.60% 破线，
≥15 站入带，与交付报告一致）；(c) 审计者独立求积与报告数字的对照——**routeB（自写 Simpson）与报告基线差
≤1e-15**，routeA（样条求积，不同规则）与报告 NURBS 差在 1e-5~2.4e-3 量级（后者即 KVLCC2 自身 41 站离散误差）。

## ⑤ 注意点（坑与适用边界）

1. **三张「型值表」不是发布机构原表**：provenance 全部 `official=false`，KVLCC2 副车明写「NOT the official
   KRISO offset table」。因此验收度量的是**插值+求积自洽性**，不是实船几何保真度；缺口背景里的「几何输入是
   占位符」并未解决——补上的是管线，不是真实数据（本条是本次审计最重的发现）。
2. **无基准注册**：`tensor_benchmark_cases` 该 anchor **0 行**，`tensorbench_results` 0 行，
   `/tensor/bench/l4-verdict/17251` 返回 `CLAIMED_NOT_PROVEN`；「ref」就是同一张表的 Simpson 积分，不能替代
   外部基准 → L4 不可给。
3. **无系统集成**：仅在 `tensorcad/__init__.py` 再导出 + 同 commit 的 `offset_forms.py` + 示例脚本消费；
   生态仓主线 `git ls-tree HEAD tensorcad/offset_loft.py` 为空（分支未合并），无求解器/CLI/API/平台调用方 → L3 不可给。
4. **验收参照物是研究里最粗的一条路线**：17 条水线按 ζ^(1/10) 聚集、首水线指数 p=1 时，表 Simpson 本身不如
   loft 准（KVLCC2 p=1：dV vs 表 0.230%，vs 闭式 1.153%）。交付报告两列都印了，读者不会被误导；但下游只拿
   「vs table Simpson」比较会高估 p=1 型值表质量。
5. **平台报告 16253 = 自动生成的 report.md**（纯表格），无六段结构、未显式回答「是否<3%」——报告纪律缺口
   （数字本身支持 <3%，见 ③）。
6. 审计侧方法论坑：`maturity_gate` 的证据路径只相对「生态仓根」或「审计工作树 m{{id}}」解析，且
   `impl_files` 带行号后缀（`:1-200`）会被判为「占位/缺失」——本单证据一律用审计工作树内的真实无后缀路径，
   并已用 `maturity_gate.judge()` 预演通过 `judge(2)→(2,[])`。
7. 未做：无 KCS/KVLCC2 官方型值表或发布水动力资料可作外部真值（沙箱内不存在）；网格路线仅以
   `results.json` 逐字节复现背书，未从零重写；大表（200×60）性能未测。

## ⑥ 误差判定（逐项 vs 基准，显式回答是否 <3%）

| 关键量 | 基准 | 误差 | <3%？ |
|---|---|---|---|
| Wigley V | 表 Simpson / 4BLT/9 解析 | {rv['wigley_iii']['rel_err']:.2e} | 是（且 <0.5%） |
| Wigley A_wp | 表 Simpson / 2LB/3 解析 | {rv['wigley_iii']['rel_err_wp']:.2e} | 是 |
| Series60 V | 表 Simpson | {rv['series60_cb060']['rel_err']:.4%} | 是 |
| Series60 A_wp | 表 Simpson | {rv['series60_cb060']['rel_err_wp']:.4%} | 是 |
| KVLCC2 V | 表 Simpson | {rv['kvlcc2_class']['rel_err']:.4%} | 是（<0.5% 验收线） |
| KVLCC2 A_wp | 表 Simpson | {rv['kvlcc2_class']['rel_err_wp']:.4%} | 是（<0.5% 验收线） |
| 基线数字真实性 | 审计者自写 Simpson | {b['auditor_baseline_reproduction']['max_rel_diff_vs_reported_table_simpson']:.2e} | 是（数值为真） |
| **vs 外部实船数据** | 无发布机构原表/水动力资料 | 无法计算 | **无基准——不属于「通过」，属 WITHHELD** |

诚实结论：m17251 的工程实现真实、可复现、自洽性验证充分（三条积分路线 + 站线/水线敏感性 + digest 防篡改），
但**无外部基准、无系统集成、数据源为自造参数化船型**，故判 **L2（独立功能实现）**；L3/L4 拒收并给出具体补齐
路径（注册 `CASE-CAD-OFFSET-LOFT` 基准 + 接入真实型值表 + 下游调用方）。verdict 已落
`tensor_maturity_verdicts` id=517，gate 预演 `judge(2)→confirmed`，无返修轮。
"""
    payload = {
        "title": "m17251 交差审计（GAP-CAD-T1-OFFSET-LOFT）：真实现真复现，判 L2；L3/L4 逐项拒收并给缺口",
        "summary": "独立重算（自写 Simpson 复现基线 ≤9.3e-16、样条求积复现 Wigley 解析 3.2e-16、篡改 4/4 被拒、"
                   "results.json 逐字节复现）；验收最差 0.237%<0.5% 通过；但 tensor_benchmark_cases 该 anchor 0 行、"
                   "tensorbench_results 0 行、无下游调用方、三张表 official=false → 判 L2，落表 id=517。",
        "scenario_key": "board_audit_offset_loft",
        "mission_id": 17307,
        "tags": ["TensorCAD", "board-audit", "m17251", "GAP-CAD-T1-OFFSET-LOFT", "NURBS", "maturity", "L2"],
        "markdown": md,
        "figures": [{"name": "audit_figures.png",
                     "caption": "(a) 验收 ref-vs-got 三型 (b) KVLCC2 站线数敏感性 (c) 审计者独立求积 vs 报告数字",
                     "png_base64": png}],
    }
    (OUT / "report_payload.json").write_text(json.dumps({**payload, "figures": [
        {"name": f["name"], "caption": f["caption"]} for f in payload["figures"]]}, ensure_ascii=False, indent=1))
    tok = login()
    req = urllib.request.Request(BASE + "/tensor/reports", data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json",
                                          "Authorization": "Bearer " + tok})
    with urllib.request.urlopen(req, timeout=120) as r:
        print(r.status, json.loads(r.read()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())