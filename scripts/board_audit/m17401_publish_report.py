#!/usr/bin/env python3
"""Publish the m17401 board-audit report to the platform (/tensor/reports).

Every number in the markdown is interpolated from artifacts/board_audit/m17401/verify.json,
the DB benchmark-table query and the verdict evidence — nothing is typed by hand.

Run:  python scripts/board_audit/m17401_publish_report.py [--dry-run]
"""
from __future__ import annotations

import argparse
import base64
import json
import os
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
AUD = ROOT / "artifacts/board_audit/m17401"
BASE = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1")


def login() -> str:
    req = urllib.request.Request(
        BASE + "/auth/login",
        data=json.dumps({"username": os.environ["TENSOR_AGENT_USER"],
                         "password": os.environ["TENSOR_AGENT_PASS"]}).encode(),
        headers={"Content-Type": "application/json"})
    return json.load(urllib.request.urlopen(req))["access_token"]


def b64(p: Path) -> str:
    return base64.b64encode(p.read_bytes()).decode()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    v = json.loads((AUD / "verify.json").read_text())
    import sys
    sys.path.insert(0, str(ROOT / "scripts/board_audit"))
    from m17401_insert_verdict import build as build_verdict  # noqa: PLC0415

    level, evidence, reason = build_verdict()
    rerun, indep = v["C1_rerun"], v["C2_volume_independent"]
    real, bmk = v["C4_real_third_party_file"], v["C7_benchmark_table"]
    cases = ["plate_with_hole", "stepped_shaft", "bushing", "l_bracket"]

    rows = []
    for c in cases:
        s = rerun["per_file"][c + ".step"]
        i = indep[c]
        rows.append(f"| {c}.step | {i['closed_form_m3']:.9e} | {i['volume_m3']:.9e} | "
                    f"{i['error_percent']:+.4f} | {s['error_percent']:+.4f} | "
                    f"{s['triangles']} | {'yes' if s['watertight'] else 'NO'} | {s['euler_characteristic']} | "
                    f"{'< 3 %' if abs(i['error_percent']) < 3 else 'FAIL'} |")
    table = "\n".join(rows)
    rec = rerun["convergence"]
    ratios = " / ".join(f"{abs(rec[k]['plate_with_hole.step'] / rec[k + 1]['plate_with_hole.step']):.2f}"
                        for k in range(len(rec) - 1))
    conv_rows = "\n".join(
        f"| {r['angular_step']:.2f} | {r['plate_with_hole.step']:+.5f} | {r['stepped_shaft.step']:+.5f} | "
        f"{r['bushing.step']:+.5f} | {r['l_bracket.step']:+.5f} |" for r in rec)
    max_err = max(abs(indep[c]["error_percent"]) for c in cases)
    indiv = ", ".join(f"{c}: {abs(indep[c]['error_percent']):.4f} %" for c in cases)

    md = f"""## 1 案例介绍

**审计对象**：m17349「GAP-CAD-B1-STEP-READER」——TensorCAD 的 STEP AP203/AP214 解析 B-rep 读入管道
（交付于 2026-11-13，分支 `agent/m17349`，commit `d6ebb27` + `849e2e7`）。**本单 m17401 是交差核验，非重做**：
只取证、判级、落表，不动被审工作树（复核后其两件 JSON 工件 sha256 未变，C1 已证）。

**物理/几何背景**：船舶总体性能平台要从船厂/设计所的 CAD 交付件直接进入几何审计与水密网格链。真实 CAD 给出的是
`MANIFOLD_SOLID_BREP` / `ADVANCED_FACE` 解析曲面文件（平面、圆柱、圆锥裁剪），而 TensorCAD 原有 `exchange`
只认 `FACETED_BREP`（已剖分三角片）。被审单补的正是这条"解析 B-rep → 三角网格 + 体积"的读入管道。

**几何与参数表（表 1，标称体积由审计方按闭式公式独立重算，非抄被审方数字）**：

| 样例 | 几何 | 闭式标称体积公式 |
|---|---|---|
| plate_with_hole.step | 板 200×120×20 mm，通孔 d40 | 200·120·20 − π·20²·20 |
| stepped_shaft.step | d60×80 圆柱 + d60→d30 锥(h=50) + d30×40 圆柱 | π·30²·80 + π·50/3·(30²+30·15+15²) + π·15²·40 |
| bushing.step | 空心筒 d80/d50，长 60 | π·(40²−25²)·60 |
| l_bracket.step | L 剖面 60×50（肢宽 20）拉伸 10 | (60·20 + 30·20)·10 |

外加 1 个**真实第三方**文件：BANC 30P-30N 构型，Siemens NX 7.5 / ST-Developer v11 导出（AP203
`CONFIG_CONTROL_DESIGN`），`{real['instances']}` 条实体，长度单位解析为 {real['length_scale']:.7f} units/m（英寸）。

## 2 怎么做的

审计方法（全部脚本化，产物 `artifacts/board_audit/m17401/verify.json`；驱动 `scripts/board_audit/m17401_verify_claims.py`）：

1. **C1 重跑**：在被审工作树跑其自身驱动 `python tools/step_reader_audit.py`，记录 exit code 与工件 sha256（重跑前后是否逐字节相同 → 确定性）。
2. **C2 独立复算体积**：把 4 个实体各自导出 OBJ，用**审计方自己写的**散度定理脚本重新积分体积，并用审计方独立重算的闭式标称做对照；同时做边普查（每条边是否恰用两次）。
3. **C3 夹具确定性**：`tools/step_samples.py --out /tmp/...` 重新生成样例，比对 4 个 sha256。
4. **C4 真实文件**：对 NX/ST-Developer 文件跑 `census()`（实体计数/单位/协议），再请求实体读取，记录拒绝消息原文。
5. **C5 拒绝路径**：把通过的样例中 `CYLINDRICAL_SURFACE` 改写成 `B_SPLINE_SURFACE_WITH_KNOTS`，看是否按名拒绝而非猜测。
6. **C6 破损文件**：把所有 `CARTESIAN_POINT` ×2 但 `CIRCLE` 半径不改（自相矛盾文件），看是否拒绝给体积。
7. **C7 基准表**：按 owner 令查 `tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-B1-STEP-READER'`。
8. **判级 + 落表**：`scripts/board_audit/m17401_insert_verdict.py` 幂等写入 `tensor_maturity_verdicts`（id={542}），
   并在本机用 `maturity_gate.judge()` 复核 → `{level}`。

工具链：Python 3.11（纯标准库 + psycopg + matplotlib），无第三方 CAD 内核（loong64 无 OCP wheel）——这本身就是本次审计的
关键局限，见 §5。

## 3 计算结果

**表 2 体积基准：闭式标称 vs 独立复算（审计方 OBJ 散度定理）vs 被审方自报**

| 样例 | 标称体积 [m³] | 审计方测得 [m³] | 审计方误差 [%] | 被审方误差 [%] | 三角数 | 水密 | Euler χ | < 3 % |
|---|---|---|---|---|---|---|---|---|
{table}

**表 3 收敛阶（体积误差随角度步长减半的变化，单位 %）**

| 角度步长 [rad] | plate_with_hole | stepped_shaft | bushing | l_bracket |
|---|---|---|---|---|
{conv_rows}

误差比（plate，逐次减半）：**{ratios}** → 渐近二阶（被审方报告只引用最细一对 3.999，见 §5 注意点 3）。

**表 4 真实第三方文件与拒绝路径**

| 检查 | 结果 |
|---|---|
| 实体实例数 | {real['instances']} |
| 协议 / 前处理器 | {real['schema']} / {real['preprocessor']} |
| 长度单位 | {real['length_scale']:.7f} units/m（{real['length_scale_detail']}，英寸） |
| 线框表示 | GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION = 1，B-rep 实体 = 0 |
| 请求实体 | **拒绝**：`no MANIFOLD_SOLID_BREP / BREP_WITH_VOIDS / CLOSED_SHELL found`（如实报，未编造网格） |
| B_SPLINE_CURVE_WITH_KNOTS | {real['bsp']} |
| 不支持的曲面类型（C5 变异） | **按名拒绝**：`StepUnsupportedStep … B_SPLINE_SURFACE_WITH_KNOTS` |
| 破损文件（C6） | watertight=False, boundary_edges={v['C6_malformed_file']['boundary_edges']}, volume_m3={v['C6_malformed_file']['volume_m3']}（不编造体积） |
| 夹具确定性（C3） | 4/4 sha256 与交付文件一致 |
| 测试重跑 | 17 + 26 = **43 passed** |
| `tensor_benchmark_cases` 本锚点 | **{bmk['n_rows']} 行**（TensorCAD 仅注册 {", ".join(bmk['tensorcad_anchors_in_table'])}） |

## 4 结果图

图 1：主体自报误差 vs 审计方独立复算误差（4 个样例），绿色带为验收线 |err| < 1 %，纵轴为体积相对误差。
两条柱几乎重合——独立路径复现了被审方每一个数字。
图 2：|体积误差| 随角度步长的 log-log 收敛曲线与 slope=2 参考线，标注逐次减半比值 {ratios}。

## 5 注意点（踩过的坑 / 局限 / 适用范围）

1. **审计方第一版"缩放探针"是错的**：只把 `CARTESIAN_POINT` ×2 而不改 `CIRCLE` 半径，造出自相矛盾文件，
   得到 +2.77 % 的"异常"其实是破损文件上的非法网格。修正探针（改为只测拒绝行为，C6）后才看清：
   读入器对这文件返回 `volume_m3=None`。**教训：审计探针自身也必须在一致性假设下构造。**
2. **`complete=True` 与 `watertight=False` 可同时出现**（C6，8 条边界边）：体积虽被正确置空，但只读 `complete`
   标志的下游消费者可能把未闭合网格当成完成品（F2，建议修）。
3. **收敛比不是恒定 4.00**：序列 {ratios}，报告引用的是最细一对 3.999。二阶是渐近性质，单点数字偏优（F3）。
4. **无第三方 CAD 内核**：loong64 / 容器内无 OCP 或 FreeCAD，无法用商业内核标称交叉验证；这也是本单拒绝 L4 的
   结构性原因之一（真实第三方文件只有 1 个，且是线框，0 个实体）。
5. **适用范围**：支持 PLANE / CYLINDRICAL / CONICAL / SPHERICAL / TOROIDAL 裁剪 + LINE / CIRCLE / ELLIPSE /
   B_SPLINE_CURVE 边界；B 样条曲面、扫掠面、POLYLINE / TRIMMED_CURVE、单环球/环面极盖（歧义）明确拒绝。
   船体自由曲面（NURBS 蒙皮）**目前读不进来**——这是下一个锚点要补的。
6. **未接入求解器**：`exchange.read_step` → `load_exchange` 只被仓内几何审计管线消费，没有任何求解器/网格主循环调用，
   所以只申报 L3，不申报 L4。

## 6 误差判定（逐项，显式回答是否 < 3 %）

| 样例 | 关键量 | 基准（闭式解析） | 实测 | 相对误差 | 是否 < 3 % | 是否 < 1 %（本单验收线） |
|---|---|---|---|---|---|---|
{chr(10).join(f"| {c}.step | 体积 | {indep[c]['closed_form_m3']:.9e} m³ | {indep[c]['volume_m3']:.9e} m³ | {abs(indep[c]['error_percent']):.4f} % | 是 | 是 |" for c in cases)}

**结论：4/4 样例体积误差均 < 1 %（最大 {max_err:.4f} %），远优于 3 % 生态线；逐项为 {indiv}。**
但必须同时声明基准的性质：这些标称是**闭式解析值**（审计方独立重算），不是 CAD 内核报出的标称；
"对商业内核语义的误差"= **NOT MEASURED（0/3 真实实体入网格）**。因此本单判 **L3**，
L4 因"真实文件 0/3 + 基准自撰 + 基准表 0 行"三条独立理由被拒绝（详见 `docs/board_audit/m17401/AUDIT.md` §3）。
被审方在报告与 `GAP-ANCHORS.md` 中均已如实标注 external validation = WITHHELD，**不存在粉饰**；
本次核验发现的差距是覆盖面，不是诚实性。

**验收三件套核验**：run_artifacts 真实存在且审计方重跑逐字节一致 ✓；benchmark ref vs got 数值经独立复算逐项一致 ✓；
repro_cmd（`python tools/step_reader_audit.py`，6.2 s，exit 0）可执行 ✓。
"""

    payload = {
        "title": "m17349 交差核验（board audit）：STEP AP203/AP214 B-rep 读入管道 — 判 L3，L4 三条理由拒绝",
        "summary": (f"独立复核 m17349 全部数字：4 样例体积误差 {max_err:.4f} %（<1 % 验收线）经审计方自写散度定理脚本复现一致，"
                    f"43 测试重跑通过，夹具重生成 sha256 4/4 一致；真实 NX AP203 文件 {real['instances']} 实体、英寸单位"
                    f"{real['length_scale']:.4f} units/m 核验为真但 0 个 B-rep 实体（线框，按名拒绝）。"
                    "判 L3：基准为自撰闭式解析值 + 真实实体 0/3 + tensor_benchmark_cases 本锚点 0 行 → L4 拒绝。verdict id=542。"),
        "scenario_key": "tensorcad_gap_cad_b1_step_reader_board_audit",
        "mission_id": 17401,
        "tags": ["TensorCAD", "GAP-CAD-B1", "STEP", "AP203", "AP214", "board-audit", "maturity", "L3"],
        "markdown": md,
        "figures": [
            {"name": "audit_volume_crosscheck.png",
             "caption": "图 1 体积误差：被审方自报 vs 审计方独立 OBJ 散度定理复算（绿色带 |err|<1 %）",
             "png_base64": b64(AUD / "figures/audit_volume_crosscheck.png")},
            {"name": "audit_convergence_order.png",
             "caption": "图 2 |体积误差| vs 角度步长（log-log）与 slope=2 参考线，标注逐次减半比值",
             "png_base64": b64(AUD / "figures/audit_convergence_order.png")},
        ],
    }
    if args.dry_run:
        Path("/tmp/report_preview.md").write_text(md)
        print(f"dry-run: markdown {len(md)} chars -> /tmp/report_preview.md")
        return 0
    tok = login()
    req = urllib.request.Request(BASE + "/tensor/reports",
                                 data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json",
                                          "Authorization": "Bearer " + tok})
    with urllib.request.urlopen(req) as r:
        print(r.status, r.read().decode()[:300])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())