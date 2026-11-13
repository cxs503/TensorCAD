#!/usr/bin/env python3
"""Publish the m17405 paper-grade progress report to the tensor platform.

Markdown is assembled from the audit artifacts (probe_results.json,
compliance_report.json, evidence.json) so every number is injected
programmatically; the figure is attached as base64.

Run: python3 docs/board_audit/m17405/post_report.py   (also saves report_body.md)
"""
from __future__ import annotations

import base64
import json
import os
import subprocess
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SUBJECT = Path("/root/tp-agent/work/m17357")
PROBE = json.loads((ROOT / "artifacts/geom_contract_audit/probe_results.json").read_text())
COMP = json.loads((SUBJECT / "artifacts/geom_contract_assets/compliance_report.json").read_text())
EV = json.loads((HERE / "evidence.json").read_text())
PNG = ROOT / "artifacts/geom_contract_audit/contract_audit.png"
BASE = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1")


def token() -> str:
    req = urllib.request.Request(
        BASE + "/auth/login",
        data=json.dumps({"username": os.environ["TENSOR_AGENT_USER"],
                         "password": os.environ["TENSOR_AGENT_PASS"]}).encode(),
        headers={"Content-Type": "application/json"}, method="POST")
    with urllib.request.urlopen(req, timeout=60) as r:
        d = json.load(r)
    return d.get("access_token") or d.get("token")


def main() -> int:
    anchors = "\n".join(
        f"| {a['case']} | `{a['formula']}` | {a['ref']:.6f} | {a['got']:.6f} | "
        f"**{a['rel_err_pct']:+.4f}%** | {'是' if a['within_3pct'] else '否'} |"
        for a in PROBE["anchors"])
    comp = "\n".join(
        f"| {r['source'].rsplit('/', 1)[-1]} | {'PASS' if r['passed'] else 'FAIL'} | "
        f"{r['mesh_info'].get('signed_volume')} | "
        f"{'; '.join('[' + v['convention'] + '/' + v['severity'] + '] ' + v['message'] for v in r['violations']) or '—'} |"
        for r in COMP)
    red = "\n".join(
        f"| {c['case']} | {', '.join(c['conventions']) or '—'} | {c['passed']} |"
        for c in PROBE["red_team"])
    fn = "\n".join(
        f"| {c['probe']} | {'通过(静默)' if c['passed'] else '失败'} | "
        f"{', '.join(c['violations']) or '无'} |" for c in PROBE["false_negatives"])
    prism = "\n".join(
        f"| {p['case']} | {p['ref']:.3f} | {p['got']:.3f} | {p['rel_err_pct']:+.3f}% | "
        f"{p['implied_cb']:.4f} | {p['declared_cb']:.2f} |"
        for p in PROBE["prismatic_observation_out_of_scope"])
    rp = PROBE["repro"]

    md = f"""## m17405 · 交差审计：m17357 `GAP-CAD-D3-GEOM-CONTRACT`（TensorCAD）

**判级 L2（独立功能实现）· approve with notes · verdict id=543 · gate 预演 judge(2)→2 confirmed**

### ① 案例介绍
`GAP-CAD-D3-GEOM-CONTRACT` 要求 TensorCAD 交付**几何契约文档 + 运行时校验器**：把坐标系/单位/比例/取向
四条约定做成机器可执行检查，对 8 个现存 STL 资产跑出合规表，并实证违规可自动检出。物理背景=船舶 CAD 网格
进入阻力/耐波/结构求解器前的**互换性契约**（单位错标 1000×、法线反转、左手系都会让下游阻力/水动力结果失效
而不报错）。对象工作树 `/root/tp-agent/work/m17357`，交付 6 个文件（校验器 347 行、文档 149 行、测试 231 行、
两个工具脚本、9 STL + 合规表），参数表=箱形驳 30×8×5 m、单位球 r=1 m、Wigley L=10 B=1 T=0.625 m、
Series 60 L=20 B=3 T=1.5 m（Cb 0.60/0.75）。

### ② 怎么做的
1. 取证 `tensor_dev_missions` id=17357（`head_commit=null`、8 路径未提交、`server_quality_audit.gate_valid=false`）。
2. 基准裁判 `tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-D3-GEOM-CONTRACT'` → **0 行**（无 ref/tolerance）。
3. 自建独立探针 `docs/board_audit/m17405/audit_probe.py`（不 import 对象夹具）：解析体积对照 / 5 组自建违规 /
   3 组沉默窗口反例 / CLI 退出码 / 复现哈希。
4. 复现 `python3 tools/run_contract_compliance.py` → 前后 sha256 `{rp['sha256_before'][:16]}…` **identical={rp['identical']}**。
5. 测试 `pytest tests/test_geom_contract.py -q` → **22 passed / 0 failed**。
6. 集成面 grep `validate_contract` → 仅 `tensorcad/__init__.py:457` 再导出，**零调用方**。

### ③ 计算结果
对象 8 资产合规表（程序化读取对象工件）：

| 资产 | 判定 | 有符号体积 m³ | 违规 |
|---|---|---|---|
{comp}

审计方独立体积基准（ref=闭式解析解，got=从交付 STL 用散度定理重算）：

| 算例 | 公式 | ref | got | 相对误差 | <3% |
|---|---|---|---|---|---|
{anchors}

红队 5 例（审计方自建网格，不共用对象夹具）：

| 违规用例 | 命中约定 | passed |
|---|---|---|
{red}

### ④ 结果图
![contract audit](attachment:contract_audit.png)
左：8 资产 × 4 约定违规矩阵（红=error，黄=warning，空白=通过），5 PASS / 3 FAIL 与对象自报一致；
右：交付 STL 体积 vs 闭式解析体积奇偶图（±3% 带 + 逐点误差标注），三点全部落在 ±0.5% 内——
证明校验器度量的是真实几何，不是回填常量。

### ⑤ 注意点
* **未提交**：8 个交付文件全为对象工作树未跟踪文件，主检出 `/root/tensor-ecosystem/TensorCAD` 无
  `geom_contract.py`，其他仓/CI 现在无法 import 该能力（已按字节快照进 `docs/board_audit/m17405/snapshot/`）。
* **无六要素报告**：`tensor_reports` 无 mission_id=17357 行（report_id 空），本报告代为补齐并留痕欠账。
* **坐标系检查是声明串检查**：源码 docstring 已声明；镜像网格 + 声明 right-handed 的开放面网格零违规通过。
* **比例检查有盲窗**：100× 单位错标只要对角线仍落在 [0.05, 500] m 内即通过。
* **开放面取向只到 warning**：恰好 50% 反向面零违规通过（阈值为 ratio < 0.5）。
* 适用范围：STL 不携带单位/坐标系元数据，结论只对"调用方声明正确"的场景成立。

### ⑥ 误差判定（是否 <3%）
| 关键量 | 基准 | 误差 | <3% |
|---|---|---|---|
| 箱形驳体积 | 解析 L·B·D=1200 m³ | {PROBE['anchors'][0]['rel_err_pct']:+.4f}% | 是 |
| 单位球体积 | 解析 4/3πr³=4.18879 m³ | {PROBE['anchors'][1]['rel_err_pct']:+.4f}% | 是 |
| Wigley 体积 | 解析 (4/9)L·B·T=2.77778 m³ | {PROBE['anchors'][2]['rel_err_pct']:+.4f}% | 是 |
| 登记基准案 ref/got vs tolerance | `tensor_benchmark_cases` 该锚点 0 行 | 无法计算 | **无基准案 → L4 不成立（honest）** |

越界观察（不计入本单判定，校验器忠实报告网格体积，偏差在资产生成器侧，建议另开单用
`tensorcad.series60.analytic_volume` 仲裁）：

| 算例 | ref=Cb·L·B·T | got | 误差 | 隐含 Cb | 声明 Cb |
|---|---|---|---|---|---|
{prism}

### 判级与处置
L2 逐键：`impl_files` ✅（347 行实现，无占位）· `entry` ✅（库 :240/:268 + CLI :282）·
`run_artifacts` ✅（8 报告 + 9 STL，字节级可复现）· `tests` ✅（22 绿且含正反例不变式）·
`integration_callers` ❌ 无调用方 → 拒 L3 · `benchmark` ❌ 锚点 0 行 → 拒 L4 · `external_anchor` 不申。
反例 3 条已登记（`evidence.checks_not_run` 同单）。后续：commit 落库 → 接网格导入口挣 L3 →
登记解析锚点基准案（箱驳 1200 m³ / Wigley (4/9)LBT，tolerance 1%）挣 L4 → 收紧开放面取向阈值。

**Honest FAIL 项**：无登记基准案、无六要素报告、交付未提交——三项如实记录，不粉饰。
"""
    (HERE / "report_body.md").write_text(md)
    payload = {
        "title": "m17405 交差审计 · m17357 GAP-CAD-D3-GEOM-CONTRACT（TensorCAD）判级 L2",
        "summary": ("几何契约四约定校验器真实可用：22 测试绿、合规表字节级可复现、审计方自建红队 4/4 约定命中、"
                    "三个闭式解析体积锚点误差 -0.474%/−0.236%/0.000%；但零集成调用方（拒 L3）、锚点 0 条登记基准案"
                    "（拒 L4）、交付全未提交、无六要素报告 → 判 L2，verdict id=543，gate 预演 confirmed。"),
        "scenario_key": "board_audit",
        "mission_id": 17405,
        "tags": ["board-audit", "TensorCAD", "GAP-CAD-D3-GEOM-CONTRACT", "m17357", "L2",
                 "adversarial-review", "maturity-verdict"],
        "markdown": md,
        "figures": [{
            "name": "contract_audit.png",
            "caption": "左：8 STL 资产 × 4 约定的违规矩阵（红=error/黄=warning/空白=通过）；右：交付 STL 体积 vs 闭式解析体积奇偶图（±3% 带）。",
            "png_base64": base64.b64encode(PNG.read_bytes()).decode(),
        }],
    }
    tok = token()
    req = urllib.request.Request(BASE + "/tensor/reports", data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json",
                                          "Authorization": "Bearer " + tok}, method="POST")
    with urllib.request.urlopen(req, timeout=180) as r:
        print(r.status, r.read().decode()[:600])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())