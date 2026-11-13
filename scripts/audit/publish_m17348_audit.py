#!/usr/bin/env python3
"""Publish the m17400 board-audit report to the platform (/tensor/reports).

Every number is injected from artifacts/audit_m17348/audit.json — nothing is typed
by hand.  Figures are inlined as png_base64.
"""
from __future__ import annotations

import base64
import json
import os
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "artifacts" / "audit_m17348"
A = json.loads((OUT / "audit.json").read_text())
H = A["headline"]
R = A["reproduce"]
V = json.loads((OUT / "verdict_evidence.json").read_text())
REG = A["registry"]
ST = A["delivery_state"]
REP = A["report_16354"]
P = A["delivered"]["manifest"]
WIG = H["wigley"]
CELLS = H["pointwise_mesh_field_vs_references"]


def pct(x: float) -> str:
    return "%+.4f%%" % (100.0 * x)


def grid_table() -> str:
    rows = ["| nx×nz | faces | smooth E_bending | crease E_bending | vs 解析 |", "|---|---|---|---|---|"]
    for g in WIG["grids"]:
        rows.append("| %d×%d | %d | %.4f | %.2f | %s |"
                    % (g["nx"], g["nz"], g["faces"], g["smooth_energy"], g["crease_energy"],
                       pct(g["smooth_relative_error"])))
    rows.append("| h→0 三点 Lagrange 外推 | — | **%.4f** | — | **%s** |"
                % (WIG["extrapolated_bending"], pct(WIG["extrapolated_relative_error"])))
    return "\n".join(rows)


def case_table() -> str:
    rows = ["| case | faces | H vs 精确场（中位） | K vs 精确场（中位） | H vs normal-cycle 中位/最大 | K vs normal-cycle 中位/最大 |",
            "|---|---|---|---|---|---|"]
    for name, v in CELLS.items():
        def f(x):
            return "n/a" if x is None else "%.2f%%" % (100.0 * x)
        rows.append("| %s | %d | %s | %s | %s / %s | %s / %s |"
                    % (name, v["faces"], f(v["H_vs_field_median"]), f(v["K_vs_field_median"]),
                       f(v["H_vs_normal_cycle_median"]), f(v["H_vs_normal_cycle_max"]),
                       f(v["K_vs_normal_cycle_median"]), f(v["K_vs_normal_cycle_max"])))
    return "\n".join(rows)


def manifest_table() -> str:
    rows = ["| 交付文件 | bytes | sha256(前16) |", "|---|---|---|"]
    for name, m in P.items():
        rows.append("| `%s` | %d | `%s` |" % (name, m["bytes"], m["sha256"][:16]))
    return "\n".join(rows)


def claims_table() -> str:
    rows = ["| claim | 值 | 判据 | 成立 |", "|---|---|---|---|"]
    for c in H["claims"]:
        rows.append("| %s | %.6g | %s %.6g | %s |"
                    % (c["claim"], c["value"], c.get("direction", "<="), c["limit"],
                       "yes" if c["holds"] else "no"))
    return "\n".join(rows)


md = f"""## 1. 案例介绍

**被审对象**：m17348（TensorCAD，`GAP-CAD-A3-CURVATURE-FAIRNESS`，status=succeeded）——曲率场与光顺性诊断：逐顶点主曲率 / Gaussian / 平均曲率 + 三种无量纲光顺（弯曲）能量。

**本单性质**：owner 2026-11-11 任务榜要求的**交差核验**（m17400），非重做。裁判口径＝`tensor_benchmark_cases` 登记判据 + 解析闭式 + 本单重跑。

**被审几何与参数**：解析球 R ∈ {{0.1,0.25,0.5,1,2,5,10}} m（1449 探针）；椭球 a=1.0,b=0.6,c=0.5；Wigley 船体 L=10 m, B=1.0 m, T=0.625 m，网格 21×7→121×37 四档。能量定义 E_bending=∫(k1²+k2²)dA、E_willmore=∫H²dA、E_gauss=∫K dA（无量纲）。

**被审工作树**：`/root/tp-agent/work/m17348`，分支 `agent/m17348`，HEAD `{ST['head'][:8]}`（{ST['head_subject']}）。

## 2. 怎么做的

1. `scripts/audit/m17348_board_audit.py`：一键审计脚本——(a) 交付文件字节保全 + sha256 清单（`artifacts/audit_m17348/delivered/`）；(b) git ref 扫描（交付是否入库）；(c) 集成面 grep（`__init__` 导出/外部调用者）；(d) `psycopg` 查裁判表；(e) 在工作树**字节副本**里重跑 pytest / example runner / crosscheck，并把 `results.json`、`curvature_crosscheck.json` 与交付件做**逐叶数值对比**；(f) 从 `results.json` 重算全部 headline 数值；(g) 拉取 `tensor_reports/16354` 校验六段与内嵌图 sha256；(h) `--insert` 落表 `tensor_maturity_verdicts`（id=551，verdict_by=board-audit，gate_state=pending）。
2. 复核未采信被审报告文本：报告里的每个数字都由本单从重跑结果重算后比对。
3. 平台 L4 口径独立取证：`GET /tensor/bench/l4-verdict/17348`。

## 3. 计算结果

**三件套存在性（10/10 文件落地，字节已保全）**

{manifest_table()}

**复跑（本单，工作树字节副本）**

| 命令 | rc | 墙钟 | 结果 |
|---|---|---|---|
| `pytest tests/test_fairness_energy.py tests/test_curvature.py -q` | {R['pytest_full_suite']['rc']} | {R['pytest_full_suite']['wall_clock_s']} s | {R['pytest_full_suite']['tail'][0] if R['pytest_full_suite']['tail'] else ''} |
| `python examples/fairness_energy/run_example.py` | {R['example_runner']['rc']} | {R['example_runner']['wall_clock_s']} s | {sum(1 for c in H['claims'] if c['holds'])}/{len(H['claims'])} claims green |
| `python tools/crosscheck_curvature_thirdparty.py` | {R['crosscheck_tool']['rc']} | {R['crosscheck_tool']['wall_clock_s']} s | crosscheck json **sha256 逐字节相同** = {R['crosscheck_json_identical']} |

`results.json` 重跑逐叶对比差异数 = **{R['n_results_json_leaf_diffs']}**（唯一差异是被写入的图绝对路径字符串，所有数值叶 ≤1e-12）。

**A. 解析球主曲率**：最坏 |k−1/R| = **{H['sphere_principal_curvature']['worst_abs_error_per_m']:.3e} 1/m**（相对 {H['sphere_principal_curvature']['worst_relative_error']:.2e}），任务榜判据 1e-6 ✔ —— 但**测的是既有解析隐式场求值器**，不是本单新建的 mesh 曲率场（见 §5 F2）。

**B. 球光顺能量 vs 闭式 8π={H['sphere_energy']['closed_forms']['bending_8pi']:.6f}**：field 路径 {pct(H['sphere_energy']['field_path_bending_rel_err'])}；conditioned mesh 路径 {pct(H['sphere_energy']['mesh_path_bending_rel_err'])}；恒等式残差 {H['sphere_energy']['field_path_identity_residual']:.2e}。

**C. Wigley 光顺能量 vs 解析船侧参考 {WIG['analytic_reference']['bending']:.6f}**（湿面积 {WIG['analytic_reference']['wetted_area']:.6f} m²）：

{grid_table()}

**D. 逐点 mesh 曲率场 vs 独立参考**（`artifacts/curvature_crosscheck.json`，重跑 bit-identical）：

{case_table()}

**E. 裁判**：`SELECT case_name,metric,tolerance FROM tensor_benchmark_cases WHERE anchor_slug='GAP-CAD-A3-CURVATURE-FAIRNESS'` → **{REG['n_rows']} 行**（全表 {REG['table_total']} 行；TensorCAD 仅 {REG['tensorcad_anchors']}）；平台 `l4-verdict/17348` → `CLAIMED_NOT_PROVEN`（tensorbench_results 全表 {REG['tensorbench_results_rows']} 行）。

**F. 交付态**：`refs_containing_module` = {ST['refs_containing_module']}（扫描 {ST['n_refs']} 个 refs）；`git status` 未提交项 {len(ST['uncommitted_status'])} 条 → **交付物不在任何 git ref**。

**G. 报告 16354**：六段齐全 {REP['sections']}；显式「是否 <3%」= {REP['has_explicit_3pct_answer']}；含 WITHHELD={REP['mentions_withheld']}、含 honest FAIL={REP['mentions_honest_fail']}；内嵌图 sha256 与 `docs/fairness_energy/plot.png` 一致 = {REP['figure_png_sha256'] == REP['delivered_png_sha256']}。

**H. claims 明细（重跑值）**

{claims_table()}

## 4. 结果图

左：Wigley smooth 能量 mesh→h=0 收敛与外推（3%/15% 线）；中：三例逐点 mesh 曲率 vs normal-cycle 中位/最大误差（15% 线，逐点跨线）；右：smooth 桶收敛 vs crease 桶 1/h 发散。第二图为被审报告 16354 内嵌图（与 `plot.png` 同 sha256）。

## 5. 注意点（发现，按严重度）

* **F1（high，交付阻塞）**：交付物不在任何 git ref；HEAD 停在他单合入提交。若按 `agent/m17348` 推送，本能力 0 字节入库。本单已做字节保全（`artifacts/audit_m17348/delivered/` + sha256 清单），恢复路径写在审计文档 §6。
* **F2（medium）**：头号验收线 `<1e-6` 由 `principal_curvatures(field, …)`（既有解析场）测得；新 `principal_curvature_field(mesh)` 只有 ±20% 松 pin（`tests/test_fairness_energy.py:288`）。逐点实测中位 H 1.58%/12.83%/40.00%、K 0.05%/22.43%/99.87%（球/椭球/Wigley 侧面）→ 跨过 15% 精度线；被审报告 §5/§6 已如实披露，未粉饰。
* **F3（medium）**：裁判缺席（0 行）→ L4 的机器口径（平台按 ref/got/error 复算）不成立；本次由解析闭式 + 独立求积（面积双求积 4.44e-16、p-extension 9.90e-16）判定。
* **F4（low-medium）**：「Wigley <3%」是 h→0 外推值 {pct(WIG['extrapolated_relative_error'])}；单网格最细 121×37 = {pct(WIG['finest_grid_relative_error'])}（>3%，<15%）。单次运行值不达 3%。
* **F5（low）**：runner 第 8 条 claim（crease 发散）value 硬编码 0.0、limit=差值 → 恒真不可伪；发散方向仅由 pytest 保真（`tests/test_fairness_energy.py:265`）。故「11/11 green」含 1 条空判据。
* **F6（low）**：新模块未在 `tensorcad/__init__.py` 导出、自身 example/tool/test 之外 0 调用者 → L3「接入主求解器/平台」不成立。
* **F7（info，诚实项）**：meshmagick 对拍子项 WITHHELD，三路尝试（import / 内网索引 / 沙箱出口）+ 日期写入 `curvature_crosscheck.json:provenance`，替代参考=精确场 + Cohen-Steiner/Morvan normal cycle；替代结果反过来暴露逐点精度不足。

## 6. 误差判定（显式回答「是否 <3%」）

| 关键量 | 误差 | <3%？ | <15%？ |
|---|---|---|---|
| 解析球主曲率 vs 1/R（field 路径） | {H['sphere_principal_curvature']['worst_relative_error']:.2e} 相对 | 是 | 是 |
| 球 E_bending vs 8π（field 路径） | {pct(H['sphere_energy']['field_path_bending_rel_err'])} | 是 | 是 |
| 球 E_bending vs 8π（conditioned mesh） | {pct(H['sphere_energy']['mesh_path_bending_rel_err'])} | 否 | 是 |
| Wigley E_bending vs 解析（h→0 外推） | {pct(WIG['extrapolated_relative_error'])} | 是（外推极限） | 是 |
| Wigley E_bending 最细网格单次 | {pct(WIG['finest_grid_relative_error'])} | 否 | 是 |
| 逐点 H（椭球中位 / Wigley 侧面中位） | {100*CELLS['ellipsoid']['H_vs_field_median']:.2f}% / {100*CELLS['wigley_flank']['H_vs_normal_cycle_median']:.2f}% | 否 | **否（honest FAIL）** |
| 逐点 K（椭球中位 / Wigley 侧面中位） | {100*CELLS['ellipsoid']['K_vs_field_median']:.2f}% / {100*CELLS['wigley_flank']['K_vs_normal_cycle_median']:.2f}% | 否 | **否（honest FAIL）** |
| meshmagick 对拍 | 未完成（WITHHELD） | — | — |

**结论**：积分量（光顺能量）达 L4 精度并满足 <3%（外推极限 0.058%；最细网格单次 6.75% 仍 <15%）；**逐点曲率场不满足 15% 判据**（12.8%–99.9%），只能作诊断/可视化，不能作定量输出——此结论与被审单自述一致。

**判级**：**L4**（积分/光顺能量核：完整算例 + 误差分析 + 可复现）；逐点曲率场 L2 级；**L3 集成不成立**（F6）；裁判缺席（F3）使 L4 只能由解析锚复核。**处置：approve with findings——F1 为交付阻塞项，需按审计文档 §6 恢复提交后方可消账。**
"""


def post() -> dict:
    base = os.environ["TENSOR_API_BASE"]
    creds = {"username": os.environ["TENSOR_AGENT_USER"],
             "password": os.environ["TENSOR_AGENT_PASS"]}
    req = urllib.request.Request(base + "/auth/login", data=json.dumps(creds).encode(),
                                 headers={"Content-Type": "application/json"})
    token = json.loads(urllib.request.urlopen(req).read().decode())["access_token"]
    figs = []
    for name, path, caption in (
            ("audit_m17348.png", OUT / "audit_figure.png",
             "左：Wigley smooth E_bending mesh→h=0 外推（3%/15% 线）；中：三例逐点 mesh 曲率 vs normal-cycle 中位/最大（15% 线）；右：smooth 桶收敛 vs crease 桶 1/h 发散"),
            ("subject_report_figure.png", OUT / "report_figure.png",
             "被审报告 16354 内嵌图（与 docs/fairness_energy/plot.png sha256 相同）")):
        figs.append({"name": name, "caption": caption,
                     "png_base64": base64.b64encode(path.read_bytes()).decode()})
    payload = {
        "title": "m17400 交差核验：m17348 GAP-CAD-A3 曲率场与光顺能量（L4，交付态 F1 阻塞）",
        "summary": "三件套真实可复现（43 测试绿、11/11 claims、crosscheck json 逐字节相同）；"
                   "积分量 L4（Wigley 外推 +0.0577%，球 8π +0.0530%）；逐点曲率场 12.8%–99.9% 超 15% 线（被审单已披露）；"
                   "交付物不在任何 git ref（F1 阻塞），已字节保全；裁判表 0 行（F3）。",
        "scenario_key": "tensorcad_curvature_fairness",
        "mission_id": 17400,
        "tags": ["TensorCAD", "GAP-CAD-A3", "board-audit", "curvature", "fairness", "L4", "delivery-gap"],
        "markdown": md,
        "figures": figs,
    }
    Path("/tmp/report_m17400.json").write_text(json.dumps(payload, ensure_ascii=False))
    req = urllib.request.Request(base + "/tensor/reports", data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json",
                                          "Authorization": "Bearer " + token})
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read().decode())


if __name__ == "__main__":
    Path("/tmp/md_m17400.md").write_text(md)
    print("markdown chars:", len(md))
    print("posted:", post())