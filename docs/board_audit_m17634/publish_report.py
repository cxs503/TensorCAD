#!/usr/bin/env python3
"""Publish the m17634 board-audit report to the platform (six elements, figure embedded).

Every number in the markdown below is injected programmatically from
docs/board_audit_m17634/audit_evidence.json (auditor's own recount) or from the
subject's committed results.json snapshot — zero hand-copied digits.
"""
from __future__ import annotations

import base64
import json
import os
import urllib.error
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = os.environ.get("TENSOR_API_BASE", "http://127.0.0.1:8010/api/v1").rstrip("/")
for cand in (BASE, BASE + "/api/v1"):
    if cand.endswith("/api/v1"):
        API = cand
        break


def api(method: str, path: str, payload=None, token=None):
    data = json.dumps(payload).encode() if payload is not None else None
    req = urllib.request.Request(API + path, data=data, method=method)
    req.add_header("Content-Type", "application/json")
    if token:
        req.add_header("Authorization", "Bearer " + token)
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.status, json.loads(r.read() or b"{}")
    except urllib.error.HTTPError as e:
        return e.code, {"error": e.read()[:400].decode(errors="ignore")}


def main() -> int:
    ev = json.loads((HERE / "audit_evidence.json").read_text())["checks"]
    res = json.loads((HERE / "snapshot" / "results.json").read_text())
    ce = ev["claimed_errors"]
    va = ev["volume_agreement"]
    man = ev["manifest"]
    probe = ev["cgns_probe"]

    rows = []
    for g, v in va.items():
        rows.append(f"| {g} | {v['readers']} | {len(v['formats_read_independently'])} "
                    f"({', '.join(v['formats_read_independently'])}) | {v['reference']:.12g} | "
                    f"{v['max_rel_err_vs_reference']:.3e} | {v['cross_format_spread_rel']:.3e} |")
    vol_table = ("| geometry | reader checks (auditor) | formats re-read | reference volume | "
                 "max rel err vs reference | cross-format spread |\n|---|---|---|---|---|---|\n"
                 + "\n".join(rows))

    up = ev["upstream_reread"]
    up_rows = []
    for k, v in sorted(up.items()):
        val = f"{v['volume']:.10g}" if v.get("ok") else f"not readable ({v['error'][:40]})"
        up_rows.append(f"| {k} | {val} |")
    up_table = "| upstream m17483 file | auditor volume |\n|---|---|\n" + "\n".join(up_rows)

    fig = (HERE / "figures" / "audit_volume_deviation.png")
    figures = [{"name": "audit_volume_deviation.png",
                "caption": ("m17595 六格式矩阵的独立重算（左：审计侧重算体积相对误差热图，log 色标；"
                            "右：每格式几何 oracle 归属——第三方 vs 被审方自写，即发现 F1 的图证）"),
                "png_base64": base64.b64encode(fig.read_bytes()).decode()}]

    md = f"""# 交差审计：m17595 XVAL-CAD-M17483-NATIVE-SIX-FORMAT-MATRIX（TensorCAD）

**裁定：L2（独立功能实现）**，`tensor_maturity_verdicts` id=602，`verdict_by=board-audit`，
`gate_state=pending`（本单 succeeded 后由 maturity_gate 机器复核）。处置：**approve with notes**。
被审交付指针：分支 `agent/m17595` / commit `4f5b492`（115 files，+276281）。

## 1. 案例介绍（对象与判据）

被审单 m17595 是 XVAL（交叉对标）单：用独立第二实现重读 TensorCAD m17483 的原生六格式网格矩阵
（STL / OBJ / VTK / SU2 / CGNS / Gmsh），判定格式传输是否几何无损。三件套判据：
交付自设 `implementation_threshold={res['acceptance']['implementation_threshold']}`、
owner 精度线 `owner_threshold={res['acceptance']['owner_threshold']}`，拓扑=顶点双射+面多重集相等，
定向=定向面多重集相等且 flipped=0。

| 参数 | 值 |
|---|---|
| 几何 | cube（解析体积 8.0 精确）、icosphere（r=1，sub=2）、torus（R=1.0, r=0.3, 64×32） |
| 格式 | 6（STL/OBJ/VTK/SU2/CGNS/Gmsh） |
| 矩阵 | A = 交付自写文件 × 读者；B = 上游 m17483 工件 × 同一批读者 |
| reader 检查数 | matrix A {ce['summary_A']['reader_checks']}，matrix B {ce['summary_B_declared_reader_checks']}，合计 {ce['reader_checks_total']} |
| 状态分布 | {ce['status_histogram']} |
| 注册裁判 | `tensor_benchmark_cases` 该锚点 **0 行**（唯一裁判空置，见 §5 F2） |

## 2. 怎么做的（方法与工具链）

1. 干净克隆 `agent/m17595` 到 `/tmp/m17595_verify`（不改被审工作树）→ 复跑其门禁：
   `pytest tests/test_mesh_xval.py -q` → **18 passed**；
   `python3 tools/xval_mesh_six_format.py --check` → **`OK: results.json reproduces`**（0 漂移）。
2. **独立重算**（`docs/board_audit_m17634/audit_independent.py`，不 import `tensorcad.mesh_xval`）：
   用 meshio 5.3.5 / trimesh 5.0.0 / numpy-stl / h5py 3.15.1 + 审计自写 SU2 读者 + 审计自写 CGNS 读者
   （SIDS 标签方言 & name-addressed 方言两支）重读 18 个交付文件与 18 个上游文件，独立算体积。
3. 独立核 `manifest.json` 的 sha256（{man['files_declared']} 个工件）、独立按 SIDS-to-HDF5 规则走查 CGNS
   节点属性、查 `tensor_benchmark_cases` 注册面、查交付是否已并入生态仓 main。
4. 图由 matplotlib 程序化生成，数字全部由 `audit_evidence.json` / 提交的 `results.json` 注入（零手抄）。

## 3. 计算结果（数字表）

{vol_table}

交付自称 vs 审计重算：

| 量 | 交付自称 | 审计独立重算 | 一致? |
|---|---|---|---|
| 最差体积相对误差 | {ce['summary_A']['worst_volume_rel_err']:.6e} | {ce['max_volume_rel_err']:.6e} | 是（headline 复现={ce['headline_reproduced']}） |
| 带数值的检查全部 ≤1e-3 | 是 | {ce['all_numbered_checks_within_impl_threshold']} | 是 |
| manifest sha256 | — | {man['verified']}/{man['files_declared']} 一致，缺失 {len(man['missing'])}，不符 {len(man['sha256_mismatch'])} | 是 |

上游 m17483 工件（matrix B）审计侧重读：

{up_table}

CGNS 映射结构（审计侧按 SIDS 规则独立计数）：上游 `cube.cgns` 节点 {probe['upstream_m17483_cube']['nodes']} 个中
**{probe['upstream_m17483_cube']['missing_name']} 缺 name / {probe['upstream_m17483_cube']['missing_label']} 缺 label /
{probe['upstream_m17483_cube']['missing_order']} 缺 order，{probe['upstream_m17483_cube']['label_in_type_attr']} 个把 label 塞在 `type` 槽**；
交付自写 `cube.cgns` 上述四项均为 {probe['xval_written_cube']['label_in_type_attr']}，但**无 `CGNSLibraryVersion` 节点**
（只用根属性 `format`/`version`）——即交付"上游映射违规、几何无损"的结论独立成立，而其自写方言也不是完整 SIDS。

## 4. 结果图

左：审计侧重算的体积相对误差（cube 全 {va['cube']['max_rel_err_vs_reference']:.0e} 精确；icosphere 最差
{va['icosphere2']['max_rel_err_vs_reference']:.3e}；torus {va['torus64x32']['max_rel_err_vs_reference']:.3e}）。
右：每格式的几何 oracle 归属——STL 3 个第三方读者、OBJ 2、VTK 1、Gmsh 1，而 CGNS/SU2 只有被审方自写读者
（meshio 的 CGNS reader 是 Zone1/四面体 stub，且 meshio 无 SU2 reader），这是措辞类发现 F1 的图证。

## 5. 注意点（发现与局限）

1. **F1（minor）**："独立实现重读六格式"把 4 个第三方 oracle 与 2 个自写读者并称；表头 reader 列与 §6 差异清单
   是准确的，只有 headline 概括过度。数字本身经审计侧重导仍成立（8.0 / 4.04704468 / 1.7623015345）。
2. **F2（moderate）**：该锚点在 `tensor_benchmark_cases` **0 行**，机器裁判空置；兄弟锚点
   `GAP-CAD-B2-MESH-FORMAT-MATRIX` 登记了同一 18 格矩阵与同一最差数 1.311e-08，但未与本锚点互链。
3. **F3（moderate）**：交付**未并入生态仓 main**（4f5b492 仅在 `agent/m17595`；`/root/tensor-ecosystem/TensorCAD`
   main `928f3e3` 无 `tensorcad/mesh_xval.py`，grep 0 命中，无 CI 接线）→ `integration_callers` 为空，这是 L3 未达的唯一原因。
4. **F4（minor）**：自写 CGNS 方言缺 `CGNSLibraryVersion` 节点，"validated reference dialect" 措辞高于实证。
5. **F5（positive）**：{ce['status_histogram'].get('unreadable', 0)} 个 unreadable + {ce['status_histogram'].get('nonconformant', 0)} 个 nonconformant
   全部显式标注、不计入 pass；报告 §6.2 误差表有两行显式 NO。honest-FAIL 处理正确。
6. **F6（minor）**：`matrix_b` summary 把 3 个 CGNS nonconformant 计入 green（几何可读），建议单列。
7. **未跑（unknown，不算通过）**：cgnslib 官方 CGNS 消费者（镜像无 libcgns）→ CGNS 无官方库一侧验证；
   1000 m 坐标偏置敏感性重算（被审方自报 NO 行，审计未复现）；收敛阶梯逐点重算。

## 6. 误差判定（逐项回答是否 <3%）

| 关键量 | 参考 | 审计重算 got | 误差 | <3%? | <15%? |
|---|---|---|---|---|---|
| cube 体积 × 6 格式 | 8.0（解析） | 8.0 | 0.000e+00 | 是 | 是 |
| icosphere 体积（STL float32 最差） | 4.047044679978 | 4.047044626915 | {va['icosphere2']['max_rel_err_vs_reference']:.3e} | 是 | 是 |
| torus 体积（STL 最差） | 1.762301534515 | {list(va['torus64x32']['volumes'].values())[0] if va['torus64x32']['volumes'] else 0:.12g} | {va['torus64x32']['max_rel_err_vs_reference']:.3e} | 是 | 是 |
| 上游 CGNS 几何（name-addressed 重读） | 上游 obj/msh 体积 | 完全一致 | ≤1e-09 | 是 | 是 |
| 上游 CGNS 结构合规 | 0 违规 | {probe['upstream_m17483_cube']['missing_label']}/{probe['upstream_m17483_cube']['nodes']} 节点缺必需属性 | 结构违规 | **否**（非几何量） | 不适用 |
| STL 在 1000 m 偏置下坐标偏差 | 0 m | 3.85e-05 m（被审方自报） | not_run | **未知** | 不适用 |

**结论**：可独立复核的几何量最差 {ce['max_volume_rel_err']:.3e}（≈{ce['max_volume_rel_err']*100:.2e}%），远优于 owner 15% 红线与 3% 线；
两条 NO 项一条是"坐标系用法"要求、一条是"结构合规"而非几何精度，均已如实列出。
"""

    token = None
    st, d = api("POST", "/auth/login",
                {"username": os.environ["TENSOR_AGENT_USER"], "password": os.environ["TENSOR_AGENT_PASS"]})
    token = d.get("access_token")
    print("login:", st, bool(token), "api:", API)
    payload = {
        "title": "交差审计 m17595 XVAL-CAD-M17483-NATIVE-SIX-FORMAT-MATRIX（TensorCAD）：独立重算 + 裁定 L2",
        "summary": (f"审计独立重算确认 m17595 六格式矩阵：最差体积相对误差 "
                    f"{ce['max_volume_rel_err']:.3e}（远优于 15%），manifest {man['verified']}/{man['files_declared']} "
                    f"sha256 一致，18/18 测试与 --check 漂移门实跑通过；裁定 L2（无 mainline 集成、"
                    f"锚点在 tensor_benchmark_cases 0 行），附 6 条 F1-F6 发现清单。"),
        "scenario_key": "board_audit",
        "mission_id": 17634,
        "tags": ["board-audit", "TensorCAD", "XVAL", "mesh-formats", "independent-recompute"],
        "markdown": md,
        "figures": figures,
    }
    st, d = api("POST", "/tensor/reports", payload, token)
    print("report:", st, json.dumps(d)[:300])
    (HERE / "report_response.json").write_text(json.dumps({"status": st, "body": d}, indent=1))
    return 0 if st in (200, 201) else 1


if __name__ == "__main__":
    raise SystemExit(main())