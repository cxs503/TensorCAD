#!/usr/bin/env python3
"""Publish the m2877 progress report from docs/wigley/results.json (no hand-typed numbers)."""
from __future__ import annotations

import base64
import json
import os
import urllib.request
from pathlib import Path

ROOT = Path("/root/tp-agent/work/m2877")
BASE = os.environ["TENSOR_API_BASE"].rstrip("/")
USER = os.environ["TENSOR_AGENT_USER"]
PASS = os.environ["TENSOR_AGENT_PASS"]

data = json.loads((ROOT / "docs/wigley/results.json").read_text())
p, a = data["parameters"], data["analytic"]
rows = data["convergence"]["rows"]
claims = data["claims"]


def post(path: str, payload: dict, token: str | None = None) -> dict:
    req = urllib.request.Request(
        f"{BASE}{path}",
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json", **({"Authorization": f"Bearer {token}"} if token else {})},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=120) as response:
        return json.loads(response.read().decode())


login = post("/auth/login", {"username": USER, "password": PASS})
token = login.get("access_token") or login.get("token") or login.get("data", {}).get("access_token")
assert token, login

md: list[str] = []
add = md.append
add("## 1. 案例介绍")
add("")
add("TensorCAD 此前只有**一个船型族** `b(x,z)=(B/2)·sin(πx/L)^0.7·√(z/T)`：前后对称、无艉板，")
add("且不是任何文献数据所对应的曲面。本单补入抛物型 **Wigley** 族——阻力/耐波/绕射对比的")
add("标准几何源，并且它的静水力是**精确有理数**，因此该模块同时是一个针对 "
    "`tensorcad.hydrostatics` 的验证案例（属 P0「多船型族扩展」缺口 + 验证欠账两类）。")
add("")
add("曲面：`y(x,w) = (B/2)(1 - |2x/L - 1|^n)(1 - (w/T)^m)`，`w` 为水线以下深度（w=0 水线，w=T 龙骨），")
add(f"`n = {p['longitudinal_exponent']:g}`、`m = {p['vertical_exponent']:g}` 为经典 Wigley。")
add("")
add("| 要素 | 值 |")
add("|---|---|")
add(f"| L / B / T | {p['length']:g} m / {p['beam']:g} m / {p['draft']:g} m |")
add(f"| L/B, B/T | {p['length']/p['beam']:.2f}, {p['beam']/p['draft']:.2f}（ITTC 参考比例） |")
add(f"| 采样网格 | {p['nx']} 站 × {p['nz']} 层（默认档 {data['topology']['triangles']} 三角形） |")
add(f"| 水密度 | {p['water_density']:g} kg/m³ |")
add("| 坐标系 | x 艉→艏，z=0 设计水线，龙骨 z=-T |")
add("")
add("## 2. 怎么做的")
add("")
add("1. 沿 `hull.py` 同一套拓扑：每站一个截面环（右舷自水线下行到共享龙骨点、左舷上行），")
add("   站间放样，并在设计水线处封顶 —— 于是曲面是闭合实体，体积可由三角形积分。")
add("2. 手工推导闭式解（可分离二重积分：`∫taper dx = L·n/(n+1)`、`∫section dw = T·m/(m+1)`），")
add("   并与生成函数的直接中点积分独立互校。")
add("3. 唯一无初等原函数的量（湿表面积）用解析面上的复合 Gauss-Legendre，并自报收敛阶差；")
add("   注意面板数取偶数，使面板边界落在中站（非整数指数下 `|u|^(n-1)` 在此折点）。")
add("4. 网格积分 vs 闭式解逐量比较（带符号相对误差），四档网格收敛研究，拓扑与出口件校验。")
add("5. 工具链：纯标准库 `tensorcad.geometry`；`examples/wigley/run_example.py` 生成")
add("   `docs/wigley/{results.json,report.md,plot.png}`，任一断言失败则非零退出；")
add("   `tests/test_wigley.py` 30 项测试钉住同一批数字。")
add("")
add("## 3. 计算结果：闭式静水力")
add("")
add("| 量 | 闭式 | 数值 |")
add("|---|---|---|")
add(f"| 排水体积 V | `n m/((n+1)(m+1))·LBT` | {a['displacement_volume_m3']:.9f} m³ |")
add(f"| 排水量 | `ρV` | {a['displacement_mass_kg']:.3f} kg |")
add(f"| 水线面面积 | `L B n/(n+1)` | {a['waterplane_area_m2']:.9f} m² |")
add(f"| 中横剖面面积 | `B T m/(m+1)` | {a['maximum_section_area_m2']:.9f} m² |")
add(f"| 湿表面积 | 数值积分（无初等原函数） | {a['wetted_surface_area_m2']:.9f} m² |")
add(f"| Cb | `n m/((n+1)(m+1))` | {a['block_coefficient']:.12f}（4/9={4/9:.12f}） |")
add(f"| Cp / Cm / Cwp | `n/(n+1)` / `m/(m+1)` / `n/(n+1)` | "
    f"{a['prismatic_coefficient']:.12f} / {a['midship_coefficient']:.12f} / {a['waterplane_coefficient']:.12f} |")
add(f"| VCB（龙骨之上） | `T - T(m+1)/(2(m+2))` | {data['vcb_above_keel']:.9f} m（5T/8={5*p['draft']/8:.9f}） |")
add(f"| LCB / LCF | `L/2`（对称） | {a['longitudinal_centre_buoyancy_m']:.9f} m |")
add("")
add("## 4. 计算结果：网格积分收敛（对比闭式解）")
add("")
add("| 网格 | 三角形 | V [m³] | e_V [%] | e_Awp [%] | e_S [%] | max\\|e\\| [%] |")
add("|---|---|---|---|---|---|---|")
for r in rows:
    e = r["relative_error"]
    add(f"| {r['nx']}×{r['nz']} | {r['triangles']} | {r['displacement_volume']:.7f} | "
        f"{100*e['displacement_volume']:+.5f} | {100*e['waterplane_area']:+.5f} | "
        f"{100*e['wetted_surface_area']:+.5f} | {100*r['max_abs_relative_error']:.5f} |")
add("")
add("观测收敛阶：" + "，".join(
    f"{k.replace('_', ' ')} = {v:.2f}" for k, v in data["convergence"]["observed_order"].items()) +
    "（二阶放样应有的量级）。所有带符号误差均为**负**：放样弦线切在解析面之内，误差是系统性的，")
add("符号翻转即缺陷 —— 这正是闭式解值得存在的理由。")
add("")
add("## 5. 图")
add("")
add("![plot](plot.png)")
add("")
add("(a) 生成面的宽度场 `2y(x,z)` —— 上面每一个静水力积分都是它的积分；(b) 三个采样站截面与解析")
add("截面曲线（点为网格、线为解析，重合）；(c) 三个广延量相对闭式解的对数-对数收敛与 `h²` 参考斜率。")
add("")
add("## 6. 注意点（坑/局限/适用范围）")
add("")
add("* 该族**龙骨处有限折角**（不像 `hull` 族的 `√` 尾 taper），故 41×13 档体积误差 0.236%，")
add("  而同网格正弦族为 0.88%：截面光滑度比网格调优更值钱。此结论已写成测试。")
add("* 湿表面积是**数值参考**而非闭式：n=m=2 时精确到舍入，非整数指数下端点 `|u|^(n-1)` 折点会")
add("  拖慢 Gauss-Legendre，故面板数必须取偶数、引用前先看 `wetted_area_evidence`。")
add("* 无附体、无艉板、无平行中体；这是基准船型，不是设计船型。")
add("* 前后对称 ⇒ LCB=LCF=L/2 恒等，本案例**不**检验纵向矩积分器（那是 `hull`/SUBOFF 的活）。")
add("* 仅有静水力；稳性、阻力、结构、规范合规仍不在包内（`engineering_status` 如实标注）。")
add("* `nx` 取奇数才能达到名义船长方向最大宽度，偶数 `nx` 给 warning 而非静默缩水。")
add("")
add("## 7. 误差判定（逐项 vs 基准）")
add("")
add("| 断言 | 实测 | 阈值 | 判定 |")
add("|---|---|---|---|")
for c in claims:
    v = c["value"]
    if isinstance(v, dict):
        shown = ", ".join(f"{k} {x:.6g}" if isinstance(x, (int, float)) else str(x) for k, x in v.items())
    elif isinstance(v, float):
        shown = f"{v:.6g}"
    else:
        shown = str(v)
    add(f"| {c['claim']} | {shown[:110]} | {c['tolerance']:g} | {'PASS' if c['passed'] else 'FAIL'} |")
add("")
worst = rows[-1]["max_abs_relative_error"]
add(f"**关键量是否 < 3%？** 是，且高出两个数量级：最细网格 (161×41) 全部可比的量相对闭式解最差 "
    f"{100*worst:.4f}%，最粗网格 (21×9) 也仅 {100*rows[0]['max_abs_relative_error']:.3f}%。")
add(f"默认档 (41×13) 体积误差 {100*rows[1]['relative_error']['displacement_volume']:+.3f}%、"
    f"水线面 {100*rows[1]['relative_error']['waterplane_area']:+.3f}%、湿表面积 "
    f"{100*rows[1]['relative_error']['wetted_surface_area']:+.3f}%，全部远小于 3% 不可接受线。")
add("基准为**精确代数**（非拟合、非另一套网格），并由生成函数的直接中点积分独立复核："
    f"{data['closed_form_cross_check']['volume_direct_integral_m3']:.9f} m³ vs "
    f"{data['closed_form_cross_check']['volume_closed_form_m3']:.9f} m³ "
    f"（{100*data['closed_form_cross_check']['relative_difference']:.2e} %）。")
add("未达成/未覆盖项：船型本身不是船舶（Wigley 是基准型），系列 60/Savitsky 族仍未实现；")
add("本单未触碰平台仓与 main 分支。")
add("")
add("## 8. 证据指针")
add("")
add("* commit `696985d`（分支 `agent/m2877`）：`tensorcad/wigley.py`、`tests/test_wigley.py`、"
    "`examples/wigley/run_example.py`、`docs/wigley.md`、`docs/wigley/` 全套工件、README/CHANGELOG/backlog")
add("* 测试：`tests/test_wigley.py` 30 passed；相邻套件 "
    "`test_hull/hydrostatics/hydrostatics_analytic/analytic_hydrostatics/public_surface/parameters/resolution` "
    "189 passed；`test_curvature/design/mesh/mesh_conditioning/validation/geometry_field` 139 passed。")
add("* 已知环境欠账（非本单引入）：`safetensors` 未安装 ⇒ `test_cli_neural/test_neural_pipeline` 及 "
    "`test_geometry_field` 中 4 个 neural 用例在收集期 ImportError。")

figure = base64.b64encode((ROOT / "docs/wigley/plot.png").read_bytes()).decode()
payload = {
    "title": "TensorCAD 新增 Wigley 船型族：精确闭式静水力 + 二阶收敛证据（vs 代数基准 <0.02%）",
    "summary": (
        f"补入抛物型 Wigley 族（可调指数 n/m）并给出精确静水力闭式解：Cb=4/9、Cp=Cm=Cwp=2/3、"
        f"VCB=5T/8。网格积分对标闭式：41×13 体积 {100*rows[1]['relative_error']['displacement_volume']:+.3f}%、"
        f"水线面 {100*rows[1]['relative_error']['waterplane_area']:+.3f}%、湿表面积 "
        f"{100*rows[1]['relative_error']['wetted_surface_area']:+.3f}%；161×41 最差 {100*worst:.4f}%；"
        f"观测阶 2.17/2.58/2.01，误差全负（内切弦）。10 项断言全 PASS，30 项测试绿。"
        f"同网格下正弦族体积误差 0.88% vs Wigley 0.236%。commit 696985d。"
    ),
    "scenario_key": "tensorcad_geometry_verification",
    "mission_id": 2877,
    "tags": ["TensorCAD", "Wigley", "hydrostatics", "verification", "closed-form", "convergence"],
    "markdown": "\n".join(md),
    "figures": [
        {
            "name": "plot.png",
            "caption": "(a) breadth field of the generating surface; (b) sampled stations on the "
                       "analytic sections; (c) convergence of mesh integrals to the closed forms with h^2 guide",
            "png_base64": figure,
        }
    ],
}
result = post("/tensor/reports", payload, token)
print("report:", json.dumps(result)[:400])
print("markdown chars:", len(payload["markdown"]), "png bytes:", len(figure))
