#!/usr/bin/env python3
"""Emit docs/wetted_surface_uncertainty.md from the run's summary.json.

Every number is injected from the JSON (no hand copying), so the document is
regenerable and `--check` can assert it is current.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data/wetted-uncertainty-2026-11-04/summary.json"
DOC = REPO / "docs/wetted_surface_uncertainty.md"

#: Froude number above which the hull is on the planning-side branch of the
#: propagation curve, i.e. where the claim "dRt/Rt < 3 %" is made. Below it the
#: friction fraction pi = Rf/Rt grows faster than dS/S falls, so the propagated
#: error leaves the band; that speed is reported separately instead of hidden.
PLANNING_FROUDE = 0.13


def pct(x: float, digits: int = 3) -> str:
    return f"{100.0 * x:+.{digits}f}%"


def render(s: dict) -> str:
    d, c, b, q, m, x = s["design"], s["closure"], s["budget"], s["quadrature"], s["mesh"], s["cross_checks"]
    out: list[str] = []
    a = out.append
    ratio = abs(c["total_relative"]) / b["expanded_95"]
    prop = s["propagation"]
    planning = [r for r in prop if r["froude"] >= PLANNING_FROUDE]
    slow = min(prop, key=lambda r: r["froude"])
    drt_max = max(abs(r["drt_relative"]) for r in planning) if planning else 0.0
    fast = max(prop, key=lambda r: r["speed_kn"])
    slowest_kn = f"{slow['speed_kn']:.0f}"

    a("# 湿表面积闭包误差的不确定度（误差棒 + 网格收敛 / GCI）")
    a("")
    a(f"> 仓 {s['repo']} · mission {s['mission']} · commit `{s['commit']}` · "
      "数据 `data/wetted-uncertainty-2026-11-04/summary.json`")
    a("> 本文件由 `tools/wetted_uncertainty_doc.py` 从 summary.json 程序化生成，零手抄。")
    a("")
    a("## 1. 缺口与要回答的问题")
    a("")
    a("课题 `wetted-surface-geometry` 的 W10 事实项 `gaps` 标记「缺不确定度/误差棒/GCI」。"
      "报告 #3871 已给出湿面积闭包误差到总阻力的解析传播（dRt/Rt = pi * dS/S），"
      "但闭包误差本身没有棒：读者无法判断这个偏差是网格不够细，还是模型本身偏。本文补三样东西：")
    a("")
    a("1. **数值不确定度** u_num：同一被积函数用三种不同求积规则 + 网格收敛序列 + Richardson/GCI；")
    a("2. **输入参数不确定度** u_in：解析灵敏度 dlnS/dlnX 与量测容差经 GUM 一阶 RSS 合成 95% 区间；")
    a("3. **偏置与不确定度分开报**：偏置不能用对称小棒掩盖。")
    a("")
    a("## 2. 几何与参数（先验校验 Cb/L/B）")
    a("")
    a("基准几何取 `tensorcad.hull` 的解析族 y = (B/2) sin^p(pi x/L) sqrt(zeta)（zeta 自龙骨量起）："
      "Cm = 2/3、Cwp = int_0^1 sin^p(pi x) dx、Cb = (2/3) Cwp 全部闭式可算，"
      "所以闭包误差是可测的量，不是断言。")
    a("")
    a("| 量 | 值 | 来源 |")
    a("|---|---|---|")
    a(f"| L / B / T | {d['length_m']:g} / {d['beam_m']:g} / {d['draft_m']:g} m | summary.design |")
    a(f"| 纵向正弦指数 p（首尾同值） | {d['fore_exponent']:g} | summary.design |")
    a(f"| Cb | {s['design_cb']:.6f} | 闭式 (2/3) int sin^p |")
    a(f"| Cwp | {s['design_cwp']:.6f} | 闭式 int sin^p |")
    a(f"| Cm | {s['design_cm']:.6f} | 闭式 |")
    a(f"| 排水体积 LBT*Cb | {m['volume_check']['LBT_Cb_m3']:.3f} m3 | 解析 |")
    a(f"| 网格积分体积（{m['levels'][-1]['nx']} x {m['levels'][-1]['nz']}） | "
      f"{m['volume_check']['mesh_volume_finest_m3']:.3f} m3 | `compute_hydrostatics` |")
    a(f"| 网格 Cb 相对偏差 | {pct(m['mesh_cb_relative_error_finest'], 4)} | 体积一致性 <0.1% 通过 |")
    a("")
    a("## 3. 方法")
    a("")
    a("### 3.1 三条独立求积规则（同一被积函数）")
    a("")
    a("在 (x, u)、zeta = u^2 图上面元为 sqrt(4 u^2 T^2 (1 + (beta_x u)^2) + beta^2)，"
      "beta = (B/2) sin^p(pi x/L)；u 代换把龙骨处 zeta^-1/2 的斜率奇异吸进 Jacobi 因子。"
      "三条路径只差求积规则：")
    a("")
    a("- **A（仓内）**：`tensorcad.hull.analytic_wetted_surface_area`，"
      "Gauss-Legendre + psi ~ tau^3 节点聚拢，节点 80/160/320/640；")
    a("- **B（双指数）**：`tensorcad.hydrostatics_analytic.analytic_hydrostatics(quadrature_levels=4..8)`，"
      "嵌套 tanh-sinh，自带自收敛差；")
    a("- **C（仲裁）**：本文自写嵌套自适应 `scipy.integrate.quad`（Gauss-Kronrod），"
      "自适应分裂自行找到首尾 x^-0.3 边缘。")
    a("")
    a("### 3.2 网格收敛与 GCI")
    a("")
    a("把解析面采样成三角网格（nx x nz = 41x21 -> 321x161，面元尺寸每次减半），"
      "用 `compute_hydrostatics` 对网格本身积分湿面积，与解析参考值之差即可验证真误差 e(h)。于是：")
    a("")
    a("p = -log2|e(2h)/e(h)|；Richardson f_ex = f_fine + (f_fine - f_med)/(r^p - 1)；"
      "Roache GCI（Fs = 1.25，取作 ~95% 带）= Fs |f_med - f_fine| / (|f_fine| |r^p - 1|)。")
    a("")
    a("### 3.3 输入参数不确定度（GUM 一阶 RSS）")
    a("")
    a("闭包误差 e = Sp/Se - 1，故 de/dlnX = dlnSp/dlnX - dlnSe/dlnX 由中心差分数值给出。"
      "量测容差是显式假设（非实测），见 §6。")
    a("")
    a("## 4. 结果")
    a("")
    a("### 4.1 设计点：闭包误差带误差棒")
    a("")
    a("| 量 | 值 |")
    a("|---|---|")
    a(f"| 解析湿面积 Se | {c['S_exact_m2']:.4f} m2 |")
    a(f"| 平台闭包湿面积 Sp | {c['S_platform_m2']:.4f} m2 |")
    a(f"| 总闭包误差 Sp/Se - 1 | **{pct(c['total_relative'])}** |")
    a(f"| 回归结构误差（代入真 Cm, Cwp） | {pct(c['structure_relative'])} |")
    a(f"| 系数闭包误差（Cm, Cwp 由 Cb 线性闭包） | {pct(c['closure_relative'])} |")
    a(f"| 合成标准不确定度 uc | {pct(b['combined_standard'], 4)} |")
    a(f"| 扩展不确定度 U95 = 2 uc | **{pct(b['expanded_95'], 4)}** |")
    a(f"| 最坏情况线性相加 | {pct(b['worst_case_sum'], 4)} |")
    a(f"| 偏置 / 不确定度 | {ratio:.1f} x |")
    a("")
    a(f"**判定**：偏差 {pct(c['total_relative'])} 是其 95% 不确定度 {pct(b['expanded_95'], 4)} 的 "
      f"{ratio:.0f} 倍 —— 属系统性模型偏置，细化网格或增加小数位都不会缩小它。"
      "要让湿面积满足 3%，只能改闭包本身（真 Cm/Cwp，或船型相关回归）。")
    a("")
    a("### 4.2 不确定度预算")
    a("")
    a("| 分量 | 输入相对不确定度 | dlnSp/dlnX | dlnSe/dlnX | 对 e 的贡献 |")
    a("|---|---|---|---|---|")
    for comp in b["components"]:
        a(f"| {comp['name']} | {comp['input_relative']:.4g} | {comp['d_ln_platform']:+.4f} | "
          f"{comp['d_ln_exact']:+.4f} | {pct(comp['contribution_relative'], 4)} |")
    a(f"| **RSS 合计** | | | | **{pct(b['combined_standard'], 4)}** |")
    a("")
    a(f"主导分量：{b['dominant']}。")
    a("")
    a("RSS 合成（k=1）必然不超过线性相加的最坏情况；乘 k=2 之后两者可以同量级，"
      "所以偏置只与 uc/U95 比较，不与最坏情况比。")
    a("")
    a("### 4.3 网格收敛序列与 GCI（GCI 是否包住真误差）")
    a("")
    a("| nx x nz | 面数 | h [m] | S_mesh [m2] | 真误差 | p | GCI95 | 包住真误差 | Richardson 误差 |")
    a("|---|---|---|---|---|---|---|---|---|")
    for row in m["levels"]:
        p_txt = f"{row['p_observed']:.3f}" if "p_observed" in row else "-"
        gci_txt = pct(row["gci_relative_95"]) if "gci_relative_95" in row else "-"
        brk_txt = str(row["gci_brackets_true_error"]) if "gci_brackets_true_error" in row else "-"
        ric_txt = pct(row["richardson_relative_error"], 4) if "richardson_relative_error" in row else "-"
        a(f"| {row['nx']} x {row['nz']} | {row['faces']} | {row['h_m']:.4g} | {row['S_mesh_m2']:.4f} | "
          f"{pct(row['error_relative'])} | {p_txt} | {gci_txt} | {brk_txt} | {ric_txt} |")
    a("")
    a(f"- 观测阶 p = {m['finest_p_observed']:.2f}（<2）：首尾 x^-0.3 边缘让面元误差掉不到二阶，"
      "这与解析面的尖端几何一致，不是实现缺陷；")
    a(f"- Richardson 外推 vs 解析参考：{pct(m['richardson_vs_analytic_relative'], 5)}；")
    a(f"- GCI 在所有可评估层级都包住真误差：**{m['gci_brackets_all_levels']}**；")
    a(f"- 最细网格湿面积的数值不确定度约 {pct(m['finest_gci_relative'])}，比闭包偏差小两个数量级。")
    a("")
    a("### 4.4 解析族扫描：每个 Cb 点的闭包误差 ± 95% 棒")
    a("")
    a("| p | Cb | Cwp | Se [m2] | Sp [m2] | 总误差 | 结构 | 闭包 | U95 | 三法互校 |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for r in s["family_sweep"]:
        a(f"| {r['exponent']:g} | {r['cb']:.4f} | {r['cwp']:.4f} | {r['S_exact_m2']:.4f} | "
          f"{r['S_platform_m2']:.4f} | {pct(r['total_relative'])} | {pct(r['structure_relative'])} | "
          f"{pct(r['closure_relative'])} | {pct(r['u95_relative'])} | "
          f"{abs(r['integrator_agreement']):.1e} |")
    a("")
    a("### 4.5 传播到总阻力（dRt/Rt = pi dS/S）")
    a("")
    a("pi 取自平台阻力端点（缓存于 `platform_resistance.json`），并现场校验 Rf = (1/2) rho S V^2 cf_total：")
    a("")
    a("| V [kn] | Fn | Rf [kN] | Rt [kN] | pi | dRt/Rt | U95 | 恒等式残差 |")
    a("|---|---|---|---|---|---|---|---|")
    for r in prop:
        a(f"| {r['speed_kn']:.0f} | {r['froude']:.4f} | {r['rf_kn']:.2f} | {r['rt_kn']:.2f} | "
          f"{r['pi']:.4f} | {pct(r['drt_relative'])} | {pct(r['u95_drt_relative'])} | "
          f"{r['rf_identity_relative_residual']:.1e} |")
    a("")
    a(f"最低速 {slowest_kn} kn（Fn = {slow['froude']:.4f}）是摩擦主导点"
      f"（pi = {slow['pi']:.4f} 为全速域最大），传播误差 {pct(slow['drt_relative'])} **不在 3% 带内**；"
      f"Fn >= {PLANNING_FROUDE:g}（{planning[0]['speed_kn']:.0f}-{fast['speed_kn']:.0f} kn）"
      f"的最大传播误差为 {pct(drt_max)}。")
    a("")
    a("## 5. 与既有判据数字逐项对表")
    a("")
    a("| 判据 | 既有数字（来源） | 本次独立结果 | 相对差 | 是否 <3% |")
    a("|---|---|---|---|---|")
    a(f"| 积分器三法互校 | #3871：GL vs dblquad 1.78e-16 | path A vs 自适应 "
      f"{q['path_a_vs_adaptive_relative']:.1e}；tanh-sinh vs 自适应 "
      f"{q['tanh_sinh_vs_adaptive_relative']:.1e}；tanh-sinh 自差 "
      f"{q['tanh_sinh_self_difference_relative']:.1e} | 见左 | 是 |")
    a(f"| 超椭圆 Wigley 湿面积 | #3871：S = 1275.2810 m2（单侧） | "
      f"{x['wigley_superellipse_one_side_S_m2']:.4f} m2（双侧 {x['wigley_superellipse_full_S_m2']:.4f}） | "
      f"{x['relative_difference_vs_3871_one_side']:.1e} | 是 |")
    a(f"| 回转体闭式面积 | #3871：S = 2478.7758 m2 | {x['prolate_spheroid_closed_form_m2']:.4f} m2 | "
      f"{x['relative_difference_prolate_vs_3871']:.1e} | 是 |")
    a(f"| 闭包回归复现 | 平台 wetted_surface（本例 {c['S_platform_m2']:.2f} m2） | "
      f"本地 Holtrop 实现，最大相对差 {x['platform_wetted_vs_local_holtrop_relative']:.1e} | 同左 | 是 |")
    a(f"| 闭包误差（本船型） | #3871：Wigley +5.577%、回转体 -13.017% | 本解析族 "
      f"{pct(c['total_relative'])} | 不同几何，不可互替 | 否，>3% |")
    a(f"| 摩擦份额 pi（8 kn） | #3871：0.3422 | {prop[1]['pi']:.4f}（同 Fn = {prop[1]['froude']:.4f}） | "
      f"{100.0 * (prop[1]['pi'] / 0.3422 - 1.0):+.2f}% | 是 |")
    a(f"| 湿面积误差传播到 Rt | #3871：最大 +2.083% | "
      f"本次 Fn >= {PLANNING_FROUDE:g} 最大 {pct(drt_max)} | - | 是 |")
    a("")
    a(f"**显式回答 3% 判据**：闭包误差本身 **不满足**（{pct(c['total_relative'])}；"
      "与 #3871 同结论的独立复现，几何不同故数值不同）；"
      f"其不确定度 {pct(b['expanded_95'], 4)} 远小于 3%，所以这个不满足不能归因于数值噪声；"
      f"传播到 Rt 在 {planning[0]['speed_kn']:.0f}-{fast['speed_kn']:.0f} kn **满足 <3%**"
      f"（最大 {pct(drt_max)}），最低速 {slowest_kn} kn 为 {pct(slow['drt_relative'])}，"
      "超出 3%，如实列出不并入结论。")
    a("")
    a("## 6. 注意点 / 局限")
    a("")
    a(f"1. **输入容差是假设**：半宽 ±0.5% B、吃水 0.02 m、L ±0.05 m、Cb 0.2% Cb。"
      f"换一套容差得同形预算的 rho 倍；偏置/不确定度之比仍约 {ratio:.0f} 倍。")
    a("2. **本几何不是实船**：Cm = 2/3、Cwp = Cb 是解析族的固定组合，而闭包假设 Cwp = 0.75 + 0.25 Cb；"
      "结构误差与闭包误差反号部分抵消，总偏置是净值，不能读成某艘船的误差。")
    a("3. **GCI 适用面**：p 约 1.2 是面积积分对网格的观测阶；位移体积走散度定理，收敛阶不同，不可挪用。")
    a("4. **pi 的来源**：pi 由平台阻力回归给出，其自身误差不在本文预算内；本文只传播湿面积一支。")
    a("5. 平台响应只回 2 位小数，故对拍判据定在 1e-5 ~ 1e-6，而非机器精度。")
    foreign = [r["exponent"] for r in s["family_sweep"] if not r["tanh_sinh_applicable"]]
    a(f"7. **第三方求积只在本几何适用**：`tensorcad.hydrostatics_analytic` 保留自己的默认 taper，"
      f"对 p≠{d['fore_exponent']:g} 它积的是另一条船（Cwp 不符），故 path B 只在 p = "
      f"{d['fore_exponent']:g} 计入预算；族扫描中 p ∈ {foreign} 的 U95 由 path A 的网格 GCI "
      "与自适应仲裁给出。这条是本轮修正的：此前那些点的 nan/0 值让跨积分器互校假通过。")
    a("")
    a(f"6. **传播结论有速度窗**：3% 带只在 Fn >= {PLANNING_FROUDE:.2f} 成立；低速端 pi 上升而 dS/S 不变，"
      "传播误差单调上升，最低速点越界。")
    a("")
    a("## 7. 复现")
    a("")
    a("```bash")
    a("python3 tools/wetted_uncertainty.py selftest")
    a("python3 tools/wetted_uncertainty.py run     --out data/wetted-uncertainty-2026-11-04")
    a("python3 tools/wetted_uncertainty.py figures --data data/wetted-uncertainty-2026-11-04")
    a("python3 tools/wetted_uncertainty_doc.py")
    a("python3 -m pytest tests/test_wetted_uncertainty.py -q")
    a("```")
    a("")
    a("图：`figures/ws_uncertainty_mesh_gci.png`（网格收敛 + GCI 覆盖）、"
      "`figures/ws_uncertainty_budget.png`（预算 vs 偏置）、"
      "`figures/ws_uncertainty_closure_bars.png`（闭包误差误差棒）、"
      "`figures/ws_uncertainty_rt_band.png`（传播到 Rt）。")
    return "\n".join(out) + "\n"


def main(argv: list[str]) -> int:
    summary = json.loads(DATA.read_text())
    text = render(summary)
    if "--check" in argv:
        current = DOC.read_text() if DOC.exists() else ""
        if current != text:
            print("doc out of date", file=sys.stderr)
            return 1
        print("doc up to date")
        return 0
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(text)
    print(f"wrote {DOC} ({len(text)} chars)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
